#include "cosmology.hpp"
#include <math.h>
#include <gsl/gsl_integration.h>
/**This file contains code to compute the Hubble function at arbitrary redshift. Most of the code integrates
 * Omega_Nu for massive neutrinos, taken from the Gadget-3 patches..**/

#define  LIGHTCGS 2.99792e10  //Speed of light in cm/s
#define BOLEVK 8.61734e-5        /*The Boltzmann constant in units of eV/K*/
/* Stefan-Boltzmann constant in cgs units*/
#define STEFAN_BOLTZMANN 5.670373e-5
#define  T_CMB0      2.7255	/* present-day CMB temperature, from Fixsen 2009 */
#define HBAR    6.582119e-16  /*hbar in units of eV s*/
/* Note there is a slight correction from 4/11
 * due to the neutrinos being slightly coupled at e+- annihilation.
 * See Mangano et al 2005 (hep-ph/0506164)
 *The correction is (3.046/3)^(1/4), for N_eff = 3.046 */
#define TNU     (T_CMB0*pow(4/11.,1/3.)*1.00381)              /* Neutrino + antineutrino background temperature in Kelvin */

/* We need to include radiation to get the early-time growth factor right,
 * which comes into the Zel'dovich approximation.*/
/* Omega_g = 4 \sigma_B T_{CMB}^4 8 \pi G / (3 c^3 H^2)*/
#define OMEGAG (4*STEFAN_BOLTZMANN*8*M_PI*GRAVITY/(3*LIGHTCGS*LIGHTCGS*LIGHTCGS*HUBBLE*HUBBLE*HubbleParam*HubbleParam)*pow(T_CMB0,4))
/*Neutrinos are included in the radiation*/
/*For massless neutrinos, rho_nu/rho_g = 7/8 (T_nu/T_cmb)^4 *N_eff, but we absorbed N_eff into T_nu above*/
#define OMEGANU (OMEGAG*7/8.*pow(TNU/T_CMB0,4)*3)
//Neutrino mass splittings
//Neglect the smaller of the two mass splittings: 7.54e-5 is close enough to zero. #define M21
#define M32 2.43e-3 //Particle data group: +- 0.06 eV

//Hubble H(z) / H0 in units of 1/T.
double Cosmology::Hubble(double a)
{
        //Begin with matter, curvature and lambda.
        double hubble_a = Omega / (a * a * a) + (1 - Omega - OmegaLambda) / (a * a) + OmegaLambda;
        //Add the radiation
        hubble_a += OmegaR(a);
        //If neutrinos are massless, add them too.
        if(MNu == 0)
                hubble_a += OMEGANU/(a*a*a*a);
        /*Otherwise add massive neutrinos, possibly slightly relativistic, to the evolution*/
        else
                hubble_a += OmegaNu(a) - OmegaNu(1)/ (a * a * a);
        return HUBBLE * sqrt(hubble_a);
}

/* Return the total matter density in neutrinos.
* rho_nu and friends are not externally callable*/
double Cosmology::OmegaNu(double a)
{
        double rhonu;
        if (MNu >= M32) {
            double split = M32;
            //Sign of the mass splitting depends on hierarchy: inverted means the heavy state is degenerate.
            if (InvertedHierarchy)
                split *= -1;
            //This one is the degenerate state with two neutrinos in it.
            //Solve a quadratic equation in total mass to get:
            double M2 = (2*MNu - sqrt(MNu*MNu + 3*split))/3;
            double M3 = MNu - 2*M2;
            rhonu =2*OmegaNu_single(a,M2)+OmegaNu_single(a,M3);
        }
        else {
            //For smaller masses, just assume that only one neutrino is massive.
            rhonu = OmegaNu_single(a,MNu)+2*OmegaNu_single(a,0);
        }
        return rhonu;
}

double Cosmology::OmegaR(double a)
{
    return OMEGAG/(a*a*a*a);
}
/*Note q carries units of eV/c. kT/c has units of eV/c.
 * M_nu has units of eV  Here c=1. */
double rho_nu_int(double q, void * params)
{
        double amnu = *((double *)params);
        double epsilon = sqrt(q*q+amnu*amnu);
        double f0 = 1./(exp(q/(BOLEVK*TNU))+1);
        return q*q*epsilon*f0;
}

/* Tables for rho_nu: stores precomputed values between
 * simulation start and a M_nu = 20 kT_nu*/
#define GSL_VAL 200
/*Get the conversion factor to go from (eV/c)^4 to g/cm^3
 * for a **single** neutrino species. */
double get_rho_nu_conversion()
{
        /*q has units of eV/c, so rho_nu now has units of (eV/c)^4*/
        double convert=4*M_PI*2; /* The factor of two is for antineutrinos*/
        /*rho_nu_val now has units of eV^4*/
        /*To get units of density, divide by (c*hbar)**3 in eV s and cm/s */
        const double chbar=1./(2*M_PI*LIGHTCGS*HBAR);
        convert*=(chbar*chbar*chbar);
        /*Now has units of (eV)/(cm^3)*/
        /* 1 eV = 1.60217646 × 10-12 g cm^2 s^(-2) */
        /* So 1eV/c^2 = 1.7826909604927859e-33 g*/
        /*So this is rho_nu_val in g /cm^3*/
        convert*=(1.60217646e-12/LIGHTCGS/LIGHTCGS);
        return convert;
}

/*Integrated value of rho_nu for regions between relativistic and non-relativistic*/
double integrate_rho_nu(double a,double mnu)
{
     double abserr;
     double amnu = a*mnu;
     gsl_function F;
     gsl_integration_workspace * w = gsl_integration_workspace_alloc (GSL_VAL);
     F.function = &rho_nu_int;
     F.params = &amnu;
     double rhonu;
     gsl_integration_qag (&F, 0, 500*BOLEVK*TNU,0 , 1e-9,GSL_VAL,6,w,&rhonu, &abserr);
     rhonu=rhonu/pow(a,4)*get_rho_nu_conversion();
     gsl_integration_workspace_free (w);
     return rhonu;
}

/* Value of kT/aM_nu on which to switch from the
 * analytic expansion to the numerical integration*/
#define NU_SW 50
//1.878 82(24) x 10-29 h02 g/cm3 = 1.053 94(13) x 104 h02 eV/cm3
/*Finds the physical density in neutrinos for a single neutrino species*/
double rho_nu(double a,double mnu)
{
        double rho_nu_val;
        double amnu=a*mnu;
        const double kT=BOLEVK*TNU;
        const double kTamnu2=(kT*kT/amnu/amnu);
        /*Do it analytically if we are in a regime where we can
         * The next term is 141682 (kT/amnu)^8.
         * At kT/amnu = 8, higher terms are larger and the series stops converging.
         * Don't go lower than 50 here. */
        if(NU_SW*NU_SW*kTamnu2 < 1){
            /*Heavily non-relativistic*/
            /*The constants are Riemann zetas: 3,5,7 respectively*/
            rho_nu_val=mnu*pow(kT/a,3)*(1.5*1.202056903159594+kTamnu2*45./4.*1.0369277551433704+2835./32.*kTamnu2*kTamnu2*1.0083492773819229+80325/32.*kTamnu2*kTamnu2*kTamnu2*1.0020083928260826)*get_rho_nu_conversion();
        }
        else if(amnu < 1e-6*kT){
            /*Heavily relativistic: we could be more accurate here,
             * but in practice this will only be called for massless neutrinos, so don't bother.*/
            rho_nu_val=7*pow(M_PI*kT/a,4)/120.*get_rho_nu_conversion();
        }
        else{
            rho_nu_val = integrate_rho_nu(a,mnu);
        }
        return rho_nu_val;
}

/* Return the matter density in a single neutrino species.
* Not externally callable*/
double Cosmology::OmegaNu_single(double a,double mnu)
{
        double rhonu;
        rhonu=rho_nu(a,mnu);
        rhonu /= (3* HUBBLE* HUBBLE / (8 * M_PI * GRAVITY));
        rhonu /= HubbleParam*HubbleParam;
        return rhonu;
}


double Cosmology::GrowthFactor(double astart, double aend)
{
  return growth(aend) / growth(astart);
}

double growth_int(double a, void * param)
{
    Cosmology * d_this = (Cosmology *) param;
    return 1./pow(a* (d_this->Hubble(a)),3);
}

/** The growth function is given as a 2nd order DE in Peacock 1999, Cosmological Physics.
 * Here we use an integral form: D1(a) = H integral(1/a^3 H^3 da)
 * Note there is a free proportionality constant, which Eisenstein takes to be 5/2 Omega_m
 * so that D(a) -> a as a-> 0
 * See astro-ph/9709054
 */
double Cosmology::growth(double a)
{
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (GSL_VAL);
  double hubble_a;
  double result,abserr;
  gsl_function F;
  F.function = &growth_int;
  //Just pass the whole structure as the params pointer, as GSL won't let us make the integrand a member function
  F.params = this;
  hubble_a = Hubble(a);
  gsl_integration_qag (&F, 0, a, 0, 1e-4,GSL_VAL,GSL_INTEG_GAUSS61, w,&result, &abserr);
//   printf("gsl_integration_qng in growth. Result %g, error: %g, intervals: %lu\n",result, abserr,w->size);
  gsl_integration_workspace_free (w);
  return 5./2*Omega*hubble_a * result;
}

/*Note q carries units of eV/c. kT/c has units of eV/c.
 * M_nu has units of eV  Here c=1.
 * This has 1/epsilon instead of epsilon as in rho_nu_int*/
double rho_nuprime_int(double q, void * params)
{
        double amnu = *((double *)params);
        double epsilon = sqrt(q*q+amnu*amnu);
        double f0 = 1./(exp(q/(BOLEVK*TNU))+1);
        return q*q*f0/epsilon;
}

/* Return the total matter density in neutrinos.
* rho_nu and friends are not externally callable*/
double Cosmology::OmegaNuPrimed(double a)
{
        double rhonu;
        if (MNu >= M32) {
            double split = M32;
            //Sign of the mass splitting depends on hierarchy: inverted means the heavy state is degenerate.
            if (InvertedHierarchy)
                split *= -1;
            //This one is the degenerate state with two neutrinos in it.
            //Solve a quadratic equation in total mass to get:
            double M2 = (2*MNu - sqrt(MNu*MNu + 3*split))/3;
            double M3 = MNu - 2*M2;
            rhonu =2*OmegaNuPrimed_single(a,M2)+OmegaNuPrimed_single(a,M3);
        }
        else {
            //For smaller masses, just assume that only one neutrino is massive.
            rhonu = OmegaNuPrimed_single(a,MNu);
        }
        return rhonu;
}

/** Do the integral of the derivative of OmegaNu wrt log a, which is
 *  -4 OmegaNu(a) + M_nu/a^2 q^2 f0(q) / epsilon da
 *  Neglect mass splitting for this one.
 */
double Cosmology::OmegaNuPrimed_single(double a, double mnu)
{
     double abserr;
     double amnu = a*mnu;
     gsl_function F;
     gsl_integration_workspace * w = gsl_integration_workspace_alloc (GSL_VAL);
     F.function = &rho_nuprime_int;
     F.params = &amnu;
     double rhonu;
     gsl_integration_qag (&F, 0, 500*BOLEVK*TNU,0 , 1e-9,GSL_VAL,6,w,&rhonu, &abserr);
     rhonu=-4*OmegaNu(a)+mnu*rhonu/pow(a,2)*get_rho_nu_conversion();
     gsl_integration_workspace_free (w);
     return rhonu;
}

/*
 * This is the Zeldovich approximation prefactor, f1 = d ln D1 / dlna.
 * By differentiating the growth function we get:
 * f = H' / H + 5/2 Omega /(a^2 H^2 D1)
 * and 2H H' = H0^2 * (-3 Omega/a^3 - 4 Omega_R/a^4 - 2 OmegaK/a^2)
 * The derivative of the neutrino density is solved numerically
 */
double Cosmology::F_Omega(double a)
{
  double Hprime = -3 * Omega/(a*a*a) -4 * OmegaR(a) - 2* (1-Omega-OmegaLambda)/(a*a);
  //Add the derivative of the neutrino mass to H'
  //d Omega_nu /da = - 4 Omega_nu + (derivative function)
  //With massive neutrinos OmegaNu is added twice
  if(MNu > 0) {
      Hprime += -3*OmegaNu(1)/(a*a*a);
      Hprime += OmegaNuPrimed(a);
  }
  else
      Hprime += -4*OMEGANU/(a*a*a*a);
  double HH = Hubble(a);
  double ff = HUBBLE*HUBBLE*Hprime/HH/HH/2. + 5.*Omega/(2*a*a*HH*HH*growth(a));
  return ff;
}

/* The 2LPT prefactor, f2 = d ln D2/dlna.
 * This is an approximation rather than the exact result, and not strictly valid for closed universes.
 */
double Cosmology::F2_Omega(double a)
{
  double omega_a;
  omega_a = Omega / (Omega + a * (1 - Omega - OmegaLambda) + a * a * a * OmegaLambda);
  if (1 - Omega - OmegaLambda < 0)
      return 2 * pow(omega_a, 4./7.);
  else
      return 2 * pow(omega_a, 6./11.);
}

