#include <math.h>
#include "allvars.h"
// #include "proto.h"
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
#define OMEGAG (4*STEFAN_BOLTZMANN*8*M_PI*GRAVITY/(3*LIGHTCGS*LIGHTCGS*LIGHTCGS*HUBBLE*HUBBLE)*pow(T_CMB0,4)/HubbleParam/HubbleParam)
/*Neutrinos are included in the radiation*/
/*For massless neutrinos, rho_nu/rho_g = 7/8 (T_nu/T_cmb)^4 *N_eff, but we absorbed N_eff into T_nu above*/
#define OMEGANU (OMEGAG*7/8.*pow(TNU/T_CMB0,4)*3)
//Neutrino mass splittings
//Neglect the smaller of the two mass splittings: 7.54e-5 is close enough to zero. #define M21
#define M32 2.43e-3 //Particle data group: +- 0.06 eV

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
double OmegaNu_single(double a,double mnu)
{
        double rhonu;
        rhonu=rho_nu(a,mnu);
        rhonu /= (3* HUBBLE* HUBBLE / (8 * M_PI * GRAVITY));
        rhonu /= HubbleParam*HubbleParam;
        return rhonu;
}

/* Return the total matter density in neutrinos.
 * rho_nu and friends are not externally callable*/
double OmegaNu(double a)
{
        double split = M32;
        //Sign of the mass splitting depends on hierarchy: inverted means the heavy state is degenerate.
        if (InvertedHierarchy)
            split *= -1;
        //This one is the degenerate state with two neutrinos in it.
        //Solve a quadratic equation in total mass to get:
        double M2 = (2*MNu - sqrt(MNu*MNu + 3*split))/3;
        double M3 = MNu - 2*M2;
        double rhonu =2*OmegaNu_single(a,M2)+OmegaNu_single(a,M3);
        return rhonu;
}

//The Hubble H(z) / H0. Thus dimensionless.
double Hubble_A(double a, double Omega, double OmegaLambda)
{
  //Begin with matter, curvature and lambda.
  double hubble_a = Omega / (a * a * a) + (1 - Omega - OmegaLambda) / (a * a) + OmegaLambda;
  //Add the radiation
  hubble_a += OMEGAG/(a*a*a*a);
  //If neutrinos are massless, add them too.
  if(MNu == 0)
        hubble_a += OMEGANU/(a*a*a*a);
  /*Otherwise add massive neutrinos, possibly slightly relativistic, to the evolution*/
  else
        hubble_a += OmegaNu(a) - OmegaNu(0)/ (a * a * a);
  return HUBBLE * UnitTime_in_s * sqrt(hubble_a);
}
