#include "cosmology.hpp"
#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>
#include "physconst.h"
#include <valarray>
#include <cassert>

/**This file contains code to compute the Hubble function at arbitrary redshift. Most of the code integrates
 * Omega_Nu for massive neutrinos, taken from the Gadget-3 patches..**/

/* We need to include radiation to get the early-time growth factor right,
 * which comes into the Zel'dovich approximation.*/
/* Omega_g = 4 \sigma_B T_{CMB}^4 8 \pi G / (3 c^3 H^2)*/
#define OMEGAG (4*STEFAN_BOLTZMANN*8*M_PI*GRAVITY/(3*LIGHTCGS*LIGHTCGS*LIGHTCGS*HUBBLE*HUBBLE*HubbleParam*HubbleParam)*pow(T_CMB0,4))
//Neutrino mass splittings
#define M21 7.53e-5 //Particle data group 2016: +- 0.18e-5 eV2
#define M32n 2.44e-3 //Particle data group: +- 0.06e-3 eV2
#define M32i 2.51e-3 //Particle data group: +- 0.06e-3 eV2
//(this is different for normal and inverted hierarchy)

/* Compute the three neutrino masses,
 * including the neutrino mass splittings
 * from oscillation experiments.
 * mnu is the total neutrino mass in eV.
 * Hierarchy is -1 for inverted (two heavy),
 * 1 for normal (two light) and 0 for degenerate.*/
std::valarray<double> NuPartMasses(double mnu, int Hierarchy)
{
    std::valarray<double> numasses(mnu/3.,3);
    //Hierarchy == 0 is 3 degenerate neutrinos
    if(Hierarchy == 0)
        return numasses;

    const double M32 = Hierarchy > 0 ? M32n : -M32i;
    //If the total mass is below that allowed by the hierarchy,
    //assume one active neutrino.
    if(mnu < ((Hierarchy > 0) ? sqrt(M32n) + sqrt(M21) : 2*sqrt(M32i) -sqrt(M21))) {
        return std::valarray<double> ({mnu, 0, 0});
    }
    //DD is the summed masses of the two closest neutrinos
    double DD1 = 4 * mnu/3. - 2/3.*sqrt(mnu*mnu + 3*M32 + 1.5*M21);
    //Last term was neglected initially. This should be very well converged.
    double DD = 4 * mnu/3. - 2/3.*sqrt(mnu*mnu + 3*M32 + 1.5*M21 + 0.75*M21*M21/DD1/DD1);
    assert(std::isfinite(DD));
    assert(fabs(DD1/DD -1) < 2e-2);
    numasses[0] = mnu - DD;
    numasses[1] = 0.5*(DD + M21/DD);
    numasses[2] = 0.5*(DD - M21/DD);
    assert(numasses[2] > 0);
    return numasses;
}

// This function converts the dimensionless units used in the integral to dimensionful units.
// Unit scaling velocity for neutrinos:
// This is an arbitrary rescaling of the unit system in the Fermi-Dirac kernel so we can integrate dimensionless quantities.
// The true thing to integrate is:
// q^2 /(e^(q c / kT) + 1 ) dq between 0 and q.
// So we choose x = (q c / kT_0) and integrate between 0 and x_0.
// The units are restored by multiplying the resulting x by kT/c for q
// To get a v we then use q = a m v/c^2
// to get:   v/c =x kT/(m a)
//NOTE: this m is the mass of a SINGLE neutrino species, not the sum of neutrinos!
double Cosmology::NU_V0(const double redshift, const double UnitVelocity_in_cm_per_s)
{
    auto masses = NuPartMasses(MNu, Hierarchy);
    //Normal or degenerate hierarchy: one heavy species, others neglected.
    double singlemass = masses[0];
    //Inverted hierarchy: two heavy species, other two neglected.
    if(Hierarchy < 0)
        singlemass = (masses[1] + masses[2])/2.;
    double NU_V0 = BOLEVK*TNU/singlemass * (1+ redshift);
    if(MNu == 0)
        NU_V0 = 1;
    return NU_V0 * sqrt(1+redshift) * (LIGHTCGS / UnitVelocity_in_cm_per_s);
}

//Hubble H(z) / H0 in units of 1/T.
double Cosmology::Hubble(double a)
{
        //Begin with curvature and lambda.
        double hubble_a = OmegaMatter(a);
        //curvature and lambda
        hubble_a += (1 - Omega - OmegaLambda) / (a * a) + OmegaLambda;
        //Add the radiation, including neutrinos
        hubble_a += OmegaR(a);
        return HUBBLE * sqrt(hubble_a);
}

/* Return the total matter density in neutrinos.
* rho_nu and friends are not externally callable*/
double Cosmology::OmegaNu(double a)
{
        double rhonu = 0;
        auto numasses = NuPartMasses(MNu, Hierarchy);
        for(int i=0; i<3;i++){
            rhonu += OmegaNu_single(a,numasses[i]);
        }
        return rhonu;
}

double Cosmology::OmegaMatter(double a)
{
    double OmegaM = Omega;
    /* Do not include neutrinos
     * in the source term of the growth function
     * in growth_ode, unless they are already not
     * free-streaming. This is not right if
     * our box overlaps with their free-streaming scale.
     * In that case the growth factor will be scale-dependent.
     * In practice the box will either be larger
     * than the horizon, which is not accurate anyway, or the neutrino
     * mass will be larger than current constraints allow,
     * so we don't worry about it too much.*/
    if(MNu > 0) {
       OmegaM -= OmegaNu(1);
       if (!NeutrinosFreeStreaming)
          OmegaM += OmegaNu(a);
    }
    return OmegaM/(a * a * a);
}

double Cosmology::OmegaR(double a)
{
    if(NoRadiation)
        return 0;
    double omr = OMEGAG/(a*a*a*a);
    //Neutrinos are included in the radiation,
    //unless they are massive and not free-streaming
    if(MNu == 0 || NeutrinosFreeStreaming)
      omr += OmegaNu(a);
    return omr;
}

bool Cosmology::SetNeutrinoFreeStream(double box, double v_th, double a)
{
    /*Neutrino free-streaming length from Lesgourgues & Pastor*/
    double nufs = 2 * M_PI * sqrt(2/3.) * v_th / Hubble(a)/a;
    NeutrinosFreeStreaming = nufs > box/4.;
    if(!NeutrinosFreeStreaming)
        fprintf(stderr,"WARNING: Neutrino free-streaming at %g boxes. May be inaccurate as scale-dependent growth is not modelled.\n",nufs/box);
    return NeutrinosFreeStreaming;
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
  return growth(aend, NULL) / growth(astart, NULL);
}

int growth_ode(double a, const double yy[], double dyda[], void * param)
{
    Cosmology * d_this = (Cosmology *) param;
    const double hub = d_this->Hubble(a)/HUBBLE;
    dyda[0] = yy[1]/pow(a,3)/hub;
    /*We want neutrinos to be non-relativistic, as only this part gravitates*/
    dyda[1] = yy[0] * 1.5 * a * d_this->OmegaMatter(a) / hub;
    return GSL_SUCCESS;
}

/** The growth function is given as a 2nd order DE in Peacock 1999, Cosmological Physics.
 * D'' + a'/a D' - 1.5 * (a'/a)^2 D = 0
 * 1/a (a D')' - 1.5 (a'/a)^2 D
 * where ' is d/d tau = a^2 H d/da
 * Define F = a^3 H dD/da
 * and we have: dF/da = 1.5 a H D
 */
double Cosmology::growth(double a, double * dDda)
{
  gsl_odeiv2_system FF;
  FF.function = &growth_ode;
  FF.jacobian = NULL;
  FF.dimension = 2;
  //Just pass the whole structure as the params pointer, as GSL won't let us make the integrand a member function
  FF.params = this;
  gsl_odeiv2_driver * drive = gsl_odeiv2_driver_alloc_standard_new(&FF,gsl_odeiv2_step_rkf45, 1e-5, 1e-8,1e-8,1,1);
   /* We start early to avoid lambda.*/
  double curtime = 1e-5;
  /* Initial velocity chosen so that D = Omegar + 3/2 Omega_m a,
   * the solution for a matter/radiation universe.*
   * Note the normalisation of D is arbitrary
   * and never seen outside this function.*/
  double yinit[2] = {OmegaR(curtime) + 1.5 * OmegaMatter(curtime)*curtime, pow(curtime,3)*Hubble(curtime)/HUBBLE * 1.5 * OmegaMatter(curtime)};
  int stat = gsl_odeiv2_driver_apply(drive, &curtime,a, yinit);
  if (stat != GSL_SUCCESS) {
      printf("gsl_odeiv in growth: %d. Result at %g is %g %g\n",stat, curtime, yinit[0], yinit[1]);
  }
  gsl_odeiv2_driver_free(drive);
  /*Store derivative of D if needed.*/
  if(dDda) {
      *dDda = yinit[1]/pow(a,3)/(Hubble(a)/HUBBLE);
  }
  return yinit[0];
}

/*
 * This is the Zeldovich approximation prefactor, 
 * f1 = d ln D1 / dlna = a / D (dD/da)
 */
double Cosmology::F_Omega(double a)
{
    double dD1da=0;
    double D1 = growth(a, &dD1da);
    return a / D1 * dD1da;
}

/* The 2LPT prefactor, f2 = d ln D2/dlna.
 * This is an approximation rather than the exact result, 
 * and not strictly valid for closed universes.
 * The approximation should be good enough since 
 * the 2lpt term is generally small.
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

