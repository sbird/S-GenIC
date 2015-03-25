#include "proto.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_rng.h>

//Amplitude of the random velocity for WDM
double WDM_V0(const double redshift, const double WDM_PartMass_in_kev, const double Omega_CDM, const double HubbleParam, const double UnitVelocity_in_cm_per_s)
{
        //Not actually sure where this equation comes from: the fiducial values are from Bode, Ostriker & Turok 2001.
        double WDM_V0 = 0.012 * (1 + redshift) * pow(Omega_CDM / 0.3, 1.0 / 3) * pow(HubbleParam / 0.65, 2.0 / 3) * pow(1.0 /WDM_PartMass_in_kev,4.0 / 3);
        printf("\nWarm dark matter rms velocity dispersion at starting redshift = %g km/sec\n\n",3.59714 * WDM_V0);
        WDM_V0 *= 1.0e5 / UnitVelocity_in_cm_per_s;
        /* convert from peculiar velocity to gadget's cosmological velocity */
        WDM_V0 *= sqrt(1 + redshift);
        return WDM_V0;
}

#ifdef NEUTRINOS

#define  T_CMB0      2.7255	/* present-day CMB temperature, from Fixsen 2009 */
/* Note there is a slight correction from 4/11
 * due to the neutrinos being slightly coupled at e+- annihilation.
 * See Mangano et al 2005 (hep-ph/0506164)
 *The correction is (3.046/3)^(1/4), for N_eff = 3.046 */
#define TNU     (T_CMB0*pow(4/11.,1/3.)*1.00381)              /* Neutrino + antineutrino background temperature in Kelvin */
#define BOLEVK 8.61734e-5        /*The Boltzmann constant in units of eV/K*/
#define  LIGHT           2.9979e10 /*Speed of light in cm/s*/

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
double NU_V0(const double redshift, const double NU_PartMass_in_ev, const double UnitVelocity_in_cm_per_s)
{
    double NU_V0 = BOLEVK*TNU/(NU_PartMass_in_ev/3.) * (1+ redshift)* (LIGHT / UnitVelocity_in_cm_per_s);
    NU_V0*=sqrt(1+redshift);
    return NU_V0;
}

#endif //NEUTRINOS

//Fermi-Dirac kernel for below
inline double fermi_dirac_kernel(double x, void * param)
{
  return x * x / (exp(x) + 1);
}

//Initialise the probability tables
FermiDiracVel::FermiDiracVel(double v_amp): m_vamp(v_amp)
{
    //Allocate random number generator
    g_rng = gsl_rng_alloc(gsl_rng_mt19937);
    /*These functions are so smooth that we don't need much space*/
    gsl_integration_workspace * w = gsl_integration_workspace_alloc (100);
    double abserr;
    gsl_function F;
    F.function = &fermi_dirac_kernel;
    F.params = NULL;
    for(int i = 0; i < LENGTH_FERMI_DIRAC_TABLE; i++) {
        fermi_dirac_vel[i] = MAX_FERMI_DIRAC * i / (LENGTH_FERMI_DIRAC_TABLE - 1.0);
        gsl_integration_qag (&F, 0, fermi_dirac_vel[i], 0, 1e-6,100,GSL_INTEG_GAUSS61, w,&(fermi_dirac_cumprob[i]), &abserr);
    //       printf("gsl_integration_qng in fermi_dirac_init_nu. Result %g, error: %g, intervals: %lu\n",fermi_dirac_cumprob[i], abserr,w->size);
    }
    gsl_integration_workspace_free (w);

    //Normalise total integral to unity
    for(int i = 0; i < LENGTH_FERMI_DIRAC_TABLE; i++)
        fermi_dirac_cumprob[i] /= fermi_dirac_cumprob[LENGTH_FERMI_DIRAC_TABLE - 1];
}

//Given a probability p, find a dimensionless X s.t. P(x< X) = p.
double FermiDiracVel::get_fermi_dirac_vel(double p)
{
    int binlow = 0;
    int binhigh = LENGTH_FERMI_DIRAC_TABLE - 1;
    int i=0;

    while(binhigh - binlow > 1)
    {
        i = (binhigh + binlow) / 2;
        if(p >= fermi_dirac_cumprob[i + 1])
            binlow = i;
        else
            binhigh = i;
    }

    const double u = (p - fermi_dirac_cumprob[i-1]) / (fermi_dirac_cumprob[i] - fermi_dirac_cumprob[i-1]);
    return fermi_dirac_vel[i-1] * (1 - u) + fermi_dirac_vel[i] * u;
}

//Add a randomly generated thermal speed to a 3-velocity
void FermiDiracVel::add_thermal_speeds(float *vel)
{
    double v, phi, theta, vx, vy, vz;
    const double p = gsl_rng_uniform (g_rng);
    //m_vamp multiples by the dimensional factor to get a velocity again.
    v = m_vamp * get_fermi_dirac_vel(p);

    //Random phase
    phi = 2 * M_PI * gsl_rng_uniform (g_rng);
    theta = acos(2 * gsl_rng_uniform (g_rng) - 1);

    vx = v * sin(theta) * cos(phi);
    vy = v * sin(theta) * sin(phi);
    vz = v * cos(theta);

    vel[0] += vx;
    vel[1] += vy;
    vel[2] += vz;
}


