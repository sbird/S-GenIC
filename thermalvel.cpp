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

//Amplitude of the random velocity for neutrinos
double NU_V0(const double redshift, const double NU_PartMass_in_ev, const double UnitVelocity_in_cm_per_s)
{
   /* I have this equation for the mean neutrino thermal velocity:

    v = 150 (1+z) \frac{1 eV}{M_\nu} km/s

    which is eq. 94 (page 41) of http://arxiv.org/abs/astro-ph/0603494

    This comes from using

    v = <p> / m \approx c \frac{3 k_B T_\nu^0}{a M_\nu}

    and T_\nu^0 = (4/11)^(1/3)T_\gamma^0
    where T_\gamma^0 = 2.726 K, the CMB temperature.
    and using k_B = 8.61734 eV/K the Boltzmann constant.

    When I plug those in and account for the slight change in neutrino temperature
    due to non-instantaneous decoupling (hep-ph/0506164) ( so T_\nu^0 = 1.00381 * (4/11)^(1/3)T_\gamma^0)
    I get 151.

    BUT: the above has <a p> = 3 T_\nu^0
    I say that this comes from this integral (making q = a p):

    < q > = int_0^\infty q n(q) dq / int_0^\infty n(q) dq

    ie, it is the value of q averaged over all states.
    For a non-relativistic Maxwell-Boltzmann distribution,
    n(q) = a M q^2 e^-(q/T) and <q> = 3 T.

    But for a Fermi-Dirac distribution (like neutrinos)
    n(q) = a M q^2 / ( e^(q/T) + 1 ) and the integral gives:

    <q> = 3 zeta(4)/zeta(3) (7/8)/(3/4)  T ~ 3.15 T

    where zeta(n) is the Riemann zeta function and I am applying the property of the Fermi-Dirac integral:

    I(p) = int_0^\infty dx x^(p-1) /(e^x+1)

     = (1-2^(1-p)) zeta(p) p!

    for integer p*/
    double NU_V0 = BOLEVK*TNU/NU_PartMass_in_ev * (1+ redshift)* (LIGHT / UnitVelocity_in_cm_per_s) * 3 * (7/8.) * (3./4)*pow(M_PI, 4)/90./1.202057;
    printf("\nNeutrino rms vel. dispersion %g (km/s)\n\n",NU_V0);
    NU_V0*=sqrt(1+redshift);
    return NU_V0;
}
#endif //NEUTRINOS

#define MAX_FERMI_DIRAC          20.0

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
    gsl_integration_workspace * w = gsl_integration_workspace_alloc (10);
    double abserr;
    gsl_function F;
    F.function = &fermi_dirac_kernel;
    F.params = NULL;
    for(int i = 0; i < LENGTH_FERMI_DIRAC_TABLE; i++) {
        fermi_dirac_vel[i] = MAX_FERMI_DIRAC * i / (LENGTH_FERMI_DIRAC_TABLE - 1.0);
        gsl_integration_qags (&F, 0, fermi_dirac_vel[i], 0, 1e-4,10,w,&(fermi_dirac_cumprob[i]), &abserr);
    //       printf("gsl_integration_qng in fermi_dirac_init_nu. Result %g, error: %g, intervals: %lu\n",fermi_dirac_cumprob[i], abserr,w->size);
    }
    gsl_integration_workspace_free (w);

    for(int i = 0; i < LENGTH_FERMI_DIRAC_TABLE; i++)
        fermi_dirac_cumprob[i] /= fermi_dirac_cumprob[LENGTH_FERMI_DIRAC_TABLE - 1];
}

//Generate a random velocity by drawing from the tabulated relationship between
//the cumulative probability function and a uniform distribution
double FermiDiracVel::get_fermi_dirac_vel(void)
{
    const double p = gsl_rng_uniform (g_rng);

    int binlow = 0;
    int binhigh = LENGTH_FERMI_DIRAC_TABLE - 2;
    int i=0;

    while(binhigh - binlow > 1)
    {
        i = (binhigh + binlow) / 2;
        if(p > fermi_dirac_cumprob[i + 1])
            binlow = i;
        else
            binhigh = i;
    }

    const double u = (p - fermi_dirac_cumprob[i]) / (fermi_dirac_cumprob[i + 1] - fermi_dirac_cumprob[i]);

    return fermi_dirac_vel[i] * (1 - u) + fermi_dirac_vel[i + 1] * u;
}

//Add a randomly generated thermal speed to a 3-velocity
void FermiDiracVel::add_thermal_speeds(float *vel)
{
    double v, phi, theta, vx, vy, vz;

    //Random amplitude of velocity
    v = m_vamp * get_fermi_dirac_vel();

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


