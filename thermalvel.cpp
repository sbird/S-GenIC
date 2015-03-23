#include "proto.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_rng.h>

//Amplitude of the random velocity for WDM
double WDM_V0(double redshift, double WDM_PartMass_in_kev)
{
        double WDM_V0 = 0.012 * (1 + Redshift) * pow((Omega - OmegaBaryon) / 0.3, 1.0 / 3) * pow(HubbleParam / 0.65, 2.0 / 3) * pow(1.0 /WDM_PartMass_in_kev,4.0 / 3);
        printf("\nWarm dark matter rms velocity dispersion at starting redshift = %g km/sec\n\n",3.59714 * WDM_V0);
        WDM_V0 *= 1.0e5 / UnitVelocity_in_cm_per_s;
        /* convert from peculiar velocity to gadget's cosmological velocity */
        WDM_V0 *= sqrt(1 + Redshift);
        return WDM_V0;
}

#ifdef NEUTRINOS
//Amplitude of the random velocity for neutrinos
double NU_V0(double redshift, double NU_PartMass_in_ev)
{
    double NU_V0 = 150.0 * (1.0e5 / UnitVelocity_in_cm_per_s) * (1 + redshift) * (1.0 / NU_PartMass_in_ev);
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


