#include "thermalvel.hpp"
#include "physconst.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_rng.h>
#include <cassert>

//Class defined in thermalvel.cpp for adding FermiDirac velocities
#define LENGTH_FERMI_DIRAC_TABLE 2000

//Amplitude of the random velocity for WDM
double WDM_V0(const double redshift, const double WDM_PartMass_in_kev, const double Omega_CDM, const double HubbleParam, const double UnitVelocity_in_cm_per_s)
{
        //Not actually sure where this equation comes from: the fiducial values are from Bode, Ostriker & Turok 2001.
        double WDM_V0 = 0.012 * (1 + redshift) * pow(Omega_CDM / 0.3, 1.0 / 3) * pow(HubbleParam / 0.65, 2.0 / 3) * pow(1.0 /WDM_PartMass_in_kev,4.0 / 3);
        WDM_V0 *= 1.0e5 / UnitVelocity_in_cm_per_s;
        /* convert from peculiar velocity to gadget's cosmological velocity */
        WDM_V0 *= sqrt(1 + redshift);
        return WDM_V0;
}

//Fermi-Dirac kernel for below
double fermi_dirac_kernel(double x, void * params)
{
  return x * x / (exp(x) + 1);
}

//Initialise the probability tables
FermiDiracVel::FermiDiracVel(const double v_amp, const double max_fd,const double min_fd) : m_vamp(v_amp), min_fd(min_fd), max_fd(std::min(MAX_FERMI_DIRAC, max_fd))
{
    assert(max_fd > min_fd);
    //Allocate random number generator
    g_rng = gsl_rng_alloc(gsl_rng_mt19937);

    double fermi_dirac_vel[LENGTH_FERMI_DIRAC_TABLE];
    double fermi_dirac_cumprob[LENGTH_FERMI_DIRAC_TABLE];

    /*These functions are so smooth that we don't need much space*/
    gsl_integration_workspace * w = gsl_integration_workspace_alloc (100);
    double abserr;
    gsl_function F;
    F.function = &fermi_dirac_kernel;
    F.params = NULL;
    for(int i = 0; i < LENGTH_FERMI_DIRAC_TABLE; i++) {
        fermi_dirac_vel[i] = min_fd+(max_fd-min_fd)* i / (LENGTH_FERMI_DIRAC_TABLE - 1.0);
        gsl_integration_qag (&F, min_fd, fermi_dirac_vel[i], 0, 1e-6,100,GSL_INTEG_GAUSS61, w,&(fermi_dirac_cumprob[i]), &abserr);
    //       printf("gsl_integration_qng in fermi_dirac_init_nu. Result %g, error: %g, intervals: %lu\n",fermi_dirac_cumprob[i], abserr,w->size);
    }
    //Save the largest cum. probability, pre-normalisation,
    //divided by the total F-D probability (which is 3 Zeta(3)/2 ~ 1.8 if MAX_FERMI_DIRAC is large enough
    double total_fd;
    gsl_integration_qag (&F, 0, MAX_FERMI_DIRAC, 0, 1e-6,100,GSL_INTEG_GAUSS61, w,&(total_fd), &abserr);
    assert(total_fd > 1.8);

    gsl_integration_workspace_free (w);

    total_frac = fermi_dirac_cumprob[LENGTH_FERMI_DIRAC_TABLE-1]/total_fd;
    //Normalise total integral to unity
    for(int i = 0; i < LENGTH_FERMI_DIRAC_TABLE; i++)
        fermi_dirac_cumprob[i] /= fermi_dirac_cumprob[LENGTH_FERMI_DIRAC_TABLE - 1];

    //Initialise the table
    fd_table = new gsl_spline_wrapper(fermi_dirac_cumprob, fermi_dirac_vel, LENGTH_FERMI_DIRAC_TABLE);
}

//Given a probability p, find a velocity v s.t. P( < v) = p.
double FermiDiracVel::get_fermi_dirac_vel(double p)
{
    return m_vamp*fd_table->eval(p);
}

//Add a randomly generated thermal speed in v_amp*(min_fd, max_fd) to a 3-velocity
std::valarray<float> FermiDiracVel::get_thermal_speeds()
{
    const double p = gsl_rng_uniform (g_rng);
    //m_vamp multiples by the dimensional factor to get a velocity again.
    const double v = get_fermi_dirac_vel(p);
    std::valarray<float> vel(0.,3);

    //Random phase
    const double phi = 2 * M_PI * gsl_rng_uniform (g_rng);
    const double theta = acos(2 * gsl_rng_uniform (g_rng) - 1);

    vel[0] = v * sin(theta) * cos(phi);
    vel[1] = v * sin(theta) * sin(phi);
    vel[2] = v * cos(theta);

    return vel;
}


