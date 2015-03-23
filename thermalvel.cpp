#include "proto.h"
#include <gsl/gsl_integration.h>

/*  Here comes the stuff to compute the thermal WDM velocity distribution */

double fermi_dirac_vel[LENGTH_FERMI_DIRAC_TABLE];
double fermi_dirac_cumprob[LENGTH_FERMI_DIRAC_TABLE];

double WDM_V0 = 0;

inline double fermi_dirac_kernel(double x, void * param)
{
  return x * x / (exp(x) + 1);
}

void fermi_dirac_init(void)
{
  int i;

  gsl_integration_workspace * w = gsl_integration_workspace_alloc (10);
  double abserr;
  gsl_function F;
  F.function = &fermi_dirac_kernel;
  F.params = NULL;
  for(i = 0; i < LENGTH_FERMI_DIRAC_TABLE; i++)
    {
      fermi_dirac_vel[i] = MAX_FERMI_DIRAC * i / (LENGTH_FERMI_DIRAC_TABLE - 1.0);
      gsl_integration_qags (&F, 0, fermi_dirac_vel[i], 0, 1e-4,10,w,&(fermi_dirac_cumprob[i]), &abserr);
//       printf("gsl_integration_qng in fermi_dirac_init. Result %g, error: %g, intervals: %lu\n",fermi_dirac_cumprob[i], abserr,w->size);
    }
  gsl_integration_workspace_free (w);

  for(i = 0; i < LENGTH_FERMI_DIRAC_TABLE; i++)
    fermi_dirac_cumprob[i] /= fermi_dirac_cumprob[LENGTH_FERMI_DIRAC_TABLE - 1];

  WDM_V0 = 0.012 * (1 + Redshift) * pow((Omega - OmegaBaryon) / 0.3, 1.0 / 3) * pow(HubbleParam / 0.65,
										    2.0 / 3) * pow(1.0 /
												   WDM_PartMass_in_kev,
												   4.0 / 3);

    printf("\nWarm dark matter rms velocity dispersion at starting redshift = %g km/sec\n\n",
	   3.59714 * WDM_V0);

  WDM_V0 *= 1.0e5 / UnitVelocity_in_cm_per_s;

  /* convert from peculiar velocity to gadget's cosmological velocity */
  WDM_V0 *= sqrt(1 + Redshift);
}



double get_fermi_dirac_vel(void)
{
  int i;
  double p, u;

  p = drand48();
  i = 0;

  while(i < LENGTH_FERMI_DIRAC_TABLE - 2)
    if(p > fermi_dirac_cumprob[i + 1])
      i++;
    else
      break;

  u = (p - fermi_dirac_cumprob[i]) / (fermi_dirac_cumprob[i + 1] - fermi_dirac_cumprob[i]);

  return fermi_dirac_vel[i] * (1 - u) + fermi_dirac_vel[i + 1] * u;
}



void add_WDM_thermal_speeds(float *vel)
{
  double v, phi, theta, vx, vy, vz;

  if(WDM_V0 == 0)
    fermi_dirac_init();

  v = WDM_V0 * get_fermi_dirac_vel();

  phi = 2 * M_PI * drand48();
  theta = acos(2 * drand48() - 1);

  vx = v * sin(theta) * cos(phi);
  vy = v * sin(theta) * sin(phi);
  vz = v * cos(theta);

  vel[0] += vx;
  vel[1] += vy;
  vel[2] += vz;
}

#ifdef NEUTRINOS

FermiDiracVelNu::FermiDiracVelNu(double redshift, double NU_PartMass_in_ev):
NU_V0(150.0 * (1.0e5 / UnitVelocity_in_cm_per_s) * (1 + redshift) * (1.0 / NU_PartMass_in_ev)*sqrt(1+redshift))
{
    /*These functions are so smooth that we don't need much space*/
    gsl_integration_workspace * w = gsl_integration_workspace_alloc (10);
    double abserr;
    gsl_function F;
    F.function = &fermi_dirac_kernel;
    F.params = NULL;
    for(int i = 0; i < LENGTH_FERMI_DIRAC_TABLE; i++) {
        fermi_dirac_vel_nu[i] = MAX_FERMI_DIRAC * i / (LENGTH_FERMI_DIRAC_TABLE - 1.0);
        gsl_integration_qags (&F, 0, fermi_dirac_vel_nu[i], 0, 1e-4,10,w,&(fermi_dirac_cumprob_nu[i]), &abserr);
    //       printf("gsl_integration_qng in fermi_dirac_init_nu. Result %g, error: %g, intervals: %lu\n",fermi_dirac_cumprob[i], abserr,w->size);
    }
    gsl_integration_workspace_free (w);

    for(int i = 0; i < LENGTH_FERMI_DIRAC_TABLE; i++)
        fermi_dirac_cumprob_nu[i] /= fermi_dirac_cumprob_nu[LENGTH_FERMI_DIRAC_TABLE - 1];

//             NU_V0 = 150.0 * (1.0e5 / UnitVelocity_in_cm_per_s) * (1 + redshift) * (1.0 / NU_PartMass_in_ev);
    /* convert from peculiar velocity to gadget's cosmological velocity */
//             NU_V0 *= sqrt(1 + redshift);

    printf("\nNeutrino rms vel. dispersion %g (km/s)\n\n",NU_V0);

}

double FermiDiracVelNu::get_fermi_dirac_vel_nu(void)
{
    int i;
    double p, u;
    int binlow = 0;
    int binhigh = LENGTH_FERMI_DIRAC_TABLE - 2;


    p = drand48();
    i = 0;

    while(binhigh - binlow > 1)
        {
        i = (binhigh + binlow) / 2;
        if(p > fermi_dirac_cumprob_nu[i + 1])
            binlow = i;
        else
            binhigh = i;
        }

    u = (p - fermi_dirac_cumprob_nu[i]) / (fermi_dirac_cumprob_nu[i + 1] - fermi_dirac_cumprob_nu[i]);

    return fermi_dirac_vel_nu[i] * (1 - u) + fermi_dirac_vel_nu[i + 1] * u;
}

void FermiDiracVelNu::add_NU_thermal_speeds(float *vel)
{
    double v, phi, theta, vx, vy, vz;

    v = NU_V0 * get_fermi_dirac_vel_nu();

    phi = 2 * M_PI * drand48();
    theta = acos(2 * drand48() - 1);

    vx = v * sin(theta) * cos(phi);
    vy = v * sin(theta) * sin(phi);
    vz = v * cos(theta);

    vel[0] += vx;
    vel[1] += vy;
    vel[2] += vz;
}

#endif

