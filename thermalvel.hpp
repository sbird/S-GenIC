#ifndef THERMALVEL_H
#define THERMALVEL_H

#include <gsl/gsl_rng.h>

//Class defined in thermalvel.cpp for adding FermiDirac velocities
#define LENGTH_FERMI_DIRAC_TABLE 2000
#define MAX_FERMI_DIRAC          20.0

/** Class to store the added thermal velocities of particles obeying Fermi Dirac statistics, like neutrinos */
class FermiDiracVel
{
    public:
        //Single parameter is the amplitude of the random velocities. All the physics is in here.
        FermiDiracVel(double v_amp);
        void add_thermal_speeds(float *vel);
        double get_fermi_dirac_vel(double p);
        ~FermiDiracVel()
        {
            gsl_rng_free(g_rng);
        }
    protected:
        double fermi_dirac_vel[LENGTH_FERMI_DIRAC_TABLE];
        double fermi_dirac_cumprob[LENGTH_FERMI_DIRAC_TABLE];
        const double m_vamp;
        gsl_rng * g_rng;
};

//Amplitude of the random velocity for WDM
double WDM_V0(const double redshift, const double WDM_PartMass_in_kev, const double Omega_CDM, const double HubbleParam, const double UnitVelocity_in_cm_per_s);

#ifdef NEUTRINOS
//Amplitude of the random velocity for neutrinos
double NU_V0(const double redshift, const double NU_PartMass_in_ev, const double UnitVelocity_in_cm_per_s);
#endif //NEUTRINOS

#endif
