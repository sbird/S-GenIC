#ifndef THERMALVEL_H
#define THERMALVEL_H

#include <gsl/gsl_rng.h>
#include "gsl_spline_wrapper.hpp"
#include <valarray>

#define MAX_FERMI_DIRAC          17.0

/** Class to store the added thermal velocities of particles obeying Fermi Dirac statistics, like neutrinos */
class FermiDiracVel
{
    public:
        //Single parameter is the amplitude of the random velocities. All the physics is in here.
        //max_fd and min_fd give the maximum and minimum velocities to integrate over.
        //Note these values are dimensionless
        FermiDiracVel(const double v_amp, const double max_fd=MAX_FERMI_DIRAC, const double min_fd=0);
        //Add a randomly generated thermal speed in v_amp*(min_fd, max_fd) to a 3-velocity
        std::valarray<float> get_thermal_speeds();
        //Given a probability p, find a velocity v s.t. P( < v) = p.
        double get_fermi_dirac_vel(double p);
        ~FermiDiracVel()
        {
            gsl_rng_free(g_rng);
            delete fd_table;
        }
        //Total fraction of the Fermi-Dirac distribution between max_fd and min_fd
        double total_frac;
    protected:
        gsl_spline_wrapper * fd_table;
        const double m_vamp;
        const double min_fd, max_fd;
        gsl_rng * g_rng;
};

//Amplitude of the random velocity for WDM
double WDM_V0(const double redshift, const double WDM_PartMass_in_kev, const double Omega_CDM, const double HubbleParam, const double UnitVelocity_in_cm_per_s);

#endif
