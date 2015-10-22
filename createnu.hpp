#ifndef CREATENU_H
#define CREATENU_H

#include "power.hpp"
#include "part_data.hpp"
#include "thermalvel.hpp"
#include <hdf5.h>
#include <hdf5_hl.h>

//Class to load the neutrino power spectra
class PowerSpec_NuTabulated: public PowerSpec
{
    public:
    PowerSpec_NuTabulated(const std::string& FileWithInputSpectrum, double atime);
    virtual double power(double k, int Type);
    virtual ~PowerSpec_NuTabulated(){
        delete[] pvals;
        delete [] kvals;
        gsl_interp_free(power_interp);
        gsl_interp_accel_free(power_accel);
    }
    private:
        int NPowerTable;
        gsl_interp * power_interp;
        gsl_interp_accel * power_accel;
        double *kvals;
        double *pvals;
};

lpt_data generate_neutrino_particles(std::string SimSpecFile, part_grid& Pgrid, size_t Nmesh, double Box, int Seed, double atime);

//Copied with minor modifications from GadgetReader
uint32_t WriteBlock(const std::string& BlockName, hid_t group, int type, void *data, hid_t d_type, int ncols, uint32_t npart, uint32_t np_write, uint32_t begin);

/* Write the data to the HDF5 file.
 * Arguments:
 * SnapFile - file name to open
 * P - particle data to write to the file
 * startPart - particle to start writing at
 * NNeutrinos - number of neutrinos to write to this file
 * FirstID - ID to start writing IDs at
 * vel_prefac - Zeldovich prefactor for velocities
 * nupartmass - mass of a single (N-body) neutrino particle in internal gadget units
 */
int64_t write_neutrino_data(const std::string & SnapFile, lpt_data& outdata, part_grid& Pgrid, FermiDiracVel& nuvels, size_t startPart, uint32_t NNeutrinos,size_t NNuTotal, int64_t FirstId, const double vel_prefac, const double nupartmass, const double Box);



#endif
