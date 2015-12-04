#ifndef CREATENU_H
#define CREATENU_H

#include "part_data.hpp"
#include "thermalvel.hpp"
#include <hdf5.h>
#include <hdf5_hl.h>

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
int64_t write_neutrino_data(const std::string & SnapFile, part_grid& Pgrid, FermiDiracVel& nuvels, size_t startPart, uint32_t NNeutrinos,size_t NNuTotal, int64_t FirstId, const double nupartmass, const double Box);



#endif
