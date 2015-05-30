#include "createnu.hpp"
#include <iostream>
#include <sstream>
#include "displacement.hpp"
//For getopt
#include <hdf5.h>
#include <hdf5_hl.h>
#include <fstream>


PowerSpec_NuTabulated::PowerSpec_NuTabulated(const std::string & FileWithInputSpectrum, double atime)
{
    //Note there is no need to convert units as we both power spectrum and snapshot come from the same place
    //Open power spectrum file
    std::ifstream file;
    file.open(FileWithInputSpectrum.c_str());
    if(!file.is_open() ) {
      std::cerr<<"Can't read input power in file "<<FileWithInputSpectrum<<std::endl;
      throw std::bad_exception();
    }
    NPowerTable = 0;
    std::string line;
    //Read first two lines with metadata
    {
        std::getline(file, line);
        std::stringstream tokeniser;
        double newatime;
        tokeniser.str(line);
        tokeniser >> newatime;
        if (fabs(newatime/atime-1) > 0.01){
            std::cerr<<"Power spectrum file "<<FileWithInputSpectrum<<" is at a="<<newatime<<" not a="<<atime<<std::endl;
            throw std::bad_exception();
        }
        std::getline(file, line);
        tokeniser.clear();
        tokeniser.str(line);
        tokeniser >> NPowerTable;
    }
    std::cerr<<"Found "<<NPowerTable<<" rows in power spectrum "<<FileWithInputSpectrum<<std::endl;
    if (NPowerTable == 0) {
      throw std::bad_alloc();
    }
    pvals = new double[NPowerTable];
    kvals = new double [NPowerTable];
    /* define matter array */
    for(int count=0; count < NPowerTable; count++) {
        std::getline(file, line);
        std::stringstream tokeniser;
        double kmat, pmat;
        tokeniser.str(line);
        tokeniser >> kmat;
        tokeniser >> pmat;
        if (!tokeniser.fail()) {
            kvals[count] = log10(kmat);
            pvals[count] = log10(pmat);
        }
    }
    power_interp = gsl_interp_alloc(gsl_interp_cspline,NPowerTable);
    power_accel = gsl_interp_accel_alloc();
    gsl_interp_init(power_interp,kvals, pvals,NPowerTable);
}

double PowerSpec_NuTabulated::power(double k, int Type)
{
    double logk = log10(k);
    if (logk < kvals[0] || logk > kvals[NPowerTable-1])
        return 0;
    double logP = gsl_interp_eval(power_interp, kvals, pvals, logk, power_accel);
    double Delta2 = pow(10.0, logP);
    return Delta2;
}


part_data generate_neutrino_particles(std::string GlassFile, std::string SimSpecFile, size_t NNeutrinos, size_t Nmesh, double Box, int Seed, double atime)
{
    //For the neutrinos we don't want this; the modes that we have are the modes that exist
    bool RayleighScatter = false;

    //Load the neutrino power spectrum from the powerspec file
    PowerSpec_NuTabulated PSpec(SimSpecFile, atime);

    //Open the glass file: this should be a regular grid and will provide the phase information for the neutrinos
    //FIXME: Get phase information from the snapshot somehow?
    std::cout<<"Initialising pre-IC file "<< GlassFile<<std::endl;
    GadgetReader::GSnap glass(GlassFile);
    printf("Nmesh = %lu Nsample = %lu\n",Nmesh,NNeutrinos);
    DisplacementFields displace(Nmesh, Seed, Box, false);
    //Output is neutrino particles
    int GlassTileFac = NNeutrinos/pow(glass.GetNpart(2),1./3.);
    part_data P(glass, 2, GlassTileFac, Box, false);
    displace.displacement_fields(2, P, &PSpec, RayleighScatter);
    return P;
}

//Copied with minor modifications from GadgetReader
uint32_t WriteBlock(const std::string& BlockName, hid_t group, int type, void *data, hid_t d_type, int ncols, uint32_t npart, uint32_t np_write, uint32_t begin)
{
        herr_t herr;
        hsize_t size[2];
        int rank=1;
        size[1] = ncols;
        if (size[1] > 1) {
                rank = 2;
        }
        /* I don't totally understand why the below works (it is not clear to me from the documentation).
         * I gleaned it from a posting to the HDF5 mailing list and a related stack overflow thread here:
         * http://stackoverflow.com/questions/24883461/hdf5-updating-a-cell-in-a-table-of-integers
         * http://lists.hdfgroup.org/pipermail/hdf-forum_lists.hdfgroup.org/2014-July/007966.html
         * The important thing seems to be that we have a dataspace for the whole array and create a hyperslab on that dataspace.
         * Then we need another dataspace with the size of the stuff we want to write.*/
        //Create a dataspace corresponding the the whole extent in the file
        size[0] = npart;
        hid_t full_space_id = H5Screate_simple(rank, size, NULL);
        //If this is the first write, create the dataset
        if (begin==0) {
            //Create a hyperslab that we will write to
            hid_t herr = H5Dcreate(group,BlockName.c_str(),d_type, full_space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            if (herr < 0) {
                std::cerr << "Could not create dataset: " <<BlockName<<std::endl;
                throw std::bad_exception();
            }

            H5Dclose(herr);
        }
        hid_t dset = H5Dopen(group,BlockName.c_str(),H5P_DEFAULT);
        if (dset < 0)
            return dset;
        size[0] = np_write;
        hid_t space_id = H5Screate_simple(rank, size, NULL);
        hsize_t begins[2]={begin,0};
        //Select the hyperslab of elements we are about to write to
        H5Sselect_hyperslab(full_space_id, H5S_SELECT_SET, begins, NULL, size, NULL);
        /* Write to the dataset */
        herr = H5Dwrite(dset, d_type, space_id, full_space_id, H5P_DEFAULT, data);
        H5Dclose(dset);
        H5Sclose(space_id);
        H5Sclose(full_space_id);
        if (herr < 0) {
            std::cerr << "Could not write dataset: " <<BlockName<<std::endl;
            throw std::bad_exception();
        }
        return np_write;
}

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
int64_t write_neutrino_data(const std::string & SnapFile, part_data& P, FermiDiracVel& nuvels, size_t startPart, uint32_t NNeutrinos,size_t NNuTotal, int64_t FirstId, const double vel_prefac, const double nupartmass)
{
    //Open the HDF5 file
    hid_t handle = H5Fopen(SnapFile.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    //Create the group
    hid_t group = H5Gcreate(handle,"PartType2",H5P_DEFAULT, H5P_DEFAULT,H5P_DEFAULT);
    if(group < 0)
        return group;
    //We now want to create Coordinates, ParticleIDs and Velocities and write particles to them
    //Make a buffer of about 100 MB
    const uint32_t blockwrite = (1<<22);
    float * buffer = new float[3*blockwrite];
    //Coordinates first : note no ++ below
    for (uint32_t curpart=0; curpart < NNeutrinos; ) {
        // Buffer initialization.
        for (uint32_t i = 0; i < std::min(blockwrite, NNeutrinos-curpart); ++i)
            for(int k = 0; k < 3; k++)
                buffer[3 * i + k] = P.Pos(startPart+i,k) + P.Vel(startPart+i,k);
        //Do the writing
        curpart += WriteBlock("Coordinates", group, 2, buffer, H5T_NATIVE_FLOAT, 3, NNeutrinos, std::min(blockwrite, NNeutrinos-curpart), curpart);
    }
    //Velocities next : note no ++ below
    for (uint32_t curpart=0; curpart < NNeutrinos; ) {
        // Buffer initialization
        for (uint32_t i = 0; i < std::min(blockwrite, NNeutrinos-curpart); ++i)
            for(int k = 0; k < 3; k++) {
                buffer[3 * i + k] = vel_prefac * P.Vel(startPart+i,k);
                //Add thermal velocities
                nuvels.add_thermal_speeds(&buffer[3 * i]);
            }
        //Do the writing
        curpart += WriteBlock("Velocities", group, 2, buffer, H5T_NATIVE_FLOAT, 3, NNeutrinos, std::min(blockwrite, NNeutrinos-curpart), curpart);
    }
    delete[] buffer;
    int64_t* idbuffer = new int64_t[blockwrite];
    //Now do Particle IDs
    for (uint32_t curpart=0; curpart < NNeutrinos; ) {
        // Buffer initialization.
        for (uint32_t i = 0; i < std::min(blockwrite, NNeutrinos-curpart); ++i)
                idbuffer[i] = FirstId+i+curpart;
        //Do the writing
        curpart += WriteBlock("ParticleIDs", group, 2, idbuffer, H5T_NATIVE_LLONG, 1, NNeutrinos, std::min(blockwrite, NNeutrinos-curpart), curpart);
    }


    H5Gclose(group);
    //Set the neutrino particle mass
    double masses[N_TYPE];
    H5LTget_attribute_double(handle,"Header","MassTable", masses);
    masses[2] = nupartmass;
    H5LTset_attribute_double(handle, "Header", "MassTable", masses, N_TYPE);
    //Set the neutrino particle numbers
    uint32_t npart[N_TYPE];
    H5LTget_attribute_uint(handle,"Header","NumPart_ThisFile", npart);
    npart[2] = NNeutrinos;
    H5LTset_attribute_uint(handle, "Header", "NumPart_ThisFile", npart, N_TYPE);
    H5LTget_attribute_uint(handle,"Header","NumPart_Total", npart);
    npart[2] = NNuTotal - ((uint64_t)(NNuTotal >> 32) << 32);
    H5LTset_attribute_uint(handle, "Header", "NumPart_Total", npart, N_TYPE);
    H5LTget_attribute_uint(handle,"Header","NumPart_Total_HighWord", npart);
    npart[2] = (NNuTotal >> 32);
    H5LTset_attribute_uint(handle, "Header", "NumPart_Total_HighWord", npart, N_TYPE);
    H5Fclose(handle);
    return NNeutrinos;
}

