#include <iostream>
#include <sstream>
#include <iomanip>
#include "createnu.hpp"
#include "cosmology.hpp"
//For getopt
#include <unistd.h>
#include <hdf5.h>
#include <hdf5_hl.h>

#include <gadgetreader.hpp>
//Check whether a filename exists
inline bool file_exists (const std::string& name)
{
    if (FILE *file = fopen(name.c_str(), "r")) {
        fclose(file);
        return true;
    } else {
        return false;
    }
}

int main(int argc, char **argv)
{
    //Get the input variables
    int c;
    int snapnum=-1;
    char * BaseDir = NULL;
    char * GlassFile = NULL;
    int numfiles = 0;
    double NUmass=0.;
    const double seed = 23;
    //Default units
    const double UnitLength_in_cm = 3.085678e21;
    const double UnitVelocity_in_cm_per_s = 1e5;
    const double UnitMass_in_g = 1.989e43;

    while ((c = getopt (argc, argv, "s:g:n:m:")) != -1) {
        switch (c)
        {
        case 's':
            BaseDir = optarg;
            break;
        case 'n':
            snapnum = atoi(optarg);
            break;
        case 'g':
            GlassFile = optarg;
            break;
        case 'm':
            NUmass = atof(optarg);
            break;
        default:
            printf("-s : Directory with simulation snapshots\n -g : glassfile for neutrino particles\n -n : Snapshot output number\n -m : total neutrino mass in eV\n");
            exit(0);
        }
    }
    if (BaseDir == NULL || GlassFile == NULL || snapnum == -1 || NUmass < 0.005) {
          printf("-s : Directory with simulation snapshots\n -g : glassfile for neutrino particles\n -n : Snapshot output number\n -m : total neutrino mass in eV\n");
          exit(0);
    }
    std::ostringstream formatter;
    formatter << std::setfill('0')<<std::setw(3)<<snapnum;
    //Contains the formatter snapshot number
    const std::string snapnum_f( formatter.str() );
    formatter.str("");
    formatter << std::string(BaseDir)<<"/snapdir_"<<snapnum_f;
    //Contains the snapshot directory
    const std::string snapdir( formatter.str() );
    //Contains the power spectrum file
    formatter.str("");
    formatter << std::string(BaseDir)<<"/powerspec_nu_"<< snapnum_f <<".txt";
    const std::string powerfile ( formatter.str() );
    //Open the first file of the simulation snapshots and read metadata from there
    formatter.str("");
    formatter << snapdir;
    //Check for single file snapshot sets
    formatter << "/snap_"<<snapnum_f <<".hdf5";
    std::string SnapFile ( formatter.str() );
    if ( file_exists(SnapFile) ) {
        numfiles = 1;
    }
    else {
        formatter.str("");
        formatter << snapdir << "/snap_"<<snapnum_f<<".0.hdf5";
        SnapFile = formatter.str();
        if (!file_exists(SnapFile)) {
            std::cout <<"Neither " << SnapFile << "nor the single-file version exist. Your snapshot set must be HDF5."<<std::endl;
            exit(1);
        }
    }
    double Box, HubbleParam, Omega0, OmegaLambda, atime;
    size_t Npart[N_TYPE];

    {
        uint32_t npart[N_TYPE];
        hid_t handle = H5Fopen(SnapFile.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
        hid_t hdf_group =H5Gopen(handle,"/Header",H5P_DEFAULT);
        /*Get the total number of particles*/
        H5LTget_attribute(hdf_group,".","NumPart_Total",H5T_NATIVE_INT, &npart);
        for(int i = 0; i< N_TYPE; i++)
                Npart[i]=npart[i];
        H5LTget_attribute(hdf_group,".","NumPart_Total_HighWord",H5T_NATIVE_INT, &npart);
        for(int i = 0; i< N_TYPE; i++)
                Npart[i]+=(1L<<32)*npart[i];
        H5LTget_attribute_double(hdf_group,".","BoxSize", &Box);
        H5LTget_attribute_double(hdf_group,".","BoxSize", &Box);
        H5LTget_attribute_int(hdf_group,".","NumFilesPerSnapshot",&numfiles);
        H5LTget_attribute_double(hdf_group,".","HubbleParam", &HubbleParam);
        H5LTget_attribute_double(hdf_group,".","Omega0", &Omega0);
        H5LTget_attribute_double(hdf_group,".","OmegaLambda", &OmegaLambda);
        H5LTget_attribute_double(hdf_group,".","Time",&atime);
        H5Gclose(hdf_group);
        H5Fclose(handle);
    }
    //FIXME: The output will have as many neutrino particles as there are CDM particles.
    //Ultimately we want to reduce this
    const size_t Nsample = round(pow(Npart[1], 1./3));
    const size_t NNeutrinos = Npart[1];
    const size_t Nmesh = Nsample;

    //Note no hierarchy right now.
    Cosmology cosmo(HubbleParam, Omega0, OmegaLambda, NUmass, false);
    //Find the desired mass of each neutrino particle
    // OmegaNu is dimensionless. Box is in kpc/h, H in h/s, G in cm^-3 g^-1 s^-2
    // so this is in g kpc^3 cm^-3
    //So the conversion factor is (kpc/cm)^3 * g/Msun)
    const double nupartmass = cosmo.OmegaNu(1.) * pow(Box, 3) * 3 * pow(HUBBLE, 2) / (8 * M_PI * GRAVITY) * pow(UnitLength_in_cm, 3) / UnitMass_in_g / NNeutrinos;
    const double hubble_a = cosmo.Hubble(atime) * UnitLength_in_cm / UnitVelocity_in_cm_per_s;
    const double vel_prefac = atime * hubble_a * cosmo.F_Omega(atime) /sqrt(atime);
    std::cout<<"Omega Nu is "<<cosmo.OmegaNu(1.)<<" at z=0, so particle mass is "<<nupartmass<<std::endl;
    //This does the FFT
    part_data P  = generate_neutrino_particles(std::string(GlassFile), powerfile, Nsample, Nmesh, Box, seed, atime);
    //Choose a high ID number
    int64_t FirstId = (Npart[0]+Npart[1]+Npart[3]+Npart[4]+Npart[5])*8;
    const double v_th = NU_V0(1./atime-1, NUmass, UnitVelocity_in_cm_per_s);
    FermiDiracVel nuvels (v_th);
    //We need to open each snapshot file in turn and write neutrinos to it until we run out.
    uint32_t Nthisfile = NNeutrinos/numfiles;
    size_t startPart = write_neutrino_data(SnapFile, P, nuvels, 0, Nthisfile, NNeutrinos,FirstId, vel_prefac, nupartmass);
    for(int i=1; i<numfiles; ++i)
    {
        //If this is the last file, write all remaining particles
        if (i == numfiles -1) {
            Nthisfile = NNeutrinos - startPart;
        }
        formatter.str("");
        formatter << snapdir << "/snap_"<<snapnum_f<<"."<<i<<".hdf5";
        SnapFile = formatter.str();
        //Base this on the routine in save.cpp
        startPart += write_neutrino_data(SnapFile, P, nuvels, startPart, Nthisfile, NNeutrinos,FirstId+startPart, vel_prefac, nupartmass);
    }
    return 0;
}

