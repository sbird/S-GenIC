#include <iostream>
#include <sstream>
#include <string>
#include <iomanip>
#include <math.h>
#include <stdlib.h>
#include "part_data.hpp"
#include "power.hpp"
#include "displacement.hpp"
#include "cosmology.hpp"
//For getopt
#include <unistd.h>
#include <hdf5.h>
#include <hdf5_hl.h>

#include <gadgetreader.hpp>

//Class to load the neutrino power spectra
class PowerSpec_NuTabulated: public PowerSpec
{
    public:
    PowerSpec_NuTabulated(const std::string& FileWithInputSpectrum);
    virtual double power(double k, int Type);
    virtual ~PowerSpec_NuTabulated(){}
    private:
        int NPowerTable;
        double scale;
        struct pow_matter *PowerMatter;
};

PowerSpec_NuTabulated::PowerSpec_NuTabulated(const std::string & FileWithInputSpectrum)
{
    //Note there is no need to convert units as we both power spectrum and snapshot come from the same place
    FILE *fd;
    if(!(fd = fopen(FileWithInputSpectrum.c_str(), "r"))) {
      printf("can't read input in file '%s'\n", FileWithInputSpectrum.c_str());
      exit(17);
    }
    NPowerTable = 0;
    while(true) {
      double ktab, Pktab;
      /* read TOTAL matter power spectrum from internal integrator*/
      if(fscanf(fd, " %lg %lg ", &ktab, &Pktab) != 2)
          break;
      NPowerTable++;
    }
    rewind(fd);
    if (NPowerTable == 0) {
      printf("Empty power spectrum file\n");
      exit(17);
    }
    printf("found %d rows in input SPECTRUM table\n", NPowerTable);
    PowerMatter = (pow_matter *) malloc(NPowerTable * sizeof(struct pow_matter));
    if (!PowerMatter){
        fprintf(stderr, "Could not allocate memory for table");
        exit(23);
    }
    /* define matter array */
    for (int count=0; count < NPowerTable; count++ ) {
        double kmat, pmat;
        if(fscanf(fd, " %lg %lg", &kmat, &pmat) != 2)
            break;
        PowerMatter[count].kmat = log10(kmat);
        PowerMatter[count].pmat = log10(pmat);
    }
    fclose(fd);
}

double PowerSpec_NuTabulated::power(double k, int Type)
{
  double logk, logD,  kold, u, dlogk, Delta2;
  int binlow, binhigh, binmid;

  kold = k;

  k *= scale;	/* convert to h/Mpc */

  logk = log10(k);

  if(logk < PowerMatter[0].kmat || logk > PowerMatter[NPowerTable - 1].kmat)
    return 0;

  binlow = 0;
  binhigh = NPowerTable - 1;

  while(binhigh - binlow > 1)
    {
      binmid = (binhigh + binlow) / 2;
      if(logk < PowerMatter[binmid].kmat)
      	binhigh = binmid;
      else
      	binlow = binmid;
    }

  dlogk = PowerMatter[binhigh].kmat - PowerMatter[binlow].kmat;

  if(dlogk == 0)
    exit(777);

  u = (logk - PowerMatter[binlow].kmat) / dlogk;

  logD = (1 - u) * PowerMatter[binlow].pmat + u * PowerMatter[binhigh].pmat;

  Delta2 = pow(10.0, logD);

  return Delta2 / (4 * M_PI * kold * kold * kold);
}


part_data generate_neutrino_particles(std::string GlassFile, std::string SimSpecFile, size_t NNeutrinos, size_t Nmesh, double Box, int Seed)
{
    //For the neutrinos we don't want this; the modes that we have are the modes that exist
    bool RayleighScatter = false;
    //Always initialise in a sphere
    bool SphereMode = true;

    //Load the neutrino power spectrum from the powerspec file
    PowerSpec_NuTabulated PSpec(SimSpecFile);

    //Open the glass file: this should be a regular grid and will provide the phase information for the neutrinos
    //FIXME: Get phase information from the snapshot somehow?
    std::cout<<"Initialising pre-IC file "<< GlassFile<<std::endl;
    GadgetReader::GSnap glass(GlassFile);
    printf("Nmesh = %lu Nsample = %lu\n",Nmesh,NNeutrinos);
    DisplacementFields displace(Nmesh, NNeutrinos, Seed, Box);
    //Output is neutrino particles
    int GlassTileFac = NNeutrinos/glass.GetNpart(2);
    part_data P(glass, 2, GlassTileFac, Box);
    displace.displacement_fields(2, NNeutrinos, P, &PSpec, SphereMode, RayleighScatter);
    return P;
}

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
        //If this is the first write, create the dataset
        if (begin==0) {
            size[0] = npart;
            hid_t space_id = H5Screate_simple(rank, size, NULL);
            H5Dcreate(group,BlockName.c_str(),d_type, space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            H5Sclose(space_id);
        }
        hid_t dset = H5Dopen(group,BlockName.c_str(),H5P_DEFAULT);
        if (dset < 0)
            return dset;
        size[0] = np_write;
        hsize_t begins[2]={begin,0};
        //Create a hyperslab that we will write to
        hid_t space_id = H5Screate_simple(rank, size, NULL);
        //Select the hyperslab of elements we are about to write to
        H5Sselect_hyperslab(space_id, H5S_SELECT_SET, begins, NULL, size, NULL);
        /* Write to the dataset */
        herr = H5Dwrite(dset, d_type, H5S_ALL, space_id, H5P_DEFAULT, data);
        H5Dclose(dset);
        if (herr < 0)
            return herr;
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
int64_t write_neutrino_data(const std::string & SnapFile, part_data P, size_t startPart, uint32_t NNeutrinos,int64_t FirstId, const double vel_prefac, const double nupartmass)
{
    //Open the HDF5 file
    hid_t handle = H5Fopen(SnapFile.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    //Create the group
    hid_t group = H5Gcreate(handle,"PartType2",H5P_DEFAULT, H5P_DEFAULT,H5P_DEFAULT);
    if(group < 0)
        return group;
    //We now want to create Coordinates, ParticleIDs and Velocities and write particles to them
    //Make a buffer of about 100 MB
    const uint32_t blockwrite = 3*(1<<26);
    float * buffer = new float[blockwrite];
    //Coordinates first : note no ++ below
    for (uint32_t curpart=0; curpart < NNeutrinos; ) {
        // Buffer initialization.
        for (uint32_t i = 0; i < blockwrite; ++i)
            for(int k = 0; k < 3; k++)
                buffer[3 * i + k] = P.Pos(startPart+i,k) + P.Vel(startPart+i,k);
        //Do the writing
        curpart += WriteBlock("Coordinates", group, 2, buffer, H5T_NATIVE_FLOAT, 3, NNeutrinos, blockwrite, curpart);
    }
    //Velocities next : note no ++ below
    for (uint32_t curpart=0; curpart < NNeutrinos; ) {
        // Buffer initialization
        for (uint32_t i = 0; i < blockwrite; ++i)
            for(int k = 0; k < 3; k++)
                buffer[3 * i + k] = vel_prefac * P.Vel(startPart+i,k);
        //Do the writing
        curpart += WriteBlock("Velocities", group, 2, buffer, H5T_NATIVE_FLOAT, 3, NNeutrinos, blockwrite, curpart);
    }
    delete[] buffer;
    H5Gclose(group);
    //Set the neutrino particle mass
    double masses[N_TYPE];
    H5LTget_attribute_double(handle,"Header","MassTable", masses);
    masses[2] = nupartmass;
    H5LTset_attribute_double(handle, "Header", "MassTable", masses, N_TYPE);
    H5Fclose(handle);
    return NNeutrinos;
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
            printf("-s : Directory with simulation snapshots\n -g : glassfile for neutrino particles\n -n : Snapshot output number\n -m : neutrino mass in eV\n");
            exit(0);
        }
    }
    if (BaseDir == NULL || GlassFile == NULL || snapnum == -1 || NUmass < 0.005) {
          printf("-s : Directory with simulation snapshots\n -g : glassfile for neutrino particles\n -n : Snapshot output number\n -m : neutrino mass in eV\n");
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
    size_t Nmesh, NNeutrinos;
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
        //FIXME: The output will have as many neutrino particles as there are CDM particles.
        //Ultimately we want to reduce this
        NNeutrinos = pow(Npart[1], 1./3);
        Nmesh = 2*NNeutrinos;
        H5LTget_attribute_double(hdf_group,".","BoxSize", &Box);
        H5LTget_attribute_double(hdf_group,".","BoxSize", &Box);
        H5LTget_attribute_int(hdf_group,".","NumFilesPerSnapshot",&numfiles);
        H5LTget_attribute_double(hdf_group,".","HubbleParam", &HubbleParam);
        H5LTget_attribute_double(hdf_group,".","Omega0", &Omega0);
        H5LTget_attribute_double(hdf_group,".","OmegaLambda", &OmegaLambda);
        H5LTget_attribute_double(hdf_group,".","Time",&atime);
    }
    //Note no hierarchy right now.
    Cosmology cosmo(HubbleParam, Omega0, OmegaLambda, NUmass, false);
    //Find the desired mass of each neutrino particle
//       omega =
//     masstot / (All.BoxSize * All.BoxSize * All.BoxSize) / (3 * All.Hubble * All.Hubble / (8 * M_PI * All.G));
    // OmegaNu is dimensionless. Box is in kpc, H in h/s, G in cm^-3 g^-1 s^-2
    // so this is in g kpc^3 cm^-3
    //So the conversion factor is
    const double nupartmass = cosmo.OmegaNu(atime) * pow(Box, 3) * 3 * pow(HUBBLE*HubbleParam, 2) / (8 * M_PI * GRAVITY) * pow(UnitLength_in_cm, 3) / UnitMass_in_g;
    const double hubble_a = cosmo.Hubble(atime) * UnitLength_in_cm / UnitVelocity_in_cm_per_s;
    const double vel_prefac = atime * hubble_a * cosmo.F_Omega(atime) /sqrt(atime);
    //This does the FFT
    part_data P  = generate_neutrino_particles(std::string(GlassFile), powerfile, NNeutrinos, Nmesh, Box, seed);
    //Choose a high ID number
    int64_t FirstId = NNeutrinos*8;
    //We need to open each snapshot file in turn and write neutrinos to it until we run out.
    size_t startPart = write_neutrino_data(SnapFile, P, 0, NNeutrinos,FirstId, vel_prefac, nupartmass);
    for(int i=1; i<numfiles; ++i)
    {
        formatter.str("");
        formatter << snapdir << "/snap_"<<snapnum_f<<"."<<i<<".hdf5";
        SnapFile = formatter.str();
        //Base this on the routine in save.cpp
        startPart += write_neutrino_data(SnapFile, P, startPart, NNeutrinos,FirstId+startPart, vel_prefac, nupartmass);
    }
    return 0;
}
