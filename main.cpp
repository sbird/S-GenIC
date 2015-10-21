#include "allvars.h"
#include "proto.h"
#include "part_data.hpp"
#include "cosmology.hpp"
#include "power.hpp"
#include "displacement.hpp"
#include "read_param.hpp"
#include <cassert>

#include <gadgetwriter.hpp>

int main(int argc, char **argv)
{
  int type;
  std::valarray<int64_t> npart((int64_t)0,(size_t)N_TYPE);
  int64_t FirstId=1;

  if(argc < 2)
    {
	  fprintf(stdout, "\nParameters are missing.\n");
	  fprintf(stdout, "Call with <ParameterFile>\n\n");
      exit(0);
    }
  /*Make sure stdout is line buffered even when not 
   * printing to a terminal but, eg, perl*/
  setlinebuf(stdout);
  SpbConfigParser config(argv[1]);
  std::map<std::string, ValueTuple> configoptions;
  //Cosmological parameters
  configoptions["Omega"] = std::make_tuple((void *) &Omega, FloatType, "");
  configoptions["OmegaLambda"] = std::make_tuple((void *) &OmegaLambda, FloatType, "");
  configoptions["OmegaBaryon"] = std::make_tuple((void *) &OmegaBaryon, FloatType, "");
  configoptions["OmegaDM_2ndSpecies"] = std::make_tuple((void *) &OmegaDM_2ndSpecies, FloatType, "");
  configoptions["HubbleParam"] = std::make_tuple((void *) &HubbleParam, FloatType, "");
  //Which output format should we use. 3 is HDF5, 2 is Gadget 2
  configoptions["ICFormat"] = std::make_tuple((void *) &ICFormat, IntType, "3");
  //How many output files to use in the set
  configoptions["NumFiles"] = std::make_tuple((void *) &NumFiles, IntType, "");
  //Parameters of the simulation box
  configoptions["Redshift"] = std::make_tuple((void *) &Redshift, FloatType, "");
  configoptions["Box"] = std::make_tuple((void *) &Box, FloatType, "");
  //Numerical FFT parameters
  configoptions["Nmesh"] = std::make_tuple((void *) &Nmesh, IntType, "");
  //Unused unless CAMB spectrum
  configoptions["FileWithInputSpectrum"] = std::make_tuple((void *) &FileWithInputSpectrum, StringType, "");
  configoptions["FileWithTransfer"] = std::make_tuple((void *) &FileWithTransfer, StringType, "");
  configoptions["OutputDir"] = std::make_tuple((void *) &OutputDir, StringType, "");
  configoptions["FileBase"] = std::make_tuple((void *) &FileBase, StringType, "");
  //Random number seed
  configoptions["Seed"] = std::make_tuple((void *) &Seed, IntType, "");
  //Various boolean flags
  configoptions["ReNormalizeInputSpectrum"] = std::make_tuple((void *) &ReNormalizeInputSpectrum, IntType, "0");
  configoptions["RayleighScatter"] = std::make_tuple((void *) &RayleighScatter, IntType, "1");
  //Power spectrum to use. Default to CAMB
  configoptions["WhichSpectrum"] = std::make_tuple((void *) &WhichSpectrum, IntType, "2");
  configoptions["TWOLPT"] = std::make_tuple((void *) &twolpt, IntType, "1");
  //Unit system
  configoptions["UnitVelocity_in_cm_per_s"] = std::make_tuple((void *) &UnitVelocity_in_cm_per_s, FloatType, "1e5");
  configoptions["UnitLength_in_cm"] = std::make_tuple((void *) &UnitLength_in_cm, FloatType, "3.085678e21");
  configoptions["UnitMass_in_g"] = std::make_tuple((void *) &UnitMass_in_g, FloatType, "1.989e43");
  configoptions["InputSpectrum_UnitLength_in_cm"] = std::make_tuple((void *) &InputSpectrum_UnitLength_in_cm, FloatType, "3.085678e24");
  //WDM options
  configoptions["WDM_On"] = std::make_tuple((void *) &WDM_On, IntType, "0");
  configoptions["WDM_Vtherm_On"] = std::make_tuple((void *) &WDM_Vtherm_On, IntType, "0");
  configoptions["WDM_PartMass_in_kev"] = std::make_tuple((void *) &WDM_PartMass_in_kev, FloatType, "0");
  //Neutrino options
  configoptions["NU_On"] = std::make_tuple((void *) &NU_On, IntType, "0");
  configoptions["NU_Vtherm_On"] = std::make_tuple((void *) &NU_Vtherm_On, IntType, "1");
  configoptions["NU_KSPACE"] = std::make_tuple((void *) &neutrinos_ks, IntType, "1");
  configoptions["NU_PartMass_in_ev"] = std::make_tuple((void *) &NU_PartMass_in_ev, FloatType, "0");
  //Parameter for the Efstathiou power spectrum. Generally does nothing.
  configoptions["ShapeGamma"] = std::make_tuple((void *) &ShapeGamma, FloatType, "0.201");
  //Needed if ReNormaliseInputSpectrum is on. Otherwise unused
  configoptions["Sigma8"] = std::make_tuple((void *) &Sigma8, FloatType, "0.8");
  configoptions["PrimordialIndex"] = std::make_tuple((void *) &PrimordialIndex, FloatType, "1.");
  //Number of particles desired
  int CbRtNpart[6]={0,0,0,0,0,0};
  configoptions["NBaryon"] = std::make_tuple((void *) &CbRtNpart[0], IntType, "0");
  configoptions["NCDM"] = std::make_tuple((void *) &CbRtNpart[1], IntType, "0");
#ifdef NEUTRINOS
  configoptions["NNeutrino"] = std::make_tuple((void *) &CbRtNpart[2], IntType, "0");
#endif
  config.parameter_parser(configoptions);
  set_units();

  npart[BARYON_TYPE] = static_cast<int64_t>(CbRtNpart[0])*static_cast<int64_t>(CbRtNpart[0])*CbRtNpart[0];
  npart[DM_TYPE] = static_cast<int64_t>(CbRtNpart[1])*static_cast<int64_t>(CbRtNpart[1])*CbRtNpart[1];
  npart[NEUTRINO_TYPE] = static_cast<int64_t>(CbRtNpart[2])*CbRtNpart[2]*CbRtNpart[2];
  printf("Particle numbers: %ld %ld %ld\n",npart[BARYON_TYPE], npart[DM_TYPE], npart[NEUTRINO_TYPE]);
  assert(npart[BARYON_TYPE] > 0 || npart[DM_TYPE] > 0 || npart[NEUTRINO_TYPE] > 0);

  if (Nmesh % 2 != 0){
    printf("Nmesh must be even or correct output is not guaranteed.\n");
    exit(1);
  }

  DisplacementFields displace(Nmesh, Seed, Box, twolpt);
  /*Set particle numbers*/
  if(npart.sum() == 0)
          exit(1);
  /*Set the global variable saying there is no gas in the glassfile,
   * so that the OmegaBaryon should be added to the DM.*/
  if (npart[BARYON_TYPE] == 0)
          no_gas = 1;
  else
          no_gas = 0;
  //Initialise a power spectrum
  PowerSpec * PSpec;
  switch(WhichSpectrum)
  {
      case 1:
            PSpec = new PowerSpec_EH(HubbleParam, Omega, OmegaBaryon, UnitLength_in_cm);
            break;
      case 2:
            PSpec = new PowerSpec_Tabulated(FileWithTransfer, FileWithInputSpectrum, Omega, OmegaLambda, OmegaBaryon, OmegaDM_2ndSpecies,InputSpectrum_UnitLength_in_cm, UnitLength_in_cm, no_gas, neutrinos_ks);
            break;
      default:
            PSpec = new PowerSpec_Efstathiou(ShapeGamma, UnitLength_in_cm);
  }

  //If normalisation or WDM are on, decorate the base power spectrum
  //to do that
  if (ReNormalizeInputSpectrum) {
      Cosmology cosmo(HubbleParam, Omega, OmegaLambda, NU_PartMass_in_ev, InvertedHierarchy);
      PSpec = new NormalizedPowerSpec(PSpec, Sigma8, PrimordialIndex, cosmo.GrowthFactor(InitTime, 1.0), UnitLength_in_cm);
  }
  if(WDM_On)
      PSpec = new WDMPowerSpec(PSpec, WDM_PartMass_in_kev, Omega, OmegaBaryon, HubbleParam, UnitLength_in_cm);

  std::string extension("");
  if (ICFormat > 3 || ICFormat < 2) {
          fprintf(stdout, "Only Gadget format 2 or HDF5 output is supported.\n");
          exit(1);
  }
  if (ICFormat == 3){
      printf("Outputting HDF5 ICs\n");
      extension=".hdf5";
  }

#ifdef NEUTRINO_PAIRS
  npart[NEUTRINO_TYPE] *= 2;
#endif //NEUTRINO_PAIRS
  GadgetWriter::GWriteSnap osnap(std::string(OutputDir)+std::string("/")+std::string(FileBase)+extension, npart,NumFiles, sizeof(id_type));
  /*Write headers*/
  gadget_header header = generate_header(npart, Omega, OmegaBaryon, OmegaDM_2ndSpecies, OmegaLambda, HubbleParam, Box, InitTime, UnitMass_in_g, UnitLength_in_cm, neutrinos_ks);

  //Generate regular particle grid
  part_grid Pgrid(CbRtNpart, header.mass, Box);

  if(osnap.WriteHeaders(header))
          FatalError(23);

  for(type=0; type<N_TYPE;type++){
      if(npart[type] == 0)
              continue;
      try{
          bool tlptpart = twolpt;
#ifdef NEUTRINOS
          if (type==2)
              tlptpart = false;
#endif
          lpt_data outdata = displace.displacement_fields(type, Pgrid, PSpec, RayleighScatter);
          FirstId = write_particle_data(osnap, type,outdata, Pgrid, FirstId, tlptpart);
      }
      catch (std::bad_alloc& ba)
      {
         size_t mem = sizeof(float)*pow(Pgrid.GetNumPart(type),3);
         if(twolpt)
#ifdef NEUTRINOS
            if (type != 2)
#endif
                mem*=2;
         fprintf(stderr, "Could not allocate %ld MB for particle velocities\n", mem);
         FatalError(24);
      }
#ifdef PRINT_SPEC
      print_spec(type);
#endif
  }

    delete PSpec;
  printf("Initial scale factor = %g\n", InitTime);

  return 0;
}
