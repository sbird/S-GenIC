#include "allvars.h"
#include "proto.h"
#include "part_data.hpp"
#include "cosmology.hpp"
#include "power.hpp"

#include <gadgetreader.hpp>
#include <gadgetwriter.hpp>

int main(int argc, char **argv)
{
  int type;
  std::valarray<int64_t> npart(N_TYPE);
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
  read_parameterfile(argv[1]);
  set_units();

  printf("Nmesh = %lu Nsample = %lu\n",Nmesh,Nsample);
  if (Nmesh % 2 != 0){
    printf("Nmesh must be even or correct output is not guaranteed.\n");
    exit(1);
  }
  initialize_ffts();
  printf("Initialising pre-IC file '%s'\n",GlassFile);
  GadgetReader::GSnap snap(GlassFile);
  /*Set particle numbers*/
  for(type = 0; type < N_TYPE; type++)
    npart[type] = snap.GetNpart(type) * GlassTileFac * GlassTileFac * GlassTileFac;
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
      Cosmology cosmo(HubbleParam, Omega, OmegaLambda, MNu, InvertedHierarchy);
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
  if(osnap.WriteHeaders(generate_header(npart)))
          FatalError(23);

  for(type=0; type<N_TYPE;type++){
      int64_t NumPart = 0;
      if(npart[type] == 0)
              continue;
      try{
        part_data P(snap, type, GlassTileFac, Box);
        NumPart = P.GetNumPart();
        displacement_fields(type, NumPart, P, PSpec, Nmesh, RayleighScatter);
        FirstId = write_particle_data(osnap, type,P, NumPart,FirstId);
      }
      catch (std::bad_alloc& ba)
      {
         size_t mem = sizeof(float)*snap.GetNpart(type)*3*GlassTileFac*GlassTileFac*GlassTileFac/1024/1024;
#ifdef TWOLPT
#ifdef NEUTRINOS
         if (type != 2)
#endif
             mem*=2;
#endif
         fprintf(stderr, "Could not allocate %ld MB for particle velocities\n", mem);
         FatalError(24);
      }
#ifdef PRINT_SPEC
      print_spec(type);
#endif
  }

  fftwf_free(Disp);
  fftwf_destroy_plan(Inverse_plan);
#ifdef TWOLPT
  /* Free  */
  fftwf_free(twosrc);
  fftwf_destroy_plan(Forward_plan2);
  for(int i=0;i<3;i++){
        fftwf_free(cdigrad[i]);
        fftwf_destroy_plan(Inverse_plan_grad[i]);
  }
#endif

  printf("Initial scale factor = %g\n", InitTime);

  return 0;
}
