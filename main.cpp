#include "part_data.hpp"
#include "cosmology.hpp"
#include "power.hpp"
#include "displacement.hpp"
#include "read_param.hpp"
#include "thermalvel.hpp"
#include "save.hpp"
#include <stdexcept>
#include <cassert>

int main(int argc, char **argv)
{
  if(argc < 2)
    {
	  fprintf(stdout, "\nParameters are missing.\n");
	  fprintf(stdout, "Call with <ParameterFile>\n\n");
      exit(0);
    }
  /*Make sure stdout is line buffered even when not
   * printing to a terminal but, eg, perl*/
  setlinebuf(stdout);
  //May throw std::runtime_error
  try {
  //Read the config file
  SpbConfigParser config(argv[1]);
  //Cosmological parameters
  const auto Omega = config.PopValue<double>("Omega");
  const auto OmegaLambda = config.PopValue<double>("OmegaLambda");
  const auto OmegaBaryon = config.PopValue<double>("OmegaBaryon");
  const auto HubbleParam = config.PopValue<double>("HubbleParam");
  //Which output format should we use. 4 is bigfile, 3 is HDF5, 2 is Gadget 2. Default to 3.
  const auto ICFormat = config.PopValue<int>("ICFormat", 3);
  //How many output files to use in the set
  const auto NumFiles = config.PopValue<int>("NumFiles");
  //Parameters of the simulation box
  const auto Box = config.PopValue<double>("Box");
  const auto Redshift = config.PopValue<double>("Redshift");
  const double InitTime = 1 / (1 + Redshift);
  //Size of FFT
  const size_t Nmesh = config.PopValue<int>("Nmesh");
  //Unused unless CAMB spectrum
  const auto FileWithInputSpectrum = config.PopValue<std::string>("FileWithInputSpectrum");
  const auto FileWithTransfer = config.PopValue<std::string>("FileWithTransfer");
  const auto InputSpectrum_UnitLength_in_cm = config.PopValue<double>("InputSpectrum_UnitLength_in_cm", 3.085678e24);
  //Output filenames
  const auto OutputDir = config.PopValue<std::string>("OutputDir");
  const auto FileBase = config.PopValue<std::string>("FileBase");
  //Random number seed
  const auto Seed = config.PopValue<int>("Seed");
  //Various boolean flags
  const auto ReNormalizeInputSpectrum = config.PopValue<bool>("ReNormalizeInputSpectrum", false);
  const auto RayleighScatter = config.PopValue<bool>("RayleighScatter", true);
  //Power spectrum to use. Default to CAMB
  const auto WhichSpectrum = config.PopValue<int>("WhichSpectrum", 2);
  //Is twolpt on?
  const auto twolpt = config.PopValue<bool>("TWOLPT",true);
  //Unit system
  const auto UnitLength_in_cm = config.PopValue<double>("UnitLength_in_cm", 3.085678e21);
  const auto UnitVelocity_in_cm_per_s = config.PopValue<double>("UnitVelocity_in_cm_per_s", 1e5);
  const auto UnitMass_in_g = config.PopValue<double>("UnitMass_in_g", 1.989e43);
  const double UnitTime_in_s = UnitLength_in_cm / UnitVelocity_in_cm_per_s;
  //WDM options
  bool WDM_Vtherm_On = config.PopValue<bool>("WDM_Vtherm_On",false);
  const auto WDM_PartMass_in_kev = config.PopValue<double>("WDM_PartMass_in_kev", 0);
  //Neutrino options
  //Add thermal velocities to type 2 particles
  bool NU_Vtherm_On = config.PopValue<bool>("NU_Vtherm_On",true);
  //This triggers the use of neutrinos via an altered transfer function
  //Should be on only if you are faking neutrinos by combining them with the dark matter,
  //and changing the transfer function, which is a terrible way of simulating neutrinos. So leave it off.
  const auto combined_neutrinos = config.PopValue<bool>("NU_in_DM",false);
  //Changes whether we have two heavy and one light neutrino or one heavy two light.
  //0 is 3 degenerate species, 1 is normal (two light species) and -1 is inverted (two heavy)
  const auto Hierarchy = config.PopValue<int>("Hierarchy",0);
  //If enabled, turn off radiation in the initial velocities
  const auto NoRadiation = config.PopValue<bool>("NoRadiation",false);
  //Total neutrino mass
  const auto NU_PartMass_in_ev = config.PopValue<double>("NU_PartMass_in_ev",0);
  //Parameter for the Efstathiou power spectrum. Generally does nothing.
  const auto ShapeGamma = config.PopValue<double>("ShapeGamma",0.201);
  //Needed if ReNormaliseInputSpectrum is on. Otherwise unused
  const auto PrimordialIndex = config.PopValue<double>("PrimordialIndex",1.);
  const auto Sigma8 = config.PopValue<double>("Sigma8",0.8);
  //Number of particles desired
  std::valarray<int64_t> npart((int64_t)0,(size_t)N_TYPE);
  int64_t CbRtNpart[6] = {0};
  CbRtNpart[BARYON_TYPE] = config.PopValue<int>("NBaryon", 0);
  CbRtNpart[DM_TYPE] = config.PopValue<int>("NCDM", 0);
  CbRtNpart[NEUTRINO_TYPE] = config.PopValue<int>("NNeutrino", 0);
  for(int type=0; type<N_TYPES; ++type)
      npart[type] = static_cast<int64_t>(CbRtNpart[type])*CbRtNpart[type]*CbRtNpart[type];

  //Make a cosmology
  Cosmology cosmo(HubbleParam, Omega, OmegaLambda, NU_PartMass_in_ev, Hierarchy, NoRadiation);
  //Init structure for neutrino velocities
  //Maximum velocity to sample Fermi-Dirac from, in km/s at z=0
  const double v_th = cosmo.NU_V0(Redshift, UnitVelocity_in_cm_per_s);
  //Convert physical km/s at z=0 in an unperturbed universe to internal gadget (comoving) velocity units at starting redshift.
  double vnumax = config.PopValue<double>("Max_nuvel", -1)*pow((1+Redshift),1.5)*(UnitVelocity_in_cm_per_s/1e5);
  if (vnumax < 0 || vnumax > v_th*MAX_FERMI_DIRAC)
      vnumax = v_th*MAX_FERMI_DIRAC;

  printf("Particle numbers: %ld %ld %ld\n",npart[BARYON_TYPE], npart[DM_TYPE], npart[NEUTRINO_TYPE]);
  assert(npart[BARYON_TYPE] > 0 || npart[DM_TYPE] > 0 || npart[NEUTRINO_TYPE] > 0);

  if (Nmesh % 2 != 0){
    printf("Nmesh must be even or correct output is not guaranteed.\n");
    exit(1);
  }

  std::vector<std::string> unexpected = config.GetRemainingKeys();
  if(unexpected.size() > 0){
    std::cerr<<"Config file contained the following unexpected keys:"<<std::endl;
    for(auto unex: unexpected)
        std::cerr<<unex<<std::endl;
    exit(1);
  }
  DisplacementFields displace(Nmesh, Seed, Box, twolpt);
  /*Set particle numbers*/
  if(npart.sum() == 0)
          exit(1);
  //Set the flag that checks whether neutrinos are free-streaming
  cosmo.SetNeutrinoFreeStream(Box*UnitLength_in_cm, v_th * UnitVelocity_in_cm_per_s, InitTime);
  //Initialise a power spectrum
  PowerSpec * PSpec;
  switch(WhichSpectrum)
  {
      case 1:
            PSpec = new PowerSpec_EH(HubbleParam, Omega, OmegaBaryon, UnitLength_in_cm);
            break;
      case 2:
            PSpec = new PowerSpec_Tabulated(FileWithTransfer, FileWithInputSpectrum, Omega, OmegaLambda, OmegaBaryon, cosmo.OmegaNu(1),InputSpectrum_UnitLength_in_cm, UnitLength_in_cm, !npart[BARYON_TYPE], combined_neutrinos);
            break;
      default:
            PSpec = new PowerSpec_Efstathiou(ShapeGamma, UnitLength_in_cm);
  }

  //If normalisation or WDM are on, decorate the base power spectrum
  //to do that
  if (ReNormalizeInputSpectrum) {
      double Dplus = cosmo.GrowthFactor(InitTime, 1.0);
      printf("Growth factor to z=0: %g \n", Dplus);
      PSpec = new NormalizedPowerSpec(PSpec, Sigma8, PrimordialIndex, Dplus, UnitLength_in_cm);
  }
  if(WDM_Vtherm_On)
      PSpec = new WDMPowerSpec(PSpec, WDM_PartMass_in_kev, Omega, OmegaBaryon, HubbleParam, UnitLength_in_cm);

  std::string extension("");
  if (ICFormat > 4 || ICFormat < 2) {
          fprintf(stderr, "Supported ICFormats:\n 2: Gadget 2 format files\n 3: HDF5\n 4: BigFile\n");
          exit(1);
  }
  if (ICFormat == 3){
      printf("Outputting HDF5 ICs\n");
      extension=".hdf5";
  }
  GadgetWriter::GWriteBaseSnap *osnap;
#ifdef HAVE_BGFL
  if(ICFormat == 4) {
     osnap = new GadgetWriter::GWriteBigSnap(OutputDir+std::string("/")+FileBase+extension, npart, NumFiles);
  }
  else
#endif
  osnap = new GadgetWriter::GWriteSnap(OutputDir+std::string("/")+FileBase+extension, npart,NumFiles, sizeof(id_type));
  assert(osnap);
  /*Write headers*/
  gadget_header header = generate_header(npart, Omega, OmegaBaryon, cosmo.OmegaNu(1), OmegaLambda, HubbleParam, Box, InitTime, UnitMass_in_g, UnitLength_in_cm, UnitVelocity_in_cm_per_s, combined_neutrinos);

  if ( vnumax/v_th < MAX_FERMI_DIRAC ) {
      FermiDiracVel therm_vels (v_th, vnumax/v_th);
      header.mass[NEUTRINO_TYPE]*=therm_vels.total_frac;
      printf("F-D velocity scale: %g (km/s). Max particle vel: %g (km/s). Fraction of mass in particles: %g\n",v_th/sqrt(1+Redshift), vnumax/sqrt(1+Redshift), therm_vels.total_frac);
  }

  //Generate regular particle grid
  part_grid Pgrid(CbRtNpart, header.mass, Box);

  if(osnap->WriteHeaders(header)) {
      fprintf(stderr, "Could not write headers to snapshot\n");
      exit(1);
  }

  int64_t FirstId=1;
  //Compute the factors to go from velocity to displacement
  const double hubble_a = cosmo.Hubble(InitTime)*UnitTime_in_s;
  const double vel_prefac = InitTime * hubble_a * cosmo.F_Omega(InitTime) /sqrt(InitTime);
  //Only used if twolpt is on
  //This is slightly approximate: we are assuming that D2 ~ -3/7 Da^2 Omega_m^{-1/143} (Bouchet, F 1995, A&A 296)
  const double vel_prefac2 = -3./7.*pow(Omega, -1./143)*InitTime * hubble_a * cosmo.F2_Omega(InitTime) /sqrt(InitTime);
  printf("vel_prefac= %g  hubble_a=%g fom=%g Omega=%g \n", vel_prefac, hubble_a, cosmo.F_Omega(InitTime), Omega);

  for(int type=0; type<N_TYPE;type++){
      if(npart[type] == 0)
              continue;
      FermiDiracVel * therm_vels = NULL;
      //For WDM thermal velocities
      if(WDM_Vtherm_On && type == DM_TYPE){
        const double wdm_vth = WDM_V0(Redshift, WDM_PartMass_in_kev, Omega-OmegaBaryon, HubbleParam, UnitVelocity_in_cm_per_s);
        therm_vels = new FermiDiracVel (wdm_vth);
        printf("\nWarm dark matter rms velocity dispersion at starting redshift = %g km/sec\n\n",3.59714 * wdm_vth);
      }
      //Neutrino thermal velocities
      if(NU_Vtherm_On && type == NEUTRINO_TYPE) {
          therm_vels = new FermiDiracVel (v_th, vnumax/v_th);
          printf("\nNeutrino rms vel. dispersion %g (km/s)\n\n",v_th/sqrt(1+Redshift));
      }
      //Only compute zeldovich displacements if we have all the neutrinos in particles.
      if(type != 2 || !therm_vels || therm_vels->total_frac > 0.95) {
        lpt_data outdata = displace.displacement_fields(type, Pgrid, PSpec, RayleighScatter);
        outdata.SetVelPrefac(vel_prefac, vel_prefac2);
        FirstId = write_particle_data(*osnap, type,&outdata, Pgrid, therm_vels, FirstId);
      }
      //Otherwise pass empty values for the velocity fields.
      else
        FirstId = write_particle_data(*osnap, type,NULL, Pgrid, therm_vels, FirstId);
      delete therm_vels;
  }

    delete osnap;
    delete PSpec;
  printf("Initial scale factor = %g\n", InitTime);

  }
  catch(std::ios_base::failure& e) {
      std::cerr<<e.what();
      return 4;
  }
  catch(std::runtime_error& e) {
      std::cerr<<e.what();
      return 1;
  }
  catch(std::invalid_argument& e) {
      std::cerr<<e.what();
      return 2;
  }
  catch(std::domain_error& e) {
      std::cerr<<e.what();
      return 3;
  }
  return 0;
}
