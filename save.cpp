#include "allvars.h"
#include "proto.h"
#include "part_data.hpp"
#include "thermalvel.hpp"
#include <gadgetwriter.hpp>

using namespace GadgetWriter;
using namespace std;

#define BUFFER 48

int64_t write_particle_data(GWriteSnap & snap, int type, lpt_data& outdata, part_grid& Pgrid, const double vel_prefac, const double vel_prefac2, int64_t FirstId, const bool twolpt)
{
  const int64_t NumPart = outdata.GetNumPart();
  float *block;
  id_type *blockid;
  int64_t written=0, pc;
  std::string posstr("POS ");
  std::string velstr("VEL ");
  std::string idstr("ID  ");
  std::string ustr("U   ");
  //HDF5 has different strings
  if(ICFormat == 3){
      posstr = "Coordinates";
      velstr = "Velocities";
      idstr = "ParticleIDs";
      ustr = "InternalEnergy";
  }
#ifdef NEUTRINOS
  //Init structure for neutrino velocities
  const double v_th = NU_V0(Redshift, NU_PartMass_in_ev, UnitVelocity_in_cm_per_s);
  FermiDiracVel nuvels (v_th);
  if(NU_On == 1 && NU_Vtherm_On == 1) {
    printf("\nNeutrino rms vel. dispersion %g (km/s)\n\n",v_th/sqrt(1+Redshift));
  }
#endif //NEUTRINOS
  //For WDM thermal velocities
  const double wdm_vth = WDM_V0(Redshift, WDM_PartMass_in_kev, Omega-OmegaBaryon, HubbleParam, UnitVelocity_in_cm_per_s);
  FermiDiracVel WDMvels (wdm_vth);
  if(WDM_On == 1 && WDM_Vtherm_On == 1)
        printf("\nWarm dark matter rms velocity dispersion at starting redshift = %g km/sec\n\n",3.59714 * wdm_vth);

  printf("\nWriting IC-file\n");
  const int64_t blockmaxlen = BUFFER * 1024 * 1024;
  if(!(block = (float *) malloc(blockmaxlen*3*sizeof(float))))
    {
      printf("failed to allocate memory for `block' (%ld MB).\n", 3*sizeof(float)*blockmaxlen/1024/1024);
      exit(24);
    }

  /*We are about to write the POS block*/
  /* Add displacement to Lagrangian coordinates, and multiply velocities by correct factor when writing VEL*/
  pc = 0;
  for(int64_t i = 0; i < NumPart; i++){
      for(int k = 0; k < 3; k++){
          block[3 * pc + k] = Pgrid.Pos(i,k, type) + outdata.Vel(i,k);
          if(twolpt)
                block[3 * pc + k] -= 3./7. * outdata.Vel2(i,k);
          block[3 * pc + k] = periodic_wrap(block[3*pc+k], Box);
      }
      pc++;

#ifdef NEUTRINO_PAIRS
      /*Add an extra copy of the position vector for the double neutrino*/
      if(type == NEUTRINO_TYPE) {
	  for(int k = 0; k < 3; k++)
	    block[3 * pc + k] = periodic_wrap(Pgrid.Pos(i,k,type) + outdata.Vel(i,k), Box);
	  pc++;
      }
#endif //NEUTRINO_PAIRS

      if(pc > blockmaxlen){
	  if(snap.WriteBlocks(posstr,type, block, pc,written) != pc)
                  FatalError(2);
          written+=pc;
	  pc = 0;
	}
  }
  if(pc > 0)
	  if(snap.WriteBlocks(posstr,type, block, pc,written) != pc)
                  FatalError(2);
  /*Done writing POS block*/
  written=0;
  pc = 0;

  /* write velocities: sizes are the same as for positions */
  for(int64_t i = 0; i < NumPart; i++) {
      for(int k = 0; k < 3; k++){
          block[3 * pc + k] = vel_prefac*outdata.Vel(i,k);
          if(twolpt)
                block[3 * pc + k] += vel_prefac2* outdata.Vel2(i,k);
      }

      //Add thermal velocities
      if(WDM_On == 1 && WDM_Vtherm_On == 1 && type == 1)
          WDMvels.add_thermal_speeds(&block[3 * pc]);
#ifdef NEUTRINOS
      if(NU_On == 1 && NU_Vtherm_On == 1 && type == 2) {
#ifdef NEUTRINO_PAIRS
          float vtherm[3];
          for(int k = 0; k < 3; k++)
              vtherm[k] = 0;
          nuvels.add_thermal_speeds(vtherm);
          for(int k = 0; k < 3; k++)
              block[3 * pc + k] = vel_prefac*outdata.Vel(i,k) + vtherm[k];
          pc++;
          for(int k = 0; k < 3; k++)
              block[3 * pc + k] = vel_prefac*outdata.Vel(i,k) - vtherm[k];
#else
          nuvels.add_thermal_speeds(&block[3 * pc]);
#endif //NEUTRINO_PAIRS
      }
#endif //NEUTRINOS
      pc++;

      if(pc > blockmaxlen){
          if(snap.WriteBlocks(velstr,type, block, pc,written) != pc)
                  FatalError(2);
          written+=pc;
          pc = 0;
      }
  }
  if(pc > 0)
	  if(snap.WriteBlocks(velstr,type, block, pc,written) != pc)
                  FatalError(2);

  /* write particle ID */
  written=0;
  pc = 0;
  blockid = (id_type *) block;
  for(int64_t i = 0; i < NumPart; i++) {
      blockid[pc] = i+FirstId;
      pc++;

#ifdef NEUTRINO_PAIRS 
      if(type == 2) {
	  blockid[pc] = i+FirstId+NumPart;
	  pc++;
	}
#endif //NEUTRINO_PAIRS

      if(pc > blockmaxlen){
	  if(snap.WriteBlocks(idstr,type, block, pc,written) != pc)
                  FatalError(2);
          written+=pc;
	  pc = 0;
      }
  }
  if(pc > 0)
	  if(snap.WriteBlocks(idstr,type, block, pc,written) != pc)
                  FatalError(2);

  /*Done writing IDs*/

  /* write zero temperatures if needed */
  if(type== 0) {
      written=0;
      pc = 0;
      for(int64_t i = 0; i < NumPart; i++) {
          block[pc] = 0;
          pc++;
          if(pc > blockmaxlen){
              if(snap.WriteBlocks(ustr,type, block, pc,written) != pc)
                  FatalError(2);
              written+=pc;
              pc = 0;
          }
      }
      if(pc > 0)
         if(snap.WriteBlocks(ustr,type, block, pc,written) != pc)
                  FatalError(2);
  }
  /*Done writing temperatures*/

  free(block);
  printf("Finished writing IC file.\n");
  FirstId+=NumPart;
#ifdef NEUTRINO_PAIRS
  if(type==2)
          FirstId+=NumPart;
#endif
  return FirstId;
}

double periodic_wrap(double x, double box)
{
  x = fmod(x, box);

  if (x < 0)
    x += box;

  return x;
}

