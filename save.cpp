#include "allvars.h"
#include "proto.h"
#include <gadgetwriter.hpp>

using namespace GadgetWriter;
using namespace std;

#define BUFFER 48

int64_t write_particle_data(GWriteSnap & snap, int type, part_data& P, int64_t NumPart, int64_t FirstId)
{
  size_t bytes;
  float *block;
  id_type *blockid;
  int blockmaxlen;
  int64_t written=0;
  int i,k,pc;
  const double hubble_a = Hubble * sqrt(Omega / pow(InitTime, 3) + (1 - Omega - OmegaLambda) / pow(InitTime, 2) + OmegaLambda);
  const double vel_prefac = InitTime * hubble_a * F_Omega(InitTime) /sqrt(InitTime);
#ifdef TWOLPT
  const double vel_prefac2 = InitTime * hubble_a * F2_Omega(InitTime) /sqrt(InitTime);
#endif
  printf("vel_prefac= %g  hubble_a=%g fom=%g Omega=%g \n", vel_prefac, hubble_a, F_Omega(InitTime), Omega);

    
  printf("\nWriting IC-file\n");

  if(!(block = (float *) malloc(bytes = BUFFER * 1024 * 1024)))
    {
      printf("failed to allocate memory for `block' (%g bytes).\n", static_cast<double>(bytes));
      exit(24);
    }

  blockmaxlen = bytes / (3 * sizeof(float));

  /*We are about to write the POS block*/
  /* Add displacement to Lagrangian coordinates, and multiply velocities by correct factor when writing VEL*/
  for(i = 0, pc = 0; i < NumPart; i++){
      for(k = 0; k < 3; k++){
	  block[3 * pc + k] = P.Pos(i,k) + P.Vel(i,k);
#ifdef TWOLPT
#ifdef NEUTRINOS
    if(type !=2)
#endif
	  block[3 * pc + k] -= 3./7. * P.Vel2(i,k);
#endif
	  block[3 * pc + k] = periodic_wrap(block[3*pc+k]);
      }
      pc++;

#ifdef NEUTRINO_PAIRS
      /*Add an extra copy of the position vector for the double neutrino*/
      if(type == NEUTRINO_TYPE) {
	  for(k = 0; k < 3; k++)
	    block[3 * pc + k] = periodic_wrap(P.Pos(i,k) + P.Vel(i,k));
	  pc++;
      }
#endif //NEUTRINO_PAIRS

      if(pc > blockmaxlen){
	  if(snap.WriteBlocks(string("POS "),type, block, pc,written) != pc)
                  FatalError(2);
          written+=pc;
	  pc = 0;
	}
  }
  if(pc > 0)
	  if(snap.WriteBlocks(string("POS "),type, block, pc,written) != pc)
                  FatalError(2);
  /*Done writing POS block*/
  written=0;

  /* write velocities: sizes are the same as for positions */
  for(i = 0, pc = 0; i < NumPart; i++)
    {
      for(k = 0; k < 3; k++){
	block[3 * pc + k] = vel_prefac*P.Vel(i,k);
#ifdef TWOLPT
#ifdef NEUTRINOS
    if(type !=2)
#endif
        block[3 * pc + k] -= 3./7. *vel_prefac2* P.Vel2(i,k);
#endif
      }

      if(WDM_On == 1 && WDM_Vtherm_On == 1 && type == 1)
	add_WDM_thermal_speeds(&block[3 * pc]);
#ifdef NEUTRINOS

#ifdef NEUTRINO_PAIRS
      if(NU_On == 1 && NU_Vtherm_On == 1 && type == 2) {
	  float vtherm[3];
	  for(k = 0; k < 3; k++)
	    vtherm[k] = 0;
	  add_NU_thermal_speeds(vtherm);
	  for(k = 0; k < 3; k++)
	    block[3 * pc + k] = vel_prefac*P.Vel(i,k) + vtherm[k];
	  pc++;
	  for(k = 0; k < 3; k++)
	    block[3 * pc + k] = vel_prefac*P.Vel(i,k) - vtherm[k];
	}
#else
      if(NU_On == 1 && NU_Vtherm_On == 1 && type == 2)
	add_NU_thermal_speeds(&block[3 * pc]);
#endif //NEUTRINO_PAIRS

#endif //NEUTRINOS
      pc++;

      if(pc > blockmaxlen){
	  if(snap.WriteBlocks(string("VEL "),type, block, pc,written) != pc)
                  FatalError(2);
          written+=pc;
	  pc = 0;
	}
    }
  if(pc > 0)
	  if(snap.WriteBlocks(string("VEL "),type, block, pc,written) != pc)
                  FatalError(2);

  /* write particle ID */
  written=0;
  blockid = (id_type *) block;
  
  for(i = 0, pc = 0; i < NumPart; i++){
      blockid[pc] = i+FirstId;
      pc++;

#ifdef NEUTRINO_PAIRS 
      if(type == 2) {
	  blockid[pc] = i+FirstId+NumPart;
	  pc++;
	}
#endif //NEUTRINO_PAIRS

      if(pc > blockmaxlen){
	  if(snap.WriteBlocks(string("ID  "),type, block, pc,written) != pc)
                  FatalError(2);
          written+=pc;
	  pc = 0;
      }
  }
  if(pc > 0)
	  if(snap.WriteBlocks(string("ID  "),type, block, pc,written) != pc)
                  FatalError(2);

  /*Done writing IDs*/

  /* write zero temperatures if needed */
  if(type== 0) {
          written=0;

      for(i = 0, pc = 0; i < NumPart; i++){
	  block[pc] = 0;
	  pc++;
          if(pc > blockmaxlen){
	  if(snap.WriteBlocks(string("U   "),type, block, pc,written) != pc)
                  FatalError(2);
              written+=pc;
              pc = 0;
          }
      }
      if(pc > 0)
	  if(snap.WriteBlocks(string("U   "),type, block, pc,written) != pc)
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

double periodic_wrap(double x)
{
  x = fmod(x, Box);

  if (x < 0)
    x += Box;

  return x;
}

