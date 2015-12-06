#include "proto.h"
#include "part_data.hpp"
#include "thermalvel.hpp"
#include <exception>
#include <gadgetwriter.hpp>

using namespace GadgetWriter;
using namespace std;

#define BUFFER 48

int64_t write_particle_data(GWriteSnap & snap, int type, lpt_data& outdata, part_grid& Pgrid, FermiDiracVel *therm_vels, const double vel_prefac, const double vel_prefac2, int64_t FirstId, const bool twolpt)
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
  if(snap.GetFormat() == 3){
      posstr = "Coordinates";
      velstr = "Velocities";
      idstr = "ParticleIDs";
      ustr = "InternalEnergy";
  }
  printf("\nWriting IC-file\n");
  const int64_t blockmaxlen = BUFFER * 1024 * 1024;
  if(!(block = (float *) malloc(blockmaxlen*3*sizeof(float))))
    {
      fprintf(stderr, "failed to allocate memory for `block' (%ld MB).\n", 3*sizeof(float)*blockmaxlen/1024/1024);
      throw std::bad_alloc();
    }

  /*We are about to write the POS block*/
  /* Add displacement to Lagrangian coordinates, and multiply velocities by correct factor when writing VEL*/
  pc = 0;
  for(int64_t i = 0; i < NumPart; i++){
      for(int k = 0; k < 3; k++){
          block[3 * pc + k] = Pgrid.Pos(i,k, type) + outdata.Vel(i,k);
          if(twolpt)
                block[3 * pc + k] -= 3./7. * outdata.Vel2(i,k);
          block[3 * pc + k] = periodic_wrap(block[3*pc+k], Pgrid.GetBox());
      }
      pc++;

#ifdef NEUTRINO_PAIRS
      /*Add an extra copy of the position vector for the double neutrino*/
      if(type == NEUTRINO_TYPE) {
	  for(int k = 0; k < 3; k++)
	    block[3 * pc + k] = periodic_wrap(Pgrid.Pos(i,k,type) + outdata.Vel(i,k), Pgrid.GetBox());
	  pc++;
      }
#endif //NEUTRINO_PAIRS
      if(pc > blockmaxlen){
	    if(snap.WriteBlocks(posstr,type, block, pc,written) != pc) {
            fprintf(stderr, "Could not write position data at particle %ld", i);
            exit(2);
        }
        written+=pc;
	    pc = 0;
	  }
  }
  if(pc > 0 && snap.WriteBlocks(posstr,type, block, pc,written) != pc) {
        fprintf(stderr, "Could not write final position data");
        exit(2);
  }
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
      if(therm_vels) {
#ifdef NEUTRINO_PAIRS
          float vtherm[3];
          for(int k = 0; k < 3; k++)
              vtherm[k] = 0;
          therm_vels->add_thermal_speeds(vtherm);
          for(int k = 0; k < 3; k++)
              block[3 * pc + k] = vel_prefac*outdata.Vel(i,k) + vtherm[k];
          pc++;
          for(int k = 0; k < 3; k++)
              block[3 * pc + k] = vel_prefac*outdata.Vel(i,k) - vtherm[k];
#else
          therm_vels->add_thermal_speeds(&block[3 * pc]);
#endif //NEUTRINO_PAIRS
      }
      pc++;

      if(pc > blockmaxlen){
          if(snap.WriteBlocks(velstr,type, block, pc,written) != pc) {
            fprintf(stderr, "Could not write velocity data at particle %ld", i);
            exit(2);
          }
          written+=pc;
          pc = 0;
      }
  }
  if(pc > 0 && snap.WriteBlocks(velstr,type, block, pc,written) != pc) {
    fprintf(stderr, "Could not write final velocity data");
    exit(2);
  }

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
          if(snap.WriteBlocks(idstr,type, block, pc,written) != pc) {
                fprintf(stderr, "Could not write id data at particle %ld", i);
                exit(2);
          }
          written+=pc;
	      pc = 0;
      }
  }
  if(pc > 0)
	  if(snap.WriteBlocks(idstr,type, block, pc,written) != pc) {
        fprintf(stderr, "Could not write final id data");
        exit(2);
      }

  /*Done writing IDs*/
  /* write zero temperatures if needed */
  if(type== 0) {
      written=0;
      pc = 0;
      for(int64_t i = 0; i < NumPart; i++) {
          block[pc] = 0;
          pc++;
          if(pc > blockmaxlen){
              if(snap.WriteBlocks(ustr,type, block, pc,written) != pc) {
                fprintf(stderr, "Could not write zero temp data at particle %ld", i);
                exit(2);
              }
              written+=pc;
              pc = 0;
          }
      }
      if(pc > 0)
         if(snap.WriteBlocks(ustr,type, block, pc,written) != pc) {
            fprintf(stderr, "Could not write final temperature data");
            exit(2);
         }
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

