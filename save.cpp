#include "allvars.h"
#include "proto.h"
#include <gadgetwriter.hpp>

using namespace GadgetWriter;
using namespace std;

#ifdef NO64BITID
        typedef int32_t id_type;
#else
        typedef int64_t id_type;
#endif //NO64BITID

void write_particle_data(int type)
{
  size_t bytes;
  float *block;
  id_type *blockid;
  int blockmaxlen;
  int64_t maxidlen,written=0;
  std::valarray<int64_t> npart(N_TYPE);
  int i,k,pc;
    
  printf("\nWriting IC-file\n");
  
  #define BUFFER 10
  /*Set particle numbers*/
  for(i = 0; i < N_TYPE; i++)
    npart[i] = header1.npartTotal[i] * GlassTileFac * GlassTileFac * GlassTileFac;

#ifdef NEUTRINO_PAIRS
  npart[NEUTRINO_TYPE] *= 2;
#endif //NEUTRINO_PAIRS

  GWriteSnap snap(string(OutputDir)+string("/")+string(FileBase), npart,NumFiles, sizeof(id_type));

  /*Write header*/
  snap.WriteHeaders(generate_header());

  if(!(block = (float *) malloc(bytes = BUFFER * 1024 * 1024)))
    {
      printf("failed to allocate memory for `block' (%g bytes).\n", (double) bytes);
      exit(24);
    }

  blockmaxlen = bytes / (3 * sizeof(float));

  /*We are about to write the POS block*/
  for(i = 0, pc = 0; i < NumPart; i++){
      for(k = 0; k < 3; k++)
	  block[3 * pc + k] = P[i].Pos[k];
      pc++;

#ifdef NEUTRINO_PAIRS
      /*Add an extra copy of the position vector for the double neutrino*/
      if(type == NEUTRINO_TYPE) {
	  for(k = 0; k < 3; k++)
	    block[3 * pc + k] = P[i].Pos[k];
	  pc++;
      }
#endif //NEUTRINO_PAIRS

      if(pc > blockmaxlen){
	  snap.WriteBlocks(string("POS "),type, block, pc,written);
          written+=pc;
	  pc = 0;
	}
  }
  if(pc > 0)
	  snap.WriteBlocks(string("POS "),type, block, pc,written);
  /*Done writing POS block*/
  written=0;

  /* write velocities: sizes are the same as for positions */
  for(i = 0, pc = 0; i < NumPart; i++)
    {
      for(k = 0; k < 3; k++)
	block[3 * pc + k] = P[i].Vel[k];

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
	    block[3 * pc + k] = P[i].Vel[k] + vtherm[k];
	  pc++;
	  for(k = 0; k < 3; k++)
	    block[3 * pc + k] = P[i].Vel[k] - vtherm[k];
	}
#else
      if(NU_On == 1 && NU_Vtherm_On == 1 && type == 2)
	add_NU_thermal_speeds(&block[3 * pc]);
#endif //NEUTRINO_PAIRS

#endif //NEUTRINOS
      pc++;

      if(pc > blockmaxlen){
	  snap.WriteBlocks(string("VEL "),type, block, pc,written);
          written+=pc;
	  pc = 0;
	}
    }
  if(pc > 0)
	  snap.WriteBlocks(string("VEL "),type, block, pc,written);

  /* write particle ID */
  written=0;
  blockid = (id_type *) block;
  maxidlen = bytes / (sizeof(id_type));
  
  for(i = 0, pc = 0; i < NumPart; i++){
      blockid[pc] = i;
      pc++;

#ifdef NEUTRINO_PAIRS /* Extra neutrinos get an ID after all other particles*/
      if(type == 2)
	{
	  blockid[pc] += header.npartTotal[0] + header.npartTotal[1] + header.npartTotal[2] +
	    header.npartTotal[3] + header.npartTotal[4] + header.npartTotal[5];
	  pc++;
	}
#endif //NEUTRINO_PAIRS

      if(pc > blockmaxlen){
	  snap.WriteBlocks(string("ID  "),type, block, pc,written);
          written+=pc;
	  pc = 0;
      }
  }
  if(pc > 0)
	  snap.WriteBlocks(string("ID  "),type, block, pc,written);

  /*Done writing IDs*/

  /* write zero temperatures if needed */
  if(npart[0]) {
          written=0;

      for(i = 0, pc = 0; i < npart[0]; i++){
	  block[pc] = 0;
	  pc++;
          if(pc > blockmaxlen){
              snap.WriteBlocks(string("U   "),type, block, pc,written);
              written+=pc;
              pc = 0;
          }
      }
      if(pc > 0)
	  snap.WriteBlocks(string("U   "),type, block, pc,written);
  }
  /*Done writing temperatures*/

  free(block);
  printf("Finished writing IC file.\n");
  return;
}

gadget_header generate_header()
{
  gadget_header header;
  double scale = 3 * Hubble * Hubble / (8 * M_PI * G) * pow(Box,3);
  /*Set masses*/
  for(int i = 0; i < N_TYPE; i++)
      header.mass[i] = 0;

  if(header1.npartTotal[BARYON_TYPE])
    header.mass[BARYON_TYPE] = (OmegaBaryon) * scale / (header.npartTotal[BARYON_TYPE]);

  if(header1.npartTotal[DM_TYPE])
    header.mass[DM_TYPE] = (Omega - OmegaBaryon - OmegaDM_2ndSpecies) * scale / (header.npartTotal[DM_TYPE]);

  if(header1.npartTotal[NEUTRINO_TYPE]){
    header.mass[NEUTRINO_TYPE] = (OmegaDM_2ndSpecies) * scale / (header.npartTotal[NEUTRINO_TYPE]);
#ifdef NEUTRINO_PAIRS
    header.mass[NEUTRINO_TYPE] /= 2;
#endif //NEUTRINO_PAIRS
  }

  header.time = InitTime;
  header.redshift = 1.0 / InitTime - 1;

  header.BoxSize = Box;
  header.Omega0 = Omega;
  header.OmegaLambda = OmegaLambda;
  header.HubbleParam = HubbleParam;
  /*Various flags; Most set by gadget later*/
  header.flag_sfr = 0;
  header.flag_feedback = 0;
  header.flag_cooling = 0;
  header.flag_stellarage = 0;
  header.flag_metals = 0;
  header.flag_entropy_instead_u=0;
  header.flag_doubleprecision=0;
  header.flag_ic_info=1;        
  header.lpt_scalingfactor=1;  
  return header;
}

