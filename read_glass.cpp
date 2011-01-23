#include <math.h>
#include <stdlib.h>
#include <gadgetreader.hpp>
#include "allvars.h"
#include "proto.h"


int64_t read_glass(GadgetReader::GSnap& snap, int type, int GlassTileFac, struct part_data *& P)
{
  int i, j, k, n;
  float *pos = NULL;
  size_t bytes;
  int64_t NumPart;
  int64_t Nglass = snap.GetNpart(type);
  double GlassBox=snap.GetHeader().BoxSize;

  printf("Reading pre-IC file '%s'. Nglass=%ld.\n",snap.GetFileName().c_str(), Nglass);
  if(!(pos = (float *) malloc(sizeof(float) * Nglass * 3))){
          fprintf(stderr,"failed to allocate %g Mbyte for glass file\n", sizeof(float) * Nglass * 3.0 / (1024.0 * 1024.0));
		  FatalError(112);
  }
  /*Read all POS data for this type*/
  if(snap.GetBlock("POS ",pos,Nglass,0,(1<<N_TYPE)-1-(1<<type)) != Nglass){
          fprintf(stderr, "Error reading particle data\n");
          free(pos);
          exit(113);
  }

  NumPart = Nglass*GlassTileFac*GlassTileFac*GlassTileFac;

  if(!(P = (struct part_data *) malloc(bytes = sizeof(struct part_data) * NumPart))){
	printf("Failed to allocate %g Mbyte (%ld particles)\n", bytes / (1024.0 * 1024.0),NumPart);
	exit(9891);
  }

  for(i = 0; i < GlassTileFac; i++)
    for(j = 0; j < GlassTileFac; j++)
      for(k = 0; k < GlassTileFac; k++){
        for(n = 0; n < Nglass; n++){
            P[n].Pos[0] = pos[3 * n] / GlassBox * (Box / GlassTileFac) + i * (Box / GlassTileFac);
            P[n].Pos[1]= pos[3 * n + 1] / GlassBox * (Box / GlassTileFac) + j * (Box / GlassTileFac);
            P[n].Pos[2]  = pos[3 * n + 2] / GlassBox * (Box / GlassTileFac) + k * (Box / GlassTileFac);
        }
      }

  free(pos);
  printf("Type %d has %ld particles\n",type,NumPart);
  return NumPart;
}

