#include <math.h>
#include <stdlib.h>
#include <gadgetreader.hpp>
#include "allvars.h"
#include "proto.h"


int64_t read_glass(GadgetReader::GSnap& snap, int type, int GlassTileFac, struct part_data *& P)
{
  float *pos = NULL;
  size_t bytes;
  int64_t NumPart;
  const int64_t Nglass = snap.GetNpart(type);
  const double GlassBox=snap.GetHeader().BoxSize;
  const int G2Tile = GlassTileFac*GlassTileFac;
  const int G3Tile = G2Tile*GlassTileFac;

  if(!(pos = (float *) malloc(sizeof(float) * Nglass * 3))){
          fprintf(stderr,"failed to allocate %ld Mbyte for glass file\n", sizeof(float) * Nglass * 3 / (1024 * 1024));
		  FatalError(112);
  }
  /*Read all POS data for this type*/
  if(snap.GetBlock("POS ",pos,Nglass,0,(1<<N_TYPE)-1-(1<<type)) != Nglass){
          fprintf(stderr, "Error reading particle data\n");
          free(pos);
          exit(113);
  }

  NumPart = Nglass*G3Tile;
  bytes = sizeof(struct part_data) * NumPart;
  if(!(P = (struct part_data *) malloc(bytes))){
	printf("Failed to allocate %ld Mbyte (%ld particles)\n", bytes / (1024 * 1024),NumPart);
	exit(9891);
  }
  else
          printf("Type %d allocated %ld MB for %ld particles\n",type,bytes/1024/1024, NumPart);
  
#pragma omp parallel
  {
  
  #pragma omp for schedule(static)
  for(int i = 0; i < GlassTileFac; i++)
    for(int j = 0; j < GlassTileFac; j++)
      for(int k = 0; k < GlassTileFac; k++){
        for(int n = 0; n < Nglass; n++){
            size_t index = n+k*Nglass+j*Nglass*GlassTileFac+i*G2Tile*Nglass;
            P[index].Pos[0] = pos[3 * n] / GlassBox * (Box / GlassTileFac) + i * (Box / GlassTileFac);
            P[index].Pos[1]= pos[3 * n + 1] / GlassBox * (Box / GlassTileFac) + j * (Box / GlassTileFac);
            P[index].Pos[2]  = pos[3 * n + 2] / GlassBox * (Box / GlassTileFac) + k * (Box / GlassTileFac);
        }
      }
  }

  free(pos);
  return NumPart;
}

