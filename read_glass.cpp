#include <math.h>
#include <stdlib.h>
#include <gadgetreader.hpp>
#include "allvars.h"
#include "proto.h"



void read_glass(char *fname)
{
  int i, j, k, n, m, count, type;
  float *pos = NULL;
  float x, y, z;
  size_t bytes;
  int num, numfiles, skip, nlocal;
  GadgetReader::GSnap snap(fname);
  GadgetReader::gadget_header header1=snap.GetHeader();

  Nglass = 0;
  for(k = 0; k < 6; k++)
        Nglass += header1.npartTotal[k];
  printf("Reading pre-IC file '%s'. Nglass=%d.\n", fname, Nglass);
  if(!(pos = (float *) malloc(sizeof(float) * Nglass * 3))){
          fprintf(stderr,"failed to allocate %g Mbyte for glass file\n", sizeof(float) * Nglass * 3.0 / (1024.0 * 1024.0));
		  FatalError(112);
  }
  /*Read all POS data*/
  if(snap.GetBlock("POS ",pos,Nglass,0,0) != Nglass){
          fprintf(stderr, "Error reading particle data\n");
          free(pos);
          FatalError(113);
  }

#if defined(MULTICOMPONENTGLASSFILE) && defined(DIFFERENT_TRANSFER_FUNC)
  MinType = 7;
  MaxType = -2;
  for(type = 0; type < 6; type++)
    if(header1.npartTotal[type])
      {
	if(MinType > type)
	  MinType = type;

	if(MaxType < type)
	  MaxType = type;
      }
#endif

  NumPart = Nglass*GlassTileFac*GlassTileFac*GlassTileFac;
  printf("\nTotal number of particles  = %ld\n",NumPart);

  if(!(P = (struct part_data *) malloc(bytes = sizeof(struct part_data) * NumPart))){
	printf("Failed to allocate %g Mbyte (%d particles)\n", bytes / (1024.0 * 1024.0),NumPart);
	FatalError(9891);
  }

  count = 0;

  IDStart = 1;

  for(i = 0; i < GlassTileFac; i++)
    for(j = 0; j < GlassTileFac; j++)
      for(k = 0; k < GlassTileFac; k++)
	{
	  for(type = 0, n = 0; type < 6; type++)
	    {
	      for(m = 0; m < header1.npartTotal[type]; m++, n++)
		{
		  x = pos[3 * n] / header1.BoxSize * (Box / GlassTileFac) + i * (Box / GlassTileFac);

		      y = pos[3 * n + 1] / header1.BoxSize * (Box / GlassTileFac) + j * (Box / GlassTileFac);
		      z = pos[3 * n + 2] / header1.BoxSize * (Box / GlassTileFac) + k * (Box / GlassTileFac);

		      P[count].Pos[0] = x;
		      P[count].Pos[1] = y;
		      P[count].Pos[2] = z;
#ifdef  MULTICOMPONENTGLASSFILE
		      P[count].Type = type;
#endif
		      P[count].ID = IDStart;

		      count++;

		  IDStart++;
		}
	    }
	}

  if(count != NumPart)
    {
      printf("fatal mismatch (%d %d)\n", count, NumPart);
      FatalError(1);
    }

  free(pos);
  return;
}

