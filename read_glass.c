#include <math.h>
#include <stdlib.h>
#include <mpi.h>
#include "allvars.h"
#include "proto.h"




void read_glass(char *fname)
{
  int i, j, k, n, m, slab, count, type;
  unsigned int dummy, dummy2;
  float *pos = 0;
  float x, y, z;
  FILE *fd = 0;
  size_t bytes;
  int *npart_Task;
  int num, numfiles, skip, nlocal;
  char buf[500];

#define SKIP {my_fread(&dummy, sizeof(int), 1, fd);}
#define SKIP2 {my_fread(&dummy2, sizeof(int), 1, fd);}

  if(ThisTask == 0)
    {
      printf("\nreading Lagrangian glass file...\n");
      fflush(stdout);

      numfiles = find_files(fname);

      for(num = 0, skip = 0; num < numfiles; num++)
	{
	  if(numfiles > 1)
	    sprintf(buf, "%s.%d", fname, num);
	  else
	    sprintf(buf, "%s", fname);

	  if(!(fd = fopen(buf, "r")))
	    {
	      printf("can't open file `%s' for reading glass file.\n", buf);
	      FatalError(1);
	    }

	  SKIP;
	  my_fread(&header1, sizeof(header1), 1, fd);
	  SKIP2;

	  if(dummy != sizeof(header1) || dummy2 != sizeof(header1))
	    {
	      printf("incorrect header size!\n");
	      FatalError(2);
	    }

	  nlocal = 0;

	  for(k = 0; k < 6; k++)
	    nlocal += header1.npart[k];

	  printf("reading '%s' with %d particles\n", fname, nlocal);

	  if(num == 0)
	    {
	      Nglass = 0;

	      for(k = 0; k < 6; k++)
		Nglass += header1.npartTotal[k];

	      printf("\nNglass= %d\n\n", Nglass);
	      pos = (float *) malloc(sizeof(float) * Nglass * 3);

	      if(!(pos))
		{
		  printf("failed to allocate %g Mbyte on Task %d for glass file\n",
			 sizeof(float) * Nglass * 3.0 / (1024.0 * 1024.0), ThisTask);
		  FatalError(112);
		}
	    }

	  SKIP;
	  my_fread(&pos[3 * skip], sizeof(float), 3 * nlocal, fd);
	  SKIP2;

	  if(dummy != sizeof(float) * 3 * nlocal || dummy2 != sizeof(float) * 3 * nlocal)
	    {
	      printf("incorrect block structure in positions block!\n");
	      FatalError(3);
	    }
	  skip += nlocal;

	  fclose(fd);
	}
    }

  MPI_Bcast(&Nglass, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&header1, sizeof(header1), MPI_BYTE, 0, MPI_COMM_WORLD);

  if(ThisTask != 0)
    {
      pos = (float *) malloc(sizeof(float) * Nglass * 3);

      if(!(pos))
	{
	  printf("failed to allocate %g Mbyte on Task %d for glass file\n",
		 sizeof(float) * Nglass * 3.0 / (1024.0 * 1024.0), ThisTask);
	  FatalError(112);
	}
    }

  MPI_Bcast(&pos[0], sizeof(float) * Nglass * 3, MPI_BYTE, 0, MPI_COMM_WORLD);


  npart_Task = malloc(sizeof(int) * NTask);

  for(i = 0; i < NTask; i++)
    npart_Task[i] = 0;

#if defined(MULTICOMPONENTGLASSFILE) && defined(DIFFERENT_TRANSFER_FUNC)
  MinType = 7;
  MaxType = -2;
  for(type = 0; type < 6; type++)
    if(header1.npartTotal[type])
      {
	if(MinType > type - 1)
	  MinType = type - 1;

	if(MaxType < type - 1)
	  MaxType = type - 1;
      }
#endif



  for(i = 0; i < GlassTileFac; i++)
    for(j = 0; j < GlassTileFac; j++)
      for(k = 0; k < GlassTileFac; k++)
	{
	  for(type = 0, n = 0; type < 6; type++)
	    {
	      for(m = 0; m < header1.npartTotal[type]; m++, n++)
		{
		  x = pos[3 * n] / header1.BoxSize * (Box / GlassTileFac) + i * (Box / GlassTileFac);

		  slab = x / Box * Nmesh;
		  if(slab >= Nmesh)
		    slab = Nmesh - 1;

		  npart_Task[Slab_to_task[slab]] += 1;
		}
	    }
	}

  TotNumPart = 0;		/* note: This is a 64 bit integer */
  NTaskWithN = 0;

  NumPart = npart_Task[ThisTask];

  for(i = 0; i < NTask; i++)
    {
      TotNumPart += npart_Task[i];
      if(npart_Task[i] > 0)
	NTaskWithN++;
    }


  if(ThisTask == 0)
    {
      for(i = 0; i < NTask; i++)
	printf("%d particles on task=%d  (slabs=%d)\n", npart_Task[i], i, Local_nx_table[i]);

      printf("\nTotal number of particles  = %d%09d\n\n",
	     (int) (TotNumPart / 1000000000), (int) (TotNumPart % 1000000000));

      fflush(stdout);
    }


  free(npart_Task);


  if(NumPart)
    {
      P = (struct part_data *) malloc(bytes = sizeof(struct part_data) * NumPart);

      if(!(P))
	{
	  printf("failed to allocate %g Mbyte (%d particles) on Task %d\n", bytes / (1024.0 * 1024.0),
		 NumPart, ThisTask);
	  FatalError(9891);
	}
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

		  slab = x / Box * Nmesh;
		  if(slab >= Nmesh)
		    slab = Nmesh - 1;

		  if(Slab_to_task[slab] == ThisTask)
		    {
		      y = pos[3 * n + 1] / header1.BoxSize * (Box / GlassTileFac) + j * (Box / GlassTileFac);
		      z = pos[3 * n + 2] / header1.BoxSize * (Box / GlassTileFac) + k * (Box / GlassTileFac);

		      P[count].Pos[0] = x;
		      P[count].Pos[1] = y;
		      P[count].Pos[2] = z;
#ifdef  MULTICOMPONENTGLASSFILE
		      P[count].Type = type - 1;
#endif
		      P[count].ID = IDStart;

		      count++;
		    }

		  IDStart++;
		}
	    }
	}

  if(count != NumPart)
    {
      printf("fatal mismatch (%d %d) on Task %d\n", count, NumPart, ThisTask);
      FatalError(1);
    }

  free(pos);
}


int find_files(char *fname)
{
  FILE *fd;
  char buf[200], buf1[200];
  int dummy;

  sprintf(buf, "%s.%d", fname, 0);
  sprintf(buf1, "%s", fname);

  if((fd = fopen(buf, "r")))
    {
      my_fread(&dummy, sizeof(dummy), 1, fd);
      my_fread(&header, sizeof(header), 1, fd);
      my_fread(&dummy, sizeof(dummy), 1, fd);

      fclose(fd);

      return header.num_files;
    }

  if((fd = fopen(buf1, "r")))
    {
      my_fread(&dummy, sizeof(dummy), 1, fd);
      my_fread(&header, sizeof(header), 1, fd);
      my_fread(&dummy, sizeof(dummy), 1, fd);

      fclose(fd);
      header.num_files = 1;

      return header.num_files;
    }

  FatalError(121);
  return 0;
}
