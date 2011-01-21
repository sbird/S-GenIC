#include <math.h>
#include <stdlib.h>

#include "allvars.h"
#include "proto.h"

void write_particle_data(void)
{
    printf("\nWriting IC-file\n");
    
#define BUFFER 10
  size_t bytes;
  float *block;
#ifdef NO64BITID
  int32_t *blockid;
#else
  int64_t *blockid;
#endif
  int blkheadsize = sizeof(int) + 4 * sizeof(char);
  int nextblock;
  int blockmaxlen; 
  int64_t maxidlen;
  int4byte dummy;
  FILE *fd;
  char buf[300];
  int64_t i, k, pc;


  if(NumPart == 0)
    return;

    sprintf(buf, "%s/%s", OutputDir, FileBase);

  if(!(fd = fopen(buf, "w"))){
      printf("Error. Can't write in file '%s'\n", buf);
      exit(10);
  }

  for(i = 0; i < 6; i++)
    {
      header.npart[i] = 0;
      header.npartTotal[i] = 0;
      header.mass[i] = 0;
    }

  qsort(P, NumPart, sizeof(struct part_data), compare_type);	/* sort particles by type, because that's how they should be stored in a gadget binary file */

  for(i = 0; i < 3; i++)
    header.npartTotal[i] = header1.npartTotal[i] * GlassTileFac * GlassTileFac * GlassTileFac;

  for(i = 0; i < NumPart; i++)
    header.npart[P[i].Type]++;

  if(header.npartTotal[0])
    header.mass[0] =
      (OmegaBaryon) * 3 * Hubble * Hubble / (8 * PI * G) * pow(Box, 3) / (header.npartTotal[0]);

  if(header.npartTotal[1])
    header.mass[1] =
      (Omega - OmegaBaryon - OmegaDM_2ndSpecies) * 3 * Hubble * Hubble / (8 * PI * G) * pow(Box,3) /
      (header.npartTotal[1]);

  if(header.npartTotal[2])
    header.mass[2] =
      (OmegaDM_2ndSpecies) * 3 * Hubble * Hubble / (8 * PI * G) * pow(Box, 3) / (header.npartTotal[2]);


#ifdef NEUTRINO_PAIRS
  header.npart[2] *= 2;
  header.npartTotal[2] *= 2;
  header.mass[2] /= 2;
#endif //NEUTRINO_PAIRS

  header.time = InitTime;
  header.redshift = 1.0 / InitTime - 1;

  header.flag_sfr = 0;
  header.flag_feedback = 0;
  header.flag_cooling = 0;
  header.flag_stellarage = 0;
  header.flag_metals = 0;

  /*FIXME*/
  header.num_files = 1;

  header.BoxSize = Box;
  header.Omega0 = Omega;
  header.OmegaLambda = OmegaLambda;
  header.HubbleParam = HubbleParam;

  header.flag_stellarage = 0;
  header.flag_metals = 0;
  header.flag_entropy_instead_u=0;	/*!< flags that IC-file contains entropy instead of u */
  header.flag_doubleprecision=0;	/*!< flags that snapshot contains double-precision instead of single precision */

  header.flag_ic_info=1;             /*!< flag to inform whether IC files are generated with ordinary Zeldovich approximation,*/
  header.lpt_scalingfactor=1;      /*!< scaling factor for 2lpt initial conditions */
  dummy = sizeof(header);
#ifdef FORMAT_TWO
      /*Write format 2 header header*/
      my_fwrite(&blkheadsize,sizeof(dummy),1,fd);
      my_fwrite("HEAD", sizeof(char), 4, fd);
      nextblock = dummy + 2 * sizeof(int);
      my_fwrite(&nextblock, sizeof(int), 1, fd);
      my_fwrite(&blkheadsize,sizeof(dummy),1,fd);
#endif //FORMAT_TWO
  my_fwrite(&dummy, sizeof(dummy), 1, fd);
  my_fwrite(&header, sizeof(header), 1, fd);
  my_fwrite(&dummy, sizeof(dummy), 1, fd);



  if(!(block = malloc(bytes = BUFFER * 1024 * 1024)))
    {
      printf("failed to allocate memory for `block' (%g bytes).\n", (double) bytes);
      exit(24);
    }

  blockmaxlen = bytes / (3 * sizeof(float));

#ifdef NO64BITID
  blockid = (int32_t *) block;
  maxidlen = bytes / (sizeof(int32_t));
#else
  blockid = (int64_t *) block;
  maxidlen = bytes / (sizeof(int64_t));
#endif //NO64BITID

  /* write coordinates */
  dummy = sizeof(float) * 3 * NumPart;
#ifdef NEUTRINO_PAIRS
  dummy =
    sizeof(float) * 3 * (header.npart[0] + header.npart[1] + header.npart[2] + header.npart[3] +
			 header.npart[4] + header.npart[5]);
#endif //NEUTRINO_PAIRS
#ifdef FORMAT_TWO
          /*Write position header*/
	  blkheadsize = sizeof(int) + 4 * sizeof(char);
      	  my_fwrite(&blkheadsize,sizeof(int),1,fd);
      	  my_fwrite("POS ", sizeof(char), 4, fd);
	  nextblock=dummy+2*sizeof(int);
	  my_fwrite(&nextblock, sizeof(int), 1, fd);
	  my_fwrite(&blkheadsize, sizeof(int), 1, fd);
	  /*Done writing position header*/
#endif //FORMAT_TWO
  my_fwrite(&dummy, sizeof(dummy), 1, fd);
  for(i = 0, pc = 0; i < NumPart; i++)
    {
      for(k = 0; k < 3; k++)
	{
	  block[3 * pc + k] = P[i].Pos[k];
	}
      pc++;

#ifdef NEUTRINO_PAIRS
      if(P[i].Type == 2)
	{
	  for(k = 0; k < 3; k++)
	    block[3 * pc + k] = P[i].Pos[k];
	  pc++;
	}
#endif //NEUTRINO_PAIRS

      if(pc >= (blockmaxlen - 1))
	{
	  my_fwrite(block, sizeof(float), 3 * pc, fd);
	  pc = 0;
	}
    }
  if(pc > 0)
    my_fwrite(block, sizeof(float), 3 * pc, fd);
  my_fwrite(&dummy, sizeof(dummy), 1, fd);



  /* write velocities */
  dummy = sizeof(float) * 3 * NumPart;
#ifdef NEUTRINO_PAIRS
  dummy =
    sizeof(float) * 3 * (header.npart[0] + header.npart[1] + header.npart[2] + header.npart[3] +
			 header.npart[4] + header.npart[5]);
#endif //NEUTRINO_PAIRS
#ifdef FORMAT_TWO
          /*Write velocity header*/
	  blkheadsize = sizeof(int) + 4 * sizeof(char);
      	  my_fwrite(&blkheadsize,sizeof(int),1,fd);
      	  my_fwrite("VEL ", sizeof(char), 4, fd);
	  nextblock=dummy+2*sizeof(int);
	  my_fwrite(&nextblock, sizeof(int), 1, fd);
	  my_fwrite(&blkheadsize, sizeof(int), 1, fd);
	  /*Done writing velocity header*/
#endif //FORMAT_TWO
  my_fwrite(&dummy, sizeof(dummy), 1, fd);
  for(i = 0, pc = 0; i < NumPart; i++)
    {
      for(k = 0; k < 3; k++)
	block[3 * pc + k] = P[i].Vel[k];

      if(WDM_On == 1 && WDM_Vtherm_On == 1 && P[i].Type == 1)
	add_WDM_thermal_speeds(&block[3 * pc]);
#ifdef NEUTRINOS

#ifdef NEUTRINO_PAIRS
      if(NU_On == 1 && NU_Vtherm_On == 1 && P[i].Type == 2)
	{
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
      if(NU_On == 1 && NU_Vtherm_On == 1 && P[i].Type == 2)
	add_NU_thermal_speeds(&block[3 * pc]);
#endif //NEUTRINO_PAIRS

#endif //NEUTRINOS

      pc++;

      if(pc >= (blockmaxlen - 1))
	{
	  my_fwrite(block, sizeof(float), 3 * pc, fd);
	  pc = 0;
	}
    }
  if(pc > 0)
    my_fwrite(block, sizeof(float), 3 * pc, fd);
  my_fwrite(&dummy, sizeof(dummy), 1, fd);


  /* write particle ID */
#ifdef NO64BITID
  dummy = sizeof(int32_t) * NumPart;
#else
  dummy = sizeof(int64_t) * NumPart;
#endif //NO64BITID
#ifdef NEUTRINO_PAIRS
  dummy =
    sizeof(int) * (header.npart[0] + header.npart[1] + header.npart[2] + header.npart[3] + header.npart[4] +
		   header.npart[5]);
#endif //NEUTRINO_PAIRS
#ifdef FORMAT_TWO
          /*Write ID header*/
	  blkheadsize = sizeof(int) + 4 * sizeof(char);
      	  my_fwrite(&blkheadsize,sizeof(int),1,fd);
      	  my_fwrite("ID  ", sizeof(char), 4, fd);
	  nextblock=dummy+2*sizeof(int);
	  my_fwrite(&nextblock, sizeof(int), 1, fd);
	  my_fwrite(&blkheadsize, sizeof(int), 1, fd);
	  /*Done writing ID header*/
#endif //FORMAT_TWO
  my_fwrite(&dummy, sizeof(dummy), 1, fd);
  for(i = 0, pc = 0; i < NumPart; i++)
    {
      blockid[pc] = i;
      pc++;

#ifdef NEUTRINO_PAIRS
      if(P[i].Type == 2)
	{
	  blockid[pc] += header.npartTotal[0] + header.npartTotal[1] + header.npartTotal[2] +
	    header.npartTotal[3] + header.npartTotal[4] + header.npartTotal[5];
	  pc++;
	}
#endif //NEUTRINO_PAIRS

      if(pc >= (maxidlen - 1))
	{
#ifdef NO64BITID
	  my_fwrite(blockid, sizeof(int32_t), pc, fd);
#else //NO64BITID
	  my_fwrite(blockid, sizeof(int64_t), pc, fd);
#endif //NO64BITID
	  pc = 0;
	}
    }
  if(pc > 0)
    {
#ifdef NO64BITID
      my_fwrite(blockid, sizeof(int32_t), pc, fd);
#else //NO64BITID
      my_fwrite(blockid, sizeof(int64_t), pc, fd);
#endif
    }

  my_fwrite(&dummy, sizeof(dummy), 1, fd);

  /* write zero temperatures if needed */
  if(header.npart[0])
    {
      dummy = sizeof(float) * header.npart[0];
#ifdef FORMAT_TWO
          /*Write temperature header*/
	  blkheadsize = sizeof(int) + 4 * sizeof(char);
      	  my_fwrite(&blkheadsize,sizeof(int),1,fd);
      	  my_fwrite("U   ", sizeof(char), 4, fd);
	  nextblock=dummy+2*sizeof(int);
	  my_fwrite(&nextblock, sizeof(int), 1, fd);
	  my_fwrite(&blkheadsize, sizeof(int), 1, fd);
	  /*Done writing temperature header*/
#endif //FORMAT_TWO
      my_fwrite(&dummy, sizeof(dummy), 1, fd);

      for(i = 0, pc = 0; i < header.npart[0]; i++)
	{
	  block[pc] = 0;

	  pc++;

	  if(pc == blockmaxlen)
	    {
	      my_fwrite(block, sizeof(float), pc, fd);
	      pc = 0;
	    }
	}
      if(pc > 0)
	my_fwrite(block, sizeof(float), pc, fd);
      my_fwrite(&dummy, sizeof(dummy), 1, fd);
    }

  free(block);

  fclose(fd);
  printf("Finished writing IC file.\n");
  return;
}


/* This catches I/O errors occuring for my_fwrite(). In this case we better stop.
 */
size_t my_fwrite(void *ptr, size_t size, size_t nmemb, FILE * stream)
{
  size_t nwritten;

  if((nwritten = fwrite(ptr, size, nmemb, stream)) != nmemb)
    {
      printf("I/O error (fwrite) has occured.\n");
      fflush(stdout);
      exit(777);
    }
  return nwritten;
}

int compare_type(const void *a, const void *b)
{
  if(((struct part_data *) a)->Type < (((struct part_data *) b)->Type))
    return -1;

  if(((struct part_data *) a)->Type > (((struct part_data *) b)->Type))
    return +1;

  return 0;
}
