#include <math.h>
#include <stdlib.h>

#include "allvars.h"
#include "proto.h"

#ifdef NO64BITID
        typedef int32_t id_type;
#else
        typedef int64_t id_type;
#endif //NO64BITID

void write_particle_data(void)
{
  size_t bytes;
  float *block;
  id_type *blockid;
  int blockmaxlen;
  int64_t maxidlen;
  int4byte dummy;
  FILE *fd;
  char buf[300];
  int64_t i, k, pc;
  gadget_header header;
    
  printf("\nWriting IC-file\n");
    
  #define BUFFER 10

  if(NumPart == 0)
    return;

  sprintf(buf, "%s/%s", OutputDir, FileBase);

  if(!(fd = fopen(buf, "w"))){
      printf("Error. Can't write in file '%s'\n", buf);
      exit(10);
  }

  /* sort particles by type, because that's how they should be stored in a gadget binary file */
  qsort(P, NumPart, sizeof(struct part_data), compare_type);	
  /*Write header*/
  /*FIXME: Multiple files*/
  header = generate_header(0,1);
  write_block(fd, "HEAD", &header,sizeof(header));

  if(!(block = malloc(bytes = BUFFER * 1024 * 1024)))
    {
      printf("failed to allocate memory for `block' (%g bytes).\n", (double) bytes);
      exit(24);
    }

  blockmaxlen = bytes / (3 * sizeof(float));

  /* write coordinates */
  dummy = sizeof(float) * 3 * NumPart;
#ifdef NEUTRINO_PAIRS
  dummy =
    sizeof(float) * 3 * (header.npart[0] + header.npart[1] + header.npart[2] + header.npart[3] +
			 header.npart[4] + header.npart[5]);
#endif //NEUTRINO_PAIRS

  /*We are about to write the POS block*/
  write_block_header(fd, "POS ",dummy);


  for(i = 0, pc = 0; i < NumPart; i++){
      for(k = 0; k < 3; k++)
	  block[3 * pc + k] = P[i].Pos[k];
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
  /*Done writing POS block*/
  write_block_footer(fd, "POS ",dummy);



  /* write velocities */
  dummy = sizeof(float) * 3 * NumPart;
#ifdef NEUTRINO_PAIRS
  dummy =
    sizeof(float) * 3 * (header.npart[0] + header.npart[1] + header.npart[2] + header.npart[3] +
			 header.npart[4] + header.npart[5]);
#endif //NEUTRINO_PAIRS

  write_block_header(fd, "VEL ",dummy);

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

  write_block_footer(fd, "VEL ",dummy);


  /* write particle ID */

  blockid = (id_type *) block;
  maxidlen = bytes / (sizeof(id_type));
  dummy = sizeof(id_type) * NumPart;
#ifdef NEUTRINO_PAIRS
  dummy =
    sizeof(int) * (header.npart[0] + header.npart[1] + header.npart[2] + header.npart[3] + header.npart[4] +
		   header.npart[5]);
#endif //NEUTRINO_PAIRS
  write_block_header(fd, "ID  ",dummy);

  for(i = 0, pc = 0; i < NumPart; i++){
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

      if(pc >= (maxidlen - 1)){
	  my_fwrite(blockid, sizeof(id_type), pc, fd);
	  pc = 0;
      }
  }
  if(pc > 0)
      my_fwrite(blockid, sizeof(id_type), pc, fd);

  write_block_footer(fd, "ID  ",dummy);
  /*Done writing IDs*/

  /* write zero temperatures if needed */
  if(header.npart[0])
    {
      dummy = sizeof(float) * header.npart[0];

      write_block_header(fd, "U   ",dummy);

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

      write_block_footer(fd, "U   ",dummy);
    }
  /*Done writing temperatures*/

  free(block);

  fclose(fd);
  printf("Finished writing IC file.\n");
  return;
}

gadget_header generate_header(int file, int num_files)
{
  gadget_header header;
  double scale = 3 * Hubble * Hubble / (8 * M_PI * G) * pow(Box,3);
  int i;
  /*Set particle numbers*/
  for(i = 0; i < N_TYPE; i++){
    uint64_t npart_t = header1.npartTotal[i] * GlassTileFac * GlassTileFac * GlassTileFac;
#ifdef NEUTRINO_PAIRS
    if(i == NEUTRINO_TYPE)
        npart_t *= 2;
#endif //NEUTRINO_PAIRS
    if(npart_t > ((int64_t)1<<32) ) {
        header.NallHW[i] = ( npart_t >> 32);
        header.npartTotal[i] = npart_t - ((int64_t)header.NallHW[i] << 32);
    }
    /* Note this is the canonical decision as to how many particles each file gets: 
     * everything else just reads this value.*/
    header.npart[i]=npart_t/num_files;
    /* If num_files is not a factor of npartTotal, give 
     * n extra particles to the final file .*/
    if(file == num_files -1) 
            header.npart[i] += npart_t % num_files;
  }
  /*Set masses*/
  for(i = 0; i < N_TYPE; i++)
      header.mass[i] = 0;

  if(header.npartTotal[BARYON_TYPE])
    header.mass[BARYON_TYPE] = (OmegaBaryon) * scale / (header.npartTotal[BARYON_TYPE]);

  if(header.npartTotal[DM_TYPE])
    header.mass[DM_TYPE] = (Omega - OmegaBaryon - OmegaDM_2ndSpecies) * scale / (header.npartTotal[DM_TYPE]);

  if(header.npartTotal[NEUTRINO_TYPE]){
    header.mass[NEUTRINO_TYPE] = (OmegaDM_2ndSpecies) * scale / (header.npartTotal[NEUTRINO_TYPE]);
#ifdef NEUTRINO_PAIRS
    header.mass[NEUTRINO_TYPE] /= 2;
#endif //NEUTRINO_PAIRS
  }

  header.num_files = num_files;

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

void write_block(FILE * fd, char * name, void * block, int blocksize)
{
  write_block_header(fd, name, blocksize);
  my_fwrite(block, blocksize, 1, fd);
  write_block_footer(fd, name, blocksize);
}

void write_block_header(FILE * fd, char * name, int blocksize)
{
#ifdef FORMAT_TWO
      /*This is the block header record, which we want for format two files*/
      int blkheadsize = sizeof(int) + 4 * sizeof(char);
      int nextblock = blocksize + 2 * sizeof(int);
      /*Write format 2 header header*/
      my_fwrite(&blkheadsize,sizeof(int),1,fd);
      my_fwrite(name, sizeof(char), 4, fd);
      my_fwrite(&nextblock, sizeof(int), 1, fd);
      my_fwrite(&blkheadsize,sizeof(int),1,fd);
#endif //FORMAT_TWO
      /*This is the record size, which we want for all files*/
      my_fwrite(&blocksize, sizeof(int), 1, fd);
      return;
}

void write_block_footer(FILE * fd, char * name, int blocksize)
{
      /*This is the record size, which we want for all files*/
      my_fwrite(&blocksize, sizeof(int), 1, fd);
      return;
}

