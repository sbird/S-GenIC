#include "proto.h"
#include "allvars.h"
#include <omp.h>

unsigned int * initialize_rng(int Seed)
{
  int i,j;

  gsl_rng * random_generator = gsl_rng_alloc(gsl_rng_ranlxd1);
  gsl_rng_set(random_generator, Seed);
  unsigned int * seedtable;

  if(!(seedtable =(unsigned int *) malloc(Nmesh * Nmesh * sizeof(unsigned int))))
    FatalError(4);

  for(i = 0; i < Nmesh / 2; i++)
    {
      for(j = 0; j < i; j++)
	seedtable[i * Nmesh + j] = 0x7fffffff * gsl_rng_uniform(random_generator);

      for(j = 0; j < i + 1; j++)
	seedtable[j * Nmesh + i] = 0x7fffffff * gsl_rng_uniform(random_generator);

      for(j = 0; j < i; j++)
	seedtable[(Nmesh - 1 - i) * Nmesh + j] = 0x7fffffff * gsl_rng_uniform(random_generator);

      for(j = 0; j < i + 1; j++)
	seedtable[(Nmesh - 1 - j) * Nmesh + i] = 0x7fffffff * gsl_rng_uniform(random_generator);

      for(j = 0; j < i; j++)
	seedtable[i * Nmesh + (Nmesh - 1 - j)] = 0x7fffffff * gsl_rng_uniform(random_generator);

      for(j = 0; j < i + 1; j++)
	seedtable[j * Nmesh + (Nmesh - 1 - i)] = 0x7fffffff * gsl_rng_uniform(random_generator);

      for(j = 0; j < i; j++)
	seedtable[(Nmesh - 1 - i) * Nmesh + (Nmesh - 1 - j)] = 0x7fffffff * gsl_rng_uniform(random_generator);

      for(j = 0; j < i + 1; j++)
	seedtable[(Nmesh - 1 - j) * Nmesh + (Nmesh - 1 - i)] = 0x7fffffff * gsl_rng_uniform(random_generator);
    }

  gsl_rng_free(random_generator);
  return seedtable;

}

void set_units(void)		/* ... set some units */
{
  UnitTime_in_s = UnitLength_in_cm / UnitVelocity_in_cm_per_s;

  G = GRAVITY / pow(UnitLength_in_cm, 3) * UnitMass_in_g * pow(UnitTime_in_s, 2);
  Hubble = HUBBLE * UnitTime_in_s;
}

void initialize_ffts(void)
{
  size_t bytes = sizeof(float) * 2*Nmesh*Nmesh*(Nmesh/2+1);
  Disp = (float *) fftwf_malloc(bytes);
#ifdef TWOLPT
  twosrc = (float *) fftwf_malloc(bytes);
  disp2 = (float *) fftwf_malloc(bytes);
  cdisp2 = (fftwf_complex *) disp2;	/* transformed array */
 for(int i=0; i<3; i++){ 
  cdigrad[i] = (fftwf_complex *) malloc(bytes);
  digrad[i] = (float *) cdigrad[i];
 }
 cdigrad_0=cdigrad[0];
 digrad_0=digrad[0];
  /*Check memory allocation*/
#endif
      
  if(Disp)
        printf("Nmesh = %d. Allocated %ld MB for FFTs\n",Nmesh,  bytes / (1024 * 1024));
  else{
      fprintf(stderr, "Nmesh = %d. Failed to allocate %ld MB for FFTs\n",Nmesh, bytes / (1024 * 1024));
      FatalError(1);
  }
  Cdata = (fftwf_complex *) Disp;	/* transformed array */

  if(!fftwf_init_threads()){
  		  fprintf(stderr,"Error initialising fftw threads\n");
  		  exit(1);
  }
  fftwf_plan_with_nthreads(omp_get_num_procs());
  Inverse_plan = fftwf_plan_dft_c2r_3d(Nmesh, Nmesh, Nmesh,Cdata,Disp, FFTW_ESTIMATE);
#ifdef TWOLPT
  Forward_plan_grad = fftwf_plan_dft_r2c_3d(Nmesh, Nmesh, Nmesh,digrad_0,cdigrad_0, FFTW_ESTIMATE);
  Inverse_plan2 = fftwf_plan_dft_c2r_3d(Nmesh, Nmesh, Nmesh,cdisp2,disp2, FFTW_ESTIMATE);
#endif
  return;
}

gadget_header generate_header(std::valarray<int64_t> & npart)
{
  gadget_header header;
  double scale = 3 * Hubble * Hubble / (8 * M_PI * G) * pow(Box,3);
  /*Set masses*/
  for(int i = 0; i < N_TYPE; i++)
      header.mass[i] = 0;
  /*Don't forget to set masses correctly when the CDM actually incorporates other species*/
  double OmegaCDM = Omega;
  if(npart[BARYON_TYPE]){
    header.mass[BARYON_TYPE] = OmegaBaryon * scale / npart[BARYON_TYPE];
    OmegaCDM -=OmegaBaryon;
  }

  if(npart[NEUTRINO_TYPE]){
    header.mass[NEUTRINO_TYPE] = OmegaDM_2ndSpecies * scale / npart[NEUTRINO_TYPE];
#ifdef NEUTRINO_PAIRS
    header.mass[NEUTRINO_TYPE] /= 2;
#endif //NEUTRINO_PAIRS
    OmegaCDM -=OmegaDM_2ndSpecies;
  }
  /*Subtract off the fake neutrinos, as their masses are
   * incorporated within the mesh*/
  else if (neutrinos_ks){
          OmegaCDM-=OmegaDM_2ndSpecies;
  }

  if(npart[DM_TYPE])
    header.mass[DM_TYPE] = OmegaCDM * scale / npart[DM_TYPE];

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
