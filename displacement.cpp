#include "displacement.hpp"

#include <math.h>
#include <stdlib.h>
#include <fftw3.h>
#include <gsl/gsl_rng.h>
#include <exception>
//For memset
#include <string.h>
#include <omp.h>

#include "power.hpp"
#include "part_data.hpp"
// #include "proto.h"

/**Initialise the memory for the FFTs*/
DisplacementFields::DisplacementFields(size_t Nmesh, size_t Nsample, int Seed, double Box, bool twolpt): Nmesh(Nmesh), Nsample(Nsample), twolpt(twolpt), Seed(Seed), Box(Box)
{
  size_t bytes = sizeof(float) * 2*Nmesh*Nmesh*(Nmesh/2+1);
  Disp = (float *) fftwf_malloc(bytes);
  if(Disp)
        printf("Nmesh = %lu. Allocated %lu MB for FFTs\n",Nmesh,  bytes / (1024 * 1024));
  else{
      fprintf(stderr, "Nmesh = %lu. Failed to allocate %lu MB for FFTs\n",Nmesh, bytes / (1024 * 1024));
      throw std::bad_alloc();
  }
  Cdata = (fftwf_complex *) Disp;	/* transformed array */

  if (twolpt) {
        twosrc = (float *) fftwf_malloc(bytes);
        ctwosrc = (fftwf_complex *) twosrc;
        for(int i=0; i<3; i++){
            cdigrad[i] = (fftwf_complex *) malloc(bytes);
            digrad[i] = (float *) cdigrad[i];
        }
        /*Check memory allocation*/
        if(cdigrad[0] && cdigrad[1] && cdigrad[2] && twosrc)
                printf("Allocated %lu MB for 2LPT term\n",4*bytes / (1024 * 1024));
        else{
            fprintf(stderr, "Failed to allocate %lu MB for 2LPT term\n",4*bytes / (1024 * 1024));
            throw std::bad_alloc();
        }
  }

  if(!fftwf_init_threads()){
  		  fprintf(stderr,"Error initialising fftw threads\n");
  		  exit(1);
  }
  fftwf_plan_with_nthreads(omp_get_num_procs());
  Inverse_plan = fftwf_plan_dft_c2r_3d(Nmesh, Nmesh, Nmesh,Cdata,Disp, FFTW_ESTIMATE);
  if(twolpt) {
    Forward_plan2 = fftwf_plan_dft_r2c_3d(Nmesh, Nmesh, Nmesh,twosrc,ctwosrc, FFTW_ESTIMATE);
    for(int i=0; i<3; i++){
            Inverse_plan_grad[i] = fftwf_plan_dft_c2r_3d(Nmesh, Nmesh, Nmesh,cdigrad[i],digrad[i], FFTW_ESTIMATE);
    }
  }
}

//Free the memory in the FFTs
DisplacementFields::~DisplacementFields()
{
    fftwf_free(Disp);
    fftwf_destroy_plan(Inverse_plan);
    if (twolpt) {
        /* Free  */
        fftwf_free(twosrc);
        fftwf_destroy_plan(Forward_plan2);
        for(int i=0;i<3;i++){
            fftwf_free(cdigrad[i]);
            fftwf_destroy_plan(Inverse_plan_grad[i]);
        }
    }
}


/**Little macro to work the storage order of the FFT.*/
inline int KVAL(const size_t n, const size_t Nmesh)
{
    return (n < Nmesh / 2 ? n : n - Nmesh);
}


#ifdef CORRECT_CIC
/* do deconvolution of CIC interpolation */
double invwindow(int kx,int ky,int kz,int n)
{
	double iwx=1.0,iwy=1.0,iwz=1.0;
        if(!n)
                return 0;
	if(kx){
		iwx=M_PI*kx/static_cast<float>(n);
		iwx=iwx/sin(iwx);
        }
	if(ky){
		iwy=M_PI*ky/static_cast<float>(n);
		iwy=iwy/sin(iwy);
        }
	if(kz){
		iwz=M_PI*kz/static_cast<float>(n);
		iwz=iwz/sin(iwz);
        }
	return pow(iwx*iwy*iwz,2);
}
#endif

/**Initialise a table of seeds for the random number generator*/
unsigned int * initialize_rng(int Seed, size_t Nmesh)
{
  gsl_rng * random_generator = gsl_rng_alloc(gsl_rng_ranlxd1);
  gsl_rng_set(random_generator, Seed);
  unsigned int * seedtable;

  if(!(seedtable =(unsigned int *) malloc(Nmesh * Nmesh * sizeof(unsigned int))))
    throw std::bad_alloc();

  for(size_t i = 0; i < Nmesh / 2; i++)
    {
      size_t j;
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

/**Function to compute Zeldovich displacement fields using a double Fourier transform*/
void DisplacementFields::displacement_fields(const int type, part_data& P, PowerSpec * PSpec, bool SphereMode, bool RayleighScatter)
{
  const double fac = pow(2 * M_PI / Box, 1.5);
  //Re-initialize every time this is called, so each particle type has the same phases
  unsigned int *seedtable = initialize_rng(Seed, Nmesh);
  double maxdisp;
  const size_t fftsize = 2*Nmesh*Nmesh*(Nmesh/2+1);

  double maxdisp2=0;
  if(twolpt)
    memset(twosrc, 0, fftsize*sizeof(float));

  for(int axes = 0; axes < 3; axes++) {
	  printf("Starting Zeldovich axis %d.\n", axes);
	  /* first, clean the array */
      memset(Disp, 0, fftsize*sizeof(float));
      #pragma omp parallel
	  {
          gsl_rng * random_generator = gsl_rng_alloc(gsl_rng_ranlxd1);
	      #pragma omp for
	      for(size_t i = 0; i < Nmesh; i++) {
		    for(size_t j = 0; j < Nmesh; j++) {
		      gsl_rng_set(random_generator, seedtable[i * Nmesh + j]);

		      for(size_t k = 0; k < Nmesh / 2; k++) {
                          double kvec[3], kmag, kmag2, p_of_k;
                          double delta, phase, ampl;
			  phase = gsl_rng_uniform(random_generator) * 2 * M_PI;
			  do
			    ampl = gsl_rng_uniform(random_generator);
			  while(ampl == 0);

			  if(i == Nmesh / 2 || j == Nmesh / 2 || k == Nmesh / 2)
			    continue;
			  if(i == 0 && j == 0 && k == 0)
			    continue;

			  kvec[0] = KVAL(i, Nmesh) * 2 * M_PI / Box;
			  kvec[1] = KVAL(j, Nmesh) * 2 * M_PI / Box;
			  kvec[2] = KVAL(k, Nmesh) * 2 * M_PI / Box;

			  kmag2 = kvec[0] * kvec[0] + kvec[1] * kvec[1] + kvec[2] * kvec[2];
			  kmag = sqrt(kmag2);

                          /* select a sphere in k-space */
			  if(SphereMode){
			      if(kmag * Box / (2 * M_PI) > Nsample / 2)
                                      continue;
                          }
                          /*Or a box*/
			  else {
			      if(fabs(kvec[0]) * Box / (2 * M_PI) > Nsample / 2)
				continue;
			      if(fabs(kvec[1]) * Box / (2 * M_PI) > Nsample / 2)
				continue;
			      if(fabs(kvec[2]) * Box / (2 * M_PI) > Nsample / 2)
				continue;
			  }

			  p_of_k = PSpec->power(kmag, type);

			  // printf(" k %d %g %g \n",Type,kmag,p_of_k);
                          if(RayleighScatter)
			        p_of_k *= -log(ampl);

			  delta = fac * sqrt(p_of_k) ;
                          /* scale back to starting redshift */
                          /*If we are using the CAMB P(k), Dplus=1.
                            * fac=(2π/Box)^1.5*/
#ifdef CORRECT_CIC
              /* do deconvolution of CIC interpolation */
              delta *= invwindow(i,j,k,Nmesh);
#endif
              if(k > 0) {
                  size_t index = (i * Nmesh + j) * (Nmesh / 2 + 1) + k;
                  (Cdata[index])[0] = -kvec[axes] / kmag2 * delta * sin(phase);
                  (Cdata[index])[1] = kvec[axes] / kmag2 * delta * cos(phase);
                }
              else    /* k=0 plane needs special treatment */
                {
                  if(i == 0)
                {
                  if(j >= Nmesh / 2)
                    continue;
                  else
                    {
                      int jj = Nmesh - j;    /* note: i=k=0 => j!=0 */
                      size_t index = j * (Nmesh / 2 + 1);
                      (Cdata[index])[0] = -kvec[axes] / kmag2 * delta * sin(phase);
                      (Cdata[index])[1] =  kvec[axes] / kmag2 * delta * cos(phase);
                      index = jj * (Nmesh / 2 + 1);
                      (Cdata[index])[0] = -kvec[axes] / kmag2 * delta * sin(phase);
                      (Cdata[index])[1] = -kvec[axes] / kmag2 * delta * cos(phase);
                    }
                }
                  else    /* here comes i!=0 */
                {
                  if(i >= Nmesh / 2)
                    continue;
                  else
                    {
                      int ii = (Nmesh - i) % Nmesh;
                      int jj = (Nmesh - j) % Nmesh;
                      size_t index = (i * Nmesh + j) * (Nmesh / 2 + 1);
                      (Cdata[index])[0] = -kvec[axes] / kmag2 * delta * sin(phase);
                      (Cdata[index])[1] =  kvec[axes] / kmag2 * delta * cos(phase);
                      index = (ii * Nmesh + jj) * (Nmesh / 2 + 1);
                      (Cdata[index])[0] = -kvec[axes] / kmag2 * delta * sin(phase);
                      (Cdata[index])[1] = -kvec[axes] / kmag2 * delta * cos(phase);
                    }
                }
                }
            }

          }
        }
        gsl_rng_free(random_generator);
        #pragma omp barrier
     } //omp_parallel
     if (twolpt ) {
            /* At this point, Cdata contains the complex Zeldovich displacement for this axis */

            /* Compute displacement gradient
            * do disp(0,0), disp(0,1), disp(0,2), disp(1,1), disp(1,2), disp(2,2) only as vector symmetric*/
            for(int ax=2;ax>=axes; ax--){
        #ifdef NEUTRINOS
                if(type == 2)
                    break;
        #endif
                #pragma omp parallel for
                for(size_t i = 0; i < Nmesh; i++) {
                    for(size_t j = 0; j < Nmesh; j++) {
                    for(size_t k = 0; k <= Nmesh / 2; k++){
                        double kvec[3];
                        size_t coord = (i * Nmesh + j) * (Nmesh / 2 + 1) + k;
                        kvec[0] = KVAL(i, Nmesh) * 2 * M_PI / Box;
                        kvec[1] = KVAL(j, Nmesh) * 2 * M_PI / Box;
                        kvec[2] = KVAL(k, Nmesh) * 2 * M_PI / Box;
                        /*Note that unlike Scoccimaro et al, we do not have
                        * memory to waste, so we only do one axis at a time */
                        /* Derivatives of ZA displacement  */
                        /* d(dis_i)/d(q_j)  -> sqrt(-1) k_j dis_i */
                        (cdigrad[axes][coord])[0] = -(Cdata[coord])[1] * kvec[ax]; /* disp0,0 */
                        (cdigrad[axes][coord])[1] = (Cdata[coord])[0] * kvec[ax];
                    }
                    }
                }
                /*At this point, cdigrad[i] contains FT(phi,ii). For grad^2 phi, want the FT. */
        //           printf("Finding gradient FT component (%d,%d)\n",ax,axes);
                fftwf_execute(Inverse_plan_grad[axes]);	/** FFT of cdigrad[axes] **/

                /* Compute second order source and store it in twosrc*/
                if(ax != axes)
                    #pragma omp parallel for
                    for(size_t i = 0; i < Nmesh; i++)
                        for(size_t j = 0; j < Nmesh; j++)
                        for(size_t k = 0; k < Nmesh; k++){
                            size_t coord = (i * Nmesh + j) * (2 * (Nmesh / 2 + 1)) + k;
                            twosrc[coord] -= digrad[axes][coord]*digrad[axes][coord];
                        }
            }
     } //end twolpt
     fftwf_execute(Inverse_plan);	/** FFT of the Zeldovich displacements **/
     /* read-out Zeldovich displacements into P.Vel*/
     maxdisp=displacement_read_out(1, P, axes);
    }

    if (twolpt) {
#ifdef NEUTRINOS
    if(type != 2){
#endif
      /* So now digrad[axes] contains phi,ii and twosrc contains  sum_(i>j)(- phi,ij^2)
       * We want to now compute phi,ii^(2), the laplacian of the 2LPT term, in twosrc */
      #pragma omp parallel for
      for(size_t i = 0; i < Nmesh; i++)
        for(size_t j = 0; j < Nmesh; j++)
          for(size_t k = 0; k < Nmesh; k++){
            size_t co = (i * Nmesh + j) * (2 * (Nmesh / 2 + 1)) + k;
            twosrc[co] += digrad[0][co]*digrad[1][co]+digrad[0][co]*digrad[2][co]+digrad[1][co]*digrad[2][co];
          }
      fftwf_execute(Forward_plan2);	/** FFT of twosrc**/
      for(int axes=0; axes< 3; axes++){
              printf("Starting 2LPT term, axis %d\n",axes);
              /* Reuse the memory used earlier for ZA field */
              memset(Disp, 0, fftsize*sizeof(float));
              /* Solve Poisson eq. and calculate 2nd order displacements */
              (Cdata[0])[0] = (Cdata[0])[1] = 0.0;
              #pragma omp parallel for
              for(size_t i = 0; i < Nmesh; i++)
                for(size_t j = 0; j < Nmesh; j++)
                    for(size_t k = 1; k <= Nmesh / 2; k++){
                      double kvec[3],kmag2;
                      size_t coord = (i * Nmesh + j) * (Nmesh / 2 + 1) + k;
                      kvec[0] = KVAL(i, Nmesh) * 2 * M_PI / Box;
                      kvec[1] = KVAL(j, Nmesh) * 2 * M_PI / Box;
                      kvec[2] = KVAL(k, Nmesh) * 2 * M_PI / Box;

                      kmag2 = kvec[0] * kvec[0] + kvec[1] * kvec[1] + kvec[2] * kvec[2];
                      /* cdisp2 = source * k / (sqrt(-1) k^2) */
                      (Cdata[coord])[0] = (ctwosrc[coord])[1] * kvec[axes] / kmag2;
                      (Cdata[coord])[1] = -(ctwosrc[coord])[0] * kvec[axes] / kmag2;
#ifdef CORRECT_CIC
                      /* do deconvolution of CIC interpolation */
                      double smth= invwindow(i,j,k,Nmesh);
                      (Cdata[coord])[0] *= smth;
                      (Cdata[coord])[1] *= smth;
#endif
               }

              /* Cdata now contains the FFT of the 2LPT term */
              fftwf_execute(Inverse_plan);	/** FFT of Cdata**/
              /* read-out displacements */
              maxdisp2=displacement_read_out(2, P, axes);
          }
#ifdef NEUTRINOS
    } //type !=2
#endif
    } //end twolpt

  printf("\nMaximum Zeldovich displacement: %g kpc/h, in units of the part-spacing= %g\n",
         maxdisp, maxdisp / (Box / Nmesh));
  printf("\nMaximum 2LPT displacement: %g kpc/h, in units of the part-spacing= %g\n",maxdisp2, maxdisp2 / (Box / Nmesh));
  free(seedtable);
  return;
}

double DisplacementFields::displacement_read_out(const int order, part_data& P, const int axes)
{
   double maxx=0;
   const double Nmesh3 = pow(1.*Nmesh, 3);
   const int64_t NumPart = P.GetNumPart();
   #pragma omp parallel
   {
      double maxdisp=0;
      #pragma omp for
      for(int n = 0; n < NumPart; n++)
        {
                  double dis;
                  double f1, f2, f3, f4, f5, f6, f7, f8;
                  double u[3];
                  size_t i[3], ii[3];
                  for(int q=0;q<3;q++){
                      u[q] = P.Pos(n,q) / Box * Nmesh;
                      i[q] = static_cast<size_t>(u[q]);
                      if(i[q] == Nmesh)
                          i[q]--;
                      u[q] -= i[q];
                      ii[q] = i[q]+1;
                      if(ii[q] >= Nmesh)
                          ii[q] -= Nmesh;
                  }

		  f1 = (1 - u[0]) * (1 - u[1]) * (1 - u[2]);
		  f2 = (1 - u[0]) * (1 - u[1]) * (u[2]);
		  f3 = (1 - u[0]) * (u[1]) * (1 - u[2]);
		  f4 = (1 - u[0]) * (u[1]) * (u[2]);
		  f5 = (u[0]) * (1 - u[1]) * (1 - u[2]);
		  f6 = (u[0]) * (1 - u[1]) * (u[2]);
		  f7 = (u[0]) * (u[1]) * (1 - u[2]);
		  f8 = (u[0]) * (u[1]) * (u[2]);

		  dis = Disp[(i[0] * Nmesh + i[1]) * (2 * (Nmesh / 2 + 1)) + i[2]] * f1 +
		    Disp[(i[0] * Nmesh + i[1]) * (2 * (Nmesh / 2 + 1)) + ii[2]] * f2 +
		    Disp[(i[0] * Nmesh + ii[1]) * (2 * (Nmesh / 2 + 1)) + i[2]] * f3 +
		    Disp[(i[0] * Nmesh + ii[1]) * (2 * (Nmesh / 2 + 1)) + ii[2]] * f4 +
		    Disp[(ii[0] * Nmesh + i[1]) * (2 * (Nmesh / 2 + 1)) + i[2]] * f5 +
		    Disp[(ii[0] * Nmesh + i[1]) * (2 * (Nmesh / 2 + 1)) + ii[2]] * f6 +
		    Disp[(ii[0] * Nmesh + ii[1]) * (2 * (Nmesh / 2 + 1)) + i[2]] * f7 +
		    Disp[(ii[0] * Nmesh + ii[1]) * (2 * (Nmesh / 2 + 1)) + ii[2]] * f8;
          /*Read out the 2lpt velocity if this is
           * being called from the 2lpt part of the code*/
          if(order == 2){
              dis /= Nmesh3;
              P.Set2Vel(dis, n,axes);
          }
          else
              P.SetVel(dis, n,axes);
          if(dis > maxdisp)
            maxdisp = dis;
        }
        #pragma omp critical
        {
              if (maxdisp > maxx)
                      maxx=maxdisp;
        }
   } //end omp_parallel
   return maxx;
}
