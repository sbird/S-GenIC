#include <math.h>
#include <stdlib.h>
#include <fftw3.h>
#include <gsl/gsl_rng.h>

#include "allvars.h"
#include "proto.h"
#include "part_data.hpp"

#include <gadgetreader.hpp>
#include <gadgetwriter.hpp>

int main(int argc, char **argv)
{
  int type;
  std::valarray<int64_t> npart(N_TYPE);
  int64_t FirstId=0;

  if(argc < 2)
    {
	  fprintf(stdout, "\nParameters are missing.\n");
	  fprintf(stdout, "Call with <ParameterFile>\n\n");
      exit(0);
    }
  /*Make sure stdout is line buffered even when not 
   * printing to a terminal but, eg, perl*/
  setlinebuf(stdout);
  read_parameterfile(argv[1]);
  set_units();

  printf("Nmesh = %d Nsample = %d\n",Nmesh,Nsample);
  initialize_ffts();
  printf("Initialising pre-IC file '%s'\n",GlassFile);
  GadgetReader::GSnap snap(GlassFile);
  /*Set particle numbers*/
  for(type = 0; type < N_TYPE; type++)
    npart[type] = snap.GetNpart(type) * GlassTileFac * GlassTileFac * GlassTileFac;
  if(npart.sum() == 0)
          exit(1);
  /*Set the global variable saying there is no gas in the glassfile,
   * so that the OmegaBaryon should be added to the DM.*/
  if (npart[BARYON_TYPE] == 0)
          no_gas = 1;
  else
          no_gas = 0;
  /*We need to initialise the power spectrum here so that no_gas is set*/
  initialize_powerspectrum();

#ifdef NEUTRINO_PAIRS
  npart[NEUTRINO_TYPE] *= 2;
#endif //NEUTRINO_PAIRS
  GadgetWriter::GWriteSnap osnap(std::string(OutputDir)+std::string("/")+std::string(FileBase), npart,NumFiles, sizeof(id_type));
  /*Write headers*/
  if(osnap.WriteHeaders(generate_header(npart)))
          FatalError(23);

  for(type=0; type<N_TYPE;type++){
      int64_t NumPart = 0;
      if(npart[type] == 0)
              continue;
      part_data P(snap, type, GlassTileFac);
      NumPart = P.GetNumPart();
      displacement_fields(type, NumPart, P, Nmesh);
      FirstId = write_particle_data(osnap, type,P, NumPart,FirstId);
#ifdef PRINT_SPEC
      print_spec(type);
#endif
  }

  fftwf_free(Disp);
  fftwf_destroy_plan(Inverse_plan);

  printf("Initial scale factor = %g\n", InitTime);

  return 0;
}

/**Little macro to work the storage order of the FFT.*/
#define KVAL(n) ((n)< Nmesh/2 ? (n) : ((n)-Nmesh))

void displacement_fields(const int type, const int64_t NumPart, part_data& P, const int Nmesh)
{
  const double fac = pow(2 * M_PI / Box, 1.5);
  const unsigned int *seedtable = initialize_rng(Seed);
  double mindisp=0, maxdisp=0;

/*I really think this is not right; Omega should be specified as total matter density, not cdm matter*/
/*  if(neutrinos_ks)
    Omega = Omega + OmegaDM_2ndSpecies;*/
#ifdef TWOLPT
  /* the final term converts to Gadget velocity */
      for(size_t i = 0; i < ((size_t) 2*Nmesh*Nmesh)*(Nmesh/2+1); i++){
              twosrc[i]=0;
      }
#endif

      for(int axes = 0; axes < 3; axes++) {
	  printf("Starting axis %d.\n", axes);

          #pragma omp parallel
	  {
          gsl_rng * random_generator = gsl_rng_alloc(gsl_rng_ranlxd1);
	  /* first, clean the array */
	  #pragma omp for 
	  for(size_t i = 0; i < ((size_t) 2*Nmesh*Nmesh)*(Nmesh/2+1); i++)
		  Disp[i] = 0;

	  #pragma omp for 
	  for(int i = 0; i < Nmesh; i++) {
		  for(int j = 0; j < Nmesh; j++) {
		      gsl_rng_set(random_generator, seedtable[i * Nmesh + j]);

		      for(int k = 0; k < Nmesh / 2; k++) {
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

			  kvec[0] = KVAL(i) * 2 * M_PI / Box;
			  kvec[1] = KVAL(j) * 2 * M_PI / Box;
			  kvec[2] = KVAL(k) * 2 * M_PI / Box;

			  kmag2 = kvec[0] * kvec[0] + kvec[1] * kvec[1] + kvec[2] * kvec[2];
			  kmag = sqrt(kmag2);

                          /* select a sphere in k-space */
			  if(SphereMode == 1){
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

			  p_of_k = PowerSpec(kmag, type);

			  // printf(" k %d %g %g \n",Type,kmag,p_of_k);
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
				  (Cdata[index])[0] =
				    -kvec[axes] / kmag2 * delta * sin(phase);
				  (Cdata[index])[1] =
				    kvec[axes] / kmag2 * delta * cos(phase);
			    }
			  else	/* k=0 plane needs special treatment */
			    {
			      if(i == 0)
				{
				  if(j >= Nmesh / 2)
				    continue;
				  else
				    {
					  int jj = Nmesh - j;	/* note: i=k=0 => j!=0 */
                                          size_t index = j * (Nmesh / 2 + 1);
					  (Cdata[index])[0] =
					    -kvec[axes] / kmag2 * delta * sin(phase);
					  (Cdata[index])[1] =
					    kvec[axes] / kmag2 * delta * cos(phase);

                                          index = jj * (Nmesh / 2 + 1);
					  (Cdata[index])[0] =
					    -kvec[axes] / kmag2 * delta * sin(phase);
					  (Cdata[index])[1] =
					    -kvec[axes] / kmag2 * delta * cos(phase);
				    }
				}
			      else	/* here comes i!=0 */
				{
				  if(i >= Nmesh / 2)
				    continue;
				  else
				    {
				      int jj = (Nmesh - j) % Nmesh;
                                      size_t index = (i * Nmesh + j) * (Nmesh / 2 + 1);
					  (Cdata[index])[0] =
					    -kvec[axes] / kmag2 * delta * sin(phase);
					  (Cdata[index])[1] =
					    kvec[axes] / kmag2 * delta * cos(phase);
                                      index = (i * Nmesh + jj) * (Nmesh / 2 + 1);
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
#ifdef TWOLPT
      /* At this point, cdisp[axes] contains the complex Zeldovich displacement */

      printf("Done Zeldovich.\n");
      
      /* Compute displacement gradient */

        /*This is the index of the gradient vector; 
         * do disp(0,0), disp(0,1), disp(0,2), disp(1,1), disp(1,2), disp(2,2) only as vector symmetric*/
      for(int ax=2;ax>=axes; ax--){ 
              for(int i = 0; i < Nmesh; i++)
        	for(int j = 0; j < Nmesh; j++)
        	  for(int k = 0; k <= Nmesh / 2; k++)
        	    {
                      double kvec[3];
        	      size_t coord = (i * Nmesh + j) * (Nmesh / 2 + 1) + k;
                      kvec[0] = KVAL(i) * 2 * M_PI / Box;
                      kvec[1] = KVAL(j) * 2 * M_PI / Box;
                      kvec[2] = KVAL(k) * 2 * M_PI / Box;
        	      /*Note that unlike Scoccimaro et al, we do not have 
                       * memory to waste, so we only do one axis at a time */ 
        	      /* Derivatives of ZA displacement  */
        	      /* d(dis_i)/d(q_j)  -> sqrt(-1) k_j dis_i */
        	      (cdigrad[axes][coord])[0] = -(Cdata[coord])[1] * kvec[ax]; /* disp0,0 */
        	      (cdigrad[axes][coord])[1] = (Cdata[coord])[0] * kvec[ax];
        	    }

       printf("Fourier transforming displacement gradient...%d",ax);
       cdigrad_0=cdigrad[axes];
       digrad_0=digrad[axes];
       fftwf_execute(Forward_plan_grad);	/** FFT **/

      /* Compute second order source and store it in digrad[3]*/

      for(int i = 0; i < Nmesh; i++)
	for(int j = 0; j < Nmesh; j++)
	  for(int k = 0; k < Nmesh; k++)
	    {
	      size_t coord = (i * Nmesh + j) * (2 * (Nmesh / 2 + 1)) + k;
                // g_00(g_11+g_22) + g_11*g_22 -g_01*g_01-g_02*g_02-g12_g12
              if(ax != axes)
        	      twosrc[coord] -= digrad[axes][coord]*digrad[axes][coord];
//		digrad[0][coord]*(digrad[3][coord]+digrad[5][coord])+digrad[3][coord]*digrad[5][coord]
//               -digrad[1][coord]*digrad[1][coord]-digrad[2][coord]*digrad[2][coord]-digrad[4][coord]*digrad[4][coord];
	    }
      }
#endif
	  fftwf_execute(Inverse_plan);	/** FFT **/

	  /* read-out displacements */
        #pragma omp parallel 
          {
	  #pragma omp for 
	  for(int n = 0; n < NumPart; n++)
	    {
                  double dis;
                  double f1, f2, f3, f4, f5, f6, f7, f8;
                  double u[3];
                  int64_t i[3], ii[3];
                  for(int q=0;q<3;q++){
        		  u[q] = P.Pos(n,q) / Box * Nmesh;
                          i[q] = static_cast<int64_t>(u[q]);
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
		  P.SetVel(dis, n,axes);
		  if(dis > maxdisp)
		    maxdisp = dis;
                  if(dis <mindisp)
                    mindisp=dis;
	    }
        } //end omp_parallel
	}
#ifdef TWOLPT
      digrad_0 = twosrc;
      cdigrad_0 = (fftwf_complex *) twosrc;
      fftwf_execute(Forward_plan_grad);	/** FFT of twosrc**/
      for(int axes=0; axes< 3; axes++){
      printf("Fourier transforming second order source...");
      
      /* The memory allocated for cdigrad[0], [1], and [2] will be used for 2nd order displacements */
      /* Freeing the rest. cdigrad[3] still has 2nd order displacement source, free later */

      /* Solve Poisson eq. and calculate 2nd order displacements */

      for(int i = 0; i < Nmesh; i++)
	for(int j = 0; j < Nmesh; j++)
	  for(int k = 0; k <= Nmesh / 2; k++)
	    {
              double kvec[3],kmag2;
	      size_t coord = (i * Nmesh + j) * (Nmesh / 2 + 1) + k;
              kvec[0] = KVAL(i) * 2 * M_PI / Box;
              kvec[1] = KVAL(j) * 2 * M_PI / Box;
              kvec[2] = KVAL(k) * 2 * M_PI / Box;

	      kmag2 = kvec[0] * kvec[0] + kvec[1] * kvec[1] + kvec[2] * kvec[2];
#ifdef CORRECT_CIC
	      /* do deconvolution of CIC interpolation */
	      double smth= invwindow(i,j,k,Nmesh);
#endif

	      /* cdisp2 = source * k / (sqrt(-1) k^2) */
		  if(kmag2 > 0.0) 
		    {
		      (cdisp2[coord])[0] = (((fftwf_complex *)twosrc)[coord])[1] * kvec[axes] / kmag2;
		      (cdisp2[coord])[1] = -(((fftwf_complex *)twosrc)[coord])[0] * kvec[axes] / kmag2;
		    }
		  else 
                      (cdisp2[coord])[0] = (cdisp2[coord])[1] = 0.0;
#ifdef CORRECT_CIC
		  (cdisp2[coord])[0] *= smth;  
                  (cdisp2[coord])[1] *= smth;
#endif
	    }
      
      /* Free cdigrad[3] */
      free(twosrc);
      /* Now, both cdisp, and cdisp2 have the ZA and 2nd order displacements */
	  fftwf_execute(Inverse_plan2);	/** FFT **/
	  /* read-out displacements */
        #pragma omp parallel 
          {
	  #pragma omp for 
	  for(int n = 0; n < NumPart; n++)
	    {
                  double dis2;
                  double f1, f2, f3, f4, f5, f6, f7, f8;
                  double u[3];
                  int64_t i[3], ii[3];
                  for(int q=0;q<3;q++){
        		  u[q] = P.Pos(n,q) / Box * Nmesh;
                          i[q] = static_cast<int64_t>(u[q]);
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

		  dis2 = disp2[(i[0] * Nmesh + i[1]) * (2 * (Nmesh / 2 + 1)) + i[2]] * f1 +
		    disp2[(i[0] * Nmesh + i[1]) * (2 * (Nmesh / 2 + 1)) + ii[2]] * f2 +
		    disp2[(i[0] * Nmesh + ii[1]) * (2 * (Nmesh / 2 + 1)) + i[2]] * f3 +
		    disp2[(i[0] * Nmesh + ii[1]) * (2 * (Nmesh / 2 + 1)) + ii[2]] * f4 +
		    disp2[(ii[0] * Nmesh + i[1]) * (2 * (Nmesh / 2 + 1)) + i[2]] * f5 +
		    disp2[(ii[0] * Nmesh + i[1]) * (2 * (Nmesh / 2 + 1)) + ii[2]] * f6 +
		    disp2[(ii[0] * Nmesh + ii[1]) * (2 * (Nmesh / 2 + 1)) + i[2]] * f7 +
		    disp2[(ii[0] * Nmesh + ii[1]) * (2 * (Nmesh / 2 + 1)) + ii[2]] * f8;
		  dis2 /= ((float) Nmesh)*Nmesh*Nmesh;
		  P.Set2Vel(dis2, n,axes);
	    }
        } //end omp_parallel
	}
#endif
      
  printf("\nMaximum displacement: %g kpc/h, in units of the part-spacing= %g\n",
         maxdisp, maxdisp / (Box / Nmesh));
  printf("Minimum displacement: %g kpc/h, in units of the part-spacing= %g\n",
         mindisp, mindisp / (Box / Nmesh));
  return;
}


int FatalError(int errnum)
{
  printf("FatalError called with number=%d\n", errnum);
  fflush(stdout);
  exit(0);
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

