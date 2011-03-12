#include <math.h>
#include <stdlib.h>
#include <fftw3.h>
#include <gsl/gsl_rng.h>

#include "allvars.h"
#include "proto.h"

#include <gadgetreader.hpp>
#include <gadgetwriter.hpp>

int main(int argc, char **argv)
{
  int type;
  std::valarray<int64_t> npart(N_TYPE);
  int64_t FirstId=0;
  struct part_data *P=NULL;

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
      NumPart = read_glass(snap, type, GlassTileFac, P);
      displacement_fields(type, NumPart, P, Nmesh);
      FirstId = write_particle_data(osnap, type,P, NumPart,FirstId);
      free(P);
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

void displacement_fields(const int type, const int64_t NumPart, struct part_data* P, const int Nmesh)
{
  const double fac = pow(2 * M_PI / Box, 1.5);
  const double hubble_a = Hubble * sqrt(Omega / pow(InitTime, 3) + (1 - Omega - OmegaLambda) / pow(InitTime, 2) + OmegaLambda);
  const unsigned int *seedtable = initialize_rng(Seed);
  const double vel_prefac = InitTime * hubble_a * F_Omega(InitTime) /sqrt(InitTime);
  double mindisp=0, maxdisp=0;

/*I really think this is not right; Omega should be specified as total matter density, not cdm matter*/
/*  if(neutrinos_ks)
    Omega = Omega + OmegaDM_2ndSpecies;*/

  /* the final term converts to Gadget velocity */

  printf("vel_prefac= %g  hubble_a=%g fom=%g Omega=%g \n", vel_prefac, hubble_a, F_Omega(InitTime), Omega);

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
			  // p_of_k *= -log(ampl);

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
        		  u[q] = P[n].Pos[q] / Box * Nmesh;
                          i[q] = (int64_t) u[q];
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

		  P[n].Vel[axes] = dis;
		  if(dis > maxdisp)
		    maxdisp = dis;
                  if(dis <mindisp)
                    mindisp=dis;
	    }
        } //end omp_parallel
	}
  /* now add displacement to Lagrangian coordinates, and multiply velocities by correct factor */
  #pragma omp parallel 
  {
  #pragma omp for 
  for(int64_t n = 0; n < NumPart; n++)
    {
      for(int axes = 0; axes < 3; axes++)
	{
	  P[n].Pos[axes] = periodic_wrap(P[n].Pos[axes] + P[n].Vel[axes]);
	  P[n].Vel[axes] *= vel_prefac;
	}
    }
  }


  printf("\nMaximum displacement: %g kpc/h, in units of the part-spacing= %g\n",
         maxdisp, maxdisp / (Box / Nmesh));
  printf("Minimum displacement: %g kpc/h, in units of the part-spacing= %g\n",
         mindisp, mindisp / (Box / Nmesh));
  return;
}

double periodic_wrap(double x)
{
  while(x >= Box)
    x -= Box;

  while(x < 0)
    x += Box;

  return x;
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
	double iwx,iwy,iwz;
        if(!n)
                return 0;
	if(!kx)
		iwx=1.0;
	else
		iwx=M_PI*kx/(n*sin(M_PI*kx/(float)n));
	if(!ky)
		iwy=1.0;
	else
		iwy=M_PI*ky/(n*sin(M_PI*ky/(float)n));
	if(!kz)
		iwz=1.0;
	else
		iwz=M_PI*kz/(n*sin(M_PI*kz/(float)n));
	return pow(iwx*iwy*iwz,2);
}
#endif

