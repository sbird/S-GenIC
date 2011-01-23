#include <math.h>
#include <stdlib.h>
#include <omp.h>
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

  read_parameterfile(argv[1]);
  NumFiles = 4;
  set_units();

  initialize_powerspectrum();
  printf("Dplus initial redshift =%g  \n\n", Dplus); 

  initialize_ffts();
  printf("Initialising pre-IC file '%s'\n",GlassFile);
  GadgetReader::GSnap snap(GlassFile);
  /*Set particle numbers*/
  for(type = 0; type < N_TYPE; type++)
    npart[type] = snap.GetNpart(type) * GlassTileFac * GlassTileFac * GlassTileFac;

#ifdef NEUTRINO_PAIRS
  npart[NEUTRINO_TYPE] *= 2;
#endif //NEUTRINO_PAIRS

  GadgetWriter::GWriteSnap osnap(std::string(OutputDir)+std::string("/")+std::string(FileBase), npart,NumFiles, sizeof(id_type));
  /*Write headers*/
  osnap.WriteHeaders(generate_header());

  for(type=0; type<N_TYPE;type++){
      int64_t NumPart = 0;
      if(npart[type] == 0)
              continue;
      NumPart = read_glass(snap, type, GlassTileFac, P);
      displacement_fields(type, NumPart, P);
      FirstId = write_particle_data(osnap, type,P, NumPart,FirstId);
      free(P);
  }

  fftwf_free(Disp);
  fftwf_destroy_plan(Inverse_plan);

  printf("Initial scale factor = %g\n", InitTime);

/*   print_spec(); */

  return 0;
}

void initialize_rng(gsl_rng * random_generator, unsigned int*& seedtable)
{
  int i,j;

  gsl_rng_set(random_generator, Seed);

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


}

void displacement_fields(int type, int64_t NumPart, struct part_data* P)
{
  int i, j, k, ii, jj, axes;
  int n;
  double fac, vel_prefac;
  double kvec[3], kmag, kmag2, p_of_k;
  double delta, phase, ampl, hubble_a;
  double mindisp=0, maxdisp=0;
  gsl_rng * random_generator = gsl_rng_alloc(gsl_rng_ranlxd1);
  unsigned int *seedtable=NULL;
  initialize_rng(random_generator, seedtable);

#ifdef CORRECT_CIC
  double fx, fy, fz, ff, smth;
#endif

/*I really think this is not right; Omega should be specified as total matter density, not cdm matter*/
/*  if(neutrinos_ks)
    Omega = Omega + OmegaDM_2ndSpecies;*/

  hubble_a = Hubble * sqrt(Omega / pow(InitTime, 3) + (1 - Omega - OmegaLambda) / pow(InitTime, 2) + OmegaLambda);

  vel_prefac = InitTime * hubble_a * F_Omega(InitTime);

  vel_prefac /= sqrt(InitTime);	/* converts to Gadget velocity */

  printf("vel_prefac= %g  hubble_a=%g fom=%g Omega=%g \n", vel_prefac, hubble_a, F_Omega(InitTime), Omega);

  fac = pow(2 * PI / Box, 1.5);


      for(axes = 0; axes < 3; axes++) {
	  printf("Starting axis %d.\n", axes);

	  /* first, clean the array */
	  for(i = 0; i < Nmesh; i++)
	    for(j = 0; j < Nmesh; j++)
	      for(k = 0; k <= Nmesh / 2; k++) {
		  (Cdata[(i * Nmesh + j) * (Nmesh / 2 + 1) + k])[0] = 0;
		  (Cdata[(i * Nmesh + j) * (Nmesh / 2 + 1) + k])[1] = 0;
	      }

	  for(i = 0; i < Nmesh; i++) {
	      ii = Nmesh - i;
	      if(ii == Nmesh)
		ii = 0;
		  for(j = 0; j < Nmesh; j++)
		    {
		      gsl_rng_set(random_generator, seedtable[i * Nmesh + j]);

		      for(k = 0; k < Nmesh / 2; k++)
			{
			  phase = gsl_rng_uniform(random_generator) * 2 * PI;
			  do
			    ampl = gsl_rng_uniform(random_generator);
			  while(ampl == 0);

			  if(i == Nmesh / 2 || j == Nmesh / 2 || k == Nmesh / 2)
			    continue;
			  if(i == 0 && j == 0 && k == 0)
			    continue;

			  if(i < Nmesh / 2)
			    kvec[0] = i * 2 * PI / Box;
			  else
			    kvec[0] = -(Nmesh - i) * 2 * PI / Box;

			  if(j < Nmesh / 2)
			    kvec[1] = j * 2 * PI / Box;
			  else
			    kvec[1] = -(Nmesh - j) * 2 * PI / Box;

			  if(k < Nmesh / 2)
			    kvec[2] = k * 2 * PI / Box;
			  else
			    kvec[2] = -(Nmesh - k) * 2 * PI / Box;

			  kmag2 = kvec[0] * kvec[0] + kvec[1] * kvec[1] + kvec[2] * kvec[2];
			  kmag = sqrt(kmag2);

			  if(SphereMode == 1)
			    {
			      if(kmag * Box / (2 * PI) > Nsample / 2)	/* select a sphere in k-space */
				continue;
			    }
			  else
			    {
			      if(fabs(kvec[0]) * Box / (2 * PI) > Nsample / 2)
				continue;
			      if(fabs(kvec[1]) * Box / (2 * PI) > Nsample / 2)
				continue;
			      if(fabs(kvec[2]) * Box / (2 * PI) > Nsample / 2)
				continue;
			    }

			  p_of_k = PowerSpec(kmag, type);

			  // printf(" k %d %g %g \n",Type,kmag,p_of_k);
			  // p_of_k *= -log(ampl);

			  delta = fac * sqrt(p_of_k) / Dplus;	/* scale back to starting redshift */
           /*If we are using the CAMB P(k), Dplus=1.
            * fac=(2π/Box)^1.5*/

#ifdef CORRECT_CIC
			  /* do deconvolution of CIC interpolation */
			  fx = fy = fz = 1;
			  if(kvec[0] != 0)
			    {
			      fx = (kvec[0] * Box / 2) / Nmesh;
			      fx = sin(fx) / fx;
			    }
			  if(kvec[1] != 0)
			    {
			      fy = (kvec[1] * Box / 2) / Nmesh;
			      fy = sin(fy) / fy;
			    }
			  if(kvec[2] != 0)
			    {
			      fz = (kvec[2] * Box / 2) / Nmesh;
			      fz = sin(fz) / fz;
			    }
			  ff = 1 / (fx * fy * fz);
			  smth = ff * ff;

			  delta *= smth;
			  /* end deconvolution */
#endif
			  if(k > 0)
			    {
				  (Cdata[(i * Nmesh + j) * (Nmesh / 2 + 1) + k])[0] =
				    -kvec[axes] / kmag2 * delta * sin(phase);
				  (Cdata[(i * Nmesh + j) * (Nmesh / 2 + 1) + k])[1] =
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
					  jj = Nmesh - j;	/* note: j!=0 surely holds at this point */

					  (Cdata[(i * Nmesh + j) * (Nmesh / 2 + 1) + k])[0] =
					    -kvec[axes] / kmag2 * delta * sin(phase);
					  (Cdata[(i * Nmesh + j) * (Nmesh / 2 + 1) + k])[1] =
					    kvec[axes] / kmag2 * delta * cos(phase);

					  (Cdata[(i * Nmesh + jj) * (Nmesh / 2 + 1) + k])[0] =
					    -kvec[axes] / kmag2 * delta * sin(phase);
					  (Cdata[(i * Nmesh + jj) * (Nmesh / 2 + 1) + k])[1] =
					    -kvec[axes] / kmag2 * delta * cos(phase);
				    }
				}
			      else	/* here comes i!=0 */
				{
				  if(i >= Nmesh / 2)
				    continue;
				  else
				    {
				      ii = Nmesh - i;
				      if(ii == Nmesh)
					ii = 0;
				      jj = Nmesh - j;
				      if(jj == Nmesh)
					jj = 0;

					  (Cdata[(i * Nmesh + j) * (Nmesh / 2 + 1)])[0] =
					    -kvec[axes] / kmag2 * delta * sin(phase);
					  (Cdata[(i * Nmesh + j) * (Nmesh / 2 + 1)])[1] =
					    kvec[axes] / kmag2 * delta * cos(phase);
					  (Cdata[(i * Nmesh + jj) * (Nmesh / 2 + 1) +
					        k])[0] = -kvec[axes] / kmag2 * delta * sin(phase);
					  (Cdata[(i * Nmesh + jj) * (Nmesh / 2 + 1) +
						k])[1] = -kvec[axes] / kmag2 * delta * cos(phase);
				    }
				}
			    }
			}
		    }
	    }

	  fftwf_execute(Inverse_plan);	/** FFT **/

	  /* read-out displacements */

	  for(n = 0; n < NumPart; n++)
	    {
                  double dis;
                  double f1, f2, f3, f4, f5, f6, f7, f8;
                  double u[3];
                  int i[3], ii[3];
                  int q;
                  for(q=0;q<3;q++){
        		  u[q] = P[n].Pos[q] / Box * Nmesh;
                          i[q] = (int) u[q];
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
	}

  /* now add displacement to Lagrangian coordinates, and multiply velocities by correct factor */
  for(n = 0; n < NumPart; n++)
    {
      for(axes = 0; axes < 3; axes++)
	{
	  P[n].Pos[axes] += P[n].Vel[axes];
	  P[n].Vel[axes] *= vel_prefac;
	  P[n].Pos[axes] = periodic_wrap(P[n].Pos[axes]);
	}
    }

  gsl_rng_free(random_generator);

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


void set_units(void)		/* ... set some units */
{
  UnitTime_in_s = UnitLength_in_cm / UnitVelocity_in_cm_per_s;

  G = GRAVITY / pow(UnitLength_in_cm, 3) * UnitMass_in_g * pow(UnitTime_in_s, 2);
  Hubble = HUBBLE * UnitTime_in_s;
}



void initialize_ffts(void)
{
  size_t bytes;

  Disp = (float *) fftwf_malloc(bytes = sizeof(float) * (2*Nmesh*Nmesh*(Nmesh/2+1)));
  if(Disp)
        printf("\nAllocated %g MB for FFTs\n", bytes / (1024.0 * 1024.0));
  else{
      fprintf(stderr, "Failed to allocate %g Mbyte\n", bytes / (1024.0 * 1024.0));
      FatalError(1);
  }
  Cdata = (fftwf_complex *) Disp;	/* transformed array */

  if(!fftwf_init_threads()){
  		  fprintf(stderr,"Error initialising fftw threads\n");
  		  exit(1);
  }
  fftwf_plan_with_nthreads(omp_get_num_procs());
  Inverse_plan = fftwf_plan_dft_c2r_3d(Nmesh, Nmesh, Nmesh,Cdata,Disp, FFTW_ESTIMATE);
  return;
}


int FatalError(int errnum)
{
  printf("FatalError called with number=%d\n", errnum);
  fflush(stdout);
  exit(0);
}


static double A, B, alpha, beta, V, gf;

double fnl(double x)		/* Peacock & Dodds formula */
{
  return x * pow((1 + B * beta * x + pow(A * x, alpha * beta)) /
		 (1 + pow(pow(A * x, alpha) * gf * gf * gf / (V * sqrt(x)), beta)), 1 / beta);
}

void print_spec(int type)
{
  double k, knl, po, dl, dnl, neff, kf, kstart, kend, po2, po1, DDD;
  char buf[1000];
  FILE *fd;

      sprintf(buf, "%s/inputspec_%s.txt", OutputDir, FileBase);

      fd = fopen(buf, "w");

      gf = GrowthFactor(0.001, 1.0) / (1.0 / 0.001);

      DDD = GrowthFactor(1.0 / (Redshift + 1), 1.0);

      fprintf(fd, "%12g %12g\n", Redshift, DDD);	/* print actual starting redshift and 
							   linear growth factor for this cosmology */

      kstart = 2 * PI / (1000.0 * (3.085678e24 / UnitLength_in_cm));	/* 1000 Mpc/h */
      kend = 2 * PI / (0.001 * (3.085678e24 / UnitLength_in_cm));	/* 0.001 Mpc/h */

      for(k = kstart; k < kend; k *= 1.025)
	{
	  po = PowerSpec(k, type);
          //printf(" po k %g %g\n ",k,po);
	  dl = 4.0 * PI * k * k * k * po;

	  kf = 0.5;

	  po2 = PowerSpec(1.001 * k * kf, type);
	  po1 = PowerSpec(k * kf, type);

	  if(po != 0 && po1 != 0 && po2 != 0)
	    {
	      neff = (log(po2) - log(po1)) / (log(1.001 * k * kf) - log(k * kf));

	      if(1 + neff / 3 > 0)
		{
		  A = 0.482 * pow(1 + neff / 3, -0.947);
		  B = 0.226 * pow(1 + neff / 3, -1.778);
		  alpha = 3.310 * pow(1 + neff / 3, -0.244);
		  beta = 0.862 * pow(1 + neff / 3, -0.287);
		  V = 11.55 * pow(1 + neff / 3, -0.423) * 1.2;

		  dnl = fnl(dl);
		  knl = k * pow(1 + dnl, 1.0 / 3);
		}
	      else
		{
		  dnl = 0;
		  knl = 0;
		}
	    }
	  else
	    {
	      dnl = 0;
	      knl = 0;
	    }

	  fprintf(fd, "%12g %12g %12g  %12g %12g\n", k,po, dl, knl, dnl);
	}
      fclose(fd);
}

gadget_header generate_header()
{
  gadget_header header;
  double scale = 3 * Hubble * Hubble / (8 * M_PI * G) * pow(Box,3);
  /*Set masses*/
  for(int i = 0; i < N_TYPE; i++)
      header.mass[i] = 0;

  if(header1.npartTotal[BARYON_TYPE])
    header.mass[BARYON_TYPE] = (OmegaBaryon) * scale / (header.npartTotal[BARYON_TYPE]);

  if(header1.npartTotal[DM_TYPE])
    header.mass[DM_TYPE] = (Omega - OmegaBaryon - OmegaDM_2ndSpecies) * scale / (header.npartTotal[DM_TYPE]);

  if(header1.npartTotal[NEUTRINO_TYPE]){
    header.mass[NEUTRINO_TYPE] = (OmegaDM_2ndSpecies) * scale / (header.npartTotal[NEUTRINO_TYPE]);
#ifdef NEUTRINO_PAIRS
    header.mass[NEUTRINO_TYPE] /= 2;
#endif //NEUTRINO_PAIRS
  }

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

