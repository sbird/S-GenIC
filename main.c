#include <math.h>
#include <stdlib.h>
#ifdef DOUBLEPRECISION_FFTW
   #include <drfftw_mpi.h>
#else
   #include <srfftw_mpi.h>
#endif
#include <mpi.h>
#include <gsl/gsl_rng.h>

#include "allvars.h"
#include "proto.h"


int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
  MPI_Comm_size(MPI_COMM_WORLD, &NTask);

  if(argc < 2)
    {
      if(ThisTask == 0)
	{
	  fprintf(stdout, "\nParameters are missing.\n");
	  fprintf(stdout, "Call with <ParameterFile>\n\n");
	}
      MPI_Finalize();
      exit(0);
    }

  read_parameterfile(argv[1]);

  set_units();

  initialize_powerspectrum();

  initialize_ffts();

  read_glass(GlassFile);

  displacement_fields();

  write_particle_data();

  if(NumPart)
    free(P);

  free_ffts();


  if(ThisTask == 0)
    {
      printf("\nIC's generated.\n\n");
      printf("Initial scale factor = %g\n", InitTime);
      printf("\n");
    }

  MPI_Barrier(MPI_COMM_WORLD);
/*   print_spec(); */

  MPI_Finalize();		/* clean up & finalize MPI */
  exit(0);
}





void displacement_fields(void)
{
  MPI_Request request;
  MPI_Status status;
  gsl_rng *random_generator;
  int i, j, k, ii, jj, kk, axes;
  int n;
  int sendTask, recvTask;
  double fac, vel_prefac;
  double kvec[3], kmag, kmag2, p_of_k;
  double delta, phase, ampl, hubble_a;
  double u, v, w;
  double f1, f2, f3, f4, f5, f6, f7, f8;
  double dis, maxdisp, max_disp_glob;
  unsigned int *seedtable;

#ifdef CORRECT_CIC
  double fx, fy, fz, ff, smth;
#endif

  if(ThisTask == 0)
    {
      printf("\nstart computing displacement fields...\n");
      fflush(stdout);
    }

  hubble_a =
    Hubble * sqrt(Omega / pow(InitTime, 3) + (1 - Omega - OmegaLambda) / pow(InitTime, 2) + OmegaLambda);

  vel_prefac = InitTime * hubble_a * F_Omega(InitTime);

  vel_prefac /= sqrt(InitTime);	/* converts to Gadget velocity */

  if(ThisTask == 0)
    printf("vel_prefac= %g  hubble_a=%g fom=%g \n", vel_prefac, hubble_a, F_Omega(InitTime));

  fac = pow(2 * PI / Box, 1.5);

  maxdisp = 0;

  random_generator = gsl_rng_alloc(gsl_rng_ranlxd1);

  gsl_rng_set(random_generator, Seed);

  if(!(seedtable = malloc(Nmesh * Nmesh * sizeof(unsigned int))))
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


#if defined(MULTICOMPONENTGLASSFILE) && defined(DIFFERENT_TRANSFER_FUNC)
  for(Type = MinType; Type <= MaxType; Type++)
#endif
    {
#if defined(MULTICOMPONENTGLASSFILE) && defined(DIFFERENT_TRANSFER_FUNC)
      if(ThisTask==0)
          fprintf(stderr, "Starting type %d\n",Type);
#endif
      for(axes = 0; axes < 3; axes++)
	{
	  if(ThisTask == 0)
	    {
	      printf("\nstarting axes=%d...\n", axes);
	      fflush(stdout);
	    }

	  /* first, clean the array */
	  for(i = 0; i < Local_nx; i++)
	    for(j = 0; j < Nmesh; j++)
	      for(k = 0; k <= Nmesh / 2; k++)
		{
		  Cdata[(i * Nmesh + j) * (Nmesh / 2 + 1) + k].re = 0;
		  Cdata[(i * Nmesh + j) * (Nmesh / 2 + 1) + k].im = 0;
		}

	  for(i = 0; i < Nmesh; i++)
	    {
	      ii = Nmesh - i;
	      if(ii == Nmesh)
		ii = 0;
	      if((i >= Local_x_start && i < (Local_x_start + Local_nx)) ||
		 (ii >= Local_x_start && ii < (Local_x_start + Local_nx)))
		{
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

			  p_of_k = PowerSpec(kmag);

			  p_of_k *= -log(ampl);

			  delta = fac * sqrt(p_of_k) / Dplus;	/* scale back to starting redshift */

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
			      if(i >= Local_x_start && i < (Local_x_start + Local_nx))
				{
				  Cdata[((i - Local_x_start) * Nmesh + j) * (Nmesh / 2 + 1) + k].re =
				    -kvec[axes] / kmag2 * delta * sin(phase);
				  Cdata[((i - Local_x_start) * Nmesh + j) * (Nmesh / 2 + 1) + k].im =
				    kvec[axes] / kmag2 * delta * cos(phase);
				}
			    }
			  else	/* k=0 plane needs special treatment */
			    {
			      if(i == 0)
				{
				  if(j >= Nmesh / 2)
				    continue;
				  else
				    {
				      if(i >= Local_x_start && i < (Local_x_start + Local_nx))
					{
					  jj = Nmesh - j;	/* note: j!=0 surely holds at this point */

					  Cdata[((i - Local_x_start) * Nmesh + j) * (Nmesh / 2 + 1) + k].re =
					    -kvec[axes] / kmag2 * delta * sin(phase);
					  Cdata[((i - Local_x_start) * Nmesh + j) * (Nmesh / 2 + 1) + k].im =
					    kvec[axes] / kmag2 * delta * cos(phase);

					  Cdata[((i - Local_x_start) * Nmesh + jj) * (Nmesh / 2 + 1) + k].re =
					    -kvec[axes] / kmag2 * delta * sin(phase);
					  Cdata[((i - Local_x_start) * Nmesh + jj) * (Nmesh / 2 + 1) + k].im =
					    -kvec[axes] / kmag2 * delta * cos(phase);
					}
				    }
				}
			      else	/* here comes i!=0 : conjugate can be on other processor! */
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

				      if(i >= Local_x_start && i < (Local_x_start + Local_nx))
					{
					  Cdata[((i - Local_x_start) * Nmesh + j) * (Nmesh / 2 + 1) + k].re =
					    -kvec[axes] / kmag2 * delta * sin(phase);
					  Cdata[((i - Local_x_start) * Nmesh + j) * (Nmesh / 2 + 1) + k].im =
					    kvec[axes] / kmag2 * delta * cos(phase);
					}

				      if(ii >= Local_x_start && ii < (Local_x_start + Local_nx))
					{
					  Cdata[((ii - Local_x_start) * Nmesh + jj) * (Nmesh / 2 + 1) +
						k].re = -kvec[axes] / kmag2 * delta * sin(phase);
					  Cdata[((ii - Local_x_start) * Nmesh + jj) * (Nmesh / 2 + 1) +
						k].im = -kvec[axes] / kmag2 * delta * cos(phase);
					}
				    }
				}
			    }
			}
		    }
		}
	    }


	  rfftwnd_mpi(Inverse_plan, 1, Disp, Workspace, FFTW_NORMAL_ORDER);	/** FFT **/

	  /* now get the plane on the right side from neighbour on the right, 
	     and send the left plane */

	  recvTask = ThisTask;
	  do
	    {
	      recvTask--;
	      if(recvTask < 0)
		recvTask = NTask - 1;
	    }
	  while(Local_nx_table[recvTask] == 0);

	  sendTask = ThisTask;
	  do
	    {
	      sendTask++;
	      if(sendTask >= NTask)
		sendTask = 0;
	    }
	  while(Local_nx_table[sendTask] == 0);

	  /* use non-blocking send */

	  if(Local_nx > 0)
	    {
	      MPI_Isend(&Disp[0],
			sizeof(fftw_real) * Nmesh * (2 * (Nmesh / 2 + 1)),
			MPI_BYTE, recvTask, 10, MPI_COMM_WORLD, &request);

	      MPI_Recv(&Disp[(Local_nx * Nmesh) * (2 * (Nmesh / 2 + 1))],
		       sizeof(fftw_real) * Nmesh * (2 * (Nmesh / 2 + 1)),
		       MPI_BYTE, sendTask, 10, MPI_COMM_WORLD, &status);

	      MPI_Wait(&request, &status);
	    }


	  /* read-out displacements */

	  for(n = 0; n < NumPart; n++)
	    {
#if defined(MULTICOMPONENTGLASSFILE) && defined(DIFFERENT_TRANSFER_FUNC)
	      if(P[n].Type == Type)
#endif
		{
		  u = P[n].Pos[0] / Box * Nmesh;
		  v = P[n].Pos[1] / Box * Nmesh;
		  w = P[n].Pos[2] / Box * Nmesh;

		  i = (int) u;
		  j = (int) v;
		  k = (int) w;

		  if(i == (Local_x_start + Local_nx))
		    i = (Local_x_start + Local_nx) - 1;
		  if(i < Local_x_start)
		    i = Local_x_start;
		  if(j == Nmesh)
		    j = Nmesh - 1;
		  if(k == Nmesh)
		    k = Nmesh - 1;

		  u -= i;
		  v -= j;
		  w -= k;

		  i -= Local_x_start;
		  ii = i + 1;
		  jj = j + 1;
		  kk = k + 1;

		  if(jj >= Nmesh)
		    jj -= Nmesh;
		  if(kk >= Nmesh)
		    kk -= Nmesh;

		  f1 = (1 - u) * (1 - v) * (1 - w);
		  f2 = (1 - u) * (1 - v) * (w);
		  f3 = (1 - u) * (v) * (1 - w);
		  f4 = (1 - u) * (v) * (w);
		  f5 = (u) * (1 - v) * (1 - w);
		  f6 = (u) * (1 - v) * (w);
		  f7 = (u) * (v) * (1 - w);
		  f8 = (u) * (v) * (w);

		  dis = Disp[(i * Nmesh + j) * (2 * (Nmesh / 2 + 1)) + k] * f1 +
		    Disp[(i * Nmesh + j) * (2 * (Nmesh / 2 + 1)) + kk] * f2 +
		    Disp[(i * Nmesh + jj) * (2 * (Nmesh / 2 + 1)) + k] * f3 +
		    Disp[(i * Nmesh + jj) * (2 * (Nmesh / 2 + 1)) + kk] * f4 +
		    Disp[(ii * Nmesh + j) * (2 * (Nmesh / 2 + 1)) + k] * f5 +
		    Disp[(ii * Nmesh + j) * (2 * (Nmesh / 2 + 1)) + kk] * f6 +
		    Disp[(ii * Nmesh + jj) * (2 * (Nmesh / 2 + 1)) + k] * f7 +
		    Disp[(ii * Nmesh + jj) * (2 * (Nmesh / 2 + 1)) + kk] * f8;
                  if(isnan(dis)){
                   fprintf(stderr, "i=%d j=%d k=%d Nmesh=%d f1=%g f2=%g f3=%g f4=%g f5=%g f6=%g f7=%g f8=%g\n",i,j,k,Nmesh,f1,f2,f3,f4,f5,f6,f7,f8);
                   fprintf(stderr, "%g %g %g %g %g %g %g %g\n",Disp[(i * Nmesh + j) * (2 * (Nmesh / 2 + 1)) + k],Disp[(i * Nmesh + j) * (2 * (Nmesh / 2 + 1)) + kk],
		    Disp[(i * Nmesh + jj) * (2 * (Nmesh / 2 + 1)) + k],
		    Disp[(i * Nmesh + jj) * (2 * (Nmesh / 2 + 1)) + kk],
		    Disp[(ii * Nmesh + j) * (2 * (Nmesh / 2 + 1)) + k], 
		    Disp[(ii * Nmesh + j) * (2 * (Nmesh / 2 + 1)) + kk],
		    Disp[(ii * Nmesh + jj) * (2 * (Nmesh / 2 + 1)) + k],
		    Disp[(ii * Nmesh + jj) * (2 * (Nmesh / 2 + 1)) + kk]);
                   FatalError(2);
                  }

		  P[n].Vel[axes] = dis;
		  if(dis > maxdisp)
		    maxdisp = dis;
		}
	    }
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

  MPI_Reduce(&maxdisp, &max_disp_glob, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      printf("\nMaximum displacement: %g kpc/h, in units of the part-spacing= %g\n",
	     max_disp_glob, max_disp_glob / (Box / Nmesh));
    }
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
  int total_size, i, additional;
  int local_ny_after_transpose, local_y_start_after_transpose;
  int *slab_to_task_local;
  size_t bytes;


  Inverse_plan = rfftw3d_mpi_create_plan(MPI_COMM_WORLD,
					 Nmesh, Nmesh, Nmesh, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE);

  rfftwnd_mpi_local_sizes(Inverse_plan, &Local_nx, &Local_x_start,
			  &local_ny_after_transpose, &local_y_start_after_transpose, &total_size);

  Local_nx_table = malloc(sizeof(int) * NTask);
  MPI_Allgather(&Local_nx, 1, MPI_INT, Local_nx_table, 1, MPI_INT, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      for(i = 0; i < NTask; i++)
	printf("Task=%d Local_nx=%d\n", i, Local_nx_table[i]);
      fflush(stdout);
    }


  Slab_to_task = malloc(sizeof(int) * Nmesh);
  slab_to_task_local = malloc(sizeof(int) * Nmesh);

  for(i = 0; i < Nmesh; i++)
    slab_to_task_local[i] = 0;

  for(i = 0; i < Local_nx; i++)
    slab_to_task_local[Local_x_start + i] = ThisTask;

  MPI_Allreduce(slab_to_task_local, Slab_to_task, Nmesh, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  free(slab_to_task_local);



  additional = (Nmesh) * (2 * (Nmesh / 2 + 1));	/* additional plane on the right side */

  Disp = (fftw_real *) malloc(bytes = sizeof(fftw_real) * (total_size + additional));
  Workspace = (fftw_real *) malloc(bytes += sizeof(fftw_real) * total_size);

  if(Disp && Workspace)
    {
      if(ThisTask == 0)
	printf("\nallocated %g Mbyte on Task %d for FFT's\n", bytes / (1024.0 * 1024.0), ThisTask);
    }
  else
    {
      printf("failed to allocate %g Mbyte on Task %d\n", bytes / (1024.0 * 1024.0), ThisTask);
      printf("bailing out.\n");
      FatalError(1);
    }

  Cdata = (fftw_complex *) Disp;	/* transformed array */
}



void free_ffts(void)
{
  free(Workspace);
  free(Disp);
  free(Slab_to_task);
  rfftwnd_mpi_destroy_plan(Inverse_plan);
}


int FatalError(int errnum)
{
  printf("FatalError called with number=%d\n", errnum);
  fflush(stdout);
  MPI_Abort(MPI_COMM_WORLD, errnum);
  exit(0);
}




static double A, B, alpha, beta, V, gf;

double fnl(double x)		/* Peacock & Dodds formula */
{
  return x * pow((1 + B * beta * x + pow(A * x, alpha * beta)) /
		 (1 + pow(pow(A * x, alpha) * gf * gf * gf / (V * sqrt(x)), beta)), 1 / beta);
}

void print_spec(void)
{
  double k, knl, po, dl, dnl, neff, kf, kstart, kend, po2, po1, DDD;
  char buf[1000];
  FILE *fd;

  if(ThisTask == 0)
    {
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
	  po = PowerSpec(k);
	  dl = 4.0 * PI * k * k * k * po;

	  kf = 0.5;

	  po2 = PowerSpec(1.001 * k * kf);
	  po1 = PowerSpec(k * kf);

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
}
