#include <math.h>
#include "allvars.h"
#include "proto.h"
#include <map>
#include <gsl/gsl_integration.h>

static double R8;

static double AA, BB, CC;
static double nu;
static double Norm;

static int NPowerTable;

/*Structure for matter power table*/
static struct pow_table
{
  double logk, logD,logDb;
  double logD2nd;
  double logDtot;
}
 *PowerTable;

static struct pow_matter
{ 
  double kmat,pmat;
 }
*PowerMatter;

double PowerSpec(double k, int Type)
{
  double power, alpha, Tf;
  switch (WhichSpectrum)
    {
    case 1:
      power = PowerSpec_EH(k);
      break;
    case 2:
#ifndef NEUTRINOS
      if (Type == 2)
        power = PowerSpec_DM_2ndSpecies(k);
#endif
      power = PowerSpec_Tabulated(k,Type);
      break;
    default:
      power = PowerSpec_Efstathiou(k);
      break;
    }

  if(WDM_On == 1)
    {
      /* Eqn. (A9) in Bode, Ostriker & Turok (2001), assuming gX=1.5  */
      alpha =
	0.048 * pow((Omega - OmegaBaryon) / 0.4, 0.15) * pow(HubbleParam / 0.65,
							     1.3) * pow(1.0 / WDM_PartMass_in_kev, 1.15);
      Tf = pow(1 + pow(alpha * k * (3.085678e24 / UnitLength_in_cm), 2 * 1.2), -5.0 / 1.2);
      power *= Tf * Tf;
    }

  if(WhichSpectrum < 2)
    power *= pow(k, PrimordialIndex - 1.0);

  return power/ (Dplus*Dplus);
}


double PowerSpec_DM_2ndSpecies(double k)
{
  /* at the moment, we simply call the Eistenstein & Hu spectrum
   * for the second DM species, but this could be replaced with
   * something more physical, say for neutrinos
   */

  double power;

  power = Norm * k * pow(tk_eh(k), 2);

  return power;
}

void read_power_table(void)
{
  FILE *fd;
  char buf[500];
  double k;

  sprintf(buf, FileWithTransfer);
  if(!(fd = fopen(buf, "r")))
    {
      printf("can't read input TRANSFER in file '%s'\n", buf);
      exit(17);
    }

  NPowerTable = 0;
  do
    {
      double T_cdm, T_b, dummy, T_nu, T_tot;

      /* read transfer function file from CAMB */
      if(fscanf(fd, " %lg %lg %lg %lg %lg %lg %lg", &k, &T_cdm, &T_b, &dummy, &dummy, &T_nu, &T_tot) == 7)
	NPowerTable++;
      else
	break;
    }
  while(1);
  fclose(fd);
      printf("found %d rows in input TRANSFER table\n", NPowerTable);
      fflush(stdout);



sprintf(buf, FileWithInputSpectrum);
  if(!(fd = fopen(buf, "r")))
    {
      printf("can't read input SPECTRUM in file '%s'\n", buf);
      exit(17);
    }
  NPowerTable = 0;
  do
    {
      double ktab, Pktab;
      /* read TOTAL matter power spectrum from CAMB*/
      if(fscanf(fd, " %lg %lg ", &ktab, &Pktab) == 2)
	NPowerTable++;
      else
	break;
    }
  while(1);
  fclose(fd);
  printf("found %d rows in input SPECTRUM table\n", NPowerTable);
  fflush(stdout);
  PowerTable = (pow_table *) malloc(NPowerTable * sizeof(struct pow_table));  
  PowerMatter = (pow_matter *) malloc(NPowerTable * sizeof(struct pow_matter)); 
  /* define matter array */
  sprintf(buf, FileWithInputSpectrum);
  if(!(fd = fopen(buf, "r")))
    {
      printf("can't read input SPECTRUM in file '%s'\n", buf);
      exit(18);
    }

  NPowerTable = 0;
  do
    {
      double kmat, pmat;
      if(fscanf(fd, " %lg %lg", &kmat, &pmat) == 2)
	{
	  PowerMatter[NPowerTable].kmat = kmat;
	  PowerMatter[NPowerTable].pmat = pmat;
	  NPowerTable++;
	}
      else
	break;
    }
  while(1);

  fclose(fd);

  qsort(PowerMatter, NPowerTable, sizeof(struct pow_matter), compare_logk);

   sprintf(buf, FileWithTransfer);
   if(!(fd = fopen(buf, "r")))
    {
      printf("can't read input spectrum in file '%s'\n", buf);
      exit(18);
    }

  NPowerTable = 0;
  do
    {
      double T_b, dummy, T_nu, T_tot, T_cdm;
      double delta_cdm, delta_nu, delta_tot, delta_b;

      if(fscanf(fd, " %lg %lg %lg %lg %lg %lg %lg", &k, &T_cdm, &T_b, &dummy, &dummy, &T_nu, &T_tot) == 7)
	{
          if(fabs(PowerMatter[NPowerTable].kmat - k) > 0.01 *k){
                  fprintf(stderr, "Error: Input spectrum row %d has k=%g, transfer has k=%g.\n",NPowerTable,PowerMatter[NPowerTable].kmat,k);
                  fprintf(stderr, "Remember you need transfer_k_interp_matterpower = F in CAMB\n");
                  exit(47);
          }
	  PowerTable[NPowerTable].logk = log10(k);

          /* The dark matter may incorporate other particle types as well, eg,
           * fake neutrinos, or baryons.
           * NOTE that CAMB defines T_tot = (T_CDM M_CDM + T_b M_b +T_nu M_nu) / (M_CDM + M_b +M_nu)
           * HOWEVER, when relativistic (ie, in the early universe), neutrinos will redshift
           * like radiation. The CAMB transfer function takes this into account when calculating T_tot,
           * so instead of summing the transfer functions, use the total transfer function and
           * optionally subtract the baryons.*/
          if(neutrinos_ks){
                T_cdm = T_tot;
                /*If we have separate gas particles, subtract
                 * the baryon transfer function */
                if(!no_gas){
                    T_cdm = (Omega*T_tot - T_b*OmegaBaryon)/(Omega-OmegaBaryon);
                }
          }
          /*Add baryons to the CDM if there are no gas particles*/
          else if(no_gas){
                    T_cdm = (Omega*T_tot - T_nu*OmegaDM_2ndSpecies)/(Omega-OmegaDM_2ndSpecies);
//                 T_cdm = (T_cdm*(Omega-OmegaDM_2ndSpecies - OmegaBaryon) + T_b*OmegaBaryon)/(Omega-OmegaDM_2ndSpecies);
          }
          /*This should be equivalent to the above in almost all cases,
           * but perhaps we want to see the effect of changing only the ICs in CAMB for the neutrinos.*/
          if(no_gas && OmegaDM_2ndSpecies == 0){
                  T_cdm = T_tot;
          }
	  /* obtain P(k) from transfer function ratios like suggested by JL*/
	  /*NOTE for this to work CAMB's transfer_k_interp_matterpower must be off!!*/
	  delta_b   = k * k * k * pow(T_b/T_tot,2)* PowerMatter[NPowerTable].pmat/(2*M_PI*M_PI);
	  delta_cdm = k * k * k * pow(T_cdm/T_tot,2)* PowerMatter[NPowerTable].pmat/(2*M_PI*M_PI);
	  delta_nu = k * k * k * pow(T_nu/T_tot,2) * PowerMatter[NPowerTable].pmat/(2*M_PI*M_PI); 
	  delta_tot = k * k * k * PowerMatter[NPowerTable].pmat/(2*M_PI*M_PI);

          PowerTable[NPowerTable].logD = log10(delta_cdm);
	  // printf("NT,d_cdm,log10(d_cdm),k %d %g %g %g \n",NPowerTable,delta_cdm,log10(delta_cdm),k);
	  PowerTable[NPowerTable].logD2nd = log10(delta_nu);
          PowerTable[NPowerTable].logDtot = log10(delta_tot);
	  PowerTable[NPowerTable].logDb = log10(delta_b);

	  NPowerTable++;
	}
      else
	break;
    }
  while(1);

  fclose(fd);

  qsort(PowerTable, NPowerTable, sizeof(struct pow_table), compare_logk);
}



int compare_logk(const void *a, const void *b)
{
  if(((struct pow_table *) a)->logk < (((struct pow_table *) b)->logk))
    return -1;

  if(((struct pow_table *) a)->logk > (((struct pow_table *) b)->logk))
    return +1;

  return 0;
}

void initialize_powerspectrum(void)
{
  double res;

  InitTime = 1 / (1 + Redshift);

  AA = 6.4 / ShapeGamma * (3.085678e24 / UnitLength_in_cm);
  BB = 3.0 / ShapeGamma * (3.085678e24 / UnitLength_in_cm);
  CC = 1.7 / ShapeGamma * (3.085678e24 / UnitLength_in_cm);
  nu = 1.13;

  R8 = 8 * (3.085678e24 / UnitLength_in_cm);	/* 8 Mpc/h */
  Norm = 1.0;
  Dplus=1.0;

  if(WhichSpectrum == 2)
    read_power_table();
  if(ReNormalizeInputSpectrum){
    res = TopHatSigma2(R8);
    if(WhichSpectrum == 2){
      printf("\nNormalization of spectrum in file:  Sigma8 = %g\n", sqrt(res));
    }

    Norm = Sigma8 * Sigma8 / res;

    printf("Normalization adjusted to  Sigma8=%g   (Normfac=%g)\n\n", Sigma8, Norm);
    Dplus = GrowthFactor(InitTime, 1.0);
    printf("Dplus initial redshift =%g  \n\n", Dplus);
  }
}


double PowerSpec_Tabulated(double k, int Type)
{
  double logk, logD, P, kold, u, dlogk, Delta2;
  int binlow, binhigh, binmid;

  kold = k;

  k *= (InputSpectrum_UnitLength_in_cm / UnitLength_in_cm);	/* convert to h/Mpc */

  logk = log10(k);

  if(logk < PowerTable[0].logk || logk > PowerTable[NPowerTable - 1].logk)
    return 0;

  binlow = 0;
  binhigh = NPowerTable - 1;

  while(binhigh - binlow > 1)
    {
      binmid = (binhigh + binlow) / 2;
      if(logk < PowerTable[binmid].logk)
      	binhigh = binmid;
      else
      	binlow = binmid;
    }

  dlogk = PowerTable[binhigh].logk - PowerTable[binlow].logk;

  if(dlogk == 0)
    exit(777);

  u = (logk - PowerTable[binlow].logk) / dlogk;
  
  /*Choose which value to return based on Type*/
  switch (Type){
          case 0:
                logD = (1 - u) * PowerTable[binlow].logDb + u * PowerTable[binhigh].logDb;
                break;
          case 1:
                logD = (1 - u) * PowerTable[binlow].logD + u * PowerTable[binhigh].logD;
                break;
          case 2:
                logD = (1 - u) * PowerTable[binlow].logD2nd + u * PowerTable[binhigh].logD2nd;
                break;
          default:
                logD = (1 - u) * PowerTable[binlow].logDtot + u * PowerTable[binhigh].logDtot;
  }

  Delta2 = pow(10.0, logD);

  P = Norm * Delta2 / (4 * M_PI * kold * kold * kold);

  //  printf(" k,P,u,logD,dlogk,binlow,binhigh,PowerTable[binlow].logD,PowerTable[binhigh].logD %g %g %g %g %g %g %d %d %g %g\n",k,P,u,logD,Delta2,dlogk,binlow,binhigh,PowerTable[binlow].logD,PowerTable[binhigh].logD);
  return P;
}


double PowerSpec_Efstathiou(double k)
{
  return Norm * k / pow(1 + pow(AA * k + pow(BB * k, 1.5) + CC * CC * k * k, nu), 2 / nu);
}


double PowerSpec_EH(double k)	/* Eisenstein & Hu */
{
  return Norm * k * pow(tk_eh(k), 2);
}


double tk_eh(double k)		/* from Martin White */
{
  double q, theta, ommh2, a, s, gamma, L0, C0;
  double tmp;
  double omegam, ombh2, hubble;

  /* other input parameters */
  hubble = HubbleParam;

  omegam = Omega;
  ombh2 = OmegaBaryon * HubbleParam * HubbleParam;

  if(OmegaBaryon == 0)
    ombh2 = 0.04 * HubbleParam * HubbleParam;

  k *= (3.085678e24 / UnitLength_in_cm);	/* convert to h/Mpc */

  theta = 2.728 / 2.7;
  ommh2 = omegam * hubble * hubble;
  s = 44.5 * log(9.83 / ommh2) / sqrt(1. + 10. * exp(0.75 * log(ombh2))) * hubble;
  a = 1. - 0.328 * log(431. * ommh2) * ombh2 / ommh2
    + 0.380 * log(22.3 * ommh2) * (ombh2 / ommh2) * (ombh2 / ommh2);
  gamma = a + (1. - a) / (1. + exp(4 * log(0.43 * k * s)));
  gamma *= omegam * hubble;
  q = k * theta * theta / gamma;
  L0 = log(2. * exp(1.) + 1.8 * q);
  C0 = 14.2 + 731. / (1. + 62.5 * q);
  tmp = L0 / (L0 + C0 * q * q);
  return (tmp);
} 

double TopHatSigma2(double R)
{
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
  double result,abserr;
  gsl_function F;
  F.function = &sigma2_int;
  F.params = &R;

  /* note: 500/R is here chosen as integration boundary (infinity) */
  gsl_integration_qags (&F, 0, 500./R, 0, 1e-4,1000,w,&result, &abserr);
//   printf("gsl_integration_qng in TopHatSigma2. Result %g, error: %g, intervals: %lu\n",result, abserr,w->size);
  gsl_integration_workspace_free (w);
  return result;
}

double sigma2_int(double k, void * params)
{
  double r = *(double *) params;
  double kr, kr2, w, x;
  kr = r * k;
  kr2 = kr * kr;

  /*Series expansion; actually good until kr~1*/
  if(kr < 1e-2)
      w = 1./3. - kr2/30. +kr2*kr2/840.;
  else
      w = 3 * (sin(kr) / kr - cos(kr)) / kr2;
  x = 4 * M_PI * k * k * w * w * PowerSpec(k,100);

  return x;
}


double GrowthFactor(double astart, double aend)
{
  return growth(aend) / growth(astart);
}

double growth(double a)
{
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (10);
  double hubble_a, Omegan=Omega;
  double result,abserr;
  gsl_function F;
  F.function = &growth_int;
  F.params = NULL;
  if(neutrinos_ks){
        Omegan = Omega + OmegaDM_2ndSpecies;
        printf("\n Omegan %g\n\n",Omegan);
  }
  hubble_a = sqrt(Omegan / (a * a * a) + (1 - Omegan - OmegaLambda) / (a * a) + OmegaLambda);
  gsl_integration_qags (&F, 0, a, 0, 1e-4,10,w,&result, &abserr);
//   printf("gsl_integration_qng in growth. Result %g, error: %g, intervals: %lu\n",result, abserr,w->size);
  gsl_integration_workspace_free (w);
  return hubble_a * result;
}


double growth_int(double a, void * param)
{
  double Omegan=Omega;
  if(neutrinos_ks)
        Omegan +=  OmegaDM_2ndSpecies;
  return pow(a / (Omegan + (1 - Omegan - OmegaLambda) * a + OmegaLambda * a * a * a), 1.5);
}


double F_Omega(double a)
{
  double omega_a,Omegan=Omega;
  if(neutrinos_ks)
        Omegan +=  OmegaDM_2ndSpecies;
  omega_a = Omegan / (Omegan + a * (1 - Omegan - OmegaLambda) + a * a * a * OmegaLambda);
  return pow(omega_a, 0.6);
}

double F2_Omega(double a)
{
  double omega_a;

  omega_a = Omega / (Omega + a * (1 - Omega - OmegaLambda) + a * a * a * OmegaLambda);

  return 2 * pow(omega_a, 4./7.);
}


/*  Here comes the stuff to compute the thermal WDM velocity distribution */


#define LENGTH_FERMI_DIRAC_TABLE 2000
#define MAX_FERMI_DIRAC          20.0

double fermi_dirac_vel[LENGTH_FERMI_DIRAC_TABLE];
double fermi_dirac_cumprob[LENGTH_FERMI_DIRAC_TABLE];

double WDM_V0 = 0;

double fermi_dirac_kernel(double x, void * param)
{
  return x * x / (exp(x) + 1);
}

void fermi_dirac_init(void)
{
  int i;

  gsl_integration_workspace * w = gsl_integration_workspace_alloc (10);
  double abserr;
  gsl_function F;
  F.function = &fermi_dirac_kernel;
  F.params = NULL;
  for(i = 0; i < LENGTH_FERMI_DIRAC_TABLE; i++)
    {
      fermi_dirac_vel[i] = MAX_FERMI_DIRAC * i / (LENGTH_FERMI_DIRAC_TABLE - 1.0);
      gsl_integration_qags (&F, 0, fermi_dirac_vel[i], 0, 1e-4,10,w,&(fermi_dirac_cumprob[i]), &abserr);
//       printf("gsl_integration_qng in fermi_dirac_init. Result %g, error: %g, intervals: %lu\n",fermi_dirac_cumprob[i], abserr,w->size);
    }
  gsl_integration_workspace_free (w);

  for(i = 0; i < LENGTH_FERMI_DIRAC_TABLE; i++)
    fermi_dirac_cumprob[i] /= fermi_dirac_cumprob[LENGTH_FERMI_DIRAC_TABLE - 1];

  WDM_V0 = 0.012 * (1 + Redshift) * pow((Omega - OmegaBaryon) / 0.3, 1.0 / 3) * pow(HubbleParam / 0.65,
										    2.0 / 3) * pow(1.0 /
												   WDM_PartMass_in_kev,
												   4.0 / 3);

    printf("\nWarm dark matter rms velocity dispersion at starting redshift = %g km/sec\n\n",
	   3.59714 * WDM_V0);

  WDM_V0 *= 1.0e5 / UnitVelocity_in_cm_per_s;

  /* convert from peculiar velocity to gadget's cosmological velocity */
  WDM_V0 *= sqrt(1 + Redshift);
}



double get_fermi_dirac_vel(void)
{
  int i;
  double p, u;

  p = drand48();
  i = 0;

  while(i < LENGTH_FERMI_DIRAC_TABLE - 2)
    if(p > fermi_dirac_cumprob[i + 1])
      i++;
    else
      break;

  u = (p - fermi_dirac_cumprob[i]) / (fermi_dirac_cumprob[i + 1] - fermi_dirac_cumprob[i]);

  return fermi_dirac_vel[i] * (1 - u) + fermi_dirac_vel[i + 1] * u;
}



void add_WDM_thermal_speeds(float *vel)
{
  double v, phi, theta, vx, vy, vz;

  if(WDM_V0 == 0)
    fermi_dirac_init();

  v = WDM_V0 * get_fermi_dirac_vel();

  phi = 2 * M_PI * drand48();
  theta = acos(2 * drand48() - 1);

  vx = v * sin(theta) * cos(phi);
  vy = v * sin(theta) * sin(phi);
  vz = v * cos(theta);

  vel[0] += vx;
  vel[1] += vy;
  vel[2] += vz;
}


#ifdef NEUTRINOS

double fermi_dirac_vel_nu[LENGTH_FERMI_DIRAC_TABLE];
double fermi_dirac_cumprob_nu[LENGTH_FERMI_DIRAC_TABLE];

double NU_V0 = 0;


void fermi_dirac_init_nu(void)
{
  int i;
  /*These functions are so smooth that we don't need much space*/
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (10);
  double abserr;
  gsl_function F;
  F.function = &fermi_dirac_kernel;
  F.params = NULL;
  for(i = 0; i < LENGTH_FERMI_DIRAC_TABLE; i++)
    {
      fermi_dirac_vel_nu[i] = MAX_FERMI_DIRAC * i / (LENGTH_FERMI_DIRAC_TABLE - 1.0);
      gsl_integration_qags (&F, 0, fermi_dirac_vel_nu[i], 0, 1e-4,10,w,&(fermi_dirac_cumprob_nu[i]), &abserr);
//       printf("gsl_integration_qng in fermi_dirac_init_nu. Result %g, error: %g, intervals: %lu\n",fermi_dirac_cumprob[i], abserr,w->size);
    }
  gsl_integration_workspace_free (w);

  for(i = 0; i < LENGTH_FERMI_DIRAC_TABLE; i++)
    fermi_dirac_cumprob_nu[i] /= fermi_dirac_cumprob_nu[LENGTH_FERMI_DIRAC_TABLE - 1];

  NU_V0 = 150.0 * (1.0e5 / UnitVelocity_in_cm_per_s) * (1 + Redshift) * (1.0 / NU_PartMass_in_ev);

  printf("\nNeutrino rms vel. dispersion %g (int vel units)\n\n",NU_V0);

  /* convert from peculiar velocity to gadget's cosmological velocity */
  NU_V0 *= sqrt(1 + Redshift);
}

double get_fermi_dirac_vel_nu(void)
{
  int i;
  double p, u;
  int binlow = 0;
  int binhigh = LENGTH_FERMI_DIRAC_TABLE - 2;


  p = drand48();
  i = 0;

  while(binhigh - binlow > 1)
    {
      i = (binhigh + binlow) / 2;
      if(p > fermi_dirac_cumprob_nu[i + 1])
        binlow = i;
      else
        binhigh = i;
    }

  u = (p - fermi_dirac_cumprob_nu[i]) / (fermi_dirac_cumprob_nu[i + 1] - fermi_dirac_cumprob_nu[i]);

  return fermi_dirac_vel_nu[i] * (1 - u) + fermi_dirac_vel_nu[i + 1] * u;
}

void add_NU_thermal_speeds(float *vel)
{
  double v, phi, theta, vx, vy, vz;

  if(NU_V0 == 0)
    fermi_dirac_init_nu();

  v = NU_V0 * get_fermi_dirac_vel_nu();

  phi = 2 * M_PI * drand48();
  theta = acos(2 * drand48() - 1);

  vx = v * sin(theta) * cos(phi);
  vy = v * sin(theta) * sin(phi);
  vz = v * cos(theta);

  vel[0] += vx;
  vel[1] += vy;
  vel[2] += vz;
}

#endif

