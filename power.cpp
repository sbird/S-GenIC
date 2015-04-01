#include "power.hpp"
#include <math.h>
#include <cstdio>
//For the normalized power spectrum
#include <gsl/gsl_integration.h>

//Define the various functions for the tabulated power spectrum
int compare_logk(const void *a, const void *b)
{
  if(((struct pow_table *) a)->logk < (((struct pow_table *) b)->logk))
    return -1;

  if(((struct pow_table *) a)->logk > (((struct pow_table *) b)->logk))
    return +1;

  return 0;
}


PowerSpec_Tabulated::PowerSpec_Tabulated(char * FileWithTransfer, char * FileWithInputSpectrum, double Omega, double OmegaLambda, double OmegaBaryon, double OmegaNu,
                        double InputSpectrum_UnitLength_in_cm, double UnitLength_in_cm, bool no_gas, bool neutrinos_ks)
{
  //Set up conversion factor between internal units and CAMB units
  scale = (InputSpectrum_UnitLength_in_cm / UnitLength_in_cm);
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
      double T_b, dummy, T_nu, T_nu2,T_tot, T_cdm;
      double delta_cdm, delta_nu, delta_nu_2nd,delta_tot, delta_b;

      if(fscanf(fd, " %lg %lg %lg %lg %lg %lg %lg", &k, &T_cdm, &T_b, &dummy, &T_nu2, &T_nu, &T_tot) == 7)
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
           * optionally subtract the baryons.
           * We want to incorporate neutrinos into the dark matter if we are changing the transfer function,
           * but not if we are using the full kspace method. */
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
//                     T_cdm = (Omega*T_tot - T_nu*OmegaDM_2ndSpecies)/(Omega-OmegaDM_2ndSpecies);
                T_cdm = (T_cdm*(Omega-OmegaNu - OmegaBaryon) + T_b*OmegaBaryon)/(Omega-OmegaNu);
          }
          /*This should be equivalent to the above in almost all cases,
           * but perhaps we want to see the effect of changing only the ICs in CAMB for the neutrinos.*/
          if(no_gas && OmegaNu == 0){
                  T_cdm = T_tot;
          }
	  /* obtain P(k) from transfer function ratios like suggested by JL*/
	  /*NOTE for this to work CAMB's transfer_k_interp_matterpower must be off!!*/
	  delta_b   = k * k * k * pow(T_b/T_tot,2)* PowerMatter[NPowerTable].pmat/(2*M_PI*M_PI);
	  delta_cdm = k * k * k * pow(T_cdm/T_tot,2)* PowerMatter[NPowerTable].pmat/(2*M_PI*M_PI);
	  delta_nu = k * k * k * pow(T_nu/T_tot,2) * PowerMatter[NPowerTable].pmat/(2*M_PI*M_PI); 
	  delta_nu_2nd = k * k * k * pow(T_nu2/T_tot,2) * PowerMatter[NPowerTable].pmat/(2*M_PI*M_PI); 
	  delta_tot = k * k * k * PowerMatter[NPowerTable].pmat/(2*M_PI*M_PI);

          PowerTable[NPowerTable].logD = log10(delta_cdm);
	  // printf("NT,d_cdm,log10(d_cdm),k %d %g %g %g \n",NPowerTable,delta_cdm,log10(delta_cdm),k);
	  PowerTable[NPowerTable].logD2nd = log10(delta_nu);
	  PowerTable[NPowerTable].logD3rd = log10(delta_nu_2nd);
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

double PowerSpec_Tabulated::power(double k, int Type)
{
  double logk, logD, P, kold, u, dlogk, Delta2;
  int binlow, binhigh, binmid;

  kold = k;

  k *= scale;	/* convert to h/Mpc */

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
          case 3:
                logD = (1 - u) * PowerTable[binlow].logD3rd + u * PowerTable[binhigh].logD3rd;
                break;
          default:
                logD = (1 - u) * PowerTable[binlow].logDtot + u * PowerTable[binhigh].logDtot;
  }

  Delta2 = pow(10.0, logD);

  P = Delta2 / (4 * M_PI * kold * kold * kold);
  //  printf(" k,P,u,logD,dlogk,binlow,binhigh,PowerTable[binlow].logD,PowerTable[binhigh].logD %g %g %g %g %g %g %d %d %g %g\n",k,P,u,logD,Delta2,dlogk,binlow,binhigh,PowerTable[binlow].logD,PowerTable[binhigh].logD);
  return P;
}

//Definitions for the normalized power spectrum
//These mostly compute sigma_8 via the GSL
NormalizedPowerSpec::NormalizedPowerSpec(PowerSpec * PSpec, double Sigma8, double PrimordialIndex, double Dplus, double UnitLength_in_cm): PSpec(PSpec), PrimordialIndex(PrimordialIndex), Dplus(Dplus)
{
        R8 = 8 * (3.085678e24 / UnitLength_in_cm);	/* 8 Mpc/h */
        //Uses R8
        Norm = Sigma8 * Sigma8 / TopHatSigma2();
        printf("Normalization adjusted to  Sigma8=%g   (Normfac=%g)\n\n", Sigma8, Norm);
//         Cosmology cosmo(HubbleParam, Omega, OmegaLambda, MNu, InvertedHierarchy);
//         Dplus = cosmo.GrowthFactor(InitTime, 1.0);
        printf("Dplus initial redshift =%g  \n\n", Dplus);
}

double sigma2_int(double k, void * params)
{
  NormalizedPowerSpec * PSpec = (NormalizedPowerSpec *) params;
  double r = PSpec->R8;
  double kr, kr2, w, x;
  kr = r * k;
  kr2 = kr * kr;

  /*Series expansion; actually good until kr~1*/
  if(kr < 1e-2)
      w = 1./3. - kr2/30. +kr2*kr2/840.;
  else
      w = 3 * (sin(kr) / kr - cos(kr)) / kr2;
  x = 4 * M_PI * k * k * w * w * PSpec->power(k,100);
  return x;
}

double NormalizedPowerSpec::TopHatSigma2()
{
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
  double result,abserr;
  gsl_function F;
  F.function = &sigma2_int;
  F.params = this;

  /* note: 500/R is here chosen as integration boundary (infinity) */
  gsl_integration_qags (&F, 0, 500./this->R8, 0, 1e-4,1000,w,&result, &abserr);
//   printf("gsl_integration_qng in TopHatSigma2. Result %g, error: %g, intervals: %lu\n",result, abserr,w->size);
  gsl_integration_workspace_free (w);
  return result;
}

