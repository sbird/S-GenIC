#include <math.h>
#include "allvars.h"
#include "proto.h"
#include <map>

static double R8;
static double r_tophat;

static double AA, BB, CC;
static double nu;
static double Norm;

/*This stores the conversion between the tables and the internal units. By default 1e-3*/
static double kctog;
/*This prints the CAMB transfer function*/
double tk_CAMB(double k, int Type);
/*Power spectra*/
double PowerSpec_CAMB(double k, int Type);
static int NPowerTable;
#define APRIM 2.43e-9
#define PIVOT_SCALE (0.05*kctog/HubbleParam)

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

 /*Structure for transfer table*/
struct trans_row{
	double T_CDM;
	double T_b;
	double T_g;
	double T_r;
	double T_n;
	double T_t;
};

static std::map<double, struct trans_row> transfer_tables;

double PowerSpec(double k, int Type)
{
  double power, alpha, Tf;
  switch (WhichSpectrum)
    {
    case 1:
      power = PowerSpec_EH(k);
      break;
    case 2:
      if (Type == 2)
#ifdef NEUTRINOS
        power = PowerSpec_Tabulated2nd(k);
#else
        power = PowerSpec_DM_2ndSpecies(k);
#endif
      else if(Type == 0)
          power = PowerSpec_Tabulated_b(k);
      else
          power = PowerSpec_Tabulated(k);
      break;
#ifdef DIFFERENT_TRANSFER_FUNC
    case 3:
      /*This is my reimplementation of Spectrum 2. */
  /*ADD THE FACTOR OF (2π)^3 to convert from CAMB conventions to GADGET conventions!!*/
      power = PowerSpec_CAMB(k,Type)/pow(2*M_PI,3);
      break;
#endif
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

#if defined(DIFFERENT_TRANSFER_FUNC)
/* Type 100 recomputes normalization over the total matter p(k) */
  if(Type == 100 && WhichSpectrum ==2) 
	power = PowerSpec_TOTAL(k);
#endif 
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


/*Read the transfer tables from CAMB*/
void read_transfer_table(void)
{
  /* Transfer file format:  
   * k/h Delta_CDM/k2 Delta_b/k2 Delta_g/k2 Delata_r/k2 Delta_nu/k2 Delta_tot/k2
   *     CDM, baryon,photon,massless neutrino, massive neutrinos, and total (massive)*/
	FILE *trans;
	char buf[200];
	struct trans_row tmp_row;
        double tmp_k;
        kctog=UnitLength_in_cm/InputSpectrum_UnitLength_in_cm;
	/*Open file*/
	sprintf(buf,FileWithTransfer);
	if(!(trans=fopen(buf,"r"))) {
		fprintf(stderr,"Can't open transfer file %s! You probably forgot to create it!\n",buf);
		exit(17);
	} 
	/*Work out how many lines in file*/
	if(!transfer_tables.empty()) {
		fprintf(stderr,"read_CAMB_tables has been called more than once! Failing.");
		exit(19);
	}
	/*Read line by line.*/
	while(fscanf(trans," %lg %lg %lg %lg %lg %lg %lg\n",&tmp_k,&tmp_row.T_CDM,&tmp_row.T_b, &tmp_row.T_g,&tmp_row.T_r,&tmp_row.T_n,&tmp_row.T_t)==7)
	{
		/* k needs to go from (h/Mpc) units to internal Gadget units (h/kpc) by default.
                 * kctog is by default 1e-3 */
		tmp_k *= kctog;
		/*Append line to table.*/
		transfer_tables[tmp_k]=tmp_row;
		if(feof(trans))
			break;
	}
	if(ferror(trans)) {
		fprintf(stderr,"Error reading transfer file: %d",errno);
		exit(21);
	}
        printf("Found %d rows in input CAMB transfer file\n",(int)transfer_tables.size());
	fclose(trans);
	/*The CAMB T_f/k_c is in units of Mpc^2! NOTE NO h!*/
	double tctog=(HubbleParam*HubbleParam)/(kctog*kctog);
        std::map<double, struct trans_row>::iterator it;
	for(it=transfer_tables.begin(); it != transfer_tables.end(); ++it) {
	/*The transfer function should be normalized to about 1 on large scales.*/
		((*it).second).T_CDM *= tctog;
		((*it).second).T_g *= tctog;
		((*it).second).T_b *= tctog;
		((*it).second).T_r *= tctog;
		((*it).second).T_n *= tctog;
		((*it).second).T_t *= tctog;
	}
	return;
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
      double T_cdmnew, T_b, dummy, T_nu, T_tot, T_cdm;
      double delta_cdm, delta_nu, delta_tot, delta_b;

      if(fscanf(fd, " %lg %lg %lg %lg %lg %lg %lg", &k, &T_cdm, &T_b, &dummy, &dummy, &T_nu, &T_tot) == 7)
	{
	  PowerTable[NPowerTable].logk = log10(k);

	  /* obtain P(k) from transfer function ratios like suggested by JL*/
	  
	  delta_b   = k * k * k * pow(T_b/T_tot,2)* PowerMatter[NPowerTable].pmat/(2*M_PI*M_PI);
	  delta_cdm = k * k * k * pow(T_cdm/T_tot,2)* PowerMatter[NPowerTable].pmat/(2*M_PI*M_PI);
	  delta_nu = k * k * k * pow(T_nu/T_tot,2) * PowerMatter[NPowerTable].pmat/(2*M_PI*M_PI); 
	  delta_tot = k * k * k * PowerMatter[NPowerTable].pmat/(2*M_PI*M_PI);

	  /* delta_nu =  k * k * k * PowerMatter[NPowerTable].pmat/(2*M_PI*M_PI);
	     delta_cdm =  k * k * k * PowerMatter[NPowerTable].pmat/(2*M_PI*M_PI); */
	  // if there are no baryons and no neutrinos then delta_cdm=delta_tot
	  if (OmegaBaryon == 0  && OmegaDM_2ndSpecies == 0)
            {
	      delta_cdm = delta_tot;

	      }  

	  // if there are no baryons but there is a nu component fake the cdm to have the baryon contribution as well
	  if (OmegaBaryon == 0  && OmegaDM_2ndSpecies != 0)
            {
	      T_cdmnew = (0.05*T_b+(Omega-0.05-OmegaDM_2ndSpecies)*T_cdm)/(Omega-OmegaDM_2ndSpecies);
              delta_cdm = k * k * k * pow(T_cdmnew/T_tot,2)* PowerMatter[NPowerTable].pmat/(2*M_PI*M_PI);
	      // printf(" check table %g %g %g \n",T_cdm,T_cdmnew,delta_cdm);
	    }

          if(neutrinos_ks){
	    // if there are no baryons but there is a nu component fake the cdm to have the baryon contribution as well
	    if (OmegaBaryon == 0  && OmegaDM_2ndSpecies != 0)
              {
	        T_cdmnew = (0.05*T_b+(Omega-0.05-OmegaDM_2ndSpecies)*T_cdm+OmegaDM_2ndSpecies*T_nu)/(Omega-OmegaDM_2ndSpecies-0.05);
                delta_cdm = k * k * k * pow(T_cdmnew/T_tot,2)* PowerMatter[NPowerTable].pmat/(2*M_PI*M_PI);
	        // printf(" check table %g %g %g \n",T_cdm,T_cdmnew,delta_cdm);
	      }
	    // if there are  baryons and a fake nu component fake the cdm to have the nu contribution in as well
            if (OmegaBaryon != 0  && OmegaDM_2ndSpecies != 0)
              {
                T_cdmnew = (OmegaDM_2ndSpecies*T_nu+(Omega-0.05-OmegaDM_2ndSpecies)*T_cdm)/(Omega-0.05);
                delta_cdm = k * k * k * pow(T_cdmnew/T_tot,2)* PowerMatter[NPowerTable].pmat/(2*M_PI*M_PI);
                // printf(" check table %g %g %g \n",T_cdm,T_cdmnew,delta_cdm);
              }
          }

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
  if(WhichSpectrum > 2)
    read_transfer_table();
  if(ReNormalizeInputSpectrum){
    res = TopHatSigma2(R8);
    if(WhichSpectrum == 2){
      printf("\nNormalization of spectrum in file:  Sigma8 = %g\n", sqrt(res));
    }

    Norm = Sigma8 * Sigma8 / res;

    if(WhichSpectrum == 2)
      printf("Normalization adjusted to  Sigma8=%g   (Normfac=%g)\n\n", Sigma8, Norm);
      Dplus = GrowthFactor(InitTime, 1.0);
      printf("Dplus initial redshift =%g  \n\n", Dplus); 
  }
}


double PowerSpec_Tabulated(double k)
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

  logD = (1 - u) * PowerTable[binlow].logD + u * PowerTable[binhigh].logD;

  Delta2 = pow(10.0, logD);

  P = Norm * Delta2 / (4 * M_PI * kold * kold * kold);

  //  printf(" k,P,u,logD,dlogk,binlow,binhigh,PowerTable[binlow].logD,PowerTable[binhigh].logD %g %g %g %g %g %g %d %d %g %g\n",k,P,u,logD,Delta2,dlogk,binlow,binhigh,PowerTable[binlow].logD,PowerTable[binhigh].logD);
  return P;
}


double PowerSpec_Tabulated_b(double k)
{
  double logk, logD, P, kold, u, dlogk, Delta2;
  int binlow, binhigh, binmid;

  kold = k;

  k *= (InputSpectrum_UnitLength_in_cm / UnitLength_in_cm);     /* convert to h/Mpc */

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

  logD = (1 - u) * PowerTable[binlow].logDb + u * PowerTable[binhigh].logDb;

  Delta2 = pow(10.0, logD);

  P = Norm * Delta2 / (4 * M_PI * kold * kold * kold);

return P;
}



#ifdef NEUTRINOS

double PowerSpec_Tabulated2nd(double k)
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

  logD = (1 - u) * PowerTable[binlow].logD2nd + u * PowerTable[binhigh].logD2nd;

  Delta2 = pow(10.0, logD);

  P = Norm * Delta2 / (4 * M_PI * kold * kold * kold);

  return P;
}
#endif



double PowerSpec_TOTAL(double k)
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

  logD = (1 - u) * PowerTable[binlow].logDtot + u * PowerTable[binhigh].logDtot;

  Delta2 = pow(10.0, logD);

  P = Norm * Delta2 / (4 * M_PI * kold * kold * kold);

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


double PowerSpec_CAMB(double k, int Type)
{
	return APRIM*2*M_PI*M_PI*k*pow(k/PIVOT_SCALE,PrimordialIndex-1.0)*pow(tk_CAMB(k, Type),2);
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

/*Return interpolated value of transfer function from table*/
double tk_CAMB(double k, int Type)
{
        std::map<double, struct trans_row>::iterator lit, uit;
        struct trans_row * lrow, * urow;
	double T1,T2,k1,k2;
	double tkout;
	if(transfer_tables.empty())
	{
		fprintf(stderr, "Some kind of error; tables not initialized!\n");
		exit(18);
	}
	uit=transfer_tables.upper_bound(k);

	/*No power outside of our boundaries.*/
	if(uit == transfer_tables.end() || uit == transfer_tables.begin())
		return 0;
        lit = uit;
        --lit;
        lrow=&((*lit).second);
        urow=&((*uit).second);

	/*Linear interpolation. Different transfer functions used for baryons and DM*/
#if defined(DIFFERENT_TRANSFER_FUNC)
        /* CDM: The dark matter may incorporate other particle types as well, eg,
         * fake neutrinos, or baryons.
         * NOTE that CAMB defines T_tot = (T_CDM M_CDM + T_b M_b +T_nu M_nu) / (M_CDM + M_b +M_nu)
         * HOWEVER, when relativistic (ie, in the early universe), neutrinos will redshift
         * like radiation. The CAMB transfer function takes this into account when calculating T_tot,
         * so instead of summing the transfer functions, use the total transfer function and
         * optionally subtract the baryons.*/
        if(Type==1){
                /* If we have fake neutrinos, use the total matter power spectrum. */
                if(neutrinos_ks){
                        T1 = lrow->T_t;
                        T2 = urow->T_t;
                        /*If we have separate gas particles, subtract
                         * the baryon transfer function */
                        if(!no_gas){
                                T1-=lrow->T_b*OmegaBaryon/Omega;
                                T2-=urow->T_b*OmegaBaryon/Omega;
                        }
                }
                /*No fake neutrinos*/
                else{
                        const double OmegaCDM = Omega -OmegaDM_2ndSpecies -OmegaBaryon;
                        double Omega_t = OmegaCDM;
                        T1=lrow->T_CDM*OmegaCDM;
                        T2=urow->T_CDM*OmegaCDM;
                        /* If we have no gas, add T_b M_b to the CDM */
                        if(no_gas){
                                T1 += lrow->T_b*OmegaBaryon;
                                T2 += urow->T_b*OmegaBaryon;
                                Omega_t +=OmegaBaryon;
                        }
                        T1/=Omega_t;
                        T2/=Omega_t;
                }
        }
        /* Baryons */
        else if(Type==0){
        	T1=(*lrow).T_b;
        	T2=(*urow).T_b;
        }
#ifdef NEUTRINOS
        /* This loads the massive neutrino type*/
        else if(Type == 2){
        	T1=(*lrow).T_n;
        	T2=(*urow).T_n;
        }
#endif //NEUTRINOS
        else{
#endif //DIFFERENT_TRANSFER_FUNCTION
        	T1=(*lrow).T_t;
        	T2=(*urow).T_t;
#if defined(DIFFERENT_TRANSFER_FUNC)
        }
#endif

	k1=(*lit).first;
	k2=(*uit).first;
	//Do it in log space!
	tkout=exp((log(T2)*(log(k)-log(k1))+log(T1)*(log(k2)-log(k)))/(log(k2)-log(k1)));
	return tkout;
}

double TopHatSigma2(double R)
{
  r_tophat = R;

  return qromb(sigma2_int, 0, 500.0 * 1 / R);	/* note: 500/R is here chosen as 
						   integration boundary (infinity) */
}


double sigma2_int(double k)
{
  double kr, kr3, kr2, w, x;

  kr = r_tophat * k;
  kr2 = kr * kr;
  kr3 = kr2 * kr;

  if(kr < 1e-8)
    return 0;

  w = 3 * (sin(kr) / kr3 - cos(kr) / kr2);
  x = 4 * PI * k * k * w * w * PowerSpec(k,100);

  return x;
}


double GrowthFactor(double astart, double aend)
{
  return growth(aend) / growth(astart);
}


double growth(double a)
{
  double hubble_a, Omegan=Omega;

  if(neutrinos_ks){
        Omegan = Omega + OmegaDM_2ndSpecies;
        printf("\n Omegan %g\n\n",Omegan);
  }
  hubble_a = sqrt(Omegan / (a * a * a) + (1 - Omegan - OmegaLambda) / (a * a) + OmegaLambda);
  return hubble_a * qromb(growth_int, 0, a);
}


double growth_int(double a)
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

/*  Here comes the stuff to compute the thermal WDM velocity distribution */


#define LENGTH_FERMI_DIRAC_TABLE 2000
#define MAX_FERMI_DIRAC          20.0

double fermi_dirac_vel[LENGTH_FERMI_DIRAC_TABLE];
double fermi_dirac_cumprob[LENGTH_FERMI_DIRAC_TABLE];

double WDM_V0 = 0;

double fermi_dirac_kernel(double x)
{
  return x * x / (exp(x) + 1);
}

void fermi_dirac_init(void)
{
  int i;

  for(i = 0; i < LENGTH_FERMI_DIRAC_TABLE; i++)
    {
      fermi_dirac_vel[i] = MAX_FERMI_DIRAC * i / (LENGTH_FERMI_DIRAC_TABLE - 1.0);
      fermi_dirac_cumprob[i] = qromb(fermi_dirac_kernel, 0, fermi_dirac_vel[i]);
    }

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

  for(i = 0; i < LENGTH_FERMI_DIRAC_TABLE; i++)
    {
      fermi_dirac_vel_nu[i] = MAX_FERMI_DIRAC * i / (LENGTH_FERMI_DIRAC_TABLE - 1.0);
      fermi_dirac_cumprob_nu[i] = qromb(fermi_dirac_kernel, 0, fermi_dirac_vel_nu[i]);
    }

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

  p = drand48();
  i = 0;

  while(i < LENGTH_FERMI_DIRAC_TABLE - 2)
    if(p > fermi_dirac_cumprob_nu[i + 1])
      i++;
    else
      break;

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

