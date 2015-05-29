#include "power.hpp"
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <limits>
#include <cassert>
//For the normalized power spectrum
#include <gsl/gsl_integration.h>

PowerSpec_Tabulated::PowerSpec_Tabulated(const char * FileWithTransfer, const char * FileWithInputSpectrum, double Omega, double OmegaLambda, double OmegaBaryon, double OmegaNu,
                        double InputSpectrum_UnitLength_in_cm, double UnitLength_in_cm, bool no_gas, bool neutrinos_ks)
{
  //Set up conversion factor between internal units and CAMB units
  scale = (InputSpectrum_UnitLength_in_cm / UnitLength_in_cm);
  std::vector<double> kmatter;
  std::vector<double> pmatter;

  std::fstream matpow;
  matpow.open(FileWithInputSpectrum, std::fstream::in);
  if ( ! matpow.is_open() ) {
      std::cerr<<"Can't open matter power spectrum in file: "<< FileWithInputSpectrum<<std::endl;
      exit(17);
  }

  /* define matter array */
  while ( matpow.good() ) {
    double ktmp,ptmp;
    matpow >> ktmp;
    matpow >> ptmp;
    if (!matpow.good() )
        break;
    kmatter.push_back(ktmp);
    pmatter.push_back(ptmp);
  }
  std::cerr<<"Found "<<pmatter.size()<<" rows in input spectrum table\n"<<std::endl;
  matpow.close();

  std::fstream transfer;
  transfer.open(FileWithTransfer, std::fstream::in);
  if ( ! transfer.is_open() ) {
      std::cerr<<"Can't open transfer function in file: "<< FileWithTransfer<<std::endl;
      exit(17);
  }

  //Temporary arrays which can be resized so we don't have to read the file twice
  std::vector<double> ktransfer_vector;
  std::vector<double> transfer_vector[N_TYPES];

  while ( transfer.good() ) {
      double T_b, dummy, T_nu, T_nu2,T_tot, T_cdm;
      double k;
      double delta_cdm, delta_nu, delta_nu_2nd,delta_tot, delta_b;
      //Read the row.
      transfer >> k;
      transfer >> T_cdm;
      transfer >> T_b;
      //radiation
      transfer >> dummy;
      transfer >> T_nu2;
      transfer >> T_nu;
      transfer >> T_tot;
      //Ignore the rest of the line
      transfer.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
      if ( ! transfer.good() )
          break;
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
//                 T_cdm = (Omega*T_tot - T_nu*OmegaDM_2ndSpecies)/(Omega-OmegaDM_2ndSpecies);
            T_cdm = (T_cdm*(Omega-OmegaNu - OmegaBaryon) + T_b*OmegaBaryon)/(Omega-OmegaNu);
      }
      /*This should be equivalent to the above in almost all cases,
       * but perhaps we want to see the effect of changing only the ICs in CAMB for the neutrinos.*/
      if(no_gas && OmegaNu == 0){
              T_cdm = T_tot;
      }
      /* obtain P(k) from transfer function ratios like suggested by JL*/
      /*NOTE for this to work CAMB's transfer_k_interp_matterpower must be off!!*/
      size_t cursize = ktransfer_vector.size();
      delta_b   = pow(T_b/T_tot,2) * pmatter[cursize]/(2*M_PI*M_PI);
      delta_cdm = pow(T_cdm/T_tot,2)* pmatter[cursize]/(2*M_PI*M_PI);
      delta_nu = pow(T_nu/T_tot,2) * pmatter[cursize]/(2*M_PI*M_PI);
      delta_nu_2nd = pow(T_nu2/T_tot,2) * pmatter[cursize]/(2*M_PI*M_PI);
      delta_tot = pmatter[cursize]/(2*M_PI*M_PI);
      //Assign row to structures
      ktransfer_vector.push_back(k);
      transfer_vector[0].push_back(delta_b);
      transfer_vector[1].push_back(delta_cdm);
      transfer_vector[2].push_back(delta_nu);
      transfer_vector[3].push_back(delta_nu_2nd);
      transfer_vector[4].push_back(delta_tot);

  }
  transfer.close();

  //Do some basic checks on the read input files.
  if (ktransfer_vector.size() != pmatter.size() ) {
      std::cerr << "Transfer function was "<<ktransfer_vector.size()<<" lines long, but input power spectrum was "<<pmatter.size()<<std::endl;
      exit(47);
  }
  //Allocate arrays for the gsl interpolation (which must be raw arrays, unfortunately)
  //Then allocate the interpolation objects
  ktransfer_table = new double[ktransfer_vector.size()];
  //Copy over the k values
  for (size_t i=0; i<ktransfer_vector.size(); i++) {
        ktransfer_table[i] = log10(ktransfer_vector[i]);
        if(fabs(kmatter[i]- ktransfer_vector[i]) > 0.01 * kmatter[i]) {
              fprintf(stderr, "Error: Input spectrum row %ld has k=%g, transfer has k=%g.\n",i,kmatter[i],ktransfer_vector[i]);
              std::cerr<<"Remember you need transfer_k_interp_matterpower = F in CAMB"<<std::endl;
              exit(47);
        }
        if(i >0 && ktransfer_table[i] <= ktransfer_table[i-1]) {
              std::cerr<<"Transfer table is not increasing in k: i = "<<i<<", k = "<<ktransfer_table[i]<<" <= "<<ktransfer_table[i-1]<<std::endl;
              exit(48);
        }
  }
  for(int type=0; type < N_TYPES; type++){
      //Check everything is the same size
      assert( ktransfer_vector.size() == transfer_vector[type].size());
      transfer_table[type] = new double[transfer_vector[type].size()];
      //Copy data from the temporary vectors to the statically sized C-arrays
      for (size_t i=0; i<transfer_vector[type].size(); i++) {
          //Make sure we do not take log(0)
          transfer_table[type][i] = log10(std::max(transfer_vector[type][i],1e-90));
      }
      //Set up the interpolation structures
      trans_interp[type] = gsl_interp_alloc(gsl_interp_cspline,transfer_vector[type].size());
      trans_interp_accel[type] = gsl_interp_accel_alloc();
      //Initialise the interpolator; same k values for each one, different T values.
      //Size asserted to be the same.
      gsl_interp_init(trans_interp[type],ktransfer_table, transfer_table[type],transfer_vector[type].size());
  }
  //Store the size of the arrays
  NPowerTable = ktransfer_vector.size();
}

double PowerSpec_Tabulated::power(double k, int Type)
{

  double logk = log10(k*scale);

  if(logk < ktransfer_table[0] || logk > ktransfer_table[NPowerTable - 1])
    return 0;

  //Clamp the requested type.
  if (Type > N_TYPES-1 || Type < 0)
      Type = N_TYPES-1;

  double logD = gsl_interp_eval(trans_interp[Type], ktransfer_table, transfer_table[Type], logk, trans_interp_accel[Type]);

  double Delta2 = pow(10.0, logD);

  double P = Delta2 * scale * scale * scale / (4 * M_PI );
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

