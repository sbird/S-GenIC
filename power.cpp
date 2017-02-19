#include "power.hpp"
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <cstddef>
#include <limits>
#include <cassert>
//For the normalized power spectrum
#include <gsl/gsl_integration.h>

PowerSpec_Tabulated::PowerSpec_Tabulated(const std::string& FileWithTransfer, const std::string & FileWithInputSpectrum, double Omega, double OmegaLambda, double OmegaBaryon, double OmegaNu,
                        double InputSpectrum_UnitLength_in_cm, double UnitLength_in_cm, bool no_gas, bool combined_neutrinos)
{
  //Set up conversion factor between internal units and CAMB units
  scale = (InputSpectrum_UnitLength_in_cm / UnitLength_in_cm);
  std::vector<double> kmatter_vector;
  std::vector<double> pmatter_vector;

  std::fstream matpow;
  matpow.open(FileWithInputSpectrum, std::fstream::in);
  if ( ! matpow.is_open() ) {
      std::cerr<<"Can't open matter power spectrum in file: "<< FileWithInputSpectrum<<std::endl;
      exit(17);
  }

  /* define matter array */
  while ( matpow.good() ) {
    double ktmp,ptmp;
    //Discard comments
    if (matpow.get() == '#')
        matpow.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
    else
        matpow.unget();
    matpow >> ktmp;
    matpow >> ptmp;
    if (!matpow.good() )
        break;
    kmatter_vector.push_back(ktmp);
    pmatter_vector.push_back(ptmp);
  }
  std::cerr<<"Found "<<pmatter_vector.size()<<" rows in input spectrum table: "<<FileWithInputSpectrum<<std::endl;
  matpow.close();

  //Allocate arrays for the gsl interpolation (which must be raw arrays, unfortunately)
  //Then allocate the interpolation objects
  kmatter_table = new double[kmatter_vector.size()];
  pmatter_table = new double[pmatter_vector.size()];
  assert(kmatter_vector.size() == pmatter_vector.size());
  //Copy over the k values and convert Fourier convention
  for (size_t i=0; i<kmatter_vector.size(); i++) {
        kmatter_table[i] = log10(kmatter_vector[i]);
        pmatter_table[i] = log10(pmatter_vector[i]/pow(2*M_PI,3));
        if(i >0 && kmatter_table[i] <= kmatter_table[i-1]) {
              std::cerr<<"Matter power table is not increasing in k: i = "<<i<<", k = "<<kmatter_table[i]<<" <= "<<kmatter_table[i-1]<<std::endl;
              exit(48);
        }
  }
  //Set up the interpolation structures for the matter power
  pmat_interp = gsl_interp_alloc(gsl_interp_cspline,kmatter_vector.size());
  pmat_interp_accel = gsl_interp_accel_alloc();
  //Initialise the interpolator; same k values for each one, different T values.
  gsl_interp_init(pmat_interp,kmatter_table, pmatter_table,pmatter_vector.size());
  NPowerTable = pmatter_vector.size();

  std::fstream transfer;
  transfer.open(FileWithTransfer, std::fstream::in);
  if ( ! transfer.is_open() ) {
      std::cerr<<"Can't open transfer function in file: "<< FileWithTransfer<<std::endl;
      exit(17);
  }

  //Temporary arrays which can be resized so we don't have to read the file twice
  std::vector<double> ktransfer_vector;
  std::vector<double> transfer_vector[N_TYPES_TAB];

  while ( transfer.good() ) {
      double T_b, dummy, T_nu, T_nu2,T_tot, T_cdm, T_dmb;
      double k;
      //Discard comments
      if (transfer.get() == '#')
          transfer.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
      else
          transfer.unget();
      //Read the row.
      transfer >> k;
      transfer >> T_cdm;
      transfer >> T_b;
      //radiation
      transfer >> dummy;
      //massless neutrinos
      transfer >> T_nu2;
      //massive neutrinos
      transfer >> T_nu;
      transfer >> T_tot;
      //DM+baryons
      transfer >> T_dmb;
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
      if(combined_neutrinos){
            T_cdm = T_tot;
            /*If we have separate gas particles, subtract
             * the baryon transfer function */
            if(!no_gas){
                T_cdm = (Omega*T_tot - T_b*OmegaBaryon)/(Omega-OmegaBaryon);
            }
      }
      /*Add baryons to the CDM if there are no gas particles*/
      else if(no_gas){
            T_cdm = T_dmb;
      }
      /*This should be equivalent to the above in almost all cases,
       * but perhaps we want to see the effect of changing only the ICs in CAMB for the neutrinos.*/
      if(no_gas && OmegaNu == 0){
              T_cdm = T_tot;
      }
      //Assign transfer functions to structures as T_s/T_t^2
      ktransfer_vector.push_back(k);
      transfer_vector[0].push_back(pow(T_b/T_tot,2));
      transfer_vector[1].push_back(pow(T_cdm/T_tot,2));
      transfer_vector[2].push_back(pow(T_nu/T_tot,2));
      transfer_vector[3].push_back(pow(T_nu2/T_tot,2));
  }
  transfer.close();

  //Allocate arrays for the gsl interpolation (which must be raw arrays, unfortunately)
  //Then allocate the interpolation objects
  ktransfer_table = new double[ktransfer_vector.size()];
  //Copy over the k values
  for (size_t i=0; i<ktransfer_vector.size(); i++) {
        ktransfer_table[i] = log10(ktransfer_vector[i]);
        if(i >0 && ktransfer_table[i] <= ktransfer_table[i-1]) {
              std::cerr<<"Transfer table is not increasing in k: i = "<<i<<", k = "<<ktransfer_table[i]<<" <= "<<ktransfer_table[i-1]<<std::endl;
              exit(48);
        }
  }
  for(int type=0; type < N_TYPES_TAB; type++){
      //Check everything is the same size
      assert( ktransfer_vector.size() == transfer_vector[type].size() );
      transfer_table[type] = new double[transfer_vector[type].size()];
      //Copy data from the temporary vectors to the statically sized C-arrays
      for (size_t i=0; i<transfer_vector[type].size(); i++) {
          //Make sure we do not take log(0)
          transfer_table[type][i] = transfer_vector[type][i];
      }
      //Set up the interpolation structures
      trans_interp[type] = gsl_interp_alloc(gsl_interp_cspline,transfer_vector[type].size());
      trans_interp_accel[type] = gsl_interp_accel_alloc();
      //Initialise the interpolator; same k values for each one, different T values.
      //Size asserted to be the same.
      gsl_interp_init(trans_interp[type],ktransfer_table, transfer_table[type],transfer_vector[type].size());
  }
  //Store the size of the arrays
  NTransferTable = ktransfer_vector.size();
}

PowerSpec_Tabulated::~PowerSpec_Tabulated()
{
        //Free memory for power table
        delete[] kmatter_table;
        delete[] pmatter_table;
        gsl_interp_free(pmat_interp);
        gsl_interp_accel_free(pmat_interp_accel);
        //Free memory for transfer table
        delete[] ktransfer_table;
        for (int i=0; i<N_TYPES_TAB;i++) {
            delete[] transfer_table[i];
            gsl_interp_free(trans_interp[i]);
            gsl_interp_accel_free(trans_interp_accel[i]);
        }
}

double PowerSpec_Tabulated::power(double k, int Type)
{
  double logk = log10(k*scale);
  double transfer;

  if(logk < ktransfer_table[0] || logk > ktransfer_table[NTransferTable - 1])
    return 0;

  if(logk < kmatter_table[0] || logk > kmatter_table[NPowerTable - 1])
    return 0;

  //If a type is requested that isn't defined, assume we want the total matter power.
  if (Type > N_TYPES_TAB-1 || Type < 0)
      transfer = 1;
  else
      transfer = gsl_interp_eval(trans_interp[Type], ktransfer_table, transfer_table[Type], logk, trans_interp_accel[Type]);
  /* Sometimes, due to numerical roundoff in CAMB, the massive neutrino transfer function will become
   * slightly negative. Set it to close to zero in this case.*/
  if(transfer < 0 && Type==2 && transfer > -1e-6){
      transfer = 1e-14;
  }

  double logP = gsl_interp_eval(pmat_interp, kmatter_table, pmatter_table, logk, pmat_interp_accel);
  double P = pow(10.0, logP) * pow(scale, 3) * transfer;
  assert(P >= 0);
  return P;
}

//Definitions for the normalized power spectrum
//These mostly compute sigma_8 via the GSL
NormalizedPowerSpec::NormalizedPowerSpec(PowerSpec * PSpec, double Sigma8, double PrimordialIndex, double Dplus, double UnitLength_in_cm): Dplus(Dplus), PSpec(PSpec), PrimordialIndex(PrimordialIndex), Norm(1.)
{
        /* 8 Mpc/h */
        R8 = 8 * (3.085678e24 / UnitLength_in_cm);
        //Uses R8
        Norm = Sigma8 * Sigma8 / TopHatSigma2();
        printf("Normalization adjusted by %g so  Sigma8 = %g\n", Norm, Sigma8);
}

double sigma2_int(double k, void * params)
{
  NormalizedPowerSpec * PSpec = (NormalizedPowerSpec *) params;
  double r = PSpec->R8;
  double kr, kr2, w, x;
  kr = r * k;
  kr2 = kr * kr;
  /*Power spectrum is at starting redshift, we want it at z=0*/
  double Dscale = pow(PSpec->Dplus,2);

  /*Series expansion; actually good until kr~1*/
  if(kr < 1e-3)
      w = 1./3. - kr2/30. +kr2*kr2/840.;
  else
      w = 3 * (sin(kr) / kr - cos(kr)) / kr2;
  x = 4 * M_PI * k * k * w * w * PSpec->power(k,100) * Dscale;
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

