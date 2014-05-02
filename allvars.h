#ifndef __ALLVARS_H
#define __ALLVARS_H

#include <fftw3.h>

#define  GRAVITY     6.672e-8
#define  HUBBLE      3.2407789e-18   /* in h/sec */


#include <gadgetheader.h>

typedef int32_t int4byte;
typedef uint32_t uint4byte;

extern int      WhichSpectrum;

extern int NumFiles;

//These should be 64 bit, to avoid Nmesh^3 overflowing.
//Make sure they are unsigned as well, as unsigned overflow
//is defined, but signed overflow is not.
extern size_t      Nmesh, Nsample;

extern int      SphereMode;

extern char     GlassFile[500]; 
extern char     FileWithInputSpectrum[500];
extern  char     FileWithTransfer[500];
extern size_t      GlassTileFac;

extern double   Box;
extern int Seed;
extern int RayleighScatter;

extern double InitTime;
extern double Redshift;

extern char OutputDir[1000], FileBase[1000];

extern fftwf_plan Inverse_plan;
extern float        *Disp;
extern fftwf_complex     *Cdata;

#ifdef TWOLPT
  extern fftwf_plan Forward_plan2;
  extern fftwf_plan Inverse_plan_grad[3];
  extern fftwf_plan Inverse_plan2;
  extern fftwf_complex *cdisp2; /* 2nd order displacements */
  extern float *disp2;
  extern float *twosrc;
  extern fftwf_complex *ctwosrc;
  extern fftwf_complex *(cdigrad[3]);
  extern float *(digrad[3]);
#endif


extern double UnitTime_in_s, UnitLength_in_cm, UnitMass_in_g, UnitVelocity_in_cm_per_s;
extern double InputSpectrum_UnitLength_in_cm;
extern double G, Hubble;

extern double Omega, OmegaLambda, OmegaDM_2ndSpecies, Sigma8;
extern double OmegaBaryon, HubbleParam;
extern double PrimordialIndex;
extern double ShapeGamma;

extern double Dplus; /* growth factor */

#ifdef NO64BITID
        typedef int32_t id_type;
#else
        typedef int64_t id_type;
#endif //NO64BITID

extern int    ReNormalizeInputSpectrum;

extern int    WDM_On;
extern int    WDM_Vtherm_On;
extern double WDM_PartMass_in_kev;

extern int    neutrinos_ks;
/* no_gas is set to true if there are no gas particles in the glass file.
 * This is so we can set OmegaBaryon to always be the value fed to CAMB*/
extern int    no_gas;
extern int    NU_On;
extern int    NU_Vtherm_On;
extern double NU_PartMass_in_ev;

#endif /* __ALLVARS_H*/

