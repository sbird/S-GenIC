#include "allvars.h"


int WhichSpectrum;

int SphereMode;

int Nmesh, Nsample;

char GlassFile[500];
char FileWithInputSpectrum[500];
char FileWithTransfer[500];

int GlassTileFac;

double Box;
int Seed;

int NumFiles;

double InitTime;
double Redshift;

char OutputDir[1000], FileBase[1000];
int NumFilesWrittenInParallel;

fftwf_plan Inverse_plan;
float *Disp;
fftwf_complex *Cdata;

#ifdef TWOLPT
  fftwf_plan Forward_plan2;
  fftwf_plan Inverse_plan_grad[3];
  float *twosrc;
  fftwf_complex *ctwosrc;
  fftwf_complex *(cdigrad[3]);
  float *(digrad[3]);
#endif


double UnitTime_in_s, UnitLength_in_cm, UnitMass_in_g, UnitVelocity_in_cm_per_s;
double InputSpectrum_UnitLength_in_cm;
double G, Hubble;

double Omega, OmegaLambda, OmegaDM_2ndSpecies, Sigma8;
double OmegaBaryon, HubbleParam;
double ShapeGamma;
double PrimordialIndex;
double Dplus;			/* growth factor */

int ReNormalizeInputSpectrum;

int WDM_On;
int WDM_Vtherm_On;
double WDM_PartMass_in_kev;


int NU_On;
int neutrinos_ks;
int no_gas;
int NU_Vtherm_On;
double NU_PartMass_in_ev;

