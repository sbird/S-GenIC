#include "allvars.h"

int WhichSpectrum;

int SphereMode;

size_t Nmesh, Nsample;

int ICFormat;

char GlassFile[500];
char FileWithInputSpectrum[500];
char FileWithTransfer[500];

size_t GlassTileFac;

double MNu;
int InvertedHierarchy;
double Box;
int Seed;
int RayleighScatter;

int NumFiles;

double InitTime;
double Redshift;

char OutputDir[1000], FileBase[1000];

double UnitTime_in_s, UnitLength_in_cm, UnitMass_in_g, UnitVelocity_in_cm_per_s;
double InputSpectrum_UnitLength_in_cm;

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
/*This triggers the use of neutrinos via an altered transfer function*/
int neutrinos_ks;
int no_gas;
int NU_Vtherm_On;
double NU_PartMass_in_ev;

