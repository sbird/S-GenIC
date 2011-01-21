#include "allvars.h"


gadget_header header1, header;

int WhichSpectrum;


int SphereMode;

FILE *FdTmp, *FdTmpInput;

int Nmesh, Nsample;

char GlassFile[500];
char FileWithInputSpectrum[500];
char FileWithTransfer[500];

int GlassTileFac;

double Box;
int Seed;

int64_t NumPart;

struct part_data *P;

int Nglass;

double InitTime;
double Redshift;
double MassTable[6];

/*Parameters for spline knots*/
#ifdef SPLINE
int NumKnots;
char KnotValues[400];
char KnotPositions[400];
#endif

char OutputDir[1000], FileBase[1000];
int NumFilesWrittenInParallel;

fftwf_plan Inverse_plan;
float *Disp;
fftwf_complex *Cdata;


double UnitTime_in_s, UnitLength_in_cm, UnitMass_in_g, UnitVelocity_in_cm_per_s;
double InputSpectrum_UnitLength_in_cm;
double G, Hubble;
double RhoCrit;

double Omega, OmegaLambda, OmegaDM_2ndSpecies, Sigma8;
double OmegaBaryon, HubbleParam;
double ShapeGamma;
double PrimordialIndex;
double Dplus;			/* growth factor */

#ifdef DIFFERENT_TRANSFER_FUNC
int Type, MinType, MaxType;
#endif

int ReNormalizeInputSpectrum;

int WDM_On;
int WDM_Vtherm_On;
double WDM_PartMass_in_kev;


int NU_On;
int neutrinos_ks;
int NU_Vtherm_On;
double NU_PartMass_in_ev;

