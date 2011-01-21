#ifndef __ALLVARS_H
#define __ALLVARS_H

#include <fftw3.h>

#define  PI          3.14159265358979323846 
#define  GRAVITY     6.672e-8
#define  HUBBLE      3.2407789e-18   /* in h/sec */


#ifdef T3E
typedef short int int4byte;	/* Note: int has 8 Bytes on the T3E ! */
typedef unsigned short int uint4byte;	/* Note: int has 8 Bytes on the T3E ! */
#else
typedef int int4byte;
typedef unsigned int uint4byte;
#endif

#include <gadgetheader.h>

extern gadget_header header, header1;


extern int      Nglass;
extern int      WhichSpectrum;


extern FILE     *FdTmp, *FdTmpInput;

extern int      Nmesh, Nsample;

extern int      SphereMode;

extern long long IDStart;


extern char     GlassFile[500]; 
extern char     FileWithInputSpectrum[500];
extern  char     FileWithTransfer[500];
extern int      GlassTileFac; 

extern double   Box;
extern int Seed;

extern long long NumPart;

/*Parameters for spline knots*/
#ifdef SPLINE
extern int NumKnots;
extern char KnotValues[400];
extern char KnotPositions[400];
#endif


extern struct part_data 
{
  float Pos[3];
  float Vel[3];
#ifdef  MULTICOMPONENTGLASSFILE                      
  int   Type;
#endif
  long long ID;
} *P;


extern double InitTime;
extern double Redshift;
extern double MassTable[6];


extern char OutputDir[1000], FileBase[1000];
extern int  NumFilesWrittenInParallel;

extern int  IdStart;

extern fftwf_plan Inverse_plan;
extern float        *Disp;
extern fftwf_complex     *Cdata;


extern double UnitTime_in_s, UnitLength_in_cm, UnitMass_in_g, UnitVelocity_in_cm_per_s;
extern double InputSpectrum_UnitLength_in_cm;
extern double G, Hubble;
extern double RhoCrit;

extern double Omega, OmegaLambda, OmegaDM_2ndSpecies, Sigma8;
extern double OmegaBaryon, HubbleParam;
extern double PrimordialIndex;
extern double ShapeGamma;

extern double Dplus; /* growth factor */


#ifdef DIFFERENT_TRANSFER_FUNC
extern int Type, MinType, MaxType;
#endif

extern int    ReNormalizeInputSpectrum;

extern int    WDM_On;
extern int    WDM_Vtherm_On;
extern double WDM_PartMass_in_kev;

extern int    neutrinos_ks;
extern int    NU_On;
extern int    NU_Vtherm_On;
extern double NU_PartMass_in_ev;

#endif /* __ALLVARS_H*/

