#ifndef __ALLVARS_H
#define __ALLVARS_H

#include <stdint.h>
#include <stddef.h>

extern int      WhichSpectrum;

extern int NumFiles;

//These should be 64 bit, to avoid Nmesh^3 overflowing.
//Make sure they are unsigned as well, as unsigned overflow
//is defined, but signed overflow is not.
extern size_t      Nmesh;

extern int ICFormat;

extern int twolpt;

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

extern double UnitTime_in_s, UnitLength_in_cm, UnitMass_in_g, UnitVelocity_in_cm_per_s;
extern double InputSpectrum_UnitLength_in_cm;

extern double Omega, OmegaLambda, OmegaDM_2ndSpecies, Sigma8;
extern double OmegaBaryon, HubbleParam;
extern double PrimordialIndex;
extern double ShapeGamma;

#ifdef NO64BITID
        typedef int32_t id_type;
#else
        typedef int64_t id_type;
#endif //NO64BITID

extern int    ReNormalizeInputSpectrum;

extern int    WDM_On;
extern int    WDM_Vtherm_On;
extern double WDM_PartMass_in_kev;

extern int    NU_On;
extern int    NU_Vtherm_On;
//Total neutrino mass
extern double NU_PartMass_in_ev;
//1 if inverted hierarchy, 0 otherwise
extern int InvertedHierarchy;


#endif /* __ALLVARS_H*/

