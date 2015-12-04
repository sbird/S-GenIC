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

extern char     FileWithInputSpectrum[500];
extern  char     FileWithTransfer[500];

extern double   Box;
extern int Seed;
extern int RayleighScatter;

extern double InitTime;
extern double Redshift;

extern char OutputDir[1000], FileBase[1000];

extern double UnitTime_in_s, UnitLength_in_cm, UnitMass_in_g, UnitVelocity_in_cm_per_s;
extern double InputSpectrum_UnitLength_in_cm;

#ifdef NO64BITID
        typedef int32_t id_type;
#else
        typedef int64_t id_type;
#endif //NO64BITID

#endif /* __ALLVARS_H*/

