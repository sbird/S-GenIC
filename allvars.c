#include "allvars.h"

int WhichSpectrum;

size_t Nmesh;

int ICFormat;

char FileWithInputSpectrum[500];
char FileWithTransfer[500];

double Box;
int Seed;
int RayleighScatter;

int twolpt;

int NumFiles;

double InitTime;
double Redshift;

char OutputDir[1000], FileBase[1000];

double UnitTime_in_s, UnitLength_in_cm, UnitMass_in_g, UnitVelocity_in_cm_per_s;
double InputSpectrum_UnitLength_in_cm;

void set_units(void)		/* ... set some units */
{
    InitTime = 1 / (1 + Redshift);
    UnitTime_in_s = UnitLength_in_cm / UnitVelocity_in_cm_per_s;
}

