#ifndef __PROTO_H
#define __PROTO_H

#include <gsl/gsl_rng.h>

#ifdef __cplusplus
#include <gadgetreader.hpp>
#include <gadgetwriter.hpp>

void   print_spec(int type);
int    FatalError(int errnum);
void displacement_fields(int type, int64_t NumPart, struct part_data* P);
void   initialize_ffts(void);
void initialize_rng(gsl_rng * random_generator, unsigned int*& seedtable);
void   set_units(void);
double fnl(double x);

double periodic_wrap(double x);

int64_t read_glass(GadgetReader::GSnap& snap, int type, int GlassTileFac, struct part_data *& P);
int64_t write_particle_data(GadgetWriter::GWriteSnap & snap, int type, struct part_data * P, int64_t NumPart, int64_t FirstId);

gadget_header generate_header();


extern "C" {
#endif

void   read_power_table(void);
int    compare_logk(const void *a, const void *b);
double PowerSpec(double kmag, int Type);
double PowerSpec_Efstathiou(double k);
double PowerSpec_EH(double k);
double PowerSpec_Tabulated(double k);
double PowerSpec_TOTAL(double k);
double PowerSpec_DM_2ndSpecies(double k);
double PowerSpec_Tabulated2nd(double k);
double PowerSpec_Tabulated_b(double k);

void   initialize_powerspectrum(void);
double GrowthFactor(double astart, double aend);
double growth(double a);
double growth_int(double);
double qromb(double (*func)(double), double a, double b);
double sigma2_int(double k);
double TopHatSigma2(double R);
double F_Omega(double a);

void  read_parameterfile(char *fname);
double tk_eh(double k);


void add_WDM_thermal_speeds(float *vel);
void add_NU_thermal_speeds(float *vel);
double get_fermi_dirac_vel_nu(void);
void fermi_dirac_init_nu(void);
#ifdef __cplusplus
}
#endif

#endif //__PROTO_H
