
#include <gsl/gsl_rng.h>

void   print_spec(void);
int    FatalError(int errnum);
void   displacement_fields(int type);
void   initialize_ffts(void);
void   set_units(void);
void   assemble_particles(void);
void   free_ffts(void);
double fnl(double x);

void   assemble_grid(void);
double periodic_wrap(double x);

void  read_glass(char *fname,int type);

void  write_particle_data(int type);
gadget_header generate_header();


#ifdef __cplusplus
extern "C" {
#endif

void   read_power_table(void);
int    compare_logk(const void *a, const void *b);
double PowerSpec(double kmag);
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
