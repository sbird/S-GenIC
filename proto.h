#ifndef __PROTO_H
#define __PROTO_H

#include <gsl/gsl_rng.h>

#ifdef __cplusplus

#include <valarray>
#include <gadgetheader.h>

//Forward type declarations to save header include
namespace GadgetWriter{
    class GWriteSnap;
}
class part_data;

#ifdef PRINT_SPEC
void   print_spec(int type);
#endif
int    FatalError(int errnum);
void displacement_fields(const int type, const int64_t NumPart, part_data& P, const size_t Nmesh, bool RayleighScatter);
void   initialize_ffts(void);
unsigned int * initialize_rng(int Seed);
void   set_units(void);
double fnl(double x);

double displacement_read_out(float * Disp, const int order, const int64_t NumPart, part_data& P, const size_t Nmesh, const int axes);

double periodic_wrap(double x, double box);

int64_t write_particle_data(GadgetWriter::GWriteSnap & snap, int type, part_data&  P, int64_t NumPart, int64_t FirstId);

gadget_header generate_header(std::valarray<int64_t> & npart);

//Class defined in thermalvel.cpp for adding FermiDirac velocities
#define LENGTH_FERMI_DIRAC_TABLE 2000
#define MAX_FERMI_DIRAC          20.0

class FermiDiracVel
{
    public:
        //Single parameter is the amplitude of the random velocities. All the physics is in here.
        FermiDiracVel(double v_amp);
        void add_thermal_speeds(float *vel);
        double get_fermi_dirac_vel(double p);
        ~FermiDiracVel()
        {
            gsl_rng_free(g_rng);
        }
    protected:
        double fermi_dirac_vel[LENGTH_FERMI_DIRAC_TABLE];
        double fermi_dirac_cumprob[LENGTH_FERMI_DIRAC_TABLE];
        const double m_vamp;
        gsl_rng * g_rng;
};

//Amplitude of the random velocity for WDM
double WDM_V0(const double redshift, const double WDM_PartMass_in_kev, const double Omega_CDM, const double HubbleParam, const double UnitVelocity_in_cm_per_s);

#ifdef NEUTRINOS
//Amplitude of the random velocity for neutrinos
double NU_V0(const double redshift, const double NU_PartMass_in_ev, const double UnitVelocity_in_cm_per_s);
#endif //NEUTRINOS

double Hubble_A(double a, double Omega, double OmegaLambda);

extern "C" {
#endif

void   read_power_table(void);
int    compare_logk(const void *a, const void *b);
double PowerSpec(double kmag, int Type);
double PowerSpec_Efstathiou(double k);
double PowerSpec_EH(double k);
double PowerSpec_Tabulated(double k,int Type);
double PowerSpec_DM_2ndSpecies(double k);

void   initialize_powerspectrum(void);
double GrowthFactor(double astart, double aend);
double growth(double a);
double growth_int(double, void *);
double sigma2_int(double k, void * params);
double TopHatSigma2(double R);
double F_Omega(double a);
double F2_Omega(double a);

void  read_parameterfile(char *fname);
double tk_eh(double k);

#ifdef __cplusplus
}
#endif

#endif //__PROTO_H
