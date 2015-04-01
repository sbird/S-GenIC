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
class PowerSpec;

#ifdef PRINT_SPEC
void   print_spec(int type);
#endif
int    FatalError(int errnum);

void displacement_fields(const int type, const int64_t NumPart, part_data& P, PowerSpec * PSpec, const size_t Nmesh, bool RayleighScatter);
void   initialize_ffts(void);
unsigned int * initialize_rng(int Seed);
void   set_units(void);
double fnl(double x);

double displacement_read_out(float * Disp, const int order, const int64_t NumPart, part_data& P, const size_t Nmesh, const int axes);

double periodic_wrap(double x, double box);

int64_t write_particle_data(GadgetWriter::GWriteSnap & snap, int type, part_data&  P, int64_t NumPart, int64_t FirstId);

gadget_header generate_header(std::valarray<int64_t> & npart);

extern "C" {
#endif

void  read_parameterfile(char *fname);

#ifdef __cplusplus
}
#endif

#endif //__PROTO_H
