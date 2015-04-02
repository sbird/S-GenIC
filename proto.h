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

void   set_units(void);

double periodic_wrap(double x, double box);

int64_t write_particle_data(GadgetWriter::GWriteSnap & snap, int type, part_data&  P, int64_t NumPart, int64_t FirstId, bool twolpt);

gadget_header generate_header(std::valarray<int64_t> & npart);

extern "C" {
#endif

void  read_parameterfile(char *fname);

#ifdef __cplusplus
}
#endif

#endif //__PROTO_H
