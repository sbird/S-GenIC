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
class part_grid;
class lpt_data;

#ifdef PRINT_SPEC
void   print_spec(int type);
#endif
int    FatalError(int errnum);

double periodic_wrap(double x, double box);

int64_t write_particle_data(GadgetWriter::GWriteSnap & snap, int type, lpt_data& outdata, part_grid&  P, const double vel_prefac, const double vel_prefac2, int64_t FirstId, const bool twolpt);

gadget_header generate_header(std::valarray<int64_t> & npart, double Omega, double OmegaBaryon, double OmegaNuPart, double OmegaLambda, double HubbleParam, double Box, double InitTime, double UnitMass_in_g, double UnitLength_in_cm, bool combined_neutrinos);
extern "C" {
#endif

void   set_units(void);

#ifdef __cplusplus
}
#endif

#endif //__PROTO_H
