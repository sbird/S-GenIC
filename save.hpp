#ifndef __SAVE_HPP
#define __SAVE_HPP

#include <valarray>
#include <gadgetheader.h>
#include <gadgetwriter.hpp>
#include "part_data.hpp"
#include "thermalvel.hpp"

#ifdef NO64BITID
        typedef int32_t id_type;
#else
        typedef int64_t id_type;
#endif //NO64BITID

double periodic_wrap(double x, double box);

int64_t write_particle_data(GadgetWriter::GWriteSnap & snap, int type, lpt_data * outdata, part_grid&  P, FermiDiracVel *therm_vels, int64_t FirstId);

gadget_header generate_header(std::valarray<int64_t> & npart, double Omega, double OmegaBaryon, double OmegaNuPart, double OmegaLambda, double HubbleParam, double Box, double InitTime, double UnitMass_in_g, double UnitLength_in_cm, double UnitVelocity_in_cm_per_s, bool combined_neutrinos);

#endif //__SAVE_HPP
