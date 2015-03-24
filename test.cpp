/* Copyright (c) 2010, Simeon Bird <spb41@cam.ac.uk>
 *
 * Permission to use, copy, modify, and/or distribute this software for any
 * purpose with or without fee is hereby granted, provided that the above
 * copyright notice and this permission notice appear in all copies.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
 * WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
 * ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 * WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
 * ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
 * OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE. */
#define BOOST_TEST_DYN_LINK

/** \file
 * Test suite using boost::test*/
#define BOOST_TEST_MODULE NGENIC
#include "allvars.h"
#include "proto.h"
#include "part_data.hpp"
#include <math.h>
#include <boost/test/unit_test.hpp>
#include <boost/test/test_tools.hpp>

#define FLOATS_NEAR_TO(x,y) \
        BOOST_CHECK_MESSAGE( fabs((x) - (y)) <= std::max<float>(fabs(x),fabs(y))/1e5,(x)<<" is not close to "<<(y))

#ifdef PRINT_SPEC
BOOST_AUTO_TEST_CASE(check_print_spec)
{
        const int dims=5;
        float field[2*dims*dims*(dims/2+1)]={0};
        float pos[30];
        for(int i=0; i<30; i++)
                pos[i]=i/3.;
        fieldize(10, dims,field,10,10,pos,1);
        //Check off-diagonal elements are zero
        FLOATS_NEAR_TO(field[0],10.7638874);
        BOOST_CHECK_EQUAL(field[3],0);
        BOOST_CHECK_EQUAL(field[20],0);
        BOOST_CHECK_EQUAL(field[125],0);
        //Check on-diagonals
        FLOATS_NEAR_TO(field[124],2.08333);

}
#endif

// BOOST_AUTO_TEST_CASE(check_displacement_fields)
// {
//
// void displacement_fields(const int type, const int64_t NumPart, part_data& P, const int Nmesh, bool RayleighScatter);
// }

BOOST_AUTO_TEST_CASE(check_generate_header)
{
gadget_header generate_header(std::valarray<int64_t> & npart);


}

// BOOST_AUTO_TEST_CASE(check_displacement_read_out)
// {
// double displacement_read_out(float * Disp, const int order, const int64_t NumPart, part_data& P, const int Nmesh, const int axes);
//
//
// }

BOOST_AUTO_TEST_CASE(check_periodic_wrap)
{
    FLOATS_NEAR_TO(periodic_wrap(15000, 25000), 15000);
    FLOATS_NEAR_TO(periodic_wrap(65000, 25000), 15000);
    FLOATS_NEAR_TO(periodic_wrap(-10000, 25000), 15000);
    FLOATS_NEAR_TO(periodic_wrap(-60000, 25000), 15000);
}

// BOOST_AUTO_TEST_CASE(check_writing_data)
// {
// int64_t write_particle_data(GadgetWriter::GWriteSnap & snap, int type, part_data&  P, int64_t NumPart, int64_t FirstId);
//
// }

// BOOST_AUTO_TEST_CASE(check_power_spec_camb)
// {
// double PowerSpec_Tabulated(double k, int Type);
//
// }

BOOST_AUTO_TEST_CASE(check_fermi_vel)
{

//     double NU_V0(const double redshift, const double NU_PartMass_in_ev, const double UnitVelocity_in_cm_per_s)

//     FLOATS_NEAR_TO(NU_V0(0, 1, 1e5), );
    //Check has units of velocity
    FLOATS_NEAR_TO(NU_V0(0, 1, 1e3), 100*NU_V0(0, 1, 1e5));
    //Check scales linearly with neutrino mass
    FLOATS_NEAR_TO(NU_V0(0, 0.1, 1e5), 10*NU_V0(0, 1, 1e5));
    //Check scales as z^3/2 (due to gadgets cosmological velocity unit)
    FLOATS_NEAR_TO(pow(0.5, 1.5)*NU_V0(1, 1, 1e5), NU_V0(0, 1, 1e5));
    //Check it is correct (roughly). This is Q*i(4)/i(3)
    //where i(p) = zeta(p) Gamma(p) (1-2^(1-p))
    //and Q = (4/11)^1/3* 2.7255* 1.00381* 8.61734e-5 * 2.99792e5
    FLOATS_NEAR_TO(NU_V0(0,1, 1e5), 158.98);


//     FermiDiracVel nuvels(100);
//     //Remember to reseed the rng...
//     float vel[3]={0,0,0};
//     nuvels.add_thermal_speeds(vel);
//     FLOATS_NEAR_TO(vel, ???);
//     //Check scaling is linear
//     //Remember to reseed the rng...
//     FermiDiracVel nuvels2(200);
//     float vel2[3]={0,0,0};
//     nuvels2.add_thermal_speeds(vel2);
//     FLOATS_NEAR_TO(vel*2[1] == vel[1]);
//     //Check some statistical properties (max, min, mean)
}
