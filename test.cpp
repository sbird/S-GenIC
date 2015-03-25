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
#define FLOATS_CLOSE_TO(x,y) \
      BOOST_CHECK_MESSAGE( fabs((x) - (y)) <= std::max<float>(fabs(x),fabs(y))/1e3,(x)<<" is not close to "<<(y))


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

class TestFermiDirac: public FermiDiracVel
{
    public:
       TestFermiDirac(int vel):FermiDiracVel(vel){}
       //For testing
       void reseed(int seed)
       {
           gsl_rng_set(g_rng, seed);
       }
       double get_fdv(int i)
       {
           if(i < LENGTH_FERMI_DIRAC_TABLE && i >= 0)
                return fermi_dirac_vel[i];
           return 0;
       }
       double get_fdvc(int i)
       {
           if(i < LENGTH_FERMI_DIRAC_TABLE && i >= 0)
                return fermi_dirac_cumprob[i];
           return 0;
       }
};

BOOST_AUTO_TEST_CASE(check_fermi_vel)
{
    //Check has units of velocity
    FLOATS_NEAR_TO(NU_V0(0, 1, 1e3), 100*NU_V0(0, 1, 1e5));
    //Check scales linearly with neutrino mass
    FLOATS_NEAR_TO(NU_V0(0, 0.1, 1e5), 10*NU_V0(0, 1, 1e5));
    //Check scales as z^3/2 (due to gadgets cosmological velocity unit)
    FLOATS_NEAR_TO(pow(0.5, 1.5)*NU_V0(1, 1, 1e5), NU_V0(0, 1, 1e5));
    //Check it is correct (roughly). This is
    //(4/11)^1/3* 2.7255* 1.00381* 8.61734e-5 * 2.99792e5/(M_nu/3)
    FLOATS_NEAR_TO(NU_V0(0,1, 1e5), 151.344);

    //Seed table with velocity of 100 km/s
    TestFermiDirac nuvels(100);
    nuvels.reseed(23);
    //Check that the probability table makes sense
    FLOATS_NEAR_TO(nuvels.get_fdvc(0), 0);
    FLOATS_NEAR_TO(nuvels.get_fdvc(LENGTH_FERMI_DIRAC_TABLE-2), 1);
    //Check that the probability table makes sense
    FLOATS_NEAR_TO(nuvels.get_fdv(0), 0);
    FLOATS_NEAR_TO(nuvels.get_fdv(LENGTH_FERMI_DIRAC_TABLE-1), MAX_FERMI_DIRAC);
    FLOATS_NEAR_TO(nuvels.get_fdv(LENGTH_FERMI_DIRAC_TABLE/2), MAX_FERMI_DIRAC/2.*(1.+1./LENGTH_FERMI_DIRAC_TABLE));
    //Test getting the distribution
    FLOATS_NEAR_TO(nuvels.get_fermi_dirac_vel(0), 0);
    FLOATS_NEAR_TO(nuvels.get_fermi_dirac_vel(1), MAX_FERMI_DIRAC);
    //Number computed by python guessing
    FLOATS_CLOSE_TO(nuvels.get_fermi_dirac_vel(0.5), 2.839);
//     for(int i = 0; i < LENGTH_FERMI_DIRAC_TABLE/2; i+=20)
//         printf("%g %g\n", nuvels.get_fdv(i), nuvels.get_fdvc(i));
    //Remember to reseed the rng...
    float vel[3]={0,0,0};
    nuvels.add_thermal_speeds(vel);
    //Check scaling is linear
    //Remember to reseed the rng...
    TestFermiDirac nuvels2(200);
    nuvels2.reseed(23);
    float vel2[3]={0,0,0};
    nuvels2.add_thermal_speeds(vel2);
    FLOATS_NEAR_TO(sqrt(vel2[0]*vel2[0]+vel2[1]*vel2[1]+vel2[2]*vel2[2]), 2*sqrt(vel[0]*vel[0]+vel[1]*vel[1]+vel[2]*vel[2]));
    for(int i=0; i<3; i++)
        FLOATS_NEAR_TO(vel2[i], 2*vel[i]);
    //Check some statistical properties (max, min, mean)
    double mean=0;
    double max = 0;
    double min = 1e10;
    int nsample;
    for (nsample=0; nsample < 10000; nsample++)
    {
        float vvel[3]={0,0,0};
        nuvels.add_thermal_speeds(vvel);
        double v2 =sqrt(vvel[0]*vvel[0]+vvel[1]*vvel[1]+vvel[2]*vvel[2]);
        max = std::max(v2, max);
        min = std::min(v2, min);
        mean+=v2;
    }
    mean/=nsample;
    //Mean should be roughly 3*zeta(4)/zeta(3)*7/8/(3/4)* m_vamp
    FLOATS_CLOSE_TO(mean, 3*pow(M_PI,4)/90./1.202057*(7./8)/(3/4.)*100);
    BOOST_CHECK_MESSAGE( min > 0,min<<" less than zero");
    BOOST_CHECK_MESSAGE( max < MAX_FERMI_DIRAC*100,max<<" greater than upper bound");

}
