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
#include "part_data.hpp"
#include "cosmology.hpp"
#include "thermalvel.hpp"
#include "power.hpp"
#include "save.hpp"
#include <math.h>
#include <boost/test/unit_test.hpp>
#include <boost/test/test_tools.hpp>
#include <gsl/gsl_sf_hyperg.h>

// BOOST_AUTO_TEST_CASE(check_displacement_fields)
// {
//
// void displacement_fields(const int type, const int64_t NumPart, part_data& P, const int Nmesh, bool RayleighScatter);
// }

BOOST_AUTO_TEST_CASE(check_generate_header)
{
    //Check the header is correctly generated
    std::valarray<int64_t> npart(0L,6);
    npart[1] = pow(512,3);
    gadget_header header = generate_header(npart, 0.2793, 0.0463, 0.00986007458598642, 0.7207, 0.7, 512000, 0.01, 1.989e43, 3.085678e21, 1e5, false);
    BOOST_CHECK_EQUAL(header.BoxSize, 512000);
    BOOST_CHECK_EQUAL(header.npartTotal[1], pow(512,3));
    //This is the most non-trivial quantity
    BOOST_CHECK_CLOSE_FRACTION(header.mass[1],7.4759756,1e-6);
    //Does it have the right dimensions? If we change mass units, is it proportional?
    gadget_header header2 = generate_header(npart, 0.2793, 0.0463, 0.00986007458598642, 0.7207, 0.7, 512000, 0.01, 1.989e33, 3.085678e21,1e5, false);
    BOOST_CHECK_CLOSE(header.mass[1], header2.mass[1]/1e10,1e-6);
    //It should not depend on length units.
    gadget_header header3 = generate_header(npart, 0.2793, 0.0463, 0.00986007458598642, 0.7207, 0.7, 512, 0.01, 1.989e43, 3.085678e24, 1e5, false);
    BOOST_CHECK_CLOSE(header.mass[1], header3.mass[1],1e-6);
}

// BOOST_AUTO_TEST_CASE(check_displacement_read_out)
// {
// double displacement_read_out(float * Disp, const int order, const int64_t NumPart, part_data& P, const int Nmesh, const int axes);
//
//
// }

BOOST_AUTO_TEST_CASE(check_periodic_wrap)
{
    BOOST_CHECK_EQUAL(periodic_wrap(15000, 25000), 15000);
    BOOST_CHECK_EQUAL(periodic_wrap(65000, 25000), 15000);
    BOOST_CHECK_EQUAL(periodic_wrap(-10000, 25000), 15000);
    BOOST_CHECK_EQUAL(periodic_wrap(-60000, 25000), 15000);
}

// BOOST_AUTO_TEST_CASE(check_writing_data)
// {
// int64_t write_particle_data(GadgetWriter::GWriteSnap & snap, int type, part_data&  P, int64_t NumPart, int64_t FirstId);
//
// }

BOOST_AUTO_TEST_CASE(check_power_spec_camb)
{
    std::string trans = "testdata/cambv1_transfer_99.dat";
    std::string matpow = "testdata/cambv1_matterpow_99.dat";
    PowerSpec_Tabulated pspec(trans.c_str(), matpow.c_str(), 0.222+0.0449, 0.7331, 0.0449,0.,3.085678e24, 3.085678e21, false, false);
    BOOST_CHECK_EQUAL(pspec.size(),335);
    std::string trans2 = "testdata/cambv2_transfer_99.dat";
    std::string matpow2 = "testdata/cambv2_matterpow_99.dat";
    PowerSpec_Tabulated pspec2(trans2.c_str(), matpow2.c_str(), 0.233+0.0463, 0.7331, 0.0463,0.,3.085678e24, 3.085678e21, false, false);
    BOOST_CHECK_EQUAL(pspec2.size(),335);
    //Check that the tabulated power spectrum gives the right answer
    //First check ranges: these should both be out of range.
    //Should be the same k as in the file (but /10^3 for Mpc -> kpc)
    BOOST_CHECK_EQUAL(pspec.power(9.8e-9,6), 0.);
    BOOST_CHECK_EQUAL(pspec.power(2.1e-1,6), 0.);
    //Now check total power: k divided by 10^3,
    //Conversion for P(k) is 10^9/(2pi)^3
    BOOST_CHECK_CLOSE(pspec.power(0.11239E-01/1e3,6), 0.53091E+01*1e9/(pow(2*M_PI,3)),1e-6);
    BOOST_CHECK_CLOSE(pspec.power(0.10022E+01/1e3,6), 0.13121E-01*1e9/(pow(2*M_PI,3)),1e-6);
    //Check that it gives reasonable results when interpolating
    for (int k = 1; k < 100; k++) {
        double newk = 0.10022E+01/1e3+ k*(0.10362E+01-0.10022E+01)/1e3/100;
        BOOST_CHECK_LT(pspec.power(newk,6), pspec.power(0.10022E+01/1e3,6));
        BOOST_CHECK_GT(pspec.power(newk,6), pspec.power(0.10362E+01/1e3,6));
    }
    //Now check transfer functions: ratio of total to species should be ratio of T_s to T_tot squared: large scales where T~ 1
    BOOST_CHECK_CLOSE(pspec.power(0.208152E-02/1e3,0)/pspec.power(0.208152E-02/1e3,6), pow(0.255362E+06/0.255697E+06,2),1e-6);
    //CDM
    BOOST_CHECK_CLOSE(pspec.power(0.208152E-02/1e3,1)/pspec.power(0.208152E-02/1e3,6), pow(0.255765E+06/0.255697E+06,2),1e-6);
    //massless NU
    BOOST_CHECK_CLOSE(pspec.power(0.208152E-02/1e3,3)/pspec.power(0.208152E-02/1e3,6), pow(0.274272E+06/0.255697E+06,2),1e-6);
    //Massive nu (which are zero here)
    BOOST_CHECK_SMALL(pspec.power(0.208152E-02/1e3,2),1e-6);
    //Small scales where there are differences
    //T_tot=0.255697E+06
    //Baryons
    BOOST_CHECK_CLOSE(pspec2.power(0.111030E+00/1e3,0)/pspec2.power(0.111030E+00/1e3,6), pow(0.200504E+05/0.277947E+05,2),1e-6);
    //CDM
    BOOST_CHECK_CLOSE(pspec2.power(0.111030E+00/1e3,1)/pspec2.power(0.111030E+00/1e3,6), pow(0.293336E+05/0.277947E+05,2),1e-6);
    //Check that normalisation works.
    PowerSpec_Tabulated *pspec_ptr = new PowerSpec_Tabulated(trans.c_str(), matpow.c_str(), 0.222+0.0449, 0.7331, 0.0449,0.,3.085678e24, 3.085678e21, false, false);
    NormalizedPowerSpec normspec(pspec_ptr, 0.918785/80., 1.0, 80, 3.085678e21);
    BOOST_CHECK_CLOSE(normspec.power(0.208152E-02/1e3,1)*80*80, pspec.power(0.208152E-02/1e3,1),1e-3);
}

class TestFermiDirac: public FermiDiracVel
{
    public:
       TestFermiDirac(int vel):FermiDiracVel(vel){}
       //For testing
       void reseed(int seed)
       {
           gsl_rng_set(g_rng, seed);
       }
//        double get_fdv(int i)
//        {
//            if(i < LENGTH_FERMI_DIRAC_TABLE && i >= 0)
//                 return fermi_dirac_vel[i];
//            return 0;
//        }
//        double get_fdvc(int i)
//        {
//            if(i < LENGTH_FERMI_DIRAC_TABLE && i >= 0)
//                 return fermi_dirac_cumprob[i];
//            return 0;
//        }
};

BOOST_AUTO_TEST_CASE(check_fermi_vel)
{
    //Make a cosmology
    Cosmology cosmo(0.7,0.3, 0.7, 1, 0, 0);
    //Check has units of velocity
    BOOST_CHECK_CLOSE(cosmo.NU_V0(0, 1e3), 100*cosmo.NU_V0(0, 1e5),1e-6);
    //Make a cosmology
    Cosmology cosmolow(0.7,0.3, 0.7, 0.1, 0, 0);
    //Check scales linearly with neutrino mass
    BOOST_CHECK_CLOSE(cosmolow.NU_V0(0, 1e5), 10*cosmo.NU_V0(0, 1e5),1e-6);
    //Check scales as z^3/2 (due to gadgets cosmological velocity unit)
//     BOOST_CHECK_CLOSE(pow(0.5, 1.5)*NU_V0(1, 1, 1e5), NU_V0(0, 1, 1e5),1e-6);
    //Check it is correct (roughly). This is
    //(4/11)^1/3* 2.7255* 1.00381* 8.61734e-5 * 2.99792e5/(M_nu/3)
    BOOST_CHECK_CLOSE(cosmo.NU_V0(0, 1e5), 151.265,1e-4);

    //Seed table with velocity of 100 km/s
    TestFermiDirac nuvels(100);
    nuvels.reseed(23);
    //Check that the probability table makes sense
//     BOOST_CHECK_CLOSE(nuvels.get_fdvc(0), 0,1e-6);
//     BOOST_CHECK_CLOSE(nuvels.get_fdvc(LENGTH_FERMI_DIRAC_TABLE-2), 1,1e-6);
    //Check that the probability table makes sense
//     BOOST_CHECK_CLOSE(nuvels.get_fdv(0), 0,1e-6);
//     BOOST_CHECK_CLOSE(nuvels.get_fdv(LENGTH_FERMI_DIRAC_TABLE-1), MAX_FERMI_DIRAC,1e-6);
//     BOOST_CHECK_CLOSE(nuvels.get_fdv(LENGTH_FERMI_DIRAC_TABLE/2), MAX_FERMI_DIRAC/2.*(1.+1./LENGTH_FERMI_DIRAC_TABLE),3e-5);
    //Test getting the distribution
    BOOST_CHECK_CLOSE(nuvels.get_fermi_dirac_vel(0), 0,1e-6);
    BOOST_CHECK_CLOSE(nuvels.get_fermi_dirac_vel(1), 100*MAX_FERMI_DIRAC,1e-3);
    //Number verified by mathematica
    BOOST_CHECK_CLOSE(nuvels.get_fermi_dirac_vel(0.5), 100*2.839075,1e-3);
//     for(int i = 0; i < LENGTH_FERMI_DIRAC_TABLE/2; i+=20)
//         printf("%g %g\n", nuvels.get_fdv(i), nuvels.get_fdvc(i));
    //Remember to reseed the rng...
    std::valarray<float> vel = nuvels.get_thermal_speeds();
    //Check scaling is linear
    //Remember to reseed the rng...
    TestFermiDirac nuvels2(200);
    nuvels2.reseed(23);
    std::valarray<float> vel2 = nuvels2.get_thermal_speeds();
    BOOST_CHECK_CLOSE(sqrt(vel2[0]*vel2[0]+vel2[1]*vel2[1]+vel2[2]*vel2[2]), 2*sqrt(vel[0]*vel[0]+vel[1]*vel[1]+vel[2]*vel[2]),1e-6);
    for(int i=0; i<3; i++)
        BOOST_CHECK_CLOSE(vel2[i], 2*vel[i],1e-6);
    //Check some statistical properties (max, min, mean)
    double mean=0;
    double max = 0;
    double min = 1e10;
    int nsample;
    for (nsample=0; nsample < 10000; nsample++)
    {
        std::valarray<float> vvel = nuvels.get_thermal_speeds();
        double v2 =sqrt(vvel[0]*vvel[0]+vvel[1]*vvel[1]+vvel[2]*vvel[2]);
        max = std::max(v2, max);
        min = std::min(v2, min);
        mean+=v2;
    }
    mean/=nsample;
    //Mean should be roughly 3*zeta(4)/zeta(3)*7/8/(3/4)* m_vamp
    BOOST_CHECK_CLOSE(mean, 3*pow(M_PI,4)/90./1.202057*(7./8)/(3/4.)*100,1e-1);
    BOOST_CHECK_MESSAGE( min > 0,min<<" less than zero");
    BOOST_CHECK_MESSAGE( max < MAX_FERMI_DIRAC*100,max<<" greater than upper bound");

}

BOOST_AUTO_TEST_CASE(check_cosmology)
{
    //Check that we get the right scalings for total matter domination.
    //Cosmology(double HubbleParam, double Omega, double OmegaLambda, double MNu, int Hierarchy, bool NoRadiation)
    Cosmology cosmo(0.7, 1., 0., 0., 0,false);
    BOOST_CHECK_CLOSE(cosmo.Hubble(1), HUBBLE,1e-2);
    BOOST_CHECK_CLOSE(cosmo.Hubble(0.1), cosmo.Hubble(1)/pow(0.1,3/2.),1e-1);
    BOOST_CHECK_CLOSE(cosmo.GrowthFactor(0.5,1), 2,2e-2);
    //Check that massless neutrinos work
    BOOST_CHECK_CLOSE(cosmo.OmegaNu(1), (cosmo.OmegaR(1)-cosmo.OmegaNu(1))*7./8.*pow(pow(4/11.,1/3.)*1.00328,4)*3,2e-3);
    BOOST_CHECK_CLOSE(cosmo.OmegaNu(0.01), cosmo.OmegaNu(1)/pow(0.01,4),1e-6);
    //Check that the velocity correction d ln D1/d lna is constant
    BOOST_CHECK_CLOSE(1.0, cosmo.F_Omega(1.5),1e-1);
    BOOST_CHECK_CLOSE(1.0, cosmo.F_Omega(2),1e-2);
    //Check radiation against exact solution from gr-qc/0504089
    double omegar = cosmo.OmegaR(1);
    auto radgrow = [omegar](double aa){ return omegar + 1.5 * 1. * aa; };
    BOOST_CHECK_CLOSE(cosmo.GrowthFactor(0.05,1), radgrow(1.)/radgrow(0.05),1e-3);
    BOOST_CHECK_CLOSE(cosmo.GrowthFactor(0.001,0.01), radgrow(0.01)/radgrow(0.001),1e-2);

    //Check against exact solutions from gr-qc/0504089: No radiation!
    //Note that the GSL hyperg needs the last argument to be < 1
    double omegal = 0.5;
    double omegam = 0.5;
    Cosmology cosmo3(0.7, omegam, omegal, 0., 0,true);
    //Omega_L + Omega_M = 1 => D+ ~ Gauss hypergeometric function
    auto growth = [omegam, omegal](double aa) { return aa * gsl_sf_hyperg_2F1(1./3, 1, 11./6, -omegal/omegam*pow(aa,3));};
    //Check growth factor during matter domination
    BOOST_CHECK_CLOSE(cosmo3.GrowthFactor(0.5,1), growth(1.)/growth(0.5),1e-3);
    BOOST_CHECK_CLOSE(cosmo3.GrowthFactor(0.15,0.3), growth(0.3)/growth(0.15),1e-3);
    BOOST_CHECK_CLOSE(cosmo3.GrowthFactor(0.01,1.), growth(1)/growth(0.01),1e-3);
    BOOST_CHECK_CLOSE(0.01*log(cosmo3.GrowthFactor(0.01-1e-5,0.01+1e-5))/2e-5, cosmo3.F_Omega(0.01),1e-3);
    //Massive neutrinos
    Cosmology nuc(0.7, 0.3, 0.7, 1.0, 0,false);
    BOOST_CHECK_CLOSE(0.01*log(nuc.GrowthFactor(0.01-1e-5,0.01+1e-5))/2e-5, nuc.F_Omega(0.01),1e-3);
    BOOST_CHECK_CLOSE(nuc.OmegaNu(0.5), nuc.OmegaNu(1.)/0.125,1e-3);
    BOOST_CHECK_CLOSE(nuc.OmegaNu(1.), 1.0/93.14/0.7/0.7,1e-2);
    BOOST_CHECK_CLOSE(nuc.OmegaNu(0.00001)*pow(0.00001,4), nuc.OmegaNu(0.00002)*pow(0.00002,4),1e-2);
}

//Neutrino mass spectrum allowed by oscillation experiments.
std::valarray<double> NuPartMasses(double mnu, int Hierarchy);
#define M21 7.53e-5 //Particle data group 2016: +- 0.18e-5 eV2
#define M32n 2.44e-3 //Particle data group: +- 0.06e-3 eV2
#define M32i 2.51e-3 //Particle data group: +- 0.06e-3 eV2

BOOST_AUTO_TEST_CASE(check_numass)
{
    auto numass = NuPartMasses(0.3,0);
    for(auto nn : numass)
        BOOST_CHECK_CLOSE(nn,0.1,1e-3);
    numass = NuPartMasses(0.3,1);
    /*Check the original inequalities are satisfied*/
    BOOST_CHECK_CLOSE(numass[0]+numass[1]+numass[2], 0.3,1e-3);
    BOOST_CHECK_CLOSE(numass[0]*numass[0] - numass[1]*numass[1], M32n,1e-3);
    BOOST_CHECK_CLOSE(numass[1]*numass[1] - numass[2]*numass[2], M21,1e-3);
    numass = NuPartMasses(0.08,1);
    BOOST_CHECK_CLOSE(numass[0]+numass[1]+numass[2], 0.08,1e-3);
    BOOST_CHECK_CLOSE(numass[0]*numass[0] - numass[1]*numass[1], M32n,1e-3);
    BOOST_CHECK_CLOSE(numass[1]*numass[1] - numass[2]*numass[2], M21,1e-3);
    numass = NuPartMasses(0.11,-1);
    BOOST_CHECK_CLOSE(numass[0]+numass[1]+numass[2], 0.11,1e-3);
    BOOST_CHECK_CLOSE(numass[0]*numass[0] - numass[1]*numass[1], -M32i,1e-3);
    BOOST_CHECK_CLOSE(numass[1]*numass[1] - numass[2]*numass[2], M21,1e-3);
}

