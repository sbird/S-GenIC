#define BOOST_TEST_DYN_LINK

/** \file
 * Test suite using boost::test*/
#define BOOST_TEST_MODULE CREATENU
#include "createnu.hpp"
#include <boost/test/unit_test.hpp>
#include <boost/test/test_tools.hpp>

BOOST_AUTO_TEST_CASE(check_power_spec_nu)
{
    //Load the neutrino power spectrum from the powerspec file
    PowerSpec_NuTabulated pspec(std::string("testdata/powerspec_nu_005.txt"), 0.5);

    //Check that the tabulated power spectrum gives the right answer
    //First check ranges: these should both be out of range.
    //Should be the same k as in the file (but /10^3 for Mpc -> kpc)
    BOOST_CHECK_EQUAL(pspec.power(1.2e-5,6), 0.);
    BOOST_CHECK_EQUAL(pspec.power(0.00612,6), 0.);
    //Now check total power: k divided by 10^3,
    //Conversion for P(k) is 10^9/(2pi)^3
    BOOST_CHECK_CLOSE(pspec.power(3.44476e-05,6), 5.83134e+09,1e-6);
    BOOST_CHECK_CLOSE(pspec.power(4.15371e-05,6), 4.92876e+09,1e-6);
    //Check that it gives reasonable results when interpolating
    for (int k = 1; k < 100; k++) {
        double newk = 3.44476e-05+ k*(3.59104e-05-3.44476e-05)/100;
        BOOST_CHECK_LT(pspec.power(newk,6), pspec.power(3.44476e-05,6));
        BOOST_CHECK_GT(pspec.power(newk,6), pspec.power(3.59104e-05,6));
    }

}
