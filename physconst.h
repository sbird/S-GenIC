#ifndef PHYS_CONSTANTS_H
#define PHYS_CONSTANTS_H
/*Define a few physical constants*/
#define  LIGHTCGS 2.99792e10  //Speed of light in cm/s
#define BOLEVK 8.61734e-5        /*The Boltzmann constant in units of eV/K*/
/* Stefan-Boltzmann constant in cgs units*/
#define STEFAN_BOLTZMANN 5.670373e-5
#define  T_CMB0      2.7255	/* present-day CMB temperature, from Fixsen 2009 */
#define HBAR    6.582119e-16  /*hbar in units of eV s*/
/* Note there is a slight correction from 4/11
 * due to the neutrinos being slightly coupled at e+- annihilation.
 * See Mangano et al 2005 (hep-ph/0506164)
 *The correction is (3.046/3)^(1/4), for N_eff = 3.046 */
#define TNU     (T_CMB0*pow(4/11.,1/3.)*1.00381)              /* Neutrino + antineutrino background temperature in Kelvin */

#define  GRAVITY     6.672e-8 // cm^3 g^-1 s^-2
#define  HUBBLE      3.2407789e-18   /* in h/sec */

#endif
