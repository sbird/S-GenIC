#ifndef __COSMOLOGY_H
#define __COSMOLOGY_H

#include "physconst.h"
#include <valarray>

/** Class to store various cosmological parameters. Any function that relies on Omega or OmegaLambda should ultimately be in here */
class Cosmology
{
    public:
        //Note: Omega here should be the matter density in the background at z=0.
        //This always includes neutrinos and baryons.
        Cosmology(double HubbleParam, double Omega, double OmegaLambda, double MNu, int Hierarchy, bool NoRadiation): HubbleParam(HubbleParam), Omega(Omega), OmegaLambda(OmegaLambda), MNu(MNu), Hierarchy(Hierarchy), NoRadiation(NoRadiation)
        {}
        //Hubble H(z) / h0. Returns in units of 1/s, not internal units.
        double Hubble(double a);
        /* Return the total matter density in all neutrino species.*/
        double OmegaNu(double a);
        //Radiation density
        double OmegaR(double a);
        double GrowthFactor(double astart, double aend);
        double F_Omega(double a);
        double F2_Omega(double a);
        double OmegaNuPrimed(double a);
        //The matter density
        double OmegaMatter(double a);
    private:
        double growth(double a, double *dDda);
        /* Return the matter density in a single neutrino species.
        * Not externally callable*/
        double OmegaNu_single(double a,double mnu);
        double OmegaNuPrimed_single(double a,double mnu);
        /*The derivative wrt a of the total matter density in all neutrino species*/
        double HubbleParam;
        double Omega;
        double OmegaLambda;
        double MNu;
        int Hierarchy;
        bool NoRadiation;
};

#endif
