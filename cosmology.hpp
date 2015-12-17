#ifndef __COSMOLOGY_H
#define __COSMOLOGY_H

#include "physconst.h"

/** Class to store various cosmological parameters. Any function that relies on Omega or OmegaLambda should ultimately be in here */
class Cosmology
{
    public:
        //Note: Omega here should be the matter density in the background at z=0.
        //This always includes neutrinos and baryons.
        Cosmology(double HubbleParam, double Omega, double OmegaLambda, double MNu, bool InvertedHierarchy, bool Radiation=true, bool ApproximateGrowth=false): HubbleParam(HubbleParam), Omega(Omega), OmegaLambda(OmegaLambda), MNu(MNu), InvertedHierarchy(InvertedHierarchy), Radiation(Radiation), ApproximateGrowth(ApproximateGrowth)
        {}
        //Hubble H(z) / h0. Returns in units of 1/s, not internal units.
        double Hubble(double a);
        /* Return the total matter density in all neutrino species.*/
        double OmegaNu(double a);
        //Radiation density
        double OmegaR(double a);
        double GrowthFactor(double astart, double aend);
        double growth(double a);
        double F_Omega(double a);
        double F2_Omega(double a);
        double OmegaNuPrimed(double a);
        //The matter density
        double OmegaMatter(double a);
    private:
        /* Return the matter density in a single neutrino species.
        * Not externally callable*/
        double OmegaNu_single(double a,double mnu);
        double OmegaNuPrimed_single(double a,double mnu);
        /*The derivative wrt a of the total matter density in all neutrino species*/
        double HubbleParam;
        double Omega;
        double OmegaLambda;
        double MNu;
        bool InvertedHierarchy;
        //Should radiation be enabled
        bool Radiation;
        //Should we use an approximation to F_Omega
        bool ApproximateGrowth;
};

#endif
