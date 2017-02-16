#ifndef POWERSPEC_H
#define POWERSPEC_H
//For pow
#include <cmath>
#include <string>
#include <vector>
#include <gsl/gsl_interp.h>


//Abstract base class to allow calling of different types of power spectrum
class PowerSpec
{
    public:
        virtual ~PowerSpec() {}
    public:
        /** Returns the power spectrum.
         * Units are L^3. L is N-GenIC internal units (usually kpc/h)
         * Fourier convention differs from CAMB by 1/(2 pi^3).
         * k should be 1/L in  and is not log.
         * Type specifies the type of particle and follows Gadget conventions. In practice it is ignored for everything apart
         * from PowerSpec_Tabulated, in which case:
         * 0 -> baryons. 1-> CDM, 2-> Neutrinos 3-> massless neutrinos 4 (or higher) -> total power
         */
       virtual double power(double k, int Type)=0;
};

//Number of types of particle defined in the tabulated power spectrum.
#define N_TYPES_TAB 4

//Derived class for a tabulated power spectrum
class PowerSpec_Tabulated: public PowerSpec
{
    public:
    PowerSpec_Tabulated(const std::string& FileWithTransfer, const std::string& FileWithInputSpectrum, double Omega, double OmegaLambda, double OmegaBaryon, double OmegaNu,
                        double InputSpectrum_UnitLength_in_cm, double UnitLength_in_cm, bool no_gas, bool combined_neutrinos);
    virtual double power(double k, int Type);
    virtual ~PowerSpec_Tabulated();
    size_t size()
    {
        return NPowerTable;
    }
    private:
        double scale;
        size_t NPowerTable, NTransferTable;
        //Data tables for the transfer function
        double * ktransfer_table;
        double * transfer_table[N_TYPES_TAB];
        gsl_interp * trans_interp[N_TYPES_TAB];
        gsl_interp_accel * trans_interp_accel[N_TYPES_TAB];
        //Data table for the power spectrum
        double * kmatter_table;
        double * pmatter_table;
        gsl_interp * pmat_interp;
        gsl_interp_accel * pmat_interp_accel;
};

//Derived class for the Efstathiou power spectrum
class PowerSpec_Efstathiou: public PowerSpec
{
    public:
        PowerSpec_Efstathiou(double ShapeGamma, double UnitLength_in_cm)
        {
            AA = 6.4 / ShapeGamma *  (3.085678e24 / UnitLength_in_cm);
            BB = 3.0 / ShapeGamma *  (3.085678e24 / UnitLength_in_cm);
            CC = 1.7 / ShapeGamma *  (3.085678e24 / UnitLength_in_cm);
            nu = 1.13;
        }
        virtual double power(double k, int Type)
        {
            return k / pow(1 + pow(AA * k + pow(BB * k, 1.5) + CC * CC * k * k, nu), 2 / nu);
        }
        virtual ~PowerSpec_Efstathiou() {}
    private:
        double AA, BB, CC;
        double nu;
};

//Derived class for the Eisenstein & Hu power spectrum
class PowerSpec_EH: public PowerSpec
{
    public:
        PowerSpec_EH(double HubbleParam, double Omega, double OmegaBaryon, double UnitLength_in_cm): omegam(Omega), hubble(HubbleParam), UnitLength_in_cm(UnitLength_in_cm)
        {
            /* Set up Omega Baryon*/
            ombh2 = OmegaBaryon * HubbleParam * HubbleParam;

            if(OmegaBaryon == 0)
                ombh2 = 0.04 * HubbleParam * HubbleParam;
        }
        virtual ~PowerSpec_EH() {}
        virtual double power(double k, int Type)
        {
              return k * pow(tk_eh(k), 2);
        }
    private:
        double tk_eh(double k)		/* from Martin White */
        {
            double q, theta, ommh2, a, s, gamma, L0, C0;
            double tmp;

            k *= (3.085678e24 / UnitLength_in_cm);	/* convert to h/Mpc */

            theta = 2.728 / 2.7;
            ommh2 = omegam * hubble * hubble;
            s = 44.5 * log(9.83 / ommh2) / sqrt(1. + 10. * exp(0.75 * log(ombh2))) * hubble;
            a = 1. - 0.328 * log(431. * ommh2) * ombh2 / ommh2
                + 0.380 * log(22.3 * ommh2) * (ombh2 / ommh2) * (ombh2 / ommh2);
            gamma = a + (1. - a) / (1. + exp(4 * log(0.43 * k * s)));
            gamma *= omegam * hubble;
            q = k * theta * theta / gamma;
            L0 = log(2. * exp(1.) + 1.8 * q);
            C0 = 14.2 + 731. / (1. + 62.5 * q);
            tmp = L0 / (L0 + C0 * q * q);
            return (tmp);
        }
        double omegam, ombh2, hubble,UnitLength_in_cm;
};

//A decorator on the power spec class
//That normalizes the input power spectrum
class NormalizedPowerSpec: public PowerSpec
{
    public:
        NormalizedPowerSpec(PowerSpec * PSpec, double Sigma8, double PrimordialIndex, double Dplus, double UnitLength_in_cm);
        virtual double power(double k, int Type)
        {
            double power = Norm * PSpec->power(k, Type);
            power *= pow(k, PrimordialIndex - 1.0);
            return power/ (Dplus*Dplus);
        }
        ~NormalizedPowerSpec()
        {
            delete PSpec;
        }
        //public as is used in the integrator function, which is not in this class
        double R8, Dplus;
    private:
        PowerSpec * PSpec;
        //Note: Uses the R8 parameter below and calls power(k, Type)
        double TopHatSigma2();
        double PrimordialIndex, Norm;
};

//A decorator for adding the effect of warm dark matter
class WDMPowerSpec: public PowerSpec
{
    public:
        WDMPowerSpec(PowerSpec * PSpec, double WDM_PartMass_in_kev, double Omega, double OmegaBaryon, double HubbleParam, double UnitLength_in_cm): PSpec(PSpec), UnitLength_in_cm(UnitLength_in_cm)
        {
            /* Eqn. (A9) in Bode, Ostriker & Turok (2001), assuming gX=1.5  */
            alpha =  0.048 * pow((Omega - OmegaBaryon) / 0.4, 0.15) * pow(HubbleParam / 0.65,1.3) * pow(1.0 / WDM_PartMass_in_kev, 1.15);
        }
        ~WDMPowerSpec()
        {
            delete PSpec;
        }
        virtual double power(double k, int Type)
        {
            double WDMTf = pow(1 + pow(alpha * k * (3.085678e24 / UnitLength_in_cm), 2 * 1.2), -5.0 / 1.2);
            double power = PSpec->power(k, Type);
            return power*WDMTf*WDMTf;
        }
    private:
        PowerSpec * PSpec;
        double alpha,UnitLength_in_cm;
};

#endif
