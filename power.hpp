#ifndef POWERSPEC_H
#define POWERSPEC_H
//For pow
#include <cmath>

//Abstract base class to allow calling of different types of power spectrum
class PowerSpec
{
    public:
        virtual double power(double k, int Type)=0;
};

/*Structure for matter power table*/
struct pow_table
{
  double logk, logD,logDb;
  double logD2nd;
  /*For a second neutrino species: type 3*/
  double logD3rd;
  double logDtot;
};

struct pow_matter
{
  double kmat,pmat;
 };

//Derived class for a tabulated power spectrum
class PowerSpec_Tabulated: public PowerSpec
{
    public:
    PowerSpec_Tabulated(char * FileWithTransfer, char * FileWithInputSpectrum, double Omega, double OmegaLambda, double OmegaBaryon, double OmegaNu,
                        double InputSpectrum_UnitLength_in_cm, double UnitLength_in_cm, bool no_gas, bool neutrinos_ks);
    virtual double power(double k, int Type);
    private:
        int NPowerTable;
        double scale;
        struct pow_table *PowerTable;
        struct pow_matter *PowerMatter;
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
        //public as is used in the integrator function, which is not in this class
        double R8;
    private:
        PowerSpec * PSpec;
        //Note: Uses the R8 parameter below and calls power(k, Type)
        double TopHatSigma2();
        double Norm;
        double PrimordialIndex, Dplus;
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
