#ifdef PRINT_SPEC
#include <string>
#include "cosmology.hpp"
#include "power.hpp"

double fnl(double x, const double neff, const double gf)		/* Peacock & Dodds formula */
{
    const double A = 0.482 * pow(1 + neff / 3, -0.947);
    const double B = 0.226 * pow(1 + neff / 3, -1.778);
    const double alpha = 3.310 * pow(1 + neff / 3, -0.244);
    const double beta = 0.862 * pow(1 + neff / 3, -0.287);
    const double V = 11.55 * pow(1 + neff / 3, -0.423) * 1.2;
    return x * pow((1 + B * beta * x + pow(A * x, alpha * beta)) /
        (1 + pow(pow(A * x, alpha) * gf * gf * gf / (V * sqrt(x)), beta)), 1 / beta);
}

void print_spec(int type, PowerSpec * PSpec, Cosmology & cosmo, std::string& filename, double Redshift, double UnitLength_in_cm)
{
      FILE * fd = fopen(filename.c_str(), "w");

      const double gf = cosmo.GrowthFactor(0.001, 1.0) / (1.0 / 0.001);

      const double DDD = cosmo.GrowthFactor(1.0 / (Redshift + 1), 1.0);

      fprintf(fd, "%12g %12g\n", Redshift, DDD);	/* print actual starting redshift and 
							   linear growth factor for this cosmology */

      const double kstart = 2 * M_PI / (1000.0 * (3.085678e24 / UnitLength_in_cm));	/* 1000 Mpc/h */
      const double kend = 2 * M_PI / (0.001 * (3.085678e24 / UnitLength_in_cm));	/* 0.001 Mpc/h */

      for(double k = kstart; k < kend; k *= 1.025)
	{
	  const double po = PSpec->power(k, type);
          //printf(" po k %g %g\n ",k,po);
	  const double dl = 4.0 * M_PI * k * k * k * po;
          double knl=0, dnl=0;

	  const double kf = 0.5;

	  const double po2 = PSpec->power(1.001 * k * kf, type);
	  const double po1 = PSpec->power(k * kf, type);

	  if(po != 0 && po1 != 0 && po2 != 0)
	    {
	      const double neff = (log(po2) - log(po1)) / (log(1.001 * k * kf) - log(k * kf));
	      if(1 + neff / 3 > 0)
		{
		  dnl = fnl(dl, neff, gf);
		  knl = k * pow(1 + dnl, 1.0 / 3);
		}
            }
	  fprintf(fd, "%12g %12g %12g  %12g %12g\n", k,po, dl, knl, dnl);
	}
      fclose(fd);
}

#endif
