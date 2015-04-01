#include "proto.h"
#include "gadgetheader.h"
#include "cosmology.hpp"
#include <fftw3.h>

int FatalError(int errnum)
{
  printf("FatalError called with number=%d\n", errnum);
  fflush(stdout);
  exit(0);
}

#include "allvars.h"

void set_units(void)		/* ... set some units */
{
    InitTime = 1 / (1 + Redshift);
    UnitTime_in_s = UnitLength_in_cm / UnitVelocity_in_cm_per_s;
}


gadget_header generate_header(std::valarray<int64_t> & npart)
{
  gadget_header header;
  //No factor of h^2 because mass is in 10^10 M_sun/h
  double scale = UnitMass_in_g / pow(UnitLength_in_cm, 3) * 3* HUBBLE* HUBBLE / (8 * M_PI * GRAVITY) * pow(Box,3);
  /*Set masses*/
  for(int i = 0; i < N_TYPE; i++)
      header.mass[i] = 0;
  /*Don't forget to set masses correctly when the CDM actually incorporates other species*/
  double OmegaCDM = Omega;
  if(npart[BARYON_TYPE]){
    header.mass[BARYON_TYPE] = OmegaBaryon * scale / npart[BARYON_TYPE];
    OmegaCDM -=OmegaBaryon;
  }

  if(npart[NEUTRINO_TYPE]){
    header.mass[NEUTRINO_TYPE] = OmegaDM_2ndSpecies * scale / npart[NEUTRINO_TYPE];
#ifdef NEUTRINO_PAIRS
    header.mass[NEUTRINO_TYPE] /= 2;
#endif //NEUTRINO_PAIRS
    OmegaCDM -=OmegaDM_2ndSpecies;
  }
  /*For the "edit the transfer function" neutrino simulation method, we would *not* want to do this.
   * For true kspace neutrinos, we do. */
  else if (!neutrinos_ks){
          OmegaCDM-=OmegaDM_2ndSpecies;
  }

  if(npart[DM_TYPE])
    header.mass[DM_TYPE] = OmegaCDM * scale / npart[DM_TYPE];

  for(int i=0; i< N_TYPE; ++i){
    header.NallHW[i] = ( npart[i] >> 32);
    header.npartTotal[i] = npart[i] - ((uint64_t)header.NallHW[i] << 32);
  }


  header.time = InitTime;
  header.redshift = 1.0 / InitTime - 1;

  header.BoxSize = Box;
  header.Omega0 = Omega;
  header.OmegaLambda = OmegaLambda;
  header.HubbleParam = HubbleParam;
  /*Various flags; Most set by gadget later*/
  header.flag_sfr = 0;
  header.flag_feedback = 0;
  header.flag_cooling = 0;
  header.flag_stellarage = 0;
  header.flag_metals = 0;
  header.flag_entropy_instead_u=0;
  header.flag_doubleprecision=0;
  header.flag_ic_info=1;
  header.lpt_scalingfactor=1;
  return header;
}
