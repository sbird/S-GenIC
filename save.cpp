#include "save.hpp"
#include <exception>
#include <iostream>
#include <string.h>
#include <cassert>
#include "physconst.h"

using namespace GadgetWriter;
using namespace std;

#define BUFFER 48

gadget_header generate_header(std::valarray<int64_t> & npart, double Omega, double OmegaBaryon, double OmegaNuPart, double OmegaLambda, double HubbleParam, double Box, double InitTime, double UnitMass_in_g, double UnitLength_in_cm, double UnitVelocity_in_cm_per_s, bool combined_neutrinos)
{
  gadget_header header = {};
  //No factor of h^2 because mass is in 10^10 M_sun/h
  double scale = 3* HUBBLE* HUBBLE / (UnitMass_in_g / pow(UnitLength_in_cm, 3) * 8 * M_PI * GRAVITY) * pow(Box,3);
  /*Set masses*/
  for(int i = 0; i < N_TYPE; i++)
      header.mass[i] = 0;
  /*Don't forget to set masses correctly when the CDM actually incorporates other species*/
  double OmegaCDM = Omega;
  if(npart[BARYON_TYPE]) {
    header.mass[BARYON_TYPE] = OmegaBaryon * scale / npart[BARYON_TYPE];
    OmegaCDM -= OmegaBaryon;
  }

  if(npart[NEUTRINO_TYPE]){
    header.mass[NEUTRINO_TYPE] = OmegaNuPart * scale / npart[NEUTRINO_TYPE];
#ifdef NEUTRINO_PAIRS
    header.mass[NEUTRINO_TYPE] /= 2;
#endif //NEUTRINO_PAIRS
    OmegaCDM -=OmegaNuPart;
  }
  /*For the "edit the transfer function" neutrino simulation method, we would *not* want to do this.
   * For true kspace neutrinos, we do. */
  else if (!combined_neutrinos){
          OmegaCDM-=OmegaNuPart;
  }

  if(npart[DM_TYPE])
    header.mass[DM_TYPE] = OmegaCDM * scale / npart[DM_TYPE];

  for(int i=0; i< N_TYPE; ++i){
    header.NallHW[i] = ( npart[i] >> 32);
    header.npartTotal[i] = npart[i] - ((uint64_t)header.NallHW[i] << 32);
    header.npart[i] = npart[i];
  }

  header.time = InitTime;
  header.redshift = 1.0 / InitTime - 1;

  header.BoxSize = Box;
  header.Omega0 = Omega;
  header.OmegaB = OmegaBaryon;
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
  header.UnitLength_in_cm = UnitLength_in_cm;
  header.UnitMass_in_g = UnitMass_in_g;
  header.UnitVelocity_in_cm_per_s = UnitVelocity_in_cm_per_s;
  return header;
}

template <typename T> class BufferedWrite
{
    public:
        BufferedWrite(GWriteBaseSnap& snap, int64_t NumPart, int ItemsPart, const std::string groupstring) :
            snap(snap), NumPart(NumPart), ItemsPart(ItemsPart), groupstring(groupstring), blockmaxlen(BUFFER * 1024 * 1024)
        {
            if(!(block = (T *) malloc(blockmaxlen*ItemsPart*sizeof(T))))
                throw std::ios_base::failure("Failed to allocate "+std::to_string(ItemsPart*sizeof(float)*blockmaxlen/1024/1024)+" MB for write buffer");
        }
        ~BufferedWrite()
        {
            free(block);
        }
        int64_t do_write(int type, void * data, int64_t np_write, int64_t begin)
        {
            int64_t retval;
            try{
                retval = dynamic_cast<GWriteSnap&>(snap).WriteBlocks(groupstring, type, data, np_write, begin);
            } catch (const std::bad_cast & e) {
#ifdef HAVE_BGFL
                try {
                retval = dynamic_cast<GWriteBigSnap&>(snap).WriteBlocks(groupstring, type, data, np_write, begin, dtype().c_str(), ItemsPart);
                } catch(const std::bad_cast & e) {
#endif
                std::cout << e.what() << '\n';
                exit(1);
#ifdef HAVE_BGFL
                }
#endif
            }
            return retval;
        }
        int writeparticles(int type)
        {
            int64_t written=0, pc = 0;
            for(int64_t i = 0; i < NumPart; i++){
                for(int k = 0; k < ItemsPart; k++){
                    assert(ItemsPart*pc + k < blockmaxlen*ItemsPart);
                    block[ItemsPart * pc + k] = setter(i,k,type);
                }
                pc++;
#ifdef NEUTRINO_PAIRS
                /*Add an extra copy of the position vector for the double neutrino*/
                if(type == NEUTRINO_TYPE) {
                    for(int k = 0; k < ItemsPart; k++)
                      block[ItemsPart * pc + k] = setter(i, k, type);
                    pc++;
                }
#endif //NEUTRINO_PAIRS
                if(pc >= blockmaxlen){
                  if(do_write(type, block, pc,written) != pc)
                      throw std::ios_base::failure("Could not write data at particle "+std::to_string(i));
                  written+=pc;
                  pc = 0;
                }
            }
            if(pc > 0 && do_write(type, block, pc,written) != pc)
                  throw std::ios_base::failure("Could not write final data for: "+groupstring);
            return written;
       }
    protected:
       virtual double setter(int i, int k, int type) = 0;
    private:
        std::string dtype(void) {
            throw std::runtime_error("Need to specialise dtype for class");
        }
        T * block;
        GWriteBaseSnap& snap;
        const int64_t NumPart;
        const int ItemsPart;
        const std::string groupstring;
        const int64_t blockmaxlen;
};

template <> std::string BufferedWrite<float>::dtype(void) {
    return "f4";
};
template <> std::string BufferedWrite<double>::dtype(void) {
    return "f8";
};
template <> std::string BufferedWrite<id_type>::dtype(void) {
        return "i"+std::to_string(sizeof(id_type));
};

class PosBufferedWrite : public BufferedWrite<float>
{
    public:
        PosBufferedWrite(GWriteBaseSnap& snap, int64_t NumPart, lpt_data * outdata, part_grid & Pgrid) :
            BufferedWrite(snap, NumPart, 3, name(snap.GetFormat())),
            Pgrid(Pgrid), outdata(outdata)
            {}
    private:
        std::string name(int format)
        {
            switch(format)
            {
                case 3:
                    return "Coordinates";
                case 4:
                    return "Position";
                default:
                    return "POS ";
            }
        }
        virtual double setter(int i, int k, int type)
        {
          double value = Pgrid.Pos(i,k, type);
          if(outdata)
            value += outdata->GetDisp(i,k);
          value = periodic_wrap(value, Pgrid.GetBox());
          return value;
        }
        part_grid & Pgrid;
        lpt_data * outdata;
};


class VelBufferedWrite : public BufferedWrite<float>
{
    public:
        VelBufferedWrite(GWriteBaseSnap& snap, int64_t NumPart, FermiDiracVel * therm_vels, lpt_data * outdata) :
            BufferedWrite(snap, NumPart, 3, name(snap.GetFormat())),
            therm_vels(therm_vels), outdata(outdata)
            {
              memset(vtherm, 0, 3);
            }
    private:
        std::string name(int format)
        {
            switch(format)
            {
                case 3:
                    return "Velocities";
                case 4:
                    return "Velocity";
                default:
                    return "VEL ";
            }
        }
        virtual double setter(int i, int k, int type)
        {
          if(k == 0)
              get_new_therm_vels();
          assert(k >= 0 || k <= 2);
          double value = vtherm[k];
          if(outdata)
            value += outdata->GetVel(i,k);
          return value;
        }
          void get_new_therm_vels()
          {
              memset(vtherm, 0, 3);
              if(!therm_vels)
                  return;
              therm_vels->add_thermal_speeds(vtherm);
          }
        FermiDiracVel *therm_vels;
        lpt_data * outdata;
        float vtherm[3];
};

class IDBufferedWrite : public BufferedWrite<id_type>
{
    public:
        IDBufferedWrite(GWriteBaseSnap& snap, int64_t NumPart, int64_t FirstId) :
            BufferedWrite(snap, NumPart, 1, name(snap.GetFormat())),
#ifdef NEUTRINO_PAIRS
            sw(0),
#endif
            FirstId(FirstId)
            {
            }
    private:
        std::string name(int format)
        {
            switch(format)
            {
                case 3:
                    return "ParticleIDs";
                case 4:
                    return "ID";
                default:
                    return "ID  ";
            }
        }
        virtual double setter(int i, int k, int type)
        {
#ifdef NEUTRINO_PAIRS
            if(type == NEUTRINO_TYPE) {
            sw != sw;
            return 2*(i+ FirstId) + sw;
            }
            else
#else
            return i + FirstId;
#endif
        }
#ifdef NEUTRINO_PAIRS
            bool sw;
#endif
            const int64_t FirstId;
};

/*Class to write zero energies*/
class EnergyBufferedWrite : public BufferedWrite<float>
{
    public:
    EnergyBufferedWrite(GWriteBaseSnap& snap, int64_t NumPart) :
        BufferedWrite(snap, NumPart, 1, snap.GetFormat() > 2 ? "InternalEnergy" : "U   ")
        {}
    private:
    virtual double setter(int i, int k, int type)
    {
        return 0;
    }
};

int64_t write_particle_data(GWriteBaseSnap& snap, int type, lpt_data * outdata, part_grid& Pgrid, FermiDiracVel *therm_vels, int64_t FirstId)
{
  const int64_t NumPart = Pgrid.GetNumPart(type)*Pgrid.GetNumPart(type)*Pgrid.GetNumPart(type);
  printf("\nWriting IC-file\n");
  {
    PosBufferedWrite pp(snap, NumPart, outdata, Pgrid);
    pp.writeparticles(type);
  }
  {
    VelBufferedWrite pp(snap, NumPart, therm_vels, outdata);
    pp.writeparticles(type);
  }
  {
    IDBufferedWrite pp(snap, NumPart, FirstId);
    pp.writeparticles(type);
  }
  {
    EnergyBufferedWrite pp(snap, NumPart);
    pp.writeparticles(type);
  }
  printf("Finished writing IC file.\n");
  FirstId+=NumPart;
#ifdef NEUTRINO_PAIRS
  if(type==2)
          FirstId+=NumPart;
#endif
  return FirstId;
}

double periodic_wrap(double x, double box)
{
  x = fmod(x, box);

  if (x < 0)
    x += box;

  return x;
}

