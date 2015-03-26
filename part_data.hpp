#ifndef __PART_DATA_HPP
#define __PART_DATA_HPP

#include <math.h>
#include <stdlib.h>
#include <valarray>
#include <gadgetreader.hpp>

class part_data{
    
    public:
        part_data(GadgetReader::GSnap& snap, int type, int GlassTileFac, double Box):
        Nglass(snap.GetNpart(type)), IGlassBox(1./snap.GetHeader().BoxSize), NglassGTile(Nglass*GlassTileFac),
        NglassG2Tile(Nglass*GlassTileFac*GlassTileFac), BoxGTile(Box/GlassTileFac),
#ifdef TWOLPT
#ifdef NEUTRINOS
        /*Do not allocate the 2lpt array for neutrinos*/
        Vel_2_data((type != 2)*3*Nglass*GlassTileFac*GlassTileFac*GlassTileFac),
#else
        Vel_2_data(3*Nglass*GlassTileFac*GlassTileFac*GlassTileFac),
#endif
#endif
        Vel_data(3*Nglass*GlassTileFac*GlassTileFac*GlassTileFac)
        {
          pos = NULL;
          
          if(!(pos = (float *) malloc(sizeof(float) * Nglass * 3))){
                  fprintf(stderr,"failed to allocate %ld Mbyte for glass file\n", sizeof(float) * Nglass * 3 / (1024 * 1024));
        		  exit(112);
          }
          /*Read all POS data for this type*/
          if(snap.GetBlock("POS ",pos,Nglass,0,(1<<N_TYPE)-1-(1<<type)) != Nglass){
                  fprintf(stderr, "Error reading particle data\n");
                  free(pos);
                  exit(113);
          }
#ifndef TWOLPT
          printf("Type %d allocated %ld MB for %ld particles\n",type,sizeof(float)*Vel_data.size()/1024/1024, GetNumPart());
#else
          printf("Type %d allocated %ld MB for %ld particles\n",type,sizeof(float)*(Vel_data.size()+Vel_2_data.size())/1024/1024, GetNumPart());
#endif
          return;
        }
        
         ~part_data(){
            if(pos)
                free(pos);
            return;
         }

        double Pos(size_t index, int axis){
           /* We want to return something equivalent to:
            * size_t index = n+k*Nglass+j*Nglass*GlassTileFac+i*G2Tile*Nglass;
            * n == index % Nglass.
            * k = (index %Nglass*GlassTileFac) /Nglass
            * j = (index % Nglass*G2Tile) / (Nglass*GlassTileFac)
            * i = index / G2Tile*Nglass 
            * P[index].Pos[0] = pos[3 * n] / GlassBox * (Box / GlassTileFac) + i * (Box / GlassTileFac);
            * P[index].Pos[1]= pos[3 * n + 1] / GlassBox * (Box / GlassTileFac) + j * (Box / GlassTileFac);
            * P[index].Pos[2]  = pos[3 * n + 2] / GlassBox * (Box / GlassTileFac) + k * (Box / GlassTileFac);
            */
            size_t n = index % Nglass;
            int i;
            switch (axis){
                case 0:
                    i = index / NglassG2Tile;
                    break;
                case 1:
                    i = index % NglassG2Tile / NglassGTile;
                    break;
                case 2:
                    i = index % NglassGTile / Nglass;
                    break;
                default: 
                    i = 0;
            }
            return pos[3 * n + axis] * IGlassBox * BoxGTile + i*BoxGTile;
        }
        
        inline double Vel(size_t index, int axis){
                return Vel_data[3*index+axis];
        }
        
        inline void SetVel(double Vel_in, size_t index, int axis){
                Vel_data[3*index+axis] = Vel_in;
                return;
        }

#ifdef TWOLPT
        inline double Vel2(size_t index, int axis){
                return Vel_2_data[3*index+axis];
        }
        
        inline void Set2Vel(double Vel_in, size_t index, int axis){
                Vel_2_data[3*index+axis] = Vel_in;
                return;
        }
#endif


        inline int64_t GetNumPart(){
                return Vel_data.size()/3;
        }

    private:
        const int64_t Nglass;
        const double IGlassBox;
        const int64_t NglassGTile;
        const int64_t NglassG2Tile;
        const double BoxGTile;
        float *pos;
#ifdef TWOLPT
        std::valarray <float> Vel_2_data;
#endif
        std::valarray <float> Vel_data;

};

#endif

