#ifndef __PART_DATA_HPP
#define __PART_DATA_HPP

#include <math.h>
#include <stdlib.h>
#include <valarray>
#include <iostream>
#include <gadgetreader.hpp>

#define FLOAT_TYPE float
#define N_TYPES 6

//Class to generate a regular grid of particle positions on demand.
class part_grid
{
    public:
        /* Construct the particle grid parameters.
         * NumPart should be the number of particles per side; the cube root of the total number of particles.
         * Masses are the mass ratios of the different particle species, for making sure no two types share the same position
         * Box is the boxsize.*/
        part_grid(const std::valarray<size_t>& NumPart, const std::valarray<double>& Masses, double Box);
        /* Get the position of a particle at some index on a grid.
         * We want to return:
         * size_t index = k + j* NumPart + i * NumPart**2
         * Pos[0] = i /NumPart * pspace[type]
         * Pos[1] = j /NumPart * pspace[type]
         * Pos[2] = k /NumPart * pspace[type]
         * In other words, k is fast, j is medium, i is slow.
         */
        double Pos(size_t index, int axis, int type);
    private:
        const std::valarray<size_t> NumPart;
        double pspace[N_TYPES];
        double shift[N_TYPES][3];
};

//Class to store the Zeldovich and 2LPT displacements for a particle type.
class lpt_data
{
    public:
        lpt_data(const size_t NumPart, bool twolpt=true): twolpt(twolpt),Vel_data(3*NumPart),Vel_2_data(twolpt*3*NumPart)
        {}
        //Note Vel_2_data takes zero space if !twolpt
        //Note these functions are not checked. In particular calling the Vel2 functions will be bad if !twolpt
        inline double Vel(size_t index, int axis){
                return Vel_data[3*index+axis];
        }
        inline void SetVel(FLOAT_TYPE Vel_in, size_t index, int axis){
                Vel_data[3*index+axis] = Vel_in;
                return;
        }
        inline double Vel2(size_t index, int axis){
                return Vel_2_data[3*index+axis];
        }
        inline void Set2Vel(FLOAT_TYPE Vel_in, size_t index, int axis){
                Vel_2_data[3*index+axis] = Vel_in;
                return;
        }
        inline size_t GetNumPart(){
                return Vel_data.size()/3;
        }
    private:
        const bool twolpt;
        std::valarray <FLOAT_TYPE> Vel_data;
        std::valarray <FLOAT_TYPE> Vel_2_data;
};

class part_data{
    public:
        part_data(GadgetReader::GSnap& snap, int type, int GlassTileFac, double Box, bool twolpt=true): twolpt(twolpt),
        Nglass(snap.GetNpart(type)), IGlassBox(1./snap.GetHeader().BoxSize), NglassGTile(Nglass*GlassTileFac),
        NglassG2Tile(Nglass*GlassTileFac*GlassTileFac), BoxGTile(Box/GlassTileFac),
        pos(new float[Nglass * 3]),
        Vel_data(3*Nglass*GlassTileFac*GlassTileFac*GlassTileFac),
        //Note this is zero space if !twolpt
        Vel_2_data(twolpt*3*Nglass*GlassTileFac*GlassTileFac*GlassTileFac)
        {
          /*Read all POS data for this type*/
          if(snap.GetBlock("POS ",pos,Nglass,0,(1<<N_TYPE)-1-(1<<type)) != Nglass) {
                  std::cerr<<"Error reading particle data"<<std::endl;
                  throw std::bad_alloc();
          }
          std::cout<<"Type "<<type<<" allocated "<<sizeof(float)*(Vel_data.size()+Vel_2_data.size())/1024/1024<<" MB for "<<GetNumPart()<<" particles"<<std::endl;
          return;
        }
        
         ~part_data(){
                delete[] pos;
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
        
        //Note these functions are not checked. In particular calling the Vel2 functions will be bad if !twolpt
        inline double Vel(size_t index, int axis){
                return Vel_data[3*index+axis];
        }
        
        inline void SetVel(double Vel_in, size_t index, int axis){
                Vel_data[3*index+axis] = Vel_in;
                return;
        }

        inline double Vel2(size_t index, int axis){
                return Vel_2_data[3*index+axis];
        }
        
        inline void Set2Vel(double Vel_in, size_t index, int axis){
                Vel_2_data[3*index+axis] = Vel_in;
                return;
        }

        inline int64_t GetNumPart(){
                return Vel_data.size()/3;
        }

    private:
        const bool twolpt;
        const int64_t Nglass;
        const double IGlassBox;
        const int64_t NglassGTile;
        const int64_t NglassG2Tile;
        const double BoxGTile;
        float *pos;
        std::valarray <float> Vel_data;
        std::valarray <float> Vel_2_data;
};

#endif

