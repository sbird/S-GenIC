#ifndef __PART_DATA_HPP
#define __PART_DATA_HPP

#include <valarray>

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
        part_grid(const int64_t NumPart[], const double Masses[], const double Box);
        /* Get the position of a particle at some index on a grid.
         * We want to return:
         * size_t index = k + j* NumPart + i * NumPart**2
         * Pos[0] = i /NumPart * pspace[type]
         * Pos[1] = j /NumPart * pspace[type]
         * Pos[2] = k /NumPart * pspace[type]
         * In other words, k is fast, j is medium, i is slow.
         */
        double Pos(size_t index, int axis, int type);
        //Get number of particles
        int64_t GetNumPart(int type);
        //Get box size
        inline double GetBox(){
            return Box;
        }
    private:
        const std::valarray<int64_t> NumPart;
        const double Box;
        double pspace[N_TYPES];
        double shift[N_TYPES][3];
};

//Class to store the Zeldovich and 2LPT displacements for a particle type.
class lpt_data
{
    public:
        //Again, NumPart is the cube root of the particle number!
        lpt_data(const int NumPart, bool twolpt=true): vel_prefac(0), vel_prefac2(0), twolpt(twolpt),Vel_data(3*static_cast<size_t>(NumPart)*NumPart*NumPart),Vel_2_data(twolpt*3*static_cast<size_t>(NumPart)*NumPart*NumPart)
        {}
        //Note Vel_2_data takes zero space if !twolpt
        inline void SetVel(FLOAT_TYPE Vel_in, size_t index, int axis){
                Vel_data[3*index+axis] = Vel_in;
                return;
        }
        inline void Set2Vel(FLOAT_TYPE Vel_in, size_t index, int axis){
                Vel_2_data[3*index+axis] = Vel_in;
                return;
        }
        inline double GetVel(size_t index, int axis)
        {
            double vel = vel_prefac*Vel_data[3*index+axis];
            if(twolpt)
                vel += vel_prefac2*Vel_2_data[3*index+axis];
            return vel;
        }
        inline double GetDisp(size_t index, int axis)
        {
            double disp = Vel_data[3*index+axis];
            if(twolpt)
                disp -= 3./7.*Vel_2_data[3*index+axis];
            return disp;
        }
        inline size_t GetNumPart(){
                return Vel_data.size()/3;
        }
        void SetVelPrefac(double vel_prefac_i, double vel_prefac2_i)
        {
            vel_prefac = vel_prefac_i;
            vel_prefac2 = vel_prefac2_i;
        }
    private:
        double vel_prefac;
        double vel_prefac2;
        const bool twolpt;
        std::valarray <FLOAT_TYPE> Vel_data;
        std::valarray <FLOAT_TYPE> Vel_2_data;
};

#endif

