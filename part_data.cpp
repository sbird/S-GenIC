#include "part_data.hpp"
#include <cassert>
//memset
#include <string.h>

part_grid::part_grid(const int NumPart_i[], const double Masses[], const double Box): NumPart(NumPart_i, (size_t)N_TYPES), Box(Box)
{
    assert(NumPart.size() == N_TYPES);
    //Mean interparticle spacing
    for(int i=0; i<N_TYPES; i++) {
        pspace[i] = NumPart[i] ? Box/NumPart[i] : 0;
    }
    //Pair up particle types
    int primary[N_TYPES/2] = {-1,-1,-1};
    int pair[N_TYPES/2] = {-1,-1,-1};
    //We want to pair each present particle type with another present type.
    //If there are an odd number of types, the last one is not paired, and the shift is zero below
    int paircnt = 0;
    for(int i=0; i < N_TYPES; i++){
        if(NumPart[i]){
           if(primary[paircnt] == -1)
              primary[paircnt] = i;
           else {
              pair[paircnt] = i;
              paircnt++;
           }
        }
    }
    //Zero the shifts
    memset(shift, 0, 3*N_TYPES*sizeof(double));
    //We shift each pair of particle types out a little,
    //so they don't have the same positions.
    //Each type will move along a different axis: this of course only works if N_TYPES <=6.
    for(int i=0; i<paircnt; i++) {
        //Move them such that M_i M_j == M_j M_i and the total displacement is less than unity.
        //Make sure to pair particles that exist together.
        double totmass = Masses[primary[i]]+Masses[pair[i]];
        assert(totmass > 0);
        shift[primary[i]][i] = -0.5*Masses[pair[i]]/totmass;
        shift[pair[i]][i] = 0.5*Masses[primary[i]]/totmass;
    }
}

/* Get the position of a particle at some index on a grid.
 * We want to return:
 * size_t index = k + j* NumPart + i * NumPart**2
 * Pos[0] = i /NumPart * pspace[type]
 * Pos[1] = j /NumPart * pspace[type]
 * Pos[2] = k /NumPart * pspace[type]
 * In other words, k is fast, j is medium, i is slow.
 */
double part_grid::Pos(size_t index, int axis, int type)
{
    //Check that we are in bounds first.
    assert(index < (size_t)NumPart[type]*NumPart[type]*NumPart[type]);
    assert(axis < 3 && axis >= 0);

    //Compute the index for this dimension
    size_t i;
    switch (axis){
        case 0:
            i = index / NumPart[type]/NumPart[type];
            break;
        case 1:
            i = index % (NumPart[type]*NumPart[type]) / NumPart[type];
            break;
        case 2:
            i = index % NumPart[type];
            break;
        default:
            i = 0;
    }
    assert(pspace[type] > 0);
    double pos = i*pspace[type] + shift[type][axis];
    return pos;
}

int part_grid::GetNumPart(int type)
{
    return NumPart[type];
}
