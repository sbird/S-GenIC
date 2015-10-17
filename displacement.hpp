#ifndef DISPLACEMENT_H
#define DISPLACEMENT_H

#include <fftw3.h>
#include <stdint.h>

class part_data;
class PowerSpec;

class DisplacementFields
{
    public:
        DisplacementFields(size_t Nmesh, int Seed, double Box, bool twolpt=true);
        void displacement_fields(const int type, part_data& P, PowerSpec * PSpec, bool SphereMode=true, bool RayleighScatter=true);
        ~DisplacementFields();
    private:
        //Takes the computed Fourier displacement field and reads them into part_data P
        double displacement_read_out(const int order, part_data& P, const int axes);
        //Size of the Fourier grid
        const size_t Nmesh;
        //Do we need a twolpt term?
        const bool twolpt;
        //Seed for the random number generator.
        //Constant for all particle types, so they are coherent.
        const int Seed;
        //Box size
        const double Box;
        //FFT variables for Zeldovich
        fftw_plan Inverse_plan;
        double *Disp;
        //This will always be a cast of Disp
        fftw_complex *Cdata;
        //Pointers for 2LPT term
        fftwf_plan Forward_plan2;
        fftwf_plan Inverse_plan_grad[3];
        float *twosrc;
        fftwf_complex *ctwosrc;
        fftwf_complex *(cdigrad[3]);
        float *(digrad[3]);
};

#endif
