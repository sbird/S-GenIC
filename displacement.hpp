#ifndef DISPLACEMENT_H
#define DISPLACEMENT_H

#include <fftw3.h>
#include <stdint.h>

class part_data;
class PowerSpec;

class DisplacementFields
{
    public:
        DisplacementFields(size_t Nmesh, size_t Nsample, int Seed, double Box, bool twolpt=true);
        void displacement_fields(const int type, const int64_t NumPart, part_data& P, PowerSpec * PSpec, bool SphereMode, bool RayleighScatter=true);
        ~DisplacementFields();
    private:
        //Takes the computed Fourier displacement field and reads them into part_data P
        double displacement_read_out(const int order, const int64_t NumPart, part_data& P, const int axes);
        //Size of the Fourier grid
        const size_t Nmesh;
        //Nsample is the number of particles
        const size_t Nsample;
        //Do we need a twolpt term?
        const bool twolpt;
        //Seed for the random number generator.
        //Constant for all particle types, so they are coherent.
        const int Seed;
        //Box size
        const double Box;
        //FFT variables for Zeldovich
        fftwf_plan Inverse_plan;
        float *Disp;
        //This will always be a cast of Disp
        fftwf_complex *Cdata;
        //Pointers for 2LPT term
        fftwf_plan Forward_plan2;
        fftwf_plan Inverse_plan_grad[3];
        float *twosrc;
        fftwf_complex *ctwosrc;
        fftwf_complex *(cdigrad[3]);
        float *(digrad[3]);
};

#endif
