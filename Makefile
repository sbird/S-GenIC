EXEC   = N-GenIC

OBJS   = main.o power.o allvars.o save.o read_param.o  read_glass.o  \
         nrsrc/nrutil.o nrsrc/qromb.o nrsrc/polint.o nrsrc/trapzd.o

INCL   = allvars.h proto.h  nrsrc/nrutil.h  Makefile



#OPT   +=  -DPRODUCEGAS   # Set this to automatically produce gas particles 
                         # for a single DM species in the input file by interleaved by a half a grid spacing


OPT   +=  -DMULTICOMPONENTGLASSFILE  # set this if the initial glass file contains multiple components

OPT   +=  -DDIFFERENT_TRANSFER_FUNC  # set this if you want to implement a transfer function that depends on
                                     # particle type

OPT    +=  -DNO64BITID    # switch this on if you want normal 32-bit IDs
#OPT   +=  -DCORRECT_CIC  # only switch this on if particles are homogenously distributed over mesh cells (say glass)

OPT   +=  -DNEUTRINOS  # this will produce a second component as slight neutrinos (needs to be in initial glass)
#OPT   +=  -DNEUTRINO_PAIRS  # this will produce an additional partner for every neutrino with opposite thermal velocities


OPTIONS =  $(OPT)

#SYSTYPE="Stella"
#SYSTYPE="Regatta"
#SYSTYPE="RZG_LinuxCluster"
#SYSTYPE="RZG_LinuxCluster-gcc"
#SYSTYPE="Solaris"

# module load mvapich2-1.2-sdr-intel/11.0
#SYSTYPE="OPA-Cluster64-Intel"

# module load mvapich2-1.2-sdr-gnu/4.1.2
#SYSTYPE="OPA-Cluster64-Gnu"



#FFTW_INCL = -I/usr/common/pdsoft/include
#FFTW_LIBS = -L/usr/common/pdsoft/lib


CC       =   mpicc        # sets the C-compiler (default)
OPTIMIZE =   -O2 -w -Wall    # optimization and warning flags (default)
#MPICHLIB =  -lmpich




ifeq ($(SYSTYPE),"Stella")
CC       =  mpicc
OPTIMIZE =  -O3 -Wall
GSL_INCL =  -I/home/schaye/libs/include
GSL_LIBS =  -L/home/schaye/libs/lib -static
FFTW_INCL=  -I/home/schaye/libs/include
FFTW_LIBS=  -L/home/schaye/libs/lib
MPICHLIB =
endif


ifeq ($(SYSTYPE),"OPA-Cluster64-Intel")  # module load mvapich2-1.2-sdr-intel/11.0
CC       =   mpicc
ifeq (SOFTDOUBLEDOUBLE,$(findstring SOFTDOUBLEDOUBLE,$(OPT))) 
CC       =   mpicxx
endif
OPTIMIZE =  -O3 -g -Wall -m64
GSL_INCL =  -I/afs/rzg/home/v/vrs/Libs/opteron64/mvapich2-1.2-sdr-intel-11.0/include
GSL_LIBS =  -L/afs/rzg/home/v/vrs/Libs/opteron64/mvapich2-1.2-sdr-intel-11.0/lib  -Xlinker -R -Xlinker /afs/rzg/home/v/vrs/Libs/opteron64/mvapich2-1.2-sdr-intel-11.0/lib
FFTW_INCL=  -I/afs/rzg/home/v/vrs/Libs/opteron64/mvapich2-1.2-sdr-intel-11.0/include
FFTW_LIBS=  -L/afs/rzg/home/v/vrs/Libs/opteron64/mvapich2-1.2-sdr-intel-11.0/lib
MPICHLIB =
HDF5INCL =  
HDF5LIB  =  
OPT      +=  -DNOCALLSOFSYSTEM
endif


ifeq ($(SYSTYPE),"OPA-Cluster64-Gnu")  # module load mvapich2-1.2-sdr-gnu/4.1.2
CC       =   mpicc
ifeq (SOFTDOUBLEDOUBLE,$(findstring SOFTDOUBLEDOUBLE,$(OPT))) 
CC       =   mpicxx
endif
OPTIMIZE =  -O3 -g -Wall -m64
GSL_INCL =  -I/afs/rzg/home/v/vrs/Libs/opteron64/mvapich2-1.2-sdr-gnu-4.1.2/include
GSL_LIBS =  -L/afs/rzg/home/v/vrs/Libs/opteron64/mvapich2-1.2-sdr-gnu-4.1.2/lib  -Xlinker -R -Xlinker /afs/rzg/home/v/vrs/Libs/opteron64/mvapich2-1.2-sdr-gnu-4.1.2/lib
FFTW_INCL=  -I/afs/rzg/home/v/vrs/Libs/opteron64/mvapich2-1.2-sdr-gnu-4.1.2/include
FFTW_LIBS=  -L/afs/rzg/home/v/vrs/Libs/opteron64/mvapich2-1.2-sdr-gnu-4.1.2/lib
MPICHLIB =
HDF5INCL =  
HDF5LIB  =  
OPT      +=  -DNOCALLSOFSYSTEM
endif

ifeq ($(SYSTYPE),"Regatta")
CC       =   mpcc_r   
OPTIMIZE =   -O5 -qstrict -qipa -q64
GSL_INCL = -I/afs/ipp-garching.mpg.de/u/vrs/gsl_psi64/include
GSL_LIBS = -L/afs/ipp-garching.mpg.de/u/vrs/gsl_psi64/lib                
FFTW_INCL= -I/afs/ipp-garching.mpg.de/u/vrs/fftw_psi64/include
FFTW_LIBS= -L/afs/ipp-garching.mpg.de/u/vrs/fftw_psi64/lib  -q64 -qipa
MPICHLIB =
endif

ifeq ($(SYSTYPE),"RZG_LinuxCluster")
CC       =   mpicci   
OPTIMIZE =   -O3 
GSL_INCL = -I/afs/ipp-garching.mpg.de/u/vrs/gsl_linux/include
GSL_LIBS = -L/afs/ipp-garching.mpg.de/u/vrs/gsl_linux/lib                -static
FFTW_INCL= -I/afs/ipp-garching.mpg.de/u/vrs/fftw_linux/include
FFTW_LIBS= -L/afs/ipp-garching.mpg.de/u/vrs/fftw_linux/lib
endif

ifeq ($(SYSTYPE),"RZG_LinuxCluster-gcc")
CC       =   /afs/ipp-garching.mpg.de/u/vrs/bin/mpiccf   
OPTIMIZE =   -O3 
GSL_INCL = -I/afs/ipp-garching.mpg.de/u/vrs/gsl_linux_gcc3.2/include
GSL_LIBS = -L/afs/ipp-garching.mpg.de/u/vrs/gsl_linux_gcc3.2/lib
FFTW_INCL= -I/afs/ipp-garching.mpg.de/u/vrs/fftw_linux_gcc3.2/include
FFTW_LIBS= -L/afs/ipp-garching.mpg.de/u/vrs/fftw_linux_gcc3.2/lib  
endif

ifeq ($(SYSTYPE),"Solaris")
CC       =   mpcc   # sets the C-compiler
OPTIMIZE =   -i -fast -xvector -xarch=v9b -xchip=ultra3 -xcache=64/32/4:8192/512/1 -I/opt/local/include

GSL_INCL = -I/opt/local/include/gsl
GSL_LIBS = -L/opt/local/lib/sparcv9               
FFTW_INCL= -I/opt/local/include
FFTW_LIBS= -L/opt/local/lib/sparcv9
endif



FFTW_LIB =  $(FFTW_LIBS) -ldrfftw_mpi -ldfftw_mpi -ldrfftw -ldfftw

LIBS   =   -lm  $(MPICHLIB)  $(FFTW_LIB)  $(GSL_LIBS)  -lgsl -lgslcblas

ifeq ($(SYSTYPE),"Solaris")
LIBS   =   -R/opt/local/lib/sparcv9 -lm  -lmpi   $(GSL_LIBS) -lgsl -lgslcblas  $(FFTW_LIB)
endif



CFLAGS =   $(OPTIONS)  $(OPTIMIZE)  $(FFTW_INCL) $(GSL_INCL)

$(EXEC): $(OBJS) 
	$(CC) $(OPTIMIZE) $(OBJS) $(LIBS)   -o  $(EXEC)  

$(OBJS): $(INCL) 


.PHONY : clean
clean:
	rm -f $(OBJS) $(EXEC)



