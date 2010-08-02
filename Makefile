SYSTYPE="#SYSTYPE#"
#SYSTYPE="Regatta"
#SYSTYPE="JUMP"
#SYSTYPE="LOUHI"
#SYSTYPE="MPA"
#SYSTYPE="Stella"
#SYSTYPE="RZG_LinuxCluster"
#SYSTYPE="RZG_LinuxCluster-gcc"
#SYSTYPE="Solaris"
#SYSTYPE="OPA-Cluster64"


EXEC   = N-GenIC

OBJS   = main.o power.o allvars.o save.o read_param.o  read_glass.o  \
         nrsrc/nrutil.o nrsrc/qromb.o nrsrc/polint.o nrsrc/trapzd.o cubspl.o

INCL   = allvars.h proto.h  nrsrc/nrutil.h  Makefile



#OPT   +=  -DPRODUCEGAS   # Set this to automatically produce gas particles 
                         # for a single DM species in the input file by interleaved by a half a grid spacing


OPT   += -DMULTICOMPONENTGLASSFILE # set this if the initial glass file contains multiple components

OPT   += -DDIFFERENT_TRANSFER_FUNC  # set this if you want to implement a transfer function that depends on particle type. Or for tk_CAMB to work.
OPT	+= -DFORMAT_TWO  #Set this if you want to output IC files in format 2.												

OPT   += -DNO64BITID # switch this on if you want normal 32-bit IDs
#OPT   +=  -DCORRECT_CIC  # only switch this on if particles are homogenously distributed over mesh cells (say glass)
OPT	+= -DDOUBLEPRECISION_FFTW


OPTIONS = $(OPT)


#gcc
#CC       =  mpicc  # sets the C-compiler (default)
#FC			=  mpif90
#OPTIMIZE =  -O3 -Wall -Wno-strict-aliasing   # optimization and warning flags (default)
#MPICHLIB = # -lmpich
#FFTW_INCL=  -I/data/store/spb41/apps/fftw/include
#FFTW_LIBS=  -L/data/store/spb41/apps/fftw/lib
#

CC       =  mpicc  # sets the C-compiler (default)
FC			=  mpif90
OPTIMIZE =  -O2 -Wall
MPICHLIB = # -lmpich
FFTW_INCL=  
FFTW_LIBS=  

ifeq ($(SYSTYPE),"OPA-Cluster64")
CC       =   mpiccg   
OPTIMIZE =  -O2 -Wall -m64
GSL_INCL =  -I/afs/rzg/bc-b/vrs/opteron64/include
GSL_LIBS =  -L/afs/rzg/bc-b/vrs/opteron64/lib  -Wl,"-R /afs/rzg/bc-b/vrs/opteron64/lib"
FFTW_INCL=  -I/afs/rzg/bc-b/vrs/opteron64/include
FFTW_LIBS=  -L/afs/rzg/bc-b/vrs/opteron64/lib 
endif

ifeq ($(SYSTYPE),"Stella")
CC       =  mpicc
OPTIMIZE =  -O3 -Wall
GSL_INCL =  -I/home/schaye/libs/include
GSL_LIBS =  -L/home/schaye/libs/lib -static
FFTW_INCL=  -I/home/schaye/libs/include
FFTW_LIBS=  -L/home/schaye/libs/lib
MPICHLIB =
endif

ifeq ($(SYSTYPE),"MPA")
CC       =  mpicc
OPTIMIZE =  -O3 -Wall
FFTW_INCL = -I/usr/common/pdsoft/include
FFTW_LIBS = -L/usr/common/pdsoft/lib
GSL_INCL =  -I/usr/common/pdsoft/include
GSL_LIBS =  -L/usr/common/pdsoft/lib
endif

ifeq ($(SYSTYPE),"Regatta")
CC       =   mpcc_r   
OPTIMIZE =   -O5 -qstrict -qipa -q64
GSL_INCL = -I/afs/ipp-garching.mpg.de/u/vrs/gsl_psi64/include -I/deisa/rzg/home/rzg00rzg/rzg0schm/DIR/psi/include
GSL_LIBS = -L/afs/ipp-garching.mpg.de/u/vrs/gsl_psi64/lib -L/deisa/rzg/home/rzg00rzg/rzg0schm/DIR/psi/lib 
FFTW_INCL= -I/afs/ipp-garching.mpg.de/u/vrs/fftw_psi64/include -I/usr/local/fftw/LOCALinclude
FFTW_LIBS= -L/afs/ipp-garching.mpg.de/u/vrs/fftw_psi64/lib  -q64 -qipa -L/usr/local/fftw/LOCALlib
MPICHLIB =
endif

ifeq ($(SYSTYPE),"JUMP")
CC       =   mpcc_r
OPTIMIZE =   -O5 -qstrict -qipa -q64
GSL_INCL = -I/afs/ipp-garching.mpg.de/u/vrs/gsl_psi64/include -I/deisa/rzg/home/rzg00rzg/rzg0schm/DIR/jump/include
GSL_LIBS = -L/afs/ipp-garching.mpg.de/u/vrs/gsl_psi64/lib -L/deisa/rzg/home/rzg00rzg/rzg0schm/DIR/jump/lib
FFTW_INCL= -I/afs/ipp-garching.mpg.de/u/vrs/fftw_psi64/include -I/usr/local/fftw/LOCALinclude
FFTW_LIBS= -L/afs/ipp-garching.mpg.de/u/vrs/fftw_psi64/lib  -q64 -qipa -L/usr/local/fftw/LOCALlib
MPICHLIB =
endif

ifeq ($(SYSTYPE),"NBENCH")
CC       =  #MPI_CC#
OPTIMIZE =  #CFLAGS#
GSL_INCL = -I#OUTDIR#/../../../run/GSL/include
GSL_LIBS = -L#OUTDIR#/../../../run/GSL/lib
FFTW_INCL= #FFT_INC#
FFTW_LIBS= #FFT_DIR#
OPT     += #FFTW_TYPE#
MPICHLIB =
endif

ifeq ($(SYSTYPE),"LOUHI")
CC       =  cc -fastsse -Mipa=fast,inline
OPTIMIZE =  -O3
GSL_INCL = -I/deisa/rzg/home/rzg00rzg/rzg0schm/DIR/louhi/include
GSL_LIBS = -L/deisa/rzg/home/rzg00rzg/rzg0schm/DIR/louhi/lib
FFTW_INCL= -I/opt/fftw/2.1.5/cnos/include
FFTW_LIBS= -L/opt/fftw/2.1.5/cnos/lib
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
CC       =   mpiccg   
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

ifeq ($(SYSTYPE),"CINECA-BCX")
CC       =  #MPI_CC#
OPTIMIZE =  #CFLAGS#
GSL_INCL = -I#OUTDIR#/../../../run/GSL/include
GSL_LIBS = -L#OUTDIR#/../../../run/GSL/lib
FFTW_INCL= #FFT_INC#
FFTW_LIBS= #FFT_DIR#
OPT     += #FFTW_TYPE#
MPICHLIB =
endif

ifeq (NOTYPEPREFIX_FFTW,$(findstring NOTYPEPREFIX_FFTW,$(OPT)))    # fftw installed with type prefix?
  FFTW_LIB = $(FFTW_LIBS) -lrfftw_mpi -lfftw_mpi -lrfftw -lfftw
else
ifeq (DOUBLEPRECISION_FFTW,$(findstring DOUBLEPRECISION_FFTW,$(OPT)))
  FFTW_LIB = $(FFTW_LIBS) -ldrfftw_mpi -ldfftw_mpi -ldrfftw -ldfftw
else
  FFTW_LIB = $(FFTW_LIBS) -lsrfftw_mpi -lsfftw_mpi -lsrfftw -lsfftw
endif
endif

LIBS   =   -lm  $(MPICHLIB)  $(FFTW_LIB)  $(GSL_LIBS)  -lgsl -lgslcblas

ifeq ($(SYSTYPE),"Solaris")
LIBS   =   -R/opt/local/lib/sparcv9 -lm  -lmpi   $(GSL_LIBS) -lgsl -lgslcblas  $(FFTW_LIB)
endif



CFLAGS = $(OPTIONS) $(OPTIMIZE) $(FFTW_INCL) $(GSL_INCL)

$(EXEC): $(OBJS) 
	$(CC) $(OPTIMIZE) $(OBJS) $(LIBS)   -o  $(EXEC)  

cubspl.o: cubspl.f
	$(FC) -c cubspl.f

$(OBJS): $(INCL)


.PHONY : clean
clean:
	rm -f $(OBJS) $(EXEC)

