#Change this to where you installed GadgetReader
GREAD=${CURDIR}/../GadgetReader
#Use this for extra library directories, eg for FFTW or the GSL.
LIBDIR =

#OPT   += -DNO64BITID # switch this on if you want normal 32-bit IDs
#OPT   +=  -DCORRECT_CIC  # only switch this on if particles are homogenously distributed over mesh cells (say glass)
OPT   +=  -DNEUTRINOS  # this will make type 2 be neutrinos instead of a second DM component
#OPT   +=  -DNEUTRINO_PAIRS  # this will produce an additional partner for every neutrino with opposite thermal velocities
#OPT   += -DPRINT_SPEC #Use this to print out the spectrum (with non-Gaussianity) after calculating ICs.

ifeq ($(CC),cc)
  ICC:=$(shell which icc --tty-only 2>&1)
  #Can we find icc?
  ifeq (/icc,$(findstring /icc,${ICC}))
     CC = icc -vec_report0
     CXX = icpc
  else
     GCC:=$(shell which gcc --tty-only 2>&1)
     #Can we find gcc?
     ifeq (/gcc,$(findstring /gcc,${GCC}))
        CC = gcc
        CXX = g++
     endif
  endif
endif

LFLAGS += $(LIBDIR) -lfftw3f_threads -lfftw3f -lgsl -lgslcblas -lpthread -lrgad -lwgad -L${GREAD} -Wl,-rpath,$(GREAD)
CFLAGS += -I${GREAD} ${OPT}
#PRO = -pg
#Are we using gcc or icc?
ifeq (icc,$(findstring icc,${CC}))
  CFLAGS +=-O2 -g -c -w1 -openmp
  LINK +=${CXX} -openmp
else
  CFLAGS +=-O2 -g -c -Wall -fopenmp
  LINK +=${CXX} -openmp $(PRO)
  LFLAGS += -lm -lgomp
endif

CXXFLAGS +=${CFLAGS}
EXEC   = N-GenIC

OBJS   = main.o power.o allvars.o save.o read_param.o  read_glass.o  \
	 initialise.o print_spec.o \
	nrsrc/qromb.o nrsrc/nrutil.o nrsrc/polint.o nrsrc/trapzd.o 

INCL   = allvars.h proto.h Makefile

.PHONY : clean all


all: $(EXEC)

nrsrc/%.o: nrsrc/%.c nrsrc/nrutil.h 
%.o: %.cpp $(INCL)

allvars.o: allvars.c allvars.h Makefile
read_param.o: read_param.c $(INCL)

$(EXEC): $(OBJS) 
	${LINK} ${LFLAGS} $^ -o  $@  

clean:
	rm -f $(OBJS) $(EXEC)

