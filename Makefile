#Change this to where you installed GadgetReader
GREAD=${CURDIR}/../GadgetReader
#Use this for extra library directories, eg for FFTW or the GSL.
LIBDIR =

#OPT   += -DNO64BITID # switch this on if you want normal 32-bit IDs
OPT   +=  -DCORRECT_CIC  # only switch this on if particles are homogenously distributed over mesh cells (say glass)
OPT   +=  -DNEUTRINOS  # this will make type 2 be neutrinos instead of a second DM component
#OPT   +=  -DNEUTRINO_PAIRS  # this will produce an additional partner for every neutrino with opposite thermal velocities
#OPT   += -DPRINT_SPEC #Use this to print out the spectrum (with non-Gaussianity) after calculating ICs.
HDF_LIB = -lhdf5 -lhdf5_hl

LFLAGS += $(LIBDIR) -lfftw3f_threads -lfftw3f -lfftw3_threads -lfftw3 -lgsl -lgslcblas -lpthread -lrgad ${HDF_LIB} -lwgad -L${GREAD} -Wl,-rpath,$(GREAD)
CFLAGS += -I${GREAD} ${OPT}
#PRO = -fprofile-generate
#PRO = -fprofile-use -fprofile-correction
#Are we using gcc or icc?
ifeq (icc,$(findstring icc,${CC}))
  CFLAGS +=-O2 -g -c -w1 -openmp
  LINK +=${CXX} -openmp
else
  CFLAGS +=-O2 -ffast-math -g -c -Wall -fopenmp $(PRO)
  LINK +=${CXX} $(PRO)
  LFLAGS += -lm -lgomp
endif

CXXFLAGS +=${CFLAGS} -std=gnu++11
EXEC   = N-GenIC

OBJS   = power.o allvars.o save.o read_param.o \
	 initialise.o print_spec.o thermalvel.o cosmology.o displacement.o

INCL   = allvars.h proto.h part_data.hpp thermalvel.hpp power.hpp read_param.hpp displacement.hpp Makefile

.PHONY : clean all test

all: $(EXEC)

createnu: createnu.o $(OBJS)
	${LINK}${LFLAGS} -lhdf5_hl -lhdf5 $^ -o $@

test: btest
	./$^

doc: Doxyfile main.cpp ${INCL}
	doxygen $<

%.o: %.cpp $(INCL)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $< -o $@

allvars.o: allvars.c allvars.h Makefile

$(EXEC): main.o $(OBJS)
	${LINK} ${LFLAGS} $^ -o  $@  

btest: test.o ${OBJS}
	${LINK} ${LFLAGS} -lboost_unit_test_framework $^ -o  $@
clean:
	rm -f $(OBJS) $(EXEC)

