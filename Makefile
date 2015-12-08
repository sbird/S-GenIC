#Change this to where you installed GadgetReader
GREAD=${CURDIR}/../GadgetReader
#Use this for extra library directories, eg for FFTW or the GSL.
LIBDIR =

#OPT   += -DNO64BITID # switch this on if you want normal 32-bit IDs
# This corrects for the CIC window function of the particles.
# It is inaccurate for a regular grid of particles because a regular grid has *extra* power at the grid scale,
# which gets convolved with the power spectrum and over-corrects.
# However, you might nevertheless want it on for a small fourier mesh, as then the loss of power becomes important.
# The most accurate option is to use a regular grid, Nmesh >= 2 * Npart and turn this off.
#OPT   +=  -DCORRECT_CIC  
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

OBJS   = power.o save.o read_param.o \
	 print_spec.o thermalvel.o cosmology.o displacement.o part_data.o gsl_spline_wrapper.o

INCL   = save.hpp part_data.hpp thermalvel.hpp power.hpp read_param.hpp displacement.hpp Makefile

.PHONY : clean all test

all: $(EXEC)

test: btest
	./$^

doc: Doxyfile main.cpp ${INCL}
	doxygen $<

%.o: %.cpp %.hpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $< -o $@

main.o: main.cpp ${INCL}
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $< -o $@

test.o: test.cpp ${INCL}
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $< -o $@

displacement.o: displacement.cpp displacement.hpp part_data.hpp power.hpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $< -o $@

cosmology.o: cosmology.cpp cosmology.hpp physconst.h
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $< -o $@

print_spec.o: print_spec.cpp cosmology.hpp power.hpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $< -o $@

save.o: save.cpp save.hpp part_data.hpp physconst.h thermalvel.hpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $< -o $@

thermalvel.o: thermalvel.cpp thermalvel.hpp gsl_spline_wrapper.hpp physconst.h
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $< -o $@

$(EXEC): main.o $(OBJS)
	${LINK} ${LFLAGS} $^ -o  $@

btest: test.o ${OBJS}
	${LINK} ${LFLAGS} -lboost_unit_test_framework $^ -o  $@
clean:
	rm -f $(OBJS) $(EXEC)

