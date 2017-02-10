#Change this to where you installed GadgetReader
GREAD=${CURDIR}/GadgetReader
BIGFILE=${GREAD}/bigfile/src/
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
OPT += -DHAVE_BGFL  #Use this if you have bigfile support compiled in
	#Check for a pkgconfig; if one exists we are probably debian.
ifeq ($(shell pkg-config --exists hdf5 && echo 1),1)
	HDF_LIB = $(shell pkg-config --libs hdf5) -lhdf5_hl
	HDF_INC = $(shell pkg-config --cflags hdf5)
else
	HDF_LIB = -lhdf5 -lhdf5_hl
	HDF_INC =
endif

LFLAGS += $(LIBDIR) -lfftw3f_threads -lfftw3f -lfftw3_threads -lfftw3 -lgsl -lgslcblas -lpthread ${HDF_LIB} -lwgad -L${GREAD} -Wl,-rpath,$(GREAD),--no-add-needed,--as-needed -L${BIGFILE} -lbigfile
CFLAGS += -I${GREAD} -I${BIGFILE} ${OPT} ${HDF_INC}
#PRO = -fprofile-generate
#PRO = -fprofile-use -fprofile-correction
#Are we using gcc or icc?
ifeq (icc,$(findstring icc,${CC}))
  CFLAGS +=-O2 -g -c -w1 -openmp
  LINK = $(CXX) -openmp
else
  CFLAGS +=-O2 -ffast-math -g -c -Wall -fopenmp $(PRO)
  LINK = $(CXX) $(PRO)
  LFLAGS += -lm -lgomp
endif

CXXFLAGS +=${CFLAGS} -std=gnu++11
EXEC   = N-GenIC

OBJS   = power.o save.o read_param.o \
	 thermalvel.o cosmology.o displacement.o part_data.o gsl_spline_wrapper.o

INCL   = save.hpp part_data.hpp thermalvel.hpp power.hpp read_param.hpp displacement.hpp Makefile

.PHONY : clean all test

all: $(EXEC)

test: btest
	./$^

doc: Doxyfile main.cpp ${INCL}
	doxygen $<

$(BIGFILE)/bigfile.h:
	./boostrap.sh

$(GREAD)/libwgad.so: $(BIGFILE)/bigfile.h $(GREAD)/gadgetwriter.hpp
	cd $(GREAD); VPATH=$(GREAD) make $(@F)
%.o: %.cpp %.hpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $< -o $@

main.o: main.cpp ${INCL}
	$(LINK) $(CXXFLAGS) $(CPPFLAGS) $< -o $@

test.o: test.cpp ${INCL}
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $< -o $@

displacement.o: displacement.cpp displacement.hpp part_data.hpp power.hpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $< -o $@

cosmology.o: cosmology.cpp cosmology.hpp physconst.h
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $< -o $@

save.o: save.cpp save.hpp part_data.hpp physconst.h thermalvel.hpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $< -o $@

thermalvel.o: thermalvel.cpp thermalvel.hpp gsl_spline_wrapper.hpp physconst.h
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $< -o $@

$(EXEC): main.o $(OBJS) $(GREAD)/libwgad.so
	${LINK} $^ ${LFLAGS} -o  $@

btest: test.o ${OBJS}
	${LINK} $^ ${LFLAGS} -lboost_unit_test_framework -o  $@
clean:
	rm -f $(OBJS) $(EXEC)

