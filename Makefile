SYSTYPE="#SYSTYPE#"

EXEC   = N-GenIC

OBJS   = main.o power.o allvars.o save.o read_param.o  read_glass.o  \
         nrsrc/nrutil.o nrsrc/qromb.o nrsrc/polint.o nrsrc/trapzd.o #cubspl.o

INCL   = allvars.h proto.h  nrsrc/nrutil.h  Makefile

OPT   += -DDIFFERENT_TRANSFER_FUNC  # set this if you want to implement a transfer function that depends on particle type. Or for tk_CAMB to work.
OPT	+= -DFORMAT_TWO  #Set this if you want to output IC files in format 2.												

#OPT   += -DNO64BITID # switch this on if you want normal 32-bit IDs
#OPT   +=  -DCORRECT_CIC  # only switch this on if particles are homogenously distributed over mesh cells (say glass)

OPT   +=  -DNEUTRINOS  # this will produce a second component as slight neutrinos (needs to be in initial glass)
#OPT   +=  -DNEUTRINO_PAIRS  # this will produce an additional partner for every neutrino with opposite thermal velocities

OPTIONS = $(OPT)

CC       =  gcc  # sets the C-compiler (default)
FC			=  gfortran
OPTIMIZE =  -O2 -g -fopenmp -Wall
MPICHLIB = 
FFTW_INCL=  
FFTW_LIBS=  

FFTW_LIB =  $(FFTW_LIBS) -lfftw3f_threads -lfftw3f -lpthread -lgomp 
LIBS   =   -lm   $(FFTW_LIB)  $(GSL_LIBS)  -lgsl -lgslcblas -L../GadgetReader -lrgad -lwgad

ifeq ($(SYSTYPE),"Solaris")
LIBS   =   -R/opt/local/lib/sparcv9 -lm  -lmpi   $(GSL_LIBS) -lgsl -lgslcblas  $(FFTW_LIB)
endif

CFLAGS = $(OPTIONS) $(OPTIMIZE) $(FFTW_INCL) $(GSL_INCL) -I../GadgetReader

CXXFLAGS = -I../GadgetReader $(OPTIMIZE) $(OPTIONS)


$(EXEC): $(OBJS) 
	$(CXX) $(OPTIMIZE) $(OBJS) $(LIBS)   -o  $(EXEC)  

cubspl.o: cubspl.f
	$(FC) -c cubspl.f

$(OBJS): $(INCL)


.PHONY : clean
clean:
	rm -f $(OBJS) $(EXEC)

