# compiler
FC = gfortran -fopenmp -fdefault-integer-8#-Wall -fcheck=all -pg

# compile flags
FCFLAGS = -c
# link flags
FLFLAGS =
LIB = -lblas -llapack -lstdc++ -lrt -ldl  -lz
# debug flags
DEBUGFLAGS = -g -Og -ffpe-trap=invalid,zero,overflow,underflow -fbacktrace

PROGS=CI_gen # ao2mo lapack_diag iHj iHj_2 prun randSCI

all : CI_gen clean
#randSCI ao2mo iHj iHj_2 prun clean

debug : CI_gen_debug

iomod.o : iomod.f90
	$(FC) $(FCFLAGS) $<

readin.o : readin.f90
	$(FC) $(FCFLAGS) $<

tools.o : tools.f90
	$(FC) $(FCFLAGS) $<

var.o : var.f90
	$(FC) $(FCFLAGS) $<

rbm.o : rbm.f90 tools.o 
	$(FC) $(FCFLAGS) $<

CI_gen : CI_gen.f90 iomod.o readin.o tools.o var.o rbm.o 
	$(FC) $^ $(LIB) -o $@

CI_gen_debug : CI_gen.f90 iomod.o readin.o tools.o var.o rbm.o 
	$(FC) $(DEBUGFLAGS) $^ $(LIB) -o CI_gen

clean : 
	rm  *.o *.mod

install: $(PROG)
	mv $(PROGS) ../bin/.
