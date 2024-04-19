
# Fortran compiler
FC=gfortran

all: potts_LW.o
potts_LW.o: main.f90
	$(FC) $^ r1279.f90 ran2.f  -o $@
