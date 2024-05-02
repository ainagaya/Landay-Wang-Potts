# Fortran compiler
FC=gfortran

#CFLAGS = -g -fcheck=all -Wall -Wextra -pedantic -std=f2008 -fbacktrace -ffpe-trap=zero,overflow,underflow -finit-real=nan -finit-integer=-9999

all: potts_LW.o thermodynamics.o

potts_LW.o: src/main.f90
	$(FC) $(CFLAGS) $^ src/r1279.f90 src/ran2.f -o $@

thermodynamics.o: src/thermodynamics.f90
	$(FC) $(CFLAGS) $^ -o $@

run: potts_LW.o
	./$^

thermo: thermodynamics.o
	./$^
	
.PHONY:	plot1
plot1: L_40_q_10/ln_n_density_q10_L40norm.dat
	python3 plot_ln_density.py $^
	
.PHONY:	plot2
plot2: L_40_q_10/res_ising_q10_L40.dat
	python3 plot_thermodynamics.py $^


.PHONY.: clean
clean:
	@rm -f *.o *.mod MD.exe
	@echo "Modules, objects and executable removed!" 

.PHONY.: cleandata
cleandata:
	@rm -f *.dat
	@echo "Data files removed!"
