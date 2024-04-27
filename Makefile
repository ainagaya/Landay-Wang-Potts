# Fortran compiler
FC=gfortran

CFLAGS = -g -fcheck=all -Wall -Wextra -pedantic -std=f2008 -fbacktrace -ffpe-trap=zero,overflow,underflow -finit-real=nan -finit-integer=-9999

all: potts_LW.o thermodynamics.o

potts_LW.o: main.f90
	$(FC) $(CFLAGS) $^ r1279.f90 ran2.f -o $@

thermodynamics.o: thermodynamics.f90
	$(FC) $(CFLAGS) $^ -o $@

run: potts_LW.o
	./$^

thermo: thermodynamics.o
	./$^
	
.PHONY:	plot1
plot1: ln_n_density_q_2_L10.dat
	python3 plot_ln_density.py $^
	
.PHONY:	plot2
plot2: Results_Ising/res_ising_q_2_L20_niter_1000000_beta.dat
	python3 plot_thermodynamics.py $^


.PHONY.: clean
clean:
	@rm -f *.o *.mod MD.exe
	@echo "Modules, objects and executable removed!" 

.PHONY.: cleandata
cleandata:
	@rm -f *.dat
	@echo "Data files removed!"