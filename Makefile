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
plot1: 
	python3 src/plot_ln_density.py 
	
.PHONY:	plot2
plot2: 
	python3 src/plot_thermodynamics.py

.PHONY:	plot3
plot3: 
	python3 src/plot_prob_energy_density.py


.PHONY.: clean
clean:
	@rm -f *.o *.mod MD.exe
	@echo "Modules, objects and executable removed!" 

.PHONY.: cleandata
cleandata:
	@rm -f *.dat
	@echo "Data files removed!"
