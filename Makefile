# Fortran compiler
FC=gfortran

all: potts_LW.o thermodynamics.o

potts_LW.o: main.f90
	$(FC) $^ r1279.f90 ran2.f -o $@

thermodynamics.o: thermodynamics.f90
	$(FC) $^ -o $@

run: potts_LW.o
	./$^

thermo: thermodynamics.o
	./$^
	
.PHONY:	plot1
plot1: Results_Ising/ln_n_density_q_2_L20_niter_1000000norm.dat
	python3 plot_ln_density.py $^
	
.PHONY:	plot2
plot2: Results_Ising/res_ising_q_2_L20_niter_1000000_beta.dat
	python3 plot_thermodynamics.py $^
