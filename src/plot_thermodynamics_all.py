import argparse
import f90nml

import matplotlib.pyplot as plt

plt.style.use('ggplot')

# Parse command line arguments
namelist_file="LWparams.nml"

# Read the namelist file
namelist = f90nml.read(namelist_file)

# Get the filename parameter from the namelist
q = namelist['LWparams']['q']
#L = namelist['LWparams']['L']

# INTERNAL ENERGY
for L in [10, 20, 30, 40]:
    filename="res_q" + str(q) + "_L" + str(L) + ".dat"

    # Read the data from the file
    data = []
    with open(filename, 'r') as file:
        for line in file:
            beta, internal_energy, free_energy, entropy, specific_heat = line.strip().split()
            data.append((beta, internal_energy, free_energy, entropy, specific_heat))


    # Separate the energy and density values
    betas = [float(d[0]) for d in data]
    print("")
    internal_energies = [float(d[1]) for d in data]
    print("")

    ################ 1/BETA 
    # Plot the data with smaller points
    plt.plot([1/beta for beta in betas], internal_energies, '-', label=f'L={L}', markersize=1)
    plt.xlim([0.66, 0.72])
    plt.xlabel('T')
    plt.ylabel('U(L,T)/N')

plt.legend()
plt.savefig("all_thermo" + "internal_energy.png")

plt.clf()

# FREE ENERGY
for L in [10, 20, 30, 40]:
    filename="res_q" + str(q) + "_L" + str(L) + ".dat"

    # Read the data from the file
    data = []
    with open(filename, 'r') as file:
        for line in file:
            beta, internal_energy, free_energy, entropy, specific_heat = line.strip().split()
            data.append((beta, internal_energy, free_energy, entropy, specific_heat))

    # Separate the energy and density values
    betas = [float(d[0]) for d in data]
    print("")
    free_energies = [float(d[2]) for d in data]
    print("")

    ################ 1/BETA 
    # Plot the data with smaller points
    plt.plot([1/beta for beta in betas], free_energies, '-', label=f'L={L}', markersize=1)
    plt.xlim([0.6, 0.8])
    plt.ylim([-2.2, -1.9])
    plt.xlabel('T')
    plt.ylabel('F(L,T)/N')

plt.legend()
plt.savefig("all_thermo" + "free_energy.png")

plt.clf()

# Entropy
for L in [10, 20, 30, 40]:
    filename="res_q" + str(q) + "_L" + str(L) + ".dat"

    # Read the data from the file
    data = []
    with open(filename, 'r') as file:
        for line in file:
            beta, internal_energy, free_energy, entropy, specific_heat = line.strip().split()
            data.append((beta, internal_energy, free_energy, entropy, specific_heat))

    # Separate the energy and density values
    betas = [float(d[0]) for d in data]
    entropies = [float(d[3]) for d in data]

    ################ 1/BETA 
    # Plot the data with smaller points
    plt.plot([1/beta for beta in betas], entropies, '-', label=f'L={L}', markersize=1)
    plt.xlim([0.6, 0.8])
    plt.xlabel('T')
    plt.ylabel('S(L,T)/N')

plt.legend()
plt.savefig("all_thermo" + "entropy.png")

plt.clf()

# speciFIC HEAT
for L in [10, 20, 30, 40]:
    filename="res_q" + str(q) + "_L" + str(L) + ".dat"

    # Read the data from the file
    data = []
    with open(filename, 'r') as file:
        for line in file:
            beta, internal_energy, free_energy, entropy, specific_heat = line.strip().split()
            data.append((beta, internal_energy, free_energy, entropy, specific_heat))

    # Separate the energy and density values
    betas = [float(d[0]) for d in data]

    specific_heats = [float(d[4]) for d in data]

    ################ 1/BETA 
    # Plot the data with smaller points
    plt.plot([1/beta for beta in betas], specific_heats, '-', label=f'L={L}', markersize=1)
    plt.xlim([0.7, 0.74])
    plt.xlabel('T')
    plt.ylabel('C(L,T)/N')
#    plt.xticks(ticks=[0.701, 0.702, 0.703, 0.704])

plt.legend()
plt.savefig("all_thermo" + "spec_heat.png")

plt.clf()



