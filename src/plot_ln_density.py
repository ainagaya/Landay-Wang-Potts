import argparse
import f90nml

import matplotlib.pyplot as plt

# Parse command line arguments
namelist_file="LWparams.nml"

# Read the namelist file
namelist = f90nml.read(namelist_file)

# Get the filename parameter from the namelist
q = namelist['LWparams']['q']
L = namelist['LWparams']['L']

filename="ln_n_density_q" + str(q) + "_L" + str(L) + "norm.dat"

# Read the data from the file
data = []
with open(filename, 'r') as file:
    for line in file:
        energy, density= line.strip().split()
        if float(density) != 0:
            data.append((float(energy), float(density)))

# Separate the energy and density values
energies = [d[0] for d in data]
densities = [d[1] for d in data]  # Filter non-zero densities

# Plot the data with points
plt.plot(energies, densities, 'o', markersize=2)
plt.xlabel('E/N')
plt.ylabel('ln_g_E')
plt.title('Energy vs Density (Non-zero values)')
plt.savefig(filename + ".png")
plt.clf()
