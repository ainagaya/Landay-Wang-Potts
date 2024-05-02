import argparse

import matplotlib.pyplot as plt
# Parse command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('filename', help='Input file name')
args = parser.parse_args()

# Read the data from the file
data = []
with open(args.filename, 'r') as file:
    for line in file:
        T, density, energy= line.strip().split()
        if float(density) != 0:
            data.append((float(T), float(density), float(energy)))

# Separate the data by beta values
T_data = {}
for T, energy_density, energy in data:
    if T not in T_data:
        T_data[T] = []
    T_data[T].append((energy_density, energy))

# Plot the data for each beta
for T, data_points in T_data.items():
    energy_densities = [float(d[0]) for d in data_points]
    energies = [float(d[1]) for d in data_points]
    
    # Plot the data with smaller points
    plt.plot([energy for energy in energies], energy_densities, 'o', label='T = {T}', markersize=1)
  #  plt.xlim([0, 8])
    plt.xlabel('Energy/N')
    plt.ylabel('g(e)exp(-betaE)')
    plt.title(f'g(e)exp(-betaE) vs Energy/N (T = {T})')
    plt.legend()
    plt.savefig(args.filename + f"energy_dens{T}.png")
    plt.clf()
#    plt.show()