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
        energy, density, caca = line.strip().split()
        if float(density) != 0:
            data.append((float(energy), float(density)))

# Separate the energy and density values
energies = [d[0] for d in data]
print(energies)
densities = [d[1] for d in data]  # Filter non-zero densities
print(densities)

# Plot the data with points
plt.plot(energies, densities, 'o')
plt.xlabel('E/N')
plt.ylabel('ln_g_E')
plt.title('Energy vs Density (Non-zero values)')
plt.savefig(args.filename + ".png")
plt.show()
