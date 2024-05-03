import f90nml

import matplotlib.pyplot as plt

#plt.style.use('ggplot')

# Parse command line arguments
namelist_file="LWparams.nml"

# Read the namelist file
namelist = f90nml.read(namelist_file)

# Get the filename parameter from the namelist
q = namelist['LWparams']['q']
L = namelist['LWparams']['L']

filename="energy_density_q" + str(q) + "_L" + str(L) + ".dat"

# Read the data from the file
data = []
with open(filename, 'r') as file:
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
    plt.plot([energy for energy in energies], energy_densities, '-', label=f'T = {T:.2f}', linewidth=0.7)
  #  plt.xlim([0, 8])


    
#    plt.clf()

plt.xlabel('E/N')
plt.ylabel('n(E)exp(-betaE)')
#plt.title(f'g(e)exp(-betaE) vs Energy/N (T = {T:.5f})')
plt.legend()
plt.savefig(filename + "_peaks.png")