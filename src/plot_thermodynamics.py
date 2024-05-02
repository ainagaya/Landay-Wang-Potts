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
        print(line)
        beta, internal_energy, free_energy, entropy, specific_heat = line.strip().split()
        data.append((beta, internal_energy, free_energy, entropy, specific_heat))


# Read the data from the ferdi file
data_control = []
with open("ferdi.D", 'r') as file:
    for line in file:
        print(line)
        beta, F, ueng, act, S, cht, zl = line.strip().split()
        data_control.append((beta, F, ueng, act, S, cht, zl))


betas_control = [float(d[0]) for d in data_control]
F_control = [float(d[1]) for d in data_control]
ueng_control = [(float(d[2])-2)/2 for d in data_control]
act_control = [float(d[3]) for d in data_control]
S_control = [float(d[4]) for d in data_control]
cht_control = [float(d[5]) for d in data_control]
zl_control = [float(d[6]) for d in data_control]


# Separate the energy and density values
betas = [float(d[0]) for d in data]
print(betas)
print("")
internal_energies = [float(d[1]) for d in data]
print(internal_energies)
print("")
free_energies = [float(d[2]) for d in data]
print(free_energies)
print("")
entropies = [float(d[3]) for d in data]
print(entropies)
print("")
specific_heats = [float(d[4]) for d in data]
print(specific_heats)
print("")




################ 1/BETA 



# Plot the data with smaller points
plt.plot([1/beta for beta in betas], internal_energies, 'o', label='Data', markersize=1)
plt.plot([0.5*1/beta for beta in betas_control], ueng_control, '.', label='Control', markersize=1)
plt.xlim([0, 8])
plt.xlabel('1/beta')
plt.ylabel('internal energy')
plt.title('Internal Energy vs Beta')
plt.legend()
plt.savefig(args.filename + "int_energy.png")
plt.show()

# Plot the data with points
plt.plot([1/beta for beta in betas], free_energies, 'o', label='Data', markersize=1)
plt.plot([0.5*1/beta for beta in betas_control], F_control, '.', label='Control', markersize=1)
plt.xlim([0, 8])
plt.xlabel('1/beta')
plt.ylabel('free energy')
plt.title('Free Energy vs Beta')
plt.legend()
plt.savefig(args.filename + "free_energy.png")
plt.show()

# Plot the data with points
plt.plot([1/beta for beta in betas], entropies, 'o', label='Data', markersize=1)
plt.plot([0.5*1/beta for beta in betas_control], S_control, '.', label='Control', markersize=1)
plt.xlim([0, 8])
plt.xlabel('1/beta')
plt.ylabel('entropy')
plt.title('Entropy vs Beta')
plt.legend()
plt.savefig(args.filename + "entropy.png")
plt.show()

# Plot the data with points
plt.plot([1/beta for beta in betas], specific_heats, 'o', label='Data', markersize=1)
plt.plot([0.5*1/beta for beta in betas_control], cht_control, '.', label='Control', markersize=1)
plt.xlim([0, 8])
plt.xlabel('1/beta')
plt.ylabel('Specific Heat')
plt.title('Specific heat vs Beta')
plt.legend()
plt.savefig(args.filename + "spec_heat.png")
plt.show()

################### BETA

# Plot the data with smaller points
plt.plot([beta for beta in betas], internal_energies, 'o', label='Data', markersize=1)
plt.plot([2*beta for beta in betas_control], ueng_control, '.', label='Ising',  markersize=1)
#plt.xlim([0, 8])
plt.xlabel('beta')
plt.ylabel('internal energy')
plt.title('Internal Energy vs Beta')
plt.legend()
plt.savefig(args.filename + "int_energy_beta.png")
plt.show()

# Plot the data with points
plt.plot([beta for beta in betas], free_energies, 'o', label='Data', markersize=1)
plt.plot([2*beta for beta in betas_control], F_control, '.', label='Ising', markersize=1)
#plt.xlim([0, 8])
plt.xlabel('beta')
plt.ylabel('free energy')
plt.title('Free Energy vs Beta')
plt.legend()
plt.savefig(args.filename + "free_energy_beta.png")
plt.show()

# Plot the data with points
plt.plot([beta for beta in betas], entropies, 'o', label='Data', markersize=1)
plt.plot([2*beta for beta in betas_control], S_control, '.', label='Ising', markersize=1)
#plt.xlim([0, 8])
plt.xlabel('beta')
plt.ylabel('entropy')
plt.title('Entropy vs Beta')
plt.legend()
plt.savefig(args.filename + "entropy_beta.png")
plt.show()

# Plot the data with points
plt.plot([beta for beta in betas], specific_heats, 'o', label='Data', markersize=1)
plt.plot([2*beta for beta in betas_control], cht_control, '.', label='Ising', markersize=1)
#plt.xlim([0, 8])
plt.xlabel('beta')
plt.ylabel('Specific Heat')
plt.title('Specific heat vs Beta')
plt.legend()
plt.savefig(args.filename + "spec_heat_beta.png")
plt.show()