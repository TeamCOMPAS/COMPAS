import numpy as np

""" This file creates a grid of stars for the HRD plot """

# use solar and low metallicity
metallicities = [0.0142, 0.001]

# circular binaries with a separation so wide that they are effectively single
eccentricity = 0.0
a = 1e20

with open("grid.txt", "w") as f:
    masses = np.round(np.logspace(np.log10(0.5), np.log10(150.0), 500), 3)

    grid_lines = ["--initial-mass-1 {} --initial-mass-2 {} --metallicity {} --eccentricity {} --semi-major-axis {} \n".format(masses[i], 0.1, metallicities[j], eccentricity, a) for j in range(len(metallicities)) for i in range(len(masses))]
    f.writelines(grid_lines)
