import numpy as np

""" This file creates a grid of stars for the HRD plot """

metallicities = [0.0001, 0.001, 0.01, 0.01416]

# make a small grid for all metallicities
with open("grid.txt", "w") as f:
    masses = np.round(np.logspace(np.log10(0.1), np.log10(150.0), 500), 3)

    grid_lines = [f"--initial-mass {masses[i]} --metallicity {metallicities[j]} \n" for j in range(len(metallicities)) for i in range(len(masses))]
    f.writelines(grid_lines)

# make a small grid just for solar
with open("rapid_grid.txt", "w") as f:
    masses = np.round(np.logspace(np.log10(0.1), np.log10(150.0), 500), 3)

    grid_lines = [f"--initial-mass {masses[i]} --metallicity {0.0001} \n" for i in range(len(masses))]
    f.writelines(grid_lines)

# make a dense grid of solar (mostly dense for NSs and low mass BHs)
with open("MM20_grid.txt", "w") as f:
    low_masses = np.round(np.logspace(np.log10(0.1), np.log10(8.0), 200), 4)
    med_masses = np.round(np.logspace(np.log10(8.0), np.log10(50.0), 4600), 4)
    high_masses = np.round(np.logspace(np.log10(50.0), np.log10(150.0), 200), 4)
    masses = np.concatenate((low_masses, med_masses, high_masses))

    grid_lines = [f"--initial-mass {masses[i]} --metallicity {0.01416} \n" for i in range(len(masses))]
    f.writelines(grid_lines)
