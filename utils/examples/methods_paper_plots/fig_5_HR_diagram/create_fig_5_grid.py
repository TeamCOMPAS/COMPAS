""" This file creates a grid of stars for the HRD plot """

# use solar and low metallicity
metallicities = [0.0142, 0.001]

# circular binaries with a separation so wide that they are effectively single
eccentricity = 0.0
a = 1e20

with open("grid.txt", "w") as f:
    # mass range: basically logspace from 0.5 to 50 but with cleaner numbers
    masses = [0.5, 0.65, 0.80, 1.0, 1.3, 1.6, 2.0, 2.5, 3.2, 4.0, 5.2, 6.5, 8.2, 10.3, 13.0, 16.5, 21, 26, 33, 45, 65, 90, 115, 150]

    grid_lines = ["--initial-mass-1 {} --initial-mass-2 {} --metallicity {} --eccentricity {} --semi-major-axis {} \n".format(masses[i], 0.1, metallicities[j], eccentricity, a) for j in range(len(metallicities)) for i in range(len(masses))]
    f.writelines(grid_lines)
