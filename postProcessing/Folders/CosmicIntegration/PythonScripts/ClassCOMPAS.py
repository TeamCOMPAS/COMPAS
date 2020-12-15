#!/usr/bin/env python3
import numpy as np
import h5py as h5
import os
import astropy.units as u
import astropy.constants as c

class COMPASData(object):
    """ A Class for reading and masking a COMPAS data file """

    def __init__(self, path, filename="COMPASOutput.h5", m1_min=5 * u.Msun, m1_max=150 * u.Msun, m2_min=0.1 * u.Msun, fbin=0.7):
        """ 
            Initialise the Class to prepare for other functions 
        
            Args:
                path     --> [string] Path to the COMPAS file that contains the output
                filename --> [string] Name of the COMPAS file
                m1_min   --> [float]  Minimum primary mass sampled by COMPAS
                m1_max   --> [float]  Maximum primary mass sampled by COMPAS
                m2_min   --> [float]  Minimum secondary mass sampled by COMPAS
                fbin     --> [float]  Binary fraction used by COMPAS
        """
        # file details for getting data
        self.compas_file = path + filename

        # mass cuts and binary fraction for total star forming mass
        self.m1_min = m1_min
        self.m1_max = m1_max
        self.m2_min = m2_min
        self.fbin = fbin
        self.n_systems = None
        self.mass_evolved_per_binary = None

        self.metallicities_allsystems = None

        # parameters for all DCOs
        self.metallicities = None
        self.delay_times = None
        self.primary_masses = None
        self.secondary_masses = None
        self.sw_weights = None

        # masks that give different DCOs
        self.DCO_masks = {
            "ALL": None,
            "BHBH": None,
            "BHNS": None,
            "NSNS": None
        }

    def get_COMPAS_variables(self, hdf5_file, var_names):
        """ 
            Get a variable or variables from a COMPAS file

            Args:
                hdf5_file --> [string]                Name of HDF5 subfile (e.g. "DoubleCompactObjects")
                var_names --> [string or string list] A variable name or list of variables names to return

            Returns:
                var_list  --> [list of lists]         A list of variables (or a single variable if only one name supplied)
        """
        # open the COMPAS file
        with h5.File(self.compas_file, "r") as compas_file:
            # if the list is only a string (i.e. one variable) then don't return a list
            if isinstance(var_names, str):
                return compas_file[hdf5_file][var_names][...].squeeze()
            # else return each variable in a list
            else:
                return [compas_file[hdf5_file][var_name][...].squeeze() for var_name in var_names]

    def set_DCO_masks(self, dco_types=["ALL", "BHBH", "BHNS", "NSNS"], merges_in_hubble_time=True, no_RLOF_after_CEE=True, pessimistic_CEE=True):
        """
            Create masks for different DCO types

            Args:
                dco_types             --> which DCO types to create masks for (default is all of them)
                merges_in_hubble_time --> whether to mask binaries that don't merge in a Hubble time
                no_RLOF_after_CEE     --> whether to mask binaries that have immediate RLOF after a CCE
                pessimistic_CEE       --> whether to mask binaries that go through Optimistic CE scenario
        """
        # get the appropriate variables from the COMPAS file
        stellar_type_1, stellar_type_2, hubble_flag, dco_seeds = \
            self.get_COMPAS_variables("DoubleCompactObjects", ["Stellar_Type(1)", "Stellar_Type(2)", "Merges_Hubble_Time", "SEED"])

        # if user wants to mask on Hubble time use the flag, otherwise just set all to True
        hubble_mask = hubble_flag if merges_in_hubble_time else True

        # mask on stellar types (where 14=BH and 13=NS), BHNS can be BHNS or NSBH
        type_masks = {
            "ALL": True,
            "BHBH": np.logical_and(stellar_type_1 == 14, stellar_type_2 == 14),
            "BHNS": np.logical_or(np.logical_and(stellar_type_1 == 14, stellar_type_2 == 13), np.logical_and(stellar_type_1 == 13, stellar_type_2 == 14)),
            "NSNS": np.logical_and(stellar_type_1 == 13, stellar_type_2 == 13),
        }

        # if the user wants to make RLOF or optimistic CEs
        if no_RLOF_after_CEE or pessimistic_CEE:

            # get the flags and unique seeds from the Common Envelopes file
            ce_seeds = self.get_COMPAS_variables("CommonEnvelopes", "SEED")
            dco_from_ce = np.in1d(ce_seeds, dco_seeds)
            dco_ce_seeds = ce_seeds[dco_from_ce]

            # if masking on RLOF, get flag and match seeds to dco seeds
            if no_RLOF_after_CEE:
                rlof_flag = self.get_COMPAS_variables("CommonEnvelopes", "Immediate_RLOF>CE")[dco_from_ce]
                rlof_seeds = np.unique(dco_ce_seeds[rlof_flag])
                rlof_mask = np.logical_not(np.in1d(dco_seeds, rlof_seeds))
            else:
                rlof_mask = True

            # if masking on pessimistic CE, get flag and match seeds to dco seeds
            if pessimistic_CEE:
                pessimistic_flag = self.get_COMPAS_variables("CommonEnvelopes", "Optimistic_CE")[dco_from_ce]
                pessimistic_seeds = np.unique(dco_ce_seeds[pessimistic_flag])
                pessimistic_mask = np.logical_not(np.in1d(dco_seeds, pessimistic_seeds))
            else:
                pessimistic_mask = True
        else:
            rlof_mask = True
            pessimistic_mask = True

        # create a mask for each dco type supplied
        for dco_type in dco_types:
            self.DCO_masks[dco_type] = type_masks[dco_type] * hubble_mask * rlof_mask * pessimistic_mask

    def set_COMPAS_data(self):
        """ 
            Set the essential COMPAS data from the file
        """
        # get the primary
        primary_masses, secondary_masses, formation_times, coalescence_times = \
            self.get_COMPAS_variables("DoubleCompactObjects", ["Mass(1)", "Mass(1)", "Time", "Coalescence_Time"])
        self.primary_masses, self.secondary_masses = primary_masses * u.Msun, secondary_masses * u.Msun
        self.delay_times = (formation_times + coalescence_times) * u.Myr

        dco_seeds = self.get_COMPAS_variables("DoubleCompactObjects", "SEED")
        systems_seeds, systems_metallicites = self.get_COMPAS_variables("SystemParameters", ["SEED", "Metallicity@ZAMS(1)"])
        self.n_systems = len(systems_seeds)
        self.metallicities_allsystems = systems_metallicites

        systems_to_dco = np.in1d(systems_seeds, dco_seeds)
        self.metallicities = systems_metallicites[systems_to_dco]

        self.sw_weights = np.ones(len(self.primary_masses))

    def set_sw_weights(self, column_name):
        """ Set STROOPWAFEL adaptive sampling weights given a column name in the DoubleCompactObjects file """
        if column_name is not None:
            self.sw_weights = self.get_COMPAS_variables("DoubleCompactObjects", column_name)
        
    def chirp_masses(self):
        """ Calculate the chirp masses of the data """
        if self.primary_masses is None or self.secondary_masses is None:
            print("You need to set data first!")
            return
        return (self.primary_masses * self.secondary_masses)**(3/5) / (self.primary_masses + self.secondary_masses)**(1/5)

    def etas(self):
        """ Calculate the symmetric mass ratios of the data """
        if self.primary_masses is None or self.secondary_masses is None:
            print("You need to set data first!")
            return
        return (self.primary_masses*self.secondary_masses)/(self.primary_masses+self.secondary_masses)**2

    def find_star_forming_mass_per_binary_sampling(self, m1=0.01, m2=0.08, m3=0.5, m4=200.0, a12=0.3, a23=1.3, a34=2.3,
            primary_mass_inverse_CDF=None, mass_ratio_inverse_CDF=None, SAMPLES=20000000):
        """
            Calculate the star forming mass evolved for each binary in the file.
            This function does this by sampling from the IMF and mass ratio distributions

            Args:
                mi                       --> [float]    masses at which to transition the slope of the IMF (ignored if primary_mass_inverse_CDF is not None)
                aij                      --> [float]    slope of the IMF between mi and mj (ignored if primary_mass_inverse_CDF is not None)
                primary_mass_inverse_CDF --> [function] a function that computes the inverse CDF functoin for the primary mass distribution
                                                        this defaults to the Kroupa IMF (which can be varied using mi, aij)
                mass_ratio_inverse_CDF   --> [function] a function that computes the inverse CDF function for the mass ratio distribution
                                                        this defaults to assuming a uniform mass ratio on [0, 1]
                SAMPLES                  --> [int]      number of samples to draw when creating a mock universe
        """
        # if primary mass inverse CDF is None, assume the Kroupa IMF
        if primary_mass_inverse_CDF is None:
            primary_mass_inverse_CDF = lambda U: inverse_CDF_IMF(U, m1=m1, m2=m2, m3=m3, m4=m4, a12=a12, a23=a23, a34=a34)

        # if mass ratio inverse CDF function is None, assume uniform
        if mass_ratio_inverse_CDF is None:
            mass_ratio_inverse_CDF = lambda q: q

        # randomly sample a large number of masses from IMF, mass ratios from supplied function, binary for boolean
        primary_mass = primary_mass_inverse_CDF(np.random.rand(SAMPLES)) * u.Msun
        mass_ratio = mass_ratio_inverse_CDF(np.random.rand(SAMPLES))
        binary = np.random.rand(SAMPLES)

        # only fbin fraction of stars have a secondary (in a binary)
        binary_mask = binary < self.fbin
        
        # assign each a random secondary mass, default 0 because single stars have m2=0 (surprisingly :P)
        secondary_mass = np.zeros(SAMPLES) * u.Msun
        secondary_mass[binary_mask] = primary_mass[binary_mask] * mass_ratio[binary_mask]

        # find the total mass of the whole population
        total_mass = np.sum(primary_mass) + np.sum(secondary_mass)

        # apply the COMPAS cuts on primary and secondary mass
        primary_mask = np.logical_and(primary_mass >= self.m1_min, primary_mass <= self.m1_max)
        secondary_mask = secondary_mass > self.m2_min
        full_mask = np.logical_and(primary_mask, secondary_mask)
        
        # find the total mass with COMPAS cuts
        total_mass_COMPAS = np.sum(primary_mass[full_mask]) + np.sum(secondary_mass[full_mask])

        # use the totals to find the ratio and return the average mass as well
        f_mass_sampled = total_mass_COMPAS / total_mass
        average_mass_COMPAS = total_mass_COMPAS / len(primary_mass[full_mask])
        
        # find the average star forming mass evolved per binary in the Universe
        self.mass_evolved_per_binary = average_mass_COMPAS / f_mass_sampled

# ============================================== #
# Initial Mass Function PDF, CDF and inverse CDF #
# ============================================== #

def IMF(m, m1=0.01, m2=0.08, m3=0.5, m4=200.0, a12=0.3, a23=1.3, a34=2.3):
    """ 
        Calculate the fraction of stellar mass between m and m + dm for a three part broken power law.
        Default values follow Kroupa (2001)
            zeta(m) ~ m^(-a_ij)
        
        Args:
            m       --> [float, list of floats] mass or masses at which to evaluate
            mi      --> [float]                 masses at which to transition the slope
            aij     --> [float]                 slope of the IMF between mi and mj
            
        Returns:
            zeta(m) --> [float, list of floats] value or values of the IMF at m
    """
    # calculate normalisation constants that ensure the IMF is continuous
    b1 = 1 / ( 
                (m2**(1 - a12) - m1**(1 - a12)) / (1 - a12) \
                + m2**(-(a12 - a23)) * (m3**(1 - a23) - m2**(1 - a23)) / (1 - a23) \
                + m2**(-(a12 - a23)) * m3**(-(a23 - a34)) * (m4**(1 - a34) - m3**(1 - a34)) / (1 - a34)
                )
    b2 = b1 * m2**(-(a12 - a23))
    b3 = b2 * m3**(-(a23 - a34))

    # evaluate IMF either at a point or for a list of points
    if isinstance(m, float):
        if m < m1:
            return 0
        elif m < m2:
            return b1 * m**(-a12)
        elif m < m3:
            return b2 * m**(-a23)
        elif m < m4:
            return b3 * m**(-a34)
        else:
            return 0
    else:
        imf_vals = np.zeros(len(m))
        imf_vals[np.logical_and(m >= m1, m < m2)] = b1 * m[np.logical_and(m >= m1, m < m2)]**(-a12)
        imf_vals[np.logical_and(m >= m2, m < m3)] = b2 * m[np.logical_and(m >= m2, m < m3)]**(-a23)
        imf_vals[np.logical_and(m >= m3, m < m4)] = b3 * m[np.logical_and(m >= m3, m < m4)]**(-a34)
        return imf_vals

def CDF_IMF(m, m1=0.01, m2=0.08, m3=0.5, m4=200.0, a12=0.3, a23=1.3, a34=2.3):
    """
        Calculate the fraction of stellar mass between 0 and m for a three part broken power law.
        Default values follow Kroupa (2001)
            F(m) ~ int_0^m zeta(m) dm
        
        Args:
            m       --> [float, list of floats] mass or masses at which to evaluate
            mi      --> [float]                 masses at which to transition the slope
            aij     --> [float]                 slope of the IMF between mi and mj
            
        Returns:
            zeta(m) --> [float, list of floats] value or values of the IMF at m

        NOTE: this is implemented recursively, probably not the most efficient if you're using this
                intensively but I'm not and it looks prettier so I'm being lazy ¯\_(ツ)_/¯ 
    """

    # calculate normalisation constants that ensure the IMF is continuous
    b1 = 1 / ( 
                (m2**(1 - a12) - m1**(1 - a12)) / (1 - a12) \
                + m2**(-(a12 - a23)) * (m3**(1 - a23) - m2**(1 - a23)) / (1 - a23) \
                + m2**(-(a12 - a23)) * m3**(-(a23 - a34)) * (m4**(1 - a34) - m3**(1 - a34)) / (1 - a34)
                )
    b2 = b1 * m2**(-(a12 - a23))
    b3 = b2 * m3**(-(a23 - a34))

    if isinstance(m, float):
        if m <= m1:
            return 0
        elif m <= m2:
            return b1 / (1 - a12) * (m**(1 - a12) - m1**(1 - a12))
        elif m <= m3:
            return CDF_IMF(m2) + b2 / (1 - a23) * (m**(1 - a23) - m2**(1 - a23))
        elif m <= m4:
            return CDF_IMF(m3) + b3 / (1 - a34) * (m**(1 - a34) - m3**(1 - a34))
        else:
            return 0
    else:
        CDF = np.zeros(len(m))
        CDF[np.logical_and(m >= m1, m < m2)] = b1 / (1 - a12) * (m[np.logical_and(m >= m1, m < m2)]**(1 - a12) - m1**(1 - a12))
        CDF[np.logical_and(m >= m2, m < m3)] = CDF_IMF(m2) + b2 / (1 - a23) * (m[np.logical_and(m >= m2, m < m3)]**(1 - a23) - m2**(1 - a23))
        CDF[np.logical_and(m >= m3, m < m4)] = CDF_IMF(m3) + b3 / (1 - a34) * (m[np.logical_and(m >= m3, m < m4)]**(1 - a34) - m3**(1 - a34))
        CDF[m >= m4] = np.ones(len(m[m >= m4]))
        return CDF

def inverse_CDF_IMF(U, m1=0.01, m2=0.08, m3=0.5, m4=200, a12=0.3, a23=1.3, a34=2.3):
    """ 
        Calculate the inverse CDF for a three part broken power law.
        Default values follow Kroupa (2001)
        
        Args:
            U       --> [float, list of floats] A uniform random variable on [0, 1]
            mi      --> [float]                 masses at which to transition the slope
            aij     --> [float]                 slope of the IMF between mi and mj
            
        Returns:
            zeta(m) --> [float, list of floats] value or values of the IMF at m

        NOTE: this is implemented recursively, probably not the most efficient if you're using this intensively but I'm not so I'm being lazy ¯\_(ツ)_/¯ 
    """
    # calculate normalisation constants that ensure the IMF is continuous
    b1 = 1 / ( 
                (m2**(1 - a12) - m1**(1 - a12)) / (1 - a12) \
                + m2**(-(a12 - a23)) * (m3**(1 - a23) - m2**(1 - a23)) / (1 - a23) \
                + m2**(-(a12 - a23)) * m3**(-(a23 - a34)) * (m4**(1 - a34) - m3**(1 - a34)) / (1 - a34)
                )
    b2 = b1 * m2**(-(a12 - a23))
    b3 = b2 * m3**(-(a23 - a34))

    # find the probabilities at which the gradient changes
    F1, F2, F3, F4 = CDF_IMF(np.array([m1, m2, m3, m4]), m1=0.01, m2=0.08, m3=0.5, m4=200, a12=0.3, a23=1.3, a34=2.3)

    masses = np.zeros(len(U))
    masses[np.logical_and(U > F1, U <= F2)] = np.power((1 - a12) / b1 * (U[np.logical_and(U > F1, U <= F2)] - F1) + m1**(1 - a12), 1 / (1 - a12))
    masses[np.logical_and(U > F2, U <= F3)] = np.power((1 - a23) / b2 * (U[np.logical_and(U > F2, U <= F3)] - F2) + m2**(1 - a23), 1 / (1 - a23))
    masses[np.logical_and(U > F3, U <= F4)] = np.power((1 - a34) / b3 * (U[np.logical_and(U > F3, U <= F4)] - F3) + m3**(1 - a34), 1 / (1 - a34))
    return masses
