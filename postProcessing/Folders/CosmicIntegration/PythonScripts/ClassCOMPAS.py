#!/usr/bin/env python3
import numpy as np
import h5py as h5
import os
import totalMassEvolvedPerZ as MPZ
import astropy.units as u

class COMPASData(object):
    def __init__(
        self,
        path=None,
        fileName="COMPAS_output.h5",
        lazyData=True,
        Mlower=None,
        Mupper=None,
        m2_min=None,
        binaryFraction=None,
        suppress_reminder=False,
    ):
        self.path = path
        self.fileName = fileName
        if self.path is None:
            print("Just to double check you create instance of ClassCOMPAS without path/Data")
        elif not os.path.isfile(path + fileName):
            raise ValueError(
                "h5 file not found. Wrong path given? %s" % (path + fileName)
            )

        # Crucial values to be able to calculate MSSFR
        self.metallicityGrid = None
        self.metallicitySystems = None
        self.delayTimes = None  # Myr
        # Crucial values I need for selection effects
        self.mass1 = None  # Msun
        self.mass2 = None  # Msun
        self.DCOmask = None
        self.allTypesMask = None
        self.BBHmask = None
        self.DNSmask = None
        self.BHNSmask = None
        self.initialZ = None
        self.sw_weights = None
        self.n_systems = None

        # Additional arrays that might be nice to store
        # to more quickly make some plots.
        # If you need more memory might help a tiny bit to not do
        self.lazyData = lazyData
        self.mChirp = None  # Msun
        self.q = None
        self.optimisticmask = None

        # Needed to recover true solar mass evolved
        self.Mlower = Mlower  # Msun
        self.Mupper = Mupper  # Msun
        self.m2_min = m2_min # Msun
        self.binaryFraction = binaryFraction
        self.totalMassEvolvedPerZ = None  # Msun
        self.mass_evolved_per_binary = None # Msun

        if not suppress_reminder:
            print("ClassCOMPAS: Remember to self.setGridAndMassEvolved() [optional]")
            print("                    then self.setCOMPASDCOmask()")
            print("                    then self.setCOMPASData()")

    def setCOMPASDCOmask(
        self, types="BBH", withinHubbleTime=True, pessimistic=True, noRLOFafterCEE=True
    ):
        # By default, we mask for BBHs that merge within a Hubble time, assumming
        # the pessimistic CEE prescription (HG donors cannot survive a CEE) and
        # not allowing immediate RLOF post-CEE
        stellar_type_1, stellar_type_2, hubble_flag, dco_seeds = \
            self.get_COMPAS_variables("DoubleCompactObjects", ["Stellar_Type(1)", "Stellar_Type(2)", "Merges_Hubble_Time", "SEED"])

        # if user wants to mask on Hubble time use the flag, otherwise just set all to True
        hubble_mask = hubble_flag if withinHubbleTime else True

        # mask on stellar types (where 14=BH and 13=NS), BHNS can be BHNS or NSBH
        type_masks = {
            "all": True,
            "BBH": np.logical_and(stellar_type_1 == 14, stellar_type_2 == 14),
            "BHNS": np.logical_or(np.logical_and(stellar_type_1 == 14, stellar_type_2 == 13), np.logical_and(stellar_type_1 == 13, stellar_type_2 == 14)),
            "BNS": np.logical_and(stellar_type_1 == 13, stellar_type_2 == 13),
        }

        # if the user wants to make RLOF or optimistic CEs
        if noRLOFafterCEE or pessimistic:

            # get the flags and unique seeds from the Common Envelopes file
            ce_seeds = self.get_COMPAS_variables("CommonEnvelopes", "SEED")
            dco_from_ce = np.in1d(ce_seeds, dco_seeds)
            dco_ce_seeds = ce_seeds[dco_from_ce]

            # if masking on RLOF, get flag and match seeds to dco seeds
            if noRLOFafterCEE:
                rlof_flag = self.get_COMPAS_variables("CommonEnvelopes", "Immediate_RLOF>CE")[dco_from_ce]
                rlof_seeds = np.unique(dco_ce_seeds[rlof_flag])
                rlof_mask = np.logical_not(np.in1d(dco_seeds, rlof_seeds))
            else:
                rlof_mask = True

            # if masking on pessimistic CE, get flag and match seeds to dco seeds
            if pessimistic:
                pessimistic_flag = self.get_COMPAS_variables("CommonEnvelopes", "Optimistic_CE")[dco_from_ce]
                pessimistic_seeds = np.unique(dco_ce_seeds[pessimistic_flag])
                pessimistic_mask = np.logical_not(np.in1d(dco_seeds, pessimistic_seeds))
            else:
                pessimistic_mask = True
        else:
            rlof_mask = True
            pessimistic_mask = True

        # create a mask for each dco type supplied
        self.DCOmask = type_masks[types] * hubble_mask * rlof_mask * pessimistic_mask
        self.BBHmask = type_masks["BBH"] * hubble_mask * rlof_mask * pessimistic_mask
        self.BHNSmask = type_masks["BHNS"] * hubble_mask * rlof_mask * pessimistic_mask
        self.DNSmask = type_masks["BNS"] * hubble_mask * rlof_mask * pessimistic_mask
        self.allTypesMask = type_masks["all"] * hubble_mask * rlof_mask * pessimistic_mask
        self.optimisticmask = pessimistic_mask

    def setGridAndMassEvolved(self):
        # The COMPAS simulation does not evolve all stars
        # give me the correction factor for the total mass evolved
        # I assume each metallicity has the same limits, and does correction
        # factor, but the total mass evolved might be different.
        # This does not change when we change types and other masks this is
        # general to the entire simulation so calculate once
        _, self.totalMassEvolvedPerZ = MPZ.totalMassEvolvedPerZ(
            path=self.path,
            fileName=self.fileName,
            Mlower=self.Mlower,
            Mupper=self.Mupper,
            binaryFraction=self.binaryFraction,
        )
        # Want to recover entire metallicity grid, assume that every metallicity
        # evolved shows in all systems again should not change within same run
        # so dont redo if we reset the data
        Data = h5.File(self.path + self.fileName, "r")
        if self.initialZ is None:
            self.initialZ = Data["SystemParameters"]["Metallicity@ZAMS_1"][()]
        self.metallicityGrid = np.unique(self.initialZ)
        Data.close()

    def setCOMPASData(self):
        
        primary_masses, secondary_masses, formation_times, coalescence_times, dco_seeds = \
            self.get_COMPAS_variables("DoubleCompactObjects", ["Mass(1)", "Mass(2)", "Time", "Coalescence_Time", "SEED"])

        initial_seeds, initial_Z = self.get_COMPAS_variables("SystemParameters", ["SEED", "Metallicity@ZAMS(1)"])

        # Get metallicity grid of DCOs
        self.seedsDCO = dco_seeds[self.DCOmask]
        if self.initialZ is None:
            self.initialZ = initial_Z
        maskMetallicity = np.in1d(initial_seeds, self.seedsDCO)
        self.metallicitySystems = self.initialZ[maskMetallicity]
        self.n_systems = len(initial_seeds)

        self.delayTimes = np.add(formation_times[self.DCOmask], coalescence_times[self.DCOmask])
        self.mass1 = primary_masses[self.DCOmask]
        self.mass2 = secondary_masses[self.DCOmask]

        # Stuff of data I dont need for integral
        # but I might be to laze to read in myself
        # and often use. Might turn it of for memory efficiency
        if self.lazyData:
            self.q = np.divide(self.mass2, self.mass1)
            boolq = self.mass2 > self.mass1
            self.q[boolq] = np.divide(self.mass1[boolq], self.mass2[boolq])
            self.mChirp = np.divide(
                (np.multiply(self.mass2, self.mass1) ** (3.0 / 5.0)),
                (np.add(self.mass2, self.mass1) ** (1.0 / 5.0)),
            )
            self.Hubble = self.get_COMPAS_variables("DoubleCompactObjects", "Merges_Hubble_Time")[self.DCOmask]

    def recalculateTrueSolarMassEvolved(self, Mlower, Mupper, binaryFraction):
        # Possibility to test assumptions of True solar mass evolved
        self.Mlower = Mlower
        self.Mupper = Mupper
        self.binaryFraction = binaryFraction
        _, self.totalMassEvolvedPerZ = MPZ.totalMassEvolvedPerZ(
            pathCOMPASh5=self.path,
            Mlower=self.Mlower,
            Mupper=self.Mupper,
            binaryFraction=self.binaryFraction,
        )

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
        with h5.File(self.path + self.fileName, "r") as compas_file:
            # if the list is only a string (i.e. one variable) then don't return a list
            if isinstance(var_names, str):
                return compas_file[hdf5_file][var_names][...].squeeze()
            # else return each variable in a list
            else:
                return [compas_file[hdf5_file][var_name][...].squeeze() for var_name in var_names]

    def set_sw_weights(self, column_name):
        """ Set STROOPWAFEL adaptive sampling weights given a column name in the DoubleCompactObjects file """
        if column_name is not None:
            self.sw_weights = self.get_COMPAS_variables("DoubleCompactObjects", column_name)[self.DCOmask]

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
        binary_mask = binary < self.binaryFraction
        
        # assign each a random secondary mass, default 0 because single stars have m2=0 (surprisingly :P)
        secondary_mass = np.zeros(SAMPLES) * u.Msun
        secondary_mass[binary_mask] = primary_mass[binary_mask] * mass_ratio[binary_mask]

        # find the total mass of the whole population
        total_mass = np.sum(primary_mass) + np.sum(secondary_mass)

        # apply the COMPAS cuts on primary and secondary mass
        primary_mask = np.logical_and(primary_mass >= self.Mlower, primary_mass <= self.Mupper)
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