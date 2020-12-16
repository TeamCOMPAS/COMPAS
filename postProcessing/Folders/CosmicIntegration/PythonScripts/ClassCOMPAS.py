#!/usr/bin/env python3
import numpy as np
import h5py as h5
import os
import totalMassEvolvedPerZ as MPZ

class COMPASData(object):
    def __init__(
        self,
        path=None,
        fileName="COMPAS_output.h5",
        lazyData=True,
        Mlower=None,
        Mupper=None,
        binaryFraction=None,
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
        self.binaryFraction = binaryFraction
        self.totalMassEvolvedPerZ = None  # Msun
        self.mass_evolved_per_binary = None # Msun

        print("ClassCOMPAS: Remember to self.setGridAndMassEvolved() [optional]")
        print("                   then  self.setCOMPASDCOmask()")
        print("                   then  self.setCOMPASData()")

    def setCOMPASDCOmask(
        self, types="BBH", withinHubbleTime=True, pessimistic=True, noRLOFafterCEE=True
    ):
        # By default, we mask for BBHs that merge within a Hubble time, assumming
        # the pessimistic CEE prescription (HG donors cannot survive a CEE) and
        # not allowing immediate RLOF post-CEE
        Data = h5.File(self.path + self.fileName, "r")
        fDCO = Data["DoubleCompactObjects"]
        fCEE = Data["CommonEnvelopes"]

        # Masks for DCO type
        maskBBH = (fDCO["Stellar_Type_1"][()] == 14) & (fDCO["Stellar_Type_2"][()] == 14)
        maskDNS = (fDCO["Stellar_Type_1"][()] == 13) & (fDCO["Stellar_Type_2"][()] == 13)
        maskBHNS = ((fDCO["Stellar_Type_1"][()] == 14) & (fDCO["Stellar_Type_2"][()] == 13)) \
                    |((fDCO["Stellar_Type_1"][()] == 13) & (fDCO["Stellar_Type_2"][()] == 14))
        maskAllTypes = maskBBH | maskDNS | maskBHNS

        if types == "BBH":
            maskTypes = maskBBH
        elif types == "BNS":
            maskTypes = maskDNS
        elif types == "BHNS":
            maskTypes = maskBHNS
        elif types == "all":
            maskTypes = maskAllTypes
        else:
            raise ValueError("type=%s not one of 'BBH', 'BNS', 'BHNS', 'all'" % (types))

        # Mask DCOs merging within Hubble time
        if withinHubbleTime:
            maskHubble = fDCO["Merges_Hubble_Time"][()] == True
        else:
            maskHubble = np.ones(len(fDCO["Merges_Hubble_Time"][()]), dtype=bool)

        # Masks related to CEEs
        CEEseeds = fCEE["SEED"][()]
        DCOseeds = fDCO["SEED"][()]

        if (pessimistic or noRLOFafterCEE):
            maskDCOCEEs = np.in1d(CEEseeds, DCOseeds) # Mask for CEEs involved in forming DCOs
            DCOCEEseeds = CEEseeds[maskDCOCEEs] # Seeds of CEEs involved in forming DCO

        # Mask for DCOs formed assumming pessimistic CEE
        if pessimistic:
            optimisticFlagForDCOCEEs = fCEE["Optimistic_CE"][()][maskDCOCEEs]
            # Seeds of DCOs that have optimistic CEE
            optimisticCEEseeds = np.unique(DCOCEEseeds[optimisticFlagForDCOCEEs])
            maskOptimistic = np.in1d(DCOseeds,optimisticCEEseeds)
            maskPessimistic = np.logical_not(maskOptimistic)
        else:
            maskPessimistic = np.ones(len(fDCO["ID"][()]), dtype=bool)

        # Mask for DCOs formed without RLOF immediately after CEE
        if noRLOFafterCEE:
            immediateRLOFflagforDCOCEEs = fCEE["Immediate_RLOF>CE"][()][maskDCOCEEs]
            # Seeds of DCOs that have immediate RLOF post-CEE
            immediateRLOFCEEseeds = np.unique(DCOCEEseeds[immediateRLOFflagforDCOCEEs])
            maskRLOFafterCEE = np.in1d(DCOseeds,immediateRLOFCEEseeds)
            maskNoRLOFafterCEE = np.logical_not(maskRLOFafterCEE)
        else:
            maskNoRLOFafterCEE = np.ones(len(fDCO["ID"][()]), dtype=bool)

        # Combine all the masks
        self.DCOmask = maskTypes & maskHubble & maskPessimistic & maskNoRLOFafterCEE
        self.BBHmask = maskBBH & maskHubble & maskPessimistic & maskNoRLOFafterCEE
        self.DNSmask = maskDNS & maskHubble & maskPessimistic & maskNoRLOFafterCEE
        self.BHNSmask = maskBHNS & maskHubble & maskPessimistic & maskNoRLOFafterCEE
        self.allTypesMask = maskAllTypes & maskHubble & maskPessimistic & maskNoRLOFafterCEE
        self.optimisticmask = maskPessimistic
        Data.close()

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
        Data = h5.File(self.path + self.fileName, "r")
        fDCO = Data["DoubleCompactObjects"]
        # sorry not the prettiest line is a boolean slice of seeds
        # this only works because seeds in systems file and DCO file are printed
        # in same order

        # Get metallicity grid of DCOs
        self.seedsDCO = fDCO["SEED"][()][self.DCOmask]
        initialSeeds = Data["SystemParameters"]["SEED"][()]
        if self.initialZ is None:
            self.initialZ = Data["SystemParameters"]["Metallicity@ZAMS_1"][()]
        maskMetallicity = np.in1d(initialSeeds, self.seedsDCO)
        self.metallicitySystems = self.initialZ[maskMetallicity]

        self.delayTimes = np.add(
            fDCO["Time"][()][self.DCOmask], fDCO["Coalescence_Time"][()][self.DCOmask]
        )
        self.mass1 = fDCO["Mass_1"][()][self.DCOmask]
        self.mass2 = fDCO["Mass_2"][()][self.DCOmask]

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
            self.Hubble = fDCO["Merges_Hubble_Time"][...].squeeze()[self.DCOmask]

        Data.close()

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
        secondary_mask = secondary_mass > 0.1 * u.Msun
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