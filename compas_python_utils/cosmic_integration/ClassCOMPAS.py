#!/usr/bin/env python3
import numpy as np
import h5py as h5
import os
from . import totalMassEvolvedPerZ as MPZ


class COMPASData(object):
    def __init__(
        self,
        path=None,
        lazyData=True,
        Mlower=None,
        Mupper=None,
        m2_min=None,
        binaryFraction=None,
        suppress_reminder=False,
    ):
        self.path = path
        if self.path is None:
            print("Template COMPASData object created with no data path")
        elif not os.path.isfile(path):
            raise ValueError( "h5 file not found. Wrong path given? {}".format(path))

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
        self.CHE_mask = None
        self.CHE_BBHmask = None
        self.NonCHE_BBHmask = None
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
            print("ClassCOMPAS: Remember to self.setCOMPASDCOmask()")
            print("                    then self.setCOMPASData()")
            print("          and optionally self.setGridAndMassEvolved() if using a metallicity grid")

    def setCOMPASDCOmask(
        self, types="BBH", withinHubbleTime=True, pessimistic=True, noRLOFafterCEE=True
    ):
        # By default, we mask for BBHs that merge within a Hubble time, assuming
        # the pessimistic CEE prescription (HG donors cannot survive a CEE) and
        # not allowing immediate RLOF post-CEE
        
        stellar_type_1, stellar_type_2, hubble_flag, dco_seeds = \
            self.get_COMPAS_variables("BSE_Double_Compact_Objects", ["Stellar_Type(1)", "Stellar_Type(2)", "Merges_Hubble_Time", "SEED"])
        dco_seeds = dco_seeds.flatten()

        if types == "CHE_BBH" or types == "NON_CHE_BBH":
            stellar_type_1_zams, stellar_type_2_zams, che_ms_1, che_ms_2, sys_seeds = \
                self.get_COMPAS_variables("BSE_System_Parameters", ["Stellar_Type@ZAMS(1)", "Stellar_Type@ZAMS(2)", "CH_on_MS(1)", "CH_on_MS(2)", "SEED"])
          
            che_mask  = np.logical_and.reduce((stellar_type_1_zams == 16, stellar_type_2_zams == 16, che_ms_1 == True, che_ms_2 == True))
            che_seeds = sys_seeds[()][che_mask]

        self.CHE_mask = np.in1d(dco_seeds, che_seeds) if types == "CHE_BBH" or types == "NON_CHE_BBH" else np.repeat(False, len(dco_seeds))

        # if user wants to mask on Hubble time use the flag, otherwise just set all to True, use astype(bool) to set masks to bool type
        hubble_mask = hubble_flag.astype(bool) if withinHubbleTime else np.repeat(True, len(dco_seeds))

        # mask on stellar types (where 14=BH and 13=NS), BHNS can be BHNS or NSBH
        type_masks = {
            "all": np.repeat(True, len(dco_seeds)),
            "BBH": np.logical_and(stellar_type_1 == 14, stellar_type_2 == 14),
            "BHNS": np.logical_or(np.logical_and(stellar_type_1 == 14, stellar_type_2 == 13), np.logical_and(stellar_type_1 == 13, stellar_type_2 == 14)),
            "BNS": np.logical_and(stellar_type_1 == 13, stellar_type_2 == 13),
        }
        type_masks["CHE_BBH"]     = np.logical_and(self.CHE_mask, type_masks["BBH"]) if types == "CHE_BBH" else np.repeat(False, len(dco_seeds))
        type_masks["NON_CHE_BBH"] = np.logical_and(np.logical_not(self.CHE_mask), type_masks["BBH"]) if types == "NON_CHE_BBH" else np.repeat(True, len(dco_seeds))

        # if the user wants to make RLOF or optimistic CEs
        if noRLOFafterCEE or pessimistic:

            # get the flags and unique seeds from the Common Envelopes file
            ce_seeds = self.get_COMPAS_variables("BSE_Common_Envelopes", "SEED")
            dco_from_ce = np.in1d(ce_seeds, dco_seeds)
            dco_ce_seeds = ce_seeds[dco_from_ce]

            # if masking on RLOF, get flag and match seeds to dco seeds
            if noRLOFafterCEE:
                rlof_flag = self.get_COMPAS_variables("BSE_Common_Envelopes", "Immediate_RLOF>CE")[dco_from_ce].astype(bool)
                rlof_seeds = np.unique(dco_ce_seeds[rlof_flag])
                rlof_mask = np.logical_not(np.in1d(dco_seeds, rlof_seeds))
            else:
                rlof_mask = np.repeat(True, len(dco_seeds))

            # if masking on pessimistic CE, get flag and match seeds to dco seeds
            if pessimistic:
                pessimistic_flag = self.get_COMPAS_variables("BSE_Common_Envelopes", "Optimistic_CE")[dco_from_ce].astype(bool)
                pessimistic_seeds = np.unique(dco_ce_seeds[pessimistic_flag])
                pessimistic_mask = np.logical_not(np.in1d(dco_seeds, pessimistic_seeds))
            else:
                pessimistic_mask = np.repeat(True, len(dco_seeds))
        else:
            rlof_mask = np.repeat(True, len(dco_seeds))
            pessimistic_mask = np.repeat(True, len(dco_seeds))

        # create a mask for each dco type supplied
        self.DCOmask = type_masks[types] * hubble_mask * rlof_mask * pessimistic_mask
        self.BBHmask = type_masks["BBH"] * hubble_mask * rlof_mask * pessimistic_mask
        self.BHNSmask = type_masks["BHNS"] * hubble_mask * rlof_mask * pessimistic_mask
        self.DNSmask = type_masks["BNS"] * hubble_mask * rlof_mask * pessimistic_mask
        self.CHE_BBHmask = type_masks["CHE_BBH"] * hubble_mask * rlof_mask * pessimistic_mask
        self.NonCHE_BBHmask = type_masks["NON_CHE_BBH"] * hubble_mask * rlof_mask * pessimistic_mask
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
            Mlower=self.Mlower,
            Mupper=self.Mupper,
            binaryFraction=self.binaryFraction,
        )
        # Want to recover entire metallicity grid, assume that every metallicity
        # evolved shows in all systems again should not change within same run
        # so dont redo if we reset the data
        Data = h5.File(self.path, "r")
        if self.initialZ is None:
            self.initialZ = Data["BSE_System_Parameters"]["Metallicity@ZAMS(1)"][()]
        self.metallicityGrid = np.unique(self.initialZ)
        Data.close()

    def setCOMPASData(self):
        
        primary_masses, secondary_masses, formation_times, coalescence_times, dco_seeds = \
            self.get_COMPAS_variables("BSE_Double_Compact_Objects", ["Mass(1)", "Mass(2)", "Time", "Coalescence_Time", "SEED"])

        initial_seeds, initial_Z = self.get_COMPAS_variables("BSE_System_Parameters", ["SEED", "Metallicity@ZAMS(1)"])

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
            self.Hubble = self.get_COMPAS_variables("BSE_Double_Compact_Objects", "Merges_Hubble_Time")[self.DCOmask]

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
                hdf5_file --> [string]                Name of HDF5 subfile (e.g. "BSE_Double_Compact_Objects")
                var_names --> [string or string list] A variable name or list of variables names to return

            Returns:
                var_list  --> [list of lists]         A list of variables (or a single variable if only one name supplied)
        """
        # open the COMPAS file
        with h5.File(self.path, "r") as compas_file:
            # if the list is only a string (i.e. one variable) then don't return a list
            if isinstance(var_names, str):
                return compas_file[hdf5_file][var_names][...].squeeze()
            # else return each variable in a list
            else:
                return [compas_file[hdf5_file][var_name][...].squeeze() for var_name in var_names]

    def set_sw_weights(self, column_name):
        """ Set STROOPWAFEL adaptive sampling weights given a column name in the BSE_Double_Compact_Objects file """
        if column_name is not None:
            self.sw_weights = self.get_COMPAS_variables("BSE_Double_Compact_Objects", column_name)[self.DCOmask]

    def find_star_forming_mass_per_binary_sampling(self):
        self.mass_evolved_per_binary = MPZ.analytical_star_forming_mass_per_binary_using_kroupa_imf(
            m1_max=self.Mupper.value,
            m1_min=self.Mlower.value,
            m2_min=self.m2_min.value,
            fbin=self.binaryFraction,
        )


