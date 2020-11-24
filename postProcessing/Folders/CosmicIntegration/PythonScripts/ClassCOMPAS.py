#!/usr/bin/env python3
import numpy as np
import h5py as h5
import os
import astropy.units as u
import astropy.constants as c

class COMPASData(object):
    """ A Class for reading and masking a COMPAS data file """

    def __init__(self, path, filename="COMPAS_output.h5", m1_min=5 * u.Msun, m1_max=150 * u.Msun, m2_min=0.1 * u.Msun, fbin=0.7):
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

    def set_DCO_masks(self, dco_types=["ALL", "BHBH", "BHNS", "NSNS"], merges_in_hubble_time=True, no_RLOF_after_CEE=True, pessimistic_CEE=True,
                            rlof_mask=None, pessimistic_mask=None):
        """
            Create masks for different DCO types

            Args:
                dco_types             --> which DCO types to create masks for (default is all of them)
                merges_in_hubble_time --> whether to mask binaries that don't merge in a Hubble time
                no_RLOF_after_CEE     --> whether to mask binaries that have immediate RLOF after a CCE
                pessimistic_CEE       --> whether to mask binaries that go through Optimistic CE scenario
                rlof_mask             --> a mask for RLOF after CEE (default is taken from CommonEnvelopes file, this is for if you don't have CE file)
                pessimistic_mask      --> a mask for Optimstic CE scenarios (default is taken from CommonEnvelopes file, this is for if you don't have CE file)
        """
        # get the appropriate variables from the COMPAS file
        stellar_type_1, stellar_type_2, hubble_flag, dco_seeds = \
            self.get_COMPAS_variables("DoubleCompactObjects", ["Stellar_Type_1", "Stellar_Type_2", "Merges_Hubble_Time", "SEED"])

        # if user wants to mask on Hubble time use the flag, otherwise just set all to True
        hubble_mask = hubble_flag if merges_in_hubble_time else True

        # mask on stellar types (where 14=BH and 13=NS), BHNS can be BHNS or NSBH
        type_masks = {
            "ALL": True,
            "BHBH": np.logical_and(stellar_type_1 == 14, stellar_type_2 == 14),
            "BHNS": np.logical_or(np.logical_and(stellar_type_1 == 14, stellar_type_2 == 13), np.logical_and(stellar_type_1 == 13, stellar_type_2 == 14)),
            "NSNS": np.logical_and(stellar_type_1 == 13, stellar_type_2 == 13),
        }

        # if user wants to mask RLOF or pessimistic and they haven't supplied masks then we need to create them
        if (no_RLOF_after_CEE or pessimistic_CEE) and (rlof_mask is None or pessimistic_mask is None):
            # Try to get the flags and unique seeds from the Common Envelopes file
            try:
                ce_seeds, rlof_flag, pessimistic_flag = self.get_COMPAS_variables("CommonEnvelopes", ["SEED", "Immediate_RLOF>CE", "Optimistic_CE"])
            except:
                ce_seeds, rlof_flag, pessimistic_flag = self.get_COMPAS_variables("SystemParameters", ["SEED", "Immediate_RLOF>CE", "Optimistic_CE"])

            # match the seeds to DCO seeds and only take corresponding flags
            dco_from_ce = np.in1d(ce_seeds, dco_seeds)
            rlof_flag, pessimistic_flag = rlof_flag[dco_from_ce], pessimistic_flag[dco_from_ce]
            
            # if user wants to mask on immediate RLOF use the inverse of the flag
            rlof_mask = np.logical_not(rlof_flag)

            # if user wants to mask on pessimistic CEE use the inverse of the flag
            pessimistic_mask = np.logical_not(pessimistic_flag)

        # if the user doesn't want to mask then just set to true
        if not no_RLOF_after_CEE:
            rlof_mask = True
        if not pessimistic_CEE:
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
            self.get_COMPAS_variables("DoubleCompactObjects", ["Mass_1", "Mass_2", "Time", "Coalescence_Time"])
        self.primary_masses, self.secondary_masses = primary_masses * u.Msun, secondary_masses * u.Msun
        self.delay_times = (formation_times + coalescence_times) * u.Myr

        dco_seeds = self.get_COMPAS_variables("DoubleCompactObjects", "SEED")
        systems_seeds, systems_metallicites = self.get_COMPAS_variables("SystemParameters", ["SEED", "Metallicity@ZAMS_1"])
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

    def find_star_forming_mass_per_binary_sampling(self):
        # this is a placeholder, will update with proper function
        self.mass_evolved_per_binary = 115