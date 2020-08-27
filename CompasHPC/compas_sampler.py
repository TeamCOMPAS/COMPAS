#!/usr/bin/env python

import os, sys
import pandas as pd
import shutil
import time
import numpy as np
from stroopwafel import sw, classes, prior, sampler, distributions, constants, utils
import argparse

########################################################
### 
###         USER GUIDE:
###
### To see all command line options for this sampler, run 
###    python sw_interface.py --help
###
### User can set these 'on the fly' on the command line, 
### or change the defaults below for more permanent changes.
###
### Below, User should set which parameters to sample and 
### specify their distributions.
###
### Then, User should set physical prescriptions and constants.
### To see explanations for all these parameters, run
###    src/COMPAS --help 
###
########################################################


# Optional user-supplied arguments
parser=argparse.ArgumentParser()
parser.add_argument('--num_systems', help = 'Total number of systems', type = int, default = 1000)
parser.add_argument('--num_cores', help = 'Number of cores to run in parallel', type = int, default = 2)
parser.add_argument('--num_per_core', help = 'Number of systems to generate in one core', type = int, default = 10)
parser.add_argument('--debug', help = 'If debug of COMPAS is to be printed', type = bool, default = False)
parser.add_argument('--mc_only', help = 'If run in MC simulation mode only', type = bool, default = True)
parser.add_argument('--run_on_helios', help = 'If we are running on helios (or other slurm) nodes', type = bool, default = False)
parser.add_argument('--output_filename', help = 'Output filename', default = 'samples.csv')
namespace, _ = parser.parse_known_args()


########################################################
###
### STROOPWAFEL SAMPLING FUNCTIONS
###
### User should specify which parameters to sample from which distributions
###
########################################################

def create_dimensions():
    """
    This Function that will create all the dimensions for stroopwafel, a dimension is basically one of the variables you want to sample
    Invoke the Dimension class to create objects for each variable. Look at the Dimension class definition in classes.py for more.
    It takes the name of the dimension, its max and min value. 
    The Sampler class will tell how to sample this dimension. Similarly, prior tells it how it calculates the prior. You can find more of these in their respective modules
    OUT:
        As Output, this should return a list containing all the instances of Dimension class.
    """
    m1 = classes.Dimension('Mass_1', 40, 50, sampler.kroupa, prior.kroupa)
    q = classes.Dimension('q', 0.2, 1, sampler.uniform, prior.uniform, should_print = False)
    a = classes.Dimension('Separation', .01, 100, sampler.flat_in_log, prior.flat_in_log)
    #kick_velocity_random_1 = classes.Dimension('Kick_Velocity_Random_1', 0, 1, sampler.uniform, prior.uniform)
    #kick_theta_1 = classes.Dimension('Kick_Theta_1', -np.pi / 2, np.pi / 2, sampler.uniform_in_cosine, prior.uniform_in_cosine)
    #kick_phi_1 = classes.Dimension('Kick_Phi_1', 0, 2 * np.pi, sampler.uniform, prior.uniform)
    #kick_velocity_random_2 = classes.Dimension('Kick_Velocity_Random_2', 0, 1, sampler.uniform, prior.uniform)
    #kick_theta_2 = classes.Dimension('Kick_Theta_2', -np.pi / 2, np.pi / 2, sampler.uniform_in_cosine, prior.uniform_in_cosine)
    #kick_phi_2 = classes.Dimension('Kick_Phi_2', 0, 2 * np.pi, sampler.uniform, prior.uniform)
    #return [m1, q, a, kick_velocity_random_1, kick_theta_1, kick_phi_1, kick_velocity_random_2, kick_theta_2, kick_phi_2]
    return [m1, q, a]

def update_properties(locations, dimensions):
    """
    This function is not mandatory, it is required only if you have some dependent variable. 
    For example, if you want to sample Mass_1 and q, then Mass_2 is a dependent variable which is product of the two.
    Similarly, you can assume that Metallicity_2 will always be equal to Metallicity_1
    IN:
        locations (list(Location)) : A list containing objects of Location class in classes.py. 
        You can play with them and update whatever fields you like or add more in the property (which is a dictionary)
    OUT: Not Required
    """
    m1 = dimensions[0]
    q = dimensions[1]
    for location in locations:
        location.properties['Mass_2'] = location.dimensions[m1] * location.dimensions[q]
        location.properties['Metallicity_2'] = location.properties['Metallicity_1'] = constants.METALLICITY_SOL
        location.properties['Eccentricity'] = 0
        #location.properties['Kick_Mean_Anomaly_1'] = np.random.uniform(0, 2 * np.pi, 1)[0]
        #location.properties['Kick_Mean_Anomaly_2'] = np.random.uniform(0, 2 * np.pi, 1)[0]



########################################################
###
### COMPAS PARAMETER SELECTION
###
### User should choose which prescriptions 
###   and constants to run with COMPAS
###
########################################################


class pythonProgramOptions:
    """
    A class to store and access COMPAS program options in python
    """

    #-- Define variables
    enable_warnings = False                                     # option to enable/disable warning messages
    populationPrinting = False

    randomSeedFileName = 'randomSeed.txt'
    if os.path.isfile(randomSeedFileName):
        random_seed = int(np.loadtxt(randomSeedFileName))
    else:
        random_seed = 0 # If you want a random seed, use: np.random.randint(2,2**63-1)

    shareSeeds = False

    single_star = False
    single_star_mass_steps = 10
    single_star_mass_min   = 1.0
    single_star_mass_max   = 75.0    

    use_mass_loss = True
    mass_transfer = True
    detailed_output = False                             # WARNING: this creates a data heavy file
    evolve_unbound_systems = False
    quiet = False

    metallicity = 0.0142                                # Solar metallicity Asplund+2010

    allow_rlof_at_birth = False;                                            # allow binaries that have one or both stars in RLOF at birth to evolve?
    allow_touching_at_birth = False;                                        # allow binaries that have stars touching at birth to evolve?

    chemically_homogeneous_evolution = 'NONE'                               # chemically homogeneous evolution.  Options are 'NONE', 'OPTIMISTIC' and 'PESSIMISTIC'

    common_envelope_alpha = 1.0
    common_envelope_lambda = 0.1                # Only if using 'LAMBDA_FIXED'
    common_envelope_lambda_prescription = 'LAMBDA_NANJING'  # Xu & Li 2010
    common_envelope_slope_Kruckow = -5.0/6.0
    stellar_zeta_prescription = 'SOBERMAN'
    common_envelope_revised_energy_formalism = False
    common_envelope_maximum_donor_mass_revised_energy_formalism = 2.0
    common_envelope_recombination_energy_density = 1.5E13
    common_envelope_alpha_thermal = 1.0                     # lambda = alpha_th*lambda_b + (1-alpha_th)*lambda_g
    common_envelope_lambda_multiplier = 1.0                 # Multiply common envelope lambda by some constant
    common_envelope_allow_main_sequence_survive = True      # Allow main sequence stars to survive CE. Was previously False by default
    common_envelope_mass_accretion_prescription = 'ZERO'
    common_envelope_mass_accretion_min = 0.04           # For 'MACLEOD+2014' [Msol]
    common_envelope_mass_accretion_max = 0.10           # For 'MACLEOD+2014' [Msol]
    envelope_state_prescription = 'LEGACY'

    mass_loss_prescription = 'VINK'
    luminous_blue_variable_multiplier = 1.5
    wolf_rayet_multiplier = 1.0

    circularise_binary_during_mass_transfer = False
    angular_momentum_conservation_during_circularisation = False
    mass_transfer_angular_momentum_loss_prescription = 'ISOTROPIC'
    mass_transfer_accretion_efficiency_prescription = 'THERMAL'
    mass_transfer_fa = 0.5  # Only if using mass_transfer_accretion_efficiency_prescription = 'FIXED'
    mass_transfer_jloss = 1.0   # Only if using mass_transfer_angular_momentum_loss_prescription = 'FIXED'
    mass_transfer_rejuvenation_prescription = 'STARTRACK'
    mass_transfer_thermal_limit_accretor= 'CFACTOR'
    mass_transfer_thermal_limit_C= 10.0
    eddington_accretion_factor = 1    #multiplication Factor for eddington accretion onto NS&BH

    #-- Stability criteria for case BB/BC mass transfer (for BNS project)
    case_bb_stability_prescription = 'ALWAYS_STABLE'
    zeta_Main_Sequence = 2.0
    zeta_Radiative_Envelope_Giant = 6.5

    maximum_evolution_time = 13700.0                    # Maximum physical time a system can be evolved [Myrs]
    maximum_number_timesteps = 99999

    initial_mass_function = 'KROUPA'
    initial_mass_min = 5.0                              # Use 1.0 for LRNe, 5.0 for DCOs  [Msol]
    initial_mass_max = 150.0                            # Stellar tracks extrapolated above 50 Msol (Hurley+2000) [Msol]

    initial_mass_power = 0.0

    semi_major_axis_distribution = 'FLATINLOG'
    semi_major_axis_min = 0.01                          #  [AU]
    semi_major_axis_max = 1000.0                        #  [AU]

    mass_ratio_distribution = 'FLAT'
    mass_ratio_min = 0.0
    mass_ratio_max = 1.0

    minimum_secondary_mass = 0.1                        # Brown dwarf limit  [Msol]

    eccentricity_distribution = 'ZERO'
    eccentricity_min = 0.0
    eccentricity_max = 1.0

    pulsar_birth_magnetic_field_distribution = 'ZERO'
    pulsar_birth_magnetic_field_min = 11.0              # [log10(B/G)]
    pulsar_birth_magnetic_field_max = 13.0              # [log10(B/G)]

    pulsar_birth_spin_period_distribution = "ZERO"
    pulsar_birth_spin_period_min = 10.0                 # [ms]
    pulsar_birth_spin_period_max = 100.0                # [ms]

    pulsar_magnetic_field_decay_timescale = 1000.0      # [Myrs]
    pulsar_magnetic_field_decay_massscale = 0.025       # [Msol]
    pulsar_minimum_magnetic_field = 8.0                 # [log10(B/G)]

    evolvePulsars = False

    rotational_velocity_distribution = 'ZERO'

    neutron_star_equation_of_state = 'SSE'

    orbital_period_min = 1.1
    orbital_period_max = 1000

    remnant_mass_prescription = 'FRYER2012'
    fryer_supernova_engine = 'DELAYED'
    black_hole_kicks = 'FALLBACK'
    kick_magnitude_distribution = 'MAXWELLIAN'

    kick_magnitude_sigma_CCSN_NS = 265.0                 #  [km/s]
    kick_magnitude_sigma_CCSN_BH = 265.0                 #  [km/s]
    kick_magnitude_sigma_ECSN = 30.0                     #  [km/s]
    kick_magnitude_sigma_USSN = 30.0                     #  [km/s]

    fix_dimensionless_kick_magnitude = -1
    kick_direction = 'ISOTROPIC'
    kick_direction_power = 0.0
    kick_scaling_factor = 1.0
    kick_magnitude_maximum = -1.0

    pair_instability_supernovae = True
    PISN_lower_limit = 60.0                             # Minimum core mass for PISN [Msol]
    PISN_upper_limit = 135.0                            # Maximum core mass for PISN [Msol]
    pulsation_pair_instability = True
    PPI_lower_limit = 35.0                              # Minimum core mass for PPI [Msol]
    PPI_upper_limit = 60.0                              # Maximum core mass for PPI [Msol]

    pulsational_pair_instability_prescription = 'MARCHANT'

    maximum_neutron_star_mass = 2.5  #  [Msol]

    log_level         = 0
    log_classes       = []

    debug_level       = 0
    debug_classes     = []

    logfile_definitions = None

    logfile_name_prefix = None
    logfile_delimiter   = 'COMMA'

    # set the logfile names here
    #
    # set to None (e.g. logfile_BSE_supernovae = None) to use the default filename
    # set to a string (e.g. logfile_BSE_supernovae = 'mySNfilename') to use that string as the filename 
    # set to empty string (e.g. logfile_BSE_supernovae = '""') to disable logging for that file (the file will not be created)

    logfile_BSE_be_binaries = None
    logfile_BSE_common_envelopes = None
    logfile_BSE_detailed_output = None
    logfile_BSE_double_compact_objects = None
    logfile_BSE_pulsar_evolution = None
    logfile_BSE_supernovae = None
    logfile_BSE_system_parameters = None

    debug_to_file  = False
    errors_to_file = False





    ###################################################################
    ###################################################################
    ###                                                             ###
    ###    User should not modify anything below this line!         ###
    ###                                                             ###
    ###################################################################
    ###################################################################

    def booleanChoices(self):
        booleanChoices = [
            self.enable_warnings,
            self.single_star,
            self.use_mass_loss,
            self.mass_transfer,
            self.detailed_output,
            self.evolve_unbound_systems,
            self.populationPrinting,
            self.circularise_binary_during_mass_transfer,
            self.angular_momentum_conservation_during_circularisation,
            self.pair_instability_supernovae,
            self.pulsation_pair_instability,
            self.quiet,
            self.common_envelope_allow_main_sequence_survive,
            self.evolvePulsars,
            self.debug_to_file,
            self.errors_to_file,
            self.allow_rlof_at_birth,
            self.allow_touching_at_birth
        ]

        return booleanChoices

    def booleanCommands(self):
        booleanCommands = [
            '--enable-warnings',
            '--single-star',
            '--use-mass-loss',
            '--massTransfer',
            '--detailedOutput',
            '--evolve-unbound-systems',
            '--populationDataPrinting',
            '--circulariseBinaryDuringMassTransfer',
            '--angularMomentumConservationDuringCircularisation',
            '--pair-instability-supernovae',
            '--pulsational-pair-instability',
            '--quiet',
            '--common-envelope-allow-main-sequence-survive',
            '--evolve-pulsars',
            '--debug-to-file',
            '--errors-to-file',
            '--allow-rlof-at-birth',
            '--allow-touching-at-birth'
        ]

        return booleanCommands

    def numericalChoices(self):
        numericalChoices = [
            #self.number_of_binaries,
            self.metallicity,
            self.common_envelope_alpha,
            self.common_envelope_lambda,
            self.common_envelope_slope_Kruckow,
            self.common_envelope_alpha_thermal,
            self.common_envelope_lambda_multiplier,
            self.luminous_blue_variable_multiplier,
            self.wolf_rayet_multiplier,
            self.mass_transfer_fa,
            self.mass_transfer_jloss,
            self.maximum_evolution_time,
            self.maximum_number_timesteps,
            self.initial_mass_min,
            self.initial_mass_max,
            self.initial_mass_power,
            self.semi_major_axis_min,
            self.semi_major_axis_max,
            self.mass_ratio_min,
            self.mass_ratio_max,
            self.minimum_secondary_mass,
            self.eccentricity_min,
            self.eccentricity_max,
            self.pulsar_birth_magnetic_field_min,
            self.pulsar_birth_magnetic_field_max,
            self.pulsar_birth_spin_period_min,
            self.pulsar_birth_spin_period_max,
            self.pulsar_magnetic_field_decay_timescale,
            self.pulsar_magnetic_field_decay_massscale,
            self.pulsar_minimum_magnetic_field,
            self.orbital_period_min,
            self.orbital_period_max,
            self.kick_magnitude_sigma_CCSN_NS,
            self.kick_magnitude_sigma_CCSN_BH,
            self.fix_dimensionless_kick_magnitude,
            self.kick_direction_power,
            self.random_seed,
            self.mass_transfer_thermal_limit_C,
            self.eddington_accretion_factor,
            self.PISN_lower_limit,
            self.PISN_upper_limit,
            self.PPI_lower_limit,
            self.PPI_upper_limit,
            self.maximum_neutron_star_mass,
            self.kick_magnitude_sigma_ECSN,
            self.kick_magnitude_sigma_USSN,
            self.kick_scaling_factor,
            self.common_envelope_maximum_donor_mass_revised_energy_formalism,
            self.common_envelope_recombination_energy_density,
            self.common_envelope_mass_accretion_max,
            self.common_envelope_mass_accretion_min,
            self.zeta_Main_Sequence,
            self.zeta_Radiative_Envelope_Giant,
            self.kick_magnitude_maximum,
            self.log_level,
            self.debug_level,
            self.single_star_mass_steps,
            self.single_star_mass_min,
            self.single_star_mass_max
        ]

        return numericalChoices

    def numericalCommands(self):
        numericalCommands = [
            #'--number-of-binaries',
            '--metallicity',
            '--common-envelope-alpha',
            '--common-envelope-lambda',
            '--common-envelope-slope-Kruckow',
            '--common-envelope-alpha-thermal',
            '--common-envelope-lambda-multiplier',
            '--luminous-blue-variable-multiplier',
            '--wolf-rayet-multiplier',
            '--mass-transfer-fa',
            '--mass-transfer-jloss',
            '--maximum-evolution-time',
            '--maximum-number-timestep-iterations',
            '--initial-mass-min',
            '--initial-mass-max',
            '--initial-mass-power',
            '--semi-major-axis-min',
            '--semi-major-axis-max',
            '--mass-ratio-min',
            '--mass-ratio-max',
            '--minimum-secondary-mass',
            '--eccentricity-min',
            '--eccentricity-max',
            '--pulsar-birth-magnetic-field-distribution-min',
            '--pulsar-birth-magnetic-field-distribution-max',
            '--pulsar-birth-spin-period-distribution-min',
            '--pulsar-birth-spin-period-distribution-max',
            '--pulsar-magnetic-field-decay-timescale',
            '--pulsar-magnetic-field-decay-massscale',
            '--pulsar-minimum-magnetic-field',
            '--orbital-period-min',
            '--orbital-period-max',
            '--kick-magnitude-sigma-CCSN-NS',
            '--kick-magnitude-sigma-CCSN-BH',
            '--fix-dimensionless-kick-magnitude',
            '--kick-direction-power',
            '--random-seed',
            '--mass-transfer-thermal-limit-C',
            '--eddington-accretion-factor',
            '--PISN-lower-limit',
            '--PISN-upper-limit','--PPI-lower-limit',
            '--PPI-upper-limit',
            '--maximum-neutron-star-mass',
            '--kick-magnitude-sigma-ECSN',
            '--kick-magnitude-sigma-USSN',
            '--kick-scaling-factor',
            '--maximum-mass-donor-Nandez-Ivanova',
            '--common-envelope-recombination-energy-density',
            '--common-envelope-mass-accretion-max',
            '--common-envelope-mass-accretion-min',
            '--zeta-main-sequence',
            '--zeta-radiative-envelope-giant',
            '--kick-magnitude-max',
            '--log-level',
            '--debug-level',
            '--single-star-mass-steps',
            '--single-star-mass-min',
            '--single-star-mass-max'
        ]

        return numericalCommands

    def stringChoices(self):
        stringChoices = [
            self.case_bb_stability_prescription,
            self.chemically_homogeneous_evolution,
            self.mass_loss_prescription,
            self.mass_transfer_angular_momentum_loss_prescription,
            self.mass_transfer_accretion_efficiency_prescription,
            self.mass_transfer_rejuvenation_prescription,
            self.initial_mass_function,
            self.semi_major_axis_distribution,
            self.mass_ratio_distribution,
            self.eccentricity_distribution,
            self.rotational_velocity_distribution,
            self.remnant_mass_prescription,
            self.fryer_supernova_engine,
            self.black_hole_kicks,
            self.kick_magnitude_distribution,
            self.kick_direction,
            #self.output,
            #self.output_container,
            self.common_envelope_lambda_prescription,
            self.stellar_zeta_prescription,
            self.mass_transfer_thermal_limit_accretor,
            self.pulsational_pair_instability_prescription,
            self.neutron_star_equation_of_state,
            self.pulsar_birth_magnetic_field_distribution,
            self.pulsar_birth_spin_period_distribution,
            self.common_envelope_mass_accretion_prescription,
            self.envelope_state_prescription,
            self.logfile_name_prefix,
            self.logfile_delimiter,
            self.logfile_definitions,
            #self.grid_filename,
            self.logfile_BSE_be_binaries,
            self.logfile_BSE_common_envelopes,
            self.logfile_BSE_detailed_output,
            self.logfile_BSE_double_compact_objects,
            self.logfile_BSE_pulsar_evolution,
            self.logfile_BSE_supernovae,
            self.logfile_BSE_system_parameters
        ]

        return stringChoices

    def stringCommands(self):
        stringCommands = [
            '--case-bb-stability-prescription',
            '--chemically-homogeneous-evolution',
            '--mass-loss-prescription',
            '--mass-transfer-angular-momentum-loss-prescription',
            '--mass-transfer-accretion-efficiency-prescription',
            '--mass-transfer-rejuvenation-prescription',
            '--initial-mass-function',
            '--semi-major-axis-distribution',
            '--mass-ratio-distribution',
            '--eccentricity-distribution',
            '--rotational-velocity-distribution',
            '--remnant-mass-prescription',
            '--fryer-supernova-engine',
            '--black-hole-kicks',
            '--kick-magnitude-distribution',
            '--kick-direction',
            #'--outputPath',
            #'--output-container',
            '--common-envelope-lambda-prescription',
            '--stellar-zeta-prescription',
            '--mass-transfer-thermal-limit-accretor',
            '--pulsational-pair-instability-prescription',
            '--neutron-star-equation-of-state',
            '--pulsar-birth-magnetic-field-distribution',
            '--pulsar-birth-spin-period-distribution',
            '--common-envelope-mass-accretion-prescription',
            '--envelope-state-prescription',
            '--logfile-name-prefix',
            '--logfile-delimiter',
            '--logfile-definitions',
            #'--grid',
            '--logfile-BSE-be-binaries',
            '--logfile-BSE-common-envelopes',
            '--logfile-BSE-detailed-output',
            '--logfile-BSE-double-compact-objects',
            '--logfile-BSE-pulsar-evolution',
            '--logfile-BSE-supernovae',
            '--logfile-BSE-system-parameters'
        ]

        return stringCommands

    def listChoices(self):
        listChoices = [
            self.log_classes,
            self.debug_classes
        ]

        return listChoices

    def listCommands(self):
        listCommands = [
            '--log-classes',
            '--debug-classes'
        ]

        return listCommands


def specifyCommandLineOptions(programOptions):
    """
    This function generates a string or strings for the terminal command to run COMPAS.
    This function is intended to be modified by the user, so that they may swap out constant values for functions etc.
    Options not to be included in the command line should be set to pythons None (except booleans, which should be set to False)

    Parameters
    -----------
    programOptions : pythonProgramOptions
        Contains program options

    Returns
    --------
    commands : str or list of strs
    """
    booleanChoices = programOptions.booleanChoices()
    booleanCommands = programOptions.booleanCommands()

    numericalChoices = programOptions.numericalChoices()
    numericalCommands = programOptions.numericalCommands()

    stringChoices = programOptions.stringChoices()
    stringCommands = programOptions.stringCommands()

    listChoices = programOptions.listChoices()
    listCommands = programOptions.listCommands()

    #if programOptions.hyperparameterGrid == True:
    #    command = hyperparameterGridCommand(programOptions.compas_executable,booleanChoices,booleanCommands,numericalChoices,numericalCommands,stringChoices,stringCommands,listChoices,listCommands,programOptions.shareSeeds)
    #elif programOptions.hyperparameterList == True:
    #    if programOptions.hyperparameterGrid == True:
    #        raise ValueError("You can't have both a list and a grid!")
    #    command = hyperparameterListCommand(programOptions.compas_executable,booleanChoices,booleanCommands,numericalChoices,numericalCommands,stringChoices,stringCommands,listChoices,listCommands,programOptions.shareSeeds)
    #else:
    #    command = [generateCommandLineOptions(programOptions.compas_executable,booleanChoices,booleanCommands,numericalChoices,numericalCommands,stringChoices,stringCommands,listChoices,listCommands)]

    #command = [generateCommandLineOptions(programOptions.compas_executable,booleanChoices,booleanCommands,numericalChoices,numericalCommands,stringChoices,stringCommands,listChoices,listCommands)]
    #return command

    return generateCommandLineOptions(booleanChoices,booleanCommands,numericalChoices,numericalCommands,stringChoices,stringCommands,listChoices,listCommands)

def generateCommandLineOptions(booleanChoices,booleanCommands,numericalChoices,numericalCommands,stringChoices,stringCommands,listChoices,listCommands):

    nBoolean = len(booleanChoices)
    assert len(booleanCommands) == nBoolean

    nNumerical = len(numericalChoices)
    assert len(numericalCommands) == nNumerical

    nString = len(stringChoices)
    assert len(stringCommands) == nString

    nList = len(listChoices)
    assert len(listCommands) == nList

    # RTW
    #command = [compas_executable] 
    command = []

    for i in range(nBoolean):
        if booleanChoices[i] == True:
            command.append(booleanCommands[i])

    for i in range(nNumerical):
        if not numericalChoices[i] == None:
            command.append(numericalCommands[i])
            command.append(str(numericalChoices[i]))

    for i in range(nString):
        if not stringChoices[i] == None:
            command.append(stringCommands[i])
            command.append(stringChoices[i])

    for i in range(nList):
        if listChoices[i]:
            command.append(listCommands[i]) 
            command.append(' '.join(map(str, listChoices[i])))

    return command


##### More stroopwafel functions


def configure_code_run(batch):
    """
    This function tells stroopwafel what program to run, along with its arguments.
    IN:
        batch(dict): This is a dictionary which stores some information about one of the runs. It has an number key which stores the unique id of the run
            It also has a subprocess which will run under the key process. Rest, it depends on the user. User is free to store any information they might need later 
            for each batch run in this dictionary. For example, here I have stored the 'output_container' and 'grid_filename' so that I can read them during discovery of interesting systems below
    OUT:
        compas_args (list(String)) : This defines what will run. It should point to the executable file along with the arguments.
        Additionally one must also store the grid_filename in the batch so that the grid file is created
    """
    batch_num = batch['number']
    grid_filename = output_folder + '/grid_' + str(batch_num) + '.csv'
    output_container = 'batch_' + str(batch_num)
    # RTW
    #compas_args = [compas_executable, "--grid", grid_filename, '--outputPath', output_folder, '--logfile-delimiter', 'COMMA', '--output-container', output_container, '--random-seed', np.random.randint(2, 2**63 - 1)]
    compas_args = [compas_executable, '--grid', grid_filename, '--outputPath', output_folder,  '--output-container', output_container, ]
            #'--outputPath', output_folder, '--logfile-delimiter', 'COMMA', '--output-container', output_container, '--random-seed', np.random.randint(2, 2**63 - 1)]
    #for params in extra_params:
        #compas_args.extend(params.split("="))
    compas_args.extend(extra_params)
    batch['grid_filename'] = grid_filename
    batch['output_container'] = output_container
    return compas_args

def interesting_systems(batch):
    """
    This is a mandatory function, it tells stroopwafel what an interesting system is. User is free to define whatever looks interesting to them.
    IN:
        batch (dict): As input you will be given the current batch which just finished its execution. You can take in all the keys you defined in the configure_code_run method above
    OUT:
        Number of interesting systems
        In the below example, I define all the NSs as interesting, so I read the files, get the SEED from the system_params file and define the key is_hit in the end for all interesting systems 
    """
    # RTW hack
    return 0 #len(batch)
    try:
        folder = os.path.join(output_folder, batch['output_container'])
        shutil.move(batch['grid_filename'], folder + '/grid_' + str(batch['number']) + '.csv')
        system_parameters = pd.read_csv(folder + '/BSE_System_Parameters.csv', skiprows = 2)
        system_parameters.rename(columns = lambda x: x.strip(), inplace = True)
        seeds = system_parameters['SEED']
        for index, sample in enumerate(batch['samples']):
            seed = seeds[index]
            sample.properties['SEED'] = seed
            sample.properties['is_hit'] = 0
            sample.properties['batch'] = batch['number']
        double_compact_objects = pd.read_csv(folder + '/BSE_Double_Compact_Objects.csv', skiprows = 2)
        double_compact_objects.rename(columns = lambda x: x.strip(), inplace = True)
        #Generally, this is the line you would want to change.
        dns = double_compact_objects[np.logical_and(double_compact_objects['Merges_Hubble_Time'] == 1, \
            np.logical_and(double_compact_objects['Stellar_Type_1'] == 14, double_compact_objects['Stellar_Type_2'] == 14))]
        interesting_systems_seeds = set(dns['SEED'])
        for sample in batch['samples']:
            if sample.properties['SEED'] in interesting_systems_seeds:
                sample.properties['is_hit'] = 1
        return len(dns)
    except IOError as error:
        return 0

def selection_effects(sw):
    """
    This is not a mandatory function, it was written to support selection effects
    Fills in selection effects for each of the distributions
    IN:
        sw (Stroopwafel) : Stroopwafel object
    """
    if hasattr(sw, 'adapted_distributions'):
        biased_masses = []
        rows = []
        for distribution in sw.adapted_distributions:
            folder = os.path.join(output_folder, 'batch_' + str(int(distribution.mean.properties['batch'])))
            dco_file = pd.read_csv(folder + '/BSE_Double_Compact_Objects.csv', skiprows = 2)
            dco_file.rename(columns = lambda x: x.strip(), inplace = True)
            row = dco_file.loc[dco_file['SEED'] == distribution.mean.properties['SEED']]
            rows.append([row.iloc[0]['Mass_1'], row.iloc[0]['Mass_2']])
            biased_masses.append(np.power(max(rows[-1]), 2.2))
        # update the weights
        mean = np.mean(biased_masses)
        for index, distribution in enumerate(sw.adapted_distributions):
            distribution.biased_weight = np.power(max(rows[index]), 2.2) / mean

def rejected_systems(locations, dimensions):
    """
    This method takes a list of locations and marks the systems which can be
    rejected by the prior distribution
    IN:
        locations (List(Location)): list of location to inspect and mark rejected
    OUT:
        num_rejected (int): number of systems which can be rejected
    """
    m1 = dimensions[0]
    q = dimensions[1]
    a = dimensions[2]
    mass_1 = [location.dimensions[m1] for location in locations]
    mass_2 = [location.properties['Mass_2'] for location in locations]
    metallicity_1 = [location.properties['Metallicity_1'] for location in locations]
    metallicity_2 = [location.properties['Metallicity_2'] for location in locations]
    eccentricity = [location.properties['Eccentricity'] for location in locations]
    num_rejected = 0
    for index, location in enumerate(locations):
        radius_1 = utils.get_zams_radius(mass_1[index], metallicity_1[index])
        radius_2 = utils.get_zams_radius(mass_2[index], metallicity_2[index])
        roche_lobe_tracker_1 = radius_1 / (location.dimensions[a] * (1 - eccentricity[index]) * utils.calculate_roche_lobe_radius(mass_1[index], mass_2[index]))
        roche_lobe_tracker_2 = radius_2 / (location.dimensions[a] * (1 - eccentricity[index]) * utils.calculate_roche_lobe_radius(mass_2[index], mass_1[index]))
        location.properties['is_rejected'] = 0
        if (mass_2[index] < constants.MINIMUM_SECONDARY_MASS) or (location.dimensions[a] <= (radius_1 + radius_2)) \
        or roche_lobe_tracker_1 > 1 or roche_lobe_tracker_2 > 1:
            location.properties['is_rejected'] = 1
            num_rejected += 1
    return num_rejected


programOptions = pythonProgramOptions()
extra_params = specifyCommandLineOptions(programOptions) # overwrite any further options

if __name__ == '__main__':
    start_time = time.time()

    # Check if we are using python 3
    python_version = sys.version_info[0]
    print("python_version =", python_version)


    #Define the parameters to the constructor of stroopwafel
    TOTAL_NUM_SYSTEMS = namespace.num_systems #total number of systems you want in the end
    NUM_CPU_CORES = namespace.num_cores #Number of cpu cores you want to run in parellel
    NUM_SYSTEMS_PER_RUN = namespace.num_per_core #Number of systems generated by each of run on each cpu core
    debug = namespace.debug #If True, will print the logs given by the external program (like COMPAS)
    run_on_helios = namespace.run_on_helios #If True, it will run on a clustered system helios, rather than your pc
    mc_only = namespace.mc_only # If you dont want to do the refinement phase and just do random mc exploration
    output_filename = namespace.output_filename #The name of the output file
    compas_executable = os.path.join(os.environ.get('COMPAS_ROOT_DIR'), 'src/COMPAS') # Location of the executable
    output_folder =  os.path.join(os.getcwd(), 'output') # Folder you want to receieve outputs, here the current working directory, but you can specify anywhere

    if os.path.exists(output_folder):
        command = input ("The output folder already exists. If you continue, I will remove all its content. Press (Y/N)\n")
        if (command == 'Y'):
            shutil.rmtree(output_folder)
        else:
            exit()
    os.makedirs(output_folder)

    # STEP 2 : Create an instance of the Stroopwafel class
    sw_object = sw.Stroopwafel(TOTAL_NUM_SYSTEMS, NUM_CPU_CORES, NUM_SYSTEMS_PER_RUN, output_folder, output_filename, debug = debug, run_on_helios = run_on_helios, mc_only = mc_only)


    #STEP 3: Initialize the stroopwafel object with the user defined functions and create dimensions and initial distribution
    dimensions = create_dimensions()
    sw_object.initialize(dimensions, interesting_systems, configure_code_run, rejected_systems, update_properties_method = update_properties)

    intial_pdf = distributions.InitialDistribution(dimensions)
    #STEP 4: Run the 4 phases of stroopwafel
    sw_object.explore(intial_pdf) #Pass in the initial distribution for exploration phase
    #sw_object.adapt(n_dimensional_distribution_type = distributions.Gaussian) #Adaptaion phase, tell stroopwafel what kind of distribution you would like to create instrumental distributions
    ## Do selection effects
    #selection_effects(sw)
    #sw_object.refine() #Stroopwafel will draw samples from the adapted distributions
    #sw_object.postprocess(distributions.Gaussian, only_hits = False) #Run it to create weights, if you want only hits in the output, then make only_hits = True

    end_time = time.time()
    print ("Total running time = %d seconds" %(end_time - start_time))
