import numpy as np
import subprocess
import sys
import os
import pickle
import itertools 
from subprocess import call

# Check if we are using python 3
python_version = sys.version_info[0]
print("python_version =", python_version)

class pythonProgramOptions:
    """
    A class to store and access COMPAS program options in python
    """

    # Do './COMPAS --help' to see all options
    #-- Define variables
    compas_executable_override = os.environ.get('COMPAS_EXECUTABLE_PATH')
    print('compas_executable_override', compas_executable_override)
    
    if (compas_executable_override is None):
        git_directory = os.environ.get('COMPAS_ROOT_DIR')
        compas_executable = os.path.join(git_directory, 'src/COMPAS')
    else:
        compas_executable = compas_executable_override

    enable_warnings = False                                     # option to enable/disable warning messages

    number_of_binaries = 10  #number of binaries per batch
    populationPrinting = False

    randomSeedFileName = 'randomSeed.txt'
    if os.path.isfile(randomSeedFileName):
        random_seed = int(np.loadtxt(randomSeedFileName))
    else:
        random_seed = 0 # If you want a random seed, use: np.random.randint(2,2**63-1)

    compas_logs_output_override = os.environ.get('COMPAS_LOGS_OUTPUT_DIR_PATH')
    
    if (compas_logs_output_override is None):
        output = os.getcwd()
        output_container = None                 # names the directory to be created and in which log files are created.  Default in COMPAS is "COMPAS_Output"
    else:
        output = compas_logs_output_override
        output_container = None

    #-- option to make a grid of hyperparameter values at which to produce populations.
    #-- If this is set to true, it will divide the number_of_binaries parameter equally
    #-- amoungst the grid points (as closely as possible). See the hyperparameterGrid method below
    #-- for more details. If this is set to True, some hyperparameter values defined in this method'gridOutputs/'+str(i)
    #-- will be overwritten
    hyperparameterGrid = False
    hyperparameterList = False
    shareSeeds = False

    single_star = False
    single_star_mass_steps = 10
    single_star_mass_min   = 1.0
    single_star_mass_max   = 75.0    

    grid_filename = None

    use_mass_loss = True
    mass_transfer = True
    detailed_output = False                             # WARNING: this creates a data heavy file
    RLOFPrinting = True
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
    kick_velocity_distribution = 'MAXWELLIAN'

    kick_velocity_sigma_CCSN_NS = 265.0                 #  [km/s]
    kick_velocity_sigma_CCSN_BH = 265.0                 #  [km/s]
    kick_velocity_sigma_ECSN = 30.0                     #  [km/s]
    kick_velocity_sigma_USSN = 30.0                     #  [km/s]

    fix_dimensionless_kick_velocity = -1
    kick_direction = 'ISOTROPIC'
    kick_direction_power = 0.0
    kick_scaling_factor = 1.0
    kick_velocity_maximum = -1.0

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
    logfile_BSE_rlof_parameters = None
    logfile_BSE_pulsar_evolution = None
    logfile_BSE_supernovae = None
    logfile_BSE_system_parameters = None

    debug_to_file  = False
    errors_to_file = False

    def booleanChoices(self):
        booleanChoices = [
            self.enable_warnings,
            self.single_star,
            self.use_mass_loss,
            self.mass_transfer,
            self.detailed_output,
            self.evolve_unbound_systems,
            self.populationPrinting,
            self.RLOFPrinting,
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
            '--RLOFPrinting',
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
            self.number_of_binaries,
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
            self.kick_velocity_sigma_CCSN_NS,
            self.kick_velocity_sigma_CCSN_BH,
            self.fix_dimensionless_kick_velocity,
            self.kick_direction_power,
            self.random_seed,
            self.mass_transfer_thermal_limit_C,
            self.eddington_accretion_factor,
            self.PISN_lower_limit,
            self.PISN_upper_limit,
            self.PPI_lower_limit,
            self.PPI_upper_limit,
            self.maximum_neutron_star_mass,
            self.kick_velocity_sigma_ECSN,
            self.kick_velocity_sigma_USSN,
            self.kick_scaling_factor,
            self.common_envelope_maximum_donor_mass_revised_energy_formalism,
            self.common_envelope_recombination_energy_density,
            self.common_envelope_mass_accretion_max,
            self.common_envelope_mass_accretion_min,
            self.zeta_Main_Sequence,
            self.zeta_Radiative_Envelope_Giant,
            self.kick_velocity_maximum,
            self.log_level,
            self.debug_level,
            self.single_star_mass_steps,
            self.single_star_mass_min,
            self.single_star_mass_max
        ]

        return numericalChoices

    def numericalCommands(self):
        numericalCommands = [
            '--number-of-binaries',
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
            '--kick-velocity-sigma-CCSN-NS',
            '--kick-velocity-sigma-CCSN-BH',
            '--fix-dimensionless-kick-velocity',
            '--kick-direction-power',
            '--random-seed',
            '--mass-transfer-thermal-limit-C',
            '--eddington-accretion-factor',
            '--PISN-lower-limit',
            '--PISN-upper-limit','--PPI-lower-limit',
            '--PPI-upper-limit',
            '--maximum-neutron-star-mass',
            '--kick-velocity-sigma-ECSN',
            '--kick-velocity-sigma-USSN',
            '--kick-scaling-factor',
            '--maximum-mass-donor-Nandez-Ivanova',
            '--common-envelope-recombination-energy-density',
            '--common-envelope-mass-accretion-max',
            '--common-envelope-mass-accretion-min',
            '--zeta-main-sequence',
            '--zeta-radiative-envelope-giant',
            '--kick-velocity-max',
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
            self.kick_velocity_distribution,
            self.kick_direction,
            self.output,
            self.output_container,
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
            self.grid_filename,
            self.logfile_BSE_be_binaries,
            self.logfile_BSE_common_envelopes,
            self.logfile_BSE_detailed_output,
            self.logfile_BSE_double_compact_objects,
            self.logfile_BSE_pulsar_evolution,
            self.logfile_BSE_rlof_parameters,
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
            '--kick-velocity-distribution',
            '--kick-direction',
            '--outputPath',
            '--output-container',
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
            '--grid',
            '--logfile-BSE-be-binaries',
            '--logfile-BSE-common-envelopes',
            '--logfile-BSE-detailed-output',
            '--logfile-BSE-double-compact-objects',
            '--logfile-BSE-pulsar-evolution',
            '--logfile-BSE-rlof-parameters',
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

    if programOptions.hyperparameterGrid == True:
        command = hyperparameterGridCommand(programOptions.compas_executable,booleanChoices,booleanCommands,numericalChoices,numericalCommands,stringChoices,stringCommands,listChoices,listCommands,programOptions.shareSeeds)
    elif programOptions.hyperparameterList == True:
        if programOptions.hyperparameterGrid == True:
            raise ValueError("You can't have both a list and a grid!")
        command = hyperparameterListCommand(programOptions.compas_executable,booleanChoices,booleanCommands,numericalChoices,numericalCommands,stringChoices,stringCommands,listChoices,listCommands,programOptions.shareSeeds)
    else:
        command = [generateCommandLineOptions(programOptions.compas_executable,booleanChoices,booleanCommands,numericalChoices,numericalCommands,stringChoices,stringCommands,listChoices,listCommands)]

    #command = [generateCommandLineOptions(programOptions.compas_executable,booleanChoices,booleanCommands,numericalChoices,numericalCommands,stringChoices,stringCommands,listChoices,listCommands)]

    return command

def generateCommandLineOptions(compas_executable,booleanChoices,booleanCommands,numericalChoices,numericalCommands,stringChoices,stringCommands,listChoices,listCommands):

    nBoolean = len(booleanChoices)
    assert len(booleanCommands) == nBoolean

    nNumerical = len(numericalChoices)
    assert len(numericalCommands) == nNumerical

    nString = len(stringChoices)
    assert len(stringCommands) == nString

    nList = len(listChoices)
    assert len(listCommands) == nList

    command = compas_executable + ' '

    for i in range(nBoolean):

        if booleanChoices[i] == True:

            command += booleanCommands[i] + ' '

    for i in range(nNumerical):

        if not numericalChoices[i] == None:

            command += numericalCommands[i] + ' ' + str(numericalChoices[i]) + ' '

    for i in range(nString):

        if not stringChoices[i] == None:

            command += stringCommands[i] + ' ' + stringChoices[i] + ' '

    for i in range(nList):

        if listChoices[i]:

            command += listCommands[i] + ' ' + ' '.join(map(str, listChoices[i]))

    return command

def hyperparameterGridCommand(compas_executable,booleanChoices,booleanCommands,numericalChoices,numericalCommands,stringChoices,stringCommands,listChoices,listCommands,shareSeeds):
    """This function allows for a range of hyperparameter values to be specified in a single run, if the hyperparameterGrid boolean is set to True in the
    specifyCommandLineOptions() function.
    This works by constructing nested output directories in the current working directory, and running a population at each combination of parameter values.
    nBinaries from specifyCommandLineOptions() is divided equally amoungst these.
    The user should follow the pattern in adding items to the commandsAndValues dictionary, the code will then handle production of
    the output folders and return a command line command to run all of the populations back to back
    """
    
    # Load up the dictionary from gridRun.py
    with open('pickledGrid.pkl', 'rb') as pg:
        commandsAndValues = pickle.load(pg)
    
    # set up lists for recursion
    if python_version >= 3:
        keys = list(commandsAndValues.keys())
    else:
        keys = commandsAndValues.keys()
    valuesLists = []
    nSimulations = 1
    for key in keys:
        nSimulations *= len(commandsAndValues[key])
        valuesLists.append(commandsAndValues[key])
    # Make folders for the messy outputs
    outPaths = []
    for i in range(nSimulations):
        path = 'gridOutputs/output-'+str(i)
        outPaths.append(path)
    #edit number of binaries per population
    for index,command in enumerate(numericalCommands):
        if command == '--number-of-binaries':
            break
    nBinariesPerSimulation = numericalChoices[index]/nSimulations

    numericalChoices[index] = int(nBinariesPerSimulation)
    print("index, nBinariesPerSimulation")
    print(index, nBinariesPerSimulation)

    bashCommands = []
    # itertools.product recurses through all combinations of the lists in valuesLists
    for en,combination in enumerate(itertools.product(*valuesLists)):
        bashCommand = ''
        pathName = outPaths[en]
        for i, val in enumerate(combination):
            for index,command in enumerate(numericalCommands):
                if command == keys[i]:
                    break
            numericalChoices[index] = val
        #change the random seed if need be
        if not shareSeeds:
            for index,command in enumerate(numericalCommands):
                if command == '--random-seed':
                    break
            numericalChoices[index] += int(nBinariesPerSimulation)
        #setup output arguments
        for index,command in enumerate(stringCommands):
            if command == '--outputPath':
                break
        stringChoices[index] = pathName + '/.'
        bashCommand += generateCommandLineOptions(compas_executable,booleanChoices,booleanCommands,numericalChoices,numericalCommands,stringChoices,stringCommands,listChoices,listCommands)
        bashCommand += '; '
        bashCommands.append(bashCommand)
    return bashCommands
    

    
def hyperparameterListCommand(compas_executable,booleanChoices,booleanCommands,numericalChoices,numericalCommands,stringChoices,stringCommands,listChoices,listCommands,shareSeeds):
    """
    """
    # Load up the dictionary from gridRun.py
    with open('pickledList.pkl', 'rb') as pl:
        commandsAndValues = pickle.load(pl)
    
    # set up lists for recursion
    if python_version >= 3:
        keys = list(commandsAndValues.keys())
    else:
        keys = commandsAndValues.keys()

    #work how many things there are in the list
    nSimulations = len(commandsAndValues[keys[0]])
    print("nSimulations = ", nSimulations)
    #make the directories
    outPaths = []
    for i in range(nSimulations):
        path = 'listOutputs/output-'+str(i)
        outPaths.append(path)
    #grab all the values out of the dictionary for syntactic ease later
    valuesLists = []
    for key in keys:
        valuesLists.append(commandsAndValues[key])
    valuesLists =np.array(valuesLists).T
    #edit number of binaries per population
    for index,command in enumerate(numericalCommands):
        if command == '--number-of-binaries':
            break
    nBinariesPerSimulation = numericalChoices[index]/nSimulations
    numericalChoices[index] = int(nBinariesPerSimulation)
    print("index, nBinariesPerSimulation")
    print(index, nBinariesPerSimulation)
    bashCommands = []
    for en in range(nSimulations):
        bashCommand = ''
        combination = valuesLists[en]
        pathName = outPaths[en]
        for i, val in enumerate(combination):
            for index,command in enumerate(numericalCommands):
                if command == keys[i]:
                    break
            numericalChoices[index] = val
        #change the random seed if need be
        if not shareSeeds:
            for index,command in enumerate(numericalCommands):
                if command == '--random-seed':
                    break
            numericalChoices[index] += int(nBinariesPerSimulation)
        #setup output arguments
        for index,command in enumerate(stringCommands):
            if command == '--outputPath':
                break
        stringChoices[index] = pathName + '/.'
        bashCommand += generateCommandLineOptions(compas_executable,booleanChoices,booleanCommands,numericalChoices,numericalCommands,stringChoices,stringCommands,listChoices,listCommands)
        bashCommand += '; '
        bashCommands.append(bashCommand)
    return bashCommands
    

def runCompas(programOptions):
    """
    """

    commands = specifyCommandLineOptions(programOptions)

    for command in commands:

        print(command)

        call(command,shell=True)

    return 0


if __name__ == "__main__":

    #-- Get the program options
    programOptions = pythonProgramOptions()

    runCompas(programOptions)

