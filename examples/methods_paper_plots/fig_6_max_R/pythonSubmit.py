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

    # environment variable COMPAS_EXECUTABLE_PATH is used for docker runs
    # if COMPAS_EXECUTABLE_PATH is not set (== None) we assume this is an
    # interactive run with python3
    # if COMPAS_EXECUTABLE_PATH is set (!= None) we assume this is a run
    # inside a docker container - we have different directories inside a 
    # docker container (src, obj, bin), and the COMPAS executable resides
    # in the bin directory (rather than the src directory)
    compas_executable_override = os.environ.get('COMPAS_EXECUTABLE_PATH')
    print('compas_executable_override', compas_executable_override)

    if (compas_executable_override is None):
        git_directory = os.environ.get('COMPAS_ROOT_DIR')
        compas_executable = os.path.join(git_directory, 'src/COMPAS')
    else:
        compas_executable = compas_executable_override

    enable_warnings = False                                     # option to enable/disable warning messages

    number_of_systems = 10  #number of systems per batch

    populationPrinting = False

    randomSeedFileName = 'randomSeed.txt'
    if os.path.isfile(randomSeedFileName):
        random_seed = int(np.loadtxt(randomSeedFileName))
    else:
        random_seed = 0 # If you want a random seed, use: np.random.randint(2,2**63-1)

    # environment variable COMPAS_LOGS_OUTPUT_DIR_PATH is used primarily for docker runs
    # if COMPAS_LOGS_OUTPUT_DIR_PATH is set (!= None) it is used as the value for the
    # --output-path option
    # if COMPAS_LOGS_OUTPUT_DIR_PATH is not set (== None) the current working directory
    # is used as the value for the --output-path option
    compas_logs_output_override = os.environ.get('COMPAS_LOGS_OUTPUT_DIR_PATH')
    if (compas_logs_output_override is None):
        output = os.getcwd()
        output_container = None                 # names the directory to be created and in which log files are created.  Default in COMPAS is "COMPAS_Output"
    else:
        output = compas_logs_output_override
        output_container = None

    # environment variable COMPAS_INPUT_DIR_PATH is used primarily for docker runs
    # if COMPAS_INPUT_DIR_PATH is set (!= None) it is prepended to input filenames
    # (such as grid_filename and logfile_definitions)
    # if COMPAS_INPUT_DIR_PATH is not set (== None) the current working directory
    # is prepended to input filenames
    compas_input_path_override = os.environ.get('COMPAS_INPUT_DIR_PATH')

    #-- option to make a grid of hyperparameter values at which to produce populations.
    #-- If this is set to true, it will divide the number_of_binaries parameter equally
    #-- amoungst the grid points (as closely as possible). See the hyperparameterGrid method below
    #-- for more details. If this is set to True, some hyperparameter values defined in this method'gridOutputs/'+str(i)
    #-- will be overwritten
    hyperparameterGrid = False
    hyperparameterList = False
    shareSeeds = False

    mode = 'BSE'                                                # evolving single stars (SSE) or binaries (BSE)?

    grid_filename = "grid.txt"                                        # grid file name (e.g. 'mygrid.txt')

    if grid_filename != None:
        if compas_input_path_override == None:
            grid_filename = os.getcwd() + '/' + grid_filename
        else:
            grid_filename = compas_input_path_override + '/' + grid_filename

    logfile_definitions = None                                  # logfile record definitions file name (e.g. 'logdefs.txt')

    if logfile_definitions != None:
        if compas_input_path_override == None:
            logfile_definitions = os.getcwd() + '/' + logfile_definitions
        else:
            logfile_definitions = compas_input_path_override + '/' + logfile_definitions

    initial_mass    = None                                      # initial mass for SSE
    initial_mass_1  = None                                      # primary initial mass for BSE
    initial_mass_2  = None                                      # secondary initial mass for BSE

    mass_ratio      = None

    eccentricity    = None                                      # eccentricity for BSE
    semi_major_axis = None                                      # semi-major axis for BSE
    orbital_period  = None                                      # orbital period for BSE


    use_mass_loss = True
    mass_transfer = True
    detailed_output = True                                     # WARNING: this creates a data heavy file
    RLOFPrinting = True
    evolve_unbound_systems = False
    quiet = False

    metallicity = None                                       # metallicity for both SSE and BSE - Solar metallicity Asplund+2010

    allow_rlof_at_birth = True                                  # allow binaries that have one or both stars in RLOF at birth to evolve?
    allow_touching_at_birth = False                             # record binaries that have stars touching at birth in output files?

    chemically_homogeneous_evolution = 'PESSIMISTIC'            # chemically homogeneous evolution.  Options are 'NONE', 'OPTIMISTIC' and 'PESSIMISTIC'

    switch_log = False

    common_envelope_alpha = 1.0
    common_envelope_lambda = 0.1                                # Only if using 'LAMBDA_FIXED'
    common_envelope_lambda_prescription = 'LAMBDA_NANJING'      # Xu & Li 2010
    common_envelope_slope_Kruckow = -5.0/6.0
    stellar_zeta_prescription = 'SOBERMAN'
    common_envelope_revised_energy_formalism = False
    common_envelope_maximum_donor_mass_revised_energy_formalism = 2.0
    common_envelope_recombination_energy_density = 1.5E13
    common_envelope_alpha_thermal = 1.0                         # lambda = alpha_th*lambda_b + (1-alpha_th)*lambda_g
    common_envelope_lambda_multiplier = 1.0                     # Multiply common envelope lambda by some constant
    common_envelope_allow_main_sequence_survive = True          # Allow main sequence stars to survive CE. Was previously False by default
    common_envelope_mass_accretion_prescription = 'ZERO'
    common_envelope_mass_accretion_min = 0.04                   # For 'MACLEOD+2014' [Msol]
    common_envelope_mass_accretion_max = 0.10                   # For 'MACLEOD+2014' [Msol]
    envelope_state_prescription = 'LEGACY'

    mass_loss_prescription = 'VINK'
    luminous_blue_variable_prescription = 'HURLEY_ADD'
    luminous_blue_variable_multiplier = 1.5
    overall_wind_mass_loss_multiplier = 1.0
    wolf_rayet_multiplier = 1.0
    cool_wind_mass_loss_multiplier = 1.0
    check_photon_tiring_limit = False

    circularise_binary_during_mass_transfer = False
    angular_momentum_conservation_during_circularisation = False
    mass_transfer_angular_momentum_loss_prescription = 'ISOTROPIC'
    mass_transfer_accretion_efficiency_prescription = 'THERMAL'
    mass_transfer_fa = 0.5                                      # Only if using mass_transfer_accretion_efficiency_prescription = 'FIXED'
    mass_transfer_jloss = 1.0                                   # Only if using mass_transfer_angular_momentum_loss_prescription = 'FIXED'
    mass_transfer_rejuvenation_prescription = 'STARTRACK'
    mass_transfer_thermal_limit_accretor= 'CFACTOR'
    mass_transfer_thermal_limit_C= 10.0
    eddington_accretion_factor = 1                              # multiplication Factor for eddington accretion onto NS&BH

    case_BB_stability_prescription = 'ALWAYS_STABLE'
    zeta_Main_Sequence = 2.0
    zeta_Radiative_Envelope_Giant = 6.5

    maximum_evolution_time = 13700.0                            # Maximum physical time a system can be evolved [Myrs]
    maximum_number_timesteps = 99999
    timestep_multiplier = 0.3                                   # Optional multiplier relative to default time step duration

    initial_mass_function = 'KROUPA'
    initial_mass_min = 5.0                                      # Use 1.0 for LRNe, 5.0 for DCOs  [Msol]
    initial_mass_max = 150.0                                    # Stellar tracks extrapolated above 50 Msol (Hurley+2000) [Msol]

    initial_mass_power = 0.0

    semi_major_axis_distribution = 'FLATINLOG'
    semi_major_axis_min = 0.01                                  # [AU]
    semi_major_axis_max = 1000.0                                # [AU]

    orbital_period_distribution = 'FLATINLOG'
    orbital_period_min = 1.1                                    # [days]
    orbital_period_max = 1000                                   # [days]

    mass_ratio_distribution = 'FLAT'
    mass_ratio_min = 0.01
    mass_ratio_max = 1.0

    minimum_secondary_mass = 0.1                                # Brown dwarf limit  [Msol]

    eccentricity_distribution = 'ZERO'
    eccentricity_min = 0.0
    eccentricity_max = 1.0

    metallicity_distribution = 'ZSOLAR'
    metallicity_min = 0.0001
    metallicity_max = 0.03

    pulsar_birth_magnetic_field_distribution = 'ZERO'
    pulsar_birth_magnetic_field_min = 11.0                      # [log10(B/G)]
    pulsar_birth_magnetic_field_max = 13.0                      # [log10(B/G)]

    pulsar_birth_spin_period_distribution = "ZERO"
    pulsar_birth_spin_period_min = 10.0                         # [ms]
    pulsar_birth_spin_period_max = 100.0                        # [ms]

    pulsar_magnetic_field_decay_timescale = 1000.0              # [Myr]
    pulsar_magnetic_field_decay_massscale = 0.025               # [Msol]
    pulsar_minimum_magnetic_field = 8.0                         # [log10(B/G)]

    evolvePulsars = False

    rotational_velocity_distribution = 'ZERO'

    neutron_star_equation_of_state = 'SSE'

    neutrino_mass_loss_BH_formation = "FIXED_MASS"              # "FIXED_FRACTION"
    neutrino_mass_loss_BH_formation_value = 0.1                 # Either fraction or mass (Msol) to lose
    
    remnant_mass_prescription   = 'FRYER2012'                   #
    fryer_supernova_engine      = 'DELAYED'
    black_hole_kicks            = 'FALLBACK'
    kick_magnitude_distribution = 'MAXWELLIAN'

    kick_magnitude_sigma_CCSN_NS = 265.0                        #  [km/s]
    kick_magnitude_sigma_CCSN_BH = 265.0                        #  [km/s]
    kick_magnitude_sigma_ECSN    = 30.0                         #  [km/s]
    kick_magnitude_sigma_USSN    = 30.0                         #  [km/s]

    fix_dimensionless_kick_magnitude = -1
    kick_direction = 'ISOTROPIC'
    kick_direction_power = 0.0
    kick_scaling_factor = 1.0
    kick_magnitude_maximum = -1.0

    kick_magnitude_random   = None                              # (SSE) used to draw the kick magnitude for the star should it undergo a supernova event
    kick_magnitude          = None                              # (SSE) (drawn) kick magnitude for the star should it undergo a supernova event [km/s]

    kick_magnitude_random_1 = None                              # (BSE) used to draw the kick magnitude for the primary star should it undergo a supernova event
    kick_magnitude_1        = None                              # (BSE) (drawn) kick magnitude for the primary star should it undergo a supernova event [km/s]
    kick_theta_1            = None                              # (BSE) angle between the orbital plane and the 'z' axis of the supernova vector for the primary star should it undergo a supernova event [radians]
    kick_phi_1              = None                              # (BSE) angle between 'x' and 'y', both in the orbital plane of the supernova vector, for the primary star should it undergo a supernova event [radians]
    kick_mean_anomaly_1     = None                              # (BSE) mean anomaly at the instant of the supernova for the primary star should it undergo a supernova event - should be uniform in [0, 2pi) [radians]

    kick_magnitude_random_2 = None                              # (BSE) used to draw the kick velocity for the secondary star should it undergo a supernova event
    kick_magnitude_2        = None                              # (BSE) (drawn) kick magnitude for the secondary star should it undergo a supernova event [km/s]
    kick_theta_2            = None                              # (BSE) angle between the orbital plane and the 'z' axis of the supernova vector for the secondary star should it undergo a supernova event [radians]
    kick_phi_2              = None                              # (BSE) angle between 'x' and 'y', both in the orbital plane of the supernova vector, for the secondary star should it undergo a supernova event [radians]
    kick_mean_anomaly_2     = None                              # (BSE) mean anomaly at the instant of the supernova for the secondary star should it undergo a supernova event - should be uniform in [0, 2pi) [radians]

    muller_mandel_kick_multiplier_BH = 200.0                    # scaling prefactor for BH kicks when using the 'MULLERMANDEL' kick magnitude distribution
    muller_mandel_kick_multiplier_NS = 400.0                    # scaling prefactor for NS kicks when using the 'MULLERMANDEL' kick magnitude distribution

    pair_instability_supernovae = True
    PISN_lower_limit = 60.0                                     # Minimum core mass for PISN [Msol]
    PISN_upper_limit = 135.0                                    # Maximum core mass for PISN [Msol]

    pulsation_pair_instability = True
    PPI_lower_limit = 35.0                                      # Minimum core mass for PPI [Msol]
    PPI_upper_limit = 60.0                                      # Maximum core mass for PPI [Msol]

    pulsational_pair_instability_prescription = 'MARCHANT'

    maximum_neutron_star_mass = 2.5  #  [Msol]

    log_level           = 0
    log_classes         = []

    debug_level         = 0
    debug_classes       = []

    logfile_name_prefix = None
    logfile_type        = 'HDF5'

    hdf5_chunk_size     = 100000
    hdf5_buffer_size    = 1


    # set the logfile names here
    #
    # set to None (e.g. logfile_BSE_supernovae = None) to use the default filename
    # set to a string (e.g. logfile_BSE_supernovae = 'mySNfilename') to use that string as the filename 
    # set to empty string (e.g. logfile_BSE_supernovae = '""') to disable logging for that file (the file will not be created)
    #
    # We don't really need the 'BSE' or 'SSE' prefixes any more - they were put there because
    # prior to the implementation of the containing folder it was too hard to locate the files
    # created by a COMPAS run - especially the detailed output files.  Now that the output
    # files are created inside a containing folder for each run there is really no need for
    # the prefixes - and if we don't have the prefixes we can share some of the options
    # (e.g. specifying the supernovae filename doesn't need to have separate options for 
    # SSE and BSE - we really just need one (we only ever run in one mode or the other))
    #
    # For now though, I'll leave them as is - we can change this when (if) we decide to
    # drop the prefixes

    logfile_common_envelopes       = None
    logfile_detailed_output        = None
    logfile_double_compact_objects = None
    logfile_rlof_parameters        = None
    logfile_pulsar_evolution       = None
    logfile_supernovae             = None
    logfile_switch_log             = None
    logfile_system_parameters      = None

    debug_to_file  = False
    errors_to_file = False


    def booleanChoices(self):
        booleanChoices = [
            self.enable_warnings,
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
            self.allow_touching_at_birth,
            self.switch_log,
            self.check_photon_tiring_limit
        ]

        return booleanChoices

    def booleanCommands(self):
        booleanCommands = [
            '--enable-warnings',
            '--use-mass-loss',
            '--mass-transfer',
            '--detailed-output',
            '--evolve-unbound-systems',
            '--population-data-printing',
            '--rlof-printing',
            '--circularise-binary-during-mass-transfer',
            '--angular-momentum-conservation-during-circularisation',
            '--pair-instability-supernovae',
            '--pulsational-pair-instability',
            '--quiet',
            '--common-envelope-allow-main-sequence-survive',
            '--evolve-pulsars',
            '--debug-to-file',
            '--errors-to-file',
            '--allow-rlof-at-birth',
            '--allow-touching-at-birth',
            '--switch-log',
            '--check-photon-tiring-limit'
        ]

        return booleanCommands

    def numericalChoices(self):
        numericalChoices = [
            self.number_of_systems,
            self.initial_mass,
            self.initial_mass_1,
            self.initial_mass_2,
            self.eccentricity,
            self.semi_major_axis,
            self.orbital_period,
            self.metallicity,
            self.common_envelope_alpha,
            self.common_envelope_lambda,
            self.common_envelope_slope_Kruckow,
            self.common_envelope_alpha_thermal,
            self.common_envelope_lambda_multiplier,
            self.luminous_blue_variable_multiplier,
            self.overall_wind_mass_loss_multiplier,
            self.wolf_rayet_multiplier,
            self.cool_wind_mass_loss_multiplier,
            self.mass_transfer_fa,
            self.mass_transfer_jloss,
            self.maximum_evolution_time,
            self.maximum_number_timesteps,
            self.timestep_multiplier,
            self.initial_mass_min,
            self.initial_mass_max,
            self.initial_mass_power,
            self.semi_major_axis_min,
            self.semi_major_axis_max,
            self.mass_ratio,
            self.mass_ratio_min,
            self.mass_ratio_max,
            self.minimum_secondary_mass,
            self.eccentricity_min,
            self.eccentricity_max,
            self.metallicity_min,
            self.metallicity_max,
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
            self.kick_magnitude_random,
            self.kick_magnitude,
            self.kick_magnitude_random_1,
            self.kick_magnitude_1,
            self.kick_theta_1,
            self.kick_phi_1,
            self.kick_mean_anomaly_1,
            self.kick_magnitude_random_2,
            self.kick_magnitude_2,
            self.kick_theta_2,
            self.kick_phi_2,
            self.kick_mean_anomaly_2,
            self.muller_mandel_kick_multiplier_BH,
            self.muller_mandel_kick_multiplier_NS,
            self.log_level,
            self.debug_level,
            self.hdf5_chunk_size,
            self.hdf5_buffer_size,
            self.neutrino_mass_loss_BH_formation_value
        ]

        return numericalChoices

    def numericalCommands(self):
        numericalCommands = [
            '--number-of-systems',
            '--initial-mass',
            '--initial-mass-1',
            '--initial-mass-2',
            '--eccentricity',
            '--semi-major-axis',
            '--orbital-period',
            '--metallicity',
            '--common-envelope-alpha',
            '--common-envelope-lambda',
            '--common-envelope-slope-kruckow',
            '--common-envelope-alpha-thermal',
            '--common-envelope-lambda-multiplier',
            '--luminous-blue-variable-multiplier',
            '--overall-wind-mass-loss-multiplier',
            '--wolf-rayet-multiplier',
            '--cool-wind-mass-loss-multiplier',
            '--mass-transfer-fa',
            '--mass-transfer-jloss',
            '--maximum-evolution-time',
            '--maximum-number-timestep-iterations',
            '--timestep-multiplier',
            '--initial-mass-min',
            '--initial-mass-max',
            '--initial-mass-power',
            '--semi-major-axis-min',
            '--semi-major-axis-max',
            '--mass-ratio',
            '--mass-ratio-min',
            '--mass-ratio-max',
            '--minimum-secondary-mass',
            '--eccentricity-min',
            '--eccentricity-max',
            '--metallicity-min',
            '--metallicity-max',
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
            '--pisn-lower-limit',
            '--pisn-upper-limit',
            '--ppi-lower-limit',
            '--ppi-upper-limit',
            '--maximum-neutron-star-mass',
            '--kick-magnitude-sigma-ECSN',
            '--kick-magnitude-sigma-USSN',
            '--kick-scaling-factor',
            '--maximum-mass-donor-nandez-ivanova',
            '--common-envelope-recombination-energy-density',
            '--common-envelope-mass-accretion-max',
            '--common-envelope-mass-accretion-min',
            '--zeta-main-sequence',
            '--zeta-radiative-envelope-giant',
            '--kick-magnitude-max',
            '--kick-magnitude-random',
            '--kick-magnitude',
            '--kick-magnitude-random-1',
            '--kick-magnitude-1',
            '--kick-theta-1',
            '--kick-phi-1',
            '--kick-mean-anomaly-1',
            '--kick-magnitude-random-2',
            '--kick-magnitude-2',
            '--kick-theta-2',
            '--kick-phi-2',
            '--kick-mean-anomaly-2',
            '--muller-mandel-kick-multiplier-BH',
            '--muller-mandel-kick-multiplier-NS',
            '--log-level',
            '--debug-level',
            '--hdf5-chunk-size',
            '--hdf5-buffer-size',
            '--neutrino-mass-loss-BH-formation-value'
        ]

        return numericalCommands

    def stringChoices(self):
        stringChoices = [
            self.mode,
            self.case_BB_stability_prescription,
            self.chemically_homogeneous_evolution,
            self.luminous_blue_variable_prescription,
            self.mass_loss_prescription,
            self.mass_transfer_angular_momentum_loss_prescription,
            self.mass_transfer_accretion_efficiency_prescription,
            self.mass_transfer_rejuvenation_prescription,
            self.initial_mass_function,
            self.semi_major_axis_distribution,
            self.orbital_period_distribution,
            self.mass_ratio_distribution,
            self.eccentricity_distribution,
            self.metallicity_distribution,
            self.rotational_velocity_distribution,
            self.remnant_mass_prescription,
            self.fryer_supernova_engine,
            self.black_hole_kicks,
            self.kick_magnitude_distribution,
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
            self.logfile_type,
            self.logfile_definitions,
            self.grid_filename,
            self.logfile_common_envelopes,
            self.logfile_detailed_output,
            self.logfile_double_compact_objects,
            self.logfile_pulsar_evolution,
            self.logfile_rlof_parameters,
            self.logfile_supernovae,
            self.logfile_switch_log,
            self.logfile_system_parameters,
            self.neutrino_mass_loss_BH_formation
        ]

        return stringChoices

    def stringCommands(self):
        stringCommands = [
            '--mode',
            '--case-BB-stability-prescription',
            '--chemically-homogeneous-evolution',
            '--luminous-blue-variable-prescription',
            '--mass-loss-prescription',
            '--mass-transfer-angular-momentum-loss-prescription',
            '--mass-transfer-accretion-efficiency-prescription',
            '--mass-transfer-rejuvenation-prescription',
            '--initial-mass-function',
            '--semi-major-axis-distribution',
            '--orbital-period-distribution',
            '--mass-ratio-distribution',
            '--eccentricity-distribution',
            '--metallicity-distribution',
            '--rotational-velocity-distribution',
            '--remnant-mass-prescription',
            '--fryer-supernova-engine',
            '--black-hole-kicks',
            '--kick-magnitude-distribution',
            '--kick-direction',
            '--output-path',
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
            '--logfile-type',
            '--logfile-definitions',
            '--grid',
            '--logfile-common-envelopes',
            '--logfile-detailed-output',
            '--logfile-double-compact-objects',
            '--logfile-pulsar-evolution',
            '--logfile-rlof-parameters',
            '--logfile-supernovae',
            '--logfile-switch-log',
            '--logfile-system-parameters',
            '--neutrino-mass-loss-BH-formation'
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


    def generateCommandLineOptionsDict(self):
        """
        This function generates a dictionary mapping COMPAS options to their specified 
        values (or empty strings for boolean options). These can be combined into a string
        and run directly as a terminal command, or passed to the stroopwafel interface
        where some of them may be overwritten. Options not to be included in the command 
        line should be set to pythons None (except booleans, which should be set to False)
    
        Parameters
        -----------
        self : pythonProgramOptions
            Contains program options
    
        Returns
        --------
        commands : str or list of strs
        """
        booleanChoices = self.booleanChoices()
        booleanCommands = self.booleanCommands()
        nBoolean = len(booleanChoices)
        assert len(booleanCommands) == nBoolean
    
        numericalChoices = self.numericalChoices()
        numericalCommands = self.numericalCommands()
        nNumerical = len(numericalChoices)
        assert len(numericalCommands) == nNumerical
    
        stringChoices = self.stringChoices()
        stringCommands = self.stringCommands()
        nString = len(stringChoices)
        assert len(stringCommands) == nString
    
        listChoices = self.listChoices()
        listCommands = self.listCommands()
        nList = len(listChoices)
        assert len(listCommands) == nList


        ### Collect all options into a dictionary mapping option name to option value

        command = {'compas_executable' : self.compas_executable}
    
        for i in range(nBoolean):
            if booleanChoices[i] == True:
                command.update({booleanCommands[i] : ''})
    
        for i in range(nNumerical):
            if not numericalChoices[i] == None:
                command.update({numericalCommands[i] : str(numericalChoices[i])})
    
        for i in range(nString):
            if not stringChoices[i] == None:
                command.update({stringCommands[i] : stringChoices[i]})
    
        for i in range(nList):
            if listChoices[i]:
                command.update({listCommands[i] : ' '.join(map(str,listChoices[i]))})
    
        return command


def combineCommandLineOptionsDictIntoShellCommand(commandOptions):
    """
    Write out the compas input parameters into a shell string.
    Ensure the Compas executable is first, and not repeated.
    Options are non-ordered.
    """

    shellCommand = commandOptions['compas_executable']
    del commandOptions['compas_executable'] 
    for key, val in commandOptions.items():
        shellCommand += ' ' + key + ' ' + val

    return shellCommand


if __name__ == "__main__":

    #-- Get the program options
    programOptions = pythonProgramOptions()
    commandOptions = programOptions.generateCommandLineOptionsDict()

    #-- Convert options into a shell string
    shellCommand = combineCommandLineOptionsDictIntoShellCommand(commandOptions)

    #-- Run exectute COMPAS shell string
    print(shellCommand)
    call(shellCommand,shell=True)

