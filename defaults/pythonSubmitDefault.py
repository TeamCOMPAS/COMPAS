import numpy as np
import subprocess
import sys
import os
from subprocess import call

class pythonProgramOptions:
    """
    A class to store and access COMPAS program options in python
    """

    # Do './COMPAS --help' to see all options
    #-- Define variables
    git_directory = os.environ.get('COMPAS_ROOT_DIR')
    compas_executable = os.path.join(git_directory, 'src/COMPAS')
    number_of_binaries = 10  #number of binaries per batch
    populationPrinting = False

    randomSeedFileName = 'randomSeed.txt'
    if os.path.isfile(randomSeedFileName):
        random_seed = int(np.loadtxt(randomSeedFileName))
    else:
        random_seed = 0 # If you want a randome seed, use: np.random.randint(2,2**63-1)

    output = os.getcwd()
    output_container = None                 # names the directory to be created and in which log files are created.  Default in COMPAS is "COMPAS_Output"

    #-- set inidividual system parameters
    single_star = False

    grid_filename = None

    use_mass_loss = True
    mass_transfer = True
    # post_newtonian_evolution = True
    detailed_output = False                         # WARNING: this creates a data heavy file
    only_double_compact_objects = False             # Delete when STROOPWAFEL fully implemented
    evolve_unbound_systems = False
    lambda_calculation_every_timestep = False
    zeta_calculation_every_timestep = False
    quiet = False

    metallicity = 0.0142                            # Solar metallicity Asplund+2010

    allow_rlof_at_birth = False;                                            # allow binaries that have one or both stars in RLOF at birth to evolve?
    allow_touching_at_birth = False;                                        # allow binaries that have stars touching at birth to evolve?

    chemically_homogeneous_evolution = 'NONE'                               # chemically homogeneous evolution.  Options are 'NONE', 'OPTIMISTIC' and 'PESSIMISTIC'

    common_envelope_alpha = 1.0
    common_envelope_lambda = 0.1                # Only if using 'LAMBDA_FIXED'
    common_envelope_lambda_prescription = 'LAMBDA_NANJING'  # Xu & Li 2010
    common_envelope_slope_Kruckow = -5.0/6.0
    common_envelope_zeta_prescription = 'SOBERMAN'
    common_envelope_revised_energy_formalism = False
    common_envelope_maximum_donor_mass_revised_energy_formalism = 2.0
    common_envelope_recombination_energy_density = 1.5E13
    common_envelope_alpha_thermal = 1.0                     # lambda = alpha_th*lambda_b + (1-alpha_th)*lambda_g
    common_envelope_lambda_multiplier = 1.0                 # Multiply common envelope lambda by some constant
    common_envelope_allow_main_sequence_survive = True      # Allow main sequence stars to survive CE. Was previously False by default
    common_envelope_mass_accretion_prescription = 'ZERO'
    common_envelope_mass_accretion_min = 0.04           # For 'MACLEOD+2014' [Msol]
    common_envelope_mass_accretion_max = 0.10           # For 'MACLEOD+2014' [Msol]

    tides_prescription = 'NONE'
    mass_loss_prescription = 'VINK'
    luminous_blue_variable_multiplier = 1.5
    wolf_rayet_multiplier = 1.0

    circularise_binary_during_mass_transfer = False
    angular_momentum_conservation_during_circularisation = False
    mass_transfer_prescription = 'DEMINK'           # Remove after cleaning MT function
    mass_transfer_angular_momentum_loss_prescription = 'ISOTROPIC'
    mass_transfer_accretion_efficiency_prescription = 'THERMAL'
    mass_transfer_fa = 0.5  # Only if using mass_transfer_accretion_efficiency_prescription = 'FIXED'
    mass_transfer_jloss = 1.0   # Only if using mass_transfer_angular_momentum_loss_prescription = 'FIXED'
    mass_transfer_rejuvenation_prescription = 'STARTRACK'
    mass_transfer_thermal_limit_accretor= 'CFACTOR'
    mass_transfer_thermal_limit_C= 10.0
    eddington_accretion_factor = 1    #multiplication Factor for eddington accretion onto NS&BH

    #-- Stability criteria for case BB/BC mass transfer (for BNS project)
    force_case_BB_BC_stability = True               # Check redundancy
    always_stable_case_BB_BC = True             # Check redundancy
    zeta_Main_Sequence = 2.0
    zeta_Hertzsprung_Gap = 6.5

    #-- Critical mass ratios for MT. 0.0 for always stable, < 0.0 to disable
    critical_mass_ratio_MS_low_mass_non_degenerate_accretor = -1.0 #1.44
    critical_mass_ratio_MS_low_mass_degenerate_accretor = -1.0 #1.0
    critical_mass_ratio_MS_high_mass_non_degenerate_accretor = -1.0 #0.625
    critical_mass_ratio_MS_high_mass_degenerate_accretor = -1.0 #0.0
    critical_mass_ratio_HG_non_degenerate_accretor = -1.0 #0.25
    critical_mass_ratio_HG_degenerate_accretor = -1.0 #0.21
    critical_mass_ratio_giant_non_degenerate_accretor = -1.0 #0.0
    critical_mass_ratio_giant_degenerate_accretor = -1.0 #0.87
    critical_mass_ratio_helium_MS_non_degenerate_accretor = -1.0 #0.625
    critical_mass_ratio_helium_MS_degenerate_accretor = -1.0 #0.0
    critical_mass_ratio_helium_HG_non_degenerate_accretor = -1.0 #0.25
    critical_mass_ratio_helium_HG_degenerate_accretor = -1.0 #0.21
    critical_mass_ratio_helium_giant_non_degenerate_accretor = -1.0 #1.28
    critical_mass_ratio_helium_giant_degenerate_accretor = -1.0 #0.87
    critical_mass_ratio_white_dwarf_non_degenerate_accretor = -1.0 #0.0
    critical_mass_ratio_white_dwarf_degenerate_accretor = -1.0 #1.6


    BeBinaries  = False
    maximum_evolution_time = 13700.0
    maximum_number_iterations = 99999


    #  STROOPWAFEL algorithm for COMPAS cf Broekgaarden+18
    #  For doumentation see the COMPAS/COMPAS/AdaptiveImportanceSampling folder on gitlab
    AIS_exploratory_phase =  False  # If TRUE COMPAS will run the STROOPWAFEL/Adaptive Importance Sampling algorithm
    #  parameters to select the DCOs of interest
    AIS_DCOtype = 'BBH'         # select "ALL", "BBH", "BHNS" or "BNS"
    AIS_Hubble = True           # focus on DCOs that merge within a Hubble time
    AIS_RLOF = True             # Focus on DCOs that have RLOF flag True
    AIS_Pessimistic = False         # Focus on only Pessimistic binaries

    kappa_gaussians = 2  # scaling factor for width of Gaussian distributions |  default = 2
    #  Run Refinement phase. This is automatically set to True after end exploratory phase:
    AIS_refinement_phase = False   # Only set True in case you want to reproduce binaries refinement phase

    initial_mass_function = 'KROUPA'
    initial_mass_min = 5.0          # Use 1.0 for LRNe, 5.0 for DCOs  [Msol]
    initial_mass_max = 150.0        # Stellar tracks extrapolated above 50 Msol (Hurley+2000) [Msol]

    initial_mass_power = 0.0

    semi_major_axis_distribution = 'FLATINLOG'
    semi_major_axis_min = 0.01  #  [AU]
    semi_major_axis_max = 1000.0  #  [AU]

    spin_distribution = 'ZERO'
    spin_assumption = 'BOTHALIGNED'
    spin_mag_min = 0.0
    spin_mag_max = 1.0

    mass_ratio_distribution = 'FLAT'
    mass_ratio_min = 0.0
    mass_ratio_max = 1.0

    minimum_secondary_mass = 0.1 # Brown dwarf limit  [Msol]

    eccentricity_distribution = 'ZERO'
    eccentricity_min = 0.0
    eccentricity_max = 1.0

    pulsar_birth_magnetic_field_distribution = 'ZERO'
    pulsar_birth_magnetic_field_min = 11.0
    pulsar_birth_magnetic_field_max = 13.0

    pulsar_birth_spin_period_distribution = "ZERO"
    pulsar_birth_spin_period_min = 10.0
    pulsar_birth_spin_period_max = 100.0

    pulsar_magnetic_field_decay_timescale = 1000.0
    pulsar_magnetic_field_decay_massscale = 0.025
    pulsar_minimum_magnetic_field = 8.0

    evolvePulsars = False

    rotational_velocity_distribution = 'ZERO'

    neutron_star_equation_of_state = 'SSE'

    orbital_period_min = 1.1
    orbital_period_max = 1000

    sample_kick_velocity_sigma = False
    sample_kick_velocity_sigma_min = 0.0
    sample_kick_velocity_sigma_max = 400.0

    sample_kick_direction_power = False
    sample_kick_direction_power_min = -10.0
    sample_kick_direction_power_max = 10.0

    sample_common_envelope_alpha = False
    sample_common_envelope_alpha_min = 0.0
    sample_common_envelope_alpha_max = 5.0

    sample_wolf_rayet_multiplier = False
    sample_wolf_rayet_multiplier_min = 0.0
    sample_wolf_rayet_multiplier_max = 10.0

    sample_luminous_blue_variable_multiplier = False
    sample_luminous_blue_variable_multiplier_min = 0.0
    sample_luminous_blue_variable_multiplier_max = 10.0

    remnant_mass_prescription = 'FRYER2012'
    fryer_supernova_engine = 'DELAYED'
    black_hole_kicks = 'FALLBACK'
    kick_velocity_distribution = 'MAXWELLIAN'

    kick_velocity_sigma_CCSN_NS = 265.0  #  [km/s]
    kick_velocity_sigma_CCSN_BH = 265.0  #  [km/s]
    kick_velocity_sigma_ECSN = 30.0  #  [km/s]
    kick_velocity_sigma_USSN = 30.0  #  [km/s]

    fix_dimensionless_kick_velocity = -1
    kick_direction = 'ISOTROPIC'
    kick_direction_power = 0.0
    kick_scaling_factor = 1.0
    kick_velocity_maximum = -1.0

    pair_instability_supernovae = True
    PISN_lower_limit = 60.0     # Minimum core mass for PISN [Msol]
    PISN_upper_limit = 135.0    # Maximum core mass for PISN [Msol]
    pulsation_pair_instability = True
    PPI_lower_limit = 35.0      # Minimum core mass for PPI [Msol]
    PPI_upper_limit = 60.0      # Maximum core mass for PPI [Msol]

    pulsational_pair_instability_prescription = 'MARCHANT'

    maximum_neutron_star_mass = 2.5  #  [Msol]

    log_level         = 0
    log_classes       = []

    debug_level       = 0
    debug_classes     = []

    logfile_definitions = None

    logfile_name_prefix = 'Compas_Log_'
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

    single_star_mass_steps = 10
    single_star_mass_min   = 1.0
    single_star_mass_max   = 75.0


    # read in nBatces for STROOPWAFEL if running Adaptive Sampling (AIS algorithm)  # you should not change this
    if AIS_exploratory_phase == True:
        compashpc_directory =  git_directory + "/CompasHPC/"
        sys.path.append(compashpc_directory)
        from compas_hpc_input import nBatches
        nbatches_used = nBatches
    else:
        nbatches_used = -1


    def booleanChoices(self):
        booleanChoices = [
            self.single_star,
            self.use_mass_loss,
            self.mass_transfer,
#            self.post_newtonian_evolution,
            self.detailed_output,
            self.only_double_compact_objects,
            self.evolve_unbound_systems,
            self.sample_kick_velocity_sigma,
            self.sample_kick_direction_power,
            self.sample_common_envelope_alpha,
            self.sample_wolf_rayet_multiplier,
            self.sample_luminous_blue_variable_multiplier,
            self.populationPrinting,
            self.lambda_calculation_every_timestep,
            self.zeta_calculation_every_timestep,
            self.circularise_binary_during_mass_transfer,
            self.force_case_BB_BC_stability,
            self.always_stable_case_BB_BC,
            self.angular_momentum_conservation_during_circularisation,
            self.BeBinaries,
            self.AIS_exploratory_phase,
            self.AIS_Hubble,
            self.AIS_RLOF,
            self.AIS_Pessimistic,
            self.AIS_refinement_phase,
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
            '--single-star',
            '--use-mass-loss',
            '--massTransfer',
#            '--PNEcc',
            '--detailedOutput',
            '--only-double-compact-objects',
            '--evolve-unbound-systems',
            '--sample-kick-velocity-sigma',
            '--sample-kick-direction-power',
            '--sample-common-envelope-alpha',
            '--sample-wolf-rayet-multiplier',
            '--sample-luminous-blue-variable-multiplier',
            '--populationDataPrinting',
            '--lambda-calculation-every-timeStep',
            '--zeta-Calculation-Every-Time-Step',
            '--circulariseBinaryDuringMassTransfer',
            '--forceCaseBBBCStabilityFlag',
            '--alwaysStableCaseBBBCFlag',
            '--angularMomentumConservationDuringCircularisation',
            '--BeBinaries',
            '--AIS-exploratory-phase',
            '--AIS-Hubble',
            '--AIS-RLOF',
            '--AIS-Pessimistic',
            '--AIS-refinement-phase',
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
            self.maximum_number_iterations,
            self.kappa_gaussians,
            self.initial_mass_min,
            self.initial_mass_max,
            self.initial_mass_power,
            self.semi_major_axis_min,
            self.semi_major_axis_max,
            self.spin_mag_min,
            self.spin_mag_max,
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
            self.sample_kick_velocity_sigma_min,
            self.sample_kick_velocity_sigma_max,
            self.sample_kick_direction_power_min,
            self.sample_kick_direction_power_max,
            self.sample_common_envelope_alpha_min,
            self.sample_common_envelope_alpha_max,
            self.sample_wolf_rayet_multiplier_min,
            self.sample_wolf_rayet_multiplier_max,
            self.sample_luminous_blue_variable_multiplier_min,
            self.sample_luminous_blue_variable_multiplier_max,
            self.random_seed,
            self.critical_mass_ratio_MS_low_mass_non_degenerate_accretor,
            self.critical_mass_ratio_MS_low_mass_degenerate_accretor,
            self.critical_mass_ratio_MS_high_mass_non_degenerate_accretor,
            self.critical_mass_ratio_MS_high_mass_degenerate_accretor,
            self.critical_mass_ratio_HG_non_degenerate_accretor,
            self.critical_mass_ratio_HG_degenerate_accretor,
            self.critical_mass_ratio_giant_non_degenerate_accretor,
            self.critical_mass_ratio_giant_degenerate_accretor,
            self.critical_mass_ratio_helium_MS_non_degenerate_accretor,
            self.critical_mass_ratio_helium_MS_degenerate_accretor,
            self.critical_mass_ratio_helium_HG_non_degenerate_accretor,
            self.critical_mass_ratio_helium_HG_degenerate_accretor,
            self.critical_mass_ratio_helium_giant_non_degenerate_accretor,
            self.critical_mass_ratio_helium_giant_degenerate_accretor,
            self.critical_mass_ratio_white_dwarf_non_degenerate_accretor,
            self.critical_mass_ratio_white_dwarf_degenerate_accretor,
            self.mass_transfer_thermal_limit_C,
            self.eddington_accretion_factor,
            self.PISN_lower_limit,
            self.PISN_upper_limit,
            self.PPI_lower_limit,
            self.PPI_upper_limit,
            self.maximum_neutron_star_mass,
            self.nbatches_used,
            self.kick_velocity_sigma_ECSN,
            self.kick_velocity_sigma_USSN,
            self.kick_scaling_factor,
            self.common_envelope_maximum_donor_mass_revised_energy_formalism,
            self.common_envelope_recombination_energy_density,
            self.common_envelope_mass_accretion_max,
            self.common_envelope_mass_accretion_min,
            self.zeta_Main_Sequence,
            self.zeta_Hertzsprung_Gap,
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
            '--maximum-number-iterations',
            '--kappa-gaussians',
            '--initial-mass-min',
            '--initial-mass-max',
            '--initial-mass-power',
            '--semi-major-axis-min',
            '--semi-major-axis-max',
            '--spin-mag-min',
            '--spin-mag-max',
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
            '--sample-kick-velocity-sigma-min',
            '--sample-kick-velocity-sigma-max',
            '--sample-kick-direction-power-min',
            '--sample-kick-direction-power-max',
            '--sample-common-envelope-alpha-min',
            '--sample-common-envelope-alpha-max',
            '--sample-wolf-rayet-multiplier-min',
            '--sample-wolf-rayet-multiplier-max',
            '--sample-luminous-blue-variable-multiplier-min',
            '--sample-luminous-blue-variable-multiplier-max',
            '--random-seed',
            '--critical-mass-ratio-MS-low-mass-non-degenerate-accretor',
            '--critical-mass-ratio-MS-low-mass-degenerate-accretor',
            '--critical-mass-ratio-MS-high-mass-non-degenerate-accretor',
            '--critical-mass-ratio-MS-high-mass-degenerate-accretor',
            '--critical-mass-ratio-HG-non-degenerate-accretor',
            '--critical-mass-ratio-HG-degenerate-accretor',
            '--critical-mass-ratio-giant-non-degenerate-accretor',
            '--critical-mass-ratio-giant-degenerate-accretor',
            '--critical-mass-ratio-helium-MS-non-degenerate-accretor',
            '--critical-mass-ratio-helium-MS-degenerate-accretor',
            '--critical-mass-ratio-helium-HG-non-degenerate-accretor',
            '--critical-mass-ratio-helium-HG-degenerate-accretor',
            '--critical-mass-ratio-helium-giant-non-degenerate-accretor',
            '--critical-mass-ratio-helium-giant-degenerate-accretor',
            '--critical-mass-ratio-white-dwarf-non-degenerate-accretor',
            '--critical-mass-ratio-white-dwarf-degenerate-accretor',
            '--mass-transfer-thermal-limit-C',
            '--eddington-accretion-factor',
            '--PISN-lower-limit',
            '--PISN-upper-limit','--PPI-lower-limit',
            '--PPI-upper-limit',
            '--maximum-neutron-star-mass',
            '--nbatches-used',
            '--kick-velocity-sigma-ECSN',
            '--kick-velocity-sigma-USSN',
            '--kick-scaling-factor',
            '--maximum-mass-donor-Nandez-Ivanova',
            '--common-envelope-recombination-energy-density',
            '--common-envelope-mass-accretion-max',
            '--common-envelope-mass-accretion-min',
            '--zeta-main-sequence',
            '--zeta-hertzsprung-gap',
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
            self.chemically_homogeneous_evolution,
            self.tides_prescription,
            self.mass_loss_prescription,
            self.mass_transfer_prescription,
            self.mass_transfer_angular_momentum_loss_prescription,
            self.mass_transfer_accretion_efficiency_prescription,
            self.mass_transfer_rejuvenation_prescription,
            self.AIS_DCOtype,
            self.initial_mass_function,
            self.semi_major_axis_distribution,
            self.spin_distribution,
            self.spin_assumption,
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
            self.common_envelope_zeta_prescription,
            self.mass_transfer_thermal_limit_accretor,
            self.pulsational_pair_instability_prescription,
            self.neutron_star_equation_of_state,
            self.pulsar_birth_magnetic_field_distribution,
            self.pulsar_birth_spin_period_distribution,
            self.common_envelope_mass_accretion_prescription,
            self.logfile_name_prefix,
            self.logfile_delimiter,
            self.logfile_definitions,
            self.grid_filename,
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
            '--chemically-homogeneous-evolution',
            '--tides-prescription',
            '--mass-loss-prescription',
            '--mass-transfer-prescription',
            '--mass-transfer-angular-momentum-loss-prescription',
            '--mass-transfer-accretion-efficiency-prescription',
            '--mass-transfer-rejuvenation-prescription',
            '--AIS-DCOtype',
            '--initial-mass-function',
            '--semi-major-axis-distribution',
            '--spin-distribution',
            '--spin-assumption',
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
            '--common-envelope-zeta-prescription',
            '--mass-transfer-thermal-limit-accretor',
            '--pulsational-pair-instability-prescription',
            '--neutron-star-equation-of-state',
            '--pulsar-birth-magnetic-field-distribution',
            '--pulsar-birth-spin-period-distribution',
            '--common-envelope-mass-accretion-prescription',
            '--logfile-name-prefix',
            '--logfile-delimiter',
            '--logfile-definitions',
            '--grid',
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

    command = [generateCommandLineOptions(programOptions.compas_executable,booleanChoices,booleanCommands,numericalChoices,numericalCommands,stringChoices,stringCommands,listChoices,listCommands)]

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

