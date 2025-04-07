COMPAS errors
=============

Following is a list of COMPAS error numbers, corresponding symbolic name, and meaning:

0. NONE |br|
   No error
#. AMBIGUOUS_REMNANT_MASS_PRESCRIPTION |br|
   Insufficient information to prescribe remnant mass
#. ARGUMENT_RANGE_PARMS_EXPECTED_FP |br|
   Expected a floating point number for range start and increment for option
#. ARGUMENT_RANGE_COUNT_EXPECTED_ULINT |br|
   Expected an unsigned long integer for range count for option
#. ARGUMENT_RANGE_PARMS_EXPECTED_INT |br|
   Expected an integer for range parameters for option
#. ARGUMENT_RANGE_PARMS_EXPECTED_LFP |br|
   Expected a long double number for range start and increment for option
#. ARGUMENT_RANGE_PARMS_EXPECTED_LINT |br|
   Expected an long integer for range parameters for option
#. ARGUMENT_RANGE_PARMS_EXPECTED_ULINT |br|
   Expected an unsigned long integer for range parameters for option
#. ARGUMENT_RANGE_NOT_SUPPORTED |br|
   Argument range not supported for option
#. ARGUMENT_RANGE_NUM_PARMS |br|
   Argument range requires exactly three parameters
#. ARGUMENT_SET_EXPECTED_BOOL |br|
   All parameters of argument set must be boolean for option
#. ARGUMENT_SET_EXPECTED_NUMERIC |br|
   All parameters of argument set must be numeric for option
#. ARGUMENT_SET_NOT_SUPPORTED |br|
   Argument set not supported for option
#. BAD_LOGFILE_RECORD_SPECIFICATIONS |br|
   Logfile record specifications error
#. BINARY_EVOLUTION_STOPPED |br|
   Evolution of current binary stopped
#. BINARY_SIMULATION_STOPPED |br|
   Binaries simulation stopped
#. BOOST_OPTION_CMDLINE |br|
   Failed to initialise Boost options descriptions for commandline options
#. BOOST_OPTION_GRIDLINE |br|
   Failed to initialise Boost options descriptions for grid line options
#. BOOST_OPTION_INTERNAL_ERROR |br|
   Internal error: Boost vm, option
#. DIRECTORY_NOT_EMPTY |br|
   Directory not empty
#. EMPTY_FILE |br|
   File is empty
#. EMPTY_FILENAME |br|
   Filename is an empty string
#. ERROR |br|
   Unpsecified error
#. ERROR_PROCESSING_CMDLINE_OPTIONS |br|
   An error occurred while processing commandline options
#. ERROR_PROCESSING_GRIDLINE_OPTIONS |br|
   An error occurred while processing grid file options
#. EXPECTED_3D_VECTOR |br|
   Expected a vector of size 3
#. EXPECTED_ASSIGNMENT_OPERATOR |br|
   Expected assignment operator: one of { '=', '-=', '+=' }
#. EXPECTED_BINARY_PROPERTY |br|
   Expected binary logfile property: one of { (STAR_1|STAR_2|SUPERNOVA|COMPANION|BINARY)_PROPERTY, PROGRAM_OPTION }
#. EXPECTED_COMMA_OR_CLOSE_BRACE |br|
   Expected a comma ',' or close brace '}'
#. EXPECTED_INTEGER |br|
   Expected an integer
#. EXPECTED_LOGFILE_RECORD_NAME |br|
   Expected logfile record specifier
#. EXPECTED_NON_NEGATIVE_INTEGER |br|
   Expected an integer >= 0
#. EXPECTED_OPEN_BRACE |br|
   Expected open brace '{'
#. EXPECTED_POSITIVE_INTEGER |br|
   Expected an integer > 0
#. EXPECTED_PROPERTY_SPECIFIER |br|
   Expected a property specifier or close brace '}'
#. EXPECTED_SN_EVENT |br|
   Expected a supernova event
#. EXPECTED_STELLAR_PROPERTY |br|
   Expected stellar logfile property: one of { STAR_PROPERTY, PROGRAM_OPTION }
#. FILE_DOES_NOT_EXIST |br|
   File does not exist
#. FILE_NOT_CLOSED |br|
   Error closing file - file not closed
#. FILE_OPEN_ERROR |br|
   Error opening file
#. FILE_READ_ERROR |br|
   Error reading from file - data not read
#. FILE_WRITE_ERROR |br|
   Error writing to file - data not written
#. FLOATING_POINT_ERROR |br|
   Unspecified floating-point error
#. FLOATING_POINT_DIVBYZERO |br|
   Floating-point divide-by-zero
#. FLOATING_POINT_INVALID_ARGUMENT |br|
   Floating-point invalid argument
#. FLOATING_POINT_OVERFLOW |br|
   Floating-point overflow
#. FLOATING_POINT_UNDERFLOW |br|
   Floating-point underflow
#. GRID_OPTIONS_ERROR |br|
   Grid File Options error
#. HIGH_TEFF_WINDS |br|
   Winds being used at high temperature
#. INDEX_OUT_OF_RANGE |br|
   Index out of range
#. INVALID_DATA_TYPE |br|
   Invalid data type
#. INVALID_INITIAL_ATTRIBUTES |br|
   Initial attributes are not valid - evolution not possible
#. INVALID_MASS_TRANSFER_DONOR |br|
   Mass transfer from NS, BH, or Massless Remnant
#. INVALID_TYPE_EDDINGTON_RATE |br|
   Invalid stellar type for Eddington critical rate calculation
#. INVALID_TYPE_MT_MASS_RATIO |br|
   Invalid stellar type for mass ratio calculation
#. INVALID_TYPE_MT_THERMAL_TIMESCALE |br|
   Invalid stellar type for thermal timescale calculation
#. INVALID_TYPE_ZETA_CALCULATION |br|
   Invalid stellar type for Zeta calculation
#. INVALID_VALUE_FOR_BOOLEAN_OPTION |br|
   Invalid value specified for BOOLEAN option
#. INVALID_VALUE_IN_FILE |br|
   Invalid value in file
#. LAMBDA_NOT_POSITIVE |br|
   Lambda <= 0.0
#. LOW_GAMMA |br|
   Very massive prescription being extrapolated to low gamma (<0.5)
#. LOW_TEFF_WINDS |br|
   Winds being used at low temperature
#. MAXIMUM_MASS_LOST |br|
   Maximum mass lost during mass loss calculations
#. MISSING_VALUE |br|
   Missing value
#. MISSING_RIGHT_BRACKET |br|
   Missing ']'
#. NO_CONVERGENCE |br|
   No convergence
#. NO_LAMBDA_DEWI |br|
   Dewi lambda calculation not supported for stellar type
#. NO_LAMBDA_NANJING |br|
   Nanjing lambda calculation not supported for stellar type
#. NO_REAL_ROOTS |br|
   No real roots
#. NO_TIMESTEPS_READ |br|
   No user timesteps read
#. NOT_INITIALISED |br|
   Object not initialised
#. OPTION_NOT_SUPPORTED_IN_GRID_FILE |br|
   Option not supported in grid file
#. OUT_OF_BOUNDS |br|
   Value out of bounds
#. PROGRAM_OPTIONS_ERROR |br|
   Commandline Options error
#. RADIUS_NOT_POSITIVE |br|
   Radius <= 0.0
#. REVERT_FAILED |br|
   Revert to previous state failed
#. ROOT_FINDER_FAILED |br|
   Exception encountered in root finder
#. STELLAR_EVOLUTION_STOPPED |br|
   Evolution of current star stopped
#. STELLAR_SIMULATION_STOPPED |br|
   Stellar simulation stopped
#. STEPS_UP |br|
   Allowed evolution timesteps exceeded
#. SWITCH_NOT_TAKEN |br|
   Switch to new stellar type not performed
#. TIMESTEP_BELOW_MINIMUM |br|
   Timestep below minimum - timestep taken
#. TIMESTEPS_EXHAUSTED |br|
   Provided timesteps exhausted, but evolution not complete
#. TIMESTEPS_NOT_CONSUMED |br|
   Evolution complete, but provided timesteps not consumed
#. TIMES_UP |br|
   Allowed evolution time exceeded
#. TOO_MANY_MASS0_ITERATIONS |br|
   Reached maximum number of iterations when looking for effective initial mass Mass_0 to match desired stellar core of HG star following case A mass transfer
#. TOO_MANY_MASS0_TRIES |br|
   Reached maximum number of tries when looking for effective initial mass Mass_0 to match desired stellar core of HG star following case A mass transfer
#. TOO_MANY_OMEGA_ITERATIONS |br|
   Reached maximum number of iterations when looking for omega when circularising and synchronising for tides
#. TOO_MANY_OMEGA_TRIES |br|
   Reached maximum number of tries when looking for omega when circularising and synchronising for tides
#. TOO_MANY_PULSAR_SPIN_ITERATIONS |br|
   Reached maximum number of iterations calculating the pulsar birth spin period
#. TOO_MANY_PULSAR_MAG_ITERATIONS |br|
   Reached maximum number of iterations calculating the pulsar birth magnetic field
#. TOO_MANY_REMNANT_MASS_ITERATIONS |br|
   Reached maximum number of iterations when calcuating remnant mass (MULLERMANDEL)
#. TOO_MANY_RETRIES |br|
   Too many retries
#. TOO_MANY_RLOF_ITERATIONS |br|
   Reached maximum number of iterations when fitting star inside Roche Lobe in RLOF
#. TOO_MANY_RLOF_TRIES |br|
   Reached maximum number of tries when fitting star inside Roche Lobe in RLOF
#. TOO_MANY_TIMESTEPS_IN_TIMESTEPS_FILE |br|
   Number of timesteps in timestpes file exceeds maximum timesteps
#. UNABLE_TO_CREATE_DIRECTORY |br|
   Unable to create directory
#. UNABLE_TO_REMOVE_DIRECTORY |br|
   Unable to remove directory
#. UNEXPECTED_ACCRETION_REGIME |br|
   Unexpected accretion regime
#. UNEXPECTED_BINARY_PROPERTY |br|
   Unexpected binary property
#. UNEXPECTED_BINARY_PROPERTY_TYPE |br|
   Unexpected binary property type
#. UNEXPECTED_END_OF_FILE |br|
   Unexpected end of file
#. UNEXPECTED_LOG_FILE_TYPE |br|
   Unexpected log file type
#. UNEXPECTED_PROGRAM_OPTION |br|
   Unexpected program option
#. UNEXPECTED_PROPERTY |br|
   Unexpected property
#. UNEXPECTED_PROPERTY_TYPE |br|
   Unexpected property type
#. UNEXPECTED_SN_EVENT |br|
   Unexpected supernova event in this context
#. UNEXPECTED_STELLAR_PROPERTY |br|
   Unexpected stellar property
#. UNEXPECTED_STELLAR_PROPERTY_TYPE |br|
   Unexpected stellar property type
#. UNEXPECTED_STELLAR_TYP |br|
   Unexpected stellar type
#. UNEXPECTED_ZETA_PRESCRIPTION |br|
   Unexpected stellar zeta prescription
#. UNHANDLED_EXCEPTION |br|
   Unhandled exception
#. UNKNOWN_A_DISTRIBUTION |br|
   Unknown semi-major-axis distribution
#. UNKNOWN_ACCRETION_REGIME |br|
   Unknown accretion regime
#. UNKNOWN_BH_KICK_MODE |br|
   Unknown black hole kicks mode
#. UNKNOWN_BINARY_PROPERTY |br|
   Unknown binary property
#. UNKNOWN_CASE_BB_STABILITY_PRESCRIPTION |br|
   Unknown case BB/BC mass transfer stability prescription
#. UNKNOWN_CE_ACCRETION_PRESCRIPTION |br|
   Unknown common envelope accretion prescription
#. UNKNOWN_CE_FORMALISM |br|
   Unknown common envelope formalism
#. UNKNOWN_CE_LAMBDA_PRESCRIPTION |br|
   Unknown common envelope lambda prescription
#. UNKNOWN_DATA_TYPE |br|
   Unknown data type
#. UNKNOWN_ENVELOPE_STATE_PRESCRIPTION |br|
   Unknown envelope state prescription
#. UNKNOWN_ENVELOPE_TYPE |br|
   Unknown envelope type
#. UNKNOWN_INITIAL_MASS_FUNCTION |br|
   Unknown initial mass function (IMF)
#. UNKNOWN_KICK_DIRECTION_DISTRIBUTION |br|
   Unknown kick direction distribution
#. UNKNOWN_KICK_MAGNITUDE_DISTRIBUTION |br|
   Unknown kick magnitude distribution
#. UNKNOWN_LBV_MASS_LOSS_PRESCRIPTION |br|
   Unknown LBV mass loss prescription
#. UNKNOWN_LOGFILE |br|
   Unknown log file
#. UNKNOWN_MT_CASE |br|
   Unknown mass transfer case
#. UNKNOWN_MT_ACCRETION_EFFICIENCY_PRESCRIPTION |br|
   Unknown mass transfer accretion efficiency prescription
#. UNKNOWN_MT_ANGULAR_MOMENTUM_LOSS_PRESCRIPTION |br|
   Unknown mass transfer angular momentum loss prescription
#. UNKNOWN_MT_REJUVENATION_PRESCRIPTION |br|
   Unknown mass transfer rejuvenation prescription
#. UNKNOWN_MT_THERMALLY_LIMITED_VARIATION |br|
   Unknown mass transfer thermally limited variation
#. UNKNOWN_MASS_LOSS_PRESCRIPTION |br|
   Unknown mass loss prescription
#. UNKNOWN_NEUTRINO_MASS_LOSS_PRESCRIPTION |br|
   Unknown neutrino mass loss prescription
#. UNKNOWN_NS_EOS |br|
   Unknown NS equation-of-state
#. UNKNOWN_OB_MASS_LOSS_PRESCRIPTION |br|
   Unknown OB mass loss prescription
#. UNKNOWN_PPI_PRESCRIPTION |br|
   Unknown pulsational pair instability prescription
#. UNKNOWN_PROGRAM_OPTION |br|
   Unknown program option
#. UNKNOWN_PROPERTY_TYPE |br|
   Unknown property type
#. UNKNOWN_PULSAR_BIRTH_MAGNETIC_FIELD_DISTRIBUTION |br|
   Unknown pulsar birth magnetic field distribution
#. UNKNOWN_PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION |br|
   Unknown pulsar birth spin period distribution
#. UNKNOWN_Q_DISTRIBUTION |br|
   Unknown q-distribution
#. UNKNOWN_QCRIT_PRESCRIPTION |br|
   Unknown QCRIT prescription
#. UNKNOWN_REMNANT_MASS_PRESCRIPTION |br|
   Unknown remnant mass prescription
#. UNKNOWN_RSG_MASS_LOSS_PRESCRIPTION |br|
   Unknown RSG mass loss prescription
#. UNKNOWN_SEMI_MAJOR_AXIS_DISTRIBUTION |br|
   Unknown semi-major axis distribution
#. UNKNOWN_SN_ENGINE |br|
   Unknown supernova engine
#. UNKNOWN_SN_EVENT |br|
   Unknown supernova event
#. UNKNOWN_STELLAR_POPULATION |br|
   Unknown stellar population
#. UNKNOWN_STELLAR_PROPERTY |br|
   Unknown stellar property
#. UNKNOWN_STELLAR_TYPE |br|
   Unknown stellar type
#. UNKNOWN_TIDES_PRESCRIPTION |br|
   Unknown tides prescription
#. UNKNOWN_VMS_MASS_LOSS_PRESCRIPTION |br|
   Unknown VMS mass loss prescription
#. UNKNOWN_VROT_PRESCRIPTION |br|
   Unknown rotational velocity prescription
#. UNKNOWN_ZETA_PRESCRIPTION |br|
   Unknown stellar ZETA prescription
#. WARNING |br|
   Unspecified warning
#. WHITE_DWARF_TOO_MASSIVE |br|
   This white dwarf exceeds the Chandrasekhar mass limit
#. UNKNOWN_WR_MASS_LOSS_PRESCRIPTION |br|
   Unknown WR mass loss prescription
