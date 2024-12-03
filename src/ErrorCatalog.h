#ifndef __ErrorCatalog_h__
#define __ErrorCatalog_h__


// This is the COMPAS error catalogue.  The error catalogue defines symbolic names for all COMPAS errors,
// and corresponding error strings for those errors.
//
// To add a new error, add the symbolic name to the ERROR enum class, and the corresponding error string
// to the ERROR_CATALOG map.
// 
// The key to the ERROR_CATALOG map is the symbolic name in the ERROR enum class.  The map entry is a
// tuple containing the ERROR_SCOPE associated with the error (see below), and the error string.

#include "constants.h"
#include "EnumHash.h"

#define ERR_MSG(x) std::get<1>(ERROR_CATALOG.at(x))                 // for convenience


// The ERROR_SCOPE enum class, defined below, allows developers to specify when, if at all, a particular
// error/warning should be displayed by the SHOW_WARN* and SHOW_ERROR* macros defined in ErrorsMacros.h
// (and by the ERRORS-SHowIt() function defined in Errors.h, which is used by the aforementioned macros):
//
//      NEVER                : the error/warning should never be displayed
//      ALWAYS               : the error/warning should always be displayed
//      FIRST                : the error/warning should only be displayed the first time it is encountered
//      FIRST_IN_FUNCTION    : the error/warning should only be displayed the first time it is encountered for a function
//      FIRST_IN_STELLAR_TYPE: the error/warning should only be displayed the first time it is encountered for a stellar type (see enum class STELLAR_TYPE in typedefs.h)
//      FIRST_IN_OBJECT_TYPE : the error/warning should only be displayed the first time it is encountered for an object type (see enum class OBJECT_TYPE in typedefs.h)
//      FIRST_IN_OBJECT_ID   : the error/warning should only be displayed the first time it is encountered for an object id (each object is assigned a unique object id - e.g. a star or binary, each constituent star of a binary)
//
// The THROW_ERROR* macros defined in ErrorsMacros.h are not affected by ERROR_SCOPE

enum class ERROR_SCOPE: int { NEVER, ALWAYS, FIRST, FIRST_IN_FUNCTION, FIRST_IN_STELLAR_TYPE, FIRST_IN_OBJECT_TYPE, FIRST_IN_OBJECT_ID };


// enum class ERROR
// Symbolic names for errors and warnings (error strings below in ERRORLabel map)
// Listed alphabetically (except for 'NONE' - first so ERROR = 0 = NONE)
enum class ERROR: int {
    NONE,                                                           // no error
    AMBIGUOUS_REMNANT_MASS_PRESCRIPTION,                            // remnant mass unclear from available parameters
    ARGUMENT_RANGE_COUNT_EXPECTED_ULINT,                            // expected an unsigned long integer for range count for option
    ARGUMENT_RANGE_NOT_SUPPORTED,                                   // argument range not supported for option 
    ARGUMENT_RANGE_NUM_PARMS,                                       // argument range requires exactly three parameters
    ARGUMENT_RANGE_PARMS_EXPECTED_FP,                               // expected a floating point number for range start and increment for option
    ARGUMENT_RANGE_PARMS_EXPECTED_INT,                              // expected an integer for range parameters for option
    ARGUMENT_RANGE_PARMS_EXPECTED_LFP,                              // expected a long double number for range start and increment for option
    ARGUMENT_RANGE_PARMS_EXPECTED_LINT,                             // expected an long integer for range parameters for option
    ARGUMENT_RANGE_PARMS_EXPECTED_ULINT,                            // expected an unsigned long integer for range parameters for option
    ARGUMENT_SET_EXPECTED_BOOL,                                     // all parameters of argument set must be boolean for option
    ARGUMENT_SET_EXPECTED_NUMERIC,                                  // all parameters of argument set must be numeric for option
    ARGUMENT_SET_NOT_SUPPORTED,                                     // argument set not supported for option
    BAD_LOGFILE_RECORD_SPECIFICATIONS,                              // error in logfile record specifications
    BINARY_EVOLUTION_STOPPED,                                       // evolution of current binary stopped
    BINARY_SIMULATION_STOPPED,                                      // binary simulation stopped
    BOOST_OPTION_CMDLINE,                                           // failed to initialise Boost options descriptions for commandline options
    BOOST_OPTION_GRIDLINE,                                          // failed to initialise Boost options descriptions for grid line options
    BOOST_OPTION_INTERNAL_ERROR,                                    // Boost option internal error
    DIRECTORY_NOT_EMPTY,                                            // Directory not empty
    EMPTY_FILE,                                                     // file is empty (contains no content)
    EMPTY_FILENAME,                                                 // filename is an empty string
    ERROR,                                                          // unspecified error
    ERROR_PROCESSING_CMDLINE_OPTIONS,                               // an error occurred while processing commandline options
    ERROR_PROCESSING_GRIDLINE_OPTIONS,                              // an error occurred while processing grid file options
    EXPECTED_3D_VECTOR,                                             // expected a vector of size 3
    EXPECTED_ASSIGNMENT_OPERATOR,                                   // expected assignment operator
    EXPECTED_BINARY_PROPERTY,                                       // expected a binary property (STAR_1_, STAR_2_, SUPERNOVA_, COMPANION_, or BINARY_PROPERTY)
    EXPECTED_COMMA_OR_CLOSE_BRACE,                                  // expected a comma or close brace
    EXPECTED_INTEGER,                                               // expected an integer
    EXPECTED_LOGFILE_RECORD_NAME,                                   // expected logfile record name
    EXPECTED_NON_NEGATIVE_INTEGER,                                  // expected a non-negative integer
    EXPECTED_OPEN_BRACE,                                            // expected an open brace
    EXPECTED_POSITIVE_INTEGER,                                      // expected a positive integer
    EXPECTED_PROPERTY_SPECIFIER,                                    // expected a valid property specifier
    EXPECTED_SN_EVENT,                                              // expected a supernova event
    EXPECTED_STELLAR_PROPERTY,                                      // expected a stellar property (STAR_PROPERTY)
    FILE_DOES_NOT_EXIST,                                            // file does not exist
    FILE_NOT_CLOSED,                                                // error closing file - file not closed
    FILE_OPEN_ERROR,                                                // error opening file
    FILE_READ_ERROR,                                                // error reading from file - data not read
    FILE_WRITE_ERROR,                                               // error writing to file - data not written
    FLOATING_POINT_ERROR,                                           // unspecified floating-point error
    FLOATING_POINT_DIVBYZERO,                                       // floating-point divide-by-zero
    FLOATING_POINT_INVALID_ARGUMENT,                                // floating-point invalid
    FLOATING_POINT_OVERFLOW,                                        // floating-point overflow
    FLOATING_POINT_UNDERFLOW,                                       // floating-point underflow
    GRID_OPTIONS_ERROR,                                             // grid file options error
    HIGH_TEFF_WINDS,                                                // winds being used at high temperature
    INDEX_OUT_OF_RANGE,                                             // index supplied is out of range
    INVALID_DATA_TYPE,                                              // invalid data type
    INVALID_INITIAL_ATTRIBUTES,                                     // initial values of stellar or binary attributes are not valid - can't evolve star or binary
    INVALID_MASS_TRANSFER_DONOR,                                    // mass transfer from NS, BH or Massless Remnant
    INVALID_TYPE_EDDINGTON_RATE,                                    // invalid stellar type for Eddington critical rate calculation
    INVALID_TYPE_MT_MASS_RATIO,                                     // invalid stellar type for mass ratio calculation
    INVALID_TYPE_MT_THERMAL_TIMESCALE,                              // invalid stellar type for thermal timescale calculation
    INVALID_TYPE_ZETA_CALCULATION,                                  // invalid stellar type for Zeta calculation
    INVALID_VALUE_FOR_BOOLEAN_OPTION,                               // invalid value specified for boolean option
    INVALID_VALUE_IN_FILE,                                          // invalid value in file
    LAMBDA_NOT_POSITIVE,                                            // lambda is <= 0.0 - invalid
    LOW_GAMMA,                                                      // very massive mass-loss prescription being extrapolated to low gamma (<0.5)
    LOW_TEFF_WINDS,                                                 // winds being used at low temperature
    MAXIMUM_MASS_LOST,                                              // (WARNING) maximum mass lost during mass loss calculations
    MISSING_VALUE,                                                  // missing value (e.g. for program option)
    MISSING_RIGHT_BRACKET,                                          // missing right bracket (e.g. for program option range or set specification)
    NO_CONVERGENCE,                                                 // iterative process did not converge
    NO_LAMBDA_DEWI,                                                 // Dewi lambda calculation not supported for stellar type
    NO_LAMBDA_NANJING,                                              // Nanjing lambda calculation not supported for stellar type
    NO_REAL_ROOTS,                                                  // equation has no real roots
    NO_TIMESTEPS_READ,                                              // no user timesteps read
    NOT_INITIALISED,                                                // object not initialised
    OPTION_NOT_SUPPORTED_IN_GRID_FILE,                              // option not supported in grid file
    OUT_OF_BOUNDS,                                                  // value out of bounds
    PROGRAM_OPTIONS_ERROR,                                          // program options error
    RADIUS_NOT_POSITIVE,                                            // radius is <= 0.0 - invalid
    REVERT_FAILED,                                                  // revert to previous state failed
    ROOT_FINDER_FAILED,                                             // root finder threw an exception
    STELLAR_EVOLUTION_STOPPED,                                      // evolution of current star stopped
    STELLAR_SIMULATION_STOPPED,                                     // stellar simulation stopped
    STEPS_UP,                                                       // allowed evolution timesteps exceeded
    SWITCH_NOT_TAKEN,                                               // switch to new stellar type not performed
    TIMESTEP_BELOW_MINIMUM,                                         // timestep too small - below minimum
    TIMESTEPS_EXHAUSTED,                                            // timesteps provided exhausted, but evolution not complete
    TIMESTEPS_NOT_CONSUMED,                                         // evolution complete, but provided timesteps not consumed
    TIMES_UP,                                                       // allowed evolution time exceeded
    TOO_MANY_MASS0_ITERATIONS,                                      // too many iterations in MASS0 root finder
    TOO_MANY_MASS0_TRIES,                                           // too many tries in MASS0 root finder
    TOO_MANY_OMEGA_ITERATIONS,                                      // too many iterations in OMEGA root finder
    TOO_MANY_OMEGA_TRIES,                                           // too many tries in OMEGA root finder
    TOO_MANY_PULSAR_SPIN_ITERATIONS,                                // too many iterations calculating the pulsar birth spin period
    TOO_MANY_RETRIES,                                               // generic too many retries
    TOO_MANY_RLOF_ITERATIONS,                                       // too many iterations in RLOF root finder
    TOO_MANY_RLOF_TRIES,                                            // too many tries in RLOF root finder
    TOO_MANY_TIMESTEPS_IN_TIMESTEPS_FILE,                           // too many timesteps in timesteps file (exceeds maximum)
    UNABLE_TO_CREATE_DIRECTORY,                                     // unable to create directory
    UNABLE_TO_REMOVE_DIRECTORY,                                     // unable to remove directory
    UNEXPECTED_ACCRETION_REGIME,                                    // unexpected accretion regime
    UNEXPECTED_BINARY_PROPERTY,                                     // unexpected binary property
    UNEXPECTED_BINARY_PROPERTY_TYPE,                                // unexpected binary property type
    UNEXPECTED_END_OF_FILE,                                         // unexpected end of file
    UNEXPECTED_LOG_FILE_TYPE,                                       // unexpected log file type
    UNEXPECTED_PROGRAM_OPTION,                                      // unexpected program option
    UNEXPECTED_PROPERTY,                                            // unexpected property
    UNEXPECTED_PROPERTY_TYPE,                                       // unexpected property type
    UNEXPECTED_ROOT_FINDER_FAILURE,                                 // unexpected root finder failure
    UNEXPECTED_SN_EVENT,                                            // unexpected supernova event in this context
    UNEXPECTED_STELLAR_PROPERTY,                                    // unexpected stellar property
    UNEXPECTED_STELLAR_PROPERTY_TYPE,                               // unexpected stellar property type
    UNEXPECTED_STELLAR_TYPE,                                        // unexpected stellar type
    UNHANDLED_EXCEPTION,                                            // unhandled exception
    UNKNOWN_A_DISTRIBUTION,                                         // unknown a-distribution
    UNKNOWN_ACCRETION_REGIME,                                       // unknown accretion regime
    UNKNOWN_BH_KICK_MODE,                                           // unknown black hole kick mode
    UNKNOWN_BINARY_PROPERTY,                                        // unknown binary property
    UNKNOWN_CASE_BB_STABILITY_PRESCRIPTION,                         // unknown case BB/BC mass transfer stability prescription
    UNKNOWN_CE_ACCRETION_PRESCRIPTION,                              // unknown common envelope accretion prescription
    UNKNOWN_CE_FORMALISM,                                           // unknown common envelope formalism
    UNKNOWN_CE_LAMBDA_PRESCRIPTION,                                 // unknown common envelope Lambda prescription
    UNKNOWN_DATA_TYPE,                                              // unknown data type
    UNKNOWN_ENVELOPE_STATE_PRESCRIPTION,                            // unknown envelope state prescription
    UNKNOWN_ENVELOPE_TYPE,                                          // unknown envelope type
    UNKNOWN_INITIAL_MASS_FUNCTION,                                  // unknown initial mass function
    UNKNOWN_KICK_DIRECTION_DISTRIBUTION,                            // unknown kick direction distribution
    UNKNOWN_KICK_MAGNITUDE_DISTRIBUTION,                            // unknown kick magnitude distribution
    UNKNOWN_LOGFILE,                                                // unknown log file
    UNKNOWN_LBV_MASS_LOSS_PRESCRIPTION,                             // unknown LBV mass loss prescription
    UNKNOWN_MT_ACCRETION_EFFICIENCY_PRESCRIPTION,                   // unknown mass transfer accretion efficiency prescription
    UNKNOWN_MT_ANGULAR_MOMENTUM_LOSS_PRESCRIPTION,                  // unknown mass transfer angular momentum loss prescription
    UNKNOWN_MASS_LOSS_PRESCRIPTION,                                 // unknown mass loss prescription
    UNKNOWN_MT_CASE,                                                // unknown mass transfer case
    UNKNOWN_MT_REJUVENATION_PRESCRIPTION,                           // unknown mass transfer rejuvenation prescription
    UNKNOWN_MT_THERMALLY_LIMITED_VARIATION,                         // unknown mass transfer thermally limited variation
    UNKNOWN_NEUTRINO_MASS_LOSS_PRESCRIPTION,                        // unknown neutrino mass loss prescription
    UNKNOWN_NS_EOS,                                                 // unknown NS equation-of-state
    UNKNOWN_OB_MASS_LOSS_PRESCRIPTION,                              // unknown OB mass loss prescription
    UNKNOWN_PPI_PRESCRIPTION,                                       // unknown pulsational pair instability prescription
    UNKNOWN_PROGRAM_OPTION,                                         // unknown program option
    UNKNOWN_PROPERTY_TYPE,                                          // unknown property type
    UNKNOWN_PULSAR_BIRTH_MAGNETIC_FIELD_DISTRIBUTION,               // unknown pulsar birth magnetic field distribution
    UNKNOWN_PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION,                  // unknown pulsar birth spin period distribution
    UNKNOWN_Q_DISTRIBUTION,                                         // unknown q-distribution
    UNKNOWN_QCRIT_PRESCRIPTION,                                     // Unknown QCRIT prescription
    UNKNOWN_REMNANT_MASS_PRESCRIPTION,                              // unknown remnant mass prescriptrion
    UNKNOWN_RSG_MASS_LOSS_PRESCRIPTION,                             // unknown RSG mass loss prescription
    UNKNOWN_SEMI_MAJOR_AXIS_DISTRIBUTION,                           // unknown sem-major axis distribution
    UNKNOWN_SN_ENGINE,                                              // unknown supernova engine
    UNKNOWN_SN_EVENT,                                               // unknown supernova event encountered
    UNKNOWN_STELLAR_POPULATION,                                     // unknown stellar population
    UNKNOWN_STELLAR_PROPERTY,                                       // unknown stellar property
    UNKNOWN_STELLAR_TYPE,                                           // unknown stellar type
    UNKNOWN_TIDES_PRESCRIPTION,                                     // unknown tides prescription
    UNKNOWN_VMS_MASS_LOSS_PRESCRIPTION,                             // unknown VMS mass loss prescription
    UNKNOWN_VROT_PRESCRIPTION,                                      // unknown rorational velocity prescription
    UNKNOWN_WR_MASS_LOSS_PRESCRIPTION,                              // unknown WR mass loss prescription
    UNKNOWN_ZETA_PRESCRIPTION,                                      // unknown stellar ZETA prescription
    WARNING,                                                        // unspecified warning
    WHITE_DWARF_TOO_MASSIVE,                                        // a white dwarf exceeds the Chandrasekhar mass limit

    // not alphabetical - here to keep "real" error numbers consecutive
    SUGGEST_HELP                                                    // suggest using --help
};


// message catalog
// for now we'll just define it here - one day we might want to have it in a
// file and read it in at program start - that way we can be more flexible
// (i.e. change messages without recompiling, internationalise etc.)
//
// unordered_map - key is integer message number (from enum class ERROR above)
// listed alphabetically

const COMPASUnorderedMap<ERROR, std::tuple<ERROR_SCOPE, std::string>> ERROR_CATALOG = {
    { ERROR::AMBIGUOUS_REMNANT_MASS_PRESCRIPTION,                   { ERROR_SCOPE::ALWAYS,              "Insufficient information to prescribe remnant mass" }},
    { ERROR::ARGUMENT_RANGE_PARMS_EXPECTED_FP,                      { ERROR_SCOPE::ALWAYS,              "Expected a floating point number for range start and increment for option" }},
    { ERROR::ARGUMENT_RANGE_COUNT_EXPECTED_ULINT,                   { ERROR_SCOPE::ALWAYS,              "Expected an unsigned long integer for range count for option" }},
    { ERROR::ARGUMENT_RANGE_PARMS_EXPECTED_INT,                     { ERROR_SCOPE::ALWAYS,              "Expected an integer for range parameters for option" }},
    { ERROR::ARGUMENT_RANGE_PARMS_EXPECTED_LFP,                     { ERROR_SCOPE::ALWAYS,              "Expected a long double number for range start and increment for option" }},
    { ERROR::ARGUMENT_RANGE_PARMS_EXPECTED_LINT,                    { ERROR_SCOPE::ALWAYS,              "Expected an long integer for range parameters for option" }},
    { ERROR::ARGUMENT_RANGE_PARMS_EXPECTED_ULINT,                   { ERROR_SCOPE::ALWAYS,              "Expected an unsigned long integer for range parameters for option" }},
    { ERROR::ARGUMENT_RANGE_NOT_SUPPORTED,                          { ERROR_SCOPE::ALWAYS,              "Argument range not supported for option" }},
    { ERROR::ARGUMENT_RANGE_NUM_PARMS,                              { ERROR_SCOPE::ALWAYS,              "Argument range requires exactly three parameters" }},
    { ERROR::ARGUMENT_SET_EXPECTED_BOOL,                            { ERROR_SCOPE::ALWAYS,              "All parameters of argument set must be boolean for option" }},
    { ERROR::ARGUMENT_SET_EXPECTED_NUMERIC,                         { ERROR_SCOPE::ALWAYS,              "All parameters of argument set must be numeric for option" }},
    { ERROR::ARGUMENT_SET_NOT_SUPPORTED,                            { ERROR_SCOPE::ALWAYS,              "Argument set not supported for option" }},
    { ERROR::BAD_LOGFILE_RECORD_SPECIFICATIONS,                     { ERROR_SCOPE::ALWAYS,              "Logfile record specifications error" }},
    { ERROR::BINARY_EVOLUTION_STOPPED,                              { ERROR_SCOPE::ALWAYS,              "Evolution of current binary stopped" }},
    { ERROR::BINARY_SIMULATION_STOPPED,                             { ERROR_SCOPE::ALWAYS,              "Binaries simulation stopped" }},
    { ERROR::BOOST_OPTION_CMDLINE,                                  { ERROR_SCOPE::ALWAYS,              "Failed to initialise Boost options descriptions for commandline options" }},
    { ERROR::BOOST_OPTION_GRIDLINE,                                 { ERROR_SCOPE::ALWAYS,              "Failed to initialise Boost options descriptions for grid line options" }},
    { ERROR::BOOST_OPTION_INTERNAL_ERROR,                           { ERROR_SCOPE::ALWAYS,              "Internal error: Boost vm, option" }},
    { ERROR::DIRECTORY_NOT_EMPTY,                                   { ERROR_SCOPE::ALWAYS,              "Directory not empty" }},
    { ERROR::EMPTY_FILE,                                            { ERROR_SCOPE::ALWAYS,              "File is empty" }},
    { ERROR::EMPTY_FILENAME,                                        { ERROR_SCOPE::ALWAYS,              "Filename is an empty string" }},
    { ERROR::ERROR,                                                 { ERROR_SCOPE::ALWAYS,              "Error!" }},
    { ERROR::ERROR_PROCESSING_CMDLINE_OPTIONS,                      { ERROR_SCOPE::ALWAYS,              "An error occurred while processing commandline options" }},
    { ERROR::ERROR_PROCESSING_GRIDLINE_OPTIONS,                     { ERROR_SCOPE::ALWAYS,              "An error occurred while processing grid file options" }},
    { ERROR::EXPECTED_3D_VECTOR,                                    { ERROR_SCOPE::ALWAYS,              "Expected a vector of size 3" }},
    { ERROR::EXPECTED_ASSIGNMENT_OPERATOR,                          { ERROR_SCOPE::ALWAYS,              "Expected assignment operator: one of { '=', '-=', '+=' }" }},
    { ERROR::EXPECTED_BINARY_PROPERTY,                              { ERROR_SCOPE::ALWAYS,              "Expected binary logfile property: one of { (STAR_1|STAR_2|SUPERNOVA|COMPANION|BINARY)_PROPERTY, PROGRAM_OPTION }" }},
    { ERROR::EXPECTED_COMMA_OR_CLOSE_BRACE,                         { ERROR_SCOPE::ALWAYS,              "Expected a comma ',' or close brace '}'" }},
    { ERROR::EXPECTED_INTEGER,                                      { ERROR_SCOPE::ALWAYS,              "Expected an integer" }},
    { ERROR::EXPECTED_LOGFILE_RECORD_NAME,                          { ERROR_SCOPE::ALWAYS,              "Expected logfile record specifier" }},
    { ERROR::EXPECTED_NON_NEGATIVE_INTEGER,                         { ERROR_SCOPE::ALWAYS,              "Expected an integer >= 0" }},
    { ERROR::EXPECTED_OPEN_BRACE,                                   { ERROR_SCOPE::ALWAYS,              "Expected open brace '{'" }},
    { ERROR::EXPECTED_POSITIVE_INTEGER,                             { ERROR_SCOPE::ALWAYS,              "Expected an integer > 0" }},
    { ERROR::EXPECTED_PROPERTY_SPECIFIER,                           { ERROR_SCOPE::ALWAYS,              "Expected a property specifier or close brace '}'" }},
    { ERROR::EXPECTED_SN_EVENT,                                     { ERROR_SCOPE::ALWAYS,              "Expected a supernova event" }},
    { ERROR::EXPECTED_STELLAR_PROPERTY,                             { ERROR_SCOPE::ALWAYS,              "Expected stellar logfile property: one of { STAR_PROPERTY, PROGRAM_OPTION }" }},
    { ERROR::FILE_DOES_NOT_EXIST,                                   { ERROR_SCOPE::ALWAYS,              "File does not exist" }},
    { ERROR::FILE_NOT_CLOSED,                                       { ERROR_SCOPE::ALWAYS,              "Error closing file - file not closed" }},
    { ERROR::FILE_OPEN_ERROR,                                       { ERROR_SCOPE::ALWAYS,              "Error opening file" }},
    { ERROR::FILE_READ_ERROR,                                       { ERROR_SCOPE::ALWAYS,              "Error reading from file - data not read" }},
    { ERROR::FILE_WRITE_ERROR,                                      { ERROR_SCOPE::ALWAYS,              "Error writing to file - data not written" }},
    { ERROR::FLOATING_POINT_ERROR,                                  { ERROR_SCOPE::ALWAYS,              "Unspecified floating-point error" }},
    { ERROR::FLOATING_POINT_DIVBYZERO,                              { ERROR_SCOPE::ALWAYS,              "Floating-point divide-by-zero" }},
    { ERROR::FLOATING_POINT_INVALID_ARGUMENT,                       { ERROR_SCOPE::ALWAYS,              "Floating-point invalid argument" }},
    { ERROR::FLOATING_POINT_OVERFLOW,                               { ERROR_SCOPE::ALWAYS,              "Floating-point overflow" }},
    { ERROR::FLOATING_POINT_UNDERFLOW,                              { ERROR_SCOPE::ALWAYS,              "Floating-point underflow" }},
    { ERROR::GRID_OPTIONS_ERROR,                                    { ERROR_SCOPE::ALWAYS,              "Grid File Options error" }},
    { ERROR::HIGH_TEFF_WINDS,                                       { ERROR_SCOPE::FIRST_IN_OBJECT_ID,  "Winds being used at high temperature" }},
    { ERROR::INDEX_OUT_OF_RANGE,                                    { ERROR_SCOPE::ALWAYS,              "Index out of range" }},
    { ERROR::INVALID_DATA_TYPE,                                     { ERROR_SCOPE::ALWAYS,              "Invalid data type" }},
    { ERROR::INVALID_INITIAL_ATTRIBUTES,                            { ERROR_SCOPE::ALWAYS,              "Initial attributes are not valid - evolution not possible" }},
    { ERROR::INVALID_MASS_TRANSFER_DONOR,                           { ERROR_SCOPE::ALWAYS,              "Mass transfer from NS, BH, or Massless Remnant" }},
    { ERROR::INVALID_TYPE_EDDINGTON_RATE,                           { ERROR_SCOPE::ALWAYS,              "Invalid stellar type for Eddington critical rate calculation" }},
    { ERROR::INVALID_TYPE_MT_MASS_RATIO,                            { ERROR_SCOPE::ALWAYS,              "Invalid stellar type for mass ratio calculation" }},
    { ERROR::INVALID_TYPE_MT_THERMAL_TIMESCALE,                     { ERROR_SCOPE::ALWAYS,              "Invalid stellar type for thermal timescale calculation" }},
    { ERROR::INVALID_TYPE_ZETA_CALCULATION,                         { ERROR_SCOPE::ALWAYS,              "Invalid stellar type for Zeta calculation" }},
    { ERROR::INVALID_VALUE_FOR_BOOLEAN_OPTION,                      { ERROR_SCOPE::ALWAYS,              "Invalid value specified for BOOLEAN option" }},
    { ERROR::INVALID_VALUE_IN_FILE,                                 { ERROR_SCOPE::ALWAYS,              "Invalid value in file" }},
    { ERROR::LAMBDA_NOT_POSITIVE,                                   { ERROR_SCOPE::FIRST_IN_OBJECT_ID,  "Lambda <= 0.0" }},
    { ERROR::LOW_GAMMA,                                             { ERROR_SCOPE::ALWAYS,              "Very massive prescription being extrapolated to low gamma (<0.5)" }},
    { ERROR::LOW_TEFF_WINDS,                                        { ERROR_SCOPE::ALWAYS,              "Winds being used at low temperature" }},
    { ERROR::MAXIMUM_MASS_LOST,                                     { ERROR_SCOPE::FIRST_IN_OBJECT_ID,  "Maximum mass lost during mass loss calculations" }},
    { ERROR::MISSING_VALUE,                                         { ERROR_SCOPE::ALWAYS,              "Missing value" }},
    { ERROR::MISSING_RIGHT_BRACKET,                                 { ERROR_SCOPE::ALWAYS,              "Missing ']'" }},
    { ERROR::NO_CONVERGENCE,                                        { ERROR_SCOPE::ALWAYS,              "No convergence" }},
    { ERROR::NO_LAMBDA_DEWI,                                        { ERROR_SCOPE::ALWAYS,              "Dewi lambda calculation not supported for stellar type" }},
    { ERROR::NO_LAMBDA_NANJING,                                     { ERROR_SCOPE::ALWAYS,              "Nanjing lambda calculation not supported for stellar type" }},
    { ERROR::NO_REAL_ROOTS,                                         { ERROR_SCOPE::ALWAYS,              "No real roots" }},
    { ERROR::NO_TIMESTEPS_READ,                                     { ERROR_SCOPE::ALWAYS,              "No user timesteps read" }},
    { ERROR::NONE,                                                  { ERROR_SCOPE::NEVER,               "No error" }},
    { ERROR::NOT_INITIALISED,                                       { ERROR_SCOPE::ALWAYS,              "Object not initialised" }},
    { ERROR::OPTION_NOT_SUPPORTED_IN_GRID_FILE,                     { ERROR_SCOPE::ALWAYS,              "Option not supported in grid file" }},
    { ERROR::OUT_OF_BOUNDS,                                         { ERROR_SCOPE::ALWAYS,              "Value out of bounds" }},
    { ERROR::PROGRAM_OPTIONS_ERROR,                                 { ERROR_SCOPE::ALWAYS,              "Commandline Options error" }},
    { ERROR::RADIUS_NOT_POSITIVE,                                   { ERROR_SCOPE::ALWAYS,              "Radius <= 0.0" }},
    { ERROR::REVERT_FAILED,                                         { ERROR_SCOPE::ALWAYS,              "Revert to previous state failed" }},
    { ERROR::ROOT_FINDER_FAILED,                                    { ERROR_SCOPE::ALWAYS,              "Exception encountered in root finder" }},
    { ERROR::STELLAR_EVOLUTION_STOPPED,                             { ERROR_SCOPE::ALWAYS,              "Evolution of current star stopped" }},
    { ERROR::STELLAR_SIMULATION_STOPPED,                            { ERROR_SCOPE::ALWAYS,              "Stellar simulation stopped" }},
    { ERROR::STEPS_UP,                                              { ERROR_SCOPE::ALWAYS,              "Allowed evolution timesteps exceeded" }},
    { ERROR::SUGGEST_HELP,                                          { ERROR_SCOPE::ALWAYS,              "Use option '-h' (or '--help') to see (descriptions of) available options" }},
    { ERROR::SWITCH_NOT_TAKEN,                                      { ERROR_SCOPE::ALWAYS,              "Switch to new stellar type not performed" }},
    { ERROR::TIMESTEP_BELOW_MINIMUM,                                { ERROR_SCOPE::ALWAYS,              "Timestep below minimum - timestep taken" }},
    { ERROR::TIMESTEPS_EXHAUSTED,                                   { ERROR_SCOPE::ALWAYS,              "Provided timesteps exhausted, but evolution not complete" }},
    { ERROR::TIMESTEPS_NOT_CONSUMED,                                { ERROR_SCOPE::ALWAYS,              "Evolution complete, but provided timesteps not consumed" }},
    { ERROR::TIMES_UP,                                              { ERROR_SCOPE::ALWAYS,              "Allowed evolution time exceeded" }},
    { ERROR::TOO_MANY_MASS0_ITERATIONS,                             { ERROR_SCOPE::ALWAYS,              "Reached maximum number of iterations when looking for effective initial mass Mass_0 to match desired stellar core of HG star following case A mass transfer" }},
    { ERROR::TOO_MANY_MASS0_TRIES,                                  { ERROR_SCOPE::ALWAYS,              "Reached maximum number of tries when looking for effective initial mass Mass_0 to match desired stellar core of HG star following case A mass transfer" }},
    { ERROR::TOO_MANY_OMEGA_ITERATIONS,                             { ERROR_SCOPE::ALWAYS,              "Reached maximum number of iterations when looking for omega when circularising and synchronising for tides" }},
    { ERROR::TOO_MANY_OMEGA_TRIES,                                  { ERROR_SCOPE::ALWAYS,              "Reached maximum number of tries when looking for omega when circularising and synchronising for tides" }},
    { ERROR::TOO_MANY_PULSAR_SPIN_ITERATIONS,                       { ERROR_SCOPE::ALWAYS,              "Reached maximum number of iterations calculating the pulsar birth spin period" }},
    { ERROR::TOO_MANY_RETRIES,                                      { ERROR_SCOPE::ALWAYS,              "Too many retries" }},
    { ERROR::TOO_MANY_RLOF_ITERATIONS,                              { ERROR_SCOPE::ALWAYS,              "Reached maximum number of iterations when fitting star inside Roche Lobe in RLOF" }},
    { ERROR::TOO_MANY_RLOF_TRIES,                                   { ERROR_SCOPE::ALWAYS,              "Reached maximum number of tries when fitting star inside Roche Lobe in RLOF" }},
    { ERROR::TOO_MANY_TIMESTEPS_IN_TIMESTEPS_FILE,                  { ERROR_SCOPE::ALWAYS,              "Number of timesteps in timestpes file exceeds maximum timesteps" }},
    { ERROR::UNABLE_TO_CREATE_DIRECTORY,                            { ERROR_SCOPE::ALWAYS,              "Unable to create directory" }},
    { ERROR::UNABLE_TO_REMOVE_DIRECTORY,                            { ERROR_SCOPE::ALWAYS,              "Unable to remove directory" }},
    { ERROR::UNEXPECTED_ACCRETION_REGIME,                           { ERROR_SCOPE::ALWAYS,              "Unexpected accretion regime" }},
    { ERROR::UNEXPECTED_BINARY_PROPERTY,                            { ERROR_SCOPE::ALWAYS,              "Unexpected binary property" }},
    { ERROR::UNEXPECTED_BINARY_PROPERTY_TYPE,                       { ERROR_SCOPE::ALWAYS,              "Unexpected binary property type" }},
    { ERROR::UNEXPECTED_END_OF_FILE,                                { ERROR_SCOPE::ALWAYS,              "Unexpected end of file" }},
    { ERROR::UNEXPECTED_LOG_FILE_TYPE,                              { ERROR_SCOPE::ALWAYS,              "Unexpected log file type" }},
    { ERROR::UNEXPECTED_PROGRAM_OPTION,                             { ERROR_SCOPE::ALWAYS,              "Unexpected program option" }},
    { ERROR::UNEXPECTED_PROPERTY,                                   { ERROR_SCOPE::ALWAYS,              "Unexpected property" }},
    { ERROR::UNEXPECTED_PROPERTY_TYPE,                              { ERROR_SCOPE::ALWAYS,              "Unexpected property type" }},
    { ERROR::UNEXPECTED_ROOT_FINDER_FAILURE,                        { ERROR_SCOPE::ALWAYS,              "Unexpected root finder failure" }},
    { ERROR::UNEXPECTED_SN_EVENT,                                   { ERROR_SCOPE::ALWAYS,              "Unexpected supernova event in this context" }},
    { ERROR::UNEXPECTED_STELLAR_PROPERTY,                           { ERROR_SCOPE::ALWAYS,              "Unexpected stellar property" }},
    { ERROR::UNEXPECTED_STELLAR_PROPERTY_TYPE,                      { ERROR_SCOPE::ALWAYS,              "Unexpected stellar property type" }},
    { ERROR::UNEXPECTED_STELLAR_TYPE,                               { ERROR_SCOPE::ALWAYS,              "Unexpected stellar type" }},
    { ERROR::UNHANDLED_EXCEPTION,                                   { ERROR_SCOPE::ALWAYS,              "Unhandled exception" }},
    { ERROR::UNKNOWN_A_DISTRIBUTION,                                { ERROR_SCOPE::ALWAYS,              "Unknown semi-major-axis distribution" }},
    { ERROR::UNKNOWN_ACCRETION_REGIME,                              { ERROR_SCOPE::ALWAYS,              "Unknown accretion regime" }},
    { ERROR::UNKNOWN_BH_KICK_MODE,                                  { ERROR_SCOPE::ALWAYS,              "Unknown black hole kicks mode" }},
    { ERROR::UNKNOWN_BINARY_PROPERTY,                               { ERROR_SCOPE::ALWAYS,              "Unknown binary property" }},
    { ERROR::UNKNOWN_CASE_BB_STABILITY_PRESCRIPTION,                { ERROR_SCOPE::ALWAYS,              "Unknown case BB/BC mass transfer stability prescription" }},
    { ERROR::UNKNOWN_CE_ACCRETION_PRESCRIPTION,                     { ERROR_SCOPE::ALWAYS,              "Unknown common envelope accretion prescription" }},
    { ERROR::UNKNOWN_CE_FORMALISM,                                  { ERROR_SCOPE::ALWAYS,              "Unknown common envelope formalism" }},
    { ERROR::UNKNOWN_CE_LAMBDA_PRESCRIPTION,                        { ERROR_SCOPE::ALWAYS,              "Unknown common envelope lambda prescription" }},
    { ERROR::UNKNOWN_DATA_TYPE,                                     { ERROR_SCOPE::ALWAYS,              "Unknown data type" }},
    { ERROR::UNKNOWN_ENVELOPE_STATE_PRESCRIPTION,                   { ERROR_SCOPE::ALWAYS,              "Unknown envelope state prescription" }},
    { ERROR::UNKNOWN_ENVELOPE_TYPE,                                 { ERROR_SCOPE::ALWAYS,              "Unknown envelope type" }},
    { ERROR::UNKNOWN_INITIAL_MASS_FUNCTION,                         { ERROR_SCOPE::ALWAYS,              "Unknown initial mass function (IMF)" }},
    { ERROR::UNKNOWN_KICK_DIRECTION_DISTRIBUTION,                   { ERROR_SCOPE::ALWAYS,              "Unknown kick direction distribution" }},
    { ERROR::UNKNOWN_KICK_MAGNITUDE_DISTRIBUTION,                   { ERROR_SCOPE::ALWAYS,              "Unknown kick magnitude distribution" }},
    { ERROR::UNKNOWN_LBV_MASS_LOSS_PRESCRIPTION,                    { ERROR_SCOPE::ALWAYS,              "Unknown LBV mass loss prescription" }},
    { ERROR::UNKNOWN_LOGFILE,                                       { ERROR_SCOPE::ALWAYS,              "Unknown log file" }},
    { ERROR::UNKNOWN_MT_CASE,                                       { ERROR_SCOPE::ALWAYS,              "Unknown mass transfer case" }},
    { ERROR::UNKNOWN_MT_ACCRETION_EFFICIENCY_PRESCRIPTION,          { ERROR_SCOPE::ALWAYS,              "Unknown mass transfer accretion efficiency prescription" }},
    { ERROR::UNKNOWN_MT_ANGULAR_MOMENTUM_LOSS_PRESCRIPTION,         { ERROR_SCOPE::ALWAYS,              "Unknown mass transfer angular momentum loss prescription" }},
    { ERROR::UNKNOWN_MT_REJUVENATION_PRESCRIPTION,                  { ERROR_SCOPE::ALWAYS,              "Unknown mass transfer rejuvenation prescription" }},
    { ERROR::UNKNOWN_MT_THERMALLY_LIMITED_VARIATION,                { ERROR_SCOPE::ALWAYS,              "Unknown mass transfer thermally limited variation" }},
    { ERROR::UNKNOWN_MASS_LOSS_PRESCRIPTION,                        { ERROR_SCOPE::ALWAYS,              "Unknown mass loss prescription" }},
    { ERROR::UNKNOWN_NEUTRINO_MASS_LOSS_PRESCRIPTION,               { ERROR_SCOPE::ALWAYS,              "Unknown neutrino mass loss prescription" }},
    { ERROR::UNKNOWN_NS_EOS,                                        { ERROR_SCOPE::ALWAYS,              "Unknown NS equation-of-state" }},
    { ERROR::UNKNOWN_OB_MASS_LOSS_PRESCRIPTION,                     { ERROR_SCOPE::ALWAYS,              "Unknown OB mass loss prescription" }},
    { ERROR::UNKNOWN_PPI_PRESCRIPTION,                              { ERROR_SCOPE::ALWAYS,              "Unknown pulsational pair instability prescription" }},
    { ERROR::UNKNOWN_PROGRAM_OPTION,                                { ERROR_SCOPE::ALWAYS,              "Unknown program option" }},
    { ERROR::UNKNOWN_PROPERTY_TYPE,                                 { ERROR_SCOPE::ALWAYS,              "Unknown property type" }},
    { ERROR::UNKNOWN_PULSAR_BIRTH_MAGNETIC_FIELD_DISTRIBUTION,      { ERROR_SCOPE::ALWAYS,              "Unknown pulsar birth magnetic field distribution" }},
    { ERROR::UNKNOWN_PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION,         { ERROR_SCOPE::ALWAYS,              "Unknown pulsar birth spin period distribution" }},
    { ERROR::UNKNOWN_Q_DISTRIBUTION,                                { ERROR_SCOPE::ALWAYS,              "Unknown q-distribution" }},
    { ERROR::UNKNOWN_QCRIT_PRESCRIPTION,                            { ERROR_SCOPE::ALWAYS,              "Unknown QCRIT prescription" }},
    { ERROR::UNKNOWN_REMNANT_MASS_PRESCRIPTION,                     { ERROR_SCOPE::ALWAYS,              "Unknown remnant mass prescription" }},
    { ERROR::UNKNOWN_RSG_MASS_LOSS_PRESCRIPTION,                    { ERROR_SCOPE::ALWAYS,              "Unknown RSG mass loss prescription" }},
    { ERROR::UNKNOWN_SEMI_MAJOR_AXIS_DISTRIBUTION,                  { ERROR_SCOPE::ALWAYS,              "Unknown semi-major axis distribution" }},
    { ERROR::UNKNOWN_SN_ENGINE,                                     { ERROR_SCOPE::ALWAYS,              "Unknown supernova engine" }},
    { ERROR::UNKNOWN_SN_EVENT,                                      { ERROR_SCOPE::ALWAYS,              "Unknown supernova event" }},
    { ERROR::UNKNOWN_STELLAR_POPULATION,                            { ERROR_SCOPE::ALWAYS,              "Unknown stellar population" }},
    { ERROR::UNKNOWN_STELLAR_PROPERTY,                              { ERROR_SCOPE::ALWAYS,              "Unknown stellar property" }},
    { ERROR::UNKNOWN_STELLAR_TYPE,                                  { ERROR_SCOPE::ALWAYS,              "Unknown stellar type" }},
    { ERROR::UNKNOWN_TIDES_PRESCRIPTION,                            { ERROR_SCOPE::ALWAYS,              "Unknown tides prescription" }},
    { ERROR::UNKNOWN_VMS_MASS_LOSS_PRESCRIPTION,                    { ERROR_SCOPE::ALWAYS,              "Unknown VMS mass loss prescription" }},
    { ERROR::UNKNOWN_VROT_PRESCRIPTION,                             { ERROR_SCOPE::ALWAYS,              "Unknown rotational velocity prescription" }},
    { ERROR::UNKNOWN_ZETA_PRESCRIPTION,                             { ERROR_SCOPE::ALWAYS,              "Unknown stellar ZETA prescription" }},
    { ERROR::WARNING,                                               { ERROR_SCOPE::ALWAYS,              "Warning!" }},
    { ERROR::WHITE_DWARF_TOO_MASSIVE,                               { ERROR_SCOPE::ALWAYS,              "This white dwarf exceeds the Chandrasekhar mass limit" }},
    { ERROR::UNKNOWN_WR_MASS_LOSS_PRESCRIPTION,                     { ERROR_SCOPE::ALWAYS,              "Unknown WR mass loss prescription" }}
};


#endif // __ErrorCatalog_h__
