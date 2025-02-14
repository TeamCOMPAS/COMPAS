#ifndef __LogTypedefs_h__
#define __LogTypedefs_h__


// This is where developer-defined types that pertain directly to the COMPAS looging functionality
// (including the definitions of the default record composition for the various log files) are
// defined.  Non-logging tyedefs are listed in typedefs.h


#include "constants.h"


// logfile file types
enum class LOGFILETYPE: int { HDF5, CSV, TSV, TXT };                        // need this declared here so can declare the constant...
const COMPASUnorderedMap<LOGFILETYPE, std::string> LOGFILETYPELabel = {     // file types
    { LOGFILETYPE::HDF5, "HDF5" },
    { LOGFILETYPE::CSV,  "CSV" },
    { LOGFILETYPE::TSV,  "TSV" },
    { LOGFILETYPE::TXT,  "TXT" }
};
const COMPASUnorderedMap<LOGFILETYPE, std::string> LOGFILETYPEFileExt = {   // file extensions
    { LOGFILETYPE::HDF5, "h5" },
    { LOGFILETYPE::CSV,  "csv" },
    { LOGFILETYPE::TSV,  "tsv" },
    { LOGFILETYPE::TXT,  "txt" }
};


// enum class TYPENAME
// Symbolic names for variable typenames (for printing)
enum class TYPENAME: int {
    NONE,
    BOOL,
    SHORTINT,
    INT,
    LONGINT,
    LONGLONGINT,
    USHORTINT,
    UINT,
    ULONGINT,
    ULONGLONGINT,
    FLOAT,
    DOUBLE,
    LONGDOUBLE,
    STRING,
    OBJECT_ID,
    ERROR,
    STELLAR_TYPE,
    MT_CASE,
    MT_TRACKING,
    MASS_TRANSFER_TIMESCALE,
    SN_EVENT,
    SN_STATE,
    STRING_VECTOR,
    EVOLUTION_STATUS
};
// labels (long and short) for typenames
// unordered_map - key is integer typename (from enum class TYPENAME above)
const COMPASUnorderedMap<TYPENAME, STR_STR> TYPENAME_LABEL = {
    { TYPENAME::NONE,             { "NONE",                   "NONE"           }},
    { TYPENAME::BOOL,             { "BOOL",                   "BOOL"           }},
    { TYPENAME::SHORTINT,         { "SHORT INT",              "INT"            }},
    { TYPENAME::INT,              { "INT",                    "INT"            }},
    { TYPENAME::LONGINT,          { "LONG_INT",               "INT"            }},
    { TYPENAME::LONGLONGINT,      { "LONG_LONG_INT",          "INT"            }},
    { TYPENAME::USHORTINT,        { "UNSIGNED_SHORT_INT",     "INT"            }},
    { TYPENAME::UINT,             { "UNSIGNED_INT",           "INT"            }},
    { TYPENAME::ULONGINT,         { "UNSIGNED_LONG_INT",      "INT"            }},
    { TYPENAME::ULONGLONGINT,     { "UNSIGNED_LONG_LONG_INT", "INT"            }},
    { TYPENAME::FLOAT,            { "FLOAT",                  "FLOAT"          }},
    { TYPENAME::DOUBLE,           { "DOUBLE",                 "FLOAT"          }},
    { TYPENAME::LONGDOUBLE,       { "LONG_DOUBLE",            "FLOAT"          }},
    { TYPENAME::STRING,           { "STRING",                 "STRING"         }},
    { TYPENAME::OBJECT_ID,        { "OBJECT_ID",              "INT"            }},
    { TYPENAME::ERROR,            { "ERROR",                  "INT"            }},
    { TYPENAME::STELLAR_TYPE,     { "STELLAR_TYPE",           "INT"            }},
    { TYPENAME::MT_CASE,          { "MT_CASE",                "INT"            }},
    { TYPENAME::MT_TRACKING,      { "MT_TRACKING",            "INT"            }},
    { TYPENAME::MASS_TRANSFER_TIMESCALE,    { "MASS_TRANSFER_TIMESCALE", "INT" }},
    { TYPENAME::SN_EVENT,         { "SN_EVENT",               "INT"            }},
    { TYPENAME::SN_STATE,         { "SN_STATE",               "INT"            }},
    { TYPENAME::STRING_VECTOR,    { "STRING_VECTOR",          "VECTOR<STRING>" }},
    { TYPENAME::EVOLUTION_STATUS, { "EVOLUTION_STATUS",       "INT"            }}
};


// (convenience) initializer list for INT data types
const std::initializer_list<TYPENAME> INT_TYPES = {
    TYPENAME::SHORTINT,
    TYPENAME::INT,
    TYPENAME::LONGINT,
    TYPENAME::LONGLONGINT,
    TYPENAME::USHORTINT,
    TYPENAME::UINT,
    TYPENAME::ULONGINT,
    TYPENAME::ULONGLONGINT
};

// (convenience) initializer list for FLOAT data types
const std::initializer_list<TYPENAME> FLOAT_TYPES = {
    TYPENAME::FLOAT,
    TYPENAME::DOUBLE,
    TYPENAME::LONGDOUBLE
};


// enum class STRING_QUALIFIER
// Qualifier for typename STRING: FIXED_LENGTH or VARIABLE_LENGTH (used for printing to HDF5 files)
enum class STRING_QUALIFIER: int { NONE, FIXED_LENGTH, VARIABLE_LENGTH };


//  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//  !!!                                                                             !!!
//  !!!   Do not change the following enum classes and maps unless you are adding   !!!
//  !!!   or deleting a new property (or changing the name of an existing property  !!!
//  !!!   for some reason)                                                          !!!
//  !!!                                                                             !!!
//  !!!   The STAR_PROPERTIES #define below defines the STELLAR variables allowed   !!!
//  !!!   for logfile record definition - if a property is not on the list it       !!!
//  !!!   cannot be selected for inclusion in a logfile via the                     !!!
//  !!!   --logfile-definitions option.                                             !!!                                            
//  !!!                                                                             !!!
//  !!!   The enum classes STAR_PROPERTY, BINARY_PROPERTY, and PROGRAM_OPTION       !!!
//  !!!   defined below defines the STELLAR and BINARY properties, and the          !!!
//  !!!   PROGRAM_OPTIONs allowed for logfile record definition - if a property is  !!!
//  !!!   not on those lists it cannot be selected for inclusion in a logfile via   !!!
//  !!!   --logfile-definitions option.                                             !!!                                            
//  !!!                                                                             !!!
//  !!!   *NOTE*                                                                    !!!
//  !!!   The following enum classes anad maps are not where header strings should  !!!
//  !!!   be changed!  These classes and maps are a lookup facility for the logfile !!!
//  !!!   definitions file parser.                                                  !!!
//  !!!                                                                             !!!
//  !!!   Header strings are in the following maps, and should be changed there:    !!!
//  !!!                                                                             !!!
//  !!!   std::map<ANY_STAR_PROPERTY, PROPERTY_DETAILS> ANY_STAR_PROPERTY_DETAIL    !!!
//  !!!   std::map<ANY_STAR_PROPERTY, PROPERTY_DETAILS> BINARY_PROPERTY_DETAIL      !!!
//  !!!                                                                             !!!
//  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

// The #define below defines the STELLAR variables allowed for logfile record definition
// Add new property names here to make the property available for specification in a logfile
// #define is used so that the same list of variables can be used for the various stellar property enum classes (see below)
#define STAR_PROPERTIES                              \
    NONE,                                            \
    AGE,                                             \
    ANGULAR_MOMENTUM,                                \
    BINDING_ENERGY_AT_COMMON_ENVELOPE,               \
    BINDING_ENERGY_FIXED,                            \
    BINDING_ENERGY_NANJING,                          \
    BINDING_ENERGY_PRE_COMMON_ENVELOPE,              \
    BINDING_ENERGY_LOVERIDGE,                        \
    BINDING_ENERGY_LOVERIDGE_WINDS,                  \
    BINDING_ENERGY_KRUCKOW,                          \
    CHEMICALLY_HOMOGENEOUS_MAIN_SEQUENCE,            \
    CO_CORE_MASS,                                    \
    CO_CORE_MASS_AT_COMMON_ENVELOPE,                 \
    CO_CORE_MASS_AT_COMPACT_OBJECT_FORMATION,        \
    CORE_MASS,                                       \
    CORE_MASS_AT_COMMON_ENVELOPE,                    \
    CORE_MASS_AT_COMPACT_OBJECT_FORMATION,           \
    DRAWN_KICK_MAGNITUDE,                            \
    DOMINANT_MASS_LOSS_RATE,                         \
    DT,                                              \
    DYNAMICAL_TIMESCALE,                             \
    DYNAMICAL_TIMESCALE_POST_COMMON_ENVELOPE,        \
    DYNAMICAL_TIMESCALE_PRE_COMMON_ENVELOPE,         \
    ECCENTRIC_ANOMALY,                               \
    ENV_MASS,                                        \
    ERROR,                                           \
    EVOL_STATUS,                                     \
    EXPERIENCED_AIC,                                 \
    EXPERIENCED_CCSN,                                \
    EXPERIENCED_HeSD,                                \
    EXPERIENCED_ECSN,                                \
    EXPERIENCED_PISN,                                \
    EXPERIENCED_PPISN,                               \
    EXPERIENCED_RLOF,                                \
    EXPERIENCED_SNIA,                                \
    EXPERIENCED_SN_TYPE,                             \
    EXPERIENCED_USSN,                                \
    FALLBACK_FRACTION,                               \
    HE_CORE_MASS,                                    \
    HE_CORE_MASS_AT_COMMON_ENVELOPE,                 \
    HE_CORE_MASS_AT_COMPACT_OBJECT_FORMATION,        \
    HELIUM_ABUNDANCE_CORE,                           \
    HELIUM_ABUNDANCE_SURFACE,                        \
    HYDROGEN_ABUNDANCE_CORE,                         \
    HYDROGEN_ABUNDANCE_SURFACE,                      \
    ID,                                              \
    INITIAL_HELIUM_ABUNDANCE,                        \
    INITIAL_HYDROGEN_ABUNDANCE,                      \
    INITIAL_STELLAR_TYPE,                            \
    INITIAL_STELLAR_TYPE_NAME,                       \
    IS_AIC,                                          \
    IS_CCSN,                                         \
    IS_HeSD,                                         \
    IS_ECSN,                                         \
    IS_HYDROGEN_POOR,                                \
    IS_PISN,                                         \
    IS_PPISN,                                        \
    IS_RLOF,                                         \
    IS_SNIA,                                         \
    IS_USSN,                                         \
    KICK_MAGNITUDE,                                  \
    LAMBDA_AT_COMMON_ENVELOPE,                       \
    LAMBDA_DEWI,                                     \
    LAMBDA_FIXED,                                    \
    LAMBDA_KRUCKOW,                                  \
    LAMBDA_KRUCKOW_BOTTOM,                           \
    LAMBDA_KRUCKOW_MIDDLE,                           \
    LAMBDA_KRUCKOW_TOP,                              \
    LAMBDA_LOVERIDGE,                                \
    LAMBDA_LOVERIDGE_WINDS,                          \
    LAMBDA_NANJING,                                  \
    LBV_PHASE_FLAG,                                  \
    LUMINOSITY,                                      \
    LUMINOSITY_POST_COMMON_ENVELOPE,                 \
    LUMINOSITY_PRE_COMMON_ENVELOPE,                  \
    MASS,                                            \
    MASS_0,                                          \
    MASS_LOSS_DIFF,                                  \
    MASS_TRANSFER_DIFF,                              \
    MASS_TRANSFER_DONOR_HISTORY,                     \
    MDOT,                                            \
    MEAN_ANOMALY,                                    \
    METALLICITY,                                     \
    MOMENT_OF_INERTIA,                               \
    MZAMS,                                           \
    OMEGA,                                           \
    OMEGA_BREAK,                                     \
    OMEGA_ZAMS,                                      \
    ORBITAL_ENERGY_POST_SUPERNOVA,                   \
    ORBITAL_ENERGY_PRE_SUPERNOVA,                    \
    PULSAR_MAGNETIC_FIELD,                           \
    PULSAR_SPIN_DOWN_RATE,                           \
    PULSAR_BIRTH_PERIOD,                             \
    PULSAR_BIRTH_SPIN_DOWN_RATE,                     \
    PULSAR_SPIN_FREQUENCY,                           \
    PULSAR_SPIN_PERIOD,                              \
    RADIAL_EXPANSION_TIMESCALE,                      \
    RADIAL_EXPANSION_TIMESCALE_POST_COMMON_ENVELOPE, \
    RADIAL_EXPANSION_TIMESCALE_PRE_COMMON_ENVELOPE,  \
    RADIUS,                                          \
    RANDOM_SEED,                                     \
    RECYCLED_NEUTRON_STAR,                           \
    RLOF_ONTO_NS,                                    \
    ROCKET_KICK_MAGNITUDE,                           \
    ROCKET_KICK_PHI,                                 \
    ROCKET_KICK_THETA,                               \
    RZAMS,                                           \
    SN_TYPE,                                         \
    SPEED,                                           \
    STELLAR_TYPE,                                    \
    STELLAR_TYPE_NAME,                               \
    STELLAR_TYPE_PREV,                               \
    STELLAR_TYPE_PREV_NAME,                          \
    SUPERNOVA_KICK_MAGNITUDE_RANDOM_NUMBER,          \
    SUPERNOVA_PHI,                                   \
    SUPERNOVA_THETA,                                 \
    TEMPERATURE,                                     \
    TEMPERATURE_POST_COMMON_ENVELOPE,                \
    TEMPERATURE_PRE_COMMON_ENVELOPE,                 \
    THERMAL_TIMESCALE,                               \
    THERMAL_TIMESCALE_POST_COMMON_ENVELOPE,          \
    THERMAL_TIMESCALE_PRE_COMMON_ENVELOPE,           \
    TIME,                                            \
    TIMESCALE_MS,                                    \
    TOTAL_MASS_AT_COMPACT_OBJECT_FORMATION,          \
    TRUE_ANOMALY,                                    \
    TZAMS,                                           \
    ZETA_HURLEY,                                     \
    ZETA_HURLEY_HE,                                  \
    ZETA_SOBERMAN,                                   \
    ZETA_SOBERMAN_HE


// enum class STAR_PROPERTY
// Symbolic names for variables of an individual star that can be selected for printing
// STAR_PROPERTY refers to an individual star of type BaseStar for SSE (differences are where the data comes from, and the column header)
enum class STAR_PROPERTY: int { STAR_PROPERTIES };

// map STAR PROPERTY to string identifying the property
// for lookup by the printing functions
// this map serves as the lookup for: STAR_PROPERTY, STAR_1_PROPERTY, STAR_2_PROPERTY, SUPERNOVA_PROPERTY, COMPANION_PROPERTY and ANY_STAR_PROPERTY
//
// Properties only need to be here if they are required to be available for printing in the logfiles.
// All keys present here should also be in the STAR_PROPERTIES #define and ANY_STAR_PROPERTY_DETAIL
const COMPASUnorderedMap<STAR_PROPERTY, std::string> STAR_PROPERTY_LABEL = {
    { STAR_PROPERTY::NONE,                                            "NONE" },
    { STAR_PROPERTY::AGE,                                             "AGE" },
    { STAR_PROPERTY::ANGULAR_MOMENTUM,                                "ANGULAR_MOMENTUM" },
    { STAR_PROPERTY::BINDING_ENERGY_AT_COMMON_ENVELOPE,               "BINDING_ENERGY_AT_COMMON_ENVELOPE" },
    { STAR_PROPERTY::BINDING_ENERGY_FIXED,                            "BINDING_ENERGY_FIXED" },
    { STAR_PROPERTY::BINDING_ENERGY_NANJING,                          "BINDING_ENERGY_NANJING" },
    { STAR_PROPERTY::BINDING_ENERGY_PRE_COMMON_ENVELOPE,              "BINDING_ENERGY_PRE_COMMON_ENVELOPE" },
    { STAR_PROPERTY::BINDING_ENERGY_LOVERIDGE,                        "BINDING_ENERGY_LOVERIDGE" },
    { STAR_PROPERTY::BINDING_ENERGY_LOVERIDGE_WINDS,                  "BINDING_ENERGY_LOVERIDGE_WINDS" },
    { STAR_PROPERTY::BINDING_ENERGY_KRUCKOW,                          "BINDING_ENERGY_KRUCKOW" },
    { STAR_PROPERTY::CHEMICALLY_HOMOGENEOUS_MAIN_SEQUENCE,            "CHEMICALLY_HOMOGENEOUS_MAIN_SEQUENCE" },
    { STAR_PROPERTY::CO_CORE_MASS,                                    "CO_CORE_MASS" },
    { STAR_PROPERTY::CO_CORE_MASS_AT_COMMON_ENVELOPE,                 "CO_CORE_MASS_AT_COMMON_ENVELOPE" },
    { STAR_PROPERTY::CO_CORE_MASS_AT_COMPACT_OBJECT_FORMATION,        "CO_CORE_MASS_AT_COMPACT_OBJECT_FORMATION" },
    { STAR_PROPERTY::CORE_MASS,                                       "CORE_MASS" },
    { STAR_PROPERTY::CORE_MASS_AT_COMMON_ENVELOPE,                    "CORE_MASS_AT_COMMON_ENVELOPE" },
    { STAR_PROPERTY::CORE_MASS_AT_COMPACT_OBJECT_FORMATION,           "CORE_MASS_AT_COMPACT_OBJECT_FORMATION" },
    { STAR_PROPERTY::DRAWN_KICK_MAGNITUDE,                            "DRAWN_KICK_MAGNITUDE" },
    { STAR_PROPERTY::DOMINANT_MASS_LOSS_RATE,                         "DOMINANT_MASS_LOSS_RATE"},
    { STAR_PROPERTY::DT,                                              "DT" },
    { STAR_PROPERTY::DYNAMICAL_TIMESCALE,                             "DYNAMICAL_TIMESCALE" },
    { STAR_PROPERTY::DYNAMICAL_TIMESCALE_POST_COMMON_ENVELOPE,        "DYNAMICAL_TIMESCALE_POST_COMMON_ENVELOPE" },
    { STAR_PROPERTY::DYNAMICAL_TIMESCALE_PRE_COMMON_ENVELOPE,         "DYNAMICAL_TIMESCALE_PRE_COMMON_ENVELOPE" },
    { STAR_PROPERTY::ECCENTRIC_ANOMALY,                               "ECCENTRIC_ANOMALY" },
    { STAR_PROPERTY::ENV_MASS,                                        "ENV_MASS" },
    { STAR_PROPERTY::ERROR,                                           "ERROR" },
    { STAR_PROPERTY::EVOL_STATUS,                                     "EVOL_STATUS" },
    { STAR_PROPERTY::EXPERIENCED_AIC,                                 "EXPERIENCED_AIC" },
    { STAR_PROPERTY::EXPERIENCED_CCSN,                                "EXPERIENCED_CCSN" },
    { STAR_PROPERTY::EXPERIENCED_HeSD,                                "EXPERIENCED_HeSD" },
    { STAR_PROPERTY::EXPERIENCED_ECSN,                                "EXPERIENCED_ECSN" },
    { STAR_PROPERTY::EXPERIENCED_PISN,                                "EXPERIENCED_PISN" },
    { STAR_PROPERTY::EXPERIENCED_PPISN,                               "EXPERIENCED_PPISN" },
    { STAR_PROPERTY::EXPERIENCED_RLOF,                                "EXPERIENCED_RLOF" },
    { STAR_PROPERTY::EXPERIENCED_SNIA,                                "EXPERIENCED_SNIA" },
    { STAR_PROPERTY::EXPERIENCED_SN_TYPE,                             "EXPERIENCED_SN_TYPE" },
    { STAR_PROPERTY::EXPERIENCED_USSN,                                "EXPERIENCED_USSN" },
    { STAR_PROPERTY::FALLBACK_FRACTION,                               "FALLBACK_FRACTION" },
    { STAR_PROPERTY::HE_CORE_MASS,                                    "HE_CORE_MASS" },
    { STAR_PROPERTY::HE_CORE_MASS_AT_COMMON_ENVELOPE,                 "HE_CORE_MASS_AT_COMMON_ENVELOPE" },
    { STAR_PROPERTY::HE_CORE_MASS_AT_COMPACT_OBJECT_FORMATION,        "HE_CORE_MASS_AT_COMPACT_OBJECT_FORMATION" },
    { STAR_PROPERTY::HELIUM_ABUNDANCE_CORE,                           "HELIUM_ABUNDANCE_CORE" },
    { STAR_PROPERTY::HELIUM_ABUNDANCE_SURFACE,                        "HELIUM_ABUNDANCE_SURFACE" },
    { STAR_PROPERTY::HYDROGEN_ABUNDANCE_CORE,                         "HYDROGEN_ABUNDANCE_CORE" },
    { STAR_PROPERTY::HYDROGEN_ABUNDANCE_SURFACE,                      "HYDROGEN_ABUNDANCE_SURFACE" },
    { STAR_PROPERTY::ID,                                              "ID" },
    { STAR_PROPERTY::INITIAL_HELIUM_ABUNDANCE,                        "INITIAL_HELIUM_ABUNDANCE" },
    { STAR_PROPERTY::INITIAL_HYDROGEN_ABUNDANCE,                      "INITIAL_HYDROGEN_ABUNDANCE" },
    { STAR_PROPERTY::INITIAL_STELLAR_TYPE,                            "INITIAL_STELLAR_TYPE" },
    { STAR_PROPERTY::INITIAL_STELLAR_TYPE,                            "INITIAL_STELLAR_TYPE" },
    { STAR_PROPERTY::INITIAL_STELLAR_TYPE_NAME,                       "INITIAL_STELLAR_TYPE_NAME" },
    { STAR_PROPERTY::IS_AIC,                                          "IS_AIC" },
    { STAR_PROPERTY::IS_CCSN,                                         "IS_CCSN" },
    { STAR_PROPERTY::IS_HeSD,                                         "IS_HeSD" },
    { STAR_PROPERTY::IS_ECSN,                                         "IS_ECSN" },
    { STAR_PROPERTY::IS_HYDROGEN_POOR,                                "IS_HYDROGEN_POOR" },
    { STAR_PROPERTY::IS_PISN,                                         "IS_PISN" },
    { STAR_PROPERTY::IS_PPISN,                                        "IS_PPISN" },
    { STAR_PROPERTY::IS_RLOF,                                         "IS_RLOF" },
    { STAR_PROPERTY::IS_SNIA,                                         "IS_SNIA" },
    { STAR_PROPERTY::IS_USSN,                                         "IS_USSN" },
    { STAR_PROPERTY::KICK_MAGNITUDE,                                  "KICK_MAGNITUDE" },
    { STAR_PROPERTY::LAMBDA_AT_COMMON_ENVELOPE,                       "LAMBDA_AT_COMMON_ENVELOPE" },
    { STAR_PROPERTY::LAMBDA_DEWI,                                     "LAMBDA_DEWI" },
    { STAR_PROPERTY::LAMBDA_FIXED,                                    "LAMBDA_FIXED" },
    { STAR_PROPERTY::LAMBDA_KRUCKOW,                                  "LAMBDA_KRUCKOW" },
    { STAR_PROPERTY::LAMBDA_KRUCKOW_BOTTOM,                           "LAMBDA_KRUCKOW_BOTTOM" },
    { STAR_PROPERTY::LAMBDA_KRUCKOW_MIDDLE,                           "LAMBDA_KRUCKOW_MIDDLE" },
    { STAR_PROPERTY::LAMBDA_KRUCKOW_TOP,                              "LAMBDA_KRUCKOW_TOP" },
    { STAR_PROPERTY::LAMBDA_LOVERIDGE,                                "LAMBDA_LOVERIDGE" },
    { STAR_PROPERTY::LAMBDA_LOVERIDGE_WINDS,                          "LAMBDA_LOVERIDGE_WINDS" },
    { STAR_PROPERTY::LAMBDA_NANJING,                                  "LAMBDA_NANJING" },
    { STAR_PROPERTY::LBV_PHASE_FLAG,                                  "LBV_PHASE_FLAG" },
    { STAR_PROPERTY::LUMINOSITY,                                      "LUMINOSITY" },
    { STAR_PROPERTY::LUMINOSITY_POST_COMMON_ENVELOPE,                 "LUMINOSITY_POST_COMMON_ENVELOPE" },
    { STAR_PROPERTY::LUMINOSITY_PRE_COMMON_ENVELOPE,                  "LUMINOSITY_PRE_COMMON_ENVELOPE" },
    { STAR_PROPERTY::MASS,                                            "MASS" },
    { STAR_PROPERTY::MASS_0,                                          "MASS_0" },
    { STAR_PROPERTY::MASS_LOSS_DIFF,                                  "MASS_LOSS_DIFF" },
    { STAR_PROPERTY::MASS_TRANSFER_DIFF,                              "MASS_TRANSFER_DIFF" },
    { STAR_PROPERTY::MASS_TRANSFER_DONOR_HISTORY,                     "MASS_TRANSFER_DONOR_HISTORY" },
    { STAR_PROPERTY::MDOT,                                            "MDOT" },
    { STAR_PROPERTY::MEAN_ANOMALY,                                    "MEAN_ANOMALY" },
    { STAR_PROPERTY::METALLICITY,                                     "METALLICITY" },
    { STAR_PROPERTY::MOMENT_OF_INERTIA,                               "MOMENT_OF_INERTIA"},
    { STAR_PROPERTY::MZAMS,                                           "MZAMS" },
    { STAR_PROPERTY::OMEGA,                                           "OMEGA" },
    { STAR_PROPERTY::OMEGA_BREAK,                                     "OMEGA_BREAK" },
    { STAR_PROPERTY::OMEGA_ZAMS,                                      "OMEGA_ZAMS" },
    { STAR_PROPERTY::ORBITAL_ENERGY_POST_SUPERNOVA,                   "ORBITAL_ENERGY_POST_SUPERNOVA" },
    { STAR_PROPERTY::ORBITAL_ENERGY_PRE_SUPERNOVA,                    "ORBITAL_ENERGY_PRE_SUPERNOVA" },
    { STAR_PROPERTY::PULSAR_MAGNETIC_FIELD,                           "PULSAR_MAGNETIC_FIELD" },
    { STAR_PROPERTY::PULSAR_SPIN_DOWN_RATE,                           "PULSAR_SPIN_DOWN_RATE" },
    { STAR_PROPERTY::PULSAR_BIRTH_PERIOD,                             "PULSAR_BIRTH_PERIOD" },
    { STAR_PROPERTY::PULSAR_BIRTH_SPIN_DOWN_RATE,                     "PULSAR_BIRTH_SPIN_DOWN_RATE" },
    { STAR_PROPERTY::PULSAR_SPIN_FREQUENCY,                           "PULSAR_SPIN_FREQUENCY" },
    { STAR_PROPERTY::PULSAR_SPIN_PERIOD,                              "PULSAR_SPIN_PERIOD" },
    { STAR_PROPERTY::RADIAL_EXPANSION_TIMESCALE,                      "RADIAL_EXPANSION_TIMESCALE" },
    { STAR_PROPERTY::RADIAL_EXPANSION_TIMESCALE_POST_COMMON_ENVELOPE, "RADIAL_EXPANSION_TIMESCALE_POST_COMMON_ENVELOPE" },
    { STAR_PROPERTY::RADIAL_EXPANSION_TIMESCALE_PRE_COMMON_ENVELOPE,  "RADIAL_EXPANSION_TIMESCALE_PRE_COMMON_ENVELOPE" },
    { STAR_PROPERTY::RADIUS,                                          "RADIUS" },
    { STAR_PROPERTY::RANDOM_SEED,                                     "RANDOM_SEED" },
    { STAR_PROPERTY::RECYCLED_NEUTRON_STAR,                           "RECYCLED_NEUTRON_STAR" },
    { STAR_PROPERTY::RLOF_ONTO_NS,                                    "RLOF_ONTO_NS" },
    { STAR_PROPERTY::ROCKET_KICK_MAGNITUDE,                           "ROCKET_KICK_MAGNITUDE" },
    { STAR_PROPERTY::ROCKET_KICK_PHI,                                 "ROCKET_KICK_PHI" },
    { STAR_PROPERTY::ROCKET_KICK_THETA,                               "ROCKET_KICK_THETA" },
    { STAR_PROPERTY::RZAMS,                                           "RZAMS" },
    { STAR_PROPERTY::SN_TYPE,                                         "SN_TYPE" },
    { STAR_PROPERTY::SPEED,                                           "SPEED" },
    { STAR_PROPERTY::STELLAR_TYPE,                                    "STELLAR_TYPE" },
    { STAR_PROPERTY::STELLAR_TYPE_NAME,                               "STELLAR_TYPE_NAME" },
    { STAR_PROPERTY::STELLAR_TYPE_PREV,                               "STELLAR_TYPE_PREV" },
    { STAR_PROPERTY::STELLAR_TYPE_PREV_NAME,                          "STELLAR_TYPE_PREV_NAME" },
    { STAR_PROPERTY::SUPERNOVA_KICK_MAGNITUDE_RANDOM_NUMBER,          "SUPERNOVA_KICK_MAGNITUDE_RANDOM_NUMBER" },
    { STAR_PROPERTY::SUPERNOVA_PHI,                                   "SUPERNOVA_PHI" },
    { STAR_PROPERTY::SUPERNOVA_THETA,                                 "SUPERNOVA_THETA" },
    { STAR_PROPERTY::TEMPERATURE,                                     "TEMPERATURE" },
    { STAR_PROPERTY::TEMPERATURE_POST_COMMON_ENVELOPE,                "TEMPERATURE_POST_COMMON_ENVELOPE" },
    { STAR_PROPERTY::TEMPERATURE_PRE_COMMON_ENVELOPE,                 "TEMPERATURE_PRE_COMMON_ENVELOPE" },
    { STAR_PROPERTY::THERMAL_TIMESCALE,                               "THERMAL_TIMESCALE" },
    { STAR_PROPERTY::THERMAL_TIMESCALE_POST_COMMON_ENVELOPE,          "THERMAL_TIMESCALE_POST_COMMON_ENVELOPE" },
    { STAR_PROPERTY::THERMAL_TIMESCALE_PRE_COMMON_ENVELOPE,           "THERMAL_TIMESCALE_PRE_COMMON_ENVELOPE" },
    { STAR_PROPERTY::TIME,                                            "TIME" },
    { STAR_PROPERTY::TIMESCALE_MS,                                    "TIMESCALE_MS" },
    { STAR_PROPERTY::TOTAL_MASS_AT_COMPACT_OBJECT_FORMATION,          "TOTAL_MASS_AT_COMPACT_OBJECT_FORMATION" },
    { STAR_PROPERTY::TRUE_ANOMALY,                                    "TRUE_ANOMALY" },
    { STAR_PROPERTY::TZAMS,                                           "TZAMS" },
    { STAR_PROPERTY::ZETA_HURLEY,                                     "ZETA_HURLEY" },
    { STAR_PROPERTY::ZETA_HURLEY_HE,                                  "ZETA_HURLEY_HE" },
    { STAR_PROPERTY::ZETA_SOBERMAN,                                   "ZETA_SOBERMAN" },
    { STAR_PROPERTY::ZETA_SOBERMAN_HE,                                "ZETA_SOBERMAN_HE" }
};


// enum class STAR_1_PROPERTY
// Symbolic names for variables of an individual star that can be selected for printing
// STAR_1_PROPERTY refers to star 1 (of type BinaryConstituentStar) of a binary for BSE (differences are where the data comes from, and the column header)
enum class STAR_1_PROPERTY: int { STAR_PROPERTIES };


// enum class STAR_2_PROPERTY
// Symbolic names for variables of an individual star that can be selected for printing
// STAR_2_PROPERTY refers to star 2 (of type BinaryConstituentStar) of a binary for BSE (differences are where the data comes from, and the column header)
enum class STAR_2_PROPERTY: int { STAR_PROPERTIES };


// enum class SUPERNOVA_PROPERTY
// Symbolic names for variables of an individual star that can be selected for printing
// SUPERNOVA_PROPERTY refers to the supernova star (of type BinaryConstituentStar) of a binary where one star has experienced a SN event for BSE (differences are where the data comes from, and the column header)
enum class SUPERNOVA_PROPERTY: int { STAR_PROPERTIES };


// enum class COMPANION_PROPERTY
// Symbolic names for variables of an individual star that can be selected for printing
// COMPANION_PROPERTY refers to the companion star (of type BinaryConstituentStar) of a binary where one star has experienced a SN event for BSE (differences are where the data comes from, and the column header)
enum class COMPANION_PROPERTY: int { STAR_PROPERTIES };


// enum class ANY_STAR_PROPERTY
// Symbolic names for variables of an individual star that can be selected for printing
// ANY_STAR_PROPERTY refers to the any individual star
enum class ANY_STAR_PROPERTY: int { STAR_PROPERTIES };


// enum class BINARY_PROPERTY
// Symbolic names for variables of binary stars that can be selected for printing
// BINARY_PROPERTY refers to a binary star of type BaseBinaryStar) for BSE
//
// Properties only need to be here if they are required to be available for 
// printing in the logfiles - all keys present here should also be in BINARY_PROPERTY_LABEL
// and BINARY_PROPERTY_DETAIL
enum class BINARY_PROPERTY: int {
    NONE,
    CIRCULARIZATION_TIMESCALE,
    COMMON_ENVELOPE_AT_LEAST_ONCE,
    COMMON_ENVELOPE_EVENT_COUNT,
    UNBOUND,
    DOUBLE_CORE_COMMON_ENVELOPE,
    DT,
    ECCENTRICITY,
    ECCENTRICITY_AT_DCO_FORMATION,
    ECCENTRICITY_INITIAL,
    ECCENTRICITY_POST_COMMON_ENVELOPE,
    ECCENTRICITY_PRE_SUPERNOVA,
    ECCENTRICITY_PRE_COMMON_ENVELOPE,
    ERROR,
    EVOL_STATUS,
    ID,
    IMMEDIATE_RLOF_POST_COMMON_ENVELOPE,
    MASS_1_POST_COMMON_ENVELOPE,
    MASS_1_PRE_COMMON_ENVELOPE,
    MASS_2_POST_COMMON_ENVELOPE,
    MASS_2_PRE_COMMON_ENVELOPE,
    MASS_ENV_1,
    MASS_ENV_2,
    MASSES_EQUILIBRATED,
    MASSES_EQUILIBRATED_AT_BIRTH,
    MASS_TRANSFER_TRACKER_HISTORY,
    MERGES_IN_HUBBLE_TIME,
    OPTIMISTIC_COMMON_ENVELOPE,
    ORBITAL_ANGULAR_VELOCITY,
    ORBITAL_VELOCITY_PRE_SUPERNOVA,
    RADIUS_1_POST_COMMON_ENVELOPE,
    RADIUS_1_PRE_COMMON_ENVELOPE,
    RADIUS_2_POST_COMMON_ENVELOPE,
    RADIUS_2_PRE_COMMON_ENVELOPE,
    RANDOM_SEED,
    RLOF_ACCRETION_EFFICIENCY,
    RLOF_MASS_LOSS_RATE,
    RLOF_MASS_TRANSFER_TIMESCALE,
    RLOF_POST_MT_COMMON_ENVELOPE,
    RLOF_POST_MT_ECCENTRICITY,
    RLOF_POST_MT_EVENT_COUNTER,
    RLOF_POST_MT_ID,
    RLOF_POST_MT_SEMI_MAJOR_AXIS,
    RLOF_POST_MT_STAR1_MASS,
    RLOF_POST_MT_STAR2_MASS,
    RLOF_POST_MT_STAR1_RADIUS,
    RLOF_POST_MT_STAR2_RADIUS,
    RLOF_POST_MT_STAR1_RLOF,
    RLOF_POST_MT_STAR2_RLOF,
    RLOF_POST_MT_STAR1_STELLAR_TYPE,
    RLOF_POST_MT_STAR1_STELLAR_TYPE_NAME,
    RLOF_POST_MT_STAR2_STELLAR_TYPE,
    RLOF_POST_MT_STAR2_STELLAR_TYPE_NAME,
    RLOF_POST_STEP_STAR_TO_ROCHE_LOBE_RADIUS_RATIO_1,
    RLOF_POST_STEP_STAR_TO_ROCHE_LOBE_RADIUS_RATIO_2,
    RLOF_PRE_MT_ECCENTRICITY,
    RLOF_PRE_MT_SEMI_MAJOR_AXIS,
    RLOF_PRE_MT_STAR1_MASS,
    RLOF_PRE_MT_STAR2_MASS,
    RLOF_PRE_MT_STAR1_RADIUS,
    RLOF_PRE_MT_STAR2_RADIUS,
    RLOF_PRE_MT_STAR1_RLOF,
    RLOF_PRE_MT_STAR2_RLOF,
    RLOF_PRE_MT_STAR1_STELLAR_TYPE,
    RLOF_PRE_MT_STAR1_STELLAR_TYPE_NAME,
    RLOF_PRE_MT_STAR2_STELLAR_TYPE,
    RLOF_PRE_MT_STAR2_STELLAR_TYPE_NAME,
    RLOF_PRE_STEP_STAR_TO_ROCHE_LOBE_RADIUS_RATIO_1,
    RLOF_PRE_STEP_STAR_TO_ROCHE_LOBE_RADIUS_RATIO_2,
    RLOF_SECONDARY_POST_COMMON_ENVELOPE,
    RLOF_TIME_POST_MT,
    RLOF_TIME_PRE_MT,
    ROCHE_LOBE_RADIUS_1,
    ROCHE_LOBE_RADIUS_2,
    ROCHE_LOBE_RADIUS_1_POST_COMMON_ENVELOPE,
    ROCHE_LOBE_RADIUS_2_POST_COMMON_ENVELOPE,
    ROCHE_LOBE_RADIUS_1_PRE_COMMON_ENVELOPE,
    ROCHE_LOBE_RADIUS_2_PRE_COMMON_ENVELOPE,
    STAR_TO_ROCHE_LOBE_RADIUS_RATIO_1,
    STAR_TO_ROCHE_LOBE_RADIUS_RATIO_2,
    SEMI_MAJOR_AXIS_AT_DCO_FORMATION,
    SEMI_MAJOR_AXIS_INITIAL,
    SEMI_MAJOR_AXIS_POST_COMMON_ENVELOPE,
    SEMI_MAJOR_AXIS_PRE_SUPERNOVA,
    SEMI_MAJOR_AXIS_PRE_SUPERNOVA_RSOL,
    SEMI_MAJOR_AXIS_PRE_COMMON_ENVELOPE,
    SEMI_MAJOR_AXIS,
    SEMI_MAJOR_AXIS_RSOL,
    SIMULTANEOUS_RLOF,
    STABLE_RLOF_POST_COMMON_ENVELOPE,
    STELLAR_MERGER,
    STELLAR_MERGER_AT_BIRTH,
    STELLAR_TYPE_1_POST_COMMON_ENVELOPE,
    STELLAR_TYPE_1_PRE_COMMON_ENVELOPE,
    STELLAR_TYPE_2_POST_COMMON_ENVELOPE,
    STELLAR_TYPE_2_PRE_COMMON_ENVELOPE,
    STELLAR_TYPE_NAME_1_POST_COMMON_ENVELOPE,
    STELLAR_TYPE_NAME_1_PRE_COMMON_ENVELOPE,
    STELLAR_TYPE_NAME_2_POST_COMMON_ENVELOPE,
    STELLAR_TYPE_NAME_2_PRE_COMMON_ENVELOPE,
    SUPERNOVA_ORBIT_INCLINATION_ANGLE,
    SUPERNOVA_ORBIT_INCLINATION_VECTOR_X,
    SUPERNOVA_ORBIT_INCLINATION_VECTOR_Y,
    SUPERNOVA_ORBIT_INCLINATION_VECTOR_Z,
    SUPERNOVA_STATE,
    SYNCHRONIZATION_TIMESCALE,
    SYSTEMIC_SPEED,
    SYSTEMIC_VELOCITY_X,
    SYSTEMIC_VELOCITY_Y,
    SYSTEMIC_VELOCITY_Z,
    TIME,
    TIME_TO_COALESCENCE,
    TOTAL_ANGULAR_MOMENTUM,
    TOTAL_ENERGY,
    ZETA_LOBE,
    ZETA_STAR
};
// map BINARY_PROPERTY to string identifying the property
// for lookup by the printing functions
//
// Property names only need to be here if they are required to be available for 
// printing in the logfiles - all keys present here should also be in BINARY_PROPERTY
// and BINARY_PROPERTY_DETAIL
const COMPASUnorderedMap<BINARY_PROPERTY, std::string> BINARY_PROPERTY_LABEL = {
    { BINARY_PROPERTY::NONE,                                               "NONE" },
    { BINARY_PROPERTY::CIRCULARIZATION_TIMESCALE,                          "CIRCULARIZATION_TIMESCALE" },
    { BINARY_PROPERTY::COMMON_ENVELOPE_AT_LEAST_ONCE,                      "COMMON_ENVELOPE_AT_LEAST_ONCE" },
    { BINARY_PROPERTY::COMMON_ENVELOPE_EVENT_COUNT,                        "COMMON_ENVELOPE_EVENT_COUNT" },
    { BINARY_PROPERTY::UNBOUND,                                            "UNBOUND" },
    { BINARY_PROPERTY::DOUBLE_CORE_COMMON_ENVELOPE,                        "DOUBLE_CORE_COMMON_ENVELOPE" },
    { BINARY_PROPERTY::DT,                                                 "DT" },
    { BINARY_PROPERTY::ECCENTRICITY,                                       "ECCENTRICITY" },
    { BINARY_PROPERTY::ECCENTRICITY_AT_DCO_FORMATION,                      "ECCENTRICITY_AT_DCO_FORMATION" },
    { BINARY_PROPERTY::ECCENTRICITY_INITIAL,                               "ECCENTRICITY_INITIAL" },
    { BINARY_PROPERTY::ECCENTRICITY_POST_COMMON_ENVELOPE,                  "ECCENTRICITY_POST_COMMON_ENVELOPE" },
    { BINARY_PROPERTY::ECCENTRICITY_PRE_SUPERNOVA,                         "ECCENTRICITY_PRE_SUPERNOVA" },
    { BINARY_PROPERTY::ECCENTRICITY_PRE_COMMON_ENVELOPE,                   "ECCENTRICITY_PRE_COMMON_ENVELOPE" },
    { BINARY_PROPERTY::ERROR,                                              "ERROR" },
    { BINARY_PROPERTY::EVOL_STATUS,                                        "EVOL_STATUS" },
    { BINARY_PROPERTY::ID,                                                 "ID" },
    { BINARY_PROPERTY::IMMEDIATE_RLOF_POST_COMMON_ENVELOPE,                "IMMEDIATE_RLOF_POST_COMMON_ENVELOPE" },
    { BINARY_PROPERTY::MASS_1_POST_COMMON_ENVELOPE,                        "MASS_1_POST_COMMON_ENVELOPE" },
    { BINARY_PROPERTY::MASS_1_PRE_COMMON_ENVELOPE,                         "MASS_1_PRE_COMMON_ENVELOPE" },
    { BINARY_PROPERTY::MASS_2_POST_COMMON_ENVELOPE,                        "MASS_2_POST_COMMON_ENVELOPE" },
    { BINARY_PROPERTY::MASS_2_PRE_COMMON_ENVELOPE,                         "MASS_2_PRE_COMMON_ENVELOPE" },
    { BINARY_PROPERTY::MASS_ENV_1,                                         "MASS_ENV_1" },
    { BINARY_PROPERTY::MASS_ENV_2,                                         "MASS_ENV_2" },
    { BINARY_PROPERTY::MASSES_EQUILIBRATED,                                "MASSES_EQUILIBRATED" },
    { BINARY_PROPERTY::MASSES_EQUILIBRATED_AT_BIRTH,                       "MASSES_EQUILIBRATED_AT_BIRTH" },
    { BINARY_PROPERTY::MASS_TRANSFER_TRACKER_HISTORY,                      "MASS_TRANSFER_TRACKER_HISTORY" },
    { BINARY_PROPERTY::MERGES_IN_HUBBLE_TIME,                              "MERGES_IN_HUBBLE_TIME" },
    { BINARY_PROPERTY::OPTIMISTIC_COMMON_ENVELOPE,                         "OPTIMISTIC_COMMON_ENVELOPE" },
    { BINARY_PROPERTY::ORBITAL_ANGULAR_VELOCITY,                           "ORBITAL_ANGULAR_VELOCITY" },
    { BINARY_PROPERTY::ORBITAL_VELOCITY_PRE_SUPERNOVA,                     "ORBITAL_VELOCITY_PRE_SUPERNOVA" },
    { BINARY_PROPERTY::RADIUS_1_POST_COMMON_ENVELOPE,                      "RADIUS_1_POST_COMMON_ENVELOPE" },
    { BINARY_PROPERTY::RADIUS_1_PRE_COMMON_ENVELOPE,                       "RADIUS_1_PRE_COMMON_ENVELOPE" },
    { BINARY_PROPERTY::RADIUS_2_POST_COMMON_ENVELOPE,                      "RADIUS_2_POST_COMMON_ENVELOPE" },
    { BINARY_PROPERTY::RADIUS_2_PRE_COMMON_ENVELOPE,                       "RADIUS_2_PRE_COMMON_ENVELOPE" },
    { BINARY_PROPERTY::RANDOM_SEED,                                        "RANDOM_SEED" },
    { BINARY_PROPERTY::RLOF_ACCRETION_EFFICIENCY,                          "RLOF_ACCRETION_EFFICIENCY"},
    { BINARY_PROPERTY::RLOF_MASS_LOSS_RATE,                                "RLOF_MASS_LOSS_RATE"},
    { BINARY_PROPERTY::RLOF_MASS_TRANSFER_TIMESCALE,                                "RLOF_MASS_TRANSFER_TIMESCALE"},
    { BINARY_PROPERTY::RLOF_POST_MT_COMMON_ENVELOPE,                       "RLOF_POST_MT_COMMON_ENVELOPE" },
    { BINARY_PROPERTY::RLOF_POST_MT_ECCENTRICITY,                          "RLOF_POST_MT_ECCENTRICITY" },
    { BINARY_PROPERTY::RLOF_POST_MT_EVENT_COUNTER,                         "RLOF_POST_MT_EVENT_COUNTER" },
    { BINARY_PROPERTY::RLOF_POST_MT_ID,                                    "RLOF_POST_MT_ID" },
    { BINARY_PROPERTY::RLOF_POST_MT_SEMI_MAJOR_AXIS,                       "RLOF_POST_MT_SEMI_MAJOR_AXIS" },
    { BINARY_PROPERTY::RLOF_POST_MT_STAR1_MASS,                            "RLOF_POST_MT_STAR1_MASS" },
    { BINARY_PROPERTY::RLOF_POST_MT_STAR2_MASS,                            "RLOF_POST_MT_STAR2_MASS" },
    { BINARY_PROPERTY::RLOF_POST_MT_STAR1_RADIUS,                          "RLOF_POST_MT_STAR1_RADIUS" },
    { BINARY_PROPERTY::RLOF_POST_MT_STAR2_RADIUS,                          "RLOF_POST_MT_STAR2_RADIUS" },
    { BINARY_PROPERTY::RLOF_POST_MT_STAR1_RLOF,                            "RLOF_POST_MT_STAR1_RLOF" },
    { BINARY_PROPERTY::RLOF_POST_MT_STAR2_RLOF,                            "RLOF_POST_MT_STAR2_RLOF" },
    { BINARY_PROPERTY::RLOF_POST_MT_STAR1_STELLAR_TYPE,                    "RLOF_POST_MT_STAR1_STELLAR_TYPE" },
    { BINARY_PROPERTY::RLOF_POST_MT_STAR1_STELLAR_TYPE_NAME,               "RLOF_POST_MT_STAR1_STELLAR_TYPE_NAME" },
    { BINARY_PROPERTY::RLOF_POST_MT_STAR2_STELLAR_TYPE,                    "RLOF_POST_MT_STAR2_STELLAR_TYPE" },
    { BINARY_PROPERTY::RLOF_POST_MT_STAR2_STELLAR_TYPE_NAME,               "RLOF_POST_MT_STAR2_STELLAR_TYPE_NAME" },
    { BINARY_PROPERTY::RLOF_POST_STEP_STAR_TO_ROCHE_LOBE_RADIUS_RATIO_1,   "RLOF_POST_STEP_STAR_TO_ROCHE_LOBE_RADIUS_RATIO_1" },
    { BINARY_PROPERTY::RLOF_POST_STEP_STAR_TO_ROCHE_LOBE_RADIUS_RATIO_2,   "RLOF_POST_STEP_STAR_TO_ROCHE_LOBE_RADIUS_RATIO_2" },
    { BINARY_PROPERTY::RLOF_PRE_MT_ECCENTRICITY,                           "RLOF_PRE_MT_ECCENTRICITY" },
    { BINARY_PROPERTY::RLOF_PRE_MT_SEMI_MAJOR_AXIS,                        "RLOF_PRE_MT_SEMI_MAJOR_AXIS" },
    { BINARY_PROPERTY::RLOF_PRE_MT_STAR1_MASS,                             "RLOF_PRE_MT_STAR1_MASS" },
    { BINARY_PROPERTY::RLOF_PRE_MT_STAR2_MASS,                             "RLOF_PRE_MT_STAR2_MASS" },
    { BINARY_PROPERTY::RLOF_PRE_MT_STAR1_RADIUS,                           "RLOF_PRE_MT_STAR1_RADIUS" },
    { BINARY_PROPERTY::RLOF_PRE_MT_STAR2_RADIUS,                           "RLOF_PRE_MT_STAR2_RADIUS" },
    { BINARY_PROPERTY::RLOF_PRE_MT_STAR1_RLOF,                             "RLOF_PRE_MT_STAR1_RLOF" },
    { BINARY_PROPERTY::RLOF_PRE_MT_STAR2_RLOF,                             "RLOF_PRE_MT_STAR2_RLOF" },
    { BINARY_PROPERTY::RLOF_PRE_MT_STAR1_STELLAR_TYPE,                     "RLOF_PRE_MT_STAR1_STELLAR_TYPE" },
    { BINARY_PROPERTY::RLOF_PRE_MT_STAR1_STELLAR_TYPE_NAME,                "RLOF_PRE_MT_STAR1_STELLAR_TYPE_NAME" },
    { BINARY_PROPERTY::RLOF_PRE_MT_STAR2_STELLAR_TYPE,                     "RLOF_PRE_MT_STAR2_STELLAR_TYPE" },
    { BINARY_PROPERTY::RLOF_PRE_MT_STAR2_STELLAR_TYPE_NAME,                "RLOF_PRE_MT_STAR2_STELLAR_TYPE_NAME" },
    { BINARY_PROPERTY::RLOF_PRE_STEP_STAR_TO_ROCHE_LOBE_RADIUS_RATIO_1,    "RLOF_PRE_STEP_STAR_TO_ROCHE_LOBE_RADIUS_RATIO_1" },
    { BINARY_PROPERTY::RLOF_PRE_STEP_STAR_TO_ROCHE_LOBE_RADIUS_RATIO_2,    "RLOF_PRE_STEP_STAR_TO_ROCHE_LOBE_RADIUS_RATIO_2" },
    { BINARY_PROPERTY::RLOF_SECONDARY_POST_COMMON_ENVELOPE,                "RLOF_SECONDARY_POST_COMMON_ENVELOPE" },
    { BINARY_PROPERTY::RLOF_TIME_POST_MT,                                  "RLOF_TIME_POST_MT" },
    { BINARY_PROPERTY::RLOF_TIME_PRE_MT,                                   "RLOF_TIME_PRE_MT" },
    { BINARY_PROPERTY::ROCHE_LOBE_RADIUS_1,                                "ROCHE_LOBE_RADIUS_1" },
    { BINARY_PROPERTY::ROCHE_LOBE_RADIUS_1_POST_COMMON_ENVELOPE,           "ROCHE_LOBE_RADIUS_1_POST_COMMON_ENVELOPE" },
    { BINARY_PROPERTY::ROCHE_LOBE_RADIUS_1_PRE_COMMON_ENVELOPE,            "ROCHE_LOBE_RADIUS_1_PRE_COMMON_ENVELOPE" },
    { BINARY_PROPERTY::ROCHE_LOBE_RADIUS_2,                                "ROCHE_LOBE_RADIUS_2" },
    { BINARY_PROPERTY::ROCHE_LOBE_RADIUS_2_POST_COMMON_ENVELOPE,           "ROCHE_LOBE_RADIUS_2_POST_COMMON_ENVELOPE" },
    { BINARY_PROPERTY::ROCHE_LOBE_RADIUS_2_PRE_COMMON_ENVELOPE,            "ROCHE_LOBE_RADIUS_2_PRE_COMMON_ENVELOPE" },
    { BINARY_PROPERTY::STAR_TO_ROCHE_LOBE_RADIUS_RATIO_1,                  "STAR_TO_ROCHE_LOBE_RADIUS_RATIO_1" },
    { BINARY_PROPERTY::STAR_TO_ROCHE_LOBE_RADIUS_RATIO_2,                  "STAR_TO_ROCHE_LOBE_RADIUS_RATIO_2" },
    { BINARY_PROPERTY::SEMI_MAJOR_AXIS_AT_DCO_FORMATION,                   "SEMI_MAJOR_AXIS_AT_DCO_FORMATION" },
    { BINARY_PROPERTY::SEMI_MAJOR_AXIS_INITIAL,                            "SEMI_MAJOR_AXIS_INITIAL" },
    { BINARY_PROPERTY::SEMI_MAJOR_AXIS_POST_COMMON_ENVELOPE,               "SEMI_MAJOR_AXIS_POST_COMMON_ENVELOPE" },
    { BINARY_PROPERTY::SEMI_MAJOR_AXIS_PRE_SUPERNOVA,                      "SEMI_MAJOR_AXIS_PRE_SUPERNOVA" },
    { BINARY_PROPERTY::SEMI_MAJOR_AXIS_PRE_SUPERNOVA_RSOL,                 "SEMI_MAJOR_AXIS_PRE_SUPERNOVA_RSOL" },
    { BINARY_PROPERTY::SEMI_MAJOR_AXIS_PRE_COMMON_ENVELOPE,                "SEMI_MAJOR_AXIS_PRE_COMMON_ENVELOPE" },
    { BINARY_PROPERTY::SEMI_MAJOR_AXIS,                                    "SEMI_MAJOR_AXIS" },
    { BINARY_PROPERTY::SEMI_MAJOR_AXIS_RSOL,                               "SEMI_MAJOR_AXIS_RSOL" },
    { BINARY_PROPERTY::SIMULTANEOUS_RLOF,                                  "SIMULTANEOUS_RLOF" },
    { BINARY_PROPERTY::STABLE_RLOF_POST_COMMON_ENVELOPE,                   "STABLE_RLOF_POST_COMMON_ENVELOPE" },
    { BINARY_PROPERTY::STELLAR_MERGER,                                     "STELLAR_MERGER" },
    { BINARY_PROPERTY::STELLAR_MERGER_AT_BIRTH,                            "STELLAR_MERGER_AT_BIRTH" },
    { BINARY_PROPERTY::STELLAR_TYPE_1_POST_COMMON_ENVELOPE,                "STELLAR_TYPE_1_POST_COMMON_ENVELOPE" },
    { BINARY_PROPERTY::STELLAR_TYPE_1_PRE_COMMON_ENVELOPE,                 "STELLAR_TYPE_1_PRE_COMMON_ENVELOPE" },
    { BINARY_PROPERTY::STELLAR_TYPE_2_POST_COMMON_ENVELOPE,                "STELLAR_TYPE_2_POST_COMMON_ENVELOPE" },
    { BINARY_PROPERTY::STELLAR_TYPE_2_PRE_COMMON_ENVELOPE,                 "STELLAR_TYPE_2_PRE_COMMON_ENVELOPE" },
    { BINARY_PROPERTY::STELLAR_TYPE_NAME_1_POST_COMMON_ENVELOPE,           "STELLAR_TYPE_NAME_1_POST_COMMON_ENVELOPE" },
    { BINARY_PROPERTY::STELLAR_TYPE_NAME_1_PRE_COMMON_ENVELOPE,            "STELLAR_TYPE_NAME_1_PRE_COMMON_ENVELOPE" },
    { BINARY_PROPERTY::STELLAR_TYPE_NAME_2_POST_COMMON_ENVELOPE,           "STELLAR_TYPE_NAME_2_POST_COMMON_ENVELOPE" },
    { BINARY_PROPERTY::STELLAR_TYPE_NAME_2_PRE_COMMON_ENVELOPE,            "STELLAR_TYPE_NAME_2_PRE_COMMON_ENVELOPE" },
    { BINARY_PROPERTY::SUPERNOVA_ORBIT_INCLINATION_ANGLE,                  "SUPERNOVA_ORBIT_INCLINATION_ANGLE" },
    { BINARY_PROPERTY::SUPERNOVA_ORBIT_INCLINATION_VECTOR_X,               "SUPERNOVA_ORBIT_INCLINATION_VECTOR_X" },
    { BINARY_PROPERTY::SUPERNOVA_ORBIT_INCLINATION_VECTOR_Y,               "SUPERNOVA_ORBIT_INCLINATION_VECTOR_Y" },
    { BINARY_PROPERTY::SUPERNOVA_ORBIT_INCLINATION_VECTOR_Z,               "SUPERNOVA_ORBIT_INCLINATION_VECTOR_Z" },
    { BINARY_PROPERTY::SUPERNOVA_STATE,                                    "SUPERNOVA_STATE" },
    { BINARY_PROPERTY::SYNCHRONIZATION_TIMESCALE,                          "SYNCHRONIZATION_TIMESCALE" },
    { BINARY_PROPERTY::SYSTEMIC_SPEED,                                     "SYSTEMIC_SPEED" },
    { BINARY_PROPERTY::SYSTEMIC_VELOCITY_X,                                "SYSTEMIC_VELOCITY_X" },
    { BINARY_PROPERTY::SYSTEMIC_VELOCITY_Y,                                "SYSTEMIC_VELOCITY_Y" },
    { BINARY_PROPERTY::SYSTEMIC_VELOCITY_Z,                                "SYSTEMIC_VELOCITY_Z" },
    { BINARY_PROPERTY::TIME,                                               "TIME" },
    { BINARY_PROPERTY::TIME_TO_COALESCENCE,                                "TIME_TO_COALESCENCE" },
    { BINARY_PROPERTY::TOTAL_ANGULAR_MOMENTUM,                             "TOTAL_ANGULAR_MOMENTUM" },
    { BINARY_PROPERTY::TOTAL_ENERGY,                                       "TOTAL_ENERGY" },
    { BINARY_PROPERTY::ZETA_LOBE,                                          "ZETA_LOBE" },
    { BINARY_PROPERTY::ZETA_STAR,                                          "ZETA_STAR" }
};


// enum class PROGRAM_OPTION
// Symbolic names for program option values
//
// Options only need to be here if they are required to be
// available for printing in the logfiles
enum class PROGRAM_OPTION: int {

    NONE,

    ADD_OPTIONS_TO_SYSPARMS,
    ALLOW_NON_STRIPPED_ECSN,
    ALLOW_IMMEDIATE_RLOF_POST_CE_TO_SURVIVE_COMMON_ENVELOPE,
    ALLOW_MS_STAR_TO_SURVIVE_COMMON_ENVELOPE,
    ALLOW_RADIATIVE_ENVELOPE_STAR_TO_SURVIVE_COMMON_ENVELOPE,
    ALLOW_RLOF_AT_BIRTH,
    ALLOW_TOUCHING_AT_BIRTH,
    ANG_MOM_CONSERVATION_DURING_CIRCULARISATION,

    BLACK_HOLE_KICKS_MODE,
    
    CASE_BB_STABILITY_PRESCRIPTION,
    
    CHECK_PHOTON_TIRING_LIMIT,

    CHE_MODE,

    CIRCULARISE_BINARY_DURING_MT,

    COMMON_ENVELOPE_ALPHA,
    COMMON_ENVELOPE_ALPHA_THERMAL,
    COMMON_ENVELOPE_FORMALISM,
    COMMON_ENVELOPE_LAMBDA,
    COMMON_ENVELOPE_LAMBDA_MULTIPLIER,
    COMMON_ENVELOPE_LAMBDA_PRESCRIPTION,
    COMMON_ENVELOPE_MASS_ACCRETION_CONSTANT,
    COMMON_ENVELOPE_MASS_ACCRETION_MAX,
    COMMON_ENVELOPE_MASS_ACCRETION_MIN,
    COMMON_ENVELOPE_MASS_ACCRETION_PRESCRIPTION,
    COMMON_ENVELOPE_RECOMBINATION_ENERGY_DENSITY,
    COMMON_ENVELOPE_SLOPE_KRUCKOW,

    CONVECTIVE_ENVELOPE_TEMPERATURE_THRESHOLD,

    COOL_WIND_MASS_LOSS_MULTIPLIER,

    ECCENTRICITY,
    ECCENTRICITY_DISTRIBUTION,
    ECCENTRICITY_DISTRIBUTION_MAX,
    ECCENTRICITY_DISTRIBUTION_MIN,
    EDDINGTON_ACCRETION_FACTOR,
    ENABLE_ROTATIONALLY_ENHANCED_MASS_LOSS,
    ENHANCE_CHE_LIFETIMES_LUMINOSITIES,
    ENVELOPE_STATE_PRESCRIPTION,
    EVOLUTION_MODE,

    FP_ERROR_MODE,

    FRYER_SUPERNOVA_ENGINE,

    FRYER22_FMIX,
    FRYER22_MCRIT,

    INITIAL_MASS,
    INITIAL_MASS_1,
    INITIAL_MASS_2,

    INITIAL_MASS_FUNCTION,
    INITIAL_MASS_FUNCTION_MAX,
    INITIAL_MASS_FUNCTION_MIN,
    INITIAL_MASS_FUNCTIONPOWER,

    KICK_DIRECTION_DISTRIBUTION,
    KICK_DIRECTION_POWER,
    KICK_SCALING_FACTOR,
    KICK_MAGNITUDE_DISTRIBUTION,
    KICK_MAGNITUDE_DISTRIBUTION_MAXIMUM,

    KICK_MAGNITUDE_DISTRIBUTION_SIGMA_CCSN_BH,
    KICK_MAGNITUDE_DISTRIBUTION_SIGMA_CCSN_NS,
    KICK_MAGNITUDE_DISTRIBUTION_SIGMA_FOR_ECSN,
    KICK_MAGNITUDE_DISTRIBUTION_SIGMA_FOR_USSN,

    KICK_MAGNITUDE,
    KICK_MAGNITUDE_1,
    KICK_MAGNITUDE_2,

    KICK_MAGNITUDE_RANDOM,
    KICK_MAGNITUDE_RANDOM_1,
    KICK_MAGNITUDE_RANDOM_2,

    KICK_MEAN_ANOMALY_1,
    KICK_MEAN_ANOMALY_2,
    KICK_PHI_1,
    KICK_PHI_2,
    KICK_THETA_1,
    KICK_THETA_2,

    LBV_FACTOR,
    LBV_MASS_LOSS_PRESCRIPTION,
    MASS_LOSS_PRESCRIPTION,

    MASS_RATIO,
    MASS_RATIO_DISTRIBUTION,
    MASS_RATIO_DISTRIBUTION_MAX,
    MASS_RATIO_DISTRIBUTION_MIN,

    MAXIMUM_EVOLUTION_TIME,
    MAXIMUM_DONOR_MASS,
    MAXIMUM_NEUTRON_STAR_MASS,
    MAXIMUM_TIMESTEPS,

    MCBUR1,

    METALLICITY,
    METALLICITY_DISTRIBUTION,
    METALLICITY_DISTRIBUTION_MAX,
    METALLICITY_DISTRIBUTION_MIN,

    MINIMUM_MASS_SECONDARY,

    MT_ACCRETION_EFFICIENCY_PRESCRIPTION,
    MT_ANG_MOM_LOSS_PRESCRIPTION,
    MT_THERMAL_LIMIT_C,

    MT_CRIT_MR_MS_LOW_MASS_DEGENERATE_ACCRETOR,
    MT_CRIT_MR_MS_LOW_MASS_NON_DEGENERATE_ACCRETOR,
    MT_CRIT_MR_MS_HIGH_MASS_DEGENERATE_ACCRETOR,
    MT_CRIT_MR_MS_HIGH_MASS_NON_DEGENERATE_ACCRETOR,
    MT_CRIT_MR_GIANT_DEGENERATE_ACCRETOR,
    MT_CRIT_MR_GIANT_NON_DEGENERATE_ACCRETOR,
    MT_CRIT_MR_HG_DEGENERATE_ACCRETOR,
    MT_CRIT_MR_HG_NON_DEGENERATE_ACCRETOR,
    MT_CRIT_MR_HE_GIANT_DEGENERATE_ACCRETOR,
    MT_CRIT_MR_HE_GIANT_NON_DEGENERATE_ACCRETOR,
    MT_CRIT_MR_HE_HG_DEGENERATE_ACCRETOR,
    MT_CRIT_MR_HE_HG_NON_DEGENERATE_ACCRETOR,
    MT_CRIT_MR_HE_MS_DEGENERATE_ACCRETOR,
    MT_CRIT_MR_HE_MS_NON_DEGENERATE_ACCRETOR,
    MT_CRIT_MR_WD_DEGENERATE_ACCRETOR,
    MT_CRIT_MR_WD_NONDEGENERATE_ACCRETOR,

    MT_FRACTION_ACCRETED,
    MT_JLOSS,
    MT_JLOSS_MACLEOD_LINEAR_FRACTION_DEGEN,
    MT_JLOSS_MACLEOD_LINEAR_FRACTION_NON_DEGEN,
    MT_REJUVENATION_PRESCRIPTION,
    MT_THERMALLY_LIMITED_VARIATION,

    MULLER_MANDEL_KICK_MULTIPLIER_BH,
    MULLER_MANDEL_KICK_MULTIPLIER_NS,
    MULLER_MANDEL_SIGMA_KICK,

    NEUTRINO_MASS_LOSS_ASSUMPTION_BH,
    NEUTRINO_MASS_LOSS_VALUE_BH,

    NOTES,

    NS_EOS,

    ORBITAL_PERIOD,
    ORBITAL_PERIOD_DISTRIBUTION,
    ORBITAL_PERIOD_DISTRIBUTION_MAX,
    ORBITAL_PERIOD_DISTRIBUTION_MIN,

    OVERALL_WIND_MASS_LOSS_MULTIPLIER,

    PISN_LOWER_LIMIT,
    PISN_UPPER_LIMIT,

    PPI_LOWER_LIMIT,
    PPI_PRESCRIPTION,
    PPI_UPPER_LIMIT,
    PPI_CO_CORE_SHIFT_HENDRIKS,

    PULSAR_MAGNETIC_FIELD_DISTRIBUTION,
    PULSAR_MAGNETIC_FIELD_DISTRIBUTION_MAX,
    PULSAR_MAGNETIC_FIELD_DISTRIBUTION_MIN,

    PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION,
    PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION_MAX,
    PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION_MIN,

    PULSAR_MAGNETIC_FIELD_DECAY_MASS_SCALE,
    PULSAR_MAGNETIC_FIELD_DECAY_TIME_SCALE,

    PULSAR_MINIMUM_MAGNETIC_FIELD,

    QCRIT_PRESCRIPTION,

    RANDOM_SEED,
    RANDOM_SEED_CMDLINE,

    REMNANT_MASS_PRESCRIPTION,

    ROCKET_KICK_MAGNITUDE_1,
    ROCKET_KICK_MAGNITUDE_2,
    ROCKET_KICK_PHI_1,
    ROCKET_KICK_PHI_2,
    ROCKET_KICK_THETA_1,
    ROCKET_KICK_THETA_2,

    ROTATIONAL_VELOCITY_DISTRIBUTION,
    ROTATIONAL_FREQUENCY,
    ROTATIONAL_FREQUENCY_1,
    ROTATIONAL_FREQUENCY_2,
    
    SCALE_CHE_MASS_LOSS_SURF_HE_ABUNDANCE,
    SCALE_TERMINAL_WIND_VEL_METALLICITY_POWER,
    SEMI_MAJOR_AXIS,
    SEMI_MAJOR_AXIS_DISTRIBUTION,
    SEMI_MAJOR_AXIS_DISTRIBUTION_MAX,
    SEMI_MAJOR_AXIS_DISTRIBUTION_MIN,

    STELLAR_ZETA_PRESCRIPTION,

    TIDES_PRESCRIPTION,

    WR_FACTOR,

    ZETA_ADIABATIC_ARBITRARY,
    ZETA_MS,
    ZETA_RADIATIVE_ENVELOPE_GIANT
};
// map PROGRAM_OPTION to string identifying the property
// for lookup by the printing functions
//
// Options only need to be here if they are required to be
// available for printing in the logfiles - all keys present here
// should also be in PROGRAM_OPTION_DETAIL ( except PROGRAM_OPTION::NOTES)
const COMPASUnorderedMap<PROGRAM_OPTION, std::string> PROGRAM_OPTION_LABEL = {

    { PROGRAM_OPTION::NONE,                                             "NONE" },

    { PROGRAM_OPTION::ALLOW_NON_STRIPPED_ECSN,                          "ALLOW_NON_STRIPPED_ECSN" },
    { PROGRAM_OPTION::ADD_OPTIONS_TO_SYSPARMS,                          "ADD_OPTIONS_TO_SYSPARMS" },
    { PROGRAM_OPTION::ALLOW_MS_STAR_TO_SURVIVE_COMMON_ENVELOPE,         "ALLOW_MS_STAR_TO_SURVIVE_COMMON_ENVELOPE" },
    { PROGRAM_OPTION::ALLOW_RADIATIVE_ENVELOPE_STAR_TO_SURVIVE_COMMON_ENVELOPE, "ALLOW_RADIATIVE_ENVELOPE_STAR_TO_SURVIVE_COMMON_ENVELOPE" },
    { PROGRAM_OPTION::ALLOW_IMMEDIATE_RLOF_POST_CE_TO_SURVIVE_COMMON_ENVELOPE,  "ALLOW_IMMEDIATE_RLOF_POST_CE_TO_SURVIVE_COMMON_ENVELOPE" },
    { PROGRAM_OPTION::ALLOW_RLOF_AT_BIRTH,                              "ALLOW_RLOF_AT_BIRTH" },
    { PROGRAM_OPTION::ALLOW_TOUCHING_AT_BIRTH,                          "ALLOW_TOUCHING_AT_BIRTH" },
    { PROGRAM_OPTION::ANG_MOM_CONSERVATION_DURING_CIRCULARISATION,      "ANG_MOM_CONSERVATION_DURING_CIRCULARISATION" },

    { PROGRAM_OPTION::BLACK_HOLE_KICKS_MODE,                            "BLACK_HOLE_KICKS_MODE" },
    
    { PROGRAM_OPTION::CASE_BB_STABILITY_PRESCRIPTION,                   "CASE_BB_STABILITY_PRESCRIPTION" },
    
    { PROGRAM_OPTION::CHECK_PHOTON_TIRING_LIMIT,                        "CHECK_PHOTON_TIRING_LIMIT" },

    { PROGRAM_OPTION::CHE_MODE,                                         "CHE_MODE" },

    { PROGRAM_OPTION::CIRCULARISE_BINARY_DURING_MT,                     "CIRCULARISE_BINARY_DURING_MT" },

    { PROGRAM_OPTION::COMMON_ENVELOPE_ALPHA,                            "COMMON_ENVELOPE_ALPHA" },
    { PROGRAM_OPTION::COMMON_ENVELOPE_ALPHA_THERMAL,                    "COMMON_ENVELOPE_ALPHA_THERMAL" },
    { PROGRAM_OPTION::COMMON_ENVELOPE_FORMALISM,                        "COMMON_ENVELOPE_FORMALISM" },
    { PROGRAM_OPTION::COMMON_ENVELOPE_LAMBDA,                           "COMMON_ENVELOPE_LAMBDA" },
    { PROGRAM_OPTION::COMMON_ENVELOPE_LAMBDA_MULTIPLIER,                "COMMON_ENVELOPE_LAMBDA_MULTIPLIER" },
    { PROGRAM_OPTION::COMMON_ENVELOPE_LAMBDA_PRESCRIPTION,              "COMMON_ENVELOPE_LAMBDA_PRESCRIPTION" },
    { PROGRAM_OPTION::COMMON_ENVELOPE_MASS_ACCRETION_CONSTANT,          "COMMON_ENVELOPE_MASS_ACCRETION_CONSTANT" },
    { PROGRAM_OPTION::COMMON_ENVELOPE_MASS_ACCRETION_MAX,               "COMMON_ENVELOPE_MASS_ACCRETION_MAX" },
    { PROGRAM_OPTION::COMMON_ENVELOPE_MASS_ACCRETION_MIN,               "COMMON_ENVELOPE_MASS_ACCRETION_MIN" },
    { PROGRAM_OPTION::COMMON_ENVELOPE_MASS_ACCRETION_PRESCRIPTION,      "COMMON_ENVELOPE_MASS_ACCRETION_PRESCRIPTION" },
    { PROGRAM_OPTION::COMMON_ENVELOPE_RECOMBINATION_ENERGY_DENSITY,     "COMMON_ENVELOPE_RECOMBINATION_ENERGY_DENSITY" },
    { PROGRAM_OPTION::COMMON_ENVELOPE_SLOPE_KRUCKOW,                    "COMMON_ENVELOPE_SLOPE_KRUCKOW" },

    { PROGRAM_OPTION::COOL_WIND_MASS_LOSS_MULTIPLIER,                   "COOL_WIND_MASS_LOSS_MULTIPLIER" },

    { PROGRAM_OPTION::ECCENTRICITY,                                     "ECCENTRICITY" },
    { PROGRAM_OPTION::ECCENTRICITY_DISTRIBUTION,                        "ECCENTRICITY_DISTRIBUTION" },
    { PROGRAM_OPTION::ECCENTRICITY_DISTRIBUTION_MAX,                    "ECCENTRICITY_DISTRIBUTION_MAX" },
    { PROGRAM_OPTION::ECCENTRICITY_DISTRIBUTION_MIN,                    "ECCENTRICITY_DISTRIBUTION_MIN" },
    { PROGRAM_OPTION::EDDINGTON_ACCRETION_FACTOR,                       "EDDINGTON_ACCRETION_FACTOR" },
    { PROGRAM_OPTION::ENABLE_ROTATIONALLY_ENHANCED_MASS_LOSS,           "ENABLE_ROTATIONALLY_ENHANCED_MASS_LOSS" },
    { PROGRAM_OPTION::ENHANCE_CHE_LIFETIMES_LUMINOSITIES,               "ENHANCE_CHE_LIFETIMES_LUMINOSITIES" }, 
    { PROGRAM_OPTION::ENVELOPE_STATE_PRESCRIPTION,                      "ENVELOPE_STATE_PRESCRIPTION" },
    { PROGRAM_OPTION::EVOLUTION_MODE,                                   "EVOLUTION_MODE" },

    { PROGRAM_OPTION::FP_ERROR_MODE,                                    "FP_ERROR_MODE" },

    { PROGRAM_OPTION::FRYER_SUPERNOVA_ENGINE,                           "FRYER_SUPERNOVA_ENGINE" },

    { PROGRAM_OPTION::FRYER22_FMIX,                                     "FRYER22_FMIX" },
    { PROGRAM_OPTION::FRYER22_MCRIT,                                    "FRYER22_MCRIT" },

    { PROGRAM_OPTION::INITIAL_MASS,                                     "INITIAL_MASS" },
    { PROGRAM_OPTION::INITIAL_MASS_1,                                   "INITIAL_MASS_1" },
    { PROGRAM_OPTION::INITIAL_MASS_2,                                   "INITIAL_MASS_2" },

    { PROGRAM_OPTION::INITIAL_MASS_FUNCTION,                            "INITIAL_MASS_FUNCTION" },
    { PROGRAM_OPTION::INITIAL_MASS_FUNCTION_MAX,                        "INITIAL_MASS_FUNCTION_MAX" },
    { PROGRAM_OPTION::INITIAL_MASS_FUNCTION_MIN,                        "INITIAL_MASS_FUNCTION_MIN" },
    { PROGRAM_OPTION::INITIAL_MASS_FUNCTIONPOWER,                       "INITIAL_MASS_FUNCTIONPOWER" },
    { PROGRAM_OPTION::KICK_DIRECTION_DISTRIBUTION,                      "KICK_DIRECTION_DISTRIBUTION" },
    { PROGRAM_OPTION::KICK_DIRECTION_POWER,                             "KICK_DIRECTION_POWER" },
    { PROGRAM_OPTION::KICK_SCALING_FACTOR,                              "KICK_SCALING_FACTOR" },
    { PROGRAM_OPTION::KICK_MAGNITUDE_DISTRIBUTION,                      "KICK_MAGNITUDE_DISTRIBUTION" },
    { PROGRAM_OPTION::KICK_MAGNITUDE_DISTRIBUTION_MAXIMUM,              "KICK_MAGNITUDE_DISTRIBUTION_MAXIMUM" },

    { PROGRAM_OPTION::KICK_MAGNITUDE_DISTRIBUTION_SIGMA_CCSN_BH,        "KICK_MAGNITUDE_DISTRIBUTION_SIGMA_CCSN_BH" },
    { PROGRAM_OPTION::KICK_MAGNITUDE_DISTRIBUTION_SIGMA_CCSN_NS,        "KICK_MAGNITUDE_DISTRIBUTION_SIGMA_CCSN_NS" },
    { PROGRAM_OPTION::KICK_MAGNITUDE_DISTRIBUTION_SIGMA_FOR_ECSN,       "KICK_MAGNITUDE_DISTRIBUTION_SIGMA_FOR_ECSN" },
    { PROGRAM_OPTION::KICK_MAGNITUDE_DISTRIBUTION_SIGMA_FOR_USSN,       "KICK_MAGNITUDE_DISTRIBUTION_SIGMA_FOR_USSN" },

    { PROGRAM_OPTION::KICK_MAGNITUDE,                                   "KICK_MAGNITUDE" },
    { PROGRAM_OPTION::KICK_MAGNITUDE_1,                                 "KICK_MAGNITUDE_1" },
    { PROGRAM_OPTION::KICK_MAGNITUDE_2,                                 "KICK_MAGNITUDE_2" },

    { PROGRAM_OPTION::KICK_MAGNITUDE_RANDOM,                            "KICK_MAGNITUDE_RANDOM" },
    { PROGRAM_OPTION::KICK_MAGNITUDE_RANDOM_1,                          "KICK_MAGNITUDE_RANDOM_1" },
    { PROGRAM_OPTION::KICK_MAGNITUDE_RANDOM_2,                          "KICK_MAGNITUDE_RANDOM_2" },

    { PROGRAM_OPTION::KICK_MEAN_ANOMALY_1,                              "KICK_MEAN_ANOMALY_1" },
    { PROGRAM_OPTION::KICK_MEAN_ANOMALY_2,                              "KICK_MEAN_ANOMALY_2" },
    { PROGRAM_OPTION::KICK_PHI_1,                                       "KICK_PHI_1" },
    { PROGRAM_OPTION::KICK_PHI_2,                                       "KICK_PHI_2" },
    { PROGRAM_OPTION::KICK_THETA_1,                                     "KICK_THETA_1" },
    { PROGRAM_OPTION::KICK_THETA_2,                                     "KICK_THETA_2" },

    { PROGRAM_OPTION::LBV_FACTOR,                                       "LBV_FACTOR" },
    { PROGRAM_OPTION::LBV_MASS_LOSS_PRESCRIPTION,                       "LBV_MASS_LOSS_PRESCRIPTION" },
    { PROGRAM_OPTION::MASS_LOSS_PRESCRIPTION,                           "MASS_LOSS_PRESCRIPTION" },

    { PROGRAM_OPTION::MASS_RATIO,                                       "MASS_RATIO" },
    { PROGRAM_OPTION::MASS_RATIO_DISTRIBUTION,                          "MASS_RATIO_DISTRIBUTION" },
    { PROGRAM_OPTION::MASS_RATIO_DISTRIBUTION_MAX,                      "MASS_RATIO_DISTRIBUTION_MAX" },
    { PROGRAM_OPTION::MASS_RATIO_DISTRIBUTION_MIN,                      "MASS_RATIO_DISTRIBUTION_MIN" },

    { PROGRAM_OPTION::MAXIMUM_EVOLUTION_TIME,                           "MAXIMUM_EVOLUTION_TIME" },
    { PROGRAM_OPTION::MAXIMUM_DONOR_MASS,                               "MAXIMUM_DONOR_MASS" },
    { PROGRAM_OPTION::MAXIMUM_NEUTRON_STAR_MASS,                        "MAXIMUM_NEUTRON_STAR_MASS" },
    { PROGRAM_OPTION::MAXIMUM_TIMESTEPS,                                "MAXIMUM_TIMESTEPS" },

    { PROGRAM_OPTION::MCBUR1,                                           "MCBUR1" },

    { PROGRAM_OPTION::METALLICITY,                                      "METALLICITY" },
    { PROGRAM_OPTION::METALLICITY_DISTRIBUTION,                         "METALLICITY_DISTRIBUTION" },
    { PROGRAM_OPTION::METALLICITY_DISTRIBUTION_MAX,                     "METALLICITY_DISTRIBUTION_MAX" },
    { PROGRAM_OPTION::METALLICITY_DISTRIBUTION_MIN,                     "METALLICITY_DISTRIBUTION_MIN" },

    { PROGRAM_OPTION::MINIMUM_MASS_SECONDARY,                           "MINIMUM_MASS_SECONDARY" },

    { PROGRAM_OPTION::MT_ACCRETION_EFFICIENCY_PRESCRIPTION,             "MT_ACCRETION_EFFICIENCY_PRESCRIPTION" },
    { PROGRAM_OPTION::MT_ANG_MOM_LOSS_PRESCRIPTION,                     "MT_ANG_MOM_LOSS_PRESCRIPTION" },
    { PROGRAM_OPTION::MT_THERMAL_LIMIT_C,                               "MT_THERMAL_LIMIT_C" },

    { PROGRAM_OPTION::MT_CRIT_MR_MS_LOW_MASS_DEGENERATE_ACCRETOR,       "MT_CRIT_MR_MS_LOW_MASS_DEGENERATE_ACCRETOR" },
    { PROGRAM_OPTION::MT_CRIT_MR_MS_LOW_MASS_NON_DEGENERATE_ACCRETOR,   "MT_CRIT_MR_MS_LOW_MASS_NON_DEGENERATE_ACCRETOR" },
    { PROGRAM_OPTION::MT_CRIT_MR_MS_HIGH_MASS_DEGENERATE_ACCRETOR,      "MT_CRIT_MR_MS_HIGH_MASS_DEGENERATE_ACCRETOR" },
    { PROGRAM_OPTION::MT_CRIT_MR_MS_HIGH_MASS_NON_DEGENERATE_ACCRETOR,  "MT_CRIT_MR_MS_HIGH_MASS_NON_DEGENERATE_ACCRETOR" },
    { PROGRAM_OPTION::MT_CRIT_MR_GIANT_DEGENERATE_ACCRETOR,             "MT_CRIT_MR_GIANT_DEGENERATE_ACCRETOR" },
    { PROGRAM_OPTION::MT_CRIT_MR_GIANT_NON_DEGENERATE_ACCRETOR,         "MT_CRIT_MR_GIANT_NON_DEGENERATE_ACCRETOR" },
    { PROGRAM_OPTION::MT_CRIT_MR_HG_DEGENERATE_ACCRETOR,                "MT_CRIT_MR_HG_DEGENERATE_ACCRETOR" },
    { PROGRAM_OPTION::MT_CRIT_MR_HG_NON_DEGENERATE_ACCRETOR,            "MT_CRIT_MR_HG_NON_DEGENERATE_ACCRETOR" },
    { PROGRAM_OPTION::MT_CRIT_MR_HE_GIANT_DEGENERATE_ACCRETOR,          "MT_CRIT_MR_HE_GIANT_DEGENERATE_ACCRETOR" },
    { PROGRAM_OPTION::MT_CRIT_MR_HE_GIANT_NON_DEGENERATE_ACCRETOR,      "MT_CRIT_MR_HE_GIANT_NON_DEGENERATE_ACCRETOR" },
    { PROGRAM_OPTION::MT_CRIT_MR_HE_HG_DEGENERATE_ACCRETOR,             "MT_CRIT_MR_HE_HG_DEGENERATE_ACCRETOR" },
    { PROGRAM_OPTION::MT_CRIT_MR_HE_HG_NON_DEGENERATE_ACCRETOR,         "MT_CRIT_MR_HE_HG_NON_DEGENERATE_ACCRETOR" },
    { PROGRAM_OPTION::MT_CRIT_MR_HE_MS_DEGENERATE_ACCRETOR,             "MT_CRIT_MR_HE_MS_DEGENERATE_ACCRETOR" },
    { PROGRAM_OPTION::MT_CRIT_MR_HE_MS_NON_DEGENERATE_ACCRETOR,         "MT_CRIT_MR_HE_MS_NON_DEGENERATE_ACCRETOR" },
    { PROGRAM_OPTION::MT_CRIT_MR_WD_DEGENERATE_ACCRETOR,                "MT_CRIT_MR_WD_DEGENERATE_ACCRETOR" },
    { PROGRAM_OPTION::MT_CRIT_MR_WD_NONDEGENERATE_ACCRETOR,             "MT_CRIT_MR_WD_NONDEGENERATE_ACCRETOR" },

    { PROGRAM_OPTION::MT_FRACTION_ACCRETED,                             "MT_FRACTION_ACCRETED" },
    { PROGRAM_OPTION::MT_JLOSS,                                         "MT_JLOSS" },
    { PROGRAM_OPTION::MT_JLOSS_MACLEOD_LINEAR_FRACTION_DEGEN,           "MT_JLOSS_MACLEOD_LINEAR_FRACTION_DEGEN" },
    { PROGRAM_OPTION::MT_JLOSS_MACLEOD_LINEAR_FRACTION_NON_DEGEN,       "MT_JLOSS_MACLEOD_LINEAR_FRACTION_NON_DEGEN" },
    { PROGRAM_OPTION::MT_REJUVENATION_PRESCRIPTION,                     "MT_REJUVENATION_PRESCRIPTION" },
    { PROGRAM_OPTION::MT_THERMALLY_LIMITED_VARIATION,                   "MT_THERMALLY_LIMITED_VARIATION" },

    { PROGRAM_OPTION::MULLER_MANDEL_KICK_MULTIPLIER_BH,                 "MULLER_MANDEL_KICK_MULTIPLIER_BH" },
    { PROGRAM_OPTION::MULLER_MANDEL_KICK_MULTIPLIER_NS,                 "MULLER_MANDEL_KICK_MULTIPLIER_NS" },
    { PROGRAM_OPTION::MULLER_MANDEL_SIGMA_KICK,                         "MULLER_MANDEL_SIGMA_KICK" },

    { PROGRAM_OPTION::NEUTRINO_MASS_LOSS_ASSUMPTION_BH,                 "NEUTRINO_MASS_LOSS_ASSUMPTION_BH" },
    { PROGRAM_OPTION::NEUTRINO_MASS_LOSS_VALUE_BH,                      "NEUTRINO_MASS_LOSS_VALUE_BH" },

    { PROGRAM_OPTION::NOTES,                                            "NOTES" },

    { PROGRAM_OPTION::NS_EOS,                                           "NS_EOS" },

    { PROGRAM_OPTION::ORBITAL_PERIOD,                                   "ORBITAL_PERIOD" },
    { PROGRAM_OPTION::ORBITAL_PERIOD_DISTRIBUTION,                      "ORBITAL_PERIOD_DISTRIBUTION" },
    { PROGRAM_OPTION::ORBITAL_PERIOD_DISTRIBUTION_MAX,                  "ORBITAL_PERIOD_DISTRIBUTION_MAX" },
    { PROGRAM_OPTION::ORBITAL_PERIOD_DISTRIBUTION_MIN,                  "ORBITAL_PERIOD_DISTRIBUTION_MIN" },

    { PROGRAM_OPTION::OVERALL_WIND_MASS_LOSS_MULTIPLIER,                "OVERALL_WIND_MASS_LOSS_MULTIPLIER" },

    { PROGRAM_OPTION::PISN_LOWER_LIMIT,                                 "PISN_LOWER_LIMIT" },
    { PROGRAM_OPTION::PISN_UPPER_LIMIT,                                 "PISN_UPPER_LIMIT" },

    { PROGRAM_OPTION::PPI_LOWER_LIMIT,                                  "PPI_LOWER_LIMIT" },
    { PROGRAM_OPTION::PPI_PRESCRIPTION,                                 "PPI_PRESCRIPTION" },
    { PROGRAM_OPTION::PPI_UPPER_LIMIT,                                  "PPI_UPPER_LIMIT" },
    { PROGRAM_OPTION::PPI_CO_CORE_SHIFT_HENDRIKS,                       "PPI_CO_CORE_SHIFT_HENDRIKS"},

    { PROGRAM_OPTION::PULSAR_MAGNETIC_FIELD_DISTRIBUTION,               "PULSAR_MAGNETIC_FIELD_DISTRIBUTION" },
    { PROGRAM_OPTION::PULSAR_MAGNETIC_FIELD_DISTRIBUTION_MAX,           "PULSAR_MAGNETIC_FIELD_DISTRIBUTION_MAX" },
    { PROGRAM_OPTION::PULSAR_MAGNETIC_FIELD_DISTRIBUTION_MIN,           "PULSAR_MAGNETIC_FIELD_DISTRIBUTION_MIN" },

    { PROGRAM_OPTION::PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION,            "PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION" },
    { PROGRAM_OPTION::PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION_MAX,        "PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION_MAX" },
    { PROGRAM_OPTION::PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION_MIN,        "PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION_MIN" },

    { PROGRAM_OPTION::PULSAR_MAGNETIC_FIELD_DECAY_MASS_SCALE,           "PULSAR_MAGNETIC_FIELD_DECAY_MASS_SCALE" },
    { PROGRAM_OPTION::PULSAR_MAGNETIC_FIELD_DECAY_TIME_SCALE,           "PULSAR_MAGNETIC_FIELD_DECAY_TIME_SCALE" },

    { PROGRAM_OPTION::PULSAR_MINIMUM_MAGNETIC_FIELD,                    "PULSAR_MINIMUM_MAGNETIC_FIELD" },

    { PROGRAM_OPTION::QCRIT_PRESCRIPTION,                               "QCRIT_PRESCRIPTION" },

    { PROGRAM_OPTION::RANDOM_SEED,                                      "RANDOM_SEED" },
    { PROGRAM_OPTION::RANDOM_SEED_CMDLINE,                              "RANDOM_SEED_CMDLINE" },

    { PROGRAM_OPTION::REMNANT_MASS_PRESCRIPTION,                        "REMNANT_MASS_PRESCRIPTION" },

    { PROGRAM_OPTION::ROCKET_KICK_MAGNITUDE_1,                          "ROCKET_KICK_MAGNITUDE_1" },
    { PROGRAM_OPTION::ROCKET_KICK_MAGNITUDE_2,                          "ROCKET_KICK_MAGNITUDE_2" },
    { PROGRAM_OPTION::ROCKET_KICK_PHI_1,                                "ROCKET_KICK_PHI_1" },
    { PROGRAM_OPTION::ROCKET_KICK_PHI_2,                                "ROCKET_KICK_PHI_2" },
    { PROGRAM_OPTION::ROCKET_KICK_THETA_1,                              "ROCKET_KICK_THETA_1" },
    { PROGRAM_OPTION::ROCKET_KICK_THETA_2,                              "ROCKET_KICK_THETA_2" },

    { PROGRAM_OPTION::ROTATIONAL_VELOCITY_DISTRIBUTION,                 "ROTATIONAL_VELOCITY_DISTRIBUTION" },
    { PROGRAM_OPTION::ROTATIONAL_FREQUENCY,                             "ROTATIONAL_FREQUENCY" },
    { PROGRAM_OPTION::ROTATIONAL_FREQUENCY_1,                           "ROTATIONAL_FREQUENCY_1" },
    { PROGRAM_OPTION::ROTATIONAL_FREQUENCY_2,                           "ROTATIONAL_FREQUENCY_2" },
   
    { PROGRAM_OPTION::SCALE_CHE_MASS_LOSS_SURF_HE_ABUNDANCE,            "SCALE_CHE_MASS_LOSS_SURF_HE_ABUNDANCE" }, 
    { PROGRAM_OPTION::SCALE_TERMINAL_WIND_VEL_METALLICITY_POWER,        "SCALE_TERMINAL_WIND_VEL_METALLICITY_POWER" }, 
    { PROGRAM_OPTION::SEMI_MAJOR_AXIS,                                  "SEMI_MAJOR_AXIS" },
    { PROGRAM_OPTION::SEMI_MAJOR_AXIS_DISTRIBUTION,                     "SEMI_MAJOR_AXIS_DISTRIBUTION" },
    { PROGRAM_OPTION::SEMI_MAJOR_AXIS_DISTRIBUTION_MAX,                 "SEMI_MAJOR_AXIS_DISTRIBUTION_MAX" },
    { PROGRAM_OPTION::SEMI_MAJOR_AXIS_DISTRIBUTION_MIN,                 "SEMI_MAJOR_AXIS_DISTRIBUTION_MIN" },

    { PROGRAM_OPTION::STELLAR_ZETA_PRESCRIPTION,                        "STELLAR_ZETA_PRESCRIPTION" },

    { PROGRAM_OPTION::TIDES_PRESCRIPTION,                               "TIDES_PRESCRIPTION" },

    { PROGRAM_OPTION::WR_FACTOR,                                        "WR_FACTOR" },

    { PROGRAM_OPTION::ZETA_ADIABATIC_ARBITRARY,                         "ZETA_ADIABATIC_ARBITRARY" },
    { PROGRAM_OPTION::ZETA_MS,                                          "ZETA_MS" },
    { PROGRAM_OPTION::ZETA_RADIATIVE_ENVELOPE_GIANT,                    "ZETA_RADIATIVE_ENVELOPE_GIANT" }
};


enum class ANY_PROPERTY_TYPE: int { T_STAR_1_PROPERTY, T_STAR_2_PROPERTY, T_SUPERNOVA_PROPERTY, T_COMPANION_PROPERTY, T_STAR_PROPERTY, T_BINARY_PROPERTY, T_PROGRAM_OPTION };

typedef boost::variant<STAR_1_PROPERTY, STAR_2_PROPERTY, SUPERNOVA_PROPERTY, COMPANION_PROPERTY, STAR_PROPERTY, BINARY_PROPERTY, PROGRAM_OPTION> T_ANY_PROPERTY;
typedef boost::variant<STAR_1_PROPERTY, STAR_2_PROPERTY, SUPERNOVA_PROPERTY, COMPANION_PROPERTY, STAR_PROPERTY                                 > T_ANY_STAR_PROPERTY;
typedef boost::variant<STAR_1_PROPERTY, STAR_2_PROPERTY, SUPERNOVA_PROPERTY, COMPANION_PROPERTY                                                > T_ANY_BINARY_CONSTITUENT_PROPERTY;
typedef boost::variant<STAR_1_PROPERTY, STAR_2_PROPERTY, SUPERNOVA_PROPERTY, COMPANION_PROPERTY,                BINARY_PROPERTY                > T_ANY_BINARY_PROPERTY;


typedef std::vector<T_ANY_PROPERTY> ANY_PROPERTY_VECTOR;


class VariantPropertyType: public boost::static_visitor<ANY_PROPERTY_TYPE> {
public:
  ANY_PROPERTY_TYPE operator()(STAR_PROPERTY   prop)    const { return ANY_PROPERTY_TYPE::T_STAR_PROPERTY;   }
  ANY_PROPERTY_TYPE operator()(STAR_1_PROPERTY prop)    const { return ANY_PROPERTY_TYPE::T_STAR_1_PROPERTY; }
  ANY_PROPERTY_TYPE operator()(STAR_2_PROPERTY prop)    const { return ANY_PROPERTY_TYPE::T_STAR_2_PROPERTY; }
  ANY_PROPERTY_TYPE operator()(SUPERNOVA_PROPERTY prop) const { return ANY_PROPERTY_TYPE::T_SUPERNOVA_PROPERTY; }
  ANY_PROPERTY_TYPE operator()(COMPANION_PROPERTY prop) const { return ANY_PROPERTY_TYPE::T_COMPANION_PROPERTY; }
  ANY_PROPERTY_TYPE operator()(BINARY_PROPERTY prop)    const { return ANY_PROPERTY_TYPE::T_BINARY_PROPERTY; }
  ANY_PROPERTY_TYPE operator()(PROGRAM_OPTION prop)     const { return ANY_PROPERTY_TYPE::T_PROGRAM_OPTION; }
};


// Property types
enum class PROPERTY_TYPE: int { NONE, STAR_PROPERTY, STAR_1_PROPERTY, STAR_2_PROPERTY, SUPERNOVA_PROPERTY, COMPANION_PROPERTY, ANY_STAR_PROPERTY, BINARY_PROPERTY, PROGRAM_OPTION };
const COMPASUnorderedMap<PROPERTY_TYPE, std::string> PROPERTY_TYPE_LABEL = {
    { PROPERTY_TYPE::NONE,               "" },
    { PROPERTY_TYPE::STAR_PROPERTY,      "STAR_PROPERTY" },
    { PROPERTY_TYPE::STAR_1_PROPERTY,    "STAR_1_PROPERTY" },
    { PROPERTY_TYPE::STAR_2_PROPERTY,    "STAR_2_PROPERTY" },
    { PROPERTY_TYPE::SUPERNOVA_PROPERTY, "SUPERNOVA_PROPERTY" },
    { PROPERTY_TYPE::COMPANION_PROPERTY, "COMPANION_PROPERTY" },
    { PROPERTY_TYPE::ANY_STAR_PROPERTY,  "ANY_STAR_PROPERTY" },
    { PROPERTY_TYPE::BINARY_PROPERTY,    "BINARY_PROPERTY" },
    { PROPERTY_TYPE::PROGRAM_OPTION,     "PROGRAM_OPTION" }
};


// typedef for property details
// the property details provide information about the property to assist in printing the property value
// the property_details tuple contains:
//
//   <data type, header string, units string, output field width, output field precision>
//
// where:
//
//    data type       is, as the name suggest, the data type of the property, from the TYPENAME enum class (above)
//    header string   is the string to be printed as the column header for the property
//    units string    is the string to be printed as the units header for the property
//    field width     is the printf() field width of the property (does not apply to HDF5 files; meaning varies per the data type, refer to printf() documentation)
//    field precision is the printf() field precision of the property (does not apply to HDF5 files; meaning varies per the data type, refer to printf() documentation)
typedef std::tuple<TYPENAME, std::string, std::string, int, int> PROPERTY_DETAILS;

// map ANY_STAR_PROPERTY_DETAIL
// Records the details of STELLAR properties.  The STELLAR properties are those that pertain
// to individual stars, whether they be a single star being evolved for SSE, or one of the
// constituent stars being evolved as part of a binary for BSE
//
// Properties only need to be here if they are required to be available for printing in 
// the logfiles - all keys present here should also be in the STAR_PROPERTIES #define and
// STAR_PROPERTIES_LABEL
const std::map<ANY_STAR_PROPERTY, PROPERTY_DETAILS> ANY_STAR_PROPERTY_DETAIL = {
    { ANY_STAR_PROPERTY::AGE,                                               { TYPENAME::DOUBLE,           "Age",                             "Myr",              24, 15}},
    { ANY_STAR_PROPERTY::ANGULAR_MOMENTUM,                                  { TYPENAME::DOUBLE,           "Ang_Momentum",                    "Msol AU^2 yr^-1",  24, 15}},
    { ANY_STAR_PROPERTY::BINDING_ENERGY_AT_COMMON_ENVELOPE,                 { TYPENAME::DOUBLE,           "Binding_Energy@CE",               "erg",              24, 15}},
    { ANY_STAR_PROPERTY::BINDING_ENERGY_FIXED,                              { TYPENAME::DOUBLE,           "BE_Fixed",                        "erg",              24, 15}},
    { ANY_STAR_PROPERTY::BINDING_ENERGY_NANJING,                            { TYPENAME::DOUBLE,           "BE_Nanjing",                      "erg",              24, 15}},
    { ANY_STAR_PROPERTY::BINDING_ENERGY_PRE_COMMON_ENVELOPE,                { TYPENAME::DOUBLE,           "Binding_Energy<CE",               "erg",              24, 15}},
    { ANY_STAR_PROPERTY::BINDING_ENERGY_LOVERIDGE,                          { TYPENAME::DOUBLE,           "BE_Loveridge",                    "erg",              24, 15}},
    { ANY_STAR_PROPERTY::BINDING_ENERGY_LOVERIDGE_WINDS,                    { TYPENAME::DOUBLE,           "BE_Loveridge_Winds",              "erg",              24, 15}},
    { ANY_STAR_PROPERTY::BINDING_ENERGY_KRUCKOW,                            { TYPENAME::DOUBLE,           "BE_Kruckow",                      "erg",              24, 15}},
    { ANY_STAR_PROPERTY::CHEMICALLY_HOMOGENEOUS_MAIN_SEQUENCE,              { TYPENAME::BOOL,             "CH_on_MS",                        "State",             0, 0 }},
    { ANY_STAR_PROPERTY::CO_CORE_MASS,                                      { TYPENAME::DOUBLE,           "Mass_CO_Core",                    "Msol",             24, 15}},
    { ANY_STAR_PROPERTY::CO_CORE_MASS_AT_COMMON_ENVELOPE,                   { TYPENAME::DOUBLE,           "Mass_CO_Core@CE",                 "Msol",             24, 15}},
    { ANY_STAR_PROPERTY::CO_CORE_MASS_AT_COMPACT_OBJECT_FORMATION,          { TYPENAME::DOUBLE,           "Mass_CO_Core@CO",                 "Msol",             24, 15}},
    { ANY_STAR_PROPERTY::CORE_MASS,                                         { TYPENAME::DOUBLE,           "Mass_Core",                       "Msol",             24, 15}},
    { ANY_STAR_PROPERTY::CORE_MASS_AT_COMMON_ENVELOPE,                      { TYPENAME::DOUBLE,           "Mass_Core@CE",                    "Msol",             24, 15}},
    { ANY_STAR_PROPERTY::CORE_MASS_AT_COMPACT_OBJECT_FORMATION,             { TYPENAME::DOUBLE,           "Mass_Core@CO",                    "Msol",             24, 15}},
    { ANY_STAR_PROPERTY::DRAWN_KICK_MAGNITUDE,                              { TYPENAME::DOUBLE,           "Drawn_Kick_Magnitude",            "kms^-1",           24, 15}},
    { ANY_STAR_PROPERTY::DOMINANT_MASS_LOSS_RATE,                           { TYPENAME::INT,              "Dominant_Mass_Loss_Rate",         "-",                 4, 1 }},
    { ANY_STAR_PROPERTY::DT,                                                { TYPENAME::DOUBLE,           "dT",                              "Myr",              24, 15}},
    { ANY_STAR_PROPERTY::DYNAMICAL_TIMESCALE,                               { TYPENAME::DOUBLE,           "Tau_Dynamical",                   "Myr",              24, 15}},
    { ANY_STAR_PROPERTY::DYNAMICAL_TIMESCALE_POST_COMMON_ENVELOPE,          { TYPENAME::DOUBLE,           "Tau_Dynamical>CE",                "Myr",              24, 15}},
    { ANY_STAR_PROPERTY::DYNAMICAL_TIMESCALE_PRE_COMMON_ENVELOPE,           { TYPENAME::DOUBLE,           "Tau_Dynamical<CE",                "Myr",              24, 15}},
    { ANY_STAR_PROPERTY::ECCENTRIC_ANOMALY,                                 { TYPENAME::DOUBLE,           "Eccentric_Anomaly",               "-",                24, 15}},
    { ANY_STAR_PROPERTY::ENV_MASS,                                          { TYPENAME::DOUBLE,           "Mass_Env",                        "Msol",             24, 15}},
    { ANY_STAR_PROPERTY::ERROR,                                             { TYPENAME::ERROR,            "Error",                           "-",                 4, 1 }},
    { ANY_STAR_PROPERTY::EVOL_STATUS,                                       { TYPENAME::EVOLUTION_STATUS, "Evolution_Status",                "-",                 4, 1 }},
    { ANY_STAR_PROPERTY::EXPERIENCED_AIC,                                   { TYPENAME::BOOL,             "Experienced_AIC",                 "Event",             0, 0 }},
    { ANY_STAR_PROPERTY::EXPERIENCED_CCSN,                                  { TYPENAME::BOOL,             "Experienced_CCSN",                "Event",             0, 0 }},
    { ANY_STAR_PROPERTY::EXPERIENCED_HeSD,                                  { TYPENAME::BOOL,             "Experienced_HeSD",                "Event",             0, 0 }},
    { ANY_STAR_PROPERTY::EXPERIENCED_ECSN,                                  { TYPENAME::BOOL,             "Experienced_ECSN",                "Event",             0, 0 }},
    { ANY_STAR_PROPERTY::EXPERIENCED_PISN,                                  { TYPENAME::BOOL,             "Experienced_PISN",                "Event",             0, 0 }},
    { ANY_STAR_PROPERTY::EXPERIENCED_PPISN,                                 { TYPENAME::BOOL,             "Experienced_PPISN",               "Event",             0, 0 }},
    { ANY_STAR_PROPERTY::EXPERIENCED_RLOF,                                  { TYPENAME::BOOL,             "Experienced_RLOF",                "Event",             0, 0 }},
    { ANY_STAR_PROPERTY::EXPERIENCED_SNIA,                                  { TYPENAME::BOOL,             "Experienced_SNIA",                "Event",             0, 0 }},
    { ANY_STAR_PROPERTY::EXPERIENCED_SN_TYPE,                               { TYPENAME::SN_EVENT,         "Experienced_SN_Type",             "-",                 4, 1 }},
    { ANY_STAR_PROPERTY::EXPERIENCED_USSN,                                  { TYPENAME::BOOL,             "Experienced_USSN",                "Event",             0, 0 }},
    { ANY_STAR_PROPERTY::FALLBACK_FRACTION,                                 { TYPENAME::DOUBLE,           "Fallback_Fraction",               "-",                24, 15}},
    { ANY_STAR_PROPERTY::HE_CORE_MASS,                                      { TYPENAME::DOUBLE,           "Mass_He_Core",                    "Msol",             24, 15}},
    { ANY_STAR_PROPERTY::HE_CORE_MASS_AT_COMMON_ENVELOPE,                   { TYPENAME::DOUBLE,           "Mass_He_Core@CE",                 "Msol",             24, 15}},
    { ANY_STAR_PROPERTY::HE_CORE_MASS_AT_COMPACT_OBJECT_FORMATION,          { TYPENAME::DOUBLE,           "Mass_He_Core@CO",                 "Msol",             24, 15}},
    { ANY_STAR_PROPERTY::HELIUM_ABUNDANCE_SURFACE,                          { TYPENAME::DOUBLE,           "Helium_Abundance_Surface",        "-",                24, 15}},
    { ANY_STAR_PROPERTY::HELIUM_ABUNDANCE_CORE,                             { TYPENAME::DOUBLE,           "Helium_Abundance_Core",           "-",                24, 15}},
    { ANY_STAR_PROPERTY::HYDROGEN_ABUNDANCE_SURFACE,                        { TYPENAME::DOUBLE,           "Hydrogen_Abundance_Surface",      "-",                24, 15}},
    { ANY_STAR_PROPERTY::HYDROGEN_ABUNDANCE_CORE,                           { TYPENAME::DOUBLE,           "Hydrogen_Abundance_Core",         "-",                24, 15}},
    { ANY_STAR_PROPERTY::ID,                                                { TYPENAME::OBJECT_ID,        "ID",                              "-",                12, 1 }},
    { ANY_STAR_PROPERTY::INITIAL_HELIUM_ABUNDANCE,                          { TYPENAME::DOUBLE,           "Helium_Abundance@ZAMS",           "-",                24, 15}},
    { ANY_STAR_PROPERTY::INITIAL_HYDROGEN_ABUNDANCE,                        { TYPENAME::DOUBLE,           "Hydrogen_Abundance@ZAMS",         "-",                24, 15}},
    { ANY_STAR_PROPERTY::INITIAL_STELLAR_TYPE,                              { TYPENAME::STELLAR_TYPE,     "Stellar_Type@ZAMS",               "-",                 4, 1 }},
    { ANY_STAR_PROPERTY::INITIAL_STELLAR_TYPE_NAME,                         { TYPENAME::STRING,           "Stellar_Type@ZAMS",               "-",                42, 1 }},
    { ANY_STAR_PROPERTY::IS_AIC,                                            { TYPENAME::BOOL,             "AIC",                             "State",             0, 0 }},
    { ANY_STAR_PROPERTY::IS_CCSN,                                           { TYPENAME::BOOL,             "CCSN",                            "State",             0, 0 }},
    { ANY_STAR_PROPERTY::IS_HeSD,                                           { TYPENAME::BOOL,             "HeSD",                            "State",             0, 0 }},
    { ANY_STAR_PROPERTY::IS_ECSN,                                           { TYPENAME::BOOL,             "ECSN",                            "State",             0, 0 }},
    { ANY_STAR_PROPERTY::IS_HYDROGEN_POOR,                                  { TYPENAME::BOOL,             "Is_Hydrogen_Poor",                "State",             0, 0 }},
    { ANY_STAR_PROPERTY::IS_PISN,                                           { TYPENAME::BOOL,             "PISN",                            "State",             0, 0 }},
    { ANY_STAR_PROPERTY::IS_PPISN,                                          { TYPENAME::BOOL,             "PPISN",                           "State",             0, 0 }},
    { ANY_STAR_PROPERTY::IS_RLOF,                                           { TYPENAME::BOOL,             "RLOF",                            "State",             0, 0 }},
    { ANY_STAR_PROPERTY::IS_SNIA,                                           { TYPENAME::BOOL,             "SNIA",                            "State",             0, 0 }},
    { ANY_STAR_PROPERTY::IS_USSN,                                           { TYPENAME::BOOL,             "USSN",                            "State",             0, 0 }},
    { ANY_STAR_PROPERTY::KICK_MAGNITUDE,                                    { TYPENAME::DOUBLE,           "Applied_Kick_Magnitude",          "kms^-1",           24, 15}},
    { ANY_STAR_PROPERTY::LAMBDA_AT_COMMON_ENVELOPE,                         { TYPENAME::DOUBLE,           "Lambda@CE",                       "-",                24, 15}},
    { ANY_STAR_PROPERTY::LAMBDA_DEWI,                                       { TYPENAME::DOUBLE,           "Lambda_Dewi",                     "-",                24, 15}},
    { ANY_STAR_PROPERTY::LAMBDA_FIXED,                                      { TYPENAME::DOUBLE,           "Lambda_Fixed",                    "-",                24, 15}},
    { ANY_STAR_PROPERTY::LAMBDA_KRUCKOW,                                    { TYPENAME::DOUBLE,           "Lambda_Kruckow",                  "-",                24, 15}},
    { ANY_STAR_PROPERTY::LAMBDA_KRUCKOW_BOTTOM,                             { TYPENAME::DOUBLE,           "Lambda_Kruckow_Bottom",           "-",                24, 15}},
    { ANY_STAR_PROPERTY::LAMBDA_KRUCKOW_MIDDLE,                             { TYPENAME::DOUBLE,           "Lambda_Kruckow_Middle",           "-",                24, 15}},
    { ANY_STAR_PROPERTY::LAMBDA_KRUCKOW_TOP,                                { TYPENAME::DOUBLE,           "Lambda_Kruckow_Top",              "-",                24, 15}},
    { ANY_STAR_PROPERTY::LAMBDA_LOVERIDGE,                                  { TYPENAME::DOUBLE,           "Lambda_Loveridge",                "-",                24, 15}},
    { ANY_STAR_PROPERTY::LAMBDA_LOVERIDGE_WINDS,                            { TYPENAME::DOUBLE,           "Lambda_Loveridge_Winds",          "-",                24, 15}},
    { ANY_STAR_PROPERTY::LAMBDA_NANJING,                                    { TYPENAME::DOUBLE,           "Lambda_Nanjing",                  "-",                24, 15}},
    { ANY_STAR_PROPERTY::LBV_PHASE_FLAG,                                    { TYPENAME::BOOL,             "LBV_Phase_Flag",                  "Event",             0, 0 }},
    { ANY_STAR_PROPERTY::LUMINOSITY,                                        { TYPENAME::DOUBLE,           "Luminosity",                      "Lsol",             24, 15}},
    { ANY_STAR_PROPERTY::LUMINOSITY_POST_COMMON_ENVELOPE,                   { TYPENAME::DOUBLE,           "Luminosity>CE",                   "Lsol",             24, 15}},
    { ANY_STAR_PROPERTY::LUMINOSITY_PRE_COMMON_ENVELOPE,                    { TYPENAME::DOUBLE,           "Luminosity<CE",                   "Lsol",             24, 15}},
    { ANY_STAR_PROPERTY::MASS,                                              { TYPENAME::DOUBLE,           "Mass",                            "Msol",             24, 15}},
    { ANY_STAR_PROPERTY::MASS_0,                                            { TYPENAME::DOUBLE,           "Mass_0",                          "Msol",             24, 15}},
    { ANY_STAR_PROPERTY::MASS_LOSS_DIFF,                                    { TYPENAME::DOUBLE,           "dmWinds",                         "Msol",             24, 15}},
    { ANY_STAR_PROPERTY::MASS_TRANSFER_DIFF,                                { TYPENAME::DOUBLE,           "dmMT",                            "Msol",             24, 15}},
    { ANY_STAR_PROPERTY::MASS_TRANSFER_DONOR_HISTORY,                       { TYPENAME::STRING,           "MT_Donor_Hist",                   "-",                16, 1 }}, 
    { ANY_STAR_PROPERTY::MDOT,                                              { TYPENAME::DOUBLE,           "Mdot",                            "Msol yr^-1",       24, 15}},
    { ANY_STAR_PROPERTY::METALLICITY,                                       { TYPENAME::DOUBLE,           "Metallicity@ZAMS",                "-",                24, 15}},
    { ANY_STAR_PROPERTY::MOMENT_OF_INERTIA,                                 { TYPENAME::DOUBLE,           "Moment_Of_Inertia",               "Msol Rsol^2",      24, 15}},
    { ANY_STAR_PROPERTY::MZAMS,                                             { TYPENAME::DOUBLE,           "Mass@ZAMS",                       "Msol",             24, 15}},
    { ANY_STAR_PROPERTY::OMEGA,                                             { TYPENAME::DOUBLE,           "Omega",                           "Hz",               24, 15}},
    { ANY_STAR_PROPERTY::OMEGA_BREAK,                                       { TYPENAME::DOUBLE,           "Omega_Break",                     "Hz",               24, 15}},
    { ANY_STAR_PROPERTY::OMEGA_ZAMS,                                        { TYPENAME::DOUBLE,           "Omega@ZAMS",                      "Hz",               24, 15}},
    { ANY_STAR_PROPERTY::ORBITAL_ENERGY_POST_SUPERNOVA,                     { TYPENAME::DOUBLE,           "Orbital_Energy>SN",               "Msol^2AU^-1",      24, 15}},
    { ANY_STAR_PROPERTY::ORBITAL_ENERGY_PRE_SUPERNOVA,                      { TYPENAME::DOUBLE,           "Orbital_Energy<SN",               "Msol^2AU^-1",      24, 15}},
    { ANY_STAR_PROPERTY::PULSAR_MAGNETIC_FIELD,                             { TYPENAME::DOUBLE,           "Pulsar_Mag_Field",                "Tesla",            24, 15}},
    { ANY_STAR_PROPERTY::PULSAR_SPIN_DOWN_RATE,                             { TYPENAME::DOUBLE,           "Pulsar_Spin_Down",                "rad/s^2",          24, 15}},
    { ANY_STAR_PROPERTY::PULSAR_BIRTH_PERIOD,                               { TYPENAME::DOUBLE,           "Pulsar_Birth_Period",             "s",                24, 15}},
    { ANY_STAR_PROPERTY::PULSAR_BIRTH_SPIN_DOWN_RATE,                       { TYPENAME::DOUBLE,           "Pulsar_Birth_Spin_Down",          "s/s",              24, 15}},
    { ANY_STAR_PROPERTY::PULSAR_SPIN_FREQUENCY,                             { TYPENAME::DOUBLE,           "Pulsar_Spin_Freq",                "rad/s",            24, 15}},
    { ANY_STAR_PROPERTY::PULSAR_SPIN_PERIOD,                                { TYPENAME::DOUBLE,           "Pulsar_Spin_Period",              "ms",               24, 15}},
    { ANY_STAR_PROPERTY::RADIAL_EXPANSION_TIMESCALE,                        { TYPENAME::DOUBLE,           "Tau_Radial",                      "Myr",              24, 15}},
    { ANY_STAR_PROPERTY::RADIAL_EXPANSION_TIMESCALE_POST_COMMON_ENVELOPE,   { TYPENAME::DOUBLE,           "Tau_Radial>CE",                   "Myr",              24, 15}},
    { ANY_STAR_PROPERTY::RADIAL_EXPANSION_TIMESCALE_PRE_COMMON_ENVELOPE,    { TYPENAME::DOUBLE,           "Tau_Radial<CE",                   "Myr",              24, 15}},
    { ANY_STAR_PROPERTY::RADIUS,                                            { TYPENAME::DOUBLE,           "Radius",                          "Rsol",             24, 15}},
    { ANY_STAR_PROPERTY::RANDOM_SEED,                                       { TYPENAME::ULONGINT,         "SEED",                            "-",                12, 1 }},
    { ANY_STAR_PROPERTY::RECYCLED_NEUTRON_STAR,                             { TYPENAME::BOOL,             "Recycled_NS",                     "Event",             0, 0 }},
    { ANY_STAR_PROPERTY::RLOF_ONTO_NS,                                      { TYPENAME::BOOL,             "RLOF->NS",                        "Event",             0, 0 }},
    { ANY_STAR_PROPERTY::ROCKET_KICK_MAGNITUDE,                             { TYPENAME::DOUBLE,           "Rocket_Kick_Magnitude",           "kms^-1",           24, 15}},
    { ANY_STAR_PROPERTY::ROCKET_KICK_PHI,                                   { TYPENAME::DOUBLE,           "Rocket_Kick_Phi",                 "-",                24, 15}},
    { ANY_STAR_PROPERTY::ROCKET_KICK_THETA,                                 { TYPENAME::DOUBLE,           "Rocket_Kick_Theta",               "-",                24, 15}},
    { ANY_STAR_PROPERTY::RZAMS,                                             { TYPENAME::DOUBLE,           "Radius@ZAMS",                     "Rsol",             24, 15}},
    { ANY_STAR_PROPERTY::SN_TYPE,                                           { TYPENAME::SN_EVENT,         "SN_Type",                         "-",                 4, 1 }},
    { ANY_STAR_PROPERTY::SPEED,                                             { TYPENAME::DOUBLE,           "ComponentSpeed",                  "kms^-1",           24, 15}},
    { ANY_STAR_PROPERTY::STELLAR_TYPE,                                      { TYPENAME::STELLAR_TYPE,     "Stellar_Type",                    "-",                 4, 1 }},
    { ANY_STAR_PROPERTY::STELLAR_TYPE_NAME,                                 { TYPENAME::STRING,           "Stellar_Type",                    "-",                42, 1 }},
    { ANY_STAR_PROPERTY::STELLAR_TYPE_PREV,                                 { TYPENAME::STELLAR_TYPE,     "Stellar_Type_Prev",               "-",                 4, 1 }},
    { ANY_STAR_PROPERTY::STELLAR_TYPE_PREV_NAME,                            { TYPENAME::STRING,           "Stellar_Type_Prev",               "-",                42, 1 }},
    { ANY_STAR_PROPERTY::SUPERNOVA_KICK_MAGNITUDE_RANDOM_NUMBER,            { TYPENAME::DOUBLE,           "SN_Kick_Magnitude_Random_Number", "-",                24, 15}},
    { ANY_STAR_PROPERTY::MEAN_ANOMALY,                                      { TYPENAME::DOUBLE,           "SN_Kick_Mean_Anomaly",            "-",                24, 15}},
    { ANY_STAR_PROPERTY::SUPERNOVA_PHI,                                     { TYPENAME::DOUBLE,           "SN_Kick_Phi",                     "-",                24, 15}},
    { ANY_STAR_PROPERTY::SUPERNOVA_THETA,                                   { TYPENAME::DOUBLE,           "SN_Kick_Theta",                   "-",                24, 15}},
    { ANY_STAR_PROPERTY::TEMPERATURE,                                       { TYPENAME::DOUBLE,           "Teff",                            "K",                24, 15}},
    { ANY_STAR_PROPERTY::TEMPERATURE_POST_COMMON_ENVELOPE,                  { TYPENAME::DOUBLE,           "Teff>CE",                         "K",                24, 15}},
    { ANY_STAR_PROPERTY::TEMPERATURE_PRE_COMMON_ENVELOPE,                   { TYPENAME::DOUBLE,           "Teff<CE",                         "K",                24, 15}},
    { ANY_STAR_PROPERTY::THERMAL_TIMESCALE,                                 { TYPENAME::DOUBLE,           "Tau_Thermal",                     "Myr",              24, 15}},
    { ANY_STAR_PROPERTY::THERMAL_TIMESCALE_POST_COMMON_ENVELOPE,            { TYPENAME::DOUBLE,           "Tau_Thermal>CE",                  "Myr",              24, 15}},
    { ANY_STAR_PROPERTY::THERMAL_TIMESCALE_PRE_COMMON_ENVELOPE,             { TYPENAME::DOUBLE,           "Tau_Thermal<CE",                  "Myr",              24, 15}},
    { ANY_STAR_PROPERTY::TIME,                                              { TYPENAME::DOUBLE,           "Time",                            "Myr",              24, 15 }},
    { ANY_STAR_PROPERTY::TIMESCALE_MS,                                      { TYPENAME::DOUBLE,           "tMS",                             "Myr",              24, 15}},
    { ANY_STAR_PROPERTY::TOTAL_MASS_AT_COMPACT_OBJECT_FORMATION,            { TYPENAME::DOUBLE,           "Mass_Total@CO",                   "Msol",             24, 15}},
    { ANY_STAR_PROPERTY::TRUE_ANOMALY,                                      { TYPENAME::DOUBLE,           "True_Anomaly(psi)",               "-",                24, 15}},
    { ANY_STAR_PROPERTY::TZAMS,                                             { TYPENAME::DOUBLE,           "Teff@ZAMS",                       "K",                24, 15}},
    { ANY_STAR_PROPERTY::ZETA_HURLEY,                                       { TYPENAME::DOUBLE,           "Zeta_Hurley",                     "-",                24, 15}},
    { ANY_STAR_PROPERTY::ZETA_HURLEY_HE,                                    { TYPENAME::DOUBLE,           "Zeta_Hurley_He",                  "-",                24, 15}},
    { ANY_STAR_PROPERTY::ZETA_SOBERMAN,                                     { TYPENAME::DOUBLE,           "Zeta_Soberman",                   "-",                24, 15}},
    { ANY_STAR_PROPERTY::ZETA_SOBERMAN_HE,                                  { TYPENAME::DOUBLE,           "Zeta_Soberman_He",                "-",                24, 15}}
};

// map BINARY_PROPERTY_DETAIL
// Records the details of BINARY properties.  The BINARY properties are those that pertain
// to exclusively to a binary star - not the constituent stars that make up the binary
//
// Properties only need to be here if they are required to be available for printing in 
// the logfiles - all keys present here should also be in BINARY_PROPERTY and BINARY_PROPERTY_LABEL
const std::map<BINARY_PROPERTY, PROPERTY_DETAILS> BINARY_PROPERTY_DETAIL = {
    { BINARY_PROPERTY::CIRCULARIZATION_TIMESCALE,                           { TYPENAME::DOUBLE,           "Tau_Circ",                  "Myr",              24, 15}},
    { BINARY_PROPERTY::COMMON_ENVELOPE_AT_LEAST_ONCE,                       { TYPENAME::BOOL,             "CEE",                       "Event",             0, 0 }},
    { BINARY_PROPERTY::COMMON_ENVELOPE_EVENT_COUNT,                         { TYPENAME::UINT,             "CE_Event_Counter",          "Count",             6, 1 }},
    { BINARY_PROPERTY::DOUBLE_CORE_COMMON_ENVELOPE,                         { TYPENAME::BOOL,             "Double_Core_CE",            "Event",             0, 0 }},
    { BINARY_PROPERTY::DT,                                                  { TYPENAME::DOUBLE,           "dT",                        "Myr",              24, 15}},
    { BINARY_PROPERTY::ECCENTRICITY,                                        { TYPENAME::DOUBLE,           "Eccentricity",              "-",                24, 15}},
    { BINARY_PROPERTY::ECCENTRICITY_AT_DCO_FORMATION,                       { TYPENAME::DOUBLE,           "Eccentricity@DCO",          "-",                24, 15}},
    { BINARY_PROPERTY::ECCENTRICITY_INITIAL,                                { TYPENAME::DOUBLE,           "Eccentricity@ZAMS",         "-",                24, 15}},
    { BINARY_PROPERTY::ECCENTRICITY_POST_COMMON_ENVELOPE,                   { TYPENAME::DOUBLE,           "Eccentricity>CE",           "-",                24, 15}},
    { BINARY_PROPERTY::ECCENTRICITY_PRE_SUPERNOVA,                          { TYPENAME::DOUBLE,           "Eccentricity<SN",           "-",                24, 15}},
    { BINARY_PROPERTY::ECCENTRICITY_PRE_COMMON_ENVELOPE,                    { TYPENAME::DOUBLE,           "Eccentricity<CE",           "-",                24, 15}},
    { BINARY_PROPERTY::ERROR,                                               { TYPENAME::ERROR,            "Error",                     "-",                 4, 1 }},
    { BINARY_PROPERTY::EVOL_STATUS,                                         { TYPENAME::EVOLUTION_STATUS, "Evolution_Status",          "-",                 4, 1 }},
    { BINARY_PROPERTY::ID,                                                  { TYPENAME::OBJECT_ID,        "ID",                        "-",                12, 1 }},
    { BINARY_PROPERTY::IMMEDIATE_RLOF_POST_COMMON_ENVELOPE,                 { TYPENAME::BOOL,             "Immediate_RLOF>CE",         "Event",             0, 0 }},
    { BINARY_PROPERTY::MASS_1_POST_COMMON_ENVELOPE,                         { TYPENAME::DOUBLE,           "Mass(1)>CE",                "Msol",             24, 15}},
    { BINARY_PROPERTY::MASS_1_PRE_COMMON_ENVELOPE,                          { TYPENAME::DOUBLE,           "Mass(1)<CE",                "Msol",             24, 15}},
    { BINARY_PROPERTY::MASS_2_POST_COMMON_ENVELOPE,                         { TYPENAME::DOUBLE,           "Mass(2)>CE",                "Msol",             24, 15}},
    { BINARY_PROPERTY::MASS_2_PRE_COMMON_ENVELOPE,                          { TYPENAME::DOUBLE,           "Mass(2)<CE",                "Msol",             24, 15}},
    { BINARY_PROPERTY::MASS_ENV_1,                                          { TYPENAME::DOUBLE,           "Mass_Env(1)",               "Msol",             24, 15}},
    { BINARY_PROPERTY::MASS_ENV_2,                                          { TYPENAME::DOUBLE,           "Mass_Env(2)",               "Msol",             24, 15}},
    { BINARY_PROPERTY::MASSES_EQUILIBRATED,                                 { TYPENAME::BOOL,             "Equilibrated",              "Event",             0, 0 }},
    { BINARY_PROPERTY::MASSES_EQUILIBRATED_AT_BIRTH,                        { TYPENAME::BOOL,             "Equilibrated_At_Birth",     "Event",             0, 0 }},
    { BINARY_PROPERTY::MASS_TRANSFER_TRACKER_HISTORY,                       { TYPENAME::MT_TRACKING,      "MT_History",                "-",                 4, 1 }},
    { BINARY_PROPERTY::MERGES_IN_HUBBLE_TIME,                               { TYPENAME::BOOL,             "Merges_Hubble_Time",        "State",             0, 0 }},
    { BINARY_PROPERTY::OPTIMISTIC_COMMON_ENVELOPE,                          { TYPENAME::BOOL,             "Optimistic_CE",             "State",             0, 0 }},
    { BINARY_PROPERTY::ORBITAL_ANGULAR_VELOCITY,                            { TYPENAME::DOUBLE,           "Orbital_Angular_Velocity",  "kms^-1",           24, 15}},
    { BINARY_PROPERTY::ORBITAL_VELOCITY_PRE_SUPERNOVA,                      { TYPENAME::DOUBLE,           "Orb_Velocity<SN",           "kms^-1",           24, 15}},
    { BINARY_PROPERTY::RADIUS_1_POST_COMMON_ENVELOPE,                       { TYPENAME::DOUBLE,           "Radius(1)>CE",              "Rsol",             24, 15}},
    { BINARY_PROPERTY::RADIUS_1_PRE_COMMON_ENVELOPE,                        { TYPENAME::DOUBLE,           "Radius(1)<CE",              "Rsol",             24, 15}},
    { BINARY_PROPERTY::RADIUS_2_POST_COMMON_ENVELOPE,                       { TYPENAME::DOUBLE,           "Radius(2)>CE",              "Rsol",             24, 15}},
    { BINARY_PROPERTY::RADIUS_2_PRE_COMMON_ENVELOPE,                        { TYPENAME::DOUBLE,           "Radius(2)<CE",              "Rsol",             24, 15}},
    { BINARY_PROPERTY::RANDOM_SEED,                                         { TYPENAME::ULONGINT,         "SEED",                      "-",                12, 1 }},
    { BINARY_PROPERTY::RLOF_ACCRETION_EFFICIENCY,                           { TYPENAME::DOUBLE,           "Beta",                      "-",                24, 15}},
    { BINARY_PROPERTY::RLOF_MASS_LOSS_RATE,                                 { TYPENAME::DOUBLE,           "MassTransferRateDonor",     "Msol/Myr",         24, 15}},
    { BINARY_PROPERTY::RLOF_MASS_TRANSFER_TIMESCALE,                        {     TYPENAME::MASS_TRANSFER_TIMESCALE,  "MassTransferTimescale",    "-",          4, 1}},
    { BINARY_PROPERTY::RLOF_POST_MT_COMMON_ENVELOPE,                        { TYPENAME::BOOL,             "CEE>MT",                    "State",             0, 0 }},
    { BINARY_PROPERTY::RLOF_POST_MT_ECCENTRICITY,                           { TYPENAME::DOUBLE,           "Eccentricity>MT",           "-",                24, 15}},
    { BINARY_PROPERTY::RLOF_POST_MT_EVENT_COUNTER,                          { TYPENAME::UINT,             "MT_Event_Counter",          "Count",             6, 1 }},
    { BINARY_PROPERTY::RLOF_POST_MT_ID,                                     { TYPENAME::OBJECT_ID,        "ID>MT",                     "-",                12, 1 }},
    { BINARY_PROPERTY::RLOF_POST_MT_SEMI_MAJOR_AXIS,                        { TYPENAME::DOUBLE,           "SemiMajorAxis>MT",          "Rsol",             24, 15}},
    { BINARY_PROPERTY::RLOF_POST_MT_STAR1_MASS,                             { TYPENAME::DOUBLE,           "Mass(1)>MT",                "Msol",             24, 15}},
    { BINARY_PROPERTY::RLOF_POST_MT_STAR2_MASS,                             { TYPENAME::DOUBLE,           "Mass(2)>MT",                "Msol",             24, 15}},
    { BINARY_PROPERTY::RLOF_POST_MT_STAR1_RADIUS,                           { TYPENAME::DOUBLE,           "Radius(1)>MT",              "Rsol",             24, 15}},
    { BINARY_PROPERTY::RLOF_POST_MT_STAR2_RADIUS,                           { TYPENAME::DOUBLE,           "Radius(2)>MT",              "Rsol",             24, 15}},
    { BINARY_PROPERTY::RLOF_POST_MT_STAR1_RLOF,                             { TYPENAME::BOOL,             "RLOF(1)>MT",                "State",             0, 0 }},
    { BINARY_PROPERTY::RLOF_POST_MT_STAR2_RLOF,                             { TYPENAME::BOOL,             "RLOF(2)>MT",                "State",             0, 0 }},
    { BINARY_PROPERTY::RLOF_POST_MT_STAR1_STELLAR_TYPE,                     { TYPENAME::STELLAR_TYPE,     "Stellar_Type(1)>MT",         "-",                4, 1 }},
    { BINARY_PROPERTY::RLOF_POST_MT_STAR1_STELLAR_TYPE_NAME,                { TYPENAME::STRING,           "Stellar_Type(1)>MT",         "-",               42, 1 }},
    { BINARY_PROPERTY::RLOF_POST_MT_STAR2_STELLAR_TYPE,                     { TYPENAME::STELLAR_TYPE,     "Stellar_Type(2)>MT",         "-",                4, 1 }},
    { BINARY_PROPERTY::RLOF_POST_MT_STAR2_STELLAR_TYPE_NAME,                { TYPENAME::STRING,           "Stellar_Type(2)>MT",         "-",               42, 1 }},
    { BINARY_PROPERTY::RLOF_POST_STEP_STAR_TO_ROCHE_LOBE_RADIUS_RATIO_1,    { TYPENAME::DOUBLE,           "Radius(1)|RL>step",          "-",               24, 15}},
    { BINARY_PROPERTY::RLOF_POST_STEP_STAR_TO_ROCHE_LOBE_RADIUS_RATIO_2,    { TYPENAME::DOUBLE,           "Radius(2)|RL>step",          "-",               24, 15}},
    { BINARY_PROPERTY::RLOF_PRE_MT_ECCENTRICITY,                            { TYPENAME::DOUBLE,           "Eccentricity<MT",            "-",               24, 15}},
    { BINARY_PROPERTY::RLOF_PRE_MT_SEMI_MAJOR_AXIS,                         { TYPENAME::DOUBLE,           "SemiMajorAxis<MT",           "Rsol",            24, 15}},
    { BINARY_PROPERTY::RLOF_PRE_MT_STAR1_MASS,                              { TYPENAME::DOUBLE,           "Mass(1)<MT",                 "Msol",            24, 15}},
    { BINARY_PROPERTY::RLOF_PRE_MT_STAR2_MASS,                              { TYPENAME::DOUBLE,           "Mass(2)<MT",                 "Msol",            24, 15}},
    { BINARY_PROPERTY::RLOF_PRE_MT_STAR1_RADIUS,                            { TYPENAME::DOUBLE,           "Radius(1)<MT",               "Rsol",            24, 15}},
    { BINARY_PROPERTY::RLOF_PRE_MT_STAR2_RADIUS,                            { TYPENAME::DOUBLE,           "Radius(2)<MT",               "Rsol",            24, 15}},
    { BINARY_PROPERTY::RLOF_PRE_MT_STAR1_RLOF,                              { TYPENAME::BOOL,             "RLOF(1)<MT",                 "Event",            0, 0 }},
    { BINARY_PROPERTY::RLOF_PRE_MT_STAR2_RLOF,                              { TYPENAME::BOOL,             "RLOF(2)<MT",                 "Event",            0, 0 }},
    { BINARY_PROPERTY::RLOF_PRE_MT_STAR1_STELLAR_TYPE,                      { TYPENAME::STELLAR_TYPE,     "Stellar_Type(1)<MT",         "-",                4, 1 }},
    { BINARY_PROPERTY::RLOF_PRE_MT_STAR1_STELLAR_TYPE_NAME,                 { TYPENAME::STRING,           "Stellar_Type(1)<MT",         "-",               42, 1 }},
    { BINARY_PROPERTY::RLOF_PRE_MT_STAR2_STELLAR_TYPE,                      { TYPENAME::STELLAR_TYPE,     "Stellar_Type(2)<MT",         "-",                4, 1 }},
    { BINARY_PROPERTY::RLOF_PRE_MT_STAR2_STELLAR_TYPE_NAME,                 { TYPENAME::STRING,           "Stellar_Type(2)<MT",         "-",               42, 1 }},
    { BINARY_PROPERTY::RLOF_PRE_STEP_STAR_TO_ROCHE_LOBE_RADIUS_RATIO_1,     { TYPENAME::DOUBLE,           "Radius(1)|RL<step",          "-",               24, 15}},
    { BINARY_PROPERTY::RLOF_PRE_STEP_STAR_TO_ROCHE_LOBE_RADIUS_RATIO_2,     { TYPENAME::DOUBLE,           "Radius(2)|RL<step",          "-",               24, 15}},
    { BINARY_PROPERTY::RLOF_SECONDARY_POST_COMMON_ENVELOPE,                 { TYPENAME::BOOL,             "RLOF_Secondary>CE",          "Event",            0, 0 }},
    { BINARY_PROPERTY::RLOF_TIME_POST_MT,                                   { TYPENAME::DOUBLE,           "Time>MT",                    "Myr",             24, 15}},
    { BINARY_PROPERTY::RLOF_TIME_PRE_MT,                                    { TYPENAME::DOUBLE,           "Time<MT",                    "Myr",             24, 15}},
    { BINARY_PROPERTY::ROCHE_LOBE_RADIUS_1,                                 { TYPENAME::DOUBLE,           "RocheLobe(1)",               "Rsol",            24, 15}},
    { BINARY_PROPERTY::ROCHE_LOBE_RADIUS_1_POST_COMMON_ENVELOPE,            { TYPENAME::DOUBLE,           "RocheLobe(1)>CE",            "Rsol",            24, 15}},
    { BINARY_PROPERTY::ROCHE_LOBE_RADIUS_1_PRE_COMMON_ENVELOPE,             { TYPENAME::DOUBLE,           "RocheLobe(1)<CE",            "Rsol",            24, 15}},
    { BINARY_PROPERTY::ROCHE_LOBE_RADIUS_2,                                 { TYPENAME::DOUBLE,           "RocheLobe(2)",               "Rsol",            24, 15}},
    { BINARY_PROPERTY::ROCHE_LOBE_RADIUS_2_POST_COMMON_ENVELOPE,            { TYPENAME::DOUBLE,           "RocheLobe(2)>CE",            "Rsol",            24, 15}},
    { BINARY_PROPERTY::ROCHE_LOBE_RADIUS_2_PRE_COMMON_ENVELOPE,             { TYPENAME::DOUBLE,           "RocheLobe(2)<CE",            "Rsol",            24, 15}},
    { BINARY_PROPERTY::STAR_TO_ROCHE_LOBE_RADIUS_RATIO_1,                   { TYPENAME::DOUBLE,           "Radius(1)|RL",               "-",               24, 15}},
    { BINARY_PROPERTY::STAR_TO_ROCHE_LOBE_RADIUS_RATIO_2,                   { TYPENAME::DOUBLE,           "Radius(2)|RL",               "-",               24, 15}},
    { BINARY_PROPERTY::SEMI_MAJOR_AXIS_AT_DCO_FORMATION,                    { TYPENAME::DOUBLE,           "SemiMajorAxis@DCO",          "AU",              24, 15}},
    { BINARY_PROPERTY::SEMI_MAJOR_AXIS_INITIAL,                             { TYPENAME::DOUBLE,           "SemiMajorAxis@ZAMS",         "AU",              24, 15}},
    { BINARY_PROPERTY::SEMI_MAJOR_AXIS_POST_COMMON_ENVELOPE,                { TYPENAME::DOUBLE,           "SemiMajorAxis>CE",           "Rsol",            24, 15}},
    { BINARY_PROPERTY::SEMI_MAJOR_AXIS_PRE_COMMON_ENVELOPE,                 { TYPENAME::DOUBLE,           "SemiMajorAxis<CE",           "Rsol",            24, 15}},
    { BINARY_PROPERTY::SEMI_MAJOR_AXIS_PRE_SUPERNOVA,                       { TYPENAME::DOUBLE,           "SemiMajorAxis<SN",           "AU",              24, 15}},
    { BINARY_PROPERTY::SEMI_MAJOR_AXIS_PRE_SUPERNOVA_RSOL,                  { TYPENAME::DOUBLE,           "SemiMajorAxis<SN",           "Rsol",            24, 15}},
    { BINARY_PROPERTY::SEMI_MAJOR_AXIS,                                     { TYPENAME::DOUBLE,           "SemiMajorAxis",              "AU",              24, 15}},
    { BINARY_PROPERTY::SEMI_MAJOR_AXIS_RSOL,                                { TYPENAME::DOUBLE,           "SemiMajorAxis",              "Rsol",            24, 15}},
    { BINARY_PROPERTY::SIMULTANEOUS_RLOF,                                   { TYPENAME::BOOL,             "Simultaneous_RLOF",          "Event",            0, 0 }},
    { BINARY_PROPERTY::STABLE_RLOF_POST_COMMON_ENVELOPE,                    { TYPENAME::BOOL,             "Stable_RLOF>CE",             "State",            0, 0 }},
    { BINARY_PROPERTY::STELLAR_MERGER,                                      { TYPENAME::BOOL,             "Merger",                     "Event",            0, 0 }},
    { BINARY_PROPERTY::STELLAR_MERGER_AT_BIRTH,                             { TYPENAME::BOOL,             "Merger_At_Birth",            "Event",            0, 0 }},
    { BINARY_PROPERTY::STELLAR_TYPE_1_POST_COMMON_ENVELOPE,                 { TYPENAME::STELLAR_TYPE,     "Stellar_Type(1)>CE",         "-",                4, 1 }},
    { BINARY_PROPERTY::STELLAR_TYPE_1_PRE_COMMON_ENVELOPE,                  { TYPENAME::STELLAR_TYPE,     "Stellar_Type(1)<CE",         "-",                4, 1 }},
    { BINARY_PROPERTY::STELLAR_TYPE_2_POST_COMMON_ENVELOPE,                 { TYPENAME::STELLAR_TYPE,     "Stellar_Type(2)>CE",         "-",                4, 1 }},
    { BINARY_PROPERTY::STELLAR_TYPE_2_PRE_COMMON_ENVELOPE,                  { TYPENAME::STELLAR_TYPE,     "Stellar_Type(2)<CE",         "-",                4, 1 }},
    { BINARY_PROPERTY::STELLAR_TYPE_NAME_1_POST_COMMON_ENVELOPE,            { TYPENAME::STRING,           "Stellar_Type(1)>CE",         "-",               42, 1 }},
    { BINARY_PROPERTY::STELLAR_TYPE_NAME_1_PRE_COMMON_ENVELOPE,             { TYPENAME::STRING,           "Stellar_Type(1)<CE",         "-",               42, 1 }},
    { BINARY_PROPERTY::STELLAR_TYPE_NAME_2_POST_COMMON_ENVELOPE,            { TYPENAME::STRING,           "Stellar_Type(2)>CE",         "-",               42, 1 }},
    { BINARY_PROPERTY::STELLAR_TYPE_NAME_2_PRE_COMMON_ENVELOPE,             { TYPENAME::STRING,           "Stellar_Type(2)<CE",         "-",               42, 1 }},
    { BINARY_PROPERTY::SUPERNOVA_ORBIT_INCLINATION_ANGLE,                   { TYPENAME::DOUBLE,           "SN_Orbit_Inclination_Angle", "-",               24, 15}},   
    { BINARY_PROPERTY::SUPERNOVA_ORBIT_INCLINATION_VECTOR_X,                { TYPENAME::DOUBLE,           "Orbital_AM_Vector>SN_X",     "-",               24, 15}},   
    { BINARY_PROPERTY::SUPERNOVA_ORBIT_INCLINATION_VECTOR_Y,                { TYPENAME::DOUBLE,           "Orbital_AM_Vector>SN_Y",     "-",               24, 15}},   
    { BINARY_PROPERTY::SUPERNOVA_ORBIT_INCLINATION_VECTOR_Z,                { TYPENAME::DOUBLE,           "Orbital_AM_Vector>SN_Z",     "-",               24, 15}},   
    { BINARY_PROPERTY::SUPERNOVA_STATE,                                     { TYPENAME::SN_STATE,         "Supernova_State",            "State",            4, 1 }},
    { BINARY_PROPERTY::SYNCHRONIZATION_TIMESCALE,                           { TYPENAME::DOUBLE,           "Tau_Sync",                   "Myr",             24, 15}},
    { BINARY_PROPERTY::SYSTEMIC_SPEED,                                      { TYPENAME::DOUBLE,           "SystemicSpeed",              "kms^-1",          24, 15}},
    { BINARY_PROPERTY::SYSTEMIC_VELOCITY_X,                                 { TYPENAME::DOUBLE,           "SystemicVelocityX",          "kms^-1",          24, 15}},
    { BINARY_PROPERTY::SYSTEMIC_VELOCITY_Y,                                 { TYPENAME::DOUBLE,           "SystemicVelocityY",          "kms^-1",          24, 15}},
    { BINARY_PROPERTY::SYSTEMIC_VELOCITY_Z,                                 { TYPENAME::DOUBLE,           "SystemicVelocityZ",          "kms^-1",          24, 15}},
    { BINARY_PROPERTY::TIME,                                                { TYPENAME::DOUBLE,           "Time",                       "Myr",             24, 15}},
    { BINARY_PROPERTY::TIME_TO_COALESCENCE,                                 { TYPENAME::DOUBLE,           "Coalescence_Time",           "Myr",             24, 15}},
    { BINARY_PROPERTY::TOTAL_ANGULAR_MOMENTUM,                              { TYPENAME::DOUBLE,           "Ang_Momentum_Total",         "Msol AU^2 yr^-1", 24, 15}},
    { BINARY_PROPERTY::TOTAL_ENERGY,                                        { TYPENAME::DOUBLE,           "Energy_Total",               "Msol AU^2 yr^-2", 24, 15}},
    { BINARY_PROPERTY::UNBOUND,                                             { TYPENAME::BOOL,             "Unbound",                    "State",            0, 0 }},
    { BINARY_PROPERTY::ZETA_LOBE,                                           { TYPENAME::DOUBLE,           "Zeta_Lobe",                  "-",               24, 15}},
    { BINARY_PROPERTY::ZETA_STAR,                                           { TYPENAME::DOUBLE,           "Zeta_Star",                  "-",               24, 15}}
};

// map PROGRAM_OPTION_DETAIL
// Records the details of PROGRAM_OPTION properties.
//
// Options only need to be here if they are required to be available for printing in 
// the logfiles - all keys present here should also be in PROGRAM_OPTION and PROGRAM_OPTION_LABEL
//
// Note that header strings here should be prefixed with "PO_" to differentiate them from stellar/binary
// properties of the same name
const std::map<PROGRAM_OPTION, PROPERTY_DETAILS> PROGRAM_OPTION_DETAIL = {

    { PROGRAM_OPTION::ADD_OPTIONS_TO_SYSPARMS,                                  { TYPENAME::INT,        "PO_Add_Options_To_SysParms",                "-",          4, 1 }},
    { PROGRAM_OPTION::ALLOW_NON_STRIPPED_ECSN,                                  { TYPENAME::BOOL,       "PO_Allow_Non_Stripped_ECSN",                "Flag",       0, 0 }},
    { PROGRAM_OPTION::ALLOW_MS_STAR_TO_SURVIVE_COMMON_ENVELOPE,                 { TYPENAME::BOOL,       "PO_Allow_MS_To_Survive_CE",                 "Flag",       0, 0 }},
    { PROGRAM_OPTION::ALLOW_RADIATIVE_ENVELOPE_STAR_TO_SURVIVE_COMMON_ENVELOPE, { TYPENAME::BOOL,       "PO_Allow_Radiative_Envelope_To_Survive_CE", "Flag",       0, 0 }},
    { PROGRAM_OPTION::ALLOW_IMMEDIATE_RLOF_POST_CE_TO_SURVIVE_COMMON_ENVELOPE,  { TYPENAME::BOOL,       "PO_Allow_Immediate_RLOF>CE_To_Survive_CE",  "Flag",       0, 0 }},
    { PROGRAM_OPTION::ALLOW_RLOF_AT_BIRTH,                                      { TYPENAME::BOOL,       "PO_Allow_RLOF@Birth",                       "Flag",       0, 0 }},
    { PROGRAM_OPTION::ALLOW_TOUCHING_AT_BIRTH,                                  { TYPENAME::BOOL,       "PO_Allow_Touching@Birth",                   "Flag",       0, 0 }},
    { PROGRAM_OPTION::ANG_MOM_CONSERVATION_DURING_CIRCULARISATION,              { TYPENAME::BOOL,       "PO_Conserve_AngMom@Circ",                   "Flag",       0, 0 }},

    { PROGRAM_OPTION::BLACK_HOLE_KICKS_MODE,                                    { TYPENAME::INT,        "PO_BH_Kicks_Mode",                          "-",          4, 1 }},
    
    { PROGRAM_OPTION::CASE_BB_STABILITY_PRESCRIPTION,                           { TYPENAME::INT,        "PO_BB_Mass_xFer_Stblty_Prscrptn",           "-",          4, 1 }},
    
    { PROGRAM_OPTION::CHECK_PHOTON_TIRING_LIMIT,                                { TYPENAME::BOOL,       "PO_Check_Photon_Tiring_Limit",              "Flag",       0, 0 }},

    { PROGRAM_OPTION::CHE_MODE,                                                 { TYPENAME::INT,        "PO_CHE_Mode",                               "-",          4, 1 }},
      
    { PROGRAM_OPTION::CIRCULARISE_BINARY_DURING_MT,                             { TYPENAME::BOOL,       "PO_Circularise@MT",                         "Flag",       0, 0 }},

    { PROGRAM_OPTION::COMMON_ENVELOPE_ALPHA,                                    { TYPENAME::DOUBLE,     "PO_CE_Alpha",                               "-",         24, 15}},
    { PROGRAM_OPTION::COMMON_ENVELOPE_ALPHA_THERMAL,                            { TYPENAME::DOUBLE,     "PO_CE_Alpha_Thermal",                       "-",         24, 15}},
    { PROGRAM_OPTION::COMMON_ENVELOPE_FORMALISM,                                { TYPENAME::INT,        "PO_CE_Formalism",                           "-",          4, 1 }},
    { PROGRAM_OPTION::COMMON_ENVELOPE_LAMBDA,                                   { TYPENAME::DOUBLE,     "PO_CE_Lambda",                              "-",         24, 15}},
    { PROGRAM_OPTION::COMMON_ENVELOPE_LAMBDA_MULTIPLIER,                        { TYPENAME::DOUBLE,     "PO_CE_Lambda_Multiplier",                   "-",         24, 15}},
    { PROGRAM_OPTION::COMMON_ENVELOPE_LAMBDA_PRESCRIPTION,                      { TYPENAME::INT,        "PO_CE_Lambda_Prscrptn",                     "-",          4, 1 }},
    { PROGRAM_OPTION::COMMON_ENVELOPE_MASS_ACCRETION_CONSTANT,                  { TYPENAME::DOUBLE,     "PO_CE_Mass_Accr_Constant",                  "-",         24, 15}},
    { PROGRAM_OPTION::COMMON_ENVELOPE_MASS_ACCRETION_MAX,                       { TYPENAME::DOUBLE,     "PO_CE_Mass_Accr_Max",                       "Msol",      24, 15}},
    { PROGRAM_OPTION::COMMON_ENVELOPE_MASS_ACCRETION_MIN,                       { TYPENAME::DOUBLE,     "PO_CE_Mass_Accr_Min",                       "Msol",      24, 15}},
    { PROGRAM_OPTION::COMMON_ENVELOPE_MASS_ACCRETION_PRESCRIPTION,              { TYPENAME::INT,        "PO_CE_Mass_Accr_Prscrptn",                  "-",          4, 1 }},
    { PROGRAM_OPTION::COMMON_ENVELOPE_RECOMBINATION_ENERGY_DENSITY,             { TYPENAME::DOUBLE,     "PO_CE_Recomb_Enrgy_Dnsty",                  "erg g^-1",  24, 15}},
    { PROGRAM_OPTION::COMMON_ENVELOPE_SLOPE_KRUCKOW,                            { TYPENAME::DOUBLE,     "PO_CE_Slope_Kruckow",                       "-",         24, 15}},

    { PROGRAM_OPTION::COOL_WIND_MASS_LOSS_MULTIPLIER,                           { TYPENAME::DOUBLE,     "PO_Cool_WindMassLoss_Multipl",              "-",         24, 15}},

    { PROGRAM_OPTION::ECCENTRICITY,                                             { TYPENAME::DOUBLE,     "PO_Eccentricity",                           "-",         24, 15}},
    { PROGRAM_OPTION::ECCENTRICITY_DISTRIBUTION,                                { TYPENAME::INT,        "PO_Eccentricity_Dstrbtn",                   "-",          4, 1 }},
    { PROGRAM_OPTION::ECCENTRICITY_DISTRIBUTION_MAX,                            { TYPENAME::DOUBLE,     "PO_Eccentricity_Dstrbtn_Max",               "-",         24, 15}},
    { PROGRAM_OPTION::ECCENTRICITY_DISTRIBUTION_MIN,                            { TYPENAME::DOUBLE,     "PO_Eccentricity_Dstrbtn_Min",               "-",         24, 15}},
    { PROGRAM_OPTION::EDDINGTON_ACCRETION_FACTOR,                               { TYPENAME::DOUBLE,     "PO_Eddington_Accr_Factor",                  "-",         24, 15}},
    { PROGRAM_OPTION::ENABLE_ROTATIONALLY_ENHANCED_MASS_LOSS,                   { TYPENAME::BOOL,       "PO_Enable_Rotationally_Enhanced_Mass_Loss", "Flag",       0,  0}},
    { PROGRAM_OPTION::ENHANCE_CHE_LIFETIMES_LUMINOSITIES,                       { TYPENAME::BOOL,       "PO_Enhance_CHE_lifetimes_luminosities",     "Flag",       0,  0}},
    { PROGRAM_OPTION::ENVELOPE_STATE_PRESCRIPTION,                              { TYPENAME::INT,        "PO_Envelope_State_Prscrptn",                "-",          4, 1 }},
    { PROGRAM_OPTION::EVOLUTION_MODE,                                           { TYPENAME::INT,        "PO_Evolution_Mode",                         "Mode",       4, 1 }},

    { PROGRAM_OPTION::FP_ERROR_MODE,                                            { TYPENAME::INT,        "PO_FP_Error_Mode",                          "Mode",       4, 1 }},

    { PROGRAM_OPTION::FRYER_SUPERNOVA_ENGINE,                                   { TYPENAME::INT,        "PO_Fryer_SN_Engine",                        "-",          4, 1 }},

    { PROGRAM_OPTION::FRYER22_FMIX,                                             { TYPENAME::DOUBLE,     "PO_Fryer22_mixing_fraction",                "-",         24, 15}},
    { PROGRAM_OPTION::FRYER22_MCRIT,                                            { TYPENAME::DOUBLE,     "PO_Fryer22_crit_COcore_Mass",               "Msol",      24, 15}},

    { PROGRAM_OPTION::INITIAL_MASS,                                             { TYPENAME::DOUBLE,     "PO_Initial_Mass",                           "Msol",      24, 15}},
    { PROGRAM_OPTION::INITIAL_MASS_1,                                           { TYPENAME::DOUBLE,     "PO_Initial_Mass(1)",                        "Msol",      24, 15}},
    { PROGRAM_OPTION::INITIAL_MASS_2,                                           { TYPENAME::DOUBLE,     "PO_Initial_Mass(2)",                        "Msol",      24, 15}},

    { PROGRAM_OPTION::INITIAL_MASS_FUNCTION,                                    { TYPENAME::INT,        "PO_Initial_Mass_Function",                  "-",          4, 1 }},
    { PROGRAM_OPTION::INITIAL_MASS_FUNCTION_MAX,                                { TYPENAME::DOUBLE,     "PO_Initial_Mass_Func_Max",                  "Msol",      24, 15}},
    { PROGRAM_OPTION::INITIAL_MASS_FUNCTION_MIN,                                { TYPENAME::DOUBLE,     "PO_Initial_Mass_Func_Min",                  "Msol",      24, 15}},
    { PROGRAM_OPTION::INITIAL_MASS_FUNCTIONPOWER,                               { TYPENAME::DOUBLE,     "PO_Initial_Mass_Func_Power",                "-",         24, 15}},

    { PROGRAM_OPTION::KICK_DIRECTION_DISTRIBUTION,                              { TYPENAME::INT,        "PO_Kick_Direction_Dstrbtn",                 "-",          4, 1 }},
    { PROGRAM_OPTION::KICK_DIRECTION_POWER,                                     { TYPENAME::DOUBLE,     "PO_Kick_Direction_Power",                   "-",         24, 15}},
    { PROGRAM_OPTION::KICK_SCALING_FACTOR,                                      { TYPENAME::DOUBLE,     "PO_Kick_Scaling_Factor",                    "-",         24, 15}},
    { PROGRAM_OPTION::KICK_MAGNITUDE_DISTRIBUTION,                              { TYPENAME::INT,        "PO_Kick_Magnitude_Dstrbtn",                 "-",          4, 1 }},
    { PROGRAM_OPTION::KICK_MAGNITUDE_DISTRIBUTION_MAXIMUM,                      { TYPENAME::DOUBLE,     "PO_Kick_Magnitude_Dstrbtn_Max",             "-",         24, 15}},

    { PROGRAM_OPTION::KICK_MAGNITUDE_DISTRIBUTION_SIGMA_CCSN_BH,                { TYPENAME::DOUBLE,     "PO_Sigma_Kick_CCSN_BH",                     "kms^-1",    24, 15}},
    { PROGRAM_OPTION::KICK_MAGNITUDE_DISTRIBUTION_SIGMA_CCSN_NS,                { TYPENAME::DOUBLE,     "PO_Sigma_Kick_CCSN_NS",                     "kms^-1",    24, 15}},
    { PROGRAM_OPTION::KICK_MAGNITUDE_DISTRIBUTION_SIGMA_FOR_ECSN,               { TYPENAME::DOUBLE,     "PO_Sigma_Kick_ECSN",                        "kms^-1",    24, 15}},
    { PROGRAM_OPTION::KICK_MAGNITUDE_DISTRIBUTION_SIGMA_FOR_USSN,               { TYPENAME::DOUBLE,     "PO_Sigma_Kick_USSN",                        "kms^-1",    24, 15}},

    { PROGRAM_OPTION::KICK_MAGNITUDE,                                           { TYPENAME::DOUBLE,     "PO_Kick_Magnitude",                         "kms^-1",    24, 15}},
    { PROGRAM_OPTION::KICK_MAGNITUDE_1,                                         { TYPENAME::DOUBLE,     "PO_Kick_Magnitude(1)",                      "kms^-1",    24, 15}},
    { PROGRAM_OPTION::KICK_MAGNITUDE_2,                                         { TYPENAME::DOUBLE,     "PO_Kick_Magnitude(2)",                      "kms^-1",    24, 15}},

    { PROGRAM_OPTION::KICK_MAGNITUDE_RANDOM,                                    { TYPENAME::DOUBLE,     "PO_Kick_Magnitude_Random",                  "-",         24, 15}},
    { PROGRAM_OPTION::KICK_MAGNITUDE_RANDOM_1,                                  { TYPENAME::DOUBLE,     "PO_Kick_Magnitude_Random(1)",               "-",         24, 15}},
    { PROGRAM_OPTION::KICK_MAGNITUDE_RANDOM_2,                                  { TYPENAME::DOUBLE,     "PO_Kick_Magnitude_Random(2)",               "-",         24, 15}},

    { PROGRAM_OPTION::KICK_MEAN_ANOMALY_1,                                      { TYPENAME::DOUBLE,     "PO_Kick_Mean_Anomaly(1)",                   "-",         24, 15}},
    { PROGRAM_OPTION::KICK_MEAN_ANOMALY_2,                                      { TYPENAME::DOUBLE,     "PO_Kick_Mean_Anomaly(2)",                   "-",         24, 15}},
    { PROGRAM_OPTION::KICK_PHI_1,                                               { TYPENAME::DOUBLE,     "PO_Kick_Phi(1)",                            "-",         24, 15}},
    { PROGRAM_OPTION::KICK_PHI_2,                                               { TYPENAME::DOUBLE,     "PO_Kick_Phi(2)",                            "-",         24, 15}},
    { PROGRAM_OPTION::KICK_THETA_1,                                             { TYPENAME::DOUBLE,     "PO_Kick_Theta(1)",                          "-",         24, 15}},
    { PROGRAM_OPTION::KICK_THETA_2,                                             { TYPENAME::DOUBLE,     "PO_Kick_Theta(2)",                          "-",         24, 15}},

    { PROGRAM_OPTION::LBV_FACTOR,                                               { TYPENAME::DOUBLE,     "PO_LBV_Factor",                             "-",         24, 15}},
    { PROGRAM_OPTION::LBV_MASS_LOSS_PRESCRIPTION,                               { TYPENAME::INT,        "PO_LBV_Mass_Loss_Prscrptn",                 "-",          4, 1 }},

    { PROGRAM_OPTION::MASS_LOSS_PRESCRIPTION,                                   { TYPENAME::INT,        "PO_Mass_Loss_Prscrptn",                     "-",          4, 1 }},

    { PROGRAM_OPTION::MASS_RATIO,                                               { TYPENAME::DOUBLE,     "PO_Mass_Ratio",                             "-",         24, 15}},
    { PROGRAM_OPTION::MASS_RATIO_DISTRIBUTION,                                  { TYPENAME::INT,        "PO_Mass_Ratio_Dstrbtn",                     "-",          4, 1 }},
    { PROGRAM_OPTION::MASS_RATIO_DISTRIBUTION_MAX,                              { TYPENAME::DOUBLE,     "PO_Mass_Ratio_Dstrbtn_Max",                 "-",         24, 15}},
    { PROGRAM_OPTION::MASS_RATIO_DISTRIBUTION_MIN,                              { TYPENAME::DOUBLE,     "PO_Mass_Ratio_Dstrbtn_Min",                 "-",         24, 15}},

    { PROGRAM_OPTION::MAXIMUM_EVOLUTION_TIME,                                   { TYPENAME::DOUBLE,     "PO_Max_Evolution_Time",                     "Myr",       24, 15}},
    { PROGRAM_OPTION::MAXIMUM_DONOR_MASS,                                       { TYPENAME::DOUBLE,     "PO_Max_Donor_Mass",                         "Msol",      24, 15}},
    { PROGRAM_OPTION::MAXIMUM_NEUTRON_STAR_MASS,                                { TYPENAME::DOUBLE,     "PO_Max_NS_Mass",                            "Msol",      24, 15}},
    { PROGRAM_OPTION::MAXIMUM_TIMESTEPS,                                        { TYPENAME::ULONGINT,   "PO_Max_Timesteps",                          "Count",     10, 1 }},

    { PROGRAM_OPTION::MCBUR1,                                                   { TYPENAME::DOUBLE,     "PO_MCBUR1",                                 "Msol",      24, 15}},

    { PROGRAM_OPTION::METALLICITY,                                              { TYPENAME::DOUBLE,     "PO_Metallicity",                            "-",         24, 15}},
    { PROGRAM_OPTION::METALLICITY_DISTRIBUTION,                                 { TYPENAME::INT,        "PO_Metallicity_Dstrbtn",                    "-",          4, 1 }},
    { PROGRAM_OPTION::METALLICITY_DISTRIBUTION_MAX,                             { TYPENAME::DOUBLE,     "PO_Metallicity_Dstrbtn_Max",                "-",         24, 15}},
    { PROGRAM_OPTION::METALLICITY_DISTRIBUTION_MIN,                             { TYPENAME::DOUBLE,     "PO_Metallicity_Dstrbtn_Min",                "-",         24, 15}},

    { PROGRAM_OPTION::MINIMUM_MASS_SECONDARY,                                   { TYPENAME::DOUBLE,     "PO_Min_Secondary_Mass",                     "Msol",      24, 15}},

    { PROGRAM_OPTION::MT_ACCRETION_EFFICIENCY_PRESCRIPTION,                     { TYPENAME::INT,        "PO_MT_Acc_Efficiency_Prscrptn",             "-",          4, 1 }},
    { PROGRAM_OPTION::MT_ANG_MOM_LOSS_PRESCRIPTION,                             { TYPENAME::INT,        "PO_MT_AngMom_Loss_Prscrptn",                "-",          4, 1 }},
    { PROGRAM_OPTION::MT_THERMAL_LIMIT_C,                                       { TYPENAME::DOUBLE,     "PO_MT_Thermal_Limit_C",                     "-",         24, 15}},

    { PROGRAM_OPTION::MT_CRIT_MR_MS_LOW_MASS_DEGENERATE_ACCRETOR,               { TYPENAME::DOUBLE,     "PO_MT_Crit_MR_MS_Low_Mass_Deg_Acc",         "-",         24, 15}},
    { PROGRAM_OPTION::MT_CRIT_MR_MS_LOW_MASS_NON_DEGENERATE_ACCRETOR,           { TYPENAME::DOUBLE,     "PO_MT_Crit_MR_MS_Low_Mass_NonDeg_Acc",      "-",         24, 15}},
    { PROGRAM_OPTION::MT_CRIT_MR_MS_HIGH_MASS_DEGENERATE_ACCRETOR,              { TYPENAME::DOUBLE,     "PO_MT_Crit_MR_MS_High_Mass_Deg_Acc",        "-",         24, 15}},
    { PROGRAM_OPTION::MT_CRIT_MR_MS_HIGH_MASS_NON_DEGENERATE_ACCRETOR,          { TYPENAME::DOUBLE,     "PO_MT_Crit_MR_MS_High_Mass_NonDeg_Acc",     "-",         24, 15}},
    { PROGRAM_OPTION::MT_CRIT_MR_GIANT_DEGENERATE_ACCRETOR,                     { TYPENAME::DOUBLE,     "PO_MT_Crit_MR_Giant_Deg_Acc",               "-",         24, 15}},
    { PROGRAM_OPTION::MT_CRIT_MR_GIANT_NON_DEGENERATE_ACCRETOR,                 { TYPENAME::DOUBLE,     "PO_MT_Crit_MR_Giant_NonDeg_Acc",            "-",         24, 15}},
    { PROGRAM_OPTION::MT_CRIT_MR_HG_DEGENERATE_ACCRETOR,                        { TYPENAME::DOUBLE,     "PO_MT_Crit_MR_HG_Deg_Acc",                  "-",         24, 15}},

    { PROGRAM_OPTION::MT_CRIT_MR_HG_NON_DEGENERATE_ACCRETOR,                    { TYPENAME::DOUBLE,     "PO_MT_Crit_MR_HG_NonDeg_Acc",               "-",         24, 15}},
    { PROGRAM_OPTION::MT_CRIT_MR_HE_GIANT_DEGENERATE_ACCRETOR,                  { TYPENAME::DOUBLE,     "PO_MT_Crit_MR_HE_Giant_Deg_Acc",            "-",         24, 15}},
    { PROGRAM_OPTION::MT_CRIT_MR_HE_GIANT_NON_DEGENERATE_ACCRETOR,              { TYPENAME::DOUBLE,     "PO_MT_Crit_MR_HE_Giant_NonDeg_Acc",         "-",         24, 15}},
    { PROGRAM_OPTION::MT_CRIT_MR_HE_HG_DEGENERATE_ACCRETOR,                     { TYPENAME::DOUBLE,     "PO_MT_Crit_MR_HE_HG_Deg_Acc",               "-",         24, 15}},
    { PROGRAM_OPTION::MT_CRIT_MR_HE_HG_NON_DEGENERATE_ACCRETOR,                 { TYPENAME::DOUBLE,     "PO_MT_Crit_MR_HE_HG_NonDeg_Acc",            "-",         24, 15}},
    { PROGRAM_OPTION::MT_CRIT_MR_HE_MS_DEGENERATE_ACCRETOR,                     { TYPENAME::DOUBLE,     "PO_MT_Crit_MR_HE_MS_Deg_Acc",               "-",         24, 15}},
    { PROGRAM_OPTION::MT_CRIT_MR_HE_MS_NON_DEGENERATE_ACCRETOR,                 { TYPENAME::DOUBLE,     "PO_MT_Crit_MR_HE_MS_NonDeg_Acc",            "-",         24, 15}},
    { PROGRAM_OPTION::MT_CRIT_MR_WD_DEGENERATE_ACCRETOR,                        { TYPENAME::DOUBLE,     "PO_MT_Crit_MR_WD_Deg_Acc",                  "-",         24, 15}},
    { PROGRAM_OPTION::MT_CRIT_MR_WD_NONDEGENERATE_ACCRETOR,                     { TYPENAME::DOUBLE,     "PO_MT_Crit_MR_WD_NonDeg_Acc",               "-",         24, 15}},
    
    { PROGRAM_OPTION::MT_FRACTION_ACCRETED,                                     { TYPENAME::DOUBLE,     "PO_MT_Fraction_Accreted",                   "-",         24, 15}},
    { PROGRAM_OPTION::MT_JLOSS,                                                 { TYPENAME::DOUBLE,     "PO_MT_JLoss",                               "-",         24, 15}},
    { PROGRAM_OPTION::MT_JLOSS_MACLEOD_LINEAR_FRACTION_DEGEN,                   { TYPENAME::DOUBLE,     "PO_MT_JLoss_Macleod_Linear_Frac_Degen",     "-",         24, 15}},
    { PROGRAM_OPTION::MT_JLOSS_MACLEOD_LINEAR_FRACTION_NON_DEGEN,               { TYPENAME::DOUBLE,     "PO_MT_JLoss_Macleod_Linear_Frac_Non_Degen", "-",         24, 15}},
    { PROGRAM_OPTION::MT_REJUVENATION_PRESCRIPTION,                             { TYPENAME::INT,        "PO_MT_Rejuvenation_Prscrptn",               "-",          4, 1 }},
    { PROGRAM_OPTION::MT_THERMALLY_LIMITED_VARIATION,                           { TYPENAME::INT,        "PO_MT_Thermally_Lmtd_Variation",            "-",          4, 1 }},

    { PROGRAM_OPTION::MULLER_MANDEL_KICK_MULTIPLIER_BH,                         { TYPENAME::DOUBLE,     "PO_MM_Kick_Multiplier_BH",                  "-",         24, 15}},
    { PROGRAM_OPTION::MULLER_MANDEL_KICK_MULTIPLIER_NS,                         { TYPENAME::DOUBLE,     "PO_MM_Kick_Multiplier_NS",                  "-",         24, 15}},
    { PROGRAM_OPTION::MULLER_MANDEL_SIGMA_KICK,                                 { TYPENAME::DOUBLE,     "PO_MM_Sigma_Kick",                          "-",         24, 15}},
    
    { PROGRAM_OPTION::NEUTRINO_MASS_LOSS_ASSUMPTION_BH,                         { TYPENAME::INT,        "PO_Neutrino_Mass_Loss_Assmptn",             "-",          4, 1 }},
    { PROGRAM_OPTION::NEUTRINO_MASS_LOSS_VALUE_BH,                              { TYPENAME::DOUBLE,     "PO_Neutrino_Mass_Loss_Value",               "-",         24, 15}},

    { PROGRAM_OPTION::NS_EOS,                                                   { TYPENAME::INT,        "PO_NS_EOS",                                 "-",          4, 1 }},

    { PROGRAM_OPTION::ORBITAL_PERIOD,                                           { TYPENAME::DOUBLE,     "PO_Orbital_Period",                         "days",      24, 15}},
    { PROGRAM_OPTION::ORBITAL_PERIOD_DISTRIBUTION,                              { TYPENAME::INT,        "PO_Orbital_Period_Dstrbtn",                 "-",          4, 1 }},
    { PROGRAM_OPTION::ORBITAL_PERIOD_DISTRIBUTION_MAX,                          { TYPENAME::DOUBLE,     "PO_Orbital_Period_Max",                     "days",      24, 15}},
    { PROGRAM_OPTION::ORBITAL_PERIOD_DISTRIBUTION_MIN,                          { TYPENAME::DOUBLE,     "PO_Orbital_Period_Min",                     "days",      24, 15}},

    { PROGRAM_OPTION::OVERALL_WIND_MASS_LOSS_MULTIPLIER,                        { TYPENAME::DOUBLE,     "PO_Overall_WindMassLoss_Multipl",           "-",         24, 15}},

    { PROGRAM_OPTION::PISN_LOWER_LIMIT,                                         { TYPENAME::DOUBLE,     "PO_PISN_Lower_Limit",                       "Msol",      24, 15}},
    { PROGRAM_OPTION::PISN_UPPER_LIMIT,                                         { TYPENAME::DOUBLE,     "PO_PISN_Upper_Limit",                       "Msol",      24, 15}},

    { PROGRAM_OPTION::PPI_LOWER_LIMIT,                                          { TYPENAME::DOUBLE,     "PO_PPI_Lower_Limit",                        "Msol",      24, 15}},
    { PROGRAM_OPTION::PPI_PRESCRIPTION,                                         { TYPENAME::INT,        "PO_PPI_Prscrptn",                           "-",          4, 1 }},
    { PROGRAM_OPTION::PPI_UPPER_LIMIT,                                          { TYPENAME::DOUBLE,     "PO_PPI_Upper_Limit",                        "Msol",      24, 15}},
    { PROGRAM_OPTION::PPI_CO_CORE_SHIFT_HENDRIKS,                               { TYPENAME::DOUBLE,     "PO_PPI_CO_CORE_SHIFT_HENDRIKS",             "Msol",      24, 15}},

    { PROGRAM_OPTION::PULSAR_MAGNETIC_FIELD_DISTRIBUTION,                       { TYPENAME::INT,        "PO_Pulsar_Mag_Field_Dstrbtn",               "-",          4, 1 }},
    { PROGRAM_OPTION::PULSAR_MAGNETIC_FIELD_DISTRIBUTION_MAX,                   { TYPENAME::DOUBLE,     "PO_Pulsar_Mag_Field_Dstrbtn_Max",           "AU",        24, 15}},
    { PROGRAM_OPTION::PULSAR_MAGNETIC_FIELD_DISTRIBUTION_MIN,                   { TYPENAME::DOUBLE,     "PO_Pulsar_Mag_Field_Dstrbtn_Min",           "AU",        24, 15}},

    { PROGRAM_OPTION::PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION,                    { TYPENAME::INT,        "PO_Pulsar_Spin_Period_Dstrbtn",             "-",          4, 1 }},
    { PROGRAM_OPTION::PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION_MAX,                { TYPENAME::DOUBLE,     "PO_Pulsar_Spin_Period_Dstrbtn_Max",         "AU",        24, 15}},
    { PROGRAM_OPTION::PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION_MIN,                { TYPENAME::DOUBLE,     "PO_Pulsar_Spin_Period_Dstrbtn_Min",         "AU",        24, 15}},

    { PROGRAM_OPTION::PULSAR_MAGNETIC_FIELD_DECAY_MASS_SCALE,                   { TYPENAME::DOUBLE,     "PO_Pulsar_Mag_Field_Decay_mScale",          "Msol",      24, 15}},
    { PROGRAM_OPTION::PULSAR_MAGNETIC_FIELD_DECAY_TIME_SCALE,                   { TYPENAME::DOUBLE,     "PO_Pulsar_Mag_Field_Decay_tScale",          "Myr",       24, 15}},

    { PROGRAM_OPTION::PULSAR_MINIMUM_MAGNETIC_FIELD,                            { TYPENAME::DOUBLE,     "PO_Pulsar_Minimum_Mag_Field",               "Gauss",     24, 15}},

    { PROGRAM_OPTION::QCRIT_PRESCRIPTION,                                       { TYPENAME::INT,        "PO_qCrit_Prescription",                     "-",          4, 1 }},

    { PROGRAM_OPTION::RANDOM_SEED,                                              { TYPENAME::ULONGINT,   "PO_SEED(OPTION)",                           "-",         12, 1 }},
    { PROGRAM_OPTION::RANDOM_SEED_CMDLINE,                                      { TYPENAME::ULONGINT,   "PO_SEED(CMDLINE)",                          "-",         12, 1 }},

    { PROGRAM_OPTION::REMNANT_MASS_PRESCRIPTION,                                { TYPENAME::INT,        "PO_Remnant_Mass_Prscrptn",                  "-",          4, 1 }},

    { PROGRAM_OPTION::ROCKET_KICK_MAGNITUDE_1,                                  { TYPENAME::DOUBLE,     "PO_Rocket_Kick_Magnitude(1)",               "kms^-1",    24, 15}},
    { PROGRAM_OPTION::ROCKET_KICK_MAGNITUDE_2,                                  { TYPENAME::DOUBLE,     "PO_Rocket_Kick_Magnitude(2)",               "kms^-1",    24, 15}},
    { PROGRAM_OPTION::ROCKET_KICK_PHI_1,                                        { TYPENAME::DOUBLE,     "PO_Rocket_Kick_Phi(1)",                     "-",         24, 15}},
    { PROGRAM_OPTION::ROCKET_KICK_PHI_2,                                        { TYPENAME::DOUBLE,     "PO_Rocket_Kick_Phi(2)",                     "-",         24, 15}},
    { PROGRAM_OPTION::ROCKET_KICK_THETA_1,                                      { TYPENAME::DOUBLE,     "PO_Rocket_Kick_Theta(1)",                   "-",         24, 15}},
    { PROGRAM_OPTION::ROCKET_KICK_THETA_2,                                      { TYPENAME::DOUBLE,     "PO_Rocket_Kick_Theta(2)",                   "-",         24, 15}},

    { PROGRAM_OPTION::ROTATIONAL_VELOCITY_DISTRIBUTION,                         { TYPENAME::INT,        "PO_Rotational_Velocity_Dstrbtn",            "-",          4, 1 }},
    { PROGRAM_OPTION::ROTATIONAL_FREQUENCY,                                     { TYPENAME::DOUBLE,     "PO_Rotational_Frequency",                   "Hz",        24, 15}},
    { PROGRAM_OPTION::ROTATIONAL_FREQUENCY_1,                                   { TYPENAME::DOUBLE,     "PO_Rotational_Frequency(1)",                "Hz",        24, 15}},
    { PROGRAM_OPTION::ROTATIONAL_FREQUENCY_2,                                   { TYPENAME::DOUBLE,     "PO_Rotational_Frequency(2)",                "Hz",        24, 15}},
   
    { PROGRAM_OPTION::SCALE_CHE_MASS_LOSS_SURF_HE_ABUNDANCE,                    { TYPENAME::BOOL,       "PO_Scale_CHE_Mass_Loss_Surf_He_Abundance",  "flag",       0,  0}},
    { PROGRAM_OPTION::SCALE_TERMINAL_WIND_VEL_METALLICITY_POWER,                { TYPENAME::DOUBLE,     "PO_Scale_Terminal_Wind_Vel_Metallicity_Power", "-",      24, 15}},
    { PROGRAM_OPTION::SEMI_MAJOR_AXIS,                                          { TYPENAME::DOUBLE,     "PO_Semi-Major_Axis",                        "AU",        24, 15}},
    { PROGRAM_OPTION::SEMI_MAJOR_AXIS_DISTRIBUTION,                             { TYPENAME::INT,        "PO_Semi-Major_Axis_Dstrbtn",                "-",          4, 1 }},
    { PROGRAM_OPTION::SEMI_MAJOR_AXIS_DISTRIBUTION_MAX,                         { TYPENAME::DOUBLE,     "PO_Semi-Major_Axis_Dstrbtn_Max",            "AU",        24, 15}},
    { PROGRAM_OPTION::SEMI_MAJOR_AXIS_DISTRIBUTION_MIN,                         { TYPENAME::DOUBLE,     "PO_Semi-Major_Axis_Dstrbtn_Min",            "AU",        24, 15}},

    { PROGRAM_OPTION::STELLAR_ZETA_PRESCRIPTION,                                { TYPENAME::INT,        "PO_Stellar_Zeta_Prscrptn",                  "-",          4, 1 }},

    { PROGRAM_OPTION::TIDES_PRESCRIPTION,                                       { TYPENAME::INT,        "PO_Tides_Prscrptn",                         "-",          4, 1 }},

    { PROGRAM_OPTION::WR_FACTOR,                                                { TYPENAME::DOUBLE,     "PO_WR_Factor",                              "-",         24, 15}},

    { PROGRAM_OPTION::ZETA_ADIABATIC_ARBITRARY,                                 { TYPENAME::DOUBLE,     "PO_Zeta_Adiabatic_Arbitrary",               "-",         24, 15}},
    { PROGRAM_OPTION::ZETA_MS,                                                  { TYPENAME::DOUBLE,     "PO_Zeta_Main_Sequence",                     "-",         24, 15}},
    { PROGRAM_OPTION::ZETA_RADIATIVE_ENVELOPE_GIANT,                            { TYPENAME::DOUBLE,     "PO_Zeta_Radiative_Envelope_Giant",          "-",         24, 15}}
};


// enum class LOGFILE
// Symbolic names for logfiles
enum class LOGFILE: int {
    NONE,

    DEBUG_LOG,
    ERROR_LOG,

    BSE_COMMON_ENVELOPES,
    BSE_DETAILED_OUTPUT,
    BSE_DOUBLE_COMPACT_OBJECTS,
    BSE_PULSAR_EVOLUTION,
    BSE_RLOF_PARAMETERS,
    BSE_SUPERNOVAE,
    BSE_SWITCH_LOG,
    BSE_SYSTEM_PARAMETERS,

    SSE_DETAILED_OUTPUT,
    SSE_PULSAR_EVOLUTION,
    SSE_SUPERNOVAE,
    SSE_SWITCH_LOG,
    SSE_SYSTEM_PARAMETERS
};


// Logfile record types
// Note all enum classes for log record types start at 1 (and *must* start at 1)
typedef unsigned int LOGRECORDTYPE;

enum class CE_RECORD_TYPE: unsigned int {                                                                           // BSE_COMMON_ENVELOPES file record type
    DEFAULT = 1                                                                                                     // 1 - default BSE_COMMON_ENVELOPES file record type
};

enum class DCO_RECORD_TYPE: unsigned int {                                                                          // BSE_DOUBLE_COMPACT_OBJECTS file record type
    DEFAULT = 1                                                                                                     // 1 - default BSE_DOUBLE_COMPACT_OBJECTS file record type
};

enum class BSE_PULSAR_RECORD_TYPE: unsigned int {                                                                   // BSE_PULSAR_EVOLUTION file record type
    DEFAULT = 1,                                                                                                    // 1 - default BSE_PULSAR_EVOLUTION file record type
    POST_SN,                                                                                                        // 2 - record was logged immediately following a supernova event
    POST_BINARY_TIMESTEP                                                                                            // 3 - record was logged immediately following binary timestep (i.e. the evolution of the binary system for a single timestep)
};

enum class RLOF_RECORD_TYPE: unsigned int {                                                                         // BSE_RLOF_PARAMETERS file record type
    DEFAULT = 1                                                                                                     // 1 - default BSE_RLOF_PARAMETERS file record type
};

enum class BSE_DETAILED_RECORD_TYPE: unsigned int {                                                                 // BSE_DETAILED_OUTPUT file record type
    INITIAL_STATE = 1,                                                                                              //  1 - record describes the initial state of the binary
    POST_STELLAR_TIMESTEP,                                                                                          //  2 - record was logged immediately following stellar timestep (i.e. the evolution of the constituent stars for a single timestep)
    POST_BINARY_TIMESTEP,                                                                                           //  3 - record was logged immediately following binary timestep (i.e. the evolution of the binary system for a single timestep)
    TIMESTEP_COMPLETED,                                                                                             //  4 - record was logged immediately following the completion of the timestep (after all changes to the binary and components)
    FINAL_STATE,                                                                                                    //  5 - record describes the final state of the binary
    STELLAR_TYPE_CHANGE_DURING_CEE,                                                                                 //  6 - record was logged immediately following a stellar type change during a common envelope event
    STELLAR_TYPE_CHANGE_DURING_MT,                                                                                  //  7 - record was logged immediately following a stellar type change during a mass transfer event
    STELLAR_TYPE_CHANGE_DURING_MASS_RESOLUTION,                                                                     //  8 - record was logged immediately following a stellar type change during mass resolution
    STELLAR_TYPE_CHANGE_DURING_CHE_EQUILIBRATION,                                                                   //  9 - record was logged immediately following a stellar type change during mass equilibration for CHE
    POST_MT,                                                                                                        // 10 - record was logged immediately following a mass transfer event
    POST_WINDS,                                                                                                     // 11 - record was logged immediately following winds mass loss
    POST_CEE,                                                                                                       // 12 - record was logged immediately following a common envelope event
    POST_SN,                                                                                                        // 13 - record was logged immediately following a supernova event
    POST_MASS_RESOLUTION,                                                                                           // 14 - record was logged immediately following mass resolution (i.e. after winds mass loss & mass transfer complete)
    POST_MASS_RESOLUTION_MERGER,                                                                                    // 15 - record was logged immediately following a merger after mass resolution
    PRE_STELLAR_TIMESTEP                                                                                            // 16 - record was logged immediately prior to stellar timestep (i.e. the evolution of the constituent stars for a single timestep)
};

enum class SSE_DETAILED_RECORD_TYPE: unsigned int {                                                                 // SSE_DETAILED_OUTPUT file record type
    INITIAL_STATE = 1,                                                                                              //  1 - record describes the initial state of the star
    PRE_MASS_LOSS,                                                                                                  //  2 - record was logged after timestep taken, but before mass loss resolution
    POST_MASS_LOSS,                                                                                                 //  3 - record was logged after after mass loss resolution
    TIMESTEP_COMPLETED,                                                                                             //  4 - record was logged immediately following the completion of the timestep (after all changes to the star)
    FINAL_STATE                                                                                                     //  5 - record describes the final state of the star
};

enum class SSE_PULSAR_RECORD_TYPE: unsigned int {                                                                   // SSE_PULSAR_EVOLUTION file record type
    DEFAULT = 1,                                                                                                    // 1 - default SSE_PULSAR_EVOLUTION file record type
    POST_SN,                                                                                                        // 2 - record was logged immediately following a supernova event
    TIMESTEP_COMPLETED                                                                                              // 3 - record was logged immediately following the completion of the timestep (after all changes to the star)
};

enum class BSE_SN_RECORD_TYPE: unsigned int {                                                                       // BSE_SUPERNOVAE file record type
    DEFAULT = 1                                                                                                     // 1 - default BSE_SUPERNOVAE file record type
};

enum class SSE_SN_RECORD_TYPE: unsigned int {                                                                       // SSE_SUPERNOVAE file record type
    DEFAULT = 1                                                                                                     // 1 - default SSE_SUPERNOVAE file record type
};

enum class BSE_SYSPARMS_RECORD_TYPE: unsigned int {                                                                 // BSE_SYSTEM_PARAMETERS file record type
    DEFAULT = 1                                                                                                     // 1 - default BSE_SYSTEM_PARAMETERS file record type
};

enum class SSE_SYSPARMS_RECORD_TYPE: unsigned int {                                                                 // SSE_SYSTEM_PARAMETERS file record type
    DEFAULT = 1                                                                                                     // 1 - default SSE_SYSTEM_PARAMETERS file record type
};


// enum class RUN_DETAILS_REC
// symbolic names for RUN DETAILS record definitions
// For preamble/stats columns only
// Program options columns are dynamic
// these must be left as default values - their order can be changed with the caveat that the sentinel "SENTINEL" must stay at the end
// it's a bit of a hack, but it lets me iterate over the enum
enum class RUN_DETAILS_COLUMNS: int { COMPAS_VERSION,
                                      RUN_START, 
                                      RUN_END, 
                                      OBJECTS_REQUESTED,
                                      OBJECTS_CREATED,
                                      CLOCK_TIME,
                                      WALL_TIME,
                                      ACTUAL_RANDOM_SEED,
                                      SENTINEL };

const COMPASUnorderedMap<RUN_DETAILS_COLUMNS, std::tuple<std::string, TYPENAME, std::size_t>> RUN_DETAILS_DETAIL = {
    { RUN_DETAILS_COLUMNS::COMPAS_VERSION,      { "COMPAS-Version",                TYPENAME::STRING,    8 }},
    { RUN_DETAILS_COLUMNS::RUN_START,           { "Run-Start",                     TYPENAME::STRING,   24 }},
    { RUN_DETAILS_COLUMNS::RUN_END,             { "Run-End",                       TYPENAME::STRING,   24 }},
    { RUN_DETAILS_COLUMNS::OBJECTS_REQUESTED,   { "Objects-Requested",             TYPENAME::INT,       0 }},
    { RUN_DETAILS_COLUMNS::OBJECTS_CREATED,     { "Objects-Created",               TYPENAME::INT,       0 }},
    { RUN_DETAILS_COLUMNS::CLOCK_TIME,          { "Clock-Time",                    TYPENAME::DOUBLE,    0 }},
    { RUN_DETAILS_COLUMNS::WALL_TIME,           { "Wall-Time",                     TYPENAME::STRING,   10 }},
    { RUN_DETAILS_COLUMNS::ACTUAL_RANDOM_SEED,  { "Actual-Random-Seed",            TYPENAME::ULONGINT,  0 }}
};


// BSE output record definitions

// BSE_COMMON_ENVELOPES_REC
//
// Default record definition for the Common Envelopes logfile
//
const ANY_PROPERTY_VECTOR BSE_COMMON_ENVELOPES_REC = {
    BINARY_PROPERTY::RANDOM_SEED,
    BINARY_PROPERTY::TIME,
    STAR_1_PROPERTY::LAMBDA_AT_COMMON_ENVELOPE,
    STAR_2_PROPERTY::LAMBDA_AT_COMMON_ENVELOPE,
    STAR_1_PROPERTY::BINDING_ENERGY_PRE_COMMON_ENVELOPE,
    STAR_2_PROPERTY::BINDING_ENERGY_PRE_COMMON_ENVELOPE,
    BINARY_PROPERTY::ECCENTRICITY_PRE_COMMON_ENVELOPE,
    BINARY_PROPERTY::ECCENTRICITY_POST_COMMON_ENVELOPE,
    BINARY_PROPERTY::SEMI_MAJOR_AXIS_PRE_COMMON_ENVELOPE,
    BINARY_PROPERTY::SEMI_MAJOR_AXIS_POST_COMMON_ENVELOPE,
    BINARY_PROPERTY::ROCHE_LOBE_RADIUS_1_PRE_COMMON_ENVELOPE,
    BINARY_PROPERTY::ROCHE_LOBE_RADIUS_1_POST_COMMON_ENVELOPE,
    BINARY_PROPERTY::ROCHE_LOBE_RADIUS_2_PRE_COMMON_ENVELOPE,
    BINARY_PROPERTY::ROCHE_LOBE_RADIUS_2_POST_COMMON_ENVELOPE,
    BINARY_PROPERTY::MASS_1_PRE_COMMON_ENVELOPE,
    BINARY_PROPERTY::MASS_1_POST_COMMON_ENVELOPE,
    BINARY_PROPERTY::MASS_ENV_1,
    BINARY_PROPERTY::RADIUS_1_PRE_COMMON_ENVELOPE,
    BINARY_PROPERTY::RADIUS_1_POST_COMMON_ENVELOPE,
    BINARY_PROPERTY::STELLAR_TYPE_1_PRE_COMMON_ENVELOPE,
    STAR_1_PROPERTY::STELLAR_TYPE,
    STAR_1_PROPERTY::LAMBDA_FIXED,
    STAR_1_PROPERTY::LAMBDA_NANJING,
    STAR_1_PROPERTY::LAMBDA_LOVERIDGE,
    STAR_1_PROPERTY::LAMBDA_LOVERIDGE_WINDS,
    STAR_1_PROPERTY::LAMBDA_KRUCKOW,
    STAR_1_PROPERTY::BINDING_ENERGY_FIXED,
    STAR_1_PROPERTY::BINDING_ENERGY_NANJING,
    STAR_1_PROPERTY::BINDING_ENERGY_LOVERIDGE,
    STAR_1_PROPERTY::BINDING_ENERGY_LOVERIDGE_WINDS,
    STAR_1_PROPERTY::BINDING_ENERGY_KRUCKOW,
    BINARY_PROPERTY::MASS_2_PRE_COMMON_ENVELOPE,
    BINARY_PROPERTY::MASS_2_POST_COMMON_ENVELOPE,
    BINARY_PROPERTY::MASS_ENV_2,
    BINARY_PROPERTY::RADIUS_2_PRE_COMMON_ENVELOPE,
    BINARY_PROPERTY::RADIUS_2_POST_COMMON_ENVELOPE,
    BINARY_PROPERTY::STELLAR_TYPE_2_PRE_COMMON_ENVELOPE,
    STAR_2_PROPERTY::STELLAR_TYPE,
    STAR_2_PROPERTY::LAMBDA_FIXED,
    STAR_2_PROPERTY::LAMBDA_NANJING,
    STAR_2_PROPERTY::LAMBDA_LOVERIDGE,
    STAR_2_PROPERTY::LAMBDA_LOVERIDGE_WINDS,
    STAR_2_PROPERTY::LAMBDA_KRUCKOW,
    STAR_2_PROPERTY::BINDING_ENERGY_FIXED,
    STAR_2_PROPERTY::BINDING_ENERGY_NANJING,
    STAR_2_PROPERTY::BINDING_ENERGY_LOVERIDGE,
    STAR_2_PROPERTY::BINDING_ENERGY_LOVERIDGE_WINDS,
    STAR_2_PROPERTY::BINDING_ENERGY_KRUCKOW,
    BINARY_PROPERTY::MASS_TRANSFER_TRACKER_HISTORY,
    BINARY_PROPERTY::STELLAR_MERGER,
    BINARY_PROPERTY::OPTIMISTIC_COMMON_ENVELOPE,
    BINARY_PROPERTY::COMMON_ENVELOPE_EVENT_COUNT,
    BINARY_PROPERTY::DOUBLE_CORE_COMMON_ENVELOPE,
    STAR_1_PROPERTY::IS_RLOF,
    STAR_1_PROPERTY::LUMINOSITY_PRE_COMMON_ENVELOPE,
    STAR_1_PROPERTY::TEMPERATURE_PRE_COMMON_ENVELOPE,
    STAR_1_PROPERTY::DYNAMICAL_TIMESCALE_PRE_COMMON_ENVELOPE,
    STAR_1_PROPERTY::THERMAL_TIMESCALE_PRE_COMMON_ENVELOPE,
    STAR_2_PROPERTY::IS_RLOF,
    STAR_2_PROPERTY::LUMINOSITY_PRE_COMMON_ENVELOPE,
    STAR_2_PROPERTY::TEMPERATURE_PRE_COMMON_ENVELOPE,
    STAR_2_PROPERTY::DYNAMICAL_TIMESCALE_PRE_COMMON_ENVELOPE,
    STAR_2_PROPERTY::THERMAL_TIMESCALE_PRE_COMMON_ENVELOPE,
    BINARY_PROPERTY::ZETA_STAR,
    BINARY_PROPERTY::ZETA_LOBE,
    BINARY_PROPERTY::SYNCHRONIZATION_TIMESCALE,
    BINARY_PROPERTY::CIRCULARIZATION_TIMESCALE,
    STAR_1_PROPERTY::RADIAL_EXPANSION_TIMESCALE_PRE_COMMON_ENVELOPE,
    STAR_2_PROPERTY::RADIAL_EXPANSION_TIMESCALE_PRE_COMMON_ENVELOPE,
    BINARY_PROPERTY::IMMEDIATE_RLOF_POST_COMMON_ENVELOPE,
    BINARY_PROPERTY::SIMULTANEOUS_RLOF
};


// BSE_DETAILED_OUTPUT_REC
//
// Default record definition for the BSE Detailed Output logfile
//
const ANY_PROPERTY_VECTOR BSE_DETAILED_OUTPUT_REC = {
    BINARY_PROPERTY::RANDOM_SEED,
    BINARY_PROPERTY::DT,
    BINARY_PROPERTY::TIME,
    BINARY_PROPERTY::UNBOUND,
    BINARY_PROPERTY::SEMI_MAJOR_AXIS_RSOL,
    BINARY_PROPERTY::ECCENTRICITY,
    STAR_1_PROPERTY::MZAMS,
    STAR_2_PROPERTY::MZAMS,
    STAR_1_PROPERTY::MASS_0,
    STAR_2_PROPERTY::MASS_0,
    STAR_1_PROPERTY::MASS,
    STAR_2_PROPERTY::MASS,
    STAR_1_PROPERTY::ENV_MASS,
    STAR_2_PROPERTY::ENV_MASS,
    STAR_1_PROPERTY::CORE_MASS,
    STAR_2_PROPERTY::CORE_MASS,
    STAR_1_PROPERTY::HE_CORE_MASS,
    STAR_2_PROPERTY::HE_CORE_MASS,
    STAR_1_PROPERTY::CO_CORE_MASS,
    STAR_2_PROPERTY::CO_CORE_MASS,
    STAR_1_PROPERTY::RADIUS,
    STAR_2_PROPERTY::RADIUS,
    BINARY_PROPERTY::ROCHE_LOBE_RADIUS_1,
    BINARY_PROPERTY::ROCHE_LOBE_RADIUS_2,
    STAR_1_PROPERTY::OMEGA,
    STAR_2_PROPERTY::OMEGA,
    STAR_1_PROPERTY::OMEGA_BREAK,
    STAR_2_PROPERTY::OMEGA_BREAK,
    STAR_1_PROPERTY::INITIAL_STELLAR_TYPE,
    STAR_2_PROPERTY::INITIAL_STELLAR_TYPE,
    STAR_1_PROPERTY::STELLAR_TYPE,
    STAR_2_PROPERTY::STELLAR_TYPE,
    STAR_1_PROPERTY::AGE,
    STAR_2_PROPERTY::AGE,
    STAR_1_PROPERTY::LUMINOSITY,
    STAR_2_PROPERTY::LUMINOSITY,
    STAR_1_PROPERTY::TEMPERATURE,
    STAR_2_PROPERTY::TEMPERATURE,
    STAR_1_PROPERTY::ANGULAR_MOMENTUM,
    STAR_2_PROPERTY::ANGULAR_MOMENTUM,
    STAR_1_PROPERTY::DYNAMICAL_TIMESCALE,
    STAR_2_PROPERTY::DYNAMICAL_TIMESCALE,
    STAR_1_PROPERTY::THERMAL_TIMESCALE,
    STAR_2_PROPERTY::THERMAL_TIMESCALE,
    STAR_1_PROPERTY::ZETA_SOBERMAN,
    STAR_2_PROPERTY::ZETA_SOBERMAN,
    STAR_1_PROPERTY::ZETA_SOBERMAN_HE,
    STAR_2_PROPERTY::ZETA_SOBERMAN_HE,
    STAR_1_PROPERTY::ZETA_HURLEY,
    STAR_2_PROPERTY::ZETA_HURLEY,
    STAR_1_PROPERTY::ZETA_HURLEY_HE,
    STAR_2_PROPERTY::ZETA_HURLEY_HE,
    STAR_1_PROPERTY::MASS_LOSS_DIFF,
    STAR_2_PROPERTY::MASS_LOSS_DIFF,
    STAR_1_PROPERTY::DOMINANT_MASS_LOSS_RATE,
    STAR_2_PROPERTY::DOMINANT_MASS_LOSS_RATE,
    STAR_1_PROPERTY::MASS_TRANSFER_DIFF,
    STAR_2_PROPERTY::MASS_TRANSFER_DIFF,
    STAR_1_PROPERTY::MDOT,
    STAR_2_PROPERTY::MDOT,
    BINARY_PROPERTY::TOTAL_ANGULAR_MOMENTUM,
    BINARY_PROPERTY::TOTAL_ENERGY,
    STAR_1_PROPERTY::METALLICITY,
    STAR_2_PROPERTY::METALLICITY,
    BINARY_PROPERTY::MASS_TRANSFER_TRACKER_HISTORY,
    STAR_1_PROPERTY::PULSAR_MAGNETIC_FIELD,
    STAR_2_PROPERTY::PULSAR_MAGNETIC_FIELD,
    STAR_1_PROPERTY::PULSAR_SPIN_FREQUENCY,
    STAR_2_PROPERTY::PULSAR_SPIN_FREQUENCY,
    STAR_1_PROPERTY::PULSAR_SPIN_DOWN_RATE,
    STAR_2_PROPERTY::PULSAR_SPIN_DOWN_RATE,
    STAR_1_PROPERTY::PULSAR_BIRTH_PERIOD,
    STAR_2_PROPERTY::PULSAR_BIRTH_PERIOD,
    STAR_1_PROPERTY::PULSAR_BIRTH_SPIN_DOWN_RATE,
    STAR_2_PROPERTY::PULSAR_BIRTH_SPIN_DOWN_RATE,
    STAR_1_PROPERTY::RADIAL_EXPANSION_TIMESCALE,
    STAR_2_PROPERTY::RADIAL_EXPANSION_TIMESCALE,
    BINARY_PROPERTY::RLOF_MASS_LOSS_RATE,
    BINARY_PROPERTY::RLOF_ACCRETION_EFFICIENCY
};


// BSE_DOUBLE_COMPACT_OBJECT_REC
//
// Default record definition for the Double Compact Objects logfile
//
const ANY_PROPERTY_VECTOR BSE_DOUBLE_COMPACT_OBJECTS_REC = {
    BINARY_PROPERTY::RANDOM_SEED,
    BINARY_PROPERTY::SEMI_MAJOR_AXIS_AT_DCO_FORMATION, 
    BINARY_PROPERTY::ECCENTRICITY_AT_DCO_FORMATION,
    STAR_1_PROPERTY::MASS,
    STAR_1_PROPERTY::STELLAR_TYPE,
    STAR_2_PROPERTY::MASS, 
    STAR_2_PROPERTY::STELLAR_TYPE,
    BINARY_PROPERTY::TIME_TO_COALESCENCE,
    BINARY_PROPERTY::TIME,
    BINARY_PROPERTY::MERGES_IN_HUBBLE_TIME, 
    STAR_1_PROPERTY::RECYCLED_NEUTRON_STAR,  
    STAR_2_PROPERTY::RECYCLED_NEUTRON_STAR
};


// BSE_PULSAR_EVOLUTION_REC
//
// Default record definition for the BSE Pulsar Evolution logfile
//
const ANY_PROPERTY_VECTOR BSE_PULSAR_EVOLUTION_REC = {
    BINARY_PROPERTY::RANDOM_SEED,
    STAR_1_PROPERTY::MASS,
    STAR_2_PROPERTY::MASS,
    STAR_1_PROPERTY::STELLAR_TYPE,
    STAR_2_PROPERTY::STELLAR_TYPE,
    BINARY_PROPERTY::SEMI_MAJOR_AXIS_RSOL,
    BINARY_PROPERTY::MASS_TRANSFER_TRACKER_HISTORY,
    STAR_1_PROPERTY::PULSAR_MAGNETIC_FIELD,
    STAR_2_PROPERTY::PULSAR_MAGNETIC_FIELD,
    STAR_1_PROPERTY::PULSAR_SPIN_FREQUENCY,
    STAR_2_PROPERTY::PULSAR_SPIN_FREQUENCY,
    STAR_1_PROPERTY::PULSAR_SPIN_DOWN_RATE,
    STAR_2_PROPERTY::PULSAR_SPIN_DOWN_RATE,
    STAR_1_PROPERTY::PULSAR_BIRTH_PERIOD,
    STAR_2_PROPERTY::PULSAR_BIRTH_PERIOD,
    STAR_1_PROPERTY::PULSAR_BIRTH_SPIN_DOWN_RATE,
    STAR_2_PROPERTY::PULSAR_BIRTH_SPIN_DOWN_RATE,
    BINARY_PROPERTY::TIME,
    BINARY_PROPERTY::DT
};


// BSE_RLOF_PARAMETERS_REC
//
// Default record definition for the RLOF Parameters logfile
//
const ANY_PROPERTY_VECTOR BSE_RLOF_PARAMETERS_REC = {
    BINARY_PROPERTY::RANDOM_SEED,
    BINARY_PROPERTY::RLOF_TIME_POST_MT,
    BINARY_PROPERTY::RLOF_TIME_PRE_MT,
    BINARY_PROPERTY::RLOF_POST_MT_STAR1_MASS,
    BINARY_PROPERTY::RLOF_POST_MT_STAR2_MASS,
    BINARY_PROPERTY::RLOF_POST_MT_STAR1_RADIUS,
    BINARY_PROPERTY::RLOF_POST_MT_STAR2_RADIUS,
    BINARY_PROPERTY::RLOF_POST_MT_STAR1_STELLAR_TYPE,
    BINARY_PROPERTY::RLOF_POST_MT_STAR2_STELLAR_TYPE,
    BINARY_PROPERTY::RLOF_POST_MT_SEMI_MAJOR_AXIS,
    BINARY_PROPERTY::RLOF_POST_MT_ECCENTRICITY,
    BINARY_PROPERTY::RLOF_POST_MT_EVENT_COUNTER,
    BINARY_PROPERTY::RLOF_POST_MT_STAR1_RLOF,
    BINARY_PROPERTY::RLOF_POST_MT_STAR2_RLOF,
    BINARY_PROPERTY::STELLAR_MERGER,
    BINARY_PROPERTY::RLOF_POST_MT_COMMON_ENVELOPE,
    BINARY_PROPERTY::RLOF_POST_STEP_STAR_TO_ROCHE_LOBE_RADIUS_RATIO_1,
    BINARY_PROPERTY::RLOF_POST_STEP_STAR_TO_ROCHE_LOBE_RADIUS_RATIO_2,
    BINARY_PROPERTY::RLOF_PRE_MT_STAR1_MASS,
    BINARY_PROPERTY::RLOF_PRE_MT_STAR2_MASS,
    BINARY_PROPERTY::RLOF_PRE_MT_STAR1_RADIUS,
    BINARY_PROPERTY::RLOF_PRE_MT_STAR2_RADIUS,
    BINARY_PROPERTY::RLOF_PRE_MT_STAR1_STELLAR_TYPE,
    BINARY_PROPERTY::RLOF_PRE_MT_STAR2_STELLAR_TYPE,
    BINARY_PROPERTY::RLOF_PRE_MT_SEMI_MAJOR_AXIS,
    BINARY_PROPERTY::RLOF_PRE_MT_ECCENTRICITY,
    BINARY_PROPERTY::RLOF_PRE_MT_STAR1_RLOF,
    BINARY_PROPERTY::RLOF_PRE_MT_STAR2_RLOF,
    BINARY_PROPERTY::RLOF_PRE_STEP_STAR_TO_ROCHE_LOBE_RADIUS_RATIO_1,
    BINARY_PROPERTY::RLOF_PRE_STEP_STAR_TO_ROCHE_LOBE_RADIUS_RATIO_2,
    BINARY_PROPERTY::RLOF_ACCRETION_EFFICIENCY,
    BINARY_PROPERTY::RLOF_MASS_LOSS_RATE,
    BINARY_PROPERTY::RLOF_MASS_TRANSFER_TIMESCALE,
    STAR_1_PROPERTY::ZETA_SOBERMAN,
    STAR_1_PROPERTY::ZETA_SOBERMAN_HE,
    STAR_1_PROPERTY::ZETA_HURLEY,
    STAR_1_PROPERTY::ZETA_HURLEY_HE,
    STAR_2_PROPERTY::ZETA_SOBERMAN,
    STAR_2_PROPERTY::ZETA_SOBERMAN_HE,
    STAR_2_PROPERTY::ZETA_HURLEY,
    STAR_2_PROPERTY::ZETA_HURLEY_HE
};


// BSE_SUPERNOVAE_REC
//
// Default record definition for the BSE Supernovae logfile
//
const ANY_PROPERTY_VECTOR BSE_SUPERNOVAE_REC = {
    BINARY_PROPERTY::RANDOM_SEED,
    SUPERNOVA_PROPERTY::DRAWN_KICK_MAGNITUDE,
    SUPERNOVA_PROPERTY::KICK_MAGNITUDE,
    SUPERNOVA_PROPERTY::FALLBACK_FRACTION,
    BINARY_PROPERTY::ORBITAL_VELOCITY_PRE_SUPERNOVA,
    SUPERNOVA_PROPERTY::MEAN_ANOMALY,
    SUPERNOVA_PROPERTY::SUPERNOVA_THETA,
    SUPERNOVA_PROPERTY::SUPERNOVA_PHI,
    SUPERNOVA_PROPERTY::SN_TYPE,
    BINARY_PROPERTY::ECCENTRICITY_PRE_SUPERNOVA,  
    BINARY_PROPERTY::ECCENTRICITY,
    BINARY_PROPERTY::SEMI_MAJOR_AXIS_PRE_SUPERNOVA_RSOL,
    BINARY_PROPERTY::SEMI_MAJOR_AXIS_RSOL,
    BINARY_PROPERTY::TIME,
    BINARY_PROPERTY::SUPERNOVA_STATE,
    BINARY_PROPERTY::UNBOUND,
    COMPANION_PROPERTY::STELLAR_TYPE,
    SUPERNOVA_PROPERTY::STELLAR_TYPE,
    SUPERNOVA_PROPERTY::STELLAR_TYPE_PREV,
    COMPANION_PROPERTY::MASS,
    SUPERNOVA_PROPERTY::TOTAL_MASS_AT_COMPACT_OBJECT_FORMATION,
    SUPERNOVA_PROPERTY::CORE_MASS_AT_COMPACT_OBJECT_FORMATION,
    SUPERNOVA_PROPERTY::CO_CORE_MASS_AT_COMPACT_OBJECT_FORMATION,
    SUPERNOVA_PROPERTY::HE_CORE_MASS_AT_COMPACT_OBJECT_FORMATION,
    SUPERNOVA_PROPERTY::MASS,
    SUPERNOVA_PROPERTY::EXPERIENCED_RLOF,
    SUPERNOVA_PROPERTY::MASS_TRANSFER_DONOR_HISTORY,
    SUPERNOVA_PROPERTY::SPEED,
    COMPANION_PROPERTY::SPEED,
    BINARY_PROPERTY::SYSTEMIC_SPEED,
    SUPERNOVA_PROPERTY::IS_HYDROGEN_POOR,
    BINARY_PROPERTY::SUPERNOVA_ORBIT_INCLINATION_ANGLE, 
};


// BSE_SWITCH_LOG
//
// Default record definition for the BSE Switch Log logfile
//
const ANY_PROPERTY_VECTOR BSE_SWITCH_LOG_REC = {
    BINARY_PROPERTY::RANDOM_SEED,
    BINARY_PROPERTY::TIME
};


// BSE_SYSTEM_PARAMETERS_REC
//
// Default record definition for the System Parameters logfile
//
const ANY_PROPERTY_VECTOR BSE_SYSTEM_PARAMETERS_REC = {
    BINARY_PROPERTY::RANDOM_SEED,
    STAR_1_PROPERTY::MZAMS,
    STAR_2_PROPERTY::MZAMS,
    BINARY_PROPERTY::SEMI_MAJOR_AXIS_INITIAL,
    BINARY_PROPERTY::ECCENTRICITY_INITIAL,
    STAR_1_PROPERTY::SUPERNOVA_KICK_MAGNITUDE_RANDOM_NUMBER,
    STAR_2_PROPERTY::SUPERNOVA_KICK_MAGNITUDE_RANDOM_NUMBER,
    STAR_1_PROPERTY::OMEGA_ZAMS,
    STAR_2_PROPERTY::OMEGA_ZAMS,
    PROGRAM_OPTION::KICK_MAGNITUDE_DISTRIBUTION_SIGMA_CCSN_NS,
    PROGRAM_OPTION::KICK_MAGNITUDE_DISTRIBUTION_SIGMA_CCSN_BH,
    PROGRAM_OPTION::KICK_MAGNITUDE_DISTRIBUTION_SIGMA_FOR_ECSN,
    PROGRAM_OPTION::KICK_MAGNITUDE_DISTRIBUTION_SIGMA_FOR_USSN,
    PROGRAM_OPTION::LBV_FACTOR,
    PROGRAM_OPTION::WR_FACTOR,
    PROGRAM_OPTION::COMMON_ENVELOPE_ALPHA,
    STAR_1_PROPERTY::METALLICITY,
    STAR_2_PROPERTY::METALLICITY,
    BINARY_PROPERTY::UNBOUND,
    BINARY_PROPERTY::STELLAR_MERGER,
    BINARY_PROPERTY::STELLAR_MERGER_AT_BIRTH,
    BINARY_PROPERTY::MASSES_EQUILIBRATED_AT_BIRTH,
    STAR_1_PROPERTY::INITIAL_STELLAR_TYPE,
    STAR_1_PROPERTY::STELLAR_TYPE,
    STAR_2_PROPERTY::INITIAL_STELLAR_TYPE,
    STAR_2_PROPERTY::STELLAR_TYPE,
    STAR_1_PROPERTY::CHEMICALLY_HOMOGENEOUS_MAIN_SEQUENCE,
    STAR_2_PROPERTY::CHEMICALLY_HOMOGENEOUS_MAIN_SEQUENCE,
    BINARY_PROPERTY::EVOL_STATUS,
    BINARY_PROPERTY::ERROR,
    PROGRAM_OPTION::NOTES
};


// SSE output record definitions

// SSE_DETAILED_OUTPUT_REC
//
// Default record definition for the SSE Detailed Output logfile
//
const ANY_PROPERTY_VECTOR SSE_DETAILED_OUTPUT_REC = {
    STAR_PROPERTY::AGE,
    STAR_PROPERTY::DT,
    STAR_PROPERTY::TIME,
    STAR_PROPERTY::STELLAR_TYPE,
    STAR_PROPERTY::METALLICITY,
    STAR_PROPERTY::MASS_0,
    STAR_PROPERTY::MASS,
    STAR_PROPERTY::RADIUS,
    STAR_PROPERTY::RZAMS,
    STAR_PROPERTY::LUMINOSITY,
    STAR_PROPERTY::TEMPERATURE,
    STAR_PROPERTY::CORE_MASS,
    STAR_PROPERTY::CO_CORE_MASS,
    STAR_PROPERTY::HE_CORE_MASS,
    STAR_PROPERTY::MDOT,
    STAR_PROPERTY::DOMINANT_MASS_LOSS_RATE,
    STAR_PROPERTY::TIMESCALE_MS
};

// SSE_PULSAR_EVOLUTION_REC
//
// Default record definition for the SSE Pulsar Evolution logfile
//
const ANY_PROPERTY_VECTOR SSE_PULSAR_EVOLUTION_REC = {
    STAR_PROPERTY::RANDOM_SEED,
    STAR_PROPERTY::MASS,
    STAR_PROPERTY::STELLAR_TYPE,
    STAR_PROPERTY::PULSAR_MAGNETIC_FIELD,
    STAR_PROPERTY::PULSAR_SPIN_FREQUENCY,
    STAR_PROPERTY::PULSAR_SPIN_DOWN_RATE,
    STAR_PROPERTY::PULSAR_BIRTH_PERIOD,
    STAR_PROPERTY::PULSAR_BIRTH_SPIN_DOWN_RATE,
    STAR_PROPERTY::TIME,
    STAR_PROPERTY::DT
};

// SSE_SUPERNOVAE_REC
//
// Default record definition for the SSE Supernovae logfile
//
const ANY_PROPERTY_VECTOR SSE_SUPERNOVAE_REC = {
    STAR_PROPERTY::RANDOM_SEED,
    STAR_PROPERTY::DRAWN_KICK_MAGNITUDE,
    STAR_PROPERTY::KICK_MAGNITUDE,
    STAR_PROPERTY::FALLBACK_FRACTION,
    STAR_PROPERTY::MEAN_ANOMALY,				
    STAR_PROPERTY::SN_TYPE,
    STAR_PROPERTY::TOTAL_MASS_AT_COMPACT_OBJECT_FORMATION,
    STAR_PROPERTY::CO_CORE_MASS_AT_COMPACT_OBJECT_FORMATION,
    STAR_PROPERTY::MASS,
    STAR_PROPERTY::STELLAR_TYPE,
    STAR_PROPERTY::STELLAR_TYPE_PREV,
    STAR_PROPERTY::CORE_MASS_AT_COMPACT_OBJECT_FORMATION,
    STAR_PROPERTY::HE_CORE_MASS_AT_COMPACT_OBJECT_FORMATION,
    STAR_PROPERTY::TIME,
    STAR_PROPERTY::IS_HYDROGEN_POOR
};


// SSE_SWITCH_LOG
//
// Default record definition for the SSE Switch Log logfile
//
const ANY_PROPERTY_VECTOR SSE_SWITCH_LOG_REC = {
    STAR_PROPERTY::RANDOM_SEED,
    STAR_PROPERTY::TIME
};


// SSE_SYSTEM_PARAMETERS_REC
//
// Default record definition for the SSE System Parameters logfile
//
const ANY_PROPERTY_VECTOR SSE_SYSTEM_PARAMETERS_REC = {
    STAR_PROPERTY::RANDOM_SEED,
    STAR_PROPERTY::MZAMS,
    STAR_PROPERTY::RZAMS,
    STAR_PROPERTY::METALLICITY,
    STAR_PROPERTY::OMEGA_ZAMS,
    STAR_PROPERTY::INITIAL_STELLAR_TYPE,
    STAR_PROPERTY::STELLAR_TYPE,
    STAR_PROPERTY::SUPERNOVA_KICK_MAGNITUDE_RANDOM_NUMBER,
    STAR_PROPERTY::MASS,
    STAR_PROPERTY::CHEMICALLY_HOMOGENEOUS_MAIN_SEQUENCE,
    STAR_PROPERTY::EVOL_STATUS,
    STAR_PROPERTY::ERROR,
    PROGRAM_OPTION::KICK_MAGNITUDE_DISTRIBUTION_SIGMA_CCSN_NS,
    PROGRAM_OPTION::KICK_MAGNITUDE_DISTRIBUTION_SIGMA_CCSN_BH,
    PROGRAM_OPTION::KICK_MAGNITUDE_DISTRIBUTION_SIGMA_FOR_ECSN,
    PROGRAM_OPTION::KICK_MAGNITUDE_DISTRIBUTION_SIGMA_FOR_USSN,
    PROGRAM_OPTION::LBV_FACTOR,
    PROGRAM_OPTION::WR_FACTOR,
    PROGRAM_OPTION::NOTES
};

// enum class LOGFILE_CLASS
// Symbolic names for logfile types
enum class LOGFILE_CLASS: int { NONE, STELLAR, BINARY };

// descriptors for logfiles
// unordered_map - key is integer logfile (from enum class LOGFILE above)
// fields are: {default filename, record descriptor, short file name, short record name, type}
// (the short names are for logfile definitions file parsing)
typedef std::tuple<std::string, ANY_PROPERTY_VECTOR, std::string, std::string, LOGFILE_CLASS> LOGFILE_DESCRIPTOR_T;
const std::map<LOGFILE, LOGFILE_DESCRIPTOR_T> LOGFILE_DESCRIPTOR = {
    { LOGFILE::NONE,                       { "" ,                              {},                             "",                "",                    LOGFILE_CLASS::NONE}},

    { LOGFILE::DEBUG_LOG,                  { "Debug_Log",                      {},                             "",                "",                    LOGFILE_CLASS::NONE }},
    { LOGFILE::ERROR_LOG,                  { "Error_Log",                      {},                             "",                "",                    LOGFILE_CLASS::NONE }},

    { LOGFILE::BSE_COMMON_ENVELOPES,       { "BSE_Common_Envelopes",           BSE_COMMON_ENVELOPES_REC,       "BSE_CEE",         "BSE_CEE_REC",         LOGFILE_CLASS::BINARY }},
    { LOGFILE::BSE_DETAILED_OUTPUT,        { "BSE_Detailed_Output",            BSE_DETAILED_OUTPUT_REC,        "BSE_DETAILED",    "BSE_DETAILED_REC",    LOGFILE_CLASS::BINARY }},
    { LOGFILE::BSE_DOUBLE_COMPACT_OBJECTS, { "BSE_Double_Compact_Objects",     BSE_DOUBLE_COMPACT_OBJECTS_REC, "BSE_DCO",         "BSE_DCO_REC",         LOGFILE_CLASS::BINARY }},
    { LOGFILE::BSE_PULSAR_EVOLUTION,       { "BSE_Pulsar_Evolution",           BSE_PULSAR_EVOLUTION_REC,       "BSE_PULSARS",     "BSE_PULSARS_REC",     LOGFILE_CLASS::BINARY }},
    { LOGFILE::BSE_RLOF_PARAMETERS,        { "BSE_RLOF",                       BSE_RLOF_PARAMETERS_REC,        "BSE_RLOF",        "BSE_RLOF_REC",        LOGFILE_CLASS::BINARY }},
    { LOGFILE::BSE_SUPERNOVAE,             { "BSE_Supernovae",                 BSE_SUPERNOVAE_REC,             "BSE_SNE",         "BSE_SNE_REC",         LOGFILE_CLASS::BINARY }},
    { LOGFILE::BSE_SWITCH_LOG,             { "BSE_Switch_Log",                 BSE_SWITCH_LOG_REC,             "BSE_SWITCH_LOG",  "BSE_SWITCH_REC",      LOGFILE_CLASS::BINARY }},
    { LOGFILE::BSE_SYSTEM_PARAMETERS,      { "BSE_System_Parameters",          BSE_SYSTEM_PARAMETERS_REC,      "BSE_SYSPARMS",    "BSE_SYSPARMS_REC",    LOGFILE_CLASS::BINARY }},

    { LOGFILE::SSE_DETAILED_OUTPUT,        { "SSE_Detailed_Output",            SSE_DETAILED_OUTPUT_REC,        "SSE_DETAILED",    "SSE_DETAILED_REC",    LOGFILE_CLASS::STELLAR }},
    { LOGFILE::SSE_PULSAR_EVOLUTION,       { "SSE_Pulsar_Evolution",           SSE_PULSAR_EVOLUTION_REC,       "SSE_PULSARS",     "SSE_PULSARS_REC",     LOGFILE_CLASS::STELLAR }},
    { LOGFILE::SSE_SUPERNOVAE,             { "SSE_Supernovae",                 SSE_SUPERNOVAE_REC,             "SSE_SNE",         "SSE_SNE_REC",         LOGFILE_CLASS::STELLAR }},
    { LOGFILE::SSE_SWITCH_LOG,             { "SSE_Switch_Log",                 SSE_SWITCH_LOG_REC,             "SSE_SWITCH_LOG",  "SSE_SWITCH_REC",      LOGFILE_CLASS::STELLAR }},
    { LOGFILE::SSE_SYSTEM_PARAMETERS,      { "SSE_System_Parameters",          SSE_SYSTEM_PARAMETERS_REC,      "SSE_SYSPARMS",    "SSE_SYSPARMS_REC",    LOGFILE_CLASS::STELLAR }}
};

#endif // __LogTypedefs_h__
