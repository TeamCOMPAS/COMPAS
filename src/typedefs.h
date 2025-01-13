#ifndef __typedefs_h__
#define __typedefs_h__


// This is where developer-defined types are defined - except for types that pertain directly to
// the COMPAS looging functionality (including the definition of the default record composition
// for the various log files) - those are listed in LogTypedefs.h


#include "EnumHash.h"
#include "LogTypedefs.h"
#include "ErrorCatalog.h"
#include <boost/math/tools/roots.hpp>
#include <boost/numeric/odeint.hpp>


// JR: todo: clean this up and document it better


// Bitwise operators for Enum Class - |, |=, &, &=, ^, ^=, ~ only
// from http://blog.bitwigglers.org/using-enum-classes-as-type-safe-bitmasks/
#define ENABLE_BITMASK_OPERATORS(x)     \
template<>                              \
struct EnableBitMaskOperators<x> {      \
    static const bool enable = true;    \
};

template<typename Enum>  
struct EnableBitMaskOperators {
    static const bool enable = false;
};

template<typename Enum>  
typename std::enable_if<EnableBitMaskOperators<Enum>::enable, Enum>::type  
operator |(Enum lhs, Enum rhs) {
    return static_cast<Enum> (
        static_cast<typename std::underlying_type<Enum>::type>(lhs) |
        static_cast<typename std::underlying_type<Enum>::type>(rhs)
    );
}

template<typename Enum>  
typename std::enable_if<EnableBitMaskOperators<Enum>::enable, Enum>::type  
operator |=(Enum &lhs, Enum rhs) {
    lhs = static_cast<Enum> (
        static_cast<typename std::underlying_type<Enum>::type>(lhs) |
        static_cast<typename std::underlying_type<Enum>::type>(rhs)           
    );

    return lhs;
}

template<typename Enum>  
typename std::enable_if<EnableBitMaskOperators<Enum>::enable, Enum>::type 
operator &(Enum lhs, Enum rhs) {
    return static_cast<Enum> (
        static_cast<typename std::underlying_type<Enum>::type>(lhs) &
        static_cast<typename std::underlying_type<Enum>::type>(rhs)
    );
}

template<typename Enum>  
typename std::enable_if<EnableBitMaskOperators<Enum>::enable, Enum>::type 
operator &=(Enum &lhs, Enum rhs) {
    lhs = static_cast<Enum> (
        static_cast<typename std::underlying_type<Enum>::type>(lhs) &
        static_cast<typename std::underlying_type<Enum>::type>(rhs)           
    );

    return lhs;
}

template<typename Enum>  
typename std::enable_if<EnableBitMaskOperators<Enum>::enable, Enum>::type 
operator ^(Enum lhs, Enum rhs) {
    return static_cast<Enum> (
        static_cast<typename std::underlying_type<Enum>::type>(lhs) ^
        static_cast<typename std::underlying_type<Enum>::type>(rhs)
    );
}

template<typename Enum>  
typename std::enable_if<EnableBitMaskOperators<Enum>::enable, Enum>::type 
operator ^=(Enum &lhs, Enum rhs) {
    lhs = static_cast<Enum> (
        static_cast<typename std::underlying_type<Enum>::type>(lhs) ^
        static_cast<typename std::underlying_type<Enum>::type>(rhs)           
    );

    return lhs;
}

template<typename Enum>  
typename std::enable_if<EnableBitMaskOperators<Enum>::enable, Enum>::type 
operator ~(Enum rhs) {
    return static_cast<Enum> (
        ~static_cast<typename std::underlying_type<Enum>::type>(rhs)
    );
}


// enum class types
// ================
//
// listed alphabetically (with the exception of stellar types - listed first)
// categories might work, but for now alphabetically makes things easy to find
//
// some enum class have associated maps allowing lookup of description, numerical velue, etc. associated
// with the enum class entries - these/ maps are define here and listed with the enum class.
//
// the order of entries is most enum classes is not significant - where entries are not set to a specific
// integer value, the integer value for the entry is its ordinal position.  Code should not rely on these
// being/remaining in any specific order.
//
// the value of entries in an enum class must be unique
//
// some enum classes are defined in constants.h because they are needed there for constants definition
//
// some typedefs are listed after enum classes (because they may use enum classes)


// we list stellar types and associated initialiser lists first so they are
// grouped and easy to find
//
// these are symolic names for the stellar types from Hurley et al. 2000
//
// the integer value for any stellar type is its ordinal position
// note that the order of entries is not significant - the code should not rely on these being in any order
// (and so the stellar type symbols having any value or order.  e.g. it is not guaranteed that the integer
// value for stellar type NEUTRON_STAR will be > the integr value for stellar type HERTZSPRUNG_GAP)
enum class STELLAR_TYPE: int {                      // Hurley
    MS_LTE_07,                                      //   0
    MS_GT_07,                                       //   1
    HERTZSPRUNG_GAP,                                //   2
    FIRST_GIANT_BRANCH,                             //   3
    CORE_HELIUM_BURNING,                            //   4
    EARLY_ASYMPTOTIC_GIANT_BRANCH,                  //   5
    THERMALLY_PULSING_ASYMPTOTIC_GIANT_BRANCH,      //   6
    NAKED_HELIUM_STAR_MS,                           //   7
    NAKED_HELIUM_STAR_HERTZSPRUNG_GAP,              //   8
    NAKED_HELIUM_STAR_GIANT_BRANCH,                 //   9
    HELIUM_WHITE_DWARF,                             //  10
    CARBON_OXYGEN_WHITE_DWARF,                      //  11
    OXYGEN_NEON_WHITE_DWARF,                        //  12
    NEUTRON_STAR,                                   //  13
    BLACK_HOLE,                                     //  14
    MASSLESS_REMNANT,                               //  15
    CHEMICALLY_HOMOGENEOUS,                         //  16  : this is here to preserve the Hurley type numbers, but note that Hurley type number progression doesn't necessarily indicate class inheritance
    STAR,                                           //  17  : star is created this way, then switches as required (down here so stellar types consistent with Hurley et al. 2000)
    BINARY_STAR,                                    //  18  : here mainly for diagnostics
    NONE                                            //  19  : here mainly for diagnostics
};
const COMPASUnorderedMap<STELLAR_TYPE, std::string> STELLAR_TYPE_LABEL = {
    { STELLAR_TYPE::MS_LTE_07,                                 "Main_Sequence_<=_0.7" },
    { STELLAR_TYPE::MS_GT_07,                                  "Main_Sequence_>_0.7" },
    { STELLAR_TYPE::HERTZSPRUNG_GAP,                           "Hertzsprung_Gap" },
    { STELLAR_TYPE::FIRST_GIANT_BRANCH,                        "First_Giant_Branch" },
    { STELLAR_TYPE::CORE_HELIUM_BURNING,                       "Core_Helium_Burning" },
    { STELLAR_TYPE::EARLY_ASYMPTOTIC_GIANT_BRANCH,             "Early_Asymptotic_Giant_Branch" },
    { STELLAR_TYPE::THERMALLY_PULSING_ASYMPTOTIC_GIANT_BRANCH, "Thermally_Pulsing_Asymptotic_Giant_Branch" },
    { STELLAR_TYPE::NAKED_HELIUM_STAR_MS,                      "Naked_Helium_Star_MS" },
    { STELLAR_TYPE::NAKED_HELIUM_STAR_HERTZSPRUNG_GAP,         "Naked_Helium_Star_Hertzsprung_Gap" },
    { STELLAR_TYPE::NAKED_HELIUM_STAR_GIANT_BRANCH,            "Naked_Helium_Star_Giant_Branch" },
    { STELLAR_TYPE::HELIUM_WHITE_DWARF,                        "Helium_White_Dwarf" },
    { STELLAR_TYPE::CARBON_OXYGEN_WHITE_DWARF,                 "Carbon-Oxygen_White_Dwarf" },
    { STELLAR_TYPE::OXYGEN_NEON_WHITE_DWARF,                   "Oxygen-Neon_White_Dwarf" },
    { STELLAR_TYPE::NEUTRON_STAR,                              "Neutron_Star" },
    { STELLAR_TYPE::BLACK_HOLE,                                "Black_Hole" },
    { STELLAR_TYPE::MASSLESS_REMNANT,                          "Massless_Remnant" },
    { STELLAR_TYPE::CHEMICALLY_HOMOGENEOUS,                    "Chemically_Homogeneous" },
    { STELLAR_TYPE::STAR,                                      "Star" },
    { STELLAR_TYPE::BINARY_STAR,                               "Binary_Star" },
    { STELLAR_TYPE::NONE,                                      "Not_a_Star!" }
};

// stellar type list initializer
typedef std::initializer_list<STELLAR_TYPE> STELLAR_TYPE_LIST;

// (convenience) initializer list for "evolvable" stellar types
// i.e. not STAR, BINARY_STAR, or NONE
const STELLAR_TYPE_LIST EVOLVABLE_TYPES = {
    STELLAR_TYPE::MS_LTE_07,
    STELLAR_TYPE::MS_GT_07,
    STELLAR_TYPE::HERTZSPRUNG_GAP,
    STELLAR_TYPE::FIRST_GIANT_BRANCH,
    STELLAR_TYPE::CORE_HELIUM_BURNING,
    STELLAR_TYPE::EARLY_ASYMPTOTIC_GIANT_BRANCH,
    STELLAR_TYPE::THERMALLY_PULSING_ASYMPTOTIC_GIANT_BRANCH,
    STELLAR_TYPE::NAKED_HELIUM_STAR_MS,
    STELLAR_TYPE::NAKED_HELIUM_STAR_HERTZSPRUNG_GAP,
    STELLAR_TYPE::NAKED_HELIUM_STAR_GIANT_BRANCH,
    STELLAR_TYPE::HELIUM_WHITE_DWARF,
    STELLAR_TYPE::CARBON_OXYGEN_WHITE_DWARF,
    STELLAR_TYPE::OXYGEN_NEON_WHITE_DWARF,
    STELLAR_TYPE::NEUTRON_STAR,
    STELLAR_TYPE::BLACK_HOLE,
    STELLAR_TYPE::MASSLESS_REMNANT,
    STELLAR_TYPE::CHEMICALLY_HOMOGENEOUS
};

// (convenience) initializer list for MAIN SEQUENCE stars
// (does not include NAKED_HELIUM_STAR_MS)
const STELLAR_TYPE_LIST MAIN_SEQUENCE = {
    STELLAR_TYPE::MS_LTE_07,
    STELLAR_TYPE::MS_GT_07,
    STELLAR_TYPE::CHEMICALLY_HOMOGENEOUS
};

// (convenience) initializer list for ALL MAIN SEQUENCE stars
// (includes NAKED_HELIUM_STAR_MS)
const STELLAR_TYPE_LIST ALL_MAIN_SEQUENCE = {
    STELLAR_TYPE::MS_LTE_07,
    STELLAR_TYPE::MS_GT_07,
    STELLAR_TYPE::CHEMICALLY_HOMOGENEOUS,
    STELLAR_TYPE::NAKED_HELIUM_STAR_MS
};

// (convenience) initializer list for ALL HERTZSPRUNG GAP
// (includes NAKED_HELIUM_STAR_HERTZSPRUNG_GAP)
const STELLAR_TYPE_LIST ALL_HERTZSPRUNG_GAP = {
    STELLAR_TYPE::HERTZSPRUNG_GAP,
    STELLAR_TYPE::NAKED_HELIUM_STAR_HERTZSPRUNG_GAP
};

// (convenience) initializer list for non-COMPACT OBJECTS
const STELLAR_TYPE_LIST NON_COMPACT_OBJECTS = {
    STELLAR_TYPE::MS_LTE_07,
    STELLAR_TYPE::MS_GT_07,
    STELLAR_TYPE::HERTZSPRUNG_GAP,
    STELLAR_TYPE::FIRST_GIANT_BRANCH,
    STELLAR_TYPE::CORE_HELIUM_BURNING,
    STELLAR_TYPE::EARLY_ASYMPTOTIC_GIANT_BRANCH,
    STELLAR_TYPE::THERMALLY_PULSING_ASYMPTOTIC_GIANT_BRANCH,
    STELLAR_TYPE::NAKED_HELIUM_STAR_MS,
    STELLAR_TYPE::NAKED_HELIUM_STAR_HERTZSPRUNG_GAP,
    STELLAR_TYPE::NAKED_HELIUM_STAR_GIANT_BRANCH,
};

// (convenience) initializer list for COMPACT OBJECTS
const STELLAR_TYPE_LIST COMPACT_OBJECTS = {
    STELLAR_TYPE::HELIUM_WHITE_DWARF,
    STELLAR_TYPE::CARBON_OXYGEN_WHITE_DWARF,
    STELLAR_TYPE::OXYGEN_NEON_WHITE_DWARF,
    STELLAR_TYPE::NEUTRON_STAR,
    STELLAR_TYPE::BLACK_HOLE,
    STELLAR_TYPE::MASSLESS_REMNANT
};

// (convenience) initializer list for SN REMNANTS
const STELLAR_TYPE_LIST SN_REMNANTS = {
    STELLAR_TYPE::NEUTRON_STAR,
    STELLAR_TYPE::BLACK_HOLE,
    STELLAR_TYPE::MASSLESS_REMNANT
};

// (convenience) initializer list for GIANTS
const STELLAR_TYPE_LIST GIANTS = {
    STELLAR_TYPE::FIRST_GIANT_BRANCH,
    STELLAR_TYPE::CORE_HELIUM_BURNING,
    STELLAR_TYPE::EARLY_ASYMPTOTIC_GIANT_BRANCH,
    STELLAR_TYPE::THERMALLY_PULSING_ASYMPTOTIC_GIANT_BRANCH
};

// (convenience) initializer list for WHITE DWARFS
const STELLAR_TYPE_LIST WHITE_DWARFS = {
    STELLAR_TYPE::HELIUM_WHITE_DWARF,
    STELLAR_TYPE::CARBON_OXYGEN_WHITE_DWARF,
    STELLAR_TYPE::OXYGEN_NEON_WHITE_DWARF
};

// (convenience) initializer list for He rich stellar types
const STELLAR_TYPE_LIST He_RICH_TYPES = {
    STELLAR_TYPE::NAKED_HELIUM_STAR_MS,
    STELLAR_TYPE::NAKED_HELIUM_STAR_HERTZSPRUNG_GAP,
    STELLAR_TYPE::NAKED_HELIUM_STAR_GIANT_BRANCH,
    STELLAR_TYPE::HELIUM_WHITE_DWARF
};


// start of alphabetical listing

// accretion regimes
// symbolic names for WD accretion regimes
enum class ACCRETION_REGIME: int {
    NONE,   // DEPRECATED June 2024 - remove end 2024 
    ZERO,
    HELIUM_ACCUMULATION,
    HELIUM_FLASHES,
    HELIUM_STABLE_BURNING,
    HELIUM_OPT_THICK_WINDS,
    HYDROGEN_FLASHES,
    HYDROGEN_STABLE_BURNING,
    HYDROGEN_OPT_THICK_WINDS,
    HELIUM_WHITE_DWARF_HELIUM_SUB_CHANDRASEKHAR,
    HELIUM_WHITE_DWARF_HELIUM_IGNITION,
    HELIUM_WHITE_DWARF_HYDROGEN_FLASHES,
    HELIUM_WHITE_DWARF_HYDROGEN_ACCUMULATION
};
const COMPASUnorderedMap<ACCRETION_REGIME, std::string> ACCRETION_REGIME_LABEL = {
    { ACCRETION_REGIME::NONE,                                        "No accretion" },
    { ACCRETION_REGIME::ZERO,                                        "No accretion" },
    { ACCRETION_REGIME::HELIUM_ACCUMULATION,                         "Helium piles up without burning, full efficiency" },
    { ACCRETION_REGIME::HELIUM_FLASHES,                              "Helium ignites in flashes, partial accretion efficiency" },
    { ACCRETION_REGIME::HELIUM_STABLE_BURNING,                       "Helium is burnt without flashes, full efficiency" },
    { ACCRETION_REGIME::HELIUM_OPT_THICK_WINDS,                      "Helium is being accreted at a high rate, producing winds and limiting the accretion efficiency to a critical value" },
    { ACCRETION_REGIME::HYDROGEN_FLASHES,                            "Hydrogen ignites in flashes, partial accretion efficiency" },
    { ACCRETION_REGIME::HYDROGEN_STABLE_BURNING,                     "Hydrogen is burnt without flashes, full efficiency" },
    { ACCRETION_REGIME::HYDROGEN_OPT_THICK_WINDS,                    "Hydrogen is being accreted at a high rate, producing winds and limiting the accretion efficiency to a critical value" },
    { ACCRETION_REGIME::HELIUM_WHITE_DWARF_HELIUM_SUB_CHANDRASEKHAR, "Full accretion leads to transient, but it would not make enough radioactive Ni-56 to be classified as a SN Ia" },
    { ACCRETION_REGIME::HELIUM_WHITE_DWARF_HELIUM_IGNITION,          "Full accretion until material is ignited in a flash and degeneracy is lifted" },
    { ACCRETION_REGIME::HELIUM_WHITE_DWARF_HYDROGEN_FLASHES,         "Unstable hydrogen flashes lead to net accretion being zero" },
    { ACCRETION_REGIME::HELIUM_WHITE_DWARF_HYDROGEN_ACCUMULATION,    "Material piles up. Depending on the companion, it could lead to a CE episode or merger" }
};

// options to add program option columns to [BSE/SSE] SYSPARMS file
enum class ADD_OPTIONS_TO_SYSPARMS: int { ALWAYS, GRID, NEVER };
const COMPASUnorderedMap<ADD_OPTIONS_TO_SYSPARMS, std::string> ADD_OPTIONS_TO_SYSPARMS_LABEL = {
    { ADD_OPTIONS_TO_SYSPARMS::ALWAYS, "ALWAYS" },
    { ADD_OPTIONS_TO_SYSPARMS::GRID,   "GRID" },
    { ADD_OPTIONS_TO_SYSPARMS::NEVER,  "NEVER" }
};

// black hole kick options
enum class BLACK_HOLE_KICKS_MODE: int { FULL, REDUCED, ZERO, FALLBACK };
const COMPASUnorderedMap<BLACK_HOLE_KICKS_MODE, std::string> BLACK_HOLE_KICKS_MODE_LABEL = {
    { BLACK_HOLE_KICKS_MODE::FULL,     "FULL" },     // FULL kicks
    { BLACK_HOLE_KICKS_MODE::REDUCED,  "REDUCED" },  // REDUCED kicks
    { BLACK_HOLE_KICKS_MODE::ZERO,     "ZERO" },     // ZERO kicks
    { BLACK_HOLE_KICKS_MODE::FALLBACK, "FALLBACK" }  // FALLBACK kicks
};

// boost map update options for program options
// see options code for use
enum class BOOST_MAP: int { UPDATE, NO_UPDATE };

// kick magnitude distributions from Bray & Eldridge 2016,2018
enum class BRAY_ELDRIDGE_CONSTANT: int { ALPHA, BETA };
const COMPASUnorderedMap<BRAY_ELDRIDGE_CONSTANT, double> BRAY_ELDRIDGE_CONSTANT_VALUES = {
    { BRAY_ELDRIDGE_CONSTANT::ALPHA, 100.0 },
    { BRAY_ELDRIDGE_CONSTANT::BETA, -170.0 }
};

// case BB mass transfer stability prescriptions
enum class CASE_BB_STABILITY_PRESCRIPTION: int{ ALWAYS_STABLE, ALWAYS_STABLE_ONTO_NSBH, TREAT_AS_OTHER_MT, ALWAYS_UNSTABLE };
const COMPASUnorderedMap<CASE_BB_STABILITY_PRESCRIPTION, std::string> CASE_BB_STABILITY_PRESCRIPTION_LABEL = {
    { CASE_BB_STABILITY_PRESCRIPTION::ALWAYS_STABLE,           "ALWAYS_STABLE" },
    { CASE_BB_STABILITY_PRESCRIPTION::ALWAYS_STABLE_ONTO_NSBH, "ALWAYS_STABLE_ONTO_NSBH" },
    { CASE_BB_STABILITY_PRESCRIPTION::TREAT_AS_OTHER_MT,       "TREAT_AS_OTHER_MT" },
    { CASE_BB_STABILITY_PRESCRIPTION::ALWAYS_UNSTABLE,         "ALWAYS_UNSTABLE" }
};

// common envelope cccretion prescriptions
enum class CE_ACCRETION_PRESCRIPTION: int { ZERO, CONSTANT, UNIFORM, MACLEOD, CHEVALIER };
const COMPASUnorderedMap<CE_ACCRETION_PRESCRIPTION, std::string> CE_ACCRETION_PRESCRIPTION_LABEL = {
    { CE_ACCRETION_PRESCRIPTION::ZERO,      "ZERO" },
    { CE_ACCRETION_PRESCRIPTION::CONSTANT,  "CONSTANT" },
    { CE_ACCRETION_PRESCRIPTION::UNIFORM,   "UNIFORM" },
    { CE_ACCRETION_PRESCRIPTION::MACLEOD,   "MACLEOD" },
    { CE_ACCRETION_PRESCRIPTION::CHEVALIER, "CHEVALIER" }
};
    
// common envelope formalisms
enum class CE_FORMALISM: int { ENERGY, TWO_STAGE };
const COMPASUnorderedMap<CE_FORMALISM, std::string> CE_FORMALISM_LABEL = {
        { CE_FORMALISM::ENERGY,    "ENERGY" },
        { CE_FORMALISM::TWO_STAGE, "TWO_STAGE" }
};

// common envelope lambda prescriptions
enum class CE_LAMBDA_PRESCRIPTION: int { FIXED, LOVERIDGE, NANJING, KRUCKOW, DEWI };
const COMPASUnorderedMap<CE_LAMBDA_PRESCRIPTION, std::string> CE_LAMBDA_PRESCRIPTION_LABEL = {
    { CE_LAMBDA_PRESCRIPTION::FIXED,     "LAMBDA_FIXED" },
    { CE_LAMBDA_PRESCRIPTION::LOVERIDGE, "LAMBDA_LOVERIDGE" },
    { CE_LAMBDA_PRESCRIPTION::NANJING,   "LAMBDA_NANJING" },
    { CE_LAMBDA_PRESCRIPTION::KRUCKOW,   "LAMBDA_KRUCKOW" },
    { CE_LAMBDA_PRESCRIPTION::DEWI,      "LAMBDA_DEWI" }
};   

// CHE (Chemically Homogeneous Evolution) Options
enum class CHE_MODE: int { NONE, OPTIMISTIC, PESSIMISTIC };
const COMPASUnorderedMap<CHE_MODE, std::string> CHE_MODE_LABEL = {
    { CHE_MODE::NONE,        "NONE" },
    { CHE_MODE::OPTIMISTIC,  "OPTIMISTIC" },
    { CHE_MODE::PESSIMISTIC, "PESSIMISTIC" }
};

// main sequence core mass prescription
enum class CORE_MASS_PRESCRIPTION: int { NONE, MANDEL, SHIKAUCHI };
const COMPASUnorderedMap<CORE_MASS_PRESCRIPTION, std::string> CORE_MASS_PRESCRIPTION_LABEL = {
    { CORE_MASS_PRESCRIPTION::NONE,      "NONE" },
    { CORE_MASS_PRESCRIPTION::MANDEL,    "MANDEL" },
    { CORE_MASS_PRESCRIPTION::SHIKAUCHI, "SHIKAUCHI" }
};

// logfile delimiters
enum class DELIMITER: int { TAB, SPACE, COMMA };
const COMPASUnorderedMap<DELIMITER, std::string> DELIMITERLabel = {         // labels
    { DELIMITER::TAB,   "TAB" },
    { DELIMITER::SPACE, "SPACE" },
    { DELIMITER::COMMA, "COMMA" }
};
const COMPASUnorderedMap<DELIMITER, std::string> DELIMITERValue = {         // values
    { DELIMITER::TAB,   "\t" },
    { DELIMITER::SPACE, " " },
    { DELIMITER::COMMA, "," }
};

// eccentricity distributions
enum class ECCENTRICITY_DISTRIBUTION: int { ZERO, FLAT, THERMAL, GELLER_2013, DUQUENNOYMAYOR1991, SANA2012 };
const COMPASUnorderedMap<ECCENTRICITY_DISTRIBUTION, std::string> ECCENTRICITY_DISTRIBUTION_LABEL = {
    { ECCENTRICITY_DISTRIBUTION::ZERO,               "ZERO" },
    { ECCENTRICITY_DISTRIBUTION::FLAT,               "FLAT" },
    { ECCENTRICITY_DISTRIBUTION::THERMAL,            "THERMAL" },
    { ECCENTRICITY_DISTRIBUTION::GELLER_2013,        "GELLER+2013" },
    { ECCENTRICITY_DISTRIBUTION::DUQUENNOYMAYOR1991, "DUQUENNOYMAYOR1991" },
    { ECCENTRICITY_DISTRIBUTION::SANA2012,           "SANA2012"}
};

// envelope types
enum class ENVELOPE: int { RADIATIVE, CONVECTIVE, REMNANT };
const COMPASUnorderedMap<ENVELOPE, std::string> ENVELOPE_LABEL = {
    { ENVELOPE::RADIATIVE,  "RADIATIVE" },
    { ENVELOPE::CONVECTIVE, "CONVECTIVE" },
    { ENVELOPE::REMNANT,    "REMNANT" }
};

// envelope state prescriptions
enum class ENVELOPE_STATE_PRESCRIPTION: int { LEGACY, HURLEY, FIXED_TEMPERATURE };
const COMPASUnorderedMap<ENVELOPE_STATE_PRESCRIPTION, std::string> ENVELOPE_STATE_PRESCRIPTION_LABEL = {
    { ENVELOPE_STATE_PRESCRIPTION::LEGACY,            "LEGACY" },
    { ENVELOPE_STATE_PRESCRIPTION::HURLEY,            "HURLEY" },
    { ENVELOPE_STATE_PRESCRIPTION::FIXED_TEMPERATURE, "FIXED_TEMPERATURE" }
};

// evolution status constants
enum class EVOLUTION_STATUS: int {
    CONTINUE,
    DONE,
    ERROR,
    TIMES_UP,
    STEPS_UP,
    NO_TIMESTEPS_READ,
    TIMESTEPS_EXHAUSTED,
    TIMESTEPS_NOT_CONSUMED,
    SSE_ERROR,
    BINARY_ERROR,
    DCO_MERGER_TIME,
    STARS_TOUCHING,
    STELLAR_MERGER,
    STELLAR_MERGER_AT_BIRTH,
    DCO,
    WD_WD,
    MASSLESS_REMNANT,
    UNBOUND,
    NOT_STARTED,
    STARTED
};
// These descritions are deliberately succinct (as much as possible) so running status doesn't scroll off the page...
const COMPASUnorderedMap<EVOLUTION_STATUS, std::string> EVOLUTION_STATUS_LABEL = {
    { EVOLUTION_STATUS::CONTINUE,                "Continue evolution" },
    { EVOLUTION_STATUS::DONE,                    "Simulation completed" },
    { EVOLUTION_STATUS::ERROR,                   "Evolution stopped because an error occurred" },
    { EVOLUTION_STATUS::TIMES_UP,                "Allowed time exceeded" },
    { EVOLUTION_STATUS::STEPS_UP,                "Allowed timesteps exceeded" },
    { EVOLUTION_STATUS::NO_TIMESTEPS_READ,       "No user-provided timesteps read" },
    { EVOLUTION_STATUS::TIMESTEPS_EXHAUSTED,     "User-provided timesteps exhausted" },
    { EVOLUTION_STATUS::TIMESTEPS_NOT_CONSUMED,  "User-provided timesteps not consumed" },
    { EVOLUTION_STATUS::SSE_ERROR,               "SSE error for one of the constituent stars" },
    { EVOLUTION_STATUS::BINARY_ERROR,            "Error evolving binary" },
    { EVOLUTION_STATUS::DCO_MERGER_TIME,         "Time exceeded DCO merger (formation + coalescence) time" },
    { EVOLUTION_STATUS::STARS_TOUCHING,          "Stars touching" },
    { EVOLUTION_STATUS::STELLAR_MERGER,          "Stars merged" },
    { EVOLUTION_STATUS::STELLAR_MERGER_AT_BIRTH, "Stars merged at birth" },
    { EVOLUTION_STATUS::DCO,                     "DCO formed" },
    { EVOLUTION_STATUS::WD_WD,                   "Double White Dwarf formed" },
    { EVOLUTION_STATUS::MASSLESS_REMNANT,        "Massless Remnant formed" },
    { EVOLUTION_STATUS::UNBOUND,                 "Unbound binary" },
    { EVOLUTION_STATUS::NOT_STARTED,             "Simulation not started" },
    { EVOLUTION_STATUS::STARTED,                 "Simulation started" }
};

// evolution mode (SSE or BSE)
enum class EVOLUTION_MODE: int { SSE, BSE };
const COMPASUnorderedMap<EVOLUTION_MODE, std::string> EVOLUTION_MODE_LABEL = {
    { EVOLUTION_MODE::SSE, "SSE" },
    { EVOLUTION_MODE::BSE, "BSE" }
};

// floating-point error handling mode
// OFF   specifies that no floating-point error checking is performed
// ON    specifies that the current star/binary will be terminted if a floating-point error occurs
// DEBUG specifies that a stack trace will be printed and the program halted if a floating-point error occurs
enum class FP_ERROR_MODE: int { OFF, ON, DEBUG };
const COMPASUnorderedMap<FP_ERROR_MODE, std::string> FP_ERROR_MODE_LABEL = {
    { FP_ERROR_MODE::OFF,   "OFF" },
    { FP_ERROR_MODE::ON,    "ON" },
    { FP_ERROR_MODE::DEBUG, "DEBUG" }
};

// symbolic names for the Gamma Constants
// these must be left as default values - their order can be changed with the caveat that the sentinel "COUNT" must stay at the end
// it's a bit of a hack, but it lets me calculate the number of GAMMA_CONSTANTS
enum class GAMMA_CONSTANTS: int { B_GAMMA, C_GAMMA, COUNT };

// symbolic names for Giant Branch Parameters
// these must be left as default values - their order can be changed with the caveat that the sentinel "COUNT" must stay at the end
// it's a bit of a hack, but it lets me calculate the number of Timescales
enum class GBP: int {
    AH,                     // Hydrogen rate constant.  Hurley et al. 2000, p553
    AHHe,                   // Effective combined rate constant for both hydrogen and helium shell burning.  Hurley et al. 2000, eq 71
    AHe,                    // Helium rate constant.  Hurley et al. 2000, eq 68
    B,                      // Hurley et al. 2000, p552, eq38 (does this represent something physical?  If so, what?  How should this be described?)
    D,                      // Hurley et al. 2000, p552, eq38 (does this represent something physical?  If so, what?  How should this be described?)
    p,                      // Hurley et al. 2000, p552, eq38 (does this represent something physical?  If so, what?  How should this be described?)
    q,                      // Hurley et al. 2000, p552, eq38 (does this represent something physical?  If so, what?  How should this be described?)
    Lx,                     // Luminosity parameter on the first giant branch (FGB) Lx as a function of the core mass (really a function of Mx).
    Mx,                     // Crosover point of high-luminosity and low-luminosity in core mass - luminosity relation. Hurley et al. 2000, p552, eq38
    McBGB,                  // Core mass at BGB (Base of Giant Branch)
    McBAGB,                 // Core mass at BAGB (Base of Asymptotic Giant Branch).  Hurley et al. 2000, eq 66 (also see eq 75 and discussion)
    McDU,                   // Core mass at second dredge up.  Hurley et al. 2000, eq 69
    McSN,                   // Core mass at which the Asymptotic Giant Branch phase is terminated in a SN/loss of envelope

    COUNT                   // Sentinel for entry count
};

// initial mass functions
enum class INITIAL_MASS_FUNCTION: int { SALPETER, POWERLAW, UNIFORM, KROUPA };
const COMPASUnorderedMap<INITIAL_MASS_FUNCTION, std::string> INITIAL_MASS_FUNCTION_LABEL = {
    { INITIAL_MASS_FUNCTION::SALPETER, "SALPETER" },
    { INITIAL_MASS_FUNCTION::POWERLAW, "POWERLAW" },
    { INITIAL_MASS_FUNCTION::UNIFORM,  "UNIFORM" },
    { INITIAL_MASS_FUNCTION::KROUPA,   "KROUPA" }
};

// kick magnitude distributions
enum class KICK_MAGNITUDE_DISTRIBUTION: int { ZERO, FIXED, FLAT, MAXWELLIAN, BRAYELDRIDGE, MULLER2016, MULLER2016MAXWELLIAN, MULLERMANDEL};
const COMPASUnorderedMap<KICK_MAGNITUDE_DISTRIBUTION, std::string> KICK_MAGNITUDE_DISTRIBUTION_LABEL = {
    { KICK_MAGNITUDE_DISTRIBUTION::ZERO,                 "ZERO" },
    { KICK_MAGNITUDE_DISTRIBUTION::FIXED,                "FIXED" },
    { KICK_MAGNITUDE_DISTRIBUTION::FLAT,                 "FLAT" },
    { KICK_MAGNITUDE_DISTRIBUTION::MAXWELLIAN,           "MAXWELLIAN" },
    { KICK_MAGNITUDE_DISTRIBUTION::BRAYELDRIDGE,         "BRAYELDRIDGE" },
    { KICK_MAGNITUDE_DISTRIBUTION::MULLER2016,           "MULLER2016" },
    { KICK_MAGNITUDE_DISTRIBUTION::MULLER2016MAXWELLIAN, "MULLER2016MAXWELLIAN" },
    { KICK_MAGNITUDE_DISTRIBUTION::MULLERMANDEL,         "MULLERMANDEL" }
};

// kick direction distributions
enum class KICK_DIRECTION_DISTRIBUTION: int { ISOTROPIC, INPLANE, PERPENDICULAR, POWERLAW, WEDGE, POLES };
const COMPASUnorderedMap<KICK_DIRECTION_DISTRIBUTION, std::string> KICK_DIRECTION_DISTRIBUTION_LABEL = {
    { KICK_DIRECTION_DISTRIBUTION::ISOTROPIC,     "ISOTROPIC" },
    { KICK_DIRECTION_DISTRIBUTION::INPLANE,       "INPLANE" },
    { KICK_DIRECTION_DISTRIBUTION::PERPENDICULAR, "PERPENDICULAR" },
    { KICK_DIRECTION_DISTRIBUTION::POWERLAW,      "POWERLAW" },
    { KICK_DIRECTION_DISTRIBUTION::WEDGE,         "WEDGE" },
    { KICK_DIRECTION_DISTRIBUTION::POLES,         "POLES" }
};

// symbolic names for the Luminosity Constants
// these must be left as default values - their order can be changed with the caveat that the sentinel "COUNT" must stay at the end
// it's a bit of a hack, but it lets me calculate the number of L_CONSTANTS
enum class L_CONSTANTS: int { B_ALPHA_L, B_BETA_L, B_DELTA_L, COUNT };

// LBV mass loss prescriptions
enum class LBV_MASS_LOSS_PRESCRIPTION: int { NONE, ZERO, HURLEY_ADD, HURLEY, BELCZYNSKI };
const COMPASUnorderedMap<LBV_MASS_LOSS_PRESCRIPTION, std::string> LBV_MASS_LOSS_PRESCRIPTION_LABEL = {
    { LBV_MASS_LOSS_PRESCRIPTION::NONE,       "NONE" },     // DEPRECATED June 2024 - remove end 2024
    { LBV_MASS_LOSS_PRESCRIPTION::ZERO,       "ZERO" },
    { LBV_MASS_LOSS_PRESCRIPTION::HURLEY_ADD, "HURLEY_ADD" },
    { LBV_MASS_LOSS_PRESCRIPTION::HURLEY,     "HURLEY" },
    { LBV_MASS_LOSS_PRESCRIPTION::BELCZYNSKI, "BELCZYNSKI" }
};

// symbolic names for GB groups described in Loveridge et al., 2011
// these are used as indices into the loveridgeCoefficients multi-dimensional vector (described below)
enum class LOVERIDGE_GROUP: int { LMR1, LMR2, LMA, HM, RECOM };
const COMPASUnorderedMap<LOVERIDGE_GROUP, std::string> LOVERIDGE_GROUP_LABEL = {
    { LOVERIDGE_GROUP::LMR1,  "Low mass early Red Giant Branch (RGB) (before dredge-up)" },
    { LOVERIDGE_GROUP::LMR2,  "Low mass late Red Giant Branch (RGB) (after dredge-up)" },
    { LOVERIDGE_GROUP::LMA,   "Low mass Asymptotic Giant Branch (AGB)" },
    { LOVERIDGE_GROUP::HM,    "High mass" },
    { LOVERIDGE_GROUP::RECOM, "Recombination energy" }
};

// symbolic names for mass cutoffs
// these must be left as default values - their order can be changed with the caveat that the sentinel "COUNT" must stay at the end
// it's a bit of a hack, but it lets me calculate the number of Timescales
enum class MASS_CUTOFF: int {
    MHook,                  // Mass above which hook appears on MS (in Msol)
    MHeF,                   // Maximum initial mass for which helium ignites degenerately in a Helium Flash (HeF)
    MFGB,                   // Maximum initial mass for which helium ignites on the First Giant Branch (FGB)
    MCHE,                   // Mass cutoff for calculation of initial angular frequency to determine if CHE occurs

    COUNT                   // Sentinel for entry count
};

// mass loss prescriptions
enum class MASS_LOSS_PRESCRIPTION: int { NONE, ZERO, HURLEY, BELCZYNSKI2010, MERRITT2024 };
const COMPASUnorderedMap<MASS_LOSS_PRESCRIPTION, std::string> MASS_LOSS_PRESCRIPTION_LABEL = {
    { MASS_LOSS_PRESCRIPTION::NONE,           "NONE" },     // DEPRECATED June 2024 - remove end 2024
    { MASS_LOSS_PRESCRIPTION::ZERO,           "ZERO" },
    { MASS_LOSS_PRESCRIPTION::HURLEY,         "HURLEY" },
    { MASS_LOSS_PRESCRIPTION::BELCZYNSKI2010, "BELCZYNSKI2010" },
    { MASS_LOSS_PRESCRIPTION::MERRITT2024,   "MERRITT2024" }
};

// symbolic names for mass loss rate type
enum class MASS_LOSS_TYPE: int { NONE, GB, LBV, OB, RSG, VMS, WR};
const COMPASUnorderedMap<MASS_LOSS_TYPE, std::string> MASS_LOSS_TYPE_LABEL = {
    { MASS_LOSS_TYPE::NONE, "NONE" },
    { MASS_LOSS_TYPE::GB,   "GB" },
    { MASS_LOSS_TYPE::LBV,  "LBV" },
    { MASS_LOSS_TYPE::OB,   "OB" },
    { MASS_LOSS_TYPE::RSG,  "RSG" },
    { MASS_LOSS_TYPE::VMS,  "VMS" },
    { MASS_LOSS_TYPE::WR,   "WR" }
};

// mass ratio distributions
enum class MASS_RATIO_DISTRIBUTION: int { FLAT, DUQUENNOYMAYOR1991, SANA2012 };
const COMPASUnorderedMap<MASS_RATIO_DISTRIBUTION, std::string> MASS_RATIO_DISTRIBUTION_LABEL = {
    { MASS_RATIO_DISTRIBUTION::FLAT,               "FLAT" },
    { MASS_RATIO_DISTRIBUTION::DUQUENNOYMAYOR1991, "DUQUENNOYMAYOR1991" },
    { MASS_RATIO_DISTRIBUTION::SANA2012,           "SANA2012" }
};

// mass transfer timescale types
enum class MASS_TRANSFER_TIMESCALE: int { NONE, NUCLEAR, THERMAL, CE };
const COMPASUnorderedMap<MASS_TRANSFER_TIMESCALE, std::string> MASS_TRANSFER_TIMESCALE_LABEL = {
    { MASS_TRANSFER_TIMESCALE::NONE,                "NONE" },
    { MASS_TRANSFER_TIMESCALE::NUCLEAR,             "NUCLEAR" },
    { MASS_TRANSFER_TIMESCALE::THERMAL,             "THERMAL" },
    { MASS_TRANSFER_TIMESCALE::CE,                  "CE" }
};

// metallicity distributions
enum class METALLICITY_DISTRIBUTION: int { ZSOLAR, LOGUNIFORM };
const COMPASUnorderedMap<METALLICITY_DISTRIBUTION, std::string> METALLICITY_DISTRIBUTION_LABEL = {
    { METALLICITY_DISTRIBUTION::ZSOLAR,     "ZSOLAR" },
    { METALLICITY_DISTRIBUTION::LOGUNIFORM, "LOGUNIFORM" }
};

// mass transfer accretion efficiency prescriptions
enum class MT_ACCRETION_EFFICIENCY_PRESCRIPTION: int { THERMALLY_LIMITED, FIXED_FRACTION };
const COMPASUnorderedMap<MT_ACCRETION_EFFICIENCY_PRESCRIPTION, std::string> MT_ACCRETION_EFFICIENCY_PRESCRIPTION_LABEL = {
    { MT_ACCRETION_EFFICIENCY_PRESCRIPTION::THERMALLY_LIMITED, "THERMAL" },
    { MT_ACCRETION_EFFICIENCY_PRESCRIPTION::FIXED_FRACTION,    "FIXED" }
};

// mass transfer angular momentum loss prescriptions
enum class MT_ANGULAR_MOMENTUM_LOSS_PRESCRIPTION: int { JEANS, ISOTROPIC_RE_EMISSION, CIRCUMBINARY_RING, MACLEOD_LINEAR, ARBITRARY };
const COMPASUnorderedMap<MT_ANGULAR_MOMENTUM_LOSS_PRESCRIPTION, std::string> MT_ANGULAR_MOMENTUM_LOSS_PRESCRIPTION_LABEL = {
    { MT_ANGULAR_MOMENTUM_LOSS_PRESCRIPTION::JEANS,                 "JEANS" },
    { MT_ANGULAR_MOMENTUM_LOSS_PRESCRIPTION::ISOTROPIC_RE_EMISSION, "ISOTROPIC" },
    { MT_ANGULAR_MOMENTUM_LOSS_PRESCRIPTION::CIRCUMBINARY_RING,     "CIRCUMBINARY" },
    { MT_ANGULAR_MOMENTUM_LOSS_PRESCRIPTION::MACLEOD_LINEAR,        "MACLEOD_LINEAR" },
    { MT_ANGULAR_MOMENTUM_LOSS_PRESCRIPTION::ARBITRARY,             "ARBITRARY" }
};

// mass transfer cases
enum class MT_CASE: int { NONE, A, B, C, OTHER };
const COMPASUnorderedMap<MT_CASE, std::string> MT_CASE_LABEL = {
    { MT_CASE::NONE, "Mass Transfer CASE NONE: No Mass Transfer" },
    { MT_CASE::A,    "Mass Transfer CASE A" },                          // mass transfer while donor is on main sequence
    { MT_CASE::B,    "Mass Transfer CASE B" },                          // donor star is in (or evolving to) Red Giant phase
    { MT_CASE::C,    "Mass Transfer CASE C" },                          // SuperGiant phase
    { MT_CASE::OTHER,"Mass Transfer CASE OTHER: Multiple MT events" }   // default value, or multiple MT events
};

// mass transfer rejuvenation prescriptions
enum class MT_REJUVENATION_PRESCRIPTION: int { HURLEY, STARTRACK };
const COMPASUnorderedMap<MT_REJUVENATION_PRESCRIPTION, std::string> MT_REJUVENATION_PRESCRIPTION_LABEL = {
    { MT_REJUVENATION_PRESCRIPTION::HURLEY,     "HURLEY" },
    { MT_REJUVENATION_PRESCRIPTION::STARTRACK,  "STARTRACK" }
};

// mass transfer thermally limited variation options
enum class MT_THERMALLY_LIMITED_VARIATION: int { C_FACTOR, RADIUS_TO_ROCHELOBE };
const COMPASUnorderedMap<MT_THERMALLY_LIMITED_VARIATION, std::string> MT_THERMALLY_LIMITED_VARIATION_LABEL = {
    { MT_THERMALLY_LIMITED_VARIATION::C_FACTOR,            "CFACTOR" },
    { MT_THERMALLY_LIMITED_VARIATION::RADIUS_TO_ROCHELOBE, "ROCHELOBE" }
};

// mass transfer timing options for writing to BSE_RLOF file
enum class MT_TIMING: int { PRE_MT, POST_MT };
const COMPASUnorderedMap<MT_TIMING, std::string> MT_TIMING_LABEL = {
    { MT_TIMING::PRE_MT,  "PRE_MT" },
    { MT_TIMING::POST_MT, "POST_MT" }
};

// mass transfer tracking constants
enum class MT_TRACKING: int { NO_MASS_TRANSFER, STABLE_1_TO_2_SURV, STABLE_2_TO_1_SURV, CE_1_TO_2_SURV, CE_2_TO_1_SURV, CE_DOUBLE_SURV, MERGER }; 
const COMPASUnorderedMap<MT_TRACKING, std::string> MT_TRACKING_LABEL = {
    { MT_TRACKING::NO_MASS_TRANSFER,   "NO MASS TRANSFER" },
    { MT_TRACKING::STABLE_1_TO_2_SURV, "MASS TRANSFER STABLE STAR1 -> STAR2" },
    { MT_TRACKING::STABLE_2_TO_1_SURV, "MASS TRANSFER STABLE STAR2 -> STAR1" },
    { MT_TRACKING::CE_1_TO_2_SURV,     "MASS TRANSFER COMMON ENVELOPE STAR1 -> STAR2" },
    { MT_TRACKING::CE_2_TO_1_SURV,     "MASS TRANSFER COMMON ENVELOPE STAR2 -> STAR1" },
    { MT_TRACKING::CE_DOUBLE_SURV,     "MASS TRANSFER COMMON ENVELOPE DOUBLE CORE" },
    { MT_TRACKING::MERGER,             "MASS TRANSFER -> MERGER" }
};

// neutrino mass loss BH formation prescriptions
enum class NEUTRINO_MASS_LOSS_PRESCRIPTION: int { FIXED_FRACTION, FIXED_MASS };
const COMPASUnorderedMap<NEUTRINO_MASS_LOSS_PRESCRIPTION, std::string> NEUTRINO_MASS_LOSS_PRESCRIPTION_LABEL = {
    { NEUTRINO_MASS_LOSS_PRESCRIPTION::FIXED_FRACTION, "FIXED_FRACTION" },
    { NEUTRINO_MASS_LOSS_PRESCRIPTION::FIXED_MASS,     "FIXED_MASS" }
};

// neutron star equations of state
enum class NS_EOS: int { SSE, ARP3 };
const COMPASUnorderedMap<NS_EOS, std::string> NS_EOSLabel = {
    { NS_EOS::SSE,  "SSE" },
    { NS_EOS::ARP3, "ARP3" }
};

// OB (main sequence) mass loss prescriptions
enum class OB_MASS_LOSS_PRESCRIPTION: int { NONE, ZERO, VINK2001, VINK2021, BJORKLUND2022, KRTICKA2018};
const COMPASUnorderedMap<OB_MASS_LOSS_PRESCRIPTION, std::string> OB_MASS_LOSS_PRESCRIPTION_LABEL = {
    { OB_MASS_LOSS_PRESCRIPTION::NONE,          "NONE" },       // DEPRECATED June 2024 - remove end 2024
    { OB_MASS_LOSS_PRESCRIPTION::ZERO,          "ZERO" },
    { OB_MASS_LOSS_PRESCRIPTION::VINK2001,      "VINK2001" },
    { OB_MASS_LOSS_PRESCRIPTION::VINK2021,      "VINK2021" },
    { OB_MASS_LOSS_PRESCRIPTION::BJORKLUND2022, "BJORKLUND2022" },
    { OB_MASS_LOSS_PRESCRIPTION::KRTICKA2018,   "KRTICKA2018" }
};

// object persistence
// specifies whether an object is permanent or ephemaeral
// EPHEMERAL is typically used for clones so that they don't participate in logging, etc.
enum class OBJECT_PERSISTENCE: int { PERMANENT, EPHEMERAL };
const COMPASUnorderedMap<OBJECT_PERSISTENCE, std::string> OBJECT_PERSISTENCE_LABEL = {
    { OBJECT_PERSISTENCE::PERMANENT, "Permanent" },
    { OBJECT_PERSISTENCE::EPHEMERAL, "Ephemeral" }
};

// object types
// identifies the type of an object
// if BASE_STAR, check STELLAR_TYPE    
enum class OBJECT_TYPE: int { NONE, MAIN, PROFILING, UTILS, STAR, BASE_STAR, BINARY_STAR, BASE_BINARY_STAR, BINARY_CONSTITUENT_STAR };
const COMPASUnorderedMap<OBJECT_TYPE, std::string> OBJECT_TYPE_LABEL = {
    { OBJECT_TYPE::NONE,                    "Not_an_Object!" },
    { OBJECT_TYPE::MAIN,                    "Main" },
    { OBJECT_TYPE::PROFILING,               "Profiling" },
    { OBJECT_TYPE::UTILS,                   "Utils" },
    { OBJECT_TYPE::STAR,                    "Star" },
    { OBJECT_TYPE::BASE_STAR,               "BaseStar" },
    { OBJECT_TYPE::BINARY_STAR,             "BinaryStar" },
    { OBJECT_TYPE::BASE_BINARY_STAR,        "BaseBinaryStar" },
    { OBJECT_TYPE::BINARY_CONSTITUENT_STAR, "BinaryConstituentStar" }
};

// program options origin indicator (command line or gridfile line)
enum class OPTIONS_ORIGIN: int { CMDLINE, GRIDFILE };

// orbital period distributions
enum class ORBITAL_PERIOD_DISTRIBUTION: int { FLATINLOG };
const COMPASUnorderedMap<ORBITAL_PERIOD_DISTRIBUTION, std::string> ORBITAL_PERIOD_DISTRIBUTION_LABEL = {
    { ORBITAL_PERIOD_DISTRIBUTION::FLATINLOG, "FLATINLOG" },
};

// pulsational pair instability prescriptions
enum class PPI_PRESCRIPTION: int { COMPAS, STARTRACK, MARCHANT, FARMER, HENDRIKS };
const COMPASUnorderedMap<PPI_PRESCRIPTION, std::string> PPI_PRESCRIPTION_LABEL = {
    { PPI_PRESCRIPTION::COMPAS,    "COMPAS" },
    { PPI_PRESCRIPTION::STARTRACK, "STARTRACK" },
    { PPI_PRESCRIPTION::MARCHANT,  "MARCHANT" },
    { PPI_PRESCRIPTION::FARMER,    "FARMER" },
    { PPI_PRESCRIPTION::HENDRIKS,  "HENDRIKS" } 
};

// program status
enum class PROGRAM_STATUS: int { SUCCESS, CONTINUE, STOPPED, ERROR_IN_COMMAND_LINE, LOGGING_FAILED, ERROR_UNHANDLED_EXCEPTION };

// pulsar birth magnetic field distributions
enum class PULSAR_BIRTH_MAGNETIC_FIELD_DISTRIBUTION: int { ZERO, FLATINLOG, UNIFORM, LOGNORMAL };
const COMPASUnorderedMap<PULSAR_BIRTH_MAGNETIC_FIELD_DISTRIBUTION, std::string> PULSAR_BIRTH_MAGNETIC_FIELD_DISTRIBUTION_LABEL = {
    { PULSAR_BIRTH_MAGNETIC_FIELD_DISTRIBUTION::ZERO,      "ZERO" },
    { PULSAR_BIRTH_MAGNETIC_FIELD_DISTRIBUTION::FLATINLOG, "FLATINLOG" },
    { PULSAR_BIRTH_MAGNETIC_FIELD_DISTRIBUTION::UNIFORM,   "UNIFORM" },
    { PULSAR_BIRTH_MAGNETIC_FIELD_DISTRIBUTION::LOGNORMAL, "LOGNORMAL" }
};

// pulsar birth spin period distributions
enum class PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION: int { ZERO, UNIFORM, NORMAL };
const COMPASUnorderedMap<PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION, std::string> PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION_LABEL = {
    { PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION::ZERO,    "ZERO" },
    { PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION::UNIFORM, "UNIFORM" },
    { PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION::NORMAL,  "NORMAL" }
};

// symbolic names for the Radius Constants
// these must be left as default values - their order can be changed with the caveat that the sentinel "COUNT" must stay at the end
// it's a bit of a hack, but it lets me calculate the number of R_CONSTANTS
enum class R_CONSTANTS: int { B_ALPHA_R, C_ALPHA_R, B_BETA_R, C_BETA_R, B_DELTA_R, COUNT };

// remnant mass prescriptions
enum class REMNANT_MASS_PRESCRIPTION: int { HURLEY2000, BELCZYNSKI2002, FRYER2012, FRYER2022, MULLER2016, MULLERMANDEL, SCHNEIDER2020, SCHNEIDER2020ALT, MALTSEV2024};
const COMPASUnorderedMap<REMNANT_MASS_PRESCRIPTION, std::string> REMNANT_MASS_PRESCRIPTION_LABEL = {
    { REMNANT_MASS_PRESCRIPTION::HURLEY2000,       "HURLEY2000" },
    { REMNANT_MASS_PRESCRIPTION::BELCZYNSKI2002,   "BELCZYNSKI2002" },
    { REMNANT_MASS_PRESCRIPTION::FRYER2012,        "FRYER2012" },
    { REMNANT_MASS_PRESCRIPTION::FRYER2022,        "FRYER2022" },
    { REMNANT_MASS_PRESCRIPTION::MULLER2016,       "MULLER2016" },
    { REMNANT_MASS_PRESCRIPTION::MULLERMANDEL,     "MULLERMANDEL" },
    { REMNANT_MASS_PRESCRIPTION::SCHNEIDER2020,    "SCHNEIDER2020" },
    { REMNANT_MASS_PRESCRIPTION::SCHNEIDER2020ALT, "SCHNEIDER2020ALT" },
    { REMNANT_MASS_PRESCRIPTION::MALTSEV2024,      "MALTSEV2024" }
};

// rotational velocity distributions
enum class ROTATIONAL_VELOCITY_DISTRIBUTION: int { ZERO, HURLEY, VLTFLAMES };
const COMPASUnorderedMap<ROTATIONAL_VELOCITY_DISTRIBUTION, std::string> ROTATIONAL_VELOCITY_DISTRIBUTION_LABEL = {
    { ROTATIONAL_VELOCITY_DISTRIBUTION::ZERO,      "ZERO" },
    { ROTATIONAL_VELOCITY_DISTRIBUTION::HURLEY,    "HURLEY" },
    { ROTATIONAL_VELOCITY_DISTRIBUTION::VLTFLAMES, "VLTFLAMES" }
};

// semi-major axis distributions
enum class SEMI_MAJOR_AXIS_DISTRIBUTION: int { FLATINLOG, DUQUENNOYMAYOR1991, SANA2012 };
const COMPASUnorderedMap<SEMI_MAJOR_AXIS_DISTRIBUTION, std::string> SEMI_MAJOR_AXIS_DISTRIBUTION_LABEL = {
    { SEMI_MAJOR_AXIS_DISTRIBUTION::FLATINLOG,          "FLATINLOG" },
    { SEMI_MAJOR_AXIS_DISTRIBUTION::DUQUENNOYMAYOR1991, "DUQUENNOYMAYOR1991" },
    { SEMI_MAJOR_AXIS_DISTRIBUTION::SANA2012,           "SANA2012" }
};

// Fryer 2012 supernova engines
enum class SN_ENGINE: int { RAPID, DELAYED };
const COMPASUnorderedMap<SN_ENGINE, std::string> SN_ENGINE_LABEL = {
    { SN_ENGINE::RAPID,   "RAPID" },
    { SN_ENGINE::DELAYED, "DELAYED" }
};

// supernova events/states
//
// The values here for SN_EVENT are powers of 2 so that they can be used in a bit map
// and manipulated with bit-wise logical operators
//
// Ordinarily we might expect that an SN event could be only one of
//
//    NONE, CCSN, ECSN, PISN, PPISN, USSN, AIC, SNIA, or HeSD
//
// Note that the CCSN value here replaces the SN value in the legacy code
// The legacy code implemented these values as boolean flags, and the SN flag was always set when
// the USSN flag was set (but not the converse).  In the legacy code when the ECSN flag was set 
// the SN flag was not set.  In the legacy code the PISN and PPISN flags were used to track history
// and we only set for the "experienced" condition (I think).
//
// To match the legacy code usage of these flags, here the "is" and "experienced" conditions 
// ("current" and "past" SN events) are implemented as bit maps - different values can be
// ORed or ANDed into the bit map (that way the USSN and CCSN flags can be set at the same
// time - necessary for the code flow (from the legacy code) - which we should probably one
// day look at and rewrite).
//
// HeSD stands for helium-shell detonation
//
// A convenience function has been provided in utils.cpp to interpret the bit map (utils::SNEventType()).
// Given an SN_EVENT bitmap (current or past), it returns (in priority order):
//     
//    SN_EVENT::NONE    iff no bits are set
//    SN_EVENT::CCSN    iff CCSN  bit is set and USSN bit is not set
//    SN_EVENT::ECSN    iff ECSN  bit is set
//    SN_EVENT::PISN    iff PISN  bit is set
//    SN_EVENT::PPISN   iff PPISN bit is set
//    SN_EVENT::USSN    iff USSN  bit is set
//    SN_EVENT::AIC     iff AIC   bit is set
//    SN_EVENT::SNIA    iff SNIA  bit is set and HeSD bit is not set
//    SN_EVENT::HeSD    iff HeSD  bit is set
//    SN_EVENT::UNKNOWN otherwise
//
enum class SN_EVENT: int { 
    NONE         = 0, 
    CCSN         = 1, 
    ECSN         = 2, 
    PISN         = 4, 
    PPISN        = 8, 
    USSN         = 16,
    AIC          = 32,
    SNIA         = 64,
    HeSD         = 128,
    UNKNOWN      = 32768 // doesn't really matter what this is because the value is never used (as long as it's > sum of all the others)
};
const COMPASUnorderedMap<SN_EVENT, std::string> SN_EVENT_LABEL = {
    { SN_EVENT::NONE,    "No Supernova" },
    { SN_EVENT::CCSN,    "Core Collapse Supernova" },
    { SN_EVENT::ECSN,    "Electron Capture Supernova" },
    { SN_EVENT::PISN,    "Pair Instability Supernova" },
    { SN_EVENT::PPISN,   "Pulsational Pair Instability Supernova" },
    { SN_EVENT::USSN,    "Ultra Stripped Supernova" },
    { SN_EVENT::AIC,     "Accretion-Induced Collapse" }, 
    { SN_EVENT::SNIA,    "Supernova Type Ia" }, 
    { SN_EVENT::HeSD,    "Helium-shell detonation" }, 
    { SN_EVENT::UNKNOWN, "Unknown Supernova Type" }
};
ENABLE_BITMASK_OPERATORS(SN_EVENT);

// supernova states
enum class SN_STATE: int { NONE, STAR1, STAR2, BOTH };
const COMPASUnorderedMap<SN_STATE, std::string> SN_STATE_LABEL = {
    { SN_STATE::NONE,  "No Supernova" },
    { SN_STATE::STAR1, "Star1 only" },
    { SN_STATE::STAR2, "Star2 only" },
    { SN_STATE::BOTH,  "Both stars" }
};

// stellar populations
enum class STELLAR_POPULATION: int { POPULATION_I, POPULATION_II }; // JR FIX popI was 1, popII 0
const COMPASUnorderedMap<STELLAR_POPULATION, std::string> STELLAR_POPULATION_LABEL = {
    { STELLAR_POPULATION::POPULATION_I,  "POPULATION_I" },
    { STELLAR_POPULATION::POPULATION_II, "POPULATION_II" }
};

// tides prescriptions
enum class TIDES_PRESCRIPTION: int { NONE, PERFECT, KAPIL2024 };
const COMPASUnorderedMap<TIDES_PRESCRIPTION, std::string> TIDES_PRESCRIPTION_LABEL = {
    { TIDES_PRESCRIPTION::NONE,      "NONE" },
    { TIDES_PRESCRIPTION::PERFECT,   "PERFECT" },
    { TIDES_PRESCRIPTION::KAPIL2024, "KAPIL2024" }
};

// symbolic names for timescales
// these must be left as default values - their order can be changed with the caveat that the sentinel "COUNT" must stay at the end
// it's a bit of a hack, but it lets me calculate the number of Timescales
enum class TIMESCALE: int {
    tMS,                    // Main sequence
    tBGB,                   // Base of Giant Branch
    tHeI,                   // Helium ignition
    tHe,                    // Helium burning

                            // First Giant Branch (FGB)
    tinf1_FGB,              // First Giant Branch tinf1
    tinf2_FGB,              // First Giant Branch tinf2
    tMx_FGB,                // First Giant Branch t(Mx)

                            // Early Asymptotic Giant Branch (EAGB) (FAGB in Hurley's sse)
    tinf1_FAGB,             // Early Asymptotic Giant Branch tinf1
    tinf2_FAGB,             // Early Asymptotic Giant Branch tinf2
    tMx_FAGB,               // Early Asymptotic Giant Branch t(Mx)

                            // Thermally Pulsating Asymptotic Giant Branch (TPAGB) (SAGB in Hurley's sse)
    tinf1_SAGB,             // Thermally Pulsating Asymptotic Giant Branch tinf1
    tinf2_SAGB,             // Thermally Pulsating Asymptotic Giant Branch tinf2
    tMx_SAGB,               // Thermally Pulsating Asymptotic Giant Branch t(Mx)
    tP,                     // (tDU?)
                            // Helium Giant Branch
    tHeMS,                  // Naked Helium Star central helium burning lifetime (HeMs)
    tinf1_HeGB,             // Helium Giant Branch tinf1
    tinf2_HeGB,             // Helium Giant Branch tinf2
    tx_HeGB,                // Helium Giant Branch tx (is this t(Mx)?)
    tau_BL,                 // Relative duration of blue loop taubl
    tauX_BL,                // Relative start of blue loop taux
    tauY_BL,                // Relative end of blue loop tauy

    COUNT                   // Sentinel for entry count
};

// critical mass ratio prescriptions
enum class QCRIT_PRESCRIPTION: int { NONE, CLAEYS, GE, GE_IC, HURLEY_HJELLMING_WEBBINK};
const COMPASUnorderedMap<QCRIT_PRESCRIPTION, std::string> QCRIT_PRESCRIPTION_LABEL = {
    { QCRIT_PRESCRIPTION::NONE,    "NONE" },
    { QCRIT_PRESCRIPTION::CLAEYS,  "CLAEYS" },
    { QCRIT_PRESCRIPTION::GE_IC, "GE_IC" },
    { QCRIT_PRESCRIPTION::GE,    "GE" },
    { QCRIT_PRESCRIPTION::HURLEY_HJELLMING_WEBBINK, "HURLEY_HJELLMING_WEBBINK" },
};

// RSG mass loss prescriptions
enum class RSG_MASS_LOSS_PRESCRIPTION: int { NONE, ZERO, VINKSABHAHIT2023, BEASOR2020, DECIN2023, YANG2023, KEE2021, NJ90};
const COMPASUnorderedMap<RSG_MASS_LOSS_PRESCRIPTION, std::string> RSG_MASS_LOSS_PRESCRIPTION_LABEL = {
    { RSG_MASS_LOSS_PRESCRIPTION::NONE,             "NONE" },   // DEPRECATED June 2024 - remove end 2024
    { RSG_MASS_LOSS_PRESCRIPTION::ZERO,             "ZERO" },
    { RSG_MASS_LOSS_PRESCRIPTION::VINKSABHAHIT2023, "VINKSABHAHIT2023" },
    { RSG_MASS_LOSS_PRESCRIPTION::BEASOR2020,       "BEASOR2020" },
    { RSG_MASS_LOSS_PRESCRIPTION::DECIN2023,        "DECIN2023" },
    { RSG_MASS_LOSS_PRESCRIPTION::YANG2023,         "YANG2023" },
    { RSG_MASS_LOSS_PRESCRIPTION::KEE2021,          "KEE2021" },
    { RSG_MASS_LOSS_PRESCRIPTION::NJ90,             "NJ90" }
};

// VMS (very massive stars) mass loss prescriptions
enum class VMS_MASS_LOSS_PRESCRIPTION: int { NONE, ZERO, VINK2011, BESTENLEHNER2020, SABHAHIT2023};
const COMPASUnorderedMap<VMS_MASS_LOSS_PRESCRIPTION, std::string> VMS_MASS_LOSS_PRESCRIPTION_LABEL = {
    { VMS_MASS_LOSS_PRESCRIPTION::NONE,             "NONE" },   // DEPRECATED June 2024 - remove end 2024
    { VMS_MASS_LOSS_PRESCRIPTION::ZERO,             "ZERO" },
    { VMS_MASS_LOSS_PRESCRIPTION::VINK2011,         "VINK2011" },
    { VMS_MASS_LOSS_PRESCRIPTION::BESTENLEHNER2020, "BESTENLEHNER2020" },
    { VMS_MASS_LOSS_PRESCRIPTION::SABHAHIT2023,     "SABHAHIT2023" }
};

// WR mass loss prescriptions
enum class WR_MASS_LOSS_PRESCRIPTION: int { BELCZYNSKI2010, SANDERVINK2023, SHENAR2019 };
const COMPASUnorderedMap<WR_MASS_LOSS_PRESCRIPTION, std::string> WR_MASS_LOSS_PRESCRIPTION_LABEL = {
    { WR_MASS_LOSS_PRESCRIPTION::BELCZYNSKI2010, "BELCZYNSKI2010" },
    { WR_MASS_LOSS_PRESCRIPTION::SANDERVINK2023, "SANDERVINK2023" },
    { WR_MASS_LOSS_PRESCRIPTION::SHENAR2019,     "SHENAR2019" }
};

// common envelope zeta prescriptions
enum class ZETA_PRESCRIPTION: int { SOBERMAN, HURLEY, ARBITRARY };
const COMPASUnorderedMap<ZETA_PRESCRIPTION, std::string> ZETA_PRESCRIPTION_LABEL = {
    { ZETA_PRESCRIPTION::SOBERMAN,  "SOBERMAN" },
    { ZETA_PRESCRIPTION::HURLEY,    "HURLEY" },
    { ZETA_PRESCRIPTION::ARBITRARY, "ARBITRARY" }
};



// boost variant definition for allowed data types
// used for variable specification to define logfile records
typedef boost::variant<
    bool,
    short int,
    int,
    long int,
    long long int,
    unsigned short int,
    unsigned int,
    unsigned long int,
    unsigned long long int,
    float,
    double,
    long double,
    std::string,
    std::vector<std::string>,
    ERROR,
    STELLAR_TYPE,
    MT_CASE,
    MT_TRACKING,
    MASS_TRANSFER_TIMESCALE,
    SN_EVENT,
    SN_STATE,
    EVOLUTION_STATUS
> COMPAS_VARIABLE;




// common type definitions
typedef std::initializer_list<SN_EVENT> SN_EVENT_LIST;
typedef std::vector<STELLAR_TYPE>       ST_VECTOR;
typedef std::vector<COMPAS_VARIABLE>    COMPAS_VARIABLE_VECTOR;



// Option details
typedef struct OptionDetails {
    std::string optionStr;                                  // name string
    std::string valueStr;                                   // value string
    std::string sourceStr;                                  // source string (COMPAS DEFAULT or USER SUPPLIED)    
    std::string typeStr;                                    // detailed data type (e.g. UNSIGNED LONG INT)
    TYPENAME    dataType;                                   // short data type (e.g. INT, FLOAT, etc)
    std::string defaultStr;                                 // default value string
    STR_VECTOR  allowedStr;                                 // vector of allowed value strings (empty for options that don't have multiple allowed values)
} OptionDetailsT;


// Log file details
typedef struct LogfileDetails {
    int                           id;                       // logfile id
    std::string                   filename;                 // filename
    int                           recordTypes;              // bitmap of record types to be written to the file (if 0, the file is disabled)
    ANY_PROPERTY_VECTOR           recordProperties;         // list of properties (columns) to be written to the logfile
    std::vector<TYPENAME>         propertyTypes;            // the COMPAS datatypes of the properties
    std::vector<STRING_QUALIFIER> stringTypes;              // the string type (fixed or variable length) for TYPENAME::STRING datatypes
    STR_VECTOR                    hdrStrings;               // the column header strings
    STR_VECTOR                    unitsStrings;             // the column units strings
    STR_VECTOR                    typeStrings;              // the column datatype strings
    STR_VECTOR                    fmtStrings;               // format strings for the columns - how the value is formatted for printing to the logfile
    BOOL_VECTOR                   annotations;              // print flags for each annotation specified by the user (e.g. OPTIONS->NotesHdrs())
} LogfileDetailsT;


// Grid file details
typedef struct Gridfile {
    std::string   filename;                                 // filename for grid file
    ERROR         error;                                    // status - ERROR::NONE if no problem, otherwise an error number
    std::ifstream handle;                                   // the file handle

    std::streamsize startLine;                              // the first line of the grid file to process (0-based)
    std::streamsize currentLine;                            // the grid line currently being processed
    std::streamsize linesProcessed;                         // the number of grid lines processed so far in this run
    std::streamsize linesToProcess;                         // the number of grid lines to process (from start line)
} GridfileT;


// RotationalVelocityParams struct for gsl root solver
struct RotationalVelocityParams {                           // Structure containing parameter (u) for the root solving function using gsl_root_solver
    double u;                                               // Value of CDF, draw in U(0,1)
};


// KickMagnitudeParams struct for gsl root solver
struct KickMagnitudeParams {
    double y;       // Value of CDF, should be drawn as U(0,1)
    double sigma;   // sigma for kick distribution
};


// struct for supernova events:
// CCSN, ECSN, PISN, PPSIN, USSN, AIC

typedef struct SNEvents {
    SN_EVENT current;                                       // Supernova event at the current timestep: NONE if no supernova event happening
    SN_EVENT past;                                          // Supernova event at any past timestep   : NONE if no supernova event happened in any past timestep
} SNEventsT;


// supernova kick struct for both SSE and BSE options
//
// some of these are only required for binary stars, but
// easier (and more logical I think (for now, anyway) to
// keep all SN-related attributes in the same place
//
// we need to know if these values were actually specified
// by the user via options - hence the boolean values

typedef struct KickParameters {
    bool   magnitudeRandomSpecified;                        // SSE and BSE
    double magnitudeRandom;                                 // SSE and BSE

    bool   magnitudeSpecified;                              // SSE and BSE
    double magnitude;                                       // SSE and BSE

    bool   phiSpecified;                                    // BSE only
    double phi;                                             // BSE only

    bool   thetaSpecified;                                  // BSE only
    double theta;                                           // BSE only

    bool   meanAnomalySpecified;                            // BSE only
    double meanAnomaly;                                     // BSE only
} KickParameters;


// struct for supernova attributes of the base star
// some of these are only required for binary stars, but
// easier (and more logical I think (for now, anyway) to
// keep all SN-related attributes in the same place

typedef struct SupernovaDetails {                           // Holds attributes, flags - if the star went supernova

    KickParameters initialKickParameters;                   // User-supplied initial kick parameters - if present used in place of drawing randomly/from distributions
    
    double         coreMassAtCOFormation;                   // Core mass of this star when it formed a compact object
    double         COCoreMassAtCOFormation;                 // Carbon Oxygen core mass of the star when it goes supernova and forms a compact object
    double         drawnKickMagnitude;                      // Kick magnitude the system received during the supernova (km s^-1)
    double         eccentricAnomaly;                        // Eccentric anomaly at instataneous time of the SN
    SNEventsT      events;                                  // Record of supernova events undergone by the star
    double         fallbackFraction;                        // Fallback fraction during a supernova event
    double         HeCoreMassAtCOFormation;                 // Helium core mass of the star when it goes supernova and forms a compact object
    bool           isHydrogenPoor;                          // Flag to indicate if exploding star is hydrogen-poor. We consider an H-rich star all SN progenitors that have an H envelope, otherwise H-poor
    double         kickMagnitude;                           // Kick magnitude the system received during the supernova (km s^-1)
    double         kickMagnitudeRandom;                     // Random number U(0,1) for choosing the supernova kick magnitude - drawn once at star creation
    double         rocketKickMagnitude;                     // Rocket kick magnitude the system received after the supernova (km s^-1)
    double         rocketKickPhi;                           // Rocket kick azimuthal angle phi the system received after the supernova 
    double         rocketKickTheta;                         // Rocket kick polar angle theta the system received after the supernova 
    double         meanAnomaly;                             // Mean anomaly at instantaneous time of the SN - uniform in [0, 2pi]
    double         phi;                                     // Kick angle in the orbital plane, defined CCW from the radial vector pointed away from the Companion (rad) [0, 2pi)
    SN_STATE       supernovaState;                          // Indicates which star (or stars) are undergoing / have undergone a supernova event
    double         theta;                                   // Kick angle out of the orbital plane, toward the orbital angular momentum axis (rad) [-pi/2, pi/2]
    double         totalMassAtCOFormation;                  // Total mass of the star when it goes supernova and forms a compact object
    double         trueAnomaly;                             // True anomaly at instantaneous time of the SN
} SupernovaDetailsT;


// pulsar parameters (if star becomes a Neutron Star)
typedef struct PulsarDetails {
    double magneticField;                                   // Pulsar magnetic field strength (G)
    double spinPeriod;                                      // Pulsar spin period (ms)
    double spinFrequency;                                   // Pulsar spin frequency in rads per second
    double spinDownRate;                                    // Pulsar spin down rate as time derivative of spin frequency (fdot, rad s^-2)
    double birthPeriod;                                     // Pulsar birth period (s)
    double birthSpinDownRate;                               // Pulsar birth down rate as Pdot (s s^-1)
} PulsarDetailsT;


// struct for Lambdas
typedef struct Lambdas {
	double dewi;                                            // JR: todo: description?
    double fixed;                                           // Set to OPTIONS->commonEnvelopeLambda
	double kruckow;                                         // Calculated using m_Radius and OPTIONS->commonEnvelopeSlopeKruckow
	double kruckowBottom;                                   // Calculated using m_Radius and -1
	double kruckowMiddle;                                   // Ccalculated using m_Radius and -4/5
	double kruckowTop;                                      // Calculated using m_Radius and -2/3
	double loveridge;                                       // No mass loss
	double loveridgeWinds;                                  // Mass loss
	double nanjing;                                         // JR: todo: description?
} LambdasT;


// struct for Zetas
// JR: add descriptive comments
typedef struct Zetas {                                      // JR: todo: descriptions for these?
	double hurley;
	double hurleyHe;
	double nuclear;
	double soberman;
	double sobermanHe;
	double thermal;
} ZetasT;


// struct for binding energies
typedef struct BindingEnergies {
    double fixed;                                           // Calculated using lambda = OPTIONS->commonEnvelopeLambda
	double nanjing;                                         // Calculated using lambda = m_Lambdas.nanjing
	double loveridge;                                       // Calculated using lambda = m_Lambdas.loveridge
	double loveridgeWinds;                                  // Calculated using lambda = m_Lambdas.loveridgeWinds
	double kruckow;                                         // Calculated using lambda = m_Lambdas.kruckow
    double dewi;                                            // Calculated using lambda = m_Lambdas.dewi
} BindingEnergiesT;


// RLOF properties
// JR: add descriptive comments
typedef struct RLOFProperties {
    OBJECT_ID    id;

    STELLAR_TYPE stellarType1;
    STELLAR_TYPE stellarType2;

    double       mass1;
    double       mass2;

    double       radius1;
    double       radius2;

    double       starToRocheLobeRadiusRatio1;                                    
    double       starToRocheLobeRadiusRatio2;

    double       eccentricity;
    double       semiMajorAxis;

    unsigned int eventCounter;

    double       time;
    double       timePrev;

    bool         isRLOF1;
    bool         isRLOF2;

    bool         isCE;
    
    double       massLossRateFromDonor;
    double       accretionEfficiency;
    MASS_TRANSFER_TIMESCALE massTransferTimescale;

} RLOFPropertiesT;

// JR: add descriptive comments
typedef struct BinaryRLOFDetails {                          // RLOF details pertinent to binaries

    bool             experiencedRLOF;
    bool             immediateRLOFPostCEE;                  // Here for now - maybe should be in Binary CEDetails struct?       JR: todo:
    bool             isRLOF;
    bool             simultaneousRLOF;                      // Here for now - maybe should be in Binary CEDetails struct?       JR: todo:
    bool             stableRLOFPostCEE;                     // Here for now - maybe should be in Binary CEDetails struct?       JR: todo:
    RLOFPropertiesT  props1;
    RLOFPropertiesT  props2;
    RLOFPropertiesT* propsPreMT;
    RLOFPropertiesT* propsPostMT;
} BinaryRLOFDetailsT;

// JR: add descriptive comments
typedef struct StellarRLOFDetails {                         // RLOF details pertinent to individual stars
    bool isRLOF;
    bool experiencedRLOF;
    bool RLOFPostCEE;
} StellarRLOFDetailsT;

// Common Envelope properties
// JR: add descriptive comments
typedef struct BinaryCEESavedValues {
    double eccentricity;
   	double rocheLobe1to2;
	double rocheLobe2to1;
    double semiMajorAxis;
} BinaryCEESavedValuesT;

// JR: add descriptive comments
typedef struct BinaryCEDetails {                            // Common Envelope details pertinent to binaries
    BinaryCEESavedValuesT preCEE;
    BinaryCEESavedValuesT postCEE;

    bool                  CEEnow;                           // Indicates whether a common envelope event is occurring now
    unsigned int          CEEcount;                         // Common Envelope Event count
    bool                  doubleCoreCE;
    bool                  optimisticCE;
} BinaryCEDetailsT;


// JR: add descriptive comments
typedef struct StellarCEESavedValues {
    double       bindingEnergy;
    double       dynamicalTimescale;
    double       luminosity;
    double       mass;
    double       radialExpansionTimescale;
    double       radius;
    STELLAR_TYPE stellarType;
    double       temperature;
    double       thermalTimescale;
} StellarCEESavedValuesT;

// JR: add descriptive comments
typedef struct StellarCEDetails {                           // Common Envelope details pertinent to individual stars
    StellarCEESavedValuesT preCEE;
    StellarCEESavedValuesT postCEE;

    double                 bindingEnergy;
    double                 COCoreMass;
    double                 CoreMass;
    double                 HeCoreMass;
    double                 lambda;
    double                 convectiveEnvelopeMass;          // for two-stage CE formalism
    double                 radiativeIntershellMass;         // for two-stage CE formalism
    double                 convectiveEnvelopeBindingEnergy; // for two-stage CE formalism
} StellarCEDetailsT; // was CommonEnvelopeDetailsT;


// For boost ODE integrators, see https://www.boost.org/doc/libs/1_83_0/libs/numeric/odeint/doc/html/boost_numeric_odeint/tutorial/harmonic_oscillator.html
typedef DBL_VECTOR state_type;
typedef boost::numeric::odeint::runge_kutta_cash_karp54<state_type> error_stepper_type;
typedef boost::numeric::odeint::controlled_runge_kutta<error_stepper_type> controlled_stepper_type;

#endif // __typedefs_h__
