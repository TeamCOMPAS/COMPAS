#ifndef __constants_h__
#define __constants_h__

#include <cmath>
#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <fstream>

// the defaults size of the boost list that handles variant types is 20 - so only 20 variant types are allowed
// we've exceeded that number - we're at 21 currently - so the size of the boost list needs to be increased
// we have to set the size of the list before we include the boost headers - otherwise boost redefines it
// boost only recognises values of 20, 30, 40, & 50: for now we set it to 30. 
#define BOOST_MPL_CFG_NO_PREPROCESSED_HEADERS
  
#if !defined(BOOST_MPL_LIMIT_LIST_SIZE)
    #if defined(BOOST_MPL_LIST_HPP_INCLUDED)
      # error "BOOST_MPL_LIMIT_LIST_SIZE must be set to accommodate the size of COMPAS_VARIABLE_TYPE before the including the Boost headers"
    #endif

    #define BOOST_MPL_LIMIT_LIST_SIZE 30

#elif BOOST_MPL_LIMIT_LIST_SIZE < 30
    #if defined(BOOST_MPL_LIST_HPP_INCLUDED)
      # error "BOOST_MPL_LIMIT_LIST_SIZE must be set to accommodate the size of COMPAS_VARIABLE_TYPE before the including the Boost headers"
    #else
      # error "BOOST_MPL_LIMIT_LIST_SIZE value is too low"
    #endif
    
#endif

#include <boost/variant.hpp>


typedef unsigned long int                                               OBJECT_ID;                  // OBJECT_ID type

typedef std::vector<double>                                             DBL_VECTOR;
typedef std::tuple <double, double>                                     DBL_DBL;
typedef std::tuple <double, double, double>                             DBL_DBL_DBL;
typedef std::tuple<std::string, std::string, std::string, std::string>  STR_STR_STR_STR;

// Hash for Enum Class
struct EnumClassHash
{
    template <typename T>
    std::size_t operator()(T t) const
    {
        return static_cast<std::size_t>(t);
    }
};

template <typename Key>
using HashType = typename std::conditional<std::is_enum<Key>::value, EnumClassHash, std::hash<Key>>::type;

template <typename Key, typename T>
using COMPASUnorderedMap = std::unordered_map<Key, T, HashType<Key>>;



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


/*
 * Trick to allow SWITCH on literal strings
 * constexpr is (or can be) evaluated at compile-time, hence the ability to SWITCH using it.
 *
 * This function returns a hash value for a string - very small chance of collisions, but if a
 * collision happens the compiler will complain (because CASE values will be the same).
 *
 * Adapted from https://rextester.com/discussion/FGBUG88403/Switch-case-with-strings-in-C-11,
 * but variations available in many places with the right Google search terms...
 */

constexpr uint64_t StrHash(char const* p_Str, uint64_t p_Hash = 14695981039346656037ull) {
    return (*p_Str == 0) ? p_Hash : StrHash(p_Str + 1, (p_Hash * 1099511628211ull) ^ static_cast<uint64_t>(*p_Str));
}

constexpr uint64_t _(char const* p_Str) {
    return StrHash(p_Str);
}



extern OBJECT_ID globalObjectId;                                                                    // used to uniquely identify objects - used primarily for error printing

// Constants in SI
// CPLB: Use CODATA values where applicable http://physics.nist.gov/cuu/Constants/index.html

// JR: cmath (included file) provides the following pi-related constants:
//
// pi	        M_PI	    3.14159265358979323846
// pi/2         M_PI_2	    1.57079632679489661923
// pi/4         M_PI_4	    0.785398163397448309616
// 1/pi         M_1_PI	    0.318309886183790671538
// 2/pi         M_2_PI	    0.636619772367581343076
// 2/sqrt(pi)	M_2_SQRTPI	1.12837916709551257390
//
// I've added _2_PI and SQRT_M_2_PI below

#undef COMPARE_GLOBAL_TOLERANCE // define/undef this to compare floats with/without tolerance (see FLOAT_TOLERANCE_ABSOLUTE, FLOAT_TOLERANCE_RELATIVE and Compare() function)

constexpr double FLOAT_TOLERANCE_ABSOLUTE               = 0.0000005;                                                // absolute tolerance for floating-point comparisons if COMPARE_GLOBAL_TOLERANCE is defined
constexpr double FLOAT_TOLERANCE_RELATIVE               = 0.0000005;                                                // relative tolerance for floating-point comparisons if COMPARE_GLOBAL_TOLERANCE is defined

constexpr double ROOT_ABS_TOLERANCE                     = 1.0E-10;                                                  // absolute tolerance for root finder
constexpr double ROOT_REL_TOLERANCE                     = 1.0E-6;                                                   // relative tolerance for root finder


// initialisation constants
constexpr double DEFAULT_INITIAL_DOUBLE_VALUE           = 0.0;                                                      // default initial value for double variables
constexpr double DEFAULT_INITIAL_INTEGER_VALUE          = 0;                                                        // default initial value for int variables
constexpr double DEFAULT_INITIAL_ULONGINT_VALUE         = 0l;                                                       // default initial value for unsigned long int variables
constexpr double DEFAULT_INITIAL_BOOLEAN_VALUE          = false;                                                    // default initial value for bool variables


// conversion constants

// mass
constexpr double G_TO_KG                                = 1.0E-3;                                                   // convert grams to kg
constexpr double MSOL_TO_G                              = 1.98892E33;                                               // convert Solar Mass to g
constexpr double MSOL_TO_KG                             = MSOL_TO_G * G_TO_KG;                                      // convert Solar Mass to kg
constexpr double KG_TO_MSOL                             = 1.0 / MSOL_TO_KG;                                         // convert kg to Solar Mass

// length
constexpr double KM_TO_CM 					            = 1.0E5;									                // convert km to cm
constexpr double KM_TO_M                                = 1000.0;                                                   // convert km to m
constexpr double CM_TO_M                                = 1.0E-2;                                                   // convert cm to m

constexpr double RSOL_TO_KM                             = 6.957E5;                                                  // convert Solar Radius (RSOL) to km
constexpr double RSOL_TO_CM                             = 6.957E10;                                                 // convert Solar Radius (RSOL) to cm
constexpr double RSOL_TO_AU                             = 0.00465047;                                               // convert Solar Radius (RSOL) to AU

constexpr double AU_TO_CM                               = 14959787070000.0;                                         // convert Astronomical Units (AU) to cm
constexpr double AU_TO_RSOL				                = 1.0 / RSOL_TO_AU;                                         // convert Astronomical Units AU to Solar Radius RSOL
constexpr double AU_TO_KM                               = AU_TO_CM / 1.0E5;                                         // convert Astronomical Units AU to km

constexpr double KM_TO_RSOL					            = 1.0 / RSOL_TO_KM;						                    // convert km to Solar Radius (RSOL)
constexpr double KM_TO_AU                               = 1.0 / AU_TO_KM;                                           // convert km to Astronomical Units AU

// time
constexpr double DAYS_IN_QUAD                           = 1461.0;                                                   // number of days in any 4-year period
constexpr double DAYS_IN_YEAR                           = DAYS_IN_QUAD / 4.0;                                       // mean days per year, given DAYS_IN_QUAD
constexpr double SECONDS_IN_DAY                         = 24.0 * 60.0 * 60.0;                                       // number of seconds in 1 day
constexpr double SECONDS_IN_QUAD                        = DAYS_IN_QUAD * SECONDS_IN_DAY;                            // number of seconds in 1 quad
constexpr double SECONDS_IN_YEAR                        = SECONDS_IN_QUAD / 4.0;                                    // number of seconds in 1 year
constexpr double SECONDS_IN_MYR                         = SECONDS_IN_YEAR * 1.0E6;                                  // number of seconds in 1 Myr
constexpr double SECONDS_IN_MS                          = 1.0E-3;                                                   // number of seconds in 1 millisecond

constexpr double MYR_TO_YEAR                            = 1.0E6;                                                    // convert Myr to year
constexpr double YEAR_TO_MYR                            = 1.0E-6;                                                   // convert year to Myr

// energy
constexpr double JOULES_TO_ERG                          = 1.0E7;                                                    // convert Joules to Erg

// B field
constexpr double TESLA_TO_GAUSS                         = 1.0E4;					                                // convert Tesla to Gauss
constexpr double GAUSS_TO_TESLA                         = 1.0 / TESLA_TO_GAUSS;                                     // convert Gauss to Tesla

// systems
constexpr double CGS_SI                                 = G_TO_KG * CM_TO_M * CM_TO_M;                              // convert CGS to SI

// opacity
constexpr double OPACITY_CGS_TO_SI                      = 0.1;                                                      // cm^2 g^-1 to m^2 kg^-1

// constants

constexpr double _2_PI                                  = M_PI * 2.0;                                               // 2PI
constexpr double PI_2                                   = M_PI * M_PI;                                              // PI squared
constexpr double SQRT_M_2_PI                            = 0.79788456080286536;                                      // sqrt(2/PI)
constexpr double DEGREE                                 = M_PI / 180.0;                                             // 1 degree in radians

constexpr double GAMMA_E                                = 0.57721566490153286;                                      // Euler's Constant

constexpr double H0                                     = 67.8;                                                     // Hubble's Constant in km s^-1 Mpc^-1  (from plank approx 67.80±0.77) CPLB: Use WMAP value
constexpr double H0SI                                   = H0 * 1000.0 / 3.0E22;                                     // Hubble's Constant in SI units, s^-1
constexpr double HUBBLE_TIME                            = 1 / H0SI;                                                 // Hubble time in s

constexpr double G                                      = 6.67E-11;                                                 // Gravitational constant in m^3 kg^-1 s^-2 (more accurately known as G M_sol)
constexpr double G_CGS                                  = 6.6743E-8;                                                // Gravitational constant in cm^3 g^-1 s^-2
constexpr double G_AU_Msol_yr                           = 4.0 * PI_2;                                               // Gravitational constant in AU^3 Msol^-1 yr^-2
constexpr double G_km_Msol_s                            = G * 1.0E-9 / KG_TO_MSOL;                                  // Gravitational constant in km^3 Msol^-1 s^-2
constexpr double G_SOLAR_YEAR                           = 3.14E7;                                                   // Gravitational constant in Lsol Rsol yr Msol^-2 for calculating photon tiring limit

constexpr double RSOL                                   = 6.957E8;                                                  // Solar Radius (in m)
constexpr double ZSOL                                   = 0.02;                                                     // Solar Metallicity used in scalings
constexpr double LOG10_ZSOL                             = -1.69897;                                                 // log10(ZSOL) - for performance
constexpr double ZSOL_ASPLUND                           = 0.0142;                                                   // Solar Metallicity (Asplund+ 2010) used in initial condition
constexpr double TSOL                                   = 5778.0;                                                   // Solar Temperature in kelvin
constexpr double LSOL                                   = 3.844E33;                                                 // Solar Luminosity in erg/s
constexpr double LSOLW                                  = 3.844E26;                                                 // Solar luminosity (in W)

constexpr double AU                                     = 149597870700.0;                                           // 1 AU (Astronomical Unit) in metres
constexpr double KM                                     = 1000.0;                                                   // 1 km (Kilometre) in metres
constexpr double C                                      = 3.0E8;                                                    // Speed of light in m s^-1

constexpr double MU_0                                   = 4.0 * M_PI * 1.0E-7;                                      // Vacuum permeability in m kg s-2 A-2

constexpr double NEUTRINO_LOSS_FALLBACK_FACTOR          = 1.0;                                                      // Factor which accounts for mass loss in neutrino winds during a supernovae. Should be made a flag and added to pythonSubmit.py

constexpr double MC_L_C1                                = 9.20925E-5;                                               // Core Mass - Luminosity relation constant c1 (Hurley et al. 2000, eq 44)
constexpr double MC_L_C2                                = 5.402216;                                                 // Core Mass - Luminosity relation constant c2 (Hurley et al. 2000, eq 44)

constexpr double HE_RATE_CONSTANT                       = 7.66E-5;                                                  // Helium rate constant (Hurley et al. 2000, eq 68)
constexpr double HHE_RATE_CONSTANT                      = 1.27E-5;                                                  // Combined rate constant for both hydrogen and helium shell burning (Hurley et al. 2000, eq 71)

constexpr double BLACK_HOLE_LUMINOSITY                  = 1.0E-10;                                                  // Black Hole luminosity

constexpr double NEUTRON_STAR_MASS                      = 1.4;                                                      // Canonical NS mass in Msol
constexpr double NEUTRON_STAR_RADIUS                    = (1.0 / 7.0) * 1.0E-4;                                     // 10km in Rsol.  Hurley et al. 2000, just after eq 93

constexpr double MCH                                    = 1.44;                                                     // Chandrasekhar mass
constexpr double MECS                                   = 1.38;                                                     // Mass of Neutron-Star (NS) formed in electron capture supernova (ECS). From Belczysnki+2008, before eq. 3.
constexpr double MECS_REM                               = 1.26;                                                     // Gravitational mass of Neutron-Star (NS) formed in electron capture supernova (ECS). From Belczysnki+2008, eq. 3
constexpr double MASS_LOSS_ETA                          = 0.5;                                                      // Mass loss efficiency -- can be set in the code as an option easily enough
constexpr double MCBUR1HURLEY					        = 1.6;							                            // Minimum core mass at base of the AGB to avoid fully degenerate CO core formation (Hurley value, Fryer+ and Belczynski+ use 1.83)
constexpr double MCBUR2					                = 2.25;							                            // Core mass at base of the AGB above which the CO core is completely non-degenerate

constexpr double NJ_MINIMUM_LUMINOSITY                  = 4.0E3;                                                    // Minimum luminosity in Lsun needed for Nieuwenhuijzen & de Jager wind mass loss
constexpr double VINK_MASS_LOSS_MINIMUM_TEMP            = 1.25E4;                                                   // Minimum temperature in K for Vink mass loss rates to be applied
constexpr double VERY_MASSIVE_MINIMUM_MASS              = 100.0;                                                    // Minimum mass for applying Very Massive (VMS) mass rates to be applied
constexpr double RSG_MAXIMUM_TEMP                       = 8.0E3;                                                    // Upper temperature in K for Red Supergiant (RSG) mass loss rates to be applied
constexpr double VINK_MASS_LOSS_BISTABILITY_TEMP        = 2.5E4;                                                    // Temperature in K for bistability jump in Vink mass loss (assumed to be 25000K following Belczysnki+2010)
constexpr double VINK_MASS_LOSS_MAXIMUM_TEMP            = 5.0E4;                                                    // Maximum temperature in K for Vink mass loss rates to be applied (show warning above this)
constexpr double LBV_LUMINOSITY_LIMIT_STARTRACK         = 6.0E5;                                                    // STARTRACK LBV luminosity limit
constexpr double LBV_LUMINOSITY_LIMIT_VANBEVEREN        = 3.0E5;                                                    // VANBEVEREN LBV luminosity limit

constexpr double CONVECTIVE_BOUNDARY_TEMPERATURE_BELCZYNSKI = 5.37E3;                                               // Threshold temperature for the star to develop a convective envelope, in Kelvin (10^3.73 K, from Belczynski+, 2008)

constexpr double TIMESTEP_QUANTUM                       = 1.0E-12;                                                  // Timestep quantum in Myr (=31.5576 seconds, given DAYS_IN_QUAD)
constexpr double ABSOLUTE_MINIMUM_TIMESTEP              = 3.0 * TIMESTEP_QUANTUM;                                   // In Myr (=94.6728 seconds, given TIMESTEP QUANTUM)
constexpr double NUCLEAR_MINIMUM_TIMESTEP               = 1.0E6 * TIMESTEP_QUANTUM;                                 // Minimum time step for nuclear evolution in My (= 1 year = 31577600 seconds, given TIMESTEP_QUANTUM)

constexpr unsigned int ABSOLUTE_MAXIMUM_TIMESTEPS       = 1000000;                                                  // Absolute maximum number of timesteps

constexpr int    MAX_BSE_INITIAL_CONDITIONS_ITERATIONS  = 100;                                                      // Maximum loop iterations looking for initial conditions for binary systems
constexpr int    MAX_TIMESTEP_RETRIES                   = 30;                                                       // Maximum retries to find a good timestep for stellar evolution

constexpr double MAXIMUM_MASS_LOSS_FRACTION             = 0.01;                                                     // Maximum allowable mass loss - 1.0% (of mass) expressed as a fraction
constexpr double MAXIMUM_RADIAL_CHANGE                  = 0.01;                                                     // Maximum allowable radial change - 1% (of radius) expressed as a fraction
constexpr double MINIMUM_MASS_SECONDARY                 = 4.0;                                                      // Minimum mass of secondary to evolve

constexpr double MAXIMUM_MASS_TRANSFER_FRACTION_PER_STEP= 0.001;                                                    // Maximal fraction of donor mass that can be transferred in one step of stable mass transfer

constexpr double LAMBDA_NANJING_ZLIMIT                  = 0.0105;                                                   // Metallicity cutoff for Nanjing lambda calculations
constexpr double LAMBDA_NANJING_POPI_Z                  = 0.02;                                                     // Population I metallicity in Xu & Li (2010)
constexpr double LAMBDA_NANJING_POPII_Z                 = 0.001;                                                    // Population II metallicity in Xu & Li (2010)
constexpr double LAMBDA_NANJING_POPI_LOGZ               = -1.69897;                                                 // Population I log metallicity in Xu & Li (2010)
constexpr double LAMBDA_NANJING_POPII_LOGZ              = -3.0;                                                     // Population II log metallicity in Xu & Li (2010)
constexpr double LAMBDA_NANJING_MIN_MASS                = 1.0;                                                      // Minimum tabulated mass model in Xu & Li (2010)
constexpr double LAMBDA_NANJING_MAX_MASS                = 100.0;                                                    // Maximum tabulated mass model in Xu & Li (2010)

constexpr int    MAX_KEPLER_ITERATIONS                  = 1000;                                                     // Maximum number of iterations to solve Kepler's equation
constexpr double NEWTON_RAPHSON_EPSILON                 = 1.0E-5;                                                   // Accuracy for Newton-Raphson method

constexpr double EPSILON_PULSAR                         = 1.0;                                                      // JR: todo: description

constexpr double MIN_HMXRB_STAR_TO_ROCHE_LOBE_RADIUS_RATIO = 0.8;                                                   // Minimum value of stellar radius | Roche Lobe radius for visible HMXRBs

constexpr double ADAPTIVE_RLOF_FRACTION_DONOR_GUESS     = 0.001;                                                    // Fraction of donor mass to use as guess in BaseBinaryStar::MassLossToFitInsideRocheLobe()
constexpr int    ADAPTIVE_RLOF_MAX_TRIES                = 30;                                                       // Maximum number of tries in BaseBinaryStar::MassLossToFitInsideRocheLobe()
constexpr int    ADAPTIVE_RLOF_MAX_ITERATIONS           = 50;                                                       // Maximum number of root finder iterations in BaseBinaryStar::MassLossToFitInsideRocheLobe()
constexpr double ADAPTIVE_RLOF_SEARCH_FACTOR_FRAC       = 1.0;                                                      // Search size factor (fractional part) in BaseBinaryStar::MassLossToFitInsideRocheLobe() (added to 1.0)

constexpr int    ADAPTIVE_MASS0_MAX_TRIES               = 30;                                                       // Maximum number of tries in HG::Mass0ToMatchDesiredCoreMass()
constexpr int    ADAPTIVE_MASS0_MAX_ITERATIONS          = 50;                                                       // Maximum number of iterations in HG::Mass0ToMatchDesiredCoreMass()
constexpr double ADAPTIVE_MASS0_SEARCH_FACTOR_FRAC      = 1.0;                                                      // Search size factor (fractional part) in HG::Mass0ToMatchDesiredCoreMass() (added to 1.0)

constexpr int    TIDES_OMEGA_MAX_TRIES                  = 30;                                                       // Maximum number of tries in BaseBinaryStar::OmegaAfterCircularisation()
constexpr int    TIDES_OMEGA_MAX_ITERATIONS             = 50;                                                       // Maximum number of root finder iterations in BaseBinaryStar::OmegaAfterCircularisation()
constexpr double TIDES_OMEGA_SEARCH_FACTOR_FRAC         = 1.0;                                                      // Search size factor (fractional part) in BaseBinaryStar::OmegaAfterCircularisation() (added to 1.0)

constexpr double FARMER_PPISN_UPP_LIM_LIN_REGIME        = 38.0;                                                     // Maximum CO core mass to result in the linear remnant mass regime of the FARMER PPISN prescription
constexpr double FARMER_PPISN_UPP_LIM_QUAD_REGIME       = 60.0;                                                     // Maximum CO core mass to result in the quadratic remnant mass regime of the FARMER PPISN prescription
constexpr double FARMER_PPISN_UPP_LIM_INSTABILLITY      = 140.0;                                                    // Maximum CO core mass to result in PI (upper edge of PISN gap) from FARMER PPISN prescription
constexpr double STARTRACK_PPISN_HE_CORE_MASS           = 45.0;                                                     // Helium core mass remaining following PPISN as assumed in StarTrack (Belczynski et al. 2017 https://arxiv.org/abs/1607.03116)


// logging constants

enum class LOGFILETYPE: int { NONE, HDF5, CSV, TSV, TXT };                                                          // Need this declared here so can declare the constant...

const LOGFILETYPE DEFAULT_LOGFILE_TYPE                  = LOGFILETYPE::HDF5;                                        // Default logfile type
const std::string DEFAULT_OUTPUT_CONTAINER_NAME         = "COMPAS_Output";                                          // Default name for output container (directory)
const std::string DETAILED_OUTPUT_DIRECTORY_NAME        = "Detailed_Output";                                        // Name for detailed output directory within output container
const std::string RUN_DETAILS_FILE_NAME                 = "Run_Details";                                            // Name for run details output file within output container

constexpr int    HDF5_DEFAULT_CHUNK_SIZE                = 100000;                                                   // default HDF5 chunk size (number of dataset entries)
constexpr int    HDF5_DEFAULT_IO_BUFFER_SIZE            = 1;                                                        // number of HDF5 chunks to buffer for IO (per open dataset)
constexpr int    HDF5_MINIMUM_CHUNK_SIZE                = 1000;                                                     // minimum HDF5 chunk size (number of dataset entries)

// Logfile record types
// Note all enum classes for log record types start at 1 (and *must* start at 1)
typedef unsigned int LOGRECORDTYPE;

enum class BE_BINARY_RECORD_TYPE: unsigned int {                                                                    // BSE_BE_BINARIES file record type
    DEFAULT = 1                                                                                                     // 1 - default BSE_BE_BINARIES file record type
};

enum class CE_RECORD_TYPE: unsigned int {                                                                           // BSE_COMMON_ENVELOPES file record type
    DEFAULT = 1                                                                                                     // 1 - default BSE_COMMON_ENVELOPES file record type
};

enum class DCO_RECORD_TYPE: unsigned int {                                                                          // BSE_DOUBLE_COMPACT_OBJECTS file record type
    DEFAULT = 1                                                                                                     // 1 - default BSE_DOUBLE_COMPACT_OBJECTS file record type
};

enum class PULSAR_RECORD_TYPE: unsigned int {                                                                       // BSE_PULSAR_EVOLUTION file record type
    DEFAULT = 1,                                                                                                    //  1 - record describes the initial state of the binary
    POST_BINARY_TIMESTEP                                                                                            //  3 - record was logged immediately following binary timestep (i.e. the evolution of the binary system for a single timestep)
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
    POST_MASS_RESOLUTION_MERGER                                                                                     // 15 - record was logged immediately following a merger after mass resolution
};

enum class SSE_DETAILED_RECORD_TYPE: unsigned int {                                                                 // SSE_DETAILED_OUTPUT file record type
    INITIAL_STATE = 1,                                                                                              //  1 - record describes the initial state of the star
    PRE_MASS_LOSS,                                                                                                  //  2 - record was logged after timestep taken, but before mass loss resolution
    POST_MASS_LOSS,                                                                                                 //  3 - record was logged after after mass loss resolution
    TIMESTEP_COMPLETED,                                                                                             //  4 - record was logged immediately following the completion of the timestep (after all changes to the star)
    FINAL_STATE                                                                                                     //  5 - record describes the final state of the star
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

// option constraints
// Use these constant to specify constraints that should be applied to program option values
// The values specified here should be checked in Options::OptionValues::CheckAndSetOptions()
// and in any relevant sampling functions

constexpr double MINIMUM_INITIAL_MASS                   = 0.00007;                                                  // Minimum initial mass (Msol) (~theoretical minimum? How low does COMPAS actually tolerate?)
constexpr double MAXIMUM_INITIAL_MASS                   = 150.0;                                                    // Maximum initial mass (Msol) (should actually be 100Msol?)

constexpr double MINIMUM_METALLICITY                    = 0.0001;                                                   // Minimum metallicity - Hurley equations known to fail for Z < 0.0001
constexpr double MAXIMUM_METALLICITY                    = 0.03;                                                     // Maximum metallicity (~> super-metal-rich?)


// IMF constants
constexpr double SALPETER_POWER                         = -2.35;
constexpr double SALPETER_MINIMUM                       = 0.5;
constexpr double SALPETER_MAXIMUM                       = 100.0;

// Kroupa IMF is a broken power law with three slopes
constexpr double KROUPA_POWER_1                         = -0.3;
constexpr double KROUPA_POWER_2                         = -1.3;
constexpr double KROUPA_POWER_3                         = -2.3;

// Declare some values here so we don't need to repeatedly calculate them in the code

// Often require the power law exponent plus one
constexpr double KROUPA_POWER_PLUS1_1                   = 0.7;
constexpr double KROUPA_POWER_PLUS1_2                   = -0.3;
constexpr double KROUPA_POWER_PLUS1_3                   = -1.3;

constexpr double ONE_OVER_KROUPA_POWER_1_PLUS1          = 1.0 / KROUPA_POWER_PLUS1_1;
constexpr double ONE_OVER_KROUPA_POWER_2_PLUS1          = 1.0 / KROUPA_POWER_PLUS1_2;
constexpr double ONE_OVER_KROUPA_POWER_3_PLUS1          = 1.0 / KROUPA_POWER_PLUS1_3;

// There are two breaks in the Kroupa power law -- they occur here (in solar masses)
constexpr double KROUPA_BREAK_1                         = 0.08;
constexpr double KROUPA_BREAK_2                         = 0.5;

// Some values that are really constants
constexpr double KROUPA_BREAK_1_PLUS1_1                 = 0.17067228025785934;                                      // pow(KROUPA_BREAK_1, KROUPA_POWER_PLUS1_1);
constexpr double KROUPA_BREAK_1_PLUS1_2                 = 2.13340350322324179;                                      // pow(KROUPA_BREAK_1, KROUPA_POWER_PLUS1_2);
constexpr double KROUPA_BREAK_1_POWER_1_2               = 0.08;                                                     // pow(KROUPA_BREAK_1, (KROUPA_POWER_1 - KROUPA_POWER_2));

constexpr double KROUPA_BREAK_2_PLUS1_2                 = 1.23114441334491628;                                      // pow(KROUPA_BREAK_2, KROUPA_POWER_PLUS1_2);
constexpr double KROUPA_BREAK_2_PLUS1_3                 = 2.46228882668983257;                                      // pow(KROUPA_BREAK_2, KROUPA_POWER_PLUS1_3);
constexpr double KROUPA_BREAK_2_POWER_2_3               = 0.5;                                                      // pow(KROUPA_BREAK_2, (KROUPA_POWER_2 - KROUPA_POWER_3));

// Constants for the Muller and Mandel remnant mass and kick prescriptions
constexpr double MULLERMANDEL_M1                        = 2.0;	
constexpr double MULLERMANDEL_M2                        = 3.0; 
constexpr double MULLERMANDEL_M3                        = 7.0; 
constexpr double MULLERMANDEL_M4                        = 8.0; 
constexpr double MULLERMANDEL_MU1                       = 1.2;
constexpr double MULLERMANDEL_SIGMA1                    = 0.02;  
constexpr double MULLERMANDEL_MU2A                      = 1.4; 
constexpr double MULLERMANDEL_MU2B                      = 0.5;
constexpr double MULLERMANDEL_SIGMA2                    = 0.05;
constexpr double MULLERMANDEL_MU3A                      = 1.4;
constexpr double MULLERMANDEL_MU3B                      = 0.4;
constexpr double MULLERMANDEL_SIGMA3                    = 0.05;
constexpr double MULLERMANDEL_MUBH                    	= 0.8;
constexpr double MULLERMANDEL_SIGMABH                   = 0.5;
constexpr double MULLERMANDEL_MINNS                     = 1.13;
constexpr double MULLERMANDEL_MAXNS                     = 2.0;
constexpr double MULLERMANDEL_KICKNS                    = 520.0;                                                    // As calibrated by Kapil+ 2023
constexpr double MULLERMANDEL_KICKBH                    = 200.0;
constexpr double MULLERMANDEL_SIGMAKICK                 = 0.3; 

// Constants for WD evolution 

constexpr double COWD_LOG_MDOT_MIN_OFF_CENTER_IGNITION  = -5.688246139;                                             // Minimum log mass accretion rate for off center ignition in a CO WD. From Wang+ 2017. Log( 2.05 x 10^-6). 
constexpr double COWD_MASS_MIN_OFF_CENTER_IGNITION      = 1.33;                                                     // Minimum mass required for off center ignition, as shown in Wang, Podsiadlowski & Han (2017), sect 3.2.
constexpr double HEWD_HE_MDOT_CRIT                      = 2.0E-8;                                                   // Critical accretion rate for He WD accreting He-rich material. From Belczynski+ 2008, Mdot_crit2 in section 5.7.1.
constexpr double HEWD_MINIMUM_MASS_IGNITION             = 0.35;                                                     // Minimum mass for HeMS burning
constexpr double MASS_DOUBLE_DETONATION_CO              = 0.9;                                                      // Minimum mass for detonation which would yield something similar to SN Ia. Ruiter+ 2014.
constexpr double Q_HYDROGEN_BURNING                     = 6.4E18 * MSOL_TO_G / (SECONDS_IN_YEAR * LSOL);            // 6.4E18 is the energy yield of H burning in erg/g as given in Nomoto+ 2007 (2007ApJ...663.1269N)
constexpr double WD_HE_SHELL_MCRIT_DETONATION           = 0.05;                                                     // Minimum shell mass of He for detonation. Should be composed of helium (so, exclude burnt material), but not implemented yet. Ruiter+ 2014.
constexpr double WD_LOG_MT_LIMIT_PIERSANTI_RG_SS_0      = -6.84;
constexpr double WD_LOG_MT_LIMIT_PIERSANTI_RG_SS_1      = 1.349;
constexpr double WD_LOG_MT_LIMIT_PIERSANTI_SS_MF_0      = -8.115;
constexpr double WD_LOG_MT_LIMIT_PIERSANTI_SS_MF_1      = 2.29;
constexpr double WD_LOG_MT_LIMIT_PIERSANTI_SF_Dt_0      = -8.313;
constexpr double WD_LOG_MT_LIMIT_PIERSANTI_SF_Dt_1      = 1.018;
constexpr double WD_LOG_MT_LIMIT_NOMOTO_REDGIANT_0      = -8.33017155;
constexpr double WD_LOG_MT_LIMIT_NOMOTO_REDGIANT_1      = 2.88247131;
constexpr double WD_LOG_MT_LIMIT_NOMOTO_REDGIANT_2      = -0.98023471;
constexpr double WD_LOG_MT_LIMIT_NOMOTO_STABLE_0        = -9.21757267;
constexpr double WD_LOG_MT_LIMIT_NOMOTO_STABLE_1        = 3.57319872;
constexpr double WD_LOG_MT_LIMIT_NOMOTO_STABLE_2        = -1.2137735;

// object types
enum class OBJECT_TYPE: int { NONE, MAIN, PROFILING, UTILS, STAR, BASE_STAR, BINARY_STAR, BASE_BINARY_STAR, BINARY_CONSTITUENT_STAR };    //  if BASE_STAR, check STELLAR_TYPE
const COMPASUnorderedMap<OBJECT_TYPE, std::string> OBJECT_TYPE_LABEL = {
    { OBJECT_TYPE::NONE,                    "Not_an_Object!" },
    { OBJECT_TYPE::MAIN,                    "Main" },
    { OBJECT_TYPE::PROFILING,               "Profiling" },
    { OBJECT_TYPE::UTILS,                   "Utils" },
    { OBJECT_TYPE::STAR,                    "Star" },
    { OBJECT_TYPE::BASE_STAR,               "BaseStar" },
    { OBJECT_TYPE::BINARY_STAR,             "BinaryStar" },
    { OBJECT_TYPE::BASE_BINARY_STAR,        "BaseBinaryStar" },
    { OBJECT_TYPE::BINARY_CONSTITUENT_STAR, "BinaryConstituentStar" },
};


// Commandline Status constants
enum class PROGRAM_STATUS: int { SUCCESS, CONTINUE, STOPPED, ERROR_IN_COMMAND_LINE, LOGGING_FAILED, ERROR_UNHANDLED_EXCEPTION };

// Boost map update options for program options
enum class BOOST_MAP: int { UPDATE, NO_UPDATE };

// Program options origin indicator (command line or gridfile line)
enum class OPTIONS_ORIGIN: int { CMDLINE, GRIDFILE };


// enum class ERROR
// Symbolic names for errors and warnings (error strings below in ERRORLabel map)
// Listed alphabetically (except for 'NONE' - first so ERROR = 0 = NONE)
enum class ERROR: int {
    NONE,                                                           // no error
    AGE_NEGATIVE_ONCE,                                              // age is < 0.0 - invalid
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
    EMPTY_FILE,                                                     // file is empty (contains no content)
    EMPTY_FILENAME,                                                 // filename is an empty string
    ERROR,                                                          // unspecified error
    ERROR_PROCESSING_CMDLINE_OPTIONS,                               // an error occurred while processing commandline options
    ERROR_PROCESSING_GRIDLINE_OPTIONS,                              // an error occurred while processing grid file options
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
    GRID_OPTIONS_ERROR,                                             // grid file options error
    HIGH_TEFF_WINDS,                                                // winds being used at high temperature
    INVALID_DATA_TYPE,                                              // invalid data type
    INVALID_ENVELOPE_TYPE,                                          // invalid envelope type
    INVALID_INITIAL_ATTRIBUTES,                                     // initial values of stellar or binary attributes are not valid - can't evolve star or binary
    INVALID_MASS_TRANSFER_DONOR,                                    // mass transfer from NS, BH or Massless Remnant
    INVALID_RADIUS_INCREASE_ONCE,                                   // radius increased when it should have decreased (or at least remained static)
    INVALID_TYPE_EDDINGTON_RATE,                                    // invalid stellar type for Eddington critical rate calculation
    INVALID_TYPE_MT_MASS_RATIO,                                     // invalid stellar type for mass ratio calculation
    INVALID_TYPE_MT_THERMAL_TIMESCALE,                              // invalid stellar type for thermal timescale calculation
    INVALID_TYPE_ZETA_CALCULATION,                                  // invalid stellar type for Zeta calculation
    INVALID_VALUE_FOR_BOOLEAN_OPTION,                               // invalid value specified for boolean option
    INVALID_VALUE_IN_TIMESTEPS_FILE,                                // invalid value in timesteps file
    LAMBDA_NOT_POSITIVE,                                            // lambda is <= 0.0 - invalid
    LOW_GAMMA,                                                      // very massive mass-loss prescription being extrapolated to low gamma (<0.5)
    LOW_TEFF_WINDS,                                                 // winds being used at low temperature
    MASS_NOT_POSITIVE_ONCE,                                         // mass is <= 0.0 - invalid
    MAXIMUM_MASS_LOST,                                              // (WARNING) maximum mass lost during mass loss calculations
    MISSING_VALUE,                                                  // missing value (e.g. for program option)
    MISSING_RIGHT_BRACKET,                                          // missing right bracket (e.g. for program option range or set specification)
    NO_CONVERGENCE,                                                 // iterative process did not converge
    NO_LAMBDA_DEWI,                                                 // Dewi lambda calculation not supported for stellar type
    NO_LAMBDA_NANJING,                                              // Nanjing lambda calculation not supported for stellar type
    NO_REAL_ROOTS,                                                  // equation has no real roots
    NOT_INITIALISED,                                                // object not initialised
    OPTION_NOT_SUPPORTED_IN_GRID_FILE,                              // option not supported in grid file
    OUT_OF_BOUNDS,                                                  // value out of bounds
    PROGRAM_OPTIONS_ERROR,                                          // program options error
    RADIUS_NOT_POSITIVE,                                            // radius is <= 0.0 - invalid
    RADIUS_NOT_POSITIVE_ONCE,                                       // radius is <= 0.0 - invalid
    RESOLVE_SUPERNOVA_IMPROPERLY_CALLED,                            // ResolveSupernova() called, but m_Supernova->IsSNevent() is false
    ROOT_FINDER_FAILED,                                             // root finder threw an exception
    STELLAR_EVOLUTION_STOPPED,                                      // evolution of current star stopped
    STELLAR_SIMULATION_STOPPED,                                     // stellar simulation stopped
    SUGGEST_HELP,                                                   // suggest using --help
    SWITCH_NOT_TAKEN,                                               // switch to new stellar type not performed
    TIMESTEP_BELOW_MINIMUM,                                         // timestep too small - below minimum
    TIMESTEPS_EXHAUSTED,                                            // timesteps provided exhausted, but evolution not complete
    TIMESTEPS_NOT_CONSUMED,                                         // evolution complete, but provided timesteps not consumed
    TOO_MANY_MASS0_ITERATIONS,                                      // too many iterations in MASS0 root finder
    TOO_MANY_MASS0_TRIES,                                           // too many tries in MASS0 root finder
    TOO_MANY_OMEGA_ITERATIONS,                                      // too many iterations in OMEGA root finder
    TOO_MANY_OMEGA_TRIES,                                           // too many tries in OMEGA root finder
    TOO_MANY_RLOF_ITERATIONS,                                       // too many iterations in RLOF root finder
    TOO_MANY_RLOF_TRIES,                                            // too many tries in RLOF root finder
    TOO_MANY_TIMESTEPS_IN_TIMESTEPS_FILE,                           // too many timesteps in timesteps file (exceeds maximum)
    UNEXPECTED_END_OF_FILE,                                         // unexpected end of file
    UNEXPECTED_LOG_FILE_TYPE,                                       // unexpected log file type
    UNEXPECTED_SN_EVENT,                                            // unexpected supernova event in this context
    UNHANDLED_EXCEPTION,                                            // unhandled exception
    UNKNOWN_A_DISTRIBUTION,                                         // unknown a-distribution
    UNKNOWN_BH_KICK_OPTION,                                         // unknown black hole kick option (in program options)
    UNKNOWN_BINARY_PROPERTY,                                        // unknown binary property
    UNKNOWN_CASE_BB_STABILITY_PRESCRIPTION,                         // unknown case BB/BC mass transfer stability prescription
    UNKNOWN_CE_ACCRETION_PRESCRIPTION,                              // unknown common envelope accretion prescription
    UNKNOWN_CE_FORMALISM,                                           // unknown common envelope formalism
    UNKNOWN_CE_LAMBDA_PRESCRIPTION,                                 // unknown common envelope Lambda Prescription
    UNKNOWN_DATA_TYPE,                                              // unknown data type
    UNKNOWN_ENVELOPE_STATE_PRESCRIPTION,                            // unknown envelope state prescription
    UNKNOWN_INITIAL_MASS_FUNCTION,                                  // unknown initial mass function
    UNKNOWN_KICK_DIRECTION_DISTRIBUTION,                            // unknown kick direction distribution
    UNKNOWN_KICK_MAGNITUDE_DISTRIBUTION,                            // unknown kick magnitude distribution
    UNKNOWN_LOGFILE,                                                // unknown log file
    UNKNOWN_LBV_PRESCRIPTION,                                       // unknown LBV mass loss prescription
    UNKNOWN_MT_ACCRETION_EFFICIENCY_PRESCRIPTION,                   // unknown mass transfer accretion efficiency prescription
    UNKNOWN_MT_ANGULAR_MOMENTUM_LOSS_PRESCRIPTION,                  // unknown mass transfer angular momentum loss prescription
    UNKNOWN_MASS_LOSS_PRESCRIPTION,                                 // unknown mass loss prescription
    UNKNOWN_MT_CASE,                                                // unknown mass transfer case
    UNKNOWN_MT_PRESCRIPTION,                                        // unknown mass transfer prescription
    UNKNOWN_MT_REJUVENATION_PRESCRIPTION,                           // unknown mass transfer rejuvenation prescription
    UNKNOWN_MT_THERMALLY_LIMITED_VARIATION,                         // unknown mass transfer thermally limited variation
    UNKNOWN_NEUTRINO_MASS_LOSS_PRESCRIPTION,                        // unknown neutrino mass loss prescription
    UNKNOWN_NS_EOS,                                                 // unknown NS equation-of-state
    UNKNOWN_PPI_PRESCRIPTION,                                       // unknown pulsational pair instability prescription
    UNKNOWN_PROGRAM_OPTION,                                         // unknown program option
    UNKNOWN_PROPERTY_TYPE,                                          // unknown property type
    UNKNOWN_PULSAR_BIRTH_MAGNETIC_FIELD_DISTRIBUTION,               // unknown pulsar birth magnetic field distribution
    UNKNOWN_PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION,                  // unknown pulsar birth spin period distribution
    UNKNOWN_Q_DISTRIBUTION,                                         // unknown q-distribution
    UNKNOWN_QCRIT_PRESCRIPTION,                                     // unknown critical mass ratio prescription
    UNKNOWN_REMNANT_MASS_PRESCRIPTION,                              // unknown remnant mass prescription
    UNKNOWN_SN_ENGINE,                                              // unknown supernova engine
    UNKNOWN_SN_EVENT,                                               // unknown supernova event encountered
    UNKNOWN_STELLAR_PROPERTY,                                       // unknown stellar property
    UNKNOWN_STELLAR_TYPE,                                           // unknown stellar type
    UNKNOWN_VROT_PRESCRIPTION,                                      // unknown rorational velocity prescription
    UNKNOWN_ZETA_PRESCRIPTION,                                      // unknown stellar Zeta prescription
    UNSUPPORTED_ZETA_PRESCRIPTION,                                  // unsupported common envelope Zeta prescription
    UNSUPPORTED_PULSAR_BIRTH_MAGNETIC_FIELD_DISTRIBUTION,           // unsupported pulsar birth magnetic field distribution
    UNSUPPORTED_PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION,              // unsupported pulsar birth spin period distribution
    UNSUPPORTED_MT_PRESCRIPTION,                                    // unsupported mass transfer prescription
    WARNING,                                                        // unspecified warning
    WHITE_DWARF_TOO_MASSIVE                                         // a white dwarf exceeds the Chandrasekhar mass limit
};


// JR: todo: DOCUMENT these
enum class ERROR_SCOPE: int { NEVER, ALWAYS, FIRST, FIRST_IN_OBJECT_TYPE, FIRST_IN_STELLAR_TYPE, FIRST_IN_OBJECT_ID, FIRST_IN_FUNCTION };


// message catalog
// for now we'll just define it here - one day we might want to have it in a
// file and read it in at program start - that way we can be more flexible
// (i.e. change messages without recompiling, internationalise etc.)
//
// unordered_map - key is integer message number (from enum class ERROR above)
// listed alphabetically

#define ERR_MSG(x) std::get<1>(ERROR_CATALOG.at(x))      // for convenience

const COMPASUnorderedMap<ERROR, std::tuple<ERROR_SCOPE, std::string>> ERROR_CATALOG = {
    { ERROR::AGE_NEGATIVE_ONCE,                                     { ERROR_SCOPE::FIRST_IN_FUNCTION,   "Age < 0.0" }},
    { ERROR::AMBIGUOUS_REMNANT_MASS_PRESCRIPTION,                   { ERROR_SCOPE::ALWAYS,              "Insufficient information to prescribe remnant mass." }},
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
    { ERROR::EMPTY_FILE,                                            { ERROR_SCOPE::ALWAYS,              "File is empty" }},
    { ERROR::EMPTY_FILENAME,                                        { ERROR_SCOPE::ALWAYS,              "Filename is an empty string" }},
    { ERROR::ERROR,                                                 { ERROR_SCOPE::ALWAYS,              "Error!" }},
    { ERROR::ERROR_PROCESSING_CMDLINE_OPTIONS,                      { ERROR_SCOPE::ALWAYS,              "An error occurred while processing commandline options" }},
    { ERROR::ERROR_PROCESSING_GRIDLINE_OPTIONS,                     { ERROR_SCOPE::ALWAYS,              "An error occurred while processing grid file options" }},
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
    { ERROR::GRID_OPTIONS_ERROR,                                    { ERROR_SCOPE::ALWAYS,              "Grid File Options error" }},
    { ERROR::HIGH_TEFF_WINDS,                                       { ERROR_SCOPE::ALWAYS,              "Winds being used at high temperature" }},
    { ERROR::INVALID_DATA_TYPE,                                     { ERROR_SCOPE::ALWAYS,              "Invalid data type" }},
    { ERROR::INVALID_ENVELOPE_TYPE,                                 { ERROR_SCOPE::ALWAYS,              "Invalid envelope type" }},
    { ERROR::INVALID_INITIAL_ATTRIBUTES,                            { ERROR_SCOPE::ALWAYS,              "Initial attributes are not valid - evolution not possible" }},
    { ERROR::INVALID_MASS_TRANSFER_DONOR,                           { ERROR_SCOPE::ALWAYS,              "Mass transfer from NS, BH, or Massless Remnant" }},
    { ERROR::INVALID_RADIUS_INCREASE_ONCE,                          { ERROR_SCOPE::FIRST_IN_FUNCTION,   "Unexpected Radius increase" }},
    { ERROR::INVALID_TYPE_EDDINGTON_RATE,                           { ERROR_SCOPE::ALWAYS,              "Invalid stellar type for Eddington critical rate calculation" }},
    { ERROR::INVALID_TYPE_MT_MASS_RATIO,                            { ERROR_SCOPE::ALWAYS,              "Invalid stellar type for mass ratio calculation" }},
    { ERROR::INVALID_TYPE_MT_THERMAL_TIMESCALE,                     { ERROR_SCOPE::ALWAYS,              "Invalid stellar type for thermal timescale calculation" }},
    { ERROR::INVALID_TYPE_ZETA_CALCULATION,                         { ERROR_SCOPE::ALWAYS,              "Invalid stellar type for Zeta calculation" }},
    { ERROR::INVALID_VALUE_FOR_BOOLEAN_OPTION,                      { ERROR_SCOPE::ALWAYS,              "Invalid value specified for BOOLEAN option" }},
    { ERROR::INVALID_VALUE_IN_TIMESTEPS_FILE,                       { ERROR_SCOPE::ALWAYS,              "Invalid value in timesteps file" }},
    { ERROR::LAMBDA_NOT_POSITIVE,                                   { ERROR_SCOPE::ALWAYS,              "Lambda <= 0.0" }},
    { ERROR::LOW_GAMMA,                                             { ERROR_SCOPE::ALWAYS,              "Very massive prescription being extrapolated to low gamma (<0.5)" }},
    { ERROR::LOW_TEFF_WINDS,                                        { ERROR_SCOPE::ALWAYS,              "Winds being used at low temperature" }},
    { ERROR::MASS_NOT_POSITIVE_ONCE,                                { ERROR_SCOPE::FIRST_IN_FUNCTION,   "Mass <= 0.0" }},
    { ERROR::MAXIMUM_MASS_LOST,                                     { ERROR_SCOPE::ALWAYS,              "Maximum mass lost during mass loss calculations" }},
    { ERROR::MISSING_VALUE,                                         { ERROR_SCOPE::ALWAYS,              "Missing value" }},
    { ERROR::MISSING_RIGHT_BRACKET,                                 { ERROR_SCOPE::ALWAYS,              "Missing ']'" }},
    { ERROR::NO_CONVERGENCE,                                        { ERROR_SCOPE::ALWAYS,              "No convergence" }},
    { ERROR::NO_LAMBDA_DEWI,                                        { ERROR_SCOPE::ALWAYS,              "Dewi lambda calculation not supported for stellar type" }},
    { ERROR::NO_LAMBDA_NANJING,                                     { ERROR_SCOPE::ALWAYS,              "Nanjing lambda calculation not supported for stellar type" }},
    { ERROR::NO_REAL_ROOTS,                                         { ERROR_SCOPE::ALWAYS,              "No real roots" }},
    { ERROR::NONE,                                                  { ERROR_SCOPE::NEVER,               "No error" }},
    { ERROR::NOT_INITIALISED,                                       { ERROR_SCOPE::ALWAYS,              "Object not initialised" }},
    { ERROR::OPTION_NOT_SUPPORTED_IN_GRID_FILE,                     { ERROR_SCOPE::ALWAYS,              "Option not supported in grid file" }},
    { ERROR::OUT_OF_BOUNDS,                                         { ERROR_SCOPE::ALWAYS,              "Value out of bounds" }},
    { ERROR::PROGRAM_OPTIONS_ERROR,                                 { ERROR_SCOPE::ALWAYS,              "Commandline Options error" }},
    { ERROR::RADIUS_NOT_POSITIVE,                                   { ERROR_SCOPE::ALWAYS,              "Radius <= 0.0" }},
    { ERROR::RADIUS_NOT_POSITIVE_ONCE,                              { ERROR_SCOPE::FIRST_IN_FUNCTION,   "Radius <= 0.0" }},
    { ERROR::RESOLVE_SUPERNOVA_IMPROPERLY_CALLED,                   { ERROR_SCOPE::ALWAYS,              "ResolveSupernova() called, but m_Supernova->IsSNevent() is false" }},
    { ERROR::ROOT_FINDER_FAILED,                                    { ERROR_SCOPE::ALWAYS,              "Exception encountered in root finder" }},
    { ERROR::STELLAR_EVOLUTION_STOPPED,                             { ERROR_SCOPE::ALWAYS,              "Evolution of current star stopped" }},
    { ERROR::STELLAR_SIMULATION_STOPPED,                            { ERROR_SCOPE::ALWAYS,              "Stellar simulation stopped" }},
    { ERROR::SUGGEST_HELP,                                          { ERROR_SCOPE::ALWAYS,              "Use option '-h' (or '--help') to see (descriptions of) available options" }},
    { ERROR::SWITCH_NOT_TAKEN,                                      { ERROR_SCOPE::ALWAYS,              "Switch to new stellar type not performed" }},
    { ERROR::TIMESTEP_BELOW_MINIMUM,                                { ERROR_SCOPE::ALWAYS,              "Timestep below minimum - timestep taken" }},
    { ERROR::TIMESTEPS_EXHAUSTED,                                   { ERROR_SCOPE::ALWAYS,              "Provided timesteps exhausted, but evolution not complete" }},
    { ERROR::TIMESTEPS_NOT_CONSUMED,                                { ERROR_SCOPE::ALWAYS,              "Evolution complete, but provided timesteps not consumed" }},
    { ERROR::TOO_MANY_MASS0_ITERATIONS,                             { ERROR_SCOPE::ALWAYS,              "Reached maximum number of iterations when looking for effective initial mass Mass_0 to match desired stellar core of HG star following case A mass transfer" }},
    { ERROR::TOO_MANY_MASS0_TRIES,                                  { ERROR_SCOPE::ALWAYS,              "Reached maximum number of tries when looking for effective initial mass Mass_0 to match desired stellar core of HG star following case A mass transfer" }},
    { ERROR::TOO_MANY_OMEGA_ITERATIONS,                             { ERROR_SCOPE::ALWAYS,              "Reached maximum number of iterations when looking for omega when circularising and synchronising for tides" }},
    { ERROR::TOO_MANY_OMEGA_TRIES,                                  { ERROR_SCOPE::ALWAYS,              "Reached maximum number of tries when looking for omega when circularising and synchronising for tides" }},
    { ERROR::TOO_MANY_RLOF_ITERATIONS,                              { ERROR_SCOPE::ALWAYS,              "Reached maximum number of iterations when fitting star inside Roche Lobe in RLOF" }},
    { ERROR::TOO_MANY_RLOF_TRIES,                                   { ERROR_SCOPE::ALWAYS,              "Reached maximum number of tries when fitting star inside Roche Lobe in RLOF" }},
    { ERROR::TOO_MANY_TIMESTEPS_IN_TIMESTEPS_FILE,                  { ERROR_SCOPE::ALWAYS,              "Number of timesteps in timestpes file exceeds maximum timesteps" }},
    { ERROR::UNEXPECTED_END_OF_FILE,                                { ERROR_SCOPE::ALWAYS,              "Unexpected end of file" }},
    { ERROR::UNEXPECTED_LOG_FILE_TYPE,                              { ERROR_SCOPE::ALWAYS,              "Unexpected log file type" }},
    { ERROR::UNEXPECTED_SN_EVENT,                                   { ERROR_SCOPE::ALWAYS,              "Unexpected supernova event in this context" }},
    { ERROR::UNHANDLED_EXCEPTION,                                   { ERROR_SCOPE::ALWAYS,              "Unhandled exception" }},
    { ERROR::UNKNOWN_A_DISTRIBUTION,                                { ERROR_SCOPE::ALWAYS,              "Unknown semi-major-axis distribution" }},
    { ERROR::UNKNOWN_BH_KICK_OPTION,                                { ERROR_SCOPE::ALWAYS,              "Unknown black hole kick option" }},
    { ERROR::UNKNOWN_BINARY_PROPERTY,                               { ERROR_SCOPE::ALWAYS,              "Unknown binary property - property details not found" }},
    { ERROR::UNKNOWN_CASE_BB_STABILITY_PRESCRIPTION,                { ERROR_SCOPE::ALWAYS,              "Unknown case BB/BC mass transfer stability prescription" }},
    { ERROR::UNKNOWN_CE_ACCRETION_PRESCRIPTION,                     { ERROR_SCOPE::ALWAYS,              "Unknown common envelope accretion prescription" }},
    { ERROR::UNKNOWN_CE_FORMALISM,                                  { ERROR_SCOPE::ALWAYS,              "Unknown common envelope formalism" }},
    { ERROR::UNKNOWN_CE_LAMBDA_PRESCRIPTION,                        { ERROR_SCOPE::ALWAYS,              "Unknown common envelope lambda prescription" }},
    { ERROR::UNKNOWN_DATA_TYPE,                                     { ERROR_SCOPE::ALWAYS,              "Unknown data type" }},
    { ERROR::UNKNOWN_ENVELOPE_STATE_PRESCRIPTION,                   { ERROR_SCOPE::ALWAYS,              "Unknown envelope state prescription" }},
    { ERROR::UNKNOWN_INITIAL_MASS_FUNCTION,                         { ERROR_SCOPE::ALWAYS,              "Unknown initial mass function (IMF)" }},
    { ERROR::UNKNOWN_KICK_DIRECTION_DISTRIBUTION,                   { ERROR_SCOPE::ALWAYS,              "Unknown kick direction distribution" }},
    { ERROR::UNKNOWN_KICK_MAGNITUDE_DISTRIBUTION,                   { ERROR_SCOPE::ALWAYS,              "Unknown kick magnitude distribution" }},
    { ERROR::UNKNOWN_LOGFILE,                                       { ERROR_SCOPE::ALWAYS,              "Unknown log file" }},
    { ERROR::UNKNOWN_MT_CASE,                                       { ERROR_SCOPE::ALWAYS,              "Unknown mass transfer case" }},
    { ERROR::UNKNOWN_MT_PRESCRIPTION,                               { ERROR_SCOPE::ALWAYS,              "Unknown mass transfer prescription" }},
    { ERROR::UNKNOWN_MT_ACCRETION_EFFICIENCY_PRESCRIPTION,          { ERROR_SCOPE::ALWAYS,              "Unknown mass transfer accretion efficiency prescription" }},
    { ERROR::UNKNOWN_MT_ANGULAR_MOMENTUM_LOSS_PRESCRIPTION,         { ERROR_SCOPE::ALWAYS,              "Unknown mass transfer angular momentum loss prescription" }},
    { ERROR::UNKNOWN_MT_REJUVENATION_PRESCRIPTION,                  { ERROR_SCOPE::ALWAYS,              "Unknown mass transfer rejuvenation prescription" }},
    { ERROR::UNKNOWN_MT_THERMALLY_LIMITED_VARIATION,                { ERROR_SCOPE::ALWAYS,              "Unknown mass transfer thermally limited variation" }},
    { ERROR::UNKNOWN_MASS_LOSS_PRESCRIPTION,                        { ERROR_SCOPE::ALWAYS,              "Unknown mass loss prescription" }},
    { ERROR::UNKNOWN_NEUTRINO_MASS_LOSS_PRESCRIPTION,               { ERROR_SCOPE::ALWAYS,              "Unknown neutrino mass loss prescription" }},
    { ERROR::UNKNOWN_NS_EOS,                                        { ERROR_SCOPE::ALWAYS,              "Unknown NS equation-of-state" }},
    { ERROR::UNKNOWN_PPI_PRESCRIPTION,                              { ERROR_SCOPE::ALWAYS,              "Unknown pulsational pair instability prescription" }},
    { ERROR::UNKNOWN_PROGRAM_OPTION,                                { ERROR_SCOPE::ALWAYS,              "Unknown program option - property details not found" }},
    { ERROR::UNKNOWN_PROPERTY_TYPE,                                 { ERROR_SCOPE::ALWAYS,              "Unknown property type - property details not found" }},
    { ERROR::UNKNOWN_PULSAR_BIRTH_MAGNETIC_FIELD_DISTRIBUTION,      { ERROR_SCOPE::ALWAYS,              "Unknown pulsar birth magnetic field distribution" }},
    { ERROR::UNKNOWN_PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION,         { ERROR_SCOPE::ALWAYS,              "Unknown pulsar birth spin period distribution" }},
    { ERROR::UNKNOWN_Q_DISTRIBUTION,                                { ERROR_SCOPE::ALWAYS,              "Unknown q-distribution" }},
    { ERROR::UNKNOWN_QCRIT_PRESCRIPTION,                            { ERROR_SCOPE::ALWAYS,              "Unknown critical mass ratio prescription" }},
    { ERROR::UNKNOWN_REMNANT_MASS_PRESCRIPTION,                     { ERROR_SCOPE::ALWAYS,              "Unknown remnant mass prescription" }},
    { ERROR::UNKNOWN_SN_ENGINE,                                     { ERROR_SCOPE::ALWAYS,              "Unknown supernova engine" }},
    { ERROR::UNKNOWN_SN_EVENT,                                      { ERROR_SCOPE::ALWAYS,              "Unknown supernova event" }},
    { ERROR::UNKNOWN_STELLAR_PROPERTY,                              { ERROR_SCOPE::ALWAYS,              "Unknown stellar property - property details not found" }},
    { ERROR::UNKNOWN_STELLAR_TYPE,                                  { ERROR_SCOPE::ALWAYS,              "Unknown stellar type" }},
    { ERROR::UNKNOWN_VROT_PRESCRIPTION,                             { ERROR_SCOPE::ALWAYS,              "Unknown rotational velocity prescription" }},
    { ERROR::UNKNOWN_ZETA_PRESCRIPTION,                             { ERROR_SCOPE::ALWAYS,              "Unknown stellar Zeta prescription" }},
    { ERROR::UNSUPPORTED_MT_PRESCRIPTION,                           { ERROR_SCOPE::ALWAYS,              "Unsupported mass transfer prescription" }},
    { ERROR::UNSUPPORTED_PULSAR_BIRTH_MAGNETIC_FIELD_DISTRIBUTION,  { ERROR_SCOPE::ALWAYS,              "Unsupported pulsar birth magnetic field distribution" }},
    { ERROR::UNSUPPORTED_PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION,     { ERROR_SCOPE::ALWAYS,              "Unsupported pulsar birth spin period distribution" }},
    { ERROR::UNSUPPORTED_ZETA_PRESCRIPTION,                         { ERROR_SCOPE::ALWAYS,              "Unsupported stellar Zeta prescription" }},
    { ERROR::WARNING,                                               { ERROR_SCOPE::ALWAYS,              "Warning!" }},
    { ERROR::WHITE_DWARF_TOO_MASSIVE,                               { ERROR_SCOPE::ALWAYS,              "This white dwarf exceeds the Chandrasekhar mass limit" }}
};


// Evolution status constants
enum class EVOLUTION_STATUS: int {
    CONTINUE,
    DONE,
    ERROR,
    TIMES_UP,
    STEPS_UP,
    NO_TIMESTEPS,
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
    UNBOUND
};

// JR: deliberately kept these messages succinct (where I could) so running status doesn't scroll off the page...
const COMPASUnorderedMap<EVOLUTION_STATUS, std::string> EVOLUTION_STATUS_LABEL = {
    { EVOLUTION_STATUS::CONTINUE,                "Continue evolution" },
    { EVOLUTION_STATUS::DONE,                    "Simulation completed" },
    { EVOLUTION_STATUS::ERROR,                   "Evolution stopped because an error occurred" },
    { EVOLUTION_STATUS::TIMES_UP,                "Allowed time exceeded" },
    { EVOLUTION_STATUS::STEPS_UP,                "Allowed timesteps exceeded" },
    { EVOLUTION_STATUS::NO_TIMESTEPS,            "No user-provided timesteps read" },
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
    { EVOLUTION_STATUS::UNBOUND,                 "Unbound binary" }
};


// Evolution mode (SSE or BSE)
enum class EVOLUTION_MODE: int { SSE, BSE };
const COMPASUnorderedMap<EVOLUTION_MODE, std::string> EVOLUTION_MODE_LABEL = {
    { EVOLUTION_MODE::SSE, "SSE" },
    { EVOLUTION_MODE::BSE, "BSE" }
};


// user specified distributions, assumptions etc.

// Black Hole Kick Options
enum class BLACK_HOLE_KICKS: int { FULL, REDUCED, ZERO, FALLBACK };
const COMPASUnorderedMap<BLACK_HOLE_KICKS, std::string> BLACK_HOLE_KICKS_LABEL = {
    { BLACK_HOLE_KICKS::FULL,     "FULL" },     // FULL kicks
    { BLACK_HOLE_KICKS::REDUCED,  "REDUCED" },  // REDUCED kicks
    { BLACK_HOLE_KICKS::ZERO,     "ZERO" },     // ZERO kicks
    { BLACK_HOLE_KICKS::FALLBACK, "FALLBACK" }  // FALLBACK kicks
};


// Kick magnitude distribution from Bray & Eldridge 2016,2018
enum class BRAY_ELDRIDGE_CONSTANT: int { ALPHA, BETA };
const COMPASUnorderedMap<BRAY_ELDRIDGE_CONSTANT, double> BRAY_ELDRIDGE_CONSTANT_VALUES = {
    { BRAY_ELDRIDGE_CONSTANT::ALPHA, 100.0 },
    { BRAY_ELDRIDGE_CONSTANT::BETA, -170.0 }
};

enum class CASE_BB_STABILITY_PRESCRIPTION: int{ ALWAYS_STABLE, ALWAYS_STABLE_ONTO_NSBH, TREAT_AS_OTHER_MT, ALWAYS_UNSTABLE };
const COMPASUnorderedMap<CASE_BB_STABILITY_PRESCRIPTION, std::string> CASE_BB_STABILITY_PRESCRIPTION_LABEL = {
    { CASE_BB_STABILITY_PRESCRIPTION::ALWAYS_STABLE,           "ALWAYS_STABLE" },
    { CASE_BB_STABILITY_PRESCRIPTION::ALWAYS_STABLE_ONTO_NSBH, "ALWAYS_STABLE_ONTO_NSBH" },
    { CASE_BB_STABILITY_PRESCRIPTION::TREAT_AS_OTHER_MT,       "TREAT_AS_OTHER_MT" },
    { CASE_BB_STABILITY_PRESCRIPTION::ALWAYS_UNSTABLE,         "ALWAYS_UNSTABLE" }
};

// Common Envelope Accretion Prescriptions
enum class CE_ACCRETION_PRESCRIPTION: int { ZERO, CONSTANT, UNIFORM, MACLEOD, CHEVALIER };
const COMPASUnorderedMap<CE_ACCRETION_PRESCRIPTION, std::string> CE_ACCRETION_PRESCRIPTION_LABEL = {
    { CE_ACCRETION_PRESCRIPTION::ZERO,      "ZERO" },
    { CE_ACCRETION_PRESCRIPTION::CONSTANT,  "CONSTANT" },
    { CE_ACCRETION_PRESCRIPTION::UNIFORM,   "UNIFORM" },
    { CE_ACCRETION_PRESCRIPTION::MACLEOD,   "MACLEOD" },
    { CE_ACCRETION_PRESCRIPTION::CHEVALIER, "CHEVALIER" }
};

// Envelope State Prescriptions
enum class ENVELOPE_STATE_PRESCRIPTION: int { LEGACY, HURLEY, FIXED_TEMPERATURE };
const COMPASUnorderedMap<ENVELOPE_STATE_PRESCRIPTION, std::string> ENVELOPE_STATE_PRESCRIPTION_LABEL = {
    { ENVELOPE_STATE_PRESCRIPTION::LEGACY,            "LEGACY" },
    { ENVELOPE_STATE_PRESCRIPTION::HURLEY,            "HURLEY" },
    { ENVELOPE_STATE_PRESCRIPTION::FIXED_TEMPERATURE, "FIXED_TEMPERATURE" }
};


// Common Envelope Lambda Prescriptions
enum class CE_LAMBDA_PRESCRIPTION: int { FIXED, LOVERIDGE, NANJING, KRUCKOW, DEWI };
const COMPASUnorderedMap<CE_LAMBDA_PRESCRIPTION, std::string> CE_LAMBDA_PRESCRIPTION_LABEL = {
    { CE_LAMBDA_PRESCRIPTION::FIXED,     "LAMBDA_FIXED" },
    { CE_LAMBDA_PRESCRIPTION::LOVERIDGE, "LAMBDA_LOVERIDGE" },
    { CE_LAMBDA_PRESCRIPTION::NANJING,   "LAMBDA_NANJING" },
    { CE_LAMBDA_PRESCRIPTION::KRUCKOW,   "LAMBDA_KRUCKOW" },
    { CE_LAMBDA_PRESCRIPTION::DEWI,      "LAMBDA_DEWI" }
};
    
    
// Common Envelope Formalism
enum class CE_FORMALISM: int { ENERGY, TWO_STAGE };
const COMPASUnorderedMap<CE_FORMALISM, std::string> CE_FORMALISM_LABEL = {
        { CE_FORMALISM::ENERGY,    "ENERGY" },
        { CE_FORMALISM::TWO_STAGE, "TWO_STAGE" }
};


// Common envelope zeta prescription
enum class ZETA_PRESCRIPTION: int { SOBERMAN, HURLEY, ARBITRARY, NONE };
const COMPASUnorderedMap<ZETA_PRESCRIPTION, std::string> ZETA_PRESCRIPTION_LABEL = {
    { ZETA_PRESCRIPTION::SOBERMAN,  "SOBERMAN" },
    { ZETA_PRESCRIPTION::HURLEY,    "HURLEY" },
    { ZETA_PRESCRIPTION::ARBITRARY, "ARBITRARY" },
    { ZETA_PRESCRIPTION::NONE,      "NONE" },
};

// Critical Mass Ratio prescription
enum class QCRIT_PRESCRIPTION: int { NONE, CLAEYS, GE20, GE20_IC, HURLEY_HJELLMING_WEBBINK};
const COMPASUnorderedMap<QCRIT_PRESCRIPTION, std::string> QCRIT_PRESCRIPTION_LABEL = {
    { QCRIT_PRESCRIPTION::NONE,    "NONE" },
    { QCRIT_PRESCRIPTION::CLAEYS,  "CLAEYS" },
    { QCRIT_PRESCRIPTION::GE20,    "GE20" },
    { QCRIT_PRESCRIPTION::GE20_IC, "GE20_IC" },
    { QCRIT_PRESCRIPTION::HURLEY_HJELLMING_WEBBINK, "HURLEY_HJELLMING_WEBBINK" },
};


// CHE (Chemically Homogeneous Evolution) Options
enum class CHE_MODE: int { NONE, OPTIMISTIC, PESSIMISTIC };
const COMPASUnorderedMap<CHE_MODE, std::string> CHE_MODE_LABEL = {
    { CHE_MODE::NONE,        "NONE" },
    { CHE_MODE::OPTIMISTIC,  "OPTIMISTIC" },
    { CHE_MODE::PESSIMISTIC, "PESSIMISTIC" }
};


// Logfile file types
const COMPASUnorderedMap<LOGFILETYPE, std::string> LOGFILETYPELabel = {     // labels
    { LOGFILETYPE::NONE, "NONE" },
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


// Logfile delimiters
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


// Logfile modifiers

// Option to add program option columns to [BSE/SSE] SYSPARMS file
enum class ADD_OPTIONS_TO_SYSPARMS: int { ALWAYS, GRID, NEVER };
const COMPASUnorderedMap<ADD_OPTIONS_TO_SYSPARMS, std::string> ADD_OPTIONS_TO_SYSPARMS_LABEL = {
    { ADD_OPTIONS_TO_SYSPARMS::ALWAYS, "ALWAYS" },
    { ADD_OPTIONS_TO_SYSPARMS::GRID,   "GRID" },
    { ADD_OPTIONS_TO_SYSPARMS::NEVER,  "NEVER" }
};


// Eccentricity distribution
enum class ECCENTRICITY_DISTRIBUTION: int { ZERO, FLAT, THERMAL, GELLER_2013, DUQUENNOYMAYOR1991, SANA2012 };
const COMPASUnorderedMap<ECCENTRICITY_DISTRIBUTION, std::string> ECCENTRICITY_DISTRIBUTION_LABEL = {
    { ECCENTRICITY_DISTRIBUTION::ZERO,               "ZERO" },
    { ECCENTRICITY_DISTRIBUTION::FLAT,               "FLAT" },
    { ECCENTRICITY_DISTRIBUTION::THERMAL,            "THERMAL" },
    { ECCENTRICITY_DISTRIBUTION::GELLER_2013,        "GELLER+2013" },
    { ECCENTRICITY_DISTRIBUTION::DUQUENNOYMAYOR1991, "DUQUENNOYMAYOR1991" },
    { ECCENTRICITY_DISTRIBUTION::SANA2012,           "SANA2012"}
};


// Envelope types
enum class ENVELOPE: int { RADIATIVE, CONVECTIVE, REMNANT };
const COMPASUnorderedMap<ENVELOPE, std::string> ENVELOPE_LABEL = {
    { ENVELOPE::RADIATIVE,  "RADIATIVE" },
    { ENVELOPE::CONVECTIVE, "CONVECTIVE" },
    { ENVELOPE::REMNANT,    "REMNANT" }
};


// Kick magnitude distribution
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


// Kick direction distribution
enum class KICK_DIRECTION_DISTRIBUTION: int { ISOTROPIC, INPLANE, PERPENDICULAR, POWERLAW, WEDGE, POLES };
const COMPASUnorderedMap<KICK_DIRECTION_DISTRIBUTION, std::string> KICK_DIRECTION_DISTRIBUTION_LABEL = {
    { KICK_DIRECTION_DISTRIBUTION::ISOTROPIC,     "ISOTROPIC" },
    { KICK_DIRECTION_DISTRIBUTION::INPLANE,       "INPLANE" },
    { KICK_DIRECTION_DISTRIBUTION::PERPENDICULAR, "PERPENDICULAR" },
    { KICK_DIRECTION_DISTRIBUTION::POWERLAW,      "POWERLAW" },
    { KICK_DIRECTION_DISTRIBUTION::WEDGE,         "WEDGE" },
    { KICK_DIRECTION_DISTRIBUTION::POLES,         "POLES" }
};


// Initial mass function
enum class INITIAL_MASS_FUNCTION: int { SALPETER, POWERLAW, UNIFORM, KROUPA };
const COMPASUnorderedMap<INITIAL_MASS_FUNCTION, std::string> INITIAL_MASS_FUNCTION_LABEL = {
    { INITIAL_MASS_FUNCTION::SALPETER, "SALPETER" },
    { INITIAL_MASS_FUNCTION::POWERLAW, "POWERLAW" },
    { INITIAL_MASS_FUNCTION::UNIFORM,  "UNIFORM" },
    { INITIAL_MASS_FUNCTION::KROUPA,   "KROUPA" }
};

// LBV Mass loss prescriptions
enum class LBV_PRESCRIPTION: int { NONE, HURLEY_ADD, HURLEY, BELCZYNSKI };
const COMPASUnorderedMap<LBV_PRESCRIPTION, std::string> LBV_PRESCRIPTION_LABEL = {
    { LBV_PRESCRIPTION::NONE,       "NONE" },
    { LBV_PRESCRIPTION::HURLEY_ADD, "HURLEY_ADD" },
    { LBV_PRESCRIPTION::HURLEY,     "HURLEY" },
    { LBV_PRESCRIPTION::BELCZYNSKI, "BELCZYNSKI" }
};

// OB (main sequence) Mass loss prescriptions
enum class OB_MASS_LOSS: int { NONE, VINK2001, VINK2021, BJORKLUND2022, KRTICKA2018};
const COMPASUnorderedMap<OB_MASS_LOSS, std::string> OB_MASS_LOSS_LABEL = {
    { OB_MASS_LOSS::NONE,          "NONE" },
    { OB_MASS_LOSS::VINK2001,      "VINK2001" },
    { OB_MASS_LOSS::VINK2021,      "VINK2021" },
    { OB_MASS_LOSS::BJORKLUND2022, "BJORKLUND2022" },
    { OB_MASS_LOSS::KRTICKA2018,   "KRTICKA2018" },
};

// Very Massive Mass loss prescriptions
enum class VMS_MASS_LOSS: int { NONE, VINK2011, BESTENLEHNER2020, SABHAHIT2023};
const COMPASUnorderedMap<VMS_MASS_LOSS, std::string> VMS_MASS_LOSS_LABEL = {
    { VMS_MASS_LOSS::NONE,             "NONE" },
    { VMS_MASS_LOSS::VINK2011,         "VINK2011" },
    { VMS_MASS_LOSS::BESTENLEHNER2020, "BESTENLEHNER2020" },
    { VMS_MASS_LOSS::SABHAHIT2023,     "SABHAHIT2023" },
};

// RSG Mass loss prescriptions
enum class RSG_MASS_LOSS: int { NONE, VINKSABHAHIT2023, BEASOR2020, DECIN2023, YANG2023, KEE2021, NJ90};
const COMPASUnorderedMap<RSG_MASS_LOSS, std::string> RSG_MASS_LOSS_LABEL = {
    { RSG_MASS_LOSS::NONE,             "NONE" },
    { RSG_MASS_LOSS::VINKSABHAHIT2023, "VINKSABHAHIT2023" },
    { RSG_MASS_LOSS::BEASOR2020,       "BEASOR2020" },
    { RSG_MASS_LOSS::DECIN2023,        "DECIN2023" },
    { RSG_MASS_LOSS::YANG2023,         "YANG2023" },
    { RSG_MASS_LOSS::KEE2021,          "KEE2021" },
    { RSG_MASS_LOSS::NJ90,             "NJ90" },
};

// WR Mass loss prescriptions
enum class WR_MASS_LOSS: int { BELCZYNSKI2010, SANDERVINK2023, SHENAR2019 };
const COMPASUnorderedMap<WR_MASS_LOSS, std::string> WR_MASS_LOSS_LABEL = {
    { WR_MASS_LOSS::BELCZYNSKI2010, "BELCZYNSKI2010" },
    { WR_MASS_LOSS::SANDERVINK2023, "SANDERVINK2023" },
    { WR_MASS_LOSS::SHENAR2019,     "SHENAR2019" },
};

// Mass loss prescriptions
enum class MASS_LOSS_PRESCRIPTION: int { NONE, HURLEY, BELCZYNSKI2010, FLEXIBLE2023 };
const COMPASUnorderedMap<MASS_LOSS_PRESCRIPTION, std::string> MASS_LOSS_PRESCRIPTION_LABEL = {
    { MASS_LOSS_PRESCRIPTION::NONE,           "NONE" },
    { MASS_LOSS_PRESCRIPTION::HURLEY,         "HURLEY" },
    { MASS_LOSS_PRESCRIPTION::BELCZYNSKI2010, "BELCZYNSKI2010" },
    { MASS_LOSS_PRESCRIPTION::FLEXIBLE2023,   "FLEXIBLE2023"}
};


// Mass ratio distribution
enum class MASS_RATIO_DISTRIBUTION: int { FLAT, DUQUENNOYMAYOR1991, SANA2012 };
const COMPASUnorderedMap<MASS_RATIO_DISTRIBUTION, std::string> MASS_RATIO_DISTRIBUTION_LABEL = {
    { MASS_RATIO_DISTRIBUTION::FLAT,               "FLAT" },
    { MASS_RATIO_DISTRIBUTION::DUQUENNOYMAYOR1991, "DUQUENNOYMAYOR1991" },
    { MASS_RATIO_DISTRIBUTION::SANA2012,           "SANA2012" }
};


// Mass Transfer types
enum class MASS_TRANSFER: int { NONE, STABLE_TIMESCALE_NUCLEAR, STABLE_TIMESCALE_THERMAL, UNSTABLE_TIMESCALE_THERMAL, UNSTABLE_TIMESCALE_DYNAMICAL};
const COMPASUnorderedMap<MASS_TRANSFER, std::string> MASS_TRANSFER_LABEL = {
    { MASS_TRANSFER::NONE,                         "No mass transfer" },
    { MASS_TRANSFER::STABLE_TIMESCALE_NUCLEAR,     "Nuclear timescale stable mass transfer" },
    { MASS_TRANSFER::STABLE_TIMESCALE_THERMAL,     "Thermal timescale stable mass transfer" },
    { MASS_TRANSFER::UNSTABLE_TIMESCALE_THERMAL,   "Thermal timescale unstable mass transfer" },
    { MASS_TRANSFER::UNSTABLE_TIMESCALE_DYNAMICAL, "Dynamical timescale unstable mass transfer" }
};


// Mass transfer accretion efficiency prescriptions
enum class MT_ACCRETION_EFFICIENCY_PRESCRIPTION: int { THERMALLY_LIMITED, FIXED_FRACTION};
const COMPASUnorderedMap<MT_ACCRETION_EFFICIENCY_PRESCRIPTION, std::string> MT_ACCRETION_EFFICIENCY_PRESCRIPTION_LABEL = {
    { MT_ACCRETION_EFFICIENCY_PRESCRIPTION::THERMALLY_LIMITED, "THERMAL" },
    { MT_ACCRETION_EFFICIENCY_PRESCRIPTION::FIXED_FRACTION,    "FIXED" }
};


// Mass transfer angular momentum loss prescriptions
enum class MT_ANGULAR_MOMENTUM_LOSS_PRESCRIPTION: int { JEANS, ISOTROPIC_RE_EMISSION, CIRCUMBINARY_RING, MACLEOD_LINEAR, ARBITRARY };
const COMPASUnorderedMap<MT_ANGULAR_MOMENTUM_LOSS_PRESCRIPTION, std::string> MT_ANGULAR_MOMENTUM_LOSS_PRESCRIPTION_LABEL = {
    { MT_ANGULAR_MOMENTUM_LOSS_PRESCRIPTION::JEANS,                 "JEANS" },
    { MT_ANGULAR_MOMENTUM_LOSS_PRESCRIPTION::ISOTROPIC_RE_EMISSION, "ISOTROPIC" },
    { MT_ANGULAR_MOMENTUM_LOSS_PRESCRIPTION::CIRCUMBINARY_RING,     "CIRCUMBINARY" },
    { MT_ANGULAR_MOMENTUM_LOSS_PRESCRIPTION::MACLEOD_LINEAR,        "MACLEOD_LINEAR" },
    { MT_ANGULAR_MOMENTUM_LOSS_PRESCRIPTION::ARBITRARY,             "ARBITRARY" }
};


// Mass transfer cases
enum class MT_CASE: int { NONE, A, B, C, OTHER };
const COMPASUnorderedMap<MT_CASE, std::string> MT_CASE_LABEL = {
    { MT_CASE::NONE, "Mass Transfer CASE NONE: No Mass Transfer" },
    { MT_CASE::A,    "Mass Transfer CASE A" },                          // mass transfer while donor is on main sequence
    { MT_CASE::B,    "Mass Transfer CASE B" },                          // donor star is in (or evolving to) Red Giant phase
    { MT_CASE::C,    "Mass Transfer CASE C" },                          // SuperGiant phase
    { MT_CASE::OTHER,"Mass Transfer CASE OTHER: Some combination" }     // default value, or multiple MT events
};


// Mass transfer prescriptions
enum class MT_PRESCRIPTION: int { HURLEY, BELCZYNSKI, NONE };
const COMPASUnorderedMap<MT_PRESCRIPTION, std::string> MT_PRESCRIPTION_LABEL = {
    { MT_PRESCRIPTION::NONE,       "NONE" },
    { MT_PRESCRIPTION::HURLEY,     "HURLEY" },
    { MT_PRESCRIPTION::BELCZYNSKI, "BELCZYNSKI" }
};


// Mass Transfer Thermally limited Variation options
enum class MT_THERMALLY_LIMITED_VARIATION: int { C_FACTOR, RADIUS_TO_ROCHELOBE };
const COMPASUnorderedMap<MT_THERMALLY_LIMITED_VARIATION, std::string> MT_THERMALLY_LIMITED_VARIATION_LABEL = {
    { MT_THERMALLY_LIMITED_VARIATION::C_FACTOR,            "CFACTOR" },
    { MT_THERMALLY_LIMITED_VARIATION::RADIUS_TO_ROCHELOBE, "ROCHELOBE" }
};


// Mass transfer rejuvenation prescription
enum class MT_REJUVENATION_PRESCRIPTION: int { NONE, STARTRACK };
const COMPASUnorderedMap<MT_REJUVENATION_PRESCRIPTION, std::string> MT_REJUVENATION_PRESCRIPTION_LABEL = {
    { MT_REJUVENATION_PRESCRIPTION::NONE,      "NONE" },
    { MT_REJUVENATION_PRESCRIPTION::STARTRACK, "STARTRACK" }
};


// Mass transfer tracking constants
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


// Mass transfer timing options for writing to BSE_RLOF file
// Doesn't need labels...
enum class MASS_TRANSFER_TIMING: int { PRE_MT, POST_MT };

// Metallicity distribution
enum class METALLICITY_DISTRIBUTION: int { ZSOLAR, LOGUNIFORM };
const COMPASUnorderedMap<METALLICITY_DISTRIBUTION, std::string> METALLICITY_DISTRIBUTION_LABEL = {
    { METALLICITY_DISTRIBUTION::ZSOLAR,     "ZSOLAR" },
    { METALLICITY_DISTRIBUTION::LOGUNIFORM, "LOGUNIFORM" }
};

// Neutrino mass loss BH formation prescriptions
enum class NEUTRINO_MASS_LOSS_PRESCRIPTION: int { FIXED_FRACTION, FIXED_MASS };
const COMPASUnorderedMap<NEUTRINO_MASS_LOSS_PRESCRIPTION, std::string> NEUTRINO_MASS_LOSS_PRESCRIPTION_LABEL = {
    { NEUTRINO_MASS_LOSS_PRESCRIPTION::FIXED_FRACTION, "FIXED_FRACTION" },
    { NEUTRINO_MASS_LOSS_PRESCRIPTION::FIXED_MASS,     "FIXED_MASS" }
};


// Neutron Star Equations of State
enum class NS_EOS: int { SSE, ARP3 };
const COMPASUnorderedMap<NS_EOS, std::string> NS_EOSLabel = {
    { NS_EOS::SSE,  "SSE" },
    { NS_EOS::ARP3, "ARP3" }
};


// Orbital Period Distributions
enum class ORBITAL_PERIOD_DISTRIBUTION: int { FLATINLOG };
const COMPASUnorderedMap<ORBITAL_PERIOD_DISTRIBUTION, std::string> ORBITAL_PERIOD_DISTRIBUTION_LABEL = {
    { ORBITAL_PERIOD_DISTRIBUTION::FLATINLOG, "FLATINLOG" },
};


// Pulsational Pair Instability Prescriptions
enum class PPI_PRESCRIPTION: int { COMPAS, STARTRACK, MARCHANT, FARMER };
const COMPASUnorderedMap<PPI_PRESCRIPTION, std::string> PPI_PRESCRIPTION_LABEL = {
    { PPI_PRESCRIPTION::COMPAS,    "COMPAS" },
    { PPI_PRESCRIPTION::STARTRACK, "STARTRACK" },
    { PPI_PRESCRIPTION::MARCHANT,  "MARCHANT" },
    { PPI_PRESCRIPTION::FARMER,    "FARMER" }
};


// Pulsar Birth Magnetic Field Distributions
enum class PULSAR_BIRTH_MAGNETIC_FIELD_DISTRIBUTION: int { ZERO, FIXED, FLATINLOG, UNIFORM, LOGNORMAL };
const COMPASUnorderedMap<PULSAR_BIRTH_MAGNETIC_FIELD_DISTRIBUTION, std::string> PULSAR_BIRTH_MAGNETIC_FIELD_DISTRIBUTION_LABEL = {
    { PULSAR_BIRTH_MAGNETIC_FIELD_DISTRIBUTION::ZERO,      "ZERO" },
    { PULSAR_BIRTH_MAGNETIC_FIELD_DISTRIBUTION::FIXED,     "FIXED" },
    { PULSAR_BIRTH_MAGNETIC_FIELD_DISTRIBUTION::FLATINLOG, "FLATINLOG" },
    { PULSAR_BIRTH_MAGNETIC_FIELD_DISTRIBUTION::UNIFORM,   "UNIFORM" },
    { PULSAR_BIRTH_MAGNETIC_FIELD_DISTRIBUTION::LOGNORMAL, "LOGNORMAL" }
};


// Pulsar Birth Spin Period Distributions
enum class PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION: int { ZERO, FIXED, UNIFORM, NORMAL };
const COMPASUnorderedMap<PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION, std::string> PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION_LABEL = {
    { PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION::ZERO,    "ZERO" },
    { PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION::FIXED,   "FIXED" },
    { PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION::UNIFORM, "UNIFORM" },
    { PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION::NORMAL,  "NORMAL" }
};


// Remnant Mass Prescriptions
enum class REMNANT_MASS_PRESCRIPTION: int { HURLEY2000, BELCZYNSKI2002, FRYER2012, FRYER2022, MULLER2016, MULLERMANDEL, SCHNEIDER2020, SCHNEIDER2020ALT};
const COMPASUnorderedMap<REMNANT_MASS_PRESCRIPTION, std::string> REMNANT_MASS_PRESCRIPTION_LABEL = {
    { REMNANT_MASS_PRESCRIPTION::HURLEY2000,       "HURLEY2000" },
    { REMNANT_MASS_PRESCRIPTION::BELCZYNSKI2002,   "BELCZYNSKI2002" },
    { REMNANT_MASS_PRESCRIPTION::FRYER2012,        "FRYER2012" },
    { REMNANT_MASS_PRESCRIPTION::FRYER2022,        "FRYER2022" },
    { REMNANT_MASS_PRESCRIPTION::MULLER2016,       "MULLER2016" },
    { REMNANT_MASS_PRESCRIPTION::MULLERMANDEL,     "MULLERMANDEL" },
    { REMNANT_MASS_PRESCRIPTION::SCHNEIDER2020,    "SCHNEIDER2020" },
    { REMNANT_MASS_PRESCRIPTION::SCHNEIDER2020ALT, "SCHNEIDER2020ALT" }
};


// Rotational Velocity Distribution options
enum class ROTATIONAL_VELOCITY_DISTRIBUTION: int { ZERO, HURLEY, VLTFLAMES };
const COMPASUnorderedMap<ROTATIONAL_VELOCITY_DISTRIBUTION, std::string> ROTATIONAL_VELOCITY_DISTRIBUTION_LABEL = {
    { ROTATIONAL_VELOCITY_DISTRIBUTION::ZERO,      "ZERO" },
    { ROTATIONAL_VELOCITY_DISTRIBUTION::HURLEY,    "HURLEY" },
    { ROTATIONAL_VELOCITY_DISTRIBUTION::VLTFLAMES, "VLTFLAMES" }
};


// Semi-major Axis Distributions
enum class SEMI_MAJOR_AXIS_DISTRIBUTION: int { FLATINLOG, DUQUENNOYMAYOR1991, SANA2012 };
const COMPASUnorderedMap<SEMI_MAJOR_AXIS_DISTRIBUTION, std::string> SEMI_MAJOR_AXIS_DISTRIBUTION_LABEL = {
    { SEMI_MAJOR_AXIS_DISTRIBUTION::FLATINLOG,          "FLATINLOG" },
    { SEMI_MAJOR_AXIS_DISTRIBUTION::DUQUENNOYMAYOR1991, "DUQUENNOYMAYOR1991" },
    { SEMI_MAJOR_AXIS_DISTRIBUTION::SANA2012,           "SANA2012" }
};

// Supernova Engines (Fryer 2012)
enum class SN_ENGINE: int { RAPID, DELAYED };
const COMPASUnorderedMap<SN_ENGINE, std::string> SN_ENGINE_LABEL = {
    { SN_ENGINE::RAPID,   "RAPID" },
    { SN_ENGINE::DELAYED, "DELAYED" }
};


// Supernova events/states
//
// The values here for SN_EVENT are powers of 2 so that they can be used in a bit map
// and manipulated with bit-wise logical operators
//
// Ordinarily we might expect that an SN event could be only one of CCSN, ECSN, PISN, PPISN, USSN, AIC, SNIA, or HeSD
// Note that the CCSN value here replaces the SN value in the legacy code
// The legacy code implemented these values as boolean flags, and the SN flag was always set when
// the uSSN flag was set (but not the converse).  In the legacy code when the ECSN flag was set 
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
//    SN_EVENT::CCSN  iff CCSN  bit is set and USSN bit is not set
//    SN_EVENT::ECSN  iff ECSN  bit is set
//    SN_EVENT::PISN  iff PISN  bit is set
//    SN_EVENT::PPISN iff PPISN bit is set
//    SN_EVENT::USSN  iff USSN  bit is set
//    SN_EVENT::AIC   iff AIC   bit is set
//    SN_EVENT::SNIA  iff SNIA  bit is set and HeSD bit is not set
//    SN_EVENT::HeSD  iff HeSD  bit is set
//    SN_EVENT::NONE  otherwise
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
};
ENABLE_BITMASK_OPERATORS(SN_EVENT);

const COMPASUnorderedMap<SN_EVENT, std::string> SN_EVENT_LABEL = {
    { SN_EVENT::NONE,  "No Supernova" },
    { SN_EVENT::CCSN,  "Core Collapse Supernova" },
    { SN_EVENT::ECSN,  "Electron Capture Supernova" },
    { SN_EVENT::PISN,  "Pair Instability Supernova" },
    { SN_EVENT::PPISN, "Pulsational Pair Instability Supernova" },
    { SN_EVENT::USSN,  "Ultra Stripped Supernova" },
    { SN_EVENT::AIC,   "Accretion-Induced Collapse" }, 
    { SN_EVENT::SNIA,  "Supernova Type Ia" }, 
    { SN_EVENT::HeSD,  "Helium-shell detonation" }, 
};


// Supernova types
enum class SN_STATE: int { NONE, STAR1, STAR2, BOTH };
const COMPASUnorderedMap<SN_STATE, std::string> SN_STATE_LABEL = {
    { SN_STATE::NONE,  "No Supernova" },
    { SN_STATE::STAR1, "Star1 only" },
    { SN_STATE::STAR2, "Star2 only" },
    { SN_STATE::BOTH,  "Both stars" }
};

// enum class ACCRETION_REGIME
// Symbolic names for WD accretion regimes
enum class ACCRETION_REGIME: int {
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
    HELIUM_WHITE_DWARF_HYDROGEN_ACCUMULATION,
    NONE
};

const COMPASUnorderedMap<ACCRETION_REGIME, std::string> ACCRETION_REGIME_LABEL = {
    { ACCRETION_REGIME::NONE,                                        "No accretion regime" },
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


// enum class L_CONSTANTS
// symbolic names for the Luminosity Constants
// these must be left as default values - their order can be changed with the caveat that the sentinel "COUNT" must stay at the end
// it's a bit of a hack, but it lets me calculate the number of L_CONSTANTS
enum class L_CONSTANTS: int { B_ALPHA_L, B_BETA_L, B_DELTA_L, COUNT };

// enum class R_CONSTANTS
// symbolic names for the Radius Constants
// these must be left as default values - their order can be changed with the caveat that the sentinel "COUNT" must stay at the end
// it's a bit of a hack, but it lets me calculate the number of R_CONSTANTS
enum class R_CONSTANTS: int { B_ALPHA_R, C_ALPHA_R, B_BETA_R, C_BETA_R, B_DELTA_R, COUNT };


// enum class GAMMA_CONSTANTS
// symbolic names for the Gamma Constants
// these must be left as default values - their order can be changed with the caveat that the sentinel "COUNT" must stay at the end
// it's a bit of a hack, but it lets me calculate the number of GAMMA_CONSTANTS
enum class GAMMA_CONSTANTS: int { B_GAMMA, C_GAMMA, COUNT };

// For the enum classes that follow:
//
// Order or entries is not significant - the code should not rely on these being in any order
// Entry value is significant for some (as noted), but must be unique (default value is ordinal position - should leave as defaults)


// enum class STELLAR_TYPE (StellarTypes from Hurley et al. 2000)
// Symbolic names for stellar types
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
    CHEMICALLY_HOMOGENEOUS,                         //  16  JR: this is here to preserve the Hurley type numbers, but note that Hurley type number progression doesn't necessarily indicate class inheritance
    STAR,                                           //  17  JR: added this - star is created this way, then switches as required (down here so stellar types consistent with Hurley et al. 2000)
    BINARY_STAR,                                    //  18  JR: added this - mainly for diagnostics
    NONE                                            //  19  JR: added this - mainly for diagnostics
};


// labels for stellar types
// unordered_map - key is integer stellar type (from enum class STELLAR_TYPE above)
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

// (convenience) initializer list for "evolvable" stellar types
// i.e. not STAR, BINARY_STAR, or NONE
const std::initializer_list<STELLAR_TYPE> EVOLVABLE_TYPES = {
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

// (convenience) initializer list for MAIN SEQUENCE stars (does not include NAKED_HELIUM_STAR_MS)
const std::initializer_list<STELLAR_TYPE> MAIN_SEQUENCE = {
    STELLAR_TYPE::MS_LTE_07,
    STELLAR_TYPE::MS_GT_07,
    STELLAR_TYPE::CHEMICALLY_HOMOGENEOUS
};


// (convenience) initializer list for ALL MAIN SEQUENCE stars (includes NAKED_HELIUM_STAR_MS)
const std::initializer_list<STELLAR_TYPE> ALL_MAIN_SEQUENCE = {
    STELLAR_TYPE::MS_LTE_07,
    STELLAR_TYPE::MS_GT_07,
    STELLAR_TYPE::CHEMICALLY_HOMOGENEOUS,
    STELLAR_TYPE::NAKED_HELIUM_STAR_MS
};


// (convenience) initializer list for ALL HERTZSPRUNG GAP (includes NAKED_HELIUM_STAR_HERTZSPRUNG_GAP)
const std::initializer_list<STELLAR_TYPE> ALL_HERTZSPRUNG_GAP = {
    STELLAR_TYPE::HERTZSPRUNG_GAP,
    STELLAR_TYPE::NAKED_HELIUM_STAR_HERTZSPRUNG_GAP
};

// (convenience) initializer list for COMPACT OBJECTS
const std::initializer_list<STELLAR_TYPE> COMPACT_OBJECTS = {
    STELLAR_TYPE::HELIUM_WHITE_DWARF,
    STELLAR_TYPE::CARBON_OXYGEN_WHITE_DWARF,
    STELLAR_TYPE::OXYGEN_NEON_WHITE_DWARF,
    STELLAR_TYPE::NEUTRON_STAR,
    STELLAR_TYPE::BLACK_HOLE,
    STELLAR_TYPE::MASSLESS_REMNANT
};

// (convenience) initializer list for GIANTS
const std::initializer_list<STELLAR_TYPE> GIANTS = {
    STELLAR_TYPE::FIRST_GIANT_BRANCH,
    STELLAR_TYPE::CORE_HELIUM_BURNING,
    STELLAR_TYPE::EARLY_ASYMPTOTIC_GIANT_BRANCH,
    STELLAR_TYPE::THERMALLY_PULSING_ASYMPTOTIC_GIANT_BRANCH
};


// (convenience) initializer list for WHITE DWARFS
const std::initializer_list<STELLAR_TYPE> WHITE_DWARFS = {
    STELLAR_TYPE::HELIUM_WHITE_DWARF,
    STELLAR_TYPE::CARBON_OXYGEN_WHITE_DWARF,
    STELLAR_TYPE::OXYGEN_NEON_WHITE_DWARF
};


// (convenience) initializer list for He rich stellar types
const std::initializer_list<STELLAR_TYPE> He_RICH_TYPES = {
    STELLAR_TYPE::NAKED_HELIUM_STAR_MS,
    STELLAR_TYPE::NAKED_HELIUM_STAR_HERTZSPRUNG_GAP,
    STELLAR_TYPE::NAKED_HELIUM_STAR_GIANT_BRANCH,
    STELLAR_TYPE::HELIUM_WHITE_DWARF
};


// White Dwarf Effective Baryon Number
// unordered_map - key is integer stellar type (from enum class ST above)
// Hurley et al. 2000, just after eq 90
const COMPASUnorderedMap<STELLAR_TYPE, double> WD_Baryon_Number = {
    {STELLAR_TYPE::HELIUM_WHITE_DWARF,         4.0},
    {STELLAR_TYPE::CARBON_OXYGEN_WHITE_DWARF, 15.0},
    {STELLAR_TYPE::OXYGEN_NEON_WHITE_DWARF,   17.0}
};


// enum class MASS_CUTOFF
// Symbolic names for mass cutoffs
// these must be left as default values - their order can be changed with the caveat that the sentinel "COUNT" must stay at the end
// it's a bit of a hack, but it lets me calculate the number of Timescales
enum class MASS_CUTOFF: int {
    MHook,                  // Mass above which hook appears on MS (in Msol)
    MHeF,                   // Maximum initial mass for which helium ignites degenerately in a Helium Flash (HeF)
    MFGB,                   // Maximum initial mass for which helium ignites on the First Giant Branch (FGB)
    MCHE,                   // Mass cutoff for calculation of initial angular frequency to determine if CHE occurs

    COUNT                   // Sentinel for entry count
};

// enum class MASS_LOSS_TYPE
// Symbolic names for mass loss rate type
enum class MASS_LOSS_TYPE: int {
    NONE,
    OB,
    VMS,
    GB,
    RSG,
    WR,
    LBV
};


// enum class TIMESCALE
// Symbolic names for timescales
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

                            // Early Asymptotic Giant Branch (EAGB)
                            // Why, then, are the following described as "FAGB"?  First AGB?
    tinf1_FAGB,             // Early Asymptotic Giant Branch tinf1
    tinf2_FAGB,             // Early Asymptotic Giant Branch tinf2
    tMx_FAGB,               // Early Asymptotic Giant Branch t(Mx)

                            // Thermally Pulsating Asymptotic Giant Branch (TPAGB)
                            // Why, then, are the following described as "SAGB"?    Second AGB?
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


// enum class GBP (Giant Branch Parameters - from Hurley et al. 2000)
// Symbolic names for Giant Branch Parameters
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
    Lx,                     // Luminosity parameter on the first giant branch (FGB) Lx as a function of the core mass (really a function of Mx).        JR: ADDED THIS
    Mx,                     // Crosover point of high-luminosity and low-luminosity in core mass - luminosity relation. Hurley et al. 2000, p552, eq38
    McBGB,                  // Core mass at BGB (Base of Giant Branch)
    McBAGB,                 // Core mass at BAGB (Base of Asymptotic Giant Branch).  Hurley et al. 2000, eq 66 (also see eq 75 and discussion)
    McDU,                   // Core mass at second dredge up.  Hurley et al. 2000, eq 69
    McSN,                   // Core mass at which the Asymptotic Giant Branch phase is terminated in a SN/loss of envelope                              JR: ADDED THIS

    COUNT                   // Sentinel for entry count
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
    SN_EVENT,
    SN_STATE,
    STRING_VECTOR,
    EVOLUTION_STATUS
};


// labels (long and short) for typenames
// unordered_map - key is integer typename (from enum class TYPENAME above)
const COMPASUnorderedMap<TYPENAME, std::tuple<std::string, std::string>> TYPENAME_LABEL = {
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
    { TYPENAME::SN_EVENT,         { "SN_EVENT",               "INT"            }},
    { TYPENAME::SN_STATE,         { "SN_STATE",               "INT"            }},
    { TYPENAME::STRING_VECTOR,    { "STRING_VECTOR",          "VECTOR<STRING>" }},
    { TYPENAME::EVOLUTION_STATUS, { "EVOLUTION_STATUS",       "INT"            }}
};


// enum class STRING_QUALIFIER
// Qualifier for typename STRING: FIXED_LENGTH or VARIABLE_LENGTH (used for printing to HDF5 files)
enum class STRING_QUALIFIER: int { NONE, FIXED_LENGTH, VARIABLE_LENGTH };


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
    SN_EVENT,
    SN_STATE,
    EVOLUTION_STATUS
> COMPAS_VARIABLE_TYPE;


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


// The #define below defines the STELLAR variables allowed for logfile record definition
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
    ID,                                              \
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
    NUCLEAR_TIMESCALE,                               \
    NUCLEAR_TIMESCALE_POST_COMMON_ENVELOPE,          \
    NUCLEAR_TIMESCALE_PRE_COMMON_ENVELOPE,           \
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


//  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//  !!!                                                                             !!!
//  !!!   Do not change the following map unless you are adding or deleting a new   !!!
//  !!!   property (or changing the name of an existing property for some reason)   !!!
//  !!!                                                                             !!!
//  !!!             This is not where header strings should be changed!             !!!
/// !!!       This is a lookup table for the logfile definitions file parser.       !!!
//  !!!                                                                             !!!
//  !!!   Header strings are in the following maps, and should be changed there:    !!!
//  !!!                                                                             !!!
//  !!!   std::map<ANY_STAR_PROPERTY, PROPERTY_DETAILS> ANY_STAR_PROPERTY_DETAIL    !!!
//  !!!   std::map<ANY_STAR_PROPERTY, PROPERTY_DETAILS> BINARY_PROPERTY_DETAIL      !!!
//  !!!                                                                             !!!
//  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

// map STAR PROPERTY to string identifying the property
// for lookup by the printing functions
// this map serves as the lookup for: STAR_PROPERTY, STAR_1_PROPERTY, STAR_2_PROPERTY, SUPERNOVA_PROPERTY, COMPANION_PROPERTY and ANY_STAR_PROPERTY
//
// Properties only need to be here if they are required to be available for 
// printing in the logfiles - all keys present here should also be in the STAR_PROPERTIES #define
// and ANY_STAR_PROPERTY_DETAIL
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
    { STAR_PROPERTY::ID,                                              "ID" },
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
    { STAR_PROPERTY::NUCLEAR_TIMESCALE,                               "NUCLEAR_TIMESCALE" },
    { STAR_PROPERTY::NUCLEAR_TIMESCALE_POST_COMMON_ENVELOPE,          "NUCLEAR_TIMESCALE_POST_COMMON_ENVELOPE" },
    { STAR_PROPERTY::NUCLEAR_TIMESCALE_PRE_COMMON_ENVELOPE,           "NUCLEAR_TIMESCALE_PRE_COMMON_ENVELOPE" },
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
    BE_BINARY_CURRENT_COMPANION_LUMINOSITY,
    BE_BINARY_CURRENT_COMPANION_MASS,
    BE_BINARY_CURRENT_COMPANION_RADIUS,
    BE_BINARY_CURRENT_COMPANION_TEFF,
    BE_BINARY_CURRENT_DT,
    BE_BINARY_CURRENT_ECCENTRICITY,
    BE_BINARY_CURRENT_ID,
    BE_BINARY_CURRENT_NS_MASS,
    BE_BINARY_CURRENT_SEMI_MAJOR_AXIS,
    BE_BINARY_CURRENT_TOTAL_TIME,
    CIRCULARIZATION_TIMESCALE,
    COMMON_ENVELOPE_AT_LEAST_ONCE,
    COMMON_ENVELOPE_EVENT_COUNT,
    DIMENSIONLESS_KICK_MAGNITUDE,
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
    { BINARY_PROPERTY::BE_BINARY_CURRENT_COMPANION_LUMINOSITY,             "BE_BINARY_CURRENT_COMPANION_LUMINOSITY" },
    { BINARY_PROPERTY::BE_BINARY_CURRENT_COMPANION_MASS,                   "BE_BINARY_CURRENT_COMPANION_MASS" },
    { BINARY_PROPERTY::BE_BINARY_CURRENT_COMPANION_RADIUS,                 "BE_BINARY_CURRENT_COMPANION_RADIUS" },
    { BINARY_PROPERTY::BE_BINARY_CURRENT_COMPANION_TEFF,                   "BE_BINARY_CURRENT_COMPANION_TEFF" },
    { BINARY_PROPERTY::BE_BINARY_CURRENT_DT,                               "BE_BINARY_CURRENT_DT" },
    { BINARY_PROPERTY::BE_BINARY_CURRENT_ECCENTRICITY,                     "BE_BINARY_CURRENT_ECCENTRICITY" },
    { BINARY_PROPERTY::BE_BINARY_CURRENT_ID,                               "BE_BINARY_CURRENT_ID" },
    { BINARY_PROPERTY::BE_BINARY_CURRENT_NS_MASS,                          "BE_BINARY_CURRENT_NS_MASS" },
    { BINARY_PROPERTY::BE_BINARY_CURRENT_SEMI_MAJOR_AXIS,                  "BE_BINARY_CURRENT_SEMI_MAJOR_AXIS" },
    { BINARY_PROPERTY::BE_BINARY_CURRENT_TOTAL_TIME,                       "BE_BINARY_CURRENT_TOTAL_TIME" },
    { BINARY_PROPERTY::CIRCULARIZATION_TIMESCALE,                          "CIRCULARIZATION_TIMESCALE" },
    { BINARY_PROPERTY::COMMON_ENVELOPE_AT_LEAST_ONCE,                      "COMMON_ENVELOPE_AT_LEAST_ONCE" },
    { BINARY_PROPERTY::COMMON_ENVELOPE_EVENT_COUNT,                        "COMMON_ENVELOPE_EVENT_COUNT" },
    { BINARY_PROPERTY::DIMENSIONLESS_KICK_MAGNITUDE,                       "DIMENSIONLESS_KICK_MAGNITUDE" },
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

    //BE_BINARIES,

    BLACK_HOLE_KICKS,
    
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
    ENVELOPE_STATE_PRESCRIPTION,
    EVOLUTION_MODE,

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
    LBV_PRESCRIPTION,
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
    MT_JLOSS_MACLEOD_LINEAR_FRACTION,
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

    ROTATIONAL_VELOCITY_DISTRIBUTION,
    ROTATIONAL_FREQUENCY,
    ROTATIONAL_FREQUENCY_1,
    ROTATIONAL_FREQUENCY_2,

    SEMI_MAJOR_AXIS,
    SEMI_MAJOR_AXIS_DISTRIBUTION,
    SEMI_MAJOR_AXIS_DISTRIBUTION_MAX,
    SEMI_MAJOR_AXIS_DISTRIBUTION_MIN,
    SEMI_MAJOR_AXIS_DISTRIBUTION_POWER,

    STELLAR_ZETA_PRESCRIPTION,

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

    //{ PROGRAM_OPTION::BE_BINARIES,                                    "BE_BINARIES" },

    { PROGRAM_OPTION::BLACK_HOLE_KICKS,                                 "BLACK_HOLE_KICKS" },
    
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
    { PROGRAM_OPTION::ENVELOPE_STATE_PRESCRIPTION,                      "ENVELOPE_STATE_PRESCRIPTION" },
    { PROGRAM_OPTION::EVOLUTION_MODE,                                   "EVOLUTION_MODE" },

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
    { PROGRAM_OPTION::LBV_PRESCRIPTION,                                 "LBV_PRESCRIPTION" },
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
    { PROGRAM_OPTION::MT_JLOSS_MACLEOD_LINEAR_FRACTION,                 "MT_JLOSS_MACLEOD_LINEAR_FRACTION" },
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

    { PROGRAM_OPTION::ROTATIONAL_VELOCITY_DISTRIBUTION,                 "ROTATIONAL_VELOCITY_DISTRIBUTION" },
    { PROGRAM_OPTION::ROTATIONAL_FREQUENCY,                             "ROTATIONAL_FREQUENCY" },
    { PROGRAM_OPTION::ROTATIONAL_FREQUENCY_1,                           "ROTATIONAL_FREQUENCY_1" },
    { PROGRAM_OPTION::ROTATIONAL_FREQUENCY_2,                           "ROTATIONAL_FREQUENCY_2" },
   
    { PROGRAM_OPTION::SEMI_MAJOR_AXIS,                                  "SEMI_MAJOR_AXIS" },
    { PROGRAM_OPTION::SEMI_MAJOR_AXIS_DISTRIBUTION,                     "SEMI_MAJOR_AXIS_DISTRIBUTION" },
    { PROGRAM_OPTION::SEMI_MAJOR_AXIS_DISTRIBUTION_MAX,                 "SEMI_MAJOR_AXIS_DISTRIBUTION_MAX" },
    { PROGRAM_OPTION::SEMI_MAJOR_AXIS_DISTRIBUTION_MIN,                 "SEMI_MAJOR_AXIS_DISTRIBUTION_MIN" },
    { PROGRAM_OPTION::SEMI_MAJOR_AXIS_DISTRIBUTION_POWER,               "SEMI_MAJOR_AXIS_DISTRIBUTION_POWER" },

    { PROGRAM_OPTION::STELLAR_ZETA_PRESCRIPTION,                        "STELLAR_ZETA_PRESCRIPTION" },

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
    { ANY_STAR_PROPERTY::ANGULAR_MOMENTUM,                                  { TYPENAME::DOUBLE,           "Ang_Momentum",                    "Msol AU^2 yr^-1",  14, 6 }},
    { ANY_STAR_PROPERTY::BINDING_ENERGY_AT_COMMON_ENVELOPE,                 { TYPENAME::DOUBLE,           "Binding_Energy@CE",               "erg",              14, 6 }},
    { ANY_STAR_PROPERTY::BINDING_ENERGY_FIXED,                              { TYPENAME::DOUBLE,           "BE_Fixed",                        "erg",              14, 6 }},
    { ANY_STAR_PROPERTY::BINDING_ENERGY_NANJING,                            { TYPENAME::DOUBLE,           "BE_Nanjing",                      "erg",              14, 6 }},
    { ANY_STAR_PROPERTY::BINDING_ENERGY_PRE_COMMON_ENVELOPE,                { TYPENAME::DOUBLE,           "Binding_Energy<CE",               "erg",              14, 6 }},
    { ANY_STAR_PROPERTY::BINDING_ENERGY_LOVERIDGE,                          { TYPENAME::DOUBLE,           "BE_Loveridge",                    "erg",              14, 6 }},
    { ANY_STAR_PROPERTY::BINDING_ENERGY_LOVERIDGE_WINDS,                    { TYPENAME::DOUBLE,           "BE_Loveridge_Winds",              "erg",              14, 6 }},
    { ANY_STAR_PROPERTY::BINDING_ENERGY_KRUCKOW,                            { TYPENAME::DOUBLE,           "BE_Kruckow",                      "erg",              14, 6 }},
    { ANY_STAR_PROPERTY::CHEMICALLY_HOMOGENEOUS_MAIN_SEQUENCE,              { TYPENAME::BOOL,             "CH_on_MS",                        "State",             0, 0 }},
    { ANY_STAR_PROPERTY::CO_CORE_MASS,                                      { TYPENAME::DOUBLE,           "Mass_CO_Core",                    "Msol",             14, 6 }},
    { ANY_STAR_PROPERTY::CO_CORE_MASS_AT_COMMON_ENVELOPE,                   { TYPENAME::DOUBLE,           "Mass_CO_Core@CE",                 "Msol",             14, 6 }},
    { ANY_STAR_PROPERTY::CO_CORE_MASS_AT_COMPACT_OBJECT_FORMATION,          { TYPENAME::DOUBLE,           "Mass_CO_Core@CO",                 "Msol",             14, 6 }},
    { ANY_STAR_PROPERTY::CORE_MASS,                                         { TYPENAME::DOUBLE,           "Mass_Core",                       "Msol",             14, 6 }},
    { ANY_STAR_PROPERTY::CORE_MASS_AT_COMMON_ENVELOPE,                      { TYPENAME::DOUBLE,           "Mass_Core@CE",                    "Msol",             14, 6 }},
    { ANY_STAR_PROPERTY::CORE_MASS_AT_COMPACT_OBJECT_FORMATION,             { TYPENAME::DOUBLE,           "Mass_Core@CO",                    "Msol",             14, 6 }},
    { ANY_STAR_PROPERTY::DRAWN_KICK_MAGNITUDE,                              { TYPENAME::DOUBLE,           "Drawn_Kick_Magnitude",            "kms^-1",           14, 6 }},
    { ANY_STAR_PROPERTY::DOMINANT_MASS_LOSS_RATE,                           { TYPENAME::INT,              "Dominant_Mass_Loss_Rate",         "-",                 4, 1 }},
    { ANY_STAR_PROPERTY::DT,                                                { TYPENAME::DOUBLE,           "dT",                              "Myr",              24, 15}},
    { ANY_STAR_PROPERTY::DYNAMICAL_TIMESCALE,                               { TYPENAME::DOUBLE,           "Tau_Dynamical",                   "Myr",              24, 15}},
    { ANY_STAR_PROPERTY::DYNAMICAL_TIMESCALE_POST_COMMON_ENVELOPE,          { TYPENAME::DOUBLE,           "Tau_Dynamical>CE",                "Myr",              24, 15}},
    { ANY_STAR_PROPERTY::DYNAMICAL_TIMESCALE_PRE_COMMON_ENVELOPE,           { TYPENAME::DOUBLE,           "Tau_Dynamical<CE",                "Myr",              24, 15}},
    { ANY_STAR_PROPERTY::ECCENTRIC_ANOMALY,                                 { TYPENAME::DOUBLE,           "Eccentric_Anomaly",               "-",                14, 6 }},
    { ANY_STAR_PROPERTY::ENV_MASS,                                          { TYPENAME::DOUBLE,           "Mass_Env",                        "Msol",             14, 6 }},
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
    { ANY_STAR_PROPERTY::FALLBACK_FRACTION,                                 { TYPENAME::DOUBLE,           "Fallback_Fraction",               "-",                14, 6 }},
    { ANY_STAR_PROPERTY::HE_CORE_MASS,                                      { TYPENAME::DOUBLE,           "Mass_He_Core",                    "Msol",             14, 6 }},
    { ANY_STAR_PROPERTY::HE_CORE_MASS_AT_COMMON_ENVELOPE,                   { TYPENAME::DOUBLE,           "Mass_He_Core@CE",                 "Msol",             14, 6 }},
    { ANY_STAR_PROPERTY::HE_CORE_MASS_AT_COMPACT_OBJECT_FORMATION,          { TYPENAME::DOUBLE,           "Mass_He_Core@CO",                 "Msol",             14, 6 }},
    { ANY_STAR_PROPERTY::ID,                                                { TYPENAME::OBJECT_ID,        "ID",                              "-",                12, 1 }},
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
    { ANY_STAR_PROPERTY::KICK_MAGNITUDE,                                    { TYPENAME::DOUBLE,           "Applied_Kick_Magnitude",          "kms^-1",           14, 6 }},
    { ANY_STAR_PROPERTY::LAMBDA_AT_COMMON_ENVELOPE,                         { TYPENAME::DOUBLE,           "Lambda@CE",                       "-",                14, 6 }},
    { ANY_STAR_PROPERTY::LAMBDA_DEWI,                                       { TYPENAME::DOUBLE,           "Lambda_Dewi",                     "-",                14, 6 }},
    { ANY_STAR_PROPERTY::LAMBDA_FIXED,                                      { TYPENAME::DOUBLE,           "Lambda_Fixed",                    "-",                14, 6 }},
    { ANY_STAR_PROPERTY::LAMBDA_KRUCKOW,                                    { TYPENAME::DOUBLE,           "Lambda_Kruckow",                  "-",                14, 6 }},
    { ANY_STAR_PROPERTY::LAMBDA_KRUCKOW_BOTTOM,                             { TYPENAME::DOUBLE,           "Lambda_Kruckow_Bottom",           "-",                14, 6 }},
    { ANY_STAR_PROPERTY::LAMBDA_KRUCKOW_MIDDLE,                             { TYPENAME::DOUBLE,           "Lambda_Kruckow_Middle",           "-",                14, 6 }},
    { ANY_STAR_PROPERTY::LAMBDA_KRUCKOW_TOP,                                { TYPENAME::DOUBLE,           "Lambda_Kruckow_Top",              "-",                14, 6 }},
    { ANY_STAR_PROPERTY::LAMBDA_LOVERIDGE,                                  { TYPENAME::DOUBLE,           "Lambda_Loveridge",                "-",                14, 6 }},
    { ANY_STAR_PROPERTY::LAMBDA_LOVERIDGE_WINDS,                            { TYPENAME::DOUBLE,           "Lambda_Loveridge_Winds",          "-",                14, 6 }},
    { ANY_STAR_PROPERTY::LAMBDA_NANJING,                                    { TYPENAME::DOUBLE,           "Lambda_Nanjing",                  "-",                14, 6 }},
    { ANY_STAR_PROPERTY::LBV_PHASE_FLAG,                                    { TYPENAME::BOOL,             "LBV_Phase_Flag",                  "Event",             0, 0 }},
    { ANY_STAR_PROPERTY::LUMINOSITY,                                        { TYPENAME::DOUBLE,           "Luminosity",                      "Lsol",             24, 15}},
    { ANY_STAR_PROPERTY::LUMINOSITY_POST_COMMON_ENVELOPE,                   { TYPENAME::DOUBLE,           "Luminosity>CE",                   "Lsol",             14, 6 }},
    { ANY_STAR_PROPERTY::LUMINOSITY_PRE_COMMON_ENVELOPE,                    { TYPENAME::DOUBLE,           "Luminosity<CE",                   "Lsol",             14, 6 }},
    { ANY_STAR_PROPERTY::MASS,                                              { TYPENAME::DOUBLE,           "Mass",                            "Msol",             24, 15}},
    { ANY_STAR_PROPERTY::MASS_0,                                            { TYPENAME::DOUBLE,           "Mass_0",                          "Msol",             14, 6 }},
    { ANY_STAR_PROPERTY::MASS_LOSS_DIFF,                                    { TYPENAME::DOUBLE,           "dmWinds",                         "Msol",             14, 6 }},
    { ANY_STAR_PROPERTY::MASS_TRANSFER_DIFF,                                { TYPENAME::DOUBLE,           "dmMT",                            "Msol",             14, 6 }},
    { ANY_STAR_PROPERTY::MASS_TRANSFER_DONOR_HISTORY,                       { TYPENAME::STRING,           "MT_Donor_Hist",                   "-",                16, 1 }}, 
    { ANY_STAR_PROPERTY::MDOT,                                              { TYPENAME::DOUBLE,           "Mdot",                            "Msol yr^-1",       14, 6 }},
    { ANY_STAR_PROPERTY::METALLICITY,                                       { TYPENAME::DOUBLE,           "Metallicity@ZAMS",                "-",                14, 6 }},
    { ANY_STAR_PROPERTY::MOMENT_OF_INERTIA,                                 { TYPENAME::DOUBLE,           "Moment_Of_Inertia",               "Msol Rsol^2",      14, 6 }},
    { ANY_STAR_PROPERTY::MZAMS,                                             { TYPENAME::DOUBLE,           "Mass@ZAMS",                       "Msol",             14, 6 }},
    { ANY_STAR_PROPERTY::NUCLEAR_TIMESCALE,                                 { TYPENAME::DOUBLE,           "Tau_Nuclear",                     "Myr",              24, 15}},
    { ANY_STAR_PROPERTY::NUCLEAR_TIMESCALE_POST_COMMON_ENVELOPE,            { TYPENAME::DOUBLE,           "Tau_Nuclear>CE",                  "Myr",              24, 15}},
    { ANY_STAR_PROPERTY::NUCLEAR_TIMESCALE_PRE_COMMON_ENVELOPE,             { TYPENAME::DOUBLE,           "Tau_Nuclear<CE",                  "Myr",              24, 15}},
    { ANY_STAR_PROPERTY::OMEGA,                                             { TYPENAME::DOUBLE,           "Omega",                           "Hz",               14, 6 }},
    { ANY_STAR_PROPERTY::OMEGA_BREAK,                                       { TYPENAME::DOUBLE,           "Omega_Break",                     "Hz",               14, 6 }},
    { ANY_STAR_PROPERTY::OMEGA_ZAMS,                                        { TYPENAME::DOUBLE,           "Omega@ZAMS",                      "Hz",               14, 6 }},
    { ANY_STAR_PROPERTY::ORBITAL_ENERGY_POST_SUPERNOVA,                     { TYPENAME::DOUBLE,           "Orbital_Energy>SN",               "Msol^2AU^-1",      14, 6 }},
    { ANY_STAR_PROPERTY::ORBITAL_ENERGY_PRE_SUPERNOVA,                      { TYPENAME::DOUBLE,           "Orbital_Energy<SN",               "Msol^2AU^-1",      14, 6 }},
    { ANY_STAR_PROPERTY::PULSAR_MAGNETIC_FIELD,                             { TYPENAME::DOUBLE,           "Pulsar_Mag_Field",                "Tesla",            14, 6 }},
    { ANY_STAR_PROPERTY::PULSAR_SPIN_DOWN_RATE,                             { TYPENAME::DOUBLE,           "Pulsar_Spin_Down",                "rad/s^2",          14, 6 }},
    { ANY_STAR_PROPERTY::PULSAR_BIRTH_PERIOD,                               { TYPENAME::DOUBLE,           "Pulsar_Birth_Period",             "s",                14, 6 }},
    { ANY_STAR_PROPERTY::PULSAR_BIRTH_SPIN_DOWN_RATE,                       { TYPENAME::DOUBLE,           "Pulsar_Birth_Spin_Down",          "s/s",              14, 6 }},
    { ANY_STAR_PROPERTY::PULSAR_SPIN_FREQUENCY,                             { TYPENAME::DOUBLE,           "Pulsar_Spin_Freq",                "rad/s",            14, 6 }},
    { ANY_STAR_PROPERTY::PULSAR_SPIN_PERIOD,                                { TYPENAME::DOUBLE,           "Pulsar_Spin_Period",              "ms",               14, 6 }},
    { ANY_STAR_PROPERTY::RADIAL_EXPANSION_TIMESCALE,                        { TYPENAME::DOUBLE,           "Tau_Radial",                      "Myr",              24, 15}},
    { ANY_STAR_PROPERTY::RADIAL_EXPANSION_TIMESCALE_POST_COMMON_ENVELOPE,   { TYPENAME::DOUBLE,           "Tau_Radial>CE",                   "Myr",              24, 15}},
    { ANY_STAR_PROPERTY::RADIAL_EXPANSION_TIMESCALE_PRE_COMMON_ENVELOPE,    { TYPENAME::DOUBLE,           "Tau_Radial<CE",                   "Myr",              24, 15}},
    { ANY_STAR_PROPERTY::RADIUS,                                            { TYPENAME::DOUBLE,           "Radius",                          "Rsol",             24, 15}},
    { ANY_STAR_PROPERTY::RANDOM_SEED,                                       { TYPENAME::ULONGINT,         "SEED",                            "-",                12, 1 }},
    { ANY_STAR_PROPERTY::RECYCLED_NEUTRON_STAR,                             { TYPENAME::BOOL,             "Recycled_NS",                     "Event",             0, 0 }},
    { ANY_STAR_PROPERTY::RLOF_ONTO_NS,                                      { TYPENAME::BOOL,             "RLOF->NS",                        "Event",             0, 0 }},
    { ANY_STAR_PROPERTY::RZAMS,                                             { TYPENAME::DOUBLE,           "Radius@ZAMS",                     "Rsol",             14, 6 }},
    { ANY_STAR_PROPERTY::SN_TYPE,                                           { TYPENAME::SN_EVENT,         "SN_Type",                         "-",                 4, 1 }},
    { ANY_STAR_PROPERTY::SPEED,                                             { TYPENAME::DOUBLE,           "ComponentSpeed",                  "kms^-1",           14, 6 }},
    { ANY_STAR_PROPERTY::STELLAR_TYPE,                                      { TYPENAME::STELLAR_TYPE,     "Stellar_Type",                    "-",                 4, 1 }},
    { ANY_STAR_PROPERTY::STELLAR_TYPE_NAME,                                 { TYPENAME::STRING,           "Stellar_Type",                    "-",                42, 1 }},
    { ANY_STAR_PROPERTY::STELLAR_TYPE_PREV,                                 { TYPENAME::STELLAR_TYPE,     "Stellar_Type_Prev",               "-",                 4, 1 }},
    { ANY_STAR_PROPERTY::STELLAR_TYPE_PREV_NAME,                            { TYPENAME::STRING,           "Stellar_Type_Prev",               "-",                42, 1 }},
    { ANY_STAR_PROPERTY::SUPERNOVA_KICK_MAGNITUDE_RANDOM_NUMBER,            { TYPENAME::DOUBLE,           "SN_Kick_Magnitude_Random_Number", "-",                14, 6 }},
    { ANY_STAR_PROPERTY::MEAN_ANOMALY,                                      { TYPENAME::DOUBLE,           "SN_Kick_Mean_Anomaly",            "-",                14, 6 }},
    { ANY_STAR_PROPERTY::SUPERNOVA_PHI,                                     { TYPENAME::DOUBLE,           "SN_Kick_Phi",                     "-",                14, 6 }},
    { ANY_STAR_PROPERTY::SUPERNOVA_THETA,                                   { TYPENAME::DOUBLE,           "SN_Kick_Theta",                   "-",                14, 6 }},
    { ANY_STAR_PROPERTY::TEMPERATURE,                                       { TYPENAME::DOUBLE,           "Teff",                            "K",                14, 6 }},
    { ANY_STAR_PROPERTY::TEMPERATURE_POST_COMMON_ENVELOPE,                  { TYPENAME::DOUBLE,           "Teff>CE",                         "K",                14, 6 }},
    { ANY_STAR_PROPERTY::TEMPERATURE_PRE_COMMON_ENVELOPE,                   { TYPENAME::DOUBLE,           "Teff<CE",                         "K",                14, 6 }},
    { ANY_STAR_PROPERTY::THERMAL_TIMESCALE,                                 { TYPENAME::DOUBLE,           "Tau_Thermal",                     "Myr",              24, 15}},
    { ANY_STAR_PROPERTY::THERMAL_TIMESCALE_POST_COMMON_ENVELOPE,            { TYPENAME::DOUBLE,           "Tau_Thermal>CE",                  "Myr",              24, 15}},
    { ANY_STAR_PROPERTY::THERMAL_TIMESCALE_PRE_COMMON_ENVELOPE,             { TYPENAME::DOUBLE,           "Tau_Thermal<CE",                  "Myr",              24, 15}},
    { ANY_STAR_PROPERTY::TIME,                                              { TYPENAME::DOUBLE,           "Time",                            "Myr",              24, 15 }},
    { ANY_STAR_PROPERTY::TIMESCALE_MS,                                      { TYPENAME::DOUBLE,           "tMS",                             "Myr",              24, 15}},
    { ANY_STAR_PROPERTY::TOTAL_MASS_AT_COMPACT_OBJECT_FORMATION,            { TYPENAME::DOUBLE,           "Mass_Total@CO",                   "Msol",             14, 6 }},
    { ANY_STAR_PROPERTY::TRUE_ANOMALY,                                      { TYPENAME::DOUBLE,           "True_Anomaly(psi)",               "-",                14, 6 }},
    { ANY_STAR_PROPERTY::TZAMS,                                             { TYPENAME::DOUBLE,           "Teff@ZAMS",                       "K",                14, 6 }},
    { ANY_STAR_PROPERTY::ZETA_HURLEY,                                       { TYPENAME::DOUBLE,           "Zeta_Hurley",                     "-",                14, 6 }},
    { ANY_STAR_PROPERTY::ZETA_HURLEY_HE,                                    { TYPENAME::DOUBLE,           "Zeta_Hurley_He",                  "-",                14, 6 }},
    { ANY_STAR_PROPERTY::ZETA_SOBERMAN,                                     { TYPENAME::DOUBLE,           "Zeta_Soberman",                   "-",                14, 6 }},
    { ANY_STAR_PROPERTY::ZETA_SOBERMAN_HE,                                  { TYPENAME::DOUBLE,           "Zeta_Soberman_He",                "-",                14, 6 }}
};

// map BINARY_PROPERTY_DETAIL
// Records the details of BINARY properties.  The BINARY properties are those that pertain
// to exclusively to a binary star - not the constituent stars that make up the binary
//
// Properties only need to be here if they are required to be available for printing in 
// the logfiles - all keys present here should also be in BINARY_PROPERTY and BINARY_PROPERTY_LABEL
const std::map<BINARY_PROPERTY, PROPERTY_DETAILS> BINARY_PROPERTY_DETAIL = {
    { BINARY_PROPERTY::BE_BINARY_CURRENT_COMPANION_LUMINOSITY,              { TYPENAME::DOUBLE,           "Companion_Lum",             "Lsol",             14, 6 }},
    { BINARY_PROPERTY::BE_BINARY_CURRENT_COMPANION_MASS,                    { TYPENAME::DOUBLE,           "Companion_Mass",            "Msol",             14, 6 }},
    { BINARY_PROPERTY::BE_BINARY_CURRENT_COMPANION_RADIUS,                  { TYPENAME::DOUBLE,           "Companion_Radius",          "Rsol",             14, 6 }},
    { BINARY_PROPERTY::BE_BINARY_CURRENT_COMPANION_TEFF,                    { TYPENAME::DOUBLE,           "Companion_Teff",            "K",                14, 6 }},
    { BINARY_PROPERTY::BE_BINARY_CURRENT_DT,                                { TYPENAME::DOUBLE,           "dT",                        "Myr",              24, 15}},
    { BINARY_PROPERTY::BE_BINARY_CURRENT_ECCENTRICITY,                      { TYPENAME::DOUBLE,           "Eccentricity",              "-",                14, 6 }},
    { BINARY_PROPERTY::BE_BINARY_CURRENT_ID,                                { TYPENAME::OBJECT_ID,        "ID",                        "-",                12, 1 }},
    { BINARY_PROPERTY::BE_BINARY_CURRENT_NS_MASS,                           { TYPENAME::DOUBLE,           "NS_Mass",                   "Msol",             14, 6 }},
    { BINARY_PROPERTY::BE_BINARY_CURRENT_SEMI_MAJOR_AXIS,                   { TYPENAME::DOUBLE,           "SemiMajorAxis",             "Rsol",             14, 6 }},
    { BINARY_PROPERTY::BE_BINARY_CURRENT_TOTAL_TIME,                        { TYPENAME::DOUBLE,           "Total_Time",                "Myr",              24, 15}},
    { BINARY_PROPERTY::CIRCULARIZATION_TIMESCALE,                           { TYPENAME::DOUBLE,           "Tau_Circ",                  "Myr",              24, 15}},
    { BINARY_PROPERTY::COMMON_ENVELOPE_AT_LEAST_ONCE,                       { TYPENAME::BOOL,             "CEE",                       "Event",             0, 0 }},
    { BINARY_PROPERTY::COMMON_ENVELOPE_EVENT_COUNT,                         { TYPENAME::UINT,             "CE_Event_Counter",          "Count",             6, 1 }},
    { BINARY_PROPERTY::DIMENSIONLESS_KICK_MAGNITUDE,                        { TYPENAME::DOUBLE,           "Kick_Magnitude(uK)",        "-",                14, 6 }},
    { BINARY_PROPERTY::DOUBLE_CORE_COMMON_ENVELOPE,                         { TYPENAME::BOOL,             "Double_Core_CE",            "Event",             0, 0 }},
    { BINARY_PROPERTY::DT,                                                  { TYPENAME::DOUBLE,           "dT",                        "Myr",              24, 15}},
    { BINARY_PROPERTY::ECCENTRICITY,                                        { TYPENAME::DOUBLE,           "Eccentricity",              "-",                14, 6 }},
    { BINARY_PROPERTY::ECCENTRICITY_AT_DCO_FORMATION,                       { TYPENAME::DOUBLE,           "Eccentricity@DCO",          "-",                14, 6 }},
    { BINARY_PROPERTY::ECCENTRICITY_INITIAL,                                { TYPENAME::DOUBLE,           "Eccentricity@ZAMS",         "-",                14, 6 }},
    { BINARY_PROPERTY::ECCENTRICITY_POST_COMMON_ENVELOPE,                   { TYPENAME::DOUBLE,           "Eccentricity>CE",           "-",                14, 6 }},
    { BINARY_PROPERTY::ECCENTRICITY_PRE_SUPERNOVA,                          { TYPENAME::DOUBLE,           "Eccentricity<SN",           "-",                14, 6 }},
    { BINARY_PROPERTY::ECCENTRICITY_PRE_COMMON_ENVELOPE,                    { TYPENAME::DOUBLE,           "Eccentricity<CE",           "-",                14, 6 }},
    { BINARY_PROPERTY::ERROR,                                               { TYPENAME::ERROR,            "Error",                     "-",                 4, 1 }},
    { BINARY_PROPERTY::EVOL_STATUS,                                         { TYPENAME::EVOLUTION_STATUS, "Evolution_Status",          "-",                 4, 1 }},
    { BINARY_PROPERTY::ID,                                                  { TYPENAME::OBJECT_ID,        "ID",                        "-",                12, 1 }},
    { BINARY_PROPERTY::IMMEDIATE_RLOF_POST_COMMON_ENVELOPE,                 { TYPENAME::BOOL,             "Immediate_RLOF>CE",         "Event",             0, 0 }},
    { BINARY_PROPERTY::MASS_1_POST_COMMON_ENVELOPE,                         { TYPENAME::DOUBLE,           "Mass(1)>CE",                "Msol",             14, 6 }},
    { BINARY_PROPERTY::MASS_1_PRE_COMMON_ENVELOPE,                          { TYPENAME::DOUBLE,           "Mass(1)<CE",                "Msol",             14, 6 }},
    { BINARY_PROPERTY::MASS_2_POST_COMMON_ENVELOPE,                         { TYPENAME::DOUBLE,           "Mass(2)>CE",                "Msol",             14, 6 }},
    { BINARY_PROPERTY::MASS_2_PRE_COMMON_ENVELOPE,                          { TYPENAME::DOUBLE,           "Mass(2)<CE",                "Msol",             14, 6 }},
    { BINARY_PROPERTY::MASS_ENV_1,                                          { TYPENAME::DOUBLE,           "Mass_Env(1)",               "Msol",             14, 6 }},
    { BINARY_PROPERTY::MASS_ENV_2,                                          { TYPENAME::DOUBLE,           "Mass_Env(2)",               "Msol",             14, 6 }},
    { BINARY_PROPERTY::MASSES_EQUILIBRATED,                                 { TYPENAME::BOOL,             "Equilibrated",              "Event",             0, 0 }},
    { BINARY_PROPERTY::MASSES_EQUILIBRATED_AT_BIRTH,                        { TYPENAME::BOOL,             "Equilibrated_At_Birth",     "Event",             0, 0 }},
    { BINARY_PROPERTY::MASS_TRANSFER_TRACKER_HISTORY,                       { TYPENAME::MT_TRACKING,      "MT_History",                "-",                 4, 1 }},
    { BINARY_PROPERTY::MERGES_IN_HUBBLE_TIME,                               { TYPENAME::BOOL,             "Merges_Hubble_Time",        "State",             0, 0 }},
    { BINARY_PROPERTY::OPTIMISTIC_COMMON_ENVELOPE,                          { TYPENAME::BOOL,             "Optimistic_CE",             "State",             0, 0 }},
    { BINARY_PROPERTY::ORBITAL_ANGULAR_VELOCITY,                            { TYPENAME::DOUBLE,           "Orbital_Angular_Velocity",  "kms^-1",           14, 6 }},
    { BINARY_PROPERTY::ORBITAL_VELOCITY_PRE_SUPERNOVA,                      { TYPENAME::DOUBLE,           "Orb_Velocity<SN",           "kms^-1",           14, 6 }},
    { BINARY_PROPERTY::RADIUS_1_POST_COMMON_ENVELOPE,                       { TYPENAME::DOUBLE,           "Radius(1)>CE",              "Rsol",             14, 6 }},
    { BINARY_PROPERTY::RADIUS_1_PRE_COMMON_ENVELOPE,                        { TYPENAME::DOUBLE,           "Radius(1)<CE",              "Rsol",             14, 6 }},
    { BINARY_PROPERTY::RADIUS_2_POST_COMMON_ENVELOPE,                       { TYPENAME::DOUBLE,           "Radius(2)>CE",              "Rsol",             14, 6 }},
    { BINARY_PROPERTY::RADIUS_2_PRE_COMMON_ENVELOPE,                        { TYPENAME::DOUBLE,           "Radius(2)<CE",              "Rsol",             14, 6 }},
    { BINARY_PROPERTY::RANDOM_SEED,                                         { TYPENAME::ULONGINT,         "SEED",                      "-",                12, 1 }},
    { BINARY_PROPERTY::RLOF_POST_MT_COMMON_ENVELOPE,                        { TYPENAME::BOOL,             "CEE>MT",                    "State",             0, 0 }},
    { BINARY_PROPERTY::RLOF_POST_MT_ECCENTRICITY,                           { TYPENAME::DOUBLE,           "Eccentricity>MT",           "-",                14, 6 }},
    { BINARY_PROPERTY::RLOF_POST_MT_EVENT_COUNTER,                          { TYPENAME::UINT,             "MT_Event_Counter",          "Count",             6, 1 }},
    { BINARY_PROPERTY::RLOF_POST_MT_ID,                                     { TYPENAME::OBJECT_ID,        "ID>MT",                     "-",                12, 1 }},
    { BINARY_PROPERTY::RLOF_POST_MT_SEMI_MAJOR_AXIS,                        { TYPENAME::DOUBLE,           "SemiMajorAxis>MT",          "Rsol",             14, 6 }},
    { BINARY_PROPERTY::RLOF_POST_MT_STAR1_MASS,                             { TYPENAME::DOUBLE,           "Mass(1)>MT",                "Msol",             14, 6 }},
    { BINARY_PROPERTY::RLOF_POST_MT_STAR2_MASS,                             { TYPENAME::DOUBLE,           "Mass(2)>MT",                "Msol",             14, 6 }},
    { BINARY_PROPERTY::RLOF_POST_MT_STAR1_RADIUS,                           { TYPENAME::DOUBLE,           "Radius(1)>MT",              "Rsol",             14, 6 }},
    { BINARY_PROPERTY::RLOF_POST_MT_STAR2_RADIUS,                           { TYPENAME::DOUBLE,           "Radius(2)>MT",              "Rsol",             14, 6 }},
    { BINARY_PROPERTY::RLOF_POST_MT_STAR1_RLOF,                             { TYPENAME::BOOL,             "RLOF(1)>MT",                "State",             0, 0 }},
    { BINARY_PROPERTY::RLOF_POST_MT_STAR2_RLOF,                             { TYPENAME::BOOL,             "RLOF(2)>MT",                "State",             0, 0 }},
    { BINARY_PROPERTY::RLOF_POST_MT_STAR1_STELLAR_TYPE,                     { TYPENAME::STELLAR_TYPE,     "Stellar_Type(1)>MT",         "-",                4, 1 }},
    { BINARY_PROPERTY::RLOF_POST_MT_STAR1_STELLAR_TYPE_NAME,                { TYPENAME::STRING,           "Stellar_Type(1)>MT",         "-",               42, 1 }},
    { BINARY_PROPERTY::RLOF_POST_MT_STAR2_STELLAR_TYPE,                     { TYPENAME::STELLAR_TYPE,     "Stellar_Type(2)>MT",         "-",                4, 1 }},
    { BINARY_PROPERTY::RLOF_POST_MT_STAR2_STELLAR_TYPE_NAME,                { TYPENAME::STRING,           "Stellar_Type(2)>MT",         "-",               42, 1 }},
    { BINARY_PROPERTY::RLOF_POST_STEP_STAR_TO_ROCHE_LOBE_RADIUS_RATIO_1,    { TYPENAME::DOUBLE,           "Radius(1)|RL>step",          "-",               14, 6 }},
    { BINARY_PROPERTY::RLOF_POST_STEP_STAR_TO_ROCHE_LOBE_RADIUS_RATIO_2,    { TYPENAME::DOUBLE,           "Radius(2)|RL>step",          "-",               14, 6 }},
    { BINARY_PROPERTY::RLOF_PRE_MT_ECCENTRICITY,                            { TYPENAME::DOUBLE,           "Eccentricity<MT",            "-",               14, 6 }},
    { BINARY_PROPERTY::RLOF_PRE_MT_SEMI_MAJOR_AXIS,                         { TYPENAME::DOUBLE,           "SemiMajorAxis<MT",           "Rsol",            14, 6 }},
    { BINARY_PROPERTY::RLOF_PRE_MT_STAR1_MASS,                              { TYPENAME::DOUBLE,           "Mass(1)<MT",                 "Msol",            14, 6 }},
    { BINARY_PROPERTY::RLOF_PRE_MT_STAR2_MASS,                              { TYPENAME::DOUBLE,           "Mass(2)<MT",                 "Msol",            14, 6 }},
    { BINARY_PROPERTY::RLOF_PRE_MT_STAR1_RADIUS,                            { TYPENAME::DOUBLE,           "Radius(1)<MT",               "Rsol",            14, 6 }},
    { BINARY_PROPERTY::RLOF_PRE_MT_STAR2_RADIUS,                            { TYPENAME::DOUBLE,           "Radius(2)<MT",               "Rsol",            14, 6 }},
    { BINARY_PROPERTY::RLOF_PRE_MT_STAR1_RLOF,                              { TYPENAME::BOOL,             "RLOF(1)<MT",                 "Event",            0, 0 }},
    { BINARY_PROPERTY::RLOF_PRE_MT_STAR2_RLOF,                              { TYPENAME::BOOL,             "RLOF(2)<MT",                 "Event",            0, 0 }},
    { BINARY_PROPERTY::RLOF_PRE_MT_STAR1_STELLAR_TYPE,                      { TYPENAME::STELLAR_TYPE,     "Stellar_Type(1)<MT",         "-",                4, 1 }},
    { BINARY_PROPERTY::RLOF_PRE_MT_STAR1_STELLAR_TYPE_NAME,                 { TYPENAME::STRING,           "Stellar_Type(1)<MT",         "-",               42, 1 }},
    { BINARY_PROPERTY::RLOF_PRE_MT_STAR2_STELLAR_TYPE,                      { TYPENAME::STELLAR_TYPE,     "Stellar_Type(2)<MT",         "-",                4, 1 }},
    { BINARY_PROPERTY::RLOF_PRE_MT_STAR2_STELLAR_TYPE_NAME,                 { TYPENAME::STRING,           "Stellar_Type(2)<MT",         "-",               42, 1 }},
    { BINARY_PROPERTY::RLOF_PRE_STEP_STAR_TO_ROCHE_LOBE_RADIUS_RATIO_1,     { TYPENAME::DOUBLE,           "Radius(1)|RL<step",          "-",               14, 6 }},
    { BINARY_PROPERTY::RLOF_PRE_STEP_STAR_TO_ROCHE_LOBE_RADIUS_RATIO_2,     { TYPENAME::DOUBLE,           "Radius(2)|RL<step",          "-",               14, 6 }},
    { BINARY_PROPERTY::RLOF_SECONDARY_POST_COMMON_ENVELOPE,                 { TYPENAME::BOOL,             "RLOF_Secondary>CE",          "Event",            0, 0 }},
    { BINARY_PROPERTY::RLOF_TIME_POST_MT,                                   { TYPENAME::DOUBLE,           "Time>MT",                    "Myr",             24, 15}},
    { BINARY_PROPERTY::RLOF_TIME_PRE_MT,                                    { TYPENAME::DOUBLE,           "Time<MT",                    "Myr",             24, 15}},
    { BINARY_PROPERTY::ROCHE_LOBE_RADIUS_1,                                 { TYPENAME::DOUBLE,           "RocheLobe(1)",               "Rsol",            14, 6 }},
    { BINARY_PROPERTY::ROCHE_LOBE_RADIUS_1_POST_COMMON_ENVELOPE,            { TYPENAME::DOUBLE,           "RocheLobe(1)>CE",            "Rsol",            14, 6 }},
    { BINARY_PROPERTY::ROCHE_LOBE_RADIUS_1_PRE_COMMON_ENVELOPE,             { TYPENAME::DOUBLE,           "RocheLobe(1)<CE",            "Rsol",            14, 6 }},
    { BINARY_PROPERTY::ROCHE_LOBE_RADIUS_2,                                 { TYPENAME::DOUBLE,           "RocheLobe(2)",               "Rsol",            14, 6 }},
    { BINARY_PROPERTY::ROCHE_LOBE_RADIUS_2_POST_COMMON_ENVELOPE,            { TYPENAME::DOUBLE,           "RocheLobe(2)>CE",            "Rsol",            14, 6 }},
    { BINARY_PROPERTY::ROCHE_LOBE_RADIUS_2_PRE_COMMON_ENVELOPE,             { TYPENAME::DOUBLE,           "RocheLobe(2)<CE",            "Rsol",            14, 6 }},
    { BINARY_PROPERTY::STAR_TO_ROCHE_LOBE_RADIUS_RATIO_1,                   { TYPENAME::DOUBLE,           "Radius(1)|RL",               "-",               14, 6 }},
    { BINARY_PROPERTY::STAR_TO_ROCHE_LOBE_RADIUS_RATIO_2,                   { TYPENAME::DOUBLE,           "Radius(2)|RL",               "-",               14, 6 }},
    { BINARY_PROPERTY::SEMI_MAJOR_AXIS_AT_DCO_FORMATION,                    { TYPENAME::DOUBLE,           "SemiMajorAxis@DCO",          "AU",              14, 6 }},
    { BINARY_PROPERTY::SEMI_MAJOR_AXIS_INITIAL,                             { TYPENAME::DOUBLE,           "SemiMajorAxis@ZAMS",         "AU",              14, 6 }},
    { BINARY_PROPERTY::SEMI_MAJOR_AXIS_POST_COMMON_ENVELOPE,                { TYPENAME::DOUBLE,           "SemiMajorAxis>CE",           "Rsol",            14, 6 }},
    { BINARY_PROPERTY::SEMI_MAJOR_AXIS_PRE_COMMON_ENVELOPE,                 { TYPENAME::DOUBLE,           "SemiMajorAxis<CE",           "Rsol",            14, 6 }},
    { BINARY_PROPERTY::SEMI_MAJOR_AXIS_PRE_SUPERNOVA,                       { TYPENAME::DOUBLE,           "SemiMajorAxis<SN",           "AU",              14, 6 }},
    { BINARY_PROPERTY::SEMI_MAJOR_AXIS_PRE_SUPERNOVA_RSOL,                  { TYPENAME::DOUBLE,           "SemiMajorAxis<SN",           "Rsol",            14, 6 }},
    { BINARY_PROPERTY::SEMI_MAJOR_AXIS,                                     { TYPENAME::DOUBLE,           "SemiMajorAxis",              "AU",              14, 6 }},
    { BINARY_PROPERTY::SEMI_MAJOR_AXIS_RSOL,                                { TYPENAME::DOUBLE,           "SemiMajorAxis",              "Rsol",            14, 6 }},
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
    { BINARY_PROPERTY::SUPERNOVA_ORBIT_INCLINATION_ANGLE,                   { TYPENAME::DOUBLE,           "SN_Orbit_Inclination_Angle", "-",               14, 6 }},   
    { BINARY_PROPERTY::SUPERNOVA_ORBIT_INCLINATION_VECTOR_X,                { TYPENAME::DOUBLE,           "Orbital_AM_Vector>SN_X",     "-",               14, 6 }},   
    { BINARY_PROPERTY::SUPERNOVA_ORBIT_INCLINATION_VECTOR_Y,                { TYPENAME::DOUBLE,           "Orbital_AM_Vector>SN_Y",     "-",               14, 6 }},   
    { BINARY_PROPERTY::SUPERNOVA_ORBIT_INCLINATION_VECTOR_Z,                { TYPENAME::DOUBLE,           "Orbital_AM_Vector>SN_Z",     "-",               14, 6 }},   
    { BINARY_PROPERTY::SUPERNOVA_STATE,                                     { TYPENAME::SN_STATE,         "Supernova_State",            "State",            4, 1 }},
    { BINARY_PROPERTY::SYNCHRONIZATION_TIMESCALE,                           { TYPENAME::DOUBLE,           "Tau_Sync",                   "Myr",             24, 15}},
    { BINARY_PROPERTY::SYSTEMIC_SPEED,                                      { TYPENAME::DOUBLE,           "SystemicSpeed",              "kms^-1",          14, 6 }},
    { BINARY_PROPERTY::TIME,                                                { TYPENAME::DOUBLE,           "Time",                       "Myr",             24, 15}},
    { BINARY_PROPERTY::TIME_TO_COALESCENCE,                                 { TYPENAME::DOUBLE,           "Coalescence_Time",           "Myr",             24, 15}},
    { BINARY_PROPERTY::TOTAL_ANGULAR_MOMENTUM,                              { TYPENAME::DOUBLE,           "Ang_Momentum_Total",         "Msol AU^2 yr^-1", 14, 6 }},
    { BINARY_PROPERTY::TOTAL_ENERGY,                                        { TYPENAME::DOUBLE,           "Energy_Total",               "Msol AU^2 yr^-2", 14, 6 }},
    { BINARY_PROPERTY::UNBOUND,                                             { TYPENAME::BOOL,             "Unbound",                    "State",            0, 0 }},
    { BINARY_PROPERTY::ZETA_LOBE,                                           { TYPENAME::DOUBLE,           "Zeta_Lobe",                  "-",               14, 6 }},
    { BINARY_PROPERTY::ZETA_STAR,                                           { TYPENAME::DOUBLE,           "Zeta_Star",                  "-",               14, 6 }}
};

// map PROGRAM_OPTION_DETAIL
// Records the details of PROGRAM_OPTION properties.
//
// Options only need to be here if they are required to be available for printing in 
// the logfiles - all keys present here should also be in PROGRAM_OPTION and PROGRAM_OPTION_LABEL
const std::map<PROGRAM_OPTION, PROPERTY_DETAILS> PROGRAM_OPTION_DETAIL = {

    { PROGRAM_OPTION::ADD_OPTIONS_TO_SYSPARMS,                                  { TYPENAME::INT,        "Add_Options_To_SysParms",                "-",          4, 1 }},
    { PROGRAM_OPTION::ALLOW_NON_STRIPPED_ECSN,                                  { TYPENAME::BOOL,       "Allow_Non_Stripped_ECSN",                "Flag",       0, 0 }},
    { PROGRAM_OPTION::ALLOW_MS_STAR_TO_SURVIVE_COMMON_ENVELOPE,                 { TYPENAME::BOOL,       "Allow_MS_To_Survive_CE",                 "Flag",       0, 0 }},
    { PROGRAM_OPTION::ALLOW_RADIATIVE_ENVELOPE_STAR_TO_SURVIVE_COMMON_ENVELOPE, { TYPENAME::BOOL,       "Allow_Radiative_Envelope_To_Survive_CE", "Flag",       0, 0 }},
    { PROGRAM_OPTION::ALLOW_IMMEDIATE_RLOF_POST_CE_TO_SURVIVE_COMMON_ENVELOPE,  { TYPENAME::BOOL,       "Allow_Immediate_RLOF>CE_To_Survive_CE",  "Flag",       0, 0 }},
    { PROGRAM_OPTION::ALLOW_RLOF_AT_BIRTH,                                      { TYPENAME::BOOL,       "Allow_RLOF@Birth",                       "Flag",       0, 0 }},
    { PROGRAM_OPTION::ALLOW_TOUCHING_AT_BIRTH,                                  { TYPENAME::BOOL,       "Allow_Touching@Birth",                   "Flag",       0, 0 }},
    { PROGRAM_OPTION::ANG_MOM_CONSERVATION_DURING_CIRCULARISATION,              { TYPENAME::BOOL,       "Conserve_AngMom@Circ",                   "Flag",       0, 0 }},
    //{ PROGRAM_OPTION::BE_BINARIES,                                              { TYPENAME::BOOL,       "Be_Binaries",                            "Flag",       0, 0 }},

    { PROGRAM_OPTION::BLACK_HOLE_KICKS,                                         { TYPENAME::INT,        "BH_Kicks",                               "-",          4, 1 }},
    
    { PROGRAM_OPTION::CASE_BB_STABILITY_PRESCRIPTION,                           { TYPENAME::INT,        "BB_Mass_xFer_Stblty_Prscrptn",           "-",          4, 1 }},
    
    { PROGRAM_OPTION::CHECK_PHOTON_TIRING_LIMIT,                                { TYPENAME::BOOL,       "Check_Photon_Tiring_Limit",              "Flag",       0, 0 }},

    { PROGRAM_OPTION::CHE_MODE,                                                 { TYPENAME::INT,        "CHE_Mode",                               "-",          4, 1 }},
      
    { PROGRAM_OPTION::CIRCULARISE_BINARY_DURING_MT,                             { TYPENAME::BOOL,       "Circularise@MT",                         "Flag",       0, 0 }},

    { PROGRAM_OPTION::COMMON_ENVELOPE_ALPHA,                                    { TYPENAME::DOUBLE,     "CE_Alpha",                               "-",         14, 6 }},
    { PROGRAM_OPTION::COMMON_ENVELOPE_ALPHA_THERMAL,                            { TYPENAME::DOUBLE,     "CE_Alpha_Thermal",                       "-",         14, 6 }},
    { PROGRAM_OPTION::COMMON_ENVELOPE_FORMALISM,                                { TYPENAME::INT,        "CE_Formalism",                           "-",          4, 1 }},
    { PROGRAM_OPTION::COMMON_ENVELOPE_LAMBDA,                                   { TYPENAME::DOUBLE,     "CE_Lambda",                              "-",         14, 6 }},
    { PROGRAM_OPTION::COMMON_ENVELOPE_LAMBDA_MULTIPLIER,                        { TYPENAME::DOUBLE,     "CE_Lambda_Multiplier",                   "-",         14, 6 }},
    { PROGRAM_OPTION::COMMON_ENVELOPE_LAMBDA_PRESCRIPTION,                      { TYPENAME::INT,        "CE_Lambda_Prscrptn",                     "-",          4, 1 }},
    { PROGRAM_OPTION::COMMON_ENVELOPE_MASS_ACCRETION_CONSTANT,                  { TYPENAME::DOUBLE,     "CE_Mass_Accr_Constant",                  "-",         14, 6 }},
    { PROGRAM_OPTION::COMMON_ENVELOPE_MASS_ACCRETION_MAX,                       { TYPENAME::DOUBLE,     "CE_Mass_Accr_Max",                       "Msol",      14, 6 }},
    { PROGRAM_OPTION::COMMON_ENVELOPE_MASS_ACCRETION_MIN,                       { TYPENAME::DOUBLE,     "CE_Mass_Accr_Min",                       "Msol",      14, 6 }},
    { PROGRAM_OPTION::COMMON_ENVELOPE_MASS_ACCRETION_PRESCRIPTION,              { TYPENAME::INT,        "CE_Mass_Accr_Prscrptn",                  "-",          4, 1 }},
    { PROGRAM_OPTION::COMMON_ENVELOPE_RECOMBINATION_ENERGY_DENSITY,             { TYPENAME::DOUBLE,     "CE_Recomb_Enrgy_Dnsty",                  "erg g^-1",  14, 6 }},
    { PROGRAM_OPTION::COMMON_ENVELOPE_SLOPE_KRUCKOW,                            { TYPENAME::DOUBLE,     "CE_Slope_Kruckow",                       "-",         14, 6 }},

    { PROGRAM_OPTION::COOL_WIND_MASS_LOSS_MULTIPLIER,                           { TYPENAME::DOUBLE,     "Cool_WindMassLoss_Multipl",              "-",         14, 6 }},

    { PROGRAM_OPTION::ECCENTRICITY,                                             { TYPENAME::DOUBLE,     "Eccentricity",                           "-",         14, 6 }},
    { PROGRAM_OPTION::ECCENTRICITY_DISTRIBUTION,                                { TYPENAME::INT,        "Eccentricity_Dstrbtn",                   "-",          4, 1 }},
    { PROGRAM_OPTION::ECCENTRICITY_DISTRIBUTION_MAX,                            { TYPENAME::DOUBLE,     "Eccentricity_Dstrbtn_Max",               "-",         14, 6 }},
    { PROGRAM_OPTION::ECCENTRICITY_DISTRIBUTION_MIN,                            { TYPENAME::DOUBLE,     "Eccentricity_Dstrbtn_Min",               "-",         14, 6 }},
    { PROGRAM_OPTION::EDDINGTON_ACCRETION_FACTOR,                               { TYPENAME::DOUBLE,     "Eddington_Accr_Factor",                  "-",         14, 6 }},
    { PROGRAM_OPTION::ENVELOPE_STATE_PRESCRIPTION,                              { TYPENAME::INT,        "Envelope_State_Prscrptn",                "-",          4, 1 }},
    { PROGRAM_OPTION::EVOLUTION_MODE,                                           { TYPENAME::INT,        "Evolution_Mode",                         "Mode",       4, 1 }},

    { PROGRAM_OPTION::FRYER_SUPERNOVA_ENGINE,                                   { TYPENAME::INT,        "Fryer_SN_Engine",                        "-",          4, 1 }},

    { PROGRAM_OPTION::FRYER22_FMIX,                                             { TYPENAME::DOUBLE,     "Fryer22_mixing_fraction",                "-",         14, 6 }},
    { PROGRAM_OPTION::FRYER22_MCRIT,                                            { TYPENAME::DOUBLE,     "Fryer22_crit_COcore_Mass",               "Msol",      14, 6 }},

    { PROGRAM_OPTION::INITIAL_MASS,                                             { TYPENAME::DOUBLE,     "Initial_Mass",                           "Msol",      14, 6 }},
    { PROGRAM_OPTION::INITIAL_MASS_1,                                           { TYPENAME::DOUBLE,     "Initial_Mass(1)",                        "Msol",      14, 6 }},
    { PROGRAM_OPTION::INITIAL_MASS_2,                                           { TYPENAME::DOUBLE,     "Initial_Mass(2)",                        "Msol",      14, 6 }},

    { PROGRAM_OPTION::INITIAL_MASS_FUNCTION,                                    { TYPENAME::INT,        "Initial_Mass_Function",                  "-",          4, 1 }},
    { PROGRAM_OPTION::INITIAL_MASS_FUNCTION_MAX,                                { TYPENAME::DOUBLE,     "Initial_Mass_Func_Max",                  "Msol",      14, 6 }},
    { PROGRAM_OPTION::INITIAL_MASS_FUNCTION_MIN,                                { TYPENAME::DOUBLE,     "Initial_Mass_Func_Min",                  "Msol",      14, 6 }},
    { PROGRAM_OPTION::INITIAL_MASS_FUNCTIONPOWER,                               { TYPENAME::DOUBLE,     "Initial_Mass_Func_Power",                "-",         14, 6 }},

    { PROGRAM_OPTION::KICK_DIRECTION_DISTRIBUTION,                              { TYPENAME::INT,        "Kick_Direction_Dstrbtn",                 "-",          4, 1 }},
    { PROGRAM_OPTION::KICK_DIRECTION_POWER,                                     { TYPENAME::DOUBLE,     "Kick_Direction_Power",                   "-",         14, 6 }},
    { PROGRAM_OPTION::KICK_SCALING_FACTOR,                                      { TYPENAME::DOUBLE,     "Kick_Scaling_Factor",                    "-",         14, 6 }},
    { PROGRAM_OPTION::KICK_MAGNITUDE_DISTRIBUTION,                              { TYPENAME::INT,        "Kick_Magnitude_Dstrbtn",                 "-",          4, 1 }},
    { PROGRAM_OPTION::KICK_MAGNITUDE_DISTRIBUTION_MAXIMUM,                      { TYPENAME::DOUBLE,     "Kick_Magnitude_Dstrbtn_Max",             "-",         14, 6 }},

    { PROGRAM_OPTION::KICK_MAGNITUDE_DISTRIBUTION_SIGMA_CCSN_BH,                { TYPENAME::DOUBLE,     "Sigma_Kick_CCSN_BH",                     "kms^-1",    14, 6 }},
    { PROGRAM_OPTION::KICK_MAGNITUDE_DISTRIBUTION_SIGMA_CCSN_NS,                { TYPENAME::DOUBLE,     "Sigma_Kick_CCSN_NS",                     "kms^-1",    14, 6 }},
    { PROGRAM_OPTION::KICK_MAGNITUDE_DISTRIBUTION_SIGMA_FOR_ECSN,               { TYPENAME::DOUBLE,     "Sigma_Kick_ECSN",                        "kms^-1",    14, 6 }},
    { PROGRAM_OPTION::KICK_MAGNITUDE_DISTRIBUTION_SIGMA_FOR_USSN,               { TYPENAME::DOUBLE,     "Sigma_Kick_USSN",                        "kms^-1",    14, 6 }},

    { PROGRAM_OPTION::KICK_MAGNITUDE,                                           { TYPENAME::DOUBLE,     "Kick_Magnitude",                         "kms^-1",    14, 6 }},
    { PROGRAM_OPTION::KICK_MAGNITUDE_1,                                         { TYPENAME::DOUBLE,     "Kick_Magnitude(1)",                      "kms^-1",    14, 6 }},
    { PROGRAM_OPTION::KICK_MAGNITUDE_2,                                         { TYPENAME::DOUBLE,     "Kick_Magnitude(2)",                      "kms^-1",    14, 6 }},

    { PROGRAM_OPTION::KICK_MAGNITUDE_RANDOM,                                    { TYPENAME::DOUBLE,     "Kick_Magnitude_Random",                  "-",         14, 6 }},
    { PROGRAM_OPTION::KICK_MAGNITUDE_RANDOM_1,                                  { TYPENAME::DOUBLE,     "Kick_Magnitude_Random(1)",               "-",         14, 6 }},
    { PROGRAM_OPTION::KICK_MAGNITUDE_RANDOM_2,                                  { TYPENAME::DOUBLE,     "Kick_Magnitude_Random(2)",               "-",         14, 6 }},

    { PROGRAM_OPTION::KICK_MEAN_ANOMALY_1,                                      { TYPENAME::DOUBLE,     "Kick_Mean_Anomaly(1)",                   "-",         14, 6 }},
    { PROGRAM_OPTION::KICK_MEAN_ANOMALY_2,                                      { TYPENAME::DOUBLE,     "Kick_Mean_Anomaly(2)",                   "-",         14, 6 }},
    { PROGRAM_OPTION::KICK_PHI_1,                                               { TYPENAME::DOUBLE,     "Kick_Phi(1)",                            "-",         14, 6 }},
    { PROGRAM_OPTION::KICK_PHI_2,                                               { TYPENAME::DOUBLE,     "Kick_Phi(2)",                            "-",         14, 6 }},
    { PROGRAM_OPTION::KICK_THETA_1,                                             { TYPENAME::DOUBLE,     "Kick_Theta(1)",                          "-",         14, 6 }},
    { PROGRAM_OPTION::KICK_THETA_2,                                             { TYPENAME::DOUBLE,     "Kick_Theta(2)",                          "-",         14, 6 }},

    { PROGRAM_OPTION::LBV_FACTOR,                                               { TYPENAME::DOUBLE,     "LBV_Factor",                             "-",         14, 6 }},
    { PROGRAM_OPTION::LBV_PRESCRIPTION,                                         { TYPENAME::INT,        "LBV_Mass_Loss_Prscrptn",                 "-",          4, 1 }},

    { PROGRAM_OPTION::MASS_LOSS_PRESCRIPTION,                                   { TYPENAME::INT,        "Mass_Loss_Prscrptn",                     "-",          4, 1 }},

    { PROGRAM_OPTION::MASS_RATIO,                                               { TYPENAME::DOUBLE,     "Mass_Ratio",                             "-",         14, 6 }},
    { PROGRAM_OPTION::MASS_RATIO_DISTRIBUTION,                                  { TYPENAME::INT,        "Mass_Ratio_Dstrbtn",                     "-",          4, 1 }},
    { PROGRAM_OPTION::MASS_RATIO_DISTRIBUTION_MAX,                              { TYPENAME::DOUBLE,     "Mass_Ratio_Dstrbtn_Max",                 "-",         14, 6 }},
    { PROGRAM_OPTION::MASS_RATIO_DISTRIBUTION_MIN,                              { TYPENAME::DOUBLE,     "Mass_Ratio_Dstrbtn_Min",                 "-",         14, 6 }},

    { PROGRAM_OPTION::MAXIMUM_EVOLUTION_TIME,                                   { TYPENAME::DOUBLE,     "Max_Evolution_Time",                     "Myr",       24, 15}},
    { PROGRAM_OPTION::MAXIMUM_DONOR_MASS,                                       { TYPENAME::DOUBLE,     "Max_Donor_Mass",                         "Msol",      14, 6 }},
    { PROGRAM_OPTION::MAXIMUM_NEUTRON_STAR_MASS,                                { TYPENAME::DOUBLE,     "Max_NS_Mass",                            "Msol",      14, 6 }},
    { PROGRAM_OPTION::MAXIMUM_TIMESTEPS,                                        { TYPENAME::ULONGINT,   "Max_Timesteps",                          "Count",     10, 1 }},

    { PROGRAM_OPTION::MCBUR1,                                                   { TYPENAME::DOUBLE,     "MCBUR1",                                 "Msol",      14, 6 }},

    { PROGRAM_OPTION::METALLICITY,                                              { TYPENAME::DOUBLE,     "Metallicity",                            "-",         14, 6 }},
    { PROGRAM_OPTION::METALLICITY_DISTRIBUTION,                                 { TYPENAME::INT,        "Metallicity_Dstrbtn",                    "-",          4, 1 }},
    { PROGRAM_OPTION::METALLICITY_DISTRIBUTION_MAX,                             { TYPENAME::DOUBLE,     "Metallicity_Dstrbtn_Max",                "-",         14, 6 }},
    { PROGRAM_OPTION::METALLICITY_DISTRIBUTION_MIN,                             { TYPENAME::DOUBLE,     "Metallicity_Dstrbtn_Min",                "-",         14, 6 }},

    { PROGRAM_OPTION::MINIMUM_MASS_SECONDARY,                                   { TYPENAME::DOUBLE,     "Min_Secondary_Mass",                     "Msol",      14, 6 }},

    { PROGRAM_OPTION::MT_ACCRETION_EFFICIENCY_PRESCRIPTION,                     { TYPENAME::INT,        "MT_Acc_Efficiency_Prscrptn",             "-",          4, 1 }},
    { PROGRAM_OPTION::MT_ANG_MOM_LOSS_PRESCRIPTION,                             { TYPENAME::INT,        "MT_AngMom_Loss_Prscrptn",                "-",          4, 1 }},
    { PROGRAM_OPTION::MT_THERMAL_LIMIT_C,                                       { TYPENAME::DOUBLE,     "MT_Thermal_Limit_C",                     "-",         14, 6 }},

    { PROGRAM_OPTION::MT_CRIT_MR_MS_LOW_MASS_DEGENERATE_ACCRETOR,               { TYPENAME::DOUBLE,     "MT_Crit_MR_MS_Low_Mass_Deg_Acc",         "-",         14, 6 }},
    { PROGRAM_OPTION::MT_CRIT_MR_MS_LOW_MASS_NON_DEGENERATE_ACCRETOR,           { TYPENAME::DOUBLE,     "MT_Crit_MR_MS_Low_Mass_NonDeg_Acc",      "-",         14, 6 }},
    { PROGRAM_OPTION::MT_CRIT_MR_MS_HIGH_MASS_DEGENERATE_ACCRETOR,              { TYPENAME::DOUBLE,     "MT_Crit_MR_MS_High_Mass_Deg_Acc",        "-",         14, 6 }},
    { PROGRAM_OPTION::MT_CRIT_MR_MS_HIGH_MASS_NON_DEGENERATE_ACCRETOR,          { TYPENAME::DOUBLE,     "MT_Crit_MR_MS_High_Mass_NonDeg_Acc",     "-",         14, 6 }},
    { PROGRAM_OPTION::MT_CRIT_MR_GIANT_DEGENERATE_ACCRETOR,                     { TYPENAME::DOUBLE,     "MT_Crit_MR_Giant_Deg_Acc",               "-",         14, 6 }},
    { PROGRAM_OPTION::MT_CRIT_MR_GIANT_NON_DEGENERATE_ACCRETOR,                 { TYPENAME::DOUBLE,     "MT_Crit_MR_Giant_NonDeg_Acc",            "-",         14, 6 }},
    { PROGRAM_OPTION::MT_CRIT_MR_HG_DEGENERATE_ACCRETOR,                        { TYPENAME::DOUBLE,     "MT_Crit_MR_HG_Deg_Acc",                  "-",         14, 6 }},

    { PROGRAM_OPTION::MT_CRIT_MR_HG_NON_DEGENERATE_ACCRETOR,                    { TYPENAME::DOUBLE,     "MT_Crit_MR_HG_NonDeg_Acc",               "-",         14, 6 }},
    { PROGRAM_OPTION::MT_CRIT_MR_HE_GIANT_DEGENERATE_ACCRETOR,                  { TYPENAME::DOUBLE,     "MT_Crit_MR_HE_Giant_Deg_Acc",            "-",         14, 6 }},
    { PROGRAM_OPTION::MT_CRIT_MR_HE_GIANT_NON_DEGENERATE_ACCRETOR,              { TYPENAME::DOUBLE,     "MT_Crit_MR_HE_Giant_NonDeg_Acc",         "-",         14, 6 }},
    { PROGRAM_OPTION::MT_CRIT_MR_HE_HG_DEGENERATE_ACCRETOR,                     { TYPENAME::DOUBLE,     "MT_Crit_MR_HE_HG_Deg_Acc",               "-",         14, 6 }},
    { PROGRAM_OPTION::MT_CRIT_MR_HE_HG_NON_DEGENERATE_ACCRETOR,                 { TYPENAME::DOUBLE,     "MT_Crit_MR_HE_HG_NonDeg_Acc",            "-",         14, 6 }},
    { PROGRAM_OPTION::MT_CRIT_MR_HE_MS_DEGENERATE_ACCRETOR,                     { TYPENAME::DOUBLE,     "MT_Crit_MR_HE_MS_Deg_Acc",               "-",         14, 6 }},
    { PROGRAM_OPTION::MT_CRIT_MR_HE_MS_NON_DEGENERATE_ACCRETOR,                 { TYPENAME::DOUBLE,     "MT_Crit_MR_HE_MS_NonDeg_Acc",            "-",         14, 6 }},
    { PROGRAM_OPTION::MT_CRIT_MR_WD_DEGENERATE_ACCRETOR,                        { TYPENAME::DOUBLE,     "MT_Crit_MR_WD_Deg_Acc",                  "-",         14, 6 }},
    { PROGRAM_OPTION::MT_CRIT_MR_WD_NONDEGENERATE_ACCRETOR,                     { TYPENAME::DOUBLE,     "MT_Crit_MR_WD_NonDeg_Acc",               "-",         14, 6 }},
    
    { PROGRAM_OPTION::MT_FRACTION_ACCRETED,                                     { TYPENAME::DOUBLE,     "MT_Fraction_Accreted",                   "-",         14, 6 }},
    { PROGRAM_OPTION::MT_JLOSS,                                                 { TYPENAME::DOUBLE,     "MT_JLoss",                               "-",         14, 6 }},
    { PROGRAM_OPTION::MT_JLOSS_MACLEOD_LINEAR_FRACTION,                         { TYPENAME::DOUBLE,     "MT_JLoss_Macleod_Linear_Frac",           "-",         14, 6 }},
    { PROGRAM_OPTION::MT_REJUVENATION_PRESCRIPTION,                             { TYPENAME::INT,        "MT_Rejuvenation_Prscrptn",               "-",          4, 1 }},
    { PROGRAM_OPTION::MT_THERMALLY_LIMITED_VARIATION,                           { TYPENAME::INT,        "MT_Thermally_Lmtd_Variation",            "-",          4, 1 }},

    { PROGRAM_OPTION::MULLER_MANDEL_KICK_MULTIPLIER_BH,                         { TYPENAME::DOUBLE,     "MM_Kick_Multiplier_BH",                  "-",         14, 6 }},
    { PROGRAM_OPTION::MULLER_MANDEL_KICK_MULTIPLIER_NS,                         { TYPENAME::DOUBLE,     "MM_Kick_Multiplier_NS",                  "-",         14, 6 }},
    { PROGRAM_OPTION::MULLER_MANDEL_SIGMA_KICK,                                 { TYPENAME::DOUBLE,     "MM_Sigma_Kick",                          "-",         14, 6 }},
    
    { PROGRAM_OPTION::NEUTRINO_MASS_LOSS_ASSUMPTION_BH,                         { TYPENAME::INT,        "Neutrino_Mass_Loss_Assmptn",             "-",          4, 1 }},
    { PROGRAM_OPTION::NEUTRINO_MASS_LOSS_VALUE_BH,                              { TYPENAME::DOUBLE,     "Neutrino_Mass_Loss_Value",               "-",         14, 6 }},

    { PROGRAM_OPTION::NS_EOS,                                                   { TYPENAME::INT,        "NS_EOS",                                 "-",          4, 1 }},

    { PROGRAM_OPTION::ORBITAL_PERIOD,                                           { TYPENAME::DOUBLE,     "Orbital_Period",                         "days",      14, 6 }},
    { PROGRAM_OPTION::ORBITAL_PERIOD_DISTRIBUTION,                              { TYPENAME::INT,        "Orbital_Period_Dstrbtn",                 "-",          4, 1 }},
    { PROGRAM_OPTION::ORBITAL_PERIOD_DISTRIBUTION_MAX,                          { TYPENAME::DOUBLE,     "Orbital_Period_Max",                     "days",      14, 6 }},
    { PROGRAM_OPTION::ORBITAL_PERIOD_DISTRIBUTION_MIN,                          { TYPENAME::DOUBLE,     "Orbital_Period_Min",                     "days",      14, 6 }},

    { PROGRAM_OPTION::OVERALL_WIND_MASS_LOSS_MULTIPLIER,                        { TYPENAME::DOUBLE,     "Overall_WindMassLoss_Multipl",           "-",         14, 6 }},

    { PROGRAM_OPTION::PISN_LOWER_LIMIT,                                         { TYPENAME::DOUBLE,     "PISN_Lower_Limit",                       "Msol",      14, 6 }},
    { PROGRAM_OPTION::PISN_UPPER_LIMIT,                                         { TYPENAME::DOUBLE,     "PISN_Upper_Limit",                       "Msol",      14, 6 }},

    { PROGRAM_OPTION::PPI_LOWER_LIMIT,                                          { TYPENAME::DOUBLE,     "PPI_Lower_Limit",                        "Msol",      14, 6 }},
    { PROGRAM_OPTION::PPI_PRESCRIPTION,                                         { TYPENAME::INT,        "PPI_Prscrptn",                           "-",          4, 1 }},
    { PROGRAM_OPTION::PPI_UPPER_LIMIT,                                          { TYPENAME::DOUBLE,     "PPI_Upper_Limit",                        "Msol",      14, 6 }},

    { PROGRAM_OPTION::PULSAR_MAGNETIC_FIELD_DISTRIBUTION,                       { TYPENAME::INT,        "Pulsar_Mag_Field_Dstrbtn",               "-",          4, 1 }},
    { PROGRAM_OPTION::PULSAR_MAGNETIC_FIELD_DISTRIBUTION_MAX,                   { TYPENAME::DOUBLE,     "Pulsar_Mag_Field_Dstrbtn_Max",           "AU",        14, 6 }},
    { PROGRAM_OPTION::PULSAR_MAGNETIC_FIELD_DISTRIBUTION_MIN,                   { TYPENAME::DOUBLE,     "Pulsar_Mag_Field_Dstrbtn_Min",           "AU",        14, 6 }},

    { PROGRAM_OPTION::PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION,                    { TYPENAME::INT,        "Pulsar_Spin_Period_Dstrbtn",             "-",          4, 1 }},
    { PROGRAM_OPTION::PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION_MAX,                { TYPENAME::DOUBLE,     "Pulsar_Spin_Period_Dstrbtn_Max",         "AU",        14, 6 }},
    { PROGRAM_OPTION::PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION_MIN,                { TYPENAME::DOUBLE,     "Pulsar_Spin_Period_Dstrbtn_Min",         "AU",        14, 6 }},

    { PROGRAM_OPTION::PULSAR_MAGNETIC_FIELD_DECAY_MASS_SCALE,                   { TYPENAME::DOUBLE,     "Pulsar_Mag_Field_Decay_mScale",          "Msol",      14, 6 }},
    { PROGRAM_OPTION::PULSAR_MAGNETIC_FIELD_DECAY_TIME_SCALE,                   { TYPENAME::DOUBLE,     "Pulsar_Mag_Field_Decay_tScale",          "Myr",       14, 6 }},

    { PROGRAM_OPTION::PULSAR_MINIMUM_MAGNETIC_FIELD,                            { TYPENAME::DOUBLE,     "Pulsar_Minimum_Mag_Field",               "Gauss",     14, 6 }},

    { PROGRAM_OPTION::QCRIT_PRESCRIPTION,                                       { TYPENAME::INT,        "qCrit_Prescription",                     "-",          4, 1 }},

    { PROGRAM_OPTION::RANDOM_SEED,                                              { TYPENAME::ULONGINT,   "SEED(OPTION)",                           "-",         12, 1 }},
    { PROGRAM_OPTION::RANDOM_SEED_CMDLINE,                                      { TYPENAME::ULONGINT,   "SEED(CMDLINE)",                          "-",         12, 1 }},

    { PROGRAM_OPTION::REMNANT_MASS_PRESCRIPTION,                                { TYPENAME::INT,        "Remnant_Mass_Prscrptn",                  "-",          4, 1 }},

    { PROGRAM_OPTION::ROTATIONAL_VELOCITY_DISTRIBUTION,                         { TYPENAME::INT,        "Rotational_Velocity_Dstrbtn",            "-",          4, 1 }},
    { PROGRAM_OPTION::ROTATIONAL_FREQUENCY,                                     { TYPENAME::DOUBLE,     "Rotational_Frequency",                   "Hz",        14, 6 }},
    { PROGRAM_OPTION::ROTATIONAL_FREQUENCY_1,                                   { TYPENAME::DOUBLE,     "Rotational_Frequency(1)",                "Hz",        14, 6 }},
    { PROGRAM_OPTION::ROTATIONAL_FREQUENCY_2,                                   { TYPENAME::DOUBLE,     "Rotational_Frequency(2)",                "Hz",        14, 6 }},
   
    { PROGRAM_OPTION::SEMI_MAJOR_AXIS,                                          { TYPENAME::DOUBLE,     "Semi-Major_Axis",                        "AU",        14, 6 }},
    { PROGRAM_OPTION::SEMI_MAJOR_AXIS_DISTRIBUTION,                             { TYPENAME::INT,        "Semi-Major_Axis_Dstrbtn",                "-",          4, 1 }},
    { PROGRAM_OPTION::SEMI_MAJOR_AXIS_DISTRIBUTION_MAX,                         { TYPENAME::DOUBLE,     "Semi-Major_Axis_Dstrbtn_Max",            "AU",        14, 6 }},
    { PROGRAM_OPTION::SEMI_MAJOR_AXIS_DISTRIBUTION_MIN,                         { TYPENAME::DOUBLE,     "Semi-Major_Axis_Dstrbtn_Min",            "AU",        14, 6 }},
    { PROGRAM_OPTION::SEMI_MAJOR_AXIS_DISTRIBUTION_POWER,                       { TYPENAME::DOUBLE,     "Semi-Major_Axis_Dstrbtn_Power",          "-",         14, 6 }},

    { PROGRAM_OPTION::STELLAR_ZETA_PRESCRIPTION,                                { TYPENAME::INT,        "Stellar_Zeta_Prscrptn",                  "-",          4, 1 }},

    { PROGRAM_OPTION::WR_FACTOR,                                                { TYPENAME::DOUBLE,     "WR_Factor",                              "-",         14, 6 }},

    { PROGRAM_OPTION::ZETA_ADIABATIC_ARBITRARY,                                 { TYPENAME::DOUBLE,     "Zeta_Adiabatic_Arbitrary",               "-",         14, 6 }},
    { PROGRAM_OPTION::ZETA_MS,                                                  { TYPENAME::DOUBLE,     "Zeta_Main_Sequence",                     "-",         14, 6 }},
    { PROGRAM_OPTION::ZETA_RADIATIVE_ENVELOPE_GIANT,                            { TYPENAME::DOUBLE,     "Zeta_Radiative_Envelope_Giant",          "-",         14, 6 }}
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

// BSE_BE_BINARY_REC
//
// Default record definition for the BeBinaries logfile
//
const ANY_PROPERTY_VECTOR BSE_BE_BINARIES_REC = {
    BINARY_PROPERTY::BE_BINARY_CURRENT_ID,
    BINARY_PROPERTY::RANDOM_SEED,
    BINARY_PROPERTY::BE_BINARY_CURRENT_DT,
    BINARY_PROPERTY::BE_BINARY_CURRENT_TOTAL_TIME,
    BINARY_PROPERTY::BE_BINARY_CURRENT_NS_MASS,
    BINARY_PROPERTY::BE_BINARY_CURRENT_COMPANION_MASS,
    BINARY_PROPERTY::BE_BINARY_CURRENT_COMPANION_LUMINOSITY,
    BINARY_PROPERTY::BE_BINARY_CURRENT_COMPANION_TEFF,
    BINARY_PROPERTY::BE_BINARY_CURRENT_COMPANION_RADIUS,
    BINARY_PROPERTY::BE_BINARY_CURRENT_SEMI_MAJOR_AXIS,
    BINARY_PROPERTY::BE_BINARY_CURRENT_ECCENTRICITY
};


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
    STAR_1_PROPERTY::NUCLEAR_TIMESCALE_PRE_COMMON_ENVELOPE,
    STAR_2_PROPERTY::IS_RLOF,
    STAR_2_PROPERTY::LUMINOSITY_PRE_COMMON_ENVELOPE,
    STAR_2_PROPERTY::TEMPERATURE_PRE_COMMON_ENVELOPE,
    STAR_2_PROPERTY::DYNAMICAL_TIMESCALE_PRE_COMMON_ENVELOPE,
    STAR_2_PROPERTY::THERMAL_TIMESCALE_PRE_COMMON_ENVELOPE,
    STAR_2_PROPERTY::NUCLEAR_TIMESCALE_PRE_COMMON_ENVELOPE,
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
    STAR_1_PROPERTY::NUCLEAR_TIMESCALE,
    STAR_2_PROPERTY::NUCLEAR_TIMESCALE,
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
    STAR_2_PROPERTY::RADIAL_EXPANSION_TIMESCALE
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
    BINARY_PROPERTY::DIMENSIONLESS_KICK_MAGNITUDE, 
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
    PROGRAM_OPTION::KICK_MAGNITUDE_DISTRIBUTION_SIGMA_CCSN_NS,
    PROGRAM_OPTION::KICK_MAGNITUDE_DISTRIBUTION_SIGMA_CCSN_BH,
    PROGRAM_OPTION::KICK_MAGNITUDE_DISTRIBUTION_SIGMA_FOR_ECSN,
    PROGRAM_OPTION::KICK_MAGNITUDE_DISTRIBUTION_SIGMA_FOR_USSN,
    PROGRAM_OPTION::LBV_FACTOR,
    PROGRAM_OPTION::WR_FACTOR,
    PROGRAM_OPTION::NOTES
};


// enum class LOGFILE
// Symbolic names for logfiles
enum class LOGFILE: int {
    NONE,

    DEBUG_LOG,
    ERROR_LOG,

    BSE_BE_BINARIES,
    BSE_COMMON_ENVELOPES,
    BSE_DETAILED_OUTPUT,
    BSE_DOUBLE_COMPACT_OBJECTS,
    BSE_PULSAR_EVOLUTION,
    BSE_RLOF_PARAMETERS,
    BSE_SUPERNOVAE,
    BSE_SWITCH_LOG,
    BSE_SYSTEM_PARAMETERS,

    SSE_DETAILED_OUTPUT,
    SSE_SUPERNOVAE,
    SSE_SWITCH_LOG,
    SSE_SYSTEM_PARAMETERS
};


// enum class LOGFILE_TYPE
// Symbolic names for logfile types
enum class LOGFILE_TYPE: int { NONE, STELLAR, BINARY };


typedef std::tuple<std::string, ANY_PROPERTY_VECTOR, std::string, std::string, LOGFILE_TYPE> LOGFILE_DESCRIPTOR_T;

// descriptors for logfiles
// unordered_map - key is integer logfile (from enum class LOGFILE above)
// fields are: {default filename, record descriptor, short file name, short record name, type}
// (the short names are for logfile definitions file parsing)
const std::map<LOGFILE, LOGFILE_DESCRIPTOR_T> LOGFILE_DESCRIPTOR = {
    { LOGFILE::NONE,                       { "" ,                              {},                             "",                "",                    LOGFILE_TYPE::NONE}},

    { LOGFILE::DEBUG_LOG,                  { "Debug_Log",                      {},                             "",                "",                    LOGFILE_TYPE::NONE }},
    { LOGFILE::ERROR_LOG,                  { "Error_Log",                      {},                             "",                "",                    LOGFILE_TYPE::NONE }},

    { LOGFILE::BSE_BE_BINARIES,            { "BSE_BE_Binaries",                BSE_BE_BINARIES_REC,            "BSE_BE_BINARIES", "BSE_BE_BINARIES_REC", LOGFILE_TYPE::BINARY }},
    { LOGFILE::BSE_COMMON_ENVELOPES,       { "BSE_Common_Envelopes",           BSE_COMMON_ENVELOPES_REC,       "BSE_CEE",         "BSE_CEE_REC",         LOGFILE_TYPE::BINARY }},
    { LOGFILE::BSE_DETAILED_OUTPUT,        { "BSE_Detailed_Output",            BSE_DETAILED_OUTPUT_REC,        "BSE_DETAILED",    "BSE_DETAILED_REC",    LOGFILE_TYPE::BINARY }},
    { LOGFILE::BSE_DOUBLE_COMPACT_OBJECTS, { "BSE_Double_Compact_Objects",     BSE_DOUBLE_COMPACT_OBJECTS_REC, "BSE_DCO",         "BSE_DCO_REC",         LOGFILE_TYPE::BINARY }},
    { LOGFILE::BSE_PULSAR_EVOLUTION,       { "BSE_Pulsar_Evolution",           BSE_PULSAR_EVOLUTION_REC,       "BSE_PULSARS",     "BSE_PULSARS_REC",     LOGFILE_TYPE::BINARY }},
    { LOGFILE::BSE_RLOF_PARAMETERS,        { "BSE_RLOF",                       BSE_RLOF_PARAMETERS_REC,        "BSE_RLOF",        "BSE_RLOF_REC",        LOGFILE_TYPE::BINARY }},
    { LOGFILE::BSE_SUPERNOVAE,             { "BSE_Supernovae",                 BSE_SUPERNOVAE_REC,             "BSE_SNE",         "BSE_SNE_REC",         LOGFILE_TYPE::BINARY }},
    { LOGFILE::BSE_SWITCH_LOG,             { "BSE_Switch_Log",                 BSE_SWITCH_LOG_REC,             "BSE_SWITCH_LOG",  "BSE_SWITCH_REC",      LOGFILE_TYPE::BINARY }},
    { LOGFILE::BSE_SYSTEM_PARAMETERS,      { "BSE_System_Parameters",          BSE_SYSTEM_PARAMETERS_REC,      "BSE_SYSPARMS",    "BSE_SYSPARMS_REC",    LOGFILE_TYPE::BINARY }},

    { LOGFILE::SSE_DETAILED_OUTPUT,        { "SSE_Detailed_Output",            SSE_DETAILED_OUTPUT_REC,        "SSE_DETAILED",    "SSE_DETAILED_REC",    LOGFILE_TYPE::STELLAR }},
    { LOGFILE::SSE_SUPERNOVAE,             { "SSE_Supernovae",                 SSE_SUPERNOVAE_REC,             "SSE_SNE",         "SSE_SNE_REC",         LOGFILE_TYPE::STELLAR }},
    { LOGFILE::SSE_SWITCH_LOG,             { "SSE_Switch_Log",                 SSE_SWITCH_LOG_REC,             "SSE_SWITCH_LOG",  "SSE_SWITCH_REC",      LOGFILE_TYPE::STELLAR }},
    { LOGFILE::SSE_SYSTEM_PARAMETERS,      { "SSE_System_Parameters",          SSE_SYSTEM_PARAMETERS_REC,      "SSE_SYSPARMS",    "SSE_SYSPARMS_REC",    LOGFILE_TYPE::STELLAR }}
};


// double vector CHE_Coefficients
// coefficients for the calculation of initial angular frequency
// Mandel from Butler 2018
const DBL_VECTOR CHE_Coefficients = { 5.7914E-04, -1.9196E-06, -4.0602E-07, 1.0150E-08, -9.1792E-11, 2.9051E-13 };


// enum class LR_TCoeff
// Symbolic names for term coefficients for Luminosity & Radius coefficients from Tout et al. 1996
enum class LR_TCoeff: int { a, b, c, d, e };
#define A LR_TCoeff::a
#define B LR_TCoeff::b
#define C LR_TCoeff::c
#define D LR_TCoeff::d
#define E LR_TCoeff::e

// enum class L_Coeff
// Symbolic names for luminosity coefficients (from Tout et al. 1996)
enum class L_Coeff: int { ALPHA, BETA, GAMMA, DELTA, EPSILON, ZETA, ETA };
#define ALPHA   L_Coeff::ALPHA
#define BETA    L_Coeff::BETA
#define GAMMA   L_Coeff::GAMMA
#define DELTA   L_Coeff::DELTA
#define EPSILON L_Coeff::EPSILON
#define ZETA    L_Coeff::ZETA
#define ETA     L_Coeff::ETA

// L (Luminosity) coefficients
// Table 1 in Tout et al 1996
// Key to map is L_Coeff.  Map element is unordered_map of term coefficient values.
const std::map<int, COMPASUnorderedMap<LR_TCoeff, double>> L_COEFF = {
    {static_cast<int>(ALPHA),   {{A, 0.39704170}, {B,  -0.32913574}, {C,  0.34776688}, {D,  0.37470851}, {E, 0.09011915}}},
    {static_cast<int>(BETA),    {{A, 8.52762600}, {B, -24.41225973}, {C, 56.43597107}, {D, 37.06152575}, {E, 5.45624060}}},
    {static_cast<int>(GAMMA),   {{A, 0.00025546}, {B,  -0.00123461}, {C, -0.00023246}, {D,  0.00045519}, {E, 0.00016176}}},
    {static_cast<int>(DELTA),   {{A, 5.43288900}, {B,  -8.62157806}, {C, 13.44202049}, {D, 14.51584135}, {E, 3.39793084}}},
    {static_cast<int>(EPSILON), {{A, 5.56357900}, {B, -10.32345224}, {C, 19.44322980}, {D, 18.97361347}, {E, 4.16903097}}},
    {static_cast<int>(ZETA),    {{A, 0.78866060}, {B,  -2.90870942}, {C,  6.54713531}, {D,  4.05606657}, {E, 0.53287322}}},
    {static_cast<int>(ETA),     {{A, 0.00586685}, {B,  -0.01704237}, {C,  0.03872348}, {D,  0.02570041}, {E, 0.00383376}}}
};

#undef ALPHA
#undef BETA
#undef GAMMA
#undef DELTA
#undef EPSILON
#undef ZETA
#undef ETA


// enum class R_Coeff
// Symbolic names for radius coefficients from Tout et al. 1996
enum class R_Coeff: int { THETA, IOTA, KAPPA, LAMBDA, MU, NU, XI, OMICRON, PI };
#define THETA   R_Coeff::THETA
#define IOTA    R_Coeff::IOTA
#define KAPPA   R_Coeff::KAPPA
#define LAMBDA  R_Coeff::LAMBDA
#define MU      R_Coeff::MU
#define NU      R_Coeff::NU
#define XI      R_Coeff::XI
#define OMICRON R_Coeff::OMICRON
#define Pi      R_Coeff::PI

// R (Radius) coefficients
// Table 2 in Tout et al. 1996
// Key to map is L_Coeff.  Map element is unordered_map of term coefficient values.
const std::map<int, COMPASUnorderedMap<LR_TCoeff, double>> R_COEFF = {
    {static_cast<int>(THETA),   {{A,  1.71535900}, {B,  0.62246212}, {C,  -0.92557761}, {D,  -1.16996966}, {E, -0.30631491}}},
    {static_cast<int>(IOTA),    {{A,  6.59778800}, {B, -0.42450044}, {C, -12.13339427}, {D, -10.73509484}, {E, -2.51487077}}},
    {static_cast<int>(KAPPA),   {{A, 10.08855000}, {B, -7.11727086}, {C, -31.67119479}, {D, -24.24848322}, {E, -5.33608972}}},
    {static_cast<int>(LAMBDA),  {{A,  1.01249500}, {B,  0.32699690}, {C,  -0.00923418}, {D,  -0.03876858}, {E, -0.00412750}}},
    {static_cast<int>(MU),      {{A,  0.07490166}, {B,  0.02410413}, {C,   0.07233664}, {D,   0.03040467}, {E,  0.00197741}}},
    {static_cast<int>(NU),      {{A,  0.01077422}, {B,  0.00000000}, {C,   0.00000000}, {D,   0.00000000}, {E,  0.00000000}}},
    {static_cast<int>(XI),      {{A,  3.08223400}, {B,  0.94472050}, {C,  -2.15200882}, {D,  -2.49219496}, {E, -0.63848738}}},
    {static_cast<int>(OMICRON), {{A, 17.84778000}, {B, -7.45345690}, {C, -48.96066856}, {D, -40.05386135}, {E, -9.09331816}}},
    {static_cast<int>(Pi),      {{A,  0.00022582}, {B, -0.00186899}, {C,   0.00388783}, {D,   0.00142402}, {E, -0.00007671}}}
};

#undef THETA
#undef IOTA
#undef KAPPA
#undef LAMBDA
#undef MU
#undef NU
#undef XI
#undef OMICRON
#undef Pi

#undef A
#undef B
#undef C
#undef D
#undef E


// enum class AB_TCoeff
// Symbolic names for term coefficients for A & B coefficients from Hurley et al. 2000
enum class AB_TCoeff: int { ALPHA, BETA, GAMMA, ETA, MU };
#define ALPHA AB_TCoeff::ALPHA
#define BETA  AB_TCoeff::BETA
#define GAMMA AB_TCoeff::GAMMA
#define ETA   AB_TCoeff::ETA
#define MU    AB_TCoeff::MU

// A coefficients
// Table in Appendix A of Hurley et al. 2000
// Key to map is n (A(n)).  Map element is unordered_map of term coefficient values.
const std::map<int, COMPASUnorderedMap<AB_TCoeff, double>> A_COEFF = {
    { 1, {{ALPHA,  1.593890E3 }, {BETA,  2.053038E3 }, {GAMMA,  1.231226E3 }, {ETA,  2.327785E2 }, {MU,  0.000000E0 }}},
    { 2, {{ALPHA,  2.706708E3 }, {BETA,  1.483131E3 }, {GAMMA,  5.772723E2 }, {ETA,  7.411230E1 }, {MU,  0.000000E0 }}},
    { 3, {{ALPHA,  1.466143E2 }, {BETA, -1.048442E2 }, {GAMMA, -6.795374E1 }, {ETA, -1.391127E1 }, {MU,  0.000000E0 }}},
    { 4, {{ALPHA,  4.141960E-2}, {BETA,  4.564888E-2}, {GAMMA,  2.958542E-2}, {ETA,  5.571483E-3}, {MU,  0.000000E0 }}},
    { 5, {{ALPHA,  3.426349E-1}, {BETA,  0.000000E0 }, {GAMMA,  0.000000E0 }, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    { 6, {{ALPHA,  1.949814E1 }, {BETA,  1.758178E0 }, {GAMMA, -6.008212E0 }, {ETA, -4.470533E0 }, {MU,  0.000000E0 }}},
    { 7, {{ALPHA,  4.903830E0 }, {BETA,  0.000000E0 }, {GAMMA,  0.000000E0 }, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    { 8, {{ALPHA,  5.212154E-2}, {BETA,  3.166411E-2}, {GAMMA, -2.750074E-3}, {ETA, -2.271549E-3}, {MU,  0.000000E0 }}},
    { 9, {{ALPHA,  1.312179E0 }, {BETA, -3.294936E-1}, {GAMMA,  9.231860E-2}, {ETA,  2.610989E-2}, {MU,  0.000000E0 }}},
    {10, {{ALPHA,  8.073972E-1}, {BETA,  0.000000E0 }, {GAMMA,  0.000000E0 }, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},

    {11, {{ALPHA,  1.031538E0 }, {BETA, -2.434480E-1}, {GAMMA,  7.732821E0 }, {ETA,  6.460705E0 }, {MU,  1.374484E0 }}},
    {12, {{ALPHA,  1.043715E0 }, {BETA, -1.577474E0 }, {GAMMA, -5.168234E0 }, {ETA, -5.596506E0 }, {MU, -1.299394E0 }}},
    {13, {{ALPHA,  7.859573E2 }, {BETA, -8.542048E0 }, {GAMMA, -2.642511E1 }, {ETA, -9.585707E0 }, {MU,  0.000000E0 }}},
    {14, {{ALPHA,  3.858911E3 }, {BETA,  2.459681E3 }, {GAMMA, -7.630093E1 }, {ETA, -3.486057E2 }, {MU, -4.861703E1 }}},
    {15, {{ALPHA,  2.888720E2 }, {BETA,  2.952979E2 }, {GAMMA,  1.850341E2 }, {ETA,  3.797254E1 }, {MU,  0.000000E0 }}},
    {16, {{ALPHA,  7.196580E0 }, {BETA,  5.613746E-1}, {GAMMA,  3.805871E-1}, {ETA,  8.398728E-2}, {MU,  0.000000E0 }}},
    {17, {{ALPHA,  0.000000E0 }, {BETA,  0.000000E0 }, {GAMMA,  0.000000E0 }, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    {18, {{ALPHA,  2.187715E-1}, {BETA, -2.154437E0 }, {GAMMA, -3.768678E0 }, {ETA, -1.975518E0 }, {MU, -3.021475E-1}}},
    {19, {{ALPHA,  1.466440E0 }, {BETA,  1.839725E0 }, {GAMMA,  6.442199E0 }, {ETA,  4.023635E0 }, {MU,  6.957529E-1}}},
    {20, {{ALPHA,  2.652091E1 }, {BETA,  8.178458E1 }, {GAMMA,  1.156058E2 }, {ETA,  7.633811E1 }, {MU,  1.950698E1 }}},

    {21, {{ALPHA,  1.472103E0 }, {BETA, -2.947609E0 }, {GAMMA, -3.312828E0 }, {ETA, -9.945065E-1}, {MU,  0.000000E0 }}},
    {22, {{ALPHA,  3.071048E0 }, {BETA, -5.679941E0 }, {GAMMA, -9.745523E0 }, {ETA, -3.594543E0 }, {MU,  0.000000E0 }}},
    {23, {{ALPHA,  2.617890E0 }, {BETA,  1.019135E0 }, {GAMMA, -3.292551E-2}, {ETA, -7.445123E-2}, {MU,  0.000000E0 }}},
    {24, {{ALPHA,  1.075567E-2}, {BETA,  1.773287E-2}, {GAMMA,  9.610479E-3}, {ETA,  1.732469E-3}, {MU,  0.000000E0 }}},
    {25, {{ALPHA,  1.476246E0 }, {BETA,  1.899331E0 }, {GAMMA,  1.195010E0 }, {ETA,  3.035051E-1}, {MU,  0.000000E0 }}},
    {26, {{ALPHA,  5.502535E0 }, {BETA, -6.601663E-2}, {GAMMA,  9.968707E-2}, {ETA,  3.599801E-2}, {MU,  0.000000E0 }}},
    {27, {{ALPHA,  9.511033E1 }, {BETA,  6.819618E1 }, {GAMMA, -1.045625E1 }, {ETA, -1.474939E1 }, {MU,  0.000000E0 }}},
    {28, {{ALPHA,  3.113458E1 }, {BETA,  1.012033E1 }, {GAMMA, -4.650511E0 }, {ETA, -2.463185E0 }, {MU,  0.000000E0 }}},
    {29, {{ALPHA,  1.413057E0 }, {BETA,  4.578814E-1}, {GAMMA, -6.850581E-2}, {ETA, -5.588658E-2}, {MU,  0.000000E0 }}},
    {30, {{ALPHA,  3.910862E1 }, {BETA,  5.196646E1 }, {GAMMA,  2.264970E1 }, {ETA,  2.873680E0 }, {MU,  0.000000E0 }}},

    {31, {{ALPHA,  4.597479E0 }, {BETA, -2.855179E-1}, {GAMMA,  2.709724E-1}, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    {32, {{ALPHA,  6.682518E0 }, {BETA,  2.827718E-1}, {GAMMA, -7.294429E-2}, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    {33, {{ALPHA,  0.000000E0 }, {BETA,  0.000000E0 }, {GAMMA,  0.000000E0 }, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    {34, {{ALPHA,  1.910302E-1}, {BETA,  1.158624E-1}, {GAMMA,  3.348990E-2}, {ETA,  2.599706E-3}, {MU,  0.000000E0 }}},
    {35, {{ALPHA,  3.931056E-1}, {BETA,  7.277637E-2}, {GAMMA, -1.366593E-1}, {ETA, -4.508946E-2}, {MU,  0.000000E0 }}},
    {36, {{ALPHA,  3.267776E-1}, {BETA,  1.204424E-1}, {GAMMA,  9.988332E-2}, {ETA,  2.455361E-2}, {MU,  0.000000E0 }}},
    {37, {{ALPHA,  5.990212E-1}, {BETA,  5.570264E-2}, {GAMMA,  6.207626E-2}, {ETA,  1.777283E-2}, {MU,  0.000000E0 }}},
    {38, {{ALPHA,  7.330122E-1}, {BETA,  5.192827E-1}, {GAMMA,  2.316416E-1}, {ETA,  8.346941E-3}, {MU,  0.000000E0 }}},
    {39, {{ALPHA,  1.172768E0 }, {BETA, -1.209262E-1}, {GAMMA, -1.193023E-1}, {ETA, -2.859837E-2}, {MU,  0.000000E0 }}},
    {40, {{ALPHA,  3.982622E-1}, {BETA, -2.296279E-1}, {GAMMA, -2.262539E-1}, {ETA, -5.219837E-2}, {MU,  0.000000E0 }}},

    {41, {{ALPHA,  3.571038E0 }, {BETA, -2.223635E-2}, {GAMMA, -2.611794E-2}, {ETA, -6.359648E-3}, {MU,  0.000000E0 }}},
    {42, {{ALPHA,  1.984800E0 }, {BETA,  1.138600E0 }, {GAMMA,  3.564000E-1}, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    {43, {{ALPHA,  6.300000E-2}, {BETA,  4.810000E-2}, {GAMMA,  9.840000E-3}, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    {44, {{ALPHA,  1.200000E0 }, {BETA,  2.450000E0 }, {GAMMA,  0.000000E0 }, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    {45, {{ALPHA,  2.321400E-1}, {BETA,  1.828075E-3}, {GAMMA, -2.232007E-2}, {ETA, -3.378734E-3}, {MU,  0.000000E0 }}},
    {46, {{ALPHA,  1.163659E-2}, {BETA,  3.427682E-3}, {GAMMA,  1.421393E-3}, {ETA, -3.710666E-3}, {MU,  0.000000E0 }}},
    {47, {{ALPHA,  1.048020E-2}, {BETA, -1.231921E-2}, {GAMMA, -1.686860E-2}, {ETA, -4.234354E-3}, {MU,  0.000000E0 }}},
    {48, {{ALPHA,  1.555590E0 }, {BETA, -3.223927E-1}, {GAMMA, -5.197429E-1}, {ETA, -1.066441E-1}, {MU,  0.000000E0 }}},
    {49, {{ALPHA,  9.770000E-2}, {BETA, -2.310000E-1}, {GAMMA, -7.530000E-2}, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    {50, {{ALPHA,  2.400000E-1}, {BETA,  1.800000E-1}, {GAMMA,  5.950000E-1}, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},

    {51, {{ALPHA,  3.300000E-1}, {BETA,  1.320000E-1}, {GAMMA,  2.180000E-1}, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    {52, {{ALPHA,  1.106400E0 }, {BETA,  4.150000E-1}, {GAMMA,  1.800000E-1}, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    {53, {{ALPHA,  1.190000E0 }, {BETA,  3.770000E-1}, {GAMMA,  1.760000E-1}, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    {54, {{ALPHA,  3.855707E-1}, {BETA, -6.104166E-1}, {GAMMA,  5.676742E0 }, {ETA,  1.060894E1 }, {MU,  5.284014E0 }}},
    {55, {{ALPHA,  3.579064E-1}, {BETA, -6.442936E-1}, {GAMMA,  5.494644E0 }, {ETA,  1.054952E1 }, {MU,  5.280991E0 }}},
    {56, {{ALPHA,  9.587587E-1}, {BETA,  8.777464E-1}, {GAMMA,  2.017321E-1}, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    {57, {{ALPHA,  1.513500E0 }, {BETA,  3.769000E-1}, {GAMMA,  0.000000E0 }, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    {58, {{ALPHA,  4.907546E-1}, {BETA, -1.683928E-1}, {GAMMA, -3.108742E-1}, {ETA, -7.202918E-2}, {MU,  0.000000E0 }}},
    {59, {{ALPHA,  4.537070E0 }, {BETA, -4.465455E0 }, {GAMMA, -1.612690E0 }, {ETA, -1.623246E0 }, {MU,  0.000000E0 }}},
    {60, {{ALPHA,  1.796220E0 }, {BETA,  2.814020E-1}, {GAMMA,  1.423325E0 }, {ETA,  3.421036E-1}, {MU,  0.000000E0 }}},

    {61, {{ALPHA,  2.256216E0 }, {BETA,  3.773400E-1}, {GAMMA,  1.537867E0 }, {ETA,  4.396373E-1}, {MU,  0.000000E0 }}},
    {62, {{ALPHA,  8.430000E-2}, {BETA, -4.750000E-2}, {GAMMA, -3.520000E-2}, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    {63, {{ALPHA,  7.360000E-2}, {BETA,  7.490000E-2}, {GAMMA,  4.426000E-2}, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    {64, {{ALPHA,  1.360000E-1}, {BETA,  3.520000E-2}, {GAMMA,  0.000000E0 }, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    {65, {{ALPHA,  1.564231E-3}, {BETA,  1.653042E-3}, {GAMMA, -4.439786E-3}, {ETA, -4.951011E-3}, {MU, -1.216530E-3}}},
    {66, {{ALPHA,  1.477000E0 }, {BETA,  2.960000E-1}, {GAMMA,  0.000000E0 }, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    {67, {{ALPHA,  5.210157E0 }, {BETA, -4.143695E0 }, {GAMMA, -2.120870E0 }, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    {68, {{ALPHA,  1.116000E0 }, {BETA,  1.660000E-1}, {GAMMA,  0.000000E0 }, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    {69, {{ALPHA,  1.071489E0 }, {BETA, -1.164852E-1}, {GAMMA, -8.623831E-2}, {ETA, -1.582349E-2}, {MU,  0.000000E0 }}},
    {70, {{ALPHA,  7.108492E-1}, {BETA,  7.935927E-1}, {GAMMA,  3.926983E-1}, {ETA,  3.622146E-2}, {MU,  0.000000E0 }}},

    {71, {{ALPHA,  3.478514E0 }, {BETA, -2.585474E-2}, {GAMMA, -1.512955E-2}, {ETA, -2.833691E-3}, {MU,  0.000000E0 }}},
    {72, {{ALPHA,  9.132108E-1}, {BETA, -1.653695E-1}, {GAMMA,  0.000000E0 }, {ETA,  3.636784E-2}, {MU,  0.000000E0 }}},
    {73, {{ALPHA,  3.969331E-3}, {BETA,  4.539076E-3}, {GAMMA,  1.720906E-3}, {ETA,  1.897857E-4}, {MU,  0.000000E0 }}},
    {74, {{ALPHA,  1.600000E0 }, {BETA,  7.640000E-1}, {GAMMA,  3.322000E-1}, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    {75, {{ALPHA,  8.109000E-1}, {BETA, -6.282000E-1}, {GAMMA,  0.000000E0 }, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    {76, {{ALPHA,  1.192334E-2}, {BETA,  1.083057E-2}, {GAMMA,  1.230969E0 }, {ETA,  1.551656E0 }, {MU,  0.000000E0 }}},
    {77, {{ALPHA, -1.668868E-1}, {BETA,  5.818123E-1}, {GAMMA, -1.105027E1 }, {ETA, -1.668070E1 }, {MU,  0.000000E0 }}},
    {78, {{ALPHA,  7.615495E-1}, {BETA,  1.068243E-1}, {GAMMA, -2.011333E-1}, {ETA, -9.371415E-2}, {MU,  0.000000E0 }}},
    {79, {{ALPHA,  9.409838E0 }, {BETA,  1.522928E0 }, {GAMMA,  0.000000E0 }, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    {80, {{ALPHA, -2.711000E-1}, {BETA, -5.756000E-1}, {GAMMA, -8.380000E-2}, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},

    {81, {{ALPHA,  2.493000E0 }, {BETA,  1.147500E0 }, {GAMMA,  0.000000E0 }, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}}
};

// Critial mass ratios and zetas, for a grid of masses and radii 
// This is a subset of an unpublished table from Hongwei Ge, which replaces the previous table 3 from Ge et al. 2020. 
// This table contains just the mass, logR, qCritSTD_full, qCritIC_full, qCritSTD_nonc, qCritIC_nonc, where
// mass is in Msol, logR uses Rsol, and qCrit is the critical mass ratio for dynamically unstable mass transfer.
// The STD and IC suffixes distinguish their adiabatic and isentropic models, respectively (see Ge et al. 2020).
// The _full and _nonc suffixes distinguish their models using fully conservative MT and fully non-conservative MT, 
// respectively, where the fully non-conservative model assumes isotropic re-emission AM loss. This fully non-conservative
// mass loss model is new, and implemented in Li et al. 2023, though the table was not made public. However, it is useful
// to have both mass loss extremes. 
//
// In all cases, q is mAccretor/mDonor, which is inverted from the Ge et al. datatable. First entry in the tuple is the vector of 
// unique mass values, second entry is the 5-tuple of vectors for logR, qCritSTD_full, qCritIC_full, qCritSTD_nonc, qCritIC_nonc.
// Note that the radius may contract several times. These points have been removed to facilitate the interpolation, so logR is monotonic.
// Previously, this table contained the commonly used Mass-Radius relation zeta's, however these are functionality identical to the critical
// mass ratios of the same derivation, and were not used in the code anyway, so they have been replaced by the non-conservative qcrits.
const std::tuple< std::vector<double>, std::vector< std::tuple<std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>>>> GE20_QCRIT = {
    {0.100, 0.130, 0.160, 0.200, 0.220, 0.250, 0.280, 0.320, 0.360, 0.400, 0.450, 0.500, 0.560, 0.630, 0.710, 0.800, 0.890, 1.000, 1.140, 1.300, 1.439, 1.600, 1.800, 2.000, 2.040, 2.500, 3.200, 4.000, 5.000, 6.300, 8.000, 10.000, 13.000, 16.000, 20.000, 25.000, 32.000, 40.000, 50.000, 63.000, 80.000, 100.000},
    {
      {{-1.1463, -1.1653, -1.1824, -1.1951, -1.1997}, {-0.0027, -0.0018, -0.0021, 0.0020, 0.0033}, {1.5165, 1.5324, 1.5344, 1.5360, 1.5374}, {-94.8767, -109.5290, -103.1992, -390.6250, 9090.9091}, {1.0060, 1.0132, 1.0143, 1.0151, 1.0159}},
      {{-1.2868, -1.3051, -1.3191, -1.3254, -1.3707, -1.3815}, {-0.0039, -0.0005, 0.0002, -0.0001, 1.0629, 4.7864}, {1.5314, 1.5341, 1.5365, 1.5396, 0.8280, 0.2087}, {-80.0000, -141.2429, -169.2047, -152.6718, -12.8584, 219.2982}, {1.0127, 1.0142, 1.0155, 1.0171, 0.6047, 0.3257}},
      {{-1.4149, -1.4295, -1.4372, -1.4382, -1.4391, -1.5056, -1.7789, -2.1738, -2.7918, -3.9062}, {0.0026, -0.0007, 0.0038, -0.0005, 0.0021, 0.6690, 1.6660, 2.7340, 4.8944, 10.3743}, {1.5333, 1.5373, 1.5411, 1.5434, 1.4267, 1.1867, 0.5800, 0.3543, 0.1962, 0.0961}, {-847.4576, -133.5113, 847.4576, -139.6648, -26.5534, 226.7574, 0.4589, 0.2982, 0.2641, 0.2498}, {1.0137, 1.0159, 1.0180, 1.0193, 0.9546, 0.8190, 0.4468, 0.2908, 0.2640, 0.2495}},
      {{-1.5692, -1.5794, -1.5860, -1.5875, -1.5885, -1.5997, -1.6349, -1.7672, -2.1612, -2.7807, -3.9064, -6.5300, -21.5332, 20.3500, 6.5902, 3.8210, 2.7573}, {-0.0006, 0.0033, 0.0001, 0.0017, 0.0040, 0.0003, 0.0039, 0.7641, 1.2152, 1.5052, 1.7401, 2.0190, 2.4122, 2.9280, 3.8741, 6.4279, 23.6956}, {1.5383, 1.5432, 1.5480, 1.5510, 1.5334, 1.4843, 1.3805, 1.1504, 0.7969, 0.6490, 0.5597, 0.4787, 0.3965, 0.3229, 0.2405, 0.1420, 0.0415}, {-137.3626, 6250.0000, -161.2903, -53.4188, -32.7869, -28.0584, -23.0044, -15.2602, -11.3934, -11.3934, 0.4602, 0.4602, 0.3359, 0.4602, 0.4602, -11.3934, 0.4872}, {1.0165, 1.0192, 1.0218, 1.0235, 1.0139, 0.9868, 0.9292, 0.7982, 0.5864, 0.4935, 0.4357, 0.3818, 0.3251, 0.2722, 0.2121, 0.2042, 0.0164}},
      {{-1.6428, -1.6536, -1.6609, -1.6624, -1.6687, -1.6865, -1.7349, -1.9019, -2.3723, -3.1522, -4.7015, -9.2524, -205.7613, 10.0361, 4.8716, 3.2141, 2.4500, 1.9207, 1.6102}, {0.0006, 0.0016, 0.0028, 0.0004, 0.0040, 0.0025, 0.0011, 0.7885, 1.1493, 1.3542, 1.4825, 1.6400, 1.8402, 2.1054, 2.4835, 3.0780, 4.0518, 6.8669, 32.8235}, {1.5414, 1.5463, 1.5510, 1.5489, 1.5223, 1.4695, 1.3637, 1.1424, 0.8436, 0.7233, 0.6594, 0.5921, 0.5229, 0.4518, 0.3778, 0.2999, 0.2237, 0.1285, 0.0296}, {-195.6947, -305.8104, -64.1026, -34.9895, 227.2727, -26.4901, -21.4915, -14.8478, 226.7574, 0.5501, 0.5103, 0.4689, 0.4257, 0.3797, 0.3297, 0.2739, 0.2151, 0.1339, 0.0216}, {1.0182, 1.0209, 1.0235, 1.0225, 1.0079, 0.9787, 0.9198, 0.7936, 0.6154, 0.5414, 0.5018, 0.4592, 0.4140, 0.3663, 0.3149, 0.2584, 0.2002, 0.1219, 0.0187}},
      {{-1.7535, -1.7657, -1.7736, -1.7770, -1.7910, -1.8191, -1.8883, -2.1083, -2.6823, -3.6839, -5.8855, -14.1323, 31.3381, 7.6005, 4.1976, 2.9839, 2.2801, 1.8316, 1.5547, 1.3527, 1.1906, 1.0574}, {0.0030, 0.0009, 0.0032, 0.0010, 0.0032, 0.0030, 0.0139, 0.8091, 1.0735, 1.2025, 1.2681, 1.3483, 1.4592, 1.5911, 1.7669, 1.9749, 2.2775, 2.7471, 3.4229, 4.5951, 7.3470, 49.1074}, {1.5458, 1.5505, 1.5510, 1.5391, 1.5055, 1.4492, 1.3434, 1.1385, 0.9045, 0.8168, 0.7744, 0.7250, 0.6648, 0.6036, 0.5367, 0.4738, 0.4042, 0.3287, 0.2582, 0.1878, 0.1134, 0.0191}, {-1612.9032, -214.5923, -36.6569, -31.9183, -28.1373, -9.4482, -9.4482, 0.8464, 0.6682, 0.6086, 0.5830, 0.5545, 0.5197, 0.4839, 0.4435, 0.4043, 0.3586, 0.3059, 0.2533, 0.1960, 0.1290, 0.0210}, {1.0207, 1.0232, 1.0238, 1.0172, 0.9987, 0.9676, 0.9085, 0.7915, 0.6530, 0.6000, 0.5752, 0.5457, 0.5087, 0.4702, 0.4273, 0.3858, 0.3383, 0.2845, 0.2320, 0.1763, 0.1130, 0.0173}},
      {{-1.8656, -1.8790, -1.8844, -1.8954, -1.9170, -1.9561, -2.0481, -2.3232, -3.0525, -4.4456, -8.2352, -59.5948, 12.4254, 5.3559, 3.4090, 2.5589, 2.0185, 1.6864, 1.4144, 1.2344, 1.0956, 0.9880, 0.8974, 0.8230}, {-0.0005, -0.0008, 0.0014, 0.0033, -0.0004, -0.0007, 0.6118, 0.8204, 1.0378, 1.1226, 1.1581, 1.2163, 1.2921, 1.3923, 1.5137, 1.6505, 1.8389, 2.0728, 2.4535, 2.9708, 3.8015, 5.2471, 8.6728, 53.9322}, {1.5497, 1.5514, 1.5480, 1.5260, 1.4888, 1.4310, 1.3265, 1.1385, 0.9271, 0.8761, 0.8505, 0.8065, 0.7537, 0.6925, 0.6291, 0.5695, 0.5031, 0.4388, 0.3622, 0.2918, 0.2215, 0.1551, 0.0891, 0.0128}, {-137.7410, -51.7331, -34.8432, -29.9491, -26.6880, -9.4518, -9.4518, 0.8367, 0.6877, 0.6451, 0.6297, 0.6055, 0.5765, 0.5427, 0.5070, 0.4723, 0.4323, 0.3916, 0.3405, 0.2900, 0.2353, 0.1784, 0.1144, 0.1556}, {1.0229, 1.0240, 1.0222, 1.0101, 0.9896, 0.9575, 0.8991, 0.7916, 0.6640, 0.6366, 0.6226, 0.5972, 0.5660, 0.5290, 0.4900, 0.4525, 0.4096, 0.3667, 0.3137, 0.2626, 0.2088, 0.1546, 0.0955, 0.0163}},
      {{-2.0177, -2.0294, -2.0374, -2.0574, -2.0894, -2.1450, -2.2696, -2.6272, -3.5556, -5.4780, -11.7716, 57.4713, 8.5653, 4.6492, 3.1699, 2.4477, 1.9605, 1.6418, 1.4082, 1.2296, 1.1046, 0.9887, 0.9055, 0.8273, 0.7645, 0.7104, 0.6654}, {-0.0004, 0.0025, 0.0023, -0.0016, 0.0005, -0.0003, 0.6713, 0.8269, 0.9930, 1.0588, 1.0742, 1.1079, 1.1636, 1.2330, 1.3148, 1.4042, 1.5200, 1.6483, 1.8153, 2.0492, 2.3298, 2.7862, 3.3805, 4.4578, 6.3077, 10.9183, 71.2088}, {1.5510, 1.5504, 1.5379, 1.5105, 1.4710, 1.4122, 1.3114, 1.1427, 0.9820, 0.9290, 0.9186, 0.8888, 0.8413, 0.7871, 0.7303, 0.6760, 0.6155, 0.5588, 0.4982, 0.4318, 0.3713, 0.3012, 0.2405, 0.1748, 0.1169, 0.0624, 0.0078}, {-9.4322, -35.2858, -31.7662, -29.4898, -26.1780, -9.4420, -9.4420, 0.8313, 0.7141, 0.6776, 0.6706, 0.6549, 0.6298, 0.6013, 0.5711, 0.5419, 0.5084, 0.4765, 0.4410, 0.3999, 0.3603, 0.3113, 0.2655, 0.2110, 0.1573, 0.0979, 0.0184}, {1.0238, 1.0236, 1.0167, 1.0016, 0.9798, 0.9471, 0.8906, 0.7941, 0.7001, 0.6690, 0.6639, 0.6477, 0.6205, 0.5890, 0.5553, 0.5227, 0.4857, 0.4504, 0.4118, 0.3679, 0.3266, 0.2765, 0.2309, 0.1787, 0.1285, 0.0765, 0.0129}},
      {{-2.1736, -2.1797, -2.1980, -2.2267, -2.2705, -2.3436, -2.5043, -2.9491, -4.1466, -6.9171, -22.7015, 18.5977, 6.5630, 4.0555, 2.9099, 2.2509, 1.8439, 1.5714, 1.3585, 1.1921, 1.0745, 0.9642, 0.8841, 0.8163, 0.7511, 0.7033, 0.6582, 0.6188, 0.5820}, {-0.0012, 0.0010, -0.0003, 0.0016, -0.0009, 0.0007, 0.6913, 0.8272, 0.9620, 1.0039, 1.0083, 1.0326, 1.0802, 1.1342, 1.1958, 1.2694, 1.3555, 1.4536, 1.5767, 1.7365, 1.9190, 2.1985, 2.5363, 3.0077, 3.8161, 4.9282, 7.0071, 12.1095, 92.7260}, {1.5518, 1.5457, 1.5252, 1.4955, 1.4548, 1.3969, 1.3009, 1.1505, 1.0150, 0.9804, 0.9790, 0.9558, 0.9091, 0.8596, 0.8075, 0.7517, 0.6945, 0.6382, 0.5782, 0.5142, 0.4555, 0.3862, 0.3249, 0.2643, 0.1982, 0.1459, 0.0943, 0.0488, 0.0041}, {-9.3809, -9.3809, -30.1114, -28.1057, -24.3605, -9.3721, -9.3721, 0.8312, 0.7335, 0.7088, 0.7072, 0.6945, 0.6703, 0.6450, 0.6186, 0.5902, 0.5603, 0.5302, 0.4973, 0.4607, 0.4256, 0.3820, 0.3408, 0.2973, 0.2453, 0.1996, 0.1490, 0.1404, 0.1345}, {1.0244, 1.0211, 1.0098, 0.9934, 0.9708, 0.9386, 0.8847, 0.7987, 0.7200, 0.7000, 0.7003, 0.6880, 0.6619, 0.6338, 0.6039, 0.5715, 0.5377, 0.5039, 0.4670, 0.4268, 0.3888, 0.3425, 0.2997, 0.2556, 0.2047, 0.1616, 0.1146, 0.0688, 0.0105}},
      {{-2.3446, -2.3706, -2.4053, -2.4518, -2.5182, -2.6256, -2.8525, -3.4593, -5.2740, -11.1235, 112.6126, 9.2464, 4.7454, 3.2443, 2.4567, 1.9693, 1.6579, 1.4110, 1.2366, 1.0969, 0.9964, 0.9008, 0.8304, 0.7622, 0.7113, 0.6618, 0.6251, 0.5530, 0.5249}, {0.0011, 0.0038, 0.0019, -0.0016, 0.0044, 0.6707, 0.7494, 0.8706, 0.9825, 0.9958, 0.9864, 1.0040, 1.0437, 1.0893, 1.1411, 1.2029, 1.2751, 1.3678, 1.4711, 1.6028, 1.7509, 1.9731, 2.2334, 2.6452, 3.1603, 4.0428, 5.2011, 14.4142, 160.6347}, {1.5331, 1.5094, 1.4803, 1.4448, 1.4002, 1.3400, 1.2464, 1.1092, 0.9985, 0.9903, 1.0016, 0.9824, 0.9396, 0.8935, 0.8444, 0.7915, 0.7369, 0.6755, 0.6171, 0.5543, 0.4964, 0.4276, 0.3662, 0.2962, 0.2369, 0.1738, 0.1255, 0.0308, 0.0015}, {-9.3467, -9.3467, -9.3467, 226.7574, -9.3440, -9.3440, -15.9974, 0.7967, 0.7210, 0.7140, 0.7207, 0.7116, 0.6905, 0.6679, 0.6441, 0.6183, 0.5908, 0.5593, 0.5287, 0.4946, 0.4617, 0.4208, 0.3820, 0.3347, 0.2912, 0.2403, 0.1973, 0.1352, 0.0131}, {1.0143, 1.0012, 0.9852, 0.9655, 0.9406, 0.9069, 0.8538, 0.7749, 0.7103, 0.7063, 0.7143, 0.7047, 0.6812, 0.6556, 0.6280, 0.5979, 0.5663, 0.5303, 0.4953, 0.4570, 0.4206, -0.0101, 0.3349, 0.2859, 0.2422, 0.1931, 0.1518, 0.0558, 0.0076}},
      {{-2.6088, -2.6543, -2.7130, -2.7904, -2.8993, -3.0739, -3.4294, -4.3501, -7.7682, -36.2845, 14.1163, 5.8197, 3.6310, 2.6553, 2.0957, 1.7316, 1.4686, 1.2897, 1.1435, 1.0237, 0.9248, 0.8517, 0.7806, 0.7272, 0.6747, 0.6299, 0.5918, 0.5594, 0.5292, 0.5033, 0.4792}, {0.5354, 0.5963, 0.6322, 0.6634, 0.7009, 0.7515, 0.8288, 0.9394, 1.0206, 0.9895, 0.9674, 0.9788, 1.0098, 1.0485, 1.0932, 1.1454, 1.2132, 1.2804, 1.3638, 1.4694, 1.6050, 1.7576, 1.9875, 2.2571, 2.6842, 3.3120, 4.2366, 5.7211, 8.7536, 17.3591, 201.3561}, {1.4789, 1.4484, 1.4128, 1.3712, 1.3209, 1.2551, 1.1606, 1.0400, 0.9654, 0.9985, 1.0211, 1.0056, 0.9684, 0.9246, 0.8776, 0.8276, 0.7700, 0.7188, 0.6631, 0.6025, 0.5376, 0.4782, 0.4079, 0.3457, 0.2757, 0.2093, 0.1501, 0.0973, 0.0506, 0.0124, 0.0012}, {1.5549, 1.0918, 1.0400, 0.9991, 0.9544, 0.9005, 0.8300, 0.7481, 0.6987, 0.7182, 0.7332, 0.7276, 0.7106, 0.6904, 0.6687, 0.6451, 0.6171, 0.5923, 0.5645, 0.5334, 0.4987, 0.4653, 0.4238, 0.3847, 0.3373, 0.2879, 0.2399, 0.1908, 0.1407, 0.0859, 0.1258}, {0.9846, 0.9677, 0.9479, 0.9246, 0.8962, 0.8589, 0.8046, 0.7344, 0.6908, 0.7118, 0.7267, 0.7195, 0.6997, 0.6758, 0.6500, 0.6221, 0.5894, 0.5602, 0.5278, 0.4919, 0.4526, 0.4154, 0.3701, 0.3285, 0.2795, 0.2306, 0.1832, 0.1375, 0.0907, 0.0432, 0.0086}},
      {{-2.9539, -3.0264, -3.1196, -3.2412, -3.4119, -3.6796, -4.2112, -5.4999, -9.4679, -35.2361, 13.5814, 5.6750, 3.6157, 2.6013, 2.0868, 1.7087, 1.4627, 1.2588, 1.1223, 1.0087, 0.9139, 0.8341, 0.7665, 0.7154, 0.6649, 0.6214, 0.5840, 0.5518, 0.5211, 0.4956, 0.4714, 0.4502}, {0.6573, 0.6808, 0.7088, 0.7413, 0.7812, 0.8341, 0.9095, 1.0034, 1.0511, 1.0106, 0.9573, 0.9497, 0.9678, 0.9991, 1.0330, 1.0770, 1.1281, 1.1930, 1.2576, 1.3389, 1.4416, 1.5732, 1.7457, 1.9440, 2.2487, 2.6802, 3.3114, 4.2097, 6.0100, 9.1342, 21.8391, 238.6966}, {1.3979, 1.3631, 1.3231, 1.2780, 1.2248, 1.1584, 1.0723, 0.9792, 0.9392, 0.9779, 1.0317, 1.0367, 1.0117, 0.9715, 0.9311, 0.8825, 0.8317, 0.7736, 0.7219, 0.6648, 0.6031, 0.5371, 0.4678, 0.4053, 0.3334, 0.2629, 0.1972, 0.1382, 0.0801, 0.0329, 0.0079, 0.0010}, {1.0069, 0.9779, 0.9456, 0.9109, 0.8720, 0.8258, 0.7684, 0.7081, 0.6819, 0.7055, 0.7396, 0.7460, 0.7361, 0.7185, 0.7006, 0.6788, 0.6553, 0.6282, 0.6038, 0.5761, 0.5451, 0.5106, 0.4724, 0.4360, 0.3912, 0.3437, 0.2952, 0.2481, 0.1943, 0.1449, 0.0836, 0.1207}, {0.9395, 0.9201, 0.8977, 0.8722, 0.8417, 0.8037, 0.7534, 0.6985, 0.6749, 0.6992, 0.7329, 0.7378, 0.7252, 0.7040, 0.6822, 0.6557, 0.6277, 0.5952, 0.5661, 0.5334, 0.4972, 0.4576, 0.4147, 0.3748, 0.3270, 0.2780, 0.2300, 0.1822, 0.1321, 0.0820, 0.0428, 0.0118}},
      {{-3.5249, -3.6462, -3.7992, -3.9965, -4.2669, -4.6744, -5.3908, -6.8932, -14.5560, 110.2536, 9.3310, 4.7890, 3.2638, 2.4807, 1.9935, 1.6571, 1.4306, 1.2567, 1.1100, 1.0012, 0.9095, 0.8408, 0.7732, 0.7151, 0.6710, 0.6270, 0.5933, 0.5596, 0.5303, 0.5021, 0.4807}, {0.7521, 0.7792, 0.8097, 0.8447, 0.8861, 0.9371, 1.0000, 1.0635, 1.0901, 1.0210, 0.9481, 0.9295, 0.9391, 0.9601, 0.9886, 1.0237, 1.0636, 1.1150, 1.1718, 1.2349, 1.3136, 1.3991, 1.5222, 1.6825, 1.8666, 2.1493, 2.4907, 3.0438, 3.8180, 5.3598, 7.5230}, {1.2723, 1.2348, 1.1938, 1.1495, 1.1002, 1.0444, 0.9823, 0.9271, 0.9060, 0.9673, 1.0407, 1.0582, 1.0414, 1.0113, 0.9732, 0.9295, 0.8842, 0.8319, 0.7785, 0.7255, 0.6676, 0.6130, 0.5466, 0.4766, 0.4134, 0.3404, 0.2772, 0.2094, 0.1473, 0.0851, 0.0374}, {0.9002, 0.8741, 0.8466, 0.8174, 0.7855, 0.7494, 0.7099, 0.6745, 0.6616, 0.6998, 0.7459, 0.7599, 0.7548, 0.7427, 0.7266, 0.7079, 0.6879, 0.6638, 0.6398, 0.6155, 0.5881, 0.5614, 0.5277, 0.4902, 0.4543, 0.4100, 0.3686, 0.3205, 0.2735, 0.2194, 0.1759}, {0.8690, 0.8477, 0.8241, 0.7986, 0.7699, 0.7370, 0.7002, 0.6670, 0.6552, 0.6934, 0.7387, 0.7509, 0.7432, 0.7278, 0.7076, 0.6844, 0.6599, 0.6311, 0.6018, 0.5723, 0.5394, 0.5079, 0.4686, 0.4259, 0.3861, 0.3385, 0.2956, 0.2475, 0.1989, 0.1473, 0.0985}},
      {{-4.4707, -4.6525, -4.8726, -5.1427, -5.4924, -5.9974, -6.9444, -9.6637, -37.0233, 20.8247, 6.6832, 3.9766, 2.8126, 2.1930, 1.7795, 1.5081, 1.3145, 1.1524, 1.0411, 0.9349, 0.8551, 0.7867, 0.7278, 0.6769, 0.6326, 0.5940, 0.5603, 0.5309, 0.5022, 0.4776, 0.4275, 0.4236}, {0.8761, 0.9031, 0.9326, 0.9645, 0.9998, 1.0393, 1.0879, 1.1447, 1.1424, 1.0409, 0.9392, 0.9121, 0.9171, 0.9351, 0.9619, 0.9928, 1.0312, 1.0818, 1.1267, 1.1900, 1.2604, 1.3488, 1.4612, 1.6076, 1.8026, 2.0703, 2.4532, 2.9916, 3.9243, 5.5533, 226.0122, 536.9089}, {1.1142, 1.0832, 1.0508, 1.0176, 0.9836, 0.9475, 0.9067, 0.8624, 0.8630, 0.9475, 1.0492, 1.0764, 1.0639, 1.0352, 0.9965, 0.9549, 0.9082, 0.8526, 0.8064, 0.7481, 0.6910, 0.6290, 0.5624, 0.4921, 0.4188, 0.3443, 0.2705, 0.2021, 0.1296, 0.0675, 0.0010, 0.0004}, {0.7929, 0.7730, 0.7525, 0.7316, 0.7100, 0.6873, 0.6616, 0.6343, 0.6361, 0.6890, 0.7522, 0.7722, 0.7703, 0.7596, 0.7441, 0.7268, 0.7066, 0.6818, 0.6618, 0.6360, 0.6100, 0.5808, 0.5480, 0.5114, 0.4710, 0.4267, 0.3787, 0.3314, 0.2754, 0.2245, 0.0255, 0.0134}, {0.7781, 0.7597, 0.7406, 0.7210, 0.7007, 0.6791, 0.6545, 0.6278, 0.6293, 0.6818, 0.7442, 0.7623, 0.7573, 0.7429, 0.7228, 0.7010, 0.6760, 0.6459, 0.6209, 0.5890, 0.5572, 0.5220, 0.4833, 0.4410, 0.3957, 0.3477, 0.2985, 0.2501, 0.1948, 0.1409, 0.0425, 0.0425}},
      {{-5.6847, -5.8994, -6.1763, -6.5419, -7.0646, -7.9195, -9.7590, -15.2068, 100.7049, 11.6009, 5.3758, 3.5226, 2.5935, 2.0600, 1.7151, 1.4484, 1.2780, 1.1224, 1.0198, 0.9205, 0.8448, 0.7793, 0.7225, 0.6731, 0.6299, 0.5921, 0.5589, 0.5298, 0.5041, 0.4791}, {0.9990, 1.0215, 1.0473, 1.0763, 1.1100, 1.1503, 1.2014, 1.2525, 1.2314, 1.0810, 0.9378, 0.8987, 0.8965, 0.9090, 0.9293, 0.9573, 0.9862, 1.0338, 1.0709, 1.1225, 1.1796, 1.2502, 1.3388, 1.4523, 1.6010, 1.8015, 2.0826, 2.4888, 3.0455, 4.1534}, {0.9851, 0.9642, 0.9412, 0.9164, 0.8887, 0.8576, 0.8211, 0.7868, 0.7984, 0.9098, 1.0485, 1.0908, 1.0870, 1.0639, 1.0315, 0.9902, 0.9508, 0.8934, 0.8508, 0.7965, 0.7428, 0.6841, 0.6205, 0.5522, 0.4801, 0.4051, 0.3286, 0.2540, 0.1845, 0.1106}, {0.7104, 0.6973, 0.6830, 0.6675, 0.6504, 0.6311, 0.6087, 0.5877, 0.5969, 0.6674, 0.7532, 0.7819, 0.7851, 0.7778, 0.7654, 0.7490, 0.7330, 0.7076, 0.6899, 0.6669, 0.6436, 0.6174, 0.5879, 0.5547, 0.5174, 0.4760, 0.4301, 0.3812, 0.3336, 0.2754}, {0.7015, 0.6891, 0.6751, 0.6602, 0.6435, 0.6246, 0.6021, 0.5812, 0.5896, 0.6595, 0.7444, 0.7711, 0.7711, 0.7601, 0.7436, 0.7223, 0.7018, 0.6709, 0.6483, 0.6193, 0.5902, 0.5577, 0.5218, 0.4823, 0.4392, 0.3928, 0.3437, 0.2944, 0.2441, 0.1888}},
      {{-7.2464, -7.6348, -8.1880, -8.9646, -10.1647, -12.3625, -17.9308, -38.8199, 17.5131, 7.1798, 4.1571, 2.9035, 2.2436, 1.8353, 1.5528, 1.3293, 1.1757, 1.0481, 0.9482, 0.8710, 0.7958, 0.7380, 0.6876, 0.6433, 0.6044, 0.5662, 0.5364, 0.5101, 0.4842, 0.4618}, {1.1237, 1.1536, 1.1875, 1.2258, 1.2704, 1.3221, 1.3814, 1.4254, 1.3674, 1.1115, 0.9251, 0.8828, 0.8806, 0.8918, 0.9108, 0.9369, 0.9682, 1.0100, 1.0493, 1.0922, 1.1520, 1.2182, 1.3010, 1.4069, 1.5461, 1.7627, 2.0433, 2.4409, 3.0997, 4.3219}, {0.8778, 0.8550, 0.8307, 0.8043, 0.7759, 0.7447, 0.7119, 0.6890, 0.7163, 0.8821, 1.0609, 1.1070, 1.1029, 1.0804, 1.0482, 1.0076, 0.9634, 0.9103, 0.8625, 0.8146, 0.7547, 0.6961, 0.6325, 0.5642, 0.4916, 0.4056, 0.3270, 0.2513, 0.1697, 0.0951}, {0.6435, 0.6295, 0.6142, 0.5979, 0.5802, 0.5609, 0.5403, 0.5262, 0.5461, 0.6522, 0.7625, 0.7942, 0.7975, 0.7907, 0.7789, 0.7632, 0.7451, 0.7223, 0.7028, 0.6831, 0.6578, 0.6322, 0.6036, 0.5711, 0.5345, 0.4876, 0.4404, 0.3915, 0.3347, 0.2767}, {0.6367, 0.6227, 0.6077, 0.5914, 0.5737, 0.5543, 0.5336, 0.5191, 0.5378, 0.6432, 0.7523, 0.7815, 0.7814, 0.7707, 0.7547, 0.7340, 0.7109, 0.6828, 0.6578, 0.6324, 0.6003, 0.5683, 0.5328, 0.4937, 0.4508, 0.3982, 0.3483, 0.2985, 0.2398, 0.1838}},
      {{-10.0150, -10.9842, -12.4673, -14.8170, -19.4137, -30.6466, -119.6172, 70.7214, 9.8707, 5.3325, 3.4992, 2.6057, 2.0653, 1.7179, 1.4602, 1.2902, 1.1356, 1.0276, 0.9337, 0.8520, 0.7886, 0.7264, 0.6783, 0.6359, 0.5984, 0.5653, 0.5359, 0.5098, 0.4841, 0.4618, 0.4404}, {1.2989, 1.3522, 1.4054, 1.4649, 1.5340, 1.6058, 1.6753, 1.7069, 1.5939, 1.1570, 0.9157, 0.8712, 0.8663, 0.8762, 0.8937, 0.9135, 0.9420, 0.9774, 1.0114, 1.0534, 1.0991, 1.1628, 1.2339, 1.3236, 1.4405, 1.5968, 1.8168, 2.1356, 2.6413, 3.5869, 7.0522}, {0.7588, 0.7280, 0.7002, 0.6708, 0.6401, 0.6106, 0.5840, 0.5727, 0.6114, 0.8450, 1.0691, 1.1196, 1.1188, 1.0982, 1.0662, 1.0334, 0.9899, 0.9422, 0.8971, 0.8459, 0.7951, 0.7322, 0.6708, 0.6044, 0.5331, 0.4573, 0.3779, 0.2980, 0.2125, 0.1293, 0.0268}, {0.5691, 0.5498, 0.5321, 0.5134, 0.4936, 0.4747, 0.4577, 0.4505, 0.4787, 0.6308, 0.7692, 0.8033, 0.8085, 0.8026, 0.7915, 0.7791, 0.7621, 0.7417, 0.7240, 0.7035, 0.6828, 0.6564, 0.6297, 0.5995, 0.5653, 0.5262, 0.4814, 0.4319, 0.3759, 0.3124, 0.2050}, {0.5627, 0.5432, 0.5255, 0.5066, 0.4868, 0.4675, 0.4501, 0.4427, 0.4696, 0.6207, 0.7578, 0.7894, 0.7912, 0.7817, 0.7658, 0.7492, 0.7270, 0.7018, 0.6786, 0.6519, 0.6251, 0.5914, 0.5578, 0.5206, 0.4795, 0.4344, 0.3853, 0.3345, 0.2761, 0.2174, 0.1883}},
      {{-21.2540, -28.4010, -46.1467, -159.2357, 95.6938, 34.9773, 21.4316, 15.9515, 6.2085, 3.8463, 2.7874, 2.1896, 1.7941, 1.5253, 1.3309, 1.1719, 1.0532, 0.9514, 0.8700, 0.7991, 0.7436, 0.6887, 0.6459, 0.6080, 0.5742, 0.5408, 0.5145, 0.4883, 0.4656, 0.4459}, {1.7098, 1.8497, 1.9671, 2.0957, 2.2242, 2.3365, 2.3986, 2.4062, 2.0864, 1.1128, 0.8847, 0.8537, 0.8535, 0.8645, 0.8818, 0.9037, 0.9275, 0.9631, 0.9967, 1.0374, 1.0817, 1.1429, 1.2112, 1.2980, 1.4109, 1.5878, 1.8183, 2.1923, 2.7839, 3.9070}, {0.5734, 0.5290, 0.4965, 0.4652, 0.4374, 0.4156, 0.4042, 0.4026, 0.4631, 0.8760, 1.1036, 1.1382, 1.1309, 1.1074, 1.0768, 1.0387, 1.0005, 0.9499, 0.9038, 0.8521, 0.8011, 0.7377, 0.6761, 0.6091, 0.5368, 0.4496, 0.3676, 0.2770, 0.1873, 0.1053}, {0.4762, 0.4636, 0.4525, 0.4411, 0.4299, 0.4193, 0.4101, 0.4031, 0.3787, 0.6522, 0.7923, 0.8175, 0.8193, 0.8123, 0.8011, 0.7872, 0.7729, 0.7520, 0.7342, 0.7139, 0.6935, 0.6675, 0.6413, 0.6113, 0.5771, 0.5320, 0.4850, 0.4287, 0.3673, 0.3039}, {0.4760, 0.4634, 0.4523, 0.4409, 0.4297, 0.4190, 0.4099, 0.4029, 0.3710, 0.6407, 0.7791, 0.8014, 0.7995, 0.7884, 0.7733, 0.7541, 0.7348, 0.7086, 0.6849, 0.6582, 0.6315, 0.5979, 0.5644, 0.5272, 0.4858, 0.4340, 0.3834, 0.3255, 0.2635, 0.2052}},
      {{32.9489, 23.0256, 17.4034, 13.8619, 11.7371, 10.3072, 9.1726, 8.1453, 6.2861, 5.1994, 3.7629, 2.9479, 2.2932, 1.8672, 1.5731, 1.3679, 1.2031, 1.0705, 0.9672, 0.8855, 0.8142, 0.7516, 0.7027, 0.6541, 0.6160, 0.5821, 0.5482, 0.5215, 0.4949, 0.4715, 0.4512, 0.4407}, {2.4334, 2.5418, 2.6331, 2.7233, 2.7964, 2.8527, 2.8979, 2.9328, 2.9826, 3.2602, 2.7563, 1.0090, 0.8536, 0.8358, 0.8403, 0.8528, 0.8710, 0.8931, 0.9236, 0.9511, 0.9843, 1.0247, 1.0680, 1.1280, 1.1943, 1.2789, 1.4066, 1.5663, 1.8339, 2.2191, 2.9159, 4.0438}, {0.4105, 0.3929, 0.3793, 0.3667, 0.3570, 0.3500, 0.3445, 0.3404, 0.3347, 0.3061, 0.3464, 0.9652, 1.1403, 1.1579, 1.1435, 1.1174, 1.0838, 1.0449, 0.9978, 0.9558, 0.9086, 0.8559, 0.8042, 0.7403, 0.6782, 0.6107, 0.5282, 0.4488, 0.3548, 0.2630, 0.1695, 0.1034}, {0.4410, 0.4305, 0.4218, 0.4135, 0.4069, 0.4018, 0.3973, 0.3933, 0.3845, 0.3665, 0.3538, 0.7091, 0.8169, 0.8327, 0.8307, 0.8223, 0.8102, 0.7959, 0.7767, 0.7610, 0.7430, 0.7227, 0.7024, 0.6765, 0.6505, 0.6207, 0.5816, 0.5402, 0.4857, 0.4285, 0.3601, 0.2945}, {0.4408, 0.4303, 0.4215, 0.4133, 0.4066, 0.4015, 0.3971, 0.3930, 0.3842, 0.3662, 0.3534, 0.6961, 0.8016, 0.8140, 0.8081, 0.7956, 0.7789, 0.7595, 0.7352, 0.7140, 0.6901, 0.6630, 0.6362, 0.6024, 0.5689, 0.5316, 0.4845, 0.4373, 0.3794, 0.3199, 0.2566, 0.2880}},
      {{9.0522, 7.6658, 6.4570, 5.4422, 4.6781, 4.1799, 3.8043, 3.1968, 2.4843, 2.2347, 1.8375, 1.5538, 1.3456, 1.1986, 1.0629, 0.9661, 0.8810, 0.8106, 0.7563, 0.7022, 0.6598, 0.6174, 0.5840, 0.5506, 0.5241, 0.4975, 0.4741, 0.4538, 0.4363}, {2.8168, 2.9553, 3.0803, 3.1873, 3.2465, 3.2359, 3.1539, 3.8705, 0.9620, 0.8611, 0.8224, 0.8239, 0.8333, 0.8482, 0.8688, 0.8890, 0.9181, 0.9478, 0.9781, 1.0190, 1.0628, 1.1229, 1.1906, 1.2895, 1.4103, 1.6083, 1.8978, 2.3467, 3.5611}, {0.3543, 0.3377, 0.3239, 0.3130, 0.3073, 0.3083, 0.3163, 0.2576, 1.0106, 1.1294, 1.1774, 1.1672, 1.1449, 1.1162, 1.0772, 1.0413, 0.9948, 0.9491, 0.9053, 0.8510, 0.7980, 0.7330, 0.6696, 0.5915, 0.5156, 0.4227, 0.3281, 0.2329, 0.1257}, {0.4134, 0.4022, 0.3921, 0.3831, 0.3765, 0.3735, 0.3732, 0.3394, 0.7383, 0.8106, 0.8437, 0.8439, 0.8376, 0.8273, 0.8134, 0.8005, 0.7823, 0.7652, 0.7488, 0.7281, 0.7075, 0.6815, 0.6549, 0.6203, 0.5837, 0.5338, 0.4778, 0.4157, 0.3210}, {0.4130, 0.4018, 0.3917, 0.3827, 0.3762, 0.3731, 0.3728, 0.3390, 0.7239, 0.7951, 0.8249, 0.8213, 0.8110, 0.7970, 0.7775, 0.7598, 0.7362, 0.7131, 0.6910, 0.6633, 0.6359, 0.6016, 0.5675, 0.5242, 0.4806, 0.4250, 0.3661, 0.3025, 0.2226}},
      {{6.4367, 5.6405, 4.7822, 4.0725, 3.5516, 3.1953, 2.9209, 2.7959, 2.6577, 2.0635, 1.7086, 1.4645, 1.2778, 1.1354, 1.0180, 0.9291, 0.8457, 0.7846, 0.7240, 0.6755, 0.6324, 0.5981, 0.5635, 0.5328, 0.5056, 0.4814, 0.4604, 0.4413, 0.4235}, {2.9920, 3.1138, 3.2481, 3.3718, 3.4542, 3.4700, 3.3967, 3.3830, 3.9222, 0.9766, 0.8341, 0.8189, 0.8232, 0.8359, 0.8540, 0.8835, 0.9095, 0.9358, 0.9715, 1.0110, 1.0590, 1.1119, 1.1872, 1.2888, 1.4339, 1.6469, 1.9493, 2.8380, 4.0068}, {0.3332, 0.3200, 0.3068, 0.2956, 0.2886, 0.2873, 0.2936, 0.2947, 0.2540, 0.9913, 1.1581, 1.1723, 1.1569, 1.1292, 1.0939, 1.0454, 1.0012, 0.9593, 0.9063, 0.8529, 0.7937, 0.7356, 0.6639, 0.5845, 0.4964, 0.4019, 0.3071, 0.1818, 0.0948}, {0.4030, 0.3934, 0.3832, 0.3738, 0.3666, 0.3623, 0.3610, 0.3591, 0.3380, 0.7298, 0.8340, 0.8484, 0.8465, 0.8378, 0.8253, 0.8049, 0.7889, 0.7736, 0.7540, 0.7338, 0.7111, 0.6881, 0.6582, 0.6227, 0.5792, 0.5273, 0.4710, 0.3674, 0.3053}, {0.4024, 0.3928, 0.3826, 0.3733, 0.3661, 0.3618, 0.3606, 0.3586, 0.3375, 0.7136, 0.8142, 0.8247, 0.8185, 0.8052, 0.7879, 0.7627, 0.7408, 0.7199, 0.6933, 0.6661, 0.6356, 0.6052, 0.5667, 0.5226, 0.4716, 0.4145, 0.3544, 0.2655, 0.2035}},
      {{5.9007, 5.0038, 4.1549, 3.4864, 2.9610, 2.5865, 2.3734, 2.2204, 1.8567, 1.8517, 1.5759, 1.3723, 1.2158, 1.0898, 0.9913, 0.9046, 0.8280, 0.7663, 0.7105, 0.6604, 0.6159, 0.5802, 0.5483, 0.5168, 0.4918, 0.4674, 0.4562, 0.4368, 0.4200}, {2.9554, 3.0845, 3.2411, 3.4070, 3.5741, 3.6823, 3.6670, 3.9804, 3.3263, 0.9654, 0.8314, 0.8188, 0.8223, 0.8330, 0.8482, 0.8691, 0.8982, 0.9246, 0.9573, 0.9980, 1.0489, 1.1067, 1.1811, 1.2948, 1.4452, 1.6907, 2.0405, 2.6437, 3.7229}, {0.3376, 0.3234, 0.3076, 0.2924, 0.2783, 0.2704, 0.2716, 0.2500, 0.2991, 1.0007, 1.1594, 1.1699, 1.1565, 1.1318, 1.1013, 1.0629, 1.0143, 0.9709, 0.9208, 0.8636, 0.7993, 0.7347, 0.6633, 0.5746, 0.4856, 0.3809, 0.2946, 0.1982, 0.1067}, {0.4060, 0.3956, 0.3838, 0.3722, 0.3611, 0.3527, 0.3494, 0.3352, 0.3406, 0.7371, 0.8366, 0.8487, 0.8474, 0.8402, 0.8299, 0.8158, 0.7967, 0.7811, 0.7628, 0.7416, 0.7171, 0.6917, 0.6621, 0.6223, 0.5777, 0.5194, 0.4528, 0.3859, 0.3195}, {0.4056, 0.3951, 0.3833, 0.3716, 0.3604, 0.3521, 0.3489, 0.3346, 0.3399, 0.7197, 0.8155, 0.8238, 0.8186, 0.8070, 0.7922, 0.7731, 0.7484, 0.7269, 0.7018, 0.6730, 0.6400, 0.6063, 0.5680, 0.5188, 0.4670, 0.4032, 0.3411, 0.2800, 0.2168}},
      {{5.3637, 4.4871, 3.7266, 3.1102, 2.6446, 2.2754, 2.0051, 1.8811, 1.7799, 1.5688, 1.3670, 1.1996, 1.0695, 0.9594, 0.8744, 0.8037, 0.7918, 0.7368, 0.6848, 0.6412, 0.6053, 0.5725, 0.5396, 0.5129, 0.4890, 0.4709, 0.4494, 0.4310, 0.4139}, {2.8863, 3.0193, 3.1844, 3.3714, 3.5604, 3.7481, 3.8918, 4.0038, 4.3911, 3.6250, 0.8422, 0.8311, 0.8338, 0.8442, 0.8605, 0.8819, 0.8954, 0.9197, 0.9515, 0.9881, 1.0289, 1.0799, 1.1536, 1.2439, 1.3723, 1.6351, 1.9118, 2.4944, 3.6956}, {0.3460, 0.3307, 0.3135, 0.2959, 0.2800, 0.2655, 0.2554, 0.2483, 0.2259, 0.2742, 1.1383, 1.1451, 1.1314, 1.1056, 1.0724, 1.0327, 1.0141, 0.9729, 0.9230, 0.8704, 0.8167, 0.7566, 0.6817, 0.6058, 0.5206, 0.4097, 0.3160, 0.2106, 0.1016}, {0.4113, 0.4001, 0.3873, 0.3739, 0.3615, 0.3498, 0.3404, 0.3335, 0.3172, 0.3292, 0.8286, 0.8395, 0.8392, 0.8330, 0.8225, 0.8089, 0.7994, 0.7850, 0.7672, 0.7480, 0.7279, 0.7048, 0.6744, 0.6412, 0.6000, 0.5274, 0.4761, 0.4024, 0.3217}, {0.4111, 0.3999, 0.3870, 0.3736, 0.3610, 0.3491, 0.3396, 0.3327, 0.3163, 0.3285, 0.8048, 0.8113, 0.8064, 0.7948, 0.7789, 0.7596, 0.7492, 0.7290, 0.7042, 0.6778, 0.6505, 0.6195, 0.5798, 0.5383, 0.4895, 0.4170, 0.3590, 0.2913, 0.2179}},
      {{4.8513, 3.9866, 3.3513, 2.8765, 2.4572, 2.1142, 1.8466, 1.7203, 1.6012, 1.3715, 1.2305, 1.0995, 0.9969, 0.9119, 0.8400, 0.7748, 0.7228, 0.6847, 0.6376, 0.5976, 0.5611, 0.5486, 0.5204, 0.4891, 0.4705, 0.4491, 0.4297, 0.4115}, {2.8279, 2.9807, 3.1487, 3.3199, 3.5173, 3.7200, 3.9100, 4.0473, 4.4401, 3.8408, 0.8581, 0.8445, 0.8443, 0.8499, 0.8619, 0.8806, 0.9024, 0.9334, 0.9694, 1.0123, 1.0685, 1.1982, 1.2704, 1.3918, 1.5024, 1.7264, 2.2283, 3.3797}, {0.3533, 0.3352, 0.3172, 0.3008, 0.2837, 0.2679, 0.2544, 0.2455, 0.2238, 0.2585, 1.1118, 1.1216, 1.1127, 1.0948, 1.0680, 1.0316, 0.9925, 0.9464, 0.8918, 0.8327, 0.7643, 0.6698, 0.6084, 0.5219, 0.4576, 0.3627, 0.2443, 0.1167}, {0.4159, 0.4027, 0.3893, 0.3768, 0.3635, 0.3507, 0.3393, 0.3317, 0.3142, 0.3208, 0.8167, 0.8293, 0.8314, 0.8289, 0.8218, 0.8102, 0.7971, 0.7780, 0.7586, 0.7370, 0.7110, 0.6501, 0.6261, 0.5902, 0.5613, 0.5119, 0.4330, 0.3378}, {0.4157, 0.4025, 0.3891, 0.3765, 0.3631, 0.3502, 0.3386, 0.3309, 0.3134, 0.3199, 0.7907, 0.7990, 0.7966, 0.7893, 0.7770, 0.7596, 0.7406, 0.7167, 0.6896, 0.6598, 0.6248, 0.5656, 0.5334, 0.4859, 0.4491, 0.3915, 0.3152, 0.2300}},
      {{4.7578, 3.9609, 3.2433, 2.7923, 2.3922, 2.0869, 1.8207, 1.6988, 1.5675, 1.3388, 1.2052, 1.0792, 0.9746, 0.8931, 0.8191, 0.7607, 0.7089, 0.6905, 0.6803, 0.6337, 0.5970, 0.5625, 0.5334, 0.5090, 0.4876, 0.4618, 0.4436, 0.4267, 0.4067}, {2.8176, 2.9604, 3.1571, 3.3293, 3.5274, 3.7151, 3.9094, 4.0487, 4.4526, 3.8498, 0.8613, 0.8473, 0.8468, 0.8526, 0.8661, 0.8842, 0.9071, 1.0244, 1.0320, 1.0670, 1.1071, 1.1609, 1.2276, 1.3096, 1.4123, 1.6128, 1.8970, 2.4870, 4.5342}, {0.3546, 0.3375, 0.3164, 0.3000, 0.2829, 0.2684, 0.2546, 0.2455, 0.2233, 0.2578, 1.1065, 1.1167, 1.1073, 1.0890, 1.0594, 1.0247, 0.9839, 0.8666, 0.8563, 0.8108, 0.7641, 0.7084, 0.6480, 0.5863, 0.5200, 0.4203, 0.3272, 0.2201, 0.0236}, {0.4167, 0.4042, 0.3885, 0.3760, 0.3627, 0.3509, 0.3392, 0.3316, 0.3134, 0.3201, 0.8143, 0.8273, 0.8297, 0.8271, 0.8191, 0.8081, 0.7945, 0.7215, 0.7178, 0.7028, 0.6865, 0.6660, 0.6425, 0.6160, 0.5860, 0.5368, 0.4821, 0.4059, 0.2902}, {0.4166, 0.4041, 0.3883, 0.3758, 0.3624, 0.3505, 0.3386, 0.3308, 0.3127, 0.3192, 0.7878, 0.7964, 0.7940, 0.7865, 0.7729, 0.7563, 0.7365, 0.6671, 0.6617, 0.6396, 0.6164, 0.5882, 0.5568, 0.5235, 0.4864, 0.4279, 0.3693, 0.2998, 0.1505}},
      {{3.8994, 3.2369, 2.7718, 2.4351, 2.1191, 1.8578, 1.6367, 1.5307, 1.4280, 1.2257, 1.0751, 0.9772, 0.8874, 0.8126, 0.7495, 0.6950, 0.6481, 0.6282, 0.5989, 0.5562, 0.5294, 0.5048, 0.4782, 0.4571, 0.4382, 0.4200, 0.4039, 0.3887}, {2.7128, 2.8839, 3.0622, 3.2317, 3.4325, 3.6365, 3.8429, 3.9942, 4.3680, 4.6077, 3.9781, 0.9017, 0.8843, 0.8824, 0.8907, 0.9087, 0.9353, 1.0190, 1.0457, 1.1022, 1.1545, 1.2203, 1.3211, 1.4535, 1.6884, 2.2320, 3.5277, 19.2706}, {0.3684, 0.3467, 0.3264, 0.3092, 0.2911, 0.2747, 0.2596, 0.2496, 0.2282, 0.2155, 0.2487, 1.0462, 1.0570, 1.0479, 1.0252, 0.9902, 0.9457, 0.8578, 0.8226, 0.7563, 0.6877, 0.6407, 0.5629, 0.4810, 0.3805, 0.2508, 0.0230, 0.4810}, {0.4253, 0.4094, 0.3944, 0.3814, 0.3673, 0.3541, 0.3416, 0.3332, 0.3157, 0.3033, 0.3128, 0.7870, 0.8021, 0.8061, 0.8028, 0.7931, 0.7785, 0.7275, 0.7157, 0.6922, 0.6719, 0.6483, 0.6158, 0.5784, 0.5236, 0.4359, 0.3283, 0.0859}, {0.4252, 0.4093, 0.3943, 0.3813, 0.3672, 0.3540, 0.3412, 0.3328, 0.3153, 0.3025, 0.3115, 0.7563, 0.7661, 0.7645, 0.7553, 0.7392, 0.7179, 0.6661, 0.6490, 0.6161, 0.5817, 0.5565, 0.5141, 0.4675, 0.4061, 0.3213, 0.1566, -0.0101}},
      {{3.1397, 2.7022, 2.3747, 2.1055, 1.8629, 1.6588, 1.4826, 1.3895, 1.3302, 1.1709, 1.0486, 0.9530, 0.8675, 0.7967, 0.7371, 0.6850, 0.6390, 0.6003, 0.5761, 0.5437, 0.5169, 0.4906, 0.4686, 0.4510, 0.4322, 0.4132}, {2.5916, 2.7574, 2.9292, 3.1093, 3.3095, 3.5127, 3.7213, 3.8790, 4.1895, 4.4505, 4.5954, 4.6344, 3.9683, 0.9649, 0.9302, 0.9264, 0.9295, 0.9472, 1.0412, 1.0838, 1.1313, 1.1945, 1.2708, 1.3648, 1.5498, 1.9772}, {0.3857, 0.3625, 0.3412, 0.3215, 0.3020, 0.2846, 0.2685, 0.2575, 0.2385, 0.2243, 0.2164, 0.2135, 0.2480, 0.9585, 0.9835, 0.9744, 0.9558, 0.9213, 0.8281, 0.7763, 0.7235, 0.6602, 0.5946, 0.5286, 0.4230, 0.2954}, {0.4357, 0.4190, 0.4034, 0.3886, 0.3737, 0.3600, 0.3469, 0.3378, 0.3224, 0.3097, 0.3012, 0.2962, 0.3124, 0.7465, 0.7721, 0.7777, 0.7790, 0.7710, 0.7178, 0.6997, 0.6807, 0.6574, 0.6320, 0.6037, 0.5556, 0.4724}, {0.4356, 0.4189, 0.4033, 0.3886, 0.3737, 0.3599, 0.3468, 0.3377, 0.3223, 0.3094, 0.3005, 0.2950, 0.3103, 0.7080, 0.7272, 0.7262, 0.7200, 0.7047, 0.6517, 0.6262, 0.5995, 0.5667, 0.5320, 0.4957, 0.4339, 0.3513}},
      {{2.6547, 2.3497, 2.0909, 1.8718, 1.6820, 1.5113, 1.3609, 1.2912, 1.2432, 1.1117, 1.0077, 0.9169, 0.8435, 0.7805, 0.7322, 0.7268, 0.6824, 0.6416, 0.6049, 0.5741, 0.5452, 0.5363, 0.5098, 0.4856, 0.4637, 0.4445, 0.4249, 0.4100, 0.3936}, {2.4835, 2.6418, 2.8169, 3.0015, 3.1946, 3.4009, 3.6182, 3.7672, 4.0396, 4.3094, 4.4926, 4.6056, 4.6353, 4.6111, 3.8104, 1.5972, 1.0669, 1.0022, 0.9932, 0.9950, 1.0097, 1.0721, 1.1100, 1.1603, 1.2278, 1.3208, 1.4948, 1.7236, 2.3615}, {0.4025, 0.3784, 0.3548, 0.3330, 0.3129, 0.2938, 0.2762, 0.2652, 0.2473, 0.2318, 0.2222, 0.2163, 0.2138, 0.2131, 0.2568, 0.5661, 0.8515, 0.8968, 0.8921, 0.8770, 0.8483, 0.7901, 0.7430, 0.6880, 0.6251, 0.5549, 0.4594, 0.3681, 0.3681}, {0.4458, 0.4286, 0.4115, 0.3954, 0.3802, 0.3656, 0.3515, 0.3428, 0.3281, 0.3151, 0.3057, 0.2986, 0.2943, 0.2928, 0.3206, 0.4877, 0.6905, 0.7303, 0.7391, 0.7415, 0.7370, 0.7025, 0.6875, 0.6685, 0.6448, 0.6153, 0.5680, 0.5171, 0.4177}, {0.4457, 0.4285, 0.4114, 0.3953, 0.3801, 0.3654, 0.3514, 0.3426, 0.3280, 0.3150, 0.3055, 0.2981, 0.2932, 0.2908, 0.3175, 0.4534, 0.6475, 0.6802, 0.6822, 0.6780, 0.6663, 0.6309, 0.6077, 0.5797, 0.5470, 0.5092, 0.4548, 0.3989, -0.0101}},
      {{2.2998, 2.0706, 1.8627, 1.7047, 1.5467, 1.3929, 1.2750, 1.2086, 1.1493, 1.0274, 0.9280, 0.8450, 0.7765, 0.7163, 0.6677, 0.6230, 0.5846, 0.5513, 0.5225, 0.4951, 0.4719, 0.4503, 0.4302, 0.4117}, {2.3755, 2.5324, 2.7100, 2.8720, 3.0629, 3.2828, 3.4830, 3.6390, 3.9400, 4.2336, 4.4529, 4.6028, 4.6723, 4.6928, 4.7124, 1.8308, 1.1385, 1.0813, 1.0759, 1.0854, 1.2075, 1.2936, 1.4281, 1.6212}, {0.4208, 0.3947, 0.3688, 0.3480, 0.3263, 0.3044, 0.2869, 0.2746, 0.2536, 0.2359, 0.2242, 0.2168, 0.2128, 0.2103, 0.2069, 0.4795, 0.7779, 0.8058, 0.7932, 0.7655, 0.6638, 0.5966, 0.5144, 0.4232}, {0.4567, 0.4382, 0.4196, 0.4045, 0.3885, 0.3721, 0.3586, 0.3489, 0.3318, 0.3175, 0.3066, 0.2983, 0.2928, 0.2898, 0.2894, 0.4369, 0.6603, 0.6933, 0.7003, 0.7007, 0.6480, 0.6203, 0.5816, 0.5357}, {0.4566, 0.4381, 0.4196, 0.4044, 0.3884, 0.3720, 0.3584, 0.3488, 0.3317, 0.3173, 0.3064, 0.2981, 0.2921, 0.2883, 0.2864, 0.3971, 0.6093, 0.6335, 0.6319, 0.6221, 0.5636, 0.5284, 0.4830, 0.4303}},
      {{2.0241, 1.8572, 1.6921, 1.5456, 1.4151, 1.2986, 1.1833, 1.1313, 1.0813, 0.9763, 0.8896, 0.8167, 0.7552, 0.7006, 0.6551, 0.6126, 0.5776, 0.5454, 0.5187, 0.4938, 0.4726, 0.4521, 0.4386, 0.4203, 0.4037, 0.3884, 0.3870, 0.3850, 0.3810, 0.3742}, {2.2632, 2.4078, 2.5791, 2.7593, 2.9472, 3.1424, 3.3722, 3.5220, 3.7798, 4.0813, 4.3241, 4.5171, 4.6551, 4.7434, 4.7973, 4.8745, 4.9765, 1.9097, 1.2350, 1.1837, 1.1803, 1.2049, 1.3567, 1.4666, 1.6438, 1.9581, 1.9626, 2.0165, 2.2308, 3.1383}, {0.4417, 0.4151, 0.3875, 0.3622, 0.3391, 0.3180, 0.2963, 0.2837, 0.2643, 0.2447, 0.2309, 0.2209, 0.2143, 0.2099, 0.2062, 0.2007, 0.1931, 0.4471, 0.6946, 0.7109, 0.6973, 0.6641, 0.5703, 0.5022, 0.4201, 0.3242, 0.3241, 0.3099, 0.2650, 0.1620}, {0.4690, 0.4505, 0.4309, 0.4128, 0.3961, 0.3806, 0.3644, 0.3547, 0.3388, 0.3236, 0.3120, 0.3026, 0.2952, 0.2896, 0.2863, 0.2843, 0.2833, 0.4263, 0.6226, 0.6497, 0.6563, 0.6524, 0.5969, 0.5683, 0.5279, 0.4713, 0.4714, 0.4632, 0.4325, 0.3448}, {0.4689, 0.4504, 0.4308, 0.4127, 0.3960, 0.3805, 0.3642, 0.3546, 0.3386, 0.3235, 0.3118, 0.3023, 0.2949, 0.2891, 0.2850, 0.2817, 0.2786, 0.3806, 0.5633, 0.5812, 0.5796, 0.5672, 0.5095, 0.4727, 0.4255, 0.3666, 0.3673, 0.3580, 0.3266, 0.2538}},
      {{1.8042, 1.6661, 1.5397, 1.4271, 1.3165, 1.2038, 1.1036, 1.0582, 0.9979, 0.9147, 0.8476, 0.7803, 0.7252, 0.6766, 0.6369, 0.6004, 0.5660, 0.5345, 0.5073, 0.4843, 0.4653, 0.4461, 0.4297, 0.4141, 0.4057, 0.3891, 0.3739}, {2.1476, 2.2961, 2.4543, 2.6169, 2.8008, 3.0194, 3.2536, 3.4060, 3.6834, 3.9592, 4.2113, 4.4606, 4.6577, 4.8197, 4.9417, 5.0490, 5.1754, 5.3572, 5.5495, 1.9697, 1.3768, 1.3432, 1.3487, 1.3794, 1.5596, 1.8006, 2.4231}, {0.4655, 0.4353, 0.4073, 0.3819, 0.3568, 0.3309, 0.3070, 0.2932, 0.2710, 0.2521, 0.2369, 0.2236, 0.2140, 0.2067, 0.2014, 0.1962, 0.1895, 0.1801, 0.1692, 0.4167, 0.5995, 0.6004, 0.5806, 0.5455, 0.4630, 0.3710, 0.2452}, {0.4830, 0.4622, 0.4427, 0.4249, 0.4071, 0.3885, 0.3711, 0.3608, 0.3427, 0.3282, 0.3163, 0.3049, 0.2960, 0.2884, 0.2827, 0.2781, 0.2745, 0.2715, 0.2693, 0.4214, 0.5768, 0.5944, 0.5981, 0.5944, 0.5425, 0.4949, 0.4095}, {0.4829, 0.4621, 0.4426, 0.4248, 0.4069, 0.3884, 0.3709, 0.3606, 0.3424, 0.3279, 0.3159, 0.3046, 0.2956, 0.2880, 0.2821, 0.2771, 0.2723, 0.2675, 0.2623, 0.3675, 0.5102, 0.5192, 0.5140, 0.4997, 0.4465, 0.3929, 0.3142}},
      {{1.6399, 1.5300, 1.4183, 1.3131, 1.2181, 1.1304, 1.0406, 0.9920, 0.9441, 0.8745, 0.8013, 0.7449, 0.6921, 0.6463, 0.6105, 0.5756, 0.5426, 0.5156, 0.4913, 0.4675, 0.4477, 0.4298, 0.4150, 0.4003, 0.3865, 0.3731, 0.3590}, {2.0378, 2.1780, 2.3411, 2.5178, 2.7017, 2.8978, 3.1368, 3.3214, 3.5595, 3.8028, 4.1094, 4.3563, 4.5945, 4.8009, 4.9581, 5.1112, 5.2655, 5.4308, 5.6469, 5.8871, 2.0755, 1.4709, 1.4835, 1.5553, 1.6623, 2.2123, 3.4115}, {0.4905, 0.4589, 0.4269, 0.3968, 0.3697, 0.3445, 0.3181, 0.3003, 0.2800, 0.2619, 0.2422, 0.2284, 0.2165, 0.2070, 0.2003, 0.1941, 0.1871, 0.1791, 0.1691, 0.1568, 0.3808, 0.5387, 0.5216, 0.4804, 0.4210, 0.2874, 0.1564}, {0.4978, 0.4763, 0.4543, 0.4336, 0.4146, 0.3970, 0.3782, 0.3653, 0.3485, 0.3350, 0.3202, 0.3091, 0.2986, 0.2898, 0.2832, 0.2769, 0.2718, 0.2683, 0.2651, 0.2625, 0.4078, 0.5526, 0.5547, 0.5413, 0.5231, 0.4322, 0.3363}, {0.4977, 0.4762, 0.4541, 0.4334, 0.4144, 0.3967, 0.3778, 0.3648, 0.3480, 0.3344, 0.3196, 0.3084, 0.2980, 0.2891, 0.2825, 0.2760, 0.2702, 0.2652, 0.2600, 0.2538, 0.3477, 0.4775, 0.4737, 0.4545, 0.4246, 0.3393, 0.2568}},
      {{1.4828, 1.3895, 1.2994, 1.2147, 1.1287, 1.0445, 0.9680, 0.9178, 0.8861, 0.8153, 0.7534, 0.7026, 0.6541, 0.6154, 0.5791, 0.5483, 0.5183, 0.4929, 0.4725, 0.4509, 0.4308, 0.4128, 0.3952, 0.3792, 0.3641, 0.3582, 0.3464}, {1.9051, 2.0480, 2.2044, 2.3720, 2.5688, 2.7936, 3.0380, 3.2646, 3.4561, 3.7619, 4.0677, 4.3413, 4.6239, 4.8637, 5.0967, 5.3014, 5.5141, 5.7104, 5.9030, 6.1723, 6.4457, 2.4514, 1.7116, 1.9032, 2.3021, 2.7195, 4.9815}, {0.5248, 0.4881, 0.4533, 0.4210, 0.3885, 0.3568, 0.3278, 0.3044, 0.2871, 0.2634, 0.2434, 0.2280, 0.2139, 0.2033, 0.1939, 0.1862, 0.1787, 0.1710, 0.1630, 0.1520, 0.1397, 0.3090, 0.4357, 0.3837, 0.2864, 0.2196, 0.0969}, {0.5181, 0.4936, 0.4701, 0.4479, 0.4258, 0.4040, 0.3837, 0.3671, 0.3520, 0.3353, 0.3208, 0.3090, 0.2975, 0.2882, 0.2794, 0.2721, 0.2651, 0.2598, 0.2561, 0.2526, 0.2503, 0.3630, 0.4984, 0.4698, 0.4227, 0.3811, 0.2799}, {0.5181, 0.4934, 0.4698, 0.4476, 0.4253, 0.4033, 0.3829, 0.3659, 0.3507, 0.3339, 0.3194, 0.3076, 0.2961, 0.2868, 0.2780, 0.2706, 0.2635, 0.2573, 0.2520, 0.2459, 0.2395, 0.3026, 0.4172, 0.3920, 0.3393, 0.2958, 0.2184}},
      {{1.3794, 1.2984, 1.2177, 1.1401, 1.0594, 0.9788, 0.8998, 0.8559, 0.8331, 0.7661, 0.7099, 0.6634, 0.6224, 0.5838, 0.5506, 0.5218, 0.4951, 0.4720, 0.4503, 0.4302, 0.4125, 0.3958, 0.3804, 0.3689, 0.3579, 0.3470, 0.3391, 0.3311, 0.3305}, {1.7953, 1.9363, 2.0962, 2.2749, 2.4911, 2.7485, 3.0611, 3.3180, 3.5212, 3.8705, 4.2226, 4.5688, 4.9219, 5.3012, 5.6666, 6.0213, 6.3919, 6.7492, 7.1036, 7.5059, 7.9539, 8.4004, 3.7248, 2.6497, 2.9279, 4.1581, 6.9235, 14.1071, 13.7325}, {0.5563, 0.5155, 0.4758, 0.4379, 0.3994, 0.3612, 0.3235, 0.2978, 0.2799, 0.2542, 0.2324, 0.2147, 0.1992, 0.1848, 0.1727, 0.1623, 0.1526, 0.1439, 0.1344, 0.1237, 0.1122, 0.1001, 0.1868, 0.2608, 0.2361, 0.1352, 0.0705, 0.0303, 0.0309}, {0.5373, 0.5098, 0.4832, 0.4575, 0.4313, 0.4050, 0.3786, 0.3600, 0.3447, 0.3270, 0.3119, 0.2986, 0.2864, 0.2745, 0.2640, 0.2546, 0.2457, 0.2378, 0.2308, 0.2243, 0.2187, 0.2148, 0.2692, 0.3668, 0.3542, 0.2986, 0.2430, 0.2007, 0.2035}, {0.5370, 0.5093, 0.4825, 0.4566, 0.4301, 0.4035, 0.3767, 0.3579, 0.3422, 0.3245, 0.3092, 0.2961, 0.2840, 0.2722, 0.2617, 0.2524, 0.2434, 0.2352, 0.2268, 0.2180, 0.2094, 0.2008, 0.2175, 0.2981, 0.2926, 0.2337, 0.1904, 0.1561, 0.1579}},
      {{1.2837, 1.2137, 1.1393, 1.0684, 0.9934, 0.9141, 0.8397, 0.7929, 0.7709, 0.7147, 0.6661, 0.6239, 0.5870, 0.5533, 0.5244, 0.4978, 0.4747, 0.4512, 0.4336, 0.4128, 0.3975, 0.3819, 0.3681, 0.3547, 0.3438, 0.3340, 0.3244}, {1.6685, 1.8090, 1.9796, 2.1703, 2.4098, 2.7218, 3.0931, 3.4401, 3.7157, 4.1105, 4.5324, 4.9757, 5.4408, 5.9491, 6.4682, 7.0431, 7.6540, 8.4214, 9.0855, 9.9174, 10.6082, 11.3890, 12.0339, 6.9862, 4.6894, 5.8851, 13.2613}, {0.5980, 0.5509, 0.5025, 0.4573, 0.4103, 0.3616, 0.3165, 0.2832, 0.2607, 0.2354, 0.2133, 0.1942, 0.1776, 0.1623, 0.1491, 0.1363, 0.1250, 0.1133, 0.1038, 0.0917, 0.0821, 0.0714, 0.0620, 0.0870, 0.1337, 0.1122, 0.0251}, {0.5628, 0.5305, 0.4983, 0.4683, 0.4372, 0.4047, 0.3737, 0.3498, 0.3289, 0.3113, 0.2958, 0.2820, 0.2695, 0.2574, 0.2465, 0.2359, 0.2260, 0.2155, 0.2075, 0.1981, 0.1911, 0.1840, 0.1779, 0.1867, 0.2593, 0.2512, 0.1941}, {0.5621, 0.5292, 0.4969, 0.4665, 0.4347, 0.4015, 0.3701, 0.3457, 0.3238, 0.3062, 0.2913, 0.2779, 0.2656, 0.2538, 0.2431, 0.2324, 0.2225, 0.2120, 0.2035, 0.1923, 0.1833, 0.1730, 0.1635, 0.1644, 0.2073, 0.2120, 0.1387}},
      {{1.2007, 1.1397, 1.0750, 1.0051, 0.9369, 0.8640, 0.7889, 0.7378, 0.7285, 0.6799, 0.6352, 0.6002, 0.5645, 0.5338, 0.5075, 0.4825, 0.4608, 0.4424, 0.4227, 0.4041, 0.3909, 0.3762, 0.3619, 0.3498, 0.3378, 0.3279, 0.3180, 0.3089}, {1.5325, 1.6778, 1.8504, 2.0726, 2.3364, 2.6866, 3.1531, 3.6321, 3.9293, 4.3700, 4.8704, 5.3472, 5.9432, 6.5775, 7.2440, 8.0225, 8.8632, 9.7471, 10.9042, 12.1832, 13.2053, 14.5004, 15.8927, 16.9543, 15.5916, 9.2167, 14.7595, 24.6247}, {0.6496, 0.5920, 0.5350, 0.4752, 0.4190, 0.3614, 0.3048, 0.2624, 0.2401, 0.2152, 0.1931, 0.1760, 0.1588, 0.1433, 0.1300, 0.1171, 0.1055, 0.0954, 0.0844, 0.0735, 0.0652, 0.0555, 0.0455, 0.0379, 0.0364, 0.0611, 0.0205, 0.0110}, {0.5946, 0.5618, 0.5274, 0.4913, 0.4569, 0.4203, 0.3816, 0.3484, 0.3138, 0.2960, 0.2799, 0.2673, 0.2542, 0.2427, 0.2324, 0.2224, 0.2133, 0.2054, 0.1967, 0.1885, 0.1825, 0.1758, 0.1693, 0.1637, 0.1584, 0.1829, 0.1786, 0.1650}, {0.5930, 0.5593, 0.5251, 0.4879, 0.4508, 0.4112, 0.3710, 0.3389, 0.3054, 0.2882, 0.2726, 0.2607, 0.2484, 0.2373, 0.2274, 0.2177, 0.2088, 0.2009, 0.1920, 0.1827, 0.1752, 0.1661, 0.1561, 0.1473, 0.1393, 0.1440, 0.1234, 0.1168}},
      {{1.1200, 1.0667, 1.0106, 0.9477, 0.8789, 0.8124, 0.7375, 0.6842, 0.6777, 0.6346, 0.5952, 0.5618, 0.5336, 0.5071, 0.4823, 0.4591, 0.4396, 0.4212, 0.4036, 0.3878, 0.3732, 0.3604, 0.3474, 0.3362, 0.3254, 0.3168, 0.3085, 0.3006}, {1.3733, 1.5158, 1.6942, 1.9391, 2.2567, 2.6464, 3.2296, 3.9058, 4.3546, 4.8960, 5.5040, 6.1321, 6.7804, 7.5235, 8.3915, 9.4092, 10.4762, 11.7269, 13.2091, 14.7715, 16.4129, 18.1942, 20.2737, 22.1012, 20.5313, 14.8994, 18.7698, 29.5311}, {0.7218, 0.6512, 0.5809, 0.5034, 0.4268, 0.3589, 0.2895, 0.2360, 0.2077, 0.1844, 0.1640, 0.1477, 0.1340, 0.1209, 0.1082, 0.0961, 0.0858, 0.0759, 0.0664, 0.0576, 0.0487, 0.0403, 0.0315, 0.0250, 0.0219, 0.0327, 0.0143, 0.0100}, {0.6392, 0.6048, 0.5695, 0.5309, 0.4879, 0.4442, 0.3954, 0.3562, 0.2964, 0.2791, 0.2635, 0.2504, 0.2394, 0.2289, 0.2190, 0.2096, 0.2015, 0.1937, 0.1861, 0.1793, 0.1730, 0.1673, 0.1613, 0.1562, 0.1513, 0.1471, 0.1635, 0.1530}, {0.6357, 0.6005, 0.5611, 0.5143, 0.4647, 0.4205, 0.3734, 0.3360, 0.2861, 0.2693, 0.2544, 0.2422, 0.2319, 0.2220, 0.2124, 0.2031, 0.1950, 0.1870, 0.1792, 0.1713, 0.1628, 0.1541, 0.1441, 0.1356, 0.1268, 0.1263, 0.1103, 0.1115}},
      {{1.0559, 1.0045, 0.9526, 0.8945, 0.8314, 0.7611, 0.7235, 0.6855, 0.6550, 0.6259, 0.6232, 0.5872, 0.5564, 0.5277, 0.5034, 0.4781, 0.4569, 0.4391, 0.4207, 0.4035, 0.3890, 0.3734, 0.3616, 0.3482, 0.3382, 0.3270, 0.3170, 0.3062, 0.2960}, {1.2225, 1.3516, 1.5089, 1.7231, 2.0359, 2.5372, 2.8908, 3.3360, 3.7893, 4.4728, 5.2957, 5.9895, 6.7023, 7.5080, 8.3349, 9.4033, 10.5182, 11.6599, 13.1172, 14.8040, 16.5230, 18.7131, 20.6047, 23.2733, 25.6296, 28.6505, 30.4520, 22.0582, 31.1532}, {0.8063, 0.7240, 0.6354, 0.5373, 0.4445, 0.3536, 0.3093, 0.2671, 0.2343, 0.1979, 0.1616, 0.1428, 0.1278, 0.1143, 0.1035, 0.0919, 0.0818, 0.0732, 0.0642, 0.0560, 0.0492, 0.0416, 0.0353, 0.0276, 0.0219, 0.0166, 0.0123, 0.0163, 0.0086}, {0.6907, 0.6520, 0.6134, 0.5661, 0.5118, 0.4514, 0.4214, 0.3926, 0.3702, 0.3475, 0.2845, 0.2673, 0.2532, 0.2403, 0.2296, 0.2184, 0.2090, 0.2011, 0.1928, 0.1850, 0.1785, 0.1727, 0.1664, 0.1604, 0.1558, 0.1507, 0.1463, 0.1497, 0.1495}, {0.6841, 0.6364, 0.5831, 0.5238, 0.4673, 0.4115, 0.3837, 0.3578, 0.3380, 0.3173, 0.2687, 0.2527, 0.2397, 0.2279, 0.2183, 0.2080, 0.1990, 0.1912, 0.1828, 0.1748, 0.1680, 0.1598, 0.1524, 0.1428, 0.1347, 0.1260, 0.1172, 0.1134, 0.1070}},
      {{0.9971, 0.9525, 0.9013, 0.8464, 0.7847, 0.7167, 0.6791, 0.6434, 0.6112, 0.5808, 0.5774, 0.5469, 0.5182, 0.4916, 0.4691, 0.4482, 0.4288, 0.4108, 0.3943, 0.3788, 0.3651, 0.3523, 0.3402, 0.3280, 0.3178, 0.3079, 0.3004, 0.2937}, {1.0711, 1.1752, 1.3188, 1.5241, 1.8504, 2.3664, 2.7434, 3.1900, 3.6965, 4.4433, 6.1604, 6.9829, 7.9161, 8.9756, 10.0753, 11.3301, 12.7707, 14.4287, 16.3106, 18.4799, 20.7713, 23.1589, 25.9383, 29.3033, 32.4612, 34.6541, 26.3523, 25.8082}, {0.9066, 0.7986, 0.6720, 0.5525, 0.4432, 0.3426, 0.2941, 0.2529, 0.2180, 0.1784, 0.1316, 0.1161, 0.1023, 0.0905, 0.0809, 0.0717, 0.0629, 0.0547, 0.0473, 0.0407, 0.0347, 0.0286, 0.0225, 0.0168, 0.0131, 0.0097, 0.0123, 0.0104}, {0.7562, 0.7162, 0.6668, 0.6066, 0.5373, 0.4650, 0.4292, 0.3975, 0.3709, 0.3451, 0.2732, 0.2565, 0.2416, 0.2283, 0.2250, 0.2154, 0.1982, 0.1993, 0.1819, 0.1864, 0.1684, 0.1627, 0.1571, 0.1515, 0.1468, 0.1424, 0.1536, 0.1345}, {0.7335, 0.6699, 0.5946, 0.5222, 0.4552, 0.3927, 0.3620, 0.3365, 0.3155, 0.2953, 0.2673, 0.2509, 0.2242, 0.2125, 0.2029, 0.1937, 0.1849, 0.1763, 0.1683, 0.1609, 0.1537, 0.1456, 0.1368, 0.1272, 0.1201, 0.1117, 0.1095, 0.0926}},
      {{0.9408, 0.8968, 0.8508, 0.7971, 0.7372, 0.6752, 0.6363, 0.5998, 0.5677, 0.5400, 0.5127, 0.4884, 0.4658, 0.4469, 0.4273, 0.4091, 0.3941, 0.3787, 0.3659, 0.3527, 0.3404, 0.3295, 0.3203, 0.3103, 0.3008, 0.2944, 0.2881}, {0.9167, 1.0127, 1.1583, 1.4014, 1.7862, 2.3287, 2.7866, 3.3303, 3.9367, 4.6451, 7.7966, 8.8503, 10.0418, 11.2404, 12.7498, 14.4797, 16.2102, 18.3289, 20.5325, 23.1038, 25.9211, 28.9688, 31.9576, 35.5240, 38.5181, 30.3877, 26.7480}, {0.9987, 0.8312, 0.6724, 0.5322, 0.4115, 0.3147, 0.2643, 0.2238, 0.1908, 0.1578, 0.0980, 0.0863, 0.0760, 0.0679, 0.0593, 0.0511, 0.0444, 0.0382, 0.0332, 0.0278, 0.0225, 0.0178, 0.0142, 0.0112, 0.0085, 0.0102, 0.0102}, {0.8440, 0.7890, 0.7181, 0.6310, 0.5389, 0.4594, 0.4145, 0.3821, 0.3534, 0.3290, 0.2467, 0.2327, 0.2202, 0.2102, 0.2000, 0.1908, 0.1834, 0.1758, 0.1696, 0.1634, 0.1578, 0.1526, 0.1481, 0.1432, 0.1387, 0.1353, 0.1375}, {0.7767, 0.6791, 0.5843, 0.4986, 0.4233, 0.3620, 0.3298, 0.3048, 0.2850, 0.2690, 0.2464, 0.2347, 0.2245, 0.2161, 0.2070, 0.1730, 0.1655, 0.1580, 0.1517, 0.1446, 0.1367, 0.1289, 0.1221, 0.1155, 0.1080, 0.1061, 0.0948}},
      {{0.8866, 0.8477, 0.8012, 0.7491, 0.6931, 0.6603, 0.6285, 0.5899, 0.5576, 0.5258, 0.4954, 0.4759, 0.4539, 0.4332, 0.4164, 0.3987, 0.3843, 0.3711, 0.3573, 0.3454, 0.3344, 0.3236, 0.3132, 0.3034, 0.2948, 0.2890, 0.2835}, {0.7653, 0.8718, 1.0627, 1.3640, 1.8112, 2.1192, 2.4701, 3.0916, 3.6627, 4.2853, 5.2810, 8.6838, 9.9779, 11.4337, 12.9488, 14.7614, 16.6421, 18.6687, 21.2150, 23.8594, 26.4680, 29.6536, 33.4581, 37.3515, 40.6568, 35.4621, 28.5297}, {1.0455, 0.8265, 0.6237, 0.4713, 0.3550, 0.3023, 0.2589, 0.2161, 0.1846, 0.1557, 0.1244, 0.0789, 0.0686, 0.0593, 0.0521, 0.0443, 0.0383, 0.0333, 0.0285, 0.0244, 0.0203, 0.0163, 0.0128, 0.0100, 0.0079, 0.0085, 0.0096}, {0.9558, 0.8681, 0.7524, 0.6319, 0.5229, 0.4737, 0.4322, 0.3880, 0.3435, 0.3252, 0.3001, 0.2360, 0.2219, 0.2096, 0.2000, 0.1903, 0.1825, 0.1756, 0.1685, 0.1627, 0.1574, 0.1520, 0.1466, 0.1413, 0.1363, 0.1330, 0.1241}, {0.7943, 0.6669, 0.5441, 0.4486, 0.3739, 0.3395, 0.3108, 0.2837, 0.2641, 0.2468, 0.2324, 0.2307, 0.2207, 0.2111, 0.2032, 0.1942, 0.1869, 0.1803, 0.1737, 0.1675, 0.1605, 0.1532, 0.1460, 0.1392, 0.1058, 0.1037, 0.1269}},
      {{0.8388, 0.7981, 0.7531, 0.7041, 0.6473, 0.6151, 0.5845, 0.5472, 0.5140, 0.4859, 0.4621, 0.4405, 0.4272, 0.4089, 0.3942, 0.3779, 0.3647, 0.3512, 0.3401, 0.3297, 0.3189, 0.3085, 0.2996, 0.2907, 0.2847, 0.2795}, {0.6584, 0.8122, 1.0605, 1.4344, 2.0771, 2.3620, 2.9239, 3.4635, 4.1812, 5.0304, 5.8371, 6.9694, 10.6828, 12.3430, 13.8651, 15.9115, 17.9813, 20.4490, 23.0285, 25.6299, 28.8723, 32.6911, 36.7904, 41.0825, 43.7272, 30.1994}, {1.0093, 0.7201, 0.5165, 0.3781, 0.2730, 0.2296, 0.1969, 0.1635, 0.1371, 0.1156, 0.0990, 0.0818, 0.0544, 0.0461, 0.0406, 0.0344, 0.0300, 0.0259, 0.0227, 0.0194, 0.0160, 0.0127, 0.0100, 0.0078, 0.0072, 0.0091}, {1.0613, 0.9063, 0.7449, 0.6007, 0.4639, 0.4320, 0.3745, 0.3487, 0.3155, 0.2899, 0.2698, 0.2553, 0.2154, 0.2033, 0.1946, 0.1847, 0.1770, 0.1694, 0.1634, 0.1583, 0.1528, 0.1468, 0.1411, 0.1353, 0.1315, 0.1245}, {0.7664, 0.5946, 0.4670, 0.3768, 0.3061, 0.2764, 0.2547, 0.2331, 0.2163, 0.2024, 0.1923, 0.1865, 0.2007, 0.1921, 0.1858, 0.1781, 0.1722, 0.1662, 0.1613, 0.1557, 0.1493, 0.1428, 0.1370, 0.1318, 0.1301, 0.1266}},
    }
};

// B coefficients
// Table in Appendix A of Hurley et al. 2000
// Key to map is n (B(n)).  Map element is unordered_map of term coefficient values.
const std::map<int, COMPASUnorderedMap<AB_TCoeff, double>> B_COEFF = {
    { 1, {{ALPHA,  3.970000E-1}, {BETA,  2.882600E-1}, {GAMMA,  5.293000E-1}, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    { 2, {{ALPHA,  0.000000E0 }, {BETA,  0.000000E0 }, {GAMMA,  0.000000E0 }, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    { 3, {{ALPHA,  0.000000E0 }, {BETA,  0.000000E0 }, {GAMMA,  0.000000E0 }, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    { 4, {{ALPHA,  9.960283E-1}, {BETA,  8.164393E-1}, {GAMMA,  2.383830E0 }, {ETA,  2.223436E0 }, {MU,  8.638115E-1}}},
    { 5, {{ALPHA,  2.561062E-1}, {BETA,  7.072646E-2}, {GAMMA, -5.444596E-2}, {ETA, -5.798167E-2}, {MU, -1.349129E-2}}},
    { 6, {{ALPHA,  1.157338E0 }, {BETA,  1.467883E0 }, {GAMMA,  4.299661E0 }, {ETA,  3.130500E0 }, {MU,  6.992080E-1}}},
    { 7, {{ALPHA,  4.022765E-1}, {BETA,  3.050010E-1}, {GAMMA,  9.962137E-1}, {ETA,  7.914079E-1}, {MU,  1.728098E-1}}},
    { 8, {{ALPHA,  0.000000E0 }, {BETA,  0.000000E0 }, {GAMMA,  0.000000E0 }, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    { 9, {{ALPHA,  2.751631E3 }, {BETA,  3.557098E2 }, {GAMMA,  0.000000E0 }, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    {10, {{ALPHA, -3.820831E-2}, {BETA,  5.872664E-2}, {GAMMA,  0.000000E0 }, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},

    {11, {{ALPHA,  1.071738E2 }, {BETA, -8.970339E1 }, {GAMMA, -3.949739E1 }, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    {12, {{ALPHA,  7.348793E2 }, {BETA, -1.531020E2 }, {GAMMA, -3.793700E1 }, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    {13, {{ALPHA,  9.219293E0 }, {BETA, -2.005865E0 }, {GAMMA, -5.561309E-1}, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    {14, {{ALPHA,  2.917412E0 }, {BETA,  1.575290E0 }, {GAMMA,  5.751814E-1}, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    {15, {{ALPHA,  3.629118E0 }, {BETA, -9.112722E-1}, {GAMMA,  1.042291E0 }, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    {16, {{ALPHA,  4.916389E0 }, {BETA,  2.862149E0 }, {GAMMA,  7.844850E-1}, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    {17, {{ALPHA,  0.000000E0 }, {BETA,  0.000000E0 }, {GAMMA,  0.000000E0 }, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    {18, {{ALPHA,  5.496045E1 }, {BETA, -1.289968E1 }, {GAMMA,  6.385758E0 }, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    {19, {{ALPHA,  1.832694E0 }, {BETA, -5.766608E-2}, {GAMMA,  5.696128E-2}, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    {20, {{ALPHA,  1.211104E2 }, {BETA,  0.000000E0 }, {GAMMA,  0.000000E0 }, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},

    {21, {{ALPHA,  2.214088E2 }, {BETA,  2.187113E2 }, {GAMMA,  1.170177E1 }, {ETA, -2.635340E1 }, {MU,  0.000000E0 }}},
    {22, {{ALPHA,  2.063983E0 }, {BETA,  7.363827E-1}, {GAMMA,  2.654323E-1}, {ETA, -6.140719E-2}, {MU,  0.000000E0 }}},
    {23, {{ALPHA,  2.003160E0 }, {BETA,  9.388871E-1}, {GAMMA,  9.656450E-1}, {ETA,  2.362266E-1}, {MU,  0.000000E0 }}},
    {24, {{ALPHA,  1.609901E1 }, {BETA,  7.391573E0 }, {GAMMA,  2.277010E1 }, {ETA,  8.334227E0 }, {MU,  0.000000E0 }}},
    {25, {{ALPHA,  1.747500E-1}, {BETA,  6.271202E-2}, {GAMMA, -2.324229E-2}, {ETA, -1.844559E-2}, {MU,  0.000000E0 }}},
    {26, {{ALPHA,  0.000000E0 }, {BETA,  0.000000E0 }, {GAMMA,  0.000000E0 }, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    {27, {{ALPHA,  2.752869E0 }, {BETA,  2.729201E-2}, {GAMMA,  4.996927E-1}, {ETA,  2.496551E-1}, {MU,  0.000000E0 }}},
    {28, {{ALPHA,  3.518506E0 }, {BETA,  1.112440E0 }, {GAMMA, -4.556216E-1}, {ETA, -2.179426E-1}, {MU,  0.000000E0 }}},
    {29, {{ALPHA,  1.626062E2 }, {BETA, -1.168838E1 }, {GAMMA, -5.498343E0 }, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    {30, {{ALPHA,  3.336833E-1}, {BETA, -1.458043E-1}, {GAMMA, -2.011751E-2}, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},

    {31, {{ALPHA,  7.425137E1 }, {BETA,  1.790236E1 }, {GAMMA,  3.033910E1 }, {ETA,  1.018259E1 }, {MU,  0.000000E0 }}},
    {32, {{ALPHA,  9.268325E2 }, {BETA, -9.739859E1 }, {GAMMA, -7.702152E1 }, {ETA, -3.158268E1 }, {MU,  0.000000E0 }}},
    {33, {{ALPHA,  2.474401E0 }, {BETA,  3.892972E-1}, {GAMMA,  0.000000E0 }, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    {34, {{ALPHA,  1.127018E1 }, {BETA,  1.622158E0 }, {GAMMA, -1.443664E0 }, {ETA, -9.474699E-1}, {MU,  0.000000E0 }}},
    {35, {{ALPHA,  0.000000E0 }, {BETA,  0.000000E0 }, {GAMMA,  0.000000E0 }, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    {36, {{ALPHA,  1.445216E-1}, {BETA, -6.180219E-2}, {GAMMA,  3.093878E-2}, {ETA,  1.567090E-2}, {MU,  0.000000E0 }}},
    {37, {{ALPHA,  1.304129E0 }, {BETA,  1.395919E-1}, {GAMMA,  4.142455E-3}, {ETA, -9.732503E-3}, {MU,  0.000000E0 }}},
    {38, {{ALPHA,  5.114149E-1}, {BETA, -1.160850E-2}, {GAMMA,  0.000000E0 }, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    {39, {{ALPHA,  1.314955E2 }, {BETA,  2.009258E1 }, {GAMMA, -5.143082E-1}, {ETA, -1.379140E0 }, {MU,  0.000000E0 }}},
    {40, {{ALPHA,  1.823973E1 }, {BETA, -3.074559E0 }, {GAMMA, -4.307878E0 }, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},

    {41, {{ALPHA,  2.327037E0 }, {BETA,  2.403445E0 }, {GAMMA,  1.208407E0 }, {ETA,  2.087263E-1}, {MU,  0.000000E0 }}},
    {42, {{ALPHA,  1.997378E0 }, {BETA, -8.126205E-1}, {GAMMA,  0.000000E0 }, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    {43, {{ALPHA,  1.079113E-1}, {BETA,  1.762409E-2}, {GAMMA,  1.096601E-2}, {ETA,  3.058818E-3}, {MU,  0.000000E0 }}},
    {44, {{ALPHA,  2.327409E0 }, {BETA,  6.901582E-1}, {GAMMA, -2.158431E-1}, {ETA, -1.084117E-1}, {MU,  0.000000E0 }}},
    {45, {{ALPHA,  0.000000E0 }, {BETA,  0.000000E0 }, {GAMMA,  0.000000E0 }, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    {46, {{ALPHA,  2.214315E0 }, {BETA, -1.975747E0 }, {GAMMA,  0.000000E0 }, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    {47, {{ALPHA,  0.000000E0 }, {BETA,  0.000000E0 }, {GAMMA,  0.000000E0 }, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    {48, {{ALPHA,  5.072525E0 }, {BETA,  1.146189E1 }, {GAMMA,  6.961724E0 }, {ETA,  1.316965E0 }, {MU,  0.000000E0 }}},
    {49, {{ALPHA,  5.139740E0 }, {BETA,  0.000000E0 }, {GAMMA,  0.000000E0 }, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    {50, {{ALPHA,  0.000000E0 }, {BETA,  0.000000E0 }, {GAMMA,  0.000000E0 }, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},

    {51, {{ALPHA,  1.125124E0 }, {BETA,  1.306486E0 }, {GAMMA,  3.622359E0 }, {ETA,  2.601976E0 }, {MU,  3.031270E-1}}},
    {52, {{ALPHA,  3.349489E-1}, {BETA,  4.531269E-3}, {GAMMA,  1.131793E-1}, {ETA,  2.300156E-1}, {MU,  7.632745E-2}}},
    {53, {{ALPHA,  1.467794E0 }, {BETA,  2.798142E0 }, {GAMMA,  9.455580E0 }, {ETA,  8.963904E0 }, {MU,  3.339719E0 }}},
    {54, {{ALPHA,  4.658512E-1}, {BETA,  2.597451E-1}, {GAMMA,  9.048179E-1}, {ETA,  7.394505E-1}, {MU,  1.607092E-1}}},
    {55, {{ALPHA,  1.042200E0 }, {BETA,  1.315600E-1}, {GAMMA,  4.500000E-2}, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    {56, {{ALPHA,  1.110866E0 }, {BETA,  9.623856E-1}, {GAMMA,  2.735487E0 }, {ETA,  2.445602E0 }, {MU,  8.826352E-1}}},
    {57, {{ALPHA, -1.584333E-1}, {BETA, -1.728865E-1}, {GAMMA, -4.461431E-1}, {ETA, -3.925259E-1}, {MU, -1.276203E-1}}}
};

#undef ALPHA
#undef BETA
#undef GAMMA
#undef ETA
#undef MU


// C coefficients
// Key to map is n (C(n)).  Map element is unordered_map of term coefficient values.
const std::unordered_map<int, double> C_COEFF = {{1, -8.672073E-2}, {2, 9.301992E0}, {3, 4.637345E0}};


// CDF from Table 7 in Dufton et al 2013 https://arxiv.org/abs/1212.2424
// There is an assumption in the code that this function is monitonically increasing - it is now
// and should remain so if the map is modified.
const std::map<double, double> BStarRotationalVelocityCDFTable = {
    {000.0, 0.000}, {020.0, 0.046}, {040.0, 0.094}, {060.0, 0.144}, {080.0, 0.192}, {100.0, 0.239},
    {120.0, 0.253}, {140.0, 0.270}, {160.0, 0.288}, {180.0, 0.322}, {200.0, 0.377}, {220.0, 0.435},
    {240.0, 0.492}, {260.0, 0.548}, {280.0, 0.609}, {300.0, 0.674}, {320.0, 0.739}, {340.0, 0.796},
    {360.0, 0.841}, {380.0, 0.879}, {400.0, 0.912}, {420.0, 0.938}, {440.0, 0.956}, {460.0, 0.971},
    {480.0, 0.983}, {500.0, 0.990}, {520.0, 0.993}, {540.0, 0.995}, {560.0, 0.996}, {580.0, 0.997}
};


// These neutron star (NS) equations-of-state (EOS) are
// taken from the review Ozel & Freire 2016,
// Masses, Radii, and Equation of State of Neutron Stars,
// Annual Reviews of Astronomy and Astrophysics,
// https://arxiv.org/abs/1603.02698, downloaded from
// their website http://xtreme.as.arizona.edu/NeutronStars/

// For now we choose one example EOS ARP3 from
// Akmal et al 1998 https://arxiv.org/abs/nucl-th/9804027
const std::map<double, double> ARP3MassRadiusRelation = {
    {0.184 , 16.518}, {0.188 , 16.292}, {0.192 , 16.067}, {0.195 , 15.857}, {0.199 , 15.658}, {0.203 , 15.46 }, {0.207 , 15.277}, {0.212, 15.102}, {0.216, 14.933},
    {0.221 , 14.774}, {0.225 , 14.619}, {0.23  , 14.473}, {0.235 , 14.334}, {0.24  , 14.199}, {0.245 , 14.073}, {0.251 , 13.951}, {0.256, 13.834}, {0.262, 13.725},
    {0.268 , 13.618}, {0.273 , 13.52 }, {0.28  , 13.423}, {0.286 , 13.332}, {0.292 , 13.245}, {0.299 , 13.162}, {0.306 , 13.084}, {0.313, 13.009}, {0.32 , 12.94 },
    {0.327 , 12.871}, {0.335 , 12.806}, {0.342 , 12.747}, {0.35  , 12.691}, {0.358 , 12.638}, {0.366 , 12.586}, {0.374 , 12.538}, {0.383, 12.493}, {0.391, 12.451},
    {0.4   , 12.409}, {0.409 , 12.371}, {0.418 , 12.336}, {0.427 , 12.302}, {0.438 , 12.269}, {0.448 , 12.239}, {0.458 , 12.211}, {0.468, 12.184}, {0.479, 12.16 },
    {0.49  , 12.136}, {0.501 , 12.116}, {0.512 , 12.096}, {0.524 , 12.078}, {0.535 , 12.061}, {0.547 , 12.046}, {0.559 , 12.031}, {0.572, 12.018}, {0.585, 12.007},
    {0.598 , 11.997}, {0.611 , 11.987}, {0.625 , 11.979}, {0.638 , 11.972}, {0.652 , 11.966}, {0.666 , 11.96 }, {0.681 , 11.955}, {0.695, 11.952}, {0.71 , 11.949},
    {0.725 , 11.947}, {0.74  , 11.946}, {0.756 , 11.945}, {0.772 , 11.945}, {0.788 , 11.945}, {0.804 , 11.946}, {0.82  , 11.947}, {0.837, 11.949}, {0.854, 11.952},
    {0.871 , 11.955}, {0.888 , 11.957}, {0.906 , 11.961}, {0.923 , 11.964}, {0.941 , 11.968}, {0.959 , 11.972}, {0.977 , 11.977}, {0.995, 11.981}, {1.014, 11.985},
    {1.032 , 11.99 }, {1.05  , 11.994}, {1.069 , 11.999}, {1.088 , 12.004}, {1.107 , 12.009}, {1.126 , 12.013}, {1.145 , 12.018}, {1.164, 12.022}, {1.184, 12.027},
    {1.203 , 12.031}, {1.222 , 12.035}, {1.242 , 12.039}, {1.261 , 12.043}, {1.281 , 12.047}, {1.3   , 12.05 }, {1.32  , 12.053}, {1.339, 12.056}, {1.358, 12.058},
    {1.378 , 12.061}, {1.397 , 12.063}, {1.416 , 12.064}, {1.436 , 12.066}, {1.455 , 12.067}, {1.474 , 12.068}, {1.493 , 12.068}, {1.512, 12.068}, {1.531, 12.068},
    {1.549 , 12.067}, {1.568 , 12.066}, {1.586 , 12.065}, {1.604 , 12.063}, {1.623 , 12.06 }, {1.64  , 12.058}, {1.658 , 12.055}, {1.676, 12.052}, {1.693, 12.048},
    {1.71  , 12.044}, {1.727 , 12.039}, {1.744 , 12.034}, {1.761 , 12.029}, {1.777 , 12.024}, {1.793 , 12.017}, {1.809 , 12.011}, {1.825, 12.004}, {1.84 , 11.997},
    {1.856 , 11.989}, {1.871 , 11.981}, {1.886 , 11.973}, {1.9   , 11.965}, {1.915 , 11.956}, {1.929 , 11.946}, {1.943 , 11.937}, {1.956, 11.927}, {1.969, 11.916},
    {1.982 , 11.906}, {1.995 , 11.895}, {2.008 , 11.884}, {2.02  , 11.827}, {2.032 , 11.86 }, {2.044 , 11.848}, {2.056 , 11.836}, {2.067, 11.823}, {2.078, 11.81 },
    {2.089 , 11.797}, {2.099 , 11.784}, {2.109 , 11.77 }, {2.119 , 11.756}, {2.129 , 11.742}, {2.139 , 11.727}, {2.148 , 11.713}, {2.157, 11.698}, {2.166, 11.683},
    {2.174 , 11.668}, {2.182 , 11.652}, {2.19  , 11.637}, {2.198 , 11.621}, {2.206 , 11.605}, {2.213 , 11.589}, {2.221 , 11.573}, {2.227, 11.556}, {2.234, 11.54 },
    {2.241 , 11.523}, {2.247 , 11.506}, {2.253 , 11.49 }, {2.259 , 11.473}, {2.264 , 11.456}, {2.27  , 11.438}, {2.275 , 11.421}, {2.28 , 11.404}, {2.285, 11.386},
    {2.29  , 11.369}, {2.294 , 11.351}, {2.299 , 11.333}, {2.303 , 11.316}, {2.307 , 11.299}, {2.31  , 11.281}, {2.314 , 11.263}, {2.317, 11.245}, {2.321, 11.227},
    {2.324 , 11.209}, {2.327 , 11.191}, {2.33  , 11.173}, {2.332 , 11.155}, {2.335 , 11.136}, {2.337 , 11.119}, {2.339 , 11.101}, {2.342, 11.083}, {2.344, 11.065},
    {2.345 , 11.046}, {2.347 , 11.028}, {2.349 , 11.01 }, {2.35  , 10.992}, {2.352 , 10.974}, {2.353 , 10.956}, {2.354 , 10.938}, {2.355, 10.92 }, {2.356, 10.902},
    {2.3571, 10.885}, {2.3572, 10.866}, {2.3581, 10.849}, {2.3582, 10.831}, {2.3591, 10.813}, {2.3592, 10.795}, {2.3593, 10.777}, {2.361, 10.76 }, {2.362, 10.742}
};


// vector NANJING_MASSES
// Mass / Msun of stellar models computed by Xu & Li (2010), with additional unpublished 50 and 100 Msun models
const DBL_VECTOR NANJING_MASSES = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0, 20.0, 50.0, 100.0 };

// vector NANJING_MASSES_MIDPOINTS
// Mass / Msun bin edges of Xu & Li (2010) lambda prescription as implemented in StarTrack. These are the midpoints between the masses in NANJING_MASSES
const DBL_VECTOR NANJING_MASSES_MIDPOINTS = { 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 11.0, 13.0, 15.0, 18.0, 35.0, 75.0 };

// Coefficients for calculating binding and recombination energy as described in Loveridge et al. 2011
// Electronic tables, program and further information in: http://astro.ru.nl/~sluys/index.php?title=BE

// struct LoveridgeCoefficients
// m, r, alpha(m,r) used in Loveridge et al. 2011, eq 5
struct LoveridgeCoefficients {
    int    m;
    int    r;
    double alpha_mr;
};


// enum class LOVERIDGE_METALLICITY
// Symbolic names for metallicities described in Loveridge et al., 2011
// These are used as indices into the loveridgeCoefficients multi-dimensional vector (described below)
// This is a bit of a hack until I figure out how to (elegantly) iterate over enum classes...
enum class LOVERIDGE_METALLICITY: int { Z000010, Z000100, Z001000, Z001500, Z002000, Z003000, COUNT };
const std::vector<std::tuple<LOVERIDGE_METALLICITY, double>> LOVERIDGE_METALLICITY_VALUE = {
    { LOVERIDGE_METALLICITY::Z000010, 0.00010 },
    { LOVERIDGE_METALLICITY::Z000100, 0.00100 },
    { LOVERIDGE_METALLICITY::Z001000, 0.01000 },
    { LOVERIDGE_METALLICITY::Z001500, 0.01500 },
    { LOVERIDGE_METALLICITY::Z002000, 0.02000 },
    { LOVERIDGE_METALLICITY::Z003000, 0.03000 }
};


// enum class LOVERIDGE_GROUP
// Symbolic names for GB groups described in Loveridge et al., 2011
// These are used as indices into the loveridgeCoefficients multi-dimensional vector (described below)
enum class LOVERIDGE_GROUP: int { LMR1, LMR2, LMA, HM, RECOM };
const COMPASUnorderedMap<LOVERIDGE_GROUP, std::string> LOVERIDGE_GROUP_LABEL = {
    { LOVERIDGE_GROUP::LMR1,  "Low mass early Red Giant Branch (RGB) (before dredge-up)" },
    { LOVERIDGE_GROUP::LMR2,  "Low mass late Red Giant Branch (RGB) (after dredge-up)" },
    { LOVERIDGE_GROUP::LMA,   "Low mass Asymptotic Giant Branch (AGB)" },
    { LOVERIDGE_GROUP::HM,    "High mass" },
    { LOVERIDGE_GROUP::RECOM, "Recombination energy" }
};


// vector LOVERIDGE_LM_HM_CUTOFFS
// Values for the division between High Mass and Low Mass, indexed by metallicity (LOVERIDGE_METALLICITY)
// Values are given in Msol
// From Loveridge et al. 2011, table 1
const DBL_VECTOR LOVERIDGE_LM_HM_CUTOFFS = { 11.7000, 11.7000, 10.2000, 11.7000, 11.7000, 13.4000 };


// vector LOVERIDGE_LM1_LM2_CUTOFFS
// Coefficients for the division between Low Mass RGB 1 and Low Mass RGB 2, indexed by metallicity (LOVERIDGE_METALLICITY)
// From Loveridge et al. 2011, table 2 & eq 4
const std::vector<DBL_VECTOR> LOVERIDGE_LM1_LM2_CUTOFFS = {
    { 6.047230E-01,  1.397240E+00, -8.249630E-02,  1.141430E+00,  0.000000E+00 },       // Metallicity Z000010 (0.00010)
    { 4.566000E-01,  1.186600E+00,  2.390880E+00, -3.054040E+00,  1.404900E+00 },       // Metallicity Z000100 (0.00100)
    { 2.821740E-01,  1.149380E+00,  1.884450E+00, -1.082300E+00,  0.000000E+00 },       // Metallicity Z001000 (0.01000)
    { 2.518180E-01,  1.210490E+00,  1.634900E+00, -8.369100E-01,  0.000000E+00 },       // Metallicity Z001500 (0.01500)
    { 2.406370E-01,  1.089220E+00,  1.953180E+00, -1.032180E+00,  0.000000E+00 },       // Metallicity Z002000 (0.02000)
    { 2.348880E-01,  8.972940E-01,  2.519950E+00, -1.414110E+00,  0.000000E+00 }        // Metallicity Z003000 (0.03000)
};


// vector LOVERIDGE_COEFFICIENTS
// Multi-dimensional vector of struct LoveridgeCoefficients, indexed by
// metallicity (LOVERIDGE_METALLICITY) and evolutionary stage (LOVERIDGE_GROUP)
// This vector records the coefficients (m, r, alpha(m,r) used in Loveridge et al. 2001, eq 5
// Using a vector indexed by metallicity and evolutionary stage because it's faster than a map, and
// the lambdas could be calculated at every timestep.
const std::vector<std::vector<std::vector<LoveridgeCoefficients>>> LOVERIDGE_COEFFICIENTS = {
    {                                                                                   // Metallicity Z000010 (0.00010)
        {                                                                               // LMR1 (Z000010)
            { 0,     0,     1.49884369566236408389E+01},
            { 0,     1,     3.55674019216888570583E+00},
            { 0,     2,    -1.50579325323499482181E+01},
            { 0,     3,     2.74507278637946647848E+01},
            { 0,     4,    -2.40420132204742671433E+01},
            { 0,     5,     6.08559902401751795509E+00},

            { 1,     0,     4.77689517615753889146E+00},
            { 1,     1,    -3.52448257631879471319E+01},
            { 1,     2,     1.26181509166749165729E+02},
            { 1,     3,    -2.05760448139415075275E+02},
            { 1,     4,     1.66826460262414656199E+02},
            { 1,     5,    -4.13681688626000720888E+01},

            { 2,     0,    -4.03038048915965774199E+00},
            { 2,     1,     5.94984473818373302834E+01},
            { 2,     2,    -3.15891070807037635859E+02},
            { 2,     3,     5.63465810280387586317E+02},
            { 2,     4,    -4.56401018552895436642E+02},
            { 2,     5,     1.12357256780075772440E+02},

            { 3,     0,    -2.22258127619659440199E+00},
            { 3,     1,     1.21818087660567726971E+02},
            { 3,     2,    -6.06793401690043339158E+01},
            { 3,     3,    -2.73750098046708558286E+02},
            { 3,     4,     3.80968083978251002009E+02},
            { 3,     5,    -1.07210865421446229107E+02},

            { 4,     0,    -4.95448396548353997559E+01},
            { 4,     1,    -1.05676079281290000722E+02},
            { 4,     2,     3.74254532751612941865E+02},
            { 4,     3,    -2.84755814237885886087E+02},
            { 4,     4,     5.32060692031168436245E+00},
            { 4,     5,     1.94841031059088862776E+01},

            { 5,     0,     5.85417149247924371025E+01},
            { 5,     1,    -6.46713025038344397899E+01},
            { 5,     2,    -8.04514035300949643670E+01},
            { 5,     3,     1.54013846765123219029E+02},
            { 5,     4,    -6.62783052076742649206E+01},
            { 5,     5,     9.83910595056972248074E+00}
        },
        {                                                                               // LMR2 (Z000010)
            { 0,     0,     2.10206064943832124925E+01},
            { 0,     1,    -2.39940628010456791230E+01},
            { 0,     2,     3.67437434259622861532E+01},
            { 0,     3,    -2.87504026348741277275E+01},
            { 0,     4,     1.10696952815601967757E+01},
            { 0,     5,    -1.67354101724841819454E+00},

            { 1,     0,     6.24726695402092602194E+01},
            { 1,     1,    -2.25859701401090774198E+02},
            { 1,     2,     3.25693445380178616233E+02},
            { 1,     3,    -2.28906270354160255920E+02},
            { 1,     4,     7.82835291167177160787E+01},
            { 1,     5,    -1.04409269263635096081E+01},

            { 2,     0,     1.68774936141528343114E+02},
            { 2,     1,    -4.70922534725343496120E+02},
            { 2,     2,     5.20150477052292671942E+02},
            { 2,     3,    -2.82942436111233064366E+02},
            { 2,     4,     7.54607477257930696624E+01},
            { 2,     5,    -7.80062541052705249456E+00},

            { 3,     0,     1.26323501968766254322E+03},
            { 3,     1,    -5.43724065618109580100E+03},
            { 3,     2,     9.47031538171058127773E+03},
            { 3,     3,    -8.20344328990647773026E+03},
            { 3,     4,     3.48253888526251739677E+03},
            { 3,     5,    -5.75361752664876235031E+02},

            { 4,     0,     1.45320316532362594444E+04},
            { 4,     1,    -6.10692503818239565589E+04},
            { 4,     2,     9.45752483181984280236E+04},
            { 4,     3,    -6.92033750093292765087E+04},
            { 4,     4,     2.43234260768021413242E+04},
            { 4,     5,    -3.32540856427475091550E+03},

            { 5,     0,    -7.83727239733487567719E+03},
            { 5,     1,     6.87101874631883547409E+04},
            { 5,     2,    -1.42788737041162559763E+05},
            { 5,     3,     1.25369407779255416244E+05},
            { 5,     4,    -5.05985607497797464021E+04},
            { 5,     5,     7.77505329663658358186E+03}
        },
        {                                                                               // LMA (Z000010)
            { 0,     0,     1.77846204423370872973E+04},
            { 0,     1,    -1.11008631122171675088E+05},
            { 0,     2,     3.07385212080689030699E+05},
            { 0,     3,    -4.97253519789625774138E+05},
            { 0,     4,     5.20899845929651521146E+05},
            { 0,     5,    -3.69562230436008889228E+05},
            { 0,     6,     1.79995285036839370150E+05},
            { 0,     7,    -5.94766776453754428076E+04},
            { 0,     8,     1.27704226161695205519E+04},
            { 0,     9,    -1.60998313917297696207E+03},
            { 0,    10,     9.05540938508377593053E+01},

            { 1,     0,     1.27448576992469941615E+05},
            { 1,     1,    -8.26126519162579439580E+05},
            { 1,     2,     2.40616127883097669110E+06},
            { 1,     3,    -4.13147857106406055391E+06},
            { 1,     4,     4.61530595326700154692E+06},
            { 1,     5,    -3.49471896526116924360E+06},
            { 1,     6,     1.81238210394326946698E+06},
            { 1,     7,    -6.34629941238884348422E+05},
            { 1,     8,     1.43452093324876565021E+05},
            { 1,     9,    -1.88911385051759207272E+04},
            { 1,    10,     1.10039680760221062883E+03},

            { 2,     0,    -1.11280780545824207366E+06},
            { 2,     1,     6.52773804973363596946E+06},
            { 2,     2,    -1.68385778483916968107E+07},
            { 2,     3,     2.51369743430132977664E+07},
            { 2,     4,    -2.40291083278050050139E+07},
            { 2,     5,     1.53511958359827548265E+07},
            { 2,     6,    -6.62599811194045469165E+06},
            { 2,     7,     1.90266653405042551458E+06},
            { 2,     8,    -3.46290151645659178030E+05},
            { 2,     9,     3.57968178594517958118E+04},
            { 2,    10,    -1.57403299302352661471E+03},

            { 3,     0,     3.61468220664994791150E+06},
            { 3,     1,    -1.42838660574260130525E+07},
            { 3,     2,     9.89719916261141002178E+06},
            { 3,     3,     4.17630757517836764455E+07},
            { 3,     4,    -1.13186791305614486337E+08},
            { 3,     5,     1.34723130784819722176E+08},
            { 3,     6,    -9.39780759352868050337E+07},
            { 3,     7,     4.08028988015334084630E+07},
            { 3,     8,    -1.08834827876409199089E+07},
            { 3,     9,     1.63703844878768618219E+06},
            { 3,    10,    -1.06502153903636033647E+05},

            { 4,     0,    -2.18383920460389368236E+07},
            { 4,     1,     1.01377685264262214303E+08},
            { 4,     2,    -1.24736986550756111741E+08},
            { 4,     3,    -1.38097782211961090565E+08},
            { 4,     4,     5.82118970734395384789E+08},
            { 4,     5,    -7.72188668410225749016E+08},
            { 4,     6,     5.69788365736976385117E+08},
            { 4,     7,    -2.56651440166880398989E+08},
            { 4,     8,     7.03184175257203429937E+07},
            { 4,     9,    -1.07993168413460906595E+07},
            { 4,    10,     7.14464107997456681915E+05},

            { 5,     0,     6.22013266083969771862E+07},
            { 5,     1,    -3.70892432569035887718E+08},
            { 5,     2,     7.32722076455112814903E+08},
            { 5,     3,    -3.41863748162672758102E+08},
            { 5,     4,    -8.72008743590860724449E+08},
            { 5,     5,     1.70439955004952502251E+09},
            { 5,     6,    -1.44412099096650505066E+09},
            { 5,     7,     7.01467443390604257584E+08},
            { 5,     8,    -2.01846242185972064734E+08},
            { 5,     9,     3.21032091475058645010E+07},
            { 5,    10,    -2.18098983308966364712E+06},

            { 6,     0,     8.32790288301304075867E+06},
            { 6,     1,     1.67320728489836782217E+08},
            { 6,     2,    -4.90329729719172537327E+08},
            { 6,     3,    -3.57015713222805708647E+07},
            { 6,     4,     1.73606751974490427971E+09},
            { 6,     5,    -2.95511773960293579102E+09},
            { 6,     6,     2.49140620948586368561E+09},
            { 6,     7,    -1.22675662513774180412E+09},
            { 6,     8,     3.58779851358991682529E+08},
            { 6,     9,    -5.79478609330464825034E+07},
            { 6,    10,     3.99136670739177288488E+06},

            { 7,     0,    -3.09375949266583919525E+08},
            { 7,     1,     1.30927392519327545166E+09},
            { 7,     2,    -2.71972201258040809631E+09},
            { 7,     3,     4.17056345501154565811E+09},
            { 7,     4,    -5.31121138472141742706E+09},
            { 7,     5,     5.18211883778930091858E+09},
            { 7,     6,    -3.52386318383868503571E+09},
            { 7,     7,     1.57889703104470300674E+09},
            { 7,     8,    -4.41837538483527064323E+08},
            { 7,     9,     6.98762237560149580240E+07},
            { 7,    10,    -4.76660235680679418147E+06},

            { 8,     0,     5.32124163954892635345E+08},
            { 8,     1,    -2.58332422589960527420E+09},
            { 8,     2,     5.78740993511894130707E+09},
            { 8,     3,    -8.18639627587050056458E+09},
            { 8,     4,     8.33603336255734443665E+09},
            { 8,     5,    -6.37392318348361968994E+09},
            { 8,     6,     3.59443787565530967712E+09},
            { 8,     7,    -1.42294472536891078949E+09},
            { 8,     8,     3.68482395798513412476E+08},
            { 8,     9,    -5.55465651649229973555E+07},
            { 8,    10,     3.67803603909488115460E+06},

            { 9,     0,    -3.67132831968976199627E+08},
            { 9,     1,     1.87202460446496915817E+09},
            { 9,     2,    -4.32871663092040729523E+09},
            { 9,     3,     6.09032624990053653717E+09},
            { 9,     4,    -5.87735082626143074036E+09},
            { 9,     5,     4.10221144979333591461E+09},
            { 9,     6,    -2.08821944738185024261E+09},
            { 9,     7,     7.54238944109313964844E+08},
            { 9,     8,    -1.81710771632690608501E+08},
            { 9,     9,     2.59757891273681670427E+07},
            { 9,    10,    -1.65624886172299087048E+06},

            {10,     0,     9.30481611941920518875E+07},
            {10,     1,    -4.87096365846687316895E+08},
            {10,     2,     1.14762701275433015823E+09},
            {10,     3,    -1.62114284406450057030E+09},
            {10,     4,     1.53847292382848095894E+09},
            {10,     5,    -1.03373922282777369022E+09},
            {10,     6,     4.99196402871803343296E+08},
            {10,     7,    -1.70266350400250434875E+08},
            {10,     8,     3.88883041814338341355E+07},
            {10,     9,    -5.31391571130599640310E+06},
            {10,    10,     3.26864794701321807224E+05}
        },
        {                                                                               // HM (Z000010)
            { 0,     0,     3.04240080558681453113E+05},
            { 0,     1,    -7.24150360611511170864E+07},
            { 0,     2,     2.70646920133822739124E+08},
            { 0,     3,    -4.33003387311595380306E+08},
            { 0,     4,     4.03091375320128858089E+08},
            { 0,     5,    -2.49456192293578892946E+08},
            { 0,     6,     9.70010504787636250257E+07},
            { 0,     7,    -8.23838906384931178764E+05},
            { 0,     8,    -3.28722313181499056518E+07},
            { 0,     9,     2.78880451327459998429E+07},
            { 0,    10,    -1.40338722991762422025E+07},
            { 0,    11,     5.06007529781199898571E+06},
            { 0,    12,    -1.37089679529935028404E+06},
            { 0,    13,     2.68858224094756355044E+05},
            { 0,    14,    -3.36690685426766431192E+04},
            { 0,    15,     1.98948270157709225714E+03},

            { 1,     0,     5.37679127182313203812E+07},
            { 1,     1,     1.46388718178314149380E+08},
            { 1,     2,    -9.35395212179676651955E+08},
            { 1,     3,     1.55915548629785561562E+09},
            { 1,     4,    -1.34035874177675724030E+09},
            { 1,     5,     7.60583187397939205170E+08},
            { 1,     6,    -3.48267020550069510937E+08},
            { 1,     7,     1.27565538716441661119E+08},
            { 1,     8,    -2.52623368581141643226E+07},
            { 1,     9,    -2.38211218436619313434E+06},
            { 1,    10,     2.83049418278661277145E+06},
            { 1,    11,    -1.12351439162609330378E+06},
            { 1,    12,     5.12258502489926875569E+05},
            { 1,    13,    -1.87547968379564146744E+05},
            { 1,    14,     3.78142640863094275119E+04},
            { 1,    15,    -3.14819542193592815238E+03},

            { 2,     0,    -2.72150296501208603382E+08},
            { 2,     1,     2.54191144958757966757E+08},
            { 2,     2,     1.09221455545737147331E+09},
            { 2,     3,    -2.27038845949890232086E+09},
            { 2,     4,     1.74805773803573489189E+09},
            { 2,     5,    -6.91637455687782287598E+08},
            { 2,     6,     2.04606336638418078423E+08},
            { 2,     7,    -7.53842029037321358919E+07},
            { 2,     8,     1.98339755878359973431E+07},
            { 2,     9,    -2.51659290353266568854E+06},
            { 2,    10,     2.78171596346601564437E+06},
            { 2,    11,    -2.17292116712976712734E+06},
            { 2,    12,     6.62271515154829132371E+05},
            { 2,    13,    -6.85233401404126780108E+04},
            { 2,    14,    -7.24533940448010434920E+03},
            { 2,    15,     1.77091708853070326768E+03},

            { 3,     0,     5.84549742654761433601E+08},
            { 3,     1,    -1.22755266303496289253E+09},
            { 3,     2,    -1.99841301217404045165E+07},
            { 3,     3,     1.62942775427007746696E+09},
            { 3,     4,    -1.29122969498875904083E+09},
            { 3,     5,     2.81517572141482591629E+08},
            { 3,     6,     3.29046875321928299963E+07},
            { 3,     7,    -2.18449690962994331494E+06},
            { 3,     8,    -2.63601564748177304864E+06},
            { 3,     9,    -4.48942869859682396054E+06},
            { 3,    10,     2.44109624343642685562E+06},
            { 3,    11,    -4.88764745659762818832E+05},
            { 3,    12,     5.56227389265404126490E+04},
            { 3,    13,    -2.70979608753431966761E+03},
            { 3,    14,     1.47652876079566203771E+03},
            { 3,    15,    -5.79993594759073516798E+02},

            { 4,     0,    -6.80110049336398720741E+08},
            { 4,     1,     1.85585044735255384445E+09},
            { 4,     2,    -1.28215713107603287697E+09},
            { 4,     3,    -3.01044179228078424931E+08},
            { 4,     4,     5.51112592108273148537E+08},
            { 4,     5,    -7.78949714611403048038E+07},
            { 4,     6,    -3.93153402725159674883E+07},
            { 4,     7,    -1.19108830926259383559E+07},
            { 4,     8,     1.43100595134035050869E+07},
            { 4,     9,    -4.91528834089644718915E+06},
            { 4,    10,     1.68842217521909135394E+06},
            { 4,    11,    -1.75640131265274452744E+05},
            { 4,    12,    -1.40626545874571602326E+05},
            { 4,    13,     3.77159267968944986933E+04},
            { 4,    14,    -1.39635312769678466793E+03},
            { 4,    15,    -1.88574192421971709166E+01},

            { 5,     0,     4.31829985381816327572E+08},
            { 5,     1,    -1.41004882914948272705E+09},
            { 5,     2,     1.38554457511777806282E+09},
            { 5,     3,    -3.82333992940614879131E+08},
            { 5,     4,    -7.42225503488553911448E+07},
            { 5,     5,     1.56423537414295095950E+07},
            { 5,     6,     1.26731311307457205839E+06},
            { 5,     7,     8.38170131250077579170E+06},
            { 5,     8,     1.69866128198041976430E+06},
            { 5,     9,    -2.63807965752872545272E+06},
            { 5,    10,    -1.30410926100970245898E+05},
            { 5,    11,     2.79376974434808362275E+05},
            { 5,    12,    -2.47335903090722194975E+03},
            { 5,    13,    -8.09120832153944047604E+03},
            { 5,    14,    -1.18437059561102614680E+03},
            { 5,    15,     2.22406017739857361448E+02},

            { 6,     0,    -1.05728919394894704223E+08},
            { 6,     1,     4.78614212386758863926E+08},
            { 6,     2,    -5.47275303584568023682E+08},
            { 6,     3,     1.81002869935932725668E+08},
            { 6,     4,     2.37412980302353836596E+07},
            { 6,     5,    -1.49761239025783911347E+07},
            { 6,     6,     6.94548275970081239939E+06},
            { 6,     7,    -7.08765067577983532101E+06},
            { 6,     8,     6.15334952530000591651E+05},
            { 6,     9,     1.04307962681467283983E+06},
            { 6,    10,    -1.66822670094482309651E+05},
            { 6,    11,    -3.14983610603470224305E+04},
            { 6,    12,    -8.67230359049541766581E+03},
            { 6,    13,     3.35310141926492042330E+03},
            { 6,    14,     5.20457690967627286227E+02},
            { 6,    15,    -9.46886765048879937012E+01},

            { 7,     0,    -4.00641674232298433781E+07},
            { 7,     1,     3.53151559413958713412E+07},
            { 7,     2,    -3.32651288367051668465E+07},
            { 7,     3,     8.78685073832835257053E+07},
            { 7,     4,    -8.83196869217660874128E+07},
            { 7,     5,     3.53058029766001477838E+07},
            { 7,     6,    -2.53460763685103179887E+06},
            { 7,     7,    -3.27478091701123584062E+06},
            { 7,     8,     1.51544590128077217378E+06},
            { 7,     9,    -1.97898065117523074150E+05},
            { 7,    10,    -4.86252464679210970644E+04},
            { 7,    11,     3.00813760869401539821E+03},
            { 7,    12,     9.57301967925047574681E+03},
            { 7,    13,    -2.60470610160503747466E+03},
            { 7,    14,     1.71621252115364512747E+02},
            { 7,    15,    -4.79603274832587800347E+00},

            { 8,     0,     3.64528375028086900711E+07},
            { 8,     1,    -7.85125069899140596390E+07},
            { 8,     2,     7.79785535495323836803E+07},
            { 8,     3,    -5.67020461767953187227E+07},
            { 8,     4,     3.75835219815801903605E+07},
            { 8,     5,    -2.08943791429358497262E+07},
            { 8,     6,     7.66625923235765285790E+06},
            { 8,     7,    -1.06592743443907494657E+06},
            { 8,     8,    -3.14079768184977932833E+05},
            { 8,     9,     1.45893387703247601166E+05},
            { 8,    10,    -2.00583693587448651670E+04},
            { 8,    11,    -1.07287338880777724626E+03},
            { 8,    12,     1.86369211872262280849E+03},
            { 8,    13,    -5.41388893973812855620E+02},
            { 8,    14,     3.82187077215814881015E+01},
            { 8,    15,     4.76195893580739504358E+00},

            { 9,     0,    -9.92938956532538309693E+06},
            { 9,     1,     2.04695142809341102839E+07},
            { 9,     2,    -1.02214338619423378259E+07},
            { 9,     3,    -9.93348938946167938411E+06},
            { 9,     4,     1.64549253709395211190E+07},
            { 9,     5,    -9.86318830775812268257E+06},
            { 9,     6,     3.13746389613276068121E+06},
            { 9,     7,    -6.67069913326343754306E+05},
            { 9,     8,     2.16856145024514931720E+05},
            { 9,     9,    -1.08820304041782830609E+05},
            { 9,    10,     3.55738569172507050098E+04},
            { 9,    11,    -3.12188240990893382332E+03},
            { 9,    12,    -1.97305440808249090878E+03},
            { 9,    13,     7.01633632749421281005E+02},
            { 9,    14,    -8.03244389802174652004E+01},
            { 9,    15,     1.98516265048327333886E+00},

            {10,     0,     1.05106274467385723256E+06},
            {10,     1,    -1.90074146260506426916E+06},
            {10,     2,    -5.14866259700076945592E+05},
            {10,     3,     4.43767094488164782524E+06},
            {10,     4,    -5.41318389280360378325E+06},
            {10,     5,     3.27827314338852465153E+06},
            {10,     6,    -1.00315844098674086854E+06},
            {10,     7,     3.13072588745724533510E+04},
            {10,     8,     1.04236830211061693262E+05},
            {10,     9,    -4.78932228588938814937E+04},
            {10,    10,     1.31090834075820221187E+04},
            {10,    11,    -3.39479313977342053477E+03},
            {10,    12,     9.29207146820721618496E+02},
            {10,    13,    -1.86795988806825192796E+02},
            {10,    14,     1.98510067854014131683E+01},
            {10,    15,    -7.40917959154837491020E-01}
        },
        {                                                                               // RECOM (Z000010)
            { 0,     0,     1.43249838796439163957E+01},
            { 0,     1,    -1.34095062880606157307E+01},
            { 0,     2,     6.08550621146791144156E+01},
            { 0,     3,    -1.46577101984033276949E+02},
            { 0,     4,     2.10671492276514470632E+02},
            { 0,     5,    -1.92066044872448998149E+02},
            { 0,     6,     1.12528193269927839992E+02},
            { 0,     7,    -4.10571661215909387010E+01},
            { 0,     8,     8.47360078103331737509E+00},
            { 0,     9,    -7.53719632151268914555E-01},

            { 1,     0,    -4.69105628993968171159E+00},
            { 1,     1,     8.90187555533843095645E+01},
            { 1,     2,    -5.17082466752443224323E+02},
            { 1,     3,     1.49287044908038524227E+03},
            { 1,     4,    -2.42034351477910740869E+03},
            { 1,     5,     2.35977318225110639105E+03},
            { 1,     6,    -1.41747263752418666627E+03},
            { 1,     7,     5.14675357902906625895E+02},
            { 1,     8,    -1.03788126354235686222E+02},
            { 1,     9,     8.93875687034091015448E+00},

            { 2,     0,    -1.98086570192114912459E+01},
            { 2,     1,     9.02725480979547114657E+01},
            { 2,     2,     3.72962652775568585639E+02},
            { 2,     3,    -2.86443864935814417549E+03},
            { 2,     4,     6.54791371890192021965E+03},
            { 2,     5,    -7.57542293321982197085E+03},
            { 2,     6,     4.98863219830037542124E+03},
            { 2,     7,    -1.90003291421173662457E+03},
            { 2,     8,     3.91370534933327633098E+02},
            { 2,     9,    -3.38888734262211954729E+01},

            { 3,     0,    -7.61727485013823990556E+00},
            { 3,     1,    -8.99676221452425437519E+01},
            { 3,     2,     2.39429289955121817002E+02},
            { 3,     3,     2.05771307220595781473E+03},
            { 3,     4,    -7.66672108820838730026E+03},
            { 3,     5,     1.07610646452268938447E+04},
            { 3,     6,    -7.83052958018395565887E+03},
            { 3,     7,     3.14853537506878637942E+03},
            { 3,     8,    -6.66805677756379054699E+02},
            { 3,     9,     5.84483701602106719974E+01},

            { 4,     0,     6.84689749428262302899E+02},
            { 4,     1,    -6.41337423137322639377E+03},
            { 4,     2,     2.14353806903193471953E+04},
            { 4,     3,    -3.96806322513512059231E+04},
            { 4,     4,     4.68443039997480664169E+04},
            { 4,     5,    -3.66080543837553850608E+04},
            { 4,     6,     1.87025281715490273200E+04},
            { 4,     7,    -5.95491842750383239036E+03},
            { 4,     8,     1.06698781883290007499E+03},
            { 4,     9,    -8.22032529086119865269E+01},
            { 5,     0,    -1.08938704804548956417E+03},
            { 5,     1,     1.48146628563555968867E+04},
            { 5,     2,    -5.75561987341030326206E+04},
            { 5,     3,     1.11350161916194265359E+05},
            { 5,     4,    -1.26655054015904039261E+05},
            { 5,     5,     9.02577656921312009217E+04},
            { 5,     6,    -4.08346371174871092080E+04},
            { 5,     7,     1.13850781915498100716E+04},
            { 5,     8,    -1.78117510051617318823E+03},
            { 5,     9,     1.19593664869470444501E+02},

            { 6,     0,     7.71158520059956913428E+01},
            { 6,     1,    -1.14063205364834775537E+04},
            { 6,     2,     5.61224719137911379221E+04},
            { 6,     3,    -1.19398022502862586407E+05},
            { 6,     4,     1.41535388629036053317E+05},
            { 6,     5,    -1.01977152037983221817E+05},
            { 6,     6,     4.57325986781135798083E+04},
            { 6,     7,    -1.24671193954050559114E+04},
            { 6,     8,     1.88748088780042576218E+03},
            { 6,     9,    -1.21427200114438761602E+02},

            { 7,     0,     7.51859271237372013275E+02},
            { 7,     1,     2.60919535160936538887E+03},
            { 7,     2,    -2.43703099252169959072E+04},
            { 7,     3,     6.04402419973716023378E+04},
            { 7,     4,    -7.66051639442676823819E+04},
            { 7,     5,     5.70820587259594394709E+04},
            { 7,     6,    -2.60263704770619915507E+04},
            { 7,     7,     7.14161111890670417779E+03},
            { 7,     8,    -1.08138895268239821235E+03},
            { 7,     9,     6.92421978398739810245E+01},

            { 8,     0,    -4.88596016673188898949E+02},
            { 8,     1,     6.06324096291598948483E+02},
            { 8,     2,     4.18032368029227109218E+03},
            { 8,     3,    -1.41908964780536971375E+04},
            { 8,     4,     1.99120584653974874527E+04},
            { 8,     5,    -1.55743938039613858564E+04},
            { 8,     6,     7.29395282901822793065E+03},
            { 8,     7,    -2.03371603812545799883E+03},
            { 8,     8,     3.11077689089091336427E+02},
            { 8,     9,    -2.00523533590040692332E+01},

            { 9,     0,     9.23231057243035166948E+01},
            { 9,     1,    -2.49558777666313233112E+02},
            { 9,     2,    -9.91094586704329287841E+01},
            { 9,     3,     1.18411592281756065859E+03},
            { 9,     4,    -1.97229180010892332575E+03},
            { 9,     5,     1.65110538528935694558E+03},
            { 9,     6,    -8.01140239767901903178E+02},
            { 9,     7,     2.28187127507413094918E+02},
            { 9,     8,    -3.54063664898269081505E+01},
            { 9,     9,     2.30675599774240458473E+00}
        }
    },
    {                                                                                   // Metallicity Z000100 (0.00100)
        {                                                                               // LMR1 (Z000100)
            { 0,     0,     1.50346099257731538046E+01},
            { 0,     1,     2.82604496789430559289E+00},
            { 0,     2,    -9.81648621516211328242E+00},
            { 0,     3,     9.56667166719136474740E+00},
            { 0,     4,    -3.13477224724504077713E+00},

            { 1,     0,     2.07944042796325678779E+00},
            { 1,     1,    -5.06461279977239176020E+00},
            { 1,     2,     2.03394884719189121824E+01},
            { 1,     3,    -1.93735154009608017134E+01},
            { 1,     4,     5.75610929539014914980E+00},

            { 2,     0,    -1.05127268925442436398E+01},
            { 2,     1,     3.68098490829458810936E+01},
            { 2,     2,    -5.91474113982236531228E+01},
            { 2,     3,     3.16463783616280913691E+01},
            { 2,     4,    -3.65289167448958851381E+00},

            { 3,     0,     2.22020213350119419715E+01},
            { 3,     1,    -7.77940541900198496705E+01},
            { 3,     2,     1.05860178274410074550E+02},
            { 3,     3,    -4.94197227165773327329E+01},
            { 3,     4,     3.82693132506987110375E+00},

            { 4,     0,    -1.12988524726900649853E+01},
            { 4,     1,     4.13367947020142807446E+01},
            { 4,     2,    -5.60688828217752686101E+01},
            { 4,     3,     2.71194601089813360772E+01},
            { 4,     4,    -2.72130850772848109642E+00}
        },
        {                                                                               // LMR2 (Z000100)
            { 0,     0,     1.88453218130641353412E+01},
            { 0,     1,    -1.66303719796976672285E+01},
            { 0,     2,     2.63222496911415824172E+01},
            { 0,     3,    -2.13131516790098238801E+01},
            { 0,     4,     8.39818275012596870965E+00},
            { 0,     5,    -1.28795907731821102082E+00},

            { 1,     0,     4.31329741689262391446E+01},
            { 1,     1,    -1.56619568143149109574E+02},
            { 1,     2,     2.26032050176760577642E+02},
            { 1,     3,    -1.59897885791484213769E+02},
            { 1,     4,     5.55606561823480475937E+01},
            { 1,     5,    -7.60112702919777394328E+00},

            { 2,     0,     2.46987072824178248709E+01},
            { 2,     1,     4.20629378936399973554E+01},
            { 2,     2,    -2.65767178117826063044E+02},
            { 2,     3,     3.25730032378355303990E+02},
            { 2,     4,    -1.58371785651665220485E+02},
            { 2,     5,     2.75938397603546903269E+01},

            { 3,     0,    -9.91381643467949629667E+02},
            { 3,     1,     3.38104557879559251887E+03},
            { 3,     2,    -4.11533190831385491038E+03},
            { 3,     3,     2.33454985382554878015E+03},
            { 3,     4,    -6.14157763172000272789E+02},
            { 3,     5,     5.77076386143955843977E+01},

            { 4,     0,     1.43797791534433781635E+03},
            { 4,     1,    -5.28916096082312014914E+03},
            { 4,     2,     6.91057399431488920527E+03},
            { 4,     3,    -4.26364756616708109505E+03},
            { 4,     4,     1.25783448603752435702E+03},
            { 4,     5,    -1.41377167956629875789E+02},

            { 5,     0,    -5.01761058636202108119E+02},
            { 5,     1,     1.99253292955415395227E+03},
            { 5,     2,    -2.72417925456980492527E+03},
            { 5,     3,     1.74907541356862247994E+03},
            { 5,     4,    -5.38448511479686771963E+02},
            { 5,     5,     6.36637946720314573668E+01}
        },
        {                                                                               // LMA (Z000100)
            { 0,     0,     7.69284893177311052568E+05},
            { 0,     1,    -4.60964236746346391737E+06},
            { 0,     2,     9.24477105650694482028E+06},
            { 0,     3,     2.96809810085771558806E+06},
            { 0,     4,    -5.47904367331531867385E+07},
            { 0,     5,     1.38449097585996419191E+08},
            { 0,     6,    -2.05590189063783049583E+08},
            { 0,     7,     2.11734651675154030323E+08},
            { 0,     8,    -1.60435957671923041344E+08},
            { 0,     9,     9.17609013885079324245E+07},
            { 0,    10,    -3.99674217866430878639E+07},
            { 0,    11,     1.32154540634083170444E+07},
            { 0,    12,    -3.26794649031273694709E+06},
            { 0,    13,     5.86210662802961654961E+05},
            { 0,    14,    -7.21158724758123280481E+04},
            { 0,    15,     5.44588834510656215571E+03},
            { 0,    16,    -1.90445008267364357835E+02},

            { 2,     0,     2.21188523097329214215E+07},
            { 2,     1,    -1.37731912987798571587E+08},
            { 2,     2,     3.57687151988543748856E+08},
            { 2,     3,    -4.63590327267392992973E+08},
            { 2,     4,     2.01845788407205879688E+08},
            { 2,     5,     2.87663499872765600681E+08},
            { 2,     6,    -5.59995404294795393944E+08},
            { 2,     7,     4.26777807943069577217E+08},
            { 2,     8,    -1.24264635254063099623E+08},
            { 2,     9,    -6.53086944807844161987E+07},
            { 2,    10,     9.24923524751316756010E+07},
            { 2,    11,    -5.24022445798896402121E+07},
            { 2,    12,     1.83838309755647853017E+07},
            { 2,    13,    -4.26215335990773420781E+06},
            { 2,    14,     6.40775684917965554632E+05},
            { 2,    15,    -5.69637563759163822397E+04},
            { 2,    16,     2.28325195741889183410E+03},

            { 4,     0,    -5.86271602527818381786E+07},
            { 4,     1,     4.78369991991747081280E+08},
            { 4,     2,    -1.74361930468868136406E+09},
            { 4,     3,     3.76621819705243778229E+09},
            { 4,     4,    -5.37751677961362648010E+09},
            { 4,     5,     5.35077163876746082306E+09},
            { 4,     6,    -3.80998322478477764130E+09},
            { 4,     7,     1.97099712467356157303E+09},
            { 4,     8,    -7.60832823835394382477E+08},
            { 4,     9,     2.42272232218582451344E+08},
            { 4,    10,    -8.05926269997112601995E+07},
            { 4,    11,     3.11124469145792499185E+07},
            { 4,    12,    -1.09693358577220626175E+07},
            { 4,    13,     2.84038850674086064100E+06},
            { 4,    14,    -4.84199551332909730263E+05},
            { 4,    15,     4.84970404654026060598E+04},
            { 4,    16,    -2.17011384837293599048E+03},

            { 6,     0,     9.83702780789473533630E+08},
            { 6,     1,    -6.23115263299602222443E+09},
            { 6,     2,     1.74950104067897109985E+10},
            { 6,     3,    -2.85244593670344696045E+10},
            { 6,     4,     2.95578184637604980469E+10},
            { 6,     5,    -1.97934299705980529785E+10},
            { 6,     6,     8.07257493996787166595E+09},
            { 6,     7,    -1.50554167974234604836E+09},
            { 6,     8,    -1.18252537929885953665E+08},
            { 6,     9,    -3.70554775613283040002E+06},
            { 6,    10,     1.23448486340704277158E+08},
            { 6,    11,    -7.41264279707325249910E+07},
            { 6,    12,     2.09071950705411061645E+07},
            { 6,    13,    -2.96454796610304014757E+06},
            { 6,    14,     1.18142086018433197751E+05},
            { 6,    15,     1.92799522383320254448E+04},
            { 6,    16,    -1.90271167151171880505E+03},

            { 8,     0,     6.69561843877016976476E+07},
            { 8,     1,     1.56256926101374864578E+09},
            { 8,     2,    -8.75745768682143020630E+09},
            { 8,     3,     1.94779956266501197815E+10},
            { 8,     4,    -2.33126852312486190796E+10},
            { 8,     5,     1.55144418171162052155E+10},
            { 8,     6,    -4.34538685387376117706E+09},
            { 8,     7,    -1.24875014824374175072E+09},
            { 8,     8,     1.36300258540026307106E+09},
            { 8,     9,    -2.47726977709512680769E+08},
            { 8,    10,    -1.37758820791601359844E+08},
            { 8,    11,     7.42609580452102422714E+07},
            { 8,    12,    -5.21190464590375591069E+06},
            { 8,    13,    -5.56763955730174109340E+06},
            { 8,    14,     1.98916589357259985991E+06},
            { 8,    15,    -2.82567601932666904759E+05},
            { 8,    16,     1.54458629248594224919E+04},

            {10,     0,    -2.91743516387244606018E+09},
            {10,     1,     1.43474200087861461639E+10},
            {10,     2,    -2.97121547997276992798E+10},
            {10,     3,     3.29244879152800636292E+10},
            {10,     4,    -1.97049122193543739319E+10},
            {10,     5,     4.44795698088213920593E+09},
            {10,     6,     1.68826391791642093658E+09},
            {10,     7,    -1.34669317233685207367E+09},
            {10,     8,     3.03101724866939246655E+08},
            {10,     9,    -4.99869293617860823870E+07},
            {10,    10,     1.42133397069574892521E+07},
            {10,    11,     2.69422055173848085105E+07},
            {10,    12,    -3.11643383659875914454E+07},
            {10,    13,     1.40837340470498204231E+07},
            {10,    14,    -3.33681038689361186698E+06},
            {10,    15,     4.13678001847658306360E+05},
            {10,    16,    -2.12956128078225046920E+04},

            {12,     0,     3.83066159035058796406E+08},
            {12,     1,    -3.57522953798943567276E+09},
            {12,     2,     1.00626903613678798676E+10},
            {12,     3,    -1.31025349462869167328E+10},
            {12,     4,     7.85384311503630828857E+09},
            {12,     5,    -2.41011892096044540405E+08},
            {12,     6,    -2.56408374643986368179E+09},
            {12,     7,     1.39171229309356403351E+09},
            {12,     8,    -1.44009585297834306955E+08},
            {12,     9,    -8.26904684745683073997E+07},
            {12,    10,    -2.38144140223466372117E+06},
            {12,    11,     1.79516107555324509740E+07},
            {12,    12,    -1.88674114853504090570E+06},
            {12,    13,    -2.42738823937244340777E+06},
            {12,    14,     9.99039133402326726355E+05},
            {12,    15,    -1.56469409353959286818E+05},
            {12,    16,     9.15292969974886182172E+03},

            {14,     0,     3.09020185410522365570E+09},
            {14,     1,    -1.36212178749835147858E+10},
            {14,     2,     2.66881548670487785339E+10},
            {14,     3,    -3.03946970897722930908E+10},
            {14,     4,     2.17784943226391410828E+10},
            {14,     5,    -9.40680629512597274780E+09},
            {14,     6,     1.39847537299111199379E+09},
            {14,     7,     1.10993825469740843773E+09},
            {14,     8,    -8.73110848697505235672E+08},
            {14,     9,     2.50292736143364071846E+08},
            {14,    10,     1.06462010817974358797E+07},
            {14,    11,    -3.22772786734968498349E+07},
            {14,    12,     1.01494283884754441679E+07},
            {14,    13,    -1.02655365919480670709E+06},
            {14,    14,    -1.37982127262780006276E+05},
            {14,    15,     4.30180039417772932211E+04},
            {14,    16,    -3.02473938756877396372E+03},

            {16,     0,    -1.56852265694742107391E+09},
            {16,     1,     7.15821216061697673798E+09},
            {16,     2,    -1.43193921312494201660E+10},
            {16,     3,     1.61786574433440303802E+10},
            {16,     4,    -1.08238279920512008667E+10},
            {16,     5,     3.66439506765549850464E+09},
            {16,     6,     3.04446730119516968727E+08},
            {16,     7,    -9.95706614031146168709E+08},
            {16,     8,     5.34070064296932458878E+08},
            {16,     9,    -1.68596081582052469254E+08},
            {16,    10,     4.50922229884207323194E+07},
            {16,    11,    -1.64262273564462922513E+07},
            {16,    12,     6.76394867549139633775E+06},
            {16,    13,    -2.05885415643938211724E+06},
            {16,    14,     3.96860804045895056333E+05},
            {16,    15,    -4.34586026756422943436E+04},
            {16,    16,     2.06736711765395784823E+03}
        },
        {                                                                               // HM (Z000100)
            { 0,     0,     2.66920628896922789863E+04},
            { 0,     1,     5.69598144516321923584E+06},
            { 0,     2,    -3.36106590539333373308E+07},
            { 0,     3,     9.25393487313766926527E+07},
            { 0,     4,    -1.64881303894646197557E+08},
            { 0,     5,     2.12786252918928653002E+08},
            { 0,     6,    -2.02235432385967493057E+08},
            { 0,     7,     1.36366495759085386992E+08},
            { 0,     8,    -5.89415461827230826020E+07},
            { 0,     9,     1.05569686857345979661E+07},
            { 0,    10,     4.26705694591859634966E+06},
            { 0,    11,    -3.31228471217951877043E+06},
            { 0,    12,     6.44460078807694022544E+05},
            { 0,    13,     1.79727468496846966445E+05},
            { 0,    14,    -1.14602002866395123419E+05},
            { 0,    15,     8.73345725996870532981E+03},
            { 0,    16,     1.03006575184149296547E+04},
            { 0,    17,    -4.45924381157494190120E+03},
            { 0,    18,     8.75712682690710266797E+02},
            { 0,    19,    -8.97432677867381869419E+01},
            { 0,    20,     3.90904107825454749658E+00},

            { 2,     0,    -4.66175954969746712595E+06},
            { 2,     1,     1.44553042693169433624E+07},
            { 2,     2,    -1.38532479348339065909E+07},
            { 2,     3,     1.61449323716102633625E+07},
            { 2,     4,    -5.04766978015377074480E+07},
            { 2,     5,     8.52236769195135086775E+07},
            { 2,     6,    -7.52658732323246151209E+07},
            { 2,     7,     3.77911912456169873476E+07},
            { 2,     8,    -1.05537492595233544707E+07},
            { 2,     9,     1.23416606663387943991E+06},
            { 2,    10,     1.07657456006634638470E+04},
            { 2,    11,     1.43573359232634276850E+05},
            { 2,    12,    -1.42859383483003708534E+05},
            { 2,    13,     5.63552937440058449283E+04},
            { 2,    14,    -1.09635084253284767328E+04},
            { 2,    15,     1.32853574780430949431E+03},
            { 2,    16,    -4.14076353943859373885E+02},
            { 2,    17,     1.47253804236653735416E+02},
            { 2,    18,    -2.47527050904476979554E+01},
            { 2,    19,     1.91210492011224597597E+00},
            { 2,    20,    -7.05892167374512213840E-02},

            { 4,     0,     1.02998189728514533490E+07},
            { 4,     1,    -3.48644858077521771193E+07},
            { 4,     2,     2.46925453657950051129E+07},
            { 4,     3,     2.65148321726851537824E+07},
            { 4,     4,    -3.81589329059887528419E+07},
            { 4,     5,     6.93028380426791380160E+05},
            { 4,     6,     1.99752902745822705328E+07},
            { 4,     7,    -1.14559117107986155897E+07},
            { 4,     8,     2.16656442928818566725E+06},
            { 4,     9,    -2.71810134198548854329E+05},
            { 4,    10,     2.68526744210913369898E+05},
            { 4,    11,    -9.77884134982192044845E+04},
            { 4,    12,    -1.53021212986694445135E+03},
            { 4,    13,     1.19037018096100418916E+04},
            { 4,    14,    -5.85794015733223022835E+03},
            { 4,    15,     1.56346892706033872855E+03},
            { 4,    16,    -2.75046129788926521087E+02},
            { 4,    17,     1.01298210263538095433E+02},
            { 4,    18,    -3.87720363616409713359E+01},
            { 4,    19,     6.66209009612612579332E+00},
            { 4,    20,    -3.93432972239732559050E-01},

            { 6,     0,    -1.23277652244743425399E+07},
            { 6,     1,     5.31219698912802562118E+07},
            { 6,     2,    -6.70881175458954200149E+07},
            { 6,     3,     1.38345048209509141743E+07},
            { 6,     4,     2.72927585462916977704E+07},
            { 6,     5,    -1.06984214016320873052E+07},
            { 6,     6,    -1.07550707413996290416E+07},
            { 6,     7,     9.12851204341692477465E+06},
            { 6,     8,    -2.45932699648965429515E+06},
            { 6,     9,     3.79086684429724293295E+05},
            { 6,    10,    -1.15384063428619745537E+05},
            { 6,    11,     9.37208910898514295695E+03},
            { 6,    12,     5.50001460585494169209E+03},
            { 6,    13,     8.99582878278114890236E+02},
            { 6,    14,     2.17976963642488072992E+02},
            { 6,    15,    -2.90083472695319642298E+02},
            { 6,    16,    -1.06912609256554702597E+02},
            { 6,    17,     7.44274626977634596869E+01},
            { 6,    18,    -1.32865712570078340349E+01},
            { 6,    19,     1.15176452822245067864E+00},
            { 6,    20,    -7.52331325998904648644E-02},

            { 8,     0,     6.81678287605022825301E+06},
            { 8,     1,    -3.98179300438016653061E+07},
            { 8,     2,     6.65612044507603123784E+07},
            { 8,     3,    -3.76181603687587901950E+07},
            { 8,     4,    -6.56919684995194338262E+06},
            { 8,     5,     1.43847081179449260235E+07},
            { 8,     6,    -4.06374374651811597869E+06},
            { 8,     7,     2.91022107704004331026E+05},
            { 8,     8,    -5.34165458237739163451E+05},
            { 8,     9,     2.58318759598450968042E+05},
            { 8,    10,    -2.04026808616285597964E+03},
            { 8,    11,    -1.91223168630281179503E+03},
            { 8,    12,    -3.69962055497544270111E+03},
            { 8,    13,    -9.23920259961619422029E+02},
            { 8,    14,     4.11256827202935312471E+02},
            { 8,    15,     1.26735718583368523582E+02},
            { 8,    16,    -4.68765072065054724249E+01},
            { 8,    17,     6.13281037162783970729E+00},
            { 8,    18,    -8.50877614045823293942E-01},
            { 8,    19,    -1.89283136262390905280E-01},
            { 8,    20,     5.75261353234685565705E-02},

            {10,     0,    -3.55154449666948523372E+05},
            {10,     1,     1.30929190946150850505E+07},
            {10,     2,    -3.15488824688208103180E+07},
            {10,     3,     2.70812807114166915417E+07},
            {10,     4,    -6.85649736623421218246E+06},
            {10,     5,    -2.40983373484842292964E+06},
            {10,     6,     1.08043478886758117005E+06},
            {10,     7,     1.11299945145412624697E+05},
            {10,     8,     4.58449855702385902987E+04},
            {10,     9,    -5.80789547527073082165E+04},
            {10,    10,    -4.28841803810808505659E+03},
            {10,    11,     4.99567543157751151739E+03},
            {10,    12,    -5.59329152251473260549E+02},
            {10,    13,     3.66269193209130548894E+02},
            {10,    14,    -4.12453739939884940213E+01},
            {10,    15,    -1.03224939821626904290E+01},
            {10,    16,    -1.26744942581002923987E+01},
            {10,    17,     2.50242050167653262704E+00},
            {10,    18,     4.89212292141251114952E-01},
            {10,    19,     7.41773856600715748508E-04},
            {10,    20,    -2.14966363670156049293E-02},

            {12,     0,    -1.71770314663222013041E+06},
            {12,     1,     1.02925354974316677544E+06},
            {12,     2,     4.75587077003400865942E+06},
            {12,     3,    -6.94997276824742369354E+06},
            {12,     4,     2.96480201554762385786E+06},
            {12,     5,     1.05305805605570232728E+05},
            {12,     6,    -2.38596544252316496568E+05},
            {12,     7,    -5.79872394392845744733E+04},
            {12,     8,     3.01928438393556607480E+04},
            {12,     9,     1.89217764968022629546E+03},
            {12,    10,    -4.88373167221754613365E+02},
            {12,    11,     7.01114489071818979937E+02},
            {12,    12,    -3.17674857603029806796E+02},
            {12,    13,    -9.41372505874503957557E+01},
            {12,    14,     1.83240842852757204184E+01},
            {12,    15,     1.08114635042327797976E+01},
            {12,    16,    -1.76741740095267596544E+00},
            {12,    17,     6.94852470861354931664E-01},
            {12,    18,    -2.89336670178484411942E-01},
            {12,    19,    -3.66873283517264768896E-03},
            {12,    20,     7.64308669905610499340E-03},

            {14,     0,     9.98775663821288966574E+05},
            {14,     1,    -2.18451837998499209061E+06},
            {14,     2,     1.40499697486835205927E+06},
            {14,     3,     4.46811799451349797891E+04},
            {14,     4,    -2.46225148374797863653E+05},
            {14,     5,    -7.05764615857048920589E+04},
            {14,     6,     9.16282262372330296785E+04},
            {14,     7,    -1.90152534930542642542E+04},
            {14,     8,     4.91854753979303222877E+03},
            {14,     9,    -2.98215458946964145071E+03},
            {14,    10,    -1.24575514575863394384E+02},
            {14,    11,     3.07045703026585499629E+02},
            {14,    12,     2.56707666795900628642E+01},
            {14,    13,    -2.55503462019516938142E+01},
            {14,    14,    -7.83210099348742638803E-01},
            {14,    15,     2.75956632611307561831E+00},
            {14,    16,    -6.89935042869725956294E-01},
            {14,    17,    -9.58434272104521506330E-02},
            {14,    18,     3.10995999077311284509E-02},
            {14,    19,     1.11788217615811057148E-02},
            {14,    20,    -2.50872008751619425190E-03},

            {16,     0,    -2.56910131149193679448E+05},
            {16,     1,     6.59413942707619862631E+05},
            {16,     2,    -5.99113926101273740642E+05},
            {16,     3,     1.79299729173374042148E+05},
            {16,     4,     3.57433247453356452752E+04},
            {16,     5,    -1.71963241067330600345E+04},
            {16,     6,    -8.07591195715742651373E+03},
            {16,     7,     1.64482614549100276236E+03},
            {16,     8,     1.58385752216435230366E+03},
            {16,     9,    -3.77870980619565330016E+02},
            {16,    10,    -1.99848538469132250839E+01},
            {16,    11,    -2.99325126666875540593E+01},
            {16,    12,     1.39774853358364321565E+01},
            {16,    13,     2.38884594919343573594E+00},
            {16,    14,    -1.77569226381652423008E+00},
            {16,    15,     2.51906789777950501641E-01},
            {16,    16,    -9.35456748727026554668E-02},
            {16,    17,     4.71639212726435025358E-02},
            {16,    18,    -2.76687708834074250132E-05},
            {16,    19,    -3.76375892191212375881E-03},
            {16,    20,     5.43491168280245645454E-04},

            {18,     0,     3.26573367579413752537E+04},
            {18,     1,    -8.28040252641195256729E+04},
            {18,     2,     6.32083466612035263097E+04},
            {18,     3,     1.10908655478984019283E+04},
            {18,     4,    -4.24696541649424107163E+04},
            {18,     5,     2.13154570304911503626E+04},
            {18,     6,    -4.01758333934846973534E+02},
            {18,     7,    -3.02668804705343745809E+03},
            {18,     8,     9.00317782135503534846E+02},
            {18,     9,    -1.71590190988977440156E+01},
            {18,    10,    -7.58773236769787295941E+00},
            {18,    11,    -1.47039559521070568593E+01},
            {18,    12,     5.86203446077681977755E+00},
            {18,    13,    -4.35902509094680423729E-01},
            {18,    14,    -3.96411294687243842549E-01},
            {18,    15,     1.90872936553090594147E-01},
            {18,    16,    -3.63349815421271857274E-02},
            {18,    17,     4.18942750214783450613E-03},
            {18,    18,    -1.95696242731476692869E-03},
            {18,    19,     6.55720172656096495292E-04},
            {18,    20,    -6.89531819695781068207E-05},

            {20,     0,    -1.79757648979984833204E+03},
            {20,     1,     4.21453422482446967479E+03},
            {20,     2,    -1.45713195934336681603E+03},
            {20,     3,    -5.27886800091239820176E+03},
            {20,     4,     8.07499233626505065331E+03},
            {20,     5,    -5.33458308189604667859E+03},
            {20,     6,     1.84109344300976272280E+03},
            {20,     7,    -2.89109084486973813455E+02},
            {20,     8,     5.54488971681655407053E+00},
            {20,     9,    -2.55037857514060917197E+00},
            {20,    10,     3.28120362746293148248E+00},
            {20,    11,    -1.43670757857304920435E+00},
            {20,    12,     8.25401835544457118665E-01},
            {20,    13,    -3.05993525251011089239E-01},
            {20,    14,     2.31463397238021693914E-02},
            {20,    15,     2.54799453700065604844E-02},
            {20,    16,    -1.16997904598278249649E-02},
            {20,    17,     2.34703472519215026668E-03},
            {20,    18,    -1.64174625118059316908E-04},
            {20,    19,    -1.87233906354765485458E-05},
            {20,    20,     3.01646763959664224869E-06}
        },
        {                                                                               // RECOM (Z000100)
            { 0,     0,     1.29041924400426033515E+01},
            { 0,     1,     1.46247169347020422592E+00},
            { 0,     2,    -2.56099379621667111451E+00},
            { 0,     3,     1.53157682928077298889E+00},
            { 0,     4,    -3.12651080854979501744E-01},

            { 1,     0,     1.13224784437094316836E+00},
            { 1,     1,    -1.21990615381324185584E+00},
            { 1,     2,     4.03103598240615390580E+00},
            { 1,     3,    -2.98697934479165105870E+00},
            { 1,     4,     7.11864900599255223668E-01},

            { 2,     0,    -1.84989977300083419109E+00},
            { 2,     1,     2.98072414316297784609E+00},
            { 2,     2,    -4.44057678245443998577E+00},
            { 2,     3,     2.62017135972970738322E+00},
            { 2,     4,    -6.33529954639474812694E-01},

            { 3,     0,     3.17264901652564645929E+00},
            { 3,     1,    -4.83562740400180235412E+00},
            { 3,     2,     4.06335321968123874825E+00},
            { 3,     3,    -1.38974858743467533095E+00},
            { 3,     4,     2.62122104574736836113E-01},

            { 4,     0,    -2.02663383998561608124E+00},
            { 4,     1,     3.18617136072649298484E+00},
            { 4,     2,    -2.29903125672929053991E+00},
            { 4,     3,     5.58176645345621169625E-01},
            { 4,     4,    -5.91094511981828871217E-02},

            { 5,     0,     3.82436170169052624956E-01},
            { 5,     1,    -6.15679347898177353748E-01},
            { 5,     2,     4.45802375081275292779E-01},
            { 5,     3,    -1.00232559556148334567E-01},
            { 5,     4,     6.38538257034518188376E-03}
        }
    },
    {                                                                                   // Metallicity Z001000 (0.01000)
        {                                                                               // LMR1 (Z001000)
            { 0,     0,     1.61732010479488863552E+01},
            { 0,     1,    -6.07511035663581200339E+00},
            { 0,     2,     8.58287622972541264232E+00},
            { 0,     3,    -3.10585034200754250833E+00},
            { 0,     4,    -7.06748010502642909358E-01},

            { 1,     0,    -1.90469269583824840630E+00},
            { 1,     1,     4.00007719345960026658E+01},
            { 1,     2,    -8.46518605722360462096E+01},
            { 1,     3,     5.91797386241528400319E+01},
            { 1,     4,    -1.13107588613914398223E+01},

            { 2,     0,    -8.70230071618597023075E+00},
            { 2,     1,    -2.89910577605233434895E+01},
            { 2,     2,     1.24597859970597440338E+02},
            { 2,     3,    -1.16207323210657051504E+02},
            { 2,     4,     3.01232264093183879083E+01},

            { 3,     0,     1.54664218136763143008E+01},
            { 3,     1,    -1.42381111546532466150E+01},
            { 3,     2,    -5.53996704254122320776E+01},
            { 3,     3,     7.85240667339935924929E+01},
            { 3,     4,    -2.52702301461299612129E+01},

            { 4,     0,    -4.23582571635512561414E+00},
            { 4,     1,     9.21228272981120710483E+00},
            { 4,     2,     6.44366319286792332832E+00},
            { 4,     3,    -1.82412432120631002874E+01},
            { 4,     4,     7.14935319397919499806E+00}
        },
        {                                                                               // LMR2 (Z001000)
            { 0,     0,     1.59314592243044277353E+01},
            { 0,     1,    -4.60401610759330548461E+00},
            { 0,     2,     6.24185962640975233739E+00},
            { 0,     3,    -4.97550125746145077699E+00},
            { 0,     4,     1.93048999642347540728E+00},
            { 0,     5,    -2.89477547622687292339E-01},

            { 1,     0,     1.01061304060731380616E+01},
            { 1,     1,    -3.46086484084549255158E+01},
            { 1,     2,     5.27667822465917879526E+01},
            { 1,     3,    -3.86837936682783691822E+01},
            { 1,     4,     1.36020700266615026663E+01},
            { 1,     5,    -1.85012309699758370485E+00},

            { 2,     0,     2.00677620782971253277E+01},
            { 2,     1,    -2.26769827650059667690E+01},
            { 2,     2,    -1.15679690058618760418E+01},
            { 2,     3,    -5.80771237470029455530E-01},
            { 2,     4,     1.37546460492662188102E+01},
            { 2,     5,    -4.53156102738701349608E+00},

            { 3,     0,    -5.09057297132129690453E+02},
            { 3,     1,     1.15813297915374073455E+03},
            { 3,     2,    -8.67528490581213873156E+02},
            { 3,     3,     3.72884696292590604116E+02},
            { 3,     4,    -1.38863311991984033966E+02},
            { 3,     5,     2.83694073141493667833E+01},

            { 4,     0,     2.61623056004032514466E+03},
            { 4,     1,    -6.70644555846752700745E+03},
            { 4,     2,     6.13584072176503832452E+03},
            { 4,     3,    -2.70631835697935139251E+03},
            { 4,     4,     6.36762028562752561811E+02},
            { 4,     5,    -7.14561574114333950547E+01},

            { 5,     0,    -2.63489815083014082120E+03},
            { 5,     1,     7.38504520674821560533E+03},
            { 5,     2,    -7.50607780122314488835E+03},
            { 5,     3,     3.63058087806650701168E+03},
            { 5,     4,    -8.69070310941866864596E+02},
            { 5,     5,     8.58816125753537704668E+01}
        },
        {                                                                               // LMA (Z001000)
            { 0,     0,    -2.09773228588053607382E+05},
            { 0,     1,     1.43389594778738403693E+06},
            { 0,     2,    -7.75133130475684534758E+06},
            { 0,     3,     3.63601069582714810967E+07},
            { 0,     4,    -1.18813444845486640930E+08},
            { 0,     5,     2.60283723461131542921E+08},
            { 0,     6,    -3.91364485037142515182E+08},
            { 0,     7,     4.10038222587301969528E+08},
            { 0,     8,    -2.94145223515916526318E+08},
            { 0,     9,     1.30304661949904114008E+08},
            { 0,    10,    -1.72744424556804075837E+07},
            { 0,    11,    -2.00820756134503111243E+07},
            { 0,    12,     1.57965446997256800532E+07},
            { 0,    13,    -5.16442762537559401244E+06},
            { 0,    14,     2.88797206094442051835E+05},
            { 0,    15,     4.85840419889588956721E+05},
            { 0,    16,    -2.37449358671291265637E+05},
            { 0,    17,     5.89009159964685095474E+04},
            { 0,    18,    -8.73970189435798602062E+03},
            { 0,    19,     7.40523972585050728412E+02},
            { 0,    20,    -2.77978725470881755655E+01},

            { 2,     0,     1.24556972785087689757E+08},
            { 2,     1,    -1.05485290802952396870E+09},
            { 2,     2,     4.02473985730900144577E+09},
            { 2,     3,    -9.11583314025122642517E+09},
            { 2,     4,     1.35512409280278778076E+10},
            { 2,     5,    -1.37244791677059879303E+10},
            { 2,     6,     9.40152573164973068237E+09},
            { 2,     7,    -4.02953830579527378082E+09},
            { 2,     8,     7.38075188234021663666E+08},
            { 2,     9,     2.04097261853698253632E+08},
            { 2,    10,    -1.31296867969154432416E+08},
            { 2,    11,    -3.72527087734079034999E+06},
            { 2,    12,     2.13105948168305270374E+07},
            { 2,    13,    -5.90374779784860368818E+06},
            { 2,    14,    -1.58750756547301629325E+05},
            { 2,    15,     1.73802705629212403437E+05},
            { 2,    16,     1.34389908452468196629E+05},
            { 2,    17,    -8.46886743101308384212E+04},
            { 2,    18,     2.01200486051142797805E+04},
            { 2,    19,    -2.33582333582660248794E+03},
            { 2,    20,     1.10848048996528689258E+02},

            { 4,     0,     5.05357352773523867130E+08},
            { 4,     1,    -3.31354820566867399216E+09},
            { 4,     2,     9.55477016233197402954E+09},
            { 4,     3,    -1.56750020973316497803E+10},
            { 4,     4,     1.55897307926384067535E+10},
            { 4,     5,    -8.83530815627405929565E+09},
            { 4,     6,     1.63784678158393263817E+09},
            { 4,     7,     1.26310745307742381096E+09},
            { 4,     8,    -9.06822750568983674049E+08},
            { 4,     9,     1.67932512482709646225E+08},
            { 4,    10,    -1.23911428774325661361E+07},
            { 4,    11,     5.71523280062625110149E+07},
            { 4,    12,    -4.26682323901077881455E+07},
            { 4,    13,     8.90929416129896789789E+06},
            { 4,    14,     1.67166402933667809702E+06},
            { 4,    15,    -4.69427729580903833266E+05},
            { 4,    16,    -5.00381548112256336026E+05},
            { 4,    17,     3.01400804692653589882E+05},
            { 4,    18,    -7.32120482670434430474E+04},
            { 4,    19,     8.80328726870056379994E+03},
            { 4,    20,    -4.34033934541470671320E+02},

            { 6,     0,    -1.07852724296717691422E+09},
            { 6,     1,     6.84582537129478645325E+09},
            { 6,     2,    -1.98636935069379844666E+10},
            { 6,     3,     3.49891808203243942261E+10},
            { 6,     4,    -4.16472301939004364014E+10},
            { 6,     5,     3.48742082858454284668E+10},
            { 6,     6,    -2.03515912933958663940E+10},
            { 6,     7,     7.63230475433247947693E+09},
            { 6,     8,    -1.32354384835169196129E+09},
            { 6,     9,    -1.69552074669079929590E+08},
            { 6,    10,     6.97440331351650208235E+07},
            { 6,    11,     4.60175066012795269489E+07},
            { 6,    12,    -2.68479745975134037435E+07},
            { 6,    13,     2.40770592574207345024E+06},
            { 6,    14,     2.26369602661154372618E+06},
            { 6,    15,    -1.54705866766018536873E+06},
            { 6,    16,     7.80276827808928559534E+05},
            { 6,    17,    -2.90525562987123674247E+05},
            { 6,    18,     6.73388625471522245789E+04},
            { 6,    19,    -8.55506609125114846393E+03},
            { 6,    20,     4.57366957264935876992E+02},

            { 8,     0,     6.01887171994113349915E+09},
            { 8,     1,    -3.59631408617360916138E+10},
            { 8,     2,     9.21373711482781219482E+10},
            { 8,     3,    -1.30859993208645553589E+11},
            { 8,     4,     1.09869343601253295898E+11},
            { 8,     5,    -5.19391019925100326538E+10},
            { 8,     6,     9.93795400740454673767E+09},
            { 8,     7,     1.17822075295209312439E+09},
            { 8,     8,     2.08547094828250437975E+08},
            { 8,     9,    -7.77891584686658740044E+08},
            { 8,    10,     9.79231656343366652727E+07},
            { 8,    11,     1.46029283089236587286E+08},
            { 8,    12,    -5.63426838521486073732E+07},
            { 8,    13,    -4.34994308839864679612E+04},
            { 8,    14,     2.67452041592313116416E+06},
            { 8,    15,     2.48601794890810997458E+05},
            { 8,    16,    -2.88661329611176857725E+05},
            { 8,    17,     3.44704079574566349038E+04},
            { 8,    18,     4.71065540904905537900E+03},
            { 8,    19,    -1.22791766059891961049E+03},
            { 8,    20,     6.52659303871330394031E+01},

            {10,     0,     4.09839446260637617111E+09},
            {10,     1,    -1.19243749172475032806E+10},
            {10,     2,     6.44260733488419055939E+08},
            {10,     3,     4.01784270278003311157E+10},
            {10,     4,    -7.02370927844391479492E+10},
            {10,     5,     5.53928774156074829102E+10},
            {10,     6,    -1.99839479763696517944E+10},
            {10,     7,     1.21178171327214643359E+08},
            {10,     8,     2.12195873823345303535E+09},
            {10,     9,    -3.04125284478327631950E+08},
            {10,    10,    -1.36583732859478831291E+08},
            {10,    11,     1.00033640716232210398E+07},
            {10,    12,     2.11985158940840736032E+07},
            {10,    13,    -1.41918619870038777590E+07},
            {10,    14,     1.06037309161468278617E+07},
            {10,    15,    -4.77449659335303772241E+06},
            {10,    16,     7.63050677152805146761E+05},
            {10,    17,     1.56993307696408010088E+05},
            {10,    18,    -8.52750815847208868945E+04},
            {10,    19,     1.33620566100046889915E+04},
            {10,    20,    -7.54749310106248799457E+02},

            {12,     0,    -2.77498112127752447128E+09},
            {12,     1,     1.99816442713947868347E+10},
            {12,     2,    -4.77683150712246398926E+10},
            {12,     3,     5.04117814886737747192E+10},
            {12,     4,    -1.67946036119296798706E+10},
            {12,     5,    -1.38693190414022674561E+10},
            {12,     6,     1.56277095143661766052E+10},
            {12,     7,    -4.57458119846564769745E+09},
            {12,     8,    -8.68468246746188163757E+08},
            {12,     9,     7.75119369125685334206E+08},
            {12,    10,    -1.16491136086922392249E+08},
            {12,    11,     1.47421104143672753125E+07},
            {12,    12,    -1.52535673836318179965E+07},
            {12,    13,    -1.70087987937960936688E+06},
            {12,    14,     4.67476437233580090106E+06},
            {12,    15,    -1.17463899375570146367E+06},
            {12,    16,    -1.10623359788225745433E+05},
            {12,    17,     8.96177599209738109494E+04},
            {12,    18,    -1.26272792247867328115E+04},
            {12,    19,     2.57183382652074556063E+02},
            {12,    20,     4.53946488958070730746E+01},

            {14,     0,    -1.19404109194643745422E+10},
            {14,     1,     4.75026388992112045288E+10},
            {14,     2,    -8.02570569421333770752E+10},
            {14,     3,     7.41607800924110260010E+10},
            {14,     4,    -4.05307739822724685669E+10},
            {14,     5,     1.41796423462739467621E+10},
            {14,     6,    -4.56807621318898487091E+09},
            {14,     7,     1.94941044600867557526E+09},
            {14,     8,    -5.25762401084450244904E+08},
            {14,     9,    -4.20587376462206020951E+07},
            {14,    10,     1.47395545866844896227E+07},
            {14,    11,     5.00008831927335187793E+07},
            {14,    12,    -3.55091173088149204850E+07},
            {14,    13,     1.31608539435279052705E+07},
            {14,    14,    -4.59113670743147470057E+06},
            {14,    15,     1.81956762125828582793E+06},
            {14,    16,    -5.94762161201574141160E+05},
            {14,    17,     1.30033906684856090578E+05},
            {14,    18,    -1.75126293747036870627E+04},
            {14,    19,     1.30486252871918713936E+03},
            {14,    20,    -4.01670777875890010478E+01}
        },
        {                                                                               // HM (Z001000)
            { 0,     0,     2.17405636868515895912E+05},
            { 0,     1,     2.06852476160948835313E+07},
            { 0,     2,    -6.14611410296131018549E+06},
            { 0,     3,    -3.71233464909802228212E+07},
            { 0,     4,     2.59646896356047354639E+07},
            { 0,     5,     2.70401430864754803479E+07},
            { 0,     6,    -5.03161212555342614651E+07},
            { 0,     7,     3.58598481040194258094E+07},
            { 0,     8,    -1.78006177931661754847E+07},
            { 0,     9,     9.30532460730251483619E+06},
            { 0,    10,    -5.15848465171259269118E+06},
            { 0,    11,     2.19027753214642358944E+06},
            { 0,    12,    -5.81644839415979455225E+05},
            { 0,    13,     7.99240321203274652362E+04},
            { 0,    14,    -2.51621643905049268142E+03},
            { 0,    15,    -3.81475220523696918917E+02},

            { 1,     0,    -3.06765020921669416130E+07},
            { 1,     1,    -7.19287141298301517963E+07},
            { 1,     2,    -4.77435954108010306954E+07},
            { 1,     3,     4.06090998470959842205E+08},
            { 1,     4,    -4.94869498500309586525E+08},
            { 1,     5,     2.69898549761012136936E+08},
            { 1,     6,    -7.97774313283056169748E+07},
            { 1,     7,     4.26528677017149850726E+07},
            { 1,     8,    -4.68106602932079061866E+07},
            { 1,     9,     3.02614077510188445449E+07},
            { 1,    10,    -1.02789476280537806451E+07},
            { 1,    11,     1.48880611471109278500E+06},
            { 1,    12,     6.37803643199263970018E+04},
            { 1,    13,    -1.54060450772956992296E+04},
            { 1,    14,    -1.20700325670922356949E+04},
            { 1,    15,     2.24656015193673965769E+03},

            { 2,     0,     1.75066912655571192503E+08},
            { 2,     1,     6.56024659875666722655E+07},
            { 2,     2,     1.23625106866149380803E+08},
            { 2,     3,    -8.66427535415050864220E+08},
            { 2,     4,     1.04647074302398228645E+09},
            { 2,     5,    -5.43304198956838369370E+08},
            { 2,     6,     6.62721930982138589025E+07},
            { 2,     7,     6.96526135511116981506E+07},
            { 2,     8,    -3.78209700669894739985E+07},
            { 2,     9,     5.84715094197653792799E+06},
            { 2,    10,     1.18073224759165989235E+06},
            { 2,    11,    -3.39494064604910614435E+05},
            { 2,    12,    -1.30009674988549915724E+05},
            { 2,    13,     2.04151331991039405693E+04},
            { 2,    14,     1.63426397272586054896E+04},
            { 2,    15,    -3.46623073630780800158E+03},

            { 3,     0,    -4.40378749002870142460E+08},
            { 3,     1,     1.90811958049112737179E+08},
            { 3,     2,    -3.46994414631571888924E+08},
            { 3,     3,     1.03839183062096667290E+09},
            { 3,     4,    -1.09191434318746423721E+09},
            { 3,     5,     6.03376176521021485329E+08},
            { 3,     6,    -1.86538522271839350462E+08},
            { 3,     7,     3.01234692302216663957E+07},
            { 3,     8,    -2.60876296034886920825E+06},
            { 3,     9,    -5.91722311695919837803E+05},
            { 3,    10,    -1.04200386924736900255E+05},
            { 3,    11,     4.46675275850050500594E+05},
            { 3,    12,    -9.51964972779181844089E+04},
            { 3,    13,     3.58457144797106002443E+03},
            { 3,    14,    -1.11371310452127327153E+04},
            { 3,    15,     2.85845040236972772618E+03},

            { 4,     0,     5.72344695093964815140E+08},
            { 4,     1,    -5.13427680762670874596E+08},
            { 4,     2,     7.89372994735144257545E+08},
            { 4,     3,    -1.00183048339205336571E+09},
            { 4,     4,     4.85858157195736885071E+08},
            { 4,     5,    -3.68810635419888198376E+07},
            { 4,     6,    -6.77671612100235223770E+07},
            { 4,     7,     3.51244092417839094996E+07},
            { 4,     8,    -4.04172301070910505950E+06},
            { 4,     9,     4.42372389482865459286E+05},
            { 4,    10,    -1.11637046404088055715E+06},
            { 4,    11,     1.32348000868789851665E+05},
            { 4,    12,     5.77673805733840563335E+04},
            { 4,    13,     1.39889883355121310160E+03},
            { 4,    14,     1.84537734558708234545E+03},
            { 4,    15,    -1.29546033045933745598E+03},

            { 5,     0,    -3.11964641996447086334E+08},
            { 5,     1,     3.20795370867759823799E+08},
            { 5,     2,    -9.73902886296924829483E+08},
            { 5,     3,     1.11109912862654781342E+09},
            { 5,     4,    -3.31644242757628619671E+08},
            { 5,     5,    -7.20401027258160263300E+07},
            { 5,     6,     5.49449315561505854130E+07},
            { 5,     7,    -1.50211363970660679042E+07},
            { 5,     8,     4.37729396379113662988E+05},
            { 5,     9,     1.78953864822479896247E+06},
            { 5,    10,    -9.64483419783796125557E+04},
            { 5,    11,    -1.03372025982109145843E+05},
            { 5,    12,     2.27654966719349886262E+03},
            { 5,    13,    -9.21556404868084791815E+03},
            { 5,    14,     2.29694125905770533791E+03},
            { 5,    15,     2.69172643454379056038E+02},

            { 6,     0,    -1.47363671236134380102E+08},
            { 6,     1,     3.98069493547490179539E+08},
            { 6,     2,     3.80600399916491806507E+08},
            { 6,     3,    -7.77899828440285444260E+08},
            { 6,     4,     2.13971953739619493484E+08},
            { 6,     5,     7.50766140657877773046E+07},
            { 6,     6,    -1.24408092643524613231E+07},
            { 6,     7,    -1.19486997694680951536E+07},
            { 6,     8,     1.72403185042079538107E+06},
            { 6,     9,     5.82624710496194384177E+04},
            { 6,    10,     6.62807449280618020566E+04},
            { 6,    11,    -5.73858576309664931614E+03},
            { 6,    12,     7.87238279678164053621E+03},
            { 6,    13,     5.58027583166682688898E+02},
            { 6,    14,    -7.80191280094344847384E+02},
            { 6,    15,    -4.41524171396966949033E+01},

            { 7,     0,     3.70057909980202019215E+08},
            { 7,     1,    -8.43582286640240669250E+08},
            { 7,     2,     3.06515664902774870396E+08},
            { 7,     3,     2.79680007567675352097E+08},
            { 7,     4,    -1.31316622141129702330E+08},
            { 7,     5,    -4.78340494412834197283E+07},
            { 7,     6,     2.05878316230771690607E+07},
            { 7,     7,     4.84660386411372851580E+06},
            { 7,     8,    -1.28697308684668410569E+06},
            { 7,     9,    -4.26216149967036908492E+05},
            { 7,    10,     1.35944377558333886554E+05},
            { 7,    11,    -2.78481544974478456425E+04},
            { 7,    12,     6.38116677538639805789E+02},
            { 7,    13,     1.86324844948431268676E+03},
            { 7,    14,    -3.54990721623104548144E+02},
            { 7,    15,     4.62881801006983621960E+01},

            { 8,     0,    -2.73364609220364511013E+08},
            { 8,     1,     6.70469859252789497375E+08},
            { 8,     2,    -4.83480623472343206406E+08},
            { 8,     3,     4.03385611767915189266E+07},
            { 8,     4,     7.35016723277738541365E+07},
            { 8,     5,    -1.62650907308186870068E+07},
            { 8,     6,    -2.72791540988075407222E+06},
            { 8,     7,     3.87529004919559811242E+04},
            { 8,     8,     2.21605940650244883727E+05},
            { 8,     9,    -4.99119267431606931495E+04},
            { 8,    10,     3.88533899694829306100E+04},
            { 8,    11,    -4.12454335026239914441E+03},
            { 8,    12,    -5.45806706724830746680E+02},
            { 8,    13,    -7.47011879871399173680E+02},
            { 8,    14,     2.58488648387669854856E+02},
            { 8,    15,    -2.59121452073608971034E+01},

            { 9,     0,     1.00383462413229718804E+08},
            { 9,     1,    -2.74271434515541374683E+08},
            { 9,     2,     2.71244504862839102745E+08},
            { 9,     3,    -1.11115631573666274548E+08},
            { 9,     4,     4.91877530747146159410E+06},
            { 9,     5,     1.25292149804861973971E+07},
            { 9,     6,    -5.59249804659485630691E+06},
            { 9,     7,     1.28330904890557797626E+06},
            { 9,     8,    -1.04317198392582082306E+05},
            { 9,     9,    -9.62038930517230619444E+03},
            { 9,    10,    -3.97809489634518195089E+03},
            { 9,    11,    -5.57959747341926686204E+02},
            { 9,    12,     8.89096577259679747840E+02},
            { 9,    13,    -5.60425945702165773099E+00},
            { 9,    14,    -4.78179067842277234490E+01},
            { 9,    15,     5.74915010119235514452E+00},

            {10,     0,    -1.55431711502586677670E+07},
            {10,     1,     4.79521485383154302835E+07},
            {10,     2,    -6.02665872499209195375E+07},
            {10,     3,     4.05979536634513586760E+07},
            {10,     4,    -1.60619130310957916081E+07},
            {10,     5,     3.70530953097213804722E+06},
            {10,     6,    -3.78917655865714070387E+05},
            {10,     7,    -1.84983833769938682963E+04},
            {10,     8,    -6.31914341540330042335E+03},
            {10,     9,     7.89696836505760438740E+03},
            {10,    10,    -2.17729950811668868482E+03},
            {10,    11,     7.29602008514560338881E+02},
            {10,    12,    -2.38483209776185589135E+02},
            {10,    13,     2.80674740177570534172E+01},
            {10,    14,     1.83297935389278832119E+00},
            {10,    15,    -4.35700622586662367208E-01}
        },
        {                                                                               // RECOM (Z001000)
            { 0,     0,     1.29690611505892015032E+01},
            { 0,     1,     1.15133146967729094179E+00},
            { 0,     2,    -2.03151799780861574973E+00},
            { 0,     3,     1.18911797645565364689E+00},
            { 0,     4,    -2.35730192355203121979E-01},

            { 1,     0,     6.77126934573034411358E-01},
            { 1,     1,     8.12545635631524709730E-01},
            { 1,     2,     1.02571355665047869721E+00},
            { 1,     3,    -1.20203984384456630252E+00},
            { 1,     4,     3.31810706642633568286E-01},

            { 2,     0,    -2.98789230039933073613E+00},
            { 2,     1,     4.40181013672451193486E+00},
            { 2,     2,    -3.18890331653727177041E+00},
            { 2,     3,     8.83079598994844006121E-01},
            { 2,     4,    -1.51558396747922291548E-01},

            { 3,     0,     6.22545012651770335310E+00},
            { 3,     1,    -1.26234063625864720848E+01},
            { 3,     2,     8.89677911090983108977E+00},
            { 3,     3,    -2.01593702143950714856E+00},
            { 3,     4,     1.24221816733169190816E-01},

            { 4,     0,    -3.48677279439661846894E+00},
            { 4,     1,     8.37556424986163960966E+00},
            { 4,     2,    -6.53850783620592235224E+00},
            { 4,     3,     1.69898655557246569536E+00},
            { 4,     4,    -1.21208024376684805890E-01},

            { 5,     0,     4.04481505949380593101E-01},
            { 5,     1,    -1.37770027286069152161E+00},
            { 5,     2,     1.26623104358849225548E+00},
            { 5,     3,    -3.70847513213560797674E-01},
            { 5,     4,     3.04915286309964846112E-02}
        }
    },
    {                                                                                   // Metallicity Z001500 (0.01500)
        {                                                                               // LMR1 (Z001500)
            { 0,     0,     1.60726718444739766767E+01},
            { 0,     1,    -6.24715135137397137299E+00},
            { 0,     2,     1.07672344037789748938E+01},
            { 0,     3,    -7.24052976085269506257E+00},
            { 0,     4,     1.46349058730338077439E+00},

            { 1,     0,    -3.63506353750313060402E+00},
            { 1,     1,     5.36864271459196729097E+01},
            { 1,     2,    -1.16942397088811162575E+02},
            { 1,     3,     9.34035699421548457622E+01},
            { 1,     4,    -2.50874003951161199666E+01},

            { 2,     0,    -7.91754787628987521941E+00},
            { 2,     1,    -5.94653927755674089894E+01},
            { 2,     2,     2.05203056353080341978E+02},
            { 2,     3,    -1.96558176123197796414E+02},
            { 2,     4,     5.97010203155139578257E+01},

            { 3,     0,     2.27278594676717915490E+01},
            { 3,     1,    -2.40820672531560520113E+00},
            { 3,     2,    -1.20196505857006002316E+02},
            { 3,     3,     1.49881801286317909216E+02},
            { 3,     4,    -5.13533879653931180087E+01},

            { 4,     0,    -1.08895531675941441563E+01},
            { 4,     1,     1.55183170137950732226E+01},
            { 4,     2,     1.95899488178189216114E+01},
            { 4,     3,    -3.88388622741310385322E+01},
            { 4,     4,     1.51805521241851977265E+01}
        },
        {                                                                               // LMR2 (Z001500)
            { 0,     0,     1.57661900292269407942E+01},
            { 0,     1,    -3.86739001681062477545E+00},
            { 0,     2,     4.90951789020641982120E+00},
            { 0,     3,    -3.83335877824608228792E+00},
            { 0,     4,     1.46458915608047957058E+00},
            { 0,     5,    -2.16702239649122857523E-01},

            { 1,     0,     8.37711942202392734202E+00},
            { 1,     1,    -2.70768859464224540545E+01},
            { 1,     2,     3.93888011339611523454E+01},
            { 1,     3,    -2.70047525691479570753E+01},
            { 1,     4,     8.70098442072441713435E+00},
            { 1,     5,    -1.06476013712186490245E+00},

            { 2,     0,     3.06826409038962175657E+01},
            { 2,     1,    -1.09405785435920606119E+02},
            { 2,     2,     1.81370175730107348500E+02},
            { 2,     3,    -1.81704582752789093547E+02},
            { 2,     4,     9.02892160295599097708E+01},
            { 2,     5,    -1.65129464015789579889E+01},

            { 3,     0,    -3.62814204845711003600E+02},
            { 3,     1,     1.08689918143087788849E+03},
            { 3,     2,    -1.31054018838510614842E+03},
            { 3,     3,     9.90790478329023699189E+02},
            { 3,     4,    -4.33419474916350452531E+02},
            { 3,     5,     7.62679478243772734913E+01},

            { 4,     0,     1.40418643511127538659E+03},
            { 4,     1,    -3.88233656649431759433E+03},
            { 4,     2,     3.84683355224100250780E+03},
            { 4,     3,    -2.03701250178350619535E+03},
            { 4,     4,     6.48142655530930028362E+02},
            { 4,     5,    -9.67319650941358304408E+01},

            { 5,     0,    -1.22407686126082307965E+03},
            { 5,     1,     3.76978222294538090864E+03},
            { 5,     2,    -3.95869928580390433126E+03},
            { 5,     3,     1.99118186826012993151E+03},
            { 5,     4,    -5.21844465216713615519E+02},
            { 5,     5,     6.04450098185751301116E+01}
        },
        {                                                                               // LMA (Z001500)
            { 0,     0,    -5.15111813953318560380E+04},
            { 0,     1,    -1.38581324034049781039E+06},
            { 0,     2,     1.39308014554255567491E+07},
            { 0,     3,    -5.71603108284442424774E+07},
            { 0,     4,     1.37981165851579815149E+08},
            { 0,     5,    -2.21159494412106066942E+08},
            { 0,     6,     2.47918615635263592005E+08},
            { 0,     7,    -1.97267821930349111557E+08},
            { 0,     8,     1.08701847878496184945E+08},
            { 0,     9,    -3.66657088282169476151E+07},
            { 0,    10,     2.47653496909504802898E+06},
            { 0,    11,     5.12285355178847070783E+06},
            { 0,    12,    -3.49757309776654699817E+06},
            { 0,    13,     1.39230419173041009344E+06},
            { 0,    14,    -4.29369586936748528387E+05},
            { 0,    15,     1.17250418664394310326E+05},
            { 0,    16,    -2.87188898192620908958E+04},
            { 0,    17,     5.69885720592272627982E+03},
            { 0,    18,    -8.01726305735812047715E+02},
            { 0,    19,     6.87272631810855187950E+01},
            { 0,    20,    -2.67041856670744914837E+00},

            { 2,     0,    -1.45339113482951931655E+07},
            { 2,     1,     1.12667823825035437942E+08},
            { 2,     2,    -4.17006587079701125622E+08},
            { 2,     3,     9.85341220899901628494E+08},
            { 2,     4,    -1.66272967677597117424E+09},
            { 2,     5,     2.10006772395531940460E+09},
            { 2,     6,    -2.00687866295555543900E+09},
            { 2,     7,     1.43029664485314631462E+09},
            { 2,     8,    -7.29760192023199796677E+08},
            { 2,     9,     2.43147726715557068586E+08},
            { 2,    10,    -4.05304287569378465414E+07},
            { 2,    11,    -4.01306759306349034887E+05},
            { 2,    12,    -7.32600274554128758609E+05},
            { 2,    13,     1.17855480548669607379E+06},
            { 2,    14,     2.50996201334879733622E+05},
            { 2,    15,    -6.29675976741036050953E+05},
            { 2,    16,     3.31786518519565521274E+05},
            { 2,    17,    -9.37220926397035946138E+04},
            { 2,    18,     1.56585215689283104439E+04},
            { 2,    19,    -1.46922968079042675527E+03},
            { 2,    20,     6.01559001927688683509E+01},

            { 4,     0,     2.33350108446424514055E+08},
            { 4,     1,    -1.60289968489695382118E+09},
            { 4,     2,     4.81990833858025169373E+09},
            { 4,     3,    -8.16360325217743110657E+09},
            { 4,     4,     8.14584602763201999664E+09},
            { 4,     5,    -4.07331865929028415680E+09},
            { 4,     6,    -5.06252506179879188538E+08},
            { 4,     7,     2.29763469151120567322E+09},
            { 4,     8,    -1.62128710050702285767E+09},
            { 4,     9,     5.32701252654437601566E+08},
            { 4,    10,    -4.48733723450238630176E+07},
            { 4,    11,    -2.11723976412126012146E+07},
            { 4,    12,     2.32885320360103854910E+06},
            { 4,    13,     2.03951553722685994580E+06},
            { 4,    14,    -3.13630035119212174322E+05},
            { 4,    15,    -7.05223482959375251085E+04},
            { 4,    16,    -5.37405402453973874799E+04},
            { 4,    17,     4.94362966298164537875E+04},
            { 4,    18,    -1.40518968918354457855E+04},
            { 4,    19,     1.83732556677333059270E+03},
            { 4,    20,    -9.50453907494923839749E+01},

            { 6,     0,    -1.76169219592734664679E+08},
            { 6,     1,     1.42698092226201272011E+09},
            { 6,     2,    -5.28067056308736801147E+09},
            { 6,     3,     1.15336708877598838806E+10},
            { 6,     4,    -1.61229222116397609711E+10},
            { 6,     5,     1.45899208559092102051E+10},
            { 6,     6,    -8.01445678044835853577E+09},
            { 6,     7,     1.90927120293207097054E+09},
            { 6,     8,     4.88043666148680567741E+08},
            { 6,     9,    -3.87662269468736469746E+08},
            { 6,    10,    -3.33959551647167652845E+07},
            { 6,    11,     9.24862097769977152348E+07},
            { 6,    12,    -2.25346824999888129532E+07},
            { 6,    13,    -5.86599007108959276229E+06},
            { 6,    14,     3.70074268507561367005E+06},
            { 6,    15,    -2.56417061140370758949E+05},
            { 6,    16,    -2.29390710704398312373E+05},
            { 6,    17,     6.94915698487437912263E+04},
            { 6,    18,    -7.11327988466724036698E+03},
            { 6,    19,     3.57413795845535204876E+01},
            { 6,    20,     2.80444597204764924925E+01},

            { 8,     0,     1.24245267843930959702E+09},
            { 8,     1,    -6.44322780962794589996E+09},
            { 8,     2,     1.38199616798107013702E+10},
            { 8,     3,    -1.50973626273141441345E+10},
            { 8,     4,     7.46748678501209640503E+09},
            { 8,     5,     8.19630542308683276176E+08},
            { 8,     6,    -3.01758840394945383072E+09},
            { 8,     7,     1.45032130975556683540E+09},
            { 8,     8,    -3.13138358723539471626E+08},
            { 8,     9,     1.84361824896876424551E+08},
            { 8,    10,    -1.69216581016204208136E+08},
            { 8,    11,     5.96086349759788066149E+07},
            { 8,    12,     3.41292303451553406194E+06},
            { 8,    13,    -8.38948524171856231987E+06},
            { 8,    14,     2.24497033832136029378E+06},
            { 8,    15,    -2.11810850483575311955E+05},
            { 8,    16,     5.84158368683356893598E+04},
            { 8,    17,    -4.11205502318392173038E+04},
            { 8,    18,     1.17356696779575449909E+04},
            { 8,    19,    -1.52952212651557533718E+03},
            { 8,    20,     7.85072147480235713601E+01},

            {10,     0,    -1.81932512063875246048E+09},
            {10,     1,     1.17609472824396915436E+10},
            {10,     2,    -3.26657704396900253296E+10},
            {10,     3,     4.99244642700828399658E+10},
            {10,     4,    -4.43242869808500061035E+10},
            {10,     5,     2.08982330227992744446E+10},
            {10,     6,    -2.02708169583341908455E+09},
            {10,     7,    -2.88544961080271482468E+09},
            {10,     8,     1.12173715222722339630E+09},
            {10,     9,     1.71997591831600487232E+08},
            {10,    10,    -1.91453234115329504013E+08},
            {10,    11,     3.25976767552792169154E+07},
            {10,    12,    -2.81013549671907734592E+05},
            {10,    13,     1.62781335264503862709E+06},
            {10,    14,    -1.02482288018334569642E+05},
            {10,    15,    -4.78252977657809446100E+05},
            {10,    16,     1.08856413307630777126E+05},
            {10,    17,     3.74121964994652880705E+04},
            {10,    18,    -2.00568317658710147953E+04},
            {10,    19,     3.30817698426615561402E+03},
            {10,    20,    -1.97037754182892768995E+02},

            {12,     0,     3.35585658860964119434E+08},
            {12,     1,     6.04561877929950594902E+08},
            {12,     2,    -5.13752404817218589783E+09},
            {12,     3,     1.06028901045264606476E+10},
            {12,     4,    -1.25273394138480796814E+10},
            {12,     5,     1.06453110355505924225E+10},
            {12,     6,    -6.95095817367120552063E+09},
            {12,     7,     3.16043713721437406540E+09},
            {12,     8,    -7.21320971392456412315E+08},
            {12,     9,    -5.26415434883427247405E+07},
            {12,    10,     3.88067712231941595674E+07},
            {12,    11,     2.20965666929983980954E+07},
            {12,    12,    -1.00345439357733502984E+07},
            {12,    13,    -2.83703418567746179178E+06},
            {12,    14,     2.26345678725755680352E+06},
            {12,    15,    -2.51391728097161831101E+05},
            {12,    16,    -1.16427576675173724652E+05},
            {12,    17,     3.66684894778185043833E+04},
            {12,    18,    -2.41298286705503824123E+03},
            {12,    19,    -3.40666737232470040908E+02},
            {12,    20,     4.25259705893395718590E+01},

            {14,     0,    -1.88270337231409168243E+09},
            {14,     1,     6.39532526253386974335E+09},
            {14,     2,    -7.82324355639520549774E+09},
            {14,     3,     2.03781762196393775940E+09},
            {14,     4,     5.26969578325209617615E+09},
            {14,     5,    -7.26948054175921058655E+09},
            {14,     6,     4.65417136459286022186E+09},
            {14,     7,    -1.63391954770838069916E+09},
            {14,     8,     1.87346108643729388714E+08},
            {14,     9,     8.89076657232203483582E+07},
            {14,    10,    -2.42219275225403495133E+07},
            {14,    11,    -1.37909139873882699758E+07},
            {14,    12,     6.82268681301416642964E+06},
            {14,    13,     7.02051292148833046667E+05},
            {14,    14,    -6.92968335637202719226E+05},
            {14,    15,    -3.38601851769478351343E+05},
            {14,    16,     3.62546876010872249026E+05},
            {14,    17,    -1.30131463530392415123E+05},
            {14,    18,     2.47381054441249470983E+04},
            {14,    19,    -2.51595863192902470473E+03},
            {14,    20,     1.08731772056314042629E+02}
        },
        {                                                                               // HM (Z001500)
            { 0,     0,    -6.40222930769867496565E+05},
            { 0,     1,     1.10169653380964145064E+08},
            { 0,     2,    -1.83288798818755358458E+08},
            { 0,     3,     2.34685038301938652992E+08},
            { 0,     4,    -2.90868321715029954910E+08},
            { 0,     5,     2.58803902657175153494E+08},
            { 0,     6,    -1.63288883819993436337E+08},
            { 0,     7,     8.19842075092330276966E+07},
            { 0,     8,    -3.74597286079473644495E+07},
            { 0,     9,     1.73063130301320850849E+07},
            { 0,    10,    -7.45417360680353827775E+06},
            { 0,    11,     2.34362116585996653885E+06},
            { 0,    12,    -3.93361441283391148318E+05},
            { 0,    13,    -1.10781031619090754248E+02},
            { 0,    14,     1.09339428886468467681E+04},
            { 0,    15,    -1.22318915545462755290E+03},

            { 1,     0,    -2.80689302632132731378E+06},
            { 1,     1,    -8.98291909237491011620E+08},
            { 1,     2,     1.18592188686839365959E+09},
            { 1,     3,    -7.65219587451042652130E+08},
            { 1,     4,     5.70295193113373160362E+08},
            { 1,     5,    -3.92352877095889389515E+08},
            { 1,     6,     1.43709948555291593075E+08},
            { 1,     7,     3.04550315389612223953E+06},
            { 1,     8,    -3.80446096340993121266E+07},
            { 1,     9,     2.31604841598522365093E+07},
            { 1,    10,    -5.65692487392537854612E+06},
            { 1,    11,    -3.33505242560792976292E+03},
            { 1,    12,     1.18223884409625112312E+05},
            { 1,    13,     8.29620922556901932694E+04},
            { 1,    14,    -3.61397359909771475941E+04},
            { 1,    15,     3.94064567463254252289E+03},

            { 2,     0,     2.01327877490426152945E+08},
            { 2,     1,     3.10118016671073770523E+09},
            { 2,     2,    -4.22208827383246850967E+09},
            { 2,     3,     1.99760734988433146477E+09},
            { 2,     4,    -7.30736768670006632805E+08},
            { 2,     5,     4.63059271373849034309E+08},
            { 2,     6,    -2.61651098515476316214E+08},
            { 2,     7,     1.29333847253416240215E+08},
            { 2,     8,    -5.09144396244669929147E+07},
            { 2,     9,     9.08995892725536786020E+06},
            { 2,    10,    -4.72198901409551501274E+05},
            { 2,    11,     8.76993455861743539572E+05},
            { 2,    12,    -4.06936925562300020829E+05},
            { 2,    13,    -6.55970599629438947886E+03},
            { 2,    14,     3.07809196263461344643E+04},
            { 2,    15,    -4.44326151571276750474E+03},

            { 3,     0,    -1.20734569145530676842E+09},
            { 3,     1,    -5.34256358175416564941E+09},
            { 3,     2,     8.35253308946118927002E+09},
            { 3,     3,    -3.55132392736585330963E+09},
            { 3,     4,     3.64000813302139878273E+08},
            { 3,     5,     6.22636854320012405515E+07},
            { 3,     6,    -7.19534735248802900314E+07},
            { 3,     7,     3.05025745163594335318E+07},
            { 3,     8,     7.68708163315938971937E+06},
            { 3,     9,    -3.07800794374816305935E+06},
            { 3,    10,    -2.14286733640326838940E+06},
            { 3,    11,     9.30921378395954845473E+05},
            { 3,    12,    -5.63446497619864894659E+04},
            { 3,    13,     1.63174566022897852235E+03},
            { 3,    14,    -1.16542300609449248441E+04},
            { 3,    15,     2.36843518043800349915E+03},

            { 4,     0,     3.21429576530229854584E+09},
            { 4,     1,     4.39998454027949714661E+09},
            { 4,     2,    -9.63131293606875228882E+09},
            { 4,     3,     4.22524392135507345200E+09},
            { 4,     4,    -1.44018430614179670811E+08},
            { 4,     5,    -1.74848804728251099586E+08},
            { 4,     6,     3.49971447999607846141E+07},
            { 4,     7,    -2.03914002418685518205E+07},
            { 4,     8,    -1.16793722923461743630E+06},
            { 4,     9,     4.07216598363398248330E+06},
            { 4,    10,    -6.63453501301569864154E+05},
            { 4,    11,    -2.07226327457722654799E+05},
            { 4,    12,     7.70487153054169175448E+04},
            { 4,    13,     3.74511172163199034912E+03},
            { 4,    14,    -2.35924085496946099738E+03},
            { 4,    15,    -2.00074150360817185401E+02},

            { 5,     0,    -4.86051768030780982971E+09},
            { 5,     1,    -2.52283411110093832016E+08},
            { 5,     2,     6.15744921080406188965E+09},
            { 5,     3,    -2.94612427147580671310E+09},
            { 5,     4,    -1.64832248763824313879E+08},
            { 5,     5,     2.56579223709281235933E+08},
            { 5,     6,    -1.66799195149700418115E+07},
            { 5,     7,     1.29958687433429018711E+05},
            { 5,     8,    -1.39862347051455453038E+06},
            { 5,     9,     6.31508472462175530382E+05},
            { 5,    10,     1.11229199938330348232E+04},
            { 5,    11,    -1.09408809112410061061E+05},
            { 5,    12,    -3.31703747407727632890E+03},
            { 5,    13,    -7.42683722320846641196E+02},
            { 5,    14,     3.67269071514491452035E+03},
            { 5,    15,    -4.51573620652630495442E+02},

            { 6,     0,     4.54118409271736431122E+09},
            { 6,     1,    -2.90767159976795244217E+09},
            { 6,     2,    -1.58025838498863267899E+09},
            { 6,     3,     1.16398375723541212082E+09},
            { 6,     4,     1.09676718899973109365E+08},
            { 6,     5,    -1.17450609424565017223E+08},
            { 6,     6,     3.66726234801336098462E+06},
            { 6,     7,     1.73320626005022414029E+06},
            { 6,     8,    -2.36895981676343455911E+06},
            { 6,     9,     6.50776589120257762261E+05},
            { 6,    10,     9.80770387879816262284E+04},
            { 6,    11,     6.08450624306344252545E+03},
            { 6,    12,     1.44103835104128138482E+04},
            { 6,    13,    -8.31901750284624540654E+03},
            { 6,    14,    -3.29760727015796192063E+02},
            { 6,    15,     2.22378172851269994226E+02},

            { 7,     0,    -2.66951759357458686829E+09},
            { 7,     1,     2.78244247697885799408E+09},
            { 7,     2,    -5.23326483537968039513E+08},
            { 7,     3,    -1.32406331732717037201E+08},
            { 7,     4,    -7.73866632845929916948E+06},
            { 7,     5,    -1.31289335217684842646E+07},
            { 7,     6,     8.68883711662754230201E+06},
            { 7,     7,     6.03206328956170473248E+06},
            { 7,     8,    -1.68968455437217419967E+06},
            { 7,     9,    -1.77911507059704919811E+05},
            { 7,    10,    -2.13170445228394455626E+04},
            { 7,    11,    -1.98704865429771484742E+03},
            { 7,    12,     1.19148875731268594791E+03},
            { 7,    13,     2.79295327351346804790E+03},
            { 7,    14,    -2.78996535146859685028E+02},
            { 7,    15,    -4.49988682514610260910E+01},

            { 8,     0,     9.50171459722543239594E+08},
            { 8,     1,    -1.19109509072384023666E+09},
            { 8,     2,     4.45397440506503641605E+08},
            { 8,     3,    -4.46003071982323154807E+07},
            { 8,     4,    -3.80750898013454861939E+06},
            { 8,     5,     1.37090059550021272153E+07},
            { 8,     6,    -5.08256188045835960656E+06},
            { 8,     7,    -1.65700011583761894144E+06},
            { 8,     8,     4.63444510281368449796E+05},
            { 8,     9,     6.44070571055476248148E+04},
            { 8,    10,     5.85348885756574745756E+04},
            { 8,    11,    -2.07048623592222938896E+04},
            { 8,    12,    -1.63741310696323444063E+03},
            { 8,    13,     6.19764619828629179210E+02},
            { 8,    14,    -5.07010754652976132206E+01},
            { 8,    15,     1.04668911802231079378E+01},

            { 9,     0,    -1.81045244850725114346E+08},
            { 9,     1,     2.23702214978076577187E+08},
            { 9,     2,    -5.39576641499614790082E+07},
            { 9,     3,    -3.48971727291814684868E+07},
            { 9,     4,     2.05198001371338181198E+07},
            { 9,     5,    -1.21495342617135704495E+06},
            { 9,     6,    -3.16323564470732118934E+06},
            { 9,     7,     1.45230267412367230281E+06},
            { 9,     8,     4.31289414334815301117E+04},
            { 9,     9,    -1.48747640901944541838E+05},
            { 9,    10,     8.63897175520341988886E+03},
            { 9,    11,     7.41770223687096313370E+03},
            { 9,    12,    -3.19056811887575008768E+02},
            { 9,    13,    -3.21467405033000488856E+02},
            { 9,    14,     5.29898436774840320140E+01},
            { 9,    15,    -3.45842797011066638291E+00},

            {10,     0,     1.29428707331854440272E+07},
            {10,     1,    -8.40693462242027930915E+06},
            {10,     2,    -1.66950000229333285242E+07},
            {10,     3,     2.57591525679362528026E+07},
            {10,     4,    -1.54869040150702800602E+07},
            {10,     5,     4.65623680094749014825E+06},
            {10,     6,    -3.66297825326654070523E+05},
            {10,     7,    -1.74355945057241042377E+05},
            {10,     8,     1.17860134753315287526E+04},
            {10,     9,     1.88480662367682598415E+04},
            {10,    10,    -1.91529993314948478655E+03},
            {10,    11,    -1.26547072451270787496E+03},
            {10,    12,     2.41500958409452209708E+02},
            {10,    13,     1.39658146862240180042E+01},
            {10,    14,    -6.16658750062571581196E+00},
            {10,    15,     4.30482965248756621612E-01}
        },
        {                                                                               // RECOM (Z001500)
            { 0,     0,     1.30114084816216948326E+01},
            { 0,     1,     9.82927123720802486950E-01},
            { 0,     2,    -1.78831205264574788494E+00},
            { 0,     3,     1.04428283556912560037E+00},
            { 0,     4,    -2.05873948734389328186E-01},

            { 1,     0,    -3.88724524733189344405E-02},
            { 1,     1,     3.79408475484787910403E+00},
            { 1,     2,    -2.91118101352734681697E+00},
            { 1,     3,     7.77250458270365651714E-01},
            { 1,     4,    -5.93043121201380050989E-03},

            { 2,     0,    -1.19760400814100687050E+00},
            { 2,     1,    -3.70516957950696035340E+00},
            { 2,     2,     8.13081389223339456862E+00},
            { 2,     3,    -4.84039188270623732535E+00},
            { 2,     4,     8.01677869943414500575E-01},

            { 3,     0,     5.25028472760498932104E+00},
            { 3,     1,    -6.49412032630474111983E+00},
            { 3,     2,    -1.37460399330354388070E+00},
            { 3,     3,     3.58589212553425529251E+00},
            { 3,     4,    -8.38094449459952040016E-01},

            { 4,     0,    -3.35142128249576209953E+00},
            { 4,     1,     6.88602644973469768530E+00},
            { 4,     2,    -2.99992638566578184722E+00},
            { 4,     3,    -4.88499106441473263107E-01},
            { 4,     4,     2.76990054008695862908E-01},

            { 5,     0,     3.32597771543929221494E-01},
            { 5,     1,    -1.23233281888843992924E+00},
            { 5,     2,     8.33374957715134034864E-01},
            { 5,     3,    -6.82770912855755612858E-02},
            { 5,     4,    -2.81426477843974842674E-02}
        }
    },
    {                                                                                   // Metallicity Z002000 (0.02000)
        {                                                                               // LMR1 (Z002000)
            { 0,     0,     1.50305864307172658556E+01},
            { 0,     1,     4.96110465703095582235E-01},
            { 0,     2,    -9.16608282566077403608E-01},
            { 0,     3,     2.48033316991016244968E-01},

            { 1,     0,     1.77238406649739155263E+00},
            { 1,     1,    -6.52669943496246296455E-01},
            { 1,     2,     6.23616025599144307989E-01},
            { 1,     3,    -1.77703681033998805994E-01}
        },
        {                                                                               // LMR2 (Z002000)
            { 0,     0,     1.56617552673034907684E+01},
            { 0,     1,    -3.38271409511212528543E+00},
            { 0,     2,     3.99712000854973803499E+00},
            { 0,     3,    -3.02506056136319223526E+00},
            { 0,     4,     1.12651646333780397491E+00},
            { 0,     5,    -1.62914377442334834534E-01},

            { 1,     0,     7.36910246144982838956E+00},
            { 1,     1,    -2.30880577686855161801E+01},
            { 1,     2,     3.32294193400783512971E+01},
            { 1,     3,    -2.22327077943441828722E+01},
            { 1,     4,     6.85986578867425578210E+00},
            { 1,     5,    -7.85596865167521363205E-01},

            { 2,     0,     3.04288311387925105578E+01},
            { 2,     1,    -1.38697534163042433875E+02},
            { 2,     2,     2.72764619820308155340E+02},
            { 2,     3,    -2.82250602598139892052E+02},
            { 2,     4,     1.37077309352478550863E+02},
            { 2,     5,    -2.43363801583426564434E+01},

            { 3,     0,    -2.40391883810966504598E+02},
            { 3,     1,     7.20123937632566253342E+02},
            { 3,     2,    -9.56047390568349555906E+02},
            { 3,     3,     8.80126339460319968566E+02},
            { 3,     4,    -4.43060386994394036719E+02},
            { 3,     5,     8.34534841118213392974E+01},

            { 4,     0,     1.19040215556355087756E+03},
            { 4,     1,    -3.06320152050599153881E+03},
            { 4,     2,     2.78453473219604893529E+03},
            { 4,     3,    -1.43582241909097660937E+03},
            { 4,     4,     4.99720580375658244066E+02},
            { 4,     5,    -8.44179485778146130315E+01},

            { 5,     0,    -1.27279708339606577283E+03},
            { 5,     1,     3.68563309227265335721E+03},
            { 5,     2,    -3.66200350857788362191E+03},
            { 5,     3,     1.74828529177925724980E+03},
            { 5,     4,    -4.43446115665101842751E+02},
            { 5,     5,     5.16038082522078980219E+01}
        },
        {                                                                               // LMA (Z002000)
            { 0,     0,    -1.51049254914130983707E+03},
            { 0,     1,     5.90804592923560267081E+03},
            { 0,     2,    -9.78411552163013584504E+03},
            { 0,     3,     9.03281646038073631644E+03},
            { 0,     4,    -5.04733341717106668511E+03},
            { 0,     5,     1.70840439977420237483E+03},
            { 0,     6,    -3.16251664756395882705E+02},
            { 0,     7,     1.74202707579704139107E+01},
            { 0,     8,     3.67713293499369164863E+00},
            { 0,     9,    -5.09055045574688058707E-01},

            { 1,     0,     1.18642296961242536781E+05},
            { 1,     1,    -4.22882345998080738354E+05},
            { 1,     2,     6.02728568485423107632E+05},
            { 1,     3,    -4.09618748060796875507E+05},
            { 1,     4,     8.90979215331482701004E+04},
            { 1,     5,     5.95414095173196401447E+04},
            { 1,     6,    -5.18590927506359803374E+04},
            { 1,     7,     1.72331732723822569824E+04},
            { 1,     8,    -2.81277261133411593619E+03},
            { 1,     9,     1.86487606226860776815E+02},

            { 2,     0,    -1.40854420295521779917E+06},
            { 2,     1,     3.16920022432578727603E+06},
            { 2,     2,     4.61300052596119523514E+05},
            { 2,     3,    -8.61105917229603044689E+06},
            { 2,     4,     1.23556398056155946106E+07},
            { 2,     5,    -8.94036121833491511643E+06},
            { 2,     6,     3.81112897666227677837E+06},
            { 2,     7,    -9.71158692551792715676E+05},
            { 2,     8,     1.37401910653425060445E+05},
            { 2,     9,    -8.32876128728674302693E+03},

            { 3,     0,     1.72458181532576158643E+07},
            { 3,     1,    -4.52772539129179418087E+07},
            { 3,     2,     2.01148809934781342745E+07},
            { 3,     3,     6.09227038679622411728E+07},
            { 3,     4,    -1.07679655401528432965E+08},
            { 3,     5,     8.29801717072715312243E+07},
            { 3,     6,    -3.64493890278426855803E+07},
            { 3,     7,     9.44902866139182634652E+06},
            { 3,     8,    -1.35153547941456316039E+06},
            { 3,     9,     8.25342644314703647979E+04},

            { 4,     0,    -1.36692876826326847076E+08},
            { 4,     1,     4.51025336598177015781E+08},
            { 4,     2,    -5.31011502233843207359E+08},
            { 4,     3,     1.64806847028889596462E+08},
            { 4,     4,     2.15493288928440541029E+08},
            { 4,     5,    -2.68977365806154668331E+08},
            { 4,     6,     1.38776598747677832842E+08},
            { 4,     7,    -3.89913542043956518173E+07},
            { 4,     8,     5.84939274544744472951E+06},
            { 4,     9,    -3.68464314099946233910E+05},

            { 5,     0,     5.88007632627574801445E+08},
            { 5,     1,    -2.19546849066643857956E+09},
            { 5,     2,     3.29578755621284484863E+09},
            { 5,     3,    -2.44640117334454774857E+09},
            { 5,     4,     7.63562415504432678223E+08},
            { 5,     5,     1.37457055430822163820E+08},
            { 5,     6,    -2.06529556624937295914E+08},
            { 5,     7,     7.50013034409245699644E+07},
            { 5,     8,    -1.26735707259846013039E+07},
            { 5,     9,     8.54653687692407402210E+05},

            { 6,     0,    -1.38327254361235642433E+09},
            { 6,     1,     5.60614067713224506378E+09},
            { 6,     2,    -9.43342679709010314941E+09},
            { 6,     3,     8.54102479910387897491E+09},
            { 6,     4,    -4.42476607837071800232E+09},
            { 6,     5,     1.22015302616584563255E+09},
            { 6,     6,    -9.17636537403069436550E+07},
            { 6,     7,    -4.19402333332933112979E+07},
            { 6,     8,     1.20023414749641995877E+07},
            { 6,     9,    -9.83676207032694481313E+05},

            { 7,     0,     1.70433842457238912582E+09},
            { 7,     1,    -7.46564615829306602478E+09},
            { 7,     2,     1.36213117744959812164E+10},
            { 7,     3,    -1.36449069424619369507E+10},
            { 7,     4,     8.23382455328104782104E+09},
            { 7,     5,    -3.05178841976632213593E+09},
            { 7,     6,     6.64362036529520750046E+08},
            { 7,     7,    -7.09507497840798497200E+07},
            { 7,     8,     7.89097510117984027602E+05},
            { 7,     9,     3.47019468185705947690E+05},

            { 8,     0,    -9.03133463246678471565E+08},
            { 8,     1,     4.47621684993862724304E+09},
            { 8,     2,    -8.96349174100218391418E+09},
            { 8,     3,     9.77159178772765159607E+09},
            { 8,     4,    -6.46222307786082077026E+09},
            { 8,     5,     2.69448710911580371857E+09},
            { 8,     6,    -7.04153944276203036308E+08},
            { 8,     7,     1.09216364710150659084E+08},
            { 8,     8,    -8.71251345965249091387E+06},
            { 8,     9,     2.35094345904298679670E+05},

            { 9,     0,     7.47826861721657812595E+07},
            { 9,     1,    -6.76538905537033677101E+08},
            { 9,     2,     1.73090915363686656952E+09},
            { 9,     3,    -2.18473390758292007446E+09},
            { 9,     4,     1.61084727759012627602E+09},
            { 9,     5,    -7.39637661531803250313E+08},
            { 9,     6,     2.13872151698962450027E+08},
            { 9,     7,    -3.76329373731744587421E+07},
            { 9,     8,     3.63438629714747751132E+06},
            { 9,     9,    -1.44132561512217595009E+05}
        },
        {                                                                               // HM (Z002000)
            { 0,    -4,     2.53729226092276048660E+09},
            { 0,    -3,    -9.90007540851090812683E+09},
            { 0,    -2,     1.52593834113249797821E+10},
            { 0,    -1,    -1.01549528384376564026E+10},
            { 0,     0,     2.99829578010063916445E+07},
            { 0,     1,     5.78893490124735832214E+09},
            { 0,     2,    -6.14609898490778446198E+09},
            { 0,     3,     4.04278241615104389191E+09},
            { 0,     4,    -1.80654362590145325661E+09},
            { 0,     5,     4.91334348953023314476E+08},
            { 0,     6,    -5.44061721887778788805E+07},
            { 0,     7,    -6.47383359300961066037E+06},
            { 0,     8,     1.76558650155847985297E+06},
            { 0,     9,     1.00096546886151510989E+05},
            { 0,    10,    -3.39961914804177431506E+04},

            { 1,    -4,    -1.43765508369811935425E+10},
            { 1,    -3,     5.53386381611402206421E+10},
            { 1,    -2,    -8.64972441405090179443E+10},
            { 1,    -1,     6.18263781388091888428E+10},
            { 1,     0,    -1.27031002818969459534E+10},
            { 1,     1,    -1.22741674009294738770E+10},
            { 1,     2,     1.32545543793032627106E+10},
            { 1,     3,    -7.72103742705006504059E+09},
            { 1,     4,     2.99797494889310359955E+09},
            { 1,     5,    -4.51316994491591095924E+08},
            { 1,     6,    -2.12359218253416419029E+08},
            { 1,     7,     1.28106564253871873021E+08},
            { 1,     8,    -2.52927841362378038466E+07},
            { 1,     9,     1.52832625190069107339E+06},
            { 1,    10,     5.33991439763792805024E+04},

            { 2,    -4,     3.53888765741276016235E+10},
            { 2,    -3,    -1.30949004957460281372E+11},
            { 2,    -2,     2.03740123625964263916E+11},
            { 2,    -1,    -1.45485633356229827881E+11},
            { 2,     0,     3.38169397075351562500E+10},
            { 2,     1,     1.75117599864639778137E+10},
            { 2,     2,    -1.68526525631832427979E+10},
            { 2,     3,     7.96671769498775196075E+09},
            { 2,     4,    -3.11421864218404531479E+09},
            { 2,     5,     6.72407024087293863297E+08},
            { 2,     6,     1.69479333713295638561E+08},
            { 2,     7,    -1.48711963495587915182E+08},
            { 2,     8,     3.25634279867792949080E+07},
            { 2,     9,    -1.75731226257729995996E+06},
            { 2,    10,    -1.35079742770175478654E+05},

            { 3,    -4,    -5.02259172277233886719E+10},
            { 3,    -3,     1.70340081645104858398E+11},
            { 3,    -2,    -2.59178892024369842529E+11},
            { 3,    -1,     1.80456668127881774902E+11},
            { 3,     0,    -4.10109717020726623535E+10},
            { 3,     1,    -1.59603847360727577209E+10},
            { 3,     2,     1.16582282107490043640E+10},
            { 3,     3,    -2.69128104888696050644E+09},
            { 3,     4,     7.08683177906285405159E+08},
            { 3,     5,    -4.34911479088057935238E+08},
            { 3,     6,     1.01044743218226939440E+08},
            { 3,     7,     1.50179933915116582066E+07},
            { 3,     8,    -4.86460542015785723925E+06},
            { 3,     9,    -1.18190563420071476139E+06},
            { 3,    10,     2.78623554417385661509E+05},

            { 4,    -4,     4.64705908425810546875E+10},
            { 4,    -3,    -1.30875109087104644775E+11},
            { 4,    -2,     1.84384613263450225830E+11},
            { 4,    -1,    -1.14315528459969711304E+11},
            { 4,     0,     1.54155640604285240173E+10},
            { 4,     1,     1.32641401056199913025E+10},
            { 4,     2,    -5.51530197341053485870E+09},
            { 4,     3,    -1.70800161550766736269E+08},
            { 4,     4,     7.78830644995285868645E+08},
            { 4,     5,    -2.06842292132320821285E+08},
            { 4,     6,    -2.58832894926396012306E+07},
            { 4,     7,     2.22092620407946482301E+07},
            { 4,     8,    -6.95661957240117993206E+06},
            { 4,     9,     1.92690144883910729550E+06},
            { 4,    10,    -2.33277229786988784326E+05},

            { 5,    -4,    -3.03660798334242973328E+10},
            { 5,    -3,     5.92434377201234664917E+10},
            { 5,    -2,    -6.54677950676404571533E+10},
            { 5,    -1,     1.98006473040911102295E+10},
            { 5,     0,     1.80741102193409042358E+10},
            { 5,     1,    -1.24110750642199649811E+10},
            { 5,     2,     2.06936088729867076874E+09},
            { 5,     3,    -5.55200935775747060776E+08},
            { 5,     4,     3.31793995932478666306E+08},
            { 5,     5,    -3.77282162409072965384E+07},
            { 5,     6,     3.77271593664872169029E+05},
            { 5,     7,    -2.39515931142189726233E+06},
            { 5,     8,     4.78976734628147096373E+05},
            { 5,     9,    -2.24492564827198133571E+05},
            { 5,    10,     4.79792411042873136466E+04},

            { 6,    -4,     1.39020977034095058441E+10},
            { 6,    -3,    -1.24214230056445674896E+10},
            { 6,    -2,     3.40159755029327201843E+09},
            { 6,    -1,     1.62527040079761524200E+10},
            { 6,     0,    -2.40572617307258148193E+10},
            { 6,     1,     1.01576765639793376923E+10},
            { 6,     2,    -2.79329356560366511345E+08},
            { 6,     3,    -3.70225353787005186081E+08},
            { 6,     4,    -1.05583829363182082772E+08},
            { 6,     5,     4.22791698590517193079E+07},
            { 6,     6,    -2.37510881991530815139E+06},
            { 6,     7,    -2.48551971484777750447E+06},
            { 6,     8,     2.30980219498265953735E+06},
            { 6,     9,    -5.97702223240277264267E+05},
            { 6,    10,     4.54666582402248095605E+04},

            { 7,    -4,    -2.83252506215792465210E+09},
            { 7,    -3,    -5.87107593619644260406E+09},
            { 7,    -2,     9.35640963237767410278E+09},
            { 7,    -1,    -9.10900080282122421265E+09},
            { 7,     0,     9.58587735224136924744E+09},
            { 7,     1,    -5.49030619645441913605E+09},
            { 7,     2,     1.36275859178348875046E+09},
            { 7,     3,    -8.57005601597409099340E+07},
            { 7,     4,    -7.86904883546570241451E+07},
            { 7,     5,     4.99692720418420732021E+07},
            { 7,     6,    -7.09569033777875266969E+06},
            { 7,     7,    -8.55338083130266284570E+05},
            { 7,     8,    -5.61264059584917849861E+05},
            { 7,     9,     3.08073506627650989685E+05},
            { 7,    10,    -3.28288523737060459098E+04},

            { 8,    -4,    -1.46525129680925726891E+09},
            { 8,    -3,     9.31328032341473007202E+09},
            { 8,    -2,    -1.05576564874556827545E+10},
            { 8,    -1,     4.15304574768818807602E+09},
            { 8,     0,    -5.00505353118861734867E+08},
            { 8,     1,     2.38979281271329969168E+08},
            { 8,     2,    -3.43075090044980287552E+08},
            { 8,     3,     1.53399431073772042990E+08},
            { 8,     4,     6.15683430398282781243E+06},
            { 8,     5,    -2.10229696774827726185E+07},
            { 8,     6,     6.68991617223582346924E+05},
            { 8,     7,     2.39900342987180408090E+06},
            { 8,     8,    -4.43480899599777883850E+05},
            { 8,     9,    -2.42983373542564295349E+04},
            { 8,    10,     7.76394675149368777056E+03},

            { 9,    -4,     1.21311579654507946968E+09},
            { 9,    -3,    -5.24070703990287113190E+09},
            { 9,    -2,     7.06700234139572620392E+09},
            { 9,    -1,    -4.37710070889741706848E+09},
            { 9,     0,     1.25442565132006978989E+09},
            { 9,     1,    -4.98846050957186818123E+07},
            { 9,     2,    -6.82854203030566722155E+07},
            { 9,     3,     3.25620574902522787452E+07},
            { 9,     4,    -2.20218804932029694319E+07},
            { 9,     5,     9.02162431089685671031E+06},
            { 9,     6,    -2.40492887758824275807E+05},
            { 9,     7,    -8.68585806973894243129E+05},
            { 9,     8,     2.25402435254328796873E+05},
            { 9,     9,    -1.55503066541077005240E+04},
            { 9,    10,    -3.94672900362825259890E+02},

            {10,    -4,    -2.52924874342063546181E+08},
            {10,    -3,     1.08433769070899581909E+09},
            {10,    -2,    -1.75274841405356764793E+09},
            {10,    -1,     1.53857146715462470055E+09},
            {10,     0,    -8.55304471302122354507E+08},
            {10,     1,     3.31495823209429860115E+08},
            {10,     2,    -9.56754638305779546499E+07},
            {10,     3,     1.92680795062012895942E+07},
            {10,     4,    -1.01626479920090991072E+06},
            {10,     5,    -6.69255818740246933885E+05},
            {10,     6,     1.56196811545844484499E+04},
            {10,     7,     9.97457630392971332185E+04},
            {10,     8,    -2.92981386761484645831E+04},
            {10,     9,     2.94914409981776043423E+03},
            {10,    10,    -6.54138557526380992613E+01}
        },
        {                                                                               // RECOM (Z002000)
            { 0,     0,     1.34142729153601880654E+01},
            { 0,     1,    -8.33977084929894418863E-01},
            { 0,     2,     6.31803234438807592710E-01},
            { 0,     3,    -1.70217751115103232973E-01},

            { 1,     0,     1.57297096125805579980E+00},
            { 1,     1,     3.04879789211891349954E-01},
            { 1,     2,    -6.10228805816398045536E-01},
            { 1,     3,     2.31090540934666771600E-01},

            { 2,     0,    -1.41108730414326188907E+00},
            { 2,     1,     1.23004856793249794933E+00},
            { 2,     2,    -2.83350980192163703908E-01},
            { 2,     3,    -4.69140139332027000796E-02},

            { 3,     0,     5.62788156340448986192E-01},
            { 3,     1,    -6.89182559641456582433E-01},
            { 3,     2,     2.69502539954297459790E-01},
            { 3,     3,    -2.14046891159626988255E-02}
        }
    },
    {                                                                                   // Metallicity Z003000 (0.03000)
        {                                                                               // LMR1 (Z003000)
            { 0,     0,     1.54866491842634630416E+01},
            { 0,     1,    -2.35824863480211543987E+00},
            { 0,     2,     1.34711604769709447638E+00},
            { 0,     3,     2.79306201502264705994E+00},
            { 0,     4,    -2.56152597553130023655E+00},

            { 1,     0,     7.89806670922946318925E-01},
            { 1,     1,     1.42719958374201905116E+01},
            { 1,     2,    -2.47713064947769225910E+01},
            { 1,     3,     8.04278071310951681028E+00},
            { 1,     4,     3.66647228234742073028E+00},

            { 2,     0,     1.66175550737136745738E+01},
            { 2,     1,    -7.35725916220352331720E+01},
            { 2,     2,     1.08031105163848295092E+02},
            { 2,     3,    -5.83176668396733290933E+01},
            { 2,     4,     6.99282806748704732769E+00},

            { 3,     0,    -4.02973895611078916090E+01},
            { 3,     1,     1.40473419123728177738E+02},
            { 3,     2,    -1.80465312924060071964E+02},
            { 3,     3,     9.78861437573207098239E+01},
            { 3,     4,    -1.77262651464499150222E+01},

            { 4,     0,     2.12553651132979481986E+01},
            { 4,     1,    -7.07553640430039791909E+01},
            { 4,     2,     8.71360130418007798880E+01},
            { 4,     3,    -4.65547494341429839437E+01},
            { 4,     4,     9.00395775876992843223E+00}
        },
        {                                                                               // LMR2 (Z003000)
            { 0,     0,     1.56146519862202897144E+01},
            { 0,     1,    -3.25295555117365031705E+00},
            { 0,     2,     3.80218538552061113833E+00},
            { 0,     3,    -2.85992446872646643996E+00},
            { 0,     4,     1.05223371503530960247E+00},
            { 0,     5,    -1.49446662816108466476E-01},

            { 1,     0,     7.77113365437422576321E+00},
            { 1,     1,    -2.69242148162415411150E+01},
            { 1,     2,     4.21377981796297760297E+01},
            { 1,     3,    -3.05791586287080576767E+01},
            { 1,     4,     1.03058811476673746199E+01},
            { 1,     5,    -1.30770786874792355192E+00},

            { 2,     0,     4.37708329067558921111E+00},
            { 2,     1,     4.05698553767332636966E+00},
            { 2,     2,    -4.39975579696197893753E+00},
            { 2,     3,    -4.02226205989940837071E+01},
            { 2,     4,     3.97982102495747582793E+01},
            { 2,     5,    -9.68255076300243366916E+00},

            { 3,     0,    -1.14171738557163791938E+02},
            { 3,     1,     1.44654983520949940612E+02},
            { 3,     2,     1.25658014849041876460E+02},
            { 3,     3,    -7.61233422914715873731E+01},
            { 3,     4,    -4.97277512863255921616E+01},
            { 3,     5,     2.27566957849081354937E+01},

            { 4,     0,     8.52284180915415390700E+02},
            { 4,     1,    -2.11793768130346688849E+03},
            { 4,     2,     1.43778303702111816165E+03},
            { 4,     3,    -3.53485239324663382376E+02},
            { 4,     4,     5.96099116904753358881E+01},
            { 4,     5,    -1.52841095719764492600E+01},

            { 5,     0,    -7.27744971361332773085E+02},
            { 5,     1,     2.37447873535868984618E+03},
            { 5,     2,    -2.29824328800873217915E+03},
            { 5,     3,     9.63017826362198434254E+02},
            { 5,     4,    -1.94941240897338843752E+02},
            { 5,     5,     1.82013102751265911650E+01}
        },
        {                                                                               // LMA (Z003000)
            { 0,    -5,    -2.57403227906218804419E+07},
            { 0,    -4,     1.91919278581367194653E+08},
            { 0,    -3,    -6.42722772365061283112E+08},
            { 0,    -2,     1.26557774018711972237E+09},
            { 0,    -1,    -1.59807275300796556473E+09},
            { 0,     0,     1.29431773149255251884E+09},
            { 0,     1,    -5.73752193441836237907E+08},
            { 0,     2,    -3.44428973121863901615E+07},
            { 0,     3,     2.58171971437299102545E+08},
            { 0,     4,    -2.08805258131750792265E+08},
            { 0,     5,     9.92256309368441253901E+07},
            { 0,     6,    -3.16719607180158011615E+07},
            { 0,     7,     6.90780604186993371695E+06},
            { 0,     8,    -9.94125170976186404005E+05},
            { 0,     9,     8.54546806862898811232E+04},
            { 0,    10,    -3.33231922729193092891E+03},

            { 2,    -5,     1.34762385503116726875E+08},
            { 2,    -4,    -1.74399176370349168777E+09},
            { 2,    -3,     8.78154209435013771057E+09},
            { 2,    -2,    -2.48007449760757026672E+10},
            { 2,    -1,     4.53496317537334213257E+10},
            { 2,     0,    -5.77108539579762878418E+10},
            { 2,     1,     5.31677144188315505981E+10},
            { 2,     2,    -3.61960228471166305542E+10},
            { 2,     3,     1.83377540802361297607E+10},
            { 2,     4,    -6.87766287682024002075E+09},
            { 2,     5,     1.87236466119424223900E+09},
            { 2,     6,    -3.54625421447190940380E+08},
            { 2,     7,     4.26550618172946497798E+07},
            { 2,     8,    -2.50770450767832808197E+06},
            { 2,     9,    -2.72280137764262872224E+04},
            { 2,    10,     8.76822961747765839391E+03},

            { 4,    -5,     3.64113230917176294327E+09},
            { 4,    -4,    -2.19993132639651412964E+10},
            { 4,    -3,     5.51823742774377899170E+10},
            { 4,    -2,    -6.74001478202328414917E+10},
            { 4,    -1,     1.97503194124485588074E+10},
            { 4,     0,     6.55027963168653793335E+10},
            { 4,     1,    -1.20047754730772918701E+11},
            { 4,     2,     1.10973078984599426270E+11},
            { 4,     3,    -6.64446966999905242920E+10},
            { 4,     4,     2.71468286742218437195E+10},
            { 4,     5,    -7.48815090425585556030E+09},
            { 4,     6,     1.28782471484118795395E+09},
            { 4,     7,    -9.90502311621036529541E+07},
            { 4,     8,    -7.30136450831011310220E+06},
            { 4,     9,     2.22677327622929215431E+06},
            { 4,    10,    -1.42107178068599314429E+05},

            { 6,    -5,    -1.69094733263322086334E+10},
            { 6,    -4,     1.18528670819833084106E+11},
            { 6,    -3,    -3.71969875084951599121E+11},
            { 6,    -2,     6.88453453818612182617E+11},
            { 6,    -1,    -8.33162353980749511719E+11},
            { 6,     0,     6.92006751423001220703E+11},
            { 6,     1,    -4.06654644455712707520E+11},
            { 6,     2,     1.78505518636014801025E+11},
            { 6,     3,    -6.90610161810206146240E+10},
            { 6,     4,     3.03147898201140594482E+10},
            { 6,     5,    -1.43313685446296138763E+10},
            { 6,     6,     5.54512398337787055969E+09},
            { 6,     7,    -1.50349138969701838493E+09},
            { 6,     8,     2.63695637427561759949E+08},
            { 6,     9,    -2.68951201610311269760E+07},
            { 6,    10,     1.21421668571477057412E+06},

            { 8,    -5,     1.62520989450171852112E+10},
            { 8,    -4,    -1.22380441325664489746E+11},
            { 8,    -3,     4.03718140234380981445E+11},
            { 8,    -2,    -7.62315076723765625000E+11},
            { 8,    -1,     8.94726561199124023438E+11},
            { 8,     0,    -6.46755386562959716797E+11},
            { 8,     1,     2.37624479000744293213E+11},
            { 8,     2,     2.53444076942211761475E+10},
            { 8,     3,    -7.69330323678240509033E+10},
            { 8,     4,     3.81364515969743728638E+10},
            { 8,     5,    -6.76366483839152717590E+09},
            { 8,     6,    -1.57825522737156033516E+09},
            { 8,     7,     1.18057807852275967598E+09},
            { 8,     8,    -2.92124731642549633980E+08},
            { 8,     9,     3.56642162956294342875E+07},
            { 8,    10,    -1.79988696399625157937E+06},

            {10,    -5,    -5.29522123969188880920E+09},
            {10,    -4,     4.37924632110345230103E+10},
            {10,    -3,    -1.66864842737637054443E+11},
            {10,    -2,     3.77193811645369873047E+11},
            {10,    -1,    -5.45755946513809143066E+11},
            {10,     0,     5.12213537621273925781E+11},
            {10,     1,    -2.95784493345500854492E+11},
            {10,     2,     7.77097973981596527100E+10},
            {10,     3,     2.15638580927770462036E+10},
            {10,     4,    -2.78122413775628700256E+10},
            {10,     5,     1.08470967326972980499E+10},
            {10,     6,    -1.76167227877269959450E+09},
            {10,     7,    -1.19039692308903947473E+08},
            {10,     8,     1.03663300992385908961E+08},
            {10,     9,    -1.72756292958132103086E+07},
            {10,    10,     1.02701810218507179525E+06},

            {12,    -5,     1.80343253078517799377E+10},
            {12,    -4,    -1.26691173309104354858E+11},
            {12,    -3,     3.56351695791545471191E+11},
            {12,    -2,    -5.53758351468571411133E+11},
            {12,    -1,     5.36138039982130004883E+11},
            {12,     0,    -3.39866845315202819824E+11},
            {12,     1,     1.44163375900454193115E+11},
            {12,     2,    -4.49762794375830764771E+10},
            {12,     3,     1.71580814489644603729E+10},
            {12,     4,    -1.12449749792033348083E+10},
            {12,     5,     6.78490694496251773834E+09},
            {12,     6,    -2.84463375643244552612E+09},
            {12,     7,     7.96841892772628188133E+08},
            {12,     8,    -1.43851289638218492270E+08},
            {12,     9,     1.52044649291066508740E+07},
            {12,    10,    -7.16406403363849385642E+05},

            {14,    -5,     1.16520079605483184814E+11},
            {14,    -4,    -5.39465539126750244141E+11},
            {14,    -3,     1.09764790470851672363E+12},
            {14,    -2,    -1.26552677553650268555E+12},
            {14,    -1,     8.70273333877262084961E+11},
            {14,     0,    -3.15431150532989440918E+11},
            {14,     1,     9.53251283628360748291E+08},
            {14,     2,     5.57105743806749877930E+10},
            {14,     3,    -2.28806497057574310303E+10},
            {14,     4,     3.86158310274266898632E+08},
            {14,     5,     3.09330777606074094772E+09},
            {14,     6,    -1.37702519450405025482E+09},
            {14,     7,     3.09963461585896790028E+08},
            {14,     8,    -3.97062345451950877905E+07},
            {14,     9,     2.66685145412934198976E+06},
            {14,    10,    -6.71761975842174288118E+04}
        },
        {                                                                               // HM (Z003000)
            { 0,     0,     3.12154248604355193675E+06},
            { 0,     1,     1.40693554095653563738E+08},
            { 0,     2,    -3.71595276305998325348E+08},
            { 0,     3,     1.76670755437760323286E+08},
            { 0,     4,     2.02951692746982872486E+08},
            { 0,     5,    -2.92336468998508751392E+08},
            { 0,     6,     1.68110673419797122478E+08},
            { 0,     7,    -5.15604220053785145283E+07},
            { 0,     8,     4.97172830721258837730E+06},
            { 0,     9,     1.57284540400977246463E+06},
            { 0,    10,    -1.75988885345202215831E+04},
            { 0,    11,    -2.97943821613070613239E+05},
            { 0,    12,     8.61773562932996719610E+04},
            { 0,    13,    -9.04684915798255497066E+03},
            { 0,    14,     5.27682281201942714688E+02},
            { 0,    15,    -4.83809114868389968933E+01},

            { 1,     0,    -7.35394878249169588089E+07},
            { 1,     1,    -5.76316277886640667915E+08},
            { 1,     2,     2.41011072754778337479E+09},
            { 1,     3,    -2.36860622483675336838E+09},
            { 1,     4,     7.09952326698126316071E+08},
            { 1,     5,     2.10905643436775207520E+08},
            { 1,     6,    -2.58844355418170630932E+08},
            { 1,     7,     1.12699839262904345989E+08},
            { 1,     8,    -2.48447791184037849307E+07},
            { 1,     9,     3.05833750104070967063E+06},
            { 1,    10,    -2.27715413729118229821E+06},
            { 1,    11,     1.28163629808912402950E+06},
            { 1,    12,    -2.45961475246419315226E+05},
            { 1,    13,     7.71686543784326659079E+03},
            { 1,    14,     9.24986993266994090845E+02},
            { 1,    15,     1.19680958374041750858E+02},

            { 2,     0,     3.01712491925559639931E+08},
            { 2,     1,     1.35997955868716746569E+08},
            { 2,     2,    -4.33322825757982063293E+09},
            { 2,     3,     5.50388031191586112976E+09},
            { 2,     4,    -2.64093562924157762527E+09},
            { 2,     5,     6.12809232748553872108E+08},
            { 2,     6,    -9.04877006798645257950E+07},
            { 2,     7,     9.55706229663263447583E+06},
            { 2,     8,    -1.24463357289487235248E+07},
            { 2,     9,     1.14754493813732806593E+07},
            { 2,    10,    -2.62512806448984425515E+06},
            { 2,    11,    -1.38746075201159896096E+05},
            { 2,    12,    -3.06434607562738347042E+04},
            { 2,    13,     5.51760989160078097484E+04},
            { 2,    14,    -7.38420266268990508252E+03},
            { 2,    15,    -7.88663649547120257921E+01},

            { 3,     0,    -2.52314621471064984798E+08},
            { 3,     1,     1.59695673223666906357E+09},
            { 3,     2,     2.56902829277027750015E+09},
            { 3,     3,    -4.80529198174804210663E+09},
            { 3,     4,     2.24889814265354108810E+09},
            { 3,     5,    -3.74476785226586341858E+08},
            { 3,     6,     5.43274903410231973976E+06},
            { 3,     7,     4.29081287132954373956E+07},
            { 3,     8,    -2.75819523820120133460E+07},
            { 3,     9,     5.40615162239633128047E+06},
            { 3,    10,    -1.34637905180512159131E+06},
            { 3,    11,     8.50683737928409944288E+05},
            { 3,    12,    -1.24430090013735956745E+05},
            { 3,    13,    -2.52241309959800564684E+04},
            { 3,    14,     3.73821107091213934837E+03},
            { 3,    15,     3.32248890148269254041E+02},

            { 4,     0,    -5.60517200779632091522E+08},
            { 4,     1,    -1.89984697471257662773E+09},
            { 4,     2,     6.79333540273263901472E+07},
            { 4,     3,     1.90352540671050238609E+09},
            { 4,     4,    -7.50747156377030849457E+08},
            { 4,     5,    -4.72567503656732067466E+07},
            { 4,     6,     3.00467742178344316781E+07},
            { 4,     7,    -1.93436903588892333210E+06},
            { 4,     8,     2.95066336620633816347E+06},
            { 4,     9,     1.70540967107301065698E+06},
            { 4,    10,    -1.16727812393627688289E+06},
            { 4,    11,     4.34517006904071604367E+04},
            { 4,    12,     2.35401513264981185785E+04},
            { 4,    13,     6.11857844834931529476E+03},
            { 4,    14,     1.03171863933184098983E+02},
            { 4,    15,    -3.75092728367472147966E+02},

            { 5,     0,     1.24819667430232143402E+09},
            { 5,     1,     2.19967497561418384314E+08},
            { 5,     2,    -7.89014919695381075144E+07},
            { 5,     3,    -5.97878301337307453156E+08},
            { 5,     4,     1.59211666115665256977E+08},
            { 5,     5,     5.16454561360483840108E+07},
            { 5,     6,    -4.18833435321484040469E+06},
            { 5,     7,    -3.39808681376820290461E+06},
            { 5,     8,    -7.22078069451165269129E+05},
            { 5,     9,     1.99512798423410247779E+05},
            { 5,    10,     9.07387635003017785493E+04},
            { 5,    11,    -1.52441052433743388974E+04},
            { 5,    12,    -8.25416524310134082043E+03},
            { 5,    13,     2.18147893849961155865E+03},
            { 5,    14,    -6.64708537198572457783E+02},
            { 5,    15,     1.44256443204359840138E+02},

            { 6,     0,    -8.01485488430039882660E+08},
            { 6,     1,     5.38297769089608192444E+08},
            { 6,     2,    -6.06655185404267787933E+08},
            { 6,     3,     4.95703300265223979950E+08},
            { 6,     4,    -1.09494887044486209750E+08},
            { 6,     5,     1.03584531923860888928E+07},
            { 6,     6,    -5.11998278185982257128E+06},
            { 6,     7,    -3.64656091820526169613E+06},
            { 6,     8,     1.75311863341458095238E+06},
            { 6,     9,    -4.11223165191869309638E+05},
            { 6,    10,     1.88387970340288331499E+05},
            { 6,    11,    -3.47362162359701178502E+04},
            { 6,    12,     8.57684754506128047069E+03},
            { 6,    13,    -2.74638773621986092621E+03},
            { 6,    14,     9.61819122641777823901E+01},
            { 6,    15,     2.69118539619979557642E+01},

            { 7,     0,    -7.32038980764418244362E+07},
            { 7,     1,     2.00004937755068242550E+08},
            { 7,     2,     1.87085510840455114841E+08},
            { 7,     3,    -2.40693202493527889252E+08},
            { 7,     4,     2.91464472802970111370E+07},
            { 7,     5,     1.80457588328361362219E+07},
            { 7,     6,    -5.22623080352768115699E+06},
            { 7,     7,     6.93950939136536209844E+05},
            { 7,     8,     1.29864879686568258330E+06},
            { 7,     9,    -6.70193084813015302643E+05},
            { 7,    10,     9.20749454564202897018E+04},
            { 7,    11,    -2.00941969415327257593E+04},
            { 7,    12,     4.96906965956355270464E+03},
            { 7,    13,    -3.33430820278479188801E+02},
            { 7,    14,     2.60022419015300613410E+02},
            { 7,    15,    -5.26054399949432109906E+01},

            { 8,     0,     3.49899361381577372551E+08},
            { 8,     1,    -6.16363983678707957268E+08},
            { 8,     2,     2.87744997564267575741E+08},
            { 8,     3,    -5.80734918315911106765E+06},
            { 8,     4,    -7.49858021421681251377E+06},
            { 8,     5,    -1.02456923512990288436E+07},
            { 8,     6,     5.04230622996888868511E+06},
            { 8,     7,     8.30524694711205520434E+04},
            { 8,     8,    -8.61025903787584858947E+05},
            { 8,     9,     1.34900060086739453254E+05},
            { 8,    10,     8.92692328603777714306E+04},
            { 8,    11,    -2.89001283128855830000E+04},
            { 8,    12,     2.12396475300008751219E+03},
            { 8,    13,     2.95241329416678070174E+02},
            { 8,    14,    -1.65338632433199251182E+02},
            { 8,    15,     2.41591961058125228590E+01},

            { 9,     0,    -1.71260510420660376549E+08},
            { 9,     1,     3.28761719424411475658E+08},
            { 9,     2,    -2.11206060614314973354E+08},
            { 9,     3,     3.46205389136431142688E+07},
            { 9,     4,     2.04971128339806981385E+07},
            { 9,     5,    -1.15401356420353483409E+07},
            { 9,     6,     2.01258777025233558379E+06},
            { 9,     7,    -1.79581004342673142673E+05},
            { 9,     8,     2.07729967156824888662E+05},
            { 9,     9,    -9.74314817908821714809E+04},
            { 9,    10,     2.22068544668985268800E+04},
            { 9,    11,    -8.59795432872551464243E+03},
            { 9,    12,     3.31306345412375048909E+03},
            { 9,    13,    -7.04637490543811736643E+02},
            { 9,    14,     9.43370988450807317349E+01},
            { 9,    15,    -7.12454991962896322377E+00},

            {10,     0,     2.84390209118411280215E+07},
            {10,     1,    -5.95417063166125416756E+07},
            {10,     2,     4.45506571666083186865E+07},
            {10,     3,    -9.78308669950580969453E+06},
            {10,     4,    -5.62624247386193368584E+06},
            {10,     5,     4.15237079235649015754E+06},
            {10,     6,    -5.87503508675953838974E+05},
            {10,     7,    -4.26206334621473739389E+05},
            {10,     8,     2.55988526607063336996E+05},
            {10,     9,    -6.22110169340419161017E+04},
            {10,    10,     3.20104591341956165707E+03},
            {10,    11,     3.22804174274844626780E+03},
            {10,    12,    -1.24817289214048719259E+03},
            {10,    13,     2.26605355126834183466E+02},
            {10,    14,    -2.29672306377457609017E+01},
            {10,    15,     1.15281981127687060962E+00}
        },
        {                                                                               // RECOM (Z003000)
            { 0,     0,     1.34060093558913830947E+01},
            { 0,     1,    -7.67393959247841483950E-01},
            { 0,     2,     5.67829538231681030247E-01},
            { 0,     3,    -1.53301717168738610431E-01},

            { 1,     0,     1.68198984963847353313E+00},
            { 1,     1,     2.75310597472883875070E-02},
            { 1,     2,    -3.85904737173993372945E-01},
            { 1,     3,     1.74248736565002115828E-01},

            { 2,     0,    -1.42654332247593096383E+00},
            { 2,     1,     1.33518023284770448456E+00},
            { 2,     2,    -3.89425375461592448989E-01},
            { 2,     3,    -1.47490049820737292169E-02},

            { 3,     0,     5.53979839492381387345E-01},
            { 3,     1,    -6.97392793919068609831E-01},
            { 3,     2,     2.81780811937450081928E-01},
            { 3,     3,    -2.63299277947934735888E-02}
        }
    }
};



#endif // __constants_h__
