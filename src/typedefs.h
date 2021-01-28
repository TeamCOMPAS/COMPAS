#ifndef __typedefs_h__
#define __typedefs_h__

#include "constants.h"

// JR: todo: clean this up and document it better


typedef std::tuple<bool, COMPAS_VARIABLE_TYPE> COMPAS_VARIABLE;
typedef std::initializer_list<STELLAR_TYPE>    STELLAR_TYPE_LIST;
typedef std::initializer_list<SN_EVENT>        SN_EVENT_LIST;
typedef std::vector<STELLAR_TYPE>              STYPE_VECTOR;


// Log file details
typedef struct LogfileDetails {
    int                      id;
    std::string              filename;
    ANY_PROPERTY_VECTOR      recordProperties;
    std::vector<TYPENAME>    propertyTypes;
    std::vector<std::string> hdrStrings;
    std::vector<std::string> unitsStrings;
    std::vector<std::string> typeStrings;
    std::vector<std::string> fmtStrings;
} LogfileDetailsT;


// Grid file details
typedef struct Gridfile {
    std::string   filename;                                 // filename for grid file
    ERROR         error;                                    // status - ERROR::NONE if no problem, otherwise an error number
    std::ifstream handle;                                   // the file handle
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
// CCSN, ECSN, PISN, PPSIN, USSN

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

    KickParameters   initialKickParameters;                 // User-supplied initial kick parameters - if present used in place of drawing randomly/from distributions
    
    double           coreMassAtCOFormation;                 // Core mass of this star when it formed a compact object
    double           COCoreMassAtCOFormation;               // Carbon Oxygen core mass of the star when it goes supernova and forms a compact object
    double           drawnKickMagnitude;                    // Kick magnitude the system received during the supernova (km s^-1)
    double           eccentricAnomaly;                      // Eccentric anomaly at instataneous time of the SN
    SNEventsT        events;                                // Record of supernova events undergone by the star
    double           fallbackFraction;                      // Fallback fraction during a supernova event
    double           HeCoreMassAtCOFormation;               // Helium core mass of the star when it goes supernova and forms a compact objec
    bool             isHydrogenPoor;                        // Flag to indicate if exploding star is hydrogen-poor. We consider an H-rich star all SN progenitors that have an H envelope, otherwise H-poor
    double           kickMagnitude;                         // Kick magnitude the system received during the supernova (km s^-1)
    double           kickMagnitudeRandom;                   // Random number U(0,1) for choosing the supernova kick magnitude - drawn once at star creation
    double           meanAnomaly;                           // Mean anomaly at instantaneous time of the SN - uniform in [0, 2pi]
    double           phi;                                   // Kick angle in the orbital plane, defined CCW from the radial vector pointed away from the Companion (rad) [0, 2pi)
    SN_STATE         supernovaState;                        // Indicates which star (or stars) are undergoing / have undergone a supernova event
    double           theta;                                 // Kick angle out of the orbital plane, toward the orbital angular momentum axis (rad) [-pi/2, pi/2]
    double           totalMassAtCOFormation;                // Total mass of the star when it goes supernova and forms a compact object
    double           trueAnomaly;                           // True anomaly at instantaneous time of the SN
} SupernovaDetailsT;


// pulsar parameters (if star becomes a Neutron Star)
typedef struct PulsarDetails {
    double magneticField;                                   // Pulsar magnetic field strength (G)
    double spinPeriod;                                      // Pulsar spin period (ms)
    double spinFrequency;                                   // Pulsar spin frequency in rads per second
    double spinDownRate;                                    // Pulsar spin down rate (Pdot, dimensionless)
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
} BindingEnergiesT;


// RLOF properties
typedef struct RLOFProperties {
    OBJECT_ID     id;
    unsigned long randomSeed;

    STELLAR_TYPE  stellarType1;
    STELLAR_TYPE  stellarType2;

    double        mass1;
    double        mass2;

    double        radius1;
    double        radius2;

    double        eccentricity;
    double        semiMajorAxis;

    unsigned int  eventCounter;

    double        time;

    bool          isRLOF1;
    bool          isRLOF2;

    bool          isCE;

} RLOFPropertiesT;

typedef struct BinaryRLOFDetails {                          // RLOF details pertinent to binaries

    bool experiencedRLOF;
    bool immediateRLOFPostCEE;                              // Here for now - maybe should be in Binary CEDetails struct?       JR: todo:
    bool isRLOF;
    bool simultaneousRLOF;                                  // Here for now - maybe should be in Binary CEDetails struct?       JR: todo:
    bool stableRLOFPostCEE;                                 // Here for now - maybe should be in Binary CEDetails struct?       JR: todo:
    RLOFPropertiesT  props1;
    RLOFPropertiesT  props2;
    RLOFPropertiesT* currentProps;
    RLOFPropertiesT* previousProps;
} BinaryRLOFDetailsT;

typedef struct StellarRLOFDetails {                         // RLOF details pertinent to individual stars
    bool   isRLOF;
    bool   experiencedRLOF;
    bool   RLOFPostCEE;
} StellarRLOFDetailsT;


// BeBinary properties
typedef struct BeBinaryProperties {
    OBJECT_ID     id;
    unsigned long randomSeed;

    double        dt;
    double        totalTime;

    double        massNS;

    double        companionMass;
    double        companionLuminosity;
    double        companionTeff;
    double        companionRadius;

    double        semiMajorAxis;
    double        eccentricity;
} BeBinaryPropertiesT;

typedef struct BeBinaryDetails {
    BeBinaryPropertiesT  props1;
    BeBinaryPropertiesT  props2;
    BeBinaryPropertiesT* currentProps;
    BeBinaryPropertiesT* previousProps;
} BeBinaryDetailsT;


// Common Envelope properties
typedef struct BinaryCEESavedValues {
    double eccentricity;
   	double rocheLobe1to2;
	double rocheLobe2to1;
    double semiMajorAxis;
} BinaryCEESavedValuesT;

typedef struct BinaryCEDetails {                            // Common Envelope details pertinent to binaries
    BinaryCEESavedValuesT preCEE;
    BinaryCEESavedValuesT postCEE;

    bool                  CEEnow;                           // Indicates whether a common envelope event is occurring now
    unsigned int          CEEcount;                         // Common Envelope Event count
    bool                  doubleCoreCE;
    bool                  optimisticCE;
} BinaryCEDetailsT;


typedef struct StellarCEESavedValues {
    double       bindingEnergy;
    double       dynamicalTimescale;
    double       luminosity;
    double       mass;
    double       nuclearTimescale;
    double       radialExpansionTimescale;
    double       radius;
    STELLAR_TYPE stellarType;
    double       temperature;
    double       thermalTimescale;
} StellarCEESavedValuesT;

typedef struct StellarCEDetails {                      // Common Envelope details pertinent to individual stars
    StellarCEESavedValuesT preCEE;
    StellarCEESavedValuesT postCEE;

    double                 bindingEnergy;
    double                 COCoreMass;
    double                 CoreMass;
    double                 HeCoreMass;
    double                 lambda;
} StellarCEDetailsT; // was CommonEnvelopeDetailsT;


#endif // __typedefs_h__
