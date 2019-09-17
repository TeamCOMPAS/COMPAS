#ifndef __typedefs_h__
#define __typedefs_h__

#include "constants.h"

// JR: todo: clean this up and document it better


typedef std::tuple<bool, COMPAS_VARIABLE_TYPE>                                      COMPAS_VARIABLE;
typedef std::initializer_list<STELLAR_TYPE>                                         STELLAR_TYPE_LIST;
typedef std::tuple<int, std::string, ANY_PROPERTY_VECTOR, std::vector<std::string>> LOGFILE_DETAILS;



// RotationalVelocityParams struct for gsl root solver
struct RotationalVelocityParams {                                           // Structure containing parameter (u) for the root solving function using gsl_root_solver
    double u;                                                               // Value of CDF, draw in U(0,1)
};

// KickVelocityParams struct for gsl root solver
struct KickVelocityParams {
    double y;       // Value of CDF, should be drawn as U(0,1)
    double sigma;   // sigma for kick distribution
};


typedef struct SNEvents {
    SN_EVENT              now;                                              // Supernova status at the current timestep: NONE if no supernova event happening
    std::vector<SN_EVENT> past;                                             // Supernova status at any past timestep   : NONE if no supernova event happened in any past timestep
} SNEventsT;

// struct for supernova attributes of the base star
    // some of these are only required for binary stars, but
    // easier (and more logical) to keep all SN-related attributes
    // in the same place

typedef struct SupernovaDetails {                                       // Holds attributes, flags - if the star went supernova
    double           coreMassAtCOFormation;                             // Core mass of this star when it formed a compact object
    double           COCoreMassAtCOFormation;                           // Carbon Oxygen core mass of the star when it goes supernova and forms a compact object
    double           drawnKickVelocity;                                 // Kick velocity the system received during the supernova (km s^-1)
    double           eccentricAnomaly;                                  // Eccentric anomaly at instataneous time of the SN
    SNEventsT        events;
    double           fallbackFraction;                                  // Fallback fraction during a supernova event
    double           HeCoreMassAtCOFormation;                           // Helium core mass of the star when it goes supernova and forms a compact objec
    HYDROGEN_CONTENT hydrogenContent;		                            // Hydrogen content of the exploding star. We consider an H-rich star all SN progenitors that have an H envelope, otherwise H-poor
    double           kickVelocity;                                      // Kick velocity the system received during the supernova (km s^-1)
    double           meanAnomaly;                                       // Mean anomaly at instantaneous time of the SN - uniform in [0, 2pi]
    double           phi;                                               // Angle between 'x' and 'y', both in the orbital plane of supernovae vector (rad)
    SN_STATE         supernovaState;                                    // indicates which star (or stars) are undergoing / hove undergone a supernova event
    double           theta;                                             // Angle between the orbital plane and the 'z' axis of supernovae vector (rad)
    double           totalMassAtCOFormation;                            // Total mass of the star when it goes supernova and forms a compact object
    double           trueAnomaly;                                       // True anomaly at instantaneous time of the SN
    double           uRand;                                             // Random number U(0,1) for choosing the supernova kick velocity magnitude - drawn once at star creation
} SupernovaDetailsT;



// pulsar parameters (if star becomes a Neutron Star)
typedef struct PulsarDetails {
    double magneticField;                                               // Pulsar magnetic field strength (G)
    double spinPeriod;                                                  // Pulsar spin period (ms)
    double spinFrequency;                                               // Pulsar spin frequency in rads per second
    double spinDownRate;                                                // Pulsar spin down rate (Pdot, dimensionless)
} PulsarDetailsT;


// struct for Lambdas
typedef struct Lambdas {
	double dewi;
    double fixed;                                                       // set to OPTIONS->commonEnvelopeLambda
	double kruckow;                                                     // calculated using m_Radius and OPTIONS->commonEnvelopeSlopeKruckow
	double kruckowBottom;                                               // calculated using m_Radius and -1
	double kruckowMiddle;                                               // calculated using m_Radius and -4/5
	double kruckowTop;                                                  // calculated using m_Radius and -2/3
	double loveridge;                                                   // no mass loss
	double loveridgeWinds;                                              // mass loss    JR: todo: or would be if the parameter wasn't ignored...
	double nanjing;
} LambdasT;


// struct for Zetas
typedef struct Zetas {
	double hurley;
	double hurleyHe;
	double nuclear;
	double soberman;
	double sobermanHe;
	double thermal;
} ZetasT;


// struct for binding energies
typedef struct BindingEnergies {
    double fixed;                                                       // calculated using lambda = OPTIONS->commonEnvelopeLambda
	double nanjing;                                                     // calculated using lambda = m_Lambdas.nanjing
	double loveridge;                                                   // calculated using lambda = m_Lambdas.loveridge
	double loveridgeWinds;                                              // calculated using lambda = m_Lambdas.loveridgeWinds
	double kruckow;                                                     // calculated using lambda = m_Lambdas.kruckow
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

    double        separation;

    unsigned int  eventCounter;

    double        time;

    bool          isRLOF1;
    bool          isRLOF2;

    bool          isCE;

    bool          monitorMS1;
    bool          monitorMS2;
    bool          monitorHeMS1;
    bool          monitorHeMS2;
} RLOFPropertiesT;

typedef struct BinaryRLOFDetails {

    bool isRLOF;
    bool experiencedRLOF;
    bool stableRLOFPostCEE;

    RLOFPropertiesT currentProps;
    RLOFPropertiesT previousProps;
} BinaryRLOFDetailsT;

typedef struct StellarRLOFDetails {
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

    double        separation;
    double        eccentricity;
} BeBinaryPropertiesT;

typedef struct BeBinaryDetails {
    BeBinaryPropertiesT  props1;
    BeBinaryPropertiesT  props2;
    BeBinaryPropertiesT* currentProps;
    BeBinaryPropertiesT* previousProps;
} BeBinaryDetailsT;


// Common Envelope properties
typedef struct CommonEnvelopeDetails {
    double   bindingEnergy;
    double   COCoreMass;
    double   CoreMass;
    double   HeCoreMass;
    double   lambda;
} CommonEnvelopeDetailsT;


#endif // __typedefs_h__
