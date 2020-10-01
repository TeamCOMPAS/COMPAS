#ifndef __Options_H__
#define __Options_H__

#define OPTIONS Options::Instance()

#include <iostream>
#include <string>
#include <sstream>

#include <boost/algorithm/string.hpp>   // Boost string manipulation
#include <boost/program_options.hpp>    // Boost command line options tools
#include <boost/filesystem.hpp>         // Boost filesystem tools for handling paths etc.

#include "constants.h"
#include "typedefs.h"
#include "profiling.h"
#include "utils.h"
#include "Rand.h"
#include "changelog.h"

using std::string;
using std::vector;
using std::get;

namespace po = boost::program_options;

// for convenience

#define COMPLAIN(complainStr)           { std::stringstream _ss; _ss << complainStr; throw _ss.str(); }
#define COMPLAIN_IF(cond, complainStr)  { if (cond) COMPLAIN(complainStr) }



/*
 * Options Singleton
 *
 * Holds program options
 *
 * Singletons are sometimes frowned-upon, but doing it this way means
 * the program options don't need to be passed around to all and sundry.
 * I think convenience and clarity sometimes trump dogma.
 */

class Options {

public:
    
    
    // The OptionsValues class holds the values for the options.  This allows the Options class
    // to hold values for both the commandline options (the options specified by the user on the
    // commandline) and the grid file options (the options specified by the user on a grid file
    // record - on a per object (str/binary) basis).  When using grid files, the grid file options
    // take precedence over the commandline options: the grid file option values are set, then
    // options that were not specified in the grid file record are set from the command line
    // options.  Options that were not specified in the grid file record, and that were not
    // specified on the commandline, are set to the COMPAS default values.

    class OptionValues {

        friend class Options;                                                                                           // So the Options class can access members directly

        private:

            // member variables - alphabetically in groups (sort of...)

            bool                                        m_AllowRLOFAtBirth;                                             // Indicates whether binaries that have one or both stars in RLOF at birth are allowed to evolve
            bool                                        m_AllowTouchingAtBirth;                                         // Indicates whether binaries that are touching at birth are allowed to evolve

            bool                                        m_DebugToFile;                                                  // Flag used to determine whether debug statements should also be written to a log file
            bool                                        m_ErrorsToFile;                                                 // Flag used to determine whether error statements should also be written to a log file

            bool                                        m_EnableWarnings;                                               // Flag used to determine if warnings (via SHOW_WARN macros) should be displayed

	        bool                                        m_BeBinaries;													// Flag if we want to print BeBinaries (main.cpp)
            bool                                        m_EvolvePulsars;                                                // Whether to evolve pulsars or not
	        bool                                        m_EvolveUnboundSystems;							                // Option to chose if unbound systems are evolved until death or the evolution stops after the system is unbound during a SN.

            bool                                        m_DetailedOutput;                                               // Print detailed output details to file (default = false)
            bool                                        m_PopulationDataPrinting;                                       // Print certain data for small populations, but not for larger one
            bool                                        m_PrintBoolAsString;                                            // flag used to indicate that boolean properties should be printed as "TRUE" or "FALSE" (default is 1 or 0)
            bool                                        m_Quiet;                                                        // suppress some output
            bool                                        m_RlofPrinting;                                                 // RLOF printing
            bool                                        m_SwitchLog;                                                    // Print switch log details to file (default = false)

            int                                         m_nBatchesUsed;                                                 // Number of batches used, only needed for STROOPWAFEL (AIS) (default = -1, not needed)


            // Micellaneous evolution variables
            EVOLUTION_MODE                              m_EvolutionMode;                                                // Mode of evolution: SSE or BSE
            string                                      m_EvolutionModeString;

            int                                         m_ObjectsToEvolve;                                              // Number of stars (SSE) or binaries (BSE) to evolve
            bool                                        m_FixedRandomSeed;                                              // Whether to use a fixed random seed given by options.randomSeed (set to true if --random-seed is passed on command line)
            unsigned long int                           m_RandomSeed;                                                   // Random seed to use
    
            double                                      m_MaxEvolutionTime;                                             // Maximum time to evolve a binary by
            int                                         m_MaxNumberOfTimestepIterations;                                // Maximum number of timesteps to evolve binary for before giving up

            // Initial distribution variables

            double                                      m_InitialMass;                                                  // Initial mass of single star (SSE)
            double                                      m_InitialMass1;                                                 // Initial mass of primary (BSE)
            double                                      m_InitialMass2;                                                 // Initial mass of secondary (BSE)

            INITIAL_MASS_FUNCTION                       m_InitialMassFunction;                                          // Which initial mass function
            string                                      m_InitialMassFunctionString;
            double                                      m_InitialMassFunctionMin;                                       // Minimum mass to generate in Msol
            double                                      m_InitialMassFunctionMax;                                       // Maximum mass to generate in Msol
            double                                      m_InitialMassFunctionPower;                                     // single IMF power law set manually

            // Mass ratio
            MASS_RATIO_DISTRIBUTION                     m_MassRatioDistribution;                                        // Which mass ratio distribution
            string                                      m_MassRatioDistributionString;
            double                                      m_MassRatioDistributionMin;                                     // Minimum initial mass ratio when using a distribution
            double                                      m_MassRatioDistributionMax;                                     // Maximum initial mass ratio when using a distribution

            double                                      m_MinimumMassSecondary;                                         // Minimum mass of secondary to draw (in Msol)

            // Semi major axis
            double                                      m_SemiMajorAxis;                                                // Semi-major axis
            SEMI_MAJOR_AXIS_DISTRIBUTION                m_SemiMajorAxisDistribution;                                    // Which semi-major axis distribution
            string                                      m_SemiMajorAxisDistributionString;
            double                                      m_SemiMajorAxisDistributionMin;                                 // Minimum a in AU
            double                                      m_SemiMajorAxisDistributionMax;                                 // Maximum a in AU
            double                                      m_SemiMajorAxisDistributionPower;                               // Set semi-major axis distribution power law slope by hand

            // Period
            double                                      m_PeriodDistributionMin;                                        // Minimum initial period in days
            double                                      m_PeriodDistributionMax;                                        // Maximum initial period in days

            // Eccentricity
            double                                      m_Eccentricity;                                                 // Eccentricity
            ECCENTRICITY_DISTRIBUTION                   m_EccentricityDistribution;                                     // Which eccentricity distribution
            string                                      m_EccentricityDistributionString;
            double                                      m_EccentricityDistributionMin;                                  // Minimum initial eccentricity when using a distribution
            double                                      m_EccentricityDistributionMax;                                  // Maximum initial eccentricity when using a distribution

            // Kick options
            KICK_MAGNITUDE_DISTRIBUTION                 m_KickMagnitudeDistribution;                                    // Which kick magnitude distribution
            string                                      m_KickMagnitudeDistributionString;
            double                                      m_KickMagnitudeDistributionSigmaCCSN_NS;                        // Kick magnitude sigma in km s^-1 for neutron stars (default = "250" )
            double                                      m_KickMagnitudeDistributionSigmaCCSN_BH;                        // Kick magnitude sigma in km s^-1 for black holes (default = "250" )
            double                                      m_KickMagnitudeDistributionMaximum;                             // Maximum kick magnitude to draw. If negative, no maximum
	        double                                      m_KickMagnitudeDistributionSigmaForECSN;			            // Kick magnitude sigma for ECSN in km s^-1 (default = "0" )
	        double                                      m_KickMagnitudeDistributionSigmaForUSSN;			            // Kick magnitude sigma for USSN in km s^-1 (default = "20" )
	        double                                      m_KickScalingFactor;								            // Arbitrary factor for scaling kicks

            // Kick direction options
            KICK_DIRECTION_DISTRIBUTION                 m_KickDirectionDistribution;                                    // Kick direction distribution
            string                                      m_KickDirectionDistributionString;
            double                                      m_KickDirectionPower;                                           // Exponent

            // User-specified supernova parameter values
            double                                      m_KickMagnitude;                                                // Supernova kick magnitude - SSE
            double                                      m_KickMagnitude1;                                               // Supernova kick magnitude - BSE primary star
            double                                      m_KickMagnitude2;                                               // Supernova kick magnitude - BSE secondary star
            double                                      m_KickMagnitudeRandom;                                          // Random number U(0,1) for choosing the supernova kick magnitude - SSE
            double                                      m_KickMagnitudeRandom1;                                         // Random number U(0,1) for choosing the supernova kick magnitude - BSE primary star
            double                                      m_KickMagnitudeRandom2;                                         // Random number U(0,1) for choosing the supernova kick magnitude - BSE secondary star
            double                                      m_KickMeanAnomaly1;                                             // Mean anomaly at instantaneous time of the SN - uniform in [0, 2pi] - BSE primary star
            double                                      m_KickMeanAnomaly2;                                             // Mean anomaly at instantaneous time of the SN - uniform in [0, 2pi] - BSE secondary star
            double                                      m_KickPhi1;                                                     // Angle between 'x' and 'y', both in the orbital plane of supernovae vector (rad) - BSE primary star
            double                                      m_KickPhi2;                                                     // Angle between 'x' and 'y', both in the orbital plane of supernovae vector (rad) - BSE secondary star
            double                                      m_KickTheta1;                                                   // Angle between the orbital plane and the 'z' axis of supernovae vector (rad) - BSE primary star
            double                                      m_KickTheta2;                                                   // Angle between the orbital plane and the 'z' axis of supernovae vector (rad) - BSE secondar star

            // Black hole kicks
            BLACK_HOLE_KICK_OPTION                      m_BlackHoleKicksOption;                                         // Which black hole kicks mode
            string                                      m_BlackHoleKicksOptionString;

            // CHE - Chemically Homogeneous Evolution
            CHE_OPTION                                  m_CheOption;                                                    // Which Chemically Homogeneous Evolution mode
            string                                      m_CheString;

            // Supernova remnant mass
            REMNANT_MASS_PRESCRIPTION                   m_RemnantMassPrescription;                                      // Which remnant mass prescription
            string                                      m_RemnantMassPrescriptionString;

            SN_ENGINE                                   m_FryerSupernovaEngine;                                         // Which Fryer et al. supernova engine
            string                                      m_FryerSupernovaEngineString;

            NEUTRINO_MASS_LOSS_PRESCRIPTION             m_NeutrinoMassLossAssumptionBH;                                 // Which neutrino mass loss assumption for BH formation
            string                                      m_NeutrinoMassLossAssumptionBHString;
            double                                      m_NeutrinoMassLossValueBH;                                      // Value (corresponding to assumption) for neutrino mass loss for BH formation

            // Fixed uk options
            bool                                        m_UseFixedUK;                                                   // Whether to fix uk to a certain value (default is to NOT fix uk)
            double                                      m_FixedUK;                                                      // Dimensionless value to fix the kick magnitude to

            // Pair instability and pulsational pair instability mass loss
            bool                                        m_UsePairInstabilitySupernovae;                                 // Whether to use pair instability supernovae (PISN)
            double                                      m_PairInstabilityLowerLimit;                                    // Minimum core mass leading to PISN
            double                                      m_PairInstabilityUpperLimit;                                    // Maximum core mass leading to PISN

            bool                                        m_UsePulsationalPairInstability;                                // Whether to use pulsational pair instability (PPI)
            double                                      m_PulsationalPairInstabilityLowerLimit;                         // Maximum core mass leading to PPI
            double                                      m_PulsationalPairInstabilityUpperLimit;                         // Minimum core mass leading to PPI

            PPI_PRESCRIPTION                            m_PulsationalPairInstabilityPrescription;                       // Which PPI prescription
            string                                      m_PulsationalPairInstabilityPrescriptionString;

	        double                                      m_MaximumNeutronStarMass;						                // Maximum mass of a neutron star allowed, set to default in StarTrack

            // Setup default output directory and desired output directory
            string                                      m_OutputPathString;                                             // String to hold the output directory
            boost::filesystem::path                     m_DefaultOutputPath;                                            // Default output location
            boost::filesystem::path                     m_OutputPath;                                                   // Desired output location
            string                                      m_OutputContainerName;                                          // Name of output container (directory)

            // Mass loss options
            bool                                        m_UseMassLoss;                                                  // Whether to activate mass loss (default = True)

            // Can also have options for modifying strength of winds etc here

            MASS_LOSS_PRESCRIPTION                      m_MassLossPrescription;                                         // Which mass loss prescription
            string                                      m_MassLossPrescriptionString;

            double                                      m_LuminousBlueVariableFactor;                                   // Multiplicitive factor for luminous blue variable (LBV) mass loss rates
            double                                      m_WolfRayetFactor;                                              // Multiplicitive factor for Wolf-Rayet (WR) wind mass loss rates

            // Mass transfer options
            bool                                        m_UseMassTransfer;                                              // Whether to use mass transfer (default = false)
	        bool                                        m_CirculariseBinaryDuringMassTransfer;						    // Whether to circularise binary when it starts (default = false)
	        bool                                        m_AngularMomentumConservationDuringCircularisation;			    // Whether to conserve angular momentum while circularising or circularise to periastron (default = false)
	
            CASE_BB_STABILITY_PRESCRIPTION              m_CaseBBStabilityPrescription;									// Which prescription for the stability of case BB/BC mass transfer
            string                                      m_CaseBBStabilityPrescriptionString;

            double                                      m_MassTransferAdaptiveAlphaParameter;                           // Parameter used in adaptive RLOF to avoid overshoot of the solution
	        double                                      m_MaxPercentageAdaptiveMassTransfer;                            // As used in binary_c for adaptive RLOF prescription in mass transfer, according to notes.


            MT_ACCRETION_EFFICIENCY_PRESCRIPTION        m_MassTransferAccretionEfficiencyPrescription;                  // Which accretion efficiency prescription
            string                                      m_MassTransferAccretionEfficiencyPrescriptionString;

            double                                      m_MassTransferFractionAccreted;                                 // In mass transfer, ammount of mass transferred that is accreted. 1 for conservative, 0 for fully-non conservative.
            double                                      m_MassTransferCParameter;                                       // Detailed model parameter used in mass transfer
            double                                      m_EddingtonAccretionFactor;                                     // Multiplication factor for eddington accretion for NS & BH
                                                                                                                        // i.e. >1 is super-eddington
                                                                                                                        //       0. is no accretion

	        MT_THERMALLY_LIMITED_VARIATION              m_MassTransferThermallyLimitedVariation;                        // Choose how to deal with mass transfer if it is set as thermally limited.
	        string                                      m_MassTransferThermallyLimitedVariationString;

            double                                      m_MassTransferJloss;                                            // Specific angular momentum of the material leaving the system (not accreted)
            MT_ANGULAR_MOMENTUM_LOSS_PRESCRIPTION       m_MassTransferAngularMomentumLossPrescription;                  // Which mass transfer angular momentum loss prescription
            string                                      m_MassTransferAngularMomentumLossPrescriptionString;

            // Mass transfer rejuvenation prescription
            MT_REJUVENATION_PRESCRIPTION                m_MassTransferRejuvenationPrescription;                         // Which mass transfer rejuvenation prescription
            string                                      m_MassTransferRejuvenationPrescriptionString;

            // Mass transfer critical mass ratios
            bool                                        m_MassTransferCriticalMassRatioMSLowMass;                       // Whether to use critical mass ratios
            double                                      m_MassTransferCriticalMassRatioMSLowMassNonDegenerateAccretor;  // Critical mass ratio for MT from a MS low mass star
            double                                      m_MassTransferCriticalMassRatioMSLowMassDegenerateAccretor;     // Critical mass ratio for MT from a MS low mass star on to a degenerate accretor

            bool                                        m_MassTransferCriticalMassRatioMSHighMass;                      // Whether to use critical mass ratios
            double                                      m_MassTransferCriticalMassRatioMSHighMassNonDegenerateAccretor; // Critical mass ratio for MT from a MS high mass star
            double                                      m_MassTransferCriticalMassRatioMSHighMassDegenerateAccretor;    // Critical mass ratio for MT from a MS high mass star on to a degenerate accretor

            bool                                        m_MassTransferCriticalMassRatioHG;                              // Whether to use critical mass ratios
            double                                      m_MassTransferCriticalMassRatioHGNonDegenerateAccretor;         // Critical mass ratio for MT from a HG star
            double                                      m_MassTransferCriticalMassRatioHGDegenerateAccretor;            // Critical mass ratio for MT from a HG star on to a degenerate accretor

            bool                                        m_MassTransferCriticalMassRatioGiant;                           // Whether to use critical mass ratios
            double                                      m_MassTransferCriticalMassRatioGiantNonDegenerateAccretor;      // Critical mass ratio for MT from a giant
            double                                      m_MassTransferCriticalMassRatioGiantDegenerateAccretor;         // Critical mass ratio for MT from a giant on to a degenerate accretor

            bool                                        m_MassTransferCriticalMassRatioHeliumMS;                        // Whether to use critical mass ratios
            double                                      m_MassTransferCriticalMassRatioHeliumMSNonDegenerateAccretor;   // Critical mass ratio for MT from a Helium MS star
            double                                      m_MassTransferCriticalMassRatioHeliumMSDegenerateAccretor;      // Critical mass ratio for MT from a Helium MS star on to a degenerate accretor

            bool                                        m_MassTransferCriticalMassRatioHeliumHG;                        // Whether to use critical mass ratios
            double                                      m_MassTransferCriticalMassRatioHeliumHGNonDegenerateAccretor;   // Critical mass ratio for MT from a Helium HG star
            double                                      m_MassTransferCriticalMassRatioHeliumHGDegenerateAccretor;      // Critical mass ratio for MT from a Helium HG star on to a degenerate accretor

            bool                                        m_MassTransferCriticalMassRatioHeliumGiant;                     // Whether to use critical mass ratios
            double                                      m_MassTransferCriticalMassRatioHeliumGiantNonDegenerateAccretor;// Critical mass ratio for MT from a helium giant
            double                                      m_MassTransferCriticalMassRatioHeliumGiantDegenerateAccretor;   // Critical mass ratio for MT from a helium giant on to a degenerate accretor

            bool                                        m_MassTransferCriticalMassRatioWhiteDwarf;                      // Whether to use critical mass ratios
            double                                      m_MassTransferCriticalMassRatioWhiteDwarfNonDegenerateAccretor; // Critical mass ratio for MT from a white dwarf
            double                                      m_MassTransferCriticalMassRatioWhiteDwarfDegenerateAccretor;    // Critical mass ratio for MT from a white dwarf on to a degenerate accretor


            // Common Envelope options
            double                                      m_CommonEnvelopeAlpha;                                          // Common envelope efficiency alpha parameter (default = X)
            double                                      m_CommonEnvelopeLambda;                                         // Common envelope Lambda parameter (default = X)
	        double                                      m_CommonEnvelopeSlopeKruckow;									// Common envelope power factor for Kruckow fit normalized according to Kruckow+2016, Fig. 1
            double                                      m_CommonEnvelopeAlphaThermal;                                   // lambda = alpha_th*lambda_b + (1-alpha_th)*lambda_g
            double                                      m_CommonEnvelopeLambdaMultiplier;                               // Multiply common envelope lambda by some constant
            bool                                        m_AllowMainSequenceStarToSurviveCommonEnvelope;                 // Whether or not to allow a main sequence star to survive a common envelope event
    
            // Prescription for envelope state (radiative or convective)
            ENVELOPE_STATE_PRESCRIPTION                 m_EnvelopeStatePrescription;
            string                                      m_EnvelopeStatePrescriptionString;

            // Accretion during common envelope
            CE_ACCRETION_PRESCRIPTION                   m_CommonEnvelopeMassAccretionPrescription;
            string                                      m_CommonEnvelopeMassAccretionPrescriptionString;
            double                                      m_CommonEnvelopeMassAccretionMin;
            double                                      m_CommonEnvelopeMassAccretionMax;
            double                                      m_CommonEnvelopeMassAccretionConstant;

	        // Common envelope lambda prescription
	        CE_LAMBDA_PRESCRIPTION                      m_CommonEnvelopeLambdaPrescription;							    // Which prescription for CE lambda
	        string                                      m_CommonEnvelopeLambdaPrescriptionString;

	        // Common envelope Nandez and Ivanova energy formalism
	        bool                                        m_RevisedEnergyFormalismNandezIvanova;			                // Use the revised energy formalism from Nandez & Ivanova 2016 (default = false)
	        double                                      m_MaximumMassDonorNandezIvanova;								// Maximum mass allowed to use the revised energy formalism in Msol (default = 2.0)
	        double                                      m_CommonEnvelopeRecombinationEnergyDensity;					    // Factor using to calculate the binding energy depending on the mass of the envelope. (default = 1.5x10^13 erg/g)


            // Adaptive Importance Sampling options
            bool                                        m_AISexploratoryPhase;                                          // Flag if we want to run Exploratory phase of Adaptive Importance Sampling // Floor
            AIS_DCO                                     m_AISDCOtype;                                                   // Which prescription for DCO type
            string                                      m_AISDCOtypeString;
            bool                                        m_AIShubble;                                                    // Whether to exclude DCOs that not merge within Hubble
            bool                                        m_AISpessimistic;                                               // Whether to exclude Optimistic binaries
            bool                                        m_AISrefinementPhase;                                           // Flag if we want to run refinement phase of Adaptive Importance Sampling
            bool                                        m_AISrlof;                                                      // Whether to exclude binaries that have RLOFSecondaryZAMS
            double                                      m_KappaGaussians;                                               // Scaling factor for the width of the Gaussian distributions in AIS main sampling phase [should be in [0,1]]


            // Zetas
            ZETA_PRESCRIPTION                           m_StellarZetaPrescription;                                 	    // Which prescription to use for calculating stellar zetas (default = SOBERMAN)
            string                                      m_StellarZetaPrescriptionString;
	        double                                      m_ZetaAdiabaticArbitrary;
	        double                                      m_ZetaMainSequence;
            double                                      m_ZetaRadiativeEnvelopeGiant;


            // Metallicity options
            bool                                        m_FixedMetallicity;                                             // Whether user has specified a metallicity to use
            double                                      m_Metallicity;                                                  // Metallicity

            double                                      m_mCBUR1;                                                       // Minimum core mass at base of the AGB to avoid fully degenerate CO core formation


            // Neutron star equation of state
            NS_EOS                                      m_NeutronStarEquationOfState;                                   // NS EOS
            string                                      m_NeutronStarEquationOfStateString;


            // Pulsar birth magnetic field distribution string
            PULSAR_BIRTH_MAGNETIC_FIELD_DISTRIBUTION    m_PulsarBirthMagneticFieldDistribution;                         // Birth magnetic field distribution for pulsars
            string                                      m_PulsarBirthMagneticFieldDistributionString;
            double                                      m_PulsarBirthMagneticFieldDistributionMin;                      // Minimum birth magnetic field (log10 B/G)
            double                                      m_PulsarBirthMagneticFieldDistributionMax;                      // Maximum birth magnetic field (log10 B/G)

            // Pulsar birth spin period distribution string
            PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION       m_PulsarBirthSpinPeriodDistribution;                            // Birth spin period distribution for pulsars
            string                                      m_PulsarBirthSpinPeriodDistributionString;
            double                                      m_PulsarBirthSpinPeriodDistributionMin;                         // Minimum birth spin period (ms)
            double                                      m_PulsarBirthSpinPeriodDistributionMax;                         // Maximum birth spin period (ms)

            double                                      m_PulsarMagneticFieldDecayTimescale;                            // Timescale on which magnetic field decays (Myrs)
            double                                      m_PulsarMagneticFieldDecayMassscale;                            // Mass scale on which magnetic field decays during accretion (solar masses)
            double                                      m_PulsarLog10MinimumMagneticField;                              // log10 of the minimum pulsar magnetic field in Gauss


            // Rotational Velocity distribution options
            ROTATIONAL_VELOCITY_DISTRIBUTION            m_RotationalVelocityDistribution;                               // Rotational velocity distribution
            string                                      m_RotationalVelocityDistributionString;


	        // grids

            string                                      m_GridFilename;                                                 // Grid filename


            // debug and logging options

            int                                         m_DebugLevel;                                                   // Debug level - used to determine which debug statements are actually written
            vector<string>                              m_DebugClasses;                                                 // Debug classes - used to determine which debug statements are actually written

            int                                         m_LogLevel;                                                     // Logging level - used to determine which logging statements are actually written
            vector<string>                              m_LogClasses;                                                   // Logging classes - used to determine which logging statements are actually written


            // Logfiles
            string                                      m_LogfileDefinitionsFilename;                                   // Filename for the logfile record definitions
            DELIMITER                                   m_LogfileDelimiter;                                             // Field delimiter for log file records
            string                                      m_LogfileDelimiterString;                                       // Field delimiter for log file records (program option string)
            string                                      m_LogfileNamePrefix;                                            // Prefix for log file names

            string                                      m_LogfileSystemParameters;                                      // output file name: system parameters
            string                                      m_LogfileDetailedOutput;                                        // output file name: detailed output
            string                                      m_LogfileDoubleCompactObjects;                                  // output file name: double compact objects
            string                                      m_LogfileSupernovae;                                            // output file name: supernovae
            string                                      m_LogfileCommonEnvelopes;                                       // output file name: common envelopes
            string                                      m_LogfileRLOFParameters;                                        // output file name: Roche Lobe overflow
            string                                      m_LogfileBeBinaries;                                            // output file name: Be Binaries
            string                                      m_LogfilePulsarEvolution;                                       // output file name: pulsar evolution
            string                                      m_LogfileSwitchLog;                                             // output file name: switch log


            // the boost variables map
            // this holds information on the options as specified by the user

            po::variables_map m_VM;


            // member functions

            std::string CheckAndSetOptions(/*const po::variables_map p_VM*/);

            void        Initialise();

            int         OptionSpecified(std::string p_OptionString);


        public:

    };  // class OptionValues


    // complex option values are values for options that the user has entered as ranges or sets
    //
    // complex option valies are described by a tuple containing:
    //
    //     optionName       (std::string)               the name of the option
    //     complexValue     (RangeOrSetDescriptorT)     the complex option value - a RANGE or a SET
    //
    // ranges and sets are described by the RangeOrSetDescriptorT struct (see below)
    // the struct elements are described as:
    //
    //     type         (INT)                       type indicates whether the entry refers to a RANGE (type 0) or SET (type 1)
    //     parameters   (std::vector<std::string>)  a vector of strings that hold the parameters as they were supplied by the user
    //                                              for a RANGE there must be exactly 3 parameters: start, count, increment
    //                                              a SET must have at least one parameter (element); there is no maximum number of elements
    //     currPos      (INT)                       the current iterator position (the code iterates over the range or set)

    typedef struct RangeOrSetDescriptor {
        int                      type;          // type = 0 = RANGE, 1 = SET
        std::vector<std::string> parameters;    // the range or set parameters
        int                      currPos;       // current position of iterator                                              
    } RangeOrSetDescriptorT;

    typedef std::vector<std::tuple<std::string, RangeOrSetDescriptorT>> COMPLEX_OPTION_VALUES;


    // we have two structs:
    //    one for the program-level options, and 
    //    one for the evolving object-level options
    //
    // each struct contains:
    //
    //    an OptionValues object - holds the values of the options 
    //    a  COMPLEX_OPTION_VALUES object - holds the complex option values (ranges, sets)

    typedef struct OptionsDescriptor {
        OptionValues            optionValues;
        po::options_description optionDescriptions;
        COMPLEX_OPTION_VALUES   complexOptionValues;
    } OptionsDescriptorT;


// class Options

private:

    Options() {};
    Options(Options const&) = delete;
    Options& operator = (Options const&) = delete;

    static Options* m_Instance;

    std::string m_OptionsDetails;


    PROGRAM_STATUS ParseCommandLineOptions(int argc, char * argv[]);

    
    std::tuple<TYPENAME, bool, std::string, std::string> OptionAttributes(const po::variables_map p_VM, const po::variables_map::const_iterator p_IT);
    string ProgramOptionDetails(const OptionValues *p_Options, const po::variables_map p_VM);

    bool SetProgramOptions(OptionValues *p_Options, po::options_description *p_OptionsDescription);
    bool SetObjectOptions(OptionValues *p_Options, po::options_description *p_OptionsDescription);

    std::string ParseOptionValues(int p_ArgCount, char *p_ArgStrings[], OptionsDescriptorT &p_OptionsDescriptor);

int GetProgramOptionValues();

    OptionsDescriptorT m_Program;
    OptionsDescriptorT m_EvolvingObject;


    OptionsDescriptorT *m_Opts;



    GridfileT m_Gridfile = {"", ERROR::EMPTY_FILENAME};




    std::vector<std::string> m_RangeExcluded = {
        "quiet", 
        "log-level", 
        "debug-level"
    }; // add others...
    
    std::vector<std::string> m_SetExcluded = {
        "quiet"
    }; // add log file names etc... debug class, log classes


public:

    static Options* Instance();


    int             ApplyCmdLineOptions(int p_ArgCount = 0, char *p_ArgStrings[] = NULL);

    int             ApplyNextGridRecord();

    void            CloseGridFile() { m_Gridfile.handle.close(); m_Gridfile.filename = ""; m_Gridfile.error = ERROR::EMPTY_FILENAME; }

    bool            Initialise(int p_OptionCount, char *p_OptionStrings[]);
    bool            InitialiseObject(const std::string p_OptionsString);


    ERROR           OpenGridFile(const std::string p_GridFilename);
    int             OptionSpecified(std::string p_OptionString);

    COMPAS_VARIABLE OptionValue(const T_ANY_PROPERTY p_Property) const;


    // getters
    // m_Opts points to the OptionsVariables object in use (program or object)
    // will always be program if no grid file in use
    // will always be object after the first grid file record is read when a grid file is in use

    AIS_DCO                                     AIS_DCOType() const                                                     { return m_Opts->optionValues.m_AISDCOtype; }
    string                                      AIS_DCOTypeString() const                                               { return m_Opts->optionValues.m_AISDCOtypeString; }
    bool                                        AIS_ExploratoryPhase() const                                            { return m_Opts->optionValues.m_AISexploratoryPhase; }
    bool                                        AIS_Hubble() const                                                      { return m_Opts->optionValues.m_AIShubble; }
    bool                                        AIS_Pessimistic() const                                                 { return m_Opts->optionValues.m_AISpessimistic; }
    bool                                        AIS_RefinementPhase() const                                             { return m_Opts->optionValues.m_AISrefinementPhase; }
    bool                                        AIS_RLOF() const                                                        { return m_Opts->optionValues.m_AISrlof; }

    bool                                        AllowMainSequenceStarToSurviveCommonEnvelope() const                    { return m_Opts->optionValues.m_AllowMainSequenceStarToSurviveCommonEnvelope; }
    bool                                        AllowRLOFAtBirth() const                                                { return m_Opts->optionValues.m_AllowRLOFAtBirth; }
    bool                                        AllowTouchingAtBirth() const                                            { return m_Opts->optionValues.m_AllowTouchingAtBirth; }
    bool                                        AngularMomentumConservationDuringCircularisation() const                { return m_Opts->optionValues.m_AngularMomentumConservationDuringCircularisation; }

    bool                                        BeBinaries() const                                                      { return m_Opts->optionValues.m_BeBinaries; }

    BLACK_HOLE_KICK_OPTION                      BlackHoleKicksOption() const                                            { return m_Opts->optionValues.m_BlackHoleKicksOption; }

    EVOLUTION_MODE                              EvolutionMode() const                                                   { return m_Opts->optionValues.m_EvolutionMode; }
    
    CASE_BB_STABILITY_PRESCRIPTION              CaseBBStabilityPrescription() const                                     { return m_Opts->optionValues.m_CaseBBStabilityPrescription; }
    
    CHE_OPTION                                  CHE_Option() const                                                      { return m_Opts->optionValues.m_CheOption; }

    bool                                        CirculariseBinaryDuringMassTransfer() const                             { return m_Opts->optionValues.m_CirculariseBinaryDuringMassTransfer; }

    double                                      CommonEnvelopeAlpha() const                                             { return m_Opts->optionValues.m_CommonEnvelopeAlpha; }
    double                                      CommonEnvelopeAlphaThermal() const                                      { return m_Opts->optionValues.m_CommonEnvelopeAlphaThermal; }
    double                                      CommonEnvelopeLambda() const                                            { return m_Opts->optionValues.m_CommonEnvelopeLambda; }
    double                                      CommonEnvelopeLambdaMultiplier() const                                  { return m_Opts->optionValues.m_CommonEnvelopeLambdaMultiplier; }
    CE_LAMBDA_PRESCRIPTION                      CommonEnvelopeLambdaPrescription() const                                { return m_Opts->optionValues.m_CommonEnvelopeLambdaPrescription; }
    double                                      CommonEnvelopeMassAccretionConstant() const                             { return m_Opts->optionValues.m_CommonEnvelopeMassAccretionConstant; }
    double                                      CommonEnvelopeMassAccretionMax() const                                  { return m_Opts->optionValues.m_CommonEnvelopeMassAccretionMax; }
    double                                      CommonEnvelopeMassAccretionMin() const                                  { return m_Opts->optionValues.m_CommonEnvelopeMassAccretionMin; }
    CE_ACCRETION_PRESCRIPTION                   CommonEnvelopeMassAccretionPrescription() const                         { return m_Opts->optionValues.m_CommonEnvelopeMassAccretionPrescription; }
    double                                      CommonEnvelopeRecombinationEnergyDensity() const                        { return m_Opts->optionValues.m_CommonEnvelopeRecombinationEnergyDensity; }
    double                                      CommonEnvelopeSlopeKruckow() const                                      { return m_Opts->optionValues.m_CommonEnvelopeSlopeKruckow; }

    vector<string>                              DebugClasses() const                                                    { return m_Opts->optionValues.m_DebugClasses; }
    int                                         DebugLevel() const                                                      { return m_Opts->optionValues.m_DebugLevel; }
    bool                                        DebugToFile() const                                                     { return m_Opts->optionValues.m_DebugToFile; }
    bool                                        DetailedOutput() const                                                  { return m_Opts->optionValues.m_DetailedOutput; }

    bool                                        EnableWarnings() const                                                  { return m_Opts->optionValues.m_EnableWarnings; }
    bool                                        ErrorsToFile() const                                                    { return m_Opts->optionValues.m_ErrorsToFile; }
    double                                      Eccentricity() const                                                    { return m_Opts->optionValues.m_Eccentricity; }
    ECCENTRICITY_DISTRIBUTION                   EccentricityDistribution() const                                        { return m_Opts->optionValues.m_EccentricityDistribution; }
    double                                      EccentricityDistributionMax() const                                     { return m_Opts->optionValues.m_EccentricityDistributionMax; }
    double                                      EccentricityDistributionMin() const                                     { return m_Opts->optionValues.m_EccentricityDistributionMin; }
    double                                      EddingtonAccretionFactor() const                                        { return m_Opts->optionValues.m_EddingtonAccretionFactor; }
    ENVELOPE_STATE_PRESCRIPTION                 EnvelopeStatePrescription() const                                       { return m_Opts->optionValues.m_EnvelopeStatePrescription; }
    bool                                        EvolvePulsars() const                                                   { return m_Opts->optionValues.m_EvolvePulsars; }
    bool                                        EvolveUnboundSystems() const                                            { return m_Opts->optionValues.m_EvolveUnboundSystems; }

    bool                                        FixedMetallicity() const                                                { return m_Opts->optionValues.m_FixedMetallicity; }
    bool                                        FixedRandomSeed() const                                                 { return m_Opts->optionValues.m_FixedRandomSeed; }
    double                                      FixedUK() const                                                         { return m_Opts->optionValues.m_FixedUK; }                 // JR: todo: this isn't consistent naming see fixedMetallicity, fixedRandomSeed)
    SN_ENGINE                                   FryerSupernovaEngine() const                                            { return m_Opts->optionValues.m_FryerSupernovaEngine; }

    string                                      GridFilename() const                                                    { return m_Opts->optionValues.m_GridFilename; }

    double                                      InitialMass() const                                                     { return m_Opts->optionValues.m_InitialMass; }
    double                                      InitialMass1() const                                                    { return m_Opts->optionValues.m_InitialMass1; }
    double                                      InitialMass2() const                                                    { return m_Opts->optionValues.m_InitialMass2; }

    INITIAL_MASS_FUNCTION                       InitialMassFunction() const                                             { return m_Opts->optionValues.m_InitialMassFunction; }
    double                                      InitialMassFunctionMax() const                                          { return m_Opts->optionValues.m_InitialMassFunctionMax; }
    double                                      InitialMassFunctionMin() const                                          { return m_Opts->optionValues.m_InitialMassFunctionMin; }
    double                                      InitialMassFunctionPower() const                                        { return m_Opts->optionValues.m_InitialMassFunctionPower; }

    KICK_DIRECTION_DISTRIBUTION                 KickDirectionDistribution() const                                       { return m_Opts->optionValues.m_KickDirectionDistribution; }
    double                                      KickDirectionPower() const                                              { return m_Opts->optionValues.m_KickDirectionPower; }
    double                                      KickScalingFactor() const                                               { return m_Opts->optionValues.m_KickScalingFactor; }
    KICK_MAGNITUDE_DISTRIBUTION                 KickMagnitudeDistribution() const                                       { return m_Opts->optionValues.m_KickMagnitudeDistribution; }

    double                                      KickMagnitudeDistributionMaximum() const                                { return m_Opts->optionValues.m_KickMagnitudeDistributionMaximum; }

    double                                      KickMagnitudeDistributionSigmaCCSN_BH() const                           { return m_Opts->optionValues.m_KickMagnitudeDistributionSigmaCCSN_BH; }
    double                                      KickMagnitudeDistributionSigmaCCSN_NS() const                           { return m_Opts->optionValues.m_KickMagnitudeDistributionSigmaCCSN_NS; }
    double                                      KickMagnitudeDistributionSigmaForECSN() const                           { return m_Opts->optionValues.m_KickMagnitudeDistributionSigmaForECSN; }
    double                                      KickMagnitudeDistributionSigmaForUSSN() const                           { return m_Opts->optionValues.m_KickMagnitudeDistributionSigmaForUSSN; }

    double                                      KickMagnitude() const                                                   { return m_Opts->optionValues.m_KickMagnitude; }
    double                                      KickMagnitude1() const                                                  { return m_Opts->optionValues.m_KickMagnitude1; }
    double                                      KickMagnitude2() const                                                  { return m_Opts->optionValues.m_KickMagnitude2; }

    double                                      KickMagnitudeRandom() const                                             { return m_Opts->optionValues.m_KickMagnitudeRandom; }
    double                                      KickMagnitudeRandom1() const                                            { return m_Opts->optionValues.m_KickMagnitudeRandom1; }
    double                                      KickMagnitudeRandom2() const                                            { return m_Opts->optionValues.m_KickMagnitudeRandom2; }

    double                                      SN_MeanAnomaly1() const                                                 { return m_Opts->optionValues.m_KickMeanAnomaly1; }
    double                                      SN_MeanAnomaly2() const                                                 { return m_Opts->optionValues.m_KickMeanAnomaly2; }
    double                                      SN_Phi1() const                                                         { return m_Opts->optionValues.m_KickPhi1; }
    double                                      SN_Phi2() const                                                         { return m_Opts->optionValues.m_KickPhi2; }
    double                                      SN_Theta1() const                                                       { return m_Opts->optionValues.m_KickTheta1; }
    double                                      SN_Theta2() const                                                       { return m_Opts->optionValues.m_KickTheta2; }

    vector<string>                              LogClasses() const                                                      { return m_Opts->optionValues.m_LogClasses; }
    string                                      LogfileBeBinaries() const                                               { return m_Opts->optionValues.m_LogfileBeBinaries; }
    string                                      LogfileCommonEnvelopes() const                                          { return m_Opts->optionValues.m_LogfileCommonEnvelopes; }
    string                                      LogfileDefinitionsFilename() const                                      { return m_Opts->optionValues.m_LogfileDefinitionsFilename; }
    DELIMITER                                   LogfileDelimiter() const                                                { return m_Opts->optionValues.m_LogfileDelimiter; }
    string                                      LogfileDelimiterString() const                                          { return m_Opts->optionValues.m_LogfileDelimiterString; }
    string                                      LogfileDetailedOutput() const                                           { return m_Opts->optionValues.m_LogfileDetailedOutput; }
    string                                      LogfileDoubleCompactObjects() const                                     { return m_Opts->optionValues.m_LogfileDoubleCompactObjects; }
    string                                      LogfileNamePrefix() const                                               { return m_Opts->optionValues.m_LogfileNamePrefix; }
    string                                      LogfilePulsarEvolution() const                                          { return m_Opts->optionValues.m_LogfilePulsarEvolution; }
    string                                      LogfileRLOFParameters() const                                           { return m_Opts->optionValues.m_LogfileRLOFParameters; }
    string                                      LogfileSupernovae() const                                               { return m_Opts->optionValues.m_LogfileSupernovae; }
    string                                      LogfileSwitchLog() const                                                { return m_Opts->optionValues.m_LogfileSwitchLog; }
    string                                      LogfileSystemParameters() const                                         { return m_Opts->optionValues.m_LogfileSystemParameters; }
    int                                         LogLevel() const                                                        { return m_Opts->optionValues.m_LogLevel; }

    double                                      LuminousBlueVariableFactor() const                                      { return m_Opts->optionValues.m_LuminousBlueVariableFactor; }

    MASS_LOSS_PRESCRIPTION                      MassLossPrescription() const                                            { return m_Opts->optionValues.m_MassLossPrescription; }

    MASS_RATIO_DISTRIBUTION                     MassRatioDistribution() const                                           { return m_Opts->optionValues.m_MassRatioDistribution; }
    double                                      MassRatioDistributionMax() const                                        { return m_Opts->optionValues.m_MassRatioDistributionMax; }
    double                                      MassRatioDistributionMin() const                                        { return m_Opts->optionValues.m_MassRatioDistributionMin; }

    MT_ACCRETION_EFFICIENCY_PRESCRIPTION        MassTransferAccretionEfficiencyPrescription() const                     { return m_Opts->optionValues.m_MassTransferAccretionEfficiencyPrescription; }
    MT_ANGULAR_MOMENTUM_LOSS_PRESCRIPTION       MassTransferAngularMomentumLossPrescription() const                     { return m_Opts->optionValues.m_MassTransferAngularMomentumLossPrescription; }
    double                                      MassTransferCParameter() const                                          { return m_Opts->optionValues.m_MassTransferCParameter; }

    //
    bool                                        MassTransferCriticalMassRatioGiant() const                              { return m_Opts->optionValues.m_MassTransferCriticalMassRatioGiant; }
    double                                      MassTransferCriticalMassRatioGiantDegenerateAccretor() const            { return m_Opts->optionValues.m_MassTransferCriticalMassRatioGiantDegenerateAccretor; }
    double                                      MassTransferCriticalMassRatioGiantNonDegenerateAccretor() const         { return m_Opts->optionValues.m_MassTransferCriticalMassRatioGiantNonDegenerateAccretor; }
    bool                                        MassTransferCriticalMassRatioHeliumGiant() const                        { return m_Opts->optionValues.m_MassTransferCriticalMassRatioHeliumGiant; }
    double                                      MassTransferCriticalMassRatioHeliumGiantDegenerateAccretor() const      { return m_Opts->optionValues.m_MassTransferCriticalMassRatioHeliumGiantDegenerateAccretor; }
    double                                      MassTransferCriticalMassRatioHeliumGiantNonDegenerateAccretor() const   { return m_Opts->optionValues.m_MassTransferCriticalMassRatioHeliumGiantNonDegenerateAccretor; }
    bool                                        MassTransferCriticalMassRatioHeliumHG() const                           { return m_Opts->optionValues.m_MassTransferCriticalMassRatioHeliumHG; }
    double                                      MassTransferCriticalMassRatioHeliumHGDegenerateAccretor() const         { return m_Opts->optionValues.m_MassTransferCriticalMassRatioHeliumHGDegenerateAccretor; }
    double                                      MassTransferCriticalMassRatioHeliumHGNonDegenerateAccretor() const      { return m_Opts->optionValues.m_MassTransferCriticalMassRatioHeliumHGNonDegenerateAccretor; }
    bool                                        MassTransferCriticalMassRatioHeliumMS() const                           { return m_Opts->optionValues.m_MassTransferCriticalMassRatioHeliumMS; }
    double                                      MassTransferCriticalMassRatioHeliumMSDegenerateAccretor() const         { return m_Opts->optionValues.m_MassTransferCriticalMassRatioHeliumMSDegenerateAccretor; }
    double                                      MassTransferCriticalMassRatioHeliumMSNonDegenerateAccretor() const      { return m_Opts->optionValues.m_MassTransferCriticalMassRatioHeliumMSNonDegenerateAccretor; }
    bool                                        MassTransferCriticalMassRatioHG() const                                 { return m_Opts->optionValues.m_MassTransferCriticalMassRatioHG; }
    double                                      MassTransferCriticalMassRatioHGDegenerateAccretor() const               { return m_Opts->optionValues.m_MassTransferCriticalMassRatioHGDegenerateAccretor; }
    double                                      MassTransferCriticalMassRatioHGNonDegenerateAccretor() const            { return m_Opts->optionValues.m_MassTransferCriticalMassRatioHGNonDegenerateAccretor; }
    bool                                        MassTransferCriticalMassRatioMSHighMass() const                         { return m_Opts->optionValues.m_MassTransferCriticalMassRatioMSHighMass; }
    double                                      MassTransferCriticalMassRatioMSHighMassDegenerateAccretor() const       { return m_Opts->optionValues.m_MassTransferCriticalMassRatioMSHighMassDegenerateAccretor; }
    double                                      MassTransferCriticalMassRatioMSHighMassNonDegenerateAccretor() const    { return m_Opts->optionValues.m_MassTransferCriticalMassRatioMSHighMassNonDegenerateAccretor; }
    bool                                        MassTransferCriticalMassRatioMSLowMass() const                          { return m_Opts->optionValues.m_MassTransferCriticalMassRatioMSLowMass; }
    double                                      MassTransferCriticalMassRatioMSLowMassDegenerateAccretor() const        { return m_Opts->optionValues.m_MassTransferCriticalMassRatioMSLowMassDegenerateAccretor; }
    double                                      MassTransferCriticalMassRatioMSLowMassNonDegenerateAccretor() const     { return m_Opts->optionValues.m_MassTransferCriticalMassRatioMSLowMassNonDegenerateAccretor; }
    bool                                        MassTransferCriticalMassRatioWhiteDwarf() const                         { return m_Opts->optionValues.m_MassTransferCriticalMassRatioWhiteDwarf; }

    double                                      MassTransferFractionAccreted() const                                    { return m_Opts->optionValues.m_MassTransferFractionAccreted; }
    double                                      MassTransferJloss() const                                               { return m_Opts->optionValues.m_MassTransferJloss; }
    MT_REJUVENATION_PRESCRIPTION                MassTransferRejuvenationPrescription() const                            { return m_Opts->optionValues.m_MassTransferRejuvenationPrescription; }
    MT_THERMALLY_LIMITED_VARIATION              MassTransferThermallyLimitedVariation() const                           { return m_Opts->optionValues.m_MassTransferThermallyLimitedVariation; }
    double                                      MaxEvolutionTime() const                                                { return m_Opts->optionValues.m_MaxEvolutionTime; }
    int                                         MaxNumberOfTimestepIterations() const                                   { return m_Opts->optionValues.m_MaxNumberOfTimestepIterations; }
    double                                      MaxPercentageAdaptiveMassTransfer() const                               { return m_Opts->optionValues.m_MaxPercentageAdaptiveMassTransfer; }

    double                                      MCBUR1() const                                                          { return m_Opts->optionValues.m_mCBUR1; }

    double                                      Metallicity() const                                                     { return m_Opts->optionValues.m_Metallicity; }

    double                                      MinimumMassSecondary() const                                            { return m_Opts->optionValues.m_MinimumMassSecondary; }
    double                                      MaximumNeutronStarMass() const                                          { return m_Opts->optionValues.m_MaximumNeutronStarMass; }

    int                                         nBatchesUsed() const                                                    { return m_Opts->optionValues.m_nBatchesUsed; }

    NEUTRINO_MASS_LOSS_PRESCRIPTION             NeutrinoMassLossAssumptionBH() const                                    { return m_Opts->optionValues.m_NeutrinoMassLossAssumptionBH; }
    double                                      NeutrinoMassLossValueBH() const                                         { return m_Opts->optionValues.m_NeutrinoMassLossValueBH; }

    NS_EOS                                      NeutronStarEquationOfState() const                                      { return m_Opts->optionValues.m_NeutronStarEquationOfState; }

    int                                         nObjectsToEvolve() const                                                { return m_Opts->optionValues.m_ObjectsToEvolve; }
    bool                                        OptimisticCHE() const                                                   { return m_Opts->optionValues.m_CheOption == CHE_OPTION::OPTIMISTIC; }

    string                                      OptionsDetails() const                                                  { return m_OptionsDetails; }

    string                                      OutputContainerName() const                                             { return m_Opts->optionValues.m_OutputContainerName; }
    string                                      OutputPathString() const                                                { return m_Opts->optionValues.m_OutputPath.string(); }

    double                                      PairInstabilityLowerLimit() const                                       { return m_Opts->optionValues.m_PairInstabilityLowerLimit; }
    double                                      PairInstabilityUpperLimit() const                                       { return m_Opts->optionValues.m_PairInstabilityUpperLimit; }

    double                                      PeriodDistributionMax() const                                           { return m_Opts->optionValues.m_PeriodDistributionMax; }
    double                                      PeriodDistributionMin() const                                           { return m_Opts->optionValues.m_PeriodDistributionMin; }

    bool                                        PopulationDataPrinting() const                                          { return m_Opts->optionValues.m_PopulationDataPrinting; }
    bool                                        PrintBoolAsString() const                                               { return m_Opts->optionValues.m_PrintBoolAsString; }

    PULSAR_BIRTH_MAGNETIC_FIELD_DISTRIBUTION    PulsarBirthMagneticFieldDistribution() const                            { return m_Opts->optionValues.m_PulsarBirthMagneticFieldDistribution; }
    double                                      PulsarBirthMagneticFieldDistributionMax() const                         { return m_Opts->optionValues.m_PulsarBirthMagneticFieldDistributionMax; }
    double                                      PulsarBirthMagneticFieldDistributionMin() const                         { return m_Opts->optionValues.m_PulsarBirthMagneticFieldDistributionMin; }

    PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION       PulsarBirthSpinPeriodDistribution() const                               { return m_Opts->optionValues.m_PulsarBirthSpinPeriodDistribution; }
    double                                      PulsarBirthSpinPeriodDistributionMax() const                            { return m_Opts->optionValues.m_PulsarBirthSpinPeriodDistributionMax; }
    double                                      PulsarBirthSpinPeriodDistributionMin() const                            { return m_Opts->optionValues.m_PulsarBirthSpinPeriodDistributionMin; }

    double                                      PulsarLog10MinimumMagneticField() const                                 { return m_Opts->optionValues.m_PulsarLog10MinimumMagneticField; }

    double                                      PulsarMagneticFieldDecayMassscale() const                               { return m_Opts->optionValues.m_PulsarMagneticFieldDecayMassscale; }
    double                                      PulsarMagneticFieldDecayTimescale() const                               { return m_Opts->optionValues.m_PulsarMagneticFieldDecayTimescale; }

    PPI_PRESCRIPTION                            PulsationalPairInstabilityPrescription() const                          { return m_Opts->optionValues.m_PulsationalPairInstabilityPrescription; }
    double                                      PulsationalPairInstabilityLowerLimit() const                            { return m_Opts->optionValues.m_PulsationalPairInstabilityLowerLimit; }
    double                                      PulsationalPairInstabilityUpperLimit() const                            { return m_Opts->optionValues.m_PulsationalPairInstabilityUpperLimit; }

    bool                                        Quiet() const                                                           { return m_Opts->optionValues.m_Quiet; }

    unsigned long int                           RandomSeed() const                                                      { return m_Opts->optionValues.m_RandomSeed; }

    REMNANT_MASS_PRESCRIPTION                   RemnantMassPrescription() const                                         { return m_Opts->optionValues.m_RemnantMassPrescription; }
    bool                                        RLOFPrinting() const                                                    { return m_Opts->optionValues.m_RlofPrinting; }

    ROTATIONAL_VELOCITY_DISTRIBUTION            RotationalVelocityDistribution() const                                  { return m_Opts->optionValues.m_RotationalVelocityDistribution; }
   
    double                                      SemiMajorAxis() const                                                   { return m_Opts->optionValues.m_SemiMajorAxis; }
    SEMI_MAJOR_AXIS_DISTRIBUTION                SemiMajorAxisDistribution() const                                       { return m_Opts->optionValues.m_SemiMajorAxisDistribution; }
    double                                      SemiMajorAxisDistributionMax() const                                    { return m_Opts->optionValues.m_SemiMajorAxisDistributionMax; }
    double                                      SemiMajorAxisDistributionMin() const                                    { return m_Opts->optionValues.m_SemiMajorAxisDistributionMin; }
    double                                      SemiMajorAxisDistributionPower() const                                  { return m_Opts->optionValues.m_SemiMajorAxisDistributionPower; }

    bool                                        RequestedHelp() const                                                   { return m_Program.optionValues.m_VM["help"].as<bool>(); }
    bool                                        RequestedVersion() const                                                { return m_Program.optionValues.m_VM["version"].as<bool>(); }

    bool                                        SwitchLog() const                                                       { return m_Opts->optionValues.m_SwitchLog; }

    ZETA_PRESCRIPTION                           StellarZetaPrescription() const                                         { return m_Opts->optionValues.m_StellarZetaPrescription; }

    bool                                        UseFixedUK() const                                                      { return m_Opts->optionValues.m_UseFixedUK; }
    bool                                        UseMassLoss() const                                                     { return m_Opts->optionValues.m_UseMassLoss; }
    bool                                        UseMassTransfer() const                                                 { return m_Opts->optionValues.m_UseMassTransfer; }
    bool                                        UsePairInstabilitySupernovae() const                                    { return m_Opts->optionValues.m_UsePairInstabilitySupernovae; }
    bool                                        UsePulsationalPairInstability() const                                   { return m_Opts->optionValues.m_UsePulsationalPairInstability; }

    double                                      WolfRayetFactor() const                                                 { return m_Opts->optionValues.m_WolfRayetFactor; }

    double                                      ZetaRadiativeEnvelopeGiant() const                                      { return m_Opts->optionValues.m_ZetaRadiativeEnvelopeGiant; }
    double                                      ZetaMainSequence() const                                                { return m_Opts->optionValues.m_ZetaMainSequence; }
    double                                      ZetaAdiabaticArbitrary() const                                          { return m_Opts->optionValues.m_ZetaAdiabaticArbitrary; }

};


#endif // __Options_H__
