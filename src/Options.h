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


#define OPT_VALUE(optName, optValue)    m_GridLine.optionValues.m_VM[optName].defaulted() ? m_CmdLine.optionValues.optValue : m_GridLine.optionValues.optValue

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

private:

    // The m_SSEOnly vector records option strings that apply to SSE only
    // The m_BSEOnly vector records option strings that apply to BSE only
    //
    // These vectors are checked when the commandline or grid line are parsed for
    // ranges and sets.  Ranges and sets are played out, and stars/binaries evolved
    // based on the grid of options defined by any ranges and sets specified by the
    // user.
    //
    // A problem arises when a user is (say) evolving single stars (in SSE mode)
    // and (perhaps inadvertently) specifies a range or set for an option that
    // applies only to BSE.  Before ranges and sets were implemented, this was not
    // a problem - the SSE code just ignores any BSE-only options (and vice-versa).
    // But with ranges and sets implemented, we need to know whether we should play
    // out the range or set.  If a range for BSE-only only option is specified the
    // SSE code will happily ignore the option, but unless we know not to, we will
    // still play out the range of values (only for them to be ignored) - so we will
    // evolve as many stars as there are values in the range, and they will all be
    // the same (because the SSE code will ignore the BSE-only option each time 
    // through the loop while the range is playing out).
    //
    // To get around this we specify in the following vectors the names of any
    // options that are SSE only or BSE only - that way we can choose not to play
    // them out as required.
    //
    // It's not the end of the world if we forget to put some entries in these
    // vectors - the worst thing that will happen is that duplicate stars/binaries
    // will be evolved as ranges/sets of options that will be ignored are played out.
    // In that case the user has a simple remedy: don't specify ranges for BSE options
    // when evolving in SSE mode (and vice-versa).  All we're doing here with these
    // vectors is helping the user avoid duplicating stars/binaries if they specify
    // inconsistent options.


// JR NEED TO IMPLEMENT CODE FOR THIS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    std::vector<std::string> m_SSEOnly = {
        "initial-mass",
        "kick-magnitude",
        "kick-magnitude-random",
    };

    std::vector<std::string> m_BSEOnly = {
        "initial-mass-1",
        "initial-mass-2",
        "semi-major-axis", "a",

        // Floor
        /*
        "ais-dcotype",
        "ais-exploratory-phase",
        "ais-hubble",
        "ais-pessimistic",
        "ais-refinement-phase",
        "ais-rlof",
        "kappa-gaussians",
        "nbatches-used",
        */

        "allow-rlof-at-birth",
        "allow-touching-at-birth",
        "angular-momentum-conservation-during-circularisation", 

        // Serena
        //"be-binaries",

        "circularise-binary-during-mass-transfer",
        "common-envelope-allow-main-sequence-survive",

        "common-envelope-alpha", 
        "common-envelope-alpha-thermal",
        "common-envelope-lambda",
        "common-envelope-lambda-multiplier",
        "common-envelope-mass-accretion-constant",
        "common-envelope-mass-accretion-max",
        "common-envelope-mass-accretion-min",
        "common-envelope-recombination-energy-density",
        "common-envelope-slope-kruckow",

        // AVG
        /*
        "critical-mass-ratio-giant-degenerate-accretor",
        "critical-mass-ratio-giant-non-degenerate-accretor",
        "critical-mass-ratio-helium-giant-degenerate-accretor",
        "critical-mass-ratio-helium-giant-non-degenerate-accretor",
        "critical-mass-ratio-helium-hg-degenerate-accretor",
        "critical-mass-ratio-helium-hg-non-degenerate-accretor",
        "critical-mass-ratio-helium-ms-degenerate-accretor",
        "critical-mass-ratio-helium-ms-non-degenerate-accretor",
        "critical-mass-ratio-hg-degenerate-accretor",
        "critical-mass-ratio-hg-non-degenerate-accretor",
        "critical-mass-ratio-ms-high-mass-degenerate-accretor",
        "critical-mass-ratio-ms-high-mass-non-degenerate-accretor",
        "critical-mass-ratio-ms-low-mass-degenerate-accretor",
        "critical-mass-ratio-ms-low-mass-non-degenerate-accretor",
        "critical-mass-ratio-white-dwarf-degenerate-accretor",
        "critical-mass-ratio-white-dwarf-non-degenerate-accretor",
        */

        "kick-magnitude-1",
        "kick-magnitude-2",
        "kick-magnitude-random-1",
        "kick-magnitude-random-2",
        "kick-mean-anomaly-1",
        "kick-mean-anomaly-2",
        "kick-phi-1",
        "kick-phi-2",
        "kick-theta-1",
        "kick-theta-2",

        "mass-ratio-max",
        "mass-ratio-min",

        "mass-transfer-fa",
        "mass-transfer-jloss",
        "mass-transfer-thermal-limit-c",
        "maximum-mass-donor-nandez-ivanova",

        "minimum-secondary-mass",

        "orbital-period-max",
        "orbital-period-min",

        "semi-major-axis", "a",
        "semi-major-axis-max",
        "semi-major-axis-min",

        "case-bb-stability-prescription",

        "common-envelope-lambda-prescription",
        "common-envelope-mass-accretion-prescription",

        "eccentricity-distribution",

        "mass-ratio-distribution", "q",

        "mass-transfer-accretion-efficiency-prescription",
        "mass-transfer-angular-momentum-loss-prescription",
        "mass-transfer-rejuvenation-prescription",
        "mass-transfer-thermal-limit-accretor",

        "evolve-pulsars",
        "evolve-unbound-systems",

        "mass-transfer",

        "rlof-printing",

        "logfile-rlof-parameters",
        "logfile-common-envelopes",
        "logfile-double-compact-objects",
        "logfile-pulsar-evolution",
        "logfile-system-parameters",
    };

    std::vector<std::string> m_RangeExcluded = {
        "help", "h",
        "version", "v",

        // Floor
        /*
        "ais-dcotype",
        "ais-exploratory-phase",
        "ais-hubble",
        "ais-pessimistic",
        "ais-refinement-phase",
        "ais-rlof",
        "kappa-gaussians",
        "nbatches-used",
        */

        "allow-rlof-at-birth",
        "allow-touching-at-birth",
        "angular-momentum-conservation-during-circularisation",

        // Serena
        //"be-binaries",

        "circularise-binary-during-mass-transfer",
        "common-envelope-allow-main-sequence-survive",

        "evolve-pulsars",
        "evolve-unbound-systems",

        "mass-transfer",

        "pair-instability-supernovae",
        "pulsational-pair-instability",

        "revised-energy-formalism-nandez-ivanova",

        "use-mass-loss",

        "black-hole-kicks",

        "case-bb-stability-prescription",
        "chemically-homogeneous-evolution",
        "common-envelope-lambda-prescription",
        "common-envelope-mass-accretion-prescription",

        "eccentricity-distribution",
        "envelope-state-prescription",

        "fryer-supernova-engine",

        "initial-mass-function", "i",

        "kick-direction",
        "kick-magnitude-distribution", 

        "mass-loss-prescription",
        "mass-ratio-distribution", "q",

        "mass-transfer-accretion-efficiency-prescription",
        "mass-transfer-angular-momentum-loss-prescription",
        "mass-transfer-rejuvenation-prescription",
        "mass-transfer-thermal-limit-accretor",

        "neutrino-mass-loss-bh-formation",
        "neutron-star-equation-of-state",

        "pulsar-birth-magnetic-field-distribution",
        "pulsar-birth-spin-period-distribution",
        "pulsational-pair-instability-prescription",
        "remnant-mass-prescription",
        "rotational-velocity-distribution",
        "semi-major-axis-distribution",
        "stellar-zeta-prescription",
        
        "quiet", 
        "log-level", 
        "log-classes",
        "debug-level",
        "debug_classes",
        "debug-to-file",
        "detailedOutput",
        "enable-warnings",
        "errors-to-file",
        "populationDataPrinting",
        "print-bool-as-string",
        "rlof-printing",
        "switchlog",
        "grid",
        "mode",

        // Serena
        //"logfile-be-binaries",

        "logfile-rlof-parameters",
        "logfile-common-envelopes",
        "logfile-detailed-output",
        "logfile-double-compact-objects",
        "logfile-pulsar-evolution",
        "logfile-supernovae",
        "logfile-system-parameters",
        "logfile-switch-log",

        "logfile-definitions",
        "logfile-delimiter",
        "logfile-name-prefix",

        "output-container", "c",
        "outputPath", "o"
    };
    
    std::vector<std::string> m_SetExcluded = {
        "help", "h",
        "version", "v",

        // Floor
        /*
        "ais-dcotype",
        "ais-exploratory-phase",
        "ais-hubble",
        "ais-pessimistic",
        "ais-refinement-phase",
        "ais-rlof",
        "kappa-gaussians",
        "nbatches-used",
        */

        "quiet",
        "log-level", 
        "log-classes",
        "debug-level",
        "debug_classes",
        "debug-to-file",
        "detailedOutput",
        "enable-warnings",
        "errors-to-file",
        "populationDataPrinting",
        "print-bool-as-string",
        "rlof-printing",
        "switchlog",
        "grid",
        "mode",

        // Serena
        //"logfile-be-binaries",

        "logfile-rlof-parameters",
        "logfile-common-envelopes",
        "logfile-detailed-output",
        "logfile-double-compact-objects",
        "logfile-pulsar-evolution",
        "logfile-supernovae",
        "logfile-system-parameters",
        "logfile-switch-log",

        "logfile-definitions",
        "logfile-delimiter",
        "logfile-name-prefix",

        "output-container", "c",
        "outputPath", "o"
    };


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


            // Miscellaneous evolution variables
            EVOLUTION_MODE                              m_EvolutionMode;                                                // Mode of evolution: SSE or BSE
            string                                      m_EvolutionModeString;

            int                                         m_ObjectsToEvolve;                                              // Number of stars (SSE) or binaries (BSE) to evolve
            bool                                        m_FixedRandomSeed;                                              // Whether to use a fixed random seed given by options.randomSeed (set to true if --random-seed is passed on command line)
            unsigned long int                           m_RandomSeed;                                                   // Random seed to use
    
            double                                      m_MaxEvolutionTime;                                             // Maximum time to evolve a binary by
            int                                         m_MaxNumberOfTimestepIterations;                                // Maximum number of timesteps to evolve binary for before giving up
            double                                      m_TimestepMultiplier;                                           // Multiplier for time step size (<1 -- shorter timesteps, >1 -- longer timesteps)

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
            double                                      m_SemiMajorAxisDistributionPower;                               // Set semi-major axis distribution power law slope by hand     ** JR: there is no option for this....
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

            // AVG
            // Mass transfer critical mass ratios
            bool                                        m_MassTransferCriticalMassRatioMSLowMass;                       // Whether to use critical mass ratios
            double                                      m_MassTransferCriticalMassRatioMSLowMassNonDegenerateAccretor;  // Critical mass ratio for MT from a MS low mass star
            double                                      m_MassTransferCriticalMassRatioMSLowMassDegenerateAccretor;     // Critical mass ratio for MT from a MS low mass star on to a degenerate accretor

            bool                                        m_MassTransferCriticalMassRatioMSHighMass;                      // Whether to use critical mass ratios
            double                                      m_MassTransferCriticalMassRatioMSHighMassNonDegenerateAccretor; // Critical mass ratio for MT from a MS high mass star
            double                                      m_MassTransferCriticalMassRatioMSHighMassDegenerateAccretor;    // Critical mass ratio for MT from a MS high mass star on to a degenerate accretor

            bool                                        m_MassTransferCriticalMassRatioGiant;                           // Whether to use critical mass ratios
            double                                      m_MassTransferCriticalMassRatioGiantNonDegenerateAccretor;      // Critical mass ratio for MT from a giant
            double                                      m_MassTransferCriticalMassRatioGiantDegenerateAccretor;         // Critical mass ratio for MT from a giant on to a degenerate accretor

            bool                                        m_MassTransferCriticalMassRatioHG;                              // Whether to use critical mass ratios
            double                                      m_MassTransferCriticalMassRatioHGNonDegenerateAccretor;         // Critical mass ratio for MT from a HG star
            double                                      m_MassTransferCriticalMassRatioHGDegenerateAccretor;            // Critical mass ratio for MT from a HG star on to a degenerate accretor

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

            template<class T>
            void ModifyVariableMap(std::map<std::string, po::variable_value>& vm, const std::string& opt, const T& val) { 
                vm[opt].value() = boost::any(val);
            }

            int         OptionSpecified(std::string p_OptionString);

        public:

    };  // class OptionValues


    // complex option values are values for options that the user has supplied as ranges or sets
    //
    // complex option valies are described by a tuple containing:
    //
    //     optionName       (std::string)               the name of the option
    //     complexValue     (RangeOrSetDescriptorT)     the complex option value - a RANGE or a SET
    //
    // ranges and sets are described by the RangeOrSetDescriptorT struct (see below)
    // the struct elements are described as:
    //
    //     type         (INT)                           type indicates whether the entry refers to a RANGE (type 0) or SET (type 1)
    //     dataType     (TYPENAME)                      the data type of the option to which the RangeOrSetDescriptorT pertaines
    //     parameters   (std::vector<std::string>)      a vector of strings that hold the parameters as they were supplied by the user
    //                                                  for a RANGE there must be exactly 3 parameters: start, count, increment
    //                                                  a SET must have at least one parameter (element); there is no maximum number of elements
    //     rangeParms   (std::vector<RangeParameterT>)  numerical values for range parameters (see RangeParameter struct)
    //     currPos      (INT)                           the current iterator position (the code iterates over the range or set)

    enum class COMPLEX_TYPE: int {NONE, RANGE, SET};

    typedef union RangeParameter {
        double fVal;                                // FLOAT
        int    iVal;                                // INT
    } RangeParameterT; 

    typedef struct RangeOrSetDescriptor {
        COMPLEX_TYPE                 type;          // RANGE or SET
        TYPENAME                     dataType;      // the option datatype
        std::vector<std::string>     parameters;    // the range or set parameters
        std::vector<RangeParameterT> rangeParms;    // range parameters numerical values
        int                          currPos;       // current position of iterator - count for RANGE, pos for SET                                             
    } RangeOrSetDescriptorT;

    typedef std::vector<std::tuple<std::string, RangeOrSetDescriptorT>> COMPLEX_OPTION_VALUES;

    typedef std::tuple<TYPENAME, bool, std::string, std::string> ATTR;  // <dataType, defaulted, typeStr, valueStr>


    // we have two structs:
    //    one for the commandline (program-level) options, and 
    //    one for the grid file line (evolving object-level) options
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


    // member variables

    std::string        m_CmdLineOptionsDetails;        // for Run_Details file

    GridfileT          m_Gridfile = {"", ERROR::EMPTY_FILENAME};

    OptionsDescriptorT m_CmdLine;
    OptionsDescriptorT m_GridLine;


    // member functions

    bool           AddOptions(OptionValues *p_Options, po::options_description *p_OptionsDescription);

    ATTR           OptionAttributes(const po::variables_map p_VM, const po::variables_map::const_iterator p_IT);
    string         OptionDetails(const OptionsDescriptorT &p_Options);

    PROGRAM_STATUS ParseCommandLineOptions(int argc, char * argv[]);
    std::string    ParseOptionValues(int p_ArgCount, char *p_ArgStrings[], OptionsDescriptorT &p_OptionsDescriptor);


public:

    static Options* Instance();



    int             AdvanceGridLineOptionValues();
    int             AdvanceCmdLineOptionValues();
    int             ApplyNextGridLine();

    void            CloseGridFile() { m_Gridfile.handle.close(); m_Gridfile.filename = ""; m_Gridfile.error = ERROR::EMPTY_FILENAME; }
 
    bool            Initialise(int p_OptionCount, char *p_OptionStrings[]);
    bool            InitialiseEvolvingObject(const std::string p_OptionsString);


    ERROR           OpenGridFile(const std::string p_GridFilename);
    int             OptionSpecified(std::string p_OptionString);

    COMPAS_VARIABLE OptionValue(const T_ANY_PROPERTY p_Property) const;

    void            RewindGridFile() { m_Gridfile.handle.clear(); m_Gridfile.handle.seekg(0); }


    // getters

    AIS_DCO                                     AIS_DCOType() const                                                     { return m_CmdLine.optionValues.m_AISDCOtype; }
    string                                      AIS_DCOTypeString() const                                               { return m_CmdLine.optionValues.m_AISDCOtypeString; }
    bool                                        AIS_ExploratoryPhase() const                                            { return m_CmdLine.optionValues.m_AISexploratoryPhase; }
    bool                                        AIS_Hubble() const                                                      { return m_CmdLine.optionValues.m_AIShubble; }
    bool                                        AIS_Pessimistic() const                                                 { return m_CmdLine.optionValues.m_AISpessimistic; }
    bool                                        AIS_RefinementPhase() const                                             { return m_CmdLine.optionValues.m_AISrefinementPhase; }
    bool                                        AIS_RLOF() const                                                        { return m_CmdLine.optionValues.m_AISrlof; }

    bool                                        AllowMainSequenceStarToSurviveCommonEnvelope() const                    { return OPT_VALUE("common-envelope-allow-main-sequence-survive", m_AllowMainSequenceStarToSurviveCommonEnvelope); }
    bool                                        AllowRLOFAtBirth() const                                                { return OPT_VALUE("allow-rlof-at-birth", m_AllowRLOFAtBirth); }
    bool                                        AllowTouchingAtBirth() const                                            { return OPT_VALUE("allow-touching-at-birth", m_AllowTouchingAtBirth); }
    bool                                        AngularMomentumConservationDuringCircularisation() const                { return OPT_VALUE("angular-momentum-conservation-during-circularisation", m_AngularMomentumConservationDuringCircularisation); }

// Serena
    bool                                        BeBinaries() const                                                      { return OPT_VALUE("be-binaries", m_BeBinaries); }

    BLACK_HOLE_KICK_OPTION                      BlackHoleKicksOption() const                                            { return OPT_VALUE("black-hole-kicks", m_BlackHoleKicksOption); }

    EVOLUTION_MODE                              EvolutionMode() const                                                   { return m_CmdLine.optionValues.m_EvolutionMode; }
    
    CASE_BB_STABILITY_PRESCRIPTION              CaseBBStabilityPrescription() const                                     { return OPT_VALUE("case-bb-stability-prescription", m_CaseBBStabilityPrescription); }
    
    CHE_OPTION                                  CHE_Option() const                                                      { return OPT_VALUE("chemically-homogeneous-evolution", m_CheOption); }

    bool                                        CirculariseBinaryDuringMassTransfer() const                             { return OPT_VALUE("circularise-binary-during-mass-transfer", m_CirculariseBinaryDuringMassTransfer); }

    double                                      CommonEnvelopeAlpha() const                                             { return OPT_VALUE("common-envelope-alpha", m_CommonEnvelopeAlpha); }
    double                                      CommonEnvelopeAlphaThermal() const                                      { return OPT_VALUE("common-envelope-alpha-thermal", m_CommonEnvelopeAlphaThermal); }
    double                                      CommonEnvelopeLambda() const                                            { return OPT_VALUE("common-envelope-lambda", m_CommonEnvelopeLambda); }
    double                                      CommonEnvelopeLambdaMultiplier() const                                  { return OPT_VALUE("common-envelope-lambda-multiplier", m_CommonEnvelopeLambdaMultiplier); }
    CE_LAMBDA_PRESCRIPTION                      CommonEnvelopeLambdaPrescription() const                                { return OPT_VALUE("common-envelope-lambda-prescription", m_CommonEnvelopeLambdaPrescription); }
    double                                      CommonEnvelopeMassAccretionConstant() const                             { return OPT_VALUE("common-envelope-mass-accretion-constant", m_CommonEnvelopeMassAccretionConstant); }
    double                                      CommonEnvelopeMassAccretionMax() const                                  { return OPT_VALUE("common-envelope-mass-accretion-max", m_CommonEnvelopeMassAccretionMax); }
    double                                      CommonEnvelopeMassAccretionMin() const                                  { return OPT_VALUE("common-envelope-mass-accretion-min", m_CommonEnvelopeMassAccretionMin); }
    CE_ACCRETION_PRESCRIPTION                   CommonEnvelopeMassAccretionPrescription() const                         { return OPT_VALUE("common-envelope-mass-accretion-prescription", m_CommonEnvelopeMassAccretionPrescription); }
    double                                      CommonEnvelopeRecombinationEnergyDensity() const                        { return OPT_VALUE("common-envelope-recombination-energy-density", m_CommonEnvelopeRecombinationEnergyDensity); }
    double                                      CommonEnvelopeSlopeKruckow() const                                      { return OPT_VALUE("common-envelope-slope-kruckow", m_CommonEnvelopeSlopeKruckow); }

    vector<string>                              DebugClasses() const                                                    { return m_CmdLine.optionValues.m_DebugClasses; }
    int                                         DebugLevel() const                                                      { return m_CmdLine.optionValues.m_DebugLevel; }
    bool                                        DebugToFile() const                                                     { return m_CmdLine.optionValues.m_DebugToFile; }
    bool                                        DetailedOutput() const                                                  { return m_CmdLine.optionValues.m_DetailedOutput; }

    bool                                        EnableWarnings() const                                                  { return m_CmdLine.optionValues.m_EnableWarnings; }
    bool                                        ErrorsToFile() const                                                    { return m_CmdLine.optionValues.m_ErrorsToFile; }
    double                                      Eccentricity() const                                                    { return OPT_VALUE("eccentricity", m_Eccentricity); }
    ECCENTRICITY_DISTRIBUTION                   EccentricityDistribution() const                                        { return OPT_VALUE("eccentricity-distribution", m_EccentricityDistribution); }
    double                                      EccentricityDistributionMax() const                                     { return OPT_VALUE("eccentricity-distribution-max", m_EccentricityDistributionMax); }
    double                                      EccentricityDistributionMin() const                                     { return OPT_VALUE("eccentricity-distribution-min", m_EccentricityDistributionMin); }
    double                                      EddingtonAccretionFactor() const                                        { return OPT_VALUE("eddington-accretion-factor", m_EddingtonAccretionFactor); }
    ENVELOPE_STATE_PRESCRIPTION                 EnvelopeStatePrescription() const                                       { return OPT_VALUE("envelope-state-prescription", m_EnvelopeStatePrescription); }
    bool                                        EvolvePulsars() const                                                   { return m_CmdLine.optionValues.m_EvolvePulsars; }
    bool                                        EvolveUnboundSystems() const                                            { return m_CmdLine.optionValues.m_EvolveUnboundSystems; }

    bool                                        FixedRandomSeedCmdLine() const                                          { return m_CmdLine.optionValues.m_FixedRandomSeed; }
    bool                                        FixedRandomSeedGridLine() const                                         { return m_GridLine.optionValues.m_FixedRandomSeed; }
    double                                      FixedUK() const                                                         { return m_GridLine.optionValues.m_UseFixedUK || m_CmdLine.optionValues.m_FixedUK; }
    SN_ENGINE                                   FryerSupernovaEngine() const                                            { return OPT_VALUE("fryer-supernova-engine", m_FryerSupernovaEngine); }

    string                                      GridFilename() const                                                    { return m_CmdLine.optionValues.m_GridFilename; }

    double                                      InitialMass() const                                                     { return OPT_VALUE("initial-mass", m_InitialMass); }
    double                                      InitialMass1() const                                                    { return OPT_VALUE("initial-mass-1", m_InitialMass1); }
    double                                      InitialMass2() const                                                    { return OPT_VALUE("initial-mass-2", m_InitialMass2); }

    INITIAL_MASS_FUNCTION                       InitialMassFunction() const                                             { return OPT_VALUE("initial-mass-function", m_InitialMassFunction); }
    double                                      InitialMassFunctionMax() const                                          { return OPT_VALUE("initial-mass-max", m_InitialMassFunctionMax); }
    double                                      InitialMassFunctionMin() const                                          { return OPT_VALUE("initial-mass-min", m_InitialMassFunctionMin); }
    double                                      InitialMassFunctionPower() const                                        { return OPT_VALUE("initial-mass-power", m_InitialMassFunctionPower); }

    KICK_DIRECTION_DISTRIBUTION                 KickDirectionDistribution() const                                       { return OPT_VALUE("kick-direction", m_KickDirectionDistribution); }
    double                                      KickDirectionPower() const                                              { return OPT_VALUE("kick-direction-power", m_KickDirectionPower); }
    double                                      KickScalingFactor() const                                               { return OPT_VALUE("kick-scaling-factor", m_KickScalingFactor); }
    KICK_MAGNITUDE_DISTRIBUTION                 KickMagnitudeDistribution() const                                       { return OPT_VALUE("kick-magnitude-distribution", m_KickMagnitudeDistribution); }

    double                                      KickMagnitudeDistributionMaximum() const                                { return OPT_VALUE("kick-magnitude-max", m_KickMagnitudeDistributionMaximum); }

    double                                      KickMagnitudeDistributionSigmaCCSN_BH() const                           { return OPT_VALUE("kick-magnitude-sigma-ccsn-bh", m_KickMagnitudeDistributionSigmaCCSN_BH); }
    double                                      KickMagnitudeDistributionSigmaCCSN_NS() const                           { return OPT_VALUE("kick-magnitude-sigma-ccsn-ns", m_KickMagnitudeDistributionSigmaCCSN_NS); }
    double                                      KickMagnitudeDistributionSigmaForECSN() const                           { return OPT_VALUE("kick-magnitude-sigma-ecsn", m_KickMagnitudeDistributionSigmaForECSN); }
    double                                      KickMagnitudeDistributionSigmaForUSSN() const                           { return OPT_VALUE("kick-magnitude-sigma-ussn", m_KickMagnitudeDistributionSigmaForUSSN); }

    double                                      KickMagnitude() const                                                   { return OPT_VALUE("kick-magnitude", m_KickMagnitude); }
    double                                      KickMagnitude1() const                                                  { return OPT_VALUE("kick-magnitude-1", m_KickMagnitude1); }
    double                                      KickMagnitude2() const                                                  { return OPT_VALUE("kick-magnitude-2", m_KickMagnitude2); }

    double                                      KickMagnitudeRandom() const                                             { return OPT_VALUE("kick-magnitude-random", m_KickMagnitudeRandom); }
    double                                      KickMagnitudeRandom1() const                                            { return OPT_VALUE("kick-magnitude-random-1", m_KickMagnitudeRandom1); }
    double                                      KickMagnitudeRandom2() const                                            { return OPT_VALUE("kick-magnitude-random-2", m_KickMagnitudeRandom2); }

    double                                      SN_MeanAnomaly1() const                                                 { return OPT_VALUE("kick-mean-anomaly-1", m_KickMeanAnomaly1); }
    double                                      SN_MeanAnomaly2() const                                                 { return OPT_VALUE("kick-mean-anomaly-2", m_KickMeanAnomaly2); }
    double                                      SN_Phi1() const                                                         { return OPT_VALUE("kick-phi-1", m_KickPhi1); }
    double                                      SN_Phi2() const                                                         { return OPT_VALUE("kick-phi-2", m_KickPhi2); }
    double                                      SN_Theta1() const                                                       { return OPT_VALUE("kick-theta-1", m_KickTheta1); }
    double                                      SN_Theta2() const                                                       { return OPT_VALUE("kick-theta-2", m_KickTheta2); }

    vector<string>                              LogClasses() const                                                      { return m_CmdLine.optionValues.m_LogClasses; }
    string                                      LogfileBeBinaries() const                                               { return m_CmdLine.optionValues.m_LogfileBeBinaries; }
    string                                      LogfileCommonEnvelopes() const                                          { return m_CmdLine.optionValues.m_LogfileCommonEnvelopes; }
    string                                      LogfileDefinitionsFilename() const                                      { return m_CmdLine.optionValues.m_LogfileDefinitionsFilename; }
    DELIMITER                                   LogfileDelimiter() const                                                { return m_CmdLine.optionValues.m_LogfileDelimiter; }
    string                                      LogfileDelimiterString() const                                          { return m_CmdLine.optionValues.m_LogfileDelimiterString; }
    string                                      LogfileDetailedOutput() const                                           { return m_CmdLine.optionValues.m_LogfileDetailedOutput; }
    string                                      LogfileDoubleCompactObjects() const                                     { return m_CmdLine.optionValues.m_LogfileDoubleCompactObjects; }
    string                                      LogfileNamePrefix() const                                               { return m_CmdLine.optionValues.m_LogfileNamePrefix; }
    string                                      LogfilePulsarEvolution() const                                          { return m_CmdLine.optionValues.m_LogfilePulsarEvolution; }
    string                                      LogfileRLOFParameters() const                                           { return m_CmdLine.optionValues.m_LogfileRLOFParameters; }
    string                                      LogfileSupernovae() const                                               { return m_CmdLine.optionValues.m_LogfileSupernovae; }
    string                                      LogfileSwitchLog() const                                                { return m_CmdLine.optionValues.m_LogfileSwitchLog; }
    string                                      LogfileSystemParameters() const                                         { return m_CmdLine.optionValues.m_LogfileSystemParameters; }
    int                                         LogLevel() const                                                        { return m_CmdLine.optionValues.m_LogLevel; }

    double                                      LuminousBlueVariableFactor() const                                      { return OPT_VALUE("luminous-blue-variable-multiplier", m_LuminousBlueVariableFactor); }

    MASS_LOSS_PRESCRIPTION                      MassLossPrescription() const                                            { return OPT_VALUE("mass-loss-prescription", m_MassLossPrescription); }

    MASS_RATIO_DISTRIBUTION                     MassRatioDistribution() const                                           { return OPT_VALUE("mass-ratio-distribution", m_MassRatioDistribution); }
    double                                      MassRatioDistributionMax() const                                        { return OPT_VALUE("mass-ratio-max", m_MassRatioDistributionMax); }
    double                                      MassRatioDistributionMin() const                                        { return OPT_VALUE("mass-ratio-min", m_MassRatioDistributionMin); }

    MT_ACCRETION_EFFICIENCY_PRESCRIPTION        MassTransferAccretionEfficiencyPrescription() const                     { return OPT_VALUE("mass-transfer-accretion-efficiency-prescription", m_MassTransferAccretionEfficiencyPrescription); }
    MT_ANGULAR_MOMENTUM_LOSS_PRESCRIPTION       MassTransferAngularMomentumLossPrescription() const                     { return OPT_VALUE("mass-transfer-angular-momentum-loss-prescription", m_MassTransferAngularMomentumLossPrescription); }
    double                                      MassTransferCParameter() const                                          { return OPT_VALUE("mass-transfer-thermal-limit-c", m_MassTransferCParameter); }

    // AVG
    bool                                        MassTransferCriticalMassRatioMSLowMass() const                          { return m_CmdLine.optionValues.m_MassTransferCriticalMassRatioMSLowMass; }     // JR: no option implemented - always FALSE
    double                                      MassTransferCriticalMassRatioMSLowMassDegenerateAccretor() const        { return OPT_VALUE("critical-mass-ratio-ms-low-mass-degenerate-accretor", m_MassTransferCriticalMassRatioMSLowMassDegenerateAccretor); }
    double                                      MassTransferCriticalMassRatioMSLowMassNonDegenerateAccretor() const     { return OPT_VALUE("critical-mass-ratio-ms-low-mass-non-degenerate-accretor", m_MassTransferCriticalMassRatioMSLowMassNonDegenerateAccretor); }
    bool                                        MassTransferCriticalMassRatioMSHighMass() const                         { return m_CmdLine.optionValues.m_MassTransferCriticalMassRatioMSHighMass; }    // JR: no option implemented - always FALSE
    double                                      MassTransferCriticalMassRatioMSHighMassDegenerateAccretor() const       { return OPT_VALUE("critical-mass-ratio-ms-high-mass-degenerate-accretor", m_MassTransferCriticalMassRatioMSHighMassDegenerateAccretor); }
    double                                      MassTransferCriticalMassRatioMSHighMassNonDegenerateAccretor() const    { return OPT_VALUE("critical-mass-ratio-ms-high-mass-non-degenerate-accretor", m_MassTransferCriticalMassRatioMSHighMassNonDegenerateAccretor); }
    bool                                        MassTransferCriticalMassRatioGiant() const                              { return m_CmdLine.optionValues.m_MassTransferCriticalMassRatioGiant; }         // JR: no option implemented - always FALSE
    double                                      MassTransferCriticalMassRatioGiantDegenerateAccretor() const            { return OPT_VALUE("critical-mass-ratio-giant-degenerate-accretor", m_MassTransferCriticalMassRatioGiantDegenerateAccretor); }
    double                                      MassTransferCriticalMassRatioGiantNonDegenerateAccretor() const         { return OPT_VALUE("critical-mass-ratio-giant-non-degenerate-accretor", m_MassTransferCriticalMassRatioGiantNonDegenerateAccretor); }
    bool                                        MassTransferCriticalMassRatioHG() const                                 { return m_CmdLine.optionValues.m_MassTransferCriticalMassRatioHG; }            // JR: no option implemented - always FALSE
    double                                      MassTransferCriticalMassRatioHGDegenerateAccretor() const               { return OPT_VALUE("critical-mass-ratio-hg-degenerate-accretor", m_MassTransferCriticalMassRatioHGDegenerateAccretor); }
    double                                      MassTransferCriticalMassRatioHGNonDegenerateAccretor() const            { return OPT_VALUE("critical-mass-ratio-hg-non-degenerate-accretor", m_MassTransferCriticalMassRatioHGNonDegenerateAccretor); }
    bool                                        MassTransferCriticalMassRatioHeliumGiant() const                        { return m_CmdLine.optionValues.m_MassTransferCriticalMassRatioHeliumGiant; }   // JR: no option implemented - always FALSE
    double                                      MassTransferCriticalMassRatioHeliumGiantDegenerateAccretor() const      { return OPT_VALUE("critical-mass-ratio-helium-giant-degenerate-accretor", m_MassTransferCriticalMassRatioHeliumGiantDegenerateAccretor); }
    double                                      MassTransferCriticalMassRatioHeliumGiantNonDegenerateAccretor() const   { return OPT_VALUE("critical-mass-ratio-helium-giant-non-degenerate-accretor", m_MassTransferCriticalMassRatioHeliumGiantNonDegenerateAccretor); }
    bool                                        MassTransferCriticalMassRatioHeliumHG() const                           { return m_CmdLine.optionValues.m_MassTransferCriticalMassRatioHeliumHG; }      // JR: no option implemented - always FALSE
    double                                      MassTransferCriticalMassRatioHeliumHGDegenerateAccretor() const         { return OPT_VALUE("critical-mass-ratio-helium-hg-degenerate-accretor", m_MassTransferCriticalMassRatioHeliumHGDegenerateAccretor); }
    double                                      MassTransferCriticalMassRatioHeliumHGNonDegenerateAccretor() const      { return OPT_VALUE("critical-mass-ratio-helium-hg-non-degenerate-accretor", m_MassTransferCriticalMassRatioHeliumHGNonDegenerateAccretor); }
    bool                                        MassTransferCriticalMassRatioHeliumMS() const                           { return m_CmdLine.optionValues.m_MassTransferCriticalMassRatioHeliumMS; }      // JR: no option implemented - always FALSE
    double                                      MassTransferCriticalMassRatioHeliumMSDegenerateAccretor() const         { return OPT_VALUE("critical-mass-ratio-helium-ms-degenerate-accretor", m_MassTransferCriticalMassRatioHeliumMSDegenerateAccretor); }
    double                                      MassTransferCriticalMassRatioHeliumMSNonDegenerateAccretor() const      { return OPT_VALUE("critical-mass-ratio-helium-ms-non-degenerate-accretor", m_MassTransferCriticalMassRatioHeliumMSNonDegenerateAccretor); }
    bool                                        MassTransferCriticalMassRatioWhiteDwarf() const                         { return m_CmdLine.optionValues.m_MassTransferCriticalMassRatioWhiteDwarf; }    // JR: no option implemented - always FALSE
    double                                      MassTransferCriticalMassRatioWhiteDwarfDegenerateAccretor() const       { return OPT_VALUE("critical-mass-ratio-white-dwarf-degenerate-accretor", m_MassTransferCriticalMassRatioWhiteDwarfDegenerateAccretor); }
    double                                      MassTransferCriticalMassRatioWhiteDwarfNonDegenerateAccretor() const    { return OPT_VALUE("critical-mass-ratio-white-dwarf-non-degenerate-accretor", m_MassTransferCriticalMassRatioWhiteDwarfNonDegenerateAccretor); }

    double                                      MassTransferFractionAccreted() const                                    { return OPT_VALUE("mass-transfer-fa", m_MassTransferFractionAccreted); }
    double                                      MassTransferJloss() const                                               { return OPT_VALUE("mass-transfer-jloss", m_MassTransferJloss); }
    MT_REJUVENATION_PRESCRIPTION                MassTransferRejuvenationPrescription() const                            { return OPT_VALUE("mass-transfer-rejuvenation-prescription", m_MassTransferRejuvenationPrescription); }
    MT_THERMALLY_LIMITED_VARIATION              MassTransferThermallyLimitedVariation() const                           { return OPT_VALUE("mass-transfer-thermal-limit-accretor", m_MassTransferThermallyLimitedVariation); }
    double                                      MaxEvolutionTime() const                                                { return m_CmdLine.optionValues.m_MaxEvolutionTime; }
    int                                         MaxNumberOfTimestepIterations() const                                   { return m_CmdLine.optionValues.m_MaxNumberOfTimestepIterations; }

    double                                      MCBUR1() const                                                          { return OPT_VALUE("mcbur1", m_mCBUR1); }

    double                                      Metallicity() const                                                     { return OPT_VALUE("metallicity", m_Metallicity); }

    double                                      MinimumMassSecondary() const                                            { return OPT_VALUE("minimum-secondary-mass", m_MinimumMassSecondary); }
    double                                      MaximumNeutronStarMass() const                                          { return OPT_VALUE("maximum-neutron-star-mass", m_MaximumNeutronStarMass); }

    int                                         nBatchesUsed() const                                                    { return m_CmdLine.optionValues.m_nBatchesUsed; }

    NEUTRINO_MASS_LOSS_PRESCRIPTION             NeutrinoMassLossAssumptionBH() const                                    { return OPT_VALUE("neutrino-mass-loss-bh-formation", m_NeutrinoMassLossAssumptionBH); }
    double                                      NeutrinoMassLossValueBH() const                                         { return OPT_VALUE("neutrino-mass-loss-bh-formation-value", m_NeutrinoMassLossValueBH); }

    NS_EOS                                      NeutronStarEquationOfState() const                                      { return OPT_VALUE("neutron-star-equation-of-state", m_NeutronStarEquationOfState); }

    int                                         nObjectsToEvolve() const                                                { return m_CmdLine.optionValues.m_ObjectsToEvolve; }
    bool                                        OptimisticCHE() const                                                   { CHE_OPTION che = OPT_VALUE("chemically-homogeneous-evolution", m_CheOption); return che == CHE_OPTION::OPTIMISTIC; }

    string                                      CmdLineOptionsDetails() const                                           { return m_CmdLineOptionsDetails; }

    string                                      OutputContainerName() const                                             { return m_CmdLine.optionValues.m_OutputContainerName; }
    string                                      OutputPathString() const                                                { return m_CmdLine.optionValues.m_OutputPath.string(); }

    double                                      PairInstabilityLowerLimit() const                                       { return OPT_VALUE("pisn-lower-limit", m_PairInstabilityLowerLimit); }
    double                                      PairInstabilityUpperLimit() const                                       { return OPT_VALUE("pisn-upper-limit", m_PairInstabilityUpperLimit); }

    double                                      PeriodDistributionMax() const                                           { return OPT_VALUE("orbital-period-max", m_PeriodDistributionMax); }
    double                                      PeriodDistributionMin() const                                           { return OPT_VALUE("orbital-period-min", m_PeriodDistributionMin); }

    bool                                        PopulationDataPrinting() const                                          { return m_CmdLine.optionValues.m_PopulationDataPrinting; }
    bool                                        PrintBoolAsString() const                                               { return m_CmdLine.optionValues.m_PrintBoolAsString; }

    PULSAR_BIRTH_MAGNETIC_FIELD_DISTRIBUTION    PulsarBirthMagneticFieldDistribution() const                            { return OPT_VALUE("pulsar-birth-magnetic-field-distribution", m_PulsarBirthMagneticFieldDistribution); }
    double                                      PulsarBirthMagneticFieldDistributionMax() const                         { return OPT_VALUE("pulsar-birth-magnetic-field-distribution-max", m_PulsarBirthMagneticFieldDistributionMax); }
    double                                      PulsarBirthMagneticFieldDistributionMin() const                         { return OPT_VALUE("pulsar-birth-magnetic-field-distribution-min", m_PulsarBirthMagneticFieldDistributionMin); }

    PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION       PulsarBirthSpinPeriodDistribution() const                               { return OPT_VALUE("pulsar-birth-spin-period-distribution", m_PulsarBirthSpinPeriodDistribution); }
    double                                      PulsarBirthSpinPeriodDistributionMax() const                            { return OPT_VALUE("pulsar-birth-spin-period-distribution-max", m_PulsarBirthSpinPeriodDistributionMax); }
    double                                      PulsarBirthSpinPeriodDistributionMin() const                            { return OPT_VALUE("pulsar-birth-spin-period-distribution-min", m_PulsarBirthSpinPeriodDistributionMin); }

    double                                      PulsarLog10MinimumMagneticField() const                                 { return OPT_VALUE("pulsar-minimum-magnetic-field", m_PulsarLog10MinimumMagneticField); }

    double                                      PulsarMagneticFieldDecayMassscale() const                               { return OPT_VALUE("pulsar-magnetic-field-decay-massscale", m_PulsarMagneticFieldDecayMassscale); }
    double                                      PulsarMagneticFieldDecayTimescale() const                               { return OPT_VALUE("pulsar-magnetic-field-decay-timescale", m_PulsarMagneticFieldDecayTimescale); }

    PPI_PRESCRIPTION                            PulsationalPairInstabilityPrescription() const                          { return OPT_VALUE("pulsational-pair-instability-prescription", m_PulsationalPairInstabilityPrescription); }
    double                                      PulsationalPairInstabilityLowerLimit() const                            { return OPT_VALUE("ppi-lower-limit", m_PulsationalPairInstabilityLowerLimit); }
    double                                      PulsationalPairInstabilityUpperLimit() const                            { return OPT_VALUE("ppi-upper-limit", m_PulsationalPairInstabilityUpperLimit); }

    bool                                        Quiet() const                                                           { return m_CmdLine.optionValues.m_Quiet; }

    unsigned long int                           RandomSeed() const                                                      { return OPT_VALUE("random-seed", m_RandomSeed); }
    unsigned long int                           RandomSeedCmdLine() const                                               { return m_CmdLine.optionValues.m_RandomSeed; }
    unsigned long int                           RandomSeedGridLine() const                                              { return m_GridLine.optionValues.m_RandomSeed; }

    REMNANT_MASS_PRESCRIPTION                   RemnantMassPrescription() const                                         { return OPT_VALUE("remnant-mass-prescription", m_RemnantMassPrescription); }
    bool                                        RLOFPrinting() const                                                    { return m_CmdLine.optionValues.m_RlofPrinting; }

    ROTATIONAL_VELOCITY_DISTRIBUTION            RotationalVelocityDistribution() const                                  { return OPT_VALUE("rotational-velocity-distribution", m_RotationalVelocityDistribution); }
   
    double                                      SemiMajorAxis() const                                                   { return OPT_VALUE("semi-major-axis", m_SemiMajorAxis); }
    SEMI_MAJOR_AXIS_DISTRIBUTION                SemiMajorAxisDistribution() const                                       { return OPT_VALUE("semi-major-axis-distribution", m_SemiMajorAxisDistribution); }
    double                                      SemiMajorAxisDistributionMax() const                                    { return OPT_VALUE("semi-major-axis-max", m_SemiMajorAxisDistributionMax); }
    double                                      SemiMajorAxisDistributionMin() const                                    { return OPT_VALUE("semi-major-axis-min", m_SemiMajorAxisDistributionMin); }
    double                                      SemiMajorAxisDistributionPower() const                                  { return m_CmdLine.optionValues.m_SemiMajorAxisDistributionPower; }     // JR: no option implemented - always -1.0

    bool                                        RequestedHelp() const                                                   { return m_CmdLine.optionValues.m_VM["help"].as<bool>(); }
    bool                                        RequestedVersion() const                                                { return m_CmdLine.optionValues.m_VM["version"].as<bool>(); }

    bool                                        SwitchLog() const                                                       { return OPT_VALUE("switchlog", m_SwitchLog); }

    ZETA_PRESCRIPTION                           StellarZetaPrescription() const                                         { return OPT_VALUE("stellar-zeta-prescription", m_StellarZetaPrescription); }

    double                                      TimestepMultiplier() const                                              { return m_CmdLine.optionValues.m_TimestepMultiplier; }

    bool                                        UseFixedUK() const                                                      { return (m_GridLine.optionValues.m_UseFixedUK || m_CmdLine.optionValues.m_UseFixedUK); }
    bool                                        UseMassLoss() const                                                     { return OPT_VALUE("use-mass-loss", m_UseMassLoss); }
    bool                                        UseMassTransfer() const                                                 { return OPT_VALUE("mass-transfer", m_UseMassTransfer); }
    bool                                        UsePairInstabilitySupernovae() const                                    { return OPT_VALUE("pair-instability-supernovae", m_UsePairInstabilitySupernovae); }
    bool                                        UsePulsationalPairInstability() const                                   { return OPT_VALUE("pulsational-pair-instability", m_UsePulsationalPairInstability); }

    double                                      WolfRayetFactor() const                                                 { return OPT_VALUE("wolf-rayet-multiplier", m_WolfRayetFactor); }

    double                                      ZetaRadiativeEnvelopeGiant() const                                      { return OPT_VALUE("zeta-radiative-envelope-giant", m_ZetaRadiativeEnvelopeGiant); }
    double                                      ZetaMainSequence() const                                                { return OPT_VALUE("zeta-main-sequence", m_ZetaMainSequence); }
    double                                      ZetaAdiabaticArbitrary() const                                          { return OPT_VALUE("zeta-adiabatic-arbitrary", m_ZetaAdiabaticArbitrary); }

};

#endif // __Options_H__