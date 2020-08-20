#include "Options.h"

namespace po = boost::program_options;
namespace fs = boost::filesystem;

Options* Options::m_Instance = nullptr;

// this is required to set default value for boost program options of type vector<string>
namespace std
{
  std::ostream& operator<<(std::ostream &os, const std::vector<string> &vec) {    
    for (auto item : vec) os << item << " ";
    return os; 
  }
} 


Options* Options::Instance() {
    if (!m_Instance) {
        m_Instance = new Options();
    }
    return m_Instance;
}


string Options::ProgramOptionDetails(const boost::program_options::variables_map p_VM) {

    std::ostringstream ss;

    ss << "COMMAND LINE OPTIONS\n-------------------\n\n";

    for (po::variables_map::const_iterator it = p_VM.begin(); it != p_VM.end(); it++) {

        // option name
        ss << it->first << " = ";
        
        if (((boost::any)it->second.value()).empty()) ss << "<EMPTY_OPTION>\n";
        else {

            // determine if option values was supplied, or whether the default was used

            string valueSource;
            if (p_VM[it->first].defaulted() || it->second.defaulted()) valueSource = "DEFAULT_USED";
            else                                                       valueSource = "USER_SUPPLIED";

            // find data type & print value
            // handles most data types - add others if they cause problems

            bool isCharPtr = false;
            bool isStr     = false;

            // (pre)check for data type = charPtr
            try {
                boost::any_cast<const char *>(it->second.value());
                isCharPtr = true;
            } catch (const boost::bad_any_cast &) {
                isCharPtr = false;
            }

            if (!isCharPtr) {
                // (pre)check for data type = string
                try {
                    boost::any_cast<std::string>(it->second.value());
                    isStr = true;
                } catch (const boost::bad_any_cast &) {
                    isStr = false;
                }
            }

            // find other data types
            // it's not pretty, but it works

            if (isCharPtr) ss << p_VM[it->first].as<const char *>() << ", " << valueSource << ", CONST_CHAR_*";

            else if (isStr) {
                std::string tmp = p_VM[it->first].as<std::string>();
                if (tmp.size()) ss << "'" << tmp << "'";
                else            ss << "''";
                ss << ", " << valueSource << ", STRING";
            }

            else if (((boost::any)it->second.value()).type() == typeid(signed                )) ss << p_VM[it->first].as<signed                >() << ", " << valueSource << ", SIGNED";
            else if (((boost::any)it->second.value()).type() == typeid(unsigned              )) ss << p_VM[it->first].as<unsigned              >() << ", " << valueSource << ", UNSIGNED";

            else if (((boost::any)it->second.value()).type() == typeid(short                 )) ss << p_VM[it->first].as<short                 >() << ", " << valueSource << ", SHORT";
            else if (((boost::any)it->second.value()).type() == typeid(signed short          )) ss << p_VM[it->first].as<signed short          >() << ", " << valueSource << ", SIGNED_SHORT";
            else if (((boost::any)it->second.value()).type() == typeid(unsigned short        )) ss << p_VM[it->first].as<unsigned short        >() << ", " << valueSource << ", UNSIGNED_SHORT";

            else if (((boost::any)it->second.value()).type() == typeid(short int             )) ss << p_VM[it->first].as<short int             >() << ", " << valueSource << ", SHORT_INT";
            else if (((boost::any)it->second.value()).type() == typeid(signed short int      )) ss << p_VM[it->first].as<signed short int      >() << ", " << valueSource << ", SIGNED_SHORT_INT";
            else if (((boost::any)it->second.value()).type() == typeid(unsigned short int    )) ss << p_VM[it->first].as<unsigned short int    >() << ", " << valueSource << ", UNSIGNED_SHORT_INT";

            else if (((boost::any)it->second.value()).type() == typeid(int                   )) ss << p_VM[it->first].as<int                   >() << ", " << valueSource << ", INT";
            else if (((boost::any)it->second.value()).type() == typeid(signed int            )) ss << p_VM[it->first].as<signed int            >() << ", " << valueSource << ", SIGNED_INT";
            else if (((boost::any)it->second.value()).type() == typeid(unsigned int          )) ss << p_VM[it->first].as<unsigned int          >() << ", " << valueSource << ", UNSIGNED_INT";

            else if (((boost::any)it->second.value()).type() == typeid(long                  )) ss << p_VM[it->first].as<long                  >() << ", " << valueSource << ", LONG";
            else if (((boost::any)it->second.value()).type() == typeid(signed long           )) ss << p_VM[it->first].as<signed long           >() << ", " << valueSource << ", SIGNED_LONG";
            else if (((boost::any)it->second.value()).type() == typeid(unsigned long         )) ss << p_VM[it->first].as<unsigned long         >() << ", " << valueSource << ", UNSIGNED_LONG";

            else if (((boost::any)it->second.value()).type() == typeid(long int              )) ss << p_VM[it->first].as<long int              >() << ", " << valueSource << ", LONG_INT";
            else if (((boost::any)it->second.value()).type() == typeid(signed long int       )) ss << p_VM[it->first].as<signed long int       >() << ", " << valueSource << ", SIGNED_LONG_INT";
            else if (((boost::any)it->second.value()).type() == typeid(unsigned long int     )) ss << p_VM[it->first].as<unsigned long int     >() << ", " << valueSource << ", UNSIGNED_LONG_INT";

            else if (((boost::any)it->second.value()).type() == typeid(long long             )) ss << p_VM[it->first].as<long long             >() << ", " << valueSource << ", LONG_LONG";
            else if (((boost::any)it->second.value()).type() == typeid(signed long long      )) ss << p_VM[it->first].as<signed long long      >() << ", " << valueSource << ", SIGNED_LONG_LONG";
            else if (((boost::any)it->second.value()).type() == typeid(unsigned long long    )) ss << p_VM[it->first].as<unsigned long long    >() << ", " << valueSource << ", UNSIGNED_LONG_LONG";

            else if (((boost::any)it->second.value()).type() == typeid(long long int         )) ss << p_VM[it->first].as<long long int         >() << ", " << valueSource << ", LONG_LONG_INT";
            else if (((boost::any)it->second.value()).type() == typeid(signed long long int  )) ss << p_VM[it->first].as<signed long long int  >() << ", " << valueSource << ", SIGNED_LONG_LONG_INT";
            else if (((boost::any)it->second.value()).type() == typeid(unsigned long long int)) ss << p_VM[it->first].as<unsigned long long int>() << ", " << valueSource << ", UNSIGNED_LONG_LONG_INT";

            else if (((boost::any)it->second.value()).type() == typeid(float                 )) ss << p_VM[it->first].as<float                 >() << ", " << valueSource << ", FLOAT";
            else if (((boost::any)it->second.value()).type() == typeid(double                )) ss << p_VM[it->first].as<double                >() << ", " << valueSource << ", DOUBLE";
            else if (((boost::any)it->second.value()).type() == typeid(long double           )) ss << p_VM[it->first].as<long double           >() << ", " << valueSource << ", LONG_DOUBLE";

            else if (((boost::any)it->second.value()).type() == typeid(char                  )) ss << p_VM[it->first].as<char                  >() << ", " << valueSource << ", CHAR";
            else if (((boost::any)it->second.value()).type() == typeid(signed char           )) ss << p_VM[it->first].as<signed char           >() << ", " << valueSource << ", SIGNED_CHAR";
            else if (((boost::any)it->second.value()).type() == typeid(unsigned char         )) ss << p_VM[it->first].as<unsigned char         >() << ", " << valueSource << ", UNSIGNED_CHAR";

            else if (((boost::any)it->second.value()).type() == typeid(bool)) {
                bool v = p_VM[it->first].as<bool>();
                ss << (v ? "TRUE" : "FALSE") << ", " << valueSource << ", BOOL";
            } 

            else {  // Assume vector<string>
                try {
                    std::ostringstream elemsSS;
                    elemsSS << "{ ";
                    vector<string> tmp = p_VM[it->first].as<vector<string>>();
                    for (std::vector<string>::iterator elem=tmp.begin(); elem != tmp.end(); elem++) {
                        elemsSS << "'" << (*elem) << "', ";
                    }
                    string elems = elemsSS.str();
                    if (elems.size() > 2) elems.erase(elems.size() - 2);
                    else if (elems.size() == 2) elems.erase(elems.size() - 1);
                    elems += " }";
                    ss << elems << ", " << valueSource << ", VECTOR<STRING>";
                } catch (const boost::bad_any_cast &) {
                    ss << "<UNKNOWN_DATA_TYPE>, " << valueSource << ", <UNKNOWN_DATA_TYPE>";
                }
            }
        }

        ss << "\n";
    }
  
    ss << "\n\nOTHER PARAMETERS\n----------------\n\n";

    ss << "fixedMetallicity   = " << (fixedMetallicity ? "TRUE" : "FALSE") << ", CALCULATED, BOOL\n";                           // fixedMetallicity
    ss << "useFixedUK         = " << (useFixedUK ? "TRUE" : "FALSE") << ", CALCULATED, BOOL\n";                                 // useFixedUK
    ss << "outputPath         = " << outputPath.string() << ", CALCULATED, STRING\n";                                           // outputPath (fully qualified)
    ss << "fixedRandomSeed    = " << (fixedRandomSeed ? "TRUE" : "FALSE") << ", CALCULATED, BOOL\n";                            // fixedRandomSeed

    return ss.str();
}


COMMANDLINE_STATUS Options::Initialise(int argc, char *argv[]) {

    InitialiseMemberVariables();

    return CommandLineSorter(argc, argv);        // parse commandline program options
}


void Options::InitialiseMemberVariables(void) {

    // This sets all of the program options to their default values -- can be modified via the command line.

    // flags

    allowRLOFAtBirth                                                = false;                                                                            // default is to not allow binaries that have one or both stars in RLOF  at birth to evolve
    allowTouchingAtBirth                                            = false;                                                                            // default is to not allow binaries that are touching at birth to evolve

    debugToFile                                                     = false;                                                                            // default is do not log debug statements to a log file
    errorsToFile                                                    = false;                                                                            // default is do not log error messages to a log file

    singleStar                                                      = false;                                                                            // Flag to evolve a single star

	lambdaCalculationEveryTimeStep                                  = false;
	zetaCalculationEveryTimeStep                                    = false;

	beBinaries                                                      = false;
    evolvePulsars                                                   = false;                                                                            // Whether to evolve pulsars
	evolveUnboundSystems                                            = false;                                                                            // Allow unbound syetms to evolve

    detailedOutput                                                  = false;                                                                            // Detailed output
    populationDataPrinting                                          = false;                                                                            // Print certain data for small populations, but not for larger one
    printBoolAsString                                               = false;                                                                            // default is do not print bool as string
    quiet                                                           = false;                                                                            // Suppress some of the printing

    // AVG - 17/03/2020 - Floor will uncomment when tested.
    //    nBatchesUsed                                                    = -1;                                                                               // Number of batches used, for STROOPWAFEL (AIS)


    // Public population synthesis variables
    nBinaries                                                       = 10;

    fixedRandomSeed                                                 = false;                                                                            // Whether to use a fixed random seed given by options.randomSeed (set to true if --random-seed is passed on command line)
    randomSeed                                                      = 0;                                                                                // Random seed to use


    // Specify how long to evolve binaries for
    maxEvolutionTime                                                = 13700.0;                                                                          // Maximum evolution time in Myrs
    maxNumberOfTimestepIterations                                   = 99999;                                                                            // Maximum number of timesteps to evolve binary for before giving up


    // Initial mass options
    initialMassFunction                                             = INITIAL_MASS_FUNCTION::KROUPA;                                                    // Default is KROUPA
    initialMassFunctionString                                       = INITIAL_MASS_FUNCTION_LABEL.at(initialMassFunction);                              // Distribution name
    initialMassFunctionMin                                          = 8.0;                                                                              // Default minimum
    initialMassFunctionMax                                          = 100.0;                                                                            // Default maximum
    initialMassFunctionPower                                        = -2.3;                                                                             // Default power


    // Initial mass ratios
    massRatioDistribution                                           = MASS_RATIO_DISTRIBUTION::FLAT;                                                    // deafult is FLAT - most likely want Flat or SANA2012
    massRatioDistributionString                                     = MASS_RATIO_DISTRIBUTION_LABEL.at(massRatioDistribution);                          // Distribution name
    massRatioDistributionMin                                        = 0.0;                                                                              // Default minimum
    massRatioDistributionMax                                        = 1.0;                                                                              // Default maximum

    minimumMassSecondary                                            = 0.0;                                                                              // Minimum mass of secondary to draw (in Msol)


    // Initial orbit options
    semiMajorAxisDistribution                                       = SEMI_MAJOR_AXIS_DISTRIBUTION::FLATINLOG;                                          // Default is FLATINLOG - most likely want FlatInLog or SANA2012
    semiMajorAxisDistributionString                                 = SEMI_MAJOR_AXIS_DISTRIBUTION_LABEL.at(SEMI_MAJOR_AXIS_DISTRIBUTION::FLATINLOG);   // Distribution name
    semiMajorAxisDistributionMin                                    = 0.1;                                                                              // Default minimum
    semiMajorAxisDistributionMax                                    = 1000.0;                                                                           // Default maximum
    semiMajorAxisDistributionPower                                  = -1.0;                                                                             // Default power


    // Initial orbital period
    periodDistributionMin                                           = 1.1;                                                                              // Minimum initial period in days
    periodDistributionMax                                           = 1000.0;                                                                           // Maximum initial period in days

    // Eccentricity
    eccentricityDistribution                                        = ECCENTRICITY_DISTRIBUTION::ZERO;                                                  // Default is ZERO
    eccentricityDistributionString                                  = ECCENTRICITY_DISTRIBUTION_LABEL.at(eccentricityDistribution);                     // Distribution name
    eccentricityDistributionMin                                     = 0.0;                                                                              // Default minimum
    eccentricityDistributionMax                                     = 1.0;                                                                              // Default maximum

    // Kick options
    kickVelocityDistribution                                        = KICK_VELOCITY_DISTRIBUTION::MAXWELLIAN;		                                    // Which kick velocity distribution to use
    kickVelocityDistributionString                                  = KICK_VELOCITY_DISTRIBUTION_LABEL.at(kickVelocityDistribution);		            // Which kick velocity distribution to use
    kickVelocityDistributionSigmaCCSN_NS                            = 250;                                                                              // Kick velocity sigma in km s^-1 for neutron stars
    kickVelocityDistributionSigmaCCSN_BH                            = 250;                                                                              // Kick velocity sigma in km s^-1 for black holes
    kickVelocityDistributionMaximum                                 = -1.0;                                                                             // Maximum kick velocity to draw in km s^-1. Ignored if < 0
    kickVelocityDistributionSigmaForECSN                            = 30.0;                                                                             // Characteristic kick velocity for an ECSN in km s^-1
    kickVelocityDistributionSigmaForUSSN   	                        = 30.0;                                                                             // Characteristic kick velocity for an USSN in km s^-1
	kickScalingFactor						                        = 1.0;				                                                                // Arbitrary factor for scaling kicks


    // Black hole kicks
    blackHoleKicksOption                                            = BLACK_HOLE_KICK_OPTION::FALLBACK;
    blackHoleKicksString                                            = BLACK_HOLE_KICK_OPTION_LABEL.at(blackHoleKicksOption);


    // Supernova remnant mass prescription options
    remnantMassPrescription                                         = REMNANT_MASS_PRESCRIPTION::FRYER2012;
    remnantMassPrescriptionString                                   = REMNANT_MASS_PRESCRIPTION_LABEL.at(remnantMassPrescription);

    fryerSupernovaEngine                                            = SN_ENGINE::DELAYED;
    fryerSupernovaEngineString                                      = SN_ENGINE_LABEL.at(fryerSupernovaEngine);

    neutrinoMassLossAssumptionBH                                    = NEUTRINO_MASS_LOSS_PRESCRIPTION::FIXED_FRACTION;                                  // Assumption to make about neutrino mass loss for BH formation
    neutrinoMassLossAssumptionBHString                              = NEUTRINO_MASS_LOSS_PRESCRIPTION_LABEL.at(neutrinoMassLossAssumptionBH);
    neutrinoMassLossValueBH                                         = 0.1;                                                                              // Value (corresponding to assumption) for neutrino mass loss for BH formation


    // Fixed uk options
    useFixedUK                                                      = false;
    fixedUK                                                         = -1.0;


    // Chemically Homogeneous Evolution
    cheOption                                                       = CHE_OPTION::NONE;                                                                 // whether and how to apply Chemically Homogeneous Evolution
    cheString                                                       = CHE_OPTION_LABEL.at(cheOption);


    // Pair instability and pulsational pair instability mass loss
    usePairInstabilitySupernovae                                    = false;                                                                            // Whether to use pair instability supernovae (PISN)
    pairInstabilityUpperLimit                                       = 135.0;                                                                            // Maximum core mass leading to PISN (Value in Belczynski+ 2016 is 135 Msol)
    pairInstabilityLowerLimit                                       = 60.0;                                                                             // Minimum core mass leading to PISN (Value in Belczynski+ 2016 is 65 Msol)

    usePulsationalPairInstability                                   = false;                                                                            // Whether to use pulsational pair instability (PPI)
    pulsationalPairInstabilityLowerLimit                            = 35.0;                                                                             // Minimum core mass leading to PPI, Value in Belczynski+ 2016 is 45 Msol
    pulsationalPairInstabilityUpperLimit                            = 60.0;                                                                             // Maximum core mass leading to PPI, Value in Belczynski+ 2016 is 65 Msol

    pulsationalPairInstabilityPrescription                          = PPI_PRESCRIPTION::COMPAS;                                                         // Prescription for PPI to use
    pulsationalPairInstabilityPrescriptionString                    = PPI_PRESCRIPTION_LABEL.at(pulsationalPairInstabilityPrescription);                // String for which PPI prescription to use

	maximumNeutronStarMass                                          = 3.0;									                                            // Maximum mass of a neutron star allowed, value in StarTrack is 3.0
    
    MCBUR1                                                          = MCBUR1HURLEY;                                                                     // Minimum core mass at base of the AGB to avoid fully degenerate CO core formation (Hurley value, Fryer+ and Belczynski+ use 1.83)


    // Kick direction option
    kickDirectionDistribution                                       = KICK_DIRECTION_DISTRIBUTION::ISOTROPIC;                                           // Which assumption for SN kicks: Possibilities: isotropic, inplane, perpendicular, powerlaw, wedge, poles
    kickDirectionDistributionString                                 = KICK_DIRECTION_DISTRIBUTION_LABEL.at(kickDirectionDistribution);		            // Which assumption for SN kicks: Possibilities: isotropic, inplane, perpendicular, powerlaw, wedge, poles
    kickDirectionPower                                              = 0.0;                                                                              // Power law power for the "power" SN kick direction choice

    // Output path
    outputPathString                                                = ".";                                                                              // String to hold the output directory
    defaultOutputPath                                               = boost::filesystem::current_path();                                                // Default output location
    outputPath                                                      = defaultOutputPath;                                                                // Desired output location
    outputContainerName                                             = DEFAULT_OUTPUT_CONTAINER_NAME;                                                    // Output container - this is a container (directory) created at outputPath to hold all output files
    
    // Mass loss options
    useMassLoss                                                     = false;                                                                            // Whether to use mass loss

    massLossPrescription                                            = MASS_LOSS_PRESCRIPTION::VINK;
    massLossPrescriptionString                                      = MASS_LOSS_PRESCRIPTION_LABEL.at(massLossPrescription);


    // Wind mass loss multiplicitive constants
    luminousBlueVariableFactor                                      = 1.5;                                                                              // Luminous blue variable mass loss enhancement factor
    wolfRayetFactor                                                 = 1.0;                                                                              // WR winds factor


    // Mass transfer options
    useMassTransfer                                                 = true;                                                                             // Whether to use mass transfer
	circulariseBinaryDuringMassTransfer         	                = false;						                                                    // Whether to circularise binary when it starts
	forceCaseBBBCStabilityFlag                                      = false;									                                        // Whether if all case BB/BC systems are forced to be stable or unstable
	alwaysStableCaseBBBCFlag                                        = false;									                                        // Whether if case BB/BC is always stable
	angularMomentumConservationDuringCircularisation                = false;		                                                                    // Whether to conserve angular momentum while circularising or circularise to periastron


    // Options adaptive Roche Lobe Overflow prescription
    massTransferAdaptiveAlphaParameter                              = 0.5;
    maxPercentageAdaptiveMassTransfer                               = 0.01;


    // Options for mass transfer accretion efficiency
    massTransferFractionAccreted                                    = 1.0;
    massTransferCParameter                                          = 10.0;

    massTransferAccretionEfficiencyPrescription                     = MT_ACCRETION_EFFICIENCY_PRESCRIPTION::THERMALLY_LIMITED;
    massTransferAccretionEfficiencyPrescriptionString               = MT_ACCRETION_EFFICIENCY_PRESCRIPTION_LABEL.at(massTransferAccretionEfficiencyPrescription);


    // Mass transfer thermally limited options
	massTransferThermallyLimitedVariation                           = MT_THERMALLY_LIMITED_VARIATION::C_FACTOR;
	massTransferThermallyLimitedVariationString                     = MT_THERMALLY_LIMITED_VARIATION_LABEL.at(massTransferThermallyLimitedVariation);

    eddingtonAccretionFactor                                        = 1;                                                                                // Multiplication factor for eddington accretion for NS & BH
                                                                                                                                                        // (>1 is super-eddington, 0. is no accretion)

    massTransferJloss                                               = 1.0;

    // Mass transfer angular momentum loss prescription options
    massTransferAngularMomentumLossPrescription                     = MT_ANGULAR_MOMENTUM_LOSS_PRESCRIPTION::ISOTROPIC_RE_EMISSION;
    massTransferAngularMomentumLossPrescriptionString               = MT_ANGULAR_MOMENTUM_LOSS_PRESCRIPTION_LABEL.at(massTransferAngularMomentumLossPrescription);


    // Mass transfer rejuvenation prescriptions
    massTransferRejuvenationPrescription                            = MT_REJUVENATION_PRESCRIPTION::NONE;
    massTransferRejuvenationPrescriptionString                      = MT_REJUVENATION_PRESCRIPTION_LABEL.at(massTransferRejuvenationPrescription);


    // Mass transfer critical mass ratios
    massTransferCriticalMassRatioMSLowMass                          = false;                    			                                            // Whether to use critical mass ratios
    massTransferCriticalMassRatioMSLowMassNonDegenerateAccretor     = 1.44;                                                                             // Critical mass ratio for MT from a MS low mass star (Claeys+ 2014 = 1.44)
    massTransferCriticalMassRatioMSLowMassDegenerateAccretor        = 1.0;                                                                              // Critical mass ratio for MT from a MS low mass star on to a degenerate accretor (Claeys+ 2014 = 1.0)

    massTransferCriticalMassRatioMSHighMass                         = false;							                                                // Whether to use critical mass ratios
    massTransferCriticalMassRatioMSHighMassNonDegenerateAccretor    = 0.625;                                                                            // Critical mass ratio for MT from a MS high mass star (Claeys+ 2014 = 0.625)
    massTransferCriticalMassRatioMSHighMassDegenerateAccretor       = 0.0;                                                                              // Critical mass ratio for MT from a MS high mass star on to a degenerate accretor

    massTransferCriticalMassRatioHG                                 = false;                                                                            // Whether to use critical mass ratios
    massTransferCriticalMassRatioHGNonDegenerateAccretor            = 0.40;                                                                             // Critical mass ratio for MT from a HG star (Claeys+ 2014 = 0.25)
    massTransferCriticalMassRatioHGDegenerateAccretor               = 0.21;                                                                             // Critical mass ratio for MT from a HG star on to a degenerate accretor (Claeys+ 2014 = 0.21)

    massTransferCriticalMassRatioGiant                              = false;                                                                            // Whether to use critical mass ratios
    massTransferCriticalMassRatioGiantNonDegenerateAccretor         = 0.0;                                                                              // Critical mass ratio for MT from a giant
    massTransferCriticalMassRatioGiantDegenerateAccretor            = 0.87;                                                                             // Critical mass ratio for MT from a giant on to a degenerate accretor (Claeys+ 2014 = 0.81)

    massTransferCriticalMassRatioHeliumMS                           = false;	                                                                        // Whether to use critical mass ratios
    massTransferCriticalMassRatioHeliumMSNonDegenerateAccretor      = 0.625;                                                                            // Critical mass ratio for MT from a Helium MS star
    massTransferCriticalMassRatioHeliumMSDegenerateAccretor         = 0.0;                                                                              // Critical mass ratio for MT from a Helium MS star on to a degenerate accretor

    massTransferCriticalMassRatioHeliumHG                           = false;                                                                            // Whether to use critical mass ratios
    massTransferCriticalMassRatioHeliumHGNonDegenerateAccretor      = 0.25;		                                                                       	// Critical mass ratio for MT from a Helium HG star (Claeys+ 2014 = 0.25)
    massTransferCriticalMassRatioHeliumHGDegenerateAccretor         = 0.21;		                                                                        // Critical mass ratio for MT from a Helium HG star on to a degenerate accretor (Claeys+ 2014 = 0.21)

    massTransferCriticalMassRatioHeliumGiant                        = false;                                                                            // Whether to use critical mass ratios
    massTransferCriticalMassRatioHeliumGiantNonDegenerateAccretor   = 1.28;                                                                             // Critical mass ratio for MT from a Helium giant (Claeys+ 2014 = 0.25)
    massTransferCriticalMassRatioHeliumGiantDegenerateAccretor      = 0.87;                                                                             // Critical mass ratio for MT from a Helium giant on to a degenerate accretor

    massTransferCriticalMassRatioWhiteDwarf                         = false;                                                                            // Whether to use critical mass ratios
	massTransferCriticalMassRatioWhiteDwarfNonDegenerateAccretor    = 0.0;                                                                              // Critical mass ratio for MT from a White Dwarf (default = 0.0)
    massTransferCriticalMassRatioWhiteDwarfDegenerateAccretor       = 1.6;                                                                              // Critical mass ratio for MT from a White Dwarf on to a degenerate accretor (Claeys+ 2014 = 1.6)


    commonEnvelopeAlpha                                             = 1.0;                                                                              // Common envelope efficiency alpha parameter
    commonEnvelopeLambda                                            = 0.1;                                                                              // Common envelope Lambda parameter
	commonEnvelopeAlphaThermal                                      = 1.0;                                                                              // lambda = (alpha_th * lambda_b) + (1-alpha_th) * lambda_g
    commonEnvelopeLambdaMultiplier                                  = 1.0;                                                                              // Multiply common envelope lambda by some constant
    allowMainSequenceStarToSurviveCommonEnvelope                    = false;                                                                            // Whether or not to allow a main sequence star to survive a common envelope event

    // Accretion during common envelope
    commonEnvelopeMassAccretionPrescription                         = CE_ACCRETION_PRESCRIPTION::ZERO;
    commonEnvelopeMassAccretionPrescriptionString                   = CE_ACCRETION_PRESCRIPTION_LABEL.at(commonEnvelopeMassAccretionPrescription);
    
    commonEnvelopeMassAccretionMin                                  = 0.04;                                                                             // Minimum amount of mass accreted during CE in solar masses
    commonEnvelopeMassAccretionMax                                  = 0.1;                                                                              // Maximum amount of mass accreted during CE in solar masses
    commonEnvelopeMassAccretionConstant                             = 0.0;                                                                              // Constant value

    // Prescription for envelope state (radiative or convective)
    envelopeStatePrescription                                       = ENVELOPE_STATE_PRESCRIPTION::LEGACY;
    envelopeStatePrescriptionString                                 = ENVELOPE_STATE_PRESCRIPTION_LABEL.at(envelopeStatePrescription);

    
	// Common envelope lambda prescription
	commonEnvelopeLambdaPrescription                                = CE_LAMBDA_PRESCRIPTION::NANJING;                                                  // Which prescription to use for CE lambda
	commonEnvelopeLambdaPrescriptionString                          = CE_LAMBDA_PRESCRIPTION_LABEL.at(commonEnvelopeLambdaPrescription);                // String containing which prescription to use for CE lambda

	// Common envelope Nandez and Ivanova energy formalism
	revisedEnergyFormalismNandezIvanova	                            = false;						                                                    // Use the revised energy formalism from Nandez & Ivanova 2016
	maximumMassDonorNandezIvanova                                   = 2.0;								                                                // Maximum mass allowed to use the revised energy formalism in Msol
	commonEnvelopeRecombinationEnergyDensity                        = 1.5E13;					                                                        // Factor using to calculate the binding energy depending on the mass of the envelope

	// Common envelope power factor for Kruckow fit
	commonEnvelopeSlopeKruckow                                      = -4.0/5.0;								                                            // Common envelope power factor for Kruckow fit normalized according to Kruckow+2016, Fig. 1

	// Which prescription to use for calculating zetas
	stellarZetaPrescription                                  = ZETA_PRESCRIPTION::SOBERMAN;					                        // Prescription to use for calculating stellar zeta
	stellarZetaPrescriptionString                            = ZETA_PRESCRIPTION_LABEL.at(stellarZetaPrescription);				    	// String containing prescription to use for calculating stellar zetas

	zetaAdiabaticArbitrary                                          = 10000.0;                                                                          // large value, which will favour stable MT
    zetaMainSequence 	                                            = 2.0;
	zetaRadiativeEnvelopeGiant	                                    = 6.5;


    // Adaptive Importance Sampling options
    AISexploratoryPhase                                             = false;                                                                            // Flag for whether to run the AIS exploratory phase
    AISDCOtype                                                      = AIS_DCO::ALL;                                                                     // Which prescription to use for DCO type
    AISDCOtypeString                                                = AIS_DCO_LABEL.at(AIS_DCO::ALL);                                                   // String containing which type of DCOs to focus on
    AIShubble                                                       = false;                                                                            // Flag for excluding DCOs that do not merge in Hubble
    AISpessimistic                                                  = false;                                                                            // Flag for excluding DCOs that are Optmistic
    AISrefinementPhase                                              = false;                                                                            // Flag for whether to run the AIS refinement phase (step 2)
    AISrlof                                                         = false;                                                                            // Flag for excluding DCOs that RLOFSecondaryZAMS
    kappaGaussians                                                  = 2;                                                                                // scaling factor for the width of the Gaussian distributions in AIS main sampling phase

    // Metallicity options
    metallicity                                                     = ZSOL;
    fixedMetallicity                                                = true;


    // Neutron star equation of state
    neutronStarEquationOfState                                      = NS_EOS::SSE;
    neutronStarEquationOfStateString                                = NS_EOSLabel.at(NS_EOS::SSE);


    // Pulsar birth magnetic field distribution
    pulsarBirthMagneticFieldDistribution                            = PULSAR_BIRTH_MAGNETIC_FIELD_DISTRIBUTION::ZERO;
    pulsarBirthMagneticFieldDistributionString                      = PULSAR_BIRTH_MAGNETIC_FIELD_DISTRIBUTION_LABEL.at(pulsarBirthMagneticFieldDistribution);  // Which birth magnetic field distribution to use for pulsars
    pulsarBirthMagneticFieldDistributionMin                         = 11.0;                                                                             // Minimum pulsar birth magnetic field distribution
    pulsarBirthMagneticFieldDistributionMax                         = 13.0;                                                                             // Maximum pulsar birth magnetic field distribution


    // Pulsar birth spin period distribution string
    pulsarBirthSpinPeriodDistribution                               = PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION::ZERO;
    pulsarBirthSpinPeriodDistributionString                         = PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION_LABEL.at(pulsarBirthSpinPeriodDistribution);// Which birth spin period distribution to use for pulsars
    pulsarBirthSpinPeriodDistributionMin                            = 0.0;                                                                              // Minimum birth spin period (ms)
    pulsarBirthSpinPeriodDistributionMax                            = 100.0;                                                                            // Maximum birth spin period (ms)

    pulsarMagneticFieldDecayTimescale                               = 1000.0;                                                                           // Timescale on which magnetic field decays (Myrs)
    pulsarMagneticFieldDecayMassscale                               = 0.025;                                                                            // Mass scale on which magnetic field decays during accretion (solar masses)
    pulsarLog10MinimumMagneticField                                 = 8.0;                                                                              // log10 of the minimum pulsar magnetic field in Gauss


    // Rotational velocity distribution options
    rotationalVelocityDistribution                                  = ROTATIONAL_VELOCITY_DISTRIBUTION::ZERO;
    rotationalVelocityDistributionString                            = ROTATIONAL_VELOCITY_DISTRIBUTION_LABEL.at(rotationalVelocityDistribution);


	//JIM BARRETT -- 06/07/2016 -- adding options to sample over some hyperparameters
	sampleKickVelocitySigma                                         = false;
	sampleKickVelocitySigmaMin                                      = 0.0;
	sampleKickVelocitySigmaMax                                      = 400.0;

	sampleKickDirectionPower                                        = false;
	sampleKickDirectionPowerMin                                     = -10.0;
	sampleKickDirectionPowerMax                                     = 10.0;

	sampleCommonEnvelopeAlpha                                       = false;
	sampleCommonEnvelopeAlphaMin                                    = 0.0;
	sampleCommonEnvelopeAlphaMax                                    = 5.0;

	sampleWolfRayetMultiplier                                       = false;
	sampleWolfRayetMultiplierMin                                    = 0.0;
	sampleWolfRayetMultiplierMax                                    = 5.0;

	sampleLuminousBlueVariableMultiplier                            = false;
	sampleLuminousBlueVariableMultiplierMin                         = 1.0;
	sampleLuminousBlueVariableMultiplierMax                         = 12.0;


	// grids

	gridFilename                                                    = "";                                                                               // default is no grid file


    // debug and logging options

    debugLevel                                                      = 0;                                                                                // default debug level
    debugClasses.clear();                                                                                                                               // default debug classes

    logLevel                                                        = 0;                                                                                // default log level
    logClasses.clear();                                                                                                                                 // default debug classes
    
    logfileNamePrefix                                               = "";                                                                               // default prefix for all log files
    logfileDelimiter                                                = DELIMITER::TAB;                                                                   // default field delimiter for all log files
    logfileDelimiterString                                          = DELIMITERLabel.at(logfileDelimiter);                                              // delimiter name - for options

    logfileDefinitionsFilename                                      = "";                                                                               // default is no log file definitions file


    // SSE options
    logfileSSEParameters                                            = get<0>(LOGFILE_DESCRIPTOR.at(LOGFILE::SSE_PARAMETERS));                           // get filename from constants.h

    singleStarMassSteps                                             = 100;
    singleStarMassMin                                               = 5.0;
    singleStarMassMax                                               = 100.0;


    // BSE options
    logfileBSESystemParameters                                      = get<0>(LOGFILE_DESCRIPTOR.at(LOGFILE::BSE_SYSTEM_PARAMETERS));                    // get default filename from constants.h
    logfileBSEDetailedOutput                                        = get<0>(LOGFILE_DESCRIPTOR.at(LOGFILE::BSE_DETAILED_OUTPUT));                      // get default filename from constants.h
    logfileBSEDoubleCompactObjects                                  = get<0>(LOGFILE_DESCRIPTOR.at(LOGFILE::BSE_DOUBLE_COMPACT_OBJECTS));               // get default filename from constants.h
    logfileBSESupernovae                                            = get<0>(LOGFILE_DESCRIPTOR.at(LOGFILE::BSE_SUPERNOVAE));                           // get default filename from constants.h
    logfileBSECommonEnvelopes                                       = get<0>(LOGFILE_DESCRIPTOR.at(LOGFILE::BSE_COMMON_ENVELOPES));                     // get default filename from constants.h
    logfileBSEBeBinaries                                            = get<0>(LOGFILE_DESCRIPTOR.at(LOGFILE::BSE_BE_BINARIES));                          // get default filename from constants.h
    logfileBSEPulsarEvolution                                       = get<0>(LOGFILE_DESCRIPTOR.at(LOGFILE::BSE_PULSAR_EVOLUTION));                     // get default filename from constants.h
}




/*
 * Read and process command line arguments using BOOST
 *
 *
 * COMMANDLINE_STATUS CommandLineSorter(int argc, char* argv[])
 *
 * @param   [IN]    argc                        Argument count - number of strings pointed to by parameter argv
 * @param   [IN]    argv                        Array of strings (char* actually) - fisrt string is program name
 * @return                                      Status
 */
COMMANDLINE_STATUS Options::CommandLineSorter(int argc, char* argv[]) {

    COMMANDLINE_STATUS programStatus = COMMANDLINE_STATUS::CONTINUE;

    
    // create default strings for vector<string> types (too hard to do inline)

    std::ostringstream ss;

    // debug classes
    string defaultDebugClasses;
    ss << "";
    for (auto debugClass = debugClasses.begin(); debugClass != debugClasses.end(); ++debugClass) ss << *debugClass << ",";
    defaultDebugClasses = ss.str();
    if (defaultDebugClasses.size() > 0) defaultDebugClasses.erase(defaultDebugClasses.size() - 1);

    // log classes
    string defaultLogClasses;
    ss << "";
    for (auto logClass = logClasses.begin(); logClass != logClasses.end(); ++logClass) ss << *logClass << ",";
    defaultLogClasses = ss.str();
    if (defaultLogClasses.size() > 0) defaultLogClasses.erase(defaultLogClasses.size() - 1);


    // create and populate program options

    try {

        // Create program options object
        po::options_description desc("Options", 128, 64);

        // add options
        desc.add_options()

	        // JR: todo: should make these names consistent ( case, hyphenated, camelCase... )

		    // JR: todo: some of the strings below declare the default value for the option - and some of them are wrong
		    // (probably have become wrong over time).  I think we should either not show the default value, or if we do
		    // then construct the string with the default value.  The second option is more work...

            // switches
		    ("help,h",                                                      po::bool_switch(), "Print this help message")
		    ("version,v",                                                   po::bool_switch(), "Print COMPAS version string")

		    // boolean options - alphabetically

            // AVG - 17/03/2020 - Floor will uncomment when tested.
            /*
            ("AIS-exploratory-phase",                                       po::value<bool>(&AISexploratoryPhase)->default_value(AISexploratoryPhase)->implicit_value(true),                                                            ("Run exploratory phase of STROOPWAFEL (default = " + std::string(AISexploratoryPhase ? "TRUE" : "FALSE") + ")").c_str())
		    ("AIS-Hubble",                                                  po::value<bool>(&AIShubble)->default_value(AIShubble)->implicit_value(true),                                                                                ("Excluding not in Hubble time mergers selection in exploratory phase of STROOPWAFEL (default = " + std::string(AIShubble ? "TRUE" : "FALSE") + ")").c_str())
		    ("AIS-Pessimistic",                                             po::value<bool>(&AISpessimistic)->default_value(AISpessimistic)->implicit_value(true),                                                                      ("Optimistic or Pessimistic selection in exploratory phase of STROOPWAFEL (default = " + std::string(AISpessimistic ? "TRUE" : "FALSE") + ")").c_str())
		    ("AIS-refinement-phase",                                        po::value<bool>(&AISrefinementPhase)->default_value(AISrefinementPhase)->implicit_value(true),                                                              ("Run main sampling phase (step2) of STROOPWAFEL (default = " + std::string(AISrefinementPhase ? "TRUE" : "FALSE") + ")").c_str())
		    ("AIS-RLOF",                                                    po::value<bool>(&AISrlof)->default_value(AISrlof)->implicit_value(true),                                                                                    ("RLOFSecondaryZAMS selection in exploratory phase of STROOPWAFEL (default = " + std::string(AISrlof ? "TRUE" : "FALSE") + ")").c_str())
            */

		    ("allow-rlof-at-birth",                                         po::value<bool>(&allowRLOFAtBirth)->default_value(allowRLOFAtBirth)->implicit_value(true),                                                                  ("Allow binaries that have one or both stars in RLOF at birth to evolve (default = " + std::string(allowRLOFAtBirth ? "TRUE" : "FALSE") + ")").c_str())
		    ("allow-touching-at-birth",                                     po::value<bool>(&allowTouchingAtBirth)->default_value(allowTouchingAtBirth)->implicit_value(true),                                                          ("Allow binaries that are touching at birth to evolve (default = " + std::string(allowTouchingAtBirth ? "TRUE" : "FALSE") + ")").c_str())

			("alwaysStableCaseBBBCFlag",                                    po::value<bool>(&alwaysStableCaseBBBCFlag)->default_value(alwaysStableCaseBBBCFlag)->implicit_value(true),                                                  ("Choose case BB/BC mass transfer to be always stable (default = " + std::string(alwaysStableCaseBBBCFlag ? "TRUE" : "FALSE") + ")").c_str())
			("angularMomentumConservationDuringCircularisation",            po::value<bool>(&angularMomentumConservationDuringCircularisation)->default_value(angularMomentumConservationDuringCircularisation)->implicit_value(true),  ("Conserve angular momentum when binary is circularised when entering a Mass Transfer episode (default = " + std::string(angularMomentumConservationDuringCircularisation ? "TRUE" : "FALSE") + ")").c_str())
			// AVG - 17/03/2020 - Serena will uncomment when tested.
            // ("BeBinaries",                                                  po::value<bool>(&beBinaries)->default_value(beBinaries)->implicit_value(true),                                                                              ("Enable Be Binaries study (default = " + std::string(beBinaries ? "TRUE" : "FALSE") + ")").c_str())
			("circulariseBinaryDuringMassTransfer",                         po::value<bool>(&circulariseBinaryDuringMassTransfer)->default_value(circulariseBinaryDuringMassTransfer)->implicit_value(true),                            ("Circularise binary when it enters a Mass Transfer episode (default = " + std::string(circulariseBinaryDuringMassTransfer ? "TRUE" : "FALSE") + ")").c_str())
		    ("common-envelope-allow-main-sequence-survive",                 po::value<bool>(&allowMainSequenceStarToSurviveCommonEnvelope)->default_value(allowMainSequenceStarToSurviveCommonEnvelope)->implicit_value(true),          ("Allow main sequence stars to survive common envelope evolution (default = " + std::string(allowMainSequenceStarToSurviveCommonEnvelope ? "TRUE" : "FALSE") + ")").c_str())

			("debug-to-file",                                               po::value<bool>(&debugToFile)->default_value(debugToFile)->implicit_value(true),                                                                            ("Write debug statements to file (default = " + std::string(debugToFile ? "TRUE" : "FALSE") + ")").c_str())
		    ("detailedOutput",                                              po::value<bool>(&detailedOutput)->default_value(detailedOutput)->implicit_value(true),                                                                      ("Print detailed output to file (default = " + std::string(detailedOutput ? "TRUE" : "FALSE") + ")").c_str())
			("errors-to-file",                                              po::value<bool>(&errorsToFile)->default_value(errorsToFile)->implicit_value(true),                                                                          ("Write error messages to file (default = " + std::string(errorsToFile ? "TRUE" : "FALSE") + ")").c_str())

		    ("evolve-pulsars",                                              po::value<bool>(&evolvePulsars)->default_value(evolvePulsars)->implicit_value(true),                                                                        ("Evolve pulsars (default = " + std::string(evolvePulsars ? "TRUE" : "FALSE") + ")").c_str())
			("evolve-unbound-systems",                                      po::value<bool>(&evolveUnboundSystems)->default_value(evolveUnboundSystems)->implicit_value(true),                                                          ("Continue evolving stars even if the binary is disrupted (default = " + std::string(evolveUnboundSystems ? "TRUE" : "FALSE") + ")").c_str())

            ("forceCaseBBBCStabilityFlag",                                  po::value<bool>(&forceCaseBBBCStabilityFlag)->default_value(forceCaseBBBCStabilityFlag)->implicit_value(true),                                              ("Force case BB/BC mass transfer to be only stable or unstable (default = " + std::string(forceCaseBBBCStabilityFlag ? "TRUE" : "FALSE") + ")").c_str())
			("lambda-calculation-every-timeStep",                           po::value<bool>(&lambdaCalculationEveryTimeStep)->default_value(lambdaCalculationEveryTimeStep)->implicit_value(true),                                      ("Calculate all values of lambda at each timestep (default = " + std::string(lambdaCalculationEveryTimeStep ? "TRUE" : "FALSE") + ")").c_str())
   		   	("massTransfer",                                                po::value<bool>(&useMassTransfer)->default_value(useMassTransfer)->implicit_value(true),                                                                    ("Enable mass transfer (default = " + std::string(useMassTransfer ? "TRUE" : "FALSE") + ")").c_str())
		    ("pair-instability-supernovae",                                 po::value<bool>(&usePairInstabilitySupernovae)->default_value(usePairInstabilitySupernovae)->implicit_value(true),                                          ("Enable pair instability supernovae (PISN) (default = " + std::string(usePairInstabilitySupernovae ? "TRUE" : "FALSE") + ")").c_str())
            ("populationDataPrinting",                                      po::value<bool>(&populationDataPrinting)->default_value(populationDataPrinting)->implicit_value(true),                                                      ("Print details of population (default = " + std::string(populationDataPrinting ? "TRUE" : "FALSE") + ")").c_str())
		    ("print-bool-as-string",                                        po::value<bool>(&printBoolAsString)->default_value(printBoolAsString)->implicit_value(true),                                                                ("Print boolean properties as 'TRUE' or 'FALSE' (default = " + std::string(printBoolAsString ? "TRUE" : "FALSE") + ")").c_str())
		    ("pulsational-pair-instability",                                po::value<bool>(&usePulsationalPairInstability)->default_value(usePulsationalPairInstability)->implicit_value(true),                                        ("Enable mass loss due to pulsational-pair-instability (PPI) (default = " + std::string(usePulsationalPairInstability ? "TRUE" : "FALSE") + ")").c_str())
		    ("quiet",                                                       po::value<bool>(&quiet)->default_value(quiet)->implicit_value(true),                                                                                        ("Suppress printing (default = " + std::string(quiet ? "TRUE" : "FALSE") + ")").c_str())
			("revised-energy-formalism-Nandez-Ivanova",                     po::value<bool>(&revisedEnergyFormalismNandezIvanova)->default_value(revisedEnergyFormalismNandezIvanova)->implicit_value(true),                            ("Enable revised energy formalism (default = " + std::string(revisedEnergyFormalismNandezIvanova ? "TRUE" : "FALSE") + ")").c_str())

            // AVG
            /*
			("sample-common-envelope-alpha",                                po::value<bool>(&sampleCommonEnvelopeAlpha)->default_value(sampleCommonEnvelopeAlpha)->implicit_value(true),                                                ("Sample over common envelope alpha (default = " + std::string(sampleCommonEnvelopeAlpha ? "TRUE" : "FALSE") + ")").c_str())
			("sample-kick-direction-power",                                 po::value<bool>(&sampleKickDirectionPower)->default_value(sampleKickDirectionPower)->implicit_value(true),                                                  ("Sample over kick direction powerlaw exponent (default = " + std::string(sampleKickDirectionPower ? "TRUE" : "FALSE") + ")").c_str())
			("sample-kick-velocity-sigma",                                  po::value<bool>(&sampleKickVelocitySigma)->default_value(sampleKickVelocitySigma)->implicit_value(true),                                                    ("Sample over Kick Velocity Sigma (default = " + std::string(sampleKickVelocitySigma ? "TRUE" : "FALSE") + ")").c_str())
			("sample-luminous-blue-variable-multiplier",                    po::value<bool>(&sampleLuminousBlueVariableMultiplier)->default_value(sampleLuminousBlueVariableMultiplier)->implicit_value(true),                          ("Sample over multiplicative constant from LBV mass loss (default = " + std::string(sampleLuminousBlueVariableMultiplier ? "TRUE" : "FALSE") + ")").c_str())
			("sample-wolf-rayet-multiplier",                                po::value<bool>(&sampleWolfRayetMultiplier)->default_value(sampleWolfRayetMultiplier)->implicit_value(true),                                                ("Sample over WR winds multiplicative constant (default = " + std::string(sampleWolfRayetMultiplier ? "TRUE" : "FALSE") + ")").c_str())
            */

            ("single-star",                                                 po::value<bool>(&singleStar)->default_value(singleStar)->implicit_value(true),                                                                              ("Evolve single star(s) (default = " + std::string(singleStar ? "TRUE" : "FALSE") + ")").c_str())

		    ("use-mass-loss",                                               po::value<bool>(&useMassLoss)->default_value(useMassLoss)->implicit_value(true),                                                                            ("Enable mass loss (default = " + std::string(useMassLoss ? "TRUE" : "FALSE") + ")").c_str())

			("zeta-calculation-every-timestep",                             po::value<bool>(&zetaCalculationEveryTimeStep)->default_value(zetaCalculationEveryTimeStep)->implicit_value(true),                                          ("Calculate all values of MT zetas at each timestep (default = " + std::string(zetaCalculationEveryTimeStep ? "TRUE" : "FALSE") + ")").c_str())


			// numerical options - alphabetically grouped by type

			// unsigned long
		    ("random-seed",                                                 po::value<unsigned long>(&randomSeed)->default_value(randomSeed),                                                                                           ("Random seed to use (default = " + std::to_string(randomSeed) + ")").c_str())

		    // int
			("debug-level",                                                 po::value<int>(&debugLevel)->default_value(debugLevel),                                                                                                     ("Determines which print statements are displayed for debugging (default = " + std::to_string(debugLevel) + ")").c_str())
		    ("log-level",                                                   po::value<int>(&logLevel)->default_value(logLevel),                                                                                                         ("Determines which print statements are included in the logfile (default = " + std::to_string(logLevel) + ")").c_str())
		    ("maximum-number-timestep-iterations",                                   po::value<int>(&maxNumberOfTimestepIterations)->default_value(maxNumberOfTimestepIterations),                                                               ("Maximum number of timesteps to evolve binary before giving up (default = " + std::to_string(maxNumberOfTimestepIterations) + ")").c_str())
			// AVG - 17/03/2020 - Floor will uncomment when tested.
            // ("nbatches-used",                                               po::value<int>(&nBatchesUsed)->default_value(nBatchesUsed),                                                                                                 ("Number of batches used, for STROOPWAFEL (AIS), -1 = not required (default = " + std::to_string(nBatchesUsed) + ")").c_str())
			("number-of-binaries,n",                                        po::value<int>(&nBinaries)->default_value(nBinaries),                                                                                                       ("Specify the number of binaries to simulate (default = " + std::to_string(nBinaries) + ")").c_str())
			("single-star-mass-steps",                                      po::value<int>(&singleStarMassSteps)->default_value(singleStarMassSteps),                                                                                   ("Specify the number of mass steps for single star evolution (default = " + std::to_string(singleStarMassSteps) + ")").c_str())

		    // double
		    ("common-envelope-alpha",                                       po::value<double>(&commonEnvelopeAlpha)->default_value(commonEnvelopeAlpha),                                                                                ("Common Envelope efficiency alpha (default = " + std::to_string(commonEnvelopeAlpha) + ")").c_str())
		    ("common-envelope-alpha-thermal",                               po::value<double>(&commonEnvelopeAlphaThermal)->default_value(commonEnvelopeAlphaThermal),                                                                  ("Defined such that lambda = alpha_th * lambda_b + (1.0 - alpha_th) * lambda_g (default = " + std::to_string(commonEnvelopeAlphaThermal) + ")").c_str())
		    ("common-envelope-lambda",                                      po::value<double>(&commonEnvelopeLambda)->default_value(commonEnvelopeLambda),                                                                              ("Common Envelope lambda (default = " + std::to_string(commonEnvelopeLambda) + ")").c_str())
		    ("common-envelope-lambda-multiplier",                           po::value<double>(&commonEnvelopeLambdaMultiplier)->default_value(commonEnvelopeLambdaMultiplier),                                                          ("Multiply lambda by some constant (default = " + std::to_string(commonEnvelopeLambdaMultiplier) + ")").c_str())
	    	("common-envelope-mass-accretion-constant",                     po::value<double>(&commonEnvelopeMassAccretionConstant)->default_value(commonEnvelopeMassAccretionConstant),                                                ("Value of mass accreted by NS/BH during common envelope evolution if assuming all NS/BH accrete same amount of mass (common-envelope-mass-accretion-prescription CONSTANT). Ignored otherwise (default = " + std::to_string(commonEnvelopeMassAccretionConstant) + ")").c_str())
		    ("common-envelope-mass-accretion-max",                          po::value<double>(&commonEnvelopeMassAccretionMax)->default_value(commonEnvelopeMassAccretionMax),                                                          ("Maximum amount of mass accreted by NS/BHs during common envelope evolution in solar masses (default = " + std::to_string(commonEnvelopeMassAccretionMax) + ")").c_str())
		    ("common-envelope-mass-accretion-min",                          po::value<double>(&commonEnvelopeMassAccretionMin)->default_value(commonEnvelopeMassAccretionMin),                                                          ("Minimum amount of mass accreted by NS/BHs during common envelope evolution in solar masses (default = " + std::to_string(commonEnvelopeMassAccretionMin) + ")").c_str())
		    ("common-envelope-recombination-energy-density",                po::value<double>(&commonEnvelopeRecombinationEnergyDensity)->default_value(commonEnvelopeRecombinationEnergyDensity),                                      ("Recombination energy density in erg/g (default = " + std::to_string(commonEnvelopeRecombinationEnergyDensity) + ")").c_str())
		    ("common-envelope-slope-Kruckow",                               po::value<double>(&commonEnvelopeSlopeKruckow)->default_value(commonEnvelopeSlopeKruckow),                                                                  ("Common Envelope slope for Kruckow lambda (default = " + std::to_string(commonEnvelopeSlopeKruckow) + ")").c_str())

            // AVG - 17/03/2020 - Uncomment mass-ratio options when fully implemented
            /*
            ("critical-mass-ratio-giant-degenerate-accretor",               po::value<double>(&massTransferCriticalMassRatioGiantDegenerateAccretor)->default_value(massTransferCriticalMassRatioGiantDegenerateAccretor),              ("Critical mass ratio for MT from a giant star (default = " + std::to_string(massTransferCriticalMassRatioGiantDegenerateAccretor) + ") Specify both giant flags to use. 0 is always stable, <0 is disabled").c_str())
            ("critical-mass-ratio-giant-non-degenerate-accretor",           po::value<double>(&massTransferCriticalMassRatioGiantNonDegenerateAccretor)->default_value(massTransferCriticalMassRatioGiantNonDegenerateAccretor),        ("Critical mass ratio for MT from a giant star (default = " + std::to_string(massTransferCriticalMassRatioGiantNonDegenerateAccretor) + ") Specify both giant flags to use. 0 is always stable, <0 is disabled").c_str())
            ("critical-mass-ratio-helium-giant-degenerate-accretor",        po::value<double>(&massTransferCriticalMassRatioHeliumGiantDegenerateAccretor)->default_value(massTransferCriticalMassRatioHeliumGiantDegenerateAccretor),  ("Critical mass ratio for MT from a helium giant star (default = " + std::to_string(massTransferCriticalMassRatioHeliumGiantDegenerateAccretor) + ") Specify both helium giant flags to use. 0 is always stable, <0 is disabled").c_str())
            ("critical-mass-ratio-helium-giant-non-degenerate-accretor",    po::value<double>(&massTransferCriticalMassRatioHeliumGiantNonDegenerateAccretor)->default_value(massTransferCriticalMassRatioHeliumGiantNonDegenerateAccretor),    ("Critical mass ratio for MT from a helium giant star (default = " + std::to_string(massTransferCriticalMassRatioHeliumGiantNonDegenerateAccretor) + ") Specify both helium giant flags to use. 0 is always stable, <0 is disabled").c_str())
            ("critical-mass-ratio-helium-HG-degenerate-accretor",           po::value<double>(&massTransferCriticalMassRatioHeliumHGDegenerateAccretor)->default_value(massTransferCriticalMassRatioHeliumHGDegenerateAccretor),        ("Critical mass ratio for MT from a helium HG star (default = " + std::to_string(massTransferCriticalMassRatioHeliumHGDegenerateAccretor) + ") Specify both helium HG flags to use. 0 is always stable, <0 is disabled").c_str())
            ("critical-mass-ratio-helium-HG-non-degenerate-accretor",       po::value<double>(&massTransferCriticalMassRatioHeliumHGNonDegenerateAccretor)->default_value(massTransferCriticalMassRatioHeliumHGNonDegenerateAccretor),  ("Critical mass ratio for MT from a helium HG star (default = " + std::to_string(massTransferCriticalMassRatioHeliumHGNonDegenerateAccretor) + ") Specify both helium HG flags to use. 0 is always stable, <0 is disabled").c_str())
            ("critical-mass-ratio-helium-MS-degenerate-accretor",           po::value<double>(&massTransferCriticalMassRatioHeliumMSDegenerateAccretor)->default_value(massTransferCriticalMassRatioHeliumMSDegenerateAccretor),        ("Critical mass ratio for MT from a helium MS star (default = " + std::to_string(massTransferCriticalMassRatioHeliumMSDegenerateAccretor) + ") Specify both helium MS flags to use. 0 is always stable, <0 is disabled").c_str())
            ("critical-mass-ratio-helium-MS-non-degenerate-accretor",       po::value<double>(&massTransferCriticalMassRatioHeliumMSNonDegenerateAccretor)->default_value(massTransferCriticalMassRatioHeliumMSNonDegenerateAccretor),  ("Critical mass ratio for MT from a helium MS star (default = " + std::to_string(massTransferCriticalMassRatioHeliumMSNonDegenerateAccretor) + ") Specify both helium MS flags to use. 0 is always stable, <0 is disabled").c_str())
            ("critical-mass-ratio-HG-degenerate-accretor",                  po::value<double>(&massTransferCriticalMassRatioHGDegenerateAccretor)->default_value(massTransferCriticalMassRatioHGDegenerateAccretor),                    ("Critical mass ratio for MT from a HG star (default = " + std::to_string(massTransferCriticalMassRatioHGDegenerateAccretor) + ") Specify both HG flags to use. 0 is always stable, <0 is disabled").c_str())
            ("critical-mass-ratio-HG-non-degenerate-accretor",              po::value<double>(&massTransferCriticalMassRatioHGNonDegenerateAccretor)->default_value(massTransferCriticalMassRatioHGNonDegenerateAccretor),              ("Critical mass ratio for MT from a HG star (default = " + std::to_string(massTransferCriticalMassRatioHGNonDegenerateAccretor) + ") Specify both HG flags to use. 0 is always stable, <0 is disabled").c_str())
            ("critical-mass-ratio-MS-high-mass-degenerate-accretor",        po::value<double>(&massTransferCriticalMassRatioMSHighMassDegenerateAccretor)->default_value(massTransferCriticalMassRatioMSHighMassDegenerateAccretor),    ("Critical mass ratio for MT from a MS star to a degenerate accretor (default = " + std::to_string(massTransferCriticalMassRatioMSHighMassDegenerateAccretor) + " Specify both MS high mass flags to use. 0 is always stable, <0 is disabled").c_str())
            ("critical-mass-ratio-MS-high-mass-non-degenerate-accretor",    po::value<double>(&massTransferCriticalMassRatioMSHighMassNonDegenerateAccretor)->default_value(massTransferCriticalMassRatioMSHighMassNonDegenerateAccretor),  ("Critical mass ratio for MT from a MS star (default = " + std::to_string(massTransferCriticalMassRatioMSHighMassNonDegenerateAccretor) + ") Specify both MS high mass flags to use. 0 is always stable, <0 is disabled").c_str())
            ("critical-mass-ratio-MS-low-mass-degenerate-accretor",         po::value<double>(&massTransferCriticalMassRatioMSLowMassDegenerateAccretor)->default_value(massTransferCriticalMassRatioMSLowMassDegenerateAccretor),      ("Critical mass ratio for MT from a MS star to a degenerate accretor (default = " + std::to_string(massTransferCriticalMassRatioMSLowMassDegenerateAccretor) + " Specify both MS low mass flags to use. 0 is always stable, <0 is disabled").c_str())
            ("critical-mass-ratio-MS-low-mass-non-degenerate-accretor",     po::value<double>(&massTransferCriticalMassRatioMSLowMassNonDegenerateAccretor)->default_value(massTransferCriticalMassRatioMSLowMassNonDegenerateAccretor),    ("Critical mass ratio for MT from a MS star (default = " + std::to_string(massTransferCriticalMassRatioMSLowMassNonDegenerateAccretor) + ") Specify both MS low mass flags to use. 0 is always stable, <0 is disabled").c_str())
            ("critical-mass-ratio-white-dwarf-degenerate-accretor",         po::value<double>(&massTransferCriticalMassRatioWhiteDwarfDegenerateAccretor)->default_value(massTransferCriticalMassRatioWhiteDwarfDegenerateAccretor),    ("Critical mass ratio for MT from a white dwarf (default = " + std::to_string(massTransferCriticalMassRatioWhiteDwarfDegenerateAccretor) + ") Specify both white dwarf flags to use. 0 is always stable, <0 is disabled").c_str())
            ("critical-mass-ratio-white-dwarf-non-degenerate-accretor",     po::value<double>(&massTransferCriticalMassRatioWhiteDwarfNonDegenerateAccretor)->default_value(massTransferCriticalMassRatioWhiteDwarfNonDegenerateAccretor),  ("Critical mass ratio for MT from a white dwarf (default = " + std::to_string(massTransferCriticalMassRatioWhiteDwarfNonDegenerateAccretor) + ") Specify both white dwarf flags to use. 0 is always stable, <0 is disabled").c_str())
            */

		    ("eccentricity-max",                                            po::value<double>(&eccentricityDistributionMax)->default_value(eccentricityDistributionMax),                                                                ("Maximum eccentricity to generate (default = " + std::to_string(eccentricityDistributionMax) + ")").c_str())
		    ("eccentricity-min",                                            po::value<double>(&eccentricityDistributionMin)->default_value(eccentricityDistributionMin),                                                                ("Minimum eccentricity to generate (default = " + std::to_string(eccentricityDistributionMin) + ")").c_str())
			("eddington-accretion-factor",                                  po::value<double>(&eddingtonAccretionFactor)->default_value(eddingtonAccretionFactor),                                                                      ("Multiplication factor for eddington accretion for NS & BH, i.e. >1 is super-eddington and 0. is no accretion (default = " + std::to_string(eddingtonAccretionFactor) + ")").c_str())

   		    ("fix-dimensionless-kick-velocity",                             po::value<double>(&fixedUK)->default_value(fixedUK),                                                                                                        ("Fix dimensionless kick velocity uk to this value (default = " + std::to_string(fixedUK) + ", -ve values false, +ve values true)").c_str())

		    ("initial-mass-max",                                            po::value<double>(&initialMassFunctionMax)->default_value(initialMassFunctionMax),                                                                          ("Maximum mass (in Msol) to generate using given IMF (default = " + std::to_string(initialMassFunctionMax) + ")").c_str())
		    ("initial-mass-min",                                            po::value<double>(&initialMassFunctionMin)->default_value(initialMassFunctionMin),                                                                          ("Minimum mass (in Msol) to generate using given IMF (default = " + std::to_string(initialMassFunctionMin) + ")").c_str())
		    ("initial-mass-power",                                          po::value<double>(&initialMassFunctionPower)->default_value(initialMassFunctionPower),                                                                      ("Single power law power to generate primary mass using given IMF (default = " + std::to_string(initialMassFunctionPower) + ")").c_str())

            // AVG - 17/03/2020 - Floor will uncomment when tested.
		    // ("kappa-gaussians",                                             po::value<double>(&kappaGaussians)->default_value(kappaGaussians),                                                                                          ("Scaling factor for the width of the Gaussian distributions in STROOPWAFEL main sampling phase (default = " + std::to_string(kappaGaussians) + ")").c_str())
		    ("kick-direction-power",                                        po::value<double>(&kickDirectionPower)->default_value(kickDirectionPower),                                                                                  ("Power for power law kick direction distribution (default = " + std::to_string(kickDirectionPower) + " = isotropic, +ve = polar, -ve = in plane)").c_str())
			("kick-scaling-factor",                                         po::value<double>(&kickScalingFactor)->default_value(kickScalingFactor),                                                                                    ("Arbitrary factor used to scale kicks (default = " + std::to_string(kickScalingFactor) + ")").c_str())
		    ("kick-velocity-max",                                           po::value<double>(&kickVelocityDistributionMaximum)->default_value(kickVelocityDistributionMaximum),                                                        ("Maximum drawn kick velocity in km s^-1. Ignored if < 0. Must be > 0 if using kick-velocity-distribution=FLAT (default = " + std::to_string(kickVelocityDistributionMaximum) + ")").c_str())
		    ("kick-velocity-sigma-CCSN-BH",                                 po::value<double>(&kickVelocityDistributionSigmaCCSN_BH)->default_value(kickVelocityDistributionSigmaCCSN_BH),                                              ("Sigma for chosen kick velocity distribution for black holes (default = " + std::to_string(kickVelocityDistributionSigmaCCSN_BH) + " km s^-1 )").c_str())
		    ("kick-velocity-sigma-CCSN-NS",                                 po::value<double>(&kickVelocityDistributionSigmaCCSN_NS)->default_value(kickVelocityDistributionSigmaCCSN_NS),                                              ("Sigma for chosen kick velocity distribution for neutron stars (default = " + std::to_string(kickVelocityDistributionSigmaCCSN_NS) + " km s^-1 )").c_str())
			("kick-velocity-sigma-ECSN",                                    po::value<double>(&kickVelocityDistributionSigmaForECSN)->default_value(kickVelocityDistributionSigmaForECSN),                                              ("Sigma for chosen kick velocity distribution for ECSN (default = " + std::to_string(kickVelocityDistributionSigmaForECSN) + " km s^-1 )").c_str())
			("kick-velocity-sigma-USSN",                                    po::value<double>(&kickVelocityDistributionSigmaForUSSN)->default_value(kickVelocityDistributionSigmaForUSSN),                                              ("Sigma for chosen kick velocity distribution for USSN (default = " + std::to_string(kickVelocityDistributionSigmaForUSSN) + " km s^-1 )").c_str())

		    ("luminous-blue-variable-multiplier",                           po::value<double>(&luminousBlueVariableFactor)->default_value(luminousBlueVariableFactor),                                                                  ("Multiplicitive constant for LBV mass loss (default = " + std::to_string(luminousBlueVariableFactor) + ", use 10 for Mennekens & Vanbeveren 2014)").c_str())

		    ("mass-ratio-max",                                              po::value<double>(&massRatioDistributionMax)->default_value(massRatioDistributionMax),                                                                      ("Maximum mass ratio m2/m1 to generate (default = " + std::to_string(massRatioDistributionMax) + ")").c_str())
		    ("mass-ratio-min",                                              po::value<double>(&massRatioDistributionMin)->default_value(massRatioDistributionMin),                                                                      ("Minimum mass ratio m2/m1 to generate (default = " + std::to_string(massRatioDistributionMin) + ")").c_str())
		    ("mass-transfer-fa",                                            po::value<double>(&massTransferFractionAccreted)->default_value(massTransferFractionAccreted),                                                              ("Mass Transfer fraction accreted in FIXED prescription (default = " + std::to_string(massTransferFractionAccreted) + ", fully conservative)").c_str())
		    ("mass-transfer-jloss",                                         po::value<double>(&massTransferJloss)->default_value(massTransferJloss),                                                                                    ("Specific angular momentum with which the non-accreted system leaves the system (default = " + std::to_string(massTransferJloss) + ")").c_str())
			("mass-transfer-thermal-limit-C",                               po::value<double>(&massTransferCParameter)->default_value(massTransferCParameter),                                                                          ("Mass Transfer Thermal rate factor fo the accretor (default = " + std::to_string(massTransferCParameter) + ")").c_str())
		    ("maximum-evolution-time",                                      po::value<double>(&maxEvolutionTime)->default_value(maxEvolutionTime),                                                                                      ("Maximum time to evolve binaries in Myrs (default = " + std::to_string(maxEvolutionTime) + ")").c_str())
		    ("maximum-mass-donor-Nandez-Ivanova",                           po::value<double>(&maximumMassDonorNandezIvanova)->default_value(maximumMassDonorNandezIvanova),                                                            ("Maximum donor mass allowed for the revised common envelope formalism in Msol (default = " + std::to_string(maximumMassDonorNandezIvanova) + ")").c_str())
			("maximum-neutron-star-mass",                                   po::value<double>(&maximumNeutronStarMass)->default_value(maximumNeutronStarMass),                                                                          ("Maximum mass of a neutron star (default = " + std::to_string(maximumNeutronStarMass) + ")").c_str())
            ("MCBUR1",                                                      po::value<double>(&MCBUR1)->default_value(MCBUR1),                                                                          ("MCBUR1: Min core mass at BAGB to avoid fully degenerate CO core  (default = " + std::to_string(MCBUR1) + ")").c_str())
            ("metallicity,z",                                               po::value<double>(&metallicity)->default_value(metallicity),                                                                                                ("Metallicity to use (default " + std::to_string(metallicity) + " Zsol)").c_str())
		    ("minimum-secondary-mass",                                      po::value<double>(&minimumMassSecondary)->default_value(minimumMassSecondary),                                                                              ("Minimum mass of secondary to generate in Msol (default = " + std::to_string(minimumMassSecondary) + ")").c_str())

		    ("neutrino-mass-loss-bh-formation-value",                       po::value<double>(&neutrinoMassLossValueBH)->default_value(neutrinoMassLossValueBH),                                                                        ("Value corresponding to neutrino mass loss assumption (default = " + std::to_string(neutrinoMassLossValueBH) + ")").c_str())

		    ("orbital-period-max",                                          po::value<double>(&periodDistributionMax)->default_value(periodDistributionMax),                                                                            ("Maximum period in days to generate (default = " + std::to_string(periodDistributionMax) + ")").c_str())
		   	("orbital-period-min",                                          po::value<double>(&periodDistributionMin)->default_value(periodDistributionMin),                                                                            ("Minimum period in days to generate (default = " + std::to_string(periodDistributionMin) + ")").c_str())

		    ("PISN-lower-limit",                                            po::value<double>(&pairInstabilityLowerLimit)->default_value(pairInstabilityLowerLimit),                                                                    ("Minimum core mass for PISN (default = " + std::to_string(pairInstabilityLowerLimit) + ")").c_str())
		    ("PISN-upper-limit",                                            po::value<double>(&pairInstabilityUpperLimit)->default_value(pairInstabilityUpperLimit),                                                                    ("Maximum core mass for PISN (default = " + std::to_string(pairInstabilityUpperLimit) + ")").c_str())
		    ("PPI-lower-limit",                                             po::value<double>(&pulsationalPairInstabilityLowerLimit)->default_value(pulsationalPairInstabilityLowerLimit),                                              ("Minimum core mass for PPI (default = " + std::to_string(pulsationalPairInstabilityLowerLimit) + ")").c_str())
		    ("PPI-upper-limit",                                             po::value<double>(&pulsationalPairInstabilityUpperLimit)->default_value(pulsationalPairInstabilityUpperLimit),                                              ("Maximum core mass for PPI (default = " + std::to_string(pulsationalPairInstabilityUpperLimit) + ")").c_str())
		    ("pulsar-birth-magnetic-field-distribution-max",                po::value<double>(&pulsarBirthMagneticFieldDistributionMax)->default_value(pulsarBirthMagneticFieldDistributionMax),                                        ("Maximum (log10) pulsar birth magnetic field (default = " + std::to_string(pulsarBirthMagneticFieldDistributionMax) + ")").c_str())
		    ("pulsar-birth-magnetic-field-distribution-min",                po::value<double>(&pulsarBirthMagneticFieldDistributionMin)->default_value(pulsarBirthMagneticFieldDistributionMin),                                        ("Minimum (log10) pulsar birth magnetic field) (default = " + std::to_string(pulsarBirthMagneticFieldDistributionMin) + ")").c_str())
		    ("pulsar-birth-spin-period-distribution-max",                   po::value<double>(&pulsarBirthSpinPeriodDistributionMax)->default_value(pulsarBirthSpinPeriodDistributionMax),                                              ("Maximum pulsar birth spin period in ms (default = " + std::to_string(pulsarBirthSpinPeriodDistributionMax) + ")").c_str())
		    ("pulsar-birth-spin-period-distribution-min",                   po::value<double>(&pulsarBirthSpinPeriodDistributionMin)->default_value(pulsarBirthSpinPeriodDistributionMin),                                              ("Minimum pulsar birth spin period in ms (default = " + std::to_string(pulsarBirthSpinPeriodDistributionMin) + ")").c_str())
		    ("pulsar-magnetic-field-decay-massscale",                       po::value<double>(&pulsarMagneticFieldDecayMassscale)->default_value(pulsarMagneticFieldDecayMassscale),                                                    ("Mass scale on which magnetic field decays during accretion in solar masses (default = " + std::to_string(pulsarMagneticFieldDecayMassscale) + ")").c_str())
		    ("pulsar-magnetic-field-decay-timescale",                       po::value<double>(&pulsarMagneticFieldDecayTimescale)->default_value(pulsarMagneticFieldDecayTimescale),                                                    ("Timescale on which magnetic field decays in Myrs (default = " + std::to_string(pulsarMagneticFieldDecayTimescale) + ")").c_str())
		    ("pulsar-minimum-magnetic-field",                               po::value<double>(&pulsarLog10MinimumMagneticField)->default_value(pulsarLog10MinimumMagneticField),                                                        ("log10 of the minimum pulsar magnetic field in Gauss (default = " + std::to_string(pulsarLog10MinimumMagneticField) + ")").c_str())

            // AVG
            /*
			("sample-common-envelope-alpha-max",                            po::value<double>(&sampleCommonEnvelopeAlphaMax)->default_value(sampleCommonEnvelopeAlphaMax),                                                              ("Maximum for Uniform sampling over common envelope alpha (default = " + std::to_string(sampleCommonEnvelopeAlphaMax) + ")").c_str())
			("sample-common-envelope-alpha-min",                            po::value<double>(&sampleCommonEnvelopeAlphaMin)->default_value(sampleCommonEnvelopeAlphaMin),                                                              ("Minimum for Uniform sampling over common envelope alpha (default = " + std::to_string(sampleCommonEnvelopeAlphaMin) + ")").c_str())
			("sample-kick-direction-power-max",                             po::value<double>(&sampleKickDirectionPowerMax)->default_value(sampleKickDirectionPowerMax),                                                                ("Maximum for Uniform sampling over kick direction powerlaw exponent (default = " + std::to_string(sampleKickDirectionPowerMax) + ")").c_str())
			("sample-kick-direction-power-min",                             po::value<double>(&sampleKickDirectionPowerMin)->default_value(sampleKickDirectionPowerMin),                                                                ("Minimum for Uniform sampling over kick direction powerlaw exponent (default = " + std::to_string(sampleKickDirectionPowerMin) + ")").c_str())
			("sample-kick-velocity-sigma-max",                              po::value<double>(&sampleKickVelocitySigmaMax)->default_value(sampleKickVelocitySigmaMax),                                                                  ("Maximum for Uniform sampling over kick velocity sigma (default = " + std::to_string(sampleKickVelocitySigmaMax) + ")").c_str())
			("sample-kick-velocity-sigma-min",                              po::value<double>(&sampleKickVelocitySigmaMin)->default_value(sampleKickVelocitySigmaMin),                                                                  ("Minimum for Uniform sampling over kick velocity sigma (default = " + std::to_string(sampleKickVelocitySigmaMin) + ")").c_str())
			("sample-luminous-blue-variable-multiplier-max",                po::value<double>(&sampleLuminousBlueVariableMultiplierMax)->default_value(sampleLuminousBlueVariableMultiplierMax),                                        ("Maximum for Uniform sampling over multiplicative constant for LBV mass loss (default = " + std::to_string(sampleLuminousBlueVariableMultiplierMax) + ")").c_str())
			("sample-luminous-blue-variable-multiplier-min",                po::value<double>(&sampleLuminousBlueVariableMultiplierMin)->default_value(sampleLuminousBlueVariableMultiplierMin),                                        ("Minimum for Uniform sampling over multiplicative constant for LBV mass loss (default = " + std::to_string(sampleLuminousBlueVariableMultiplierMin) + ")").c_str())
			("sample-wolf-rayet-multiplier-max",                            po::value<double>(&sampleWolfRayetMultiplierMax)->default_value(sampleWolfRayetMultiplierMax),                                                              ("Maximum for Uniform sampling over multiplicative constant for WR winds (default = " + std::to_string(sampleWolfRayetMultiplierMax) + ")").c_str())
			("sample-wolf-rayet-multiplier-min",                            po::value<double>(&sampleWolfRayetMultiplierMin)->default_value(sampleWolfRayetMultiplierMin),                                                              ("Minimum for Uniform sampling over multiplicative constant for WR winds (default = " + std::to_string(sampleWolfRayetMultiplierMin) + ")").c_str())
		    */
            ("semi-major-axis-max",                                         po::value<double>(&semiMajorAxisDistributionMax)->default_value(semiMajorAxisDistributionMax),                                                              ("Maximum semi major axis in AU to generate (default = " + std::to_string(semiMajorAxisDistributionMax) + ")").c_str())
		    ("semi-major-axis-min",                                         po::value<double>(&semiMajorAxisDistributionMin)->default_value(semiMajorAxisDistributionMin),                                                              ("Minimum semi major axis in AU to generate (default = " + std::to_string(semiMajorAxisDistributionMin) + ")").c_str())
		    ("single-star-mass-max",                                        po::value<double>(&singleStarMassMax)->default_value(singleStarMassMax),                                                                                    ("Maximum mass (in Msol) for single star evolution (default = " + std::to_string(singleStarMassMax) + ")").c_str())
            ("single-star-mass-min",                                        po::value<double>(&singleStarMassMin)->default_value(singleStarMassMin),                                                                                    ("Minimum mass (in Msol) for single star evolution (default = " + std::to_string(singleStarMassMin) + ")").c_str())

		    ("wolf-rayet-multiplier",                                       po::value<double>(&wolfRayetFactor)->default_value(wolfRayetFactor),                                                                                        ("Multiplicitive constant for WR winds (default = " + std::to_string(wolfRayetFactor) + ")").c_str())

		    ("zeta-adiabatic-arbitrary",                                    po::value<double>(&zetaAdiabaticArbitrary)->default_value(zetaAdiabaticArbitrary),                                                                          ("Value of mass-radius exponent zeta adiabatic (default = " + std::to_string(zetaAdiabaticArbitrary) + ")").c_str())
		    ("zeta-radiative-envelope-giant",                               po::value<double>(&zetaRadiativeEnvelopeGiant)->default_value(zetaRadiativeEnvelopeGiant),                                                                                  ("Value of mass-radius exponent zeta for radiative envelope giants (default = " + std::to_string(zetaRadiativeEnvelopeGiant) + ")").c_str())
		    ("zeta-main-sequence",                                          po::value<double>(&zetaMainSequence)->default_value(zetaMainSequence),                                                                                      ("Value of mass-radius exponent zeta on the main sequence (default = " + std::to_string(zetaMainSequence) + ")").c_str())


		    // string options - alphabetically
            // AVG - 17/03/2020 - Floor will uncomment when tested.
            // ("AIS-DCOtype",                                                 po::value<string>(&AISDCOtypeString)->default_value(AISDCOtypeString),                                                                                      ("DCO type selection in exploratory phase of STROOPWAFEL, (options: ALL, BBH, BNS or BHNS), default = " + AISDCOtypeString + ")").c_str())

		  	("black-hole-kicks",                                            po::value<string>(&blackHoleKicksString)->default_value(blackHoleKicksString),                                                                              ("Black hole kicks relative to NS kicks (options: FULL, REDUCED, ZERO, FALLBACK), default = " + blackHoleKicksString + ")").c_str())

		  	("chemically-homogeneous-evolution",                            po::value<string>(&cheString)->default_value(cheString),                                                                                                    ("Chemically Homogeneous Evolution (options: NONE, OPTIMISTIC, PESSIMISTIC), default = " + cheString + ")").c_str())

			("common-envelope-lambda-prescription",                         po::value<string>(&commonEnvelopeLambdaPrescriptionString)->default_value(commonEnvelopeLambdaPrescriptionString),                                          ("CE lambda prescription (options: LAMBDA_FIXED, LAMBDA_LOVERIDGE, LAMBDA_NANJING, LAMBDA_KRUCKOW, LAMBDA_DEWI), default = " + commonEnvelopeLambdaPrescriptionString + ")").c_str())
		    ("common-envelope-mass-accretion-prescription",                 po::value<string>(&commonEnvelopeMassAccretionPrescriptionString)->default_value(commonEnvelopeMassAccretionPrescriptionString),                            ("Assumption about whether NS/BHs can accrete mass during common envelope evolution (options: ZERO, CONSTANT, UNIFORM, MACLEOD), default = " + commonEnvelopeMassAccretionPrescriptionString + ")").c_str())
        
            ("envelope-state-prescription",                                 po::value<string>(&envelopeStatePrescriptionString)->default_value(envelopeStatePrescriptionString),                                   ("Prescription for whether the envelope is radiative or convective (options: LEGACY, HURLEY, FIXED_TEMPERATURE), default = " + envelopeStatePrescriptionString + ")").c_str())
        
			("stellar-zeta-prescription",                           po::value<string>(&stellarZetaPrescriptionString)->default_value(stellarZetaPrescriptionString),                                              ("Prescription for stellar zeta (default = " + stellarZetaPrescriptionString + ")").c_str())

		    ("eccentricity-distribution,e",                                 po::value<string>(&eccentricityDistributionString)->default_value(eccentricityDistributionString),                                                          ("Initial eccentricity distribution, e (options: ZERO, FIXED, FLAT, THERMALISED, GELLER+2013), default = " + eccentricityDistributionString + ")").c_str())

		    ("fryer-supernova-engine",                                      po::value<string>(&fryerSupernovaEngineString)->default_value(fryerSupernovaEngineString),                                                                  ("If using Fryer et al 2012 fallback prescription, select between 'delayed' and 'rapid' engines (default = " + fryerSupernovaEngineString + ")").c_str())

            ("grid",                                                        po::value<string>(&gridFilename)->default_value(gridFilename)->implicit_value(""),                                                                          ("Grid filename (default = " + gridFilename + ")").c_str())

		    ("initial-mass-function,i",                                     po::value<string>(&initialMassFunctionString)->default_value(initialMassFunctionString),                                                                    ("Initial mass function (options: SALPETER, POWERLAW, UNIFORM, KROUPA), default = " + initialMassFunctionString + ")").c_str())

		    ("kick-direction",                                              po::value<string>(&kickDirectionDistributionString)->default_value(kickDirectionDistributionString),                                                        ("Natal kick direction distribution (options: ISOTROPIC, INPLANE, PERPENDICULAR, POWERLAW, WEDGE, POLES), default = " + kickDirectionDistributionString + ")").c_str())
		    ("kick-velocity-distribution",                                  po::value<string>(&kickVelocityDistributionString)->default_value(kickVelocityDistributionString),                                                          ("Natal kick velocity distribution (options: ZERO, FIXED, FLAT, MAXWELLIAN, BRAYELDRIDGE, MULLER2016, MULLER2016MAXWELLIAN, MULLERMANDEL), default = " + kickVelocityDistributionString + ")").c_str())

            // JR - 01/04/2020 - Serena will uncomment when tested.
            // ("logfile-BSE-be-binaries",                                     po::value<string>(&logfileBSEBeBinaries)->default_value(logfileBSEBeBinaries),                                                                              ("Filename for BSE Be Binaries logfile (default = " + logfileBSEBeBinaries + ")").c_str())
            ("logfile-BSE-common-envelopes",                                po::value<string>(&logfileBSECommonEnvelopes)->default_value(logfileBSECommonEnvelopes),                                                                    ("Filename for BSE Common Envelopes logfile (default = " + logfileBSECommonEnvelopes + ")").c_str())
            ("logfile-BSE-detailed-output",                                 po::value<string>(&logfileBSEDetailedOutput)->default_value(logfileBSEDetailedOutput),                                                                      ("Filename for BSE Detailed Output logfile (default = " + logfileBSEDetailedOutput + ")").c_str())
            ("logfile-BSE-double-compact-objects",                          po::value<string>(&logfileBSEDoubleCompactObjects)->default_value(logfileBSEDoubleCompactObjects),                                                          ("Filename for BSE Double Compact Objects logfile (default = " + logfileBSEDoubleCompactObjects + ")").c_str())
            ("logfile-BSE-pulsar-evolution",                                po::value<string>(&logfileBSEPulsarEvolution)->default_value(logfileBSEPulsarEvolution),                                                                    ("Filename for BSE Pulsar Evolution logfile (default = " + logfileBSEPulsarEvolution + ")").c_str())
            ("logfile-BSE-supernovae",                                      po::value<string>(&logfileBSESupernovae)->default_value(logfileBSESupernovae),                                                                              ("Filename for BSE Supernovae logfile (default = " + logfileBSESupernovae + ")").c_str())
            ("logfile-BSE-system-parameters",                               po::value<string>(&logfileBSESystemParameters)->default_value(logfileBSESystemParameters),                                                                  ("Filename for BSE System Parameters logfile (default = " + logfileBSESystemParameters + ")").c_str())
            ("logfile-definitions",                                         po::value<string>(&logfileDefinitionsFilename)->default_value(logfileDefinitionsFilename)->implicit_value(""),                                              ("Filename for logfile record definitions (default = " + logfileDefinitionsFilename + ")").c_str())
            ("logfile-delimiter",                                           po::value<string>(&logfileDelimiterString)->default_value(logfileDelimiterString),                                                                          ("Field delimiter for logfile records (default = " + logfileDelimiterString + ")").c_str())
            ("logfile-name-prefix",                                         po::value<string>(&logfileNamePrefix)->default_value(logfileNamePrefix)->implicit_value(""),                                                                ("Prefix for logfile names (default = " + logfileNamePrefix + ")").c_str())
            ("logfile-SSE-parameters",                                      po::value<string>(&logfileSSEParameters)->default_value(logfileSSEParameters),                                                                              ("Filename for SSE Parameters logfile (default = " + logfileSSEParameters + ")").c_str())

		    ("mass-loss-prescription",                                      po::value<string>(&massLossPrescriptionString)->default_value(massLossPrescriptionString),                                                                  ("Mass loss prescription (options: NONE, HURLEY, VINK), default = " + massLossPrescriptionString + ")").c_str())
		    ("mass-ratio-distribution,q",                                   po::value<string>(&massRatioDistributionString)->default_value(massRatioDistributionString),                                                                ("Initial mass ratio distribution for q=m2/m1 (options: FLAT, DuquennoyMayor1991, SANA2012), default = " + massRatioDistributionString + ")").c_str())
		    ("mass-transfer-accretion-efficiency-prescription",             po::value<string>(&massTransferAccretionEfficiencyPrescriptionString)->default_value(massTransferAccretionEfficiencyPrescriptionString),                    ("Mass Transfer Accretion Efficiency prescription (options: THERMAL, FIXED), default = " + massTransferAccretionEfficiencyPrescriptionString + ")").c_str())
		    ("mass-transfer-angular-momentum-loss-prescription",            po::value<string>(&massTransferAngularMomentumLossPrescriptionString)->default_value(massTransferAngularMomentumLossPrescriptionString),                    ("Mass Transfer Angular Momentum Loss prescription (options: JEANS, ISOTROPIC, CIRCUMBINARY, ARBITRARY), default = " + massTransferAngularMomentumLossPrescriptionString + ")").c_str())
		    ("mass-transfer-rejuvenation-prescription",                     po::value<string>(&massTransferRejuvenationPrescriptionString)->default_value(massTransferRejuvenationPrescriptionString),                                  ("Mass Transfer Rejuvenation prescription (options: NONE, STARTRACK), default = " + massTransferRejuvenationPrescriptionString + ")").c_str())
			("mass-transfer-thermal-limit-accretor",                        po::value<string>(&massTransferThermallyLimitedVariationString)->default_value(massTransferThermallyLimitedVariationString),                                ("Mass Transfer Thermal Accretion limit (default = " + massTransferThermallyLimitedVariationString + ")").c_str())

		    ("neutrino-mass-loss-bh-formation",                             po::value<string>(&neutrinoMassLossAssumptionBHString)->default_value(neutrinoMassLossAssumptionBHString),                                                  ("Assumption about neutrino mass loss during BH formation (options: FIXED_FRACTION, FIXED_MASS), default = " + neutrinoMassLossAssumptionBHString + ")").c_str())

		    ("neutron-star-equation-of-state",                              po::value<string>(&neutronStarEquationOfStateString)->default_value(neutronStarEquationOfStateString),                                                      ("Neutron star equation of state to use (options: SSE, ARP3), default = " + neutronStarEquationOfStateString + ")").c_str())

   		    ("output-container,c",                                          po::value<string>(&outputContainerName)->default_value(outputContainerName)->implicit_value(""),                                                            ("Container (directory) name for output files (default = " + outputContainerName + ")").c_str())
   		    ("outputPath,o",                                                po::value<string>(&outputPathString)->default_value(outputPathString)->implicit_value(""),                                                                  ("Directory for output (default = " + outputPathString + ")").c_str())

		    ("pulsar-birth-magnetic-field-distribution",                    po::value<string>(&pulsarBirthMagneticFieldDistributionString)->default_value(pulsarBirthMagneticFieldDistributionString),                                  ("Pulsar Birth Magnetic Field distribution (options: ZERO, FIXED, FLATINLOG, UNIFORM, LOGNORMAL), default = " + pulsarBirthMagneticFieldDistributionString + ")").c_str())
		    ("pulsar-birth-spin-period-distribution",                       po::value<string>(&pulsarBirthSpinPeriodDistributionString)->default_value(pulsarBirthSpinPeriodDistributionString),                                        ("Pulsar Birth Spin Period distribution (options: ZERO, FIXED, UNIFORM, NORMAL), default = " + pulsarBirthSpinPeriodDistributionString + ")").c_str())
		    ("pulsational-pair-instability-prescription",                   po::value<string>(&pulsationalPairInstabilityPrescriptionString)->default_value(pulsationalPairInstabilityPrescriptionString),                              ("Pulsational Pair Instability prescription (options: COMPAS, STARTRACK, MARCHANT), default = " + pulsationalPairInstabilityPrescriptionString + ")").c_str())

		    ("remnant-mass-prescription",                                   po::value<string>(&remnantMassPrescriptionString)->default_value(remnantMassPrescriptionString),                                                            ("Choose remnant mass prescription (options: HURLEY2000, BELCZYNSKI2002, FRYER2012, MULLER2016, MULLERMANDEL), default = " + remnantMassPrescriptionString + ")").c_str())
		    ("rotational-velocity-distribution",                            po::value<string>(&rotationalVelocityDistributionString)->default_value(rotationalVelocityDistributionString),                                              ("Initial rotational velocity distribution (options: ZERO, HURLEY, VLTFLAMES), default = " + rotationalVelocityDistributionString + ")").c_str())

		    ("semi-major-axis-distribution,a",                              po::value<string>(&semiMajorAxisDistributionString)->default_value(semiMajorAxisDistributionString),                                                        ("Initial semi-major axis distribution, a (options: FLATINLOG, CUSTOM, DuquennoyMayor1991, SANA2012), default = " + semiMajorAxisDistributionString + ")").c_str())

            // vector (list) options - alphabetically
            ("debug-classes",                                               po::value<vector<string>>(&debugClasses)->multitoken()->default_value(debugClasses),                                                                        ("Debug classes enabled (default = " + defaultDebugClasses + ")").c_str())
            ("log-classes",                                                 po::value<vector<string>>(&logClasses)->multitoken()->default_value(logClasses),                                                                            ("Logging classes enabled (default = " + defaultLogClasses + ")").c_str())
		;

        po::variables_map vm;                                                                                                           // Variables map

        try {

            po::store(po::parse_command_line(argc, argv, desc), vm);

            // --help option
            if (vm["help"].as<bool>()) {                                                                                                // user requested help?
                utils::SplashScreen();                                                                                                  // yes - show splash screen
                ANNOUNCE(desc);                                                                                                         // and help
                programStatus = COMMANDLINE_STATUS::SUCCESS;                                                                            // ok
            }

            // --version option
            if (vm["version"].as<bool>()) {                                                                                             // user requested version?
                ANNOUNCE("COMPAS v" << VERSION_STRING);                                                                                 // yes, show version string
                programStatus = COMMANDLINE_STATUS::SUCCESS;                                                                            // ok
            }


            po::notify(vm);                                                                                                             // invoke notify to assign user-input values to variables.  Throws an error, so do after help just in case there are any problems.

            fixedRandomSeed  = !vm["random-seed"].defaulted();                                                                          // use random seed if it is provided by the user
            fixedMetallicity = !vm["metallicity"].defaulted();                                                                          // determine if user supplied a metallicity value
            useFixedUK       = !vm["fix-dimensionless-kick-velocity"].defaulted() && (fixedUK >= 0.0);                                  // determine if user supplied a valid kick velocity


            // check & set prescriptions, distributions, assumptions etc. options - alphabetically

            bool found;

            // AVG - 17/03/2020 - Floor will uncomment when tested.
            /*
            if (!vm["AIS-DCOtype"].defaulted()) {                                                                                       // Adaptive Importance Sampling DCO type
                std::tie(found, AISDCOtype) = utils::GetMapKey(AISDCOtypeString, AIS_DCO_LABEL, AISDCOtype);
                COMPLAIN_IF(!found, "Unknown AIS DCO Type");
            }
            */

            if (!vm["black-hole-kicks"].defaulted()) {                                                                                  // black hole kicks option
                std::tie(found, blackHoleKicksOption) = utils::GetMapKey(blackHoleKicksString, BLACK_HOLE_KICK_OPTION_LABEL, blackHoleKicksOption);
                COMPLAIN_IF(!found, "Unknown Black Hole Kicks Option");
            }

            if (!vm["chemically-homogeneous-evolution"].defaulted()) {                                                                  // Chemically Homogeneous Evolution
                std::tie(found, cheOption) = utils::GetMapKey(cheString, CHE_OPTION_LABEL, cheOption);
                COMPLAIN_IF(!found, "Unknown Chemically Homogeneous Evolution Option");
            }

            if (!vm["common-envelope-lambda-prescription"].defaulted()) {                                                               // common envelope lambda prescription
                std::tie(found, commonEnvelopeLambdaPrescription) = utils::GetMapKey(commonEnvelopeLambdaPrescriptionString, CE_LAMBDA_PRESCRIPTION_LABEL, commonEnvelopeLambdaPrescription);
                COMPLAIN_IF(!found, "Unknown CE Lambda Prescription");
            }

            if (!vm["common-envelope-mass-accretion-prescription"].defaulted()) {                                                       // common envelope mass accretion prescription
                std::tie(found, commonEnvelopeMassAccretionPrescription) = utils::GetMapKey(commonEnvelopeMassAccretionPrescriptionString, CE_ACCRETION_PRESCRIPTION_LABEL, commonEnvelopeMassAccretionPrescription);
                COMPLAIN_IF(!found, "Unknown CE Mass Accretion Prescription");
            }
            
            if (!vm["envelope-state-prescription"].defaulted()) {                                                       // envelope state prescription
                std::tie(found, envelopeStatePrescription) = utils::GetMapKey(envelopeStatePrescriptionString, ENVELOPE_STATE_PRESCRIPTION_LABEL, envelopeStatePrescription);
                COMPLAIN_IF(!found, "Unknown Envelope State Prescription");
            }

            if (!vm["stellar-zeta-prescription"].defaulted()) {                                                                 // common envelope zeta prescription
                std::tie(found, stellarZetaPrescription) = utils::GetMapKey(stellarZetaPrescriptionString, ZETA_PRESCRIPTION_LABEL, stellarZetaPrescription);
                COMPLAIN_IF(!found, "Unknown stellar Zeta Prescription");
            }

            if (!vm["eccentricity-distribution"].defaulted()) {                                                                         // eccentricity distribution
                std::tie(found, eccentricityDistribution) = utils::GetMapKey(eccentricityDistributionString, ECCENTRICITY_DISTRIBUTION_LABEL, eccentricityDistribution);
                COMPLAIN_IF(!found, "Unknown Eccentricity Distribution");
            }

            if (!vm["fryer-supernova-engine"].defaulted()) {                                                                            // Fryer et al. 2012 supernova engine
                std::tie(found, fryerSupernovaEngine) = utils::GetMapKey(fryerSupernovaEngineString, SN_ENGINE_LABEL, fryerSupernovaEngine);
                COMPLAIN_IF(!found, "Unknown Fryer et al. Supernova Engine");
            }

            if (!vm["initial-mass-function"].defaulted()) {                                                                             // initial mass function
                std::tie(found, initialMassFunction) = utils::GetMapKey(initialMassFunctionString, INITIAL_MASS_FUNCTION_LABEL, initialMassFunction);
                COMPLAIN_IF(!found, "Unknown Initial Mass Function");
            }

            if (!vm["kick-direction"].defaulted()) {                                                                                    // kick direction
                std::tie(found, kickDirectionDistribution) = utils::GetMapKey(kickDirectionDistributionString, KICK_DIRECTION_DISTRIBUTION_LABEL, kickDirectionDistribution);
                COMPLAIN_IF(!found, "Unknown Kick Direction Distribution");
            }

            if (!vm["kick-velocity-distribution"].defaulted()) {                                                                        // kick velocity
                std::tie(found, kickVelocityDistribution) = utils::GetMapKey(kickVelocityDistributionString, KICK_VELOCITY_DISTRIBUTION_LABEL, kickVelocityDistribution);
                COMPLAIN_IF(!found, "Unknown Kick Velocity Distribution");
            }

			if (!vm["logfile-delimiter"].defaulted()) {                                                                                 // logfile field delimiter
                std::tie(found, logfileDelimiter) = utils::GetMapKey(logfileDelimiterString, DELIMITERLabel, logfileDelimiter);
                COMPLAIN_IF(!found, "Unknown Logfile Delimiter");
            }

            if (!vm["mass-loss-prescription"].defaulted()) {                                                                            // mass loss prescription
                std::tie(found, massLossPrescription) = utils::GetMapKey(massLossPrescriptionString, MASS_LOSS_PRESCRIPTION_LABEL, massLossPrescription);
                COMPLAIN_IF(!found, "Unknown Mass Loss Prescription");
            }

            if (!vm["mass-ratio-distribution"].defaulted()) {                                                                           // mass ratio distribution
                std::tie(found, massRatioDistribution) = utils::GetMapKey(massRatioDistributionString, MASS_RATIO_DISTRIBUTION_LABEL, massRatioDistribution);
                COMPLAIN_IF(!found, "Unknown Mass Ratio Distribution");
            }


            if (useMassTransfer && !vm["mass-transfer-accretion-efficiency-prescription"].defaulted()) {                                // mass transfer accretion efficiency prescription
                std::tie(found, massTransferAccretionEfficiencyPrescription) = utils::GetMapKey(massTransferAccretionEfficiencyPrescriptionString, MT_ACCRETION_EFFICIENCY_PRESCRIPTION_LABEL, massTransferAccretionEfficiencyPrescription);
                COMPLAIN_IF(!found, "Unknown Mass Transfer Angular Momentum Loss Prescription");
            }

            if (useMassTransfer && !vm["mass-transfer-angular-momentum-loss-prescription"].defaulted()) {                               // mass transfer angular momentum loss prescription
                std::tie(found, massTransferAngularMomentumLossPrescription) = utils::GetMapKey(massTransferAngularMomentumLossPrescriptionString, MT_ANGULAR_MOMENTUM_LOSS_PRESCRIPTION_LABEL, massTransferAngularMomentumLossPrescription);
                COMPLAIN_IF(!found, "Unknown Mass Transfer Angular Momentum Loss Prescription");
            }

            if (useMassTransfer && !vm["mass-transfer-rejuvenation-prescription"].defaulted()) {                                        // mass transfer rejuvenation prescription
                std::tie(found, massTransferRejuvenationPrescription) = utils::GetMapKey(massTransferRejuvenationPrescriptionString, MT_REJUVENATION_PRESCRIPTION_LABEL, massTransferRejuvenationPrescription);
                COMPLAIN_IF(!found, "Unknown Mass Transfer Rejuvenation Prescription");
            }

            if (useMassTransfer && !vm["mass-transfer-thermal-limit-accretor"].defaulted()) {                                           // mass transfer accretor thermal limit
                std::tie(found, massTransferThermallyLimitedVariation) = utils::GetMapKey(massTransferThermallyLimitedVariationString, MT_THERMALLY_LIMITED_VARIATION_LABEL, massTransferThermallyLimitedVariation);
                COMPLAIN_IF(!found, "Unknown Mass Transfer Accretor Thermal Limit");

                if (found) {                                                                                                            // if user didn't specify choice of C factor, use default based on choice of thermally limited variation
                    if (massTransferThermallyLimitedVariation == MT_THERMALLY_LIMITED_VARIATION::C_FACTOR) {
                        massTransferCParameter = vm["mass-transfer-thermal-limit-C"].defaulted() ? 10.0 : massTransferCParameter;
                    }

                    if (massTransferThermallyLimitedVariation == MT_THERMALLY_LIMITED_VARIATION::RADIUS_TO_ROCHELOBE) {
                        massTransferCParameter = vm["mass-transfer-thermal-limit-C"].defaulted() ? 1.0 : massTransferCParameter;
                    }
                }
            }

            if (!vm["neutrino-mass-loss-bh-formation"].defaulted()) {                                                                   // neutrino mass loss assumption
                std::tie(found, neutrinoMassLossAssumptionBH) = utils::GetMapKey(neutrinoMassLossAssumptionBHString, NEUTRINO_MASS_LOSS_PRESCRIPTION_LABEL, neutrinoMassLossAssumptionBH);
                COMPLAIN_IF(!found, "Unknown Neutrino Mass Loss Assumption");
            }

            if (!vm["neutron-star-equation-of-state"].defaulted()) {                                                                    // neutron star equation of state
                std::tie(found, neutronStarEquationOfState) = utils::GetMapKey(neutronStarEquationOfStateString, NS_EOSLabel, neutronStarEquationOfState);
                COMPLAIN_IF(!found, "Unknown Neutron Star Equation of State");
            }

            if (!vm["pulsar-birth-magnetic-field-distribution"].defaulted()) {                                                          // pulsar birth magnetic field distribution
                std::tie(found, pulsarBirthMagneticFieldDistribution) = utils::GetMapKey(pulsarBirthMagneticFieldDistributionString, PULSAR_BIRTH_MAGNETIC_FIELD_DISTRIBUTION_LABEL, pulsarBirthMagneticFieldDistribution);
                COMPLAIN_IF(!found, "Unknown Pulsar Birth Magnetic Field Distribution");
            }

            if (!vm["pulsar-birth-spin-period-distribution"].defaulted()) {                                                             // pulsar birth spin period distribution
                std::tie(found, pulsarBirthSpinPeriodDistribution) = utils::GetMapKey(pulsarBirthSpinPeriodDistributionString, PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION_LABEL, pulsarBirthSpinPeriodDistribution);
                COMPLAIN_IF(!found, "Unknown Pulsar Birth Spin Period Distribution");
            }

			if (!vm["pulsational-pair-instability-prescription"].defaulted()) {                                                         // pulsational pair instability prescription
                std::tie(found, pulsationalPairInstabilityPrescription) = utils::GetMapKey(pulsationalPairInstabilityPrescriptionString, PPI_PRESCRIPTION_LABEL, pulsationalPairInstabilityPrescription);
                COMPLAIN_IF(!found, "Unknown Pulsational Pair Instability Prescription");
			}

            if (!vm["remnant-mass-prescription"].defaulted()) {                                                                         // remnant mass prescription
                std::tie(found, remnantMassPrescription) = utils::GetMapKey(remnantMassPrescriptionString, REMNANT_MASS_PRESCRIPTION_LABEL, remnantMassPrescription);
                COMPLAIN_IF(!found, "Unknown Remnant Mass Prescription");
            }

            if (!vm["rotational-velocity-distribution"].defaulted()) {                                                                  // rotational velocity distribution
                std::tie(found, rotationalVelocityDistribution) = utils::GetMapKey(rotationalVelocityDistributionString, ROTATIONAL_VELOCITY_DISTRIBUTION_LABEL, rotationalVelocityDistribution);
                COMPLAIN_IF(!found, "Unknown Rotational Velocity Distribution");
            }

            if (!vm["semi-major-axis-distribution"].defaulted()) {                                                                      // semi-major axis distribution
                std::tie(found, semiMajorAxisDistribution) = utils::GetMapKey(semiMajorAxisDistributionString, SEMI_MAJOR_AXIS_DISTRIBUTION_LABEL, semiMajorAxisDistribution);
                COMPLAIN_IF(!found, "Unknown Semi-Major Axis Distribution");
            }

            // constraint/value/range checks - alphabetically (where possible)

            COMPLAIN_IF(!vm["common-envelope-alpha"].defaulted() && commonEnvelopeAlpha < 0.0, "CE alpha (--common-envelope-alpha) < 0");
            COMPLAIN_IF(!vm["common-envelope-alpha-thermal"].defaulted() && (commonEnvelopeAlphaThermal < 0.0 || commonEnvelopeAlphaThermal > 1.0), "CE alpha thermal (--common-envelope-alpha-thermal) must be between 0 and 1");
            COMPLAIN_IF(!vm["common-envelope-lambda-multiplier"].defaulted() && commonEnvelopeLambdaMultiplier < 0.0, "CE lambda multiplie (--common-envelope-lambda-multiplier < 0");
            COMPLAIN_IF(!vm["common-envelope-mass-accretion-constant"].defaulted() && commonEnvelopeMassAccretionConstant < 0.0, "CE mass accretion constant (--common-envelope-mass-accretion-constant) < 0");
            COMPLAIN_IF(!vm["common-envelope-mass-accretion-max"].defaulted() && commonEnvelopeMassAccretionMax < 0.0, "Maximum accreted mass (--common-envelope-mass-accretion-max) < 0");
            COMPLAIN_IF(!vm["common-envelope-mass-accretion-min"].defaulted() && commonEnvelopeMassAccretionMin < 0.0, "Minimum accreted mass (--common-envelope-mass-accretion-min) < 0");

            COMPLAIN_IF(debugLevel < 0, "Debug level (--debug-level) < 0");

            COMPLAIN_IF(eccentricityDistributionMin < 0.0 || eccentricityDistributionMin > 1.0, "Minimum eccentricity (--eccentricity-min) must be between 0 and 1");
            COMPLAIN_IF(eccentricityDistributionMax < 0.0 || eccentricityDistributionMax > 1.0, "Maximum eccentricity (--eccentricity-max) must be between 0 and 1");
            COMPLAIN_IF(eccentricityDistributionMax <= eccentricityDistributionMin, "Maximum eccentricity (--eccentricity-max) must be > Minimum eccentricity (--eccentricity-min)");

            COMPLAIN_IF(initialMassFunctionMin < 0.0, "Minimum initial mass (--initial-mass-min) < 0");
            COMPLAIN_IF(initialMassFunctionMax < 0.0, "Maximum initial mass (--initial-mass-max) < 0");
            COMPLAIN_IF(initialMassFunctionMax <= initialMassFunctionMin, "Maximum initial mass (--initial-mass-max) must be > Minimum initial mass (--initial-mass-min)");

            if (kickVelocityDistribution == KICK_VELOCITY_DISTRIBUTION::FLAT) {
                COMPLAIN_IF(kickVelocityDistributionMaximum <= 0.0, "User specified --kick-velocity-distribution = FLAT with Maximum kick velocity (--kick-velocity-max) <= 0.0");
            }

            COMPLAIN_IF(logLevel < 0, "Logging level (--log-level) < 0");

            COMPLAIN_IF(massRatioDistributionMin < 0.0 || massRatioDistributionMin > 1.0, "Minimum mass ratio (--mass-ratio-min) must be between 0 and 1");
            COMPLAIN_IF(massRatioDistributionMax < 0.0 || massRatioDistributionMax > 1.0, "Maximum mass ratio (--mass-ratio-max) must be between 0 and 1");
            COMPLAIN_IF(massRatioDistributionMax <= massRatioDistributionMin, "Maximum mass ratio (--mass-ratio-max) must be > Minimum mass ratio (--mass-ratio-min)");

            COMPLAIN_IF(maxEvolutionTime <= 0.0, "Maximum evolution time in Myr (--maxEvolutionTime) must be > 0");

            COMPLAIN_IF(metallicity < 0.0 || metallicity > 1.0, "Metallicity (--metallicity) should be absolute metallicity and must be between 0 and 1");

            COMPLAIN_IF(minimumMassSecondary < 0.0, "Seconday minimum mass (--minimum-secondary-mass) must be >= 0");
            COMPLAIN_IF(minimumMassSecondary > initialMassFunctionMax, "Seconday minimum mass (--minimum-secondary-mass) must be <= Maximum initial mass (--initial-mass-max)");

            if (neutrinoMassLossAssumptionBH == NEUTRINO_MASS_LOSS_PRESCRIPTION::FIXED_MASS) {
               COMPLAIN_IF(neutrinoMassLossValueBH < 0.0, "Neutrino mass loss value < 0");
            }

            COMPLAIN_IF(nBinaries <= 0, "Number of binaries requested <= 0");

            if (neutrinoMassLossAssumptionBH == NEUTRINO_MASS_LOSS_PRESCRIPTION::FIXED_FRACTION) {
               COMPLAIN_IF(neutrinoMassLossValueBH < 0.0 || neutrinoMassLossValueBH > 1.0, "Neutrino mass loss must be between 0 and 1");
            }

            if (!vm["outputPath"].defaulted()) {                                                                                        // user specified output path?
                                                                                                                                        // yes
                fs::path userPath = outputPathString;                                                                                   // user-specifed path
                if (fs::is_directory(userPath)) {                                                                                       // valid directory?
                    outputPath = userPath;                                                                                              // yes - set outputPath to user-specified path
                }
                else {                                                                                                                  // not a valid directory
                    WARN("User-specified output path is not a valid directory. Using CWD.");                                            // show warning
                    outputPath = defaultOutputPath;                                                                                     // use default path = CWD
                }
            }

            COMPLAIN_IF(periodDistributionMin < 0.0, "Minimum orbital period (--orbital-period-min) < 0");
            COMPLAIN_IF(periodDistributionMax < 0.0, "Maximum orbital period (--orbital-period-max) < 0");

            COMPLAIN_IF(!vm["pulsar-magnetic-field-decay-timescale"].defaulted() && pulsarMagneticFieldDecayTimescale <= 0.0, "Pulsar magnetic field decay timescale (--pulsar-magnetic-field-decay-timescale) <= 0");
            COMPLAIN_IF(!vm["pulsar-magnetic-field-decay-massscale"].defaulted() && pulsarMagneticFieldDecayMassscale <= 0.0, "Pulsar Magnetic field decay massscale (--pulsar-magnetic-field-decay-massscale) <= 0");

            COMPLAIN_IF(semiMajorAxisDistributionMin < 0.0, "Minimum semi-major Axis (--semi-major-axis-min) < 0");
            COMPLAIN_IF(semiMajorAxisDistributionMax < 0.0, "Maximum semi-major Axis (--semi-major-axis-max) < 0");

            COMPLAIN_IF(singleStarMassMax   <= 0.0,               "Single star mass maximum (--single-star-mass-max) <= 0");
            COMPLAIN_IF(singleStarMassSteps > 1 && (singleStarMassMax <= singleStarMassMin), "Single star mass maximum (--single-star-mass-max) <= minimum (--single-star-mass-min)");
            COMPLAIN_IF(singleStarMassMin   <= 0.0,               "Single star mass minimum (--single-star-mass-min) <= 0");
            COMPLAIN_IF(singleStarMassSteps <= 0,                 "Single star mass steps (--single-star-mass-steps) <= 0");

            m_OptionsDetails = ProgramOptionDetails(vm);                                                                                  // construct options details string for output

        }
        catch (po::error& e) {                                                                                                          // program options exception
            std::cerr << "Program Options error: " << e.what() << std::endl;
            std::cerr << desc << std::endl;
            programStatus = COMMANDLINE_STATUS::ERROR_IN_COMMAND_LINE;                                                                  // set status
        }

        catch (...) {                                                                                                                   // unhandled problem - command line error
            std::cerr << "Error in CommandLine: unknown or invalid option" << std::endl;                                                // announce error
            programStatus = COMMANDLINE_STATUS::ERROR_IN_COMMAND_LINE;                                                                  // set status
        }

    } catch (std::exception& e) {                                                                                                       // unhandled exception - command line error
        std::cerr << "Unhandled exception: " << e.what() << std::endl;                                                                  // announce error
        programStatus = COMMANDLINE_STATUS::ERROR_UNHANDLED_EXCEPTION;                                                                  // set status
    }

    return programStatus;
}


/*
 * Determine the value of the requested program option
 *
 * The property is a boost variant variable, and is one of the following types:
 *
 *      STAR_PROPERTY           - any individual star property
 *      STAR_1_PROPERTY         - property of the primary (m_Star1)
 *      STAR_2_PROPERTY         - property of the secondary (m_Star2)
 *      SUPERNOVA_PROPERTY      - property of the star that has gone supernova
 *      COMPANION_PROPERTY      - property of the companion to the supernova
 *      BINARY_PROPERTY         - property of the binary
 *      PROGRAM_OPTION          - program option
 *
 * This function handles properties of type PROGRAM_OPTION
 *
 * This is the function used to retrieve values for properties required to be printed.
 * This allows the composition of the log records to be dynamically modified - this is
 * how we allow users to specify what properties they want recorded in log files.
 *
 * The functional return is the value of the property requested.  The type of the
 * functional return is a tuple: std::tuple<bool, COMPAS_VARIABLE_TYPE>.  This type
 * is COMPAS_VARIABLE by typedef.
 *
 * The bool returned indicates whether the property value was retrieved ok: true = yes, fales = no
 * The COMPAS_VARIABLE_TYPE variable returned is a boost variant variable, the value of which is the
 * value of the underlying primitive variable.
 *
 *
 * COMPAS_VARIABLE OptionValue(const T_ANY_PROPERTY p_Property) const
 *
 * @param   [IN]    p_Property                  The property for which the value is required
 * @return                                      The value of the requested property
 */
COMPAS_VARIABLE Options::OptionValue(const T_ANY_PROPERTY p_Property) const {

    bool ok = true;                                                                                                     // default is no error

    COMPAS_VARIABLE_TYPE value;                                                                                         // default property value

    PROGRAM_OPTION property = boost::get<PROGRAM_OPTION>(p_Property);                                                   // get property
                                                                                                                        // get property value
    switch (property) {

        case PROGRAM_OPTION::KICK_VELOCITY_DISTRIBUTION_SIGMA_CCSN_BH:  value = KickVelocityDistributionSigmaCCSN_BH(); break;
        case PROGRAM_OPTION::KICK_VELOCITY_DISTRIBUTION_SIGMA_CCSN_NS:  value = KickVelocityDistributionSigmaCCSN_NS(); break;
        case PROGRAM_OPTION::KICK_VELOCITY_DISTRIBUTION_SIGMA_FOR_ECSN: value = KickVelocityDistributionSigmaForECSN(); break;
        case PROGRAM_OPTION::KICK_VELOCITY_DISTRIBUTION_SIGMA_FOR_USSN: value = KickVelocityDistributionSigmaForUSSN(); break;
        case PROGRAM_OPTION::RANDOM_SEED:                               value = RandomSeed();                           break;

        default:                                                                                                        // unknown property
            ok    = false;                                                                                              // that's not ok...
            value = "UNKNOWN";                                                                                          // default value
            std::cerr << ERR_MSG(ERROR::UNKNOWN_PROGRAM_OPTION) << std::endl;                                           // show warning (don't have logging or errors here...)
    }

    return std::make_tuple(ok, value);
}


