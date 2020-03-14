#include "Options.h"

Options* Options::m_Instance = nullptr;

// this is required to set default value for boost program options of type vector<string>
namespace std
{
  std::ostream& operator<<(std::ostream &os, const std::vector<string> &vec) 
  {    
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


    namespace po = boost::program_options;
    namespace fs = boost::filesystem;


void PrintProgramOptionDetails(const boost::program_options::variables_map vm) {
    for (po::variables_map::const_iterator it = vm.begin(); it != vm.end(); it++) {
        std::cout << it->first;
        
        if (((boost::any)it->second.value()).empty()) {
            std::cout << "(empty)";
        }
        if (vm[it->first].defaulted() || it->second.defaulted()) {
            std::cout << "(default)";
        }


        bool is_char = false;
        try {
            boost::any_cast<const char *>(it->second.value());
            is_char = true;
            std::cout << "CHAR!!\n";
        } catch (const boost::bad_any_cast &) {
            is_char = false;
        }
        bool is_str = false;
        if (!is_char) {
        try {
            boost::any_cast<std::string>(it->second.value());
            is_str = true;
        } catch (const boost::bad_any_cast &) {
            is_str = false;
        }
        }

        if (((boost::any)it->second.value()).type() == typeid(unsigned)) {
            std::cout << "(unsigned)=";
            std::cout << vm[it->first].as<unsigned>() << std::endl;
        } else if (((boost::any)it->second.value()).type() == typeid(unsigned char)) {
            std::cout << "(unsigned char)=";
            std::cout << vm[it->first].as<unsigned char>() << std::endl;
        } else if (((boost::any)it->second.value()).type() == typeid(unsigned short)) {
            std::cout << "(unsigned short)=";
            std::cout << vm[it->first].as<unsigned short>() << std::endl;
        } else if (((boost::any)it->second.value()).type() == typeid(unsigned int)) {
            std::cout << "(unsigned int)=";
            std::cout << vm[it->first].as<unsigned int>() << std::endl;
        } else if (((boost::any)it->second.value()).type() == typeid(unsigned long)) {
            std::cout << "(unsigned long)=";
            std::cout << vm[it->first].as<unsigned long>() << std::endl;
        } else if (((boost::any)it->second.value()).type() == typeid(int)) {
            std::cout << "(int)=";
            std::cout << vm[it->first].as<int>() << std::endl;
        } else if (((boost::any)it->second.value()).type() == typeid(long int)) {
            std::cout << "(long int)=";
            std::cout << vm[it->first].as<long int>() << std::endl;
        } else if (((boost::any)it->second.value()).type() == typeid(unsigned long int)) {
            std::cout << "(unsigned long int)=";
            std::cout << vm[it->first].as<unsigned long int>() << std::endl;
        } else if (((boost::any)it->second.value()).type() == typeid(bool)) {
            std::cout << "(bool)=";
            bool v = vm[it->first].as<bool>();
            std::cout << (v ? ("true (" + std::to_string(v) + ")") : ("false (" + std::to_string(v) + ")")) << std::endl;
        } else if (((boost::any)it->second.value()).type() == typeid(double)) {
            std::cout << "(double)=";
            std::cout << vm[it->first].as<double>() << std::endl;
        } else if (is_char) {
            std::cout << "(char)=";
            std::cout << vm[it->first].as<const char * >() << std::endl;
        } else if (is_str) {
            std::cout << "(string)=";
            std::string temp = vm[it->first].as<std::string>();
            if (temp.size()) {
                std::cout << temp << std::endl;
            } else {
                std::cout << "<empty string>" << std::endl;
            }
        } else { // Assumes that the only remainder is vector<string>
            try {
                std::cout << "(vector<string>)= {";
                std::vector<std::string> vect = vm[it->first].as<std::vector<std::string> >();
                for (std::vector<std::string>::iterator oit=vect.begin(); oit != vect.end(); oit++) {
                    std::cout << (*oit) << ", ";
                }
                std::cout << "}" << std::endl;
            } catch (const boost::bad_any_cast &) {
                std::cout << "UnknownType(" << ((boost::any)it->second.value()).type().name() << ")" << std::endl;
            }
        }
    }
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
    onlyDoubleCompactObjects                                        = false;                                                                            // Flag to turn on some shortcuts to only evolve systems which may form double compact objects
    PNevolution                                                     = false;

    detailedOutput                                                  = false;                                                                            // Detailed output
    populationDataPrinting                                          = false;                                                                            // Print certain data for small populations, but not for larger one
    printBoolAsString                                               = false;                                                                            // default is do not print bool as string
    quiet                                                           = false;                                                                            // Suppress some of the printing

    nBatchesUsed                                                    = -1;                                                                               // Number of batches used, for STROOPWAFEL (AIS)


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
    usePairInstabilitySupernovae                                    = true;                                                                             // Whether to use pair instability supernovae (PISN)
    pairInstabilityUpperLimit                                       = 135.0;                                                                            // Maximum core mass leading to PISN (Value in Belczynski+ 2016 is 135 Msol)
    pairInstabilityLowerLimit                                       = 60.0;                                                                             // Minimum core mass leading to PISN (Value in Belczynski+ 2016 is 65 Msol)

    usePulsationalPairInstability                                   = true;                                                                             // Whether to use pulsational pair instability (PPI)
    pulsationalPairInstabilityLowerLimit                            = 35.0;                                                                             // Minimum core mass leading to PPI, Value in Belczynski+ 2016 is 45 Msol
    pulsationalPairInstabilityUpperLimit                            = 60.0;                                                                             // Maximum core mass leading to PPI, Value in Belczynski+ 2016 is 65 Msol

    pulsationalPairInstabilityPrescription                          = PPI_PRESCRIPTION::COMPAS;                                                         // Prescription for PPI to use
    pulsationalPairInstabilityPrescriptionString                    = PPI_PRESCRIPTION_LABEL.at(pulsationalPairInstabilityPrescription);                // String for which PPI prescription to use

	maximumNeutronStarMass                                          = 3.0;									                                            // Maximum mass of a neutron star allowed, value in StarTrack is 3.0


    // Kick direction option
    kickDirectionDistribution                                       = KICK_DIRECTION_DISTRIBUTION::ISOTROPIC;                                           // Which assumption for SN kicks: Possibilities: isotropic, inplane, perpendicular, powerlaw, wedge, poles
    kickDirectionDistributionString                                 = KICK_DIRECTION_DISTRIBUTION_LABEL.at(kickDirectionDistribution);		            // Which assumption for SN kicks: Possibilities: isotropic, inplane, perpendicular, powerlaw, wedge, poles
    kickDirectionPower                                              = 0.0;                                                                              // Power law power for the "power" SN kick direction choice

    // Get default output path
    outputPathString                                                = ".";                                                                               // String to hold the output directory
    defaultOutputPath                                               = boost::filesystem::current_path();                                                // Default output location
    outputPath                                                      = defaultOutputPath;                                                                // Desired output location

    // Spin options
    spinDistribution                                                = SPIN_DISTRIBUTION::FIXED;
    spinDistributionString                                          = SPIN_DISTRIBUTION_LABEL.at(spinDistribution);
    spinDistributionMin                                             = 0.60;
    spinDistributionMax                                             = 0.98;

    spinAssumption                                                  = SPIN_ASSUMPTION::ALIGNED;
    spinAssumptionString                                            = SPIN_ASSUMPTION_LABEL.at(spinAssumption);


    // Tides options
    tidesPrescription                                               = TIDES_PRESCRIPTION::NONE;                                                         // Tides prescription that will be used by the code
    tidesPrescriptionString                                         = TIDES_PRESCRIPTION_LABEL.at(tidesPrescription);                                   // String containing which tides prescription to use


    // Mass loss options
    useMassLoss                                                     = true;                                                                             // Whether to use mass loss

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

    massTransferPrescription                                        = MT_PRESCRIPTION::DEMINK;
	massTransferPrescriptionString                                  = MT_PRESCRIPTION_LABEL.at(massTransferPrescription);


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


    // Common Envelope parameters
    commonEnvelopePrescriptionFlag                                  = COMMON_ENVELOPE_PRESCRIPTION::WEBBINK;                                            // Which common envelope prescription to use
    commonEnvelopeAlpha                                             = 1.0;                                                                              // Common envelope efficiency alpha parameter
    commonEnvelopeLambda                                            = 0.1;                                                                              // Common envelope Lambda parameter
    commonEnvelopeHertzsprungGapDonor                               = COMMON_ENVELOPE_PRESCRIPTION::OPTIMISTIC_HG;                                      // Which prescription to use for Hertzsprung gap donors in a CE
    commonEnvelopeHertzsprungGapDonorString                         = COMMON_ENVELOPE_PRESCRIPTION_LABEL.at(commonEnvelopeHertzsprungGapDonor);         // String containing which prescription to use for Hertzsprung gap donors in a CE
	commonEnvelopeAlphaThermal                                      = 1.0;                                                                              // lambda = (alpha_th * lambda_b) + (1-alpha_th) * lambda_g
    commonEnvelopeLambdaMultiplier                                  = 1.0;                                                                              // Multiply common envelope lambda by some constant
    allowMainSequenceStarToSurviveCommonEnvelope                    = false;                                                                            // Whether or not to allow a main sequence star to survive a common envelope event

    // Accretion during common envelope
    commonEnvelopeMassAccretionPrescription                         = CE_ACCRETION_PRESCRIPTION::ZERO;
    commonEnvelopeMassAccretionPrescriptionString                   = CE_ACCRETION_PRESCRIPTION_LABEL.at(commonEnvelopeMassAccretionPrescription);

    commonEnvelopeMassAccretionMin                                  = 0.04;                                                                             // Minimum amount of mass accreted during CE in solar masses
    commonEnvelopeMassAccretionMax                                  = 0.1;                                                                              // Maximum amount of mass accreted during CE in solar masses
    commonEnvelopeMassAccretionConstant                             = 0.0;                                                                              // Constant value

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
	commonEnvelopeZetaPrescription                                  = CE_ZETA_PRESCRIPTION::SOBERMAN;					                                // Which prescription to use for calculating CE zetas
	commonEnvelopeZetaPrescriptionString                            = CE_ZETA_PRESCRIPTION_LABEL.at(commonEnvelopeZetaPrescription);				    // String containing which prescription to use for calculating CE zetas

	zetaAdiabaticArbitrary                                          = 0.0;
	zetaThermalArbitrary                                            = 0.0;
    zetaMainSequence 	                                            = 2.0;
	zetaHertzsprungGap	                                            = 6.5;


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
    logClasses = {"one", "two", "three"};

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


void Options::SetToFiducialValues(void) {

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
    onlyDoubleCompactObjects                                        = false;                                                                            // Flag to turn on some shortcuts to only evolve systems which may form double compact objects
    PNevolution                                                     = false;

    detailedOutput                                                  = false;                                                                            // Detailed output
    populationDataPrinting                                          = false;                                                                            // Print certain data for small populations, but not for larger one
    printBoolAsString                                               = false;                                                                            // default is do not print bool as string
    quiet                                                           = false;                                                                            // Suppress some of the printing

    nBatchesUsed                                                    = -1;                                                                               // Number of batches used, for STROOPWAFEL (AIS)


    // Public population synthesis variables
    nBinaries                                                       = 10;

    fixedRandomSeed                                                 = true;                                                                             // Whether to use a fixed random seed given by options.randomSeed (set to true if --random-seed is passed on command line)
    randomSeed                                                      = 0;                                                                                // Random seed to use


    // Specify how long to evolve binaries for
    maxEvolutionTime                                                = 13700.0;                                                                          // Maximum evolution time in Myrs
    maxNumberOfTimestepIterations                                   = 99999;                                                                            // Maximum number of timesteps to evolve binary for before giving up


    // Initial mass options
    initialMassFunction                                             = INITIAL_MASS_FUNCTION::KROUPA;
    initialMassFunctionString                                       = INITIAL_MASS_FUNCTION_LABEL.at(initialMassFunction);
    initialMassFunctionMin                                          = 8.;
    initialMassFunctionMax                                          = 100.;
    initialMassFunctionPower                                        = 0.;


    // Initial mass ratios
    massRatioDistribution                                           = MASS_RATIO_DISTRIBUTION::FLAT;                                                    // Most likely want Flat or SANA2012
    massRatioDistributionString                                     = MASS_RATIO_DISTRIBUTION_LABEL.at(massRatioDistribution);                          // Most likely want Flat or SANA2012
    massRatioDistributionMin                                        = 0.0;
    massRatioDistributionMax                                        = 1.0;

    minimumMassSecondary                                            = 0.1;                                                                              // Minimum mass of secondary to draw (in Msol) (brown dwarf limit = 0.1)


    // Initial orbit options
    semiMajorAxisDistribution                                       = SEMI_MAJOR_AXIS_DISTRIBUTION::FLATINLOG;
    semiMajorAxisDistributionString                                 = SEMI_MAJOR_AXIS_DISTRIBUTION_LABEL.at(SEMI_MAJOR_AXIS_DISTRIBUTION::FLATINLOG);   // Most likely want FlatInLog or SANA2012
    semiMajorAxisDistributionMin                                    = 0.1;
    semiMajorAxisDistributionMax                                    = 1000.0;
    semiMajorAxisDistributionPower                                  = -1.0;


    // Initial orbital period
    periodDistributionMin                                           = 1.1;                                                                              // Minimum initial period in days
    periodDistributionMax                                           = 1000.0;                                                                           // Maximum initial period in days


    // Eccentricity
    eccentricityDistribution                                        = ECCENTRICITY_DISTRIBUTION::ZERO;                                                  // Which eccentricity distribution to use
    eccentricityDistributionString                                  = ECCENTRICITY_DISTRIBUTION_LABEL.at(eccentricityDistribution);                     // Which eccentricity distribution to use
    eccentricityDistributionMin                                     = 0.0;                                                                              // Minimum initial eccentricity to sample
    eccentricityDistributionMax                                     = 1.0;                                                                              // Maximum initial eccentricity to sample


    // Kick options
    kickVelocityDistribution                                        = KICK_VELOCITY_DISTRIBUTION::MAXWELLIAN;		                                    // Which kick velocity distribution to use
    kickVelocityDistributionString                                  = KICK_VELOCITY_DISTRIBUTION_LABEL.at(kickVelocityDistribution);		            // Which kick velocity distribution to use


    // Kick velocity options
    kickVelocityDistributionSigmaCCSN_NS                            = 250;                                                                              // Kick velocity sigma in km s^-1 for neutron stars
    kickVelocityDistributionSigmaCCSN_BH                            = 250;                                                                              // Kick velocity sigma in km s^-1 for black holes
    kickVelocityDistributionMaximum                                 = -1.0;                                                                             // Maximum kick velocity to draw in km s^-1. Ignored if < 0
    kickVelocityDistributionSigmaForECSN   	                        = 30.0;                                                                             // Characteristic kick velocity for an ECSN in km s^-1
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
    usePairInstabilitySupernovae                                    = true;                                                                             // Whether to use pair instability supernovae (PISN)
    pairInstabilityUpperLimit                                       = 135.0;                                                                            // Maximum core mass leading to PISN (Value in Belczynski+ 2016 is 135 Msol)
    pairInstabilityLowerLimit                                       = 60.0;                                                                             // Minimum core mass leading to PISN (Value in Belczynski+ 2016 is 65 Msol)

    usePulsationalPairInstability                                   = true;                                                                             // Whether to use pulsational pair instability (PPI)
    pulsationalPairInstabilityLowerLimit                            = 35.0;                                                                             // Minimum core mass leading to PPI (Value in Belczynski+ 2016 is 45 Msol)
    pulsationalPairInstabilityUpperLimit                            = 60.0;                                                                             // Maximum core mass leading to PPI (Value in Belczynski+ 2016 is 65 Msol)

    pulsationalPairInstabilityPrescription                          = PPI_PRESCRIPTION::COMPAS;                                                         // Prescription for PPI to use
    pulsationalPairInstabilityPrescriptionString                    = PPI_PRESCRIPTION_LABEL.at(pulsationalPairInstabilityPrescription);                // String for which PPI prescription to use

	maximumNeutronStarMass                                          = 3.0;								                                                // Maximum mass of a neutron star allowed, valuse in StarTrack = 3.0


    // Kick direction option
    kickDirectionDistribution                                       = KICK_DIRECTION_DISTRIBUTION::ISOTROPIC;                                           // Which assumption for SN kicks: Possibilities: isotropic, inplane, perpendicular, powerlaw, wedge, poles
    kickDirectionDistributionString                                 = KICK_DIRECTION_DISTRIBUTION_LABEL.at(kickDirectionDistribution);		            // Which assumption for SN kicks: Possibilities: isotropic, inplane, perpendicular, powerlaw, wedge, poles
    kickDirectionPower                                              = 0.0;                                                                              // Power law power for the "power" SN kick direction choice

    // Get default output path
    outputPathString                                                = ".";                                                                              // String to hold the output directory
    defaultOutputPath                                               = boost::filesystem::current_path();                                                // Default output location
    outputPath                                                      = defaultOutputPath;                                                                // Desired output location

    // Spin options
    spinDistribution                                                = SPIN_DISTRIBUTION::ZERO;
    spinDistributionString                                          = SPIN_DISTRIBUTION_LABEL.at(spinDistribution);
    spinDistributionMin                                             = 0.0;
    spinDistributionMax                                             = 1.0;

    spinAssumption                                                  = SPIN_ASSUMPTION::ALIGNED;
    spinAssumptionString                                            = SPIN_ASSUMPTION_LABEL.at(spinAssumption);


    // Tides options
    tidesPrescription                                               = TIDES_PRESCRIPTION::NONE;                                                         // Tides prescription that will be used by the code
    tidesPrescriptionString                                         = TIDES_PRESCRIPTION_LABEL.at(tidesPrescription);                                   // String containing which tides prescription to use


    // Mass loss options
    useMassLoss                                                     = true;                                                                             // Whether to use mass loss

    massLossPrescription                                            = MASS_LOSS_PRESCRIPTION::VINK;
    massLossPrescriptionString                                      = MASS_LOSS_PRESCRIPTION_LABEL.at(massLossPrescription);


    // Wind mass loss multiplicitive constants
    luminousBlueVariableFactor                                      = 1.5;                                                                              // Luminous blue variable mass loss enhancement factor
    wolfRayetFactor                                                 = 1.0;                                                                              // WR winds factor


    // Mass transfer options
    useMassTransfer                                                 = true;											                                    // Whether to use mass transfer
	circulariseBinaryDuringMassTransfer	                            = false;						                                                    // Whether to circularise binary when it starts
	forceCaseBBBCStabilityFlag                                      = false;									                                        // Whether if all case BB/BC systems are forced to be stable or unstable
	alwaysStableCaseBBBCFlag                                        = false;							                                                // Whether if case BB/BC is always stable
	angularMomentumConservationDuringCircularisation                = false;		                                                                    // Whether to conserve angular momentum while circularising or circularise to periastron

    massTransferPrescription                                        = MT_PRESCRIPTION::DEMINK;
    massTransferPrescriptionString                                  = MT_PRESCRIPTION_LABEL.at(massTransferPrescription);


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
    massTransferRejuvenationPrescription                            = MT_REJUVENATION_PRESCRIPTION::STARTRACK;
    massTransferRejuvenationPrescriptionString                      = MT_REJUVENATION_PRESCRIPTION_LABEL.at(massTransferRejuvenationPrescription);


    // Mass transfer critical mass ratios
    massTransferCriticalMassRatioMSLowMass                          = false;                    			                                            // Whether to use critical mass ratios
    massTransferCriticalMassRatioMSLowMassNonDegenerateAccretor     = 1.44;                    			                                                // Critical mass ratio for MT from a MS low mass star (Claeys+ 2014 = 1.44)
    massTransferCriticalMassRatioMSLowMassDegenerateAccretor        = 1.0;                    			                                                // Critical mass ratio for MT from a MS low mass star on to a degenerate accretor (Claeys+ 2014 = 1.0)

    massTransferCriticalMassRatioMSHighMass                         = false;	           			                                                    // Whether to use critical mass ratios
    massTransferCriticalMassRatioMSHighMassNonDegenerateAccretor    = 0.625;                 			                                                // Critical mass ratio for MT from a MS high mass star (Claeys+ 2014 = 0.625)
    massTransferCriticalMassRatioMSHighMassDegenerateAccretor       = 0.0;                      			                                            // Critical mass ratio for MT from a MS high mass star on to a degenerate accretor

    massTransferCriticalMassRatioHG                                 = false;                 			                                                // Whether to use critical mass ratios
    massTransferCriticalMassRatioHGNonDegenerateAccretor            = 0.40;                  			                                                // Critical mass ratio for MT from a HG star (Claeys+ 2014 = 0.25)
    massTransferCriticalMassRatioHGDegenerateAccretor               = 0.21;                  			                                                // Critical mass ratio for MT from a HG star on to a degenerate accretor (Claeys+ 2014 = 0.21)

    massTransferCriticalMassRatioGiant                              = false;                  			                                                // Whether to use critical mass ratios
    massTransferCriticalMassRatioGiantNonDegenerateAccretor         = 0.0;                   			                                                // Critical mass ratio for MT from a giant
    massTransferCriticalMassRatioGiantDegenerateAccretor            = 0.87;                  			                                                // Critical mass ratio for MT from a giant on to a degenerate accretor (Claeys+ 2014 = 0.81)

    massTransferCriticalMassRatioHeliumMS                           = false;	           			                                                    // Whether to use critical mass ratios
    massTransferCriticalMassRatioHeliumMSNonDegenerateAccretor      = 0.625;                			                                                // Critical mass ratio for MT from a Helium MS star
    massTransferCriticalMassRatioHeliumMSDegenerateAccretor         = 0.0;                  			                                                // Critical mass ratio for MT from a Helium MS star on to a degenerate accretor

    massTransferCriticalMassRatioHeliumHG                           = false;                  			                                                // Whether to use critical mass ratios
    massTransferCriticalMassRatioHeliumHGNonDegenerateAccretor      = 0.25;		           			                                                    // Critical mass ratio for MT from a Helium HG star (Claeys+ 2014 = 0.25)
    massTransferCriticalMassRatioHeliumHGDegenerateAccretor         = 0.21;			           			                                                // Critical mass ratio for MT from a Helium HG star on to a degenerate accretor (Claeys+ 2014 = 0.21)

    massTransferCriticalMassRatioHeliumGiant                        = false;                 			                                                // Whether to use critical mass ratios
    massTransferCriticalMassRatioHeliumGiantNonDegenerateAccretor   = 1.28;                 			                                                // Critical mass ratio for MT from a Helium giant (Claeys+ 2014 = 0.25)
    massTransferCriticalMassRatioHeliumGiantDegenerateAccretor      = 0.87;                 			                                                // Critical mass ratio for MT from a Helium giant on to a degenerate accretor

    massTransferCriticalMassRatioWhiteDwarf                         = false;                			                                                // Whether to use critical mass ratios
	massTransferCriticalMassRatioWhiteDwarfNonDegenerateAccretor    = 0.0;                  			                                                // Critical mass ratio for MT from a White Dwarf
    massTransferCriticalMassRatioWhiteDwarfDegenerateAccretor       = 1.6;                     			                                                // Critical mass ratio for MT from a White Dwarf on to a degenerate accretor (Claeys+ 2014 = 1.6)


    // Common Envelope parameters
    commonEnvelopePrescriptionFlag                                  = COMMON_ENVELOPE_PRESCRIPTION::WEBBINK;                                            // Which common envelope prescription to use
    commonEnvelopeAlpha                                             = 1.0;                                                                              // Common envelope efficiency alpha parameter
    commonEnvelopeLambda                                            = 0.1;                                                                              // Common envelope Lambda parameter
    commonEnvelopeHertzsprungGapDonor                               = COMMON_ENVELOPE_PRESCRIPTION::PESSIMISTIC_HG;                                     // Which prescription to use for Hertzsprung gap donors in a CE
    commonEnvelopeHertzsprungGapDonorString                         = COMMON_ENVELOPE_PRESCRIPTION_LABEL.at(commonEnvelopeHertzsprungGapDonor);         // String containing which prescription to use for Hertzsprung gap donors in a CE
    commonEnvelopeAlphaThermal                                      = 1.0;                                                                              // lambda = (alpha_th * lambda_b) + (1-alpha_th) * lambda_g
    commonEnvelopeLambdaMultiplier                                  = 1.0;                                                                              // Multiply common envelope lambda by some constant
    allowMainSequenceStarToSurviveCommonEnvelope                    = false;                                                                            // Whether or not to allow a main sequence star to survive a common envelope event


    // Accretion during common envelope
    commonEnvelopeMassAccretionPrescription                         = CE_ACCRETION_PRESCRIPTION::ZERO;
    commonEnvelopeMassAccretionPrescriptionString                   = CE_ACCRETION_PRESCRIPTION_LABEL.at(commonEnvelopeMassAccretionPrescription);

    commonEnvelopeMassAccretionMin                                  = 0.04;                                                                             // Minimum amount of mass accreted during CE in solar masses
    commonEnvelopeMassAccretionMax                                  = 0.1;                                                                              // Maximum amount of mass accreted during CE in solar masses
    commonEnvelopeMassAccretionConstant                             = 0.0;                                                                              // Constant value


	// Common envelope lambda prescription
	commonEnvelopeLambdaPrescription                                = CE_LAMBDA_PRESCRIPTION::NANJING;                                                  // Which prescription to use for CE lambda
	commonEnvelopeLambdaPrescriptionString                          = CE_LAMBDA_PRESCRIPTION_LABEL.at(commonEnvelopeLambdaPrescription);                // String containing which prescription to use for CE lambda


	// Common envelope Nandez and Ivanova energy formalism
	revisedEnergyFormalismNandezIvanova	                            = false;						                                                    // Use the revised energy formalism from Nandez & Ivanova 2016
	maximumMassDonorNandezIvanova                                   = 2.0;								                                                // Maximum mass allowed to use the revised energy formalism in Msol
	commonEnvelopeRecombinationEnergyDensity                        = 1.5E13;					                                                        // Factor using to calculate the bin


	// Common envelope power factor for Kruckow fit
	commonEnvelopeSlopeKruckow                                      = -2.0/3.0;								                                            // Common envelope power factor for Kruckow fit normalized according to Kruckow+2016, Fig. 1


	// Which prescription to use for calculating zetas
	commonEnvelopeZetaPrescription                                  = CE_ZETA_PRESCRIPTION::SOBERMAN;					                                // Which prescription to use for calculating CE zetas
	commonEnvelopeZetaPrescriptionString                            = CE_ZETA_PRESCRIPTION_LABEL.at(commonEnvelopeZetaPrescription);					// String containing which prescription to use for calculating CE zetas


	zetaAdiabaticArbitrary                                          = 0.0;
	zetaThermalArbitrary                                            = 0.0;
    zetaMainSequence 	                                            = 6.5;
	zetaHertzsprungGap	                                            = 2.0;


    // Adaptive Importance Sampling Exploratory phase
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
    neutronStarEquationOfStateString                                = NS_EOSLabel.at(neutronStarEquationOfState);


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

            ("AIS-exploratory-phase",                                       po::value<bool>(&AISexploratoryPhase)->default_value(AISexploratoryPhase)->implicit_value(true),                                                            ("Run exploratory phase of STROOPWAFEL (default = " + std::string(AISexploratoryPhase ? "TRUE" : "FALSE") + ")").c_str())
		    ("AIS-Hubble",                                                  po::value<bool>(&AIShubble)->default_value(AIShubble)->implicit_value(true),                                                                                ("Excluding not in Hubble time mergers selection in exploratory phase of STROOPWAFEL (default = " + std::string(AIShubble ? "TRUE" : "FALSE") + ")").c_str())
		    ("AIS-Pessimistic",                                             po::value<bool>(&AISpessimistic)->default_value(AISpessimistic)->implicit_value(true),                                                                      ("Optimistic or Pessimistic selection in exploratory phase of STROOPWAFEL (default = " + std::string(AISpessimistic ? "TRUE" : "FALSE") + ")").c_str())
		    ("AIS-refinement-phase",                                        po::value<bool>(&AISrefinementPhase)->default_value(AISrefinementPhase)->implicit_value(true),                                                              ("Run main sampling phase (step2) of STROOPWAFEL (default = " + std::string(AISrefinementPhase ? "TRUE" : "FALSE") + ")").c_str())
		    ("AIS-RLOF",                                                    po::value<bool>(&AISrlof)->default_value(AISrlof)->implicit_value(true),                                                                                    ("RLOFSecondaryZAMS selection in exploratory phase of STROOPWAFEL (default = " + std::string(AISrlof ? "TRUE" : "FALSE") + ")").c_str())

		    ("allow-rlof-at-birth",                                         po::value<bool>(&allowRLOFAtBirth)->default_value(allowRLOFAtBirth)->implicit_value(true),                                                                  ("Allow binaries that have one or both stars in RLOF at birth to evolve (default = " + std::string(allowRLOFAtBirth ? "TRUE" : "FALSE") + ")").c_str())
		    ("allow-touching-at-birth",                                     po::value<bool>(&allowTouchingAtBirth)->default_value(allowTouchingAtBirth)->implicit_value(true),                                                          ("Allow binaries that are touching at birth to evolve (default = " + std::string(allowTouchingAtBirth ? "TRUE" : "FALSE") + ")").c_str())

			("alwaysStableCaseBBBCFlag",                                    po::value<bool>(&alwaysStableCaseBBBCFlag)->default_value(alwaysStableCaseBBBCFlag)->implicit_value(true),                                                  ("Choose case BB/BC mass transfer to be always stable (default = " + std::string(alwaysStableCaseBBBCFlag ? "TRUE" : "FALSE") + ")").c_str())
			("angularMomentumConservationDuringCircularisation",            po::value<bool>(&angularMomentumConservationDuringCircularisation)->default_value(angularMomentumConservationDuringCircularisation)->implicit_value(true),  ("Conserve angular momentum when binary is circularised when entering a Mass Transfer episode (default = " + std::string(angularMomentumConservationDuringCircularisation ? "TRUE" : "FALSE") + ")").c_str())
			("BeBinaries",                                                  po::value<bool>(&beBinaries)->default_value(beBinaries)->implicit_value(true),                                                                              ("Enable Be Binaries study (default = " + std::string(beBinaries ? "TRUE" : "FALSE") + ")").c_str())
			("circulariseBinaryDuringMassTransfer",                         po::value<bool>(&circulariseBinaryDuringMassTransfer)->default_value(circulariseBinaryDuringMassTransfer)->implicit_value(true),                            ("Circularise binary when it enters a Mass Transfer episode (default = " + std::string(circulariseBinaryDuringMassTransfer ? "TRUE" : "FALSE") + ")").c_str())
		    ("common-envelope-allow-main-sequence-survive",                 po::value<bool>(&allowMainSequenceStarToSurviveCommonEnvelope)->default_value(allowMainSequenceStarToSurviveCommonEnvelope)->implicit_value(true),          ("Allow main sequence stars to survive common envelope evolution (default = " + std::string(allowMainSequenceStarToSurviveCommonEnvelope ? "TRUE" : "FALSE") + ")").c_str())

			("debug-to-file",                                               po::value<bool>(&debugToFile)->default_value(debugToFile)->implicit_value(debugToFile), "Write debug statements to file")
		    ("detailedOutput",                                              po::value<bool>(&detailedOutput)->default_value(detailedOutput)->implicit_value(detailedOutput), "Print detailed output to file")
			("errors-to-file",                                              po::value<bool>(&errorsToFile)->default_value(errorsToFile)->implicit_value(errorsToFile), "Write error messages to file")

		    ("evolve-pulsars",                                              po::value<bool>(&evolvePulsars)->default_value(evolvePulsars)->implicit_value(evolvePulsars), "Whether to evolve pulsars")
			("evolve-unbound-systems",                                      po::value<bool>(&evolveUnboundSystems)->default_value(evolveUnboundSystems)->implicit_value(evolveUnboundSystems), "Keeps evolving stars even if the binary is disrupted")

            ("forceCaseBBBCStabilityFlag",                                  po::value<bool>(&forceCaseBBBCStabilityFlag)->default_value(forceCaseBBBCStabilityFlag)->implicit_value(true), "Force case BB/BC mass transfer to be only stable or unstable (default = True)")
			("lambda-calculation-every-timeStep",                           po::value<bool>(&lambdaCalculationEveryTimeStep)->default_value(lambdaCalculationEveryTimeStep)->implicit_value(lambdaCalculationEveryTimeStep), "Calculate all values of lambda at each timestep")
   		   	("massTransfer",                                                po::value<bool>()->default_value(false)->implicit_value(), "Enable mass transfer")
		    ("only-double-compact-objects",                                 po::value<bool>(&onlyDoubleCompactObjects)->default_value(onlyDoubleCompactObjects)->implicit_value(onlyDoubleCompactObjects), "Only evolve binaries which may form double compact objects")
		    ("pair-instability-supernovae",                                 po::value<bool>(&usePairInstabilitySupernovae)->default_value(usePairInstabilitySupernovae)->implicit_value(usePairInstabilitySupernovae), "Enable pair instability supernovae (PISN)")
			("PN",                                                          po::value<bool>(&PNevolution)->default_value(PNevolution)->implicit_value(PNevolution), "Enable post-newtonian evolution")
            ("populationDataPrinting",                                      po::value<bool>(&populationDataPrinting)->default_value(populationDataPrinting)->implicit_value(populationDataPrinting), "Print details of population")
		    ("print-bool-as-string",                                        po::value<bool>(&printBoolAsString)->default_value(printBoolAsString)->implicit_value(printBoolAsString), "Print boolean properties as 'TRUE' or 'FALSE'")
		    ("pulsational-pair-instability",                                po::value<bool>(&usePulsationalPairInstability)->default_value(usePulsationalPairInstability)->implicit_value(usePulsationalPairInstability), "Enable mass loss due to pulsational-pair-instability (PPI)")
		    ("quiet",                                                       po::value<bool>(&quiet)->default_value(quiet)->implicit_value(quiet), "Suppress printing")
			("revised-energy-formalism-Nandez-Ivanova",                     po::value<bool>(&revisedEnergyFormalismNandezIvanova)->default_value(revisedEnergyFormalismNandezIvanova)->implicit_value(revisedEnergyFormalismNandezIvanova), "Enable revised energy formalism")

			("sample-common-envelope-alpha",                                po::value<bool>(&sampleCommonEnvelopeAlpha)->default_value(sampleCommonEnvelopeAlpha)->implicit_value(sampleCommonEnvelopeAlpha), "Sample over common envelope alpha")
			("sample-kick-direction-power",                                 po::value<bool>(&sampleKickDirectionPower)->default_value(sampleKickDirectionPower)->implicit_value(sampleKickDirectionPower), "Sample over kick direction powerlaw exponent")
			("sample-kick-velocity-sigma",                                  po::value<bool>(&sampleKickVelocitySigma)->default_value(sampleKickVelocitySigma)->implicit_value(sampleKickVelocitySigma), "Sample over Kick Velocity Sigma")
			("sample-luminous-blue-variable-multiplier",                    po::value<bool>(&sampleLuminousBlueVariableMultiplier)->default_value(sampleLuminousBlueVariableMultiplier)->implicit_value(sampleLuminousBlueVariableMultiplier), "Sample over multiplicative constant from LBV mass loss")
			("sample-wolf-rayet-multiplier",                                po::value<bool>(&sampleWolfRayetMultiplier)->default_value(falsampleWolfRayetMultiplierse)->implicit_value(sampleWolfRayetMultiplier), "Sample over WR winds multiplicative constant")

            ("single-star",                                                 po::value<bool>(&singleStar)->default_value(singleStar)->implicit_value(singleStar), "Evolve single star(s)")

		    ("use-mass-loss",                                               po::value<bool>(&useMassLoss)->default_value(useMassLoss)->implicit_value(useMassLoss), "Enable mass loss")

			("zeta-calculation-every-timestep",                             po::value<bool>(&zetaCalculationEveryTimeStep)->default_value(zetaCalculationEveryTimeStep)->implicit_value(zetaCalculationEveryTimeStep), "Calculate all values of MT zetas at each timestep")


			// numerical options - alphabetically by groups

			// unsigned long
		    ("random-seed",                                                 po::value<unsigned long>()->default_value(randomSeed),                              ("Random seed to use (default = " + std::to_string(randomSeed) + ")").c_str())

		    // int
			("debug-level",                                                 po::value<int>()->default_value(debugLevel),                                        ("Determines which print statements are displayed for debugging (default = " + std::to_string(debugLevel) + ")").c_str())
		    ("log-level",                                                   po::value<int>()->default_value(logLevel),                                          ("Determines which print statements are included in the logfile (default = " + std::to_string(logLevel) + ")").c_str())
		    ("maximum-number-iterations",                                   po::value<int>()->default_value(maxNumberOfTimestepIterations),                     ("Maximum number of timesteps to evolve binary before giving up (default = " + std::to_string(maxNumberOfTimestepIterations) + ")").c_str())
			("number-of-binaries,n",                                        po::value<int>()->default_value(nBinaries),                                         ("Specify the number of binaries to simulate (default = " + std::to_string(nBinaries) + ")").c_str())
			("single-star-mass-steps",                                      po::value<int>()->default_value(singleStarMassSteps),                               ("Specify the number of mass steps for single star evolution (default = " + std::to_string(singleStarMassSteps) + ")").c_str())

		    // double
		    ("common-envelope-alpha",                                       po::value<double>()->default_value(commonEnvelopeAlpha),                            ("Common Envelope efficiency alpha (default = " + std::to_string(commonEnvelopeAlpha) + ")").c_str())
		    ("common-envelope-alpha-thermal",                               po::value<double>()->default_value(commonEnvelopeAlphaThermal),                     ("Defined such that lambda = alpha_th * lambda_b + (1.0 - alpha_th) * lambda_g (default = " + std::to_string(commonEnvelopeAlphaThermal) + ")").c_str())
		    ("common-envelope-lambda",                                      po::value<double>()->default_value(commonEnvelopeLambda),                           ("Common Envelope lambda (default = " + std::to_string(commonEnvelopeLambda) + ")").c_str())
		    ("common-envelope-lambda-multiplier",                           po::value<double>()->default_value(commonEnvelopeLambdaMultiplier),                 ("Multiply lambda by some constant (default = " + std::to_string(commonEnvelopeLambdaMultiplier) + ")").c_str())
	    	("common-envelope-mass-accretion-constant",                     po::value<double>()->default_value(commonEnvelopeMassAccretionConstant),            ("Value of mass accreted by NS/BH during common envelope evolution if assuming all NS/BH accrete same amount of mass (common-envelope-mass-accretion-prescription CONSTANT). Ignored otherwise (default = " + std::to_string(commonEnvelopeMassAccretionConstant) + ")").c_str())
		    ("common-envelope-mass-accretion-max",                          po::value<double>()->default_value(commonEnvelopeMassAccretionMax),                 ("Maximum amount of mass accreted by NS/BHs during common envelope evolution in solar masses (default = " + std::to_string(commonEnvelopeMassAccretionMax) + ")").c_str())
		    ("common-envelope-mass-accretion-min",                          po::value<double>()->default_value(commonEnvelopeMassAccretionMin),                 ("Minimum amount of mass accreted by NS/BHs during common envelope evolution in solar masses (default = " + std::to_string(commonEnvelopeMassAccretionMin) + ")").c_str())
		    ("common-envelope-recombination-energy-density",                po::value<double>()->default_value(commonEnvelopeRecombinationEnergyDensity),       ("Recombination energy density in ergs/g (default = " + std::to_string(commonEnvelopeRecombinationEnergyDensity) + ")").c_str())
		    ("common-envelope-slope-Kruckow",                               po::value<double>()->default_value(commonEnvelopeSlopeKruckow),                     ("Common Envelope slope for Kruckow lambda (default = " + std::to_string(commonEnvelopeSlopeKruckow) + ")").c_str())

            ("critical-mass-ratio-giant-degenerate-accretor",               po::value<double>()->default_value(massTransferCriticalMassRatioGiantDegenerateAccretor),           ("Critical mass ratio for MT from a giant star (default = " + std::to_string(massTransferCriticalMassRatioGiantDegenerateAccretor) + " from Claeys+ 2014) Specify both giant flags to use. 0 is always stable, <0 is disabled").c_str())
            ("critical-mass-ratio-giant-non-degenerate-accretor",           po::value<double>()->default_value(massTransferCriticalMassRatioGiantNonDegenerateAccretor),        ("Critical mass ratio for MT from a giant star (default = " + std::to_string(massTransferCriticalMassRatioGiantNonDegenerateAccretor) + " from Claeys+ 2014) Specify both giant flags to use. 0 is always stable, <0 is disabled").c_str())
            ("critical-mass-ratio-helium-giant-degenerate-accretor",        po::value<double>()->default_value(massTransferCriticalMassRatioHeliumGiantDegenerateAccretor),     ("Critical mass ratio for MT from a helium giant star (default = " + std::to_string(massTransferCriticalMassRatioHeliumGiantDegenerateAccretor) + " from Claeys+ 2014) Specify both helium giant flags to use. 0 is always stable, <0 is disabled").c_str())
            ("critical-mass-ratio-helium-giant-non-degenerate-accretor",    po::value<double>()->default_value(massTransferCriticalMassRatioHeliumGiantNonDegenerateAccretor),  ("Critical mass ratio for MT from a helium giant star (default = " + std::to_string(massTransferCriticalMassRatioHeliumGiantNonDegenerateAccretor) + " from Claeys+ 2014) Specify both helium giant flags to use. 0 is always stable, <0 is disabled").c_str())
            ("critical-mass-ratio-helium-HG-degenerate-accretor",           po::value<double>()->default_value(massTransferCriticalMassRatioHeliumHGDegenerateAccretor),        ("Critical mass ratio for MT from a helium HG star (default = " + std::to_string(massTransferCriticalMassRatioHeliumHGDegenerateAccretor) + " from Claeys+ 2014) Specify both helium HG flags to use. 0 is always stable, <0 is disabled").c_str())
            ("critical-mass-ratio-helium-HG-non-degenerate-accretor",       po::value<double>()->default_value(massTransferCriticalMassRatioHeliumHGNonDegenerateAccretor),     ("Critical mass ratio for MT from a helium HG star (default = " + std::to_string(massTransferCriticalMassRatioHeliumHGNonDegenerateAccretor) + " from Claeys+ 2014) Specify both helium HG flags to use. 0 is always stable, <0 is disabled").c_str())
            ("critical-mass-ratio-helium-MS-degenerate-accretor",           po::value<double>()->default_value(massTransferCriticalMassRatioHeliumMSDegenerateAccretor),        ("Critical mass ratio for MT from a helium MS star (default = " + std::to_string(massTransferCriticalMassRatioHeliumMSDegenerateAccretor) + ") Specify both helium MS flags to use. 0 is always stable, <0 is disabled").c_str())
            ("critical-mass-ratio-helium-MS-non-degenerate-accretor",       po::value<double>()->default_value(massTransferCriticalMassRatioHeliumMSNonDegenerateAccretor),     ("Critical mass ratio for MT from a helium MS star (default = " + std::to_string(massTransferCriticalMassRatioHeliumMSNonDegenerateAccretor) + ") Specify both helium MS flags to use. 0 is always stable, <0 is disabled").c_str())
            ("critical-mass-ratio-HG-degenerate-accretor",                  po::value<double>()->default_value(massTransferCriticalMassRatioHGDegenerateAccretor),              ("Critical mass ratio for MT from a HG star (default = " + std::to_string(massTransferCriticalMassRatioHGDegenerateAccretor) + " from Claeys+ 2014) Specify both HG flags to use. 0 is always stable, <0 is disabled").c_str())
            ("critical-mass-ratio-HG-non-degenerate-accretor",              po::value<double>()->default_value(massTransferCriticalMassRatioHGNonDegenerateAccretor),           ("Critical mass ratio for MT from a HG star (default = " + std::to_string(massTransferCriticalMassRatioHGNonDegenerateAccretor) + " from de Mink+ 2013) Specify both HG flags to use. 0 is always stable, <0 is disabled").c_str())
            ("critical-mass-ratio-MS-high-mass-degenerate-accretor",        po::value<double>()->default_value(massTransferCriticalMassRatioMSHighMassDegenerateAccretor),      ("Critical mass ratio for MT from a MS star to a degenerate accretor (default = " + std::to_string(massTransferCriticalMassRatioMSHighMassDegenerateAccretor) + " from Claeys+ 2014) Specify both MS high mass flags to use. 0 is always stable, <0 is disabled").c_str())
            ("critical-mass-ratio-MS-high-mass-non-degenerate-accretor",    po::value<double>()->default_value(massTransferCriticalMassRatioMSHighMassNonDegenerateAccretor),   ("Critical mass ratio for MT from a MS star (default = " + std::to_string(massTransferCriticalMassRatioMSHighMassNonDegenerateAccretor) + ", Claeys+ 2014). Specify both MS high mass flags to use. 0 is always stable, <0 is disabled").c_str())
            ("critical-mass-ratio-MS-low-mass-degenerate-accretor",         po::value<double>()->default_value(massTransferCriticalMassRatioMSLowMassDegenerateAccretor),       ("Critical mass ratio for MT from a MS star to a degenerate accretor (default = " + std::to_string(massTransferCriticalMassRatioMSLowMassDegenerateAccretor) + " from Claeys+ 2014) Specify both MS low mass flags to use. 0 is always stable, <0 is disabled").c_str())
            ("critical-mass-ratio-MS-low-mass-non-degenerate-accretor",     po::value<double>()->default_value(massTransferCriticalMassRatioMSLowMassNonDegenerateAccretor),    ("Critical mass ratio for MT from a MS star (default = " + std::to_string(massTransferCriticalMassRatioMSLowMassNonDegenerateAccretor) + ", Claeys+ 2014). Specify both MS low mass flags to use. 0 is always stable, <0 is disabled").c_str())
            ("critical-mass-ratio-white-dwarf-degenerate-accretor",         po::value<double>()->default_value(massTransferCriticalMassRatioWhiteDwarfDegenerateAccretor),      ("Critical mass ratio for MT from a white dwarf (default = " + std::to_string(massTransferCriticalMassRatioWhiteDwarfDegenerateAccretor) + " from Claeys+ 2014) Specify both white dwarf flags to use. 0 is always stable, <0 is disabled").c_str())
            ("critical-mass-ratio-white-dwarf-non-degenerate-accretor",     po::value<double>()->default_value(massTransferCriticalMassRatioWhiteDwarfNonDegenerateAccretor),   ("Critical mass ratio for MT from a white dwarf (default = " + std::to_string(massTransferCriticalMassRatioWhiteDwarfNonDegenerateAccretor) + ") Specify both white dwarf flags to use. 0 is always stable, <0 is disabled").c_str())

		    ("eccentricity-max",                                            po::value<double>()->default_value(eccentricityDistributionMax),                                    ("Maximum eccentricity to generate (default = " + std::to_string(eccentricityDistributionMax) + ")").c_str())
		    ("eccentricity-min",                                            po::value<double>()->default_value(eccentricityDistributionMin),                                    ("Minimum eccentricity to generate (default = " + std::to_string(eccentricityDistributionMin) + ")").c_str())
			("eddington-accretion-factor",                                  po::value<double>()->default_value(eddingtonAccretionFactor),                                       ("Multiplication factor for eddington accretion for NS & BH, i.e. >1 is super-eddington and 0. is no accretion (default = " + std::to_string(eddingtonAccretionFactor) + ")").c_str())

   		    ("fix-dimensionless-kick-velocity",                             po::value<double>()->default_value(fixedUK),                                                        ("Fix dimensionless kick velocity uk to this value (default = " + std::to_string(fixedUK) + ", -ve values false, +ve values true)").c_str())

		    ("initial-mass-max",                                            po::value<double>()->default_value(initialMassFunctionMax),                                         ("Maximum mass (in Msol) to generate using given IMF (default = " + std::to_string(initialMassFunctionMax) + ")").c_str())
		    ("initial-mass-min",                                            po::value<double>()->default_value(initialMassFunctionMin),                                         ("Minimum mass (in Msol) to generate using given IMF (default = " + std::to_string(initialMassFunctionMin) + ")").c_str())
		    ("initial-mass-power",                                          po::value<double>()->default_value(initialMassFunctionPower),                                       ("Single power law power to generate primary mass using given IMF (default = " + std::to_string(initialMassFunctionPower) + ")").c_str())

		    ("kappa-gaussians",                                             po::value<double>()->default_value(kappaGaussians),                                                 ("Scaling factor for the width of the Gaussian distributions in STROOPWAFEL main sampling phase (default = " + std::to_string(kappaGaussians) + ")").c_str())
		    ("kick-direction-power",                                        po::value<double>()->default_value(kickDirectionPower),                                             ("Power for power law kick direction distribution (default = " + std::to_string(kickDirectionPower) + " = isotropic, +ve = polar, -ve = in plane)").c_str())
			("kick-scaling-factor",                                         po::value<double>()->default_value(kickScalingFactor),                                              ("Arbitrary factor used to scale kicks (default = " + std::to_string(kickScalingFactor) + ")").c_str())
		    ("kick-velocity-max",                                           po::value<double>()->default_value(kickVelocityDistributionMaximum),                                ("Maximum drawn kick velocity in km s^-1. Ignored if < 0. Must be > 0 if using kick-velocity-distribution=FLAT (default = " + std::to_string(kickVelocityDistributionMaximum) + ")").c_str())
		    ("kick-velocity-sigma-CCSN-BH",                                 po::value<double>()->default_value(kickVelocityDistributionSigmaCCSN_BH),                           ("Sigma for chosen kick velocity distribution for black holes (default = " + std::to_string(kickVelocityDistributionSigmaCCSN_BH) + " km s^-1 )").c_str())
		    ("kick-velocity-sigma-CCSN-NS",                                 po::value<double>()->default_value(kickVelocityDistributionSigmaCCSN_NS),                           ("Sigma for chosen kick velocity distribution for neutron stars (default = " + std::to_string(kickVelocityDistributionSigmaCCSN_NS) + " km s^-1 )").c_str())
			("kick-velocity-sigma-ECSN",                                    po::value<double>()->default_value(kickVelocityDistributionSigmaForECSN),                           ("Sigma for chosen kick velocity distribution for ECSN (default = " + std::to_string(kickVelocityDistributionSigmaForECSN) + " km s^-1 )").c_str())
			("kick-velocity-sigma-USSN",                                    po::value<double>()->default_value(kickVelocityDistributionSigmaForUSSN),                           ("Sigma for chosen kick velocity distribution for USSN (default = " + std::to_string(kickVelocityDistributionSigmaForUSSN) + " km s^-1 )").c_str())

		    ("luminous-blue-variable-multiplier",                           po::value<double>()->default_value(luminousBlueVariableFactor),                                     ("Multiplicitive constant for LBV mass loss (default = " + std::to_string(luminousBlueVariableFactor) + ", use 10 for Mennekens & Vanbeveren 2014)").c_str())

		    ("mass-ratio-max",                                              po::value<double>()->default_value(massRatioDistributionMax),                                       ("Maximum mass ratio m2/m1 to generate (default = " + std::to_string(massRatioDistributionMax) + ")").c_str())
		    ("mass-ratio-min",                                              po::value<double>()->default_value(massRatioDistributionMin),                                       ("Minimum mass ratio m2/m1 to generate (default = " + std::to_string(massRatioDistributionMin) + ")").c_str())
		    ("mass-transfer-fa",                                            po::value<double>()->default_value(massTransferFractionAccreted),                                   ("Mass Transfer fraction accreted (default = " + std::to_string(massTransferFractionAccreted) + ", fully conservative)").c_str())
		    ("mass-transfer-jloss",                                         po::value<double>()->default_value(massTransferJloss),                                              ("Specific angular momentum with which the non-accreted system leaves the system (default = " + std::to_string(massTransferJloss) + ")").c_str())
			("mass-transfer-thermal-limit-C",                               po::value<double>()->default_value(massTransferCParameter),                                         ("Mass Transfer Thermal rate factor fo the accretor (default = " + std::to_string(massTransferCParameter) + ", Hurley+2002)").c_str())
		    ("maximum-evolution-time",                                      po::value<double>()->default_value(maxEvolutionTime),                                               ("Maximum time to evolve binaries in Myrs (default = " + std::to_string(maxEvolutionTime) + ")").c_str())
		    ("maximum-mass-donor-Nandez-Ivanova",                           po::value<double>()->default_value(maximumMassDonorNandezIvanova),                                  ("Maximum donor mass allowed for the revised common envelope formalism in Msol (default = " + std::to_string(maximumMassDonorNandezIvanova) + ")").c_str())
			("maximum-neutron-star-mass",                                   po::value<double>()->default_value(maximumNeutronStarMass),                                         ("Maximum mass of a neutron star (default = " + std::to_string(maximumNeutronStarMass) + ", as in StarTrack)").c_str())
            ("metallicity,z",                                               po::value<double>()->default_value(metallicity),                                                    ("Metallicity to use (default " + std::to_string(metallicity) + " is Zsol)").c_str())
		    ("minimum-secondary-mass",                                      po::value<double>()->default_value(minimumMassSecondary),                                           ("Minimum mass of secondary to generate in Msol (default = " + std::to_string(minimumMassSecondary) + ")").c_str())

			("nbatches-used",                                               po::value<int>()->default_value(nBatchesUsed),                                                      ("Number of batches used, for STROOPWAFEL (AIS) (default = " + std::to_string(nBatchesUsed) + ")").c_str())

		    ("orbital-period-max",                                          po::value<double>()->default_value(periodDistributionMax),                                          ("Maximum period in days to generate (default = " + std::to_string(periodDistributionMax) + ")").c_str())
		   	("orbital-period-min",                                          po::value<double>()->default_value(periodDistributionMin),                                          ("Minimum period in days to generate (default = " + std::to_string(periodDistributionMin) + ")").c_str())

		    ("PISN-lower-limit",                                            po::value<double>()->default_value(pairInstabilityLowerLimit),                                      ("Minimum core mass for PISN (default = " + std::to_string(pairInstabilityLowerLimit) + ")").c_str())
		    ("PISN-upper-limit",                                            po::value<double>()->default_value(pairInstabilityUpperLimit),                                      ("Maximum core mass for PISN (default = " + std::to_string(pairInstabilityUpperLimit) + ")").c_str())
		    ("PPI-lower-limit",                                             po::value<double>()->default_value(pulsationalPairInstabilityLowerLimit),                           ("Minimum core mass for PPI (default = " + std::to_string(pulsationalPairInstabilityLowerLimit) + ")").c_str())
		    ("PPI-upper-limit",                                             po::value<double>()->default_value(pulsationalPairInstabilityUpperLimit),                           ("Maximum core mass for PPI (default = " + std::to_string(pulsationalPairInstabilityUpperLimit) + ")").c_str())
		    ("pulsar-birth-magnetic-field-distribution-max",                po::value<double>()->default_value(pulsarBirthMagneticFieldDistributionMax),                        ("Maximum (log10) pulsar birth magnetic field (default = " + std::to_string(pulsarBirthMagneticFieldDistributionMax) + ")").c_str())
		    ("pulsar-birth-magnetic-field-distribution-min",                po::value<double>()->default_value(pulsarBirthMagneticFieldDistributionMin),                        ("Minimum (log10) pulsar birth magnetic field) (default = " + std::to_string(pulsarBirthMagneticFieldDistributionMin) + ")").c_str())
		    ("pulsar-birth-spin-period-distribution-max",                   po::value<double>()->default_value(pulsarBirthSpinPeriodDistributionMax),                           ("Maximum pulsar birth spin period in ms (default = " + std::to_string(pulsarBirthSpinPeriodDistributionMax) + ")").c_str())
		    ("pulsar-birth-spin-period-distribution-min",                   po::value<double>()->default_value(pulsarBirthSpinPeriodDistributionMin),                           ("Minimum pulsar birth spin period in ms (default = " + std::to_string(pulsarBirthSpinPeriodDistributionMin) + ")").c_str())
		    ("pulsar-magnetic-field-decay-massscale",                       po::value<double>()->default_value(pulsarMagneticFieldDecayMassscale),                              ("Mass scale on which magnetic field decays during accretion in solar masses (default = " + std::to_string(pulsarMagneticFieldDecayMassscale) + ")").c_str())
		    ("pulsar-magnetic-field-decay-timescale",                       po::value<double>()->default_value(pulsarMagneticFieldDecayTimescale),                              ("Timescale on which magnetic field decays in Myrs (default = " + std::to_string(pulsarMagneticFieldDecayTimescale) + ")").c_str())
		    ("pulsar-minimum-magnetic-field",                               po::value<double>()->default_value(pulsarLog10MinimumMagneticField),                                ("log10 of the minimum pulsar magnetic field in Gauss (default = " + std::to_string(pulsarLog10MinimumMagneticField) + ")").c_str())

			("sample-common-envelope-alpha-max",                            po::value<double>()->default_value(sampleCommonEnvelopeAlphaMax),                                   ("Maximum for Uniform sampling over common envelope alpha (default = " + std::to_string(sampleCommonEnvelopeAlphaMax) + ")").c_str())
			("sample-common-envelope-alpha-min",                            po::value<double>()->default_value(sampleCommonEnvelopeAlphaMin),                                   ("Minimum for Uniform sampling over common envelope alpha (default = " + std::to_string(sampleCommonEnvelopeAlphaMin) + ")").c_str())
			("sample-kick-direction-power-max",                             po::value<double>()->default_value(sampleKickDirectionPowerMax),                                    ("Maximum for Uniform sampling over kick direction powerlaw exponent (default = " + std::to_string(sampleKickDirectionPowerMax) + ")").c_str())
			("sample-kick-direction-power-min",                             po::value<double>()->default_value(sampleKickDirectionPowerMin),                                    ("Minimum for Uniform sampling over kick direction powerlaw exponent (default = " + std::to_string(sampleKickDirectionPowerMin) + ")").c_str())
			("sample-kick-velocity-sigma-max",                              po::value<double>()->default_value(sampleKickVelocitySigmaMax),                                     ("Maximum for Uniform sampling over kick velocity sigma (default = " + std::to_string(sampleKickVelocitySigmaMax) + ")").c_str())
			("sample-kick-velocity-sigma-min",                              po::value<double>()->default_value(sampleKickVelocitySigmaMin),                                     ("Minimum for Uniform sampling over kick velocity sigma (default = " + std::to_string(sampleKickVelocitySigmaMin) + ")").c_str())
			("sample-luminous-blue-variable-multiplier-max",                po::value<double>()->default_value(sampleLuminousBlueVariableMultiplierMax),                        ("Maximum for Uniform sampling over multiplicative constant for LBV mass loss (default = " + std::to_string(sampleLuminousBlueVariableMultiplierMax) + ")").c_str())
			("sample-luminous-blue-variable-multiplier-min",                po::value<double>()->default_value(sampleLuminousBlueVariableMultiplierMin),                        ("Minimum for Uniform sampling over multiplicative constant for LBV mass loss (default = " + std::to_string(sampleLuminousBlueVariableMultiplierMin) + ")").c_str())
			("sample-wolf-rayet-multiplier-max",                            po::value<double>()->default_value(sampleWolfRayetMultiplierMax),                                   ("Maximum for Uniform sampling over multiplicative constant for WR winds (default = " + std::to_string(sampleWolfRayetMultiplierMax) + ")").c_str())
			("sample-wolf-rayet-multiplier-min",                            po::value<double>()->default_value(sampleWolfRayetMultiplierMin),                                   ("Minimum for Uniform sampling over multiplicative constant for WR winds (default = " + std::to_string(sampleWolfRayetMultiplierMin) + ")").c_str())
		    ("semi-major-axis-max",                                         po::value<double>()->default_value(semiMajorAxisDistributionMax),                                   ("Maximum semi major axis in AU to generate (default = " + std::to_string(semiMajorAxisDistributionMax) + ")").c_str())
		    ("semi-major-axis-min",                                         po::value<double>()->default_value(semiMajorAxisDistributionMin),                                   ("Minimum semi major axis in AU to generate (default = " + std::to_string(semiMajorAxisDistributionMin) + ")").c_str())
		    ("single-star-mass-max",                                        po::value<double>()->default_value(singleStarMassMax),                                              ("Maximum mass (in Msol) for single star evolution (default = " + std::to_string(singleStarMassMax) + ")").c_str())
            ("single-star-mass-min",                                        po::value<double>()->default_value(singleStarMassMin),                                              ("Minimum mass (in Msol) for single star evolution (default = " + std::to_string(singleStarMassMin) + ")").c_str())
		    ("spin-mag-max",                                                po::value<double>()->default_value(spinDistributionMax),                                            ("Maximum magnitude of spin (default = " + std::to_string(spinDistributionMax) + ")").c_str())
		    ("spin-mag-min",                                                po::value<double>()->default_value(spinDistributionMin),                                            ("Minimum magnitude of spin (default = " + std::to_string(spinDistributionMin) + ")").c_str())

		    ("wolf-rayet-multiplier",                                       po::value<double>()->default_value(wolfRayetFactor),                                                ("Multiplicitive constant for WR winds (default = " + std::to_string(wolfRayetFactor) + ")").c_str())

		    ("zeta-adiabatic-arbitrary",                                    po::value<double>()->default_value(zetaAdiabaticArbitrary),                                         ("Value of mass-radius exponent zeta adiabatic (default = " + std::to_string(zetaAdiabaticArbitrary) + ")").c_str())
		    ("zeta-hertzsprung-gap",                                        po::value<double>()->default_value(zetaHertzsprungGap),                                             ("Value of mass-radius exponent zeta on the hertzstrpung gap (default = " + std::to_string(zetaHertzsprungGap) + ")").c_str())
		    ("zeta-main-sequence",                                          po::value<double>()->default_value(zetaMainSequence),                                               ("Value of mass-radius exponent zeta on the main sequence (default = " + std::to_string(zetaMainSequence) + ")").c_str())


		    // string options - alphabetically
            ("AIS-DCOtype",                                                 po::value<string>()->default_value(AISDCOtypeString),                                               ("DCO type selection in exploratory phase of STROOPWAFEL, (options: ALL, BBH, BNS or BHNS), default = " + AISDCOtypeString + ")").c_str())

		  	("black-hole-kicks",                                            po::value<string>()->default_value(blackHoleKicksString),                                           ("Black hole kicks relative to NS kicks (options: FULL, REDUCED, ZERO, FALLBACK), default = " + blackHoleKicksString + ")").c_str())

		  	("chemically-homogeneous-evolution",                            po::value<string>()->default_value(cheString),                                                      ("Chemically Homogeneous Evolution (options: NONE, OPTIMISTIC, PESSIMISTIC), default = " + cheString + ")").c_str())

			("common-envelope-hertzsprung-gap-assumption",                  po::value<string>()->default_value(commonEnvelopeHertzsprungGapDonorString),                        ("Assumption to make about HG stars in CE (default = " + commonEnvelopeHertzsprungGapDonorString + ")").c_str())
			("common-envelope-lambda-prescription",                         po::value<string>()->default_value(commonEnvelopeLambdaPrescriptionString),                         ("Prescription for CE lambda (options: LAMBDA_FIXED, LAMBDA_LOVERIDGE, LAMBDA_NANJING, LAMBDA_KRUCKOW, LAMBDA_DEWI), default = " + commonEnvelopeLambdaPrescriptionString + ")").c_str())
		    ("common-envelope-mass-accretion-prescription",                 po::value<string>()->default_value(commonEnvelopeMassAccretionPrescriptionString),                  ("Assumption about whether NS/BHs can accrete mass during common envelope evolution (options: ZERO, CONSTANT, UNIFORM, MACLEOD+2014), default = " + commonEnvelopeMassAccretionPrescriptionString + ")").c_str())
			("common-envelope-zeta-prescription",                           po::value<string>()->default_value(commonEnvelopeZetaPrescriptionString),                           ("Prescription for CE zeta (default = " + commonEnvelopeZetaPrescriptionString + ")").c_str())

		    ("eccentricity-distribution,e",                                 po::value<string>()->default_value(eccentricityDistributionString),                                 ("Initial eccentricity distribution, e (options: ZERO, FIXED, FLAT, THERMALISED, GELLER+2013), default = " + eccentricityDistributionString + ")").c_str())

		    ("fryer-supernova-engine",                                      po::value<string>()->default_value(fryerSupernovaEngineString),                                     ("If using Fryer et al 2012 fallback prescription, select between 'delayed' and 'rapid' engines (default = " + fryerSupernovaEngineString + ")").c_str())

            ("grid",                                                        po::value<string>()->default_value(gridFilename)->implicit_value(""),                               ("Grid filename (default = " + gridFilename + ")").c_str())

		    ("initial-mass-function,i",                                     po::value<string>()->default_value(initialMassFunctionString),                                      ("Specify initial mass function to use (options: SALPETER, POWERLAW, UNIFORM, KROUPA), default = " + initialMassFunctionString + ")").c_str())

		    ("kick-direction",                                              po::value<string>()->default_value(kickDirectionDistributionString),                                ("Distribution for natal kick direction (options: ISOTROPIC, INPLANE, PERPENDICULAR, POWERLAW, WEDGE, POLES), default = " + kickDirectionDistributionString + ")").c_str())
		    ("kick-velocity-distribution",                                  po::value<string>()->default_value(kickVelocityDistributionString),                                 ("Natal kick velocity distribution (options: ZERO, FLAT, MAXWELLIAN, MUELLER2016, MUELLER2016MAXWELLIAN, BRAYELDRIDGE), default = " + kickVelocityDistributionString + ")").c_str())

            ("logfile-BSE-be-binaries",                                     po::value<string>()->default_value(logfileBSEBeBinaries),                                           ("Filename for BSE Be Binaries logfile (default = " + logfileBSEBeBinaries + ")").c_str())
            ("logfile-BSE-common-envelopes",                                po::value<string>()->default_value(logfileBSECommonEnvelopes),                                      ("Filename for BSE Common Envelopes logfile (default = " + logfileBSECommonEnvelopes + ")").c_str())
            ("logfile-BSE-detailed-output",                                 po::value<string>()->default_value(logfileBSEDetailedOutput),                                       ("Filename for BSE Detailed Output logfile (default = " + logfileBSEDetailedOutput + ")").c_str())
            ("logfile-BSE-double-compact-objects",                          po::value<string>()->default_value(logfileBSEDoubleCompactObjects),                                 ("Filename for BSE Double Compact Objects logfile (default = " + logfileBSEDoubleCompactObjects + ")").c_str())
            ("logfile-BSE-pulsar-evolution",                                po::value<string>()->default_value(logfileBSEPulsarEvolution),                                      ("Filename for BSE Pulsar Evolution logfile (default = " + logfileBSEPulsarEvolution + ")").c_str())
            ("logfile-BSE-supernovae",                                      po::value<string>()->default_value(logfileBSESupernovae),                                           ("Filename for BSE Supernovae logfile (default = " + logfileBSESupernovae + ")").c_str())
            ("logfile-BSE-system-parameters",                               po::value<string>()->default_value(logfileBSESystemParameters),                                     ("Filename for BSE System Parameters logfile (default = " + logfileBSESystemParameters + ")").c_str())
            ("logfile-definitions",                                         po::value<string>()->default_value(logfileDefinitionsFilename)->implicit_value(""),                 ("Filename for logfile record definitions (default = " + logfileDefinitionsFilename + ")").c_str())
            ("logfile-delimiter",                                           po::value<string>()->default_value(logfileDelimiterString),                                         ("Field delimiter for logfile records (default = " + logfileDelimiterString + ")").c_str())
            ("logfile-name-prefix",                                         po::value<string>()->default_value(logfileNamePrefix),                                              ("Prefix for logfile names (default = " + logfileNamePrefix + ")").c_str())
            ("logfile-SSE-parameters",                                      po::value<string>()->default_value(logfileSSEParameters),                                           ("Filename for SSE Parameters logfile (default = " + logfileSSEParameters + ")").c_str())

		    ("mass-loss-prescription",                                      po::value<string>()->default_value(massLossPrescriptionString),                                     ("Mass loss prescription to use (options: NONE, HURLEY, VINK), default = " + massLossPrescriptionString + ")").c_str())
		    ("mass-ratio-distribution,q",                                   po::value<string>()->default_value(massRatioDistributionString),                                    ("Initial mass ratio distribution for q=m2/m1 (options: FLAT, DuquennoyMayor1991, SANA2012), default = " + massRatioDistributionString + ")").c_str())
		    ("mass-transfer-accretion-efficiency-prescription",             po::value<string>()->default_value(massTransferAccretionEfficiencyPrescriptionString),              ("Mass Transfer Accretion Efficiency prescription to use (options: THERMAL, FIXED, CENTRIFUGAL), default = " + massTransferAngularMomentumLossPrescriptionString + ")").c_str())
		    ("mass-transfer-angular-momentum-loss-prescription",            po::value<string>()->default_value(massTransferAngularMomentumLossPrescriptionString),              ("Mass Transfer Angular Momentum Loss prescription to use (options: JEANS, ISOTROPIC, CIRCUMBINARY, ARBITRARY), default = " + massTransferAngularMomentumLossPrescriptionString + ")").c_str())
		    ("mass-transfer-prescription",                                  po::value<string>()->default_value(massTransferPrescriptionString),                                 ("Mass Transfer prescription to use (default = " + massTransferPrescriptionString + ")").c_str())
		    ("mass-transfer-rejuvenation-prescription",                     po::value<string>()->default_value(massTransferRejuvenationPrescriptionString),                     ("Mass Transfer Rejuvenation prescription to use (options: NONE, STARTRACK), default = " + massTransferRejuvenationPrescriptionString + ")").c_str())
			("mass-transfer-thermal-limit-accretor",                        po::value<string>()->default_value(massTransferThermallyLimitedVariationString),                    ("Mass Transfer Thermal Accretion limit to use (default = " + massTransferThermallyLimitedVariationString + ")").c_str())

		    ("neutron-star-equation-of-state",                              po::value<string>()->default_value(neutronStarEquationOfStateString),                               ("Specify which neutron star equation of state to use (options: SSE, ARP3), default = " + neutronStarEquationOfStateString + ")").c_str())

   		    ("outputPath,o",                                                po::value<string>()->default_value(outputPathString),                                               ("Directory for output (default = " + outputPathString + ")").c_str())

		    ("pulsar-birth-magnetic-field-distribution",                    po::value<string>()->default_value(pulsarBirthMagneticFieldDistributionString),                     ("Distribution of (log10 of) pulsar birth magnetic field in G (options: ZERO, FIXED, FLATINLOG, UNIFORM, LOGNORMAL), default = " + pulsarBirthMagneticFieldDistributionString + ")").c_str())
		    ("pulsar-birth-spin-period-distribution",                       po::value<string>()->default_value(pulsarBirthSpinPeriodDistributionString),                        ("Distribution of pulsar birth spin period in ms (options: ZERO, FIXED, UNIFORM, NORMAL), default = " + pulsarBirthSpinPeriodDistributionString + ")").c_str())
		    ("pulsational-pair-instability-prescription",                   po::value<string>()->default_value(pulsationalPairInstabilityPrescriptionString),                   ("Specify which prescription to use for pulsational pair instability (options: COMPAS, STARTRACK, MARCHANT), default = " + pulsationalPairInstabilityPrescriptionString + ")").c_str())

		    ("remnant-mass-prescription",                                   po::value<string>()->default_value(remnantMassPrescriptionString),                                  ("Choose remnant mass prescription (options: postitnote, hurley2000, belczynski2002, fryer2012, muller2016, muller2016Maxwellian), default = " + remnantMassPrescriptionString + ")").c_str())
		    ("rotational-velocity-distribution",                            po::value<string>()->default_value(rotationalVelocityDistributionString),                           ("Initial rotational velocity distribution (options: ZERO, HURLEY, VLTFLAMES), default = " + rotationalVelocityDistributionString + ")").c_str())

		    ("semi-major-axis-distribution,a",                              po::value<string>()->default_value(semiMajorAxisDistributionString),                                ("Initial semi-major axis distribution, a (options: FLATINLOG, CUSTOM, DuquennoyMayor1991, SANA2012), default = " + semiMajorAxisDistributionString + ")").c_str())

		    ("spin-assumption",                                             po::value<string>()->default_value(spinAssumptionString),                                           ("Which assumption of misalignedments to use (default = " + spinAssumptionString + ")").c_str())
		    ("spin-distribution",                                           po::value<string>()->default_value(spinDistributionString),                                         ("Which distribution of spins to use (default = " + spinDistributionString + ")").c_str())

		    ("tides-prescription",                                          po::value<string>()->default_value(tidesPrescriptionString),                                        ("Tides prescription to use (options: default = " + tidesPrescriptionString + ")").c_str())


            // vector (list) options - alphabetically
            ("debug-classes",                                               po::value<vector<string>>()->multitoken()->default_value(debugClasses),                             ("Debug classes enabled (default = " + defaultDebugClasses + ")").c_str())
            ("log-classes",                                                 po::value<vector<string>>()->multitoken()->default_value(logClasses),                               ("Logging classes enabled (default = " + defaultLogClasses + ")").c_str())
		;

        po::variables_map vm;   // Variables map

        try {

            po::store(po::parse_command_line(argc, argv, desc), vm);

            // --help option
            if (vm["help"].as<bool>()) {                                                                                                     // user requested help?
                utils::SplashScreen();                                                                                                  // yes - show splash screen
                ANNOUNCE(desc);                                                                                                         // and help
                programStatus = COMMANDLINE_STATUS::SUCCESS;                                                                            // ok
            }

            // --version option
            if (vm["version"].as<bool>()) {                                                                                                  // user requested version?
                ANNOUNCE("COMPAS v" << VERSION_STRING);                                                                                 // yes, show version string
                programStatus = COMMANDLINE_STATUS::SUCCESS;                                                                            // ok
            }

            po::notify(vm);                                                                                                             // invoke notify to assign user-input values to variables.  Throws an error, so do after help just in case there are any problems.

PrintProgramOptionDetails(vm);

            // boolean options (flags) - alphabetically (where possible - dependencies)

//            AISexploratoryPhase                             = vm["AIS-exploratory-phase"].as<bool>();                                        // exploratory phase of Adaptive Importance Sampling - Floor 24-04-2018.  Do not retain previous (default) value
//            AIShubble                                       = vm["AIS-Hubble"].as<bool>();                                                   // excluding binaries that merge outside Hubble time (exploratory phase AIS)?  Do not retain previous (default) value
            AISpessimistic                                  = vm["AIS-Pessimistic"].as<bool>();                                              // excluding binaries that are Optimistic (exploratory phase AIS)?  Do not retain previous (default) value
            AISrefinementPhase                              = vm["AIS-refinement-phase"].as<bool>();
            AISrlof                                         = vm["AIS-RLOF"].as<bool>();                                                     // excluding binaries that RLOFSecondaryZAMS (exploratory phase AIS)?  Do not retain previous (default) value

            allowTouchingAtBirth                            = vm["allow-touching-at-birth"].as<bool>();                                      // allow binaries that are touching at birth to evolve
            allowRLOFAtBirth                                = vm["allow-rlof-at-birth"].as<bool>();                                          // allow binaries that have one or both stars in RLOF at birth to evolve

            allowMainSequenceStarToSurviveCommonEnvelope    = vm["common-envelope-allow-main-sequence-survive"].as<bool>();                  // allow MS stars to survive CE event?  Do not retain previous (default) value

			beBinaries                                      = vm["BeBinaries"].as<bool>();                                                   // enable Be Binaries?  Do not retain previous (default) value

            debugToFile                                     = vm.count("debug-to-file") ? vm["debug-to-file"].as<bool>() : debugToFile;      // write debug records to file?  Retain previous (default) value if not specified

            detailedOutput                                  = vm["detailedOutput"].as<bool>();                                               // detailed output of each simulated system?  Do not retain previous (default) value

            errorsToFile                                    = vm.count("errors-to-file") ? vm["errors-to-file"].as<bool>() : errorsToFile;   // write error messages to file?  Retain previous (default) value

            evolvePulsars                                   = vm["evolve-pulsars"].as<bool>();                                               // evolve pulsars?  Do not retain previous (default) value

            evolveUnboundSystems                            = vm["evolve-unbound-systems"].as<bool>();                                       // evolve unbound systems?  Do not retain previous (default) value

            fixedMetallicity                                = vm.count("metallicity") || fixedMetallicity;     // user-specified a metallicity value?  Retain previous (default) value

            fixedRandomSeed                                 = vm.count("random-seed") || fixedRandomSeed;                                 // user-specified random seed?  Retain previous (default) value

			lambdaCalculationEveryTimeStep                  = vm["lambda-calculation-every-timeStep"].as<bool>();                            // calculate lambdas at every timestep?  Do not retain previous (default) value if not specified

            massTransferCriticalMassRatioMSLowMass          = vm.count("critical-mass-ratio-MS-low-mass-non-degenerate-accretor") &&
                                                              vm.count("critical-mass-ratio-MS-low-mass-degenerate-accretor")     &&
                                                              massTransferCriticalMassRatioMSLowMassNonDegenerateAccretor >= 0.0  &&
                                                              massTransferCriticalMassRatioMSLowMassDegenerateAccretor    >= 0.0;       // Do not retain previous (default) value

            massTransferCriticalMassRatioMSHighMass         = vm.count("critical-mass-ratio-MS-high-mass-non-degenerate-accretor") &&
                                                              vm.count("critical-mass-ratio-MS-high-mass-degenerate-accretor")     &&
                                                              massTransferCriticalMassRatioMSHighMassNonDegenerateAccretor >= 0.0  &&
                                                              massTransferCriticalMassRatioMSHighMassDegenerateAccretor    >= 0.0;      // Do not retain previous (default) value

            massTransferCriticalMassRatioHG                 = vm.count("critical-mass-ratio-HG-non-degenerate-accretor")  &&
                                                              vm.count("critical-mass-ratio-HG-degenerate-accretor")      &&
                                                              massTransferCriticalMassRatioHGNonDegenerateAccretor >= 0.0 &&
                                                              massTransferCriticalMassRatioHGDegenerateAccretor    >= 0.0;              // Do not retain previous (default) value

            massTransferCriticalMassRatioGiant              = vm.count("critical-mass-ratio-giant-non-degenerate-accretor")  &&
                                                              vm.count("critical-mass-ratio-giant-degenerate-accretor")      &&
                                                              massTransferCriticalMassRatioGiantNonDegenerateAccretor >= 0.0 &&
                                                              massTransferCriticalMassRatioGiantDegenerateAccretor    >= 0.0;           // Do not retain previous (default) value

            massTransferCriticalMassRatioHeliumMS           = vm.count("critical-mass-ratio-helium-MS-non-degenerate-accretor") &&
                                                              vm.count("critical-mass-ratio-helium-MS-degenerate-accretor")     &&
                                                              massTransferCriticalMassRatioHeliumMSNonDegenerateAccretor >= 0.0 &&
                                                              massTransferCriticalMassRatioHeliumMSDegenerateAccretor    >= 0.0;        // Do not retain previous (default) value

            massTransferCriticalMassRatioHeliumHG           = vm.count("critical-mass-ratio-helium-HG-non-degenerate-accretor") &&
                                                              vm.count("critical-mass-ratio-helium-HG-degenerate-accretor")     &&
                                                              massTransferCriticalMassRatioHeliumHGNonDegenerateAccretor >= 0.0 &&
                                                              massTransferCriticalMassRatioHeliumHGDegenerateAccretor    >= 0.0;        // Do not retain previous (default) value

            massTransferCriticalMassRatioHeliumGiant        = vm.count("critical-mass-ratio-helium-giant-non-degenerate-accretor") &&
                                                              vm.count("critical-mass-ratio-helium-giant-degenerate-accretor")     &&
                                                              massTransferCriticalMassRatioHeliumGiantNonDegenerateAccretor >= 0.0 &&
                                                              massTransferCriticalMassRatioHeliumGiantDegenerateAccretor    >= 0.0;     // Do not retain previous (default) value

            massTransferCriticalMassRatioWhiteDwarf         = vm.count("critical-mass-ratio-white-dwarf-non-degenerate-accretor") &&
                                                              vm.count("critical-mass-ratio-white-dwarf-degenerate-accretor")     &&
                                                              massTransferCriticalMassRatioWhiteDwarfNonDegenerateAccretor >= 0.0 &&
                                                              massTransferCriticalMassRatioWhiteDwarfDegenerateAccretor    >= 0.0;      // Do not retain previous (default) value

            onlyDoubleCompactObjects                        = vm["only-double-compact-objects"].as<bool>();                             // only evolve DCOs?  Do not retain previous (default) value

            PNevolution                                     = vm.count("PN") ? vm["PN"].as<bool>() : PNevolution;                       // integrate spins using PN equations?  Retain previous (default) value

            populationDataPrinting                          = vm["populationDataPrinting"].as<bool>();                                  // print certain values while running a population?  Do not retain previous (default) valued

            printBoolAsString                               = vm.count("print-bool-as-string") ? vm["print-bool-as-string"].as<bool>() : printBoolAsString; // print boolean values as "TRUE"/"FALSE"?  Retain previous (default) value

            quiet                                           = vm["quiet"].as<bool>();                                                   // verbose or quiet mode?  Do not retain previous (default) value

            revisedEnergyFormalismNandezIvanova             = vm["revised-energy-formalism-Nandez-Ivanova"].as<bool>();                 // Do not retain previous (default) value

			sampleCommonEnvelopeAlpha                       = vm["sample-common-envelope-alpha"].as<bool>();                            // sample CE alpha?  Do not retain previous (default) value
			sampleKickDirectionPower                        = vm["sample-kick-direction-power"].as<bool>();                             // sample kick direction power?  Do not retain previous (default) value
			sampleKickVelocitySigma                         = vm["sample-kick-velocity-sigma"].as<bool>();                              // sample over some hyperparameters? - JIM BARRETT -- 06/07/2016.  Do not retain previous (default) value
			sampleLuminousBlueVariableMultiplier            = vm["sample-luminous-blue-variable-multiplier"].as<bool>();                // sample LBV multiplier?  Do not retain previous (default) value
			sampleWolfRayetMultiplier                       = vm["sample-wolf-rayet-multiplier"].as<bool>();                            // sample wolf-rayet multiplier?  Do not retain previous (default) valued

            singleStar                                      = vm["single-star"].as<bool>();                                             // evolve a single star?  Do not retain previous (default) valued

            useFixedUK                                      = vm.count("fix-dimensionless-kick-velocity") ? (utils::Compare(fixedUK, 0.0) >= 0) : false;   // fix the dimensionless kick velocity?  Do not retain previous (default) value

            useMassTransfer                                 = vm.count("massTransfer") ? vm["massTransfer"].as<bool>() : useMassTransfer;   // use mass transfer?  Retain previous (default) value

            useMassLoss                                     = vm["use-mass-loss"].as<bool>();                                           // use mass loss?  Do not retain previous (default) value

            usePairInstabilitySupernovae                    = vm["pair-instability-supernovae"].as<bool>();                             // pair instability supernovae?  Do not retain previous (default) value

            usePulsationalPairInstability                   = vm["pulsational-pair-instability"].as<bool>();                            // pulsational pair instability supernovae?  Do not retain previous (default) value

			zetaCalculationEveryTimeStep                    = vm["zeta-calculation-every-timestep"].as<bool>();                        // calculate zetas at every timestep?  Do not retain previous (default) value

			if (useMassTransfer) {                                                                                                      // the following depend on the value of useMassTransfer
                if (vm.count("circulariseBinaryDuringMassTransfer")) {
					angularMomentumConservationDuringCircularisation = vm["angularMomentumConservationDuringCircularisation"].as<bool>();   // Do not retain previous (default) value
					circulariseBinaryDuringMassTransfer              = true;
                }
				else{
					circulariseBinaryDuringMassTransfer              = false;
				}

                if (vm.count("forceCaseBBBCStabilityFlag")) {
	                alwaysStableCaseBBBCFlag   = vm["alwaysStableCaseBBBCFlag"].as<bool>();                                             // Do not retain previous (default) value
					forceCaseBBBCStabilityFlag = true;
                }
				else{
					forceCaseBBBCStabilityFlag = false;
				}
            }


            // numerical options - alphabetically group by type (where possible - dependencies)

            randomSeed                               = vm.count("random-seed") ? vm["random-seed"].as<unsigned long>() : randomSeed;                                                                                            // Random Seed

            debugLevel                               = vm.count("debug-level") ? vm["debug-level"].as<int>() : debugLevel;                                                                                                      // Debug Level
            logLevel                                 = vm.count("log-level") ? vm["log-level"].as<int>() : logLevel;                                                                                                            // Logging Level
            maxNumberOfTimestepIterations            = vm.count("maximum-number-iterations") ? vm["maximum-number-iterations"].as<int>() : maxNumberOfTimestepIterations;                                                       // Maximum number of timesteps to evolve binary for before giving up
            nBatchesUsed                             = vm.count("nbatches-used") ? vm["nbatches-used"].as<int>() : nBatchesUsed;                                                                                                // Number of batches used, only needed for STROOPWAFEL (AIS)
            nBinaries                                = vm.count("number-of-binaries,n") ? vm["number-of-binaries,n"].as<int>() : nBinaries;                                                                                     // Number of binaries to simulate
            singleStarMassSteps                      = vm.count("single-star-mass-steps") ? vm["single-star-mass-steps"].as<int>() : singleStarMassSteps;                                                                       // Number of mass steps for single star evolution

            commonEnvelopeAlpha                      = vm.count("common-envelope-alpha") ? vm["common-envelope-alpha"].as<double>() : commonEnvelopeAlpha;                                                                     // Common envelope efficiency alpha parameter
            commonEnvelopeAlphaThermal               = vm.count("common-envelope-alpha-thermal") ? vm["common-envelope-alpha-thermal"].as<double>() : commonEnvelopeAlphaThermal;                                              // lambda = (alpha_th * lambda_b) + (1-alpha_th) * lambda_g
            commonEnvelopeLambda                     = vm.count("common-envelope-lambda") ? vm["common-envelope-lambda"].as<double>() : commonEnvelopeLambda;                                                                  // Common Envelope lambda
            commonEnvelopeLambdaMultiplier           = vm.count("common-envelope-lambda-multiplier") ? vm["common-envelope-lambda-multiplier"].as<double>() : commonEnvelopeLambdaMultiplier;                                  // Multiply lambda by some constant
            commonEnvelopeMassAccretionConstant      = vm.count("common-envelope-mass-accretion-constant") ? vm["common-envelope-mass-accretion-constant"].as<double>() : commonEnvelopeMassAccretionConstant;                 // Value of mass accreted by NS/BH during common envelope evolution if assuming all NS/BH accrete same amount of mass
            commonEnvelopeMassAccretionMax           = vm.count("common-envelope-mass-accretion-max") ? vm["common-envelope-mass-accretion-max"].as<double>() : commonEnvelopeMassAccretionMax;                                // Maximum amount of mass accreted by NS/BHs during common envelope evolution in solar masses
            commonEnvelopeMassAccretionMin           = vm.count("common-envelope-mass-accretion-min") ? vm["common-envelope-mass-accretion-min"].as<double>() : commonEnvelopeMassAccretionMin;                                // Minimum amount of mass accreted by NS/BHs during common envelope evolution in solar masses
            commonEnvelopeRecombinationEnergyDensity = vm.count("common-envelope-recombination-energy-density") ? vm["common-envelope-recombination-energy-density"].as<double>() : commonEnvelopeRecombinationEnergyDensity;  // Recombination energy density in ergs/g
            commonEnvelopeSlopeKruckow               = vm.count("common-envelope-slope-Kruckow") ? vm["common-envelope-slope-Kruckow"].as<double>() : commonEnvelopeSlopeKruckow;                                              // Common Envelope slope for Kruckow lambda

            massTransferCriticalMassRatioGiantDegenerateAccretor          = vm.count("critical-mass-ratio-giant-degenerate-accretor") ? vm["critical-mass-ratio-giant-degenerate-accretor"].as<double>() : massTransferCriticalMassRatioGiantDegenerateAccretor;                                // Critical mass ratio for MT from a giant on to a degenerate accretor
            massTransferCriticalMassRatioGiantNonDegenerateAccretor       = vm.count("critical-mass-ratio-giant-non-degenerate-accretor") ? vm["critical-mass-ratio-giant-non-degenerate-accretor"].as<double>() : massTransferCriticalMassRatioGiantNonDegenerateAccretor;                     // Critical mass ratio for MT from a giant on to a non degenerate accretor
            massTransferCriticalMassRatioHeliumGiantDegenerateAccretor    = vm.count("critical-mass-ratio-helium-giant-degenerate-accretor") ? vm["critical-mass-ratio-helium-giant-degenerate-accretor"].as<double>() : massTransferCriticalMassRatioHeliumGiantDegenerateAccretor;            // Critical mass ratio for MT from a Helium giant on to a degenerate accretor
            massTransferCriticalMassRatioHeliumGiantNonDegenerateAccretor = vm.count("critical-mass-ratio-helium-giant-non-degenerate-accretor") ? vm["critical-mass-ratio-helium-giant-non-degenerate-accretor"].as<double>() : massTransferCriticalMassRatioHeliumGiantNonDegenerateAccretor; // Critical mass ratio for MT from a Helium giant on to a non degenerate accretor
            massTransferCriticalMassRatioHeliumHGDegenerateAccretor       = vm.count("critical-mass-ratio-helium-HG-degenerate-accretor") ? vm["critical-mass-ratio-helium-HG-degenerate-accretor"].as<double>() : massTransferCriticalMassRatioHeliumHGDegenerateAccretor;                     // Critical mass ratio for MT from a Helium HG star on to a degenerate accretor
            massTransferCriticalMassRatioHeliumHGNonDegenerateAccretor    = vm.count("critical-mass-ratio-helium-HG-non-degenerate-accretor") ? vm["critical-mass-ratio-helium-HG-non-degenerate-accretor"].as<double>() : massTransferCriticalMassRatioHeliumHGNonDegenerateAccretor;          // Critical mass ratio for MT from a Helium HG star on to a non degenerate accretor
            massTransferCriticalMassRatioHeliumMSDegenerateAccretor       = vm.count("critical-mass-ratio-helium-MS-degenerate-accretor") ? vm["critical-mass-ratio-helium-MS-degenerate-accretor"].as<double>() : massTransferCriticalMassRatioHeliumMSDegenerateAccretor;                     // Critical mass ratio for MT from a Helium MS star on to a degenerate accretor
            massTransferCriticalMassRatioHeliumMSNonDegenerateAccretor    = vm.count("critical-mass-ratio-helium-MS-non-degenerate-accretor") ? vm["critical-mass-ratio-helium-MS-non-degenerate-accretor"].as<double>() : massTransferCriticalMassRatioHeliumMSNonDegenerateAccretor;          // Critical mass ratio for MT from a Helium MS star on to a non degenerate accretor
            massTransferCriticalMassRatioHGDegenerateAccretor             = vm.count("critical-mass-ratio-HG-degenerate-accretor") ? vm["critical-mass-ratio-HG-degenerate-accretor"].as<double>() : massTransferCriticalMassRatioHGDegenerateAccretor;                                         // Critical mass ratio for MT from a HG star on to a degenerate accretor
            massTransferCriticalMassRatioHGNonDegenerateAccretor          = vm.count("critical-mass-ratio-HG-non-degenerate-accretor") ? vm["critical-mass-ratio-HG-non-degenerate-accretor"].as<double>() : massTransferCriticalMassRatioHGNonDegenerateAccretor;                              // Critical mass ratio for MT from a HG star on to a non degenerate accretor
            massTransferCriticalMassRatioMSHighMassDegenerateAccretor     = vm.count("critical-mass-ratio-MS-high-mass-degenerate-accretor") ? vm["critical-mass-ratio-MS-high-mass-degenerate-accretor"].as<double>() : massTransferCriticalMassRatioMSHighMassDegenerateAccretor;             // Critical mass ratio for MT from a MS high mass star on to a degenerate accretor
            massTransferCriticalMassRatioMSHighMassNonDegenerateAccretor  = vm.count("critical-mass-ratio-MS-high-mass-non-degenerate-accretor") ? vm["critical-mass-ratio-MS-high-mass-non-degenerate-accretor"].as<double>() : massTransferCriticalMassRatioMSHighMassNonDegenerateAccretor;  // Critical mass ratio for MT from a MS high mass star on to a non degenerate accretor
            massTransferCriticalMassRatioMSLowMassDegenerateAccretor      = vm.count("critical-mass-ratio-MS-low-mass-degenerate-accretor") ? vm["critical-mass-ratio-MS-low-mass-degenerate-accretor"].as<double>() : massTransferCriticalMassRatioMSLowMassDegenerateAccretor;                // Critical mass ratio for MT from a MS low mass star on to a degenerate accretor
            massTransferCriticalMassRatioMSLowMassNonDegenerateAccretor   = vm.count("critical-mass-ratio-MS-low-mass-non-degenerate-accretor") ? vm["critical-mass-ratio-MS-low-mass-non-degenerate-accretor"].as<double>() : massTransferCriticalMassRatioMSLowMassNonDegenerateAccretor;     // Critical mass ratio for MT from a MS low mass star on to a non degenerate accretor
            massTransferCriticalMassRatioWhiteDwarfDegenerateAccretor     = vm.count("critical-mass-ratio-white-dwarf-degenerate-accretor") ? vm["critical-mass-ratio-white-dwarf-degenerate-accretor"].as<double>() : massTransferCriticalMassRatioWhiteDwarfDegenerateAccretor;               // Critical mass ratio for MT from a White Dwarf on to a degenerate accretor
            massTransferCriticalMassRatioWhiteDwarfNonDegenerateAccretor  = vm.count("critical-mass-ratio-white-dwarf-non-degenerate-accretor") ? vm["critical-mass-ratio-white-dwarf-non-degenerate-accretor"].as<double>() : massTransferCriticalMassRatioWhiteDwarfNonDegenerateAccretor;    // Critical mass ratio for MT from a White Dwarf on to a non degenerate accretor

            eccentricityDistributionMax             = vm.count("eccentricity-max") ? vm["eccentricity-max"].as<double>() : eccentricityDistributionMax;                                     // Maximum eccentricity to generate
            eccentricityDistributionMin             = vm.count("eccentricity-min") ? vm["eccentricity-min"].as<double>() : eccentricityDistributionMin;                                     // Minimum eccentricity to generate

            eddingtonAccretionFactor                = vm.count("eddington-accretion-factor") ? vm["eddington-accretion-factor"].as<double>() : eddingtonAccretionFactor;                    // Multiplication factor for eddington accretion for NS & BH

            fixedUK                                 = vm.count("fix-dimensionless-kick-velocity") ? vm["fix-dimensionless-kick-velocity"].as<double>() : fixedUK;                           // Value to which dimensionless kick velocity uk is fixed

            initialMassFunctionMax                  = vm.count("initial-mass-max") ? vm["initial-mass-max"].as<double>() : initialMassFunctionMax;                                          // Maximum mass (in Msol) to generate using given IMF
            initialMassFunctionMin                  = vm.count("initial-mass-min") ? vm["initial-mass-min"].as<double>() : initialMassFunctionMin;                                          // Minimum mass (in Msol) to generate using given IMF
            initialMassFunctionPower                = vm.count("initial-mass-power") ? vm["initial-mass-power"].as<double>() : initialMassFunctionPower;                                    // Single power law power to generate primary mass using given IMF

            kappaGaussians                          = vm.count("kappa-gaussians") ? vm["kappa-gaussians"].as<double>() : kappaGaussians;                                                    // Scaling factor for the width of the Gaussian distributions in STROOPWAFEL main sampling phase

            kickDirectionPower                      = vm.count("kick-direction-power") ? vm["kick-direction-power"].as<double>() : kickDirectionPower;                                      // Power law power for the "power" SN kick direction choice
            kickScalingFactor                       = vm.count("kick-scaling-factor") ? vm["kick-scaling-factor"].as<double>() : kickScalingFactor;                                         // Arbitrary factor for scaling kicks
            kickVelocityDistributionMaximum         = vm.count("kick-velocity-max") ? vm["kick-velocity-max"].as<double>() : kickVelocityDistributionMaximum;                               // Maximum kick velocity to draw in km s^-1
            kickVelocityDistributionSigmaCCSN_BH    = vm.count("kick-velocity-sigma-CCSN-BH") ? vm["kick-velocity-sigma-CCSN-BH"].as<double>() : kickVelocityDistributionSigmaCCSN_BH;      // Kick velocity sigma in km s^-1 for black holes
            kickVelocityDistributionSigmaCCSN_NS    = vm.count("kick-velocity-sigma-CCSN-NS") ? vm["kick-velocity-sigma-CCSN-NS"].as<double>() : kickVelocityDistributionSigmaCCSN_NS;      // Kick velocity sigma in km s^-1 for neutron stars
            kickVelocityDistributionSigmaForECSN    = vm.count("kick-velocity-sigma-ECSN") ? vm["kick-velocity-sigma-ECSN"].as<double>() : kickVelocityDistributionSigmaForECSN;            // Characteristic kick velocity for an ECSN in km s^-1
            kickVelocityDistributionSigmaForUSSN    = vm.count("kick-velocity-sigma-USSN") ? vm["kick-velocity-sigma-USSN"].as<double>() : kickVelocityDistributionSigmaForUSSN;            // Characteristic kick velocity for an USSN in km s^-1

            luminousBlueVariableFactor              = vm.count("luminous-blue-variable-multiplier") ? vm["luminous-blue-variable-multiplier"].as<double>() : luminousBlueVariableFactor;    // Luminous blue variable mass loss enhancement factor

            massRatioDistributionMax                = vm.count("mass-ratio-max") ? vm["mass-ratio-max"].as<double>() : massRatioDistributionMax;                                            // Maximum mass ratio m2/m1 to generate
            massRatioDistributionMin                = vm.count("mass-ratio-min") ? vm["mass-ratio-min"].as<double>() : massRatioDistributionMin;                                            // Minimum mass ratio m2/m1 to generate

            massTransferFractionAccreted            = vm.count("mass-transfer-fa") ? vm["mass-transfer-fa"].as<double>() : massTransferFractionAccreted;                                    // Mass Transfer fraction accreted

            massTransferJloss                       = vm.count("mass-transfer-jloss") ? vm["mass-transfer-jloss"].as<double>() : massTransferJloss;                                         // Specific angular momentum with which the non-accreted system leaves the system
            massTransferCParameter                  = vm.count("mass-transfer-thermal-limit-C") ? vm["mass-transfer-thermal-limit-C"].as<double>() : massTransferCParameter;                // Mass Transfer Thermal rate factor fo the accretor

            maxEvolutionTime                        = vm.count("maximum-evolution-time") ? vm["maximum-evolution-time"].as<double>() : maxEvolutionTime;                                    // Maximum evolution time in Myrs

            maximumMassDonorNandezIvanova           = vm.count("maximum-mass-donor-Nandez-Ivanova") ? vm["maximum-mass-donor-Nandez-Ivanova"].as<double>() : maximumMassDonorNandezIvanova; // Maximum mass allowed to use the revised energy formalism in Msol
            maximumNeutronStarMass                  = vm.count("maximum-neutron-star-mass") ? vm["maximum-neutron-star-mass"].as<double>() : maximumNeutronStarMass;                        // Maximum mass of a neutron star allowed

            metallicity                             = vm.count("metallicity,z") ? vm["metallicity,z"].as<double>() : metallicity;                                                           // Metallicity to use

            minimumMassSecondary                    = vm.count("minimum-secondary-mass") ? vm["minimum-secondary-mass"].as<double>() : minimumMassSecondary;                                // Minimum mass of secondary to draw (in Msol)

            periodDistributionMax                   = vm.count("orbital-period-max") ? vm["orbital-period-max"].as<double>() : periodDistributionMax;                                                                       // Maximum initial period in days
            periodDistributionMin                   = vm.count("orbital-period-min") ? vm["orbital-period-min"].as<double>() : periodDistributionMin;                                                                       // Minimum initial period in days

            pairInstabilityLowerLimit               = vm.count("PISN-lower-limit") ? vm["PISN-lower-limit"].as<double>() : pairInstabilityLowerLimit;                                                                       // Minimum core mass leading to PISN
            pairInstabilityUpperLimit               = vm.count("PISN-upper-limit") ? vm["PISN-upper-limit"].as<double>() : pairInstabilityUpperLimit;                                                                       // Maximum core mass leading to PISN
            pulsationalPairInstabilityLowerLimit    = vm.count("PPI-lower-limit") ? vm["PPI-lower-limit"].as<double>() : pulsationalPairInstabilityLowerLimit;                                                              // Minimum core mass leading to PPISN
            pulsationalPairInstabilityUpperLimit    = vm.count("PPI-upper-limit") ? vm["PPI-upper-limit"].as<double>() : pulsationalPairInstabilityUpperLimit;                                                              // Maximum core mass leading to PPISN
            pulsarBirthMagneticFieldDistributionMax = vm.count("pulsar-birth-magnetic-field-distribution-max") ? vm["pulsar-birth-magnetic-field-distribution-max"].as<double>() : pulsarBirthMagneticFieldDistributionMax; // Maximum pulsar birth magnetic field distribution
            pulsarBirthMagneticFieldDistributionMin = vm.count("pulsar-birth-magnetic-field-distribution-min") ? vm["pulsar-birth-magnetic-field-distribution-min"].as<double>() : pulsarBirthMagneticFieldDistributionMin; // Minimum pulsar birth magnetic field distribution
            pulsarBirthSpinPeriodDistributionMax    = vm.count("pulsar-birth-spin-period-distribution-max") ? vm["pulsar-birth-spin-period-distribution-max"].as<double>() : pulsarBirthSpinPeriodDistributionMax;          // Maximum birth spin period (ms)
            pulsarBirthSpinPeriodDistributionMin    = vm.count("pulsar-birth-spin-period-distribution-min") ? vm["pulsar-birth-spin-period-distribution-min"].as<double>() : pulsarBirthSpinPeriodDistributionMin;          // Minimum birth spin period (ms)
            pulsarMagneticFieldDecayMassscale       = vm.count("pulsar-magnetic-field-decay-massscale") ? vm["pulsar-magnetic-field-decay-massscale"].as<double>() : pulsarMagneticFieldDecayMassscale;                     // Mass scale on which magnetic field decays during accretion (solar masses)
            pulsarMagneticFieldDecayTimescale       = vm.count("pulsar-magnetic-field-decay-timescale") ? vm["pulsar-magnetic-field-decay-timescale"].as<double>() : pulsarMagneticFieldDecayTimescale;                     // Timescale on which magnetic field decays (Myrs)
            pulsarLog10MinimumMagneticField         = vm.count("pulsar-minimum-magnetic-field") ? vm["pulsar-minimum-magnetic-field"].as<double>() : pulsarLog10MinimumMagneticField;                                       // log10 of the minimum pulsar magnetic field in Gauss

            sampleCommonEnvelopeAlphaMax            = vm.count("sample-common-envelope-alpha-max") ? vm["sample-common-envelope-alpha-max"].as<double>() : sampleCommonEnvelopeAlphaMax;                                    // Maximum for Uniform sampling over common envelope alpha
            sampleCommonEnvelopeAlphaMin            = vm.count("sample-common-envelope-alpha-min") ? vm["sample-common-envelope-alpha-min"].as<double>() : sampleCommonEnvelopeAlphaMin;                                    // Minimum for Uniform sampling over common envelope alpha
            sampleKickDirectionPowerMax             = vm.count("sample-kick-direction-power-max") ? vm["sample-kick-direction-power-max"].as<double>() : sampleKickDirectionPowerMax;                                       // Maximum for Uniform sampling over kick direction powerlaw exponent
            sampleKickDirectionPowerMin             = vm.count("sample-kick-direction-power-min") ? vm["sample-kick-direction-power-min"].as<double>() : sampleKickDirectionPowerMin;                                       // Minimum for Uniform sampling over kick direction powerlaw exponent
            sampleKickVelocitySigmaMax              = vm.count("sample-kick-velocity-sigma-max") ? vm["sample-kick-velocity-sigma-max"].as<double>() : sampleKickVelocitySigmaMax;                                          // Maximum for Uniform sampling over kick velocity sigma
            sampleKickVelocitySigmaMin              = vm.count("sample-kick-velocity-sigma-min") ? vm["sample-kick-velocity-sigma-min"].as<double>() : sampleKickVelocitySigmaMin;                                          // Minimum for Uniform sampling over kick velocity sigma
            sampleLuminousBlueVariableMultiplierMax = vm.count("sample-luminous-blue-variable-multiplier-max") ? vm["sample-luminous-blue-variable-multiplier-max"].as<double>() : sampleLuminousBlueVariableMultiplierMax; // Maximum for Uniform sampling over multiplicative constant for LBV mass loss
            sampleLuminousBlueVariableMultiplierMin = vm.count("sample-luminous-blue-variable-multiplier-min") ? vm["sample-luminous-blue-variable-multiplier-min"].as<double>() : sampleLuminousBlueVariableMultiplierMin; // Minimum for Uniform sampling over multiplicative constant for LBV mass loss
            sampleWolfRayetMultiplierMax            = vm.count("sample-wolf-rayet-multiplier-max") ? vm["sample-wolf-rayet-multiplier-max"].as<double>() : sampleWolfRayetMultiplierMax;                                    // Maximum for Uniform sampling over multiplicative constant for WR winds
            sampleWolfRayetMultiplierMin            = vm.count("sample-wolf-rayet-multiplier-min") ? vm["sample-wolf-rayet-multiplier-min"].as<double>() : sampleWolfRayetMultiplierMin;                                    // Minimum for Uniform sampling over multiplicative constant for WR winds
            semiMajorAxisDistributionMax            = vm.count("semi-major-axis-max") ? vm["semi-major-axis-max"].as<double>() : semiMajorAxisDistributionMax;                                                              // Maximum semi major axis in AU to generate
            semiMajorAxisDistributionMin            = vm.count("semi-major-axis-min") ? vm["semi-major-axis-min"].as<double>() : semiMajorAxisDistributionMin;                                                              // Minimum semi major axis in AU to generate
            singleStarMassMax                       = vm.count("single-star-mass-max") ? vm["single-star-mass-max"].as<double>() : singleStarMassMax;                                                                       // Maximum mass (in Msol) for single star evolution
            singleStarMassMin                       = vm.count("single-star-mass-min") ? vm["single-star-mass-min"].as<double>() : singleStarMassMin;                                                                       // Minimum mass (in Msol) for single star evolution
            spinDistributionMax                     = vm.count("spin-mag-max") ? vm["spin-mag-max"].as<double>() : spinDistributionMax;                                                                                     // Maximum magnitude of spin
            spinDistributionMin                     = vm.count("spin-mag-min") ? vm["spin-mag-min"].as<double>() : spinDistributionMin;                                                                                     // Minimum magnitude of spin

            wolfRayetFactor                         = vm.count("wolf-rayet-multiplier") ? vm["wolf-rayet-multiplier"].as<double>() : wolfRayetFactor;                                                                       // Multiplicitive constant for WR winds

            zetaAdiabaticArbitrary                  = vm.count("zeta-adiabatic-arbitrary") ? vm["zeta-adiabatic-arbitrary"].as<double>() : zetaAdiabaticArbitrary;                                                          // Value of mass-radius exponent zeta adiabatic
            zetaHertzsprungGap                      = vm.count("zeta-hertzsprung-gap") ? vm["zeta-hertzsprung-gap"].as<double>() : zetaHertzsprungGap;                                                                      // Value of mass-radius exponent zeta on the hertzstrpung gap
            zetaMainSequence                        = vm.count("zeta-main-sequence") ? vm["zeta-main-sequence"].as<double>() : zetaMainSequence;                                                                            // Value of mass-radius exponent zeta on the main sequence


            // string options - alphabetically (sort of - where possible - dependencies)

            AISDCOtypeString                        = vm.count("AIS-DCOtype") ? vm["AIS-DCOtype"].as<string>() : AISDCOtypeString;                                                                                          // DCO type selection in exploratory phase of STROOPWAFEL

            blackHoleKicksString                    = vm.count("black-hole-kicks") ? vm["black-hole-kicks"].as<string>() : blackHoleKicksString;                                                                            // Black hole kicks relative to NS kicks

            cheString                               = vm.count("chemically-homogeneous-evolution") ? vm["chemically-homogeneous-evolution"].as<string>() : cheString;                                                       // Chemically Homogeneous Evolution

            commonEnvelopeHertzsprungGapDonorString         = vm.count("common-envelope-hertzsprung-gap-assumption") ? vm["common-envelope-hertzsprung-gap-assumption"].as<string>() : commonEnvelopeHertzsprungGapDonorString;                 // Assumption to make about HG stars in CE
            commonEnvelopeLambdaPrescriptionString          = vm.count("common-envelope-lambda-prescription") ? vm["common-envelope-lambda-prescription"].as<string>() : commonEnvelopeLambdaPrescriptionString;                                // Prescription for CE lambda
            commonEnvelopeMassAccretionPrescriptionString   = vm.count("common-envelope-mass-accretion-prescription") ? vm["common-envelope-mass-accretion-prescription"].as<string>() : commonEnvelopeMassAccretionPrescriptionString;         // Assumption about whether NS/BHs can accrete mass during common envelope evolution
            commonEnvelopeZetaPrescriptionString            = vm.count("common-envelope-zeta-prescription") ? vm["common-envelope-zeta-prescription"].as<string>() : commonEnvelopeZetaPrescriptionString;                                      // Prescription for CE zeta

            eccentricityDistributionString                  = vm.count("eccentricity-distribution,e") ? vm["eccentricity-distribution,e"].as<string>() : eccentricityDistributionString;                                                        // Initial eccentricity distribution

            fryerSupernovaEngineString                      = vm.count("fryer-supernova-engine") ? vm["fryer-supernova-engine"].as<string>() : fryerSupernovaEngineString;                                                                      // If using Fryer et al 2012 fallback prescription

            gridFilename                                    = vm.count("grid") ? vm["grid"].as<string>() : gridFilename;                                                                                                                        // Grid filename

            initialMassFunctionString                       = vm.count("initial-mass-function,i") ? vm["initial-mass-function,i"].as<string>() : initialMassFunctionString;                                                                     // Initial mass function to use

            kickDirectionDistributionString                 = vm.count("kick-direction") ? vm["kick-direction"].as<string>() : kickDirectionDistributionString;                                                                                 // Distribution for natal kick direction
            kickVelocityDistributionString                  = vm.count("kick-velocity-distribution") ? vm["kick-velocity-distribution"].as<string>() : kickVelocityDistributionString;                                                          // Natal kick velocity distribution

            logfileBSEBeBinaries                            = vm.count("logfile-BSE-be-binaries") ? vm["logfile-BSE-be-binaries"].as<string>() : logfileBSEBeBinaries;                                                                          // Filename for BSE Be Binaries logfile
            logfileBSECommonEnvelopes                       = vm.count("logfile-BSE-common-envelopes") ? vm["logfile-BSE-common-envelopes"].as<string>() : logfileBSECommonEnvelopes;                                                           // Filename for BSE Common Envelopes logfile
            logfileBSEDetailedOutput                        = vm.count("logfile-BSE-detailed-output") ? vm["logfile-BSE-detailed-output"].as<string>() : logfileBSEDetailedOutput;                                                              // Filename for BSE Detailed Output logfile
            logfileBSEDoubleCompactObjects                  = vm.count("logfile-BSE-double-compact-objects") ? vm["logfile-BSE-double-compact-objects"].as<string>() : logfileBSEDoubleCompactObjects;                                          // Filename for BSE Double Compact Objects logfile
            logfileBSEPulsarEvolution                       = vm.count("logfile-BSE-pulsar-evolution") ? vm["logfile-BSE-pulsar-evolution"].as<string>() : logfileBSEPulsarEvolution;                                                           // Filename for BSE Pulsar Evolution logfile
            logfileBSESupernovae                            = vm.count("logfile-BSE-supernovae") ? vm["logfile-BSE-supernovae"].as<string>() : logfileBSESupernovae;                                                                            // Filename for BSE Supernovae logfile
            logfileBSESystemParameters                      = vm.count("logfile-BSE-system-parameters") ? vm["logfile-BSE-system-parameters"].as<string>() : logfileBSESystemParameters;                                                        // Filename for BSE System Parameters logfile
            logfileDefinitionsFilename                      = vm.count("logfile-definitions") ? vm["logfile-definitions"].as<string>() : logfileDefinitionsFilename;                                                                            // Filename for logfile record definitions
            logfileDelimiterString                          = vm.count("logfile-delimiter") ? vm["logfile-delimiter"].as<string>() : logfileDelimiterString;                                                                                    // Field delimiter for logfile records
            logfileNamePrefix                               = vm.count("logfile-name-prefix") ? vm["logfile-name-prefix"].as<string>() : logfileNamePrefix;                                                                                     // Prefix for logfile names
            logfileSSEParameters                            = vm.count("logfile-SSE-parameters") ? vm["logfile-SSE-parameters"].as<string>() : logfileSSEParameters;                                                                            // Filename for SSE Parameters logfile

            massLossPrescriptionString                          = vm.count("mass-loss-prescription") ? vm["mass-loss-prescription"].as<string>() : massLossPrescriptionString;                                                                              // Mass loss prescription
            massRatioDistributionString                         = vm.count("mass-ratio-distribution,q") ? vm["mass-ratio-distribution,q"].as<string>() : massRatioDistributionString;                                                                       // Initial mass ratio distribution for q=m2/m1
            massTransferAccretionEfficiencyPrescriptionString   = vm.count("mass-transfer-accretion-efficiency-prescription") ? vm["mass-transfer-accretion-efficiency-prescription"].as<string>() : massTransferAccretionEfficiencyPrescriptionString;     // Mass Transfer Accretion Efficiency prescription
            massTransferAngularMomentumLossPrescriptionString   = vm.count("mass-transfer-angular-momentum-loss-prescription") ? vm["mass-transfer-angular-momentum-loss-prescription"].as<string>() : massTransferAngularMomentumLossPrescriptionString;   // Mass Transfer Angular Momentum Loss prescription
            massTransferPrescriptionString                      = vm.count("mass-transfer-prescription") ? vm["mass-transfer-prescription"].as<string>() : massTransferPrescriptionString;                                                                  // Mass Transfer prescription
            massTransferRejuvenationPrescriptionString          = vm.count("mass-transfer-rejuvenation-prescription") ? vm["mass-transfer-rejuvenation-prescription"].as<string>() : massTransferRejuvenationPrescriptionString;                            // Mass Transfer Rejuvenation prescription
            massTransferThermallyLimitedVariationString         = vm.count("mass-transfer-thermal-limit-accretor") ? vm["mass-transfer-thermal-limit-accretor"].as<string>() : massTransferThermallyLimitedVariationString;                                 // Mass Transfer Thermal Accretion limit

            neutronStarEquationOfStateString                    = vm.count("neutron-star-equation-of-state") ? vm["neutron-star-equation-of-state"].as<string>() : neutronStarEquationOfStateString;                                                        // Neutron star equation of state

            outputPathString                                    = vm.count("outputPath,o") ? vm["outputPath,o"].as<string>() : outputPathString;                                                                                                            // Directory for output

            pulsarBirthMagneticFieldDistributionString          = vm.count("pulsar-birth-magnetic-field-distribution") ? vm["pulsar-birth-magnetic-field-distribution"].as<string>() : pulsarBirthMagneticFieldDistributionString;                          // Distribution of (log10 of) pulsar birth magnetic field in G
            pulsarBirthSpinPeriodDistributionString             = vm.count("pulsar-birth-spin-period-distribution") ? vm["pulsar-birth-spin-period-distribution"].as<string>() : pulsarBirthSpinPeriodDistributionString;                                   // Distribution of pulsar birth spin period in ms
            pulsationalPairInstabilityPrescriptionString        = vm.count("pulsational-pair-instability-prescription") ? vm["pulsational-pair-instability-prescription"].as<string>() : pulsationalPairInstabilityPrescriptionString;                      // Pulsational pair instability distribution

            remnantMassPrescriptionString                       = vm.count("remnant-mass-prescription") ? vm["remnant-mass-prescription"].as<string>() : remnantMassPrescriptionString;                                                                     // Remnant mass prescription
            rotationalVelocityDistributionString                = vm.count("rotational-velocity-distribution") ? vm["rotational-velocity-distribution"].as<string>() : rotationalVelocityDistributionString;                                                // Initial rotational velocity distribution

            semiMajorAxisDistributionString                     = vm.count("semi-major-axis-distribution,a") ? vm["semi-major-axis-distribution,a"].as<string>() : semiMajorAxisDistributionString;                                                         // Initial semi-major axis distribution
            spinAssumptionString                                = vm.count("spin-assumption") ? vm["spin-assumption"].as<string>() : spinAssumptionString;                                                                                                  // Assumption of misalignedments
            spinDistributionString                              = vm.count("spin-distribution") ? vm["spin-distribution"].as<string>() : spinDistributionString;                                                                                            // Distribution of spins

            tidesPrescriptionString                             = vm.count("tides-prescription") ? vm["tides-prescription"].as<string>() : tidesPrescriptionString;                                                                                         // Tides prescription 


            // vector (list) options - alphabetically

            debugClasses                                        = vm.count("debug-classes") ? vm["debug-classes"].as<vector<string>>() : debugClasses;                                                                                                      // Debug classes enabled
            logClasses                                          = vm.count("log-classes") ? vm["log-classes"].as<vector<string>>() : logClasses;                                                                                                            // Logging classes enabled


            // prescriptions, distributions, assumptions etc. options - alphabetically

            bool found;

            if (vm.count("AIS-DCOtype")) {                                                                                              // Adaptive Importance Sampling DCO type
                std::tie(found, AISDCOtype) = utils::GetMapKey(AISDCOtypeString, AIS_DCO_LABEL, AISDCOtype);
                COMPLAIN_IF(!found, "Unknown AIS DCO Type");
            }

            if (vm.count("black-hole-kicks")) {                                                                                         // black hole kicks option
                std::tie(found, blackHoleKicksOption) = utils::GetMapKey(blackHoleKicksString, BLACK_HOLE_KICK_OPTION_LABEL, blackHoleKicksOption);
                COMPLAIN_IF(!found, "Unknown Black Hole Kicks Option");
            }

            if (vm.count("chemically-homogeneous-evolution")) {                                                                         // Chemically Homogeneous Evolution
                std::tie(found, cheOption) = utils::GetMapKey(cheString, CHE_OPTION_LABEL, cheOption);
                COMPLAIN_IF(!found, "Unknown Chemically Homogeneous Evolution Option");
            }

            if (vm.count("common-envelope-hertzsprung-gap-assumption")) {                                                               // common envelope hertzsprung gap assumption
                std::tie(found, commonEnvelopeHertzsprungGapDonor) = utils::GetMapKey(commonEnvelopeHertzsprungGapDonorString, COMMON_ENVELOPE_PRESCRIPTION_LABEL, commonEnvelopeHertzsprungGapDonor);
                COMPLAIN_IF(!found, "Unknown CE HG Assumption");
            }

            if (vm.count("common-envelope-lambda-prescription")) {                                                                      // common envelope lambda prescription
                std::tie(found, commonEnvelopeLambdaPrescription) = utils::GetMapKey(commonEnvelopeLambdaPrescriptionString, CE_LAMBDA_PRESCRIPTION_LABEL, commonEnvelopeLambdaPrescription);
                COMPLAIN_IF(!found, "Unknown CE Lambda Prescription");

                if (commonEnvelopeLambdaPrescription == CE_LAMBDA_PRESCRIPTION::KRUCKOW) {                                              // CE Lambda prescription = Kruckow
                    commonEnvelopeSlopeKruckow = vm.count("common-envelope-slope-Kruckow") ? commonEnvelopeSlopeKruckow : -4.0 / 5.0;   // if user didn't specify choice of the slope, use default
                }
            }

            if (vm.count("common-envelope-mass-accretion-prescription")) {                                                              // common envelope mass accretion prescription
                std::tie(found, commonEnvelopeMassAccretionPrescription) = utils::GetMapKey(commonEnvelopeMassAccretionPrescriptionString, CE_ACCRETION_PRESCRIPTION_LABEL, commonEnvelopeMassAccretionPrescription);
                COMPLAIN_IF(!found, "Unknown CE Mass Accretion Prescription");
            }

            if (vm.count("common-envelope-zeta-prescription")) {                                                                        // common envelope zeta prescription
                std::tie(found, commonEnvelopeZetaPrescription) = utils::GetMapKey(commonEnvelopeZetaPrescriptionString, CE_ZETA_PRESCRIPTION_LABEL, commonEnvelopeZetaPrescription);
                COMPLAIN_IF(!found, "Unknown CE Zeta Prescription");

                if (commonEnvelopeZetaPrescription == CE_ZETA_PRESCRIPTION::ARBITRARY) {                                                // CE Lambda prescription = Kruckow
                    zetaAdiabaticArbitrary = vm.count("zeta-adiabatic-arbitrary") ? zetaAdiabaticArbitrary : 10000.0;                   // if user didn't specify Zetas, asign large value, which will favour stable MT
                    zetaThermalArbitrary   = vm.count("zeta-thermal-arbitrary")   ? zetaThermalArbitrary   : 10000.0;                   // if user didn't specify Zetas, asign large value, which will favour stable MT
                }
            }

            if (vm.count("eccentricity-distribution")) {                                                                                // eccentricity distribution
                std::tie(found, eccentricityDistribution) = utils::GetMapKey(eccentricityDistributionString, ECCENTRICITY_DISTRIBUTION_LABEL, eccentricityDistribution);
                COMPLAIN_IF(!found, "Unknown Eccentricity Distribution");
            }

            if (vm.count("fryer-supernova-engine")) {                                                                                   // Fryer et al. 2012 supernova engine
                std::tie(found, fryerSupernovaEngine) = utils::GetMapKey(fryerSupernovaEngineString, SN_ENGINE_LABEL, fryerSupernovaEngine);
                COMPLAIN_IF(!found, "Unknown Fryer et al. Supernova Engine");
            }

            if (vm.count("initial-mass-function")) {                                                                                    // initial mass function
                std::tie(found, initialMassFunction) = utils::GetMapKey(initialMassFunctionString, INITIAL_MASS_FUNCTION_LABEL, initialMassFunction);
                COMPLAIN_IF(!found, "Unknown Initial Mass Function");
            }

            if (vm.count("kick-direction")) {                                                                                           // kick direction
                std::tie(found, kickDirectionDistribution) = utils::GetMapKey(kickDirectionDistributionString, KICK_DIRECTION_DISTRIBUTION_LABEL, kickDirectionDistribution);
                COMPLAIN_IF(!found, "Unknown Kick Direction Distribution");
            }

            if (vm.count("kick-velocity-distribution")) {                                                                               // kick velocity
                std::tie(found, kickVelocityDistribution) = utils::GetMapKey(kickVelocityDistributionString, KICK_VELOCITY_DISTRIBUTION_LABEL, kickVelocityDistribution);
                COMPLAIN_IF(!found, "Unknown Kick Velocity Distribution");
            }

			if (vm.count("logfile-delimiter")) {                                                                                        // logfile field delimiter
                std::tie(found, logfileDelimiter) = utils::GetMapKey(logfileDelimiterString, DELIMITERLabel, logfileDelimiter);
                COMPLAIN_IF(!found, "Unknown Logfile Delimiter");
            }

            if (vm.count("mass-loss-prescription")) {                                                                                   // mass loss prescription
                std::tie(found, massLossPrescription) = utils::GetMapKey(massLossPrescriptionString, MASS_LOSS_PRESCRIPTION_LABEL, massLossPrescription);
                COMPLAIN_IF(!found, "Unknown Mass Loss Prescription");
            }

            if (vm.count("mass-ratio-distribution")) {                                                                                  // mass ratio distribution
                std::tie(found, massRatioDistribution) = utils::GetMapKey(massRatioDistributionString, MASS_RATIO_DISTRIBUTION_LABEL, massRatioDistribution);
                COMPLAIN_IF(!found, "Unknown Mass Ratio Distribution");
            }

            if (useMassTransfer && vm.count("mass-transfer-prescription")) {                                                            // mass transfer prescription
                std::tie(found, massTransferPrescription) = utils::GetMapKey(massTransferPrescriptionString, MT_PRESCRIPTION_LABEL, massTransferPrescription);
                COMPLAIN_IF(!found, "Unknown Mass Transfer Prescription");
            }

            if (useMassTransfer && vm.count("mass-transfer-accretion-efficiency-prescription")) {                                       // mass transfer accretion efficiency prescription
                std::tie(found, massTransferAccretionEfficiencyPrescription) = utils::GetMapKey(massTransferAccretionEfficiencyPrescriptionString, MT_ACCRETION_EFFICIENCY_PRESCRIPTION_LABEL, massTransferAccretionEfficiencyPrescription);
                COMPLAIN_IF(!found, "Unknown Mass Transfer Angular Momentum Loss Prescription");
            }

            if (useMassTransfer && vm.count("mass-transfer-angular-momentum-loss-prescription")) {                                      // mass transfer angular momentum loss prescription
                std::tie(found, massTransferAngularMomentumLossPrescription) = utils::GetMapKey(massTransferAngularMomentumLossPrescriptionString, MT_ANGULAR_MOMENTUM_LOSS_PRESCRIPTION_LABEL, massTransferAngularMomentumLossPrescription);
                COMPLAIN_IF(!found, "Unknown Mass Transfer Angular Momentum Loss Prescription");
            }

            if (useMassTransfer && vm.count("mass-transfer-rejuvenation-prescription")) {                                               // mass transfer rejuvenation prescription
                std::tie(found, massTransferRejuvenationPrescription) = utils::GetMapKey(massTransferRejuvenationPrescriptionString, MT_REJUVENATION_PRESCRIPTION_LABEL, massTransferRejuvenationPrescription);
                COMPLAIN_IF(!found, "Unknown Mass Transfer Rejuvenation Prescription");
            }

            if (useMassTransfer && vm.count("mass-transfer-thermal-limit-accretor")) {                                                  // mass transfer accretor thermal limit
                std::tie(found, massTransferThermallyLimitedVariation) = utils::GetMapKey(massTransferThermallyLimitedVariationString, MT_THERMALLY_LIMITED_VARIATION_LABEL, massTransferThermallyLimitedVariation);
                COMPLAIN_IF(!found, "Unknown Mass Transfer Accretor Thermal Limit");

                if (massTransferThermallyLimitedVariation == MT_THERMALLY_LIMITED_VARIATION::C_FACTOR) {
                    massTransferCParameter = vm.count("mass-transfer-thermal-limit-C") ? massTransferCParameter : 10.0;                 // if user didn't specify choice of C factor, use default based on choice of thermally limited variation
                }

                if (massTransferThermallyLimitedVariation == MT_THERMALLY_LIMITED_VARIATION::RADIUS_TO_ROCHELOBE) {
                    massTransferCParameter = vm.count("mass-transfer-thermal-limit-C") ? massTransferCParameter : 1.0;                  // if user didn't specify choice of C factor, use default based on choice of thermally limited variation
                }
            }

            if (vm.count("neutrino-mass-loss-bh-formation")) {                                                                          // neutrino mass loss assumption
                std::tie(found, neutrinoMassLossAssumptionBH) = utils::GetMapKey(neutrinoMassLossAssumptionBHString, NEUTRINO_MASS_LOSS_PRESCRIPTION_LABEL, neutrinoMassLossAssumptionBH);
                COMPLAIN_IF(!found, "Unknown Neutrino Mass Loss Assumption");
            }

            if (vm.count("neutron-star-equation-of-state")) {                                                                           // neutron star equation of state
                std::tie(found, neutronStarEquationOfState) = utils::GetMapKey(neutronStarEquationOfStateString, NS_EOSLabel, neutronStarEquationOfState);
                COMPLAIN_IF(!found, "Unknown Neutron Star Equation of State");
            }

            if (vm.count("pulsar-birth-magnetic-field-distribution")) {                                                                 // pulsar birth magnetic field distribution
                std::tie(found, pulsarBirthMagneticFieldDistribution) = utils::GetMapKey(pulsarBirthMagneticFieldDistributionString, PULSAR_BIRTH_MAGNETIC_FIELD_DISTRIBUTION_LABEL, pulsarBirthMagneticFieldDistribution);
                COMPLAIN_IF(!found, "Unknown Pulsar Birth Magnetic Field Distribution");
            }

            if (vm.count("pulsar-birth-spin-period-distribution")) {                                                                    // pulsar birth spin period distribution
                std::tie(found, pulsarBirthSpinPeriodDistribution) = utils::GetMapKey(pulsarBirthSpinPeriodDistributionString, PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION_LABEL, pulsarBirthSpinPeriodDistribution);
                COMPLAIN_IF(!found, "Unknown Pulsar Birth Spin Period Distribution");
            }


			if (vm.count("pulsational-pair-instability-prescription")) {                                                                // pulsational pair instability prescription
                std::tie(found, pulsationalPairInstabilityPrescription) = utils::GetMapKey(pulsationalPairInstabilityPrescriptionString, PPI_PRESCRIPTION_LABEL, pulsationalPairInstabilityPrescription);
                COMPLAIN_IF(!found, "Unknown Pulsational Pair Instability Prescription");
			}

            if (vm.count("remnant-mass-prescription")) {                                                                                // remnant mass prescription
                std::tie(found, remnantMassPrescription) = utils::GetMapKey(remnantMassPrescriptionString, REMNANT_MASS_PRESCRIPTION_LABEL, remnantMassPrescription);
                COMPLAIN_IF(!found, "Unknown Remnant Mass Prescription");
            }

            if (vm.count("rotational-velocity-distribution")) {                                                                         // rotational velocity distribution
                std::tie(found, rotationalVelocityDistribution) = utils::GetMapKey(rotationalVelocityDistributionString, ROTATIONAL_VELOCITY_DISTRIBUTION_LABEL, rotationalVelocityDistribution);
                COMPLAIN_IF(!found, "Unknown Rotational Velocity Distribution");
            }

            if (vm.count("semi-major-axis-distribution")) {                                                                             // semi-major axis distribution
                std::tie(found, semiMajorAxisDistribution) = utils::GetMapKey(semiMajorAxisDistributionString, SEMI_MAJOR_AXIS_DISTRIBUTION_LABEL, semiMajorAxisDistribution);
                COMPLAIN_IF(!found, "Unknown Semi-Major Axis Distribution");
            }

			if (vm.count("spin-assumption")) {                                                                                          // spin assumption
                std::tie(found, spinAssumption) = utils::GetMapKey(spinAssumptionString, SPIN_ASSUMPTION_LABEL, spinAssumption);
                COMPLAIN_IF(!found, "Unknown Spin Assumption");
			}

			if (vm.count("spin-distribution")) {                                                                                        // spin distribution
                std::tie(found, spinDistribution) = utils::GetMapKey(spinDistributionString, SPIN_DISTRIBUTION_LABEL, spinDistribution);
                COMPLAIN_IF(!found, "Unknown Spin Assumption");
			}

            if (vm.count("tides-prescription")) {                                                                                       // tides prescription
                std::tie(found, tidesPrescription) = utils::GetMapKey(tidesPrescriptionString, TIDES_PRESCRIPTION_LABEL, tidesPrescription);
                COMPLAIN_IF(!found, "Unknown Tides Prescription");
            }


            // constraint/value/range checks - alphabetically (where possible)

            COMPLAIN_IF(vm.count("common-envelope-alpha") && commonEnvelopeAlpha < 0.0, "CE alpha (--common-envelope-alpha) < 0");
            COMPLAIN_IF(vm.count("common-envelope-alpha-thermal") && (commonEnvelopeAlphaThermal < 0.0 || commonEnvelopeAlphaThermal > 1.0), "CE alpha thermal (--common-envelope-alpha-thermal) must be between 0 and 1");
            COMPLAIN_IF(vm.count("common-envelope-lambda-multiplier") && commonEnvelopeLambdaMultiplier < 0.0, "CE lambda multiplie (--common-envelope-lambda-multiplier < 0");
            COMPLAIN_IF(vm.count("common-envelope-mass-accretion-constant") && commonEnvelopeMassAccretionConstant < 0.0, "CE mass accretion constant (--common-envelope-mass-accretion-constant) < 0");
            COMPLAIN_IF(vm.count("common-envelope-mass-accretion-max") && commonEnvelopeMassAccretionMax < 0.0, "Maximum accreted mass (--common-envelope-mass-accretion-max) < 0");
            COMPLAIN_IF(vm.count("common-envelope-mass-accretion-min") && commonEnvelopeMassAccretionMin < 0.0, "Minimum accreted mass (--common-envelope-mass-accretion-min) < 0");

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

            if (neutrinoMassLossAssumptionBH == NEUTRINO_MASS_LOSS_PRESCRIPTION::FIXED_MASS) {
               COMPLAIN_IF(neutrinoMassLossValueBH < 0.0, "Neutrino mass loss value < 0");                                              // JR: todo: this is not a user-specified option
            }

            COMPLAIN_IF(logLevel < 0, "Logging level (--log-level) < 0");

            COMPLAIN_IF(nBinaries <= 0, "Number of binaries requested <= 0");

            if (neutrinoMassLossAssumptionBH == NEUTRINO_MASS_LOSS_PRESCRIPTION::FIXED_FRACTION) {
               COMPLAIN_IF(neutrinoMassLossValueBH < 0.0 || neutrinoMassLossValueBH > 1.0, "Neutrino mass loss must be between 0 and 1"); // JR: todo: this is not a user-specified option
            }

            COMPLAIN_IF(massRatioDistributionMin < 0.0 || massRatioDistributionMin > 1.0, "Minimum mass ratio (--mass-ratio-min) must be between 0 and 1");
            COMPLAIN_IF(massRatioDistributionMax < 0.0 || massRatioDistributionMax > 1.0, "Maximum mass ratio (--mass-ratio-max) must be between 0 and 1");
            COMPLAIN_IF(massRatioDistributionMax <= massRatioDistributionMin, "Maximum mass ratio (--mass-ratio-max) must be > Minimum mass ratio (--mass-ratio-min)");

            COMPLAIN_IF(maxEvolutionTime <= 0.0, "Maximum evolution time in Myr (--maxEvolutionTime) must be > 0");

            COMPLAIN_IF(metallicity < 0.0 || metallicity > 1.0, "Metallicity (--metallicity) should be absolute metallicity and must be between 0 and 1");

            COMPLAIN_IF(minimumMassSecondary < 0.0, "Seconday minimum mass (--minimum-secondary-mass) must be >= 0");
            COMPLAIN_IF(minimumMassSecondary > initialMassFunctionMax, "Seconday minimum mass (--minimum-secondary-mass) must be <= Maximum initial mass (--initial-mass-max)");

            if (vm.count("outputPath") or vm.count("o")) {                                                                              // user specified output path?
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

            COMPLAIN_IF(vm.count("pulsar-magnetic-field-decay-timescale") && pulsarMagneticFieldDecayTimescale <= 0.0, "Pulsar magnetic field decay timescale (--pulsar-magnetic-field-decay-timescale) <= 0");
            COMPLAIN_IF(vm.count("pulsar-magnetic-field-decay-massscale") && pulsarMagneticFieldDecayMassscale <= 0.0, "Pulsar Magnetic field decay massscale (--pulsar-magnetic-field-decay-massscale) <= 0");

            COMPLAIN_IF(semiMajorAxisDistributionMin < 0.0, "Minimum semi-major Axis (--semi-major-axis-min) < 0");
            COMPLAIN_IF(semiMajorAxisDistributionMax < 0.0, "Maximum semi-major Axis (--semi-major-axis-max) < 0");

            COMPLAIN_IF(singleStarMassMax   <= 0.0,               "Single star mass maximum (--single-star-mass-max) <= 0");
            COMPLAIN_IF(singleStarMassMax   <= singleStarMassMin, "Single star mass maximum (--single-star-mass-max) <= minimum (--single-star-mass-min)");
            COMPLAIN_IF(singleStarMassMin   <= 0.0,               "Single star mass minimum (--single-star-mass-min) <= 0");
            COMPLAIN_IF(singleStarMassSteps <= 0,                 "Single star mass steps (--single-star-mass-steps) <= 0");

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


