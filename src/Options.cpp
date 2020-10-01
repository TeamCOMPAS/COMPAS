#include "Options.h"
#include "changelog.h"

Options* Options::m_Instance = nullptr;

namespace po  = boost::program_options;
namespace cls = po::command_line_style;
namespace fs  = boost::filesystem;


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


/*
 * Initialise option values
 * 
 * This function initialises the values of the options - this is where the defaults
 * are set.  If a user does not specify an option, either on the commandline or in
 * a grid file, the default values are taken from here.
 * 
 * Note this is a class OptionValues function.
 * 
 * 
 * void Options::OptionValues::Initialise()
 * 
 */
void Options::OptionValues::Initialise() {

    // First set all options to their default values

    // flags

    m_AllowRLOFAtBirth                                              = false;
    m_AllowTouchingAtBirth                                          = false;

    m_DebugToFile                                                   = false;
    m_ErrorsToFile                                                  = false;

    m_EnableWarnings                                                = false;

	m_BeBinaries                                                    = false;
    m_EvolvePulsars                                                 = false;
	m_EvolveUnboundSystems                                          = false;

    m_DetailedOutput                                                = false;
    m_PopulationDataPrinting                                        = false;
    m_PrintBoolAsString                                             = false;
    m_Quiet                                                         = false;
    m_RlofPrinting                                                  = false;
    m_SwitchLog                                                     = false;

    m_nBatchesUsed                                                  = -1;


    // Evolution mode: SSE or BSE
    m_EvolutionMode                                                 = EVOLUTION_MODE::BSE;
    m_EvolutionModeString                                           = EVOLUTION_MODE_LABEL.at(m_EvolutionMode);

    // Population synthesis variables
    m_ObjectsToEvolve                                               = 10;

    m_FixedRandomSeed                                               = false;                                                // TRUE if --random-seed is passed on command line
    m_RandomSeed                                                    = 0;

    // Specify how long to evolve binaries for
    m_MaxEvolutionTime                                              = 13700.0;
    m_MaxNumberOfTimestepIterations                                 = 99999;


    // Initial mass options
    m_InitialMass                                                   = 5.0;
    m_InitialMass1                                                  = 5.0;
    m_InitialMass2                                                  = 5.0;

    m_InitialMassFunction                                           = INITIAL_MASS_FUNCTION::KROUPA;
    m_InitialMassFunctionString                                     = INITIAL_MASS_FUNCTION_LABEL.at(m_InitialMassFunction);
    m_InitialMassFunctionMin                                        = 8.0;
    m_InitialMassFunctionMax                                        = 100.0; 
    m_InitialMassFunctionPower                                      = -2.3;


    // Initial mass ratios
    m_MassRatioDistribution                                         = MASS_RATIO_DISTRIBUTION::FLAT;                        // Most likely want FLAT or SANA2012
    m_MassRatioDistributionString                                   = MASS_RATIO_DISTRIBUTION_LABEL.at(m_MassRatioDistribution);
    m_MassRatioDistributionMin                                      = 0.0;
    m_MassRatioDistributionMax                                      = 1.0;

    m_MinimumMassSecondary                                          = 0.0;


    // Initial orbit options
    m_SemiMajorAxis                                                 = 0.1;
    m_SemiMajorAxisDistribution                                     = SEMI_MAJOR_AXIS_DISTRIBUTION::FLATINLOG;              // Most likely want FLATINLOG or SANA2012
    m_SemiMajorAxisDistributionString                               = SEMI_MAJOR_AXIS_DISTRIBUTION_LABEL.at(SEMI_MAJOR_AXIS_DISTRIBUTION::FLATINLOG);
    m_SemiMajorAxisDistributionMin                                  = 0.1;
    m_SemiMajorAxisDistributionMax                                  = 1000.0;
    m_SemiMajorAxisDistributionPower                                = -1.0; 


    // Initial orbital period
    m_PeriodDistributionMin                                         = 1.1;
    m_PeriodDistributionMax                                         = 1000.0;

    // Eccentricity
    m_Eccentricity                                                  = 0.0;
    m_EccentricityDistribution                                      = ECCENTRICITY_DISTRIBUTION::ZERO; 
    m_EccentricityDistributionString                                = ECCENTRICITY_DISTRIBUTION_LABEL.at(m_EccentricityDistribution);
    m_EccentricityDistributionMin                                   = 0.0;
    m_EccentricityDistributionMax                                   = 1.0;

    // Kick options
    m_KickMagnitudeDistribution                                     = KICK_MAGNITUDE_DISTRIBUTION::MAXWELLIAN;
    m_KickMagnitudeDistributionString                               = KICK_MAGNITUDE_DISTRIBUTION_LABEL.at(m_KickMagnitudeDistribution);
    m_KickMagnitudeDistributionSigmaCCSN_NS                         = 250;
    m_KickMagnitudeDistributionSigmaCCSN_BH                         = 250;
    m_KickMagnitudeDistributionMaximum                              = -1.0; 
    m_KickMagnitudeDistributionSigmaForECSN                         = 30.0;
    m_KickMagnitudeDistributionSigmaForUSSN   	                    = 30.0;
	m_KickScalingFactor						                        = 1.0;

    // Kick direction option
    m_KickDirectionDistribution                                     = KICK_DIRECTION_DISTRIBUTION::ISOTROPIC;
    m_KickDirectionDistributionString                               = KICK_DIRECTION_DISTRIBUTION_LABEL.at(m_KickDirectionDistribution);
    m_KickDirectionPower                                            = 0.0;

    // Kick magnitude
    m_KickMagnitude                                                 = 0.0;
    m_KickMagnitude1                                                = 0.0;
    m_KickMagnitude2                                                = 0.0;                               

    // Kick magnitude random number (used to draw kick magnitude if necessary)
    m_KickMagnitudeRandom                                           = RAND->Random();
    m_KickMagnitudeRandom1                                          = RAND->Random();
    m_KickMagnitudeRandom2                                          = RAND->Random();

    // Mean anomaly
    m_KickMeanAnomaly1                                              = RAND->Random(0.0, _2_PI);
    m_KickMeanAnomaly2                                              = RAND->Random(0.0, _2_PI);

    // Phi
    m_KickPhi1                                                      = 0.0;                                              // actual value set later
    m_KickPhi2                                                      = 0.0;                                              // actual value set later

    // Theta
    m_KickTheta1                                                    = 0.0;                                              // actual value set later 
    m_KickTheta2                                                    = 0.0;                                              // actual value set later

    // Black hole kicks
    m_BlackHoleKicksOption                                          = BLACK_HOLE_KICK_OPTION::FALLBACK;
    m_BlackHoleKicksOptionString                                    = BLACK_HOLE_KICK_OPTION_LABEL.at(m_BlackHoleKicksOption);


    // Chemically Homogeneous Evolution
    m_CheOption                                                     = CHE_OPTION::NONE;
    m_CheString                                                     = CHE_OPTION_LABEL.at(m_CheOption);


    // Supernova remnant mass prescription options
    m_RemnantMassPrescription                                       = REMNANT_MASS_PRESCRIPTION::FRYER2012;
    m_RemnantMassPrescriptionString                                 = REMNANT_MASS_PRESCRIPTION_LABEL.at(m_RemnantMassPrescription);

    m_FryerSupernovaEngine                                          = SN_ENGINE::DELAYED;
    m_FryerSupernovaEngineString                                    = SN_ENGINE_LABEL.at(m_FryerSupernovaEngine);

    m_NeutrinoMassLossAssumptionBH                                  = NEUTRINO_MASS_LOSS_PRESCRIPTION::FIXED_FRACTION;
    m_NeutrinoMassLossAssumptionBHString                            = NEUTRINO_MASS_LOSS_PRESCRIPTION_LABEL.at(m_NeutrinoMassLossAssumptionBH);
    m_NeutrinoMassLossValueBH                                       = 0.1;


    // Fixed uk options
    m_UseFixedUK                                                    = false;
    m_FixedUK                                                       = -1.0;


    // Pair instability and pulsational pair instability mass loss
    m_UsePairInstabilitySupernovae                                  = false;
    m_PairInstabilityLowerLimit                                     = 60.0;                                                 // Belczynski+ 2016 is 65 Msol
    m_PairInstabilityUpperLimit                                     = 135.0;                                                // Belczynski+ 2016 is 135 Msol

    m_UsePulsationalPairInstability                                 = false;
    m_PulsationalPairInstabilityLowerLimit                          = 35.0;                                                 // Belczynski+ 2016 is 45 Msol
    m_PulsationalPairInstabilityUpperLimit                          = 60.0;                                                 // Belczynski+ 2016 is 65 Msol

    m_PulsationalPairInstabilityPrescription                        = PPI_PRESCRIPTION::COMPAS;
    m_PulsationalPairInstabilityPrescriptionString                  = PPI_PRESCRIPTION_LABEL.at(m_PulsationalPairInstabilityPrescription);

	m_MaximumNeutronStarMass                                        = 3.0;                                                  // StarTrack is 3.0
    
    m_mCBUR1                                                        = MCBUR1HURLEY;                                         // MHurley value, Fryer+ and Belczynski+ use 1.83


    // Output path
    m_OutputPathString                                              = ".";
    m_DefaultOutputPath                                             = boost::filesystem::current_path();
    m_OutputPath                                                    = m_DefaultOutputPath;
    m_OutputContainerName                                           = DEFAULT_OUTPUT_CONTAINER_NAME;
    

    // Mass loss options
    m_UseMassLoss                                                   = false;

    m_MassLossPrescription                                          = MASS_LOSS_PRESCRIPTION::VINK;
    m_MassLossPrescriptionString                                    = MASS_LOSS_PRESCRIPTION_LABEL.at(m_MassLossPrescription);


    // Wind mass loss multiplicitive constants
    m_LuminousBlueVariableFactor                                    = 1.5;
    m_WolfRayetFactor                                               = 1.0;


    // Mass transfer options
    m_UseMassTransfer                                               = true;
	m_CirculariseBinaryDuringMassTransfer         	                = false;
	m_AngularMomentumConservationDuringCircularisation              = false;

    // Case BB/BC mass transfer stability prescription
    m_CaseBBStabilityPrescription                                   = CASE_BB_STABILITY_PRESCRIPTION::ALWAYS_STABLE;
    m_CaseBBStabilityPrescriptionString                             = CASE_BB_STABILITY_PRESCRIPTION_LABEL.at(m_CaseBBStabilityPrescription);

    // Options adaptive Roche Lobe Overflow prescription
    m_MassTransferAdaptiveAlphaParameter                            = 0.5;
    m_MaxPercentageAdaptiveMassTransfer                             = 0.01;

    // Options for mass transfer accretion efficiency
    m_MassTransferAccretionEfficiencyPrescription                   = MT_ACCRETION_EFFICIENCY_PRESCRIPTION::THERMALLY_LIMITED;
    m_MassTransferAccretionEfficiencyPrescriptionString             = MT_ACCRETION_EFFICIENCY_PRESCRIPTION_LABEL.at(m_MassTransferAccretionEfficiencyPrescription);

    m_MassTransferFractionAccreted                                  = 1.0;
    m_MassTransferCParameter                                        = 10.0;
    m_EddingtonAccretionFactor                                      = 1;                                                    // >1 is super-eddington, 0 is no accretion

    // Mass transfer thermally limited options
	m_MassTransferThermallyLimitedVariation                         = MT_THERMALLY_LIMITED_VARIATION::C_FACTOR;
	m_MassTransferThermallyLimitedVariationString                   = MT_THERMALLY_LIMITED_VARIATION_LABEL.at(m_MassTransferThermallyLimitedVariation);

    // Mass transfer angular momentum loss prescription options
    m_MassTransferJloss                                             = 1.0;
    m_MassTransferAngularMomentumLossPrescription                   = MT_ANGULAR_MOMENTUM_LOSS_PRESCRIPTION::ISOTROPIC_RE_EMISSION;
    m_MassTransferAngularMomentumLossPrescriptionString             = MT_ANGULAR_MOMENTUM_LOSS_PRESCRIPTION_LABEL.at(m_MassTransferAngularMomentumLossPrescription);

    // Mass transfer rejuvenation prescriptions
    m_MassTransferRejuvenationPrescription                          = MT_REJUVENATION_PRESCRIPTION::NONE;
    m_MassTransferRejuvenationPrescriptionString                    = MT_REJUVENATION_PRESCRIPTION_LABEL.at(m_MassTransferRejuvenationPrescription);

    // Mass transfer critical mass ratios
    m_MassTransferCriticalMassRatioMSLowMass                        = false;
    m_MassTransferCriticalMassRatioMSLowMassNonDegenerateAccretor   = 1.44;                                                 // Claeys+ 2014 = 1.44
    m_MassTransferCriticalMassRatioMSLowMassDegenerateAccretor      = 1.0;                                                  // Claeys+ 2014 = 1.0

    m_MassTransferCriticalMassRatioMSHighMass                       = false;
    m_MassTransferCriticalMassRatioMSHighMassNonDegenerateAccretor  = 0.625;                                                // Claeys+ 2014 = 0.625
    m_MassTransferCriticalMassRatioMSHighMassDegenerateAccretor     = 0.0;

    m_MassTransferCriticalMassRatioHG                               = false;
    m_MassTransferCriticalMassRatioHGNonDegenerateAccretor          = 0.40;                                                 // Claeys+ 2014 = 0.25
    m_MassTransferCriticalMassRatioHGDegenerateAccretor             = 0.21;                                                 // Claeys+ 2014 = 0.21

    m_MassTransferCriticalMassRatioGiant                            = false;
    m_MassTransferCriticalMassRatioGiantNonDegenerateAccretor       = 0.0;
    m_MassTransferCriticalMassRatioGiantDegenerateAccretor          = 0.87;                                                 // Claeys+ 2014 = 0.81

    m_MassTransferCriticalMassRatioHeliumMS                         = false;
    m_MassTransferCriticalMassRatioHeliumMSNonDegenerateAccretor    = 0.625;
    m_MassTransferCriticalMassRatioHeliumMSDegenerateAccretor       = 0.0;

    m_MassTransferCriticalMassRatioHeliumHG                         = false;
    m_MassTransferCriticalMassRatioHeliumHGNonDegenerateAccretor    = 0.25;                                                 // Claeys+ 2014 = 0.25
    m_MassTransferCriticalMassRatioHeliumHGDegenerateAccretor       = 0.21;                                                 // Claeys+ 2014 = 0.21

    m_MassTransferCriticalMassRatioHeliumGiant                      = false;
    m_MassTransferCriticalMassRatioHeliumGiantNonDegenerateAccretor = 1.28;                                                 // Claeys+ 2014 = 0.25
    m_MassTransferCriticalMassRatioHeliumGiantDegenerateAccretor    = 0.87;

    m_MassTransferCriticalMassRatioWhiteDwarf                       = false;
	m_MassTransferCriticalMassRatioWhiteDwarfNonDegenerateAccretor  = 0.0;
    m_MassTransferCriticalMassRatioWhiteDwarfDegenerateAccretor     = 1.6;                                                  // Claeys+ 2014 = 1.6


    // Common Envelope options
    m_CommonEnvelopeAlpha                                           = 1.0;
    m_CommonEnvelopeLambda                                          = 0.1;
	m_CommonEnvelopeSlopeKruckow                                    = -4.0 / 5.0;
	m_CommonEnvelopeAlphaThermal                                    = 1.0;
    m_CommonEnvelopeLambdaMultiplier                                = 1.0;
    m_AllowMainSequenceStarToSurviveCommonEnvelope                  = false;

    // Prescription for envelope state (radiative or convective)
    m_EnvelopeStatePrescription                                     = ENVELOPE_STATE_PRESCRIPTION::LEGACY;
    m_EnvelopeStatePrescriptionString                               = ENVELOPE_STATE_PRESCRIPTION_LABEL.at(m_EnvelopeStatePrescription);

    // Accretion during common envelope
    m_CommonEnvelopeMassAccretionPrescription                       = CE_ACCRETION_PRESCRIPTION::ZERO;
    m_CommonEnvelopeMassAccretionPrescriptionString                 = CE_ACCRETION_PRESCRIPTION_LABEL.at(m_CommonEnvelopeMassAccretionPrescription);
    
    m_CommonEnvelopeMassAccretionMin                                = 0.04;
    m_CommonEnvelopeMassAccretionMax                                = 0.1;
    m_CommonEnvelopeMassAccretionConstant                           = 0.0;

	// Common envelope lambda prescription
	m_CommonEnvelopeLambdaPrescription                              = CE_LAMBDA_PRESCRIPTION::NANJING;
	m_CommonEnvelopeLambdaPrescriptionString                        = CE_LAMBDA_PRESCRIPTION_LABEL.at(m_CommonEnvelopeLambdaPrescription);

	// Common envelope Nandez and Ivanova energy formalism
	m_RevisedEnergyFormalismNandezIvanova	                        = false;
	m_MaximumMassDonorNandezIvanova                                 = 2.0;
	m_CommonEnvelopeRecombinationEnergyDensity                      = 1.5E13;


    // Adaptive Importance Sampling options
    m_AISexploratoryPhase                                           = false;
    m_AISDCOtype                                                    = AIS_DCO::ALL;
    m_AISDCOtypeString                                              = AIS_DCO_LABEL.at(AIS_DCO::ALL);
    m_AIShubble                                                     = false;
    m_AISpessimistic                                                = false;
    m_AISrefinementPhase                                            = false;
    m_AISrlof                                                       = false;
    m_KappaGaussians                                                = 2;


	// Zetas
	m_StellarZetaPrescription                                       = ZETA_PRESCRIPTION::SOBERMAN;
	m_StellarZetaPrescriptionString                                 = ZETA_PRESCRIPTION_LABEL.at(m_StellarZetaPrescription);
	m_ZetaAdiabaticArbitrary                                        = 10000.0;                                              // large value favours stable MT
    m_ZetaMainSequence 	                                            = 2.0;
	m_ZetaRadiativeEnvelopeGiant	                                = 6.5;


    // Metallicity options
    m_Metallicity                                                   = ZSOL;
    m_FixedMetallicity                                              = true;


    // Neutron star equation of state
    m_NeutronStarEquationOfState                                    = NS_EOS::SSE;
    m_NeutronStarEquationOfStateString                              = NS_EOSLabel.at(NS_EOS::SSE);


    // Pulsar birth magnetic field distribution
    m_PulsarBirthMagneticFieldDistribution                          = PULSAR_BIRTH_MAGNETIC_FIELD_DISTRIBUTION::ZERO;
    m_PulsarBirthMagneticFieldDistributionString                    = PULSAR_BIRTH_MAGNETIC_FIELD_DISTRIBUTION_LABEL.at(m_PulsarBirthMagneticFieldDistribution);
    m_PulsarBirthMagneticFieldDistributionMin                       = 11.0;
    m_PulsarBirthMagneticFieldDistributionMax                       = 13.0;


    // Pulsar birth spin period distribution string
    m_PulsarBirthSpinPeriodDistribution                             = PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION::ZERO;
    m_PulsarBirthSpinPeriodDistributionString                       = PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION_LABEL.at(m_PulsarBirthSpinPeriodDistribution);
    m_PulsarBirthSpinPeriodDistributionMin                          = 0.0;
    m_PulsarBirthSpinPeriodDistributionMax                          = 100.0;

    m_PulsarMagneticFieldDecayTimescale                             = 1000.0;
    m_PulsarMagneticFieldDecayMassscale                             = 0.025;
    m_PulsarLog10MinimumMagneticField                               = 8.0;


    // Rotational velocity distribution options
    m_RotationalVelocityDistribution                                = ROTATIONAL_VELOCITY_DISTRIBUTION::ZERO;
    m_RotationalVelocityDistributionString                          = ROTATIONAL_VELOCITY_DISTRIBUTION_LABEL.at(m_RotationalVelocityDistribution);


	// grids

	m_GridFilename                                                  = "";


    // debug and logging options

    m_DebugLevel                                                    = 0;
    m_DebugClasses.clear();

    m_LogLevel                                                      = 0;
    m_LogClasses.clear();


    // Logfiles    
    m_LogfileDefinitionsFilename                                    = "";
    m_LogfileDelimiter                                              = DELIMITER::TAB;
    m_LogfileDelimiterString                                        = DELIMITERLabel.at(m_LogfileDelimiter);
    m_LogfileNamePrefix                                             = "";

    m_LogfileSystemParameters                                       = get<0>(LOGFILE_DESCRIPTOR.at(LOGFILE::SYSTEM_PARAMETERS));
    m_LogfileDetailedOutput                                         = (m_EvolutionMode == EVOLUTION_MODE::SSE) ? get<0>(LOGFILE_DESCRIPTOR.at(LOGFILE::SSE_DETAILED_OUTPUT)) : get<0>(LOGFILE_DESCRIPTOR.at(LOGFILE::BSE_DETAILED_OUTPUT));
    m_LogfileDoubleCompactObjects                                   = get<0>(LOGFILE_DESCRIPTOR.at(LOGFILE::DOUBLE_COMPACT_OBJECTS));
    m_LogfileSupernovae                                             = (m_EvolutionMode == EVOLUTION_MODE::SSE) ? get<0>(LOGFILE_DESCRIPTOR.at(LOGFILE::SSE_SUPERNOVAE)) : get<0>(LOGFILE_DESCRIPTOR.at(LOGFILE::BSE_SUPERNOVAE));
    m_LogfileCommonEnvelopes                                        = get<0>(LOGFILE_DESCRIPTOR.at(LOGFILE::COMMON_ENVELOPES));
    m_LogfileRLOFParameters                                         = get<0>(LOGFILE_DESCRIPTOR.at(LOGFILE::RLOF_PARAMETERS));
    m_LogfileBeBinaries                                             = get<0>(LOGFILE_DESCRIPTOR.at(LOGFILE::BE_BINARIES));
    m_LogfilePulsarEvolution                                        = get<0>(LOGFILE_DESCRIPTOR.at(LOGFILE::BSE_PULSAR_EVOLUTION)); // only BSE for now
    m_LogfileSwitchLog                                              = (m_EvolutionMode == EVOLUTION_MODE::SSE) ? get<0>(LOGFILE_DESCRIPTOR.at(LOGFILE::SSE_SWITCH_LOG)) : get<0>(LOGFILE_DESCRIPTOR.at(LOGFILE::BSE_SWITCH_LOG));


    
}
    

/*
 * Sanity check options and option values
 * 
 * We can currently sample mass, metallicity, separation, eccentricity etc. within COMPAS,
 * so those options don't need to be mandatory - but when we move all sampling out of 
 * COMPAS we will need to enforce those as mandatory options - unless we decide to leave
 * some minimal sampling inside COMPAS to allow for missing options.  We would only need
 * to leave a single distribution for each - we wouldn't want to give the user the option
 * of choosing a distribution - the functionality would only be for convenience if an
 * option was missing.
 * 
 * The boost variable map from the parsed options should already have been set before calling
 * this function.  This records, for each option, whether the user specified a value and, if 
 * so, the value specified by the user.  This function sanity checks the user specified values, 
 * sets the values if all pass the sanity checks, then sets the values of the options not 
 * specified by the user the the defaults specifed here.
 * 
 * Note this is a class OptionValues function.
 * 
 * 
 * std::string Options::OptionValues::CheckAndSetOptions(const po::variables_map p_VM)
 * 
 * @return                                      String containing an error string
 *                                              If no error occurred the return string will be the empty string 
 */
std::string Options::OptionValues::CheckAndSetOptions() {

    std::string errStr = "";                                                                                // error string

    // check & set prescriptions, distributions, assumptions etc. options - alphabetically

    try {

        bool found;

        m_FixedRandomSeed  = !m_VM["random-seed"].defaulted();                                              // use random seed if it is provided by the user
        m_FixedMetallicity = !m_VM["metallicity"].defaulted();                                              // determine if user supplied a metallicity value
        m_UseFixedUK       = !m_VM["fix-dimensionless-kick-magnitude"].defaulted() && (m_FixedUK >= 0.0);   // determine if user supplied a valid kick magnitude


        // Floor
        /*
        if (!vmCmdLine["AIS-DCOtype"].defaulted()) {                                                        // Adaptive Importance Sampling DCO type
            std::tie(found, p_OptionValues->m_AISDCOtype) = utils::GetMapKey(p_OptionValues->m_AISDCOtypeString, AIS_DCO_LABEL, p_OptionValues->m_AISDCOtype);
            return "Unknown AIS DCO Type";
        }
        */

        if (!m_VM["black-hole-kicks"].defaulted()) {                                                        // black hole kicks option
            std::tie(found, m_BlackHoleKicksOption) = utils::GetMapKey(m_BlackHoleKicksOptionString, BLACK_HOLE_KICK_OPTION_LABEL, m_BlackHoleKicksOption);
            if (!found) return "Unknown Black Hole Kicks Option";
        }

        if (!m_VM["case-bb-stability-prescription"].defaulted()) {                                          //case BB/BC mass transfer stability prescription
            std::tie(found, m_CaseBBStabilityPrescription) = utils::GetMapKey(m_CaseBBStabilityPrescriptionString, CASE_BB_STABILITY_PRESCRIPTION_LABEL, m_CaseBBStabilityPrescription);
            if (!found) return "Unknown Case BB/BC Mass Transfer Stability Prescription";
        }
           
        if (!m_VM["chemically-homogeneous-evolution"].defaulted()) {                                        // Chemically Homogeneous Evolution
            std::tie(found, m_CheOption) = utils::GetMapKey(m_CheString, CHE_OPTION_LABEL, m_CheOption);
            if (!found) return "Unknown Chemically Homogeneous Evolution Option";
        }

        if (!m_VM["common-envelope-lambda-prescription"].defaulted()) {                                     // common envelope lambda prescription
            std::tie(found, m_CommonEnvelopeLambdaPrescription) = utils::GetMapKey(m_CommonEnvelopeLambdaPrescriptionString, CE_LAMBDA_PRESCRIPTION_LABEL, m_CommonEnvelopeLambdaPrescription);
            if (!found) return "Unknown CE Lambda Prescription";
        }

        if (!m_VM["common-envelope-mass-accretion-prescription"].defaulted()) {                             // common envelope mass accretion prescription
            std::tie(found, m_CommonEnvelopeMassAccretionPrescription) = utils::GetMapKey(m_CommonEnvelopeMassAccretionPrescriptionString, CE_ACCRETION_PRESCRIPTION_LABEL, m_CommonEnvelopeMassAccretionPrescription);
            if (!found) return "Unknown CE Mass Accretion Prescription";
        }
            
        if (!m_VM["envelope-state-prescription"].defaulted()) {                                             // envelope state prescription
            std::tie(found, m_EnvelopeStatePrescription) = utils::GetMapKey(m_EnvelopeStatePrescriptionString, ENVELOPE_STATE_PRESCRIPTION_LABEL, m_EnvelopeStatePrescription);
            if (!found) return "Unknown Envelope State Prescription";
        }

        if (!m_VM["eccentricity-distribution"].defaulted()) {                                               // eccentricity distribution
            std::tie(found, m_EccentricityDistribution) = utils::GetMapKey(m_EccentricityDistributionString, ECCENTRICITY_DISTRIBUTION_LABEL, m_EccentricityDistribution);
            if (!found) return "Unknown Eccentricity Distribution";
        }

        if (!m_VM["fryer-supernova-engine"].defaulted()) {                                                  // Fryer et al. 2012 supernova engine
            std::tie(found, m_FryerSupernovaEngine) = utils::GetMapKey(m_FryerSupernovaEngineString, SN_ENGINE_LABEL, m_FryerSupernovaEngine);
            if (!found) return "Unknown Fryer et al. Supernova Engine";
        }

        if (!m_VM["initial-mass-function"].defaulted()) {                                                   // initial mass function
            std::tie(found, m_InitialMassFunction) = utils::GetMapKey(m_InitialMassFunctionString, INITIAL_MASS_FUNCTION_LABEL, m_InitialMassFunction);
            if (!found) return "Unknown Initial Mass Function";
        }

        if (!m_VM["kick-direction"].defaulted()) {                                                          // kick direction
            std::tie(found, m_KickDirectionDistribution) = utils::GetMapKey(m_KickDirectionDistributionString, KICK_DIRECTION_DISTRIBUTION_LABEL, m_KickDirectionDistribution);
            if (!found) return "Unknown Kick Direction Distribution";
        }

        if (!m_VM["kick-magnitude-distribution"].defaulted()) {                                             // kick magnitude
            std::tie(found, m_KickMagnitudeDistribution) = utils::GetMapKey(m_KickMagnitudeDistributionString, KICK_MAGNITUDE_DISTRIBUTION_LABEL, m_KickMagnitudeDistribution);
            if (!found) return "Unknown Kick Magnitude Distribution";
        }

        // set values for m_KickPhi[1/2] and m_KickTheta[1/2] here
        // we now have the kick direction distribution and kick direction power (exponent) required by the user (either default or specified)

        bool phi1Defaulted   = m_VM["kick-phi-1"].defaulted();
        bool theta1Defaulted = m_VM["kick-theta-1"].defaulted();

        if (phi1Defaulted || theta1Defaulted) {
            double phi1, theta1;
            std::tie(phi1, theta1) = utils::DrawKickDirection(m_KickDirectionDistribution, m_KickDirectionPower);
            if (phi1Defaulted  ) m_KickPhi1   = phi1;
            if (theta1Defaulted) m_KickTheta1 = theta1;
        }

        bool phi2Defaulted   = m_VM["kick-phi-2"].defaulted();
        bool theta2Defaulted = m_VM["kick-theta-2"].defaulted();

        if (phi2Defaulted || theta2Defaulted) {
            double phi2, theta2;
            std::tie(phi2, theta2) = utils::DrawKickDirection(m_KickDirectionDistribution, m_KickDirectionPower);
            if (phi2Defaulted  ) m_KickPhi2   = phi2;
            if (theta2Defaulted) m_KickTheta2 = theta2;
        }

        if (!m_VM["logfile-delimiter"].defaulted()) {                                                       // logfile field delimiter
            std::tie(found, m_LogfileDelimiter) = utils::GetMapKey(m_LogfileDelimiterString, DELIMITERLabel, m_LogfileDelimiter);
            if (!found) return "Unknown Logfile Delimiter";
        }

        if (!m_VM["mass-loss-prescription"].defaulted()) {                                                  // mass loss prescription
            std::tie(found, m_MassLossPrescription) = utils::GetMapKey(m_MassLossPrescriptionString, MASS_LOSS_PRESCRIPTION_LABEL, m_MassLossPrescription);
            if (!found) return "Unknown Mass Loss Prescription";
        }

        if (!m_VM["mass-ratio-distribution"].defaulted()) {                                                 // mass ratio distribution
            std::tie(found, m_MassRatioDistribution) = utils::GetMapKey(m_MassRatioDistributionString, MASS_RATIO_DISTRIBUTION_LABEL, m_MassRatioDistribution);
            if (!found) return "Unknown Mass Ratio Distribution";
        }

        if (m_UseMassTransfer && !m_VM["mass-transfer-accretion-efficiency-prescription"].defaulted()) {    // mass transfer accretion efficiency prescription
            std::tie(found, m_MassTransferAccretionEfficiencyPrescription) = utils::GetMapKey(m_MassTransferAccretionEfficiencyPrescriptionString, MT_ACCRETION_EFFICIENCY_PRESCRIPTION_LABEL, m_MassTransferAccretionEfficiencyPrescription);
            if (!found) return "Unknown Mass Transfer Angular Momentum Loss Prescription";
        }

        if (m_UseMassTransfer && !m_VM["mass-transfer-angular-momentum-loss-prescription"].defaulted()) {   // mass transfer angular momentum loss prescription
            std::tie(found, m_MassTransferAngularMomentumLossPrescription) = utils::GetMapKey(m_MassTransferAngularMomentumLossPrescriptionString, MT_ANGULAR_MOMENTUM_LOSS_PRESCRIPTION_LABEL, m_MassTransferAngularMomentumLossPrescription);
            if (!found) return "Unknown Mass Transfer Angular Momentum Loss Prescription";
        }

        if (m_UseMassTransfer && !m_VM["mass-transfer-rejuvenation-prescription"].defaulted()) {            // mass transfer rejuvenation prescription
            std::tie(found, m_MassTransferRejuvenationPrescription) = utils::GetMapKey(m_MassTransferRejuvenationPrescriptionString, MT_REJUVENATION_PRESCRIPTION_LABEL, m_MassTransferRejuvenationPrescription);
            if (!found) return "Unknown Mass Transfer Rejuvenation Prescription";
        }

        if (m_UseMassTransfer && !m_VM["mass-transfer-thermal-limit-accretor"].defaulted()) {               // mass transfer accretor thermal limit
            std::tie(found, m_MassTransferThermallyLimitedVariation) = utils::GetMapKey(m_MassTransferThermallyLimitedVariationString, MT_THERMALLY_LIMITED_VARIATION_LABEL, m_MassTransferThermallyLimitedVariation);
            if (!found) return "Unknown Mass Transfer Accretor Thermal Limit";

            if (m_MassTransferThermallyLimitedVariation == MT_THERMALLY_LIMITED_VARIATION::C_FACTOR) {
                m_MassTransferCParameter = m_VM["mass-transfer-thermal-limit-C"].defaulted() ? 10.0 : m_MassTransferCParameter;
            }

            if (m_MassTransferThermallyLimitedVariation == MT_THERMALLY_LIMITED_VARIATION::RADIUS_TO_ROCHELOBE) {
                m_MassTransferCParameter = m_VM["mass-transfer-thermal-limit-C"].defaulted() ? 1.0 : m_MassTransferCParameter;
            }
        }

        if (!m_VM["mode"].defaulted()) {                                                                    // mode
            std::tie(found, m_EvolutionMode) = utils::GetMapKey(m_EvolutionModeString, EVOLUTION_MODE_LABEL, m_EvolutionMode);
            if (!found) return "Unknown Mode";
        }

        if (!m_VM["neutrino-mass-loss-bh-formation"].defaulted()) {                                         // neutrino mass loss assumption
            std::tie(found, m_NeutrinoMassLossAssumptionBH) = utils::GetMapKey(m_NeutrinoMassLossAssumptionBHString, NEUTRINO_MASS_LOSS_PRESCRIPTION_LABEL, m_NeutrinoMassLossAssumptionBH);
            if (!found) return "Unknown Neutrino Mass Loss Assumption";
        }

        if (!m_VM["neutron-star-equation-of-state"].defaulted()) {                                          // neutron star equation of state
            std::tie(found, m_NeutronStarEquationOfState) = utils::GetMapKey(m_NeutronStarEquationOfStateString, NS_EOSLabel, m_NeutronStarEquationOfState);
            if (!found) return "Unknown Neutron Star Equation of State";
        }

        if (!m_VM["pulsar-birth-magnetic-field-distribution"].defaulted()) {                                // pulsar birth magnetic field distribution
            std::tie(found, m_PulsarBirthMagneticFieldDistribution) = utils::GetMapKey(m_PulsarBirthMagneticFieldDistributionString, PULSAR_BIRTH_MAGNETIC_FIELD_DISTRIBUTION_LABEL, m_PulsarBirthMagneticFieldDistribution);
            if (!found) return "Unknown Pulsar Birth Magnetic Field Distribution";
        }

        if (!m_VM["pulsar-birth-spin-period-distribution"].defaulted()) {                                   // pulsar birth spin period distribution
            std::tie(found, m_PulsarBirthSpinPeriodDistribution) = utils::GetMapKey(m_PulsarBirthSpinPeriodDistributionString, PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION_LABEL, m_PulsarBirthSpinPeriodDistribution);
            if (!found) return "Unknown Pulsar Birth Spin Period Distribution";
        }

        if (!m_VM["pulsational-pair-instability-prescription"].defaulted()) {                               // pulsational pair instability prescription
            std::tie(found, m_PulsationalPairInstabilityPrescription) = utils::GetMapKey(m_PulsationalPairInstabilityPrescriptionString, PPI_PRESCRIPTION_LABEL, m_PulsationalPairInstabilityPrescription);
            if (!found) return "Unknown Pulsational Pair Instability Prescription";
        }

        if (!m_VM["remnant-mass-prescription"].defaulted()) {                                               // remnant mass prescription
            std::tie(found, m_RemnantMassPrescription) = utils::GetMapKey(m_RemnantMassPrescriptionString, REMNANT_MASS_PRESCRIPTION_LABEL, m_RemnantMassPrescription);
            if (!found) return "Unknown Remnant Mass Prescription";
        }

        if (!m_VM["rotational-velocity-distribution"].defaulted()) {                                        // rotational velocity distribution
            std::tie(found, m_RotationalVelocityDistribution) = utils::GetMapKey(m_RotationalVelocityDistributionString, ROTATIONAL_VELOCITY_DISTRIBUTION_LABEL, m_RotationalVelocityDistribution);
            if (!found) return "Unknown Rotational Velocity Distribution";
        }

        if (!m_VM["semi-major-axis-distribution"].defaulted()) {                                            // semi-major axis distribution
            std::tie(found, m_SemiMajorAxisDistribution) = utils::GetMapKey(m_SemiMajorAxisDistributionString, SEMI_MAJOR_AXIS_DISTRIBUTION_LABEL, m_SemiMajorAxisDistribution);
            if (!found) return "Unknown Semi-Major Axis Distribution";
        }

        if (!m_VM["stellar-zeta-prescription"].defaulted()) {                                               // common envelope zeta prescription
            std::tie(found, m_StellarZetaPrescription) = utils::GetMapKey(m_StellarZetaPrescriptionString, ZETA_PRESCRIPTION_LABEL, m_StellarZetaPrescription);
            if (!found) return "Unknown stellar Zeta Prescription";
        }

        // constraint/value/range checks - alphabetically (where possible)

        if (!m_VM["common-envelope-alpha"].defaulted() && m_CommonEnvelopeAlpha < 0.0) return "CE alpha (--common-envelope-alpha) < 0";
        if (!m_VM["common-envelope-alpha-thermal"].defaulted() && (m_CommonEnvelopeAlphaThermal < 0.0 || m_CommonEnvelopeAlphaThermal > 1.0)) return "CE alpha thermal (--common-envelope-alpha-thermal) must be between 0 and 1";
        if (!m_VM["common-envelope-lambda-multiplier"].defaulted() && m_CommonEnvelopeLambdaMultiplier < 0.0) return "CE lambda multiplie (--common-envelope-lambda-multiplier < 0";
        if (!m_VM["common-envelope-mass-accretion-constant"].defaulted() && m_CommonEnvelopeMassAccretionConstant < 0.0) return "CE mass accretion constant (--common-envelope-mass-accretion-constant) < 0";
        if (!m_VM["common-envelope-mass-accretion-max"].defaulted() && m_CommonEnvelopeMassAccretionMax < 0.0) return "Maximum accreted mass (--common-envelope-mass-accretion-max) < 0";
        if (!m_VM["common-envelope-mass-accretion-min"].defaulted() && m_CommonEnvelopeMassAccretionMin < 0.0) return "Minimum accreted mass (--common-envelope-mass-accretion-min) < 0";

        if (m_DebugLevel < 0) return "Debug level (--debug-level) < 0";

        if (m_Eccentricity < 0.0 || m_Eccentricity > 1.0) return "Eccentricity (--eccentricity) must be between 0 and 1";
        if (m_EccentricityDistributionMin < 0.0 || m_EccentricityDistributionMin > 1.0) return "Minimum eccentricity (--eccentricity-min) must be between 0 and 1";
        if (m_EccentricityDistributionMax < 0.0 || m_EccentricityDistributionMax > 1.0) return "Maximum eccentricity (--eccentricity-max) must be between 0 and 1";
        if (m_EccentricityDistributionMax <= m_EccentricityDistributionMin) return "Maximum eccentricity (--eccentricity-max) must be > Minimum eccentricity (--eccentricity-min)";

        if (m_InitialMassFunctionMin < 0.0) return "Minimum initial mass (--initial-mass-min) < 0";
        if (m_InitialMassFunctionMax < 0.0) return "Maximum initial mass (--initial-mass-max) < 0";
        if (m_InitialMassFunctionMax <= m_InitialMassFunctionMin) return "Maximum initial mass (--initial-mass-max) must be > Minimum initial mass (--initial-mass-min)";

        if (m_KickMagnitudeDistribution == KICK_MAGNITUDE_DISTRIBUTION::FLAT) {
            if (m_KickMagnitudeDistributionMaximum <= 0.0) return "User specified --kick-magnitude-distribution = FLAT with Maximum kick magnitude (--kick-magnitude-max) <= 0.0";
        }

        if (m_LogLevel < 0) return "Logging level (--log-level) < 0";
 
        if (m_LuminousBlueVariableFactor < 0.0) return "LBV multiplier (--luminous-blue-variable-multiplier) < 0";

        if (m_MassRatioDistributionMin < 0.0 || m_MassRatioDistributionMin > 1.0) return "Minimum mass ratio (--mass-ratio-min) must be between 0 and 1";
        if (m_MassRatioDistributionMax < 0.0 || m_MassRatioDistributionMax > 1.0) return "Maximum mass ratio (--mass-ratio-max) must be between 0 and 1";
        if (m_MassRatioDistributionMax <= m_MassRatioDistributionMin) return "Maximum mass ratio (--mass-ratio-max) must be > Minimum mass ratio (--mass-ratio-min)";

        if (m_MaxEvolutionTime <= 0.0) return "Maximum evolution time in Myr (--maxEvolutionTime) must be > 0";

        if (m_Metallicity < 0.0 || m_Metallicity > 1.0) return "Metallicity (--metallicity) should be absolute metallicity and must be between 0 and 1";

        if (m_MinimumMassSecondary < 0.0) return "Seconday minimum mass (--minimum-secondary-mass) must be >= 0";
        if (m_MinimumMassSecondary > m_InitialMassFunctionMax) return "Seconday minimum mass (--minimum-secondary-mass) must be <= Maximum initial mass (--initial-mass-max)";

        if (m_NeutrinoMassLossAssumptionBH == NEUTRINO_MASS_LOSS_PRESCRIPTION::FIXED_MASS) {
            if (m_NeutrinoMassLossValueBH < 0.0) return "Neutrino mass loss value < 0";
        }

        if (m_ObjectsToEvolve <= 0) return m_EvolutionMode == EVOLUTION_MODE::SSE ? "Number of stars requested <= 0" : "Number of binaries requested <= 0";
    
        if (m_NeutrinoMassLossAssumptionBH == NEUTRINO_MASS_LOSS_PRESCRIPTION::FIXED_FRACTION) {
            if (m_NeutrinoMassLossValueBH < 0.0 || m_NeutrinoMassLossValueBH > 1.0) return "Neutrino mass loss must be between 0 and 1";
        }

        if (!m_VM["outputPath"].defaulted()) {                                                              // user specified output path?
                                                                                                            // yes
            fs::path userPath = m_OutputPathString;                                                         // user-specifed path
            if (fs::is_directory(userPath)) {                                                               // valid directory?
                m_OutputPath = userPath;                                                                    // yes - set outputPath to user-specified path
            }
            else {                                                                                          // not a valid directory
                m_OutputPath = m_DefaultOutputPath;                                                         // use default path = CWD
            }
        }

        if (m_PeriodDistributionMin < 0.0) return "Minimum orbital period (--orbital-period-min) < 0";
        if (m_PeriodDistributionMax < 0.0) return "Maximum orbital period (--orbital-period-max) < 0";

        if (!m_VM["pulsar-magnetic-field-decay-timescale"].defaulted() && m_PulsarMagneticFieldDecayTimescale <= 0.0) return "Pulsar magnetic field decay timescale (--pulsar-magnetic-field-decay-timescale) <= 0";
        if (!m_VM["pulsar-magnetic-field-decay-massscale"].defaulted() && m_PulsarMagneticFieldDecayMassscale <= 0.0) return "Pulsar Magnetic field decay massscale (--pulsar-magnetic-field-decay-massscale) <= 0";

        if (m_SemiMajorAxisDistributionMin < 0.0) return "Minimum semi-major Axis (--semi-major-axis-min) < 0";
        if (m_SemiMajorAxisDistributionMax < 0.0) return "Maximum semi-major Axis (--semi-major-axis-max) < 0";

        if (m_WolfRayetFactor < 0.0) return "WR multiplier (--wolf-rayet-multiplier) < 0";


    } catch (po::error& e) {                                                                                // program options exception
        errStr = e.what();
    } catch (...) {                                                                                         // unhandled exception
        errStr = "unhandled exception";
    }

    return errStr;
}


/*
 * determine if the user specified a value for the option
 * 
 * 
 * int OptionSpecified(std::string p_OptionString) 
 * 
 * @param
 * @param   [IN]    p_OptionString              String containing option name
 * @return                                      Int result:
 *                                                  -1: invalid/unknown option name
 *                                                   0: option was not specified by user - default value used
 *                                                   1: option specified by user - user specified value used
 */
int Options::OptionSpecified(std::string p_OptionString) {

    int  result = -1;                                                                           // default = invalid/unknown option
    
    try {

        if (m_EvolvingObject.optionValues.m_VM.count(p_OptionString) > 0) {                     // option exists at object level?
            result = m_EvolvingObject.optionValues.m_VM[p_OptionString].defaulted() ? 0 : 1;    // yes - set result
        }
        else {                                                                                  // option does not exist at object level
            if (m_Program.optionValues.m_VM.count(p_OptionString) > 0) {                        // option exists at program level
                result = m_Program.optionValues.m_VM[p_OptionString].defaulted() ? 0 : 1;       // yes - set result
            }
            else {                                                                              // option does not exist at program level
                result = -1;                                                                    // set result
            }    
        }
    } catch (po::error& e) {                                                                    // program options exception
        result = -1;
    } catch (...) {                                                                             // unhandled exception
        result = -1;                                                                            // set return value - invalid/unknown option
    }
    
    return result;
}


// JR: todo: For SetProgramOptions() and SetObjectOptions(), one day we should construct the list of
//           options shown in the help string from the maps in constants.h rather than have them
//           here as literal strings - too much opportunity for then to get out of sync doing it
//           this way (and more work to update the strings here every time something changes)


bool Options::SetProgramOptions(OptionValues *p_Options, po::options_description *p_OptionsDescription) {

    bool ok = true;                             // status - unless a problem occurs

    // create default strings for vector<string> types (too hard to do inline)

    std::ostringstream ss;

    // debug classes
    string defaultDebugClasses;
    ss << "";
    for (auto debugClass = p_Options->m_DebugClasses.begin(); debugClass != p_Options->m_DebugClasses.end(); ++debugClass) ss << *debugClass << ",";
    defaultDebugClasses = ss.str();
    if (defaultDebugClasses.size() > 0) defaultDebugClasses.erase(defaultDebugClasses.size() - 1);

    // log classes
    string defaultLogClasses;
    ss << "";
    for (auto logClass = p_Options->m_LogClasses.begin(); logClass != p_Options->m_LogClasses.end(); ++logClass) ss << *logClass << ",";
    defaultLogClasses = ss.str();
    if (defaultLogClasses.size() > 0) defaultLogClasses.erase(defaultLogClasses.size() - 1);


    // add options

    try {

        p_OptionsDescription->add_options()     // begin the list of options to be added - boost syntactic sugar

        // there is no good way of formatting these - the boost syntax doesn't help that much
        // there is just a boatload of options, so this function (and similar functions) are just going to be long...


        // switches

        (
            "help,h",                                                      
            po::bool_switch(), "Print this help message"
        )
        (
            "version,v",                                                   
            po::bool_switch(), "Print COMPAS version string"
        )


        // boolean options - alphabetically

        // Floor
        /*
        (
            "AIS-exploratory-phase",                                       
            po::value<bool>(&p_Options->m_AISexploratoryPhase)->default_value(p_Options->m_AISexploratoryPhase)->implicit_value(true),                                                            
            ("Run exploratory phase of STROOPWAFEL (default = " + std::string(p_Options->m_AISexploratoryPhase ? "TRUE" : "FALSE") + ")").c_str()
        )
        (
            "AIS-Hubble",                                                  
            po::value<bool>(&p_Options->m_AIShubble)->default_value(p_Options->m_AIShubble)->implicit_value(true),                                                                                
            ("Excluding not in Hubble time mergers selection in exploratory phase of STROOPWAFEL (default = " + std::string(p_Options->m_AIShubble ? "TRUE" : "FALSE") + ")").c_str()
        )
        (
            "AIS-Pessimistic",                                             
            po::value<bool>(&p_Options->m_AISpessimistic)->default_value(p_Options->m_AISpessimistic)->implicit_value(true),                                                                      
            ("Optimistic or Pessimistic selection in exploratory phase of STROOPWAFEL (default = " + std::string(p_Options->m_AISpessimistic ? "TRUE" : "FALSE") + ")").c_str()
        )
        (
            "AIS-refinement-phase",                                        
            po::value<bool>(&p_Options->m_AISrefinementPhase)->default_value(p_Options->m_AISrefinementPhase)->implicit_value(true),                                                              
            ("Run main sampling phase (step2) of STROOPWAFEL (default = " + std::string(p_Options->m_AISrefinementPhase ? "TRUE" : "FALSE") + ")").c_str()
        )
        (
            "AIS-RLOF",                                                    
            po::value<bool>(&p_Options->m_AISrlof)->default_value(p_Options->m_AISrlof)->implicit_value(true),                                                                                    
            ("RLOFSecondaryZAMS selection in exploratory phase of STROOPWAFEL (default = " + std::string(p_Options->m_AISrlof ? "TRUE" : "FALSE") + ")").c_str()
       )
        */

        (
            "debug-to-file",                                               
            po::value<bool>(&p_Options->m_DebugToFile)->default_value(p_Options->m_DebugToFile)->implicit_value(true),                                                                            
            ("Write debug statements to file (default = " + std::string(p_Options->m_DebugToFile ? "TRUE" : "FALSE") + ")").c_str()
        )
        (
            "detailedOutput",                                              
            po::value<bool>(&p_Options->m_DetailedOutput)->default_value(p_Options->m_DetailedOutput)->implicit_value(true),                                                                      
            ("Print detailed output to file (default = " + std::string(p_Options->m_DetailedOutput ? "TRUE" : "FALSE") + ")").c_str()
        )
        (
            "enable-warnings",                                             
            po::value<bool>(&p_Options->m_EnableWarnings)->default_value(p_Options->m_EnableWarnings)->implicit_value(true),                                                                      
            ("Display warning messages to stdout (default = " + std::string(p_Options->m_EnableWarnings ? "TRUE" : "FALSE") + ")").c_str()
        )
        (
            "errors-to-file",                                              
            po::value<bool>(&p_Options->m_ErrorsToFile)->default_value(p_Options->m_ErrorsToFile)->implicit_value(true),                                                                          
            ("Write error messages to file (default = " + std::string(p_Options->m_ErrorsToFile ? "TRUE" : "FALSE") + ")").c_str()
        )
        (
            "populationDataPrinting",                                      
            po::value<bool>(&p_Options->m_PopulationDataPrinting)->default_value(p_Options->m_PopulationDataPrinting)->implicit_value(true),                                                      
            ("Print details of population (default = " + std::string(p_Options->m_PopulationDataPrinting ? "TRUE" : "FALSE") + ")").c_str()
        )
        (
            "quiet",                                                       
            po::value<bool>(&p_Options->m_Quiet)->default_value(p_Options->m_Quiet)->implicit_value(true),                                                                                        
            ("Suppress printing (default = " + std::string(p_Options->m_Quiet ? "TRUE" : "FALSE") + ")").c_str()
        )
        (
            "RLOFPrinting",                                                
            po::value<bool>(&p_Options->m_RlofPrinting)->default_value(p_Options->m_RlofPrinting)->implicit_value(true),                                                                          
            ("Enable output parameters before/after RLOF (default = " + std::string(p_Options->m_RlofPrinting ? "TRUE" : "FALSE") + ")").c_str()
        )
        (
            "switchLog",                                                
            po::value<bool>(&p_Options->m_SwitchLog)->default_value(p_Options->m_SwitchLog)->implicit_value(true),                                                                          
            ("Print switch log to file (default = " + std::string(p_Options->m_SwitchLog ? "TRUE" : "FALSE") + ")").c_str()
        )


        // numerical options - alphabetically grouped by type 

        // int

        (
            "debug-level",                                                 
            po::value<int>(&p_Options->m_DebugLevel)->default_value(p_Options->m_DebugLevel),                                                                                                     
            ("Determines which print statements are displayed for debugging (default = " + std::to_string(p_Options->m_DebugLevel) + ")").c_str()
        )
        (
            "log-level",                                                   
            po::value<int>(&p_Options->m_LogLevel)->default_value(p_Options->m_LogLevel),                                                                                                         
            ("Determines which print statements are included in the logfile (default = " + std::to_string(p_Options->m_LogLevel) + ")").c_str()
        )

        // Floor
        /*
        (
            "nbatches-used",                                               
            po::value<int>(&p_Options->m_nBatchesUsed)->default_value(p_Options->m_nBatchesUsed),                                                                                                 
            ("Number of batches used, for STROOPWAFEL (AIS), -1 = not required (default = " + std::to_string(p_Options->m_nBatchesUsed) + ")").c_str()
        )
        */

        (
            "number-of-stars",                                        
            po::value<int>(&p_Options->m_ObjectsToEvolve)->default_value(p_Options->m_ObjectsToEvolve),                                                                                                       
            ("Specify the number of stars to simulate (SSE) (default = " + std::to_string(p_Options->m_ObjectsToEvolve) + ")").c_str()
        )


        // double

        // Floor
        /*
        (
            "kappa-gaussians",                                             
            po::value<double>(&p_Options->m_KappaGaussians)->default_value(p_Options->m_KappaGaussians),                                                                                          
            ("Scaling factor for the width of the Gaussian distributions in STROOPWAFEL main sampling phase (default = " + std::to_string(p_Options->m_KappaGaussians) + ")").c_str()
        )
        */


        // string options - alphabetically

        // Floor
        /*
        (
            "AIS-DCOtype",                                                 
            po::value<string>(&p_Options->m_AISDCOtypeString)->default_value(p_Options->m_AISDCOtypeString),                                                                                      
            ("DCO type selection in exploratory phase of STROOPWAFEL, (options: [ALL, BBH, BNS, BHNS], default = " + p_Options->m_AISDCOtypeString + ")").c_str()
        )
        */

        (
            "grid",                                                        
            po::value<string>(&p_Options->m_GridFilename)->default_value(p_Options->m_GridFilename)->implicit_value(""),                                                                      
            ("Grid filename (default = " + p_Options->m_GridFilename + ")").c_str()
        )

        // Serena
        /*
        (
            "logfile-be-binaries",                                     
            po::value<string>(&p_Options->m_LogfileBeBinaries)->default_value(p_Options->m_LogfileBeBinaries),                                                                              
            ("Filename for Be Binaries logfile (default = " + p_Options->m_LogfileBeBinaries + ")").c_str()
        )
        */

        (
            "logfile-rlof-parameters",                                 
            po::value<string>(&p_Options->m_LogfileRLOFParameters)->default_value(p_Options->m_LogfileRLOFParameters),                                                                      
            ("Filename for RLOF Parameters logfile ( default = " + p_Options->m_LogfileRLOFParameters + ")").c_str()
        )
        (
            "logfile-common-envelopes",                                
            po::value<string>(&p_Options->m_LogfileCommonEnvelopes)->default_value(p_Options->m_LogfileCommonEnvelopes),                                                                    
            ("Filename for Common Envelopes logfile (default = " + p_Options->m_LogfileCommonEnvelopes + ")").c_str()
        )
        (
            "logfile-detailed-output",                                 
            po::value<string>(&p_Options->m_LogfileDetailedOutput)->default_value(p_Options->m_LogfileDetailedOutput),                                                                      
            ("Filename for Detailed Output logfile (default = " + p_Options->m_LogfileDetailedOutput + ")").c_str()
        )
        (
            "logfile-double-compact-objects",                          
            po::value<string>(&p_Options->m_LogfileDoubleCompactObjects)->default_value(p_Options->m_LogfileDoubleCompactObjects),                                                          
            ("Filename for Double Compact Objects logfile (default = " + p_Options->m_LogfileDoubleCompactObjects + ")").c_str()
        )
        (
            "logfile-pulsar-evolution",                                
            po::value<string>(&p_Options->m_LogfilePulsarEvolution)->default_value(p_Options->m_LogfilePulsarEvolution),                                                                    
            ("Filename for Pulsar Evolution logfile (default = " + p_Options->m_LogfilePulsarEvolution + ")").c_str()
        )
        (
            "logfile-supernovae",                                      
            po::value<string>(&p_Options->m_LogfileSupernovae)->default_value(p_Options->m_LogfileSupernovae),                                                                              
            ("Filename for Supernovae logfile (default = " + p_Options->m_LogfileSupernovae + ")").c_str()
        )
        (
            "logfile-system-parameters",                               
            po::value<string>(&p_Options->m_LogfileSystemParameters)->default_value(p_Options->m_LogfileSystemParameters),                                                                  
            ("Filename for System Parameters logfile (default = " + p_Options->m_LogfileSystemParameters + ")").c_str()
        )
        (
            "logfile-definitions",                                         
            po::value<string>(&p_Options->m_LogfileDefinitionsFilename)->default_value(p_Options->m_LogfileDefinitionsFilename)->implicit_value(""),                                              
            ("Filename for logfile record definitions (default = " + p_Options->m_LogfileDefinitionsFilename + ")").c_str()
        )
        (
            "logfile-delimiter",                                           
            po::value<string>(&p_Options->m_LogfileDelimiterString)->default_value(p_Options->m_LogfileDelimiterString),                                                                          
            ("Field delimiter for logfile records (default = " + p_Options->m_LogfileDelimiterString + ")").c_str()
        )
        (
            "logfile-name-prefix",                                         
            po::value<string>(&p_Options->m_LogfileNamePrefix)->default_value(p_Options->m_LogfileNamePrefix)->implicit_value(""),                                                                
            ("Prefix for logfile names (default = " + p_Options->m_LogfileNamePrefix + ")").c_str()
        )
        (
            "logfile-switch-log",                                      
            po::value<string>(&p_Options->m_LogfileSwitchLog)->default_value(p_Options->m_LogfileSwitchLog),                                                                                
            ("Filename for Switch Log logfile (default = " + p_Options->m_LogfileSwitchLog + ")").c_str()
        )

        (
            "mode",                                                 
            po::value<string>(&p_Options->m_EvolutionModeString)->default_value(p_Options->m_EvolutionModeString),                                                                              
            ("Evolution mode (options: [SSE, BSE], default = " + p_Options->m_EvolutionModeString + ")").c_str()
        )

        (
            "output-container,c",                                          
            po::value<string>(&p_Options->m_OutputContainerName)->default_value(p_Options->m_OutputContainerName)->implicit_value(""),                                                            
            ("Container (directory) name for output files (default = " + p_Options->m_OutputContainerName + ")").c_str()
        )
        (
            "outputPath,o",                                                
            po::value<string>(&p_Options->m_OutputPathString)->default_value(p_Options->m_OutputPathString)->implicit_value(""),                                                                  
            ("Directory for output (default = " + p_Options->m_OutputPathString + ")").c_str()
        )


        // vector (list) options - alphabetically

        (
            "debug-classes",                                               
            po::value<vector<string>>(&p_Options->m_DebugClasses)->multitoken()->default_value(p_Options->m_DebugClasses),                                                                        
            ("Debug classes enabled (default = " + defaultDebugClasses + ")").c_str()
        )
        (
            "log-classes",                                                 
            po::value<vector<string>>(&p_Options->m_LogClasses)->multitoken()->default_value(p_Options->m_LogClasses),                                                                            
            ("Logging classes enabled (default = " + defaultLogClasses + ")").c_str()
        )
    
        ;   // end the list of options to be added

    } catch (po::error& e) {    // program options exception
        ok = false;             // set status
    } catch (...) {             // unhandled exception
        ok = false;             // set status
    }

    return ok;
}


bool Options::SetObjectOptions(OptionValues *p_Options, po::options_description *p_OptionsDescription) {

    bool ok = true;                             // status - unless a problem occurs

    // add options

    try {

        p_OptionsDescription->add_options()     // begin the list of options to be added - boost syntactic sugar

        // there is no good way of formatting these - the boost syntz doesn't help that much
        // there is just a boatload of options, so this function (and similar functions) are just going to be long...
    

        // boolean options - alphabetically 

        (
            "allow-rlof-at-birth",                                         
            po::value<bool>(&p_Options->m_AllowRLOFAtBirth)->default_value(p_Options->m_AllowRLOFAtBirth)->implicit_value(true),                                                                  
            ("Allow binaries that have one or both stars in RLOF at birth to evolve (default = " + std::string(p_Options->m_AllowRLOFAtBirth ? "TRUE" : "FALSE") + ")").c_str()
        )
        (
            "allow-touching-at-birth",                                     
            po::value<bool>(&p_Options->m_AllowTouchingAtBirth)->default_value(p_Options->m_AllowTouchingAtBirth)->implicit_value(true),                                                          
            ("Allow binaries that are touching at birth to evolve (default = " + std::string(p_Options->m_AllowTouchingAtBirth ? "TRUE" : "FALSE") + ")").c_str()
        )
        (
            "angularMomentumConservationDuringCircularisation",            
            po::value<bool>(&p_Options->m_AngularMomentumConservationDuringCircularisation)->default_value(p_Options->m_AngularMomentumConservationDuringCircularisation)->implicit_value(true),  
            ("Conserve angular momentum when binary is circularised when entering a Mass Transfer episode (default = " + std::string(p_Options->m_AngularMomentumConservationDuringCircularisation ? "TRUE" : "FALSE") + ")").c_str()
        )

        // Serena
        /* 
        (
            "BeBinaries",                                                  
            po::value<bool>(&p_Options->m_BeBinaries)->default_value(p_Options->m_BeBinaries)->implicit_value(true),                                                                              
            ("Enable Be Binaries study (default = " + std::string(p_Options->m_BeBinaries ? "TRUE" : "FALSE") + ")").c_str()
        )
        */

        (
            "circulariseBinaryDuringMassTransfer",                         
            po::value<bool>(&p_Options->m_CirculariseBinaryDuringMassTransfer)->default_value(p_Options->m_CirculariseBinaryDuringMassTransfer)->implicit_value(true),                            
            ("Circularise binary when it enters a Mass Transfer episode (default = " + std::string(p_Options->m_CirculariseBinaryDuringMassTransfer ? "TRUE" : "FALSE") + ")").c_str()
        )
        (
            "common-envelope-allow-main-sequence-survive",                 
            po::value<bool>(&p_Options->m_AllowMainSequenceStarToSurviveCommonEnvelope)->default_value(p_Options->m_AllowMainSequenceStarToSurviveCommonEnvelope)->implicit_value(true),          
            ("Allow main sequence stars to survive common envelope evolution (default = " + std::string(p_Options->m_AllowMainSequenceStarToSurviveCommonEnvelope ? "TRUE" : "FALSE") + ")").c_str()
        )
        (
            "evolve-pulsars",                                              
            po::value<bool>(&p_Options->m_EvolvePulsars)->default_value(p_Options->m_EvolvePulsars)->implicit_value(true),                                                                        
            ("Evolve pulsars (default = " + std::string(p_Options->m_EvolvePulsars ? "TRUE" : "FALSE") + ")").c_str()
        )
        (
            "evolve-unbound-systems",                                      
            po::value<bool>(&p_Options->m_EvolveUnboundSystems)->default_value(p_Options->m_EvolveUnboundSystems)->implicit_value(true),                                                          
            ("Continue evolving stars even if the binary is disrupted (default = " + std::string(p_Options->m_EvolveUnboundSystems ? "TRUE" : "FALSE") + ")").c_str()
        )
        (
            "massTransfer",                                                
            po::value<bool>(&p_Options->m_UseMassTransfer)->default_value(p_Options->m_UseMassTransfer)->implicit_value(true),                                                                    
            ("Enable mass transfer (default = " + std::string(p_Options->m_UseMassTransfer ? "TRUE" : "FALSE") + ")").c_str()
        )
        (
            "pair-instability-supernovae",                                 
            po::value<bool>(&p_Options->m_UsePairInstabilitySupernovae)->default_value(p_Options->m_UsePairInstabilitySupernovae)->implicit_value(true),                                          
            ("Enable pair instability supernovae (PISN) (default = " + std::string(p_Options->m_UsePairInstabilitySupernovae ? "TRUE" : "FALSE") + ")").c_str()
        )
        (
            "print-bool-as-string",                                        
            po::value<bool>(&p_Options->m_PrintBoolAsString)->default_value(p_Options->m_PrintBoolAsString)->implicit_value(true),                                                                
            ("Print boolean properties as 'TRUE' or 'FALSE' (default = " + std::string(p_Options->m_PrintBoolAsString ? "TRUE" : "FALSE") + ")").c_str()
        )
        (
            "pulsational-pair-instability",                                
            po::value<bool>(&p_Options->m_UsePulsationalPairInstability)->default_value(p_Options->m_UsePulsationalPairInstability)->implicit_value(true),                                        
            ("Enable mass loss due to pulsational-pair-instability (PPI) (default = " + std::string(p_Options->m_UsePulsationalPairInstability ? "TRUE" : "FALSE") + ")").c_str()
        )
        (
            "revised-energy-formalism-Nandez-Ivanova",                     
            po::value<bool>(&p_Options->m_RevisedEnergyFormalismNandezIvanova)->default_value(p_Options->m_RevisedEnergyFormalismNandezIvanova)->implicit_value(true),                            
            ("Enable revised energy formalism (default = " + std::string(p_Options->m_RevisedEnergyFormalismNandezIvanova ? "TRUE" : "FALSE") + ")").c_str()
        )

        (
            "use-mass-loss",                                               
            po::value<bool>(&p_Options->m_UseMassLoss)->default_value(p_Options->m_UseMassLoss)->implicit_value(true),                                                                            
            ("Enable mass loss (default = " + std::string(p_Options->m_UseMassLoss ? "TRUE" : "FALSE") + ")").c_str()
        )


        // numerical options - alphabetically grouped by type

        // unsigned long

        (
            "random-seed",                                                 
            po::value<unsigned long>(&p_Options->m_RandomSeed)->default_value(p_Options->m_RandomSeed),                                                                                           
            ("Random seed to use (default = " + std::to_string(p_Options->m_RandomSeed) + ")").c_str()
        )


        // int

        (
            "maximum-number-timestep-iterations",                          
            po::value<int>(&p_Options->m_MaxNumberOfTimestepIterations)->default_value(p_Options->m_MaxNumberOfTimestepIterations),                                                               
            ("Maximum number of timesteps to evolve binary before giving up (default = " + std::to_string(p_Options->m_MaxNumberOfTimestepIterations) + ")").c_str()
        )


        // double

        (
            "common-envelope-alpha",                                       
            po::value<double>(&p_Options->m_CommonEnvelopeAlpha)->default_value(p_Options->m_CommonEnvelopeAlpha),                                                                                
            ("Common Envelope efficiency alpha (default = " + std::to_string(p_Options->m_CommonEnvelopeAlpha) + ")").c_str()
        )
        (
            "common-envelope-alpha-thermal",                               
            po::value<double>(&p_Options->m_CommonEnvelopeAlphaThermal)->default_value(p_Options->m_CommonEnvelopeAlphaThermal),                                                                  
            ("Defined such that lambda = alpha_th * lambda_b + (1.0 - alpha_th) * lambda_g (default = " + std::to_string(p_Options->m_CommonEnvelopeAlphaThermal) + ")").c_str()
        )
        (
            "common-envelope-lambda",                                      
            po::value<double>(&p_Options->m_CommonEnvelopeLambda)->default_value(p_Options->m_CommonEnvelopeLambda),                                                                              
            ("Common Envelope lambda (default = " + std::to_string(p_Options->m_CommonEnvelopeLambda) + ")").c_str()
        )
        (
            "common-envelope-lambda-multiplier",                           
            po::value<double>(&p_Options->m_CommonEnvelopeLambdaMultiplier)->default_value(p_Options->m_CommonEnvelopeLambdaMultiplier),                                                          
            ("Multiply lambda by some constant (default = " + std::to_string(p_Options->m_CommonEnvelopeLambdaMultiplier) + ")").c_str()
        )
        (
            "common-envelope-mass-accretion-constant",                     
            po::value<double>(&p_Options->m_CommonEnvelopeMassAccretionConstant)->default_value(p_Options->m_CommonEnvelopeMassAccretionConstant),                                                
            ("Value of mass accreted by NS/BH during common envelope evolution if assuming all NS/BH accrete same amount of mass (common-envelope-mass-accretion-prescription CONSTANT). Ignored otherwise (default = " + std::to_string(p_Options->m_CommonEnvelopeMassAccretionConstant) + ")").c_str()
        )
        (
            "common-envelope-mass-accretion-max",                          
            po::value<double>(&p_Options->m_CommonEnvelopeMassAccretionMax)->default_value(p_Options->m_CommonEnvelopeMassAccretionMax),                                                          
            ("Maximum amount of mass accreted by NS/BHs during common envelope evolution in solar masses (default = " + std::to_string(p_Options->m_CommonEnvelopeMassAccretionMax) + ")").c_str()
        )
        (
            "common-envelope-mass-accretion-min",                          
            po::value<double>(&p_Options->m_CommonEnvelopeMassAccretionMin)->default_value(p_Options->m_CommonEnvelopeMassAccretionMin),                                                          
            ("Minimum amount of mass accreted by NS/BHs during common envelope evolution in solar masses (default = " + std::to_string(p_Options->m_CommonEnvelopeMassAccretionMin) + ")").c_str()
        )
        (
            "common-envelope-recombination-energy-density",                
            po::value<double>(&p_Options->m_CommonEnvelopeRecombinationEnergyDensity)->default_value(p_Options->m_CommonEnvelopeRecombinationEnergyDensity),                                      
            ("Recombination energy density in erg/g (default = " + std::to_string(p_Options->m_CommonEnvelopeRecombinationEnergyDensity) + ")").c_str()
        )
        (
            "common-envelope-slope-Kruckow",                               
            po::value<double>(&p_Options->m_CommonEnvelopeSlopeKruckow)->default_value(p_Options->m_CommonEnvelopeSlopeKruckow),                                                                  
            ("Common Envelope slope for Kruckow lambda (default = " + std::to_string(p_Options->m_CommonEnvelopeSlopeKruckow) + ")").c_str()
        )

        // AVG
        /*
        (
            "critical-mass-ratio-giant-degenerate-accretor",               
            po::value<double>(&p_Options->m_MassTransferCriticalMassRatioGiantDegenerateAccretor)->default_value(p_Options->m_MassTransferCriticalMassRatioGiantDegenerateAccretor),              
            ("Critical mass ratio for MT from a giant star (default = " + std::to_string(p_Options->m_MassTransferCriticalMassRatioGiantDegenerateAccretor) + ") Specify both giant flags to use. 0 is always stable, <0 is disabled").c_str()
        )
        (
            "critical-mass-ratio-giant-non-degenerate-accretor",           
            po::value<double>(&p_Options->m_MassTransferCriticalMassRatioGiantNonDegenerateAccretor)->default_value(p_Options->m_MassTransferCriticalMassRatioGiantNonDegenerateAccretor),        
            ("Critical mass ratio for MT from a giant star (default = " + std::to_string(p_Options->m_MassTransferCriticalMassRatioGiantNonDegenerateAccretor) + ") Specify both giant flags to use. 0 is always stable, <0 is disabled").c_str()
        )
        (
            "critical-mass-ratio-helium-giant-degenerate-accretor",        
            po::value<double>(&p_Options->m_MassTransferCriticalMassRatioHeliumGiantDegenerateAccretor)->default_value(p_Options->m_MassTransferCriticalMassRatioHeliumGiantDegenerateAccretor),  
            ("Critical mass ratio for MT from a helium giant star (default = " + std::to_string(p_Options->m_MassTransferCriticalMassRatioHeliumGiantDegenerateAccretor) + ") Specify both helium giant flags to use. 0 is always stable, <0 is disabled").c_str()
        )
        (
            "critical-mass-ratio-helium-giant-non-degenerate-accretor",    
            po::value<double>(&p_Options->m_MassTransferCriticalMassRatioHeliumGiantNonDegenerateAccretor)->default_value(p_Options->m_MassTransferCriticalMassRatioHeliumGiantNonDegenerateAccretor), 
            ("Critical mass ratio for MT from a helium giant star (default = " + std::to_string(p_Options->m_MassTransferCriticalMassRatioHeliumGiantNonDegenerateAccretor) + ") Specify both helium giant flags to use. 0 is always stable, <0 is disabled").c_str()
        )
        (
            "critical-mass-ratio-helium-HG-degenerate-accretor",           
            po::value<double>(&p_Options->m_MassTransferCriticalMassRatioHeliumHGDegenerateAccretor)->default_value(p_Options->m_MassTransferCriticalMassRatioHeliumHGDegenerateAccretor),        
            ("Critical mass ratio for MT from a helium HG star (default = " + std::to_string(p_Options->m_MassTransferCriticalMassRatioHeliumHGDegenerateAccretor) + ") Specify both helium HG flags to use. 0 is always stable, <0 is disabled").c_str()
        )
        (
            "critical-mass-ratio-helium-HG-non-degenerate-accretor",       
            po::value<double>(&p_Options->m_MassTransferCriticalMassRatioHeliumHGNonDegenerateAccretor)->default_value(p_Options->m_MassTransferCriticalMassRatioHeliumHGNonDegenerateAccretor),  
            ("Critical mass ratio for MT from a helium HG star (default = " + std::to_string(p_Options->m_MassTransferCriticalMassRatioHeliumHGNonDegenerateAccretor) + ") Specify both helium HG flags to use. 0 is always stable, <0 is disabled").c_str()
        )
        (
            "critical-mass-ratio-helium-MS-degenerate-accretor",           
            po::value<double>(&p_Options->m_MassTransferCriticalMassRatioHeliumMSDegenerateAccretor)->default_value(p_Options->m_MassTransferCriticalMassRatioHeliumMSDegenerateAccretor),        
            ("Critical mass ratio for MT from a helium MS star (default = " + std::to_string(p_Options->m_MassTransferCriticalMassRatioHeliumMSDegenerateAccretor) + ") Specify both helium MS flags to use. 0 is always stable, <0 is disabled").c_str()
        )
        (
            "critical-mass-ratio-helium-MS-non-degenerate-accretor",       
            po::value<double>(&p_Options->m_MassTransferCriticalMassRatioHeliumMSNonDegenerateAccretor)->default_value(p_Options->m_MassTransferCriticalMassRatioHeliumMSNonDegenerateAccretor),  
            ("Critical mass ratio for MT from a helium MS star (default = " + std::to_string(p_Options->m_MassTransferCriticalMassRatioHeliumMSNonDegenerateAccretor) + ") Specify both helium MS flags to use. 0 is always stable, <0 is disabled").c_str()
        )
        (
            "critical-mass-ratio-HG-degenerate-accretor",                  
            po::value<double>(&p_Options->m_MassTransferCriticalMassRatioHGDegenerateAccretor)->default_value(p_Options->m_MassTransferCriticalMassRatioHGDegenerateAccretor),                    
            ("Critical mass ratio for MT from a HG star (default = " + std::to_string(p_Options->m_MassTransferCriticalMassRatioHGDegenerateAccretor) + ") Specify both HG flags to use. 0 is always stable, <0 is disabled").c_str()
        )
        (
            "critical-mass-ratio-HG-non-degenerate-accretor",              
            po::value<double>(&p_Options->m_MassTransferCriticalMassRatioHGNonDegenerateAccretor)->default_value(p_Options->m_MassTransferCriticalMassRatioHGNonDegenerateAccretor),              
            ("Critical mass ratio for MT from a HG star (default = " + std::to_string(p_Options->m_MassTransferCriticalMassRatioHGNonDegenerateAccretor) + ") Specify both HG flags to use. 0 is always stable, <0 is disabled").c_str()
        )
        (
            "critical-mass-ratio-MS-high-mass-degenerate-accretor",        
            po::value<double>(&p_Options->m_MassTransferCriticalMassRatioMSHighMassDegenerateAccretor)->default_value(p_Options->m_MassTransferCriticalMassRatioMSHighMassDegenerateAccretor),    
            ("Critical mass ratio for MT from a MS star to a degenerate accretor (default = " + std::to_string(p_Options->m_MassTransferCriticalMassRatioMSHighMassDegenerateAccretor) + " Specify both MS high mass flags to use. 0 is always stable, <0 is disabled").c_str()
        )
        (
            "critical-mass-ratio-MS-high-mass-non-degenerate-accretor",    
            po::value<double>(&p_Options->m_MassTransferCriticalMassRatioMSHighMassNonDegenerateAccretor)->default_value(p_Options->m_MassTransferCriticalMassRatioMSHighMassNonDegenerateAccretor), 
            ("Critical mass ratio for MT from a MS star (default = " + std::to_string(p_Options->m_MassTransferCriticalMassRatioMSHighMassNonDegenerateAccretor) + ") Specify both MS high mass flags to use. 0 is always stable, <0 is disabled").c_str()
        )
        (
            "critical-mass-ratio-MS-low-mass-degenerate-accretor",         
            po::value<double>(&p_Options->m_MassTransferCriticalMassRatioMSLowMassDegenerateAccretor)->default_value(p_Options->m_MassTransferCriticalMassRatioMSLowMassDegenerateAccretor),      
            ("Critical mass ratio for MT from a MS star to a degenerate accretor (default = " + std::to_string(p_Options->m_MassTransferCriticalMassRatioMSLowMassDegenerateAccretor) + " Specify both MS low mass flags to use. 0 is always stable, <0 is disabled").c_str()
        )
        (
            "critical-mass-ratio-MS-low-mass-non-degenerate-accretor",     
            po::value<double>(&p_Options->m_MassTransferCriticalMassRatioMSLowMassNonDegenerateAccretor)->default_value(p_Options->m_MassTransferCriticalMassRatioMSLowMassNonDegenerateAccretor), 
            ("Critical mass ratio for MT from a MS star (default = " + std::to_string(p_Options->m_MassTransferCriticalMassRatioMSLowMassNonDegenerateAccretor) + ") Specify both MS low mass flags to use. 0 is always stable, <0 is disabled").c_str()
        )
        (
            "critical-mass-ratio-white-dwarf-degenerate-accretor",         
            po::value<double>(&p_Options->m_MassTransferCriticalMassRatioWhiteDwarfDegenerateAccretor)->default_value(p_Options->m_MassTransferCriticalMassRatioWhiteDwarfDegenerateAccretor),    
            ("Critical mass ratio for MT from a white dwarf (default = " + std::to_string(p_Options->m_MassTransferCriticalMassRatioWhiteDwarfDegenerateAccretor) + ") Specify both white dwarf flags to use. 0 is always stable, <0 is disabled").c_str()
        )
        (
            "critical-mass-ratio-white-dwarf-non-degenerate-accretor",     
            po::value<double>(&p_Options->m_MassTransferCriticalMassRatioWhiteDwarfNonDegenerateAccretor)->default_value(p_Options->m_MassTransferCriticalMassRatioWhiteDwarfNonDegenerateAccretor), 
            ("Critical mass ratio for MT from a white dwarf (default = " + std::to_string(p_Options->m_MassTransferCriticalMassRatioWhiteDwarfNonDegenerateAccretor) + ") Specify both white dwarf flags to use. 0 is always stable, <0 is disabled").c_str()
        )
        */

        (
            "eccentricity,e",                                            
            po::value<double>(&p_Options->m_Eccentricity)->default_value(p_Options->m_Eccentricity),                                                                
            ("Eccentricity, e (default = " + std::to_string(p_Options->m_Eccentricity) + ")").c_str()
        )

        (
            "eccentricity-max",                                            
            po::value<double>(&p_Options->m_EccentricityDistributionMax)->default_value(p_Options->m_EccentricityDistributionMax),                                                                
            ("Maximum eccentricity to generate (default = " + std::to_string(p_Options->m_EccentricityDistributionMax) + ")").c_str()
        )
        (
            "eccentricity-min",                                            
            po::value<double>(&p_Options->m_EccentricityDistributionMin)->default_value(p_Options->m_EccentricityDistributionMin),                                                                
            ("Minimum eccentricity to generate (default = " + std::to_string(p_Options->m_EccentricityDistributionMin) + ")").c_str()
        )
        (
            "eddington-accretion-factor",                                  
            po::value<double>(&p_Options->m_EddingtonAccretionFactor)->default_value(p_Options->m_EddingtonAccretionFactor),                                                                      
            ("Multiplication factor for eddington accretion for NS & BH, i.e. >1 is super-eddington and 0. is no accretion (default = " + std::to_string(p_Options->m_EddingtonAccretionFactor) + ")").c_str()
        )

        (
            "fix-dimensionless-kick-magnitude",                            
            po::value<double>(&p_Options->m_FixedUK)->default_value(p_Options->m_FixedUK),                                                                                                        
            ("Fix dimensionless kick magnitude uk to this value (default = " + std::to_string(p_Options->m_FixedUK) + ", -ve values false, +ve values true)").c_str()
        )

        (
            "initial-mass",                                            
            po::value<double>(&p_Options->m_InitialMass)->default_value(p_Options->m_InitialMass),                                                                          
            ("Initial mass (in Msol) for the star (SSE) (default = " + std::to_string(p_Options->m_InitialMass) + ")").c_str()
        )
        (
            "initial-mass-1",                                            
            po::value<double>(&p_Options->m_InitialMass1)->default_value(p_Options->m_InitialMass1),                                                                          
            ("Initial mass (in Msol) for the primary star (BSE) (default = " + std::to_string(p_Options->m_InitialMass1) + ")").c_str()
        )
        (
            "initial-mass-2",                                            
            po::value<double>(&p_Options->m_InitialMass2)->default_value(p_Options->m_InitialMass2),
            ("Initial mass (in Msol) for the secondary star (BSE) (default = " + std::to_string(p_Options->m_InitialMass2) + ")").c_str()
        )
        (
            "initial-mass-max",                                            
            po::value<double>(&p_Options->m_InitialMassFunctionMax)->default_value(p_Options->m_InitialMassFunctionMax),                                                                          
            ("Maximum mass (in Msol) to generate using given IMF (default = " + std::to_string(p_Options->m_InitialMassFunctionMax) + ")").c_str()
        )
        (
            "initial-mass-min",                                            
            po::value<double>(&p_Options->m_InitialMassFunctionMin)->default_value(p_Options->m_InitialMassFunctionMin),                                                                          
            ("Minimum mass (in Msol) to generate using given IMF (default = " + std::to_string(p_Options->m_InitialMassFunctionMin) + ")").c_str()
        )
        (
            "initial-mass-power",                                          
            po::value<double>(&p_Options->m_InitialMassFunctionPower)->default_value(p_Options->m_InitialMassFunctionPower),                                                                      
            ("Single power law power to generate primary mass using given IMF (default = " + std::to_string(p_Options->m_InitialMassFunctionPower) + ")").c_str()
        )

        (
            "kick-direction-power",                                        
            po::value<double>(&p_Options->m_KickDirectionPower)->default_value(p_Options->m_KickDirectionPower),                                                                                  
            ("Power for power law kick direction distribution (default = " + std::to_string(p_Options->m_KickDirectionPower) + " = isotropic, +ve = polar, -ve = in plane)").c_str()
        )
        (
            "kick-magnitude-max",                                          
            po::value<double>(&p_Options->m_KickMagnitudeDistributionMaximum)->default_value(p_Options->m_KickMagnitudeDistributionMaximum),                                                      
            ("Maximum drawn kick magnitude in km s^-1. Ignored if < 0. Must be > 0 if using kick-magnitude-distribution=FLAT (default = " + std::to_string(p_Options->m_KickMagnitudeDistributionMaximum) + ")").c_str()
        )
        (
            "kick-magnitude",                                          
            po::value<double>(&p_Options->m_KickMagnitude)->default_value(p_Options->m_KickMagnitude),                                                      
            ("The magnitude of the kick velocity the star receives during the a supernova (default = " + std::to_string(p_Options->m_KickMagnitude) + " km s^-1 )").c_str()
        )
        (
            "kick-magnitude-1",                                          
            po::value<double>(&p_Options->m_KickMagnitude1)->default_value(p_Options->m_KickMagnitude1),                                                      
            ("The magnitude of the kick velocity the primary star receives during the a supernova (default = " + std::to_string(p_Options->m_KickMagnitude1) + " km s^-1 )").c_str()
        )
        (
            "kick-magnitude-2",                                          
            po::value<double>(&p_Options->m_KickMagnitude2)->default_value(p_Options->m_KickMagnitude2),                                                      
            ("The magnitude of the kick velocity the secondary star receives during the a supernova (default = " + std::to_string(p_Options->m_KickMagnitude2) + " km s^-1 )").c_str()
        )
        (
            "kick-magnitude-random",                                          
            po::value<double>(&p_Options->m_KickMagnitudeRandom)->default_value(p_Options->m_KickMagnitudeRandom),                                                      
            "Number used to choose the kick velocity magnitude for the star during the a supernova (default = uniform random number [0.0, 1.0))"
        )
        (
            "kick-magnitude-random-1",                                          
            po::value<double>(&p_Options->m_KickMagnitudeRandom1)->default_value(p_Options->m_KickMagnitudeRandom1),                                                      
            "Number used to choose the kick velocity magnitude for the primary star during the a supernova (default = uniform random number [0.0, 1.0))"
        )
        (
            "kick-magnitude-random-2",                                          
            po::value<double>(&p_Options->m_KickMagnitudeRandom2)->default_value(p_Options->m_KickMagnitudeRandom2),                                                      
            "Number used to choose the kick velocity magnitude for the secondary during the a supernova (default = uniform random number [0.0, 1.0))"
        )
        (
            "kick-magnitude-sigma-CCSN-BH",                                
            po::value<double>(&p_Options->m_KickMagnitudeDistributionSigmaCCSN_BH)->default_value(p_Options->m_KickMagnitudeDistributionSigmaCCSN_BH),                                            
            ("Sigma for chosen kick magnitude distribution for black holes (default = " + std::to_string(p_Options->m_KickMagnitudeDistributionSigmaCCSN_BH) + " km s^-1 )").c_str()
        )
        (
            "kick-magnitude-sigma-CCSN-NS",                                
            po::value<double>(&p_Options->m_KickMagnitudeDistributionSigmaCCSN_NS)->default_value(p_Options->m_KickMagnitudeDistributionSigmaCCSN_NS),                                            
            ("Sigma for chosen kick magnitude distribution for neutron stars (default = " + std::to_string(p_Options->m_KickMagnitudeDistributionSigmaCCSN_NS) + " km s^-1 )").c_str()
        )
        (
            "kick-magnitude-sigma-ECSN",                                   
            po::value<double>(&p_Options->m_KickMagnitudeDistributionSigmaForECSN)->default_value(p_Options->m_KickMagnitudeDistributionSigmaForECSN),                                            
            ("Sigma for chosen kick magnitude distribution for ECSN (default = " + std::to_string(p_Options->m_KickMagnitudeDistributionSigmaForECSN) + " km s^-1 )").c_str()
        )
        (
            "kick-magnitude-sigma-USSN",                                   
            po::value<double>(&p_Options->m_KickMagnitudeDistributionSigmaForUSSN)->default_value(p_Options->m_KickMagnitudeDistributionSigmaForUSSN),                                            
            ("Sigma for chosen kick magnitude distribution for USSN (default = " + std::to_string(p_Options->m_KickMagnitudeDistributionSigmaForUSSN) + " km s^-1 )").c_str()
        )
        (
            "kick-mean-anomaly-1",
            po::value<double>(&p_Options->m_KickMeanAnomaly1)->default_value(p_Options->m_KickMeanAnomaly1),                                                                                  
            "Mean anomaly for the primary star at instantaneous time of the supernova (default = uniform random number [0.0, 2pi))"
        )
        (
            "kick-mean-anomaly-2",
            po::value<double>(&p_Options->m_KickMeanAnomaly2)->default_value(p_Options->m_KickMeanAnomaly2),                                                                                  
            "Mean anomaly for the secondary star at instantaneous time of the supernova (default = uniform random number [0.0, 2pi))"
        )
        (
            "kick-phi-1",
            po::value<double>(&p_Options->m_KickPhi1)->default_value(p_Options->m_KickPhi1),                                                                                  
            "Angle between 'x' and 'y', both in the orbital plane of the supernovae vector, for the primary star (default = drawn from kick direction distribution)"
        )
        (
            "kick-phi-2",
            po::value<double>(&p_Options->m_KickPhi2)->default_value(p_Options->m_KickPhi2),                                                                                  
            "Angle between 'x' and 'y', both in the orbital plane of the supernovae vector, for the secondary star (default = drawn from kick direction distribution)"
        )
        (
            "kick-scaling-factor",                                         
            po::value<double>(&p_Options->m_KickScalingFactor)->default_value(p_Options->m_KickScalingFactor),                                                                                    
            ("Arbitrary factor used to scale kicks (default = " + std::to_string(p_Options->m_KickScalingFactor) + ")").c_str()
        )
        (
            "kick-theta-1",                                        
            po::value<double>(&p_Options->m_KickTheta1)->default_value(p_Options->m_KickTheta1),                                                                                  
            "Angle between the orbital plane and the 'z' axis of the supernovae vector, for the primary star (default = drawn from kick direction distribution)"
        )
        (
            "kick-theta-2",                                        
            po::value<double>(&p_Options->m_KickTheta2)->default_value(p_Options->m_KickTheta2),                                                                                  
            "Angle between the orbital plane and the 'z' axis of the supernovae vector, for the secondary star (default = drawn from kick direction distribution)"
        )

        (
            "luminous-blue-variable-multiplier",                           
            po::value<double>(&p_Options->m_LuminousBlueVariableFactor)->default_value(p_Options->m_LuminousBlueVariableFactor),                                                                  
            ("Multiplicitive constant for LBV mass loss (default = " + std::to_string(p_Options->m_LuminousBlueVariableFactor) + ", use 10 for Mennekens & Vanbeveren 2014)").c_str()
        )

        (
            "mass-ratio-max",                                              
            po::value<double>(&p_Options->m_MassRatioDistributionMax)->default_value(p_Options->m_MassRatioDistributionMax),                                                                      
            ("Maximum mass ratio m2/m1 to generate (default = " + std::to_string(p_Options->m_MassRatioDistributionMax) + ")").c_str()
        )
        (
            "mass-ratio-min",                                              
            po::value<double>(&p_Options->m_MassRatioDistributionMin)->default_value(p_Options->m_MassRatioDistributionMin),                                                                      
            ("Minimum mass ratio m2/m1 to generate (default = " + std::to_string(p_Options->m_MassRatioDistributionMin) + ")").c_str()
        )
        (
            "mass-transfer-fa",                                            
            po::value<double>(&p_Options->m_MassTransferFractionAccreted)->default_value(p_Options->m_MassTransferFractionAccreted),                                                              
            ("Mass Transfer fraction accreted in FIXED prescription (default = " + std::to_string(p_Options->m_MassTransferFractionAccreted) + ", fully conservative)").c_str()
        )
        (
            "mass-transfer-jloss",                                         
            po::value<double>(&p_Options->m_MassTransferJloss)->default_value(p_Options->m_MassTransferJloss),                                                                                    
            ("Specific angular momentum with which the non-accreted system leaves the system (default = " + std::to_string(p_Options->m_MassTransferJloss) + ")").c_str()
        )
        (
            "mass-transfer-thermal-limit-C",                               
            po::value<double>(&p_Options->m_MassTransferCParameter)->default_value(p_Options->m_MassTransferCParameter),                                                                          
            ("Mass Transfer Thermal rate factor fo the accretor (default = " + std::to_string(p_Options->m_MassTransferCParameter) + ")").c_str()
        )
        (
            "maximum-evolution-time",                                      
            po::value<double>(&p_Options->m_MaxEvolutionTime)->default_value(p_Options->m_MaxEvolutionTime),                                                                                      
            ("Maximum time to evolve binaries in Myrs (default = " + std::to_string(p_Options->m_MaxEvolutionTime) + ")").c_str()
        )
        (
            "maximum-mass-donor-Nandez-Ivanova",                           
            po::value<double>(&p_Options->m_MaximumMassDonorNandezIvanova)->default_value(p_Options->m_MaximumMassDonorNandezIvanova),                                                            
            ("Maximum donor mass allowed for the revised common envelope formalism in Msol (default = " + std::to_string(p_Options->m_MaximumMassDonorNandezIvanova) + ")").c_str()
        )
        (
            "maximum-neutron-star-mass",                                   
            po::value<double>(&p_Options->m_MaximumNeutronStarMass)->default_value(p_Options->m_MaximumNeutronStarMass),                                                                          
            ("Maximum mass of a neutron star (default = " + std::to_string(p_Options->m_MaximumNeutronStarMass) + ")").c_str()
        )
        (
            "MCBUR1",                                                      
            po::value<double>(&p_Options->m_mCBUR1)->default_value(p_Options->m_mCBUR1),                                                                                                          
            ("MCBUR1: Min core mass at BAGB to avoid fully degenerate CO core  (default = " + std::to_string(p_Options->m_mCBUR1) + ")").c_str()
        )
        (
            "metallicity,z",                                               
            po::value<double>(&p_Options->m_Metallicity)->default_value(p_Options->m_Metallicity),                                                                                                
            ("Metallicity to use (default " + std::to_string(p_Options->m_Metallicity) + " Zsol)").c_str()
        )
        (
            "minimum-secondary-mass",                                      
            po::value<double>(&p_Options->m_MinimumMassSecondary)->default_value(p_Options->m_MinimumMassSecondary),                                                                              
            ("Minimum mass of secondary to generate in Msol (default = " + std::to_string(p_Options->m_MinimumMassSecondary) + ")").c_str()
        )

        (
            "neutrino-mass-loss-bh-formation-value",                       
            po::value<double>(&p_Options->m_NeutrinoMassLossValueBH)->default_value(p_Options->m_NeutrinoMassLossValueBH),                                                                        
            ("Value corresponding to neutrino mass loss assumption (default = " + std::to_string(p_Options->m_NeutrinoMassLossValueBH) + ")").c_str()
        )

        (
            "orbital-period-max",                                          
            po::value<double>(&p_Options->m_PeriodDistributionMax)->default_value(p_Options->m_PeriodDistributionMax),                                                                            
            ("Maximum period in days to generate (default = " + std::to_string(p_Options->m_PeriodDistributionMax) + ")").c_str()
        )
        (
            "orbital-period-min",                                          
            po::value<double>(&p_Options->m_PeriodDistributionMin)->default_value(p_Options->m_PeriodDistributionMin),                                                                            
            ("Minimum period in days to generate (default = " + std::to_string(p_Options->m_PeriodDistributionMin) + ")").c_str()
        )

        (
            "PISN-lower-limit",                                            
            po::value<double>(&p_Options->m_PairInstabilityLowerLimit)->default_value(p_Options->m_PairInstabilityLowerLimit),                                                                    
            ("Minimum core mass for PISN (default = " + std::to_string(p_Options->m_PairInstabilityLowerLimit) + ")").c_str()
        )
        (
            "PISN-upper-limit",                                            
            po::value<double>(&p_Options->m_PairInstabilityUpperLimit)->default_value(p_Options->m_PairInstabilityUpperLimit),                                                                    
            ("Maximum core mass for PISN (default = " + std::to_string(p_Options->m_PairInstabilityUpperLimit) + ")").c_str()
        )
        (
            "PPI-lower-limit",                                             
            po::value<double>(&p_Options->m_PulsationalPairInstabilityLowerLimit)->default_value(p_Options->m_PulsationalPairInstabilityLowerLimit),                                              
            ("Minimum core mass for PPI (default = " + std::to_string(p_Options->m_PulsationalPairInstabilityLowerLimit) + ")").c_str()
        )
        (
            "PPI-upper-limit",                                             
            po::value<double>(&p_Options->m_PulsationalPairInstabilityUpperLimit)->default_value(p_Options->m_PulsationalPairInstabilityUpperLimit),                                              
            ("Maximum core mass for PPI (default = " + std::to_string(p_Options->m_PulsationalPairInstabilityUpperLimit) + ")").c_str()
        )
        (
            "pulsar-birth-magnetic-field-distribution-max",                
            po::value<double>(&p_Options->m_PulsarBirthMagneticFieldDistributionMax)->default_value(p_Options->m_PulsarBirthMagneticFieldDistributionMax),                                        
            ("Maximum (log10) pulsar birth magnetic field (default = " + std::to_string(p_Options->m_PulsarBirthMagneticFieldDistributionMax) + ")").c_str()
        )
        (
            "pulsar-birth-magnetic-field-distribution-min",                
            po::value<double>(&p_Options->m_PulsarBirthMagneticFieldDistributionMin)->default_value(p_Options->m_PulsarBirthMagneticFieldDistributionMin),                                        
            ("Minimum (log10) pulsar birth magnetic field) (default = " + std::to_string(p_Options->m_PulsarBirthMagneticFieldDistributionMin) + ")").c_str()
        )
        (
            "pulsar-birth-spin-period-distribution-max",                   
            po::value<double>(&p_Options->m_PulsarBirthSpinPeriodDistributionMax)->default_value(p_Options->m_PulsarBirthSpinPeriodDistributionMax),                                              
            ("Maximum pulsar birth spin period in ms (default = " + std::to_string(p_Options->m_PulsarBirthSpinPeriodDistributionMax) + ")").c_str()
        )
        (
            "pulsar-birth-spin-period-distribution-min",                   
            po::value<double>(&p_Options->m_PulsarBirthSpinPeriodDistributionMin)->default_value(p_Options->m_PulsarBirthSpinPeriodDistributionMin),                                              
            ("Minimum pulsar birth spin period in ms (default = " + std::to_string(p_Options->m_PulsarBirthSpinPeriodDistributionMin) + ")").c_str()
        )
        (
            "pulsar-magnetic-field-decay-massscale",                       
            po::value<double>(&p_Options->m_PulsarMagneticFieldDecayMassscale)->default_value(p_Options->m_PulsarMagneticFieldDecayMassscale),                                                    
            ("Mass scale on which magnetic field decays during accretion in solar masses (default = " + std::to_string(p_Options->m_PulsarMagneticFieldDecayMassscale) + ")").c_str()
        )
        (
            "pulsar-magnetic-field-decay-timescale",                       
            po::value<double>(&p_Options->m_PulsarMagneticFieldDecayTimescale)->default_value(p_Options->m_PulsarMagneticFieldDecayTimescale),                                                    
            ("Timescale on which magnetic field decays in Myrs (default = " + std::to_string(p_Options->m_PulsarMagneticFieldDecayTimescale) + ")").c_str()
        )
        (
            "pulsar-minimum-magnetic-field",                               
            po::value<double>(&p_Options->m_PulsarLog10MinimumMagneticField)->default_value(p_Options->m_PulsarLog10MinimumMagneticField),                                                        
            ("log10 of the minimum pulsar magnetic field in Gauss (default = " + std::to_string(p_Options->m_PulsarLog10MinimumMagneticField) + ")").c_str()
        )

        (
            "semi-major-axis,a",                              
            po::value<double>(&p_Options->m_SemiMajorAxis)->default_value(p_Options->m_SemiMajorAxis),                                                        
            ("Initial semi-major axis, a (default = " + std::to_string(p_Options->m_SemiMajorAxis) + ")").c_str()
        )        
        (
            "semi-major-axis-max",                                         
            po::value<double>(&p_Options->m_SemiMajorAxisDistributionMax)->default_value(p_Options->m_SemiMajorAxisDistributionMax),                                                              
            ("Maximum semi major axis in AU to generate (default = " + std::to_string(p_Options->m_SemiMajorAxisDistributionMax) + ")").c_str()
        )
        (
            "semi-major-axis-min",                                         
            po::value<double>(&p_Options->m_SemiMajorAxisDistributionMin)->default_value(p_Options->m_SemiMajorAxisDistributionMin),                                                              
            ("Minimum semi major axis in AU to generate (default = " + std::to_string(p_Options->m_SemiMajorAxisDistributionMin) + ")").c_str()
        )

        (
            "wolf-rayet-multiplier",                                       
            po::value<double>(&p_Options->m_WolfRayetFactor)->default_value(p_Options->m_WolfRayetFactor),                                                                                        
            ("Multiplicitive constant for WR winds (default = " + std::to_string(p_Options->m_WolfRayetFactor) + ")").c_str()
        )

        (
            "zeta-adiabatic-arbitrary",                                    
            po::value<double>(&p_Options->m_ZetaAdiabaticArbitrary)->default_value(p_Options->m_ZetaAdiabaticArbitrary),                                                                          
            ("Value of mass-radius exponent zeta adiabatic (default = " + std::to_string(p_Options->m_ZetaAdiabaticArbitrary) + ")").c_str()
        )
        (
            "zeta-main-sequence",                                          
            po::value<double>(&p_Options->m_ZetaMainSequence)->default_value(p_Options->m_ZetaMainSequence),                                                                                      
            ("Value of mass-radius exponent zeta on the main sequence (default = " + std::to_string(p_Options->m_ZetaMainSequence) + ")").c_str()
        )
        (
            "zeta-radiative-envelope-giant",                               
            po::value<double>(&p_Options->m_ZetaRadiativeEnvelopeGiant)->default_value(p_Options->m_ZetaRadiativeEnvelopeGiant),                                                                  
            ("Value of mass-radius exponent zeta for radiative envelope giants (default = " + std::to_string(p_Options->m_ZetaRadiativeEnvelopeGiant) + ")").c_str()
        )


        // string options - alphabetically

        (
            "black-hole-kicks",                                            
            po::value<string>(&p_Options->m_BlackHoleKicksOptionString)->default_value(p_Options->m_BlackHoleKicksOptionString),                                                                              
            ("Black hole kicks relative to NS kicks (options: [FULL, REDUCED, ZERO, FALLBACK], default = " + p_Options->m_BlackHoleKicksOptionString + ")").c_str()
        )

        (
            "case-bb-stability-prescription",                              
            po::value<string>(&p_Options->m_CaseBBStabilityPrescriptionString)->default_value(p_Options->m_CaseBBStabilityPrescriptionString),                                                    
            ("Case BB/BC mass transfer stability prescription (options: [ALWAYS_STABLE, ALWAYS_STABLE_ONTO_NSBH, TREAT_AS_OTHER_MT, ALWAYS_UNSTABLE], default = " + p_Options->m_CaseBBStabilityPrescriptionString + ")").c_str()
        )
        (
            "chemically-homogeneous-evolution",                            
            po::value<string>(&p_Options->m_CheString)->default_value(p_Options->m_CheString),                                                                                                    
            ("Chemically Homogeneous Evolution (options: [NONE, OPTIMISTIC, PESSIMISTIC], default = " + p_Options->m_CheString + ")").c_str()
        )
        (
            "common-envelope-lambda-prescription",                         
            po::value<string>(&p_Options->m_CommonEnvelopeLambdaPrescriptionString)->default_value(p_Options->m_CommonEnvelopeLambdaPrescriptionString),                                          
            ("CE lambda prescription (options: [LAMBDA_FIXED, LAMBDA_LOVERIDGE, LAMBDA_NANJING, LAMBDA_KRUCKOW, LAMBDA_DEWI], default = " + p_Options->m_CommonEnvelopeLambdaPrescriptionString + ")").c_str()
        )
        (
            "common-envelope-mass-accretion-prescription",                 
            po::value<string>(&p_Options->m_CommonEnvelopeMassAccretionPrescriptionString)->default_value(p_Options->m_CommonEnvelopeMassAccretionPrescriptionString),                            
            ("Assumption about whether NS/BHs can accrete mass during common envelope evolution (options: [ZERO, CONSTANT, UNIFORM, MACLEOD], default = " + p_Options->m_CommonEnvelopeMassAccretionPrescriptionString + ")").c_str()
        )
        
        (
            "eccentricity-distribution",                                 
            po::value<string>(&p_Options->m_EccentricityDistributionString)->default_value(p_Options->m_EccentricityDistributionString),                                                          
            ("Initial eccentricity distribution (options: [ZERO, FIXED, FLAT, THERMALISED, GELLER+2013], default = " + p_Options->m_EccentricityDistributionString + ")").c_str()
        )
        (
            "envelope-state-prescription",                                 
            po::value<string>(&p_Options->m_EnvelopeStatePrescriptionString)->default_value(p_Options->m_EnvelopeStatePrescriptionString),                                                        
            ("Prescription for whether the envelope is radiative or convective (options: [LEGACY, HURLEY, FIXED_TEMPERATURE], default = " + p_Options->m_EnvelopeStatePrescriptionString + ")").c_str()
        )

        (
            "fryer-supernova-engine",                                      
            po::value<string>(&p_Options->m_FryerSupernovaEngineString)->default_value(p_Options->m_FryerSupernovaEngineString),                                                                  
            ("If using Fryer et al 2012 fallback prescription, select between 'delayed' and 'rapid' engines (default = " + p_Options->m_FryerSupernovaEngineString + ")").c_str()
        )

        (
            "initial-mass-function,i",                                     
            po::value<string>(&p_Options->m_InitialMassFunctionString)->default_value(p_Options->m_InitialMassFunctionString),                                                                    
            ("Initial mass function (options: [SALPETER, POWERLAW, UNIFORM, KROUPA], default = " + p_Options->m_InitialMassFunctionString + ")").c_str()
        )

        (
            "kick-direction",                                              
            po::value<string>(&p_Options->m_KickDirectionDistributionString)->default_value(p_Options->m_KickDirectionDistributionString),                                                        
            ("Natal kick direction distribution (options: [ISOTROPIC, INPLANE, PERPENDICULAR, POWERLAW, WEDGE, POLES], default = " + p_Options->m_KickDirectionDistributionString + ")").c_str()
        )
        (
            "kick-magnitude-distribution",                                 
            po::value<string>(&p_Options->m_KickMagnitudeDistributionString)->default_value(p_Options->m_KickMagnitudeDistributionString),                                                        
            ("Natal kick magnitude distribution (options: [ZERO, FIXED, FLAT, MAXWELLIAN, BRAYELDRIDGE, MULLER2016, MULLER2016MAXWELLIAN, MULLERMANDEL], default = " + p_Options->m_KickMagnitudeDistributionString + ")").c_str()
        )

        (
            "mass-loss-prescription",                                      
            po::value<string>(&p_Options->m_MassLossPrescriptionString)->default_value(p_Options->m_MassLossPrescriptionString),                                                                  
            ("Mass loss prescription (options: [NONE, HURLEY, VINK], default = " + p_Options->m_MassLossPrescriptionString + ")").c_str()
        )
        (
            "mass-ratio-distribution,q",                                   
            po::value<string>(&p_Options->m_MassRatioDistributionString)->default_value(p_Options->m_MassRatioDistributionString),                                                                
            ("Initial mass ratio distribution for q=m2/m1 (options: [FLAT, DuquennoyMayor1991, SANA2012], default = " + p_Options->m_MassRatioDistributionString + ")").c_str()
        )
        (
            "mass-transfer-accretion-efficiency-prescription",             
            po::value<string>(&p_Options->m_MassTransferAccretionEfficiencyPrescriptionString)->default_value(p_Options->m_MassTransferAccretionEfficiencyPrescriptionString),                    
            ("Mass Transfer Accretion Efficiency prescription (options: [THERMAL, FIXED], default = " + p_Options->m_MassTransferAccretionEfficiencyPrescriptionString + ")").c_str()
        )
        (
            "mass-transfer-angular-momentum-loss-prescription",            
            po::value<string>(&p_Options->m_MassTransferAngularMomentumLossPrescriptionString)->default_value(p_Options->m_MassTransferAngularMomentumLossPrescriptionString),                    
            ("Mass Transfer Angular Momentum Loss prescription (options: [JEANS, ISOTROPIC, CIRCUMBINARY, ARBITRARY], default = " + p_Options->m_MassTransferAngularMomentumLossPrescriptionString + ")").c_str()
        )
        (
            "mass-transfer-rejuvenation-prescription",                     
            po::value<string>(&p_Options->m_MassTransferRejuvenationPrescriptionString)->default_value(p_Options->m_MassTransferRejuvenationPrescriptionString),                                  
            ("Mass Transfer Rejuvenation prescription (options: [NONE, STARTRACK], default = " + p_Options->m_MassTransferRejuvenationPrescriptionString + ")").c_str()
        )
        (
            "mass-transfer-thermal-limit-accretor",                        
            po::value<string>(&p_Options->m_MassTransferThermallyLimitedVariationString)->default_value(p_Options->m_MassTransferThermallyLimitedVariationString),                                
            ("Mass Transfer Thermal Accretion limit (default = " + p_Options->m_MassTransferThermallyLimitedVariationString + ")").c_str()
        )

        (
            "neutrino-mass-loss-bh-formation",                             
            po::value<string>(&p_Options->m_NeutrinoMassLossAssumptionBHString)->default_value(p_Options->m_NeutrinoMassLossAssumptionBHString),                                                  
            ("Assumption about neutrino mass loss during BH formation (options: [FIXED_FRACTION, FIXED_MASS], default = " + p_Options->m_NeutrinoMassLossAssumptionBHString + ")").c_str()
        )
        (
            "neutron-star-equation-of-state",                              
            po::value<string>(&p_Options->m_NeutronStarEquationOfStateString)->default_value(p_Options->m_NeutronStarEquationOfStateString),                                                      
            ("Neutron star equation of state to use (options: [SSE, ARP3], default = " + p_Options->m_NeutronStarEquationOfStateString + ")").c_str()
        )

        (
            "pulsar-birth-magnetic-field-distribution",                    
            po::value<string>(&p_Options->m_PulsarBirthMagneticFieldDistributionString)->default_value(p_Options->m_PulsarBirthMagneticFieldDistributionString),                                  
            ("Pulsar Birth Magnetic Field distribution (options: [ZERO, FIXED, FLATINLOG, UNIFORM, LOGNORMAL], default = " + p_Options->m_PulsarBirthMagneticFieldDistributionString + ")").c_str()
        )
        (
            "pulsar-birth-spin-period-distribution",                       
            po::value<string>(&p_Options->m_PulsarBirthSpinPeriodDistributionString)->default_value(p_Options->m_PulsarBirthSpinPeriodDistributionString),                                        
            ("Pulsar Birth Spin Period distribution (options: [ZERO, FIXED, UNIFORM, NORMAL], default = " + p_Options->m_PulsarBirthSpinPeriodDistributionString + ")").c_str()
        )
        (
            "pulsational-pair-instability-prescription",                   
            po::value<string>(&p_Options->m_PulsationalPairInstabilityPrescriptionString)->default_value(p_Options->m_PulsationalPairInstabilityPrescriptionString),                              
            ("Pulsational Pair Instability prescription (options: [COMPAS, STARTRACK, MARCHANT], default = " + p_Options->m_PulsationalPairInstabilityPrescriptionString + ")").c_str()
        )

        (
            "remnant-mass-prescription",                                   
            po::value<string>(&p_Options->m_RemnantMassPrescriptionString)->default_value(p_Options->m_RemnantMassPrescriptionString),                                                            
            ("Choose remnant mass prescription (options: [HURLEY2000, BELCZYNSKI2002, FRYER2012, MULLER2016, MULLERMANDEL], default = " + p_Options->m_RemnantMassPrescriptionString + ")").c_str()
        )
        (
            "rotational-velocity-distribution",                            
            po::value<string>(&p_Options->m_RotationalVelocityDistributionString)->default_value(p_Options->m_RotationalVelocityDistributionString),                                              
            ("Initial rotational velocity distribution (options: [ZERO, HURLEY, VLTFLAMES], default = " + p_Options->m_RotationalVelocityDistributionString + ")").c_str()
        )

        (
            "semi-major-axis-distribution",                              
            po::value<string>(&p_Options->m_SemiMajorAxisDistributionString)->default_value(p_Options->m_SemiMajorAxisDistributionString),                                                        
            ("Initial semi-major axis distribution (options: [FLATINLOG, CUSTOM, DuquennoyMayor1991, SANA2012], default = " + p_Options->m_SemiMajorAxisDistributionString + ")").c_str()
        )        
        (
            "stellar-zeta-prescription",                                   
            po::value<string>(&p_Options->m_StellarZetaPrescriptionString)->default_value(p_Options->m_StellarZetaPrescriptionString),                                                            
            ("Prescription for stellar zeta (default = " + p_Options->m_StellarZetaPrescriptionString + ")").c_str()
        )
   
        ;   // end the list of options to be added

    } catch (po::error& e) {    // program options exception
        ok = false;             // set status
    } catch (...) {             // unhandled exception
        ok = false;             // set status
    }

    return ok;
}



/* DOCUMENTATION!!!!!!!!!!!!!!!!
 *
 * 
 * 
 * dataType is overall type (INT, FLOAT, STRING)
 * typeStr is detailed data type ("UNSIGNED LONG INT" etc)
 */
std::tuple<TYPENAME, bool, std::string, std::string> Options::OptionAttributes(const po::variables_map p_VM, const po::variables_map::const_iterator p_IT) {
            
    TYPENAME    dataType  = TYPENAME::NONE;
    std::string typeStr   = "";
    bool        defaulted = false;
    std::string valueStr  = "";

    if (((boost::any)p_IT->second.value()).empty()) return std::make_tuple(TYPENAME::NONE, true, "", "");    // empty option 

    // determine if option values was supplied, or whether the default was used

    defaulted = (p_VM[p_IT->first].defaulted() || p_IT->second.defaulted());

    // find data type and format the option value into a string
    // handles most data types - add others if they cause problems

    bool isCharPtr = false;
    bool isStr     = false;

    // (pre)check for data type = charPtr
    try {
        boost::any_cast<const char *>(p_IT->second.value());
        isCharPtr = true;
    } catch (const boost::bad_any_cast &) {
        isCharPtr = false;
    }

    if (!isCharPtr) {
        // (pre)check for data type = string
        try {
            boost::any_cast<std::string>(p_IT->second.value());
            isStr = true;
        } catch (const boost::bad_any_cast &) {
            isStr = false;
        }
    }

    // find other data types
    // it's not pretty, but it works

    if (isCharPtr) { 
        dataType = TYPENAME::NONE;                                          // not supported by COMPAS as an option data type
        typeStr  = "CONST_CHAR_*";                                          // ... but we know what type it is, and
        valueStr = p_VM[p_IT->first].as<const char *>();                    // ... we can still format the value
    }

    else if (isStr) {
        dataType = TYPENAME::STRING;
        typeStr  = "STRING";
        std::string tmp = p_VM[p_IT->first].as<std::string>();
        if (tmp.size()) valueStr = "'" + tmp + "'";
        else            valueStr = "''";
    }

    else if (((boost::any)p_IT->second.value()).type() == typeid(signed                )) { dataType = TYPENAME::INT;   typeStr = "SIGNED";                 valueStr = std::to_string(p_VM[p_IT->first].as<signed                >()); }
    else if (((boost::any)p_IT->second.value()).type() == typeid(unsigned              )) { dataType = TYPENAME::INT;   typeStr = "UNSIGNED";               valueStr = std::to_string(p_VM[p_IT->first].as<unsigned              >()); }

    else if (((boost::any)p_IT->second.value()).type() == typeid(short                 )) { dataType = TYPENAME::INT;   typeStr = "SHORT";                  valueStr = std::to_string(p_VM[p_IT->first].as<short                 >()); }
    else if (((boost::any)p_IT->second.value()).type() == typeid(signed short          )) { dataType = TYPENAME::INT;   typeStr = "SIGNED_SHORT";           valueStr = std::to_string(p_VM[p_IT->first].as<signed short          >()); }
    else if (((boost::any)p_IT->second.value()).type() == typeid(unsigned short        )) { dataType = TYPENAME::INT;   typeStr = "UNSIGNED_SHORT";         valueStr = std::to_string(p_VM[p_IT->first].as<unsigned short        >()); }

    else if (((boost::any)p_IT->second.value()).type() == typeid(short int             )) { dataType = TYPENAME::INT;   typeStr = "SHORT_INT";              valueStr = std::to_string(p_VM[p_IT->first].as<short int             >()); }
    else if (((boost::any)p_IT->second.value()).type() == typeid(signed short int      )) { dataType = TYPENAME::INT;   typeStr = "SIGNED_SHORT_INT";       valueStr = std::to_string(p_VM[p_IT->first].as<signed short int      >()); }
    else if (((boost::any)p_IT->second.value()).type() == typeid(unsigned short int    )) { dataType = TYPENAME::INT;   typeStr = "UNSIGNED_SHORT_INT";     valueStr = std::to_string(p_VM[p_IT->first].as<unsigned short int    >()); }

    else if (((boost::any)p_IT->second.value()).type() == typeid(int                   )) { dataType = TYPENAME::INT;   typeStr = "INT";                    valueStr = std::to_string(p_VM[p_IT->first].as<int                   >()); }
    else if (((boost::any)p_IT->second.value()).type() == typeid(signed int            )) { dataType = TYPENAME::INT;   typeStr = "SIGNED_INT";             valueStr = std::to_string(p_VM[p_IT->first].as<signed int            >()); }
    else if (((boost::any)p_IT->second.value()).type() == typeid(unsigned int          )) { dataType = TYPENAME::INT;   typeStr = "UNSIGNED_INT";           valueStr = std::to_string(p_VM[p_IT->first].as<unsigned int          >()); }

    else if (((boost::any)p_IT->second.value()).type() == typeid(long                  )) { dataType = TYPENAME::INT;   typeStr = "LONG";                   valueStr = std::to_string(p_VM[p_IT->first].as<long                  >()); }
    else if (((boost::any)p_IT->second.value()).type() == typeid(signed long           )) { dataType = TYPENAME::INT;   typeStr = "SIGNED_LONG";            valueStr = std::to_string(p_VM[p_IT->first].as<signed long           >()); }
    else if (((boost::any)p_IT->second.value()).type() == typeid(unsigned long         )) { dataType = TYPENAME::INT;   typeStr = "UNSIGNED_LONG";          valueStr = std::to_string(p_VM[p_IT->first].as<unsigned long         >()); }

    else if (((boost::any)p_IT->second.value()).type() == typeid(long int              )) { dataType = TYPENAME::INT;   typeStr = "LONG_INT";               valueStr = std::to_string(p_VM[p_IT->first].as<long int              >()); }
    else if (((boost::any)p_IT->second.value()).type() == typeid(signed long int       )) { dataType = TYPENAME::INT;   typeStr = "SIGNED_LONG_INT";        valueStr = std::to_string(p_VM[p_IT->first].as<signed long int       >()); }
    else if (((boost::any)p_IT->second.value()).type() == typeid(unsigned long int     )) { dataType = TYPENAME::INT;   typeStr = "UNSIGNED_LONG_INT";      valueStr = std::to_string(p_VM[p_IT->first].as<unsigned long int     >()); }

    else if (((boost::any)p_IT->second.value()).type() == typeid(long long             )) { dataType = TYPENAME::INT;   typeStr = "LONG_LONG";              valueStr = std::to_string(p_VM[p_IT->first].as<long long             >()); }
    else if (((boost::any)p_IT->second.value()).type() == typeid(signed long long      )) { dataType = TYPENAME::INT;   typeStr = "SIGNED_LONG_LONG";       valueStr = std::to_string(p_VM[p_IT->first].as<signed long long      >()); }
    else if (((boost::any)p_IT->second.value()).type() == typeid(unsigned long long    )) { dataType = TYPENAME::INT;   typeStr = "UNSIGNED_LONG_LONG";     valueStr = std::to_string(p_VM[p_IT->first].as<unsigned long long    >()); }

    else if (((boost::any)p_IT->second.value()).type() == typeid(long long int         )) { dataType = TYPENAME::INT;   typeStr = "LONG_LONG_INT";          valueStr = std::to_string(p_VM[p_IT->first].as<long long int         >()); }
    else if (((boost::any)p_IT->second.value()).type() == typeid(signed long long int  )) { dataType = TYPENAME::INT;   typeStr = "SIGNED_LONG_LONG_INT";   valueStr = std::to_string(p_VM[p_IT->first].as<signed long long int  >()); }
    else if (((boost::any)p_IT->second.value()).type() == typeid(unsigned long long int)) { dataType = TYPENAME::INT;   typeStr = "UNSIGNED_LONG_LONG_INT"; valueStr = std::to_string(p_VM[p_IT->first].as<unsigned long long int>()); }

    else if (((boost::any)p_IT->second.value()).type() == typeid(float                 )) { dataType = TYPENAME::FLOAT; typeStr = "FLOAT";                  valueStr = std::to_string(p_VM[p_IT->first].as<float                 >()); }
    else if (((boost::any)p_IT->second.value()).type() == typeid(double                )) { dataType = TYPENAME::FLOAT; typeStr = "DOUBLE";                 valueStr = std::to_string(p_VM[p_IT->first].as<double                >()); }
    else if (((boost::any)p_IT->second.value()).type() == typeid(long double           )) { dataType = TYPENAME::FLOAT; typeStr = "LONG_DOUBLE";            valueStr = std::to_string(p_VM[p_IT->first].as<long double           >()); }

    else if (((boost::any)p_IT->second.value()).type() == typeid(char                  )) { dataType = TYPENAME::INT;   typeStr = "CHAR";                   valueStr = std::to_string(p_VM[p_IT->first].as<char                  >()); }
    else if (((boost::any)p_IT->second.value()).type() == typeid(signed char           )) { dataType = TYPENAME::INT;   typeStr = "SIGNED_CHAR";            valueStr = std::to_string(p_VM[p_IT->first].as<signed char           >()); }
    else if (((boost::any)p_IT->second.value()).type() == typeid(unsigned char         )) { dataType = TYPENAME::INT;   typeStr = "UNSIGNED_CHAR";          valueStr = std::to_string(p_VM[p_IT->first].as<unsigned char         >()); }

    else if (((boost::any)p_IT->second.value()).type() == typeid(bool)) {
        dataType = TYPENAME::BOOL;
        typeStr  = "BOOL";
        valueStr = p_VM[p_IT->first].as<bool>() ? "TRUE" : "FALSE";
    } 

    else {  // Assume vector<string>
        try {
            std::ostringstream elemsSS;
            elemsSS << "{ ";
            vector<string> tmp = p_VM[p_IT->first].as<vector<string>>();
            for (std::vector<string>::iterator elem=tmp.begin(); elem != tmp.end(); elem++) {
                elemsSS << "'" << (*elem) << "', ";
            }
            string elems = elemsSS.str();
            if (elems.size() > 2) elems.erase(elems.size() - 2);
            else if (elems.size() == 2) elems.erase(elems.size() - 1);
            elems += " }";

            dataType = TYPENAME::NONE;                                                  // not supported by COMPAS as an option data type            
            typeStr  = "VECTOR<STRING>";                                                // ... but we know what type it is, and
            valueStr = elems;                                                           // ... we can still format the value
        } catch (const boost::bad_any_cast &) {
            dataType = TYPENAME::NONE;                                                  // unknown data type               
            typeStr  = "<UNKNOWN_DATA_TYPE>";
            valueStr = "<UNKNOWN_DATA_TYPE>";
        }
    }

    return std::make_tuple(dataType, defaulted, typeStr, valueStr);
}


/*
 * DOCUMENTATION HERE!
 * 
 */
string Options::ProgramOptionDetails(const OptionValues *p_Options, const po::variables_map p_VM) {
            
    TYPENAME    dataType  = TYPENAME::NONE;
    std::string typeStr   = "";
    bool        defaulted = false;
    std::string valueStr  = "";

    std::ostringstream ss;                                                                                              // output string

    ss << "COMMAND LINE OPTIONS\n-------------------\n\n";

    for (po::variables_map::const_iterator it = p_VM.begin(); it != p_VM.end(); it++) {                                 // for all options in the variable map
  
        ss << it->first << " = ";                                                                                       // add option name to output string

        std::tie(dataType, defaulted, typeStr, valueStr) = OptionAttributes(p_VM, it);                                  // get option attributes

        if (valueStr == "")                                                                                             // empty option?
            ss << "<EMPTY_OPTION>\n";                                                                                   // yes - say so
        else                                                                                                            // no
            ss << valueStr + ", " << (defaulted ? "DEFAULT_USED, " : "USER_SUPPLIED, ") << typeStr << "\n";             // add option details to output string
    }
  
    ss << "\n\nOTHER PARAMETERS\n----------------\n\n";

// JRFIX <<<<<<<<<<<<<<<<<<<<<<<<<<<< fix these!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ss << "fixedMetallicity   = " << (p_Options->m_FixedMetallicity ? "TRUE" : "FALSE") << ", CALCULATED, BOOL\n";     // fixedMetallicity
    ss << "useFixedUK         = " << (p_Options->m_UseFixedUK ? "TRUE" : "FALSE") << ", CALCULATED, BOOL\n";           // useFixedUK
    ss << "outputPath         = " << p_Options->m_OutputPath.string() << ", CALCULATED, STRING\n";                     // outputPath (fully qualified)
    ss << "fixedRandomSeed    = " << (p_Options->m_FixedRandomSeed ? "TRUE" : "FALSE") << ", CALCULATED, BOOL\n";      // fixedRandomSeed

    return ss.str();
}



            // before we give the options to boost we need to determine if the user passed
            // any ranges or sets and, if they did, handle those - boost doesn't know anything
            // about them
            //
            // a range is allowed only for numeric options (i.e. INT or FLOAT types),
            // but is not allowed for all numeric options (e.g. --log-level)
            // a set is allowed for numeric, string, and bool options - but not all
            // of them (e.g. --quiet)
            //
            // we define a vector of options excluded from the range and set constructs
            // (one vector each).  We don't need to exclude non-numeric options from range
            // here - that is done later - here we just exclude options for which range/set
            // makes no sense

std::string Options::ParseOptionValues(int p_ArgCount, char *p_ArgStrings[], OptionsDescriptorT &p_OptionsDescriptor) {

    bool error         = false;                                                                                     // for now...
    std::string errStr = "";                                                                                        // also for now...
    
    // JRFIX: we should use ERROR:: in constants.h for error strings here...!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    try {

        std::string optionName         = "";                                                                        // option name
        int type                       = -1;                                                                        // arg type (range = 0, set = 1, neither = -1
        std::vector<std::string> parms = {};                                                                        // the range or set parameters

        for (int iArg = 0; iArg < p_ArgCount; iArg++) {                                                             // for each arg string

            if (iArg == 0) continue;                                                                                // ignore the executable name

            optionName = p_ArgStrings[iArg - 1];                                                                    // get the option name for the argument we're processing
            if (optionName[0] == '-') optionName.erase(0, optionName.find_first_not_of("-"));                       // remove the "-" or "--"

            if (p_ArgStrings[iArg] != nullptr) {                                                                    // null arg?
                                                                                                                    // no
                std::string str(p_ArgStrings[iArg]);                                                                // convert char* to std::string
                str = utils::ToLower(utils::trim(str));                                                             // downshift - for comparisons

                // check for RANGE or SET
                // range is indicated by 'range[start,count,inc]', 'r[start,count,inc]', or just '[start,count,inc]'

                if ((str[0] == '[') || (str.rfind("r[", 0) == 0) || (str.rfind("range[", 0) == 0)) {                // starts with '[', 'r[' or 'range[', so...
                    error = true;                                                                                   // unless set otherwise
                    if (str[str.length()-1] == ']') {                                                               // ... needs to end with ']' to be a valid RANGE
                        type  = 0;                                                                                  // it did - so RANGE

                        // check for RANGE requested for option in range excluded list
                        if (iArg > 1) {                                                                             // range not valid for arg[1]
                            if (std::find(m_RangeExcluded.begin(), m_RangeExcluded.end(), optionName) != m_RangeExcluded.end())
                                errStr = std::string("argument range not supported for option '") + optionName + std::string("'");
                            else
                                error = false;                                                                      // we're good
                        }
                    }
                }

                if (!error) {                                                                                       // still ok?
                                                                                                                    // yes
                    // set is indicated by 'set[elem1,elem2,...,elemN]', or 's[elem1,elem2,...,elemN]'

                    if ((str.rfind("s[", 0) == 0) || (str.rfind("set[", 0) == 0)) {                                 // starts with 's[' or 'set[', so ...
                        error = true;                                                                               // unless set otherwise
                        if (str[str.length()-1] == ']') {                                                           // ... needs to end with ']' to be a valid SET
                            type  = 1;                                                                              // it did - so SET

                            // check for SET requested for option in set excluded list
                            if (iArg > 1) {                                                                         // set not valid for arg[1]
                                if (std::find(m_SetExcluded.begin(), m_SetExcluded.end(), optionName) != m_SetExcluded.end())
                                    errStr = std::string("argument set not supported for option '") + optionName + std::string("'");
                                else
                                    error = false;                                                                  // we're good
                            }
                        }
                    }
                }

                if (!error && type != -1) {                                                                         // range or set?
                                                                                                                    // yes
                    // we have what looks like a 'range' or 'set' argument
                    // for now, just stash the details away and substitute the
                    // first value for the argument so we can check parsing

                    // look for comma separated values - there will be no 
                    // spaces - the OS/shell would have complained...

                    if (str.rfind("range", 0) == 0) str.erase(0, 5);                                                // strip 'range' (range indicator) if present
                    if (str.rfind("set", 0) == 0) str.erase(0, 3);                                                  // strip 'set' (set indicator) if present
                    if (str[0] == 'r' || str[0] == 's') str.erase(0, 1);                                            // strip 'r' or 's' (range or set indicator) if present
                    str = str.substr(1, str.size() - 2);                                                            // strip enclosing brackets (must be present)

                    if (str.length() == 0 || str[str.length() - 1] == ',') error = true;                            // no values, or trailing comma is an error
                    else {

                        parms.clear();                                                                              // start empty

                        size_t start = 0;                                                                           // start position
                        size_t pos   = 0;                                                                           // current position
                        while (!error && start < str.length() && pos != string::npos) {                             // comma found before the end of the string?

                            std::string value = "";                                                                 // value

                            pos = str.find(",", start);                                                             // next comma
                                                                                                        
                            if ((pos - start) > 0) {                                                                // non-zero length string?
                                value = str.substr(start, pos - start);                                             // yes - grab it
                                parms.push_back(value);                                                             // store value

                                start = pos + 1;                                                                    // next start
                            }
                            else error = true;                                                                      // empty value - stop
                        }

                        if (!error) {                                                                               // still ok?
                                                                                                                    // yes
                            if (type == 0 && parms.size() != 3) {                                                   // ranges require exactly 3 parameters
                                error  = true;                                                                      // error
                                errStr = "argument range requires exactly three parameters"; 
                            }
                            else {
                                // if range, then we have 3 parameters (checked above)
                                // if set, we have at least one parameter (check earlier), so we're good

//                                p_OptionsDescriptor.complexOptionValues.push_back(std::make_tuple(optionName, std::make_tuple(type, parms, 0))); // store the range/set
                                RangeOrSetDescriptorT details = {type, parms, 0};
                                p_OptionsDescriptor.complexOptionValues.push_back(std::make_tuple(optionName, details)); // store the range/set

                                strncpy(p_ArgStrings[iArg], parms[0].c_str(), parms[0].length());                   // replace arg value (temporarily)
                                p_ArgStrings[iArg][parms[0].length()] = '\0';
                            }
                        }
                    }
                }          
            }
            if (error) break;                                                                                       // stop parsing if error encountered
        }

        // boost parse_command_line() expects the first arg to be the program name
        // (it thinks it is getting the values that were passed to main() from the 
        // OS/shell), so for options from a grid file we insert a dummy argument as 
        // arg[0] and set the argument count appropriately.
        //
        // if valid ranges or sets were specified by the user they've been temporariliy
        // replaced for the boost parse, but if they were not valid they've been left
        // in the argument strings that will be passed to boost - so boost will fail
        // and complain about the offending parameter (which is what we want)

        po::parsed_options const parsedOptions = po::parse_command_line(p_ArgCount, p_ArgStrings, p_OptionsDescriptor.optionDescriptions);  // parse user-supplied options
        po::store(parsedOptions, p_OptionsDescriptor.optionValues.m_VM);                                            // store parsed options into variable map

        // if we've made it this far then boost parsed the commandline arguments ok
        // if there were any ranges or sets specified by the user we can now work
        // out the data types of the options for which they (the ranges/sets) were
        // specified and sanity check them.
        //
        //
        // iterate through the specified ranges and sets to sanity check
        // need to check:
        //     - ranges have not bee specified for non-numeric options
        //     - range values are all numeric

// MAKE THIS A TYPEDEF - in .h as well!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        RangeOrSetDescriptorT details = {};                                                // initialised to please the compiler...

        for (size_t idx = 0; idx < p_OptionsDescriptor.complexOptionValues.size(); idx++) {                         // for each range or set specified

            error = false;                                                                                          // for now...

            parms.clear();                                                                                          // clear the range/set parameters




// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// FIGURE OUT HOW TO CHECK RANGE AND SET PARAMETERS AGINST VALID VALUES
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




            optionName = get<0>(p_OptionsDescriptor.complexOptionValues[idx]);                                      // the option name
            details    = get<1>(p_OptionsDescriptor.complexOptionValues[idx]);                                      // range/set details for this optionName
            type       = details.type;                                                                           // range = 0, set = 1
            parms      = details.parameters;                                                                           // range/set parameter values

            if (type == 0) {                                                                                        // range?
                po::variables_map::const_iterator it = p_OptionsDescriptor.optionValues.m_VM.find(optionName);      // yes - find the option in the boost variables map
                if (it != p_OptionsDescriptor.optionValues.m_VM.end()) {                                            // found?
                    TYPENAME dataType;                                                                              // yes
                    std::tie(dataType, std::ignore, std::ignore, std::ignore) = OptionAttributes(p_OptionsDescriptor.optionValues.m_VM, it);    // get data type
                    if (dataType != TYPENAME::INT && dataType != TYPENAME::FLOAT) {                                 // numeric?
                        error  = true;                                                                              // no - that's not ok
                        errStr = std::string("argument range not supported for option '") + optionName + std::string("'");
                    }
                    else {                                                                                          // yes - numeric
                        for (size_t ip = 0; ip < parms.size(); ip++) {                                              // for each range parameter specified
                            if (parms[ip].empty() || !std::all_of(parms[ip].begin(), parms[ip].end(), ::isdigit)) { // numeric?
                                error  = true;                                                                      // no - that's not ok
                                errStr = std::string("all parameters of argument range must be numeric for option '") + optionName  + std::string("'");
                                break;
                            }
                        }


                        if (!error) {

                                // CALCULATE FIRST AND LAST (=FIRST + (COUNT * INC)) AND RANEG CHECK AGANIST MIN & MAX FOR OPTION (MAY BE 0 AND INFINITE)
                        }


                    }
                }
                else {                                                                                              // option not found in boost variables map
                    error  = true;                                                                                  // that can't be good...
                    errStr = std::string("internal error: boost vm, option '") + optionName + std::string("'");
                }
            }



            else {      // SET

                // CHECK SET VALUES FOR OPTION ARE VALID - CF CHECKANDSET()

            }



            if (error) break;                                                                                       // stop now

        }
    
    } catch (po::error& e) {                                                                                        // program options exception
        errStr = e.what();                                                                                          // set error string
    } catch (...) {                                                                                                 // unhandled exception
        errStr = "uhandled exception";                                                                              // set errors tring
    }

    return errStr;
}


int Options::GetProgramOptionValues() {

    m_Opts = &m_Program;                        // point at program options

    // upon entry iterators for ranges and sets will be pointing at the
    // values that should be returned - so those values need to be copied
    // into place, then the iterators advanced as required



     


    return 0;

}



/*
 * Initialise options service
 * 
 * Intitialises the options service.  Constructs options objects for the program options
 * (options that are specified only on the commandline and that cannot be specified in a 
 * grid file on a per object (star/binary) basis), and the grid file options (options that
 * can be specified in a grid file on a per object (star/binary) basis).  Populates the
 * program options object from the commandline arguments passed to main() - this object
 * stays static throughout the life of the program.
 * 
 * bool Options::Initialise(int p_ArgCount, char *p_ArgStrings[])
 * 
 * @param   [IN]    p_ArgCount                  Integer number of args passed in p_ArgStrings
 * @param   [IN]    p_ArgStrings                Arg strings - 1 per option
 *                                              Note that the first arg string is ignored (expected to be program name)
 * @return                                      Boolean status (true = ok, false = problem) (this function displays any error string)
 */
bool Options::Initialise(int p_ArgCount, char *p_ArgStrings[]) {

    bool ok = true;                                                                                                         // status - unless something changes

    try {
        m_Program.optionsServed = false;                                                                                    // not yet

        m_Program.optionValues.Initialise();                                                                                // initialise option variables for program options
        m_EvolvingObject.optionValues.Initialise();                                                                         // initialise option variables for evolving object options

        // initialise program options - these are the options that are specified only on the 
        // commandline and that cannot be specified in a grid file on a per evolving object (star/binary) basis

        po::options_description programLevelOptions("Program Level Options");                                               // boost options descriptions object for program-level options
        ok = SetProgramOptions(&m_Program.optionValues, &programLevelOptions);                                              // ... populate
        if (!ok) {                                                                                                          // ok?
            COMPLAIN("failed to initialise options descriptions for program-level options");                                // no, complain - this throws an exception
        }
        else {                                                                                                              // yes, ok
            // initialise evolving object options - these are the options that can be specified in a grid file
            // on a per evolving object (star/binary) basis.  The values of options specified in a grid file
            // take precedence over the values of the same options specified on the commandline, but only for
            // the object (star/binary) corresponding to the grid file record.

            po::options_description objectLevelOptions("Star/Binary Level Options");                                        // boost options descriptions object for evolving object (star/binary) level options
            ok = SetObjectOptions(&m_EvolvingObject.optionValues, &objectLevelOptions);                                     // ... populate
            if (!ok) {                                                                                                      // ok?
                COMPLAIN("failed to initialise options descriptions for evolving object-level options");                    // no, complain - this throws an exception
            }
            else {                                                                                                          // yes, ok
                // populate the commandline options
                // this stays static throughout the life of the program
                m_Program.optionDescriptions.add(programLevelOptions).add(objectLevelOptions);                              // both program and object options are available on the commandline
    
                // we parse the option values before handing them over to boost
                // boost knows nothing about ranges and sets, so we have to handlde
                // them ourselves first
                m_Program.complexOptionValues = {};                                                                         // no ranges or sets - unless we find them in the parse
                std::string errStr = ParseOptionValues(p_ArgCount, p_ArgStrings, m_Program);                                // parse the option values - specifically for ranges and sets
                if (!errStr.empty()) {                                                                                      // parsed ok?
                    COMPLAIN(errStr);                                                                                       // no, complain - this throws an exception
                }
                else {
                    errStr = m_Program.optionValues.CheckAndSetOptions();                                                   // yes - sanity check, and set, values
                    if (!errStr.empty()) {                                                                                  // check ok?
                        COMPLAIN(errStr);                                                                                   // no, complain - this throws an exception
                    }
                    else {

                        m_OptionsDetails = ProgramOptionDetails(&m_Program.optionValues, m_Program.optionValues.m_VM);      // yes - get Run_Details contents

                        m_Opts = &m_Program;                                                                                // point at program options

                        // We now have the options the user entered at the commandline, including any ranges and/or 
                        // sets, so this is where we stop the initialisation - from here we just play out the options 
                        // that are specified by any ranges and sets via the GetProgramOptionValues() function.
                        //
                        // If the user has specified any ranges or sets we set the options to the first value in each
                        // range or set (already done by the time we get here).  Calls to GetProgramOptionValues() will
                        // then advance the option values through the ranges and sets as required - *however*, the very
                        // first call to GetProgramOptionValues() just returns the options as they are set here, so that 
                        // a call to GetProgramOptionValues() can be put in a loop, and the first evaluation of the loop
                        // will be the first star/binary - and if only 1 star/binary is required then no loop, and no call 
                        // to GetProgramOptionValues() is required (because the initial values have already been set).
                        //
                        // Note that there are analogous functikns for object (star/binary) initialisation and retrievel
                        // option values: InitialiseObject() and GetObjectOptionValues().  These functions intitialise
                        // and retrieve options specified in grid file records.
                    }
                }
            }
        }  
    } catch (po::error& e) {                                                                                                // program options exception
        std::cerr << "Program Options error: " << e.what() << std::endl;                                                    // show the problem
        std::cerr << m_Program.optionDescriptions << std::endl;                                                             // show help
        ok = false;                                                                                                         // set status
    }  catch (const std::string eStr) {                                                                                     // custom exception
        std::cerr << "Program Options error: " << eStr << std::endl;                                                        // show the problem
        std::cerr << m_Program.optionDescriptions << std::endl;                                                             // show help
        ok = false;                                                                                                         // set status
    } catch (...) {                                                                                                         // unhandled exception
        std::cerr << "Program Options error: unhandled exception" << std::endl;                                             // show the problem
        std::cerr << m_Program.optionDescriptions << std::endl;                                                             // show help
        ok = false;                                                                                                         // set status
    }

    return ok;
}


/*
 * Initialise grid file options
 * 
 * Intitialises the grid file options (options that can be specified in a grid file on a 
 * per object (star/binary) basis).  Populates the object options object (created by
 * Initialise()) from the arguments specified in the grid file - this object is updated
 * for each grid file record.
 * 
 * bool Options::InitialiseObject(int p_ArgCount, char *p_ArgStrings[])
 * 
 * @param   [IN]    p_OptionsString             String containg all options - the grid file record
 * @return                                      String containing an error string
 *                                              If no error occurred the return string will be the empty string 
 */
bool Options::InitialiseObject(const std::string p_OptionsString) {

    bool ok = true;                                                                                                         // status - unless something changes

    po::options_description gridfileOptions;                                                                                // boost options descriptions object for options available on commandline

    try {

        // parse the option string (just as the OS/shell would do)

        std::vector<std::string> parsedStrings;                                                                             // parsed option strings

        size_t start      = 0;                                                                                              // start position of parsed option string
        size_t end        = 0;                                                                                              // end position of parsed option strinf
        std::string delim = " ";                                                                                            // delimiter
        while (end != string::npos) {                                                                                       // iterate over input string
            end = p_OptionsString.find(delim, start);                                                                       // find delimiter
            std::string s = p_OptionsString.substr(start, end - start);                                                     // grab option/argument string
            parsedStrings.push_back(utils::trim(s));                                                                        // trim whitespace and store
            start = end + delim.length();                                                                                   // new start position
        }
    
        std::vector<char const*> args {"placeHolder"};                                                                      // place-holder - boost expects command name as argv[0]
        for (auto& arg : parsedStrings)                                                                                     // iterate over the parsed strings
            args.push_back(arg.c_str());                                                                                    // and grab the c_str  

        m_EvolvingObject.optionValues = m_Program.optionValues;                                                             // start with program/commandline options

        // initialise object options - these are the options that can be specified in a grid 
        // file on a per object (star/binary) basis.  The values of options specified in a grid 
        // file take precedence over the values of the same options specified on the commandline,
        // but only for the object (star/binary) corresponding to the grid file record.

        po::options_description objectLevelOptions("Star/Binary Level Options");                                            // boost options descriptions object for per object (star/binary) options
        ok = SetObjectOptions(&m_EvolvingObject.optionValues, &objectLevelOptions);                                         // ... populate
        if (!ok) {                                                                                                          // ok?
            COMPLAIN("Failed to initialise program options descriptions for ObjectOptions");                                // no, complain - this throws an exception
        }
        else {                                                                                                              // yes, ok

            // populate the grid file options
            // this changes for every record in the grid file

            m_EvolvingObject.optionDescriptions.add(objectLevelOptions);                                                    // object options only
            po::parsed_options const parsedOptions = po::parse_command_line(args.size(), args.data(), gridfileOptions);     // parse user-supplied options

            po::store(parsedOptions, m_EvolvingObject.optionValues.m_VM);                                                   // store parsed options into variable map
            po::notify(m_EvolvingObject.optionValues.m_VM);                                                                 // populate the variables with option values

            m_EvolvingObject.optionValues.CheckAndSetOptions();                                                             // check user-supplied options and set values as appropriate
    
            m_Opts = &m_EvolvingObject;   // REVISIT THIS!!!!!!!!!!!!!!!!!                                                                                  // point at object options
        }

    } catch (po::error& e) {                                                                                                // program options exception
        std::cerr << "Program Options error: " << e.what() << std::endl;
        std::cerr << gridfileOptions << std::endl;
        ok = false;                                                                                                         // set status
    } catch (...) {                                                                                                           // unhandled exception - something wrong...
        std::cerr << "Program Options error: unhandled exception" << std::endl;
        std::cerr << gridfileOptions << std::endl;
        ok = false;                                                                                                         // set status
    }
    
    return ok;
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

    bool ok = true;                                                                                                     // status - unless a problem occurs

    COMPAS_VARIABLE_TYPE value;                                                                                         // default property value

    PROGRAM_OPTION property = boost::get<PROGRAM_OPTION>(p_Property);                                                   // get property
                                                                                                                        // get property value
    switch (property) {

        case PROGRAM_OPTION::KICK_MAGNITUDE_DISTRIBUTION_SIGMA_CCSN_BH:  value = KickMagnitudeDistributionSigmaCCSN_BH(); break;
        case PROGRAM_OPTION::KICK_MAGNITUDE_DISTRIBUTION_SIGMA_CCSN_NS:  value = KickMagnitudeDistributionSigmaCCSN_NS(); break;
        case PROGRAM_OPTION::KICK_MAGNITUDE_DISTRIBUTION_SIGMA_FOR_ECSN: value = KickMagnitudeDistributionSigmaForECSN(); break;
        case PROGRAM_OPTION::KICK_MAGNITUDE_DISTRIBUTION_SIGMA_FOR_USSN: value = KickMagnitudeDistributionSigmaForUSSN(); break;
        case PROGRAM_OPTION::RANDOM_SEED:                                value = RandomSeed();                            break;

        default:                                                                                                        // unknown property
            ok    = false;                                                                                              // that's not ok...
            value = "UNKNOWN";                                                                                          // default value
            std::cerr << ERR_MSG(ERROR::UNKNOWN_PROGRAM_OPTION) << std::endl;                                           // show warning (don't have logging or errors here...)
    }

    return std::make_tuple(ok, value);
}


/*
 * Read and apply the next record in the grid file
 * 
 * The record from the grid file is read as one string, then passed
 * to InitialiseObject() for processing.
 * 
 * In InitialiseObject() the record is parsed into separate tokens ready 
 * to be passed to the boost program option parser.  Once that is done
 * the options are handed over to the boost functions for parsing and
 * detting of values.
 * 
 * pon return from this function the option values will be set to the
 * values specified by the user, or their default values - either way
 * ready for the star/binary to be evolved.
 * 
 * The grid file struct (m_Gridfile) will be used and updated by this
 * function.  The grid file name, error status, and file handle are
 * stored in the struct.  
 * 
 * 
 * int ApplyNextGridRecord()
 * 
 * @return                                      Int result:
 *                                                  -1: Error reading grid file record (error value in grid file struct)
 *                                                   0: No record to read - end of file
 *                                                   1: Grid file record read and applied ok
 */
int Options::ApplyNextGridRecord() {

    int status = -1;                                        // default status is failure

    if (m_Gridfile.handle.is_open()) {                      // file open?
                                                            // yes
        std::string record;                                 // the record read
        std::getline(m_Gridfile.handle, record);            // read the next record
        if (m_Gridfile.handle.fail()) {                     // read ok?
            if (m_Gridfile.handle.eof()) status = 0;        // eof?
            else {                                          // no
                m_Gridfile.error = ERROR::FILE_READ_ERROR;  // record error
                status = -1;                                // set status
            }
        }
        else {                                              // read ok
            status = InitialiseObject(record) ? 1 : -1;     // apply record and set status
        }
    }

    return status;
}


/*
 * Open the grid file
 *
 * The grid file is opened at the start of the simulation and stays open until the
 * simulation is complete - we just pick a record off and process the record, and
 * when we hit the end of the file the file is closed and the simulation complete.
 *
 * 
 * ERROR OpenGridFile(const std::string p_GridFilename)
 *
 * @param   [IN]        p_Filename              The filename of the Grid file
 * @return                                      ERROR indicator - will be ERROR::NONE if file opened sccessfully
 */
ERROR Options::OpenGridFile(const std::string p_GridFilename) {

    m_Gridfile.filename = p_GridFilename;                       // record filename

    if (!m_Gridfile.filename.empty()) {                         // have grid filename?
        m_Gridfile.handle.open(m_Gridfile.filename);            // yes - open the file
        if (m_Gridfile.handle.fail()) {                         // open ok?
            m_Gridfile.error = ERROR::FILE_OPEN_ERROR;          // no - record error
        }
        else m_Gridfile.error = ERROR::NONE;                    // open ok - no error
    }
    else m_Gridfile.error = ERROR::EMPTY_FILENAME;              // empty filename

    return m_Gridfile.error;
}
