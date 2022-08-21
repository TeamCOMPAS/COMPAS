/******************************************************************************************/
/*                                                                                        */
/* Instructions for adding a new option:                                                  */
/*                                                                                        */
/* 1. Decide on a string for the option - this is the string the user will use on the     */
/*    commandline or in the grid file (e.g. "random-seed")                                */
/*    The convention I've settled on is hyphenated lower case - I don't mind what         */
/*    convention we settle on, as long as it's just one.                                  */
/*                                                                                        */
/* 2. Decide on a class member variable name for the option (e.g. m_RandomSeed).          */
/*                                                                                        */
/* 3. Decide on a default value for the option.                                           */
/*                                                                                        */
/* 4. Add the class member variable to the PRIVATE area of the OptionValues class in      */
/*    Options.h.                                                                          */
/*                                                                                        */
/* 5. Add a getter for the class member variable to the PUBLIC area of the Options class  */
/*    in Options.h.                                                                       */
/*                                                                                        */
/*    Decide if the getter should always retrieve the value specified on the commandline, */
/*    or whether it should retrieve the grid line value if one was specified by the user, */
/*    and only retrieve the commandline value if the user did not specify a grid line     */
/*    value - see the OPT_VALUE macro defined in Options.h.                               */
/*                                                                                        */
/* 6. Add the class member variable initialisation (to the default value) to the          */
/*    Options::OptionValues::Initialise() function in Options.cpp.                        */
/*                                                                                        */
/* 7. Add the option to the Options::AddOptions() function in Options.cpp.                */
/*    This is where we tell Boost about the option - the option string, the variable in   */
/*    which the user-specified value should be stored, the default value to use if the    */
/*    user does not specify a value, and the (text) description of the option.            */
/*                                                                                        */
/*    Options::AddOptions() is a bit of a beast - there's no easy, short way to specify   */
/*    the details of as many options as we have.  I have tried to make it semi-readable,  */
/*    but it is, and always will be, just long...  The best we can do is keep it neat so  */
/*    it doesn't become too hard to read.                                                 */
/*                                                                                        */
/*    When adding options to Options::AddOptions(), the convention I have used is that    */
/*    the option string (e.g. "random-seed") is predominantly in lower case - that's not  */
/*    strictly required, but I think it's easier for users to remember the option names   */
/*    if they don't have to remember a mixture of upper and lower case.  I have           */
/*    configured Boost to perform a case-insensitive match, so option strings can have    */
/*    mixed case if necessary (e.g. muller-mandel-kick-multiplier-BH).                    */
/*                                                                                        */
/* 8. Add any sanity checks: constraint/range/dependency checks etc. for the new option,  */
/*    and any affected existing options, to Options::OptionValues::CheckAndSetOptions()   */
/*    in Options.cpp.  It is also here you can set any final values that, perhaps due to  */
/*    dependencies on options that had not yet been parsed, could not be set directly by  */
/*    Boost when the options were parsed (also see SetCalculatedOptionDefaults(); viz.    */
/*    m_KickPhi1 etc.).                                                                   */
/*                                                                                        */
/* 9. Add the new option to one or more of the following vectors in Options.h, as         */
/*    required:                                                                           */
/*                                                                                        */
/*        m_ShorthandAllowed: options for which shorthand notation is allowed             */
/*                                                                                        */
/*        m_GridLineExcluded: option strings that may not be specified on a grid line     */
/*                                                                                        */
/*        m_SSEOnly         : option strings that apply to SSE only                       */
/*        m_BSEOnly         : option strings that apply to BSE only                       */
/*                                                                                        */
/*        m_RangeExcluded   : option strings for which a range may not be specified       */
/*        m_SetExcluded     : option strings for which a set may not be specified         */
/*                                                                                        */
/*    Read the explanations for each of the vectors in Options.h to get a better idea of  */
/*    what they are for and where the new option should go.                               */
/*                                                                                        */
/* 10. Add the new option to the following structures in constants.h (only required if    */
/*     the option is required to be available for printing in the logfiles):              */
/*                                                                                        */
/*        - enum class PROGRAM_OPTION                                                     */
/*        - const COMPASUnorderedMap<PROGRAM_OPTION, std::string> PROGRAM_OPTION_LABEL    */
/*        - const std::map<PROGRAM_OPTION, PROPERTY_DETAILS> PROGRAM_OPTION_DETAIL        */
/*                                                                                        */
/* 11. Add the new option to Options::OptionValue() - this enables selection of the       */
/*     option value for printing in the output (log) files.  Only required if the option  */
/*     is required to be available for printing in the logfiles.                          */
/*                                                                                        */
/******************************************************************************************/


#include "Options.h"
#include "changelog.h"

Options* Options::m_Instance = nullptr;

namespace po  = boost::program_options;
namespace cls = po::command_line_style;
namespace fs  = boost::filesystem;

// this is required to set default value for boost program options of type vector<std::string>
namespace std
{
  std::ostream& operator<<(std::ostream &os, const std::vector<std::string> &vec) {    
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

// for convenience

#define COMPLAIN(complainStr)           { std::stringstream _ss; _ss << complainStr; throw _ss.str(); }
#define COMPLAIN_IF(cond, complainStr)  { if (cond) COMPLAIN(complainStr) }

#define WARNUSER(warnStr)               { std::cerr << warnStr << std::endl; }
#define WARNUSER_IF(cond, warnStr)      { if (cond) WARNUSER(warnStr) }


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

    m_Populated = false;        

    // set all options to their default values

    // flags

    m_AllowNonStrippedECSN                                          = true;
    m_AllowRLOFAtBirth                                              = true;
    m_AllowTouchingAtBirth                                          = false;

    m_DebugToFile                                                   = false;
    m_ErrorsToFile                                                  = false;

    m_EnableWarnings                                                = false;

	m_BeBinaries                                                    = false;
    m_HMXRBinaries                                                  = false;

    m_EvolvePulsars                                                 = false;
	m_EvolveUnboundSystems                                          = false;

    m_DetailedOutput                                                = false;
    m_PopulationDataPrinting                                        = false;
    m_PrintBoolAsString                                             = false;
    m_Quiet                                                         = false;
    m_RlofPrinting                                                  = true;

    m_ShortHelp                                                     = true;

    m_StoreInputFiles                                               = true;


    m_SwitchLog                                                     = false;


    // annotations
    m_Notes.clear();
    m_NotesHdrs.clear();


    // Evolution mode: SSE or BSE
    m_EvolutionMode.type                                            = EVOLUTION_MODE::BSE;
    m_EvolutionMode.typeString                                      = EVOLUTION_MODE_LABEL.at(m_EvolutionMode.type);

    // Population synthesis variables
    m_ObjectsToEvolve                                               = 10;

    m_FixedRandomSeed                                               = false;
    m_RandomSeed                                                    = 0;

    // Specify how long to evolve for
    m_MaxEvolutionTime                                              = 13700.0;
    m_MaxNumberOfTimestepIterations                                 = 99999;
    m_TimestepMultiplier                                            = 1.0;

    // Initial mass options
    m_InitialMass                                                   = 5.0;
    m_InitialMass1                                                  = 5.0;
    m_InitialMass2                                                  = 5.0;

    m_InitialMassFunction.type                                      = INITIAL_MASS_FUNCTION::KROUPA;
    m_InitialMassFunction.typeString                                = INITIAL_MASS_FUNCTION_LABEL.at(m_InitialMassFunction.type);
    m_InitialMassFunctionMin                                        = 5.0;
    m_InitialMassFunctionMax                                        = 150.0;
    m_InitialMassFunctionPower                                      = 0.0;


    // Initial mass ratio
    m_MassRatio                                                     = 1.0;
    m_MassRatioDistribution.type                                    = MASS_RATIO_DISTRIBUTION::FLAT;
    m_MassRatioDistribution.typeString                              = MASS_RATIO_DISTRIBUTION_LABEL.at(m_MassRatioDistribution.type);
    m_MassRatioDistributionMin                                      = 0.01;
    m_MassRatioDistributionMax                                      = 1.0;

    m_MinimumMassSecondary                                          = 0.1;


    // Initial orbit options
    m_SemiMajorAxis                                                 = 0.1;

    m_SemiMajorAxisDistribution.type                                = SEMI_MAJOR_AXIS_DISTRIBUTION::FLATINLOG;
    m_SemiMajorAxisDistribution.typeString                          = SEMI_MAJOR_AXIS_DISTRIBUTION_LABEL.at(m_SemiMajorAxisDistribution.type);
    m_SemiMajorAxisDistributionMin                                  = 0.01;
    m_SemiMajorAxisDistributionMax                                  = 1000.0;
    m_SemiMajorAxisDistributionPower                                = -1.0;

    // Initial orbital period
    m_OrbitalPeriod                                                 = 0.1;                                                  // Only used if user specified and semi-major axis not specified

    // There is a single distribution available for orbital period (and actually,
    // eventually, it will be the case for all initial attributes that we have a 
    // single "convenience" distribution available in the C++ code - we expect 
    // users will sample outside the C++ code (with Stroopwafel etc.) so that we
    // don't have to code and maintain everybody's favourite distribution inside
    // the C++ code).
    //
    // The orbital period distribution will only used if it is specified by the
    // user AND semi-major axis (--semi-major-axis), orbital period (--orbital-period),
    // and semi-major axis distribution (--semi-major-axis-distribution) are NOT 
    // specified by the user.
    //
    // This --orbital-period-distribution option exists even though there is no real
    // choice (there is a single distribution available) so that users can specify
    // the distribution, and so it will be used according to the rule stated above.

    m_OrbitalPeriodDistribution.type                                = ORBITAL_PERIOD_DISTRIBUTION::FLATINLOG;
    m_OrbitalPeriodDistribution.typeString                          = ORBITAL_PERIOD_DISTRIBUTION_LABEL.at(m_OrbitalPeriodDistribution.type);
    m_OrbitalPeriodDistributionMin                                  = 1.1;
    m_OrbitalPeriodDistributionMax                                  = 1000.0;

    // Eccentricity
    m_Eccentricity                                                  = 0.0;
    m_EccentricityDistribution.type                                 = ECCENTRICITY_DISTRIBUTION::ZERO; 
    m_EccentricityDistribution.typeString                           = ECCENTRICITY_DISTRIBUTION_LABEL.at(m_EccentricityDistribution.type);
    m_EccentricityDistributionMin                                   = 0.0;
    m_EccentricityDistributionMax                                   = 1.0;

    // Kick options
    m_KickMagnitudeDistribution.type                                = KICK_MAGNITUDE_DISTRIBUTION::MAXWELLIAN;
    m_KickMagnitudeDistribution.typeString                          = KICK_MAGNITUDE_DISTRIBUTION_LABEL.at(m_KickMagnitudeDistribution.type);
    m_KickMagnitudeDistributionSigmaCCSN_NS                         = 265;
    m_KickMagnitudeDistributionSigmaCCSN_BH                         = 265;
    m_KickMagnitudeDistributionMaximum                              = -1.0; 
    m_KickMagnitudeDistributionSigmaForECSN                         = 30.0;
    m_KickMagnitudeDistributionSigmaForUSSN   	                    = 30.0;
	m_KickScalingFactor						                        = 1.0;

    // Kick direction option
    m_KickDirectionDistribution.type                                = KICK_DIRECTION_DISTRIBUTION::ISOTROPIC;
    m_KickDirectionDistribution.typeString                          = KICK_DIRECTION_DISTRIBUTION_LABEL.at(m_KickDirectionDistribution.type);
    m_KickDirectionPower                                            = 0.0;

    // Kick magnitude
    m_KickMagnitude                                                 = 0.0;
    m_KickMagnitude1                                                = 0.0;
    m_KickMagnitude2                                                = 0.0;                               

    m_MullerMandelKickBH                                            = MULLERMANDEL_KICKBH;
    m_MullerMandelKickNS                                            = MULLERMANDEL_KICKNS;
    m_MullerMandelSigmaKick                                         = MULLERMANDEL_SIGMAKICK;

    // Kick magnitude random number (used to draw kick magnitude if necessary)
    m_KickMagnitudeRandom                                           = 0.0;                                                  // actual value set later
    m_KickMagnitudeRandom1                                          = 0.0;                                                  // actual value set later
    m_KickMagnitudeRandom2                                          = 0.0;                                                  // actual value set later

    // Mean anomaly
    m_KickMeanAnomaly1                                              = 0.0;                                                  // actual value set later
    m_KickMeanAnomaly2                                              = 0.0;                                                  // actual value set later

    // Phi
    m_KickPhi1                                                      = 0.0;                                                  // actual value set later
    m_KickPhi2                                                      = 0.0;                                                  // actual value set later

    // Theta
    m_KickTheta1                                                    = 0.0;                                                  // actual value set later 
    m_KickTheta2                                                    = 0.0;                                                  // actual value set later

    // Black hole kicks
    m_BlackHoleKicks.type                                           = BLACK_HOLE_KICKS::FALLBACK;
    m_BlackHoleKicks.typeString                                     = BLACK_HOLE_KICKS_LABEL.at(m_BlackHoleKicks.type);


    // Chemically Homogeneous Evolution
    m_CheMode.type                                                  = CHE_MODE::PESSIMISTIC;
    m_CheMode.typeString                                            = CHE_MODE_LABEL.at(m_CheMode.type);


    // Supernova remnant mass prescription options
    m_RemnantMassPrescription.type                                  = REMNANT_MASS_PRESCRIPTION::FRYER2012;
    m_RemnantMassPrescription.typeString                            = REMNANT_MASS_PRESCRIPTION_LABEL.at(m_RemnantMassPrescription.type);

    m_FryerSupernovaEngine.type                                     = SN_ENGINE::DELAYED;
    m_FryerSupernovaEngine.typeString                               = SN_ENGINE_LABEL.at(m_FryerSupernovaEngine.type);

    m_Fryer22fmix                                                   = 0.5; //default is similar to DELAYED engine in Fryer 2012
    m_Fryer22Mcrit                                                  = 5.75; //

    m_NeutrinoMassLossAssumptionBH.type                             = NEUTRINO_MASS_LOSS_PRESCRIPTION::FIXED_MASS;
    m_NeutrinoMassLossAssumptionBH.typeString                       = NEUTRINO_MASS_LOSS_PRESCRIPTION_LABEL.at(m_NeutrinoMassLossAssumptionBH.type);
    m_NeutrinoMassLossValueBH                                       = 0.1;


    // Fixed uk options
    m_UseFixedUK                                                    = false;
    m_FixedUK                                                       = -1.0;


    // Pair instability and pulsational pair instability mass loss
    m_UsePairInstabilitySupernovae                                  = true;
    m_PairInstabilityLowerLimit                                     = 60.0;                                                 // Belczynski+ 2016 is 65 Msol
    m_PairInstabilityUpperLimit                                     = 135.0;                                                // Belczynski+ 2016 is 135 Msol

    m_UsePulsationalPairInstability                                 = true;
    m_PulsationalPairInstabilityLowerLimit                          = 35.0;                                                 // Belczynski+ 2016 is 45 Msol
    m_PulsationalPairInstabilityUpperLimit                          = 60.0;                                                 // Belczynski+ 2016 is 65 Msol

    m_PulsationalPairInstabilityPrescription.type                   = PPI_PRESCRIPTION::MARCHANT;
    m_PulsationalPairInstabilityPrescription.typeString             = PPI_PRESCRIPTION_LABEL.at(m_PulsationalPairInstabilityPrescription.type);

	m_MaximumNeutronStarMass                                        = 2.5;                                                  // StarTrack is 3.0
    
    m_mCBUR1                                                        = MCBUR1HURLEY;                                         // MHurley value, Fryer+ and Belczynski+ use 1.83


    // Output path
    m_OutputPathString                                              = ".";
    m_DefaultOutputPath                                             = boost::filesystem::current_path();
    m_OutputPath                                                    = m_DefaultOutputPath;
    m_OutputContainerName                                           = DEFAULT_OUTPUT_CONTAINER_NAME;
    

    // Mass loss options
    m_UseMassLoss                                                   = true;
    m_CheckPhotonTiringLimit                                        = false;

    m_MassLossPrescription.type                                     = MASS_LOSS_PRESCRIPTION::VINK;
    m_MassLossPrescription.typeString                               = MASS_LOSS_PRESCRIPTION_LABEL.at(m_MassLossPrescription.type);

    m_LuminousBlueVariablePrescription.type                         = LBV_PRESCRIPTION::HURLEY_ADD;
    m_LuminousBlueVariablePrescription.typeString                   = LBV_PRESCRIPTION_LABEL.at(m_LuminousBlueVariablePrescription.type);

    // Wind mass loss multiplicitive constants
    m_CoolWindMassLossMultiplier                                    = 1.0;
    m_LuminousBlueVariableFactor                                    = 1.5;
    m_OverallWindMassLossMultiplier                                 = 1.0;
    m_WolfRayetFactor                                               = 1.0;


    // Mass transfer options
    m_UseMassTransfer                                               = true;
	m_CirculariseBinaryDuringMassTransfer         	                = true;
	m_AngularMomentumConservationDuringCircularisation              = false;
    m_RetainCoreMassDuringCaseAMassTransfer                         = false;
    m_ConvectiveEnvelopeTemperatureThreshold                        = CONVECTIVE_BOUNDARY_TEMPERATURE_BELCZYNSKI;

    // Case BB/BC mass transfer stability prescription
    m_CaseBBStabilityPrescription.type                              = CASE_BB_STABILITY_PRESCRIPTION::ALWAYS_STABLE;
    m_CaseBBStabilityPrescription.typeString                        = CASE_BB_STABILITY_PRESCRIPTION_LABEL.at(m_CaseBBStabilityPrescription.type);

    // Options for mass transfer accretion efficiency
    m_MassTransferAccretionEfficiencyPrescription.type              = MT_ACCRETION_EFFICIENCY_PRESCRIPTION::THERMALLY_LIMITED;
    m_MassTransferAccretionEfficiencyPrescription.typeString        = MT_ACCRETION_EFFICIENCY_PRESCRIPTION_LABEL.at(m_MassTransferAccretionEfficiencyPrescription.type);

    m_MassTransferFractionAccreted                                  = 0.5;
    m_MassTransferCParameter                                        = 10.0;
    m_EddingtonAccretionFactor                                      = 1;                                                    // >1 is super-eddington, 0 is no accretion

    // Mass transfer thermally limited options
	m_MassTransferThermallyLimitedVariation.type                    = MT_THERMALLY_LIMITED_VARIATION::C_FACTOR;
	m_MassTransferThermallyLimitedVariation.typeString              = MT_THERMALLY_LIMITED_VARIATION_LABEL.at(m_MassTransferThermallyLimitedVariation.type);

    // Mass transfer angular momentum loss prescription options
    m_MassTransferJloss                                             = 1.0;
    m_MassTransferJlossMacLeodLinearFraction                        = 0.5;
    m_MassTransferAngularMomentumLossPrescription.type              = MT_ANGULAR_MOMENTUM_LOSS_PRESCRIPTION::ISOTROPIC_RE_EMISSION;
    m_MassTransferAngularMomentumLossPrescription.typeString        = MT_ANGULAR_MOMENTUM_LOSS_PRESCRIPTION_LABEL.at(m_MassTransferAngularMomentumLossPrescription.type);

    // Mass transfer rejuvenation prescriptions
    m_MassTransferRejuvenationPrescription.type                     = MT_REJUVENATION_PRESCRIPTION::STARTRACK;
    m_MassTransferRejuvenationPrescription.typeString               = MT_REJUVENATION_PRESCRIPTION_LABEL.at(m_MassTransferRejuvenationPrescription.type);

    // Mass transfer critical mass ratios - defined here as (accretor mass / donor mass)
    // A value of 0.0 means the mass ratio will never be unstable - this does not guaruntee stability of the MT, just that instability is not based on the mass ratio

    m_QCritPrescription.type                                        = QCRIT_PRESCRIPTION::NONE;                             // Assume no critical mass ratio prescription
    m_QCritPrescription.typeString                                  = QCRIT_PRESCRIPTION_LABEL.at(m_QCritPrescription.type);
    m_MassTransferCriticalMassRatioMSLowMassNonDegenerateAccretor   = 1.44;                                                 // Claeys+ 2014 = 1.44
    m_MassTransferCriticalMassRatioMSLowMassDegenerateAccretor      = 1.0;                                                  // Claeys+ 2014 = 1.0

    m_MassTransferCriticalMassRatioMSHighMassNonDegenerateAccretor  = 0.625;                                                // Claeys+ 2014 = 0.625
    m_MassTransferCriticalMassRatioMSHighMassDegenerateAccretor     = 0.0;                                                  // Claeys+ 2014 = unspecified

    m_MassTransferCriticalMassRatioHGNonDegenerateAccretor          = 0.25;                                                 // Claeys+ 2014 = 0.25
    m_MassTransferCriticalMassRatioHGDegenerateAccretor             = 0.21;                                                 // Claeys+ 2014 = 0.21

    m_MassTransferCriticalMassRatioGiantNonDegenerateAccretor       = -1.0;                                                 // Value not used! Claeys+ 2014 uses an equation in mass-radius exponent and core mass (equivalent to Hurley zeta adiabatic definition), so if -1, this value is overwritten later
    m_MassTransferCriticalMassRatioGiantDegenerateAccretor          = 0.87;                                                 // Claeys+ 2014 = 0.87

    m_MassTransferCriticalMassRatioHeliumMSNonDegenerateAccretor    = 0.0;                                                  // Claeys+ 2014 = unspecified 
    m_MassTransferCriticalMassRatioHeliumMSDegenerateAccretor       = 0.0;                                                  // Claeys+ 2014 = unspecified

    m_MassTransferCriticalMassRatioHeliumHGNonDegenerateAccretor    = 0.25;                                                 // Claeys+ 2014 = 0.25
    m_MassTransferCriticalMassRatioHeliumHGDegenerateAccretor       = 0.21;                                                 // Claeys+ 2014 = 0.21

    m_MassTransferCriticalMassRatioHeliumGiantNonDegenerateAccretor = 1.28;                                                 // Claeys+ 2014 = 1.28
    m_MassTransferCriticalMassRatioHeliumGiantDegenerateAccretor    = 0.87;                                                 // Claeys+ 2014 = 0.87

	m_MassTransferCriticalMassRatioWhiteDwarfNonDegenerateAccretor  = 0.0;                                                  // Claeys+ 2014 = unspecified
    m_MassTransferCriticalMassRatioWhiteDwarfDegenerateAccretor     = 1.6;                                                  // Claeys+ 2014 = 1.6

    // Common Envelope options
    m_CommonEnvelopeAlpha                                           = 1.0;
    m_CommonEnvelopeLambda                                          = 0.1;
	m_CommonEnvelopeSlopeKruckow                                    = -5.0 / 6.0;
	m_CommonEnvelopeAlphaThermal                                    = 1.0;
    m_CommonEnvelopeLambdaMultiplier                                = 1.0;
    m_CommonEnvelopeLambdaNanjingEnhanced                           = false;
    m_CommonEnvelopeLambdaNanjingInterpolateInMass                  = false;
    m_CommonEnvelopeLambdaNanjingInterpolateInMetallicity           = false;
    m_CommonEnvelopeLambdaNanjingUseRejuvenatedMass                 = false;
    m_AllowRadiativeEnvelopeStarToSurviveCommonEnvelope             = false;
    m_AllowMainSequenceStarToSurviveCommonEnvelope                  = true;
    m_AllowImmediateRLOFpostCEToSurviveCommonEnvelope               = false;

    // Prescription for envelope state (radiative or convective)
    m_EnvelopeStatePrescription.type                                = ENVELOPE_STATE_PRESCRIPTION::LEGACY;
    m_EnvelopeStatePrescription.typeString                          = ENVELOPE_STATE_PRESCRIPTION_LABEL.at(m_EnvelopeStatePrescription.type);

    // Accretion during common envelope
    m_CommonEnvelopeMassAccretionPrescription.type                  = CE_ACCRETION_PRESCRIPTION::ZERO;
    m_CommonEnvelopeMassAccretionPrescription.typeString            = CE_ACCRETION_PRESCRIPTION_LABEL.at(m_CommonEnvelopeMassAccretionPrescription.type);
    
    m_CommonEnvelopeMassAccretionMin                                = 0.04;
    m_CommonEnvelopeMassAccretionMax                                = 0.1;
    m_CommonEnvelopeMassAccretionConstant                           = 0.0;

	// Common envelope lambda prescription
	m_CommonEnvelopeLambdaPrescription.type                         = CE_LAMBDA_PRESCRIPTION::NANJING;
	m_CommonEnvelopeLambdaPrescription.typeString                   = CE_LAMBDA_PRESCRIPTION_LABEL.at(m_CommonEnvelopeLambdaPrescription.type);

	// Common envelope Nandez and Ivanova energy formalism
	m_RevisedEnergyFormalismNandezIvanova	                        = false;
	m_MaximumMassDonorNandezIvanova                                 = 2.0;
	m_CommonEnvelopeRecombinationEnergyDensity                      = 1.5E13;


	// Zetas
	m_StellarZetaPrescription.type                                  = ZETA_PRESCRIPTION::SOBERMAN;
	m_StellarZetaPrescription.typeString                            = ZETA_PRESCRIPTION_LABEL.at(m_StellarZetaPrescription.type);

	m_ZetaAdiabaticArbitrary                                        = 10000.0;                                              // large value favours stable MT
    m_ZetaMainSequence 	                                            = 2.0;
	m_ZetaRadiativeEnvelopeGiant	                                = 6.5;


    // Metallicity options
    m_Metallicity                                                   = ZSOL_ASPLUND;
    m_MetallicityDistribution.type                                  = METALLICITY_DISTRIBUTION::ZSOLAR;
    m_MetallicityDistribution.typeString                            = METALLICITY_DISTRIBUTION_LABEL.at(m_MetallicityDistribution.type);
    m_MetallicityDistributionMin                                    = MINIMUM_METALLICITY;
    m_MetallicityDistributionMax                                    = MAXIMUM_METALLICITY;


    // Neutron star equation of state
    m_NeutronStarEquationOfState.type                               = NS_EOS::SSE;
    m_NeutronStarEquationOfState.typeString                         = NS_EOSLabel.at(m_NeutronStarEquationOfState.type);


    // Pulsar birth magnetic field distribution
    m_PulsarBirthMagneticFieldDistribution.type                     = PULSAR_BIRTH_MAGNETIC_FIELD_DISTRIBUTION::ZERO;
    m_PulsarBirthMagneticFieldDistribution.typeString               = PULSAR_BIRTH_MAGNETIC_FIELD_DISTRIBUTION_LABEL.at(m_PulsarBirthMagneticFieldDistribution.type);
    m_PulsarBirthMagneticFieldDistributionMin                       = 11.0;
    m_PulsarBirthMagneticFieldDistributionMax                       = 13.0;


    // Pulsar birth spin period distribution string
    m_PulsarBirthSpinPeriodDistribution.type                        = PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION::ZERO;
    m_PulsarBirthSpinPeriodDistribution.typeString                  = PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION_LABEL.at(m_PulsarBirthSpinPeriodDistribution.type);
    m_PulsarBirthSpinPeriodDistributionMin                          = 10.0;
    m_PulsarBirthSpinPeriodDistributionMax                          = 100.0;

    m_PulsarMagneticFieldDecayTimescale                             = 1000.0;
    m_PulsarMagneticFieldDecayMassscale                             = 0.025;
    m_PulsarLog10MinimumMagneticField                               = 8.0;


    // Rotational velocity distribution options
    m_RotationalVelocityDistribution.type                           = ROTATIONAL_VELOCITY_DISTRIBUTION::ZERO;
    m_RotationalVelocityDistribution.typeString                     = ROTATIONAL_VELOCITY_DISTRIBUTION_LABEL.at(m_RotationalVelocityDistribution.type);

    m_RotationalFrequency                                           = 0.0;
    m_RotationalFrequency1                                          = 0.0;
    m_RotationalFrequency2                                          = 0.0;


	// grids

	m_GridFilename                                                  = "";
    m_GridStartLine                                                 = 0;
    m_GridLinesToProcess                                            = std::numeric_limits<std::streamsize>::max();                  // effectively no limit - process to EOF

    // debug and logging options

    m_DebugLevel                                                    = 0;
    m_DebugClasses.clear();

    m_LogLevel                                                      = 0;
    m_LogClasses.clear();

    // Logfiles    
    m_LogfileDefinitionsFilename                                    = "";
    m_LogfileNamePrefix                                             = "";
    m_LogfileType.type                                              = LOGFILETYPE::HDF5;
    m_LogfileType.typeString                                        = LOGFILETYPELabel.at(m_LogfileType.type);

    m_LogfileBeBinaries                                             = std::get<0>(LOGFILE_DESCRIPTOR.at(LOGFILE::BSE_BE_BINARIES));
    m_LogfileCommonEnvelopes                                        = std::get<0>(LOGFILE_DESCRIPTOR.at(LOGFILE::BSE_COMMON_ENVELOPES));
    m_LogfileDetailedOutput                                         = std::get<0>(LOGFILE_DESCRIPTOR.at(LOGFILE::BSE_DETAILED_OUTPUT));  // assume BSE - get real answer when we know mode
    m_LogfileDoubleCompactObjects                                   = std::get<0>(LOGFILE_DESCRIPTOR.at(LOGFILE::BSE_DOUBLE_COMPACT_OBJECTS));
    m_LogfilePulsarEvolution                                        = std::get<0>(LOGFILE_DESCRIPTOR.at(LOGFILE::BSE_PULSAR_EVOLUTION)); // only BSE for now
    m_LogfileRLOFParameters                                         = std::get<0>(LOGFILE_DESCRIPTOR.at(LOGFILE::BSE_RLOF_PARAMETERS));
    m_LogfileSupernovae                                             = std::get<0>(LOGFILE_DESCRIPTOR.at(LOGFILE::BSE_SUPERNOVAE));       // assume BSE - get real answer when we know mode
    m_LogfileSwitchLog                                              = std::get<0>(LOGFILE_DESCRIPTOR.at(LOGFILE::BSE_SWITCH_LOG));       // assume BSE - get real answer when we know mode
    m_LogfileSystemParameters                                       = std::get<0>(LOGFILE_DESCRIPTOR.at(LOGFILE::BSE_SYSTEM_PARAMETERS));

    m_AddOptionsToSysParms.type                                     = ADD_OPTIONS_TO_SYSPARMS::GRID;
    m_AddOptionsToSysParms.typeString                               = ADD_OPTIONS_TO_SYSPARMS_LABEL.at(m_AddOptionsToSysParms.type);
    
    m_HDF5BufferSize                                                = HDF5_DEFAULT_IO_BUFFER_SIZE;
    m_HDF5ChunkSize                                                 = HDF5_DEFAULT_CHUNK_SIZE;

    po::variables_map vm;
    m_VM = vm;
}


/*
 * This function populates the Boost options_description object with all defined
 * options.  The options strings are associated with the variable to be populated
 * for each option, and the default values are defined.
 * 
 * JR: todo: One day we should construct the list of options shown in the help
 *           string from the maps in constants.h rather than have them here as 
 *           literal strings - too much opportunity for then to get out of sync 
 *           doing it this way (and more work to update the strings here every 
 *           time something changes)
 * 
 * Note that both parameters are modified here.
 * 
 * 
 * bool AddOptions(OptionValues *p_Options, po::options_description *p_OptionsDescription)
 * 
 * @param   [IN/OUT]    p_Options                   Object containing option values
 * @param   [IN/OUT]    p_OptionsDescription        options_description onject
 * @return                                          Boolean result:
 *                                                      true  if options added ok
 *                                                      false if an error occurred
 */
bool Options::AddOptions(OptionValues *p_Options, po::options_description *p_OptionsDescription) {

    bool ok = true;                             // status - unless a problem occurs

    // create default strings for std::vector<std::string> types (too hard to do inline)

    std::ostringstream ss;

    // debug classes
    std::string defaultDebugClasses;
    ss << "";
    for (auto debugClass = p_Options->m_DebugClasses.begin(); debugClass != p_Options->m_DebugClasses.end(); ++debugClass) ss << *debugClass << ",";
    defaultDebugClasses = ss.str();
    if (defaultDebugClasses.length() > 0) defaultDebugClasses.erase(defaultDebugClasses.length() - 1);

    // log classes
    std::string defaultLogClasses;
    ss << "";
    for (auto logClass = p_Options->m_LogClasses.begin(); logClass != p_Options->m_LogClasses.end(); ++logClass) ss << *logClass << ",";
    defaultLogClasses = ss.str();
    if (defaultLogClasses.length() > 0) defaultLogClasses.erase(defaultLogClasses.length() - 1);

    // annotations
    std::string defaultNotes;
    ss << "";
    for (auto note = p_Options->m_Notes.begin(); note != p_Options->m_Notes.end(); ++note) ss << *note << ",";
    defaultNotes = ss.str();
    if (defaultNotes.length() > 0) defaultNotes.erase(defaultNotes.length() - 1);

    // annotation headers
    std::string defaultNotesHdrs;
    ss << "";
    for (auto noteHdr = p_Options->m_NotesHdrs.begin(); noteHdr != p_Options->m_NotesHdrs.end(); ++noteHdr) ss << *noteHdr << ",";
    defaultNotesHdrs = ss.str();
    if (defaultNotesHdrs.length() > 0) defaultNotesHdrs.erase(defaultNotesHdrs.length() - 1);


    // add options

    try {

        p_OptionsDescription->add_options()     // begin the list of options to be added - boost syntactic sugar

        // there is no good way of formatting these - the boost syntax doesn't help that much
        // there is just a boatload of options, so this function is just going to be long...
        // the options are (kind-of) ordered by data type (i.e. bool, int, double etc.) and
        // mostly alphabetic in the data types (but might be grouped by functionality in the
        // data types if it makes more sense)


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

        (
            "allow-non-stripped-ECSN",
            po::value<bool>(&p_Options->m_AllowNonStrippedECSN)->default_value(p_Options->m_AllowNonStrippedECSN)->implicit_value(true),                                                                  
            ("Allow ECSN to occur in unstripped progenitors (default = " + std::string(p_Options->m_AllowNonStrippedECSN? "TRUE" : "FALSE") + ")").c_str()
        )
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
            "angular-momentum-conservation-during-circularisation",            
            po::value<bool>(&p_Options->m_AngularMomentumConservationDuringCircularisation)->default_value(p_Options->m_AngularMomentumConservationDuringCircularisation)->implicit_value(true),  
            ("Conserve angular momentum when binary is circularised when entering a Mass Transfer episode (default = " + std::string(p_Options->m_AngularMomentumConservationDuringCircularisation ? "TRUE" : "FALSE") + ")").c_str()
        )

        /*
        (
            "BE-binaries",                                                  
            po::value<bool>(&p_Options->m_BeBinaries)->default_value(p_Options->m_BeBinaries)->implicit_value(true),                                                                              
            ("Enable Be Binaries study (default = " + std::string(p_Options->m_BeBinaries ? "TRUE" : "FALSE") + ")").c_str()
        )
        */
        (
            "check-photon-tiring-limit",
            po::value<bool>(&p_Options->m_CheckPhotonTiringLimit)->default_value(p_Options->m_CheckPhotonTiringLimit)->implicit_value(true),                            
            ("Check the photon tiring limit hasn't been exceeded by wind mass loss (default = " + std::string(p_Options->m_CheckPhotonTiringLimit ? "TRUE" : "FALSE") + ")").c_str()
        )
        (
            "circularise-binary-during-mass-transfer",                         
            po::value<bool>(&p_Options->m_CirculariseBinaryDuringMassTransfer)->default_value(p_Options->m_CirculariseBinaryDuringMassTransfer)->implicit_value(true),                            
            ("Circularise binary when it enters a Mass Transfer episode (default = " + std::string(p_Options->m_CirculariseBinaryDuringMassTransfer ? "TRUE" : "FALSE") + ")").c_str()
        )
        (
            "common-envelope-allow-immediate-RLOF-post-CE-survive",
            po::value<bool>(&p_Options->m_AllowImmediateRLOFpostCEToSurviveCommonEnvelope)->default_value(p_Options->m_AllowImmediateRLOFpostCEToSurviveCommonEnvelope)->implicit_value(true),
            ("Allow immediate post CE RLOF to survive common envelope evolution (default = " + std::string(p_Options->m_AllowImmediateRLOFpostCEToSurviveCommonEnvelope ? "TRUE" : "FALSE") + ")").c_str()
        )
        (
            "common-envelope-allow-main-sequence-survive",                 
            po::value<bool>(&p_Options->m_AllowMainSequenceStarToSurviveCommonEnvelope)->default_value(p_Options->m_AllowMainSequenceStarToSurviveCommonEnvelope)->implicit_value(true),          
            ("Allow main sequence stars to survive common envelope evolution (default = " + std::string(p_Options->m_AllowMainSequenceStarToSurviveCommonEnvelope ? "TRUE" : "FALSE") + ")").c_str()
        )
        (
            "common-envelope-lambda-nanjing-enhanced",
            po::value<bool>(&p_Options->m_CommonEnvelopeLambdaNanjingEnhanced)->default_value(p_Options->m_CommonEnvelopeLambdaNanjingEnhanced)->implicit_value(true),
            ("Use Nanjing lambda's with enhanced extrapolation in stellar radius (default = " + std::string(p_Options->m_CommonEnvelopeLambdaNanjingEnhanced ? "TRUE" : "FALSE") + ")").c_str()
        )
        (
            "common-envelope-lambda-nanjing-interpolate-in-mass",
            po::value<bool>(&p_Options->m_CommonEnvelopeLambdaNanjingInterpolateInMass)->default_value(p_Options->m_CommonEnvelopeLambdaNanjingInterpolateInMass)->implicit_value(true),
            ("Use Nanjing lambda's with mass interpolation (only used when using enhanced Nanjing lambda's) (default = " + std::string(p_Options->m_CommonEnvelopeLambdaNanjingInterpolateInMass ? "TRUE" : "FALSE") + ")").c_str()
        )
        (
            "common-envelope-lambda-nanjing-interpolate-in-metallicity",
            po::value<bool>(&p_Options->m_CommonEnvelopeLambdaNanjingInterpolateInMetallicity)->default_value(p_Options->m_CommonEnvelopeLambdaNanjingInterpolateInMetallicity)->implicit_value(true),
            ("Use Nanjing lambda's with metallicity interpolation (only used when using enhanced Nanjing lambda's) (default = " + std::string(p_Options->m_CommonEnvelopeLambdaNanjingInterpolateInMetallicity ? "TRUE" : "FALSE") + ")").c_str()
        )
        (
            "common-envelope-lambda-nanjing-use-rejuvenated-mass",
            po::value<bool>(&p_Options->m_CommonEnvelopeLambdaNanjingUseRejuvenatedMass)->default_value(p_Options->m_CommonEnvelopeLambdaNanjingUseRejuvenatedMass)->implicit_value(true),
            ("Use rejuvenated mass to calculate Nanjing lambda's (default = " + std::string(p_Options->m_CommonEnvelopeLambdaNanjingUseRejuvenatedMass ? "TRUE" : "FALSE") + ")").c_str()
        )
        (
            "common-envelope-allow-radiative-envelope-survive",
            po::value<bool>(&p_Options->m_AllowRadiativeEnvelopeStarToSurviveCommonEnvelope)->default_value(p_Options->m_AllowRadiativeEnvelopeStarToSurviveCommonEnvelope)->implicit_value(true),
            ("Allow radiative envelope stars to survive common envelope evolution (default = " + std::string(p_Options->m_AllowRadiativeEnvelopeStarToSurviveCommonEnvelope ? "TRUE" : "FALSE") + ")").c_str()
        )
        (
            "cool-wind-mass-loss-multiplier",                           
            po::value<double>(&p_Options->m_CoolWindMassLossMultiplier)->default_value(p_Options->m_CoolWindMassLossMultiplier),                                                                  
            ("Multiplicative constant for wind mass loss of cool stars (default = " + std::to_string(p_Options->m_CoolWindMassLossMultiplier)+ ")").c_str()
        )
        (
            "debug-to-file",                                               
            po::value<bool>(&p_Options->m_DebugToFile)->default_value(p_Options->m_DebugToFile)->implicit_value(true),                                                                            
            ("Write debug statements to file (default = " + std::string(p_Options->m_DebugToFile ? "TRUE" : "FALSE") + ")").c_str()
        )
        (
            "detailed-output",                                              
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
            "mass-transfer",                                                
            po::value<bool>(&p_Options->m_UseMassTransfer)->default_value(p_Options->m_UseMassTransfer)->implicit_value(true),                                                                    
            ("Enable mass transfer (default = " + std::string(p_Options->m_UseMassTransfer ? "TRUE" : "FALSE") + ")").c_str()
        )
        (
            "pair-instability-supernovae",                                 
            po::value<bool>(&p_Options->m_UsePairInstabilitySupernovae)->default_value(p_Options->m_UsePairInstabilitySupernovae)->implicit_value(true),                                          
            ("Enable pair instability supernovae (PISN) (default = " + std::string(p_Options->m_UsePairInstabilitySupernovae ? "TRUE" : "FALSE") + ")").c_str()
        )
        (
            "population-data-printing",                                      
            po::value<bool>(&p_Options->m_PopulationDataPrinting)->default_value(p_Options->m_PopulationDataPrinting)->implicit_value(true),                                                      
            ("Print details of population (default = " + std::string(p_Options->m_PopulationDataPrinting ? "TRUE" : "FALSE") + ")").c_str()
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
            "quiet",                                                       
            po::value<bool>(&p_Options->m_Quiet)->default_value(p_Options->m_Quiet)->implicit_value(true),                                                                                        
            ("Suppress printing (default = " + std::string(p_Options->m_Quiet ? "TRUE" : "FALSE") + ")").c_str()
        )
        (
            "retain-core-mass-during-caseA-mass-transfer",
            po::value<bool>(&p_Options->m_RetainCoreMassDuringCaseAMassTransfer)->default_value(p_Options->m_RetainCoreMassDuringCaseAMassTransfer)->implicit_value(true),
            ("Retain approximate core mass of a case A donor as a minimum core at end of MS or HeMS (default = " + std::string(p_Options->m_RetainCoreMassDuringCaseAMassTransfer ? "TRUE" : "FALSE") + ")").c_str()
        )
        (
            "revised-energy-formalism-nandez-ivanova",                     
            po::value<bool>(&p_Options->m_RevisedEnergyFormalismNandezIvanova)->default_value(p_Options->m_RevisedEnergyFormalismNandezIvanova)->implicit_value(true),                            
            ("Enable revised energy formalism (default = " + std::string(p_Options->m_RevisedEnergyFormalismNandezIvanova ? "TRUE" : "FALSE") + ")").c_str()
        )
        (
            "rlof-printing",                                                
            po::value<bool>(&p_Options->m_RlofPrinting)->default_value(p_Options->m_RlofPrinting)->implicit_value(true),                                                                          
            ("Enable output parameters before/after RLOF (default = " + std::string(p_Options->m_RlofPrinting ? "TRUE" : "FALSE") + ")").c_str()
        )
        (
            "store-input-files",                                                
            po::value<bool>(&p_Options->m_StoreInputFiles)->default_value(p_Options->m_StoreInputFiles)->implicit_value(true),                                                                          
            ("Store input files in output container (default = " + std::string(p_Options->m_StoreInputFiles ? "TRUE" : "FALSE") + ")").c_str()
        )
        (
            "hmxr-binaries",
            po::value<bool>(&p_Options->m_HMXRBinaries)->default_value(p_Options->m_HMXRBinaries)->implicit_value(true),
            ("Store HMXRB candidates in BSE_RLOF output file (default = " + std::string(p_Options->m_HMXRBinaries ? "TRUE" : "FALSE") + ")").c_str()
        )
        (
            "switch-log",                                                
            po::value<bool>(&p_Options->m_SwitchLog)->default_value(p_Options->m_SwitchLog)->implicit_value(true),                                                                          
            ("Print switch log to file (default = " + std::string(p_Options->m_SwitchLog ? "TRUE" : "FALSE") + ")").c_str()
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


        // int / unsigned int

        (
            "debug-level",                                                 
            po::value<int>(&p_Options->m_DebugLevel)->default_value(p_Options->m_DebugLevel),                                                                                                     
            ("Determines which print statements are displayed for debugging (default = " + std::to_string(p_Options->m_DebugLevel) + ")").c_str()
        )
        (
            "grid-start-line",                                                 
            po::value<std::streamsize>(&p_Options->m_GridStartLine)->default_value(p_Options->m_GridStartLine),                                                                                                     
            ("Specifies which line of the grid file is processed first (0-based) (default = " + std::to_string(p_Options->m_GridStartLine) + ")").c_str()
        )
        (
            "grid-lines-to-process",                                                 
            po::value<std::streamsize>(&p_Options->m_GridLinesToProcess)->default_value(p_Options->m_GridLinesToProcess),                                                                                                     
            ("Specifies how many grid lines should be processed (from the start line - see grid-start-line) (default = " + (p_Options->m_GridLinesToProcess == std::numeric_limits<std::streamsize>::max() ? "Process to EOF" : std::to_string(p_Options->m_GridLinesToProcess)) + ")").c_str()
        )
        (
            "hdf5-chunk-size",                                                 
            po::value<int>(&p_Options->m_HDF5ChunkSize)->default_value(p_Options->m_HDF5ChunkSize),                                                                                                     
            ("HDF5 file dataset chunk size (number of dataset entries, default = " + std::to_string(p_Options->m_HDF5ChunkSize) + ")").c_str()
        )
        (
            "hdf5-buffer-size",                                                 
            po::value<int>(&p_Options->m_HDF5BufferSize)->default_value(p_Options->m_HDF5BufferSize),                                                                                                     
            ("HDF5 file dataset IO buffer size (number of chunks, default = " + std::to_string(p_Options->m_HDF5BufferSize) + ")").c_str()
        )
        (
            "log-level",                                                   
            po::value<int>(&p_Options->m_LogLevel)->default_value(p_Options->m_LogLevel),                                                                                                         
            ("Determines which print statements are included in the logfile (default = " + std::to_string(p_Options->m_LogLevel) + ")").c_str()
        )
        (
            "maximum-number-timestep-iterations",                          
            po::value<int>(&p_Options->m_MaxNumberOfTimestepIterations)->default_value(p_Options->m_MaxNumberOfTimestepIterations),                                                               
            ("Maximum number of timesteps to evolve binary before giving up (default = " + std::to_string(p_Options->m_MaxNumberOfTimestepIterations) + ")").c_str()
        )
        (
            "number-of-systems,n",                                        
            po::value<int>(&p_Options->m_ObjectsToEvolve)->default_value(p_Options->m_ObjectsToEvolve),                                                                                                       
            ("Specify the number of systems to simulate (SSE) (default = " + std::to_string(p_Options->m_ObjectsToEvolve) + ")").c_str()
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
            ("Value of mass accreted by NS/BH, in Msol, during common envelope evolution, assuming all NS/BH accrete same amount of mass (common-envelope-mass-accretion-prescription CONSTANT). Ignored otherwise (default = " + std::to_string(p_Options->m_CommonEnvelopeMassAccretionConstant) + ")").c_str()
        )
        (
            "common-envelope-mass-accretion-max",                          
            po::value<double>(&p_Options->m_CommonEnvelopeMassAccretionMax)->default_value(p_Options->m_CommonEnvelopeMassAccretionMax),                                                          
            ("Maximum amount of mass accreted by NS/BHs, in Msol, during common envelope evolution in Msol (default = " + std::to_string(p_Options->m_CommonEnvelopeMassAccretionMax) + ")").c_str()
        )
        (
            "common-envelope-mass-accretion-min",                          
            po::value<double>(&p_Options->m_CommonEnvelopeMassAccretionMin)->default_value(p_Options->m_CommonEnvelopeMassAccretionMin),                                                          
            ("Minimum amount of mass accreted by NS/BHs, in Msol, during common envelope evolution in Msol (default = " + std::to_string(p_Options->m_CommonEnvelopeMassAccretionMin) + ")").c_str()
        )
        (
            "common-envelope-recombination-energy-density",                
            po::value<double>(&p_Options->m_CommonEnvelopeRecombinationEnergyDensity)->default_value(p_Options->m_CommonEnvelopeRecombinationEnergyDensity),                                      
            ("Recombination energy density, in erg/g (default = " + std::to_string(p_Options->m_CommonEnvelopeRecombinationEnergyDensity) + ")").c_str()
        )
        (
            "common-envelope-slope-kruckow",                               
            po::value<double>(&p_Options->m_CommonEnvelopeSlopeKruckow)->default_value(p_Options->m_CommonEnvelopeSlopeKruckow),                                                                  
            ("Common Envelope slope for Kruckow lambda (default = " + std::to_string(p_Options->m_CommonEnvelopeSlopeKruckow) + ")").c_str()
        )
        (
            "convective-envelope-temperature-threshold",                               
            po::value<double>(&p_Options->m_ConvectiveEnvelopeTemperatureThreshold)->default_value(p_Options->m_ConvectiveEnvelopeTemperatureThreshold),                                                                  
            ("Temperature [K] threshold, below which the envelopes of giants are convective. Only used for --envelope-state-prescription = FIXED_TEMPERATURE, ignored otherwise. (default = " + std::to_string(p_Options->m_ConvectiveEnvelopeTemperatureThreshold) + ")").c_str()
        )

        (
            "critical-mass-ratio-giant-degenerate-accretor",
            po::value<double>(&p_Options->m_MassTransferCriticalMassRatioGiantDegenerateAccretor)->default_value(p_Options->m_MassTransferCriticalMassRatioGiantDegenerateAccretor),
            ("Critical mass ratio for MT from a giant star to a degenerate accretor (default = " + std::to_string(p_Options->m_MassTransferCriticalMassRatioGiantDegenerateAccretor) + ")\n  0 is always stable, <0 is disabled.\n  Only used for --critical-mass-ratio-prescription CLAEYS, ignored otherwise.").c_str()
        )
        (
            "critical-mass-ratio-giant-non-degenerate-accretor",
            po::value<double>(&p_Options->m_MassTransferCriticalMassRatioGiantNonDegenerateAccretor)->default_value(p_Options->m_MassTransferCriticalMassRatioGiantNonDegenerateAccretor),
            ("Critical mass ratio for MT from a giant star to a non-degenerate accretor (default = " + std::to_string(p_Options->m_MassTransferCriticalMassRatioGiantNonDegenerateAccretor) + ", which triggers a call to a function of the core mass ratio [Claeys+2014]).\n  0 is always stable, <0 is disabled.\n  Only used for --critical-mass-ratio-prescription CLAEYS, ignored otherwise.").c_str()  
        )
        (
            "critical-mass-ratio-helium-giant-degenerate-accretor",
            po::value<double>(&p_Options->m_MassTransferCriticalMassRatioHeliumGiantDegenerateAccretor)->default_value(p_Options->m_MassTransferCriticalMassRatioHeliumGiantDegenerateAccretor),
            ("Critical mass ratio for MT from a helium giant star to a degenerate accretor (default = " + std::to_string(p_Options->m_MassTransferCriticalMassRatioHeliumGiantDegenerateAccretor) + ")\n  0 is always stable, <0 is disabled.\n  Only used for --critical-mass-ratio-prescription CLAEYS, ignored otherwise.").c_str()
        )
        (
            "critical-mass-ratio-helium-giant-non-degenerate-accretor",
            po::value<double>(&p_Options->m_MassTransferCriticalMassRatioHeliumGiantNonDegenerateAccretor)->default_value(p_Options->m_MassTransferCriticalMassRatioHeliumGiantNonDegenerateAccretor),
            ("Critical mass ratio for MT from a helium giant star to a non-degenerate accretor (default = " + std::to_string(p_Options->m_MassTransferCriticalMassRatioHeliumGiantNonDegenerateAccretor) + ")\n  0 is always stable, <0 is disabled.\n  Only used for --critical-mass-ratio-prescription CLAEYS, ignored otherwise.").c_str()
        )
        (
            "critical-mass-ratio-helium-HG-degenerate-accretor",
            po::value<double>(&p_Options->m_MassTransferCriticalMassRatioHeliumHGDegenerateAccretor)->default_value(p_Options->m_MassTransferCriticalMassRatioHeliumHGDegenerateAccretor),
            ("Critical mass ratio for MT from a helium HG star to a degenerate accretor (default = " + std::to_string(p_Options->m_MassTransferCriticalMassRatioHeliumHGDegenerateAccretor) + ")\n  0 is always stable, <0 is disabled.\n  Only used for --critical-mass-ratio-prescription CLAEYS, ignored otherwise.").c_str()
        )
        (
            "critical-mass-ratio-helium-HG-non-degenerate-accretor",
            po::value<double>(&p_Options->m_MassTransferCriticalMassRatioHeliumHGNonDegenerateAccretor)->default_value(p_Options->m_MassTransferCriticalMassRatioHeliumHGNonDegenerateAccretor),
            ("Critical mass ratio for MT from a helium HG star to a non-degenerate accretor (default = " + std::to_string(p_Options->m_MassTransferCriticalMassRatioHeliumHGNonDegenerateAccretor) + ")\n  0 is always stable, <0 is disabled.\n  Only used for --critical-mass-ratio-prescription CLAEYS, ignored otherwise.").c_str()
        )
        (
            "critical-mass-ratio-helium-MS-degenerate-accretor",
            po::value<double>(&p_Options->m_MassTransferCriticalMassRatioHeliumMSDegenerateAccretor)->default_value(p_Options->m_MassTransferCriticalMassRatioHeliumMSDegenerateAccretor),
            ("Critical mass ratio for MT from a helium MS star to a degenerate accretor (default = " + std::to_string(p_Options->m_MassTransferCriticalMassRatioHeliumMSDegenerateAccretor) + ")\n  0 is always stable, <0 is disabled.\n  Only used for --critical-mass-ratio-prescription CLAEYS, ignored otherwise.").c_str()
        )
        (
            "critical-mass-ratio-helium-MS-non-degenerate-accretor",
            po::value<double>(&p_Options->m_MassTransferCriticalMassRatioHeliumMSNonDegenerateAccretor)->default_value(p_Options->m_MassTransferCriticalMassRatioHeliumMSNonDegenerateAccretor),
            ("Critical mass ratio for MT from a helium MS star to a non-degenerate accretor (default = " + std::to_string(p_Options->m_MassTransferCriticalMassRatioHeliumMSNonDegenerateAccretor) + ")\n  0 is always stable, <0 is disabled.\n  Only used for --critical-mass-ratio-prescription CLAEYS, ignored otherwise.").c_str()
        )
        (
            "critical-mass-ratio-HG-degenerate-accretor",
            po::value<double>(&p_Options->m_MassTransferCriticalMassRatioHGDegenerateAccretor)->default_value(p_Options->m_MassTransferCriticalMassRatioHGDegenerateAccretor),
            ("Critical mass ratio for MT from a HG star to a degenerate accretor (default = " + std::to_string(p_Options->m_MassTransferCriticalMassRatioHGDegenerateAccretor) + ")\n  0 is always stable, <0 is disabled.\n  Only used for --critical-mass-ratio-prescription CLAEYS, ignored otherwise.").c_str()
        )
        (
            "critical-mass-ratio-HG-non-degenerate-accretor",
            po::value<double>(&p_Options->m_MassTransferCriticalMassRatioHGNonDegenerateAccretor)->default_value(p_Options->m_MassTransferCriticalMassRatioHGNonDegenerateAccretor),
            ("Critical mass ratio for MT from a HG star to a non-degenerate accretor (default = " + std::to_string(p_Options->m_MassTransferCriticalMassRatioHGNonDegenerateAccretor) + ")\n  0 is always stable, <0 is disabled.\n  Only used for --critical-mass-ratio-prescription CLAEYS, ignored otherwise.").c_str()
        )
        (
            "critical-mass-ratio-MS-high-mass-degenerate-accretor",
            po::value<double>(&p_Options->m_MassTransferCriticalMassRatioMSHighMassDegenerateAccretor)->default_value(p_Options->m_MassTransferCriticalMassRatioMSHighMassDegenerateAccretor),
            ("Critical mass ratio for MT from a MS star to a degenerate accretor (default = " + std::to_string(p_Options->m_MassTransferCriticalMassRatioMSHighMassDegenerateAccretor) + ")\n  0 is always stable, <0 is disabled.\n  Only used for --critical-mass-ratio-prescription CLAEYS, ignored otherwise.").c_str()
        )
        (
            "critical-mass-ratio-MS-high-mass-non-degenerate-accretor",
            po::value<double>(&p_Options->m_MassTransferCriticalMassRatioMSHighMassNonDegenerateAccretor)->default_value(p_Options->m_MassTransferCriticalMassRatioMSHighMassNonDegenerateAccretor),
            ("Critical mass ratio for MT from a MS star to a non-degenerate accretor (default = " + std::to_string(p_Options->m_MassTransferCriticalMassRatioMSHighMassNonDegenerateAccretor) + ")\n  0 is always stable, <0 is disabled.\n  Only used for --critical-mass-ratio-prescription CLAEYS, ignored otherwise.").c_str()
        )
        (
            "critical-mass-ratio-MS-low-mass-degenerate-accretor",
            po::value<double>(&p_Options->m_MassTransferCriticalMassRatioMSLowMassDegenerateAccretor)->default_value(p_Options->m_MassTransferCriticalMassRatioMSLowMassDegenerateAccretor),
            ("Critical mass ratio for MT from a MS star to a degenerate accretor (default = " + std::to_string(p_Options->m_MassTransferCriticalMassRatioMSLowMassDegenerateAccretor) + ")\n  0 is always stable, <0 is disabled.\n  Only used for --critical-mass-ratio-prescription CLAEYS, ignored otherwise.").c_str()
        )
        (
            "critical-mass-ratio-MS-low-mass-non-degenerate-accretor",
            po::value<double>(&p_Options->m_MassTransferCriticalMassRatioMSLowMassNonDegenerateAccretor)->default_value(p_Options->m_MassTransferCriticalMassRatioMSLowMassNonDegenerateAccretor),
            ("Critical mass ratio for MT from a MS star to a non-degenerate accretor (default = " + std::to_string(p_Options->m_MassTransferCriticalMassRatioMSLowMassNonDegenerateAccretor) + ")\n  0 is always stable, <0 is disabled.\n  Only used for --critical-mass-ratio-prescription CLAEYS, ignored otherwise.").c_str()
        )
        (
            "critical-mass-ratio-white-dwarf-degenerate-accretor",
            po::value<double>(&p_Options->m_MassTransferCriticalMassRatioWhiteDwarfDegenerateAccretor)->default_value(p_Options->m_MassTransferCriticalMassRatioWhiteDwarfDegenerateAccretor),
            ("Critical mass ratio for MT from a white dwarf to a degenerate accretor (default = " + std::to_string(p_Options->m_MassTransferCriticalMassRatioWhiteDwarfDegenerateAccretor) + ")\n  0 is always stable, <0 is disabled.\n  Only used for --critical-mass-ratio-prescription CLAEYS, ignored otherwise.").c_str()
        )
        (
            "critical-mass-ratio-white-dwarf-non-degenerate-accretor",
            po::value<double>(&p_Options->m_MassTransferCriticalMassRatioWhiteDwarfNonDegenerateAccretor)->default_value(p_Options->m_MassTransferCriticalMassRatioWhiteDwarfNonDegenerateAccretor),
            ("Critical mass ratio for MT from a white dwarf to a non-degenerate accretor (default = " + std::to_string(p_Options->m_MassTransferCriticalMassRatioWhiteDwarfNonDegenerateAccretor) + ")\n  0 is always stable, <0 is disabled.\n  Only used for --critical-mass-ratio-prescription CLAEYS, ignored otherwise.").c_str()
        )

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
            ("Single power law power to generate primary mass using POWERLAW IMF (default = " + std::to_string(p_Options->m_InitialMassFunctionPower) + ")").c_str()
        )

        (
            "kick-direction-power",                                        
            po::value<double>(&p_Options->m_KickDirectionPower)->default_value(p_Options->m_KickDirectionPower),                                                                                  
            ("Power for power law kick direction distribution (default = " + std::to_string(p_Options->m_KickDirectionPower) + " = isotropic, +ve = polar, -ve = in plane)").c_str()
        )
        (
            "kick-magnitude",                                          
            po::value<double>(&p_Options->m_KickMagnitude)->default_value(p_Options->m_KickMagnitude),                                                      
            ("The magnitude of the kick velocity, in km/s, that the star receives during the a supernova (default = " + std::to_string(p_Options->m_KickMagnitude) + ")").c_str()
        )
        (
            "kick-magnitude-1",                                          
            po::value<double>(&p_Options->m_KickMagnitude1)->default_value(p_Options->m_KickMagnitude1),                                                      
            ("The magnitude of the kick velocity, in km/s, that the primary star receives during the a supernova (default = " + std::to_string(p_Options->m_KickMagnitude1) + ")").c_str()
        )
        (
            "kick-magnitude-2",                                          
            po::value<double>(&p_Options->m_KickMagnitude2)->default_value(p_Options->m_KickMagnitude2),                                                      
            ("The magnitude of the kick velocity, in km/s, that the secondary star receives during the a supernova (default = " + std::to_string(p_Options->m_KickMagnitude2) + ")").c_str()
        )
        (
            "kick-magnitude-max",                                          
            po::value<double>(&p_Options->m_KickMagnitudeDistributionMaximum)->default_value(p_Options->m_KickMagnitudeDistributionMaximum),                                                      
            ("Maximum drawn kick magnitude in km/s. Ignored if < 0. Must be > 0 if using kick-magnitude-distribution=FLAT (default = " + std::to_string(p_Options->m_KickMagnitudeDistributionMaximum) + ")").c_str()
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
            ("Sigma for chosen kick magnitude distribution, in km/s, for black holes (default = " + std::to_string(p_Options->m_KickMagnitudeDistributionSigmaCCSN_BH) + ")").c_str()
        )
        (
            "kick-magnitude-sigma-CCSN-NS",                                
            po::value<double>(&p_Options->m_KickMagnitudeDistributionSigmaCCSN_NS)->default_value(p_Options->m_KickMagnitudeDistributionSigmaCCSN_NS),                                            
            ("Sigma for chosen kick magnitude distribution, in km/s, for neutron stars (default = " + std::to_string(p_Options->m_KickMagnitudeDistributionSigmaCCSN_NS) + ")").c_str()
        )
        (
            "kick-magnitude-sigma-ECSN",                                   
            po::value<double>(&p_Options->m_KickMagnitudeDistributionSigmaForECSN)->default_value(p_Options->m_KickMagnitudeDistributionSigmaForECSN),                                            
            ("Sigma for chosen kick magnitude distribution, in km/s, for ECSN (default = " + std::to_string(p_Options->m_KickMagnitudeDistributionSigmaForECSN) + ")").c_str()
        )
        (
            "kick-magnitude-sigma-USSN",                                   
            po::value<double>(&p_Options->m_KickMagnitudeDistributionSigmaForUSSN)->default_value(p_Options->m_KickMagnitudeDistributionSigmaForUSSN),                                            
            ("Sigma for chosen kick magnitude distribution, in km/s, for USSN (default = " + std::to_string(p_Options->m_KickMagnitudeDistributionSigmaForUSSN) + ")").c_str()
        )
        (
            "kick-mean-anomaly-1",
            po::value<double>(&p_Options->m_KickMeanAnomaly1)->default_value(p_Options->m_KickMeanAnomaly1),                                                                                  
            "Mean anomaly, in rad, for the primary star at instantaneous time of the supernova (default = uniform random number [0.0, 2pi))"
        )
        (
            "kick-mean-anomaly-2",
            po::value<double>(&p_Options->m_KickMeanAnomaly2)->default_value(p_Options->m_KickMeanAnomaly2),                                                                                  
            "Mean anomaly, in rad, for the secondary star at instantaneous time of the supernova (default = uniform random number [0.0, 2pi))"
        )
        (
            "kick-phi-1",
            po::value<double>(&p_Options->m_KickPhi1)->default_value(p_Options->m_KickPhi1),                                                                                  
            "Planar angle, in rad, of the supernova vector, for the primary star (default = drawn from kick direction distribution)"
        )
        (
            "kick-phi-2",
            po::value<double>(&p_Options->m_KickPhi2)->default_value(p_Options->m_KickPhi2),                                                                                  
            "Planar angle, in rad, of the supernova vector, for the secondary star (default = drawn from kick direction distribution)"
        )
        (
            "kick-scaling-factor",                                         
            po::value<double>(&p_Options->m_KickScalingFactor)->default_value(p_Options->m_KickScalingFactor),                                                                                    
            ("Arbitrary factor used to scale kicks (default = " + std::to_string(p_Options->m_KickScalingFactor) + ")").c_str()
        )
        (
            "kick-theta-1",                                        
            po::value<double>(&p_Options->m_KickTheta1)->default_value(p_Options->m_KickTheta1),                                                                                  
            "Polar angle, in rad, of the supernova vector, for the primary star (default = drawn from kick direction distribution)"
        )
        (
            "kick-theta-2",                                        
            po::value<double>(&p_Options->m_KickTheta2)->default_value(p_Options->m_KickTheta2),                                                                                  
            "Polar angle, in rad, of the supernova vector, for the secondary star (default = drawn from kick direction distribution)"
        )

        (
            "luminous-blue-variable-multiplier",                           
            po::value<double>(&p_Options->m_LuminousBlueVariableFactor)->default_value(p_Options->m_LuminousBlueVariableFactor),                                                                  
            ("Multiplicitive constant for LBV mass loss (default = " + std::to_string(p_Options->m_LuminousBlueVariableFactor) + ", use 10 for Mennekens & Vanbeveren 2014)").c_str()
        )

        (
            "mass-ratio,q",                                              
            po::value<double>(&p_Options->m_MassRatio)->default_value(p_Options->m_MassRatio),                                                                      
            ("Mass ratio m2/m1 used to determine secondary mass if not specified (default = " + std::to_string(p_Options->m_MassRatio) + ")").c_str()
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
            ("Fraction of specific angular momentum which non-accreted matter removes from the system (default = " + std::to_string(p_Options->m_MassTransferJloss) + ")").c_str()
        )
        (
            "mass-transfer-jloss-macleod-linear-fraction",
            po::value<double>(&p_Options->m_MassTransferJlossMacLeodLinearFraction)->default_value(p_Options->m_MassTransferJlossMacLeodLinearFraction),                                                                                    
            ("Interpolation fraction for jloss prescription if --mass-transfer-angular-momentum-loss-prescription=MACLEOD_LINEAR. 0 is gamma_acc, 1 is gamma_L2 (default = " + std::to_string(p_Options->m_MassTransferJlossMacLeodLinearFraction) + ")").c_str()
        )
        (
            "mass-transfer-thermal-limit-C",                               
            po::value<double>(&p_Options->m_MassTransferCParameter)->default_value(p_Options->m_MassTransferCParameter),                                                                          
            ("Mass Transfer Thermal rate factor of the accretor (default = " + std::to_string(p_Options->m_MassTransferCParameter) + ")").c_str()
        )
        (
            "maximum-evolution-time",                                      
            po::value<double>(&p_Options->m_MaxEvolutionTime)->default_value(p_Options->m_MaxEvolutionTime),                                                                                      
            ("Maximum time to evolve binaries, in Myr (default = " + std::to_string(p_Options->m_MaxEvolutionTime) + ")").c_str()
        )
        (
            "maximum-mass-donor-nandez-ivanova",                           
            po::value<double>(&p_Options->m_MaximumMassDonorNandezIvanova)->default_value(p_Options->m_MaximumMassDonorNandezIvanova),                                                            
            ("Maximum donor mass, in Msol, allowed for the revised common envelope formalism in Msol (default = " + std::to_string(p_Options->m_MaximumMassDonorNandezIvanova) + ")").c_str()
        )
        (
            "maximum-neutron-star-mass",                                   
            po::value<double>(&p_Options->m_MaximumNeutronStarMass)->default_value(p_Options->m_MaximumNeutronStarMass),                                                                          
            ("Maximum mass of a neutron star, in Msol (default = " + std::to_string(p_Options->m_MaximumNeutronStarMass) + ")").c_str()
        )
        (
            "mcbur1",                                                      
            po::value<double>(&p_Options->m_mCBUR1)->default_value(p_Options->m_mCBUR1),                                                                                                          
            ("Minimum core mass at BAGB, in Msol, to avoid fully degenerate CO core  (default = " + std::to_string(p_Options->m_mCBUR1) + ")").c_str()
        )
        (
            "metallicity,z",                                               
            po::value<double>(&p_Options->m_Metallicity)->default_value(p_Options->m_Metallicity),                                                                                                
            ("Metallicity to use (default = " + std::to_string(p_Options->m_Metallicity) + ")").c_str()
        )
        (
            "metallicity-max",                                            
            po::value<double>(&p_Options->m_MetallicityDistributionMax)->default_value(p_Options->m_MetallicityDistributionMax),                                                                
            ("Maximum metallicity to generate (default = " + std::to_string(p_Options->m_MetallicityDistributionMax) + ")").c_str()
        )
        (
            "metallicity-min",                                            
            po::value<double>(&p_Options->m_MetallicityDistributionMin)->default_value(p_Options->m_MetallicityDistributionMin),                                                                
            ("Minimum metallicity to generate (default = " + std::to_string(p_Options->m_MetallicityDistributionMin) + ")").c_str()
        )
        (
            "minimum-secondary-mass",                                      
            po::value<double>(&p_Options->m_MinimumMassSecondary)->default_value(p_Options->m_MinimumMassSecondary),                                                                              
            ("Minimum mass of secondary to generate, in Msol (default = " + std::to_string(p_Options->m_MinimumMassSecondary) + ")").c_str()
        )
        (
            "muller-mandel-kick-multiplier-BH",                                        
            po::value<double>(&p_Options->m_MullerMandelKickBH)->default_value(p_Options->m_MullerMandelKickBH),                                                                                  
            ("Scaling prefactor for BH kicks when using the 'MULLERMANDEL' kick magnitude distribution (default = " + std::to_string(p_Options->m_MullerMandelKickBH) + ")").c_str()
        )
        (
            "muller-mandel-kick-multiplier-NS",                                        
            po::value<double>(&p_Options->m_MullerMandelKickNS)->default_value(p_Options->m_MullerMandelKickNS),                                                                                  
            ("Scaling prefactor for NS kicks when using the 'MULLERMANDEL' kick magnitude distribution (default = " + std::to_string(p_Options->m_MullerMandelKickNS) + ")").c_str()
        )
        (
            "muller-mandel-sigma-kick",                                        
            po::value<double>(&p_Options->m_MullerMandelSigmaKick)->default_value(p_Options->m_MullerMandelSigmaKick),                                                                                  
            ("Kick scatter when using the 'MULLERMANDEL' kick magnitude distribution (default = " + std::to_string(p_Options->m_MullerMandelSigmaKick) + ")").c_str()
        )

        (
            "neutrino-mass-loss-BH-formation-value",                       
            po::value<double>(&p_Options->m_NeutrinoMassLossValueBH)->default_value(p_Options->m_NeutrinoMassLossValueBH),                                                                        
            ("Amount of BH mass lost due to neutrinos (either fraction or fixed value, depending on --neutrino-mass-loss-BH-formation) (default = " + std::to_string(p_Options->m_NeutrinoMassLossValueBH) + ")").c_str()
        )

        (
            "orbital-period",                                          
            po::value<double>(&p_Options->m_OrbitalPeriod)->default_value(p_Options->m_OrbitalPeriod),                                                                            
            ("Initial orbital period, in days (default = " + std::to_string(p_Options->m_OrbitalPeriod) + ")").c_str()
        )
        (
            "orbital-period-max",                                          
            po::value<double>(&p_Options->m_OrbitalPeriodDistributionMax)->default_value(p_Options->m_OrbitalPeriodDistributionMax),                                                                            
            ("Maximum period, in days, to generate (default = " + std::to_string(p_Options->m_OrbitalPeriodDistributionMax) + ")").c_str()
        )
        (
            "orbital-period-min",                                          
            po::value<double>(&p_Options->m_OrbitalPeriodDistributionMin)->default_value(p_Options->m_OrbitalPeriodDistributionMin),                                                                            
            ("Minimum period, in days, to generate (default = " + std::to_string(p_Options->m_OrbitalPeriodDistributionMin) + ")").c_str()
        )

        (
            "overall-wind-mass-loss-multiplier",                           
            po::value<double>(&p_Options->m_OverallWindMassLossMultiplier)->default_value(p_Options->m_OverallWindMassLossMultiplier),                                                                  
            ("Multiplicitive constant for overall wind mass loss (default = " + std::to_string(p_Options->m_OverallWindMassLossMultiplier)+ ")").c_str()
        )

        (
            "PISN-lower-limit",                                            
            po::value<double>(&p_Options->m_PairInstabilityLowerLimit)->default_value(p_Options->m_PairInstabilityLowerLimit),                                                                    
            ("Minimum core mass for PISN, in Msol (default = " + std::to_string(p_Options->m_PairInstabilityLowerLimit) + ")").c_str()
        )
        (
            "PISN-upper-limit",                                            
            po::value<double>(&p_Options->m_PairInstabilityUpperLimit)->default_value(p_Options->m_PairInstabilityUpperLimit),                                                                    
            ("Maximum core mass for PISN, in Msol (default = " + std::to_string(p_Options->m_PairInstabilityUpperLimit) + ")").c_str()
        )
        (
            "PPI-lower-limit",                                             
            po::value<double>(&p_Options->m_PulsationalPairInstabilityLowerLimit)->default_value(p_Options->m_PulsationalPairInstabilityLowerLimit),                                              
            ("Minimum core mass for PPI, in Msol (default = " + std::to_string(p_Options->m_PulsationalPairInstabilityLowerLimit) + ")").c_str()
        )
        (
            "PPI-upper-limit",                                             
            po::value<double>(&p_Options->m_PulsationalPairInstabilityUpperLimit)->default_value(p_Options->m_PulsationalPairInstabilityUpperLimit),                                              
            ("Maximum core mass for PPI, in Msol (default = " + std::to_string(p_Options->m_PulsationalPairInstabilityUpperLimit) + ")").c_str()
        )
        (
            "pulsar-birth-magnetic-field-distribution-max",                
            po::value<double>(&p_Options->m_PulsarBirthMagneticFieldDistributionMax)->default_value(p_Options->m_PulsarBirthMagneticFieldDistributionMax),                                        
            ("Maximum pulsar birth magnetic field, in log10(Gauss) (default = " + std::to_string(p_Options->m_PulsarBirthMagneticFieldDistributionMax) + ")").c_str()
        )
        (
            "pulsar-birth-magnetic-field-distribution-min",                
            po::value<double>(&p_Options->m_PulsarBirthMagneticFieldDistributionMin)->default_value(p_Options->m_PulsarBirthMagneticFieldDistributionMin),                                        
            ("Minimum pulsar birth magnetic field, in log10(Gauss) (default = " + std::to_string(p_Options->m_PulsarBirthMagneticFieldDistributionMin) + ")").c_str()
        )
        (
            "pulsar-birth-spin-period-distribution-max",                   
            po::value<double>(&p_Options->m_PulsarBirthSpinPeriodDistributionMax)->default_value(p_Options->m_PulsarBirthSpinPeriodDistributionMax),                                              
            ("Maximum pulsar birth spin period, in ms (default = " + std::to_string(p_Options->m_PulsarBirthSpinPeriodDistributionMax) + ")").c_str()
        )
        (
            "pulsar-birth-spin-period-distribution-min",                   
            po::value<double>(&p_Options->m_PulsarBirthSpinPeriodDistributionMin)->default_value(p_Options->m_PulsarBirthSpinPeriodDistributionMin),                                              
            ("Minimum pulsar birth spin period, in ms (default = " + std::to_string(p_Options->m_PulsarBirthSpinPeriodDistributionMin) + ")").c_str()
        )
        (
            "pulsar-magnetic-field-decay-massscale",                       
            po::value<double>(&p_Options->m_PulsarMagneticFieldDecayMassscale)->default_value(p_Options->m_PulsarMagneticFieldDecayMassscale),                                                    
            ("Mass scale on which magnetic field decays during accretion, in Msol (default = " + std::to_string(p_Options->m_PulsarMagneticFieldDecayMassscale) + ")").c_str()
        )
        (
            "pulsar-magnetic-field-decay-timescale",                       
            po::value<double>(&p_Options->m_PulsarMagneticFieldDecayTimescale)->default_value(p_Options->m_PulsarMagneticFieldDecayTimescale),                                                    
            ("Timescale on which magnetic field decays, in Myrs (default = " + std::to_string(p_Options->m_PulsarMagneticFieldDecayTimescale) + ")").c_str()
        )
        (
            "pulsar-minimum-magnetic-field",                               
            po::value<double>(&p_Options->m_PulsarLog10MinimumMagneticField)->default_value(p_Options->m_PulsarLog10MinimumMagneticField),                                                        
            ("Minimum pulsar magnetic field, in log10(Gauss) (default = " + std::to_string(p_Options->m_PulsarLog10MinimumMagneticField) + ")").c_str()
        )

        (
            "rotational-frequency",                              
            po::value<double>(&p_Options->m_RotationalFrequency)->default_value(p_Options->m_RotationalFrequency),                                                        
            ("Initial rotational frequency for the star for SSE (Hz) (default = " + std::to_string(p_Options->m_RotationalFrequency) + ")").c_str()
        )        

        (
            "rotational-frequency-1",                              
            po::value<double>(&p_Options->m_RotationalFrequency1)->default_value(p_Options->m_RotationalFrequency1),                                                        
            ("Initial rotational frequency for the primary star for BSE (Hz) (default = " + std::to_string(p_Options->m_RotationalFrequency1) + ")").c_str()
        )        

        (
            "rotational-frequency-2",                              
            po::value<double>(&p_Options->m_RotationalFrequency2)->default_value(p_Options->m_RotationalFrequency2),                                                        
            ("Initial rotational frequency for the secondary star for BSE (Hz) (default = " + std::to_string(p_Options->m_RotationalFrequency2) + ")").c_str()
        )        

        (
            "semi-major-axis,a",                              
            po::value<double>(&p_Options->m_SemiMajorAxis)->default_value(p_Options->m_SemiMajorAxis),                                                        
            ("Initial semi-major axis, in AU (default = " + std::to_string(p_Options->m_SemiMajorAxis) + ")").c_str()
        )        
        (
            "semi-major-axis-max",                                         
            po::value<double>(&p_Options->m_SemiMajorAxisDistributionMax)->default_value(p_Options->m_SemiMajorAxisDistributionMax),                                                              
            ("Maximum semi-major axis, in AU, to generate (default = " + std::to_string(p_Options->m_SemiMajorAxisDistributionMax) + ")").c_str()
        )
        (
            "semi-major-axis-min",                                         
            po::value<double>(&p_Options->m_SemiMajorAxisDistributionMin)->default_value(p_Options->m_SemiMajorAxisDistributionMin),                                                              
            ("Minimum semi-major axis, in AU, to generate (default = " + std::to_string(p_Options->m_SemiMajorAxisDistributionMin) + ")").c_str()
        )

        (
            "timestep-multiplier",
            po::value<double>(&p_Options->m_TimestepMultiplier)->default_value(p_Options->m_TimestepMultiplier),
            ("Timestep multiplier for SSE and BSE (default = " + std::to_string(p_Options->m_TimestepMultiplier) + ")").c_str()
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
            "add-options-to-sysparms",                                            
            po::value<std::string>(&p_Options->m_AddOptionsToSysParms.typeString)->default_value(p_Options->m_AddOptionsToSysParms.typeString),                                                                              
            ("Add program options columns to BSE/SSE SysParms file (options: [ALWAYS, GRID, NEVER], default = " + p_Options->m_AddOptionsToSysParms.typeString + ")").c_str()
        )

        (
            "black-hole-kicks",                                            
            po::value<std::string>(&p_Options->m_BlackHoleKicks.typeString)->default_value(p_Options->m_BlackHoleKicks.typeString),                                                                              
            ("Black hole kicks relative to NS kicks (options: [FULL, REDUCED, ZERO, FALLBACK], default = " + p_Options->m_BlackHoleKicks.typeString + ")").c_str()
        )

        (
            "case-BB-stability-prescription",                              
            po::value<std::string>(&p_Options->m_CaseBBStabilityPrescription.typeString)->default_value(p_Options->m_CaseBBStabilityPrescription.typeString),                                                    
            ("Case BB/BC mass transfer stability prescription (options: [ALWAYS_STABLE, ALWAYS_STABLE_ONTO_NSBH, TREAT_AS_OTHER_MT, ALWAYS_UNSTABLE], default = " + p_Options->m_CaseBBStabilityPrescription.typeString + ")").c_str()
        )
        (
            "chemically-homogeneous-evolution",                            
            po::value<std::string>(&p_Options->m_CheMode.typeString)->default_value(p_Options->m_CheMode.typeString),                                                                                                    
            ("Chemically Homogeneous Evolution (options: [NONE, OPTIMISTIC, PESSIMISTIC], default = " + p_Options->m_CheMode.typeString + ")").c_str()
        )
        (
            "common-envelope-lambda-prescription",                         
            po::value<std::string>(&p_Options->m_CommonEnvelopeLambdaPrescription.typeString)->default_value(p_Options->m_CommonEnvelopeLambdaPrescription.typeString),                                          
            ("CE lambda prescription (options: [LAMBDA_FIXED, LAMBDA_LOVERIDGE, LAMBDA_NANJING, LAMBDA_KRUCKOW, LAMBDA_DEWI], default = " + p_Options->m_CommonEnvelopeLambdaPrescription.typeString + ")").c_str()
        )
        (
            "common-envelope-mass-accretion-prescription",                 
            po::value<std::string>(&p_Options->m_CommonEnvelopeMassAccretionPrescription.typeString)->default_value(p_Options->m_CommonEnvelopeMassAccretionPrescription.typeString),                            
            ("Assumption about whether NS/BHs can accrete mass during common envelope evolution (options: [ZERO, CONSTANT, UNIFORM, MACLEOD], default = " + p_Options->m_CommonEnvelopeMassAccretionPrescription.typeString + ")").c_str()
        )

        (
            "critical-mass-ratio-prescription",                                 
            po::value<std::string>(&p_Options->m_QCritPrescription.typeString)->default_value(p_Options->m_QCritPrescription.typeString),
            ("Prescription for which critical mass ratio prescription to use, if any (options: [NONE, CLAEYS, GE15, GE15_IC], default = " + p_Options->m_QCritPrescription.typeString + ")").c_str()
        )
        
        (
            "eccentricity-distribution",                                 
            po::value<std::string>(&p_Options->m_EccentricityDistribution.typeString)->default_value(p_Options->m_EccentricityDistribution.typeString),                                                          
            ("Initial eccentricity distribution (options: [ZERO, FLAT, THERMAL, GELLER+2013, DUQUENNOYMAYOR1991, SANA2012], default = " + p_Options->m_EccentricityDistribution.typeString + ")").c_str()
        )
        (
            "envelope-state-prescription",                                 
            po::value<std::string>(&p_Options->m_EnvelopeStatePrescription.typeString)->default_value(p_Options->m_EnvelopeStatePrescription.typeString),                                                        
            ("Prescription for whether the envelope is radiative or convective (options: [LEGACY, HURLEY, FIXED_TEMPERATURE], default = " + p_Options->m_EnvelopeStatePrescription.typeString + ")").c_str()
        )

        (
            "fryer-supernova-engine",                                      
            po::value<std::string>(&p_Options->m_FryerSupernovaEngine.typeString)->default_value(p_Options->m_FryerSupernovaEngine.typeString),                                                                  
            ("If using Fryer et al 2012 fallback prescription. (options: [DELAYED, RAPID], default = " + p_Options->m_FryerSupernovaEngine.typeString + ")").c_str()
        )
        (
            "fryer-22-fmix",                                        
            po::value<double>(&p_Options->m_Fryer22fmix)->default_value(p_Options->m_Fryer22fmix),                                                                                  
            ("paramter describing the mixing growth time when using the 'FRYER2022' remnant mass distribution (default = " + std::to_string(p_Options->m_Fryer22fmix) + ")").c_str()
        )
        (
            "fryer-22-mcrit",                                        
            po::value<double>(&p_Options->m_Fryer22Mcrit)->default_value(p_Options->m_Fryer22Mcrit),                                                                                  
            ("Critical CO core mass for black hole formation when using the 'FRYER2022' remnant mass distribution (default = " + std::to_string(p_Options->m_Fryer22Mcrit) + ")").c_str()
        )

        (
            "grid",                                                        
            po::value<std::string>(&p_Options->m_GridFilename)->default_value(p_Options->m_GridFilename)->implicit_value(""),
            ("Grid filename (default = " + p_Options->m_GridFilename + ")").c_str()
        )

        (
            "initial-mass-function,i",                                     
            po::value<std::string>(&p_Options->m_InitialMassFunction.typeString)->default_value(p_Options->m_InitialMassFunction.typeString),                                                                    
            ("Initial mass function (options: [SALPETER, POWERLAW, UNIFORM, KROUPA], default = " + p_Options->m_InitialMassFunction.typeString + ")").c_str()
        )

        (
            "kick-direction",                                              
            po::value<std::string>(&p_Options->m_KickDirectionDistribution.typeString)->default_value(p_Options->m_KickDirectionDistribution.typeString),                                                        
            ("Natal kick direction distribution (options: [ISOTROPIC, INPLANE, PERPENDICULAR, POWERLAW, WEDGE, POLES], default = " + p_Options->m_KickDirectionDistribution.typeString + ")").c_str()
        )
        (
            "kick-magnitude-distribution",                                 
            po::value<std::string>(&p_Options->m_KickMagnitudeDistribution.typeString)->default_value(p_Options->m_KickMagnitudeDistribution.typeString),                                                        
            ("Natal kick magnitude distribution (options: [ZERO, FIXED, FLAT, MAXWELLIAN, BRAYELDRIDGE, MULLER2016, MULLER2016MAXWELLIAN, MULLERMANDEL], default = " + p_Options->m_KickMagnitudeDistribution.typeString + ")").c_str()
        )

        /*
        (
            "logfile-BE-binaries",                                     
            po::value<std::string>(&p_Options->m_LogfileBeBinaries)->default_value(p_Options->m_LogfileBeBinaries),                                                                              
            ("Filename for BSE Be Binaries logfile (default = " + p_Options->m_LogfileBeBinaries + ")").c_str()
        )
        */

        (
            "logfile-rlof-parameters",                                 
            po::value<std::string>(&p_Options->m_LogfileRLOFParameters)->default_value(p_Options->m_LogfileRLOFParameters),                                                                      
            ("Filename for BSE RLOF Parameters logfile ( default = " + p_Options->m_LogfileRLOFParameters + ")").c_str()
        )
        (
            "logfile-common-envelopes",                                
            po::value<std::string>(&p_Options->m_LogfileCommonEnvelopes)->default_value(p_Options->m_LogfileCommonEnvelopes),                                                                    
            ("Filename for BSE Common Envelopes logfile (default = " + p_Options->m_LogfileCommonEnvelopes + ")").c_str()
        )
        (
            "logfile-detailed-output",                                 
            po::value<std::string>(&p_Options->m_LogfileDetailedOutput)->default_value(p_Options->m_LogfileDetailedOutput),                                                                      
            ("Filename for BSE Detailed Output logfile (default = " + p_Options->m_LogfileDetailedOutput + ")").c_str()
        )
        (
            "logfile-double-compact-objects",                          
            po::value<std::string>(&p_Options->m_LogfileDoubleCompactObjects)->default_value(p_Options->m_LogfileDoubleCompactObjects),                                                          
            ("Filename for Double Compact Objects logfile (default = " + p_Options->m_LogfileDoubleCompactObjects + ")").c_str()
        )
        (
            "logfile-pulsar-evolution",                                
            po::value<std::string>(&p_Options->m_LogfilePulsarEvolution)->default_value(p_Options->m_LogfilePulsarEvolution),                                                                    
            ("Filename for Pulsar Evolution logfile (default = " + p_Options->m_LogfilePulsarEvolution + ")").c_str()
        )
        (
            "logfile-supernovae",                                      
            po::value<std::string>(&p_Options->m_LogfileSupernovae)->default_value(p_Options->m_LogfileSupernovae),                                                                              
            ("Filename for Supernovae logfile (default = " + p_Options->m_LogfileSupernovae + ")").c_str()
        )
        (
            "logfile-system-parameters",                               
            po::value<std::string>(&p_Options->m_LogfileSystemParameters)->default_value(p_Options->m_LogfileSystemParameters),                                                                  
            ("Filename for System Parameters logfile (default = " + p_Options->m_LogfileSystemParameters + ")").c_str()
        )
        (
            "logfile-definitions",                                         
            po::value<std::string>(&p_Options->m_LogfileDefinitionsFilename)->default_value(p_Options->m_LogfileDefinitionsFilename)->implicit_value(""),                                              
            ("Filename for logfile record definitions (default = " + p_Options->m_LogfileDefinitionsFilename + ")").c_str()
        )
        (
            "logfile-name-prefix",                                         
            po::value<std::string>(&p_Options->m_LogfileNamePrefix)->default_value(p_Options->m_LogfileNamePrefix)->implicit_value(""),                                                                
            ("Prefix for logfile names (default = " + p_Options->m_LogfileNamePrefix + ")").c_str()
        )
        (
            "logfile-switch-log",                                      
            po::value<std::string>(&p_Options->m_LogfileSwitchLog)->default_value(p_Options->m_LogfileSwitchLog),                                                                                
            ("Filename for Switch Log logfile (default = " + p_Options->m_LogfileSwitchLog + ")").c_str()
        )
        (
            "logfile-type",                                           
            po::value<std::string>(&p_Options->m_LogfileType.typeString)->default_value(p_Options->m_LogfileType.typeString),                                                                          
            ("File type for logfiles (options: [HDF5, CSV, TSV, TXT], default = " + p_Options->m_LogfileType.typeString + ")").c_str()
        )

        (
            "luminous-blue-variable-prescription",                                      
            po::value<std::string>(&p_Options->m_LuminousBlueVariablePrescription.typeString)->default_value(p_Options->m_LuminousBlueVariablePrescription.typeString),                                                                  
            ("LBV Mass loss prescription (options: [NONE, HURLEY_ADD, HURLEY, BELCZYNSKI], default = " + p_Options->m_LuminousBlueVariablePrescription.typeString + ")").c_str()
        )
        (
            "mass-loss-prescription",                                      
            po::value<std::string>(&p_Options->m_MassLossPrescription.typeString)->default_value(p_Options->m_MassLossPrescription.typeString),                                                                  
            ("Mass loss prescription (options: [NONE, HURLEY, VINK], default = " + p_Options->m_MassLossPrescription.typeString + ")").c_str()
        )
        (
            "mass-ratio-distribution",                                   
            po::value<std::string>(&p_Options->m_MassRatioDistribution.typeString)->default_value(p_Options->m_MassRatioDistribution.typeString),                                                                
            ("Initial mass ratio distribution for q=m2/m1 (options: [FLAT, DUQUENNOYMAYOR1991, SANA2012], default = " + p_Options->m_MassRatioDistribution.typeString + ")").c_str()
        )
        (
            "mass-transfer-accretion-efficiency-prescription",             
            po::value<std::string>(&p_Options->m_MassTransferAccretionEfficiencyPrescription.typeString)->default_value(p_Options->m_MassTransferAccretionEfficiencyPrescription.typeString),                    
            ("Mass Transfer Accretion Efficiency prescription (options: [THERMAL, FIXED], default = " + p_Options->m_MassTransferAccretionEfficiencyPrescription.typeString + ")").c_str()
        )
        (
            "mass-transfer-angular-momentum-loss-prescription",            
            po::value<std::string>(&p_Options->m_MassTransferAngularMomentumLossPrescription.typeString)->default_value(p_Options->m_MassTransferAngularMomentumLossPrescription.typeString),                    
            ("Mass Transfer Angular Momentum Loss prescription (options: [JEANS, ISOTROPIC, CIRCUMBINARY, MACLEOD_LINEAR, ARBITRARY], default = " + p_Options->m_MassTransferAngularMomentumLossPrescription.typeString + ")").c_str()
        )
        (
            "mass-transfer-rejuvenation-prescription",                     
            po::value<std::string>(&p_Options->m_MassTransferRejuvenationPrescription.typeString)->default_value(p_Options->m_MassTransferRejuvenationPrescription.typeString),                                  
            ("Mass Transfer Rejuvenation prescription (options: [NONE, STARTRACK], default = " + p_Options->m_MassTransferRejuvenationPrescription.typeString + ")").c_str()
        )
        (
            "mass-transfer-thermal-limit-accretor",                        
            po::value<std::string>(&p_Options->m_MassTransferThermallyLimitedVariation.typeString)->default_value(p_Options->m_MassTransferThermallyLimitedVariation.typeString),                                
            ("Mass Transfer Thermal Accretion limit (options: [CFACTOR, ROCHELOBE], default = " + p_Options->m_MassTransferThermallyLimitedVariation.typeString + ")").c_str()
        )
        (
            "metallicity-distribution",                                 
            po::value<std::string>(&p_Options->m_MetallicityDistribution.typeString)->default_value(p_Options->m_MetallicityDistribution.typeString),                                                          
            ("Metallicity distribution (options: [ZSOLAR, LOGUNIFORM], default = " + p_Options->m_MetallicityDistribution.typeString + ")").c_str()
        )
        (
            "mode",                                                 
            po::value<std::string>(&p_Options->m_EvolutionMode.typeString)->default_value(p_Options->m_EvolutionMode.typeString),                                                                              
            ("Evolution mode (options: [SSE, BSE], default = " + p_Options->m_EvolutionMode.typeString + ")").c_str()
        )

        (
            "neutrino-mass-loss-BH-formation",                             
            po::value<std::string>(&p_Options->m_NeutrinoMassLossAssumptionBH.typeString)->default_value(p_Options->m_NeutrinoMassLossAssumptionBH.typeString),                                                  
            ("Assumption about neutrino mass loss during BH formation (options: [FIXED_FRACTION, FIXED_MASS], default = " + p_Options->m_NeutrinoMassLossAssumptionBH.typeString + ")").c_str()
        )
        (
            "neutron-star-equation-of-state",                              
            po::value<std::string>(&p_Options->m_NeutronStarEquationOfState.typeString)->default_value(p_Options->m_NeutronStarEquationOfState.typeString),                                                      
            ("Neutron star equation of state to use (options: [SSE, ARP3], default = " + p_Options->m_NeutronStarEquationOfState.typeString + ")").c_str()
        )

        (
            "orbital-period-distribution",                              
            po::value<std::string>(&p_Options->m_OrbitalPeriodDistribution.typeString)->default_value(p_Options->m_OrbitalPeriodDistribution.typeString),                                                        
            ("Initial orbital period distribution (options: [FLATINLOG], default = " + p_Options->m_OrbitalPeriodDistribution.typeString + ")").c_str()
        )        
        (
            "output-container,c",                                          
            po::value<std::string>(&p_Options->m_OutputContainerName)->default_value(p_Options->m_OutputContainerName)->implicit_value(""),                                                            
            ("Container (directory) name for output files (default = " + p_Options->m_OutputContainerName + ")").c_str()
        )
        (
            "output-path,o",                                                
            po::value<std::string>(&p_Options->m_OutputPathString)->default_value(p_Options->m_OutputPathString)->implicit_value(""),                                                                  
            ("Directory for output (default = " + p_Options->m_OutputPathString + ")").c_str()
        )

        (
            "pulsar-birth-magnetic-field-distribution",                    
            po::value<std::string>(&p_Options->m_PulsarBirthMagneticFieldDistribution.typeString)->default_value(p_Options->m_PulsarBirthMagneticFieldDistribution.typeString),                                  
            ("Pulsar Birth Magnetic Field distribution (options: [ZERO, FIXED, FLATINLOG, UNIFORM, LOGNORMAL], default = " + p_Options->m_PulsarBirthMagneticFieldDistribution.typeString + ")").c_str()
        )
        (
            "pulsar-birth-spin-period-distribution",                       
            po::value<std::string>(&p_Options->m_PulsarBirthSpinPeriodDistribution.typeString)->default_value(p_Options->m_PulsarBirthSpinPeriodDistribution.typeString),                                        
            ("Pulsar Birth Spin Period distribution (options: [ZERO, FIXED, UNIFORM, NORMAL], default = " + p_Options->m_PulsarBirthSpinPeriodDistribution.typeString + ")").c_str()
        )
        (
            "pulsational-pair-instability-prescription",                   
            po::value<std::string>(&p_Options->m_PulsationalPairInstabilityPrescription.typeString)->default_value(p_Options->m_PulsationalPairInstabilityPrescription.typeString),                              
            ("Pulsational Pair Instability prescription (options: [COMPAS, STARTRACK, MARCHANT, FARMER], default = " + p_Options->m_PulsationalPairInstabilityPrescription.typeString + ")").c_str()
        )

        (
            "remnant-mass-prescription",                                   
            po::value<std::string>(&p_Options->m_RemnantMassPrescription.typeString)->default_value(p_Options->m_RemnantMassPrescription.typeString),                                                            
            ("Choose remnant mass prescription (options: [HURLEY2000, BELCZYNSKI2002, FRYER2012, FRYER2022, MULLER2016, MULLERMANDEL, SCHNEIDER2020, SCHNEIDER2020ALT], default = " + p_Options->m_RemnantMassPrescription.typeString + ")").c_str()
        )
        (
            "rotational-velocity-distribution",                            
            po::value<std::string>(&p_Options->m_RotationalVelocityDistribution.typeString)->default_value(p_Options->m_RotationalVelocityDistribution.typeString),                                              
            ("Initial rotational velocity distribution (options: [ZERO, HURLEY, VLTFLAMES], default = " + p_Options->m_RotationalVelocityDistribution.typeString + ")").c_str()
        )

        (
            "semi-major-axis-distribution",                              
            po::value<std::string>(&p_Options->m_SemiMajorAxisDistribution.typeString)->default_value(p_Options->m_SemiMajorAxisDistribution.typeString),                                                        
            ("Initial semi-major axis distribution (options: [FLATINLOG, DUQUENNOYMAYOR1991, SANA2012], default = " + p_Options->m_SemiMajorAxisDistribution.typeString + ")").c_str()
        )        
        (
            "stellar-zeta-prescription",                                   
            po::value<std::string>(&p_Options->m_StellarZetaPrescription.typeString)->default_value(p_Options->m_StellarZetaPrescription.typeString),                                                            
            ("Prescription for stellar zeta (default = " + p_Options->m_StellarZetaPrescription.typeString + ")").c_str()
        )


        // vector (list) options - alphabetically

        (
            "debug-classes",                                               
            po::value<std::vector<std::string>>(&p_Options->m_DebugClasses)->multitoken()->default_value(p_Options->m_DebugClasses),                                                                        
            ("Debug classes enabled (default = " + defaultDebugClasses + ")").c_str()
        )

        (
            "log-classes",                                                 
            po::value<std::vector<std::string>>(&p_Options->m_LogClasses)->multitoken()->default_value(p_Options->m_LogClasses),                                                                            
            ("Logging classes enabled (default = " + defaultLogClasses + ")").c_str()
        )

        (
            "notes",                                                 
            po::value<std::vector<std::string>>(&p_Options->m_Notes)->multitoken()->default_value(p_Options->m_Notes),                                                                            
            ("User-specified annotations (default = " + defaultNotes + ")").c_str()
        )
        (
            "notes-hdrs",                                                 
            po::value<std::vector<std::string>>(&p_Options->m_NotesHdrs)->multitoken()->default_value(p_Options->m_NotesHdrs),                                                                            
            ("User-specified annotation header strings (default = " + defaultNotesHdrs + ")").c_str()
        )
    
        ;   // end the list of options to be added

    }
    catch (po::error& e) {      // program options exception
        ok = false;             // set status
    }
    catch (...) {               // unhandled exception
        ok = false;             // set status
    }

    return ok;
}


/*
 * Sets new values for options that are calculated or drawn from distributions
 * This is broken out into this function so that that it can be called each time
 * an options "variation" is advanced
 * 
 * Note this is a class OptionValues function.
 * 
 * 
 * std::string SetCalculatedOptionDefaults(const BOOST_MAP p_ModifyMap)
 * 
 * @param   [IN]    p_UpdateMap                 Flag indicating whether the Boost variables map should be updated
 * @return                                      String containing an error string
 *                                              If no error occurred the return string will be the empty string 
 */
std::string Options::OptionValues::SetCalculatedOptionDefaults(const BOOST_MAP p_UpdateMap) {
#define DEFAULTED(opt) m_VM[opt].defaulted()    // for convenience and readability - undefined at end of function

    std::string errStr = "";                                        // error string

    try {

        // set "default" values for magnitude random number
        // only set if the user did not specify a value

        if (DEFAULTED("kick-magnitude-random")) {
            m_KickMagnitudeRandom = RAND->Random();
            if (p_UpdateMap == BOOST_MAP::UPDATE) {
                ModifyVariableMap(m_VM, "kick-magnitude-random", m_KickMagnitudeRandom);
                po::notify(m_VM);
            }
        }
    
        if (DEFAULTED("kick-magnitude-random-1")) {
            m_KickMagnitudeRandom1 = RAND->Random();
            if (p_UpdateMap == BOOST_MAP::UPDATE) {
                ModifyVariableMap(m_VM, "kick-magnitude-random-1", m_KickMagnitudeRandom1);
                po::notify(m_VM);
            }
        }
    
        if (DEFAULTED("kick-magnitude-random-2")) {
            m_KickMagnitudeRandom2 = RAND->Random();
            if (p_UpdateMap == BOOST_MAP::UPDATE) {
                ModifyVariableMap(m_VM, "kick-magnitude-random-2", m_KickMagnitudeRandom2);
                po::notify(m_VM);
            }
        }

        // set "default" values for mean anomaly
        // only set if the user did not specify a value
    
        if (DEFAULTED("kick-mean-anomaly-1")) {
            m_KickMeanAnomaly1 = RAND->Random(0.0, _2_PI);
            if (p_UpdateMap == BOOST_MAP::UPDATE) {
                ModifyVariableMap(m_VM, "kick-mean-anomaly-1", m_KickMeanAnomaly1);
                po::notify(m_VM);
            }
        }

        if (DEFAULTED("kick-mean-anomaly-2")) {
            m_KickMeanAnomaly2 = RAND->Random(0.0, _2_PI);
            if (p_UpdateMap == BOOST_MAP::UPDATE) {
                ModifyVariableMap(m_VM, "kick-mean-anomaly-2", m_KickMeanAnomaly2);
                po::notify(m_VM);
            }
        }

        // set "default" values for m_KickPhi[1/2] and m_KickTheta[1/2]
        // we now have the kick direction distribution and kick direction 
        // power (exponent) required by the user (either default or specified)
        // only set if the user did not specify a value

        bool phi1Defaulted   = DEFAULTED("kick-phi-1");
        bool theta1Defaulted = DEFAULTED("kick-theta-1");

        if (phi1Defaulted || theta1Defaulted) {
            double phi1, theta1;
            std::tie(phi1, theta1) = utils::DrawKickDirection(m_KickDirectionDistribution.type, m_KickDirectionPower);
            if (phi1Defaulted) {
                m_KickPhi1 = phi1;
                if (p_UpdateMap == BOOST_MAP::UPDATE) {
                    ModifyVariableMap(m_VM, "kick-phi-1", m_KickPhi1);
                    po::notify(m_VM);
                }
            }

            if (theta1Defaulted) {
                m_KickTheta1 = theta1;
                if (p_UpdateMap == BOOST_MAP::UPDATE) {
                    ModifyVariableMap(m_VM, "kick-theta-1", m_KickTheta1);
                    po::notify(m_VM);
                }
            }
        }

        bool phi2Defaulted   = DEFAULTED("kick-phi-2");
        bool theta2Defaulted = DEFAULTED("kick-theta-2");

        if (phi2Defaulted || theta2Defaulted) {
            double phi2, theta2;
            std::tie(phi2, theta2) = utils::DrawKickDirection(m_KickDirectionDistribution.type, m_KickDirectionPower);
            if (phi2Defaulted) {
                m_KickPhi2 = phi2;
                if (p_UpdateMap == BOOST_MAP::UPDATE) {
                    ModifyVariableMap(m_VM, "kick-phi-2", m_KickPhi2);
                    po::notify(m_VM);
                }
            }

            if (theta2Defaulted) {
                m_KickTheta2 = theta2;
                if (p_UpdateMap == BOOST_MAP::UPDATE) {
                    ModifyVariableMap(m_VM, "kick-theta-2", m_KickTheta2);
                    po::notify(m_VM);
                }
            }
        }
    }
    catch (po::error& e) {                                          // program options exception
        errStr = e.what();                                          // set the error string
    }
    catch (const std::string eStr) {                                // custom exception
        errStr = eStr;                                              // set the error string
    }
    catch (...) {                                                   // unhandled exception
        errStr = ERR_MSG(ERROR::UNHANDLED_EXCEPTION);               // set the error string
    }

    return errStr;
#undef DEFAULTED
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
 * std::string Options::OptionValues::CheckAndSetOptions()
 * 
 * @return                                      String containing an error string
 *                                              If no error occurred the return string will be the empty string 
 */
std::string Options::OptionValues::CheckAndSetOptions() {
#define DEFAULTED(opt) m_VM[opt].defaulted()    // for convenience and readability - undefined at end of function

    std::string errStr = "";                                                                                                        // error string

    // check & set prescriptions, distributions, assumptions etc. options - alphabetically

    try {

        bool found;

        m_FixedRandomSeed  = !DEFAULTED("random-seed");                                                                             // use random seed if it is provided by the user
        m_UseFixedUK       = !DEFAULTED("fix-dimensionless-kick-magnitude") && (m_FixedUK >= 0.0);                                  // determine if user supplied a valid kick magnitude

        if (!DEFAULTED("add-options_to-sysparms")) {                                                                                // add program options to BSE/SSE sysparms
            std::tie(found, m_AddOptionsToSysParms.type) = utils::GetMapKey(m_AddOptionsToSysParms.typeString, ADD_OPTIONS_TO_SYSPARMS_LABEL, m_AddOptionsToSysParms.type);
            COMPLAIN_IF(!found, "Unknown Add Options to SysParms Option");
        }

        if (!DEFAULTED("black-hole-kicks")) {                                                                                       // black hole kicks
            std::tie(found, m_BlackHoleKicks.type) = utils::GetMapKey(m_BlackHoleKicks.typeString, BLACK_HOLE_KICKS_LABEL, m_BlackHoleKicks.type);
            COMPLAIN_IF(!found, "Unknown Black Hole Kicks Option");
        }

        if (!DEFAULTED("case-BB-stability-prescription")) {                                                                         //case BB/BC mass transfer stability prescription
            std::tie(found, m_CaseBBStabilityPrescription.type) = utils::GetMapKey(m_CaseBBStabilityPrescription.typeString, CASE_BB_STABILITY_PRESCRIPTION_LABEL, m_CaseBBStabilityPrescription.type);
            COMPLAIN_IF(!found, "Unknown Case BB/BC Mass Transfer Stability Prescription");
        }
           
        if (!DEFAULTED("chemically-homogeneous-evolution")) {                                                                       // Chemically Homogeneous Evolution
            std::tie(found, m_CheMode.type) = utils::GetMapKey(m_CheMode.typeString, CHE_MODE_LABEL, m_CheMode.type);
            COMPLAIN_IF(!found, "Unknown Chemically Homogeneous Evolution Option");
        }

        if (!DEFAULTED("common-envelope-lambda-prescription")) {                                                                    // common envelope lambda prescription
            std::tie(found, m_CommonEnvelopeLambdaPrescription.type) = utils::GetMapKey(m_CommonEnvelopeLambdaPrescription.typeString, CE_LAMBDA_PRESCRIPTION_LABEL, m_CommonEnvelopeLambdaPrescription.type);
            COMPLAIN_IF(!found, "Unknown CE Lambda Prescription");
        }

        if (!DEFAULTED("common-envelope-mass-accretion-prescription")) {                                                            // common envelope mass accretion prescription
            std::tie(found, m_CommonEnvelopeMassAccretionPrescription.type) = utils::GetMapKey(m_CommonEnvelopeMassAccretionPrescription.typeString, CE_ACCRETION_PRESCRIPTION_LABEL, m_CommonEnvelopeMassAccretionPrescription.type);
            COMPLAIN_IF(!found, "Unknown CE Mass Accretion Prescription");
        }

        if (!DEFAULTED("critical-mass-ratio-prescription")) {                                                                       // critical mass ratio prescription
            std::tie(found, m_QCritPrescription.type) = utils::GetMapKey(m_QCritPrescription.typeString, QCRIT_PRESCRIPTION_LABEL, m_QCritPrescription.type);
            COMPLAIN_IF(!found, "Unknown qCrit Prescription");
        }
            
        if (!DEFAULTED("envelope-state-prescription")) {                                                                            // envelope state prescription
            std::tie(found, m_EnvelopeStatePrescription.type) = utils::GetMapKey(m_EnvelopeStatePrescription.typeString, ENVELOPE_STATE_PRESCRIPTION_LABEL, m_EnvelopeStatePrescription.type);
            COMPLAIN_IF(!found, "Unknown Envelope State Prescription");
        }

        if (!DEFAULTED("eccentricity-distribution")) {                                                                              // eccentricity distribution
            std::tie(found, m_EccentricityDistribution.type) = utils::GetMapKey(m_EccentricityDistribution.typeString, ECCENTRICITY_DISTRIBUTION_LABEL, m_EccentricityDistribution.type);
            COMPLAIN_IF(!found, "Unknown Eccentricity Distribution");
        }

        if (!DEFAULTED("fryer-supernova-engine")) {                                                                                 // Fryer et al. 2012 supernova engine
            std::tie(found, m_FryerSupernovaEngine.type) = utils::GetMapKey(m_FryerSupernovaEngine.typeString, SN_ENGINE_LABEL, m_FryerSupernovaEngine.type);
            COMPLAIN_IF(!found, "Unknown Fryer et al. Supernova Engine");
        }

        if (!DEFAULTED("initial-mass-function")) {                                                                                  // initial mass function
            std::tie(found, m_InitialMassFunction.type) = utils::GetMapKey(m_InitialMassFunction.typeString, INITIAL_MASS_FUNCTION_LABEL, m_InitialMassFunction.type);
            COMPLAIN_IF(!found, "Unknown Initial Mass Function");
        }

        if (!DEFAULTED("kick-direction")) {                                                                                         // kick direction
            std::tie(found, m_KickDirectionDistribution.type) = utils::GetMapKey(m_KickDirectionDistribution.typeString, KICK_DIRECTION_DISTRIBUTION_LABEL, m_KickDirectionDistribution.type);
            COMPLAIN_IF(!found, "Unknown Kick Direction Distribution");
        }

        if (!DEFAULTED("kick-magnitude-distribution")) {                                                                            // kick magnitude
            std::tie(found, m_KickMagnitudeDistribution.type) = utils::GetMapKey(m_KickMagnitudeDistribution.typeString, KICK_MAGNITUDE_DISTRIBUTION_LABEL, m_KickMagnitudeDistribution.type);
            COMPLAIN_IF(!found, "Unknown Kick Magnitude Distribution");
        }

        if (!DEFAULTED("logfile-type")) {                                                                                           // logfile type
            std::tie(found, m_LogfileType.type) = utils::GetMapKey(m_LogfileType.typeString, LOGFILETYPELabel, m_LogfileType.type);
            COMPLAIN_IF(!found, "Unknown Logfile Type");
        }

        if (!DEFAULTED("luminous-blue-variable-prescription")) {                                                                    // LBV mass loss prescription
            std::tie(found, m_LuminousBlueVariablePrescription.type) = utils::GetMapKey(m_LuminousBlueVariablePrescription.typeString, LBV_PRESCRIPTION_LABEL, m_LuminousBlueVariablePrescription.type);
            COMPLAIN_IF(!found, "Unknown LBV Mass Loss Prescription");
        }

        if (!DEFAULTED("mass-loss-prescription")) {                                                                                 // mass loss prescription
            std::tie(found, m_MassLossPrescription.type) = utils::GetMapKey(m_MassLossPrescription.typeString, MASS_LOSS_PRESCRIPTION_LABEL, m_MassLossPrescription.type);
            COMPLAIN_IF(!found, "Unknown Mass Loss Prescription");
        }

        if (!DEFAULTED("mass-ratio-distribution")) {                                                                                // mass ratio distribution
            std::tie(found, m_MassRatioDistribution.type) = utils::GetMapKey(m_MassRatioDistribution.typeString, MASS_RATIO_DISTRIBUTION_LABEL, m_MassRatioDistribution.type);
            COMPLAIN_IF(!found, "Unknown Mass Ratio Distribution");
        }

        if (m_UseMassTransfer && !DEFAULTED("mass-transfer-accretion-efficiency-prescription")) {                                   // mass transfer accretion efficiency prescription
            std::tie(found, m_MassTransferAccretionEfficiencyPrescription.type) = utils::GetMapKey(m_MassTransferAccretionEfficiencyPrescription.typeString, MT_ACCRETION_EFFICIENCY_PRESCRIPTION_LABEL, m_MassTransferAccretionEfficiencyPrescription.type);
            COMPLAIN_IF(!found, "Unknown Mass Transfer Angular Momentum Loss Prescription");
        }

        if (m_UseMassTransfer && !DEFAULTED("mass-transfer-angular-momentum-loss-prescription")) {                                  // mass transfer angular momentum loss prescription
            std::tie(found, m_MassTransferAngularMomentumLossPrescription.type) = utils::GetMapKey(m_MassTransferAngularMomentumLossPrescription.typeString, MT_ANGULAR_MOMENTUM_LOSS_PRESCRIPTION_LABEL, m_MassTransferAngularMomentumLossPrescription.type);
            COMPLAIN_IF(!found, "Unknown Mass Transfer Angular Momentum Loss Prescription");
        }

        if (m_UseMassTransfer && !DEFAULTED("mass-transfer-rejuvenation-prescription")) {                                           // mass transfer rejuvenation prescription
            std::tie(found, m_MassTransferRejuvenationPrescription.type) = utils::GetMapKey(m_MassTransferRejuvenationPrescription.typeString, MT_REJUVENATION_PRESCRIPTION_LABEL, m_MassTransferRejuvenationPrescription.type);
            COMPLAIN_IF(!found, "Unknown Mass Transfer Rejuvenation Prescription");
        }

        if (m_UseMassTransfer && !DEFAULTED("mass-transfer-thermal-limit-accretor")) {                                              // mass transfer accretor thermal limit
            std::tie(found, m_MassTransferThermallyLimitedVariation.type) = utils::GetMapKey(m_MassTransferThermallyLimitedVariation.typeString, MT_THERMALLY_LIMITED_VARIATION_LABEL, m_MassTransferThermallyLimitedVariation.type);
            COMPLAIN_IF(!found, "Unknown Mass Transfer Accretor Thermal Limit");

            if (m_MassTransferThermallyLimitedVariation.type == MT_THERMALLY_LIMITED_VARIATION::C_FACTOR) {
                m_MassTransferCParameter = DEFAULTED("mass-transfer-thermal-limit-C") ? 10.0 : m_MassTransferCParameter;
            }

            if (m_MassTransferThermallyLimitedVariation.type == MT_THERMALLY_LIMITED_VARIATION::RADIUS_TO_ROCHELOBE) {
                m_MassTransferCParameter = DEFAULTED("mass-transfer-thermal-limit-C") ? 1.0 : m_MassTransferCParameter;
            }
        }

        if (!DEFAULTED("metallicity-distribution")) {                                                                               // metallicity distribution
            std::tie(found, m_MetallicityDistribution.type) = utils::GetMapKey(m_MetallicityDistribution.typeString, METALLICITY_DISTRIBUTION_LABEL, m_MetallicityDistribution.type);
            COMPLAIN_IF(!found, "Unknown Metallicity Distribution");
        }

        if (!DEFAULTED("mode")) {                                                                                                   // mode
            std::tie(found, m_EvolutionMode.type) = utils::GetMapKey(m_EvolutionMode.typeString, EVOLUTION_MODE_LABEL, m_EvolutionMode.type);
            COMPLAIN_IF(!found, "Unknown Mode");
        }

        if (!DEFAULTED("neutrino-mass-loss-BH-formation")) {                                                                        // neutrino mass loss assumption
            std::tie(found, m_NeutrinoMassLossAssumptionBH.type) = utils::GetMapKey(m_NeutrinoMassLossAssumptionBH.typeString, NEUTRINO_MASS_LOSS_PRESCRIPTION_LABEL, m_NeutrinoMassLossAssumptionBH.type);
            COMPLAIN_IF(!found, "Unknown Neutrino Mass Loss Assumption");
        }

        if (!DEFAULTED("neutron-star-equation-of-state")) {                                                                         // neutron star equation of state
            std::tie(found, m_NeutronStarEquationOfState.type) = utils::GetMapKey(m_NeutronStarEquationOfState.typeString, NS_EOSLabel, m_NeutronStarEquationOfState.type);
            COMPLAIN_IF(!found, "Unknown Neutron Star Equation of State");
        }

        if (!DEFAULTED("pulsar-birth-magnetic-field-distribution")) {                                                               // pulsar birth magnetic field distribution
            std::tie(found, m_PulsarBirthMagneticFieldDistribution.type) = utils::GetMapKey(m_PulsarBirthMagneticFieldDistribution.typeString, PULSAR_BIRTH_MAGNETIC_FIELD_DISTRIBUTION_LABEL, m_PulsarBirthMagneticFieldDistribution.type);
            COMPLAIN_IF(!found, "Unknown Pulsar Birth Magnetic Field Distribution");
        }

        if (!DEFAULTED("pulsar-birth-spin-period-distribution")) {                                                                  // pulsar birth spin period distribution
            std::tie(found, m_PulsarBirthSpinPeriodDistribution.type) = utils::GetMapKey(m_PulsarBirthSpinPeriodDistribution.typeString, PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION_LABEL, m_PulsarBirthSpinPeriodDistribution.type);
            COMPLAIN_IF(!found, "Unknown Pulsar Birth Spin Period Distribution");
        }

        if (!DEFAULTED("pulsational-pair-instability-prescription")) {                                                              // pulsational pair instability prescription
            std::tie(found, m_PulsationalPairInstabilityPrescription.type) = utils::GetMapKey(m_PulsationalPairInstabilityPrescription.typeString, PPI_PRESCRIPTION_LABEL, m_PulsationalPairInstabilityPrescription.type);
            COMPLAIN_IF(!found, "Unknown Pulsational Pair Instability Prescription");
        }

        if (!DEFAULTED("remnant-mass-prescription")) {                                                                              // remnant mass prescription
            std::tie(found, m_RemnantMassPrescription.type) = utils::GetMapKey(m_RemnantMassPrescription.typeString, REMNANT_MASS_PRESCRIPTION_LABEL, m_RemnantMassPrescription.type);
            COMPLAIN_IF(!found, "Unknown Remnant Mass Prescription");
        }

        if (!DEFAULTED("rotational-velocity-distribution")) {                                                                       // rotational velocity distribution
            std::tie(found, m_RotationalVelocityDistribution.type) = utils::GetMapKey(m_RotationalVelocityDistribution.typeString, ROTATIONAL_VELOCITY_DISTRIBUTION_LABEL, m_RotationalVelocityDistribution.type);
            COMPLAIN_IF(!found, "Unknown Rotational Velocity Distribution");
        }

        if (!DEFAULTED("semi-major-axis-distribution")) {                                                                           // semi-major axis distribution
            std::tie(found, m_SemiMajorAxisDistribution.type) = utils::GetMapKey(m_SemiMajorAxisDistribution.typeString, SEMI_MAJOR_AXIS_DISTRIBUTION_LABEL, m_SemiMajorAxisDistribution.type);
            COMPLAIN_IF(!found, "Unknown Semi-Major Axis Distribution");
        }

        if (!DEFAULTED("stellar-zeta-prescription")) {                                                                              // common envelope zeta prescription
            std::tie(found, m_StellarZetaPrescription.type) = utils::GetMapKey(m_StellarZetaPrescription.typeString, ZETA_PRESCRIPTION_LABEL, m_StellarZetaPrescription.type);
            COMPLAIN_IF(!found, "Unknown stellar Zeta Prescription");
        }

        // constraint/value/range checks - alphabetically (where possible)

        COMPLAIN_IF(m_CommonEnvelopeAlpha < 0.0, "CE alpha (--common-envelope-alpha) < 0");
        COMPLAIN_IF(m_CommonEnvelopeAlphaThermal < 0.0 || m_CommonEnvelopeAlphaThermal > 1.0, "CE alpha thermal (--common-envelope-alpha-thermal) must be between 0 and 1");
        COMPLAIN_IF(m_CommonEnvelopeLambdaMultiplier < 0.0, "CE lambda multiplie (--common-envelope-lambda-multiplier < 0");
        COMPLAIN_IF(m_CommonEnvelopeMassAccretionConstant < 0.0, "CE mass accretion constant (--common-envelope-mass-accretion-constant) < 0");
        COMPLAIN_IF(m_CommonEnvelopeMassAccretionMax < 0.0, "Maximum accreted mass (--common-envelope-mass-accretion-max) < 0");
        COMPLAIN_IF(m_CommonEnvelopeMassAccretionMin < 0.0, "Minimum accreted mass (--common-envelope-mass-accretion-min) < 0");

        COMPLAIN_IF(m_CoolWindMassLossMultiplier < 0.0, "Wind mass loss multiplier for cool stars (--cool-wind-mass-loss-multiplier) < 0.0");

        COMPLAIN_IF(m_DebugLevel < 0, "Debug level (--debug-level) < 0");

        COMPLAIN_IF(m_Eccentricity < 0.0 || m_Eccentricity > 1.0, "Eccentricity (--eccentricity) must be between 0 and 1");
        COMPLAIN_IF(m_EccentricityDistributionMin < 0.0 || m_EccentricityDistributionMin > 1.0, "Minimum eccentricity (--eccentricity-min) must be between 0 and 1");
        COMPLAIN_IF(m_EccentricityDistributionMax < 0.0 || m_EccentricityDistributionMax > 1.0, "Maximum eccentricity (--eccentricity-max) must be between 0 and 1");
        COMPLAIN_IF(m_EccentricityDistributionMax <= m_EccentricityDistributionMin, "Maximum eccentricity (--eccentricity-max) must be > Minimum eccentricity (--eccentricity-min)");

        COMPLAIN_IF(m_GridStartLine < 0, "Grid file start line (--grid-start-line) < 0");
        COMPLAIN_IF(!DEFAULTED("grid-lines-to-process") && m_GridLinesToProcess < 1, "Grid file lines to process (--grid-lines-to-process) < 1");

        COMPLAIN_IF(m_HDF5BufferSize < 1, "HDF5 IO buffer size (--hdf5-buffer-size) must be >= 1");
        COMPLAIN_IF(m_HDF5ChunkSize < HDF5_MINIMUM_CHUNK_SIZE, "HDF5 file dataset chunk size (--hdf5-chunk-size) must be >= minimum chunk size of " + std::to_string(HDF5_MINIMUM_CHUNK_SIZE));

        COMPLAIN_IF(m_InitialMass < MINIMUM_INITIAL_MASS || m_InitialMass > MAXIMUM_INITIAL_MASS, "Initial mass (--initial-mass) must be between " + std::to_string(MINIMUM_INITIAL_MASS) + " and " + std::to_string(MAXIMUM_INITIAL_MASS) + " Msol");
        COMPLAIN_IF(m_InitialMass1 < MINIMUM_INITIAL_MASS || m_InitialMass1 > MAXIMUM_INITIAL_MASS, "Primary initial mass (--initial-mass-1) must be between " + std::to_string(MINIMUM_INITIAL_MASS) + " and " + std::to_string(MAXIMUM_INITIAL_MASS) + " Msol");
        COMPLAIN_IF(m_InitialMass2 < MINIMUM_INITIAL_MASS || m_InitialMass2 > MAXIMUM_INITIAL_MASS, "Secondary initial mass (--initial-mass-2) must be between " + std::to_string(MINIMUM_INITIAL_MASS) + " and " + std::to_string(MAXIMUM_INITIAL_MASS) + " Msol");

        COMPLAIN_IF(m_InitialMassFunctionMin < MINIMUM_INITIAL_MASS, "Minimum initial mass (--initial-mass-min) must be >= " + std::to_string(MINIMUM_INITIAL_MASS) + " Msol");
        COMPLAIN_IF(m_InitialMassFunctionMax > MAXIMUM_INITIAL_MASS, "Maximum initial mass (--initial-mass-max) must be <= " + std::to_string(MAXIMUM_INITIAL_MASS) + " Msol");
        COMPLAIN_IF(m_InitialMassFunctionMax <= m_InitialMassFunctionMin, "Maximum initial mass (--initial-mass-max) must be > Minimum initial mass (--initial-mass-min)");

        if (m_KickMagnitudeDistribution.type == KICK_MAGNITUDE_DISTRIBUTION::FLAT) {
            COMPLAIN_IF(m_KickMagnitudeDistributionMaximum <= 0.0, "User specified --kick-magnitude-distribution = FLAT with Maximum kick magnitude (--kick-magnitude-max) <= 0.0");
        }

        COMPLAIN_IF(m_LogLevel < 0, "Logging level (--log-level) < 0");
 
        COMPLAIN_IF(m_LuminousBlueVariableFactor < 0.0, "LBV multiplier (--luminous-blue-variable-multiplier) < 0");

        COMPLAIN_IF(m_MassRatio <= 0.0 || m_MassRatio > 1.0, "Mass ratio (--mass-ratio) must be greater than 0 and less than or equal to 1");

        COMPLAIN_IF(m_MassRatioDistributionMin <= 0.0 || m_MassRatioDistributionMin > 1.0, "Minimum mass ratio (--mass-ratio-min) must be greater than 0 and less than or equal to 1");
        COMPLAIN_IF(m_MassRatioDistributionMax <= 0.0 || m_MassRatioDistributionMax > 1.0, "Maximum mass ratio (--mass-ratio-max) must be greater than 0 and less than or equal to 1");
        COMPLAIN_IF(m_MassRatioDistributionMax <= m_MassRatioDistributionMin, "Maximum mass ratio (--mass-ratio-max) must be > Minimum mass ratio (--mass-ratio-min)");

        COMPLAIN_IF(m_MaxEvolutionTime <= 0.0, "Maximum evolution time in Myr (--maxEvolutionTime) must be > 0");

        COMPLAIN_IF(m_Metallicity < MINIMUM_METALLICITY || m_Metallicity > MAXIMUM_METALLICITY, "Metallicity (--metallicity) should be absolute metallicity and must be between " + std::to_string(MINIMUM_METALLICITY) + " and " + std::to_string(MAXIMUM_METALLICITY));
        COMPLAIN_IF(m_MetallicityDistributionMin < MINIMUM_METALLICITY || m_MetallicityDistributionMin > MAXIMUM_METALLICITY, "Minimum metallicity (--metallicity-min) must be between " + std::to_string(MINIMUM_METALLICITY) + " and " + std::to_string(MAXIMUM_METALLICITY));
        COMPLAIN_IF(m_MetallicityDistributionMax < MINIMUM_METALLICITY || m_MetallicityDistributionMax > MAXIMUM_METALLICITY, "Maximum metallicity (--metallicity-max) must be between " + std::to_string(MINIMUM_METALLICITY) + " and " + std::to_string(MAXIMUM_METALLICITY));
        COMPLAIN_IF(m_MetallicityDistributionMax <= m_MetallicityDistributionMin, "Maximum metallicity (--metallicity-max) must be > Minimum metallicity (--metallicity-min)");

        COMPLAIN_IF(m_MinimumMassSecondary < MINIMUM_INITIAL_MASS, "Seconday minimum mass (--minimum-secondary-mass) must be >= minimum initial mass of " + std::to_string(MINIMUM_INITIAL_MASS) + " Msol");
        COMPLAIN_IF(m_MinimumMassSecondary > MAXIMUM_INITIAL_MASS, "Seconday minimum mass (--minimum-secondary-mass) must be <= maximum initial mass of " + std::to_string(MAXIMUM_INITIAL_MASS) + " Msol");

        if (m_NeutrinoMassLossAssumptionBH.type == NEUTRINO_MASS_LOSS_PRESCRIPTION::FIXED_MASS) {
            COMPLAIN_IF(m_NeutrinoMassLossValueBH < 0.0, "Neutrino mass loss value < 0");
        }

        COMPLAIN_IF(m_ObjectsToEvolve <= 0, (m_EvolutionMode.type == EVOLUTION_MODE::SSE ? "Number of stars requested <= 0" : "Number of binaries requested <= 0"));
    
        if (m_NeutrinoMassLossAssumptionBH.type == NEUTRINO_MASS_LOSS_PRESCRIPTION::FIXED_FRACTION) {
            COMPLAIN_IF(m_NeutrinoMassLossValueBH < 0.0 || m_NeutrinoMassLossValueBH > 1.0, "Neutrino mass loss must be between 0 and 1");
        }

        if (!DEFAULTED("notes")) {                                                                                                  // user specified notes?
            WARNUSER_IF(m_Notes.size() > Options::Instance()->NotesHdrs().size(), "WARNING: Annotations: more notes than headers - extra notes ignored"); // yes - check counts
        }

        if (!DEFAULTED("output-path")) {                                                                                            // user specified output path?
                                                                                                                                    // yes
            fs::path userPath = m_OutputPathString;                                                                                 // user-specifed path
            if (fs::is_directory(userPath)) {                                                                                       // valid directory?
                m_OutputPath = userPath;                                                                                            // yes - set outputPath to user-specified path
            }
            else {                                                                                                                  // not a valid directory
                m_OutputPath = m_DefaultOutputPath;                                                                                 // use default path = CWD
            }
        }

        COMPLAIN_IF(m_OrbitalPeriodDistributionMin < 0.0, "Minimum orbital period (--orbital-period-min) < 0");
        COMPLAIN_IF(m_OrbitalPeriodDistributionMax < 0.0, "Maximum orbital period (--orbital-period-max) < 0");
        COMPLAIN_IF(m_OrbitalPeriodDistributionMax <= m_OrbitalPeriodDistributionMin, "Maximum orbital period (--orbital-period-max) must be > Minimum orbital period (--orbital-period-min)");

        COMPLAIN_IF(m_OverallWindMassLossMultiplier < 0.0, "Overall wind mass loss multiplier (--overall-wind-mass-loss-multiplier) < 0.0");

        COMPLAIN_IF(!DEFAULTED("pulsar-magnetic-field-decay-timescale") && m_PulsarMagneticFieldDecayTimescale <= 0.0, "Pulsar magnetic field decay timescale (--pulsar-magnetic-field-decay-timescale) <= 0");
        COMPLAIN_IF(!DEFAULTED("pulsar-magnetic-field-decay-massscale") && m_PulsarMagneticFieldDecayMassscale <= 0.0, "Pulsar Magnetic field decay massscale (--pulsar-magnetic-field-decay-massscale) <= 0");

        COMPLAIN_IF(!DEFAULTED("rotational-frequency")  && m_RotationalFrequency < 0.0, "Rotational frequency (--rotational-frequency) < 0");
        COMPLAIN_IF(!DEFAULTED("rotational-frequency-1") && m_RotationalFrequency1 < 0.0, "Primary rotational frequency (--rotational-frequency-1) < 0");
        COMPLAIN_IF(!DEFAULTED("rotational-frequency-2") && m_RotationalFrequency2 < 0.0, "Secondary rotational frequency (--rotational-frequency-2) < 0");

        COMPLAIN_IF(m_SemiMajorAxisDistributionMin < 0.0, "Minimum semi-major Axis (--semi-major-axis-min) < 0");
        COMPLAIN_IF(m_SemiMajorAxisDistributionMax < 0.0, "Maximum semi-major Axis (--semi-major-axis-max) < 0");

        COMPLAIN_IF(m_TimestepMultiplier <= 0.0, "Timestep multiplier (--timestep-multiplier) <= 0");

        COMPLAIN_IF(m_WolfRayetFactor < 0.0, "WR multiplier (--wolf-rayet-multiplier) < 0");

        COMPLAIN_IF(!DEFAULTED("initial-mass")   && m_InitialMass  <= 0.0, "Initial mass (--initial-mass) <= 0");                   // initial mass must be > 0.0
        COMPLAIN_IF(!DEFAULTED("initial-mass-1") && m_InitialMass1 <= 0.0, "Primary initial mass (--initial-mass-1) <= 0");         // primary initial mass must be > 0.0
        COMPLAIN_IF(!DEFAULTED("initial-mass-2") && m_InitialMass2 <= 0.0, "Secondary initial mass (--initial-mass-2) <= 0");       // secondary initial mass must be > 0.0

        COMPLAIN_IF(!DEFAULTED("semi-major-axis") && m_SemiMajorAxis <= 0.0, "Semi-major axis (--semi-major-axis) <= 0");           // semi-major axis must be > 0.0
        COMPLAIN_IF(!DEFAULTED("orbital-period")  && m_OrbitalPeriod <= 0.0, "Orbital period (--orbital-period) <= 0");             // orbital period must be > 0.0

        COMPLAIN_IF(m_KickMagnitude  < 0.0, "Kick magnitude (--kick-magnitude) must be >= 0");
        COMPLAIN_IF(m_KickMagnitude1 < 0.0, "Kick magnitude (--kick-magnitude-1) must be >= 0");
        COMPLAIN_IF(m_KickMagnitude2 < 0.0, "Kick magnitude (--kick-magnitude-2) must be >= 0");

        COMPLAIN_IF(m_KickMagnitudeRandom  < 0.0 || m_KickMagnitudeRandom  >= 1.0, "Kick magnitude random (--kick-magnitude-random) must be >= 0 and < 1");
        COMPLAIN_IF(m_KickMagnitudeRandom1 < 0.0 || m_KickMagnitudeRandom1 >= 1.0, "Kick magnitude random (--kick-magnitude-random-1) must be >= 0 and < 1");
        COMPLAIN_IF(m_KickMagnitudeRandom2 < 0.0 || m_KickMagnitudeRandom2 >= 1.0, "Kick magnitude random (--kick-magnitude-random-2) must be >= 0 and < 1");

        errStr = SetCalculatedOptionDefaults(BOOST_MAP::NO_UPDATE);                                                                 // set calculated option values
    }
    catch (po::error& e) {                                                                                                          // program options exception
        errStr = e.what();                                                                                                          // set the error string
    }
    catch (const std::string eStr) {                                                                                                // custom exception
        errStr = eStr;                                                                                                              // set the error string
    }
    catch (...) {                                                                                                                   // unhandled exception
        errStr = ERR_MSG(ERROR::UNHANDLED_EXCEPTION);                                                                               // set the error string
    }

    return errStr;
#undef DEFAULTED
}


/*
 * Determine if the user specified a value for the option
 *
 * Note that this function does not check whether the option string
 * pass as p_OptionString is a valid option string - it just checks
 * whether the user specfied it, either at the grid line level, or
 * at the commandline level.
 * 
 * 
 * int OptionSpecified(std::string p_OptionString) 
 * 
 * 
 * @param   [IN]    p_OptionString              String containing option name
 * @return                                      Int result:
 *                                                   0: option was not specified by user
 *                                                   1: option specified by user
 */
int Options::OptionSpecified(const std::string p_OptionString) {

    std::string opt = p_OptionString;
    opt = utils::ToLower(utils::trim(opt));

    // check if option specified at grid line (evolving object) level?
    auto gridLineIt = std::find_if(
        m_GridLine.optionsSpecified.begin(), m_GridLine.optionsSpecified.end(), [&opt](const OPTIONSTR& e) {
            return std::get<0>(e) == opt || std::get<1>(e) == opt || std::get<2>(e) == opt || std::get<3>(e) == opt;
        }
    );
    if (gridLineIt != m_GridLine.optionsSpecified.end()) return 1;

    // check if option specified at commandline (program) level?
    auto cmdLineIt = std::find_if(
        m_CmdLine.optionsSpecified.begin(), m_CmdLine.optionsSpecified.end(), [&opt](const OPTIONSTR& e) {
            return std::get<0>(e) == opt || std::get<1>(e) == opt || std::get<2>(e) == opt || std::get<3>(e) == opt;
        }
    );
    if (cmdLineIt != m_CmdLine.optionsSpecified.end()) return 1;

    return 0;       // not specified
}


/*
 * Retrieve the attributes of an option
 *
 * The option for which the attributes are to be retreived is passed as an iterator pointing at the option in the boost
 * variables map.  Note that this function is private to the Options class and is intended for Options internal use
 * only.  External actors should use the public function Options::OptionValue() to get option values.
 * 
 * The attributes are returned as a tuple, described by typedef ATTR, containing
 * 
 *     - dataType       TYPENAME (high-level) data type of the attribute.  Will be one of {NONE, BOOL, INT, FLOAT, STRING}
 *     - defaulted      BOOL     flag to indicate if the option was specified by the user or defaulted to the defaulty value
 *     - typeStr        STRING   detailed data type returned as a string (e.g. "UNSIGNED LONG INT" etc.)
 *     - valueStr       STRING   the value of the option returned as a string (e.g. "2.3", "BSE" etc.)
 * 
 * 
 * Options::ATTR OptionAttributes(const po::variables_map p_VM, const po::variables_map::const_iterator p_IT)
 * 
 *
 * @param   [IN]    p_VM                        The boost variables map
 * @param   [IN]    p_IT                        Iterator for the boost variables map pointing to the option
 *                                              for which the attributed are required to be retrieved
 * @return                                      Tuple (type ATTR) containing the option attributes
 */
Options::ATTR Options::OptionAttributes(const po::variables_map p_VM, const po::variables_map::const_iterator p_IT) {
            
    TYPENAME    dataType  = TYPENAME::NONE;
    std::string typeStr   = "";
    bool        defaulted = false;
    std::string valueStr  = "";

    if (((boost::any)p_IT->second.value()).empty()) return std::make_tuple(TYPENAME::NONE, true, "", "");   // empty option 

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
    }
    catch (const boost::bad_any_cast &) {
        isCharPtr = false;
    }

    if (!isCharPtr) {
        // (pre)check for data type = string
        try {
            boost::any_cast<std::string>(p_IT->second.value());
            isStr = true;
        }
        catch (const boost::bad_any_cast &) {
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
        if (tmp.length()) valueStr = "'" + tmp + "'";
        else              valueStr = "''";
    }

    else if (((boost::any)p_IT->second.value()).type() == typeid(signed                )) { dataType = TYPENAME::INT;          typeStr = "SIGNED";                 valueStr = std::to_string(p_VM[p_IT->first].as<signed                >()); }
    else if (((boost::any)p_IT->second.value()).type() == typeid(unsigned              )) { dataType = TYPENAME::INT;          typeStr = "UNSIGNED";               valueStr = std::to_string(p_VM[p_IT->first].as<unsigned              >()); }

    else if (((boost::any)p_IT->second.value()).type() == typeid(short                 )) { dataType = TYPENAME::INT;          typeStr = "SHORT";                  valueStr = std::to_string(p_VM[p_IT->first].as<short                 >()); }
    else if (((boost::any)p_IT->second.value()).type() == typeid(signed short          )) { dataType = TYPENAME::INT;          typeStr = "SIGNED_SHORT";           valueStr = std::to_string(p_VM[p_IT->first].as<signed short          >()); }
    else if (((boost::any)p_IT->second.value()).type() == typeid(unsigned short        )) { dataType = TYPENAME::INT;          typeStr = "UNSIGNED_SHORT";         valueStr = std::to_string(p_VM[p_IT->first].as<unsigned short        >()); }

    else if (((boost::any)p_IT->second.value()).type() == typeid(short int             )) { dataType = TYPENAME::INT;          typeStr = "SHORT_INT";              valueStr = std::to_string(p_VM[p_IT->first].as<short int             >()); }
    else if (((boost::any)p_IT->second.value()).type() == typeid(signed short int      )) { dataType = TYPENAME::INT;          typeStr = "SIGNED_SHORT_INT";       valueStr = std::to_string(p_VM[p_IT->first].as<signed short int      >()); }
    else if (((boost::any)p_IT->second.value()).type() == typeid(unsigned short int    )) { dataType = TYPENAME::INT;          typeStr = "UNSIGNED_SHORT_INT";     valueStr = std::to_string(p_VM[p_IT->first].as<unsigned short int    >()); }

    else if (((boost::any)p_IT->second.value()).type() == typeid(int                   )) { dataType = TYPENAME::INT;          typeStr = "INT";                    valueStr = std::to_string(p_VM[p_IT->first].as<int                   >()); }
    else if (((boost::any)p_IT->second.value()).type() == typeid(signed int            )) { dataType = TYPENAME::INT;          typeStr = "SIGNED_INT";             valueStr = std::to_string(p_VM[p_IT->first].as<signed int            >()); }
    else if (((boost::any)p_IT->second.value()).type() == typeid(unsigned int          )) { dataType = TYPENAME::INT;          typeStr = "UNSIGNED_INT";           valueStr = std::to_string(p_VM[p_IT->first].as<unsigned int          >()); }

    else if (((boost::any)p_IT->second.value()).type() == typeid(long                  )) { dataType = TYPENAME::LONGINT;      typeStr = "LONG";                   valueStr = std::to_string(p_VM[p_IT->first].as<long                  >()); }
    else if (((boost::any)p_IT->second.value()).type() == typeid(signed long           )) { dataType = TYPENAME::LONGINT;      typeStr = "SIGNED_LONG";            valueStr = std::to_string(p_VM[p_IT->first].as<signed long           >()); }
    else if (((boost::any)p_IT->second.value()).type() == typeid(unsigned long         )) { dataType = TYPENAME::ULONGINT;     typeStr = "UNSIGNED_LONG";          valueStr = std::to_string(p_VM[p_IT->first].as<unsigned long         >()); }

    else if (((boost::any)p_IT->second.value()).type() == typeid(long int              )) { dataType = TYPENAME::LONGINT;      typeStr = "LONG_INT";               valueStr = std::to_string(p_VM[p_IT->first].as<long int              >()); }
    else if (((boost::any)p_IT->second.value()).type() == typeid(signed long int       )) { dataType = TYPENAME::LONGINT;      typeStr = "SIGNED_LONG_INT";        valueStr = std::to_string(p_VM[p_IT->first].as<signed long int       >()); }
    else if (((boost::any)p_IT->second.value()).type() == typeid(unsigned long int     )) { dataType = TYPENAME::ULONGINT;     typeStr = "UNSIGNED_LONG_INT";      valueStr = std::to_string(p_VM[p_IT->first].as<unsigned long int     >()); }

    else if (((boost::any)p_IT->second.value()).type() == typeid(long long             )) { dataType = TYPENAME::LONGLONGINT;  typeStr = "LONG_LONG";              valueStr = std::to_string(p_VM[p_IT->first].as<long long             >()); }
    else if (((boost::any)p_IT->second.value()).type() == typeid(signed long long      )) { dataType = TYPENAME::LONGLONGINT;  typeStr = "SIGNED_LONG_LONG";       valueStr = std::to_string(p_VM[p_IT->first].as<signed long long      >()); }
    else if (((boost::any)p_IT->second.value()).type() == typeid(unsigned long long    )) { dataType = TYPENAME::LONGLONGINT;  typeStr = "UNSIGNED_LONG_LONG";     valueStr = std::to_string(p_VM[p_IT->first].as<unsigned long long    >()); }

    else if (((boost::any)p_IT->second.value()).type() == typeid(long long int         )) { dataType = TYPENAME::LONGLONGINT;  typeStr = "LONG_LONG_INT";          valueStr = std::to_string(p_VM[p_IT->first].as<long long int         >()); }
    else if (((boost::any)p_IT->second.value()).type() == typeid(signed long long int  )) { dataType = TYPENAME::LONGLONGINT;  typeStr = "SIGNED_LONG_LONG_INT";   valueStr = std::to_string(p_VM[p_IT->first].as<signed long long int  >()); }
    else if (((boost::any)p_IT->second.value()).type() == typeid(unsigned long long int)) { dataType = TYPENAME::ULONGLONGINT; typeStr = "UNSIGNED_LONG_LONG_INT"; valueStr = std::to_string(p_VM[p_IT->first].as<unsigned long long int>()); }

    else if (((boost::any)p_IT->second.value()).type() == typeid(float                 )) { dataType = TYPENAME::FLOAT;        typeStr = "FLOAT";                  valueStr = std::to_string(p_VM[p_IT->first].as<float                 >()); }
    else if (((boost::any)p_IT->second.value()).type() == typeid(double                )) { dataType = TYPENAME::DOUBLE;       typeStr = "DOUBLE";                 valueStr = std::to_string(p_VM[p_IT->first].as<double                >()); }
    else if (((boost::any)p_IT->second.value()).type() == typeid(long double           )) { dataType = TYPENAME::LONGDOUBLE;   typeStr = "LONG_DOUBLE";            valueStr = std::to_string(p_VM[p_IT->first].as<long double           >()); }

    else if (((boost::any)p_IT->second.value()).type() == typeid(char                  )) { dataType = TYPENAME::INT;          typeStr = "CHAR";                   valueStr = std::to_string(p_VM[p_IT->first].as<char                  >()); }
    else if (((boost::any)p_IT->second.value()).type() == typeid(signed char           )) { dataType = TYPENAME::INT;          typeStr = "SIGNED_CHAR";            valueStr = std::to_string(p_VM[p_IT->first].as<signed char           >()); }
    else if (((boost::any)p_IT->second.value()).type() == typeid(unsigned char         )) { dataType = TYPENAME::INT;          typeStr = "UNSIGNED_CHAR";          valueStr = std::to_string(p_VM[p_IT->first].as<unsigned char         >()); }

    else if (((boost::any)p_IT->second.value()).type() == typeid(bool)) {
        dataType = TYPENAME::BOOL;
        typeStr  = "BOOL";
        valueStr = p_VM[p_IT->first].as<bool>() ? "TRUE" : "FALSE";
    } 

    else {  // Assume std::vector<std::string>
        try {
            std::ostringstream elemsSS;
            elemsSS << "{ ";
            std::vector<std::string> tmp = p_VM[p_IT->first].as<std::vector<std::string>>();
            for (std::vector<std::string>::iterator elem=tmp.begin(); elem != tmp.end(); elem++) {
                elemsSS << "'" << (*elem) << "', ";
            }
            std::string elems = elemsSS.str();
            if (elems.length() > 2) elems.erase(elems.length() - 2);
            else if (elems.length() == 2) elems.erase(elems.length() - 1);
            elems += " }";

            // the following options are declared as std::vector<std::string>>:
            //
            //     debug-classes
            //     log-classes
            //     notes
            //     notes-hdrs
            // 
            // The vector of strings is just formatted as a string here - with braces
            // sourrounding comma-separated values.
            //
            // We return dateType = TYPENAME::STRING, but typeStr = "VECTOR<STRING>"

            dataType = TYPENAME::STRING;  
            typeStr  = "VECTOR<STRING>";
            valueStr = elems;
        }
        catch (const boost::bad_any_cast &) {
            dataType = TYPENAME::NONE;                                                  // unknown data type               
            typeStr  = "<UNKNOWN_DATA_TYPE>";
            valueStr = "<UNKNOWN_DATA_TYPE>";
        }
    }

    return std::make_tuple(dataType, defaulted, typeStr, valueStr);
}


/*
 * Get option details for the Run_Details file
 *
 * The parameter passed is the options descriptor - the grid line options, or
 * the commandline options.  Ordinarily we would build the Run_Details contents
 * from the commandline options, but the flexibility exists to use a set of
 * grid line options (maybe one day we will want to (optionally) produce a
 * per star/binary Run_Details file)
 * 
 * 
 * std::vector<std::tuple<std::string, std::string, std::string, std::string, TYPENAME>> Options::OptionDetails(const OptionsDescriptorT &p_Options)
 * 
 * @param   [IN]    p_Options                   The options descriptor to use to build the output string
 * @return                                      Vector containing the option details for Run_Details
 */
std::vector<std::tuple<std::string, std::string, std::string, std::string, TYPENAME>> Options::OptionDetails(const OptionsDescriptorT &p_Options) {

    std::vector<std::tuple<std::string, std::string, std::string, std::string, TYPENAME>> optionDetails = {};

    TYPENAME    dataType  = TYPENAME::NONE;
    std::string typeStr   = "";
    bool        defaulted = false;
    std::string valueStr  = "";

    for (po::variables_map::const_iterator it = p_Options.optionValues.m_VM.begin(); it != p_Options.optionValues.m_VM.end(); it++) {                                   // for all options in the variable map

        std::tie(dataType, defaulted, typeStr, valueStr) = OptionAttributes(p_Options.optionValues.m_VM, it);                                                           // get option attributes

        if (valueStr == "")                                                                                                                                             // empty option?
            optionDetails.push_back(std::make_tuple(it->first, "<EMPTY_OPTION>", "<EMPTY_OPTION>", "<EMPTY_OPTION>", TYPENAME::NONE));                                  // yes - say so
        else                                                                                                                                                            // no
            optionDetails.push_back(std::make_tuple(it->first, valueStr, (defaulted ? "DEFAULT_USED" : "USER_SUPPLIED"), typeStr, dataType));                           // add option details to return vector
    }
  
    // add other (calculated) options

    optionDetails.push_back(std::make_tuple("useFixedUK", (p_Options.optionValues.m_UseFixedUK ? "TRUE" : "FALSE"), "CALCULATED", "BOOL", TYPENAME::BOOL));             // useFixedUK
    optionDetails.push_back(std::make_tuple("actual-output-path", p_Options.optionValues.m_OutputPath.string(), "CALCULATED", "STRING", TYPENAME::STRING));             // output-path
    optionDetails.push_back(std::make_tuple("fixedRandomSeed", (p_Options.optionValues.m_FixedRandomSeed ? "TRUE" : "FALSE"), "CALCULATED", "BOOL", TYPENAME::BOOL));   // fixedRandomSeed

    return optionDetails;
}


/*
 * Show available options
 * 
 * If the p_Verbose parameter is false, just print option names
 * If the p_Verbose parameter is true, print option names and descriptions
 * 
 * 
 * void PrintOptionHelp(const bool p_Verbose)
 * 
 * @param   [IN]    p_Verbose                   Boolean to indicate whether the option descriptions should
 *                                              be printed (true), or just option names (false)
 */
void Options::PrintOptionHelp(const bool p_Verbose) {

    std::cout << "Options:" << std::endl;

    for (po::variables_map::const_iterator it = m_CmdLine.optionValues.m_VM.begin(); it != m_CmdLine.optionValues.m_VM.end(); it++) {
  
        po::option_description const& opt = m_CmdLine.optionDescriptions.find(it->first, false, false, false); 

        std::string optionLongName  = opt.canonical_display_name(cls::allow_long);                          // long name ('--') prefix
        if (optionLongName[0] == '-') optionLongName.erase(0, optionLongName.find_first_not_of("-"));       // remove the "-" or "--"

        std::string optionShortName = opt.canonical_display_name(cls::allow_dash_for_short);                // short name ('-') prefix
        if (optionShortName[0] == '-') optionShortName.erase(0, optionShortName.find_first_not_of("-"));    // remove the "-" or "--"

        std::cout << "--" << optionLongName;
        if (optionLongName != optionShortName) std::cout << " [ -" << optionShortName << " ]";
        std::cout << std::endl;

        if (p_Verbose) {
            std::cout << "  " << opt.description() << std::endl;
        }
    }
}


/*
 * Returns TRUE if parameter p_TypeName is a supported COMPAS numeric datatype
 * for program options, otherwise FALSE
 *
 * The datatypes here should cover our options for now - but we might have to 
 * refine them over time
 * 
 * 
 * bool IsSupportedNumericDataType(TYPENAME p_TypeName)
 * 
 * @param   [IN]    p_TypeName                  COMPAS datatype name
 * @return                                      True if p_TypeName is a supported numeric datatype, else false
 */
bool Options::IsSupportedNumericDataType(TYPENAME p_TypeName) {

    bool supported = false;

    switch(p_TypeName) {
        case TYPENAME::INT:
        case TYPENAME::LONGINT:
        case TYPENAME::ULONGINT:
        case TYPENAME::FLOAT:
        case TYPENAME::DOUBLE:
        case TYPENAME::LONGDOUBLE:
            supported = true;
            break;
        default:
            supported = false;
    }
    return supported;
}


/*
 * Preprocess the options provided by the user - expand any shorthand notation
 * 
 * Before we parse the options provided by the user, we replace/expand any shorthand devices the user has taken
 * advantage of, so we can present the expanded form of the options to boost.
 * 
 * For example, we provide shorthand for users to specify annotations and annotation headers.  Both of these options
 * are defined as boost vector options, and would typically be specified by the user thus:
 * 
 * ./compas --notes-hdrs hdrStr1 hdrStr2 hdrStr3 --notes "note 1" "another note" "this is note 3" --option-name option-value ...
 * 
 * We allow blank notes, but they must be entered as empty strings using this method.  e.g.:
 * 
 * ./compas --notes-hdrs hdrStr1 hdrStr2 hdrStr3 --notes "note 1" "" "this is note 3" --option-name option-value ... (note 2 is blank)
 * 
 * That format could become awkward, so we provide a shorthand method for specifying vector options.  The shorthand method allows
 * users to list the comma-separated option values enclosed in square brackets "[...]", and any blank values can just be omitted.
 * 
 * e.g., the second example above could be specified as:
 * 
 * ./compas --notes-hdrs [hdrStr1,hdrStr2,hdrStr3] --notes ["note 1",,"this is note 3"] --option-name option-value ... (note 2 is omitted)
 * ./compas --notes-hdrs [hdrStr1,hdrStr2,hdrStr3] --notes ["note 1",,] --option-name option-value ... (note 2 and note 3 are omitted)
 * 
 * This function will expand this shorthand to the example shown above.  There is no checking for correctness here - we just expand any
 * shorthand necessary and pass the argument vector back - correctness checking is done elsewhere.  
 * 
 * The value for any omitted option values will be a string of length 1, with the char value NOT_PROVIDED (constant define in Options.h).
 * Since at this stage the option names and values are just strings that will be parsed by boost, we don't need to worry about data type - 
 * code processing the options can check for NOT_PROVIDED and deal with it then.  Since we may not know the maxmimum number of values
 * expected (e.g. the number of notes-hdrs specifies the maximum number of notes expected, and we may not have that number yet), we leave
 * it to later to pad out missing values beyond the last one specified here.  e.g. a specification shuch as:
 * 
 * ./compas --notes-hdrs [hdrStr1,hdrStr2,hdrStr3,hdrStr4,hdrStr5] --notes ["note 1",,"this is note 3"] --option-name option-value ...
 * 
 * has notes 2, 4 & 5 omitted - we specify note 2 as NOT_PROVIDED here, and notes 4 & 5 will be added later.
 * 
 * 
 * std::tuple<std::string, int, std::vector<std::string>> ExpandShorthandOptionValues(int p_ArgCount, char *p_ArgStrings[])
 * 
 * 
 * @param   [IN]    p_ArgCount                  The number of argument strings. (note below for p_ArgStrings)
 * @param   [IN]    p_ArgStrings                The argument strings.  The first argument string is expected
 *                                              (by boost) to be the executable name (boost expects the arguments
 *                                              to be command-line arguments passed to main())
 * @return                                      Tuple containing:
 *                                                  String containing an error string (if no error occurred the return string will be the empty string)
 *                                                  The expanded arguments:
 *                                                      - an integer indicating the number of arguments (analogous to p_ArgCount)
 *                                                      - a vector of strings containing the arguments (analogous to p_ArgStrings)
 */
std::tuple<std::string, int, std::vector<std::string>> Options::ExpandShorthandOptionValues(int p_ArgCount, char *p_ArgStrings[]) {

    std::string errStr = "";                                                                                                // for now

    std::vector<std::string> strargs = {std::string(p_ArgStrings[0])};                                                      // new args vector - command name is arg[0]

    std::string  argString  = "";                                                                                           // argument string
    std::string  optionName = "";                                                                                           // option name
    for (size_t iArg = 1; iArg < (size_t)p_ArgCount; iArg++) {                                                              // for each arg string

        argString = p_ArgStrings[iArg];                                                                                     // the argument we're processing

        if (argString[0] == '-') {                                                                                          // is it actually an option name?
            strargs.push_back(argString);                                                                                   // yes - add it to the new vector (unadulterated)

            optionName = argString;                                                                                         // get the option name for the argument we'll be processing
            optionName = utils::ToLower(utils::trim(optionName));                                                           // downshift and trim whitespace
            optionName.erase(0, optionName.find_first_not_of("-"));                                                         // remove the "-" or "--"
        }
        else {                                                                                                              // argument is not an option name - process it
            if (iArg == 1) {                                                                                                // first arg (i.e. no option name)?
                strargs.push_back(argString);                                                                               // yes - add it to the new vector (unadulterated)
            }
            else {                                                                                                          // no - not first arg - process it
                // check if option is on the shorthand allowed list, and process it if it is
                // presumably invalid options won't be on the list...
                auto it = std::find_if(m_ShorthandAllowed.begin(), m_ShorthandAllowed.end(), [&optionName](const SHORTHAND_ENTRY& e) { return std::get<0>(e) == optionName; });
                if (it == m_ShorthandAllowed.end()) {                                                                       // option in shorthand allowed list?
                    strargs.push_back(argString);                                                                           // no - add option value to the new vector
                }
                else {                                                                                                      // yes - shorthand allowed - process it
                    if (!argString.empty()) {                                                                               // null arg?
                                                                                                                            // no
                        std::string str(p_ArgStrings[iArg]);                                                                // convert char* to std::string

                        if (argString[0] != '[' || argString[argString.length()-1] != ']') {                                // starts with '[' and ends with ']'?
                            strargs.push_back(argString);                                                                   // no - not shorthand - add option value to the new vector
                        }
                        else {                                                                                              // yes - shorthand - process it
                            argString = argString.substr(1, argString.length() - 2);                                        // yes - strip enclosing brackets

                            bool defaultAllowed = std::get<1>(*it);                                                         // ok to omit values?

                            if (argString.empty()) {                                                                        // have null parameter?
                                if (defaultAllowed) {                                                                       // yes, null - defaults allowed?
                                    // we don't know how many default values to put here - how many default values
                                    // we need to specify depends on the option, and possibly other option values
                                    // (e.g. 'notes-hdrs' for 'notes').  So for now we just push a single default
                                    // value and deal with it later

                                    strargs.push_back(NOT_PROVIDED);                                                        // "not provided" indicator"
                                }
                                else {                                                                                      // no - defaults not allowed
                                    errStr = ERR_MSG(ERROR::MISSING_VALUE) + std::string(" for option '") + optionName + std::string("'"); // error
                                }
                            }
                            else {                                                                                          // non-null parameter
                                size_t start = 0;                                                                           // start position
                                size_t pos   = 0;                                                                           // current position
                                size_t idx   = 0;                                                                           // vector index of parameter
                                while (start < argString.length() && pos != std::string::npos) {                            // comma found before the end of the string?
                                                                                                                            // yes
                                    pos = argString.find(",", start);                                                       // next comma
                                    if (pos == std::string::npos) pos = argString.length();                                 // last character?
                                                                                                        
                                    if ((pos - start) > 0) {                                                                // non-zero length string?
                                        strargs.push_back(argString.substr(start, pos - start));                            // yes - grab it
                                    }
                                    else {                                                                                  // empty value
                                        if (defaultAllowed) {                                                               // defaults allowed?
                                            strargs.push_back(NOT_PROVIDED);                                                // "not provided" indicator"
                                        }
                                        else {                                                                              // no - defaults not allowed
                                            errStr = ERR_MSG(ERROR::MISSING_VALUE) + std::string(" for option '") + optionName + std::string("'"); // error
                                        }
                                    }
                                    start = pos + 1;                                                                        // next start
                                    idx++;                                                                                  // next vector index
                                }

                                if (argString[argString.length() - 1] == ',') {                                             // trailing comma in shorthand values?
                                    if (defaultAllowed) {                                                                   // defaults allowed?
                                        strargs.push_back(NOT_PROVIDED);                                                    // "not provided" indicator"
                                    }
                                    else {                                                                                  // no - defaults not allowed
                                        errStr = ERR_MSG(ERROR::MISSING_VALUE) + std::string(" for option '") + optionName + std::string("'"); // error
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        if (!errStr.empty()) break;                                                                                         // stop on error
    }

    // return arg count and arg strings to caller
    // if there was an error the strings returned may not be the complete command line, so should not be used

    return std::make_tuple(errStr, (int)strargs.size(), strargs);
}


/*
 * Parse the options provided by the user
 * 
 * We first expand any shorthand notation the user might have used (for options that allow shorthand
 * notation).
 * 
 * Before we give the options to boost we need to determine if the user passed any ranges or sets and, 
 * if they did, handle those - boost doesn't know anything about them.
 *
 * A range is allowed only for numeric options (i.e. INT or FLOAT types), but is not  allowed for all 
 * numeric options (e.g. --log-level)
 * A set is allowed for numeric, string, and bool options - but not all of them (e.g. --quiet)
 * 
 * We define a vector of options excluded from the range and set constructs (one vector each).  We don't
 * need to exclude non-numeric options from range here - that is done later - here we just exclude options
 * for which range/set makes no sense.  We have defined vectors of option names that are excluded from ranges
 * (m_RangeExcluded) and sets (m_SetExcluded).
 * 
 * 
 * std::string ParseOptionValues(int p_ArgCount, char *p_ArgStrings[], OptionsDescriptorT &p_OptionsDescriptor)
 * 
 * 
 * @param   [IN]    p_ArgCount                  The number of argument strings. (note below for p_ArgStrings)
 * @param   [IN]    p_ArgStrings                The argument strings.   The first argument string is expected
 *                                              (by boost) to be the executable name (boost expects the arguments
 *                                              to be commandline arguments passed to main())
 * @param   [IN]    p_OptionsDescriptor         Struct containing options descriptions.  This struct holds the
 *                                              boost options_description object, the option valued, and a struct
 *                                              containing the complex option values (the ranges and sets)
 * @return                                      String containing an error string
 *                                              If no error occurred the return string will be the empty string 
 */
std::string Options::ParseOptionValues(int p_ArgCount, char *p_ArgStrings[], OptionsDescriptorT &p_OptionsDescriptor) {

    std::string errStr = "";                                                                                                // initially

    int argCount;                                                                                                           // number or arg strings
    std::vector<std::string> sArgStrings;                                                                                   // arg strings - as std::strings


    //********************************************************************//
    // first expand any shorthand notation used in the options            //
    // if this returns an error, we return immediately from this function //
    //********************************************************************//

    std::tie(errStr, argCount, sArgStrings) = ExpandShorthandOptionValues(p_ArgCount, p_ArgStrings);                        // expand any shorthand option specifications

    if (!errStr.empty()) return errStr;                                                                                     // stop on error

    //********************************************************************//
    // shorthand notation expanded - proceed with parsing                 //
    //********************************************************************//


    std::vector<char const *> args {};                                                                                      // copy string vector to char * vector
    for (size_t idx = 0; idx < sArgStrings.size(); idx++) {
        args.push_back(sArgStrings[idx].c_str());
    }
                 /***** << do *not* try this at home >> *****/
    char **argStrings = const_cast<char**>(args.data());                                                                    // arg strings - as array of char*

    try {

        p_OptionsDescriptor.optionsSpecified    = {};                                                                       // initially
        p_OptionsDescriptor.complexOptionValues = {};                                                                       // initially

        std::string  optionName        = "";                                                                                // option name
        COMPLEX_TYPE type              = COMPLEX_TYPE::NONE;                                                                // complex arg type (range, set, neither/none)
        std::vector<std::string> parms = {};                                                                                // the range or set parameters

        for (size_t iArg = 1; iArg < (size_t)argCount; iArg++) {                                                            // for each arg string

            if (iArg <= 1) continue;                                                                                        // step over the executable name

            type = COMPLEX_TYPE::NONE;                                                                                      // initially

            optionName = argStrings[iArg - 1];                                                                              // get the option name for the argument we're processing
            optionName = utils::ToLower(utils::trim(optionName));                                                           // downshift and trim whitespace
            if (optionName[0] == '-') optionName.erase(0, optionName.find_first_not_of("-"));                               // remove the "-" or "--"

            if (argStrings[iArg] != nullptr) {                                                                              // null arg?
                                                                                                                            // no
                std::string str(argStrings[iArg]);                                                                          // convert char* to std::string
                str = utils::ToLower(utils::trim(str));                                                                     // downshift and trim whitespace

                // check for RANGE or SET
                // range is indicated by 'range[start,count,inc]', 'r[start,count,inc]', or just '[start,count,inc]'

                if ((str[0] == '[') || (str.rfind("r[", 0) == 0) || (str.rfind("range[", 0) == 0)) {                        // starts with '[', 'r[' or 'range[', so...
                    if (str[str.length()-1] == ']') {                                                                       // ... needs to end with ']' to be a valid RANGE
                        type = COMPLEX_TYPE::RANGE;                                                                         // it did - so RANGE

                        // check for RANGE requested for option in range excluded list
                        if (std::find(m_RangeExcluded.begin(), m_RangeExcluded.end(), optionName) != m_RangeExcluded.end()) {
                            errStr = ERR_MSG(ERROR::ARGUMENT_RANGE_NOT_SUPPORTED) + std::string(" for option '") + optionName + std::string("'");
                        }
                    }
                    else {
                        errStr = ERR_MSG(ERROR::MISSING_RIGHT_BRACKET) + std::string(" for option '") + optionName + std::string("'");
                    }
                }

                if (errStr.empty()) {                                                                                       // still ok?
                                                                                                                            // yes
                    // set is indicated by 'set[elem1,elem2,...,elemN]', or 's[elem1,elem2,...,elemN]'

                    if ((str.rfind("s[", 0) == 0) || (str.rfind("set[", 0) == 0)) {                                         // starts with 's[' or 'set[', so ...
                        if (str[str.length()-1] == ']') {                                                                   // ... needs to end with ']' to be a valid SET
                            type = COMPLEX_TYPE::SET;                                                                       // it did - so SET

                            // check for SET requested for option in set excluded list
                            if (std::find(m_SetExcluded.begin(), m_SetExcluded.end(), optionName) != m_SetExcluded.end())
                                errStr = ERR_MSG(ERROR::ARGUMENT_SET_NOT_SUPPORTED) + std::string(" for option '") + optionName + std::string("'");
                        }
                        else {
                            errStr = ERR_MSG(ERROR::MISSING_RIGHT_BRACKET) + std::string(" for option '") + optionName + std::string("'");
                        }
                    }
                }

                if (errStr.empty() && type != COMPLEX_TYPE::NONE) {                                                         // range or set?
                                                                                                                            // yes
                    // we have what looks like a 'range' or 'set' argument
                    // for now, just stash the details away and substitute the
                    // first value for the argument so we can check parsing

                    // look for comma separated values - there will be no 
                    // spaces - the OS/shell would have complained...

                    if (str.rfind("range", 0) == 0) str.erase(0, 5);                                                        // strip 'range' (range indicator) if present
                    if (str.rfind("set", 0) == 0) str.erase(0, 3);                                                          // strip 'set' (set indicator) if present
                    if (str[0] == 'r' || str[0] == 's') str.erase(0, 1);                                                    // strip 'r' or 's' (range or set indicator) if present
                    str = str.substr(1, str.length() - 2);                                                                  // strip enclosing brackets (must be present)

                    if (str.length() == 0 || str[str.length() - 1] == ',') {
                        errStr = ERR_MSG(ERROR::MISSING_VALUE) + std::string(" for option '") + optionName + std::string("'"); // no values, or trailing comma is an error
                    } 
                    else {

                        parms.clear();                                                                                      // start empty

                        size_t start = 0;                                                                                   // start position
                        size_t pos   = 0;                                                                                   // current position
                        while (errStr.empty() && start < str.length() && pos != std::string::npos) {                        // comma found before the end of the string?

                            std::string value = "";                                                                         // value

                            pos = str.find(",", start);                                                                     // next comma
                                                                                                        
                            if ((pos - start) > 0) {                                                                        // non-zero length string?
                                value = str.substr(start, pos - start);                                                     // yes - grab it
                                parms.push_back(value);                                                                     // store value

                                start = pos + 1;                                                                            // next start
                            }
                            else {                                                                                          // empty value - stop
                                errStr = ERR_MSG(ERROR::MISSING_VALUE) + std::string(" for option '") + optionName + std::string("'"); // error
                            }
                        }

                        if (errStr.empty()) {                                                                               // still ok?
                                                                                                                            // yes
                            if (type == COMPLEX_TYPE::RANGE && parms.size() != 3) {                                         // ranges require exactly 3 parameters
                                errStr = ERR_MSG(ERROR::ARGUMENT_RANGE_NUM_PARMS);                                          // error
                            }
                            else {
                                // if range, then we have 3 parameters (checked above)
                                // if set, we have at least one parameter (check earlier), so we're good

                                RangeOrSetDescriptorT details = {type, TYPENAME::NONE, parms, {}, 0};                       // dummy values for datatype and numerical parms
                                p_OptionsDescriptor.complexOptionValues.push_back(std::make_tuple(optionName, details));    // store the range/set

                                strncpy(argStrings[iArg], parms[0].c_str(), parms[0].length());                             // replace arg value (temporarily)
                                argStrings[iArg][parms[0].length()] = '\0';
                            }
                        }
                    }
                }          
            }
            if (!errStr.empty()) break;                                                                                     // stop parsing if error encountered
        }

        // if we've found errors, don't bother parsing

        if (errStr.empty()) {                                                                                               // no need if we've already flagged an error

            // boost parse_command_line() expects the first arg to be the program name
            // (it thinks it is getting the values that were passed to main() from the 
            // OS/shell), so for options from a grid file we insert a dummy argument as 
            // arg[0] and set the argument count appropriately.
            //
            // if valid ranges or sets were specified by the user they've been temporariliy
            // replaced for the boost parse, but if they were not valid they've been left
            // in the argument strings that will be passed to boost - so boost will fail
            // and complain about the offending parameter (which is what we want)

            po::parsed_options const parsedOptions = po::parse_command_line(argCount, argStrings, p_OptionsDescriptor.optionDescriptions, cls::unix_style|cls::case_insensitive); // parse user-supplied options
            po::store(parsedOptions, p_OptionsDescriptor.optionValues.m_VM);                                              // store parsed options into variable map
            po::notify(p_OptionsDescriptor.optionValues.m_VM);                                                            // populate the variables with option values

            // this is our opportunity to distinguish beteen "-h" and "--help" (if specified)
            for (auto& entry : parsedOptions.options) {
                po::option_description const& opt = p_OptionsDescriptor.optionDescriptions.find(entry.string_key, false, false, false);
                std::string originalTok = entry.original_tokens[0];
                std::string thisTok = utils::ToLower(utils::trim(originalTok));

                if (!thisTok.empty()) {
                    std::string shortOpt = utils::ToLower(opt.canonical_display_name(cls::allow_dash_for_short));
                    std::string longOpt  = utils::ToLower(opt.canonical_display_name(cls::allow_long));

                    if ((shortOpt == "-h") || (longOpt == "--help")) {
                        p_OptionsDescriptor.optionValues.m_ShortHelp = thisTok == "-h";
                    }

                    if (originalTok[0] == '-') originalTok.erase(0, originalTok.find_first_not_of("-"));                    // remove the "-" or "--"
                    if (thisTok[0]     == '-') thisTok.erase(0, thisTok.find_first_not_of("-"));                            // remove the "-" or "--"
                    if (longOpt[0]     == '-') longOpt.erase(0, longOpt.find_first_not_of("-"));                            // remove the "-" or "--"
                    if (shortOpt[0]    == '-') shortOpt.erase(0, shortOpt.find_first_not_of("-"));                          // remove the "-" or "--"

                    p_OptionsDescriptor.optionsSpecified.push_back(std::make_tuple(originalTok, thisTok, longOpt, shortOpt));        
                }
            }

            // If we've made it this far then boost parsed the command-line arguments ok.
            //
            // If there were any ranges or sets specified by the user we can now work out the 
            // data types of the options for which they (the ranges/sets) were specified and 
            // sanity check them.
            //
            // First iterate through the specified ranges and sets to sanity check - and
            // manually set the value for the option to the first value in the range or set.
            // Need to check:
            //     - ranges have not been specified for non-numeric options
            //     - range values are all numeric
            //     - values for ranges and sets match the data type of the option
            //       (i.e. INTs for integer options, FLOATs for fp options, STRINGs for string options)
            //
            // It would be preferable to range check values against valid values for specific
            // options here but, for now at least, too problematic - we'll just let the
            // evolution fail if an option value is bad (we don't want to read through every
            // record of a grid file and check...)


            // Now's a good time to pull out SSE-only/BSE-only options that we want to ignore
            // We really only need to remove the complex options - the COMPAS code naturally
            // ignores BSE options when in SSE mode and vice-versa.  The only reason we want
            // to remove the complex options is that we don't want to iterate over them and
            // produce too many stars/binaries (see explanation of m_SSEOnly and m_BSEOnly
            // in constants.h).
            //
            // First, loop through the complex options and identify the indices of the options
            // we want to remove.

            bool bseMode = utils::ToLower(m_CmdLine.optionValues.m_EvolutionMode.typeString) == "bse";                      // mode

            std::vector<size_t> removeOpts = {};                                                                            // vector of option indices to be removed
            for (size_t iOpt = 0; iOpt < p_OptionsDescriptor.complexOptionValues.size(); iOpt++) {                          // for each range or set specified
        
                std::string opt = std::get<0>(p_OptionsDescriptor.complexOptionValues[iOpt]);                               // option name

                if (bseMode) {                                                                                              // BSE?
                    if (std::find(m_SSEOnly.begin(), m_SSEOnly.end(), opt) != m_SSEOnly.end()) removeOpts.push_back(iOpt);  // remove SSE only option
                }
                else {                                                                                                      // SSE
                    if (std::find(m_BSEOnly.begin(), m_BSEOnly.end(), opt) != m_BSEOnly.end()) removeOpts.push_back(iOpt);  // remove BSE only option
                }
            }

            // We now know the indices of options we want to remove from the complex options,
            // so we can loop through and remove them.
            // This only works because the removeOpts vector is populated in the same order as
            // the complexOptionValues vector - that allows me to loop through complexOptionValues
            // in reverse order to delete elements (otherwise the index numbers would change from 
            // under me as I was deleting them)

            size_t removeCount(removeOpts.size());
            if (removeCount > 0) {                                                                                          // anything to remove?
                for (size_t iOpt = removeCount - 1; iOpt >= 0; iOpt--) {                                                    // loop through all complex options identified to be removed
                    p_OptionsDescriptor.complexOptionValues.erase(p_OptionsDescriptor.complexOptionValues.begin() + removeOpts[iOpt]);                                                    // erase it
                }
            }


            // ok, now process the list of complex options
            RangeOrSetDescriptorT details = {};
            std::string longOptionName;

            size_t count = p_OptionsDescriptor.complexOptionValues.size();                                                  // count of complex values (ranges or sets)
            for (size_t idx = 0; idx < count; idx++) {                                                                      // for each range or set specified

                optionName     = std::get<0>(p_OptionsDescriptor.complexOptionValues[idx]);                                 // the option name
                longOptionName = optionName;
                details        = std::get<1>(p_OptionsDescriptor.complexOptionValues[idx]);                                 // range/set details for this optionName
                type           = details.type;                                                                              // range or set
                parms          = details.parameters;                                                                        // range/set parameter values (as strings)

                // we want to use the long name of the option for this next bit
                // look for the option in the options specified - if it's not 
                // found the find in m_VM will fail for us...
                std::string longOptionName = optionName;
                auto thisIt = std::find_if(
                    p_OptionsDescriptor.optionsSpecified.begin(), p_OptionsDescriptor.optionsSpecified.end(), [&optionName](const OPTIONSTR& e) {
                        return std::get<0>(e) == optionName || std::get<1>(e) == optionName || std::get<2>(e) == optionName || std::get<3>(e) == optionName;
                    }
                );
                if (thisIt != p_OptionsDescriptor.optionsSpecified.end()) {
                    longOptionName = std::get<2>(*thisIt);
                    p_OptionsDescriptor.complexOptionValues[idx] = std::make_tuple(longOptionName, details);
                }

                po::variables_map::const_iterator it = p_OptionsDescriptor.optionValues.m_VM.find(longOptionName);          // yes - find the option in the boost variables map
                if (it != p_OptionsDescriptor.optionValues.m_VM.end()) {                                                    // found?
                    TYPENAME dataType = TYPENAME::NONE;                                                                     // yes
                    std::tie(dataType, std::ignore, std::ignore, std::ignore) = OptionAttributes(p_OptionsDescriptor.optionValues.m_VM, it); // data type
                    details.dataType = dataType;                                                                            // set data type

                    if (idx == (count - 1)) details.currPos = 0;                                                            // initial position for inner iterator

                    if (type == COMPLEX_TYPE::RANGE) {                                                                      // RANGE?
                        if (!IsSupportedNumericDataType(dataType)) {                                                        // yes - numeric? 
                            errStr = ERR_MSG(ERROR::ARGUMENT_RANGE_NOT_SUPPORTED) + std::string(" for option '") + optionName + std::string("'"); // no - that's not ok
                        }
                        else {                                                                                              // yes - numeric
                                                                                                                            // yes - determine numerical range parameters
                            switch (dataType) {                                                                             // which data type?

                                case TYPENAME::INT: {                                                                       // INT
                                    std::string complaint1 = ERR_MSG(ERROR::ARGUMENT_RANGE_PARMS_EXPECTED_INT) + std::string(" for option '") + optionName + std::string("'");
                                    std::string complaint2 = ERR_MSG(ERROR::ARGUMENT_RANGE_COUNT_EXPECTED_ULINT) + std::string(" for option '") + optionName + std::string("'");
                                    try {
                                        RangeParameterT tmp = {0.0L};                                                       // dummy value
                                        details.rangeParms = {tmp, tmp, tmp};                                               // create the vector

                                        size_t lastChar;
                                        details.rangeParms[0].iVal   = std::stoi(details.parameters[0], &lastChar);         // integer start
                                        COMPLAIN_IF(lastChar != details.parameters[0].length(), complaint1);                // not a valid int
                                        details.rangeParms[2].iVal   = std::stoi(details.parameters[2], &lastChar);         // integer inc
                                        COMPLAIN_IF(lastChar != details.parameters[2].length(), complaint1);                // not a valid int

                                        try {
                                            size_t lastChar;
                                            details.rangeParms[1].ulVal = std::stoul(details.parameters[1], &lastChar);     // unsigned long int count
                                            COMPLAIN_IF(lastChar != details.parameters[1].length(), complaint2);            // not a valid unsigned long int

                                            p_OptionsDescriptor.complexOptionValues[idx] = std::make_tuple(longOptionName, details); // reset values
                                        }
                                        catch (const std::out_of_range& e) {                                                // not a valid unsigned long int
                                            errStr = complaint2;
                                        }
                                        catch (const std::invalid_argument& e) {                                            // not a valid unsigned long int
                                            errStr = complaint2;
                                        }
                                    }
                                    catch (const std::out_of_range& e) {                                                    // not a valid int
                                        errStr = complaint1;
                                    }
                                    catch (const std::invalid_argument& e) {                                                // not a valid int
                                        errStr = complaint1;
                                    }
                                } break;

                                case TYPENAME::LONGINT: {                                                                   // LONG INT
                                    std::string complaint1 = ERR_MSG(ERROR::ARGUMENT_RANGE_PARMS_EXPECTED_LINT) + std::string(" for option '") + optionName + std::string("'");
                                    std::string complaint2 = ERR_MSG(ERROR::ARGUMENT_RANGE_COUNT_EXPECTED_ULINT) + std::string(" for option '") + optionName + std::string("'");
                                    try {
                                        RangeParameterT tmp = {0.0};                                                        // dummy value
                                        details.rangeParms = {tmp, tmp, tmp};                                               // create the vector

                                        size_t lastChar;
                                        details.rangeParms[0].lVal = std::stol(details.parameters[0], &lastChar);           // unsigned long int start
                                        COMPLAIN_IF(lastChar != details.parameters[0].length(), complaint1);                // not a valid long int
                                        details.rangeParms[2].lVal = std::stol(details.parameters[2], &lastChar);           // unsigned long int inc
                                        COMPLAIN_IF(lastChar != details.parameters[2].length(), complaint1);                // not a valid long int

                                        try {
                                            size_t lastChar;
                                            details.rangeParms[1].ulVal = std::stoul(details.parameters[1], &lastChar);     // unsigned long int count
                                            COMPLAIN_IF(lastChar != details.parameters[1].length(), complaint2);            // not a valid unsigned long int

                                            p_OptionsDescriptor.complexOptionValues[idx] = std::make_tuple(longOptionName, details); // reset values
                                        }
                                        catch (const std::out_of_range& e) {                                                // not a valid unsigned long int
                                            errStr = complaint2;
                                        }
                                        catch (const std::invalid_argument& e) {                                            // not a valid unsigned long int
                                            errStr = complaint2;
                                        }
                                    }
                                    catch (const std::out_of_range& e) {                                                    // not a valid long int
                                        errStr = complaint1;
                                    }
                                    catch (const std::invalid_argument& e) {                                            // not a valid long int
                                        errStr = complaint1;
                                    }
                                } break;

                                case TYPENAME::ULONGINT: {                                                                  // UNSIGNED LONG INT
                                    std::string complaint1 = ERR_MSG(ERROR::ARGUMENT_RANGE_PARMS_EXPECTED_ULINT) + std::string(" for option '") + optionName + std::string("'");
                                    std::string complaint2 = ERR_MSG(ERROR::ARGUMENT_RANGE_COUNT_EXPECTED_ULINT) + std::string(" for option '") + optionName + std::string("'");
                                    try {
                                        RangeParameterT tmp = {0.0};                                                        // dummy value
                                        details.rangeParms = {tmp, tmp, tmp};                                               // create the vector

                                        size_t lastChar;
                                        details.rangeParms[0].ulVal = std::stoul(details.parameters[0], &lastChar);         // unsigned long int start
                                        COMPLAIN_IF(lastChar != details.parameters[0].length(), complaint1);                // not a valid unsigned long int
                                        details.rangeParms[2].ulVal = std::stoul(details.parameters[2], &lastChar);         // unsigned long int inc
                                        COMPLAIN_IF(lastChar != details.parameters[2].length(), complaint1);                // not a valid unsigned long int

                                        try {
                                            size_t lastChar;
                                            details.rangeParms[1].ulVal = std::stoul(details.parameters[1], &lastChar);     // unsigned long int count
                                            COMPLAIN_IF(lastChar != details.parameters[2].length(), complaint2);            // not a valid unsigned long int

                                            p_OptionsDescriptor.complexOptionValues[idx] = std::make_tuple(longOptionName, details); // reset values
                                        }
                                        catch (const std::out_of_range& e) {                                                // not a valid unsigned long int
                                            errStr = complaint2;
                                        }
                                        catch (const std::invalid_argument& e) {                                            // not a valid unsigned long int
                                            errStr = complaint2;
                                        }
                                    }
                                    catch (const std::out_of_range& e) {                                                    // not a valid unsigned long int
                                        errStr = complaint1;
                                    }
                                    catch (const std::invalid_argument& e) {                                                // not a valid unsigned long int
                                        errStr = complaint1;
                                    }
                                } break;

                                case TYPENAME::FLOAT: {                                                                     // FLOAT
                                    std::string complaint1 = ERR_MSG(ERROR::ARGUMENT_RANGE_PARMS_EXPECTED_FP) + std::string(" for option '") + optionName + std::string("'");
                                    std::string complaint2 = ERR_MSG(ERROR::ARGUMENT_RANGE_COUNT_EXPECTED_ULINT) + std::string(" for option '") + optionName + std::string("'");
                                    try {
                                        RangeParameterT tmp = {0.0};                                                        // dummy value
                                        details.rangeParms = {tmp, tmp, tmp};                                               // create the vector

                                        size_t lastChar;
                                        details.rangeParms[0].dVal = std::stof(details.parameters[0], &lastChar);           // floating point start
                                        COMPLAIN_IF(lastChar != details.parameters[0].length(), complaint1);                // not a valid float
                                        details.rangeParms[2].dVal = std::stof(details.parameters[2], &lastChar);           // floating point inc
                                        COMPLAIN_IF(lastChar != details.parameters[2].length(), complaint1);                // not a valid float

                                        try {
                                            size_t lastChar;
                                            details.rangeParms[1].ulVal = std::stoul(details.parameters[1], &lastChar);     // unsigned long int count
                                            COMPLAIN_IF(lastChar != details.parameters[1].length(), complaint2);            // not a valid unsigned long int

                                            p_OptionsDescriptor.complexOptionValues[idx] = std::make_tuple(longOptionName, details); // reset values
                                        }
                                        catch (const std::out_of_range& e) {                                                // not a valid unsigned long int
                                            errStr = complaint2;
                                        }
                                        catch (const std::invalid_argument& e) {                                            // not a valid unsigned long int
                                            errStr = complaint2;
                                        }
                                    }
                                    catch (const std::out_of_range& e) {                                                    // not a valid floating point number
                                        errStr = complaint1;
                                    }
                                    catch (const std::invalid_argument& e) {                                                // not a valid floating point number
                                        errStr = complaint1;
                                    }
                                } break;

                                case TYPENAME::DOUBLE: {                                                                    // DOUBLE
                                    std::string complaint1 = ERR_MSG(ERROR::ARGUMENT_RANGE_PARMS_EXPECTED_FP) + std::string(" for option '") + optionName + std::string("'");
                                    std::string complaint2 = ERR_MSG(ERROR::ARGUMENT_RANGE_COUNT_EXPECTED_ULINT) + std::string(" for option '") + optionName + std::string("'");
                                    try {
                                        RangeParameterT tmp = {0.0};                                                        // dummy value
                                        details.rangeParms = {tmp, tmp, tmp};                                               // create the vector

                                        size_t lastChar;
                                        details.rangeParms[0].dVal = std::stod(details.parameters[0], &lastChar);           // floating point start
                                        COMPLAIN_IF(lastChar != details.parameters[0].length(), complaint1);                // not a valid double
                                        details.rangeParms[2].dVal = std::stod(details.parameters[2], &lastChar);           // floating point inc
                                        COMPLAIN_IF(lastChar != details.parameters[2].length(), complaint1);                // not a valid double

                                        try {
                                            size_t lastChar;
                                            details.rangeParms[1].ulVal = std::stoul(details.parameters[1], &lastChar);     // unsigned long int count
                                            COMPLAIN_IF(lastChar != details.parameters[1].length(), complaint2);            // not a valid unsigned long int

                                            p_OptionsDescriptor.complexOptionValues[idx] = std::make_tuple(longOptionName, details); // reset values
                                        }
                                        catch (const std::out_of_range& e) {                                                // not a valid unsigned long int
                                            errStr = complaint2;
                                        }
                                        catch (const std::invalid_argument& e) {                                            // not a valid unsigned long int
                                            errStr = complaint2;
                                        }
                                    }
                                    catch (const std::out_of_range& e) {                                                    // not a valid floating point number
                                        errStr = complaint1;
                                    }
                                    catch (const std::invalid_argument& e) {                                                // not a valid floating point number
                                        errStr = complaint1;
                                    }
                                } break;

                                case TYPENAME::LONGDOUBLE: {                                                                // LONG DOUBLE
                                    std::string complaint1 = ERR_MSG(ERROR::ARGUMENT_RANGE_PARMS_EXPECTED_LFP) + std::string(" for option '") + optionName + std::string("'");
                                    std::string complaint2 = ERR_MSG(ERROR::ARGUMENT_RANGE_COUNT_EXPECTED_ULINT) + std::string(" for option '") + optionName + std::string("'");
                                    try {
                                        RangeParameterT tmp = {0.0};                                                        // dummy value
                                        details.rangeParms = {tmp, tmp, tmp};                                               // create the vector

                                        size_t lastChar;
                                        details.rangeParms[0].ldVal = std::stold(details.parameters[0], &lastChar);         // long double start
                                        COMPLAIN_IF(lastChar != details.parameters[0].length(), complaint1);                // not a valid long double
                                        details.rangeParms[2].ldVal = std::stold(details.parameters[2], &lastChar);         // long double inc
                                        COMPLAIN_IF(lastChar != details.parameters[2].length(), complaint1);                // not a valid long double

                                        try {
                                            size_t lastChar;
                                            details.rangeParms[1].ulVal = std::stoul(details.parameters[1], &lastChar);     // unsigned long int count
                                            COMPLAIN_IF(lastChar != details.parameters[1].length(), complaint2);            // not a valid unsigned long int

                                            p_OptionsDescriptor.complexOptionValues[idx] = std::make_tuple(longOptionName, details); // reset values
                                        }
                                        catch (const std::out_of_range& e) {                                                // not a valid unsigned long int
                                            errStr = complaint2;
                                        }
                                        catch (const std::invalid_argument& e) {                                            // not a valid unsigned long int
                                            errStr = complaint2;
                                        }
                                    }
                                    catch (const std::out_of_range& e) {                                                    // not a valid long double number
                                        errStr = complaint1;
                                    }
                                    catch (const std::invalid_argument& e) {                                                // not a valid long double number
                                        errStr = complaint1;
                                    }
                                } break;
                                
                                default:                                                                                    // that's a problem...
                                    COMPLAIN(ERR_MSG(ERROR::INVALID_DATA_TYPE));                                            // complain
                            }
                        }
                    }
                    else {                                                                                                  // SET
                        // check for numeric/bool data types only that all set parameters are numeric/bool
                        // can't check for string data types 
                        
                        if (IsSupportedNumericDataType(dataType)) {                                                         // numeric?
                            
                            for (size_t ip = 0; ip < parms.size(); ip++) {                                                  // yes - for each set parameter specified

                                if ((dataType == TYPENAME::INT        && !utils::IsINT(parms[ip]))         ||               // INT?
                                    (dataType == TYPENAME::LONGINT    && !utils::IsLONGINT(parms[ip]))     ||               // LONG INT?
                                    (dataType == TYPENAME::ULONGINT   && !utils::IsULONGINT(parms[ip]))    ||               // UNSIGNED LONG INT?
                                    (dataType == TYPENAME::FLOAT      && !utils::IsFLOAT(parms[ip]))       ||               // FLOAT?
                                    (dataType == TYPENAME::DOUBLE     && !utils::IsDOUBLE(parms[ip]))      ||               // DOUBLE?
                                    (dataType == TYPENAME::LONGDOUBLE && !utils::IsLONGDOUBLE(parms[ip]))) {                // LONG DOUBLE?
                                    errStr = ERR_MSG(ERROR::ARGUMENT_SET_EXPECTED_NUMERIC) + std::string(" for option '") + optionName + std::string("'"); // no - that's not ok
                                    break;
                                }
                            }
                        }
                        else if (dataType == TYPENAME::BOOL) {                                                              // bool?
                            // we allow boolean set values to be specified as 1|0, TRUE|FALSE, YES|NO, ON|OFF - but not mixed
                            // i.e. all must be 0|1, or all must be TRUE|FALSE etc. - case is not significant (downshifted already)

                            size_t checks[4] = {0};                                                                         // all (Boost) boolean representations
                            for (size_t ip = 0; ip < parms.size(); ip++) {                                                  // yes - for each set parameter specified
                                int check = std::abs(utils::IsBOOL(parms[ip]));                                             // check parm for valid BOOL
                                if (check == 0) break;
                                checks[check - 1]++;
                            }

                            int validCheck = 0;                                                                             // any valid boolean?
                            bool errorhere = false;
                            for (size_t iChk = 0; iChk < 4; iChk++) {                                                       // for all (Boost) boolean representations
                                validCheck += checks[iChk];                                                                 // valid booleans
                                if (checks[iChk] != 0 && checks[iChk] != parms.size()) {                                    // all parms boolean, and consistent?
                                    errorhere = true;                                                                           // no - that's not ok
                                    break;
                                }
                            }
                            if (errorhere || validCheck == 0) errStr = ERR_MSG(ERROR::ARGUMENT_SET_EXPECTED_BOOL) + std::string(" for option '") + optionName + std::string("'");
                        }

                        if (errStr.empty()) {                                                                               // all ok?
                            p_OptionsDescriptor.complexOptionValues[idx] = std::make_tuple(longOptionName, details);        // yes - reset values
                        }
                    }
                }
                else {                                                                                                      // option not found in boost variables map
                    errStr = ERR_MSG(ERROR::BOOST_OPTION_INTERNAL_ERROR) + std::string(" for option '") + optionName + std::string("'"); // that can't be good...
                }

                if (!errStr.empty()) break;                                                                                 // stop on error
            }


            // if we've made it this far we've parsed everything, and now is the time to fix up any default values
            // for vector options.  Any vector option values that we know were not provided will have the value
            // NOT_PROVIDED, so they can just be replaced with whatever the default value should be.  However,
            // we couldn't really know until now how many values we should expect for vector options like 'notes',
            // because (in that case) the number of values expected is the number of 'notes-hdrs' specified by the
            // user, and that we didn't really know until now (and moreover, we didn't know what the command-line
            // values were until now, so couldn't default grid-line values to command-line values until now).
            //
            // but - if we've found errors, don't bother...

            if (errStr.empty()) {                                                                                           // no need if we've already flagged an error

                for (auto& elem: m_ShorthandAllowed) {                                                                      // for each option for which shorthand is allowed
                    std::string optionName     = std::get<0>(elem);                                                         // option name
                    bool        defaultAllowed = std::get<1>(elem);                                                         // omissions allowed?
                    std::string defaultString  = std::get<2>(elem);                                                         // the COMPAS default for omissions

                    if (defaultAllowed) {                                                                                   // omission allowed for this option?
                                                                                                                            // yes
                        switch (_(optionName.c_str())) {                                                                    // which option?

                            // '--notes' is the only option affected at the moment

                            case _("notes"):                                                                                // notes
                                for (size_t idx = 0; idx < NotesHdrs().size(); idx++) {                                     // for each specified notes-hdr
                                    if (idx < p_OptionsDescriptor.optionValues.m_Notes.size()) {                            // have parsed value?
                                        if (p_OptionsDescriptor.optionValues.m_Notes[idx] == NOT_PROVIDED) {                // yes - notes value provided?
                                                                                                                            // no - get default
                                            if (p_OptionsDescriptor.optionsOrigin == OPTIONS_ORIGIN::CMDLINE) {             // from command line?
                                                p_OptionsDescriptor.optionValues.m_Notes[idx] = defaultString;              // yes - use COMPAS default
                                            }
                                            else {                                                                          // no - grid file line
                                                p_OptionsDescriptor.optionValues.m_Notes[idx] = m_CmdLine.optionValues.m_Notes[idx]; // use command-line value
                                            }
                                        }
                                    }
                                    else {                                                                                  // no parsed value
                                        if (p_OptionsDescriptor.optionsOrigin == OPTIONS_ORIGIN::CMDLINE) {                 // from command line?
                                            p_OptionsDescriptor.optionValues.m_Notes.push_back(defaultString);              // yes - use COMPAS default
                                        }
                                        else {                                                                              // no - grid file line
                                            p_OptionsDescriptor.optionValues.m_Notes.push_back(m_CmdLine.optionValues.m_Notes[idx]); // use command-line value
                                        }
                                    }
                                }
                                break;

                            default:                                                                                        // default - shouldn't happen
                                break;                                                                                      // do nothing - the parse will fail
                        }
                    }
                }
            }

        }
    }
    catch (po::error& e) {                                                                                                  // program options exception
        errStr = e.what();                                                                                                  // set error string
    }
    catch (const std::string eStr) {                                                                                        // custom exception
        errStr = eStr;                                                                                                      // set error string
    }
    catch (...) {                                                                                                           // unhandled exception
        errStr = ERR_MSG(ERROR::UNHANDLED_EXCEPTION);                                                                       // set error string
    }

    return errStr;
}


/*
 * Initialise options service
 * 
 * Intitialises the options service.  Constructs options objects for the program options
 * (options that are specified only on the commandline and that cannot be specified in a 
 * grid file on a per object (star/binary) basis), and the grid file options (options that
 * can be specified in a grid file on a per object (star/binary) basis).
 * 
 * Initialises the commandline (program-level) options object, and the grid line (evolving 
 * object-level) options object.  Populates the commandline (program-level) options object
 * from the commandline arguments passed to main() - this object stays static throughout the
 * life of the program.
 * 
 * 
 * bool Options::Initialise(int p_ArgCount, char *p_ArgStrings[])
 * 
 * @param   [IN]    p_ArgCount                  Integer number of args passed in p_ArgStrings
 * @param   [IN]    p_ArgStrings                Arg strings - 1 per option
 *                                              Note that the first arg string is ignored (expected to be program name)
 * @return                                      Boolean status (true = ok, false = problem) (this function displays any error string)
 */
bool Options::Initialise(int p_ArgCount, char *p_ArgStrings[]) {

    bool ok = true;                                                                                                 // status - unless something changes

    // we need the strings in the following vectors to be downshifted for comparisons
    // we'll do it here once - just in case someone added non-lower case strings...

    for (size_t idx = 0; idx < m_GridLineExcluded.size(); idx++) m_GridLineExcluded[idx] = utils::ToLower(utils::trim(m_GridLineExcluded[idx]));
    for (size_t idx = 0; idx < m_SSEOnly.size();          idx++) m_SSEOnly[idx]          = utils::ToLower(utils::trim(m_SSEOnly[idx]));
    for (size_t idx = 0; idx < m_BSEOnly.size();          idx++) m_BSEOnly[idx]          = utils::ToLower(utils::trim(m_BSEOnly[idx]));
    for (size_t idx = 0; idx < m_RangeExcluded.size();    idx++) m_RangeExcluded[idx]    = utils::ToLower(utils::trim(m_RangeExcluded[idx]));
    for (size_t idx = 0; idx < m_SetExcluded.size();      idx++) m_SetExcluded[idx]      = utils::ToLower(utils::trim(m_SetExcluded[idx]));
    for (size_t idx = 0; idx < m_ShorthandAllowed.size(); idx--) m_ShorthandAllowed[idx] = std::make_tuple(utils::ToLower(utils::trim(std::get<0>(m_ShorthandAllowed[idx]))), std::get<1>(m_ShorthandAllowed[idx]), std::get<2>(m_ShorthandAllowed[idx]));

    try {

        m_CmdLine.optionValues.Initialise();                                                                        // initialise option variables for program-level options
        m_GridLine.optionValues.Initialise();                                                                       // initialise option variables for evolving object-level options

        po::options_description programLevelOptions("Program Options", 128);                                        // boost options descriptions object for program-level options
        ok = AddOptions(&m_CmdLine.optionValues, &programLevelOptions);                                             // ... add
        if (!ok) {                                                                                                  // ok?
            COMPLAIN(ERR_MSG(ERROR::BOOST_OPTION_CMDLINE));                                                         // no, complain - this throws an exception
        }
        else {                                                                                                      // yes, ok

            m_CmdLine.optionDescriptions.add(programLevelOptions);                                                  // commandline options - stays static throughout the life of the program
    
            // we parse the option values before handing them over to boost
            // boost knows nothing about shorthand, ranges and sets, so we have to handlde
            // them ourselves first
            m_CmdLine.complexOptionValues = {};                                                                     // no ranges or sets - unless we find them in the parse
            std::string errStr = ParseOptionValues(p_ArgCount, p_ArgStrings, m_CmdLine);                            // parse & populate the option values - specifically for ranges and sets
            if (!errStr.empty()) {                                                                                  // parsed ok?
                COMPLAIN(errStr);                                                                                   // no, complain - this throws an exception
            }
            else {
                errStr = m_CmdLine.optionValues.CheckAndSetOptions();                                               // yes - sanity check, and set, program-level values
                if (!errStr.empty()) {                                                                              // check ok?
                    COMPLAIN(errStr);                                                                               // no, complain - this throws an exception
                }
                else {

                    m_CmdLineOptionsDetails = OptionDetails(m_CmdLine);                                             // yes - get Run_Details contents

                    // initialise evolving object-level options.  The values of options specified in a grid file
                    // take precedence over the values of the same options specified on the commandline, but only
                    // for the object (star/binary) corresponding to the grid file record.
                    po::options_description objectLevelOptions("Program Options");                                  // boost options descriptions object for per object (star/binary) options
                    ok = AddOptions(&m_GridLine.optionValues, &objectLevelOptions);                                 // ... add
                    if (!ok) {                                                                                      // ok?
                        COMPLAIN(ERR_MSG(ERROR::BOOST_OPTION_GRIDLINE));                                            // no, complain - this throws an exception
                    }
                    else {                                                                                          // yes, ok
                        m_GridLine.optionDescriptions.add(objectLevelOptions);                                      // grid line options - stays static throughout the life of the program
                    }

                    // We now have the options the user entered at the commandline, including any expanded shorthand,
                    // ranges and/or sets, so this is where we stop the initialisation - from here we just play out 
                    // the options that are specified by any ranges and sets via the AdvanceCmdLineOptionValues() function.
                    //
                    // If the user has specified any ranges or sets we set the options to the first value in each
                    // range or set (already done by the time we get here).  Calls to AdvanceCmdLineOptionValues()
                    // will then advance the option values through the ranges and sets as required - *however*, the
                    // very first call to AdvanceCmdLineOptionValues() just returns the options as they are set here,
                    // so that a call to AdvanceCmdLineOptionValues() can be put in a loop, and the first evaluation
                    // of the loop will be the first star/binary - and if only 1 star/binary is required then no loop,
                    // and no call to AdvanceCmdLineOptionValues() is required (because the initial values have already
                    // been set).
                    //
                    // Note that there are analogous functions for object (star/binary) initialisation and retrievel
                    // option values: InitialiseEvolvingObject() and AdvanceGridLineOptionValues().  These functions 
                    // intitialise and retrieve options specified in grid file records.
                }
            }
        }  
    }
    catch (po::error& e) {                                                                                          // program options exception
        std::cerr << ERR_MSG(ERROR::PROGRAM_OPTIONS_ERROR) << ": " << e.what() << std::endl;                        // show the problem
        std::cerr << ERR_MSG(ERROR::SUGGEST_HELP) << std::endl;                                                     // suggest using --help
        ok = false;                                                                                                 // set status
    } 
    catch (const std::string eStr) {                                                                                // custom exception
        std::cerr << ERR_MSG(ERROR::PROGRAM_OPTIONS_ERROR) << ": " << eStr << std::endl;                            // show the problem
        std::cerr << ERR_MSG(ERROR::SUGGEST_HELP) << std::endl;                                                     // suggest using --help
        ok = false;                                                                                                 // set status
    }
    catch (...) {                                                                                                   // unhandled exception
        std::cerr << ERR_MSG(ERROR::PROGRAM_OPTIONS_ERROR) << ": " << ERR_MSG(ERROR::UNHANDLED_EXCEPTION) << std::endl; // show the problem
        std::cerr << ERR_MSG(ERROR::SUGGEST_HELP) << std::endl;                                                     // suggest using --help
        ok = false;                                                                                                 // set status
    }

    m_CmdLine.optionValues.m_Populated = ok;                                                                        // flag use

    return ok;
}


/*
 * Advances the command or grid line options to the next variation (depending
 * upon the parameters passed).  A "variation" is a combination of options 
 * defined by the option values, ranges, and sets the user specified on the
 * commandline or grid file line.
 * 
 * When the commandline or grid file line is parsed, the option values are set
 * to the initial variation of range and set values, then once that is processed
 * this function is called to advance to the next variation - the ranges and
 * sets are played out in order.
 * 
 * 
 * int AdvanceOptionVariation(OptionsDescriptorT &p_OptionsDescriptor)
 * 
 * @param   [IN]    p_OptionsDescriptor         Commandline or grid line options descriptor
 * @return                                      Int result:
 *                                                  -1: an error occurred
 *                                                   0: no more variations - all done
 *                                                   1: new variation applied - option values are set
 */
int Options::AdvanceOptionVariation(OptionsDescriptorT &p_OptionsDescriptor) {

    int retVal = 0;

    if (p_OptionsDescriptor.complexOptionValues.size() == 0) {          // more variations?
        // no - set calculated option defaults and return
        return p_OptionsDescriptor.optionValues.SetCalculatedOptionDefaults(BOOST_MAP::NO_UPDATE) == "" ? 0 : -1;
    }

    // Upon entry iterators for ranges and sets need to be advanced in order
    // to pick up the correct values to be loaded into the options.  Really
    // all we need to do is pick up the inner (fastest-counting) iterator (the
    // right-most in terms of placement on the commandline).  If we have to 
    // wrap-around to get that value, then we have to increment the next outer
    // (immediately left) iterator.
    
    // Traverse each of the complex option values and gather the option values.
    // We only need values for the options that have changing values.
    // Values have already been sanity checked by the time we get here.
    // We traverse the complex options values in reverse order because the
    // fastest change is to the right...

    bool stop  = false;                                                 // stop once we have all the values we need
    int idx = p_OptionsDescriptor.complexOptionValues.size() - 1;
    while (!stop && idx >= 0) {

        stop = true;                                                    // assume we have what we need

        std::string optionName        = std::get<0>(p_OptionsDescriptor.complexOptionValues[idx]);  
        RangeOrSetDescriptorT details = std::get<1>(p_OptionsDescriptor.complexOptionValues[idx]);
        details.currPos++;                                              // advance iterator

        if (details.type == COMPLEX_TYPE::SET) {                        // SET
            if (details.currPos >= int(details.parameters.size())) {    // currPos is set position - wrap?
                if (idx == 0) {                                         // outermost iterator?
                    retVal = 0;                                         // yes - we're done
                    break;
                }
                else {                                                  // no - wrap
                    details.currPos = 0;                                // ... back to zero
                    stop = false;                                       // ... and don't stop at this iterator
                }
            }
            std::string optionValue = details.parameters[details.currPos];  // option value (as string)

            switch (details.dataType) {                                 // which data type?
                // these should cover our options for now - but we might have to refine them over time

                case TYPENAME::INT: {                                   // INT
                    int thisVal = std::stoi(optionValue);
                    p_OptionsDescriptor.optionValues.ModifyVariableMap(p_OptionsDescriptor.optionValues.m_VM, optionName, thisVal);
                    po::notify(p_OptionsDescriptor.optionValues.m_VM);
                }   break;

                case TYPENAME::LONGINT: {                               // LONG INT
                    long int thisVal = std::stol(optionValue);
                    p_OptionsDescriptor.optionValues.ModifyVariableMap(p_OptionsDescriptor.optionValues.m_VM, optionName, thisVal);
                    po::notify(p_OptionsDescriptor.optionValues.m_VM);
                }   break;

                case TYPENAME::ULONGINT: {                              // UNSIGNED LONG INT
                    unsigned long int thisVal = std::stoul(optionValue);
                    p_OptionsDescriptor.optionValues.ModifyVariableMap(p_OptionsDescriptor.optionValues.m_VM, optionName, thisVal);
                    po::notify(p_OptionsDescriptor.optionValues.m_VM);
                }   break;

                case TYPENAME::FLOAT: {                                 // FLOAT
                    double thisVal = std::stof(optionValue);                    
                    p_OptionsDescriptor.optionValues.ModifyVariableMap(p_OptionsDescriptor.optionValues.m_VM, optionName, thisVal);
                    po::notify(p_OptionsDescriptor.optionValues.m_VM);
                }   break;

                case TYPENAME::DOUBLE: {                                // DOUBLE
                    double thisVal = std::stod(optionValue);                    
                    p_OptionsDescriptor.optionValues.ModifyVariableMap(p_OptionsDescriptor.optionValues.m_VM, optionName, thisVal);
                    po::notify(p_OptionsDescriptor.optionValues.m_VM);
                }   break;

                case TYPENAME::LONGDOUBLE: {                            // LONG DOUBLE
                    long double thisVal = std::stold(optionValue);                    
                    p_OptionsDescriptor.optionValues.ModifyVariableMap(p_OptionsDescriptor.optionValues.m_VM, optionName, thisVal);
                    po::notify(p_OptionsDescriptor.optionValues.m_VM);
                }   break;

                case TYPENAME::BOOL: {                                  // BOOL
                    bool thisVal = utils::IsBOOL(optionValue) > 0;      // already checked to be valid - just need to know if true or false
                    p_OptionsDescriptor.optionValues.ModifyVariableMap(p_OptionsDescriptor.optionValues.m_VM, optionName, thisVal);
                    po::notify(p_OptionsDescriptor.optionValues.m_VM);
                }   break;

                case TYPENAME::STRING: {                                // STRING
                    p_OptionsDescriptor.optionValues.ModifyVariableMap(p_OptionsDescriptor.optionValues.m_VM, optionName, optionValue);
                    po::notify(p_OptionsDescriptor.optionValues.m_VM);
                }   break;
                                
                default: break;                                         // already checked before we get here
            }            
        }
        else {                                                          // RANGE

            if (details.currPos >= details.rangeParms[1].iVal) {        // currPos is range count - wrap?
                if (idx == 0) {                                         // outermost iterator?
                    retVal = 0;                                         // yes - we're done
                    break;
                }
                else {                                                  // no - wrap
                    details.currPos = 0;                                // ... back to zero
                    stop = false;                                       // ... and don't stop at this iterator
                }
            }

            switch (details.dataType) {                                 // which data type?
                // these should cover our options for now - but we might have to refine them over time

                case TYPENAME::INT: {                                   // INT
                    int start = details.rangeParms[0].iVal;
                    int inc   = details.rangeParms[2].iVal;
                    int thisVal = start + (details.currPos * inc);
                    
                    p_OptionsDescriptor.optionValues.ModifyVariableMap(p_OptionsDescriptor.optionValues.m_VM, optionName, thisVal);
                    po::notify(p_OptionsDescriptor.optionValues.m_VM);
                }   break;

                case TYPENAME::LONGINT: {                               // LONG INT
                    long int start   = details.rangeParms[0].iVal;
                    long int inc     = details.rangeParms[2].iVal;
                    long int thisVal = start + (details.currPos * inc);
                    
                    p_OptionsDescriptor.optionValues.ModifyVariableMap(p_OptionsDescriptor.optionValues.m_VM, optionName, thisVal);
                    po::notify(p_OptionsDescriptor.optionValues.m_VM);
                  }   break;

                case TYPENAME::ULONGINT: {                              // UNSIGNED LONG INT
                    unsigned long int start   = details.rangeParms[0].iVal;
                    unsigned long int inc     = details.rangeParms[2].iVal;
                    unsigned long int thisVal = start + (details.currPos * inc);
                    
                    p_OptionsDescriptor.optionValues.ModifyVariableMap(p_OptionsDescriptor.optionValues.m_VM, optionName, thisVal);
                    po::notify(p_OptionsDescriptor.optionValues.m_VM);
                }   break;

                case TYPENAME::FLOAT: {                                 // FLOAT

                    double start   = details.rangeParms[0].dVal;
                    double inc     = details.rangeParms[2].dVal;
                    double thisVal = start + (details.currPos * inc);
                    
                    p_OptionsDescriptor.optionValues.ModifyVariableMap(p_OptionsDescriptor.optionValues.m_VM, optionName, thisVal);
                    po::notify(p_OptionsDescriptor.optionValues.m_VM);
                }   break;

                case TYPENAME::DOUBLE: {                                // DOUBLE

                    double start   = details.rangeParms[0].dVal;
                    double inc     = details.rangeParms[2].dVal;
                    double thisVal = start + (details.currPos * inc);
                    
                    p_OptionsDescriptor.optionValues.ModifyVariableMap(p_OptionsDescriptor.optionValues.m_VM, optionName, thisVal);
                    po::notify(p_OptionsDescriptor.optionValues.m_VM);
                }   break;

                case TYPENAME::LONGDOUBLE: {                            // LONG DOUBLE

                    long double start   = details.rangeParms[0].dVal;
                    long double inc     = details.rangeParms[2].dVal;
                    long double thisVal = start + (details.currPos * inc);
                    
                    p_OptionsDescriptor.optionValues.ModifyVariableMap(p_OptionsDescriptor.optionValues.m_VM, optionName, thisVal);
                    po::notify(p_OptionsDescriptor.optionValues.m_VM);
                }   break;

                default: break;                                         // already checked before we get here
            }
        }

        p_OptionsDescriptor.complexOptionValues[idx] = std::make_tuple(optionName, details); // reset values

        retVal = 1;                                                     // set return value
        idx--;                                                          // next (outer) iterator
    }

    std::string errStr = p_OptionsDescriptor.optionValues.CheckAndSetOptions();
    if (!errStr.empty()) std::cerr << ERR_MSG(ERROR::PROGRAM_OPTIONS_ERROR) << ": " << errStr << std::endl; // show the problem

    retVal = errStr.empty() ? retVal : -1;

    return retVal;
}


/*
 * Initialise grid file options
 * 
 * Initialises and populates the grid line (evolving object-level) options object, using the
 * options the user specified in the grid file record - this object is updated for each grid
 * line read.
 *
 * 
 * bool Options::InitialiseEvolvingObject(const std::string p_OptionsString)
 * 
 * @param   [IN]    p_OptionsString             String containing all options - the grid file record
 * @return                                      Boolean value indicating status: true = ok, false = an error occurred
 */
bool Options::InitialiseEvolvingObject(const std::string p_OptionsString) {

    bool ok = true;                                                                                                 // status - unless something changes

    try {

        m_GridLine.optionValues.Initialise();                                                                       // initialise option variables for evolving object-level options

        // parse the option string (just as the OS/shell would do)

        std::vector<std::string> parsedStrings;                                                                     // parsed option strings

        size_t start      = 0;                                                                                      // start position of parsed option string
        size_t end        = 0;                                                                                      // end position of parsed option strinf
        std::string delim = " ";                                                                                    // delimiter
        bool done         = false;
        while (!done && end != std::string::npos) {                                                                 // iterate over input string
            end = p_OptionsString.find(delim, start);                                                               // find delimiter
            std::string str = p_OptionsString.substr(start, end - start);                                           // grab option/argument string
            std::string trimmedStr = utils::trim(str);                                                              // trim whitespace
            if (trimmedStr[0] == '#') {                                                                             // comment?
                done = true;                                                                                        // yes - done with this line 
            }
            else {                                                                                                  // no - not a comment
                if (!trimmedStr.empty()) parsedStrings.push_back(trimmedStr);                                       // store if not empty string
                start = end + delim.length();                                                                       // new start position
            }
        }
    
        std::vector<char const *> args {"placeHolder"};                                                             // place-holder - boost expects command name as argv[0]

        for (auto& arg : parsedStrings)                                                                             // iterate over the parsed strings
            args.push_back(arg.c_str());                                                                            // and grab the c_str  

        // check here for excluded grid file options
        // if any are found, issue a warning and remove the option (and any
        // associated value) from the option strings
        // this should work here - I should be able to figure out what are
        // option names and what are option values - if it doesn't pan out
        // then we may need to move it to after Boost has parsed the options

        if (args.size() > 1) {                                                                                      // args present (other than the executable/placeholder)
            std::vector<int> removeArgs = {};                                                                       // vector of argument indices to be removed
            size_t iArg = 1;                                                                                        // start after the executable name/placeholder
            while (iArg < args.size()) {                                                                            // for each arg string (except the executable/placeholder)

                std::string optionName(args[iArg]);                                                                 // get the string (we'll call it the option name for now)

                // check whether the string really is an option name
                // we assume any string starting with "--" is an option name, and
                // any string starting with "-*", where '*' is an alphabetic
                // character, is an option name
                bool haveOptionName = false;                                                                        // is the string really an option name - default false
                if (optionName.length() > 1 && optionName[0] == '-') {                                              // check for '--' or '-*' (* is alpha character)
                    if (optionName[1] == '-') haveOptionName = true;                                                // '--' - option name
                    else if (isalpha(optionName[1])) haveOptionName = true;                                         // '-*', where * is alpha character - option name
                }
                iArg++;                                                                                             // next argument string

                if (haveOptionName) {                                                                               // do we think we have an option name?
                                                                                                                    // yes
                    if (optionName[0] == '-') optionName.erase(0, optionName.find_first_not_of("-"));               // remove the "-" or "--"

                    if (std::find(m_GridLineExcluded.begin(), m_GridLineExcluded.end(), optionName) != m_GridLineExcluded.end()) {  // on excluded list?
                        
                        removeArgs.push_back(iArg - 1);                                                             // yes - we need to remove it and any associated values

                        std::cerr << "WARNING: " << ERR_MSG(ERROR::OPTION_NOT_SUPPORTED_IN_GRID_FILE) << ": '" << optionName << "': ignored\n";  // show warning

                        // remove all option values for option to be removed

                        bool done = false;
                        while (!done) {

                            // we need to determine if the next argument string is a value for the
                            // option name we have, or whether it is the next option name - not all
                            // options need to specify a value (e.g. boolean switches)

                            std::string optionValue = "";
                            if (iArg < args.size()) optionValue = std::string(args[iArg]);                          // get the (potential) option value string

                            // check whether the string really is an option value
                            // as noted above, we assume any string starting with "--" is an option name,
                            // and any string starting with "-*", where '*' is an alphabetic character, 
                            // is an option name
                            bool haveOptionValue = false;                                                           // is the string really an option value - default false
                            if (!optionValue.empty()) {                                                             // empty string?
                                haveOptionValue = true;                                                             // no - we'll assume an option value, unless...
                                if (optionValue.length() > 1 && optionValue[0] == '-') {                            // check for '--' or '-*' (* is alpha character)
                                    if (optionValue[1] == '-') haveOptionValue = false;                             // '--' - option name, not value
                                    else if (isalpha(optionValue[1])) haveOptionValue = false;                      // '-*', where * is alpha character - option name, not value
                                }
                            }

                            if (haveOptionValue) {                                                                  // do we think we have an option value?
                                removeArgs.push_back(iArg);                                                         // yes - we need to remove it
                                iArg++;                                                                             // next argument string
                            }
                            else done = true;
                        }
                    }
                }
            }

            // Remove any argument strings identified as being excluded options
            // and the value associated with excluded options.
            // This only works because the removeArgs vector is populated in the
            // same order as the args vector - that allows me to loop through args
            // in reverse order to delete elements (otherwise the index numbers
            // would change from under me as I was deleting them)

            size_t removeCount(removeArgs.size());
            if (removeCount > 0) {                                                                                  // anything to remove?
                for (size_t iArg = removeCount; iArg > 0; iArg--) {                                                 // loop through all args identified to be removed
                    args.erase(args.begin() + removeArgs[iArg - 1]);                                                // erase it
                }
            }
        }

        // parse the option values before handing them over to boost
        // boost knows nothing about shorthand, ranges and sets, so we have to handle
        // them ourselves first
        m_GridLine.complexOptionValues = {};                                                                        // no ranges or sets - unless we find them in the parse
                                                    /*****  << do *not* try this at home >> *****/
        std::string errStr = ParseOptionValues(args.size(), const_cast<char**>(args.data()), m_GridLine);           // parse the option values - specifically for ranges and sets
        if (!errStr.empty()) {                                                                                      // parsed ok?
            COMPLAIN(errStr);                                                                                       // no, complain - this throws an exception
        }
        else {
            errStr = m_GridLine.optionValues.CheckAndSetOptions();                                                  // yes - sanity check, and set, evolving object-level values
            if (!errStr.empty()) {                                                                                  // check ok?
                COMPLAIN(errStr);                                                                                   // no, complain - this throws an exception
            }
       }
    }
    catch (po::error& e) {                                                                                          // program options exception
        std::cerr << ERR_MSG(ERROR::GRID_OPTIONS_ERROR) << ": " << e.what() << std::endl;                           // show the problem
        std::cerr << ERR_MSG(ERROR::SUGGEST_HELP) << std::endl;                                                     // suggest using --help
        ok = false;                                                                                                 // set status
    } 
    catch (const std::string eStr) {                                                                                // custom exception
        std::cerr << ERR_MSG(ERROR::GRID_OPTIONS_ERROR) << ": " << eStr << std::endl;                               // show the problem
        std::cerr << ERR_MSG(ERROR::SUGGEST_HELP) << std::endl;                                                     // suggest using --help
        ok = false;                                                                                                 // set status
    }
    catch (...) {                                                                                                   // unhandled exception
        std::cerr << ERR_MSG(ERROR::GRID_OPTIONS_ERROR) << ": " << ERR_MSG(ERROR::UNHANDLED_EXCEPTION) << std::endl; // show the problem
        std::cerr << ERR_MSG(ERROR::SUGGEST_HELP) << std::endl;                                                     // suggest using --help
        ok = false;                                                                                                 // set status
    }

    m_GridLine.optionValues.m_Populated = ok;                                                                       // flag use
    
    return ok;
}


/*
 * Read and apply the next record in the grid file
 * 
 * The record from the grid file is read as one string, then passed to
 * InitialiseEvolvingObject() for processing.
 * 
 * In InitialiseEvolvingObject() the record is parsed into separate tokens
 * ready to be passed to the boost program option parser.  Once that is 
 * done the options are handed over to the boost functions for parsing and
 * detting of values.
 * 
 * pon return from this function the option values will be set to the values 
 * specified by the user, or their default values - either way ready for the 
 * star/binary to be evolved.
 * 
 * The grid file struct (m_Gridfile) will be used and updated by this function.
 * The grid file name, error status, and file handle are stored in the struct.  
 * 
 * 
 * int ApplyNextGridLine()
 * 
 * @return                                      Int result:
 *                                                  -1: Error reading grid file record (error value in grid file struct)
 *                                                   0: No record to read - end of file
 *                                                   1: Grid file record read and applied ok
 */
int Options::ApplyNextGridLine() {

    int status = -1;                                                                                // default status is failure

    if (m_Gridfile.handle.is_open()) {                                                              // file open?

        if (m_Gridfile.linesProcessed >= m_Gridfile.linesToProcess) {                               // yes - have we already processed all the records the user wants processed?
            status = 0;                                                                             // yes - return EOF so processing of grid file stops
        }
        else {                                                                                      // no - process current record
            bool done = false;
            while (!done) {
                std::string record;                                                                 // the record read
                std::getline(m_Gridfile.handle, record);                                            // read the next record
                if (m_Gridfile.handle.fail()) {                                                     // read ok?
                    if (m_Gridfile.handle.eof()) {                                                  // no - eof?
                                                                                                    // yes - eof
                        if (OPTIONS->OptionSpecified("grid-lines-to-process") == 1 &&               // user specified number of grid lines to process?
                            m_Gridfile.linesProcessed < m_Gridfile.linesToProcess) {                // yes - did we process all the lines the user asked for?
                            m_Gridfile.error = ERROR::UNEXPECTED_END_OF_FILE;                       // no - not all lines user asked for were processed before EOF - record error
                            status = -1;                                                            // set error status
                        }
                        else status = 0;                                                            // set EOF status         
                    }
                    else {                                                                          // not eof - some other error
                        m_Gridfile.error = ERROR::FILE_READ_ERROR;                                  // record error
                        status = -1;                                                                // set error status
                    }
                    done = true;                                                                    // we're done
                }
                else {                                                                              // read ok
                    m_Gridfile.currentLine++;                                                       // increment line about to be processed (will be current)
                    m_Gridfile.linesProcessed++;                                                    // increment lines processed
                    record = utils::ltrim(record);                                                  // trim leading white space
                    if (!record.empty() && record[0] != '#') {                                      // blank line or comment?
                        status = InitialiseEvolvingObject(record) ? 1 : -1;                         // no - apply record and set status
                        done = true;                                                                // we're done
                    }
                }
            }
        }
    }

    return status;
}


/*
 * Seek to grid file line
 *
 * The grid file is a variable-length file (in that each line is not a fixed number of bytes),
 * so we can't use the filesystem seek functions - they rely on each line being the same number of bytes.
 * We have to just read each line until we get to the one we want - use ignore() instead of getline() 
 * because it's a little faster.  There will be some overhead in doing this, especially for very large 
 * grid file, but it shouldn't be significant, especially compared to the overall runtime.
 * (I tested this on a grid file of 1,000,000 lines, each line 85 bytes - to process just the first line of
 * the file was 0.03 CPU seconds, and to process just the last line of the file (skipping the first 999,999
 * lines) was 0.07 CPU seconds - suggesting that even skipping to the end of a grid file of several million 
 * records might only add overhead of one or two tenths of a second of CPU time to the entire run)
 * 
 * 
 * ERROR SeekToGridFileLine(const unsigned int p_Line)
 *
 * @param   [IN]        p_Line                  The line number to seek to - this will be the next line read from the file
 * @return                                      ERROR indicator - will be ERROR::NONE if seek is successful
 */
ERROR Options::SeekToGridFileLine(const unsigned int p_Line) { 

    m_Gridfile.error = ERROR::NONE;                                                         // default is no error

    if (m_Gridfile.startLine > 0) {                                                         // need to seek?
        for (unsigned int line = 0; line < m_Gridfile.startLine; line++) {                  // yes - for each line until start line
            m_Gridfile.handle.ignore(std::numeric_limits<std::streamsize>::max(), '\n');    // skip to end of record
            if (m_Gridfile.handle.fail()) {                                                 // skip ok?
                if (m_Gridfile.handle.eof()) {                                              // no - eof?
                    m_Gridfile.error = ERROR::UNEXPECTED_END_OF_FILE;                       // yes - record error
                }
                else {                                                                      // not eof
                    m_Gridfile.error = ERROR::FILE_READ_ERROR;                              // record error
                    break;                                                                  // we're done
                }
            }
            else {                                                                          // skip ok
                m_Gridfile.currentLine++;                                                   // set line about to be processed (will be current)
            }
        } 
    }
 
    return m_Gridfile.error;    
}


/*
 * Rewind the grid file
 *
 * The grid file is rewound to prepare for the next commandline options variation.
 * Here we have to seek to the start of the file, then advance to the first line 
 * the user asked to be processed.
 *
 * 
 * ERROR RewindGridFile()
 *
 * @return                                      ERROR indicator - will be ERROR::NONE if file opened successfully
 */
ERROR Options::RewindGridFile() { 
    
    m_Gridfile.handle.clear();                          // clear file errors
    m_Gridfile.handle.seekg(0);                         // go to start of file

    m_Gridfile.error          = ERROR::NONE;            // reset error (none)
    m_Gridfile.currentLine    = 0;                      // reset current lines
    m_Gridfile.linesProcessed = 0;                      // reset number of lines processed

    return SeekToGridFileLine(m_Gridfile.startLine);    // seek to start line requested by user and return error
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
 * @return                                      ERROR indicator - will be ERROR::NONE if file opened successfully
 */
ERROR Options::OpenGridFile(const std::string p_GridFilename) {

    m_Gridfile.filename = p_GridFilename;                                               // record filename

    if (!m_Gridfile.filename.empty()) {                                                 // have grid filename?
        m_Gridfile.handle.open(m_Gridfile.filename);                                    // yes - open the file
        if (m_Gridfile.handle.fail()) {                                                 // open ok?
            m_Gridfile.error = ERROR::FILE_OPEN_ERROR;                                  // no - record error
        }
        else {                                                                          // open ok - no error
            m_Gridfile.error = ERROR::NONE;                                             // record success

            m_Gridfile.startLine      = OPTIONS->GridStartLine();                       // set first line to process (0-based)
            m_Gridfile.linesToProcess = OPTIONS->GridLinesToProcess();                  // set number of lines to process (-1 = process to EOF)
            m_Gridfile.currentLine    = 0;                                              // set line about to be processed (will be current)
            m_Gridfile.linesProcessed = 0;                                              // set number of lines processed in this run

            m_Gridfile.error = SeekToGridFileLine(m_Gridfile.startLine);                // Seek to first line to be processed (if necessary)
        }
    }
    else m_Gridfile.error = ERROR::EMPTY_FILENAME;                                      // empty filename

    return m_Gridfile.error;
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

        case PROGRAM_OPTION::ADD_OPTIONS_TO_SYSPARMS                        : value = static_cast<int>(AddOptionsToSysParms());                             break;
        case PROGRAM_OPTION::ALLOW_NON_STRIPPED_ECSN                        : value = AllowNonStrippedECSN();                                               break;
        case PROGRAM_OPTION::ALLOW_MS_STAR_TO_SURVIVE_COMMON_ENVELOPE       : value = AllowMainSequenceStarToSurviveCommonEnvelope();                       break;
        case PROGRAM_OPTION::ALLOW_RADIATIVE_ENVELOPE_STAR_TO_SURVIVE_COMMON_ENVELOPE : value = AllowRadiativeEnvelopeStarToSurviveCommonEnvelope();        break;
        case PROGRAM_OPTION::ALLOW_IMMEDIATE_RLOF_POST_CE_TO_SURVIVE_COMMON_ENVELOPE  : value = AllowImmediateRLOFpostCEToSurviveCommonEnvelope();          break;
        case PROGRAM_OPTION::ALLOW_RLOF_AT_BIRTH                            : value = AllowRLOFAtBirth();                                                   break;
        case PROGRAM_OPTION::ALLOW_TOUCHING_AT_BIRTH                        : value = AllowTouchingAtBirth();                                               break;
        case PROGRAM_OPTION::ANG_MOM_CONSERVATION_DURING_CIRCULARISATION    : value = AngularMomentumConservationDuringCircularisation();                   break;

        //case PROGRAM_OPTION::BE_BINARIES                                    : value = BeBinaries();                                                         break;

        case PROGRAM_OPTION::BLACK_HOLE_KICKS                               : value = static_cast<int>(BlackHoleKicks());                                   break;
    
        case PROGRAM_OPTION::CASE_BB_STABILITY_PRESCRIPTION                 : value = static_cast<int>(CaseBBStabilityPrescription());                      break;
    
        case PROGRAM_OPTION::CHECK_PHOTON_TIRING_LIMIT                      : value = CheckPhotonTiringLimit();                                             break;

        case PROGRAM_OPTION::CHE_MODE                                       : value = static_cast<int>(CHEMode());                                          break;

        case PROGRAM_OPTION::CIRCULARISE_BINARY_DURING_MT                   : value = CirculariseBinaryDuringMassTransfer();                                break;

        case PROGRAM_OPTION::COMMON_ENVELOPE_ALPHA                          : value = CommonEnvelopeAlpha();                                                break;
        case PROGRAM_OPTION::COMMON_ENVELOPE_ALPHA_THERMAL                  : value = CommonEnvelopeAlphaThermal();                                         break;
        case PROGRAM_OPTION::COMMON_ENVELOPE_LAMBDA                         : value = CommonEnvelopeLambda();                                               break;
        case PROGRAM_OPTION::COMMON_ENVELOPE_LAMBDA_MULTIPLIER              : value = CommonEnvelopeLambdaMultiplier();                                     break;
        case PROGRAM_OPTION::COMMON_ENVELOPE_LAMBDA_PRESCRIPTION            : value = static_cast<int>(CommonEnvelopeLambdaPrescription());                 break;
        case PROGRAM_OPTION::COMMON_ENVELOPE_MASS_ACCRETION_CONSTANT        : value = CommonEnvelopeMassAccretionConstant();                                break;
        case PROGRAM_OPTION::COMMON_ENVELOPE_MASS_ACCRETION_MAX             : value = CommonEnvelopeMassAccretionMax();                                     break;
        case PROGRAM_OPTION::COMMON_ENVELOPE_MASS_ACCRETION_MIN             : value = CommonEnvelopeMassAccretionMin();                                     break;
        case PROGRAM_OPTION::COMMON_ENVELOPE_MASS_ACCRETION_PRESCRIPTION    : value = static_cast<int>(CommonEnvelopeMassAccretionPrescription());          break;
        case PROGRAM_OPTION::COMMON_ENVELOPE_RECOMBINATION_ENERGY_DENSITY   : value = CommonEnvelopeRecombinationEnergyDensity();                           break;
        case PROGRAM_OPTION::COMMON_ENVELOPE_SLOPE_KRUCKOW                  : value = CommonEnvelopeSlopeKruckow();                                         break;

        case PROGRAM_OPTION::CONVECTIVE_ENVELOPE_TEMPERATURE_THRESHOLD      : value = ConvectiveEnvelopeTemperatureThreshold();                             break;

        case PROGRAM_OPTION::COOL_WIND_MASS_LOSS_MULTIPLIER                 : value = CoolWindMassLossMultiplier();                                         break;

        case PROGRAM_OPTION::ECCENTRICITY                                   : value = Eccentricity();                                                       break;
        case PROGRAM_OPTION::ECCENTRICITY_DISTRIBUTION                      : value = static_cast<int>(EccentricityDistribution());                         break;
        case PROGRAM_OPTION::ECCENTRICITY_DISTRIBUTION_MAX                  : value = EccentricityDistributionMax();                                        break;
        case PROGRAM_OPTION::ECCENTRICITY_DISTRIBUTION_MIN                  : value = EccentricityDistributionMin();                                        break;
        case PROGRAM_OPTION::EDDINGTON_ACCRETION_FACTOR                     : value = EddingtonAccretionFactor();                                           break;
        case PROGRAM_OPTION::ENVELOPE_STATE_PRESCRIPTION                    : value = static_cast<int>(EnvelopeStatePrescription());                        break;
        case PROGRAM_OPTION::EVOLUTION_MODE                                 : value = static_cast<int>(EvolutionMode());                                    break;

        case PROGRAM_OPTION::FRYER_SUPERNOVA_ENGINE                         : value = static_cast<int>(FryerSupernovaEngine());                             break;

        case PROGRAM_OPTION::FRYER22_FMIX                                   : value = Fryer22fmix();                                       break;
        case PROGRAM_OPTION::FRYER22_MCRIT                                  : value = Fryer22Mcrit();                                       break;

        case PROGRAM_OPTION::INITIAL_MASS                                   : value = InitialMass();                                                        break;
        case PROGRAM_OPTION::INITIAL_MASS_1                                 : value = InitialMass1();                                                       break;
        case PROGRAM_OPTION::INITIAL_MASS_2                                 : value = InitialMass2();                                                       break;

        case PROGRAM_OPTION::INITIAL_MASS_FUNCTION                          : value = static_cast<int>(InitialMassFunction());                              break;
        case PROGRAM_OPTION::INITIAL_MASS_FUNCTION_MAX                      : value = InitialMassFunctionMax();                                             break;
        case PROGRAM_OPTION::INITIAL_MASS_FUNCTION_MIN                      : value = InitialMassFunctionMin();                                             break;
        case PROGRAM_OPTION::INITIAL_MASS_FUNCTIONPOWER                     : value = InitialMassFunctionPower();                                           break;

        case PROGRAM_OPTION::KICK_DIRECTION_DISTRIBUTION                    : value = static_cast<int>(KickDirectionDistribution());                        break;
        case PROGRAM_OPTION::KICK_DIRECTION_POWER                           : value = KickDirectionPower();                                                 break;
        case PROGRAM_OPTION::KICK_SCALING_FACTOR                            : value = KickScalingFactor();                                                  break;
        case PROGRAM_OPTION::KICK_MAGNITUDE_DISTRIBUTION                    : value = static_cast<int>(KickMagnitudeDistribution());                        break;
        case PROGRAM_OPTION::KICK_MAGNITUDE_DISTRIBUTION_MAXIMUM            : value = KickMagnitudeDistributionMaximum();                                   break;

        case PROGRAM_OPTION::KICK_MAGNITUDE_DISTRIBUTION_SIGMA_CCSN_BH      : value = KickMagnitudeDistributionSigmaCCSN_BH();                              break;
        case PROGRAM_OPTION::KICK_MAGNITUDE_DISTRIBUTION_SIGMA_CCSN_NS      : value = KickMagnitudeDistributionSigmaCCSN_NS();                              break;
        case PROGRAM_OPTION::KICK_MAGNITUDE_DISTRIBUTION_SIGMA_FOR_ECSN     : value = KickMagnitudeDistributionSigmaForECSN();                              break;
        case PROGRAM_OPTION::KICK_MAGNITUDE_DISTRIBUTION_SIGMA_FOR_USSN     : value = KickMagnitudeDistributionSigmaForUSSN();                              break;

        case PROGRAM_OPTION::KICK_MAGNITUDE                                 : value = KickMagnitude();                                                      break;
        case PROGRAM_OPTION::KICK_MAGNITUDE_1                               : value = KickMagnitude1();                                                     break;
        case PROGRAM_OPTION::KICK_MAGNITUDE_2                               : value = KickMagnitude2();                                                     break;

        case PROGRAM_OPTION::KICK_MAGNITUDE_RANDOM                          : value = KickMagnitudeRandom();                                                break;
        case PROGRAM_OPTION::KICK_MAGNITUDE_RANDOM_1                        : value = KickMagnitudeRandom1();                                               break;
        case PROGRAM_OPTION::KICK_MAGNITUDE_RANDOM_2                        : value = KickMagnitudeRandom2();                                               break;

        case PROGRAM_OPTION::KICK_MEAN_ANOMALY_1                            : value = SN_MeanAnomaly1();                                                    break;
        case PROGRAM_OPTION::KICK_MEAN_ANOMALY_2                            : value = SN_MeanAnomaly2();                                                    break;
        case PROGRAM_OPTION::KICK_PHI_1                                     : value = SN_Phi1();                                                            break;
        case PROGRAM_OPTION::KICK_PHI_2                                     : value = SN_Phi2();                                                            break;
        case PROGRAM_OPTION::KICK_THETA_1                                   : value = SN_Theta1();                                                          break;
        case PROGRAM_OPTION::KICK_THETA_2                                   : value = SN_Theta2();                                                          break;

        case PROGRAM_OPTION::LBV_FACTOR                                     : value = LuminousBlueVariableFactor();                                         break;
        case PROGRAM_OPTION::LBV_PRESCRIPTION                               : value = static_cast<int>(LuminousBlueVariablePrescription());                 break;

        case PROGRAM_OPTION::MASS_LOSS_PRESCRIPTION                         : value = static_cast<int>(MassLossPrescription());                             break;

        case PROGRAM_OPTION::MASS_RATIO                                     : value = MassRatio();                                                          break;                     
        case PROGRAM_OPTION::MASS_RATIO_DISTRIBUTION                        : value = static_cast<int>(MassRatioDistribution());                            break;
        case PROGRAM_OPTION::MASS_RATIO_DISTRIBUTION_MAX                    : value = MassRatioDistributionMax();                                           break;
        case PROGRAM_OPTION::MASS_RATIO_DISTRIBUTION_MIN                    : value = MassRatioDistributionMin();                                           break;

        case PROGRAM_OPTION::MAXIMUM_EVOLUTION_TIME                         : value = MaxEvolutionTime();                                                   break;
        case PROGRAM_OPTION::MAXIMUM_DONOR_MASS                             : value = MaximumDonorMass();                                                   break;
        case PROGRAM_OPTION::MAXIMUM_NEUTRON_STAR_MASS                      : value = MaximumNeutronStarMass();                                             break;
        case PROGRAM_OPTION::MAXIMUM_TIMESTEPS                              : value = MaxNumberOfTimestepIterations();                                      break;

        case PROGRAM_OPTION::MCBUR1                                         : value = MCBUR1();                                                             break;

        case PROGRAM_OPTION::METALLICITY                                    : value = Metallicity();                                                        break;
        case PROGRAM_OPTION::METALLICITY_DISTRIBUTION                       : value = static_cast<int>(MetallicityDistribution());                          break;
        case PROGRAM_OPTION::METALLICITY_DISTRIBUTION_MAX                   : value = MetallicityDistributionMax();                                         break;
        case PROGRAM_OPTION::METALLICITY_DISTRIBUTION_MIN                   : value = MetallicityDistributionMin();                                         break;

        case PROGRAM_OPTION::MINIMUM_MASS_SECONDARY                         : value = MinimumMassSecondary();                                               break;

        case PROGRAM_OPTION::MT_ACCRETION_EFFICIENCY_PRESCRIPTION           : value = static_cast<int>(MassTransferAccretionEfficiencyPrescription());      break;
        case PROGRAM_OPTION::MT_ANG_MOM_LOSS_PRESCRIPTION                   : value = static_cast<int>(MassTransferAngularMomentumLossPrescription());      break;
        case PROGRAM_OPTION::MT_THERMAL_LIMIT_C                             : value = MassTransferCParameter();                                             break;

        case PROGRAM_OPTION::MT_CRIT_MR_MS_LOW_MASS_DEGENERATE_ACCRETOR     : value = MassTransferCriticalMassRatioMSLowMassDegenerateAccretor();           break;
        case PROGRAM_OPTION::MT_CRIT_MR_MS_LOW_MASS_NON_DEGENERATE_ACCRETOR : value = MassTransferCriticalMassRatioMSLowMassNonDegenerateAccretor();        break;
        case PROGRAM_OPTION::MT_CRIT_MR_MS_HIGH_MASS_DEGENERATE_ACCRETOR    : value = MassTransferCriticalMassRatioMSHighMassDegenerateAccretor();          break;
        case PROGRAM_OPTION::MT_CRIT_MR_MS_HIGH_MASS_NON_DEGENERATE_ACCRETOR: value = MassTransferCriticalMassRatioMSHighMassNonDegenerateAccretor();       break;
        case PROGRAM_OPTION::MT_CRIT_MR_GIANT_DEGENERATE_ACCRETOR           : value = MassTransferCriticalMassRatioGiantDegenerateAccretor();               break;
        case PROGRAM_OPTION::MT_CRIT_MR_GIANT_NON_DEGENERATE_ACCRETOR       : value = MassTransferCriticalMassRatioGiantNonDegenerateAccretor();            break;
        case PROGRAM_OPTION::MT_CRIT_MR_HG_DEGENERATE_ACCRETOR              : value = MassTransferCriticalMassRatioHGDegenerateAccretor();                  break;
        case PROGRAM_OPTION::MT_CRIT_MR_HG_NON_DEGENERATE_ACCRETOR          : value = MassTransferCriticalMassRatioHGNonDegenerateAccretor();               break;
        case PROGRAM_OPTION::MT_CRIT_MR_HE_GIANT_DEGENERATE_ACCRETOR        : value = MassTransferCriticalMassRatioHeliumGiantDegenerateAccretor();         break;
        case PROGRAM_OPTION::MT_CRIT_MR_HE_GIANT_NON_DEGENERATE_ACCRETOR    : value = MassTransferCriticalMassRatioHeliumGiantNonDegenerateAccretor();      break;
        case PROGRAM_OPTION::MT_CRIT_MR_HE_HG_DEGENERATE_ACCRETOR           : value = MassTransferCriticalMassRatioHeliumHGDegenerateAccretor();            break;
        case PROGRAM_OPTION::MT_CRIT_MR_HE_HG_NON_DEGENERATE_ACCRETOR       : value = MassTransferCriticalMassRatioHeliumHGNonDegenerateAccretor();         break;
        case PROGRAM_OPTION::MT_CRIT_MR_HE_MS_DEGENERATE_ACCRETOR           : value = MassTransferCriticalMassRatioHeliumMSDegenerateAccretor();            break;
        case PROGRAM_OPTION::MT_CRIT_MR_HE_MS_NON_DEGENERATE_ACCRETOR       : value = MassTransferCriticalMassRatioHeliumMSNonDegenerateAccretor();         break;
        case PROGRAM_OPTION::MT_CRIT_MR_WD_DEGENERATE_ACCRETOR              : value = MassTransferCriticalMassRatioWhiteDwarfDegenerateAccretor();          break;
        case PROGRAM_OPTION::MT_CRIT_MR_WD_NONDEGENERATE_ACCRETOR           : value = MassTransferCriticalMassRatioWhiteDwarfNonDegenerateAccretor();       break;

        case PROGRAM_OPTION::MT_FRACTION_ACCRETED                           : value = MassTransferFractionAccreted();                                       break;
        case PROGRAM_OPTION::MT_JLOSS                                       : value = MassTransferJloss();                                                  break;
        case PROGRAM_OPTION::MT_JLOSS_MACLEOD_LINEAR_FRACTION               : value = MassTransferJlossMacLeodLinearFraction();                             break; 
        case PROGRAM_OPTION::MT_REJUVENATION_PRESCRIPTION                   : value = static_cast<int>(MassTransferRejuvenationPrescription());             break;
        case PROGRAM_OPTION::MT_THERMALLY_LIMITED_VARIATION                 : value = static_cast<int>(MassTransferThermallyLimitedVariation());            break;

        case PROGRAM_OPTION::MULLER_MANDEL_KICK_MULTIPLIER_BH               : value = MullerMandelKickMultiplierBH();                                       break;
        case PROGRAM_OPTION::MULLER_MANDEL_KICK_MULTIPLIER_NS               : value = MullerMandelKickMultiplierNS();                                       break;
        case PROGRAM_OPTION::MULLER_MANDEL_SIGMA_KICK                       : value = MullerMandelSigmaKick();                                             break;

        case PROGRAM_OPTION::NEUTRINO_MASS_LOSS_ASSUMPTION_BH               : value = static_cast<int>(NeutrinoMassLossAssumptionBH());                     break;
        case PROGRAM_OPTION::NEUTRINO_MASS_LOSS_VALUE_BH                    : value = NeutrinoMassLossValueBH();                                            break;

        case PROGRAM_OPTION::NOTES                                          : value = Notes();                                                              break;

        case PROGRAM_OPTION::NS_EOS                                         : value = static_cast<int>(NeutronStarEquationOfState());                       break;

        case PROGRAM_OPTION::ORBITAL_PERIOD                                 : value = OrbitalPeriod();                                                      break;
        case PROGRAM_OPTION::ORBITAL_PERIOD_DISTRIBUTION                    : value = static_cast<int>(OrbitalPeriodDistribution());                        break;
        case PROGRAM_OPTION::ORBITAL_PERIOD_DISTRIBUTION_MAX                : value = OrbitalPeriodDistributionMax();                                       break;
        case PROGRAM_OPTION::ORBITAL_PERIOD_DISTRIBUTION_MIN                : value = OrbitalPeriodDistributionMin();                                       break;

        case PROGRAM_OPTION::OVERALL_WIND_MASS_LOSS_MULTIPLIER              : value = OverallWindMassLossMultiplier();                                      break;

        case PROGRAM_OPTION::PISN_LOWER_LIMIT                               : value = PairInstabilityLowerLimit();                                          break;
        case PROGRAM_OPTION::PISN_UPPER_LIMIT                               : value = PairInstabilityUpperLimit();                                          break;

        case PROGRAM_OPTION::PULSAR_MAGNETIC_FIELD_DISTRIBUTION             : value = static_cast<int>(PulsarBirthMagneticFieldDistribution());             break;
        case PROGRAM_OPTION::PULSAR_MAGNETIC_FIELD_DISTRIBUTION_MAX         : value = PulsarBirthMagneticFieldDistributionMax();                            break;
        case PROGRAM_OPTION::PULSAR_MAGNETIC_FIELD_DISTRIBUTION_MIN         : value = PulsarBirthMagneticFieldDistributionMin();                            break;

        case PROGRAM_OPTION::PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION          : value = static_cast<int>(PulsarBirthSpinPeriodDistribution());                break;
        case PROGRAM_OPTION::PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION_MAX      : value = PulsarBirthSpinPeriodDistributionMax();                               break;
        case PROGRAM_OPTION::PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION_MIN      : value = PulsarBirthSpinPeriodDistributionMin();                               break;

        case PROGRAM_OPTION::PULSAR_MINIMUM_MAGNETIC_FIELD                  : value = PulsarLog10MinimumMagneticField();                                    break;

        case PROGRAM_OPTION::PULSAR_MAGNETIC_FIELD_DECAY_MASS_SCALE         : value = PulsarMagneticFieldDecayMassscale();                                  break;
        case PROGRAM_OPTION::PULSAR_MAGNETIC_FIELD_DECAY_TIME_SCALE         : value = PulsarMagneticFieldDecayTimescale();                                  break;

        case PROGRAM_OPTION::PPI_PRESCRIPTION                               : value = static_cast<int>(PulsationalPairInstabilityPrescription());           break;
        case PROGRAM_OPTION::PPI_LOWER_LIMIT                                : value = PulsationalPairInstabilityLowerLimit();                               break;
        case PROGRAM_OPTION::PPI_UPPER_LIMIT                                : value = PulsationalPairInstabilityUpperLimit();                               break;

        case PROGRAM_OPTION::QCRIT_PRESCRIPTION                             : value = static_cast<int>(QCritPrescription());                                break;

        case PROGRAM_OPTION::RANDOM_SEED                                    : value = RandomSeed();                                                         break;
        case PROGRAM_OPTION::RANDOM_SEED_CMDLINE                            : value = RandomSeedCmdLine();                                                  break;

        case PROGRAM_OPTION::REMNANT_MASS_PRESCRIPTION                      : value = static_cast<int>(RemnantMassPrescription());                          break;

        case PROGRAM_OPTION::ROTATIONAL_VELOCITY_DISTRIBUTION               : value = static_cast<int>(RotationalVelocityDistribution());                   break;

        case PROGRAM_OPTION::ROTATIONAL_FREQUENCY                           : value = RotationalFrequency();                                                break;
        case PROGRAM_OPTION::ROTATIONAL_FREQUENCY_1                         : value = RotationalFrequency1();                                               break;
        case PROGRAM_OPTION::ROTATIONAL_FREQUENCY_2                         : value = RotationalFrequency2();                                               break;
   
        case PROGRAM_OPTION::SEMI_MAJOR_AXIS                                : value = SemiMajorAxis();                                                      break;
        case PROGRAM_OPTION::SEMI_MAJOR_AXIS_DISTRIBUTION                   : value = static_cast<int>(SemiMajorAxisDistribution());                        break;
        case PROGRAM_OPTION::SEMI_MAJOR_AXIS_DISTRIBUTION_MAX               : value = SemiMajorAxisDistributionMax();                                       break;
        case PROGRAM_OPTION::SEMI_MAJOR_AXIS_DISTRIBUTION_MIN               : value = SemiMajorAxisDistributionMin();                                       break;
        case PROGRAM_OPTION::SEMI_MAJOR_AXIS_DISTRIBUTION_POWER             : value = SemiMajorAxisDistributionPower();                                     break;

        case PROGRAM_OPTION::STELLAR_ZETA_PRESCRIPTION                      : value = static_cast<int>(StellarZetaPrescription());                          break;

        case PROGRAM_OPTION::WR_FACTOR                                      : value = WolfRayetFactor();                                                    break;

        case PROGRAM_OPTION::ZETA_RADIATIVE_ENVELOPE_GIANT                  : value = ZetaRadiativeEnvelopeGiant();                                         break;
        case PROGRAM_OPTION::ZETA_MS                                        : value = ZetaMainSequence();                                                   break;
        case PROGRAM_OPTION::ZETA_ADIABATIC_ARBITRARY                       : value = ZetaAdiabaticArbitrary();                                             break;

        default:                                                                                                        // unknown property
            ok    = false;                                                                                              // that's not ok...
            value = "UNKNOWN";                                                                                          // default value
            std::cerr << ERR_MSG(ERROR::UNKNOWN_PROGRAM_OPTION) << std::endl;                                           // show warning (don't have logging or errors here...)
    }

    return std::make_tuple(ok, value);
}


/*
 * Sets the seed for the pseudo random number generator.
 * 
 * After setting the seed for the pseudo random number generator, recalculates the "calculated"
 * program option values - those that are calculated from other variables or drawn from
 * distributions.  This is required to ensure that all default option values are drawn or calculated
 * using the random seed for the current system (star or binary).  The parameter "p_OptionsSet"
 * indicates whether the command-line or gridfile-line set of option values is updated.
 * 
 * 
 * int AdvanceOptionVariation(SetRandomSeed(OptionsDescriptorT &p_OptionsDescriptor, OPTIONS_ORIGIN p_OptionSet)
 * 
 * @param   [IN]    p_RandomSeed                The random seed to use as the seed for the pseudo random number generator
 * @param   [IN]    p_OptionsSet                Indicates which set of options to update (command-line or gridfile-line)
 * @return                                      Int result:
 *                                                  -1: en error occurred
 *                                                   0: no more variations - all done
 *                                                   1: new variation applied - option values are set
 */

int Options::SetRandomSeed(const unsigned long int p_RandomSeed, const OPTIONS_ORIGIN p_OptionsSet) {

    RAND->Seed(p_RandomSeed);       // seed the pseudo random number generator

    // recalculate "calculated" option values for the relevant set of option values
    std::string err = p_OptionsSet == OPTIONS_ORIGIN::CMDLINE
                        ? m_CmdLine.optionValues.SetCalculatedOptionDefaults(BOOST_MAP::NO_UPDATE)
                        : m_GridLine.optionValues.SetCalculatedOptionDefaults(BOOST_MAP::NO_UPDATE);

    return err == "" ? 0 : -1;
}
