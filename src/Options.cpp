#include "Options.h"


Options* Options::m_Instance = nullptr;


Options* Options::Instance() {
    if (!m_Instance) {
        m_Instance = new Options();
    }
    return m_Instance;
}


COMMANDLINE_STATUS Options::Initialise(int argc, char *argv[]) {

    InitialiseMemberVariables();

    return CommandLineSorter(argc, argv);        // parse commandline program options
}


void Options::InitialiseMemberVariables(void) {

    // This sets all of the program options to their default values -- can be modified via the command line.

    // flags

    debugToFile                                                     = false;                                                                            // default is do not log debug statements to a log file
    errorsToFile                                                    = false;                                                                            // default is do not log error messages to a log file

    individualSystem                                                = false;                                                                            // Flag to evolve a specific individual system which you can specify initial parameters of
    singleStar                                                      = false;                                                                            // Flag to evolve a single star

	lambdaCalculationEveryTimeStep                                  = false;
	zetaCalculationEveryTimeStep                                    = false;

	beBinaries                                                      = false;
    evolvePulsars                                                   = false;                                                                            // Whether to evolve pulsars
	evolveUnboundSystems                                            = false;                                                                            // Allow unbound syetms to evolve
//    onlyDoubleCompactObjects                                        = false;                                                                            // Flag to turn on some shortcuts to only evolve systems which may form double compact objects

    detailedOutput                                                  = false;                                                                            // Detailed output
    populationDataPrinting                                          = false;                                                                            // Print certain data for small populations, but not for larger one
    printBoolAsString                                               = false;                                                                            // default is do not print bool as string
    quiet                                                           = false;                                                                            // Suppress some of the printing
    rlofPrinting                                                    = false;

    useImportanceSampling                                           = false;
//    useMCMC                                                         = false;

    nBatchesUsed                                                    = -1;                                                                               // nr of batches used, only needed for STROOPWAFEL (AIS) (default = -1, not needed)


	// Individual system variables
    primaryMass                                                     = 20.0;                                                                             // Initial primary mass in solar masses
    secondaryMass                                                   = 10.0;                                                                             // Initial secondary mass in solar masses

    initialPrimaryMetallicity                                       = 0.02;                                                                             // Initial metallicity of the primary
    initialSecondaryMetallicity                                     = 0.02;                                                                             // Initial metallicity of the secondary

    binarySeparation                                                = -1.0;                                                                             // Initial separation in AU
    binaryOrbitalPeriod                                             = -1.0;                                                                             // Initial orbital period in day
    binaryEccentricity                                              = 0.0;                                                                              // Initial eccentricity


    // Variables required to restart a binary/star part-way through
//    primaryStellarType                                              = -1;                                                                               // Initial primary stellar type (not yet implemented)
//    secondaryStellarType                                            = -1;                                                                               // Initial secondary stellar type (not yet implemented)

//    primaryEffectiveInitialMass                                     = primaryMass;                                                                      // Effective initial mass for the primary in solar masses (not yet implemented)
//    secondaryEffectiveInitialMass                                   = secondaryMass;                                                                    // Effective initial mass for the secondary in solar masses (not yet implemented)

//    primaryCoreMass                                                 = 0.0;                                                                              // Initial primary core mass in solar masses (not yet implemented)
//    secondaryCoreMass                                               = 0.0;                                                                              // Initial secondary core mass in solar masses (not yet implemented)

//    primaryAge                                                      = 0.0;                                                                              // Effective age for the primary star in Myrs (not yet implemented)
//    secondaryAge                                                    = 0.0;                                                                              // Effective age for the secondary star in Myrs (not yet implemented)

//    primaryRotationalVelocity                                       = 0.0;                                                                              // Initial rotational velocity of the primary (not yet implemented)
//    secondaryRotationalVelocity                                     = 0.0;                                                                              // Initial rotational velocity of the secondary (not yet implemented)


    // Public population synthesis variables
    nBinaries                                                       = 10;

    fixedRandomSeed                                                 = false;                                                                            // Whether to use a fixed random seed given by options.randomSeed (set to true if --random-seed is passed on command line) (default = false)
    randomSeed                                                      = 0;                                                                                // Random seed to use (default = 0)


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
    kickVelocityDistributionSigmaCCSN_NS                            = 250;                                                                              // Kick velocity sigma in km s^-1 for neutron stars (default = "250" )
    kickVelocityDistributionSigmaCCSN_BH                            = 250;                                                                              // Kick velocity sigma in km s^-1 for black holes (default = "250" )
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
    pairInstabilityUpperLimit                                       = 135.0;                                                                            // Maximum core mass leading to PISN (default = 135, Value in Belczynski+ 2016 is 135 Msol)
    pairInstabilityLowerLimit                                       = 60.0;                                                                             // Minimum core mass leading to PISN (default = 65,  Value in Belczynski+ 2016 is 65 Msol)

    usePulsationalPairInstability                                   = true;                                                                             // Whether to use pulsational pair instability (PPI)
    pulsationalPairInstabilityLowerLimit                            = 35.0;                                                                             // Minimum core mass leading to PPI, default = 40, Value in Belczynski+ 2016 is 45 Msol
    pulsationalPairInstabilityUpperLimit                            = 60.0;                                                                             // Maximum core mass leading to PPI, default = 65, Value in Belczynski+ 2016 is 65 Msol

    pulsationalPairInstabilityPrescription                          = PPI_PRESCRIPTION::COMPAS;                                                         // Prescription for PPI to use
    pulsationalPairInstabilityPrescriptionString                    = PPI_PRESCRIPTION_LABEL.at(pulsationalPairInstabilityPrescription);                // String for which PPI prescription to use

	maximumNeutronStarMass                                          = 3.0;									                                            // Maximum mass of a neutron star allowed, set to default in StarTrack


    // Kick direction option
    kickDirectionDistribution                                       = KICK_DIRECTION_DISTRIBUTION::ISOTROPIC;                                           // Which assumption for SN kicks: Possibilities: isotropic, inplane, perpendicular, powerlaw, wedge, poles, default = isotropic
    kickDirectionDistributionString                                 = KICK_DIRECTION_DISTRIBUTION_LABEL.at(kickDirectionDistribution);		            // Which assumption for SN kicks: Possibilities: isotropic, inplane, perpendicular, powerlaw, wedge, poles, default = isotropic
    kickDirectionPower                                              = 0.0;                                                                              // Power law power for the "power" SN kick direction choice

    // Get default output path
    outputPathString                                                = "";                                                                               // String to hold the output directory
    defaultOutputPath                                               = boost::filesystem::current_path();                                                // Default output location
    outputPath                                                      = defaultOutputPath;                                                                // Desired output location (default = CWD)

    // Spin options
    spinDistribution                                                = SPIN_DISTRIBUTION::FIXED;
    spinDistributionString                                          = SPIN_DISTRIBUTION_LABEL.at(spinDistribution);
    spinDistributionMin                                             = 0.60;
    spinDistributionMax                                             = 0.98;

    spinAssumption                                                  = SPIN_ASSUMPTION::ALIGNED;
    spinAssumptionString                                            = SPIN_ASSUMPTION_LABEL.at(spinAssumption);


    // Tides options
    tidesPrescription                                               = TIDES_PRESCRIPTION::NONE;                                                         // Tides prescription that will be used by the code
    tidesPrescriptionString                                         = TIDES_PRESCRIPTION_LABEL.at(tidesPrescription);                                   // String containing which tides prescription to use (default = "None")


    // Mass loss options
    useMassLoss                                                     = true;                                                                             // Whether to use mass loss

    massLossPrescription                                            = MASS_LOSS_PRESCRIPTION::VINK;
    massLossPrescriptionString                                      = MASS_LOSS_PRESCRIPTION_LABEL.at(massLossPrescription);


    // Wind mass loss multiplicitive constants
    luminousBlueVariableFactor                                      = 1.5;                                                                              // Luminous blue variable mass loss enhancement factor
    wolfRayetFactor                                                 = 1.0;                                                                              // WR winds factor

    // Mass transfer options
    useMassTransfer                                                 = true;                                                                             // Whether to use mass transfer (default = true)
	circulariseBinaryDuringMassTransfer         	                = false;						                                                    // Whether to circularise binary when it starts (default = false)
	forceCaseBBBCStabilityFlag                                      = true;									                                            // Whether if all case BB/BC systems are forced to be stable or unstable
	alwaysStableCaseBBBCFlag                                        = true;									                                            // Whether if case BB/BC is always stable

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
    massTransferCriticalMassRatioMSLowMassNonDegenerateAccretor     = 1.44;                                                                             // Critical mass ratio for MT from a MS low mass star (default = 1.44, Claeys+ 2014)
    massTransferCriticalMassRatioMSLowMassDegenerateAccretor        = 1.0;                                                                              // Critical mass ratio for MT from a MS low mass star on to a degenerate accretor (default = 1.0, Claeys+ 2014)

    massTransferCriticalMassRatioMSHighMass                         = false;							                                                // Whether to use critical mass ratios
    massTransferCriticalMassRatioMSHighMassNonDegenerateAccretor    = 0.625;                                                                            // Critical mass ratio for MT from a MS high mass star (default = 0.625, Claeys+ 2014)
    massTransferCriticalMassRatioMSHighMassDegenerateAccretor       = 0.0;                                                                              // Critical mass ratio for MT from a MS high mass star on to a degenerate accretor (default = 0)

    massTransferCriticalMassRatioHG                                 = false;                                                                            // Whether to use critical mass ratios
    massTransferCriticalMassRatioHGNonDegenerateAccretor            = 0.40;                                                                             // Critical mass ratio for MT from a HG star (default = 0.25, Claeys+ 2014)
    massTransferCriticalMassRatioHGDegenerateAccretor               = 0.21;                                                                             // Critical mass ratio for MT from a HG star on to a degenerate accretor (default = 0.21, Claeys+ 2014)

    massTransferCriticalMassRatioGiant                              = false;                                                                            // Whether to use critical mass ratios
    massTransferCriticalMassRatioGiantNonDegenerateAccretor         = 0.0;                                                                              // Critical mass ratio for MT from a giant (default = 0.0)
    massTransferCriticalMassRatioGiantDegenerateAccretor            = 0.87;                                                                             // Critical mass ratio for MT from a giant on to a degenerate accretor (default = 0.81, Claeys+ 2014)

    massTransferCriticalMassRatioHeliumMS                           = false;	                                                                        // Whether to use critical mass ratios
    massTransferCriticalMassRatioHeliumMSNonDegenerateAccretor      = 0.625;                                                                            // Critical mass ratio for MT from a Helium MS star (default = 0.625)
    massTransferCriticalMassRatioHeliumMSDegenerateAccretor         = 0.0;                                                                              // Critical mass ratio for MT from a Helium MS star on to a degenerate accretor (default = 0)

    massTransferCriticalMassRatioHeliumHG                           = false;                                                                            // Whether to use critical mass ratios
    massTransferCriticalMassRatioHeliumHGNonDegenerateAccretor      = 0.25;		                                                                       	// Critical mass ratio for MT from a Helium HG star (default = 0.25, de Claeys+ 2014)
    massTransferCriticalMassRatioHeliumHGDegenerateAccretor         = 0.21;		                                                                        // Critical mass ratio for MT from a Helium HG star on to a degenerate accretor (default = 0.21, Claeys+ 2014)

    massTransferCriticalMassRatioHeliumGiant                        = false;                                                                            // Whether to use critical mass ratios
    massTransferCriticalMassRatioHeliumGiantNonDegenerateAccretor   = 1.28;                                                                             // Critical mass ratio for MT from a Helium giant (default = 0.25, Claeys+ 2014)
    massTransferCriticalMassRatioHeliumGiantDegenerateAccretor      = 0.87;                                                                             // Critical mass ratio for MT from a Helium giant on to a degenerate accretor

    massTransferCriticalMassRatioWhiteDwarf                         = false;                                                                            // Whether to use critical mass ratios
	massTransferCriticalMassRatioWhiteDwarfNonDegenerateAccretor    = 0.0;                                                                              // Critical mass ratio for MT from a White Dwarf (default = 0.0)
    massTransferCriticalMassRatioWhiteDwarfDegenerateAccretor       = 1.6;                                                                              // Critical mass ratio for MT from a White Dwarf on to a degenerate accretor (default = 1.6, Claeys+ 2014)


    // Common Envelope parameters
    commonEnvelopePrescriptionFlag                                  = COMMON_ENVELOPE_PRESCRIPTION::WEBBINK;                                            // Which common envelope prescription to use
    commonEnvelopeAlpha                                             = 1.0;                                                                              // Common envelope efficiency alpha parameter (default = 1.0)
    commonEnvelopeLambda                                            = 0.1;                                                                              // Common envelope Lambda parameter (default = 0.1)
    commonEnvelopeHertzsprungGapDonor                               = COMMON_ENVELOPE_PRESCRIPTION::OPTIMISTIC_HG;                                      // Which prescription to use for Hertzsprung gap donors in a CE (default = OPTIMISTIC_HG_CE)
    commonEnvelopeHertzsprungGapDonorString                         = COMMON_ENVELOPE_PRESCRIPTION_LABEL.at(commonEnvelopeHertzsprungGapDonor);         // String containing which prescription to use for Hertzsprung gap donors in a CE (default = "OPTIMISTIC_HG_CE")
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
	commonEnvelopeLambdaPrescription                                = CE_LAMBDA_PRESCRIPTION::NANJING;                                                  // Which prescription to use for CE lambda (default = LAMBDA_NANJING)
	commonEnvelopeLambdaPrescriptionString                          = CE_LAMBDA_PRESCRIPTION_LABEL.at(commonEnvelopeLambdaPrescription);                // String containing which prescription to use for CE lambda (default = "LAMBDA_NANJING")

	// Common envelope Nandez and Ivanova energy formalism
	revisedEnergyFormalismNandezIvanova	                            = false;						                                                    // Use the revised energy formalism from Nandez & Ivanova 2016 (default = false)
	maximumMassDonorNandezIvanova                                   = 2.0;								                                                // Maximum mass allowed to use the revised energy formalism in Msol (default = 2.0)
	commonEnvelopeRecombinationEnergyDensity                        = 1.5E13;					                                                        // Factor using to calculate the binding energy depending on the mass of the envelope. (default = 1.5x10^13 ergs/g)

	// Common envelope power factor for Kruckow fit
	commonEnvelopeSlopeKruckow                                      = -4.0/5.0;								                                            // Common envelope power factor for Kruckow fit normalized according to Kruckow+2016, Fig. 1

	// Which prescription to use for calculating zetas
	commonEnvelopeZetaPrescription                                  = CE_ZETA_PRESCRIPTION::STARTRACK;					                                // Which prescription to use for calculating CE zetas (default = ZETA_ADIABATIC)
	commonEnvelopeZetaPrescriptionString                            = CE_ZETA_PRESCRIPTION_LABEL.at(commonEnvelopeZetaPrescription);				    // String containing which prescription to use for calculating CE zetas (default = STARTRACK)

	zetaAdiabaticArbitrary                                          = 0.0;
	zetaThermalArbitrary                                            = 0.0;
    zetaMainSequence 	                                            = 2.0;
	zetaHertzsprungGap	                                            = 6.5;

	// Which prescription to use for calculating zetas
	commonEnvelopeZetaPrescription                                  = CE_ZETA_PRESCRIPTION::STARTRACK;					                                // Which prescription to use for calculating CE zetas (default = ZETA_ADIABATIC)
	commonEnvelopeZetaPrescriptionString                            = CE_ZETA_PRESCRIPTION_LABEL.at(commonEnvelopeZetaPrescription);				    // String containing which prescription to use for calculating CE zetas (default = STARTRACK)


    // Afaptive Importance Sampling options
    AISexploratoryPhase                                             = false;                                                                            // Flag for whether to run the AIS exploratory phase
    AISDCOtype                                                      = AIS_DCO::ALL;                                                                     // Which prescription to use for DCO type (default = ALL)
    AISDCOtypeString                                                = AIS_DCO_LABEL.at(AIS_DCO::ALL);                                                   // String containing which type of DCOs to focus on (default = "ALL")
    AIShubble                                                       = false;                                                                            // Flag for excluding DCOs that do not merge in Hubble
    AISrlof                                                         = false;                                                                            // Flag for excluding DCOs that RLOFSecondaryZAMS
    AISpessimistic                                                  = false;                                                                            // Flag for excluding DCOs that are Optmistic
    AISrefinementPhase                                              = false;                                                                            // Flag for whether to run the AIS refinement phase (step 2)
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
    logfileBSERLOFParameters                                        = get<0>(LOGFILE_DESCRIPTOR.at(LOGFILE::BSE_RLOF_PARAMETERS));                      // get default filename from constants.h
    logfileBSEBeBinaries                                            = get<0>(LOGFILE_DESCRIPTOR.at(LOGFILE::BSE_BE_BINARIES));                          // get default filename from constants.h
    logfileBSEPulsarEvolution                                       = get<0>(LOGFILE_DESCRIPTOR.at(LOGFILE::BSE_PULSAR_EVOLUTION));                     // get default filename from constants.h
}


void Options::SetToFiducialValues(void) {

    // flags

    debugToFile                                                     = false;                                                                            // default is do not log debug statements to a log file
    errorsToFile                                                    = false;                                                                            // default is do not log error messages to a log file

    individualSystem                                                = false;                                                                            // Flag to evolve a specific individual system which you can specify initial parameters of
    singleStar                                                      = false;                                                                            // Flag to evolve a single star

	lambdaCalculationEveryTimeStep                                  = false;
	zetaCalculationEveryTimeStep                                    = false;

	beBinaries                                                      = false;
    evolvePulsars                                                   = false;                                                                            // Whether to evolve pulsars
	evolveUnboundSystems                                            = false;                                                                            // Allow unbound syetms to evolve
//    onlyDoubleCompactObjects                                        = false;                                                                            // Flag to turn on some shortcuts to only evolve systems which may form double compact objects

    detailedOutput                                                  = false;                                                                            // Detailed output
    populationDataPrinting                                          = false;                                                                            // Print certain data for small populations, but not for larger one
    printBoolAsString                                               = false;                                                                            // default is do not print bool as string
    quiet                                                           = false;                                                                            // Suppress some of the printing
    rlofPrinting                                                    = false;

    useImportanceSampling                                           = false;
//    useMCMC                                                         = false;

    nBatchesUsed                                                    = -1;                                                                               // nr of batches used, only needed for STROOPWAFEL (AIS) (default = -1, not needed)


    // Individual system variables
    primaryMass                                                     = 96.2;                                                                             // Initial primary mass in solar masses
    secondaryMass                                                   = 60.2;                                                                             // Initial secondary mass in solar masses

    initialPrimaryMetallicity                                       = 0.02;                                                                             // Initial metallicity of the primary
    initialSecondaryMetallicity                                     = 0.02;                                                                             // Initial metallicity of the secondary

    binarySeparation                                                = 11.5;                                                                             // Initial separation in AU
    binaryOrbitalPeriod                                             = -1.0;                                                                             // Initial orbital period in day
    binaryEccentricity                                              = 0.0;                                                                              // Initial eccentricity


    // Variables required to restart a binary/star halfway through
//    primaryStellarType                                              = 1;                                                                                // Initial primary stellar type (not yet implemented)
//    secondaryStellarType                                            = 1;                                                                                // Initial secondary stellar type (not yet implemented)

//    primaryEffectiveInitialMass                                     = 0.;                                                                               // Effective initial mass for the primary in solar masses (not yet implemented)
//    secondaryEffectiveInitialMass                                   = 0.;                                                                               // Effective initial mass for the secondary in solar masses (not yet implemented)

//    primaryCoreMass                                                 = 0.0;                                                                              // Initial primary core mass in solar masses (not yet implemented)
//    secondaryCoreMass                                               = 0.0;                                                                              // Initial secondary core mass in solar masses (not yet implemented)

//    primaryAge                                                      = 0.0;                                                                              // Effective age for the primary star in Myrs (not yet implemented)
//    secondaryAge                                                    = 0.0;                                                                              // Effective age for the secondary star in Myrs (not yet implemented)

//    primaryRotationalVelocity                                       = 0.0;                                                                              // Initial rotational velocity of the primary (not yet implemented)
//    secondaryRotationalVelocity                                     = 0.0;                                                                              // Initial rotational velocity of the secondary (not yet implemented)


    // Public population synthesis variables
    nBinaries                                                       = 10;

    fixedRandomSeed                                                 = true;                                                                             // Whether to use a fixed random seed given by options.randomSeed (set to true if --random-seed is passed on command line) (default = false)
    randomSeed                                                      = 0;                                                                                // Random seed to use (default = 0)


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

    minimumMassSecondary                                            = 0.1;                                                                              // Minimum mass of secondary to draw (in Msol) (default = 0.1, brown dwarf limit)


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
    eccentricityDistribution                                        = ECCENTRICITY_DISTRIBUTION::ZERO;                                                  // Which eccentricity distribution to use (default = "Zero")
    eccentricityDistributionString                                  = ECCENTRICITY_DISTRIBUTION_LABEL.at(eccentricityDistribution);                     // Which eccentricity distribution to use (default = "Zero")
    eccentricityDistributionMin                                     = 0.0;                                                                              // Minimum initial eccentricity to sample
    eccentricityDistributionMax                                     = 1.0;                                                                              // Maximum initial eccentricity to sample


    // Kick options
    kickVelocityDistribution                                        = KICK_VELOCITY_DISTRIBUTION::MAXWELLIAN;		                                    // Which kick velocity distribution to use
    kickVelocityDistributionString                                  = KICK_VELOCITY_DISTRIBUTION_LABEL.at(kickVelocityDistribution);		            // Which kick velocity distribution to use


    // Kick velocity options
    kickVelocityDistributionSigmaCCSN_NS                            = 250;                                                                              // Kick velocity sigma in km s^-1 for neutron stars (default = "250" )
    kickVelocityDistributionSigmaCCSN_BH                            = 250;                                                                              // Kick velocity sigma in km s^-1 for black holes (default = "250" )
    kickVelocityDistributionMaximum                                 = -1.0;                                                                             // Maximum kick velocity to draw in km s^-1. Ignored if < 0
    kickVelocityDistributionSigmaForECSN   	                        = 30.0;                                                                             // Characteristic kick velocity for an ECSN in km s^-1 (default = "30")
    kickVelocityDistributionSigmaForUSSN   	                        = 30.0;                                                                             // Characteristic kick velocity for an USSN in km s^-1 (default = "30")
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
    pairInstabilityUpperLimit                                       = 135.0;                                                                            // Maximum core mass leading to PISN (default = 135, Value in Belczynski+ 2016 is 135 Msol)
    pairInstabilityLowerLimit                                       = 60.0;                                                                             // Minimum core mass leading to PISN (default = 60, Value in Belczynski+ 2016 is 65 Msol)

    usePulsationalPairInstability                                   = true;                                                                             // Whether to use pulsational pair instability (PPI)
    pulsationalPairInstabilityLowerLimit                            = 35.0;                                                                             // Minimum core mass leading to PPI (default = 40, Value in Belczynski+ 2016 is 45 Msol)
    pulsationalPairInstabilityUpperLimit                            = 60.0;                                                                             // Maximum core mass leading to PPI (default = 60, Value in Belczynski+ 2016 is 65 Msol)

    pulsationalPairInstabilityPrescription                          = PPI_PRESCRIPTION::COMPAS;                                                         // Prescription for PPI to use
    pulsationalPairInstabilityPrescriptionString                    = PPI_PRESCRIPTION_LABEL.at(pulsationalPairInstabilityPrescription);                // String for which PPI prescription to use

	maximumNeutronStarMass                                          = 3.0;								                                                // Maximum mass of a neutron star allowed, set to default in StarTrack


    // Kick direction option
    kickDirectionDistribution                                       = KICK_DIRECTION_DISTRIBUTION::ISOTROPIC;                                           // Which assumption for SN kicks: Possibilities: isotropic, inplane, perpendicular, powerlaw, wedge, poles, default = isotropic
    kickDirectionDistributionString                                 = KICK_DIRECTION_DISTRIBUTION_LABEL.at(kickDirectionDistribution);		            // Which assumption for SN kicks: Possibilities: isotropic, inplane, perpendicular, powerlaw, wedge, poles, default = isotropic
    kickDirectionPower                                              = 0.0;                                                                              // Power law power for the "power" SN kick direction choice

    // Get default output path
    outputPathString                                                = ".";                                                                              // String to hold the output directory
    defaultOutputPath                                               = boost::filesystem::current_path();                                                // Default output location
    outputPath                                                      = defaultOutputPath;                                                                // Desired output location (default = CWD)

    // Spin options
    spinDistribution                                                = SPIN_DISTRIBUTION::ZERO;
    spinDistributionString                                          = SPIN_DISTRIBUTION_LABEL.at(spinDistribution);
    spinDistributionMin                                             = 0.0;
    spinDistributionMax                                             = 1.0;

    spinAssumption                                                  = SPIN_ASSUMPTION::ALIGNED;
    spinAssumptionString                                            = SPIN_ASSUMPTION_LABEL.at(spinAssumption);


    // Tides options
    tidesPrescription                                               = TIDES_PRESCRIPTION::NONE;                                                         // Tides prescription that will be used by the code
    tidesPrescriptionString                                         = TIDES_PRESCRIPTION_LABEL.at(tidesPrescription);                                   // String containing which tides prescription to use (default = "None")


    // Mass loss options
    useMassLoss                                                     = true;                                                                             // Whether to use mass loss

    massLossPrescription                                            = MASS_LOSS_PRESCRIPTION::VINK;
    massLossPrescriptionString                                      = MASS_LOSS_PRESCRIPTION_LABEL.at(massLossPrescription);


    // Wind mass loss multiplicitive constants
    luminousBlueVariableFactor                                      = 1.5;                                                                              // Luminous blue variable mass loss enhancement factor
    wolfRayetFactor                                                 = 1.0;                                                                              // WR winds factor


    // Mass transfer options
    useMassTransfer                                                 = true;											                                    // Whether to use mass transfer (default = true)
	circulariseBinaryDuringMassTransfer	                            = true;						                                                        // Whether to circularise binary when it starts (default = true)
	forceCaseBBBCStabilityFlag                                      = true;									                                            // Whether if all case BB/BC systems are forced to be stable or unstable
	alwaysStableCaseBBBCFlag                                        = true;									                                            // Whether if case BB/BC is always stable
	angularMomentumConservationDuringCircularisation                = true;		                                                                        // Whether to conserve angular momentum while circularising or circularise to periastron (default = true)


    // Mass transfer prescription options
    massTransferPrescription                                        = MT_PRESCRIPTION::DEMINK;
    massTransferPrescriptionString                                  = MT_PRESCRIPTION_LABEL.at(massTransferPrescription);


    // Options adaptive Roche Lobe Overflow prescription
    massTransferAdaptiveAlphaParameter                              = 0.5;
    maxPercentageAdaptiveMassTransfer                               = 0.01;


    // Options for mass transfer accretion efficiency
    massTransferFractionAccreted                                    = 1.0;
    massTransferCParameter                                          = 10.0;


    // Mass transfer accretion efficiency prescription options
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
    massTransferCriticalMassRatioMSLowMassNonDegenerateAccretor     = 1.44;                    			                                                // Critical mass ratio for MT from a MS low mass star (default = 1.44, Claeys+ 2014)
    massTransferCriticalMassRatioMSLowMassDegenerateAccretor        = 1.0;                    			                                                // Critical mass ratio for MT from a MS low mass star on to a degenerate accretor (default = 1.0, Claeys+ 2014)

    massTransferCriticalMassRatioMSHighMass                         = false;	           			                                                    // Whether to use critical mass ratios
    massTransferCriticalMassRatioMSHighMassNonDegenerateAccretor    = 0.625;                 			                                                // Critical mass ratio for MT from a MS high mass star (default = 0.625, Claeys+ 2014)
    massTransferCriticalMassRatioMSHighMassDegenerateAccretor       = 0.0;                      			                                            // Critical mass ratio for MT from a MS high mass star on to a degenerate accretor (default = 0)

    massTransferCriticalMassRatioHG                                 = false;                 			                                                // Whether to use critical mass ratios
    massTransferCriticalMassRatioHGNonDegenerateAccretor            = 0.40;                  			                                                // Critical mass ratio for MT from a HG star (default = 0.25, Claeys+ 2014)
    massTransferCriticalMassRatioHGDegenerateAccretor               = 0.21;                  			                                                // Critical mass ratio for MT from a HG star on to a degenerate accretor (default = 0.21, Claeys+ 2014)

    massTransferCriticalMassRatioGiant                              = false;                  			                                                // Whether to use critical mass ratios
    massTransferCriticalMassRatioGiantNonDegenerateAccretor         = 0.0;                   			                                                // Critical mass ratio for MT from a giant (default = 0.0)
    massTransferCriticalMassRatioGiantDegenerateAccretor            = 0.87;                  			                                                // Critical mass ratio for MT from a giant on to a degenerate accretor (default = 0.81, Claeys+ 2014)

    massTransferCriticalMassRatioHeliumMS                           = false;	           			                                                    // Whether to use critical mass ratios
    massTransferCriticalMassRatioHeliumMSNonDegenerateAccretor      = 0.625;                			                                                // Critical mass ratio for MT from a Helium MS star (default = 0.625)
    massTransferCriticalMassRatioHeliumMSDegenerateAccretor         = 0.0;                  			                                                // Critical mass ratio for MT from a Helium MS star on to a degenerate accretor (default = 0)

    massTransferCriticalMassRatioHeliumHG                           = false;                  			                                                // Whether to use critical mass ratios
    massTransferCriticalMassRatioHeliumHGNonDegenerateAccretor      = 0.25;		           			                                                    // Critical mass ratio for MT from a Helium HG star (default = 0.25, de Claeys+ 2014)
    massTransferCriticalMassRatioHeliumHGDegenerateAccretor         = 0.21;			           			                                                // Critical mass ratio for MT from a Helium HG star on to a degenerate accretor (default = 0.21, Claeys+ 2014)

    massTransferCriticalMassRatioHeliumGiant                        = false;                 			                                                // Whether to use critical mass ratios
    massTransferCriticalMassRatioHeliumGiantNonDegenerateAccretor   = 1.28;                 			                                                // Critical mass ratio for MT from a Helium giant (default = 0.25, Claeys+ 2014)
    massTransferCriticalMassRatioHeliumGiantDegenerateAccretor      = 0.87;                 			                                                // Critical mass ratio for MT from a Helium giant on to a degenerate accretor

    massTransferCriticalMassRatioWhiteDwarf                         = false;                			                                                // Whether to use critical mass ratios
	massTransferCriticalMassRatioWhiteDwarfNonDegenerateAccretor    = 0.0;                  			                                                // Critical mass ratio for MT from a White Dwarf (default = 0.0)
    massTransferCriticalMassRatioWhiteDwarfDegenerateAccretor       = 1.6;                     			                                                // Critical mass ratio for MT from a White Dwarf on to a degenerate accretor (default = 1.6, Claeys+ 2014)


    // Common Envelope parameters
    commonEnvelopePrescriptionFlag                                  = COMMON_ENVELOPE_PRESCRIPTION::WEBBINK;                                            // Which common envelope prescription to use
    commonEnvelopeAlpha                                             = 1.0;                                                                              // Common envelope efficiency alpha parameter (default = 1.0)
    commonEnvelopeLambda                                            = 0.1;                                                                              // Common envelope Lambda parameter (default = 0.1)
    commonEnvelopeHertzsprungGapDonor                               = COMMON_ENVELOPE_PRESCRIPTION::PESSIMISTIC_HG;                                     // Which prescription to use for Hertzsprung gap donors in a CE (default = OPTIMISTIC_HG_CE)
    commonEnvelopeHertzsprungGapDonorString                         = COMMON_ENVELOPE_PRESCRIPTION_LABEL.at(commonEnvelopeHertzsprungGapDonor);         // String containing which prescription to use for Hertzsprung gap donors in a CE (default = "OPTIMISTIC_HG_CE")
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
	commonEnvelopeLambdaPrescription                                = CE_LAMBDA_PRESCRIPTION::NANJING;                                                  // Which prescription to use for CE lambda (default = LAMBDA_NANJING)
	commonEnvelopeLambdaPrescriptionString                          = CE_LAMBDA_PRESCRIPTION_LABEL.at(commonEnvelopeLambdaPrescription);                // String containing which prescription to use for CE lambda (default = "LAMBDA_NANJING")


	// Common envelope Nandez and Ivanova energy formalism
	revisedEnergyFormalismNandezIvanova	                            = false;						                                                    // Use the revised energy formalism from Nandez & Ivanova 2016 (default = false)
	maximumMassDonorNandezIvanova                                   = 2.0;								                                                // Maximum mass allowed to use the revised energy formalism in Msol (default = 2.0)
	commonEnvelopeRecombinationEnergyDensity                        = 1.5E13;					                                                        // Factor using to calculate the bin


	// Common envelope power factor for Kruckow fit
	commonEnvelopeSlopeKruckow                                      = -2.0/3.0;								                                            // Common envelope power factor for Kruckow fit normalized according to Kruckow+2016, Fig. 1


	// Which prescription to use for calculating zetas
	commonEnvelopeZetaPrescription                                  = CE_ZETA_PRESCRIPTION::STARTRACK;					                                // Which prescription to use for calculating CE zetas (default = ZETA_ADIABATIC)
	commonEnvelopeZetaPrescriptionString                            = CE_ZETA_PRESCRIPTION_LABEL.at(commonEnvelopeZetaPrescription);					// String containing which prescription to use for calculating CE zetas (default = STARTRACK)


	zetaAdiabaticArbitrary                                          = 0.0;
	zetaThermalArbitrary                                            = 0.0;
    zetaMainSequence 	                                            = 6.5;
	zetaHertzsprungGap	                                            = 2.0;


    // Adaptive Importance Sampling Exploratory phase
    AISexploratoryPhase                                             = false;  // Floor
    AISDCOtype                                                      = AIS_DCO::ALL;                                                                     // Which prescription to use for DCO type (default = ALL)
    AISDCOtypeString                                                = AIS_DCO_LABEL.at(AIS_DCO::ALL);                                                   // String containing which type of DCOs to focus on (default = "ALL")
    AIShubble                                                       = false;                                                                            // Flag for excluding DCOs that do not merge in Hubble
    AISrlof                                                         = false;                                                                            // Flag for excluding DCOs that RLOFSecondaryZAMS
    AISpessimistic                                                  = false;                                                                            // Flag for excluding DCOs that are Optmistic
    AISrefinementPhase                                              = false;                                                                            // Flag for whether to run the AIS refinement phase (step 2)
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
    logfileBSERLOFParameters                                        = get<0>(LOGFILE_DESCRIPTOR.at(LOGFILE::BSE_RLOF_PARAMETERS));                      // get default filename from constants.h
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

    namespace po = boost::program_options;
    namespace fs = boost::filesystem;

    COMMANDLINE_STATUS programStatus = COMMANDLINE_STATUS::CONTINUE;

    try {

        // Create program options object
        po::options_description desc("Options");
        desc.add_options()

		    // JR: todo: should make these names consistent ( case, hyphenated, camelCase... )

		    // JR: todo: some of the strings below declare the default value for the option - and some of them are wrong
		    // (probably have become wrong over time).  I think we should either not show the default value, or if we do
		    // then construct the string with the default value.  The second option is more work...


		    // boolean options (flags) - alphabetically

            ("AIS-exploratory-phase",                               "Run exploratory phase of STROOPWAFEL") // Floor
		    ("AIS-Hubble",                                          "Excluding not in Hubble time mergers selection in exploratory phase of STROOPWAFEL")
		    ("AIS-Pessimistic",                                     "Optimistic or Pessimistic selection in exploratory phase of STROOPWAFEL")
		    ("AIS-refinement-phase",                                "If true: run main sampling phase (step2) of STROOPWAFEL")
		    ("AIS-RLOF",                                            "RLOFSecondaryZAMS selection in exploratory phase of STROOPWAFEL")
			("alwaysStableCaseBBBCFlag",                            "Choose case BB/BC mass transfer to be always stable (default = True)")
			("angularMomentumConservationDuringCircularisation",    "Conserve angular momentum when binary is circularised when entering a Mass Transfer episode (default = False)")
			("BeBinaries",                                          "Enable Be Binaries study")
			("circulariseBinaryDuringMassTransfer",                 "Circularise binary when it enters a Mass Transfer episode (default = False)")
		    ("common-envelope-allow-main-sequence-survive",         "Allow main sequence stars to survive common envelope evolution")
			("debug-to-file",                                       "Write debug statements to file")
		    ("detailedOutput",                                      "Print detailed output to file")
			("errors-to-file",                                      "Write error messages to file")
		    ("evolve-pulsars",                                      "Whether to evolve pulsars")
			("evolve-unbound-systems",                              "Keeps evolving stars even if the binary is disrupted")
            ("forceCaseBBBCStabilityFlag",                          "Force case BB/BC mass transfer to be only stable or unstable (default = True)")
		    ("help,h",                                              "Print this help message")
		    ("individual-system",                                   "Run an individual system")
			("lambda-calculation-every-timeStep",                   "Calculate all values of lambda at each timestep")
   		   	("massTransfer",                                        "Enable mass transfer")
//			("mcmc",                                                "Use MCMC sampling (Not yet implemented. default = false)")
//		    ("only-double-compact-objects",                         "Only evolve binaries which may form double compact objects")
		    ("pair-instability-supernovae",                         "Enable pair instability supernovae (PISN)")
            ("populationDataPrinting",                              "Print details of population")
		    ("print-bool-as-string",                                "Print boolean properties as 'TRUE' or 'FALSE'")
		    ("pulsational-pair-instability",                        "Enable mass loss due to pulsational-pair-instability (PPI)")
		    ("quiet",                                               "Suppress printing")
			("revised-energy-formalism-Nandez-Ivanova",             "Enable revised energy formalism")
            ("RLOFPrinting",                                        "Enable output parameters before/after RLOF ")
			("sample-common-envelope-alpha",                        "Sample over common envelope alpha")
			("sample-kick-direction-power",                         "Sample over kick direction powerlaw exponent")
			("sample-kick-velocity-sigma",                          "Sample over Kick Velocity Sigma")
			("sample-luminous-blue-variable-multiplier",            "Sample over multiplicative constant from LBV mass loss")
			("sample-wolf-rayet-multiplier",                        "Sample over WR winds multiplicative constant")
            ("single-star",                                         "Evolve single star(s)")
		    ("use-mass-loss",                                       "Enable mass loss")
			("zeta-calculation-every-time-Step",                    "Calculate all values of MT zetas at each timestep")


			// numerical options - alphabetically by groups

			// unsigned long
		    ("random-seed",                                                 po::value<unsigned long>(&randomSeed),                                              "Random seed to use (default = 0)")

		    // int
			("debug-level",                                                 po::value<int>(&debugLevel),                                                        "Determines which print statements are displayed for debugging")
//		    ("individual-initial-primary-type",                             po::value<int>(&primaryStellarType),                                                "Initial stellar type for primary (not yet implemented) (default ZAMS from mass)")
//		    ("individual-initial-secondary-type",                           po::value<int>(&secondaryStellarType),                                              "Initial stellar type for secondary (not yet implemented) (default ZAMS from mass)")
		    ("log-level",                                                   po::value<int>(&logLevel),                                                          "Determines which print statements are included in the logfile")
		    ("maximum-number-iterations",                                   po::value<int>(&maxNumberOfTimestepIterations),                                     "Maximum number of timesteps to evolve binary before giving up (default = 99999)")
			("number-of-binaries,n",                                        po::value<int>(&nBinaries),                                                         "Specify the number of binaries to simulate (default = 1000000)")
			("single-star-mass-steps",                                      po::value<int>(&singleStarMassSteps),                                               "Specify the number of mass steps for single star evolution (default = 100)")

		    // double
		    ("common-envelope-alpha",                                       po::value<double>(&commonEnvelopeAlpha),                                            "Common Envelope efficiency alpha (default = 1.0)")
		    ("common-envelope-alpha-thermal",                               po::value<double>(&commonEnvelopeAlphaThermal),                                     "Defined such that lambda = alpha_th * lambda_b + (1.0 - alpha_th) * lambda_g (default = 1.0)")
		    ("common-envelope-lambda",                                      po::value<double>(&commonEnvelopeLambda),                                           "Common Envelope lambda (default = 0.1)")
		    ("common-envelope-lambda-multiplier",                           po::value<double>(&commonEnvelopeLambdaMultiplier),                                 "Multiply lambda by some constant (default = 1.0)")
	    	("common-envelope-mass-accretion-constant",                     po::value<double>(&commonEnvelopeMassAccretionConstant),                            "Value of mass accreted by NS/BH during common envelope evolution if assuming all NS/BH accrete same amount of mass (common-envelope-mass-accretion-prescription CONSTANT). Ignored otherwise (default = 0)")
		    ("common-envelope-mass-accretion-max",                          po::value<double>(&commonEnvelopeMassAccretionMax),                                 "Maximum amount of mass accreted by NS/BHs during common envelope evolution in solar masses (default = 0.1)")
		    ("common-envelope-mass-accretion-min",                          po::value<double>(&commonEnvelopeMassAccretionMin),                                 "Minimum amount of mass accreted by NS/BHs during common envelope evolution in solar masses (default = 0.04)")
		    ("common-envelope-recombination-energy-density",                po::value<double>(&commonEnvelopeRecombinationEnergyDensity),                       "Recombination energy density in ergs/g (default = 1.5x10^13)")
		    ("common-envelope-slope-Kruckow",                               po::value<double>(&commonEnvelopeSlopeKruckow),                                     "Common Envelope slope for Kruckow lambda (default = -4/5)")
            ("critical-mass-ratio-giant-degenerate-accretor",               po::value<double>(&massTransferCriticalMassRatioGiantDegenerateAccretor),           "Critical mass ratio for MT from a giant star (default = 0.87 from Claeys+ 2014) Specify both giant flags to use. 0 is always stable, <0 is disabled")
            ("critical-mass-ratio-giant-non-degenerate-accretor",           po::value<double>(&massTransferCriticalMassRatioGiantNonDegenerateAccretor),        "Critical mass ratio for MT from a giant star (default = not implemented from Claeys+ 2014) Specify both giant flags to use. 0 is always stable, <0 is disabled")
            ("critical-mass-ratio-helium-giant-degenerate-accretor",        po::value<double>(&massTransferCriticalMassRatioHeliumGiantDegenerateAccretor),     "Critical mass ratio for MT from a helium giant star (default = 0.87 from Claeys+ 2014) Specify both helium giant flags to use. 0 is always stable, <0 is disabled")
            ("critical-mass-ratio-helium-giant-non-degenerate-accretor",    po::value<double>(&massTransferCriticalMassRatioHeliumGiantNonDegenerateAccretor),  "Critical mass ratio for MT from a helium giant star (default = 1.28 from Claeys+ 2014) Specify both helium giant flags to use. 0 is always stable, <0 is disabled")
            ("critical-mass-ratio-helium-HG-degenerate-accretor",           po::value<double>(&massTransferCriticalMassRatioHeliumHGDegenerateAccretor),        "Critical mass ratio for MT from a helium HG star (default = 0.21 from Claeys+ 2014) Specify both helium HG flags to use. 0 is always stable, <0 is disabled")
            ("critical-mass-ratio-helium-HG-non-degenerate-accretor",       po::value<double>(&massTransferCriticalMassRatioHeliumHGNonDegenerateAccretor),     "Critical mass ratio for MT from a helium HG star (default = 0.25 from Claeys+ 2014) Specify both helium HG flags to use. 0 is always stable, <0 is disabled")
            ("critical-mass-ratio-helium-MS-degenerate-accretor",           po::value<double>(&massTransferCriticalMassRatioHeliumMSDegenerateAccretor),        "Critical mass ratio for MT from a helium MS star (default = 0.0) Specify both helium MS flags to use. 0 is always stable, <0 is disabled")
            ("critical-mass-ratio-helium-MS-non-degenerate-accretor",       po::value<double>(&massTransferCriticalMassRatioHeliumMSNonDegenerateAccretor),     "Critical mass ratio for MT from a helium MS star (default = 0.625) Specify both helium MS flags to use. 0 is always stable, <0 is disabled")
            ("critical-mass-ratio-HG-degenerate-accretor",                  po::value<double>(&massTransferCriticalMassRatioHGDegenerateAccretor),              "Critical mass ratio for MT from a HG star (default = 0.21 from Claeys+ 2014) Specify both HG flags to use. 0 is always stable, <0 is disabled")
            ("critical-mass-ratio-HG-non-degenerate-accretor",              po::value<double>(&massTransferCriticalMassRatioHGNonDegenerateAccretor),           "Critical mass ratio for MT from a HG star (default = 0.40 from de Mink+ 2013) Specify both HG flags to use. 0 is always stable, <0 is disabled")
            ("critical-mass-ratio-MS-high-mass-degenerate-accretor",        po::value<double>(&massTransferCriticalMassRatioMSHighMassDegenerateAccretor),      "Critical mass ratio for MT from a MS star to a degenerate accretor (default = 0.0 from Claeys+ 2014) Specify both MS high mass flags to use. 0 is always stable, <0 is disabled")
            ("critical-mass-ratio-MS-high-mass-non-degenerate-accretor",    po::value<double>(&massTransferCriticalMassRatioMSHighMassNonDegenerateAccretor),   "Critical mass ratio for MT from a MS star (default = 0.625, Claeys+ 2014). Specify both MS high mass flags to use. 0 is always stable, <0 is disabled")
            ("critical-mass-ratio-MS-low-mass-degenerate-accretor",         po::value<double>(&massTransferCriticalMassRatioMSLowMassDegenerateAccretor),       "Critical mass ratio for MT from a MS star to a degenerate accretor (default = 1.0 from Claeys+ 2014) Specify both MS low mass flags to use. 0 is always stable, <0 is disabled")
            ("critical-mass-ratio-MS-low-mass-non-degenerate-accretor",     po::value<double>(&massTransferCriticalMassRatioMSLowMassNonDegenerateAccretor),    "Critical mass ratio for MT from a MS star (default = 1.44, Claeys+ 2014). Specify both MS low mass flags to use. 0 is always stable, <0 is disabled")
            ("critical-mass-ratio-white-dwarf-degenerate-accretor",         po::value<double>(&massTransferCriticalMassRatioWhiteDwarfDegenerateAccretor),      "Critical mass ratio for MT from a white dwarf (default = 1.6 from Claeys+ 2014) Specify both white dwarf flags to use. 0 is always stable, <0 is disabled")
            ("critical-mass-ratio-white-dwarf-non-degenerate-accretor",     po::value<double>(&massTransferCriticalMassRatioWhiteDwarfNonDegenerateAccretor),   "Critical mass ratio for MT from a white dwarf (default = 0.0) Specify both white dwarf flags to use. 0 is always stable, <0 is disabled")

		    ("eccentricity-max",                                            po::value<double>(&eccentricityDistributionMax),                                    "Maximum eccentricity to generate (default = 1.0)")
		    ("eccentricity-min",                                            po::value<double>(&eccentricityDistributionMin),                                    "Minimum eccentricity to generate (default = 0.0)")
			("eddington-accretion-factor",                                  po::value<double>(&eddingtonAccretionFactor),                                       "Multiplication factor for eddington accretion for NS & BH, i.e. >1 is super-eddington and 0. is no accretion")

   		    ("fix-dimensionless-kick-velocity",                             po::value<double>(&fixedUK),                                                        "Fix dimensionless kick velocity uk to this value (default = -1, -ve values false, +ve values true)")

//		    ("individual-effective-initial-primary-mass",                   po::value<double>(&primaryEffectiveInitialMass),                                    "Effective initial mass for primary in Msol (default = Mass)")
//		    ("individual-effective-initial-secondary-mass",                 po::value<double>(&secondaryEffectiveInitialMass),                                  "Effective initial mass for secondary in Msol (default = Mass)")
		    ("individual-initial-orbital-eccentricity",                     po::value<double>(&binaryEccentricity),                                             "Initial orbital eccentricity (default = 0.0)")
		    ("individual-initial-orbital-period",                           po::value<double>(&binaryOrbitalPeriod),                                            "Initial orbital period in days (default from masses and separation)")
		    ("individual-initial-orbital-separation",                       po::value<double>(&binarySeparation),                                               "Initial orbital separation in AU (default = 1.0)")
//		    ("individual-initial-primary-age",                              po::value<double>(&primaryAge),                                                     "Initial age for primary in Myrs (default = 0.0)")
//		    ("individual-initial-primary-core-mass",                        po::value<double>(&primaryCoreMass),                                                "Initial core mass for primary in Msol (default = 0.0)")
		    ("individual-initial-primary-mass",                             po::value<double>(&primaryMass),                                                    "Initial mass for primary in Msol (default = 20.0)")
		    ("individual-initial-primary-metallicity",                      po::value<double>(&initialPrimaryMetallicity),                                      "Initial metallicity for primary (default = 0.02)")
//		    ("individual-initial-primary-rotational-velocity",              po::value<double>(&primaryRotationalVelocity),                                      "Initial rotational velocity for primary (not yet implemented) (default = 0.0)")
//		    ("individual-initial-secondary-age",                            po::value<double>(&secondaryAge),                                                   "Initial age for secondary in Myrs (default = 0.0)")
//		    ("individual-initial-secondary-core-mass",                      po::value<double>(&secondaryCoreMass),                                              "Initial core mass for secondary in Msol (default = 0.0)")
		    ("individual-initial-secondary-mass",                           po::value<double>(&secondaryMass),                                                  "Initial mass for secondary in Msol (default = 10.0)")
		    ("individual-initial-secondary-metallicity",                    po::value<double>(&initialSecondaryMetallicity),                                    "Initial metallicity for secondary (default = 0.02)")
//		    ("individual-initial-secondary-rotational-velocity",            po::value<double>(&secondaryRotationalVelocity),                                    "Initial rotational velocity for secondary (not yet implemented) (default = 0.0)")
		    ("initial-mass-max",                                            po::value<double>(&initialMassFunctionMax),                                         "Maximum mass (in Msol) to generate using given IMF (default = 100)")
		    ("initial-mass-min",                                            po::value<double>(&initialMassFunctionMin),                                         "Minimum mass (in Msol) to generate using given IMF (default = 8)")
		    ("initial-mass-power",                                          po::value<double>(&initialMassFunctionPower),                                       "Single power law power to generate primary mass using given IMF")

		    ("kappa-gaussians",                                             po::value<double>(&kappaGaussians),                                                 "Scaling factor for the width of the Gaussian distributions in STROOPWAFEL main sampling phase" )
		    ("kick-direction-power",                                        po::value<double>(&kickDirectionPower),                                             "Power for power law kick direction distribution (default = 0.0 = isotropic, +ve = polar, -ve = in plane)")
			("kick-scaling-factor",                                         po::value<double>(&kickScalingFactor),                                              "Arbitrary factor used to scale kicks (default = 1.0 )")
		    ("kick-velocity-max",                                           po::value<double>(&kickVelocityDistributionMaximum),                                "Maximum drawn kick velocity in km s^-1. Ignored if < 0. Must be > 0 if using kick-velocity-distribution=FLAT")
		    ("kick-velocity-sigma-CCSN-BH",                                 po::value<double>(&kickVelocityDistributionSigmaCCSN_BH),                           "Sigma for chosen kick velocity distribution for black holes (default = 250 km s^-1 )")
		    ("kick-velocity-sigma-CCSN-NS",                                 po::value<double>(&kickVelocityDistributionSigmaCCSN_NS),                           "Sigma for chosen kick velocity distribution for neutron stars (default = 250 km s^-1 )")
			("kick-velocity-sigma-ECSN",                                    po::value<double>(&kickVelocityDistributionSigmaForECSN),                           "Sigma for chosen kick velocity distribution for ECSN (default = 0 km s^-1 )")
			("kick-velocity-sigma-USSN",                                    po::value<double>(&kickVelocityDistributionSigmaForUSSN),                           "Sigma for chosen kick velocity distribution for USSN (default = 20 km s^-1 )")

		    ("luminous-blue-variable-multiplier",                           po::value<double>(&luminousBlueVariableFactor),                                     "Multiplicitive constant for LBV mass loss (default = 1.5, use 10 for Mennekens & Vanbeveren 2014)")

		    ("mass-ratio-max",                                              po::value<double>(&massRatioDistributionMax),                                       "Maximum mass ratio m2/m1 to generate (default = 1.0)")
		    ("mass-ratio-min",                                              po::value<double>(&massRatioDistributionMin),                                       "Minimum mass ratio m2/m1 to generate (default = 0.0)")
		    ("mass-transfer-fa",                                            po::value<double>(&massTransferFractionAccreted),                                   "Mass Transfer fraction accreted (default = 1.0, fully conservative)")
		    ("mass-transfer-jloss",                                         po::value<double>(&massTransferJloss),                                              "Specific angular momentum with which the non-accreted system leaves the system (default = 1.0)")
			("mass-transfer-thermal-limit-C",                               po::value<double>(&massTransferCParameter),                                         "Mass Transfer Thermal rate factor fo the accretor (default = 10.0, Hurley+2002)")
		    ("maximum-evolution-time",                                      po::value<double>(&maxEvolutionTime),                                               "Maximum time to evolve binaries in Myrs (default = 5.0)")
		    ("maximum-mass-donor-Nandez-Ivanova",                           po::value<double>(&maximumMassDonorNandezIvanova),                                  "Maximum donor mass allowed for the revised common envelope formalism in Msol (default = 2.0)")
			("maximum-neutron-star-mass",                                   po::value<double>(&maximumNeutronStarMass),                                         "Maximum mass of a neutron star (default = 3.0, as in StarTrack)")
            ("metallicity,z",                                               po::value<double>(&metallicity),                                                    "Metallicity to use (default 0.02 is Zsol)")
		    ("minimum-secondary-mass",                                      po::value<double>(&minimumMassSecondary),                                           "Minimum mass of secondary to generate in Msol (default = 0.0)")

			("nbatches-used",                                               po::value<int>(&nBatchesUsed),                                                      "nr of batches used, only needed for STROOPWAFEL (AIS) (default = -1, not needed) sets itself automatically in pythonSubmit")

		    ("orbital-period-max",                                          po::value<double>(&periodDistributionMax),                                          "Maximum period in days to generate (default = 1000)")
		   	("orbital-period-min",                                          po::value<double>(&periodDistributionMin),                                          "Minimum period in days to generate (default = 1.1)")

		    ("PISN-lower-limit",                                            po::value<double>(&pairInstabilityLowerLimit),                                      "Minimum core mass for PISN (default = 60.0)")
		    ("PISN-upper-limit",                                            po::value<double>(&pairInstabilityUpperLimit),                                      "Maximum core mass for PISN (default = 135.0)")
		    ("PPI-lower-limit",                                             po::value<double>(&pulsationalPairInstabilityLowerLimit),                           "Minimum core mass for PPI (default = 35.0)")
		    ("PPI-upper-limit",                                             po::value<double>(&pulsationalPairInstabilityUpperLimit),                           "Maximum core mass for PPI (default = 60.0)")
		    ("pulsar-birth-magnetic-field-distribution-max",                po::value<double>(&pulsarBirthMagneticFieldDistributionMax),                        "Maximum (log10) pulsar birth magnetic field (default = 13.0)")
		    ("pulsar-birth-magnetic-field-distribution-min",                po::value<double>(&pulsarBirthMagneticFieldDistributionMin),                        "Minimum (log10) pulsar birth magnetic field) (default = 11.0)")
		    ("pulsar-birth-spin-period-distribution-max",                   po::value<double>(&pulsarBirthSpinPeriodDistributionMax),                           "Maximum pulsar birth spin period in ms (default = 100.0)")
		    ("pulsar-birth-spin-period-distribution-min",                   po::value<double>(&pulsarBirthSpinPeriodDistributionMin),                           "Minimum pulsar birth spin period in ms (default = 0.0)")
		    ("pulsar-magnetic-field-decay-massscale",                       po::value<double>(&pulsarMagneticFieldDecayMassscale),                              "Mass scale on which magnetic field decays during accretion in solar masses (default = 0.025)")
		    ("pulsar-magnetic-field-decay-timescale",                       po::value<double>(&pulsarMagneticFieldDecayTimescale),                              "Timescale on which magnetic field decays in Myrs (default = 1000.0)")
		    ("pulsar-minimum-magnetic-field",                               po::value<double>(&pulsarLog10MinimumMagneticField),                                "log10 of the minimum pulsar magnetic field in Gauss (default = 8.0)")

			("sample-common-envelope-alpha-max",                            po::value<double>(&sampleCommonEnvelopeAlphaMax),                                   "Maximum for Uniform sampling over common envelope alpha")
			("sample-common-envelope-alpha-min",                            po::value<double>(&sampleCommonEnvelopeAlphaMin),                                   "Minimum for Uniform sampling over common envelope alpha")
			("sample-kick-direction-power-max",                             po::value<double>(&sampleKickDirectionPowerMax),                                    "Maximum for Uniform sampling over kick direction powerlaw exponent")
			("sample-kick-direction-power-min",                             po::value<double>(&sampleKickDirectionPowerMin),                                    "Minimum for Uniform sampling over kick direction powerlaw exponent")
			("sample-kick-velocity-sigma-max",                              po::value<double>(&sampleKickVelocitySigmaMax),                                     "Maximum for Uniform sampling over kick velocity sigma")
			("sample-kick-velocity-sigma-min",                              po::value<double>(&sampleKickVelocitySigmaMin),                                     "Minimum for Uniform sampling over kick velocity sigma")
			("sample-luminous-blue-variable-multiplier-max",                po::value<double>(&sampleLuminousBlueVariableMultiplierMax),                        "Maximum for Uniform sampling over multiplicative constant for LBV mass loss")
			("sample-luminous-blue-variable-multiplier-min",                po::value<double>(&sampleLuminousBlueVariableMultiplierMin),                        "Minimum for Uniform sampling over multiplicative constant for LBV mass loss")
			("sample-wolf-rayet-multiplier-max",                            po::value<double>(&sampleWolfRayetMultiplierMax),                                   "Maximum for Uniform sampling over multiplicative constant for WR winds")
			("sample-wolf-rayet-multiplier-min",                            po::value<double>(&sampleWolfRayetMultiplierMin),                                   "Minimum for Uniform sampling over multiplicative constant for WR winds")
		    ("semi-major-axis-max",                                         po::value<double>(&semiMajorAxisDistributionMax),                                   "Maximum semi major axis in AU to generate (default = 1000)")
		    ("semi-major-axis-min",                                         po::value<double>(&semiMajorAxisDistributionMin),                                   "Minimum semi major axis in AU to generate (default = 0.1)")
		    ("single-star-mass-max",                                        po::value<double>(&singleStarMassMax),                                              "Maximum mass (in Msol) for single star evolution (default = 100.0)")
            ("single-star-mass-min",                                        po::value<double>(&singleStarMassMin),                                              "Minimum mass (in Msol) for single star evolution (default = 5.0)")
		    ("spin-mag-max",                                                po::value<double>(&spinDistributionMax),                                            "Maximum magnitude of spin (default = )")
		    ("spin-mag-min",                                                po::value<double>(&spinDistributionMin),                                            "Minimum magnitude of spin (default = )")

		    ("wolf-rayet-multiplier",                                       po::value<double>(&wolfRayetFactor),                                                "Multiplicitive constant for WR winds (default = 1.0)")

		    ("zeta-adiabatic-arbitrary",                                    po::value<double>(&zetaAdiabaticArbitrary),                                         "Value of mass-radius exponent zeta adiabatic")
		    ("zeta-hertzsprung-gap",                                        po::value<double>(&zetaHertzsprungGap),                                             "Value of mass-radius exponent zeta on the hertzstrpung gap (default = 6.5)")
		    ("zeta-main-sequence",                                          po::value<double>(&zetaMainSequence),                                               "Value of mass-radius exponent zeta on the main sequence (default = 2.0)")


		    // string options - alphabetically
            ("AIS-DCOtype",                                                 po::value<string>(&AISDCOtypeString),                                               "DCO type selection in exploratory phase of STROOPWAFEL, select ALL, BBH, BNS or BHNS")

		  	("black-hole-kicks",                                            po::value<string>(&blackHoleKicksString),                                           "Black hole kicks relative to NS kicks (options: FULL, REDUCED, ZERO, FALLBACK. Default = FALLBACK)")

		  	("chemically-homogeneous-evolution",                            po::value<string>(&cheString),                                                      "Chemically Homogenesous Evolution (options: NONE, OPTIMISTIC, PESSIMISTIC. Default = NONE)")

			("common-envelope-hertzsprung-gap-assumption",                  po::value<string>(&commonEnvelopeHertzsprungGapDonorString),                        "Assumption to make about HG stars in CE (default = OPTIMISTIC_HG_CE)")
			("common-envelope-lambda-prescription",                         po::value<string>(&commonEnvelopeLambdaPrescriptionString),                         "Prescription for CE lambda (options: LAMBDA_FIXED, LAMBDA_LOVERIDGE, LAMBDA_NANJING, LAMBDA_KRUCKOW, LAMBDA_DEWI. Default = LAMBDA_FIXED)")
		    ("common-envelope-mass-accretion-prescription",                 po::value<string>(&commonEnvelopeMassAccretionPrescriptionString),                  "Assumption about whether NS/BHs can accrete mass during common envelope evolution (options: ZERO, CONSTANT, UNIFORM, MACLEOD+2014 . Default = ZERO)")
			("common-envelope-zeta-prescription",                           po::value<string>(&commonEnvelopeZetaPrescriptionString),                           "Prescription for CE zeta (default = STARTRACK)")

		    ("eccentricity-distribution,e",                                 po::value<string>(&eccentricityDistributionString),                                 "Initial eccentricity distribution, e (options: ZERO, FIXED, FLAT, THERMALISED, GELLER+2013. Default = ZERO)")

		    ("fryer-supernova-engine",                                      po::value<string>(&fryerSupernovaEngineString),                                     "If using Fryer et al 2012 fallback prescription, select between 'delayed' and 'rapid' engines (default = 'delayed')")

            ("grid",                                                        po::value<string>(&gridFilename)->implicit_value(""),                               "Grid filename - SSE or BSE")

		    ("initial-mass-function,i",                                     po::value<string>(&initialMassFunctionString),                                      "Specify initial mass function to use (options: SALPETER, POWERLAW, UNIFORM, KROUPA. default = KROUPA)")

		    ("kick-direction",                                              po::value<string>(&kickDirectionDistributionString),                                "Distribution for natal kick direction (options: ISOTROPIC, INPLANE, PERPENDICULAR, POWERLAW, WEDGE, POLES. Default = ISOTROPIC)")
		    ("kick-velocity-distribution",                                  po::value<string>(&kickVelocityDistributionString),                                 "Natal kick velocity distribution (options: ZERO, FLAT, MAXWELLIAN, MUELLER2016, MUELLER2016MAXWELLIAN, BRAYELDRIDGE. Default = MAXWELLIAN)")

            ("logfile-BSE-be-binaries",                                     po::value<string>(&logfileBSEBeBinaries),                                           "Filename for BSE Be Binaries logfile")
            ("logfile-BSE-common-envelopes",                                po::value<string>(&logfileBSECommonEnvelopes),                                      "Filename for BSE Common Envelopes logfile")
            ("logfile-BSE-detailed-output",                                 po::value<string>(&logfileBSEDetailedOutput),                                       "Filename for BSE Detailed Output logfile")
            ("logfile-BSE-double-compact-objects",                          po::value<string>(&logfileBSEDoubleCompactObjects),                                 "Filename for BSE Double Compact Objects logfile")
            ("logfile-BSE-pulsar-evolution",                                po::value<string>(&logfileBSEPulsarEvolution),                                      "Filename for BSE Pulsar Evolution logfile")
            ("logfile-BSE-rlof-parameters",                                 po::value<string>(&logfileBSERLOFParameters),                                       "Filename for BSE RLOF Parameters logfile")
            ("logfile-BSE-supernovae",                                      po::value<string>(&logfileBSESupernovae),                                           "Filename for BSE Supernovae logfile")
            ("logfile-BSE-system-parameters",                               po::value<string>(&logfileBSESystemParameters),                                     "Filename for BSE System Parameters logfile")
            ("logfile-definitions",                                         po::value<string>(&logfileDefinitionsFilename)->implicit_value(""),                 "Filename for logfile record definitions")
            ("logfile-delimiter",                                           po::value<string>(&logfileDelimiterString),                                         "Field delimiter for logfile records.  Default = TAB")
            ("logfile-name-prefix",                                         po::value<string>(&logfileNamePrefix),                                              "Prefix for logfile names")
            ("logfile-SSE-parameters",                                      po::value<string>(&logfileSSEParameters),                                           "Filename for SSE Parameters logfile")

		    ("mass-loss-prescription",                                      po::value<string>(&massLossPrescriptionString),                                     "Mass loss prescription to use (options: NONE, HURLEY, VINK. Default = NONE)")
		    ("mass-ratio-distribution,q",                                   po::value<string>(&massRatioDistributionString),                                    "Initial mass ratio distribution for q=m2/m1 (options: FLAT, DuquennoyMayor1991, SANA2012. Default = FLAT)")
		    ("mass-transfer-accretion-efficiency-prescription",             po::value<string>(&massTransferAccretionEfficiencyPrescriptionString),              "Mass Transfer Accretion Efficiency prescription to use (options: THERMAL, FIXED, CENTRIFUGAL. Default = THERMAL)")
		    ("mass-transfer-angular-momentum-loss-prescription",            po::value<string>(&massTransferAngularMomentumLossPrescriptionString),              "Mass Transfer Angular Momentum Loss prescription to use (options: JEANS, ISOTROPIC, CIRCUMBINARY, ARBITRARY. Default = ISOTROPIC)")
		    ("mass-transfer-prescription",                                  po::value<string>(&massTransferPrescriptionString),                                 "Mass Transfer prescription to use (default = DEMINK)")
		    ("mass-transfer-rejuvenation-prescription",                     po::value<string>(&massTransferRejuvenationPrescriptionString),                     "Mass Transfer Rejuvenation prescription to use (options: NONE, STARTRACK. Default = NONE)")
			("mass-transfer-thermal-limit-accretor",                        po::value<string>(&massTransferThermallyLimitedVariationString),                    "Mass Transfer Thermal Accretion limit to use (default = CFACTOR)")

		    ("neutron-star-equation-of-state",                              po::value<string>(&neutronStarEquationOfStateString),                               "Specify which neutron star equation of state to use (options: SSE, ARP3 default = SSE)")

   		    ("outputPath,o",                                                po::value<string>(&outputPathString),                                               "Directory for output (default = CWD)")

		    ("pulsar-birth-magnetic-field-distribution",                    po::value<string>(&pulsarBirthMagneticFieldDistributionString),                     "Distribution of (log10 of) pulsar birth magnetic field in G (options: ZERO, FIXED, FLATINLOG, UNIFORM, LOGNORMAL. Default = ZERO)")
		    ("pulsar-birth-spin-period-distribution",                       po::value<string>(&pulsarBirthSpinPeriodDistributionString),                        "Distribution of pulsar birth spin period in ms (options: ZERO, FIXED, UNIFORM, NORMAL. Default = ZERO)")
		    ("pulsational-pair-instability-prescription",                   po::value<string>(&pulsationalPairInstabilityPrescriptionString),                   "Specify which prescription to use for pulsational pair instability (options: COMPAS, STARTRACK, MARCHANT default = COMPAS)")

		    ("remnant-mass-prescription",                                   po::value<string>(&remnantMassPrescriptionString),                                  "Choose remnant mass prescription (options: postitnote, hurley2000, belczynski2002, fryer2012, muller2016, muller2016Maxwellian. Default = fryer2012)")
		    ("rotational-velocity-distribution",                            po::value<string>(&rotationalVelocityDistributionString),                           "Initial rotational velocity distribution (options: ZERO, HURLEY, VLTFLAMES. Default = ZERO)")

		    ("semi-major-axis-distribution,a",                              po::value<string>(&semiMajorAxisDistributionString),                                "Initial semi-major axis distribution, a (options: FLATINLOG, CUSTOM, DuquennoyMayor1991, SANA2012. default = FLATINLOG)")
		    ("spin-assumption",                                             po::value<string>(&spinAssumptionString),                                           "Which assumption of misalignedments to use (default = bothAligned)")
		    ("spin-distribution",                                           po::value<string>(&spinDistributionString),                                         "Which distribution of spins to use (default = 0)")

		    ("tides-prescription",                                          po::value<string>(&tidesPrescriptionString),                                        "Tides prescription to use (options: default = None)")


            // vector (list) options - alphabetically
            ("debug-classes",                                               po::value<vector<string>>(&debugClasses)->multitoken(),                             "Debug classes enabled")
            ("log-classes",                                                 po::value<vector<string>>(&logClasses)->multitoken(),                               "Logging classes enabled")
		;


        po::variables_map vm;   // Variables map

        try {

            po::store(po::parse_command_line(argc, argv, desc), vm);

            // --help option
            if (vm.count("help")) {                                                                                                     // user requested help
                utils::SplashScreen();                                                                                                  // yes - show splash screen
                ANNOUNCE(desc);                                                                                                         // and help
                programStatus = COMMANDLINE_STATUS::SUCCESS;                                                                            // ok
            }

            po::notify(vm);                                                                                                             // invoke notify to assign user-input values to variables.  Throws an error, so do after help just in case there are any problems.


            // boolean options (flags) - alphabetically (where possible - dependencies)

            AISexploratoryPhase                             = vm.count("AIS-exploratory-phase");                                        // exploratory phase of Adaptive Importance Sampling - Floor 24-04-2018.  Do not retain previous (default) value
            AIShubble                                       = vm.count("AIS-Hubble");                                                   // excluding binaries that merge outside Hubble time (exploratory phase AIS)?  Do not retain previous (default) value
            AISpessimistic                                  = vm.count("AIS-Pessimistic");                                              // excluding binaries that are Optimistic (exploratory phase AIS)?  Do not retain previous (default) value
            AISrefinementPhase                              = vm.count("AIS-refinement-phase");
            AISrlof                                         = vm.count("AIS-RLOF");                                                     // excluding binaries that RLOFSecondaryZAMS (exploratory phase AIS)?  Do not retain previous (default) value

            allowMainSequenceStarToSurviveCommonEnvelope    = vm.count("common-envelope-allow-main-sequence-survive");                  // allow MS stars to survive CE event?  Do not retain previous (default) value

			beBinaries                                      = vm.count("BeBinaries");                                                   // enable Be Binaries?  Do not retain previous (default) value

            debugToFile                                     = vm.count("debug-to-file") ? true : debugToFile;                           // write debug records to file?  Retain previous (default) value if not specified

            detailedOutput                                  = vm.count("detailedOutput");                                               // detailed output of each simulated system?  Do not retain previous (default) value

            errorsToFile                                    = vm.count("errors-to-file") ? true : errorsToFile;                         // write error messages to file?  Retain previous (default) value

            evolvePulsars                                   = vm.count("evolve-pulsars");                                               // evolve pulsars?  Do not retain previous (default) value

            evolveUnboundSystems                            = vm.count("evolve-unbound-systems");                                       // evolve unbound systems?  Do not retain previous (default) value

            fixedMetallicity                                = vm.count("metallicity") ? true : fixedMetallicity;                        // user-specified a metallicity value?  Retain previous (default) value

            fixedRandomSeed                                 = fixedRandomSeed || vm.count("random-seed");                               // user-specified random seed?  Do not retain previous (default) value

            individualSystem                                = vm.count("individual-system");                                            // run an individual system or do population synthesis?  Do not retain previous (default) value

			lambdaCalculationEveryTimeStep                  = vm.count("lambda-calculation-every-timeStep");                            // calculate lambdas at every timestep?  Do not retain previous (default) value if not specified

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

//            onlyDoubleCompactObjects                        = vm.count("only-double-compact-objects");                                  // only evolve DCOs?  Do not retain previous (default) value

            populationDataPrinting                          = vm.count("populationDataPrinting");                                       // print certain values while running a population?  Do not retain previous (default) valued

            printBoolAsString                               = vm.count("print-bool-as-string") ? true : printBoolAsString;              // print boolean values as "TRUE"/"FALSE"?  Retain previous (default) value

            quiet                                           = vm.count("quiet");                                                        // verbose or quiet mode?  Do not retain previous (default) value

            revisedEnergyFormalismNandezIvanova             = vm.count("revised-energy-formalism-Nandez-Ivanova");                      // Do not retain previous (default) value

            rlofPrinting                                    = vm.count("RLOFPrinting");                                                 // print Roche Lobe overflow details?  Do not retain previous (default) value

			sampleCommonEnvelopeAlpha                       = vm.count("sample-common-envelope-alpha");                                 // sample CE alpha?  Do not retain previous (default) value
			sampleKickDirectionPower                        = vm.count("sample-kick-direction-power");                                  // sample kick direction power?  Do not retain previous (default) value
			sampleKickVelocitySigma                         = vm.count("sample-kick-velocity-sigma");                                   // sample over some hyperparameters? - JIM BARRETT -- 06/07/2016.  Do not retain previous (default) value
			sampleLuminousBlueVariableMultiplier            = vm.count("sample-luminous-blue-variable-multiplier");                     // sample LBV multiplier?  Do not retain previous (default) value
			sampleWolfRayetMultiplier                       = vm.count("sample-wolf-rayet-multiplier");                                 // sample wolf-rayet multiplier?  Do not retain previous (default) valued

            singleStar                                      = vm.count("single-star");                                                  // evolve a single star?  Do not retain previous (default) valued

            useFixedUK                                      = vm.count("fix-dimensionless-kick-velocity") && (utils::Compare(fixedUK, 0.0) >= 0);   // fix the dimensionless kick velocity?  Do not retain previous (default) value

			useImportanceSampling                           = vm.count("importance-sampling");                                          // Do not retain previous (default) value

            useMassTransfer                                 = vm.count("massTransfer") ? true : useMassTransfer;                        // use mass transfer?  Retain previous (default) value

            useMassLoss                                     = vm.count("use-mass-loss");                                                // use mass loss?  Do not retain previous (default) value

//			useMCMC                                         = vm.count("mcmc");			                                                // Simon Stevenson - 15/03/2018 - adding MCMC functionality.  Do not retain previous (default) value

            usePairInstabilitySupernovae                    = vm.count("pair-instability-supernovae");                                  // pair instability supernovae?  Do not retain previous (default) value

            usePulsationalPairInstability                   = vm.count("pulsational-pair-instability");                                 // pulsational pair instability supernovae?  Do not retain previous (default) value

			zetaCalculationEveryTimeStep                    = vm.count("zeta-Calculation-Every-Time-Step");                             // calculate zetas at every timestep?  Do not retain previous (default) value


			if (useMassTransfer) {                                                                                                      // the following depend on the value of useMassTransfer
                if (vm.count("circulariseBinaryDuringMassTransfer")) {
					angularMomentumConservationDuringCircularisation = vm.count("angularMomentumConservationDuringCircularisation");    // Do not retain previous (default) value
					circulariseBinaryDuringMassTransfer              = true;
                }
				else{
					circulariseBinaryDuringMassTransfer              = false;
				}

                if (vm.count("forceCaseBBBCStabilityFlag")) {
	                alwaysStableCaseBBBCFlag   = vm.count("alwaysStableCaseBBBCFlag");                                                  // Do not retain previous (default) value
					forceCaseBBBCStabilityFlag = true;
                }
				else{
					forceCaseBBBCStabilityFlag = false;
				}
            }


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

            if (individualSystem && !gridFilename.empty()) individualSystem = false;                                                    // ignore individual-system if have a grid filename

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

        default:                                                                                                        // unknown property
            ok    = false;                                                                                              // that's not ok...
            value = "UNKNOWN";                                                                                          // default value
            std::cerr << ERR_MSG(ERROR::UNKNOWN_PROGRAM_OPTION) << std::endl;                                           // show warning (don't have logging or errors here...)
    }

    return std::make_tuple(ok, value);
}


