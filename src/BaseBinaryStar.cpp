#include "BaseBinaryStar.h"
#include <fenv.h>

#include "vector3d.h"

// gsl includes
#include <gsl/gsl_poly.h>


/* Constructor
 *
 * Parameter p_Seed is the seed for the random number generator - see main.cpp for an
 * explanation of how p_Seed is derived.
 * 
 * Parameter p_Id is the id of the binary - effectively an index - which is added as
 * a suffix to the filenames of any detailed output files created.
 */


// binary is generated according to distributions specified in program options
BaseBinaryStar::BaseBinaryStar(const unsigned long int p_Seed, const long int p_Id) {

    ERROR error = ERROR::NONE;

    SetInitialValues(p_Seed, p_Id);                                                                                                     // start construction of the binary
                        
    // generate initial properties of binary
    // check that the constituent stars are not touching
    // also check m2 > m2min

    bool done                            = false;
    bool merger                          = false;
    bool rlof                            = false;
    bool secondarySmallerThanMinimumMass = false;

    // determine if any if the initial conditions are sampled
    // we consider eccentricity distribution = ECCENTRICITY_DISTRIBUTION::ZERO to be not sampled!
    // we consider metallicity distribution = METALLICITY_DISTRIBUTION::ZSOLAR to be not sampled!
    bool sampled = !OPTIONS->OptionSpecified("initial-mass-1")  ||
                   !OPTIONS->OptionSpecified("initial-mass-2")  ||
                  (!OPTIONS->OptionSpecified("metallicity")     && OPTIONS->MetallicityDistribution() != METALLICITY_DISTRIBUTION::ZSOLAR) ||
                  (!OPTIONS->OptionSpecified("semi-major-axis") && !OPTIONS->OptionSpecified("orbital-period"))                            ||
                  (!OPTIONS->OptionSpecified("eccentricity")    && OPTIONS->EccentricityDistribution() != ECCENTRICITY_DISTRIBUTION::ZERO);


    // Single stars are provided with a kick structure that specifies the values of the random
    // number to be used to generate to kick magnitude, and the actual kick magnitude specified
    // by the user via program option --kick-magnitude 
    //
    // See typedefs.h for the kick structure.
    //
    // We can't just pick up the values of the options inside Basestar.cpp because the constituents
    // of binaries get different values, so use different options. The Basestar.cpp code doesn't 
    // know if the star is a single star (SSE) or a constituent of a binary (BSE) - it only knows 
    // that it is a star - so we have to setup the kick structures here for each constituent star.

    KickParameters kickParameters1;
    kickParameters1.magnitudeRandomSpecified = OPTIONS->OptionSpecified("kick-magnitude-random-1");
    kickParameters1.magnitudeRandom          = OPTIONS->KickMagnitudeRandom1();
    kickParameters1.magnitudeSpecified       = OPTIONS->OptionSpecified("kick-magnitude-1");
    kickParameters1.magnitude                = OPTIONS->KickMagnitude1();
    kickParameters1.phiSpecified             = OPTIONS->OptionSpecified("kick-phi-1");
    kickParameters1.phi                      = OPTIONS->SN_Phi1();
    kickParameters1.thetaSpecified           = OPTIONS->OptionSpecified("kick-theta-1");
    kickParameters1.theta                    = OPTIONS->SN_Theta1();
    kickParameters1.meanAnomalySpecified     = OPTIONS->OptionSpecified("kick-mean-anomaly-1");
    kickParameters1.meanAnomaly              = OPTIONS->SN_MeanAnomaly1();

    KickParameters kickParameters2;
    kickParameters2.magnitudeRandomSpecified = OPTIONS->OptionSpecified("kick-magnitude-random-2");
    kickParameters2.magnitudeRandom          = OPTIONS->KickMagnitudeRandom2();
    kickParameters2.magnitudeSpecified       = OPTIONS->OptionSpecified("kick-magnitude-2");
    kickParameters2.magnitude                = OPTIONS->KickMagnitude2();
    kickParameters2.phiSpecified             = OPTIONS->OptionSpecified("kick-phi-2");
    kickParameters2.phi                      = OPTIONS->SN_Phi2();
    kickParameters2.thetaSpecified           = OPTIONS->OptionSpecified("kick-theta-2");
    kickParameters2.theta                    = OPTIONS->SN_Theta2();
    kickParameters2.meanAnomalySpecified     = OPTIONS->OptionSpecified("kick-mean-anomaly-2");
    kickParameters2.meanAnomaly              = OPTIONS->SN_MeanAnomaly2();

    // loop here to find initial conditions that suit our needs
    // if the user supplied all initial conditions, no loop
    // loop for a maximum of MAX_BSE_INITIAL_CONDITIONS_ITERATIONS - it hasn't (that I
    // know of) been a problem in the past, but we should have a guard on the loop so
    // that we don't loop forever - probably more important now that the user can specify
    // initial conditions (so might leave the insufficient space for the  (say) one to be
    // sampled...)

    int tries = 0;
    do {

        double mass1 = OPTIONS->OptionSpecified("initial-mass-1")                                                                       // user specified primary mass?
                        ? OPTIONS->InitialMass1()                                                                                       // yes, use it
                        : utils::SampleInitialMass(OPTIONS->InitialMassFunction(),                                                      // no - sample it 
                                                   OPTIONS->InitialMassFunctionMax(), 
                                                   OPTIONS->InitialMassFunctionMin(), 
                                                   OPTIONS->InitialMassFunctionPower());

        double mass2 = 0.0;                      
        if (OPTIONS->OptionSpecified("initial-mass-2")) {                                                                               // user specified secondary mass?
            mass2 = OPTIONS->InitialMass2();                                                                                            // yes, use it
        }
        else {                                                                                                                          // no - sample it
            // first, determine mass ratio q    
            double q = OPTIONS->OptionSpecified("mass-ratio")                                                                           // user specified mass ratio?
                        ? OPTIONS->MassRatio()                                                                                          // yes, use it
                        : utils::SampleMassRatio(OPTIONS->MassRatioDistribution(),                                                      // no - sample it
                                                 OPTIONS->MassRatioDistributionMax(), 
                                                 OPTIONS->MassRatioDistributionMin());

            mass2 = mass1 * q;                                                                                                          // calculate mass2 using mass ratio                                                                     
        }

        double metallicity = OPTIONS->OptionSpecified("metallicity")                                                                    // user specified metallicity?
                                ? OPTIONS->Metallicity()                                                                                // yes, use it
                                : utils::SampleMetallicity(OPTIONS->MetallicityDistribution(),                                          // no, sample it
                                                           OPTIONS->MetallicityDistributionMax(), 
                                                           OPTIONS->MetallicityDistributionMin());

        if (OPTIONS->OptionSpecified("semi-major-axis")) {                                                                              // user specified semi-major axis?
            m_SemiMajorAxis = OPTIONS->SemiMajorAxis();                                                                                 // yes, use it
        }
        else {                                                                                                                          // no, semi-major axis not specified
            if (OPTIONS->OptionSpecified("orbital-period")) {                                                                           // user specified orbital period?
                m_SemiMajorAxis = utils::ConvertPeriodInDaysToSemiMajorAxisInAU(mass1, mass2, OPTIONS->OrbitalPeriod());                // yes - calculate semi-major axis from period
            }
            else {                                                                                                                      // no
                if (OPTIONS->OptionSpecified("semi-major-axis-distribution") ||                                                         // user specified semi-major axis distribution, or
                   !OPTIONS->OptionSpecified("orbital-period-distribution" )) {                                                         // user did not specify oprbital period distribution
                    ERROR error;
                    std::tie(error, m_SemiMajorAxis) = utils::SampleSemiMajorAxis(OPTIONS->SemiMajorAxisDistribution(),                 // yes, sample from semi-major axis distribution (might be default), assumes Opik's law (-1.0 exponent)
                                                                                  OPTIONS->SemiMajorAxisDistributionMax(), 
                                                                                  OPTIONS->SemiMajorAxisDistributionMin(),
                                                                                  OPIKS_LAW_SEMIMAJOR_AXIS_DISTRIBUTION_POWER,
                                                                                  OPTIONS->OrbitalPeriodDistributionMax(), 
                                                                                  OPTIONS->OrbitalPeriodDistributionMin(), 
                                                                                  mass1, 
                                                                                  mass2);
                    THROW_ERROR_IF(error == ERROR::UNKNOWN_SEMI_MAJOR_AXIS_DISTRIBUTION, error, "Sampling semi-major axis");            // throw error if necessary
                    SHOW_WARN_IF(error == ERROR::NO_CONVERGENCE, error, "Sampling semi-major axis");                                    // show warning if necessary
                }
                else {                                                                                                                  // no - sample from orbital period distribution
                    double orbitalPeriod = utils::SampleOrbitalPeriod(OPTIONS->OrbitalPeriodDistribution(),                              
                                                                      OPTIONS->OrbitalPeriodDistributionMax(), 
                                                                      OPTIONS->OrbitalPeriodDistributionMin());

                    m_SemiMajorAxis = utils::ConvertPeriodInDaysToSemiMajorAxisInAU(mass1, mass2, orbitalPeriod);                       // calculate semi-major axis from period
                }
            }
        }

        m_Eccentricity = OPTIONS->OptionSpecified("eccentricity")                                                                       // user specified eccentricity?
                            ? OPTIONS->Eccentricity()                                                                                   // yes, use it
                            : utils::SampleEccentricity(OPTIONS->EccentricityDistribution(),                                            // no, sample it
                                                        OPTIONS->EccentricityDistributionMax(), 
                                                        OPTIONS->EccentricityDistributionMin());

        // binary star contains two instances of star to hold masses, radii and luminosities.
        // star 1 initially more massive
        m_Star1 = OPTIONS->OptionSpecified("rotational-frequency-1")                                                                    // user specified primary rotational frequency?
                    ? new BinaryConstituentStar(m_RandomSeed, mass1, metallicity, kickParameters1, OPTIONS->RotationalFrequency1() * SECONDS_IN_YEAR) // yes - use it (convert from Hz to cycles per year - see BaseStar::CalculateZAMSAngularFrequency())
                    : new BinaryConstituentStar(m_RandomSeed, mass1, metallicity, kickParameters1);                                     // no - let it be calculated

        m_Star2 = OPTIONS->OptionSpecified("rotational-frequency-2")                                                                    // user specified secondary rotational frequency?
                    ? new BinaryConstituentStar(m_RandomSeed, mass2, metallicity, kickParameters2, OPTIONS->RotationalFrequency2() * SECONDS_IN_YEAR) // yes - use it (convert from Hz to cycles per year - see BaseStar::CalculateZAMSAngularFrequency())
                    : new BinaryConstituentStar(m_RandomSeed, mass2, metallicity, kickParameters2);                                     // no - let it be calculated

        double starToRocheLobeRadiusRatio1 = (m_Star1->Radius() * RSOL_TO_AU) / (m_SemiMajorAxis * (1.0 - m_Eccentricity) * CalculateRocheLobeRadius_Static(mass1, mass2));
        double starToRocheLobeRadiusRatio2 = (m_Star2->Radius() * RSOL_TO_AU) / (m_SemiMajorAxis * (1.0 - m_Eccentricity) * CalculateRocheLobeRadius_Static(mass2, mass1));

        m_Flags.massesEquilibrated         = false;                                                                                     // default
        m_Flags.massesEquilibratedAtBirth  = false;                                                                                     // default

        rlof = utils::Compare(starToRocheLobeRadiusRatio1, 1.0) > 0 || utils::Compare(starToRocheLobeRadiusRatio2, 1.0) > 0;            // either star overflowing Roche Lobe?

        if (rlof && OPTIONS->AllowRLOFAtBirth()) {                                                                                      // over-contact binaries at birth allowed?    
            m_Flags.massesEquilibratedAtBirth = true;                                                                                   // record that we've equilbrated at birth

            mass1            = (mass1 + mass2) / 2.0;                                                                                   // equilibrate masses
            mass2            = mass1;                                                                                                   // ditto
            
            double M         = mass1 + mass2;
            double m1m2      = mass1 * mass2;
            m_SemiMajorAxis *= 16.0 * m1m2 * m1m2 / (M * M * M * M) * (1.0 - (m_Eccentricity * m_Eccentricity));                        // circularise; conserve angular momentum

            m_Eccentricity   = 0.0;                                                                                                     // now circular

            // create new stars with equal masses - all other ZAMS values recalculated
            delete m_Star1;
            m_Star1 = OPTIONS->OptionSpecified("rotational-frequency-1")                                                                // user specified primary rotational frequency?
                        ? new BinaryConstituentStar(m_RandomSeed, mass1, metallicity, kickParameters1, OPTIONS->RotationalFrequency1() * SECONDS_IN_YEAR) // yes - use it (convert from Hz to cycles per year - see BaseStar::CalculateZAMSAngularFrequency())
                        : new BinaryConstituentStar(m_RandomSeed, mass1, metallicity, kickParameters1);                                 // no - let it be calculated

            delete m_Star2;
            m_Star2 = OPTIONS->OptionSpecified("rotational-frequency-2")                                                                // user specified secondary rotational frequency?
                        ? new BinaryConstituentStar(m_RandomSeed, mass2, metallicity, kickParameters2, OPTIONS->RotationalFrequency2() * SECONDS_IN_YEAR) // yes - use it (convert from Hz to cycles per year - see BaseStar::CalculateZAMSAngularFrequency())
                        : new BinaryConstituentStar(m_RandomSeed, mass2, metallicity, kickParameters2);                                 // no - let it be calculated
        
            starToRocheLobeRadiusRatio1 = (m_Star1->Radius() * RSOL_TO_AU) / (m_SemiMajorAxis * CalculateRocheLobeRadius_Static(mass1, mass2)); //eccentricity already zero
            starToRocheLobeRadiusRatio2 = (m_Star2->Radius() * RSOL_TO_AU) / (m_SemiMajorAxis * CalculateRocheLobeRadius_Static(mass2, mass1));
        }

        m_Star1->SetCompanion(m_Star2);
        m_Star2->SetCompanion(m_Star1);

        merger                          = (m_SemiMajorAxis * AU_TO_RSOL) < (m_Star1->Radius() + m_Star2->Radius());
        secondarySmallerThanMinimumMass = utils::Compare(mass2, OPTIONS->MinimumMassSecondary()) < 0;

        // check whether our initial conditions are good
        // if they are - evolve the binary
        // if they are not ok:
        //    - if we sampled at least one of them, sample again
        //    - if all were user supplied, set error - Evolve() will show the error and return without evolving

        bool ok = !((!OPTIONS->AllowRLOFAtBirth() && rlof) || (!OPTIONS->AllowTouchingAtBirth() && merger) || secondarySmallerThanMinimumMass);

        done = ok;
        if (!sampled && !ok) {
            error = ERROR::INVALID_INITIAL_ATTRIBUTES;
            done = true;
        }

    } while (!done && ++tries < MAX_BSE_INITIAL_CONDITIONS_ITERATIONS);

    if (!done) error = ERROR::INVALID_INITIAL_ATTRIBUTES;                                                                               // too many iterations - bad initial conditions

    if (error != ERROR::NONE) THROW_ERROR(error);                                                                                       // throw error if necessary

    SetRemainingValues();                                                                                                               // complete the construction of the binary
}


/*
 * Initiate the construction of the binary - initial values
 *
 *
 * void SetInitialValues(const long int p_Id)
 *
 * @param   [IN]    p_Id                        Ordinal value of binary - see constructor notes above
 */
void BaseBinaryStar::SetInitialValues(const unsigned long int p_Seed, const long int p_Id) {

    m_Error             = ERROR::NONE;                                                                                                  // we can safely set this here

    m_ObjectId          = globalObjectId++;
    m_ObjectPersistence = OBJECT_PERSISTENCE::PERMANENT;
    m_RandomSeed        = p_Seed;
    m_Id                = p_Id;

    m_EvolutionStatus   = EVOLUTION_STATUS::CONTINUE;

    if (OPTIONS->PopulationDataPrinting()) {                                                                                            // user wants to see details of binary?
        SAY("Using supplied random seed " << m_RandomSeed << " for Binary Star id = " << m_ObjectId);                                   // yes - show them
    }
}


/*
 * Complete the construction of the binary - remaining values
 *
 *
 * void SetRemainingValues()
 */
void BaseBinaryStar::SetRemainingValues() {

    // Initialise other parameters
    m_SemiMajorAxisPrev           = m_SemiMajorAxis;
    m_EccentricityPrev            = m_Eccentricity;

    // initial binary parameters - kept constant as a record of the initial parameters of the binary
    m_SemiMajorAxisInitial        = m_SemiMajorAxis;
    m_EccentricityInitial         = m_Eccentricity;

    // initialise variables to hold parameters prior to supernova explosion
    m_SemiMajorAxisPreSN          = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_EccentricityPreSN           = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_OrbitalVelocityPreSN        = DEFAULT_INITIAL_DOUBLE_VALUE;

    // initialise variables to hold parameters at DCO formation
    m_SemiMajorAxisAtDCOFormation = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_EccentricityAtDCOFormation  = DEFAULT_INITIAL_DOUBLE_VALUE;

    double momentOfInertia1       = m_Star1->CalculateMomentOfInertiaAU();
    double momentOfInertia2       = m_Star2->CalculateMomentOfInertiaAU();

    m_TotalEnergy                 = CalculateTotalEnergy(m_SemiMajorAxis, m_Star1->Mass(), m_Star2->Mass(), m_Star1->Omega(), m_Star2->Omega(), momentOfInertia1, momentOfInertia2);

    m_TotalAngularMomentum        = CalculateAngularMomentum(m_SemiMajorAxis, m_Eccentricity, m_Star1->Mass(), m_Star2->Mass(), m_Star1->Omega(), m_Star2->Omega(), momentOfInertia1, momentOfInertia2);
    m_TotalAngularMomentumPrev    = m_TotalAngularMomentum;
    
    m_Omega                       = 0.0;

    if (OPTIONS->CHEMode() != CHE_MODE::NONE) {                                                                                                         // CHE enabled?

        // CHE enabled, update rotational frequency for constituent stars - assume tidally locked

        double omega = OrbitalAngularVelocity();                                                                                                        // orbital angular velocity

        m_Star1->SetOmega(omega);
        m_Star2->SetOmega(omega);

        // check for CHE
        //
        // because we've changed the rotational frequency of the constituent stars we
        // have to reset the stellar type - at this stage, based on their rotational
        // frequency at birth, they may have already been assigned one of MS_LTE_07,
        // MS_GT_07, or CHEMICALLY_HOMOGENEOUS
        //
        // here we need to change from MS_* -> CH, or from CH->MS* based on the
        // newly-assigned rotational frequencies

        // star 1
        if (utils::Compare(m_Star1->Omega(), m_Star1->OmegaCHE()) >= 0) {                                                                               // star 1 CH?
            if (m_Star1->StellarType() != STELLAR_TYPE::CHEMICALLY_HOMOGENEOUS) (void)m_Star1->SwitchTo(STELLAR_TYPE::CHEMICALLY_HOMOGENEOUS, true);    // yes, switch if not already Chemically Homogeneous
        }
        else if (m_Star1->MZAMS() <= 0.7) {                                                                                                             // no - MS - initial mass determines actual type  (don't use utils::Compare() here)
            if (m_Star1->StellarType() != STELLAR_TYPE::MS_LTE_07) (void)m_Star1->SwitchTo(STELLAR_TYPE::MS_LTE_07, true);                              // MS <= 0.7 Msol - switch if necessary
        }
        else {
            if (m_Star1->StellarType() != STELLAR_TYPE::MS_GT_07) (void)m_Star1->SwitchTo(STELLAR_TYPE::MS_GT_07, true);                                // MS > 0.7 Msol - switch if necessary
        }

        // star 2
        if (utils::Compare(m_Star2->Omega(), m_Star2->OmegaCHE()) >= 0) {                                                                               // star 2 CH?
            if (m_Star2->StellarType() != STELLAR_TYPE::CHEMICALLY_HOMOGENEOUS) (void)m_Star2->SwitchTo(STELLAR_TYPE::CHEMICALLY_HOMOGENEOUS, true);    // yes, switch if not already Chemically Homogeneous
        }
        else if (m_Star2->MZAMS() <= 0.7) {                                                                                                             // no - MS - initial mass determines actual type  (don't use utils::Compare() here)
            if (m_Star2->StellarType() != STELLAR_TYPE::MS_LTE_07) (void)m_Star2->SwitchTo(STELLAR_TYPE::MS_LTE_07, true);                              // MS <= 0.0 Msol - switch if necessary
        }
        else {
            if (m_Star2->StellarType() != STELLAR_TYPE::MS_GT_07) (void)m_Star2->SwitchTo(STELLAR_TYPE::MS_GT_07, true);                                // MS > 0.7 Msol - switch if necessary
        }

        // if both stars evolving as chemically homogeneous stars set m_Omega for binary
        if (HasTwoOf({STELLAR_TYPE::CHEMICALLY_HOMOGENEOUS})) m_Omega = omega;
    }

	double totalMass 					             = m_Star1->Mass() + m_Star2->Mass();
	double reducedMass					             = (m_Star1->Mass() * m_Star2->Mass()) / totalMass;
	m_OrbitalEnergy 			                     = CalculateOrbitalEnergy(reducedMass, totalMass, m_SemiMajorAxis);
	m_OrbitalEnergyPrev 			                 = m_OrbitalEnergy;

	m_OrbitalAngularMomentum 	                     = CalculateOrbitalAngularMomentum(m_Star1->Mass(), m_Star2->Mass(), m_SemiMajorAxis, m_Eccentricity);
	m_OrbitalAngularMomentumPrev 	                 = m_OrbitalAngularMomentum;

    m_Time                                           = DEFAULT_INITIAL_DOUBLE_VALUE;
	m_Dt                                             = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_TimePrev                                       = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_DCOFormationTime                               = DEFAULT_INITIAL_DOUBLE_VALUE;

    m_aMassLossDiff                                  = DEFAULT_INITIAL_DOUBLE_VALUE;

    m_aMassTransferDiff                              = DEFAULT_INITIAL_DOUBLE_VALUE;

	m_MassTransferTrackerHistory                     = MT_TRACKING::NO_MASS_TRANSFER;
    m_MassTransfer                                   = false;

    m_JLoss                                          = OPTIONS->MassTransferJloss();

	m_FractionAccreted                               = OPTIONS->MassTransferFractionAccreted();

    // Common Envelope
    m_CEDetails.CEEcount                             = 0;
    m_CEDetails.CEEnow                               = false;
    m_CEDetails.doubleCoreCE                         = false;
	m_CEDetails.optimisticCE                         = false;
	m_CEDetails.postCEE.eccentricity                 = DEFAULT_INITIAL_DOUBLE_VALUE;
	m_CEDetails.postCEE.rocheLobe1to2                = DEFAULT_INITIAL_DOUBLE_VALUE;
	m_CEDetails.postCEE.rocheLobe2to1                = DEFAULT_INITIAL_DOUBLE_VALUE;
	m_CEDetails.postCEE.semiMajorAxis                = DEFAULT_INITIAL_DOUBLE_VALUE;
	m_CEDetails.preCEE.eccentricity                  = DEFAULT_INITIAL_DOUBLE_VALUE;
	m_CEDetails.preCEE.rocheLobe1to2                 = DEFAULT_INITIAL_DOUBLE_VALUE;
	m_CEDetails.preCEE.rocheLobe2to1                 = DEFAULT_INITIAL_DOUBLE_VALUE;
	m_CEDetails.preCEE.semiMajorAxis                 = DEFAULT_INITIAL_DOUBLE_VALUE;

    m_Flags.stellarMerger                            = false;
    m_Flags.stellarMergerAtBirth                     = false;

	m_Mass1Final                                     = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_Mass2Final                                     = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_MassEnv1                                       = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_MassEnv2                                       = DEFAULT_INITIAL_DOUBLE_VALUE;

    m_ZetaLobe                                       = DEFAULT_INITIAL_DOUBLE_VALUE;
	m_ZetaStar	                                     = DEFAULT_INITIAL_DOUBLE_VALUE;

    // Initialise other parameters to 0
    m_CosIPrime                                      = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_IPrime                                         = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_TimeToCoalescence                              = DEFAULT_INITIAL_DOUBLE_VALUE;

    m_SupernovaState                                 = SN_STATE::NONE;

    m_Flags.mergesInHubbleTime                       = false;
    m_Unbound                                        = false;

    m_SystemicVelocity                               = Vector3d();
    m_NormalizedOrbitalAngularMomentumVector         = Vector3d();
	m_ThetaE                                         = DEFAULT_INITIAL_DOUBLE_VALUE;
	m_PhiE                                           = DEFAULT_INITIAL_DOUBLE_VALUE;
	m_PsiE                                           = DEFAULT_INITIAL_DOUBLE_VALUE;

	m_SynchronizationTimescale                       = DEFAULT_INITIAL_DOUBLE_VALUE;
	m_CircularizationTimescale                       = DEFAULT_INITIAL_DOUBLE_VALUE;

	// RLOF details
    m_RLOFDetails.experiencedRLOF                    = false;
    m_RLOFDetails.immediateRLOFPostCEE               = false;
    m_RLOFDetails.isRLOF                             = false;
    m_RLOFDetails.simultaneousRLOF                   = false;
    m_RLOFDetails.stableRLOFPostCEE                  = false;

	// RLOF details - properties 1
    m_RLOFDetails.props1.id                          = -1l;

    m_RLOFDetails.props1.stellarType1                = STELLAR_TYPE::NONE;
    m_RLOFDetails.props1.stellarType2                = STELLAR_TYPE::NONE;

    m_RLOFDetails.props1.mass1                       = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_RLOFDetails.props1.mass2                       = DEFAULT_INITIAL_DOUBLE_VALUE;

    m_RLOFDetails.props1.radius1                     = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_RLOFDetails.props1.radius2                     = DEFAULT_INITIAL_DOUBLE_VALUE;

    m_RLOFDetails.props1.starToRocheLobeRadiusRatio1 = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_RLOFDetails.props1.starToRocheLobeRadiusRatio2 = DEFAULT_INITIAL_DOUBLE_VALUE;

    m_RLOFDetails.props1.semiMajorAxis               = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_RLOFDetails.props1.eccentricity                = DEFAULT_INITIAL_DOUBLE_VALUE;
    
    m_RLOFDetails.props1.eventCounter                = DEFAULT_INITIAL_ULONGINT_VALUE;

    m_RLOFDetails.props1.time                        = DEFAULT_INITIAL_DOUBLE_VALUE;
    
    m_RLOFDetails.props1.accretionEfficiency         = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_RLOFDetails.props1.massLossRateFromDonor       = DEFAULT_INITIAL_DOUBLE_VALUE;

    m_RLOFDetails.props1.isRLOF1                     = false;
    m_RLOFDetails.props1.isRLOF2                     = false;

    m_RLOFDetails.props1.isCE                        = false;


	// RLOF details - properties 2
    m_RLOFDetails.props2.id = -1l;

    m_RLOFDetails.props2.stellarType1                = STELLAR_TYPE::NONE;
    m_RLOFDetails.props2.stellarType2                = STELLAR_TYPE::NONE;

    m_RLOFDetails.props2.mass1                       = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_RLOFDetails.props2.mass2                       = DEFAULT_INITIAL_DOUBLE_VALUE;

    m_RLOFDetails.props2.radius1                     = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_RLOFDetails.props2.radius2                     = DEFAULT_INITIAL_DOUBLE_VALUE;

    m_RLOFDetails.props2.starToRocheLobeRadiusRatio1 = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_RLOFDetails.props2.starToRocheLobeRadiusRatio2 = DEFAULT_INITIAL_DOUBLE_VALUE;


    m_RLOFDetails.props2.semiMajorAxis               = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_RLOFDetails.props2.eccentricity                = DEFAULT_INITIAL_DOUBLE_VALUE;
    
    m_RLOFDetails.props2.eventCounter                = DEFAULT_INITIAL_ULONGINT_VALUE;

    m_RLOFDetails.props2.time                        = DEFAULT_INITIAL_DOUBLE_VALUE;
    
    m_RLOFDetails.props2.accretionEfficiency         = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_RLOFDetails.props2.massLossRateFromDonor       = DEFAULT_INITIAL_DOUBLE_VALUE;

    m_RLOFDetails.props2.isRLOF1                     = false;
    m_RLOFDetails.props2.isRLOF2                     = false;

    m_RLOFDetails.props2.isCE                        = false;

    // RLOF details - pre/post-MT props pointers
    m_RLOFDetails.propsPostMT                        = &m_RLOFDetails.props1;
    m_RLOFDetails.propsPreMT                         = &m_RLOFDetails.props2;


    // BeBinary details - properties 1
//    m_BeBinaryDetails.props1.id                      = -1l;
//
//    m_BeBinaryDetails.props1.dt                      = DEFAULT_INITIAL_DOUBLE_VALUE;
//    m_BeBinaryDetails.props1.totalTime               = DEFAULT_INITIAL_DOUBLE_VALUE;
//
//    m_BeBinaryDetails.props1.massNS                  = DEFAULT_INITIAL_DOUBLE_VALUE;
//
//    m_BeBinaryDetails.props1.companionMass           = DEFAULT_INITIAL_DOUBLE_VALUE;
//    m_BeBinaryDetails.props1.companionLuminosity     = DEFAULT_INITIAL_DOUBLE_VALUE;
//    m_BeBinaryDetails.props1.companionTeff           = DEFAULT_INITIAL_DOUBLE_VALUE;
//    m_BeBinaryDetails.props1.companionRadius         = DEFAULT_INITIAL_DOUBLE_VALUE;
//
//    m_BeBinaryDetails.props1.semiMajorAxis           = DEFAULT_INITIAL_DOUBLE_VALUE;
//    m_BeBinaryDetails.props1.eccentricity            = DEFAULT_INITIAL_DOUBLE_VALUE;

    // BeBinary details - properties 2
//    m_BeBinaryDetails.props2.id                      = -1l;
//
//    m_BeBinaryDetails.props2.dt                      = DEFAULT_INITIAL_DOUBLE_VALUE;
//    m_BeBinaryDetails.props2.totalTime               = DEFAULT_INITIAL_DOUBLE_VALUE;
//
//    m_BeBinaryDetails.props2.massNS                  = DEFAULT_INITIAL_DOUBLE_VALUE;
//
//    m_BeBinaryDetails.props2.companionMass           = DEFAULT_INITIAL_DOUBLE_VALUE;
//    m_BeBinaryDetails.props2.companionLuminosity     = DEFAULT_INITIAL_DOUBLE_VALUE;
//    m_BeBinaryDetails.props2.companionTeff           = DEFAULT_INITIAL_DOUBLE_VALUE;
//    m_BeBinaryDetails.props2.companionRadius         = DEFAULT_INITIAL_DOUBLE_VALUE;
//
//    m_BeBinaryDetails.props2.semiMajorAxis           = DEFAULT_INITIAL_DOUBLE_VALUE;
//    m_BeBinaryDetails.props2.eccentricity            = DEFAULT_INITIAL_DOUBLE_VALUE;

    // BeBinary details - current/prev props pointers
//    m_BeBinaryDetails.currentProps                   = &m_BeBinaryDetails.props1;
//    m_BeBinaryDetails.previousProps                  = &m_BeBinaryDetails.props2;

    // pointers

    m_Donor                                          = nullptr;
    m_Accretor                                       = nullptr;

    m_Supernova                                      = nullptr;
    m_Companion                                      = nullptr;
}


/*
 * Determine the value of the requested property of the binary (parameter p_Property)
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
 * This function handles properties of type BINARY_PROPERTY only.
 *
 * This is the function used to retrieve values for properties required to be printed.
 * This allows the composition of the log records to be dynamically modified - this is
 * how we allow users to specify what properties they want recorded in log files.
 *
 * The functional return is the value of the property requested.
 *
 *
 * COMPAS_VARIABLE BinaryPropertyValue(const T_ANY_PROPERTY p_Property) const
 *
 * @param   [IN]    p_Property                  The property for which the value is required
 * @return                                      The value of the requested property
 */
COMPAS_VARIABLE BaseBinaryStar::BinaryPropertyValue(const T_ANY_PROPERTY p_Property) const {

    COMPAS_VARIABLE value;                                                                                              // property value

    BINARY_PROPERTY property = boost::get<BINARY_PROPERTY>(p_Property);                                                 // get the id of the property required

    switch (property) {                                                                                                 // which property?

        case BINARY_PROPERTY::CIRCULARIZATION_TIMESCALE:                            value = CircularizationTimescale();                                         break;
        case BINARY_PROPERTY::COMMON_ENVELOPE_AT_LEAST_ONCE:                        value = CEAtLeastOnce();                                                    break;
        case BINARY_PROPERTY::COMMON_ENVELOPE_EVENT_COUNT:                          value = CommonEnvelopeEventCount();                                         break;
        case BINARY_PROPERTY::UNBOUND:                                              value = Unbound();                                                          break;
        case BINARY_PROPERTY::DOUBLE_CORE_COMMON_ENVELOPE:                          value = DoubleCoreCE();                                                     break;
        case BINARY_PROPERTY::DT:                                                   value = Dt();                                                               break;
        case BINARY_PROPERTY::ECCENTRICITY:                                         value = Eccentricity();                                                     break;
        case BINARY_PROPERTY::ECCENTRICITY_AT_DCO_FORMATION:                        value = EccentricityAtDCOFormation();                                       break;
        case BINARY_PROPERTY::ECCENTRICITY_INITIAL:                                 value = EccentricityInitial();                                              break;
        case BINARY_PROPERTY::ECCENTRICITY_POST_COMMON_ENVELOPE:                    value = EccentricityPostCEE();                                              break;
        case BINARY_PROPERTY::ECCENTRICITY_PRE_SUPERNOVA:                           value = EccentricityPreSN();                                                break;
        case BINARY_PROPERTY::ECCENTRICITY_PRE_COMMON_ENVELOPE:                     value = EccentricityPreCEE();                                               break;
        case BINARY_PROPERTY::ERROR:                                                value = Error();                                                            break;
        case BINARY_PROPERTY::EVOL_STATUS:                                          value = EvolutionStatus();                                                  break;
        case BINARY_PROPERTY::ID:                                                   value = ObjectId();                                                         break;
        case BINARY_PROPERTY::IMMEDIATE_RLOF_POST_COMMON_ENVELOPE:                  value = ImmediateRLOFPostCEE();                                             break;
        case BINARY_PROPERTY::MASS_1_POST_COMMON_ENVELOPE:                          value = Mass1PostCEE();                                                     break;
        case BINARY_PROPERTY::MASS_1_PRE_COMMON_ENVELOPE:                           value = Mass1PreCEE();                                                      break;
        case BINARY_PROPERTY::MASS_2_POST_COMMON_ENVELOPE:                          value = Mass2PostCEE();                                                     break;
        case BINARY_PROPERTY::MASS_2_PRE_COMMON_ENVELOPE:                           value = Mass2PreCEE();                                                      break;
        case BINARY_PROPERTY::MASS_ENV_1:                                           value = MassEnv1();                                                         break;
        case BINARY_PROPERTY::MASS_ENV_2:                                           value = MassEnv2();                                                         break;
        case BINARY_PROPERTY::MASSES_EQUILIBRATED:                                  value = MassesEquilibrated();                                               break;
        case BINARY_PROPERTY::MASSES_EQUILIBRATED_AT_BIRTH:                         value = MassesEquilibratedAtBirth();                                        break;
        case BINARY_PROPERTY::MASS_TRANSFER_TRACKER_HISTORY:                        value = MassTransferTrackerHistory();                                       break;
        case BINARY_PROPERTY::MERGES_IN_HUBBLE_TIME:                                value = MergesInHubbleTime();                                               break;
        case BINARY_PROPERTY::OPTIMISTIC_COMMON_ENVELOPE:                           value = OptimisticCommonEnvelope();                                         break;
        case BINARY_PROPERTY::ORBITAL_ANGULAR_VELOCITY:                             value = OrbitalAngularVelocity();                                           break;
        case BINARY_PROPERTY::ORBITAL_VELOCITY_PRE_SUPERNOVA:                       value = OrbitalVelocityPreSN();                                             break;
        case BINARY_PROPERTY::RADIUS_1_POST_COMMON_ENVELOPE:                        value = Radius1PostCEE();                                                   break;
        case BINARY_PROPERTY::RADIUS_1_PRE_COMMON_ENVELOPE:                         value = Radius1PreCEE();                                                    break;
        case BINARY_PROPERTY::RADIUS_2_POST_COMMON_ENVELOPE:                        value = Radius2PostCEE();                                                   break;
        case BINARY_PROPERTY::RADIUS_2_PRE_COMMON_ENVELOPE:                         value = Radius2PreCEE();                                                    break;
        case BINARY_PROPERTY::RANDOM_SEED:                                          value = RandomSeed();                                                       break;
        case BINARY_PROPERTY::RLOF_ACCRETION_EFFICIENCY:                            value = RLOFDetails().propsPostMT->accretionEfficiency;                     break;
        case BINARY_PROPERTY::RLOF_MASS_LOSS_RATE:                                  value = RLOFDetails().propsPostMT->massLossRateFromDonor;                   break;
        case BINARY_PROPERTY::RLOF_POST_MT_COMMON_ENVELOPE:                         value = RLOFDetails().propsPostMT->isCE;                                    break;
        case BINARY_PROPERTY::RLOF_POST_MT_ECCENTRICITY:                            value = RLOFDetails().propsPostMT->eccentricity;                            break;
        case BINARY_PROPERTY::RLOF_POST_MT_EVENT_COUNTER:                           value = RLOFDetails().propsPostMT->eventCounter;                            break;
        case BINARY_PROPERTY::RLOF_POST_MT_ID:                                      value = RLOFDetails().propsPostMT->id;                                      break;
        case BINARY_PROPERTY::RLOF_POST_MT_SEMI_MAJOR_AXIS:                         value = RLOFDetails().propsPostMT->semiMajorAxis;                           break;
        case BINARY_PROPERTY::RLOF_POST_MT_STAR1_MASS:                              value = RLOFDetails().propsPostMT->mass1;                                   break;
        case BINARY_PROPERTY::RLOF_POST_MT_STAR2_MASS:                              value = RLOFDetails().propsPostMT->mass2;                                   break;
        case BINARY_PROPERTY::RLOF_POST_MT_STAR1_RADIUS:                            value = RLOFDetails().propsPostMT->radius1;                                 break;
        case BINARY_PROPERTY::RLOF_POST_MT_STAR2_RADIUS:                            value = RLOFDetails().propsPostMT->radius2;                                 break;
        case BINARY_PROPERTY::RLOF_POST_MT_STAR1_RLOF:                              value = RLOFDetails().propsPostMT->isRLOF1;                                 break;
        case BINARY_PROPERTY::RLOF_POST_MT_STAR2_RLOF:                              value = RLOFDetails().propsPostMT->isRLOF2;                                 break;
        case BINARY_PROPERTY::RLOF_POST_MT_STAR1_STELLAR_TYPE:                      value = RLOFDetails().propsPostMT->stellarType1;                            break;
        case BINARY_PROPERTY::RLOF_POST_MT_STAR1_STELLAR_TYPE_NAME:                 value = STELLAR_TYPE_LABEL.at(RLOFDetails().propsPostMT->stellarType1);     break;
        case BINARY_PROPERTY::RLOF_POST_MT_STAR2_STELLAR_TYPE:                      value = RLOFDetails().propsPostMT->stellarType2;                            break;
        case BINARY_PROPERTY::RLOF_POST_MT_STAR2_STELLAR_TYPE_NAME:                 value = STELLAR_TYPE_LABEL.at(RLOFDetails().propsPostMT->stellarType2);     break;
        case BINARY_PROPERTY::RLOF_POST_STEP_STAR_TO_ROCHE_LOBE_RADIUS_RATIO_1:     value = RLOFDetails().propsPostMT->starToRocheLobeRadiusRatio1;             break;
        case BINARY_PROPERTY::RLOF_POST_STEP_STAR_TO_ROCHE_LOBE_RADIUS_RATIO_2:     value = RLOFDetails().propsPostMT->starToRocheLobeRadiusRatio2;             break;
        case BINARY_PROPERTY::RLOF_PRE_MT_ECCENTRICITY:                             value = RLOFDetails().propsPreMT->eccentricity;                             break;
        case BINARY_PROPERTY::RLOF_PRE_MT_SEMI_MAJOR_AXIS:                          value = RLOFDetails().propsPreMT->semiMajorAxis;                            break;
        case BINARY_PROPERTY::RLOF_PRE_MT_STAR1_MASS:                               value = RLOFDetails().propsPreMT->mass1;                                    break;
        case BINARY_PROPERTY::RLOF_PRE_MT_STAR2_MASS:                               value = RLOFDetails().propsPreMT->mass2;                                    break;
        case BINARY_PROPERTY::RLOF_PRE_MT_STAR1_RADIUS:                             value = RLOFDetails().propsPreMT->radius1;                                  break;
        case BINARY_PROPERTY::RLOF_PRE_MT_STAR2_RADIUS:                             value = RLOFDetails().propsPreMT->radius2;                                  break;
        case BINARY_PROPERTY::RLOF_PRE_MT_STAR1_RLOF:                               value = RLOFDetails().propsPreMT->isRLOF1;                                  break;
        case BINARY_PROPERTY::RLOF_PRE_MT_STAR2_RLOF:                               value = RLOFDetails().propsPreMT->isRLOF2;                                  break;
        case BINARY_PROPERTY::RLOF_PRE_MT_STAR1_STELLAR_TYPE:                       value = RLOFDetails().propsPreMT->stellarType1;                             break;
        case BINARY_PROPERTY::RLOF_PRE_MT_STAR1_STELLAR_TYPE_NAME:                  value = STELLAR_TYPE_LABEL.at(RLOFDetails().propsPreMT->stellarType1);      break;
        case BINARY_PROPERTY::RLOF_PRE_MT_STAR2_STELLAR_TYPE:                       value = RLOFDetails().propsPreMT->stellarType2;                             break;
        case BINARY_PROPERTY::RLOF_PRE_MT_STAR2_STELLAR_TYPE_NAME:                  value = STELLAR_TYPE_LABEL.at(RLOFDetails().propsPreMT->stellarType2);      break;
        case BINARY_PROPERTY::RLOF_PRE_STEP_STAR_TO_ROCHE_LOBE_RADIUS_RATIO_1:      value = RLOFDetails().propsPreMT->starToRocheLobeRadiusRatio1;              break;
        case BINARY_PROPERTY::RLOF_PRE_STEP_STAR_TO_ROCHE_LOBE_RADIUS_RATIO_2:      value = RLOFDetails().propsPreMT->starToRocheLobeRadiusRatio2;              break;
        case BINARY_PROPERTY::RLOF_SECONDARY_POST_COMMON_ENVELOPE:                  value = RLOFSecondaryPostCEE();                                             break;
        case BINARY_PROPERTY::RLOF_TIME_POST_MT:                                    value = RLOFDetails().propsPreMT->time;                                     break;
        case BINARY_PROPERTY::RLOF_TIME_PRE_MT:                                     value = RLOFDetails().propsPreMT->timePrev;                                 break;
        case BINARY_PROPERTY::ROCHE_LOBE_RADIUS_1:                                  value = RocheLobeRadius1();                                                 break;
        case BINARY_PROPERTY::ROCHE_LOBE_RADIUS_1_POST_COMMON_ENVELOPE:             value = RocheLobe1to2PostCEE();                                             break;
        case BINARY_PROPERTY::ROCHE_LOBE_RADIUS_1_PRE_COMMON_ENVELOPE:              value = RocheLobe1to2PreCEE();                                              break;
        case BINARY_PROPERTY::ROCHE_LOBE_RADIUS_2:                                  value = RocheLobeRadius2();                                                 break;
        case BINARY_PROPERTY::ROCHE_LOBE_RADIUS_2_POST_COMMON_ENVELOPE:             value = RocheLobe2to1PostCEE();                                             break;
        case BINARY_PROPERTY::ROCHE_LOBE_RADIUS_2_PRE_COMMON_ENVELOPE:              value = RocheLobe2to1PreCEE();                                              break;
        case BINARY_PROPERTY::SEMI_MAJOR_AXIS_AT_DCO_FORMATION:                     value = SemiMajorAxisAtDCOFormation();                                      break;
        case BINARY_PROPERTY::SEMI_MAJOR_AXIS_INITIAL:                              value = SemiMajorAxisInitial();                                             break;
        case BINARY_PROPERTY::SEMI_MAJOR_AXIS_POST_COMMON_ENVELOPE:                 value = SemiMajorAxisPostCEE();                                             break;
        case BINARY_PROPERTY::SEMI_MAJOR_AXIS_PRE_SUPERNOVA:                        value = SemiMajorAxisPreSN();                                               break;
        case BINARY_PROPERTY::SEMI_MAJOR_AXIS_PRE_SUPERNOVA_RSOL:                   value = SemiMajorAxisPreSN() * AU_TO_RSOL;                                  break;
        case BINARY_PROPERTY::SEMI_MAJOR_AXIS_PRE_COMMON_ENVELOPE:                  value = SemiMajorAxisPreCEE();                                              break;
        case BINARY_PROPERTY::SEMI_MAJOR_AXIS:                                      value = SemiMajorAxis();                                                    break;
        case BINARY_PROPERTY::SEMI_MAJOR_AXIS_RSOL:                                 value = SemiMajorAxis() * AU_TO_RSOL;                                       break;
        case BINARY_PROPERTY::SIMULTANEOUS_RLOF:                                    value = SimultaneousRLOF();                                                 break;
        case BINARY_PROPERTY::STABLE_RLOF_POST_COMMON_ENVELOPE:                     value = StableRLOFPostCEE();                                                break;
        case BINARY_PROPERTY::STAR_TO_ROCHE_LOBE_RADIUS_RATIO_1:                    value = StarToRocheLobeRadiusRatio1();                                      break;
        case BINARY_PROPERTY::STAR_TO_ROCHE_LOBE_RADIUS_RATIO_2:                    value = StarToRocheLobeRadiusRatio2();                                      break;
        case BINARY_PROPERTY::STELLAR_MERGER:                                       value = StellarMerger();                                                    break;
        case BINARY_PROPERTY::STELLAR_MERGER_AT_BIRTH:                              value = StellarMergerAtBirth();                                             break;
        case BINARY_PROPERTY::STELLAR_TYPE_1_POST_COMMON_ENVELOPE:                  value = StellarType1PostCEE();                                              break;
        case BINARY_PROPERTY::STELLAR_TYPE_1_PRE_COMMON_ENVELOPE:                   value = StellarType1PreCEE();                                               break;
        case BINARY_PROPERTY::STELLAR_TYPE_2_POST_COMMON_ENVELOPE:                  value = StellarType2PostCEE();                                              break;
        case BINARY_PROPERTY::STELLAR_TYPE_2_PRE_COMMON_ENVELOPE:                   value = StellarType2PreCEE();                                               break;
        case BINARY_PROPERTY::STELLAR_TYPE_NAME_1_POST_COMMON_ENVELOPE:             value = STELLAR_TYPE_LABEL.at(StellarType1PostCEE());                       break;
        case BINARY_PROPERTY::STELLAR_TYPE_NAME_1_PRE_COMMON_ENVELOPE:              value = STELLAR_TYPE_LABEL.at(StellarType1PreCEE());                        break;
        case BINARY_PROPERTY::STELLAR_TYPE_NAME_2_POST_COMMON_ENVELOPE:             value = STELLAR_TYPE_LABEL.at(StellarType2PostCEE());                       break;
        case BINARY_PROPERTY::STELLAR_TYPE_NAME_2_PRE_COMMON_ENVELOPE:              value = STELLAR_TYPE_LABEL.at(StellarType2PreCEE());                        break;
        case BINARY_PROPERTY::SUPERNOVA_ORBIT_INCLINATION_ANGLE:                    value = SN_OrbitInclinationAngle();                                         break;
        case BINARY_PROPERTY::SUPERNOVA_ORBIT_INCLINATION_VECTOR_X:                 value = SN_OrbitInclinationVectorX();                                       break;
        case BINARY_PROPERTY::SUPERNOVA_ORBIT_INCLINATION_VECTOR_Y:                 value = SN_OrbitInclinationVectorY();                                       break;
        case BINARY_PROPERTY::SUPERNOVA_ORBIT_INCLINATION_VECTOR_Z:                 value = SN_OrbitInclinationVectorZ();                                       break;
        case BINARY_PROPERTY::SUPERNOVA_STATE:                                      value = SN_State();                                                         break;
        case BINARY_PROPERTY::SYNCHRONIZATION_TIMESCALE:                            value = SynchronizationTimescale();                                         break;
        case BINARY_PROPERTY::SYSTEMIC_SPEED:                                       value = SystemicSpeed();                                                    break;
        case BINARY_PROPERTY::TIME:                                                 value = Time();                                                             break;
        case BINARY_PROPERTY::TIME_TO_COALESCENCE:                                  value = TimeToCoalescence();                                                break;
        case BINARY_PROPERTY::TOTAL_ANGULAR_MOMENTUM:                               value = TotalAngularMomentum();                                             break;
        case BINARY_PROPERTY::TOTAL_ENERGY:                                         value = TotalEnergy();                                                      break;
        case BINARY_PROPERTY::ZETA_LOBE:                                            value = ZetaLobe();                                                         break;
        case BINARY_PROPERTY::ZETA_STAR:                                            value = ZetaStar();                                                         break;

        default:                                                                                                        // unexpected binary property
            // the only ways this can happen are if someone added a binary property (into BINARY_PROPERTY),
            // or allowed users to specify a binary property (via the logfile definitions file), and it isn't
            // accounted for in this code.  We should not default here, with or without a warning - this is a
            // code defect, so we flag it as an error and that will result in termination of the evolution of
            // the binary.
            // The correct fix for this is to add code for the missing property, or prevent it from being 
            // specified in the logfile definitions file.

            THROW_ERROR(ERROR::UNEXPECTED_BINARY_PROPERTY);                                                             // throw error
    }

    return value;
}


/*
 * Determine the value of the requested property of the binary (parameter p_Property)
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
 * This function handles properties of type:
 * 
 *    STAR_1_PROPERTY, STAR_2_PROPERTY, SUPERNOVA_PROPERTY, COMPANION_PROPERTY, BINARY_PROPERTY, PROGRAM_OPTION
 * 
 * only - anything else will result in an error being thrown and the evolution of the star (or binary)
 * terminated.
 * 
 * This function calls the appropriate helper function to retrieve the value.
 *
 * This is the function used to retrieve values for properties required to be printed.
 * This allows the composition of the log records to be dynamically modified - this is
 * how we allow users to specify what properties they want recorded in log files.
 *
 * The functional return is the value of the property requested. 
 *
 *
 * COMPAS_VARIABLE PropertyValue(const T_ANY_PROPERTY p_Property) const
 *
 * @param   [IN]    p_Property                  The property for which the value is required
 * @return                                      The value of the requested property
 */
COMPAS_VARIABLE BaseBinaryStar::PropertyValue(const T_ANY_PROPERTY p_Property) const {

    COMPAS_VARIABLE value;                                                                                              // property value

    switch (boost::apply_visitor(VariantPropertyType(), p_Property)) {                                                  // which property type?

        case ANY_PROPERTY_TYPE::T_BINARY_PROPERTY:                                                                      // BSE binary star property
            value = BinaryPropertyValue(p_Property);                                                                    // get the value
            break;

        case ANY_PROPERTY_TYPE::T_STAR_1_PROPERTY:                                                                      // star 1 of BSE binary star property
            if (m_Star1) value = m_Star1->StellarPropertyValue(p_Property);                                             // if have pointer to primary, get the value
            break;

        case ANY_PROPERTY_TYPE::T_STAR_2_PROPERTY:                                                                      // star 2 of BSE binary star property
            if (m_Star2) value = m_Star2->StellarPropertyValue(p_Property);                                             // if have pointer to secondary, get the value
            break;

        case ANY_PROPERTY_TYPE::T_SUPERNOVA_PROPERTY:                                                                   // supernova star of BSE binary star property
            if (m_Supernova) value = m_Supernova->StellarPropertyValue(p_Property);                                     // if have pointer to supernova, get the value
            break;

        case ANY_PROPERTY_TYPE::T_COMPANION_PROPERTY:                                                                   // companion star of BSE binary star property
            if (m_Companion) value = m_Companion->StellarPropertyValue(p_Property);                                     // if have pointer to companion, get the value
            break;

        case ANY_PROPERTY_TYPE::T_PROGRAM_OPTION:                                                                       // program option
            value = OPTIONS->OptionValue(p_Property);                                                                   // get the value
            break;

        default:                                                                                                        // unexpected binary property type
            // the only ways this can happen are if someone added a stellar type property (into ANY_PROPERTY_TYPE)
            // and it isn't accounted for in this code, or if there is a defect in the code that causes
            // this function to be called with a bad parameter.  We should not default here, with or without a
            // warning - this is a code defect, so we flag it as an error and that will result in termination of
            // the evolution of the binary.
            // The correct fix for this is to add code for the missing property type or find and fix the code defect.

            THROW_ERROR(ERROR::UNEXPECTED_BINARY_PROPERTY_TYPE);                                                        // throw error
    }

    return value;
}


/*
 * Determines if the binary contains only one star which is one of a list of stellar types passed
 *
 *
 * bool HasOnlyOneOf(STELLAR_TYPE_LIST p_List)
 *
 * @param   [IN]    p_List                      List of stellar types
 * @return                                      Boolean - true if only one of the stars of the binary is in list, false if neither or both
 */
bool BaseBinaryStar::HasOnlyOneOf(STELLAR_TYPE_LIST p_List) const {
    int matchCount = 0;
    for (auto elem: p_List) {
        if (m_Star1->StellarType() == elem) matchCount++;
        if (m_Star2->StellarType() == elem) matchCount++;
    }
    return matchCount == 1;
}


/*
 * Determines if the binary contains at least one star which is one of a list of stellar types passed
 *
 *
 * bool HasOneOf(STELLAR_TYPE_LIST p_List)
 *
 * @param   [IN]    p_List                      List of stellar types
 * @return                                      Boolean - true if one of the stars of the binary is in list, false if not
 */
bool BaseBinaryStar::HasOneOf(STELLAR_TYPE_LIST p_List) const {
    for (auto elem: p_List) {
        if ((m_Star1->StellarType() == elem) || (m_Star2->StellarType() == elem)) return true;
    }
	return false;
}


/*
 * Determines if the binary contains two stars from the list of stellar types passed
 *
 *
 * bool HasTwoOf(STELLAR_TYPE_LIST p_List)
 *
 * @param   [IN]    p_List                      List of stellar types
 * @return                                      Boolean - true if both of the stars of the binary are in list, false if not
 */
bool BaseBinaryStar::HasTwoOf(STELLAR_TYPE_LIST p_List) const {
    int matchCount = 0;
    for (auto elem: p_List) {
        if (m_Star1->StellarType() == elem) matchCount++;
        if (m_Star2->StellarType() == elem) matchCount++;
        if (matchCount > 1) return true;
    }
	return false;
}

/*
 * Determines if the binary is a high-mass XRB candidate (one compact object with a companion at >80% Roche lobe filling)
 *
 *
 * bool IsHMXRBinary()
 * @return                                      Boolean - true if the binary is a HMXRB candidate
 *
 */
bool BaseBinaryStar::IsHMXRBinary() const {
    if (HasOnlyOneOf({STELLAR_TYPE::NEUTRON_STAR, STELLAR_TYPE::BLACK_HOLE})){
        if (m_Star1->StellarType() < STELLAR_TYPE::NEUTRON_STAR && utils::Compare(StarToRocheLobeRadiusRatio1(), MIN_HMXRB_STAR_TO_ROCHE_LOBE_RADIUS_RATIO) > 0) return true;
        if (m_Star2->StellarType() < STELLAR_TYPE::NEUTRON_STAR && utils::Compare(StarToRocheLobeRadiusRatio2(), MIN_HMXRB_STAR_TO_ROCHE_LOBE_RADIUS_RATIO) > 0) return true;
    }
    return false;
}


/*
 * Write RLOF parameters to RLOF logfile if RLOF printing is enabled and at least one of the stars is in RLOF
 * and / or HMXRBs are being printed and IsHMXRBinary is true
 *
 *
 * bool PrintRLOFParameters(const RLOF_RECORD_TYPE p_RecordType)
 * 
 * @param   [IN]    p_RecordType                Record type to be written
 * @return                                      Boolean status (true = success, false = failure)
 * 
 */
bool BaseBinaryStar::PrintRLOFParameters(const RLOF_RECORD_TYPE p_RecordType) {

    bool ok = true;

    if (!OPTIONS->RLOFPrinting()) return ok;                            // do not print if printing option off

    StashRLOFProperties(MT_TIMING::POST_MT);                            // stash properties immediately post-Mass Transfer 

    if (m_Star1->IsRLOF() || m_Star2->IsRLOF()) {                       // print if either star is in RLOF
        m_RLOFDetails.propsPostMT->eventCounter += 1;                   // every time we print a MT event happened, increment counter
        ok = LOGGING->LogRLOFParameters(this, p_RecordType);            // yes - write to log file
    }

    if (OPTIONS->HMXRBinaries()) {
        if (IsHMXRBinary()) {                                           // print if star is HMXRB candidate
            ok = LOGGING->LogRLOFParameters(this, p_RecordType); 
        }
    }

    return ok;
}


/*
 * Squirrel RLOF properties away
 *
 * Various binary property values are stashed into either the m_RLOFDetails.propsPreMT or 
 * m_RLOFDetails.propsPostMT struct for use/printing later. 
 * The switch is so that pre-MT props store the binary state immediately before EvaluateBinary(),
 * to avoid recording problems when a stellar type changes twice in one timestep.
 *
 * void StashRLOFProperties()
 *
 * @param   [IN]    p_Which                     MT_TIMING (PRE_MT or POST_MT)
 */
void BaseBinaryStar::StashRLOFProperties(const MT_TIMING p_Which) {

    if (!OPTIONS->RLOFPrinting()) return;                                                       // nothing to do

    // set whether to update pre-MT or post-MT parameters depending on input argument
    RLOFPropertiesT* rlofPropertiesToReset = (p_Which == MT_TIMING::PRE_MT) ? m_RLOFDetails.propsPreMT : m_RLOFDetails.propsPostMT;

    // update properties for appropriate timestep
    rlofPropertiesToReset->id                          = m_ObjectId;
    rlofPropertiesToReset->mass1                       = m_Star1->Mass();
    rlofPropertiesToReset->mass2                       = m_Star2->Mass();
    rlofPropertiesToReset->radius1                     = m_Star1->Radius();
    rlofPropertiesToReset->radius2                     = m_Star2->Radius();
    rlofPropertiesToReset->starToRocheLobeRadiusRatio1 = StarToRocheLobeRadiusRatio1();
    rlofPropertiesToReset->starToRocheLobeRadiusRatio2 = StarToRocheLobeRadiusRatio2();
    rlofPropertiesToReset->stellarType1                = m_Star1->StellarType();
    rlofPropertiesToReset->stellarType2                = m_Star2->StellarType();
    rlofPropertiesToReset->eccentricity                = m_Eccentricity;
    rlofPropertiesToReset->semiMajorAxis               = m_SemiMajorAxis * AU_TO_RSOL;          // semi-major axis - change units to Rsol
    rlofPropertiesToReset->time                        = m_Time;
    rlofPropertiesToReset->timePrev                    = m_TimePrev;
    rlofPropertiesToReset->isRLOF1                     = m_Star1->IsRLOF();
    rlofPropertiesToReset->isRLOF2                     = m_Star2->IsRLOF();
    rlofPropertiesToReset->isCE                        = m_CEDetails.CEEnow;
    rlofPropertiesToReset->massLossRateFromDonor       = m_MassLossRateInRLOF;
    rlofPropertiesToReset->accretionEfficiency         = m_FractionAccreted;
}


/*
 * Squirrel BeBinaries properties away
 *
 * Various binary property values are stashed into the m_BeBinaryDetails.currentProps struct for use/printing later
 * The existing m_BeBinaryDetails.currentProps struct is copied to the m_BeBinaryDetails.previousProps struct first
 * (actually there is no copying - just switch pointers...)
 *
 *
 * void StashBeBinaryProperties()
 */
//void BaseBinaryStar::StashBeBinaryProperties() {
//
//    if (!OPTIONS->BeBinaries() || !IsBeBinary()) return;                                                            // nothing to do;
//
//    // switch previous<->current (preserves existing current as (new) previous)
//    BeBinaryPropertiesT* tmp        = m_BeBinaryDetails.previousProps;                                              // save pointer to existing previous props
//    m_BeBinaryDetails.previousProps = m_BeBinaryDetails.currentProps;                                               // existing current props become new previous props (values will be preserved)
//    m_BeBinaryDetails.currentProps  = tmp;                                                                          // new current props points at existing previous (values will be replaced)
//
//    // now save (new) current
//    m_BeBinaryDetails.currentProps->id            = m_ObjectId;                                                      // object id
//    m_BeBinaryDetails.currentProps->dt            = m_Dt;                                                            // timestep
//    m_BeBinaryDetails.currentProps->totalTime     = m_BeBinaryDetails.previousProps->dt + m_Dt;                      // total time - accumulate, don't just replace
//    m_BeBinaryDetails.currentProps->semiMajorAxis = m_SemiMajorAxis * AU_TO_RSOL;                                    // semi-major axis - change units to Rsol
//    m_BeBinaryDetails.currentProps->eccentricity  = m_Eccentricity;                                                  // eccentricity
//
//    BinaryConstituentStar* neutronStar   = m_Star1->IsOneOf({ STELLAR_TYPE::NEUTRON_STAR }) ? m_Star1 : m_Star2;    // pointer to neutron star
//    BinaryConstituentStar* companionStar = m_Star1->IsOneOf({ STELLAR_TYPE::NEUTRON_STAR }) ? m_Star2 : m_Star1;    // pointer to companion
//
//    m_BeBinaryDetails.currentProps->massNS              = neutronStar->Mass();                                      // neutron star mass
//    m_BeBinaryDetails.currentProps->companionMass       = companionStar->Mass();                                    // companion mass
//    m_BeBinaryDetails.currentProps->companionLuminosity = companionStar->Luminosity();                              // companion luminosity
//    m_BeBinaryDetails.currentProps->companionTeff       = companionStar->Temperature();                             // companion temperature
//    m_BeBinaryDetails.currentProps->companionRadius     = companionStar->Radius();                                  // companion radius
//}


/*
 * Calculate (or set) pre common envelope values for the binary:
 *
 *    m_CommonEnvelopeDetails.preCEE.eccentricity
 *    m_CommonEnvelopeDetails.preCEE.semiMajorAxis
 *    m_CommonEnvelopeDetails.preCEE.rocheLobe1to2
 *    m_CommonEnvelopeDetails.preCEE.rocheLobe2to1
 *
 *
 * void SetPreCEEValues(const double p_SemiMajorAxis,
 *                      const double p_Eccentricity,
 *                      const double p_RocheLobe1to2,
 *                      const double p_RocheLobe2to1)
 *
 * @param   [IN]    p_SemiMajorAxis             pre CEE semi-major axis in AU
 * @param   [IN]    p_Eccentricity              pre CEE eccentricity
 * @param   [IN]    p_RocheLobe1to2             pre CEE Roche Lobe radius in AU as seen by star1
 * @param   [IN]    p_RocheLobe2to1             pre CEE Roche Lobe radius in AU as seen by star2
 */
void BaseBinaryStar::SetPreCEEValues(const double p_SemiMajorAxis,
                                     const double p_Eccentricity,
                                     const double p_RocheLobe1to2,
                                     const double p_RocheLobe2to1) {

	m_CEDetails.preCEE.semiMajorAxis = p_SemiMajorAxis;
	m_CEDetails.preCEE.eccentricity  = p_Eccentricity;
	m_CEDetails.preCEE.rocheLobe1to2 = p_RocheLobe1to2;
	m_CEDetails.preCEE.rocheLobe2to1 = p_RocheLobe2to1;
}


/*
 * Calculate (or set) post common envelope values for the binary:
 *
 *    m_CommonEnvelopeDetails.postCEE.eccentricity
 *    m_CommonEnvelopeDetails.postCEE.semiMajorAxis
 *    m_CommonEnvelopeDetails.postCEE.rocheLobe1to2
 *    m_CommonEnvelopeDetails.postCEE.rocheLobe2to1
 *    m_RLOFDetails.immediateRLOFPostCEE
 *
 *
 * void SetPostCEEValues(const double p_SemiMajorAxis,
 *                       const double p_Eccentricity,
 *                       const double p_RocheLobe1to2,
 *                       const double p_RocheLobe2to1)
 *
 * @param   [IN]    p_SemiMajorAxis             post CEE semi-major axis in AU
 * @param   [IN]    p_Eccentricity              post CEE eccentricity
 * @param   [IN]    p_RocheLobe1to2             post CEE Roche Lobe radius in AU as seen by star1
 * @param   [IN]    p_RocheLobe2to1             post CEE Roche Lobe radius in AU as seen by star2
 */
void BaseBinaryStar::SetPostCEEValues(const double p_SemiMajorAxis,
                                      const double p_Eccentricity,
                                      const double p_RocheLobe1to2,
                                      const double p_RocheLobe2to1) {

	m_CEDetails.postCEE.semiMajorAxis = p_SemiMajorAxis;
    m_CEDetails.postCEE.eccentricity  = p_Eccentricity;
	m_CEDetails.postCEE.rocheLobe1to2 = p_RocheLobe1to2;
	m_CEDetails.postCEE.rocheLobe2to1 = p_RocheLobe2to1;

    if (utils::Compare(m_Star1->RadiusPostCEE(), m_CEDetails.postCEE.rocheLobe1to2) >= 0 ||         // Check for RLOF immediately after the CEE
        utils::Compare(m_Star2->RadiusPostCEE(), m_CEDetails.postCEE.rocheLobe2to1) >= 0) {
        m_RLOFDetails.immediateRLOFPostCEE = true;
    }
}


/*
 * Calculate the time to coalescence for a binary with arbitrary eccentricity
 *
 * Mandel 2021 https://iopscience.iop.org/article/10.3847/2515-5172/ac2d35, eq 5
 * 
 * Accurate to within 3% over the full range of initial eccentricities up to 0.99999
 * Will return time = 0.0 for eccentricities < 0.0 and >= 1.0
 *
 *
 * double CalculateTimeToCoalescence(const double p_SemiMajorAxis,
 *                                   const double p_Eccentricity,
 *                                   const double p_Mass1,
 *                                   const double p_Mass2)
 *
 * @param   [IN]    p_SemiMajorAxis             Initial semi-major axis in SI units
 * @param   [IN]    p_Eccentricity              Initial eccentricity
 * @param   [IN]    p_Mass1                     Primary mass in SI units
 * @param   [IN]    p_Mass2                     Secondary mass in SI units
 * @return                                      Time to coalescence in SI units (s): returns 0.0 if p_Eccentricity < 0 or p_Eccentricity >= 1
 */
double BaseBinaryStar::CalculateTimeToCoalescence(const double p_SemiMajorAxis,
                                                  const double p_Eccentricity,
                                                  const double p_Mass1,
                                                  const double p_Mass2) const {

    if (p_Eccentricity < 0.0 || p_Eccentricity >= 1.0) return 0.0;                                      // save some cpu cycles...

    // pow() is slow - use multiplication where possible

    // calculate time for a circular binary to merge - Mandel 2021, eq 2
    double numerator = 5.0 * C * C * C * C * C * p_SemiMajorAxis * p_SemiMajorAxis * p_SemiMajorAxis * p_SemiMajorAxis;
    double denominator = 256.0 * G * G * G * p_Mass1 * p_Mass2 * (p_Mass1 + p_Mass2);

    double tC = numerator / denominator;                                                                // time for a circular binary to merge

    if (utils::Compare(p_Eccentricity, 0.0) > 0) {                                                      // eccentricity > 0.0?
                                                                                                        // yes - not circular
        // calculate time for eccentric binary to merge - Mandel 2021, eq 5
        double e0     = p_Eccentricity;
        double e0_10  = e0 * e0 * e0 * e0 * e0 * e0 * e0 * e0 * e0 * e0;
        double e0_20  = e0_10 * e0_10;
        double e0_100 = e0_10 * e0_10 * e0_10 * e0_10 * e0_10 * e0_10 * e0_10 * e0_10 * e0_10 * e0_10;
        double f      = 1.0 - (e0 * e0);
        double f_3    = f * f * f;
    
        tC = f <= 0.0 ? 0.0 : tC * (1.0 + 0.27 * e0_10 + 0.33 * e0_20 + 0.2 * e0_100) * f_3 * std::sqrt(f);  // check f <= 0.0 just in case a rounding error hurts us
    }

    return tC;
}


/*
 * Resolve coalescence of the binary
 *
 * Calculates:
 *
 *   - time to coaslescence
 *   - whether the binary merges within hubble time
 *
 * Records details of binaries that merge within hubble time
 *
 * void ResolveCoalescence()
 */
void BaseBinaryStar::ResolveCoalescence() {

    // Calculate the time for the binary to coalesce due to emission of gravitational radiation.

    // define DCO formation to be now
    m_SemiMajorAxisAtDCOFormation = m_SemiMajorAxis;
    m_EccentricityAtDCOFormation  = m_Eccentricity;

    double tC           = CalculateTimeToCoalescence(m_SemiMajorAxis * AU, m_Eccentricity, m_Star1->Mass() * MSOL_TO_KG, m_Star2->Mass() * MSOL_TO_KG);
    m_TimeToCoalescence = (tC / SECONDS_IN_YEAR) * YEAR_TO_MYR;                                                                                 // coalescence time in Myr

    if (utils::Compare(tC, HUBBLE_TIME) < 0) {                                                                                                  // shorter than HubbleTime
        m_Flags.mergesInHubbleTime = true;
    }
    else {
        m_Flags.mergesInHubbleTime = false;
    }

    if (!IsUnbound()) {
        (void)PrintDoubleCompactObjects();                                                                                                      // print (log) double compact object details
    }
}


/*
 * Calculate the change in eccentricity based on secular equations for tidal evolution given the tidal Love number
 * Zahn, 1977, Eq. (3.7)
 *
 *
 * double BaseBinaryStar::CalculateDEccentricityTidalDt(const DBL_DBL_DBL_DBL p_ImKlm, const BinaryConstituentStar* p_Star)
 *
 * @param   [IN]    p_ImKlm                     Imaginary [(1,0), (1,2), (2,2), (3,2)] components of the potential tidal Love number of star (unitless)
 * @param   [IN]    p_Star                      Star for which impact on eccentricity is to be calculated
 * @return                                      Change in Eccentricity for binary (1/yr)
 */    
double BaseBinaryStar::CalculateDEccentricityTidalDt(const DBL_DBL_DBL_DBL p_ImKlm, const BinaryConstituentStar* p_Star) {
    
    double massStar      = p_Star->Mass();
    double radiusStar    = p_Star->Radius();
    double massCompanion = p_Star == m_Star1 ? m_Star2->Mass() : m_Star1->Mass();

    double ImK10, ImK12, ImK22, ImK32;
    std::tie(ImK10, ImK12, ImK22, ImK32) = p_ImKlm;

    double R1_AU       = radiusStar * RSOL_TO_AU;
    double R1_over_a   = R1_AU / m_SemiMajorAxis;
    double R1_over_a_8 = R1_over_a * R1_over_a * R1_over_a * R1_over_a * R1_over_a * R1_over_a * R1_over_a * R1_over_a;

    return -(3.0 / 4.0) * (m_Eccentricity / m_Omega) * (1.0 + (massCompanion / massStar)) * (G_AU_Msol_yr * massCompanion / R1_AU / R1_AU / R1_AU) * R1_over_a_8 * ((3.0 * ImK10 / 2.0) - (ImK12 / 4.0) - ImK22 + (49.0 * ImK32 / 4.0));
}


/*
 * Calculate the change in spin based on secular equations for tidal evolution given the tidal Love number
 * Zahn, 1977, Eq. (3.8)
 *
 *
 * double BaseBinaryStar::CalculateDOmegaTidalDt(const DBL_DBL_DBL_DBL p_ImKlm, const BinaryConstituentStar* p_Star)
 *
 * @param   [IN]    p_ImKlm                     Imaginary [(1,0), (1,2), (2,2), (3,2)] components of the potential tidal Love number of star (unitless)
 * @param   [IN]    p_Star                      Star for which impact on spin is to be calculated
 * @return                                      Change in Omega for star (1/yr/yr)
 */    
double BaseBinaryStar::CalculateDOmegaTidalDt(const DBL_DBL_DBL_DBL p_ImKlm, const BinaryConstituentStar* p_Star) {
 
    double MoIstar       = p_Star->CalculateMomentOfInertiaAU();
    double radiusStar    = p_Star->Radius();
    double massCompanion = p_Star == m_Star1 ? m_Star2->Mass() : m_Star1->Mass();

    double ImK10, ImK12, ImK22, ImK32;
    std::tie(ImK10, ImK12, ImK22, ImK32) = p_ImKlm;

    double R1_AU       = radiusStar * RSOL_TO_AU;
    double R1_over_a   = R1_AU / m_SemiMajorAxis;
    double R1_over_a_6 = R1_over_a * R1_over_a * R1_over_a * R1_over_a * R1_over_a * R1_over_a;

    return (3.0 / 2.0) * (1.0 / MoIstar) * (G_AU_Msol_yr * massCompanion * massCompanion / R1_AU) * R1_over_a_6 * (ImK22 + ((m_Eccentricity * m_Eccentricity) *  ((ImK12 / 4.0) - (5.0 * ImK22) + (49.0 * ImK32 / 4.0))));
}

/*
 * Calculate the change in semi-major axis based on secular equations for tidal evolution given the tidal Love number
 * Zahn, 1977, Eq. (3.6)
 *
 *
 * double BaseBinaryStar::CalculateDSemiMajorAxisTidalDt(const DBL_DBL_DBL_DBL p_ImKlm, const BinaryConstituentStar* p_Star)
 *
 * @param   [IN]    p_ImKlm                     Imaginary [(1,0), (1,2), (2,2), (3,2)] components of the potential tidal Love number of star (unitless)
 * @param   [IN]    p_Star                      Star for which impact on semi-major axis is to be calculated
 * @return                                      Change in semi-major axis for binary (AU/yr)
 */    
double BaseBinaryStar::CalculateDSemiMajorAxisTidalDt(const DBL_DBL_DBL_DBL p_ImKlm, const BinaryConstituentStar* p_Star) {
    
    double massStar      = p_Star->Mass();
    double radiusStar    = p_Star->Radius();
    double massCompanion = p_Star == m_Star1 ? m_Star2->Mass() : m_Star1->Mass();
    
    double ImK10, ImK12, ImK22, ImK32;
    std::tie(ImK10, ImK12, ImK22, ImK32) = p_ImKlm;

    double R1_AU       = radiusStar * RSOL_TO_AU;
    double R1_over_a   = R1_AU / m_SemiMajorAxis;
    double R1_over_a_7 = R1_over_a * R1_over_a * R1_over_a * R1_over_a * R1_over_a * R1_over_a * R1_over_a;

    return -(3.0 / m_Omega) * (1.0 + (massCompanion / massStar)) * (G_AU_Msol_yr * massCompanion / R1_AU / R1_AU) * R1_over_a_7 * (ImK22 + ((m_Eccentricity * m_Eccentricity) * ((3.0 * ImK10 / 4.0) + (ImK12 / 8.0) - (5.0 * ImK22) + (147.0 * ImK32 / 8.0))));
}


/*
 * Resolves supernova event - one of the stars has gone supernova!
 *
 * Assign a random supernova kick according to the user specified options and then update the orbit and velocities.
 * Vector algebra is directly based on Pfahl, Rappaport, Podsiadlowski 2002, Appendix B:
 * https://arxiv.org/abs/astro-ph/0106141 
 * The change of reference basis angles, ThetaE, PhiE, and PsiE, are the standard Euler angles (see vector3d.h)
 *
 * Note: the systemic speed is only valid for intact binaries, and component speeds are only valid for disrupted binaries.
 * 
 * Logic:
 *  
 *     if (Unbound before SN):
 *  
 *         Must be 2nd SN, only need to update starSN component velocity (rotated into previous reference frame).
 *  
 *     else: (Intact before SN)
 *  
 *         Evolve binary according to vector algebra to determine centerofmass velocity, h', e', a', and whether bound or unbound.
 *         Update binary systemic velocity (even if disrupted, just for consistency) - rotate into previous reference frame if needed.
 *   
 *         if now unbound:
 *  
 *             Set m_Unbound to True - should be the only place in the code this is done.
 *  
 *             Continue vector algebra to find v1inf and v2inf.
 *             Add these values to previous component velocities (rotated if need be) which will be the systemic velocity if this is the 2nd SN. 
 *  
 *             For unbound binary, new Euler Angles should be randomized (see vector3d.cpp).
 *  
 *         if still intact:
 *  
 *             Binary systemic velocity has already been set, so just set the component velocities to the same vector.
 *             (this is to make it easier to add just a component velocity later).
 *  
 *             For intact binary, Euler Angles must be calculated according to the vector algebra (see vector3d.h).
 *
 *
 * void ResolveSupernova()
 *
 */
void BaseBinaryStar::ResolveSupernova() {
// Functions defined in vector3d.h
// Defined here for convenience - undefined later
#define cross(x,y)        Vector3d::Cross(x, y)
#define dot(x,y)          Vector3d::Dot(x, y) 
#define angleBetween(x,y) Vector3d::AngleBetween(x, y)
#define mag               Magnitude()
#define hat               UnitVector()

    // set relevant preSN parameters 
    m_EccentricityPreSN     = m_Eccentricity;                                                 
    m_SemiMajorAxisPreSN    = m_SemiMajorAxis;                                               

    double totalMassPreSN   = m_Supernova->SN_TotalMassAtCOFormation() + m_Companion->Mass();                                   // total Mass preSN
    double reducedMassPreSN = m_Supernova->SN_TotalMassAtCOFormation() * m_Companion->Mass() / totalMassPreSN;                  // reduced Mass preSN
    m_Supernova->SetOrbitalEnergyPreSN(CalculateOrbitalEnergy(reducedMassPreSN, totalMassPreSN, m_SemiMajorAxisPreSN));         // orbital energy preSN

    // define the natal kick vector (see above for precise definitions of the angles)
    double theta             = m_Supernova->SN_Theta();                                                                         // angle out of the binary plane
    double phi               = m_Supernova->SN_Phi();                                                                           // angle in the binary plane
    Vector3d natalKickVector = m_Supernova->SN_KickMagnitude() * Vector3d(cos(theta) * cos(phi), cos(theta) * sin(phi), sin(theta));
    
    // Define the rocket kick vector - will be 0 if unused. 
    // The rocket is aligned with the NS spin axis, which by default is aligned with the pre-SN orbit (0.0, 0.0, 1.0)
    // Defined here in case the system is already unbound.
    double rocketTheta        = m_Supernova->SN_RocketKickTheta();                                                              // azimuthal angle
    double rocketPhi          = m_Supernova->SN_RocketKickPhi();                                                                // polar angle
    Vector3d rocketKickVector = m_Supernova->SN_RocketKickMagnitude() * Vector3d(sin(rocketTheta) * cos(rocketPhi), sin(rocketTheta) * sin(rocketPhi), cos(rocketTheta));

    // Check if the system is already unbound
    if (IsUnbound()) {                                                                                                          // is system already unbound?
                                                                                                                                // yes
        m_Supernova->UpdateComponentVelocity( (natalKickVector+rocketKickVector).ChangeBasis(m_ThetaE, m_PhiE, m_PsiE));        // only need to update the velocity of the star undergoing SN

        m_OrbitalVelocityPreSN = 0.0;
    }
    else {                                                                                                                      // no, not unbound - evaluate orbital changes and calculate velocities
        // Evolve SN out of binary       
        
        // Pre-SN parameters
        double semiMajorAxisPrev_km     = m_SemiMajorAxis * AU_TO_KM;                                                           // semi-Major axis in km
        double eccentricityPrev         = m_Eccentricity;                                                                       // eccentricity prior to any updates to m_Eccentricity
        double sqrt1MinusEccPrevSquared = std::sqrt(1.0 - eccentricityPrev * eccentricityPrev);                                 // useful function of eccentricity

        double m1Prev                   = m_Supernova->SN_TotalMassAtCOFormation();                                             // supernova pre-SN mass (Msol)
        double m2Prev                   = m_Companion->Mass();                                                                  // companion pre-SN mass (Msol)
        double totalMassPrev            = m1Prev + m2Prev;                                                                      // total binary pre-SN mass (Msol)
        
        // Functions of eccentric anomaly
        m_Supernova->CalculateSNAnomalies(eccentricityPrev);
        double cosEccAnomaly = cos(m_Supernova->SN_EccentricAnomaly());        
        double sinEccAnomaly = sin(m_Supernova->SN_EccentricAnomaly());

        // Derived quantities
        double aPrev   = semiMajorAxisPrev_km;
        double aPrev_2 = aPrev * aPrev;
        double aPrev_3 = aPrev_2 * aPrev;

        double omega   = std::sqrt(G_km_Msol_s * totalMassPrev / aPrev_3);                                                      // Keplerian orbital frequency (rad/s)

        Vector3d separationVectorPrev = Vector3d(aPrev * (cosEccAnomaly - eccentricityPrev), aPrev * (sinEccAnomaly) * sqrt1MinusEccPrevSquared, 0.0); // relative position vector, from m1Prev to m2Prev (km)
        double separationPrev         = separationVectorPrev.mag;                                                               // instantaneous Separation (km)
        double fact1                  = aPrev_2 * omega / separationPrev;

        Vector3d relativeVelocityVectorPrev       = Vector3d(-fact1 * sinEccAnomaly, fact1 * cosEccAnomaly * sqrt1MinusEccPrevSquared, 0.0); // relative velocity vector, in the m1Prev rest frame (km/s)
        Vector3d orbitalAngularMomentumVectorPrev = cross(separationVectorPrev, relativeVelocityVectorPrev);                    // specific orbital angular momentum vector (km^2 s^-1)
        Vector3d eccentricityVectorPrev           = cross(relativeVelocityVectorPrev, orbitalAngularMomentumVectorPrev) / 
                                                    (G_km_Msol_s * totalMassPrev) - separationVectorPrev.hat;                   // Laplace-Runge-Lenz vector (magnitude = eccentricity)

        m_OrbitalVelocityPreSN = relativeVelocityVectorPrev.mag;                                                                // pre-SN orbital velocity (km/s) 

        // Note: In the following,
        // orbitalAngularMomentumVectorPrev defines the Z-axis, 
        // eccentricityVectorPrev defines the X-axis, and
        // (orbitalAngularMomentumVectorPrev x eccentricityVectorPrev) defines the Y-axis
        
        // Apply supernova natal kick and mass loss  
        //
        // Note: the code allows for mass loss and kick in the companion 
        // (due to ablation), though we currently do not apply these.
        
        Vector3d companionRecoilVector = Vector3d(0.0, 0.0, 0.0);                                                               // km/s - The recoil of the companion due to ablation

        double m1        = m_Supernova->Mass();                                                                                 // supernova post-SN mass (Msol)
        double m2        = m_Companion->Mass();                                                                                 // companion post-SN mass (Msol)
        double totalMass = m1 + m2;                                                                                             // total binary post-SN mass (Msol)
        double fact2     = totalMassPrev * totalMass;       
        double dm1       = (m1Prev - m1);                                                                                       // mass difference of supernova (Msol)
        double dm2       = (m2Prev - m2);                                                                                       // mass difference of companion (Msol)

        Vector3d centerOfMassVelocity         = (-m2Prev * dm1 / fact2 + m1Prev * dm2 / fact2) * relativeVelocityVectorPrev + 
                                                (m1 / totalMass) * natalKickVector + (m2 / totalMass) * companionRecoilVector;  // post-SN center of mass velocity vector (km/s)

        Vector3d relativeVelocityVector       = relativeVelocityVectorPrev + (natalKickVector - companionRecoilVector);         // post-SN relative velocity vector (km/s)

        Vector3d orbitalAngularMomentumVector = cross(separationVectorPrev, relativeVelocityVector);                            // post-SN specific orbital angular momentum vector (km^2 s^-1)
        double   orbitalAngularMomentum = orbitalAngularMomentumVector.mag;                                                     // post-SN specific orbital angular momentum (km^2 s^-1)
        m_NormalizedOrbitalAngularMomentumVector = orbitalAngularMomentumVector/orbitalAngularMomentum;                         // set unit vector here to make printing out the inclination vector easier

        Vector3d eccentricityVector           = cross(relativeVelocityVector, orbitalAngularMomentumVector) / 
                                                (G_km_Msol_s * totalMass) - separationVectorPrev / separationPrev;              // post-SN Laplace-Runge-Lenz vector
        m_Eccentricity                        = eccentricityVector.mag;                                                         // post-SN eccentricity
        double eccSquared                     = m_Eccentricity * m_Eccentricity;                                                // useful function of eccentricity

        double semiMajorAxis_km               = (orbitalAngularMomentum * orbitalAngularMomentum) / (G_km_Msol_s * totalMass * (1.0 - eccSquared));                                 // post-SN semi-major axis (km)
        m_SemiMajorAxis                       = semiMajorAxis_km * KM_TO_AU;                                                    // post-SN semi-major axis (AU)

        // Note: similar to above,
        // orbitalAngularMomentumVector defines the Z'-axis, 
        // eccentricityVector defines the X'-axis, and
        // (orbitalAngularMomentumVector x eccentricityVector) defines the Y'-axis
         
        UpdateSystemicVelocity(centerOfMassVelocity.ChangeBasis(m_ThetaE, m_PhiE, m_PsiE));                                     // update the system velocity with the new center of mass velocity
        double reducedMass = m_Supernova->Mass() * m_Companion->Mass() / totalMass;                                             // reduced Mass
        m_Supernova->SetOrbitalEnergyPostSN(CalculateOrbitalEnergy(reducedMass, totalMass, m_SemiMajorAxis));                   // orbital energy

        // Split off and evaluate depending on whether the binary is now bound or unbound
	    if (utils::Compare(m_Eccentricity, 1.0) >= 0) {                                                                         // unbound?
                                                                                                                                // yes, unbound            
            m_Unbound = true;

            // Calculate the asymptotic Center of Mass velocity 
            double   relativeVelocityAtInfinity       = (G_km_Msol_s*totalMass/orbitalAngularMomentum) * std::sqrt(eccSquared - 1.0);
            Vector3d relativeVelocityVectorAtInfinity = relativeVelocityAtInfinity 
                                                        * (-1.0 * (eccentricityVector.hat / m_Eccentricity) 
                                                        + std::sqrt(1.0 - 1.0 / eccSquared) * cross(orbitalAngularMomentumVector.hat, eccentricityVector.hat));

            // Calculate the asymptotic velocities of Star1 (SN) and Star2 (CP)
            Vector3d component1VelocityVectorAtInfinity =  (m2 / totalMass) * relativeVelocityVectorAtInfinity + centerOfMassVelocity;
            Vector3d component2VelocityVectorAtInfinity = -(m1 / totalMass) * relativeVelocityVectorAtInfinity + centerOfMassVelocity;

            // Update the component velocities 
            m_Supernova->UpdateComponentVelocity(component1VelocityVectorAtInfinity.ChangeBasis(m_ThetaE, m_PhiE, m_PsiE));
            m_Companion->UpdateComponentVelocity(component2VelocityVectorAtInfinity.ChangeBasis(m_ThetaE, m_PhiE, m_PsiE));

            // Set Euler Angles 
            m_ThetaE = angleBetween(orbitalAngularMomentumVectorPrev, orbitalAngularMomentumVector);                            // angle between the angular momentum unit vectors, always well defined
            m_PhiE   = _2_PI * RAND->Random(); 
            m_PsiE   = _2_PI * RAND->Random(); 
        }
        else {                                                                                                                  // no - binary still bound

            // Set the component velocites to the system velocity. System velocity was already correctly set above.
             
            m_Supernova->UpdateComponentVelocity(centerOfMassVelocity.ChangeBasis(m_ThetaE, m_PhiE, m_PsiE));
            m_Companion->UpdateComponentVelocity(centerOfMassVelocity.ChangeBasis(m_ThetaE, m_PhiE, m_PsiE));

            // Calculate Euler angles - see ChangeBasis() in vector.cpp for details
            m_ThetaE = angleBetween(orbitalAngularMomentumVector, orbitalAngularMomentumVectorPrev); // angle between the angular momentum unit vectors, always well defined

            // If the new orbital A.M. is parallel or anti-parallel to the previous orbital A.M., 
            // then the cross product is not well-defined, and we need to account for degeneracy between eccentricity vectors.
            // Also, if either eccentricity is 0.0, then the eccentricity vector is not well defined.

            if ((utils::Compare(m_ThetaE, 0.0) == 0) &&                                                                         // orbitalAngularMomentumVectorPrev parallel to orbitalAngularMomentumVector?
                ((utils::Compare(eccentricityPrev, 0.0) > 0) && (utils::Compare(m_Eccentricity, 0.0) > 0))) {                   // ... and both eccentricityVectorPrev and eccentricityVector well-defined?

                double psiPlusPhi = angleBetween(eccentricityVector, eccentricityVectorPrev);                                   // yes - then psi + phi is constant
                m_PhiE            = _2_PI * RAND->Random();    
                m_PsiE            = psiPlusPhi - m_PhiE;
            }
            else if ((utils::Compare(m_ThetaE, M_PI) == 0) &&                                                                   // orbitalAngularMomentumVectorPrev anti-parallel to orbitalAngularMomentumVector?
                    ((utils::Compare(eccentricityPrev, 0.0) > 0) &&  (utils::Compare(m_Eccentricity, 0.0) > 0))) {              // ... and both eccentricityVectorPrev and eccentricityVector well-defined?
                                                                                                   
                double psiMinusPhi = angleBetween(eccentricityVector, eccentricityVectorPrev);                                  // yes - then psi - phi is constant
                m_PhiE             = _2_PI * RAND->Random();    
                m_PsiE             = psiMinusPhi + m_PhiE;
            }
            else {                                                                                                              // neither - the cross product of the orbit normals is well-defined
                Vector3d orbitalPivotAxis = cross(orbitalAngularMomentumVectorPrev, orbitalAngularMomentumVector);              // cross product of the orbit normals

                if (utils::Compare(eccentricityPrev, 0.0) == 0 ) {                                                              // eccentricityVectorPrev well-defined?
                    m_PhiE = _2_PI * RAND->Random();                                                                            // no - set phi random
                }
                else {                                                                                                          // yes - phi is +/- angle between eccentricityVectorPrev and orbitalPivotAxis
                    m_PhiE = utils::Compare( dot(eccentricityVectorPrev, orbitalAngularMomentumVector), 0.0) >= 0               // are eccentricityVectorPrev and orbitalAngularMomentumVector in the same hemisphere?
                        ? angleBetween(eccentricityVectorPrev, orbitalPivotAxis)                                                // yes - phi in [0,pi)
                        : -angleBetween(eccentricityVectorPrev, orbitalPivotAxis);                                              // no  - phi in [-pi,0)
                }

                if ( utils::Compare(m_Eccentricity, 0.0) == 0 ) {                                                               // is eccentricityVector well-defined?
                    m_PsiE = _2_PI * RAND->Random();                                                                            // no - set psi random 
                }                                                                                              
                else {                                                                                                          // yes - psi is +/- angle between eccentricityVector and orbitalPivotAxis
                    m_PsiE = utils::Compare( dot(eccentricityVector, orbitalAngularMomentumVectorPrev), 0.0) >= 0               // are eccentricityVector and orbitalAngularMomentumVectorPrev in the same hemisphere?
                    ? angleBetween(eccentricityVector, orbitalPivotAxis)                                                        // yes - psi in [0,pi)
                    : -angleBetween(eccentricityVector, orbitalPivotAxis);                                                      // no  - psi in [-pi,0)
                }
            }

            // Note: There is some evidence for evolution of periapsis in mass transferring binaries (see e.g Dosopoulou & Kalogera 2016, 2018). 
            // This should be investigated in more depth, but until then, we assume that the periapsis *may* evolve, and accordingly randomize
            // the angle of periapsis around the new orbital angular momentum, (i.e, Psi) - RTW 15/05/20
            m_PsiE = _2_PI * RAND->Random();
        }
        
        // account for possible neutrino rocket - see Hirai+ 2024
        if (ShouldResolveNeutrinoRocketMechanism()) {

            if (IsUnbound()) {                                                                                                  // is system unbound? 
                m_Supernova->UpdateComponentVelocity(rocketKickVector.ChangeBasis(m_ThetaE, m_PhiE, m_PsiE));                   // yes - simply update the component velocity
            }
            else {                                                                                                              // no - need to update the eccentricity and system velocity
                Vector3d eccentricityVectorPreRocket             = eccentricityVector;                                          // defined earlier
                double averageOrbitalVelocityPreRocket           = std::sqrt(-2.0 * m_OrbitalEnergy/reducedMass);               // average orbital velocity post-SN (AU/yr)
                double kGrav                                     = averageOrbitalVelocityPreRocket * averageOrbitalVelocityPreRocket * reducedMass * m_SemiMajorAxis; // AU^3 * Msol / yr^2
                Vector3d totalAmVectorPreRocket                  = orbitalAngularMomentumVector * reducedMass * KM_TO_AU * KM_TO_AU * SECONDS_IN_YEAR; // Msol * AU^2 / yr (orbitalAngularMomentumVector is the specific orbital AM)
                Vector3d amVectorNormalizedByCircularAmPreRocket = totalAmVectorPreRocket * (averageOrbitalVelocityPreRocket / kGrav); // unitless!
                double theta_rotation                            = 3.0 * rocketKickVector.mag * KM_TO_AU * SECONDS_IN_YEAR / (2.0 * averageOrbitalVelocityPreRocket); // rad - need to convert velocities to same units
                    
                // apply hPlus and hMinus support vectors
                Vector3d hPlusVector  = amVectorNormalizedByCircularAmPreRocket + eccentricityVectorPreRocket;
                Vector3d hMinusVector = amVectorNormalizedByCircularAmPreRocket - eccentricityVectorPreRocket;

                // rotate hPlus and hMinus vectors so that the thrust is parallel to the z-axis, in order to apply the rotation below
                hPlusVector  = hPlusVector.RotateVectorAboutZ( -rocketPhi).RotateVectorAboutY(-rocketTheta);
                hMinusVector = hMinusVector.RotateVectorAboutZ(-rocketPhi).RotateVectorAboutY(-rocketTheta);

                // rotate vectors about the new "z-axis" - parallel to the rocket thrust
                Vector3d hPlusVector_prime  = hPlusVector.RotateVectorAboutZ(  theta_rotation);
                Vector3d hMinusVector_prime = hMinusVector.RotateVectorAboutZ(-theta_rotation);

                // rotate new hPlus and hMinus vectors back to the original frame
                hPlusVector  = hPlusVector.RotateVectorAboutY( rocketTheta).RotateVectorAboutZ(rocketPhi);
                hMinusVector = hMinusVector.RotateVectorAboutY(rocketTheta).RotateVectorAboutZ(rocketPhi);

                // calculate post-rocket values
                Vector3d normalizedAngularMomentumVectorPostRocket = 0.5 * (hPlusVector_prime + hMinusVector_prime);
                Vector3d eccentricityVectorPostRocket              = 0.5 * (hPlusVector_prime - hMinusVector_prime);

                m_NormalizedOrbitalAngularMomentumVector = normalizedAngularMomentumVectorPostRocket ;                 
                m_Eccentricity                           = eccentricityVectorPostRocket.mag;                                                        

                UpdateSystemicVelocity(rocketKickVector.ChangeBasis(m_ThetaE, m_PhiE, m_PsiE));                            
                m_Supernova->UpdateComponentVelocity(rocketKickVector.ChangeBasis(m_ThetaE, m_PhiE, m_PsiE));
                m_Companion->UpdateComponentVelocity(rocketKickVector.ChangeBasis(m_ThetaE, m_PhiE, m_PsiE));
            }
        }

        #undef hat
        #undef mag        
        #undef angleBetween
        #undef dot
        #undef cross
    }

    // Do for all systems 

    m_IPrime    = m_ThetaE;                                                                                                     // inclination angle between preSN and postSN orbital planes 
    m_CosIPrime = cos(m_IPrime);

    (void)PrintSupernovaDetails();                                                                                              // log record to supernovae logfile
    m_Supernova->ClearCurrentSNEvent();

#undef hat
#undef mag        
#undef angleBetween
#undef dot
#undef cross
}


/*
 * Determine if one or both of the stars are undergoing a supernova event,
 * and if so resolve the event(s) by calling ResolveSupernova() for each of
 * the stars as appropriate.
 *
 * void EvaluateSupernovae
 * 
 */
void BaseBinaryStar::EvaluateSupernovae() {

    m_SupernovaState = SN_STATE::NONE;                                  // not yet determined
    
    if (m_Star1->IsSNevent()) {                                         // star1 supernova
        m_SupernovaState = SN_STATE::STAR1;                             // star1

        // resolve star1 supernova
        m_Supernova = m_Star1;                                          // supernova
        m_Companion = m_Star2;                                          // companion
        ResolveSupernova();                                             // resolve supernova
    }

    if (m_Star2->IsSNevent()) {                                         // star2 supernova                                                                                                        
        m_SupernovaState = m_SupernovaState == SN_STATE::NONE           // star1 not supernova?
                            ? SN_STATE::STAR2                           // yes - just star2
                            : SN_STATE::BOTH;                           // no - both 

        m_Supernova = m_Star2;                                          // supernova
        m_Companion = m_Star1;                                          // companion
        ResolveSupernova();                                             // resolve supernova
    }
}


/*
 * Resolve a Common Envelope Event
 *
 * The binary has entered a common envelope event. This function updates the binary parameters accordingly
 *
 * From Hurley et al. 2002, section 2.7.1:
 *
 *    Common-envelope evolution occurs either as a result of a collision between
 *    a star with a dense core (k1 {2,3,4,5,6,8,9}) or at the onset of RLOF where mass
 *    is transferred from a giant (k1 {2,3,4,5,6,8,9}) on a dynamical time-scale
 *
 *
 * void ResolveCommonEnvelopeEvent()
 */
void BaseBinaryStar::ResolveCommonEnvelopeEvent() {

    double alphaCE = OPTIONS->CommonEnvelopeAlpha();                                                                    // CE efficiency parameter

	double eccentricity      = Eccentricity();								                                            // current eccentricity (before CEE)
    double semiMajorAxisRsol = SemiMajorAxisRsol();                                                                     // current semi-major axis in default units, Rsol (before CEE)
    double periastronRsol    = PeriastronRsol();                                                                        // periastron, Rsol (before CEE)
    double rRLd1Rsol         = periastronRsol * CalculateRocheLobeRadius_Static(m_Star1->Mass(), m_Star2->Mass());      // Roche-lobe radius at periastron in Rsol at the moment where CEE begins, seen by star1
    double rRLd2Rsol         = periastronRsol * CalculateRocheLobeRadius_Static(m_Star2->Mass(), m_Star1->Mass());      // Roche-lobe radius at periastron in Rsol at the moment where CEE begins, seen by star2
    
    bool isDonorMS = false;                                                                                             // check for main sequence donor
    if (OPTIONS->AllowMainSequenceStarToSurviveCommonEnvelope()) {                                                      // allow main sequence stars to survive CEE?
        if (m_Star1->IsOneOf(ALL_MAIN_SEQUENCE)) {                                                                      // yes - star1 MS_LTE_07, MS_GT_07 or NAKED_HELIUM_STAR_MS?
            isDonorMS    = isDonorMS || m_Star1->IsRLOF();                                                              // yes - donor MS?
            m_Mass1Final = m_Star1->Mass();                                                                             // set mass
            m_MassEnv1   = 0.0;                                                                                         // no envelope
        }
        else {                                                                                                          // no, star1 not MS_LTE_07, MS_GT_07 or NAKED_HELIUM_STAR_MS
            m_Mass1Final = m_Star1->CoreMass();                                                                         // set mass
            m_MassEnv1   = m_Star1->Mass() - m_Star1->CoreMass();                                                       // and envelope
        }

        if (m_Star2->IsOneOf(ALL_MAIN_SEQUENCE)) {                                                                      // star2 MS_LTE_07, MS_GT_07 or NAKED_HELIUM_STAR_MS?
            isDonorMS    = isDonorMS || m_Star2->IsRLOF();                                                              // yes - donor MS?
            m_Mass2Final = m_Star2->Mass();                                                                             // yes - set mass
            m_MassEnv2   = 0.0;                                                                                         // no envelope
        }
        else {                                                                                                          // no, star2 not MS_LTE_07, MS_GT_07 or NAKED_HELIUM_STAR_MS
            m_Mass2Final = m_Star2->CoreMass();                                                                         // set mass
            m_MassEnv2   = m_Star2->Mass() - m_Star2->CoreMass();                                                       // and envelope
        }
    }
    else {                                                                                                              // no don't allow main sequence stars to survive CEE; should lead to stellar merger
        m_Mass1Final = m_Star1->CoreMass();                                                                             // set mass1
        m_MassEnv1   = m_Star1->Mass() - m_Star1->CoreMass();                                                           // and envelope1
        m_Mass2Final = m_Star2->CoreMass();                                                                             // set mass2
        m_MassEnv2   = m_Star2->Mass() - m_Star2->CoreMass();                                                           // and envelope2
    }

    bool envelopeFlag1 = utils::Compare(m_MassEnv1, 0.0) > 0 && utils::Compare(m_Mass1Final, 0.0) > 0;                  // star1 not massless remnant and has envelope?
    bool envelopeFlag2 = utils::Compare(m_MassEnv2, 0.0) > 0 && utils::Compare(m_Mass2Final, 0.0) > 0;                  // star2 not massless remnant and has envelope?
    m_CEDetails.doubleCoreCE = envelopeFlag1 && envelopeFlag2;

    m_CEDetails.CEEcount++;                                                                                             // increment CEE count
    m_RLOFDetails.simultaneousRLOF = m_Star1->IsRLOF() && m_Star2->IsRLOF();                                            // check for simultaneous RLOF

	m_Star1->CalculateLambdas(m_MassEnv1);                                                                              // calculate lambdas for star1
	m_Star2->CalculateLambdas(m_MassEnv2);                                                                              // calculate lambdas for star2

    m_Star1->CalculateBindingEnergies(m_Mass1Final, m_MassEnv1, m_Star1->Radius());                                     // calculate binding energies for star1 (uses lambdas)
    m_Star2->CalculateBindingEnergies(m_Mass2Final, m_MassEnv2, m_Star2->Radius());                                     // calculate binding energies for star2 (uses lambdas)

    m_Star1->CalculateCommonEnvelopeValues();                                                                           // calculate common envelope values for star1
    m_Star2->CalculateCommonEnvelopeValues();                                                                           // calculate common envelope values for star2

    double lambda1 = m_Star1->LambdaAtCEE();                                                                            // measures the envelope binding energy of star 1
    double lambda2 = m_Star2->LambdaAtCEE();                                                                            // measures the envelope binding energy of star 2

    m_Star1->SetPreCEEValues();                                                                                         // squirrel away pre CEE stellar values for star 1
    m_Star2->SetPreCEEValues();                                                                                         // squirrel away pre CEE stellar values for star 2
  	SetPreCEEValues(semiMajorAxisRsol, eccentricity, rRLd1Rsol, rRLd2Rsol);                                             // squirrel away pre CEE binary values
    
	// double common envelope phase prescription (Brown 1995) to calculate new semi-major axis
	// due to the CEE as described in Belczynsky et al. 2002, eq. (12)
    
    switch (OPTIONS->CommonEnvelopeFormalism()) {
        case CE_FORMALISM::ENERGY: {

            double k1         = m_Star1->IsOneOf(COMPACT_OBJECTS) ? 0.0 : (2.0 / (lambda1 * alphaCE)) * m_Star1->Mass() * m_MassEnv1 / m_Star1->Radius();
            double k2         = m_Star2->IsOneOf(COMPACT_OBJECTS) ? 0.0 : (2.0 / (lambda2 * alphaCE)) * m_Star2->Mass() * m_MassEnv2 / m_Star2->Radius();
            double k3         = m_Star1->Mass() * m_Star2->Mass() / periastronRsol;                                     // assumes immediate circularisation at periastron at start of CE
            double k4         = (m_Mass1Final * m_Mass2Final);
            double aFinalRsol = k4 / (k1 + k2 + k3);
            m_SemiMajorAxis   = aFinalRsol * RSOL_TO_AU;
            } break;

        case CE_FORMALISM::TWO_STAGE: {
            // two-stage common envelope, Hirai & Mandel (2022)

            double convectiveEnvelopeMass1, maxConvectiveEnvelopeMass1;
            std::tie(convectiveEnvelopeMass1, maxConvectiveEnvelopeMass1) = m_Star1->CalculateConvectiveEnvelopeMass();

            double radiativeIntershellMass1 = m_MassEnv1 - convectiveEnvelopeMass1;
            double endOfFirstStageMass1     = m_Mass1Final + radiativeIntershellMass1;

            double convectiveEnvelopeMass2, maxConvectiveEnvelopeMass2;
            std::tie(convectiveEnvelopeMass2, maxConvectiveEnvelopeMass2) = m_Star2->CalculateConvectiveEnvelopeMass();

            double radiativeIntershellMass2 = m_MassEnv2 - convectiveEnvelopeMass2;
            double endOfFirstStageMass2     = m_Mass2Final + radiativeIntershellMass2;
        
            // stage 1: convective envelope removal on a dynamical timescale; assumes lambda = lambda_He
            double lambda1    = m_Star1->CalculateConvectiveEnvelopeLambdaPicker(convectiveEnvelopeMass1, maxConvectiveEnvelopeMass1);
            double lambda2    = m_Star1->CalculateConvectiveEnvelopeLambdaPicker(convectiveEnvelopeMass2, maxConvectiveEnvelopeMass2);
        
            double k1         = m_Star1->IsOneOf(COMPACT_OBJECTS) ? 0.0 : (2.0 / (lambda1 * alphaCE)) * m_Star1->Mass() * convectiveEnvelopeMass1 / m_Star1->Radius();
            double k2         = m_Star2->IsOneOf(COMPACT_OBJECTS) ? 0.0 : (2.0 / (lambda2 * alphaCE)) * m_Star2->Mass() * convectiveEnvelopeMass2 / m_Star2->Radius();
            double k3         = m_Star1->Mass() * m_Star2->Mass() / periastronRsol;                                     // assumes immediate circularisation at periastron at start of CE
            double k4         = endOfFirstStageMass1 * endOfFirstStageMass2;

            double aFinalRsol = k4 / (k1 + k2 + k3);
            m_SemiMajorAxis   = aFinalRsol * RSOL_TO_AU;

            // stage 2: radiative envelope removal on a thermal timescale; assumed to be fully non-conservative
            if (utils::Compare(radiativeIntershellMass1, 0.0) > 0) {
                m_SemiMajorAxis = CalculateMassTransferOrbit(endOfFirstStageMass1, -radiativeIntershellMass1, *m_Star2, 0.0);
            }

            if (utils::Compare(radiativeIntershellMass2, 0.0) > 0) {
                m_SemiMajorAxis = CalculateMassTransferOrbit(endOfFirstStageMass2, -radiativeIntershellMass2, *m_Star1, 0.0);
            }
        } break;

        default:                                                                                                        // unknown prescription
            // the only way this can happen is if someone added a CE_FORMALISM
            // and it isn't accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a prescription this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the binary.
            // The correct fix for this is to add code for the missing prescription or, if the missing
            // prescription is superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNKNOWN_CE_FORMALISM);                                                                   // throw error
    }
    
	double rRLdfin1     = m_SemiMajorAxis * CalculateRocheLobeRadius_Static(m_Mass1Final, m_Mass2Final);                // Roche-lobe radius in AU after CEE, seen by star1
	double rRLdfin2     = m_SemiMajorAxis * CalculateRocheLobeRadius_Static(m_Mass2Final, m_Mass1Final);                // Roche-lobe radius in AU after CEE, seen by star2
    double rRLdfin1Rsol = rRLdfin1 * AU_TO_RSOL;                                                                        // Roche-lobe radius in Rsol after CEE, seen by star1
    double rRLdfin2Rsol = rRLdfin2 * AU_TO_RSOL;                                                                        // Roche-lobe radius in Rsol after CEE, seen by star2
    m_Eccentricity      = 0.0;                                                                                          // we assume that a common envelope event (CEE) circularises the binary

    m_Star1->ResolveCommonEnvelopeAccretion(m_Mass1Final);                                                              // update star1's mass after CE accretion
    m_Star2->ResolveCommonEnvelopeAccretion(m_Mass2Final);                                                              // update star2's mass after CE accretion

    // update stellar type after losing its envelope. Star1, Star2 or both if double CEE.

    if (isDonorMS || (!envelopeFlag1 && !envelopeFlag2)) {                                                              // stellar merger
        m_MassTransferTrackerHistory = MT_TRACKING::MERGER; 
        m_Flags.stellarMerger        = true;
    }
    else if ( (m_Star1->DetermineEnvelopeType()==ENVELOPE::RADIATIVE && !m_Star1->IsOneOf(ALL_MAIN_SEQUENCE)) ||
              (m_Star2->DetermineEnvelopeType()==ENVELOPE::RADIATIVE && !m_Star2->IsOneOf(ALL_MAIN_SEQUENCE)) ) {       // check if we have a non-MS radiative-envelope star
        if (!OPTIONS->AllowRadiativeEnvelopeStarToSurviveCommonEnvelope() && OPTIONS->CommonEnvelopeFormalism()!=CE_FORMALISM::TWO_STAGE) {                                            // stellar merger
            m_CEDetails.optimisticCE = true;
            m_MassTransferTrackerHistory = MT_TRACKING::MERGER;
            m_Flags.stellarMerger        = true;
        }
    }

	if (!m_Flags.stellarMerger) {

        STELLAR_TYPE stellarType1 = m_Star1->StellarType();                                                             // star 1 stellar type before resolving envelope loss
        STELLAR_TYPE stellarType2 = m_Star2->StellarType();                                                             // star 2 stellar type before resolving envelope loss
        
        if (envelopeFlag1) {
            m_Star1->ResolveEnvelopeLossAndSwitch();                                                                    // resolve envelope loss for star1 and switch to new stellar type
            m_MassTransferTrackerHistory = MT_TRACKING::CE_1_TO_2_SURV;
        }
        if (envelopeFlag2) {
            m_Star2->ResolveEnvelopeLossAndSwitch();                                                                    // resolve envelope loss for star1 and switch to new stellar type
            m_MassTransferTrackerHistory = MT_TRACKING::CE_2_TO_1_SURV;
        }
        if (m_CEDetails.doubleCoreCE)
            m_MassTransferTrackerHistory = MT_TRACKING::CE_DOUBLE_SURV;                                                 // record history - double CEE

        m_Star1->UpdateAttributes(0.0, 0.0, true);
        m_Star2->UpdateAttributes(0.0, 0.0, true);

        if (m_Star1->StellarType() != stellarType1 || m_Star2->StellarType() != stellarType2) {                         // stellar type change?
            (void)PrintDetailedOutput(m_Id, BSE_DETAILED_RECORD_TYPE::STELLAR_TYPE_CHANGE_DURING_CEE);                  // yes - print (log) detailed output
        }
	}

    if (utils::Compare(m_SemiMajorAxis, 0.0) <= 0 || utils::Compare(m_Star1->Radius() + m_Star2->Radius(), m_SemiMajorAxis * AU_TO_RSOL) > 0) {
        m_Flags.stellarMerger = true;
    }

    // if CHE enabled, update rotational frequency for constituent stars - assume tidally locked
    double omega = OrbitalAngularVelocity();                                                                            // orbital angular velocity
    if (OPTIONS->CHEMode() != CHE_MODE::NONE) m_Star1->SetOmega(omega);
    if (OPTIONS->CHEMode() != CHE_MODE::NONE) m_Star2->SetOmega(omega);

    // if both stars evolving as chemically homogeneous stars set m_Omega for binary
    if (HasTwoOf({STELLAR_TYPE::CHEMICALLY_HOMOGENEOUS})) m_Omega = omega;
    
    m_Star1->SetPostCEEValues();                                                                                        // squirrel away post CEE stellar values for star 1
    m_Star2->SetPostCEEValues();                                                                                        // squirrel away post CEE stellar values for star 2
    SetPostCEEValues(m_SemiMajorAxis * AU_TO_RSOL, m_Eccentricity, rRLdfin1Rsol, rRLdfin2Rsol);                         // squirrel away post CEE binary values (checks for post-CE RLOF, so should be done at end)

    if (m_RLOFDetails.immediateRLOFPostCEE == true && !OPTIONS->AllowImmediateRLOFpostCEToSurviveCommonEnvelope()) {    // is there immediate post-CE RLOF which is not allowed?
            m_MassTransferTrackerHistory = MT_TRACKING::MERGER;
            m_Flags.stellarMerger        = true;
    }

    (void)PrintCommonEnvelope();                                                                                        // print (log) common envelope details
}

/*
 * Resolve a main-sequence merger event
 *
 * Star1 will become the merger product; Star2 will become a massless remnant
 *
 * void ResolveMainSequenceMerger()
 *
 */
void BaseBinaryStar::ResolveMainSequenceMerger() {
    if (!(m_Star1->IsOneOf(MAIN_SEQUENCE) && m_Star2->IsOneOf(MAIN_SEQUENCE) && OPTIONS->EvolveMainSequenceMergers()))
        return;                                                                                 // nothing to do if does not satisfy conditions for MS merger
	
    double mass1 = m_Star1->Mass();
    double mass2 = m_Star2->Mass();
    double tau1  = m_Star1->Tau();
    double tau2  = m_Star2->Tau();

    // /*ILYA*/ temporary solution, should use TAMS core mass
    double TAMSCoreMass1 = 0.3 * mass1;
    double TAMSCoreMass2 = 0.3 * mass2;
    
    double q   = std::min(mass1 / mass2, mass2 / mass1);
    double phi = 0.3 * q / (1.0 + q) / (1.0 + q);                                               // fraction of mass lost in merger, Wang+ 2022, https://www.nature.com/articles/s41550-021-01597-5
	
    double finalMass               = (1.0 - phi) * (mass1 + mass2);
    double initialHydrogenFraction = 1.0 - utils::MESAZAMSHeliumFractionByMetallicity(m_Star1->Metallicity()) - m_Star1->Metallicity();
    double finalHydrogenMass       = finalMass * initialHydrogenFraction - tau1 * TAMSCoreMass1 * initialHydrogenFraction - tau2 * TAMSCoreMass2 * initialHydrogenFraction;
    
    m_SemiMajorAxis = std::numeric_limits<float>::infinity();                                   // set separation to infinity to avoid subsequent fake interactions with a massless companion (RLOF, CE, etc.)
    
    m_Star1->UpdateAfterMerger(finalMass, finalHydrogenMass);
    
    m_Star2->SwitchTo(STELLAR_TYPE::MASSLESS_REMNANT);
}


/*
 * Calculate the Roche Lobe radius given the input masses
 *
 * Eggleton 1983, eq 2
 *
 *
 * double CalculateRocheLobeRadius_Static(const double p_MassPrimary, const double p_MassSecondary)
 *
 * @param   [IN]    p_MassPrimary               Mass, in Msol, of the primary star
 * @param   [IN]    p_MassSecondary             Mass, in Msol, of the secondary star
 * @return                                      Radius of Roche Lobe in units of the semi-major axis a
 */
double BaseBinaryStar::CalculateRocheLobeRadius_Static(const double p_MassPrimary, const double p_MassSecondary) {
    double q         = p_MassPrimary / p_MassSecondary;
    double qCubeRoot = std::cbrt(q);                                                                                    // cube roots are expensive, only compute once
    return 0.49 / (0.6 + log(1.0 + qCubeRoot) / qCubeRoot / qCubeRoot);
}


/*
 * Calculate the fraction of specific angular momentum with which the non-accreted mass leaves the system
 *
 * This is gamma (as in Pols's notes) or jloss (as in Belczynski et al. 2008
 * which is the fraction of specific angular momentum with which the non-accreted mass leaves the system.
 * Macleod_linear comes from Willcox et al. (2022)
 *
 * Calculation is based on user-specified Angular Momentum Loss prescription
 *
 *
 * double CalculateGammaAngularMomentumLoss_Static(const double p_DonorMass, const double p_AccretorMass, const bool p_IsAccretorDegenerate)
 *
 * @param   [IN]    p_DonorMass                 The mass of the donor (Msol)
 * @param   [IN]    p_AccretorMass              The mass of the accretor (Msol)
 * @param   [IN]    p_IsAccretorDegenerate      True if the accretor is a degenerate star, false otherwise (need to know up front to keep this function static)
 * @return                                      The fraction of specific angular momentum with which the non-accreted mass leaves the system
 */
double BaseBinaryStar::CalculateGammaAngularMomentumLoss_Static(const double p_DonorMass, const double p_AccretorMass, const bool p_IsAccretorDegenerate) {

	double gamma;

	switch (OPTIONS->MassTransferAngularMomentumLossPrescription()) {                                                               // which prescription?

        case MT_ANGULAR_MOMENTUM_LOSS_PRESCRIPTION::JEANS                : gamma = p_AccretorMass / p_DonorMass; break;             // vicinity of the donor 

        case MT_ANGULAR_MOMENTUM_LOSS_PRESCRIPTION::ISOTROPIC_RE_EMISSION: gamma = p_DonorMass / p_AccretorMass; break;             // vicinity of the accretor
        
        case MT_ANGULAR_MOMENTUM_LOSS_PRESCRIPTION::ARBITRARY            : gamma = OPTIONS->MassTransferJloss(); break;

        case MT_ANGULAR_MOMENTUM_LOSS_PRESCRIPTION::CIRCUMBINARY_RING: {                                                            // based on the assumption that a_ring = 2*a, Vinciguerra+, 2020 
            double sumMasses = p_DonorMass + p_AccretorMass;
            gamma            = (M_SQRT2 * sumMasses * sumMasses) / (p_DonorMass * p_AccretorMass);
            } break;

        case MT_ANGULAR_MOMENTUM_LOSS_PRESCRIPTION::MACLEOD_LINEAR : {                                                              // linear interpolation on separation between accretor and L2 point
            // interpolate in separation between a_acc and a_L2, both normalized to units of separation a
            double q        = p_AccretorMass / p_DonorMass;
            double qPlus1   = 1.0 + q;
            double aL2      = std::sqrt(M_SQRT2);                                                                                   // roughly, coincides with CIRCUMBINARY_RING def above
            double aAcc     = 1.0 / qPlus1;
            double fMacleod = p_IsAccretorDegenerate 
                                ? OPTIONS->MassTransferJlossMacLeodLinearFractionDegen() 
                                : OPTIONS->MassTransferJlossMacLeodLinearFractionNonDegen();
            double aGamma   = aAcc + (aL2 - aAcc) * fMacleod;
            gamma           = aGamma * aGamma * qPlus1 * qPlus1 / q;
            } break;

        default:                                                                                                                    // unknown prescription
            // the only way this can happen is if someone added an MT_ANGULAR_MOMENTUM_LOSS_PRESCRIPTION
            // and it isn't accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a prescription this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the binary.
            // The correct fix for this is to add code for the missing prescription or, if the missing
            // prescription is superfluous, remove it from the option.

            THROW_ERROR_STATIC(ERROR::UNKNOWN_MT_ANGULAR_MOMENTUM_LOSS_PRESCRIPTION);                                               // throw error
    }
    return gamma;
}


/*
 * Calculate new semi-major axis due to angular momentum loss
 *
 * Pols et al. notes; Belczynski et al. 2008, eq 32, 33
 *
 *
 * double CalculateMassTransferOrbit (const double                 p_DonorMass, 
 *                                    const double                 p_DeltaMassDonor, 
 *                                          BinaryConstituentStar& p_Accretor,
 *                                    const double                 p_FractionAccreted)
 *
 * @param   [IN]    p_DonorMass                 Donor mass
 * @param   [IN]    p_DeltaMassDonor            Change in donor mass
 * @param   [IN]    p_Accretor                  Pointer to accretor
 * @param   [IN]    p_FractionAccreted          Mass fraction lost from donor accreted by accretor
 * @return                                      Semi-major axis
 */
double BaseBinaryStar::CalculateMassTransferOrbit(const double                 p_DonorMass, 
                                                  const double                 p_DeltaMassDonor, 
                                                        BinaryConstituentStar& p_Accretor,
                                                  const double                 p_FractionAccreted) {

    double semiMajorAxis = m_SemiMajorAxis;
    
    if (utils::Compare(p_DeltaMassDonor, 0.0) < 0) {    // mass loss from donor?

        controlled_stepper_type controlled_stepper;
        state_type x(1);
        x[0] = semiMajorAxis;

        // Use boost adaptive ODE solver for speed and accuracy
        struct ode {
            double p_MassDonor0, p_MassAccretor0, p_FractionAccreted;
            bool   p_IsAccretorDegenerate;
            ode(double massDonor0, double massAccretor0, double fractionAccreted, bool isAccretorDegenerate) : p_MassDonor0(massDonor0), p_MassAccretor0(massAccretor0), p_FractionAccreted(fractionAccreted), p_IsAccretorDegenerate(isAccretorDegenerate) {}

            void operator()(state_type const& x, state_type& dxdt, double p_MassChange ) const {
                double massD = p_MassDonor0 + p_MassChange;
                double massA = p_MassAccretor0 - p_MassChange * p_FractionAccreted;
                double jLoss = CalculateGammaAngularMomentumLoss_Static(massD, massA, p_IsAccretorDegenerate);
                dxdt[0]      = (-2.0 / massD) * (1.0 - (p_FractionAccreted * (massD / massA)) - ((1.0 - p_FractionAccreted) * (jLoss + 0.5) * (massD / (massA + massD)))) * x[0];
            }
        };

        integrate_adaptive(controlled_stepper, ode{ p_DonorMass, p_Accretor.Mass(), p_FractionAccreted, p_Accretor.IsDegenerate() }, x, 0.0, p_DeltaMassDonor, p_DeltaMassDonor / 1000.0);
        semiMajorAxis = x[0];
    }
    
    return semiMajorAxis;
}



/*
 * Calculate the response of the donor Roche Lobe to mass loss during mass transfer per Sluys 2013, Woods et al., 2012
 *
 * Sluys 2013, eq 60, Woods et al., 2012
 * Formula from M. Sluys notes "Binary evolution in a nutshell"
 *
 *
 * double CalculateZetaRocheLobe(const double p_jLoss, const double p_beta) const
 *
 * @param   [IN]    p_jLoss                     Specific angular momentum with which mass is lost during non-conservative mass transfer
 *                                              (Podsiadlowski et al. 1992, Beta: specific angular momentum of matter [2Pia^2/P])
 * @param   [IN]    p_beta                      Fraction of donated mass that is accreted by the accretor
 * @return                                      Roche Lobe response
 */
double BaseBinaryStar::CalculateZetaRocheLobe(const double p_jLoss, const double p_beta) const {

    double donorMass    = m_Donor->Mass();                  // donor mass
    double accretorMass = m_Accretor->Mass();               // accretor mass
    double gamma        = p_jLoss;
    double q            = donorMass / accretorMass;
    double cbrt_q       = std::cbrt(q);

    double k1 = -2.0 * (1.0 - (p_beta * q) - (1.0 - p_beta) * (gamma + 0.5) * (q / (1.0 + q)));
    double k2 = (2.0 / 3.0) - cbrt_q * (1.2 * cbrt_q + 1.0 / (1.0 + cbrt_q)) / (3.0 * (0.6 * cbrt_q * cbrt_q + log(1.0 + cbrt_q)));
    double k3 = 1.0 + (p_beta * q);

    return k1 + (k2 * k3);
}


/*
 * Calculate mass loss due to winds for each star and apply loss
 *
 * JR: todo: flesh-out this documentation
 *
 *
 * void CalculateWindsMassLoss()
 */
void BaseBinaryStar::CalculateWindsMassLoss() {

    m_aMassLossDiff = 0.0;                                                                                                      // initially - no change to orbit (semi-major axis) due to winds mass loss

    // Halt mass loss due to winds if the binary is in mass transfer and set the Mdot parameters of both stars appropriately
    if (OPTIONS->UseMassTransfer() && m_MassTransfer) {
            m_Star1->SetMassLossDiff(0.0);                                                                                      // JR would prefer to avoid a Setter for aesthetic reasons
            m_Star2->SetMassLossDiff(0.0);
            m_Star1->HaltWinds();
            m_Star2->HaltWinds();
    }
    else {
        if (OPTIONS->UseMassLoss()) {                                                                                           // mass loss enabled?

            double mWinds1 = m_Star1->CalculateMassLossValues(true);                                                            // calculate new values assuming mass loss applied
            double mWinds2 = m_Star2->CalculateMassLossValues(true);                                                            // calculate new values assuming mass loss applied

            double aWinds  = m_SemiMajorAxisPrev / (2.0 - ((m_Star1->MassPrev() + m_Star2->MassPrev()) / (mWinds1 + mWinds2))); // new semi-major axis for circularlised orbit

            m_Star1->SetMassLossDiff(mWinds1 - m_Star1->Mass());                                                                // JR: todo: find a better way?
            m_Star2->SetMassLossDiff(mWinds2 - m_Star2->Mass());                                                                // JR: todo: find a better way?

            m_aMassLossDiff = aWinds - m_SemiMajorAxisPrev;                                                                     // change to orbit (semi-major axis) due to winds mass loss
        }
    }
}


/*
 *  Check if mass transfer should happen (either star, but not both, overflowing Roche Lobe)
 *  Perform mass transfer if required and update individual stars accordingly
 *
 *  Updates class member variables
 * 
 *
 * void CalculateMassTransfer(const double p_Dt)
 *
 * @param   [IN]    p_Dt                        timestep in Myr
 */
void BaseBinaryStar::CalculateMassTransfer(const double p_Dt) {
    
    InitialiseMassTransfer();                                                                                                   // initialise - even if not using mass transfer (sets some flags we might need)
    
    if (Unbound()) return;                                                                                                      // do nothing for unbound binaries
    
    if (!OPTIONS->UseMassTransfer()) return;                                                                                    // mass transfer not enabled - nothing to do
    
    if (!m_Star1->IsRLOF() && !m_Star2->IsRLOF()) return;                                                                       // neither star is overflowing its Roche Lobe - no mass transfer - nothing to do
    
    if (OPTIONS->CHEMode() != CHE_MODE::NONE && HasTwoOf({STELLAR_TYPE::CHEMICALLY_HOMOGENEOUS}) && HasStarsTouching()) {       // CHE enabled and both stars CH?
        m_Flags.stellarMerger = true;
        return;
    }
    
    if (HasOneOf({STELLAR_TYPE::MASSLESS_REMNANT})) return;                                                                     // one of the stars is already a massless remnant, nothing to do
    
    if (m_Star1->IsRLOF() && m_Star2->IsRLOF()) {                                                                               // both stars overflowing their Roche Lobe?
        m_CEDetails.CEEnow = true;                                                                                              // yes - common envelope event - no mass transfer
        return;                                                                                                                 // and return - nothing (else) to do
    }
    
    // one, and only one, star is overflowing its Roche Lobe - resolve mass transfer
    m_Donor    = m_Star2->IsRLOF() ? m_Star2 : m_Star1;                                                                         // donor is primary unless secondary is overflowing its Roche Lobe
    m_Accretor = m_Star2->IsRLOF() ? m_Star1 : m_Star2;                                                                         // accretor is secondary unless secondary is overflowing its Roche Lobe
    
    m_Donor->UpdateMassTransferDonorHistory();                                                                                  // add event to MT history of the donor
    
    // Calculate accretion fraction if stable
    // This passes the accretor's Roche lobe radius to m_Accretor->CalculateThermalMassAcceptanceRate()
    // just in case MT_THERMALLY_LIMITED_VARIATION::RADIUS_TO_ROCHELOBE is used; otherwise, the radius input is ignored
    double accretorRLradius = CalculateRocheLobeRadius_Static(m_Accretor->Mass(), m_Donor->Mass()) * AU_TO_RSOL * m_SemiMajorAxis * (1.0 - m_Eccentricity);
    bool donorIsHeRich      = m_Donor->IsOneOf(He_RICH_TYPES);
    
    double jLoss = m_JLoss;                                                                                                     // specific angular momentum with which mass is lost during non-conservative mass transfer, current timestep
    if (OPTIONS->MassTransferAngularMomentumLossPrescription() != MT_ANGULAR_MOMENTUM_LOSS_PRESCRIPTION::ARBITRARY) {           // arbitrary angular momentum loss prescription?
        jLoss = CalculateGammaAngularMomentumLoss();                                                                            // no - re-calculate angular momentum
    }
    
    double betaThermal              = 0.0;                                                                                      // fraction of mass accreted if accretion proceeds on thermal timescale
    double betaNuclear              = 0.0;                                                                                      // fraction of mass accreted if accretion proceeds on nuclear timescale
    double donorMassLossRateThermal = m_Donor->CalculateThermalMassLossRate();
    double donorMassLossRateNuclear = m_Donor->CalculateNuclearMassLossRate();
    
    std::tie(std::ignore, betaThermal) = m_Accretor->CalculateMassAcceptanceRate(donorMassLossRateThermal,
                                                                                 m_Accretor->CalculateThermalMassAcceptanceRate(accretorRLradius),
                                                                                 donorIsHeRich);
    std::tie(std::ignore, betaNuclear) = m_Accretor->CalculateMassAcceptanceRate(donorMassLossRateNuclear,
                                                                                 m_Accretor->CalculateThermalMassAcceptanceRate(accretorRLradius),
                                                                                 donorIsHeRich);
    
    m_ZetaStar             = m_Donor->CalculateZetaAdiabatic();
    double zetaEquilibrium = m_Donor->CalculateZetaEquilibrium();
    
    m_ZetaLobe = CalculateZetaRocheLobe(jLoss, betaNuclear);                                                                    // try nuclear timescale mass transfer first
    if(m_Donor->IsOneOf(ALL_MAIN_SEQUENCE) && utils::Compare(zetaEquilibrium, m_ZetaLobe) > 0) {
        m_MassLossRateInRLOF = donorMassLossRateNuclear;
        m_FractionAccreted   = betaNuclear;
    }
    else {
        m_ZetaLobe = CalculateZetaRocheLobe(jLoss, betaThermal);
        m_MassLossRateInRLOF = donorMassLossRateThermal;
        m_FractionAccreted   = betaThermal;
    }
        
    double aInitial = m_SemiMajorAxis;                                                                                          // semi-major axis in default units, AU, current timestep
    double aFinal;                                                                                                              // semi-major axis in default units, AU, after next timestep

    // Calculate conditions for automatic (in)stability for case BB
    bool caseBBAlwaysStable           = OPTIONS->CaseBBStabilityPrescription() == CASE_BB_STABILITY_PRESCRIPTION::ALWAYS_STABLE;
    bool caseBBAlwaysUnstable         = OPTIONS->CaseBBStabilityPrescription() == CASE_BB_STABILITY_PRESCRIPTION::ALWAYS_UNSTABLE;
    bool caseBBAlwaysUnstableOntoNSBH = OPTIONS->CaseBBStabilityPrescription() == CASE_BB_STABILITY_PRESCRIPTION::ALWAYS_STABLE_ONTO_NSBH;
    bool donorIsHeHGorHeGB            = m_Donor->IsOneOf({ STELLAR_TYPE::NAKED_HELIUM_STAR_HERTZSPRUNG_GAP, STELLAR_TYPE::NAKED_HELIUM_STAR_GIANT_BRANCH });
    bool accretorIsNSorBH             = m_Accretor->IsOneOf({ STELLAR_TYPE::NEUTRON_STAR, STELLAR_TYPE::BLACK_HOLE });
    bool accretorIsWD                 = m_Accretor->IsOneOf(WHITE_DWARFS); 

    // Determine stability
    bool isUnstable = false;
    if (donorIsHeHGorHeGB && (caseBBAlwaysStable || caseBBAlwaysUnstable || (caseBBAlwaysUnstableOntoNSBH && accretorIsNSorBH))) { // determine stability based on case BB 
        isUnstable = (caseBBAlwaysUnstable || (caseBBAlwaysUnstableOntoNSBH && accretorIsNSorBH));                              // already established that donor is HeHG or HeGB - need to check if new case BB prescriptions are added
    } 
    else if (accretorIsWD && (m_Accretor->WhiteDwarfAccretionRegime() == ACCRETION_REGIME::HELIUM_WHITE_DWARF_HYDROGEN_ACCUMULATION)) { 
        isUnstable = true;
        if (!m_Donor->IsOneOf(GIANTS)) m_Flags.stellarMerger = true;
    }
    else if (OPTIONS->QCritPrescription() != QCRIT_PRESCRIPTION::NONE) {                                                        // determine stability based on critical mass ratios
        // NOTE: Critical mass ratio is defined as mAccretor/mDonor
        double qCrit = m_Donor->CalculateCriticalMassRatio(m_Accretor->IsDegenerate(), m_FractionAccreted);
        isUnstable   = utils::Compare((m_Accretor->Mass() / m_Donor->Mass()), qCrit) < 0;
    }
    else {                                                                                                                      // determine stability based on zetas
        isUnstable   = (utils::Compare(m_ZetaStar, m_ZetaLobe) < 0);
    }

    // Evaluate separately for stable / unstable MT
    if (isUnstable) {                                                                                                           // unstable Mass Transfer
         m_CEDetails.CEEnow = true;
    }
    else {                                                                                                                      // stable MT
            
        m_MassTransferTrackerHistory = m_Donor == m_Star1                                                                       // record what happened - for later printing
            ? MT_TRACKING::STABLE_1_TO_2_SURV
            : MT_TRACKING::STABLE_2_TO_1_SURV; 

        double massDiffDonor;
        double envMassDonor    = m_Donor->Mass() - m_Donor->CoreMass();
        bool isEnvelopeRemoved = false;

        if (utils::Compare(m_Donor->CoreMass(), 0) > 0 && utils::Compare(envMassDonor, 0) > 0) {                                // donor has a core and an envelope?
            massDiffDonor     = -envMassDonor;                                                                                  // yes - set donor mass loss to (negative of) the envelope mass
            isEnvelopeRemoved = true;
        }
        else {                                                                                                                  // donor has no envelope
            massDiffDonor = MassLossToFitInsideRocheLobe(this, m_Donor, m_Accretor, m_FractionAccreted);                        // use root solver to determine how much mass should be lost from the donor to allow it to fit within the Roche lobe
            
            if (massDiffDonor <= 0.0) {                                                                                         // no root found
                // if donor cannot lose mass to fit inside Roche lobe, the only viable action is to enter CE phase
                m_CEDetails.CEEnow = true;                                                                                      // flag CE
            }
            else {                                                                                                              // have required mass loss
                if (utils::Compare(m_MassLossRateInRLOF,donorMassLossRateNuclear) == 0)                                         // if transferring mass on nuclear timescale, limit mass loss amount to rate * timestep (thermal timescale MT always happens in one timestep)
                    massDiffDonor = std::min(massDiffDonor, m_MassLossRateInRLOF * m_Dt);
                massDiffDonor = -massDiffDonor;                                                                                 // set mass difference
                m_Donor->UpdateMinimumCoreMass();                                                                               // reset the minimum core mass following case A
            }
        }

        if (!m_CEDetails.CEEnow) {                                                                                              // CE flagged?
                                                                                                                                // no
            double massGainAccretor = -massDiffDonor * m_FractionAccreted;                                                      // set accretor mass gain to mass loss * conservativeness

            m_Donor->SetMassTransferDiffAndResolveWDShellChange(massDiffDonor);                                                 // set new mass of donor
            m_Accretor->SetMassTransferDiffAndResolveWDShellChange(massGainAccretor);                                           // set new mass of accretor

            aFinal              = CalculateMassTransferOrbit(m_Donor->Mass(), massDiffDonor, *m_Accretor, m_FractionAccreted);  // calculate new orbit
            m_aMassTransferDiff = aFinal - aInitial;                                                                            // set change in orbit (semi-major axis)
                                                                                                                    
            STELLAR_TYPE stellarTypeDonor = m_Donor->StellarType();                                                             // donor stellar type before resolving envelope loss
            if (isEnvelopeRemoved) m_Donor->ResolveEnvelopeLossAndSwitch();                                                     // if this was an envelope stripping episode, resolve envelope loss
            if (m_Donor->StellarType() != stellarTypeDonor) {                                                                   // stellar type change?
                (void)PrintDetailedOutput(m_Id, BSE_DETAILED_RECORD_TYPE::STELLAR_TYPE_CHANGE_DURING_MT);                       // yes - print (log) detailed output
            }
        
            // Check if this was stable mass transfer after a CEE
            if (m_CEDetails.CEEcount > 0 && !m_RLOFDetails.stableRLOFPostCEE) {
                m_RLOFDetails.stableRLOFPostCEE = m_MassTransferTrackerHistory == MT_TRACKING::STABLE_2_TO_1_SURV ||
                                                  m_MassTransferTrackerHistory == MT_TRACKING::STABLE_1_TO_2_SURV;
            }
        }
    }
    
	// Check for recycled pulsars. Not considering CEE as a way of recycling NSs.
	if (!m_CEDetails.CEEnow && m_Accretor->IsOneOf({ STELLAR_TYPE::NEUTRON_STAR })) {                                           // accretor is a neutron star
        m_Donor->SetRLOFOntoNS();                                                                                               // donor donated mass to a neutron star
        m_Accretor->SetRecycledNS();                                                                                            // accretor is (was) a recycled NS
	}
}


/*
 * Setup parameters for mass transfer/common envelope event
 *
 *
 * void InitialiseMassTransfer()
 */
void BaseBinaryStar::InitialiseMassTransfer() {

	m_MassTransferTrackerHistory = MT_TRACKING::NO_MASS_TRANSFER;	                                                            // Initiating flag, every timestep, to NO_MASS_TRANSFER. If it undergoes to MT or CEE, it should change.

    m_Star1->InitialiseMassTransfer(m_CEDetails.CEEnow, m_SemiMajorAxis, m_Eccentricity);                                       // initialise mass transfer for star1
    m_Star2->InitialiseMassTransfer(m_CEDetails.CEEnow, m_SemiMajorAxis, m_Eccentricity);                                       // initialise mass transfer for star2
    
    if (m_Star1->IsRLOF() || m_Star2->IsRLOF()) {                                                                               // either star overflowing its Roche Lobe?
                                                                                                                                // yes - mass transfer if not both CH
        if (OPTIONS->CHEMode() != CHE_MODE::NONE && HasTwoOf({STELLAR_TYPE::CHEMICALLY_HOMOGENEOUS})) {                         // CHE enabled and both stars CH?
                                                                                                                                // yes
            // equilibrate masses and circularise (check for merger is done later)

            if (utils::Compare(m_Star1->Mass(), m_Star2->Mass()) != 0) {                                                        // masses already equal?
                                                                                                                                // no - make them equal
                STELLAR_TYPE stellarType1 = m_Star1->StellarType();                                                             // star 1 stellar type before updating attributes
                STELLAR_TYPE stellarType2 = m_Star2->StellarType();                                                             // star 2 stellar type before updating attributes

                double mass = (m_Star1->Mass() + m_Star2->Mass()) / 2.0;                                                        // share mass equally
                if ((m_Star1->UpdateAttributes(mass - m_Star1->Mass(), mass - m_Star1->Mass0(), true) != stellarType1) ||       // set new mass, mass0 for star 1
                    (m_Star2->UpdateAttributes(mass - m_Star2->Mass(), mass - m_Star2->Mass0(), true) != stellarType2)) {       // set new mass, mass0 for star 2
                    (void)PrintDetailedOutput(m_Id, BSE_DETAILED_RECORD_TYPE::STELLAR_TYPE_CHANGE_DURING_CHE_EQUILIBRATION);    // print (log) detailed output if stellar type changed
                }
                m_Flags.massesEquilibrated = true;                                                                              // record that we've equilbrated
            }

            // circularise if not already
            if (utils::Compare(m_Eccentricity, 0.0) != 0) {                                                                     // eccentricity = 0.0?
                                                                                                                                // no - circularise
                // conserve angular momentum
                // use J = m1 * m2 * sqrt(G * a * (1 - e^2) / (m1 + m2))

                double M         = m_Star1->Mass() + m_Star2->Mass();
                double m1m2      = m_Star1->Mass() * m_Star2->Mass();
                m_SemiMajorAxis *= 16.0 * m1m2 * m1m2 / (M * M * M * M) * (1.0 - (m_Eccentricity * m_Eccentricity));            // circularise; conserve angular momentum
                m_Eccentricity   = 0.0;                                                                                         // now circular
            }
            
            m_Star1->InitialiseMassTransfer(m_CEDetails.CEEnow, m_SemiMajorAxis, m_Eccentricity);                               // re-initialise mass transfer for star1
            m_Star2->InitialiseMassTransfer(m_CEDetails.CEEnow, m_SemiMajorAxis, m_Eccentricity);                               // re-initialise mass transfer for star2

            m_MassTransfer     = false;                                                                                         // no mass transfer
            m_CEDetails.CEEnow = false;                                                                                         // no common envelope
        }
        else {                                                                                                                  // not both CH, so ...
		    m_MassTransfer     = true;                                                                                          // ... mass transfer
            m_CEDetails.CEEnow = false;                                                                                         // no common envelope

		    if (OPTIONS->CirculariseBinaryDuringMassTransfer()) {                                                               // circularise binary
                m_SemiMajorAxis *= OPTIONS->AngularMomentumConservationDuringCircularisation()                                  // yes - conserve angular momentum?
                                        ? (1.0 - (m_Eccentricity * m_Eccentricity))                                             // yes - conserve angular momentum
                                        : (1.0 - m_Eccentricity);                                                               // no - angular momentum not conserved, circularise at periapsis

			    m_Eccentricity = 0.0;

                m_Star1->InitialiseMassTransfer(m_CEDetails.CEEnow, m_SemiMajorAxis, m_Eccentricity);                           // re-initialise mass transfer for star1
                m_Star2->InitialiseMassTransfer(m_CEDetails.CEEnow, m_SemiMajorAxis, m_Eccentricity);                           // re-initialise mass transfer for star2
                
			    // Update previous timestep values to those of the circularised binary to serve as a baseline for future updates.
			    m_SemiMajorAxisPrev = m_SemiMajorAxis;
			    m_EccentricityPrev  = m_Eccentricity;
		    }
        }
    }
    else {
        m_MassTransfer     = false;                                                                                             // no mass transfer
        m_CEDetails.CEEnow = false;                                                                                             // no common envelope
    }

    m_aMassTransferDiff = 0.0;                                                                                                  // iniitially - no change to orbit (semi-major axis) due to mass transfer
}


/*
 * Calculate the total energy of the binary
 *
 * The energy consists of the spin kinetic energies of the two stars, the kinetic energy of the binary, and the gravitational potential energy of the binary
 *
 *
 * double CalculateTotalEnergy(const double p_SemiMajorAxis,
 *                             const double p_Star1Mass,
 *                             const double p_Star2Mass,
 *                             const double p_Star1SpinAngularVelocity,
 *                             const double p_Star2SpinAngularVelocity,
 *                             const double p_Star1MomentOfInertia,
 *                             const double p_Star2MomentOfInertia)
 *
 * @param   [IN]    p_SemiMajorAxis             Semi-major axis of the binary
 * @param   [IN]    p_Star1Mass                 Mass of star 1
 * @param   [IN]    p_Star2Mass                 Mass of star 2
 * @param   [IN]    p_Star1SpinAngularVelocity  Spin angular velocity of star 1
 * @param   [IN]    p_Star2SpinAngularVelocity  Spin angular velocity of star 2
 * @param   [IN]    p_Star1MomentOfInertia      Moment of inertia of star 1
 * @param   [IN]    p_Star2MomentOfInertia      Moment of inertia of star 2
 * @return                                      Total energy of the binary
 */
double BaseBinaryStar::CalculateTotalEnergy(const double p_SemiMajorAxis,
                                            const double p_Star1Mass,
                                            const double p_Star2Mass,
                                            const double p_Star1SpinAngularVelocity,
                                            const double p_Star2SpinAngularVelocity,
                                            const double p_Star1MomentOfInertia,
                                            const double p_Star2MomentOfInertia) const {

	double w1_2 = p_Star1SpinAngularVelocity * p_Star1SpinAngularVelocity;
	double w2_2 = p_Star2SpinAngularVelocity * p_Star2SpinAngularVelocity;

	return 0.5 * ((p_Star1MomentOfInertia * w1_2) + (p_Star2MomentOfInertia * w2_2) - (G_AU_Msol_yr * p_Star1Mass * p_Star2Mass / p_SemiMajorAxis));
}


/*
 * Calculate the angular momentum of the binary
 *
 * The angular momentum consists of the spin angular momenta of the two stars and the orbital angular momentum of the binary
 *
 *
 * double CalculateAngularMomentum(const double p_SemiMajorAxis,
 *                                 const double p_Eccentricity,
 *                                 const double p_Star1Mass,
 *                                 const double p_Star2Mass,
 *                                 const double p_Star1SpinAngularVelocity,
 *                                 const double p_Star1SpinAngularVelocity,
 *                                 const double p_Star1MomentOfInertia,
 *                                 const double p_Star2MomentOfInertia)
 *
 * @param   [IN]    p_SemiMajorAxis             Semi-major axis of the binary
 * @param   [IN]    p_Eccentricity              Eccentricity of the binary
 * @param   [IN]    p_Star1Mass                 Mass of the primary
 * @param   [IN]    p_Star2Mass                 Mass of the secondary
 * @param   [IN]    p_Star1SpinAngularVelocity  Orbital frequency of the primary
 * @param   [IN]    p_Star1SpinAngularVelocity  Orbital frequency of the secondary
 * @param   [IN]    p_Star1MomentOfInertia      Moment of inertia of the primary
 * @param   [IN]    p_Star2MomentOfInertia      Moment of inertia of the secondary
 * @return                                      Angular momentum of the binary
 */
double BaseBinaryStar::CalculateAngularMomentum(const double p_SemiMajorAxis,
                                                const double p_Eccentricity,
                                                const double p_Star1Mass,
                                                const double p_Star2Mass,
                                                const double p_Star1SpinAngularVelocity,
                                                const double p_Star2SpinAngularVelocity,
                                                const double p_Star1MomentOfInertia,
                                                const double p_Star2MomentOfInertia) const {

    double Jorb = CalculateOrbitalAngularMomentum(p_Star1Mass, p_Star2Mass, p_SemiMajorAxis, p_Eccentricity);

	return (p_Star1MomentOfInertia * p_Star1SpinAngularVelocity) + (p_Star2MomentOfInertia * p_Star2SpinAngularVelocity) + Jorb;
}


/*
 * Calculate total energy and angular momentum of the binary
 *
 * Calls CalculateTotalEnergy() and CalculateAngularMomentum()
 * Updates class member variables
 *
 *
 * void CalculateEnergyAndAngularMomentum()
 */
void BaseBinaryStar::CalculateEnergyAndAngularMomentum() {

    if (m_Star1->IsOneOf({ STELLAR_TYPE::MASSLESS_REMNANT }) || m_Star2->IsOneOf({ STELLAR_TYPE::MASSLESS_REMNANT })) return;

    // Calculate orbital energy and angular momentum
    m_OrbitalEnergyPrev          = m_OrbitalEnergy;
    m_OrbitalAngularMomentumPrev = m_OrbitalAngularMomentum;
    m_TotalAngularMomentumPrev   = m_TotalAngularMomentum;

	double totalMass             = m_Star1->Mass() + m_Star2->Mass();
	double reducedMass           = (m_Star1->Mass() * m_Star2->Mass()) / totalMass;

    m_OrbitalEnergy              = CalculateOrbitalEnergy(reducedMass, totalMass, m_SemiMajorAxis);
    m_OrbitalAngularMomentum     = CalculateOrbitalAngularMomentum(m_Star1->Mass(), m_Star2->Mass(), m_SemiMajorAxis, m_Eccentricity);

    // Calculate total energy and angular momentum using regular conservation of energy, especially useful for checking tides and rotational effects
    m_TotalEnergy                = CalculateTotalEnergy();
    m_TotalAngularMomentum       = CalculateAngularMomentum();
}


/*
 * Resolve mass changes
 *
 * Applies mass changes to both stars
 * Updates attributes of both stars in response to mass changes
 * Calculates orbital velocity and semi-major axis of binary after mass changes
 * Calculate total energy and angular momentum of binary after mass changes
 *
 *
 * void ResolveMassChanges()
 *
 */
void BaseBinaryStar::ResolveMassChanges() {

    STELLAR_TYPE stellarType1 = m_Star1->StellarTypePrev();                                             // star 1 stellar type before updating attributes
    STELLAR_TYPE stellarType2 = m_Star2->StellarTypePrev();                                             // star 2 stellar type before updating attributes

    // update mass of star1 according to mass loss and mass transfer, then update age accordingly
    (void)m_Star1->UpdateAttributes(m_Star1->MassPrev() - m_Star1->Mass() + m_Star1->MassLossDiff() + m_Star1->MassTransferDiff(), 0.0); // update mass for star1
    m_Star1->UpdateInitialMass();                                                                       // update effective initial mass of star1 (MS, HG & HeMS)
    m_Star1->UpdateAgeAfterMassLoss();                                                                  // update age of star1
    m_Star1->ApplyMassTransferRejuvenationFactor();                                                     // apply age rejuvenation factor for star1
    m_Star1->UpdateAttributes(0.0, 0.0, true);

    // rinse and repeat for star2
    (void)m_Star2->UpdateAttributes(m_Star2->MassPrev() - m_Star2->Mass() + m_Star2->MassLossDiff() + m_Star2->MassTransferDiff(), 0.0); // update mass for star2
    m_Star2->UpdateInitialMass();                                                                       // update effective initial mass of star 2 (MS, HG & HeMS)
    m_Star2->UpdateAgeAfterMassLoss();                                                                  // update age of star2
    m_Star2->ApplyMassTransferRejuvenationFactor();                                                     // apply age rejuvenation factor for star2
    m_Star2->UpdateAttributes(0.0, 0.0, true);
    
    // update binary separation, but only if semimajor axis not already infinite and binary does not contain a massless remnant
    if(!isinf(m_SemiMajorAxis) && !HasOneOf({STELLAR_TYPE::MASSLESS_REMNANT}))
        m_SemiMajorAxis = m_SemiMajorAxisPrev + m_aMassLossDiff + m_aMassTransferDiff;
    
    //Envelope ejection for convective envelope stars exceeding threshold luminosity to mass ratio: assume the entire envelope was lost on timescales long relative to the orbit
    if (m_Star1->EnvelopeJustExpelledByPulsations() || m_Star2->EnvelopeJustExpelledByPulsations()) {
        m_SemiMajorAxis /= (2.0 - ((m_Star1->MassPrev() + m_Star2->MassPrev()) / (m_Star1->Mass() + m_Star2->Mass()))); // update separation in response to pulsational mass loss
        m_Star1->ResetEnvelopeExpulsationByPulsations();
        m_Star2->ResetEnvelopeExpulsationByPulsations();
    }

    // if CHE enabled, update rotational frequency for constituent stars - assume tidally locked
    double omega = OrbitalAngularVelocity();                                                           // orbital angular velocity
    if (OPTIONS->CHEMode() != CHE_MODE::NONE) m_Star1->SetOmega(omega);
    if (OPTIONS->CHEMode() != CHE_MODE::NONE) m_Star2->SetOmega(omega);

    // if both stars evolving as chemically homogeneous stars set m_Omega for binary
    if (HasTwoOf({STELLAR_TYPE::CHEMICALLY_HOMOGENEOUS})) m_Omega = omega;

    CalculateEnergyAndAngularMomentum();                                                                // perform energy and angular momentum calculations

    if ((m_Star1->StellarType() != stellarType1) || (m_Star2->StellarType() != stellarType2)) {         // stellar type change?
        (void)PrintDetailedOutput(m_Id, BSE_DETAILED_RECORD_TYPE::STELLAR_TYPE_CHANGE_DURING_MASS_RESOLUTION); // yes - print (log) detailed output
    }
}


/*
 * Process tides if required
 *
 * 
 * void BaseBinaryStar::ProcessTides(const double p_Dt)
 *
 * @param   [in]        p_Dt                    Timestep (in Myr)
 */
void BaseBinaryStar::ProcessTides(const double p_Dt) {

    if (!m_Unbound) {                                                                                                           // binary bound?
                                                                                                                                // yes - process tides if enabled
        if (OPTIONS->TidesPrescription() != TIDES_PRESCRIPTION::NONE) {                                                         // tides enabled?

            // if m_Omega == 0.0 (should only happen on the first timestep), calculate m_Omega here
            if (utils::Compare(m_Omega, 0.0) == 0) m_Omega = OrbitalAngularVelocity();
        }

        switch (OPTIONS->TidesPrescription()) {                                                                                 // which tides prescription?
            case TIDES_PRESCRIPTION::NONE: break;                                                                               // NONE - tides not enabled - do nothing
        
            case TIDES_PRESCRIPTION::KAPIL2024: {                                                                               // KAPIL2024

                // Evolve binary semi-major axis, eccentricity, and spin of each star based on Kapil et al., 2024

                DBL_DBL_DBL_DBL ImKlm1   = m_Star1->CalculateImKlmTidal(m_Omega, m_SemiMajorAxis, m_Star2->Mass());
                DBL_DBL_DBL_DBL ImKlm2   = m_Star2->CalculateImKlmTidal(m_Omega, m_SemiMajorAxis, m_Star1->Mass());

                double DSemiMajorAxis1Dt = CalculateDSemiMajorAxisTidalDt(ImKlm1, m_Star1);                                     // change in semi-major axis from star1
                double DSemiMajorAxis2Dt = CalculateDSemiMajorAxisTidalDt(ImKlm2, m_Star2);                                     // change in semi-major axis from star2

                double DEccentricity1Dt  = CalculateDEccentricityTidalDt(ImKlm1, m_Star1);                                      // change in eccentricity from star1
                double DEccentricity2Dt  = CalculateDEccentricityTidalDt(ImKlm2, m_Star2);                                      // change in eccentricity from star2

                double DOmega1Dt         = CalculateDOmegaTidalDt(ImKlm1, m_Star1);                                             // change in spin from star1
                double DOmega2Dt         = CalculateDOmegaTidalDt(ImKlm2, m_Star2);                                             // change in spin from star2

                m_Star1->SetOmega(m_Star1->Omega() + (DOmega1Dt * p_Dt * MYR_TO_YEAR));                                         // evolve star 1 spin
                m_Star2->SetOmega(m_Star2->Omega() + (DOmega2Dt * p_Dt * MYR_TO_YEAR));                                         // evolve star 2 spin

                m_SemiMajorAxis          = m_SemiMajorAxis + ((DSemiMajorAxis1Dt + DSemiMajorAxis2Dt) * p_Dt * MYR_TO_YEAR);    // evolve separation
                m_Eccentricity           = m_Eccentricity + ((DEccentricity1Dt + DEccentricity2Dt) * p_Dt * MYR_TO_YEAR);       // evolve eccentricity 
                m_Omega                  = OrbitalAngularVelocity();                                                            // re-calculate orbital frequency
                m_TotalAngularMomentum   = CalculateAngularMomentum();                                                          // re-calculate total angular momentum

                } break;

            case TIDES_PRESCRIPTION::PERFECT: {                                                                                 // PERFECT

                // find omega assuming instantaneous synchronisation
                // use current value of m_Omega as best guess for root

                m_Omega = OmegaAfterSynchronisation(m_Star1->Mass(), m_Star2->Mass(), m_Star1->CalculateMomentOfInertiaAU(), m_Star2->CalculateMomentOfInertiaAU(), m_TotalAngularMomentum, m_Omega);

                if (m_Omega >= 0.0) {                                                                                           // root found?
                                                                                                                                // yes
                    m_Star1->SetOmega(m_Omega);                                                                                 // synchronise star 1
                    m_Star2->SetOmega(m_Omega);                                                                                 // synchronise star 2

                    m_SemiMajorAxis        = std::cbrt(G_AU_Msol_yr * (m_Star1->Mass() + m_Star2->Mass()) / m_Omega / m_Omega); // re-calculate semi-major axis
                    m_Eccentricity         = 0.0;                                                                               // circularise
                    m_TotalAngularMomentum = CalculateAngularMomentum();                                                        // re-calculate total angular momentum
                }
                else {                                                                                                          // no (real) root found

                    // no real root found - push the binary to a common envelope
                    // place the constituent star closest to RLOF at RLOF and use that to
                    // calculate semi-major axis, then use that to calculate m_Omega

                    double ratio1 = m_Star1->StarToRocheLobeRadiusRatio(m_SemiMajorAxis, m_Star1->Mass());                      // star 1 ratio of radius to Roche lobe radius
                    double ratio2 = m_Star2->StarToRocheLobeRadiusRatio(m_SemiMajorAxis, m_Star2->Mass());                      // star 2 ratio of radius to Roche lobe radius

                    double radius;
                    double mass1;
                    double mass2;
                    if (ratio1 >= ratio2) {                                                                                     // star 1 closer to RLOF than star 2 (or same)?
                        radius = m_Star1->Radius();                                                                             // yes - use star 1 to calculate semi-major axis at RLOF
                        mass1  = m_Star1->Mass();
                        mass2  = m_Star2->Mass();
                    }
                    else {                                                                                                      // no - star 2 closer to RLOF than star 1
                        radius = m_Star2->Radius();                                                                             // use star 2 to calculate semi-major axis at RLOF
                        mass1  = m_Star2->Mass();
                        mass2  = m_Star1->Mass();
                    }
            
                    m_Eccentricity  = 0.0;                                                                                      // assume circular
                    m_SemiMajorAxis = radius * RSOL_TO_AU / CalculateRocheLobeRadius_Static(mass1, mass2);                      // new semi-major axis - should tip into CE
                    m_Omega         = OrbitalAngularVelocity();                                                                 // m_Omega at new semi-major axis
                }
                } break;

            default:
                // the only way this can happen is if someone added a TIDES_PRESCRIPTION
                // and it isn't accounted for in this code.  We should not default here, with or without a warning.
                // We are here because the user chose a prescription this code doesn't account for, and that should
                // be flagged as an error and result in termination of the evolution of the binary.
                // The correct fix for this is to add code for the missing prescription or, if the missing
                // prescription is superfluous, remove it from the option.

                THROW_ERROR(ERROR::UNKNOWN_TIDES_PRESCRIPTION);                                                                 // throw error
        }
    }
}


/*
 * Calculate and emit gravitational radiation.
 *
 * This function uses Peters 1964 to approximate the effects of GW emission with two steps:
 * - Calculate the change in semi-major axis (m_SemiMajorAxis) per time given by eq 5.6.
 * - Calculate the change in eccentricity (m_Eccentricity) per time given by eq 5.7.
 * 
 * m_DaDtGW and m_DeDtGW are updated so that they can be used to calculate the timestep dynamically.
 * 
 *
 * void CalculateGravitationalRadiation()
 */
void BaseBinaryStar::CalculateGravitationalRadiation() {

    // Useful values
    double eccentricitySquared = m_Eccentricity * m_Eccentricity;
    double oneMinusESq         = 1.0 - eccentricitySquared;
    double oneMinusESq_5       = oneMinusESq * oneMinusESq * oneMinusESq * oneMinusESq * oneMinusESq;
    double G_AU_Msol_yr_3      = G_AU_Msol_yr * G_AU_Msol_yr * G_AU_Msol_yr;
    double C_AU_Yr_5           = C_AU_yr * C_AU_yr * C_AU_yr * C_AU_yr * C_AU_yr;
    double m_SemiMajorAxis_3   = m_SemiMajorAxis * m_SemiMajorAxis * m_SemiMajorAxis;
    double massAndGAndCTerm    = G_AU_Msol_yr_3 * m_Star1->Mass() * m_Star2->Mass() * (m_Star1->Mass() + m_Star2->Mass()) / C_AU_Yr_5;						// G^3 * m1 * m2(m1 + m2) / c^5 in units of Msol, AU and yr

    // Approximate rate of change in semimajor axis
    double numeratorA   = -64.0 * massAndGAndCTerm;
    double denominatorA = 5.0 * m_SemiMajorAxis_3 * std::sqrt(oneMinusESq_5 * oneMinusESq * oneMinusESq);
    m_DaDtGW            = (numeratorA / denominatorA) * (1.0 + (73.0 / 24.0) * eccentricitySquared + (37.0 / 96.0) * eccentricitySquared * eccentricitySquared) * MYR_TO_YEAR;  // units of AU Myr^-1

    // Approximate rate of change in eccentricity
    double numeratorE   = -304.0 * m_Eccentricity * massAndGAndCTerm;
    double denominatorE = 15.0 * m_SemiMajorAxis_3 * m_SemiMajorAxis * std::sqrt(oneMinusESq_5);
    m_DeDtGW            = (numeratorE / denominatorE) * (1.0 + (121.0 / 304.0) * eccentricitySquared) * YEAR_TO_MYR;									// units of Myr^-1
}


/*
 * Emit a GW based on the effects calculated by BaseBinaryStar::CalculateGravitationalRadiation().
 * 
 * This function updates the semi-major axis, eccentricity, and previous eccentricity values
 * (m_SemiMajorAxis, m_Eccentricity, m_SemiMajorAxisPrev, and m_EccentricityPrev) as a result of emitting GWs.
 * 
 *
 * void EmitGravitationalRadiation(const double p_Dt)
 *
 * @param   [IN]    p_Dt                        timestep in Myr
 */
void BaseBinaryStar::EmitGravitationalWave(const double p_Dt) {

    // Update semimajor axis
    double aNew     = m_SemiMajorAxis + (m_DaDtGW * p_Dt);
    m_SemiMajorAxis = utils::Compare(aNew, 0.0) > 0 ? aNew : 1E-20;  // if <0, set to arbitrarily small number

    // Update the eccentricity
    m_Eccentricity += m_DeDtGW * p_Dt;

    // Save values as previous timestep	
    m_SemiMajorAxisPrev = m_SemiMajorAxis;	
    m_EccentricityPrev  = m_Eccentricity;
}


/* 
 * Choose a timestep based on the parameters of the binary.
 *
 * This function will return the minimum of (i) a timestep based on the
 * orbital timescale of the binary and (ii) (if configured to emit GWs)
 * a timestep based on the magnitude of gravitational radiation.
 * 
 *
 * double ChooseTimestep(const double p_Dt)
 * 
 * @param   [IN]    p_Dt                        previous timestep in Myr
 * @return                                      new timestep in Myr
 */
double BaseBinaryStar::ChooseTimestep(const double p_Dt) {

    double newDt = std::min(m_Star1->CalculateTimestep(), m_Star2->CalculateTimestep());        			// new timestep

    if (OPTIONS->EmitGravitationalRadiation()) {                                                                        // emitting GWs?
        newDt = std::min(newDt, -1.0E-2 * m_SemiMajorAxis / m_DaDtGW);                                                  // reduce timestep if necessary to ensure that the orbital separation does not change by more than ~1% per timestep due to GW emission
    }

    newDt *= OPTIONS->TimestepMultiplier();	

    return std::max(std::round(newDt / TIMESTEP_QUANTUM) * TIMESTEP_QUANTUM, NUCLEAR_MINIMUM_TIMESTEP);                // quantised and not less than minimum
}


/*
 * Evaluate the binary system
 *
 *    - calculate any mass transfer
 *    - calculate mass loss due to winds
 *    - resolve any Common Envelope Event
 *    - resolve any Supernova Event
 *    - resolve mass changes - apply mass loss and mass transfer
 *    - resolve tidal interactions
 *    - calculate total energy and angular momentum after mass changes
 *    - update pulsar parameters
 *
 * 
 * void EvaluateBinary(const double p_Dt)
 *
 * @param   [in]        p_Dt                    Timestep (in Myr)
 */
void BaseBinaryStar::EvaluateBinary(const double p_Dt) {

    CalculateMassTransfer(p_Dt);                                                                                        // calculate mass transfer if necessary

    (void)PrintDetailedOutput(m_Id, BSE_DETAILED_RECORD_TYPE::POST_MT);                                                 // print (log) detailed output

    CalculateWindsMassLoss();                                                                                           // calculate mass loss dues to winds

    (void)PrintDetailedOutput(m_Id, BSE_DETAILED_RECORD_TYPE::POST_WINDS);                                              // print (log) detailed output

    if ((m_CEDetails.CEEnow || StellarMerger()) &&                                                                      // CEE or merger?
        !(OPTIONS->CHEMode() != CHE_MODE::NONE && HasTwoOf({STELLAR_TYPE::CHEMICALLY_HOMOGENEOUS}))
        && !HasOneOf({STELLAR_TYPE::MASSLESS_REMNANT}) ) {                                                              // yes - avoid CEE if CH+CH or one star is a massless remnant

        ResolveCommonEnvelopeEvent();                                                                                   // resolve CEE - immediate event
        (void)PrintDetailedOutput(m_Id, BSE_DETAILED_RECORD_TYPE::POST_CEE);                                            // print (log) detailed output
    }
    else if (m_Star1->IsSNevent() || m_Star2->IsSNevent()) {
        EvaluateSupernovae();                                                                                           // evaluate supernovae (both stars) - immediate event
        (void)PrintDetailedOutput(m_Id, BSE_DETAILED_RECORD_TYPE::POST_SN);                                             // print (log) detailed output
        if (HasOneOf({ STELLAR_TYPE::NEUTRON_STAR })) {
            (void)PrintPulsarEvolutionParameters(PULSAR_RECORD_TYPE::POST_SN);                                          // print (log) pulsar evolution parameters 
        }
    }
    else {
        ResolveMassChanges();                                                                                           // apply mass loss and mass transfer as necessary
        (void)PrintDetailedOutput(m_Id, BSE_DETAILED_RECORD_TYPE::POST_MASS_RESOLUTION);                                // print (log) detailed output

        if (HasStarsTouching()) {                                                                                       // if stars emerged from mass transfer as touching, it's a merger
            m_Flags.stellarMerger = true;
		
            // Set Roche lobe flags for both stars so that they show correct RLOF status
            m_Star1->SetRocheLobeFlags(m_CEDetails.CEEnow, m_SemiMajorAxis, m_Eccentricity);                            // set Roche lobe flags for star1
            m_Star2->SetRocheLobeFlags(m_CEDetails.CEEnow, m_SemiMajorAxis, m_Eccentricity);                            // set Roche lobe flags for star2
            (void)PrintDetailedOutput(m_Id, BSE_DETAILED_RECORD_TYPE::POST_MASS_RESOLUTION_MERGER);                     // print (log) detailed output
        }
    }

    if ((m_Star1->IsSNevent() || m_Star2->IsSNevent())) {
        EvaluateSupernovae();                                                                                           // evaluate supernovae (both stars) if mass changes are responsible for a supernova
        (void)PrintDetailedOutput(m_Id, BSE_DETAILED_RECORD_TYPE::POST_SN);                                             // print (log) detailed output
        if (HasOneOf({ STELLAR_TYPE::NEUTRON_STAR })) {
            (void)PrintPulsarEvolutionParameters(PULSAR_RECORD_TYPE::POST_SN);                                          // print (log) pulsar evolution parameters 
        }
    }

    CalculateEnergyAndAngularMomentum();                                                                                // perform energy and angular momentum calculations

    ProcessTides(p_Dt);                                                                                                 // process tides if required

    // assign new values to "previous" values, for following timestep
    m_EccentricityPrev  = m_Eccentricity;
    m_SemiMajorAxisPrev = m_SemiMajorAxis;

    m_Star1->UpdateMagneticFieldAndSpin(m_CEDetails.CEEnow, m_Dt * MYR_TO_YEAR * SECONDS_IN_YEAR, EPSILON_PULSAR);      // update pulsar parameters for star1
    m_Star2->UpdateMagneticFieldAndSpin(m_CEDetails.CEEnow, m_Dt * MYR_TO_YEAR * SECONDS_IN_YEAR, EPSILON_PULSAR);      // update pulsar parameters for star2
}


/*
 * Set parameters required before evolving one timestep - modify binary attributes
 *
 *
 * void EvolveOneTimestepPreamble(const double p_Dt)
 *
 * @param   [IN]    p_Dt                        Timestep
 */
void BaseBinaryStar::EvolveOneTimestepPreamble(const double p_Dt) {

    if (p_Dt > 0.0) {           // if dt > 0    (don't use utils::Compare() here)
        m_TimePrev = m_Time;    // Remember current simulation time
        m_Time    += p_Dt;      // Advance physical simulation time
        m_Dt       = p_Dt;      // Set timestep
    }
}


/*
 * Evolve the constituent stars of the binary a single timestep - timestep is provided.
 * Each individual star is aged for the same timestep
 *
 * See AgeOneTimestep() documentation in Star.cpp for details
 *
 * We catch any exceptions thrown by the SSE code here, set the binary error value as
 * necessary, and return the error value to the caller.  We don't rethrow exceptions
 * here - we expect th ecaller to examine the error return and do whatever is required
 * to manage any errors.
 * 
 *
 * ERROR EvolveOneTimestep(const double p_Dt, const int p_LogFileId)
 *
 * @param   [IN]    p_Dt                        The suggested timestep to evolve
 * @return                                      Error value
 */
ERROR BaseBinaryStar::EvolveOneTimestep(const double p_Dt) {

    EvolveOneTimestepPreamble(p_Dt);

    try {    

        m_Star1->AgeOneTimestep(p_Dt, true);                                            // Age the primary one timestep and switch to the new stellar type if necessary
        m_Star2->AgeOneTimestep(p_Dt, true);                                            // Age the secondary one timestep and switch to the new stellar type if necessary
    }

    // if we catch an error here it happened during the SSE evolution of one of the
    // constituent stars.  The error may have been displayed to the user already, but
    // the binary error value (m_Error) will not have been set - we set it here so we
    // know an error has occurred.

    catch (const std::runtime_error& e) {                                               // catch runtime exceptions
        // anything we catch here should not already have been displayed to the user,
        // so set and display the error (do not rethrow the error)
        if (std::string(e.what()) == "FPE") m_Error = ERROR::FLOATING_POINT_ERROR;      // set error value - floating-point error
        else                                m_Error = ERROR::ERROR;                     // set error value - unspecified error
        SHOW_ERROR(m_Error);                                                            // display the error
    }
    catch (int e) {                                                                     // catch errors thrown
        // anything we catch here should already have been displayed to the user,
        // so just set the error (do not rethrow the error)
        if (e != static_cast<int>(ERROR::NONE)) m_Error = static_cast<ERROR>(e);        // set error value - specified errpr
        else                                    m_Error = ERROR::ERROR;                 // set error value - unspecified error
    }
    catch (...) {                                                                       // catchall
        // anything we catch here should not already have been displayed to the user,
        // so set and display the error (do not rethrow the error)
        m_Error = ERROR::ERROR;                                                         // set error value - unspecified error
        SHOW_ERROR(m_Error);                                                            // unspecified error
    }

    return m_Error;
}


/*
 * Evolve the binary up to the maximum evolution time (and number of steps)
 *
 * The functional return is the status of the evolution (will indicate why the evolution stopped, and if an error occurred)
 *
 * JR: todo: flesh-out this documentation
 *
 *
 * EVOLUTION_STATUS Evolve()
 *
 * @return                                      Status of the evolution (EVOLUTION_STATUS)
 */
EVOLUTION_STATUS BaseBinaryStar::Evolve() {

    ERROR error = ERROR::NONE;

    EVOLUTION_STATUS evolutionStatus = EVOLUTION_STATUS::CONTINUE;

    try {

        if (HasStarsTouching()) {                                                                                                       // check if stars are touching
            m_Flags.stellarMerger        = true;
            m_Flags.stellarMergerAtBirth = true;
            evolutionStatus              = EVOLUTION_STATUS::STELLAR_MERGER_AT_BIRTH;                                                   // binary components are touching - merger at birth
        }

        (void)PrintDetailedOutput(m_Id, BSE_DETAILED_RECORD_TYPE::INITIAL_STATE);                                                       // print (log) detailed output: this is the initial state of the binary

        if (OPTIONS->PopulationDataPrinting()) {
            SAY("\nGenerating a new binary - " << m_Id);
            SAY("Binary has masses " << m_Star1->Mass() << " & " << m_Star2->Mass() << " Msol");
            SAY("Binary has initial semiMajorAxis " << m_SemiMajorAxis << " AU");
            SAY("RandomSeed " << m_RandomSeed);
        }

        if (evolutionStatus == EVOLUTION_STATUS::CONTINUE) {                                                                            // continue evolution

            // if the user provided timestep values, get them now
            bool usingProvidedTimesteps = false;                                                                                        // using user-provided timesteps?
            DBL_VECTOR timesteps;
            if (!OPTIONS->TimestepsFileName().empty()) {                                                                                // have timesteps filename?
                                                                                                                                        // yes
                std::tie(error, timesteps) = utils::ReadTimesteps(OPTIONS->TimestepsFileName());                                        // read timesteps from file
                if (error != ERROR::NONE) {                                                                                             // ok?
                    THROW_ERROR(error, ERR_MSG(ERROR::NO_TIMESTEPS_READ));                                                              // no - throw error - this is not what the user asked for
                }
                else usingProvidedTimesteps = true;                                                                                     // have user-provided timesteps
            }

            // evolve the current binary up to the maximum evolution time (and number of steps)

            double dt;                                                                                                                  // timestep
            unsigned long int stepNum = 0;                                                                                              // initialise step number

            while (evolutionStatus == EVOLUTION_STATUS::CONTINUE) {                                                                     // perform binary evolution - iterate over timesteps until told to stop

                stepNum++;                                                                                                              // increment stepNum

                // if user selects to emit GWs, calculate the effects of radiation
                //     - note that this is placed before ChooseTimestep() is called because
                //       the timestep is a function of graviational radiation
                if (OPTIONS->EmitGravitationalRadiation()) {
                    CalculateGravitationalRadiation();
                }

                if (stepNum > 1) {                                                                                                      // after the first timestep, set previous timestep
                    m_Star2->UpdatePreviousTimestepDuration();
                    m_Star1->UpdatePreviousTimestepDuration();
                }
                if (usingProvidedTimesteps) {                                                                                           // user-provided timesteps?
                    // select a timestep
                    //   - don't quantise
                    //   - don't apply timestep multiplier
                    // (we assume user wants the timesteps in the file)
                    dt = timesteps[stepNum - 1];
                }
                else {                                                                                                                  // no - not using user-provided timesteps
                    dt = ChooseTimestep(dt);
                }

                error = EvolveOneTimestep(dt);                                                                                          // evolve the binary system one timestep

                if (OPTIONS->EmitGravitationalRadiation()) {
                    EmitGravitationalWave(dt);                                                                                          // emit a graviataional wave
                }

                if (error != ERROR::NONE) {                                                                                             // SSE error for either constituent star?
                    evolutionStatus = EVOLUTION_STATUS::SSE_ERROR;                                                                      // yes - stop evolution
                }
                else {                                                                                                                  // continue evolution

                    (void)PrintDetailedOutput(m_Id, BSE_DETAILED_RECORD_TYPE::POST_STELLAR_TIMESTEP);                                   // print (log) detailed output

                    if (OPTIONS->RLOFPrinting()) StashRLOFProperties(MT_TIMING::PRE_MT);                                                // stash properties immediately pre-Mass Transfer 

                    EvaluateBinary(dt);                                                                                                 // evaluate the binary at this timestep

                    (void)PrintDetailedOutput(m_Id, BSE_DETAILED_RECORD_TYPE::POST_BINARY_TIMESTEP);                                    // print (log) detailed output
                
                    (void)PrintRLOFParameters();                                                                                        // print (log) RLOF parameters

                    // check for reasons to not continue evolution
                    if (StellarMerger() && !HasOneOf({ STELLAR_TYPE::MASSLESS_REMNANT })) {                                             // have stars merged without merger already being resolved?
                        if (m_Star1->IsOneOf(MAIN_SEQUENCE) && m_Star2->IsOneOf(MAIN_SEQUENCE) && OPTIONS->EvolveMainSequenceMergers()) // yes - both MS and evolving MS merger products?
                            ResolveMainSequenceMerger();                                                                                // yes - handle main sequence mergers gracefully; no need to change evolution status
                        else
                            evolutionStatus = EVOLUTION_STATUS::STELLAR_MERGER;                                                         // no - for now, stop evolution
                    }
                    else if (HasStarsTouching()) {                                                                                      // binary components touching? (should usually be avoided as MT or CE or merger should happen prior to this)
                        evolutionStatus = EVOLUTION_STATUS::STARS_TOUCHING;                                                             // yes - stop evolution
                    }
                    else if (IsUnbound()) {                                                                                             // binary is unbound?
                        m_Flags.mergesInHubbleTime = false;                                                                             // yes - won't merge in a Hubble time

                        if (IsDCO()) {                                                                                                  // DCO (has two COs)?
                            if (m_DCOFormationTime == DEFAULT_INITIAL_DOUBLE_VALUE) {                                                   // DCO not yet evaluated
                                m_DCOFormationTime = m_Time;                                                                            // set the DCO formation time
                            }
                        }

                        if (!OPTIONS->EvolveUnboundSystems() || IsDCO()) {                                                              // should we evolve unbound systems?
                            evolutionStatus = EVOLUTION_STATUS::UNBOUND;                                                                // no - stop evolution
                        }
                    }

                    if (evolutionStatus == EVOLUTION_STATUS::CONTINUE) {                                                                // continue evolution?

                        if (HasOneOf({ STELLAR_TYPE::NEUTRON_STAR })) {
                            (void)PrintPulsarEvolutionParameters(PULSAR_RECORD_TYPE::POST_BINARY_TIMESTEP);                             // print (log) pulsar evolution parameters 
                        }

                        //(void)PrintBeBinary();                                                                                          // print (log) BeBinary properties
                        
                        if (IsDCO() && !IsUnbound()) {                                                                                  // bound double compact object?
                            if (m_DCOFormationTime == DEFAULT_INITIAL_DOUBLE_VALUE) {                                                   // DCO not yet evaluated -- to ensure that the coalescence is only resolved once
                                ResolveCoalescence();                                                                                   // yes - resolve coalescence
                                m_DCOFormationTime = m_Time;                                                                            // set the DCO formation time
                            }

                            if (!(OPTIONS->EvolvePulsars() && HasOneOf({ STELLAR_TYPE::NEUTRON_STAR }))) {                              // evolve pulsar?
                                evolutionStatus = EVOLUTION_STATUS::DCO;                                                                // no - have DCO - stop evolving
                            }
                        }

                        // check whether to continue evolution
                        if (evolutionStatus == EVOLUTION_STATUS::CONTINUE) {                                                            // continue evolution?

                            // check for other reasons to stop evolution
                            if (IsDCO() && m_Time > (m_DCOFormationTime + m_TimeToCoalescence) && !IsUnbound()) {                       // evolution time exceeds DCO merger time?
                                evolutionStatus = EVOLUTION_STATUS::DCO_MERGER_TIME;                                                    // yes - stop evolution
                            }
                            else if (m_Time > OPTIONS->MaxEvolutionTime()) {                                                            // evolution time exceeds maximum?
                                evolutionStatus = EVOLUTION_STATUS::TIMES_UP;                                                           // yes - stop evolution
                            }
                            else if (!OPTIONS->EvolveDoubleWhiteDwarfs() && IsWDandWD()) {                                              // double WD and their evolution is not enabled?
                                evolutionStatus = EVOLUTION_STATUS::WD_WD;                                                              // yes - do not evolve double WD systems
                            }
                            else if ((HasOneOf({ STELLAR_TYPE::MASSLESS_REMNANT }) && !OPTIONS->EvolveMainSequenceMergers()) || IsMRandRemant()) {           // at least one massless remnant and not evolving MS merger products, or is MR + stellar remnant
                                evolutionStatus = EVOLUTION_STATUS::MASSLESS_REMNANT;                                                   // yes - stop evolution
                            }
                        }
                    }
                }

                (void)PrintDetailedOutput(m_Id, BSE_DETAILED_RECORD_TYPE::TIMESTEP_COMPLETED);                                          // print (log) detailed output: this is after all changes made in the timestep

                if (stepNum >= OPTIONS->MaxNumberOfTimestepIterations()) evolutionStatus = EVOLUTION_STATUS::STEPS_UP;                  // number of timesteps for evolution exceeds maximum
                else if (evolutionStatus == EVOLUTION_STATUS::CONTINUE && usingProvidedTimesteps && stepNum >= timesteps.size()) {      // using user-provided timesteps and all consumed
                    evolutionStatus = EVOLUTION_STATUS::TIMESTEPS_EXHAUSTED;                                                            // yes - set status
                    SHOW_WARN(ERROR::TIMESTEPS_EXHAUSTED);                                                                              // show warning
                }

            }

            if (usingProvidedTimesteps && timesteps.size() > stepNum) {                                                                 // all user-defined timesteps consumed?
                evolutionStatus = EVOLUTION_STATUS::TIMESTEPS_NOT_CONSUMED;                                                             // no - set status
                SHOW_WARN(ERROR::TIMESTEPS_NOT_CONSUMED);                                                                               // show warning
            }
        }

        (void)PrintDetailedOutput(m_Id, BSE_DETAILED_RECORD_TYPE::FINAL_STATE);                                                         // print (log) detailed output: this is the final state of the binary

        // if we trapped a floating-point error we set the binary's error value to indicate a
        // floating-point error occurred, but we don't terminate evolution (we can only have
        // floating-point errors trapped here if the user has not activated the floating-point
        // error instrumentation.  i.e --fp-error-mode OFF)
        // Set the error here so that users know that a floating-point error occurred, even though
        // the evolution of the binary was not terminated because an error occurred.

        if (fetestexcept(FE_DIVBYZERO) ||
            fetestexcept(FE_INVALID)   ||
            fetestexcept(FE_OVERFLOW)  ||
            fetestexcept(FE_UNDERFLOW)) m_Error = ERROR::FLOATING_POINT_ERROR;                                                     // floating-point error

            feclearexcept(FE_ALL_EXCEPT);                                                                                          // clear all FE traps
            
    }
    catch (const std::runtime_error& e) {                                                                                               // catch runtime exceptions
        // anything we catch here should not already have been displayed to the user,
        // so set the error value, display the error, and flag termination (do not rethrow the error)
        if (std::string(e.what()) == "FPE") m_Error = ERROR::FLOATING_POINT_ERROR;                                                      // floating-point error
        else                                m_Error = ERROR::ERROR;                                                                     // unspecified error
        SHOW_ERROR(m_Error);                                                                                                            // display error (don't throw here - handled by returning status)
        evolutionStatus = EVOLUTION_STATUS::ERROR;                                                                                      // evolution terminated
    }
    catch (int e) {
        // anything we catch here should already have been displayed to the user,
        // so just ensure error value is set and flag termination (do not rethrow the error)
        if (e != static_cast<int>(ERROR::NONE)) m_Error = static_cast<ERROR>(e);                                                        // specified errpr
        else                                    m_Error = ERROR::ERROR;                                                                 // unspecified error
        evolutionStatus = EVOLUTION_STATUS::ERROR;                                                                                      // evolution terminated
    }
    catch (...) {
        // anything we catch here should not already have been displayed to the user,
        // so set the error value, display the error, and flag termination (do not rethrow the error)
        m_Error = ERROR::ERROR;                                                                                                         // unspecified error
        SHOW_ERROR(m_Error);                                                                                                            // display error (don't throw here - handled by returning status)
        evolutionStatus = EVOLUTION_STATUS::ERROR;                                                                                      // evolution terminated
    }

    m_EvolutionStatus = evolutionStatus;

    (void)PrintBinarySystemParameters();                                                                                                // print (log) binary system parameters

    return evolutionStatus;
}

