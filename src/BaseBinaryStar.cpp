#include "BaseBinaryStar.h"
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
    bool sampled = OPTIONS->OptionSpecified("initial-mass-1") == 0 ||
                   OPTIONS->OptionSpecified("initial-mass-2") == 0 ||
                  (OPTIONS->OptionSpecified("metallicity") == 0 && OPTIONS->MetallicityDistribution() != METALLICITY_DISTRIBUTION::ZSOLAR) ||
                  (OPTIONS->OptionSpecified("semi-major-axis") == 0 && OPTIONS->OptionSpecified("orbital-period") == 0) ||
                  (OPTIONS->OptionSpecified("eccentricity") == 0 && OPTIONS->EccentricityDistribution() != ECCENTRICITY_DISTRIBUTION::ZERO);


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
    kickParameters1.magnitudeRandomSpecified = OPTIONS->OptionSpecified("kick-magnitude-random-1") == 1;
    kickParameters1.magnitudeRandom          = OPTIONS->KickMagnitudeRandom1();
    kickParameters1.magnitudeSpecified       = OPTIONS->OptionSpecified("kick-magnitude-1") == 1;
    kickParameters1.magnitude                = OPTIONS->KickMagnitude1();
    kickParameters1.phiSpecified             = OPTIONS->OptionSpecified("kick-phi-1") == 1;
    kickParameters1.phi                      = OPTIONS->SN_Phi1();
    kickParameters1.thetaSpecified           = OPTIONS->OptionSpecified("kick-theta-1") == 1;
    kickParameters1.theta                    = OPTIONS->SN_Theta1();
    kickParameters1.meanAnomalySpecified     = OPTIONS->OptionSpecified("kick-mean-anomaly-1") == 1;
    kickParameters1.meanAnomaly              = OPTIONS->SN_MeanAnomaly1();

    KickParameters kickParameters2;
    kickParameters2.magnitudeRandomSpecified = OPTIONS->OptionSpecified("kick-magnitude-random-2") == 1;
    kickParameters2.magnitudeRandom          = OPTIONS->KickMagnitudeRandom2();
    kickParameters2.magnitudeSpecified       = OPTIONS->OptionSpecified("kick-magnitude-2") == 1;
    kickParameters2.magnitude                = OPTIONS->KickMagnitude2();
    kickParameters2.phiSpecified             = OPTIONS->OptionSpecified("kick-phi-2") == 1;
    kickParameters2.phi                      = OPTIONS->SN_Phi2();
    kickParameters2.thetaSpecified           = OPTIONS->OptionSpecified("kick-theta-2") == 1;
    kickParameters2.theta                    = OPTIONS->SN_Theta2();
    kickParameters2.meanAnomalySpecified     = OPTIONS->OptionSpecified("kick-mean-anomaly-2") == 1;
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

        double mass1 = OPTIONS->OptionSpecified("initial-mass-1") == 1                                                                  // user specified primary mass?
                        ? OPTIONS->InitialMass1()                                                                                       // yes, use it
                        : utils::SampleInitialMass(OPTIONS->InitialMassFunction(), 
                                                   OPTIONS->InitialMassFunctionMax(), 
                                                   OPTIONS->InitialMassFunctionMin(), 
                                                   OPTIONS->InitialMassFunctionPower());                                                // no - asmple it

        double mass2 = 0.0;                      
        if (OPTIONS->OptionSpecified("initial-mass-2") == 1) {                                                                          // user specified secondary mass?
            mass2 = OPTIONS->InitialMass2();                                                                                            // yes, use it
        }
        else {                                                                                                                          // no - sample it
            // first, determine mass ratio q    
            double q = OPTIONS->OptionSpecified("mass-ratio") == 1                                                                      // user specified mass ratio?
                        ? OPTIONS->MassRatio()                                                                                          // yes, use it
                        : utils::SampleMassRatio(OPTIONS->MassRatioDistribution(),
                                                 OPTIONS->MassRatioDistributionMax(), 
                                                 OPTIONS->MassRatioDistributionMin());                                                  // no - sample it

            mass2 = mass1 * q;                                                                                                          // calculate mass2 using mass ratio                                                                     
        }

        double metallicity = OPTIONS->OptionSpecified("metallicity") == 1                                                               // user specified metallicity?
                                ? OPTIONS->Metallicity()                                                                                // yes, use it
                                : utils::SampleMetallicity(OPTIONS->MetallicityDistribution(), 
                                                           OPTIONS->MetallicityDistributionMax(), 
                                                           OPTIONS->MetallicityDistributionMin());                                      // no, sample it

        if (OPTIONS->OptionSpecified("semi-major-axis") == 1) {                                                                         // user specified semi-major axis?
            m_SemiMajorAxis = OPTIONS->SemiMajorAxis();                                                                                 // yes, use it
        }
        else {                                                                                                                          // no, semi-major axis not specified
            if (OPTIONS->OptionSpecified("orbital-period") == 1) {                                                                      // user specified orbital period?
                m_SemiMajorAxis = utils::ConvertPeriodInDaysToSemiMajorAxisInAU(mass1, mass2, OPTIONS->OrbitalPeriod());                // yes - calculate semi-major axis from period
            }
            else {                                                                                                                      // no
                if (OPTIONS->OptionSpecified("semi-major-axis-distribution") == 1 ||                                                    // user specified semi-major axis distribution, or
                    OPTIONS->OptionSpecified("orbital-period-distribution" ) == 0) {                                                    // user did not specify oprbital period distribution
                    m_SemiMajorAxis = utils::SampleSemiMajorAxis(OPTIONS->SemiMajorAxisDistribution(),                              
                                                                 OPTIONS->SemiMajorAxisDistributionMax(), 
                                                                 OPTIONS->SemiMajorAxisDistributionMin(),
                                                                 OPTIONS->SemiMajorAxisDistributionPower(), 
                                                                 OPTIONS->OrbitalPeriodDistributionMax(), 
                                                                 OPTIONS->OrbitalPeriodDistributionMin(), 
                                                                 mass1, 
                                                                 mass2);                                                                // yes, sample from semi-major axis distribution (might be default)
                }
                else {                                                                                                                  // no - sample from orbital period distribution
                    double orbitalPeriod = utils::SampleOrbitalPeriod(OPTIONS->OrbitalPeriodDistribution(),                              
                                                                      OPTIONS->OrbitalPeriodDistributionMax(), 
                                                                      OPTIONS->OrbitalPeriodDistributionMin());

                    m_SemiMajorAxis = utils::ConvertPeriodInDaysToSemiMajorAxisInAU(mass1, mass2, orbitalPeriod);                       // calculate semi-major axis from period
                }
            }
        }

        m_Eccentricity = OPTIONS->OptionSpecified("eccentricity") == 1                                                                  // user specified eccentricity?
                            ? OPTIONS->Eccentricity()                                                                                   // yes, use it
                            : utils::SampleEccentricity(OPTIONS->EccentricityDistribution(), 
                                                        OPTIONS->EccentricityDistributionMax(), 
                                                        OPTIONS->EccentricityDistributionMin());                                        // no, sample it

        // binary star contains two instances of star to hold masses, radii and luminosities.
        // star 1 initially more massive
        m_Star1 = OPTIONS->OptionSpecified("rotational-frequency-1") == 1                                                               // user specified primary rotational frequency?
                    ? new BinaryConstituentStar(m_RandomSeed, mass1, metallicity, kickParameters1, OPTIONS->RotationalFrequency1() * SECONDS_IN_YEAR) // yes - use it (convert from Hz to cycles per year - see BaseStar::CalculateZAMSAngularFrequency())
                    : new BinaryConstituentStar(m_RandomSeed, mass1, metallicity, kickParameters1);                                     // no - let it be calculated

        m_Star2 = OPTIONS->OptionSpecified("rotational-frequency-2") == 1                                                               // user specified secondary rotational frequency?
                    ? new BinaryConstituentStar(m_RandomSeed, mass2, metallicity, kickParameters2, OPTIONS->RotationalFrequency2() * SECONDS_IN_YEAR) // yes - use it (convert from Hz to cycles per year - see BaseStar::CalculateZAMSAngularFrequency())
                    : new BinaryConstituentStar(m_RandomSeed, mass2, metallicity, kickParameters2);                                     // no - let it be calculated

        double starToRocheLobeRadiusRatio1 = (m_Star1->Radius() * RSOL_TO_AU) / (m_SemiMajorAxis * (1.0 - m_Eccentricity) * CalculateRocheLobeRadius_Static(mass1, mass2));
        double starToRocheLobeRadiusRatio2 = (m_Star2->Radius() * RSOL_TO_AU) / (m_SemiMajorAxis * (1.0 - m_Eccentricity) * CalculateRocheLobeRadius_Static(mass2, mass1));

        m_Flags.massesEquilibrated        = false;                                                                                      // default
        m_Flags.massesEquilibratedAtBirth = false;                                                                                      // default

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
            m_Star1 = OPTIONS->OptionSpecified("rotational-frequency-1") == 1                                                           // user specified primary rotational frequency?
                        ? new BinaryConstituentStar(m_RandomSeed, mass1, metallicity, kickParameters1, OPTIONS->RotationalFrequency1() * SECONDS_IN_YEAR) // yes - use it (convert from Hz to cycles per year - see BaseStar::CalculateZAMSAngularFrequency())
                        : new BinaryConstituentStar(m_RandomSeed, mass1, metallicity, kickParameters1);                                 // no - let it be calculated

            delete m_Star2;
            m_Star2 = OPTIONS->OptionSpecified("rotational-frequency-2") == 1                                                           // user specified secondary rotational frequency?
                        ? new BinaryConstituentStar(m_RandomSeed, mass2, metallicity, kickParameters2, OPTIONS->RotationalFrequency2() * SECONDS_IN_YEAR) // yes - use it (convert from Hz to cycles per year - see BaseStar::CalculateZAMSAngularFrequency())
                        : new BinaryConstituentStar(m_RandomSeed, mass2, metallicity, kickParameters2);                                 // no - let it be calculated
        
            starToRocheLobeRadiusRatio1 = (m_Star1->Radius() * RSOL_TO_AU) / (m_SemiMajorAxis * CalculateRocheLobeRadius_Static(mass1, mass2));   //eccentricity already zero
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
            m_Error = ERROR::INVALID_INITIAL_ATTRIBUTES;
            done = true;
        }

    } while (!done && ++tries < MAX_BSE_INITIAL_CONDITIONS_ITERATIONS);

    if (!done) m_Error = ERROR::INVALID_INITIAL_ATTRIBUTES;                                                                             // too many iterations - bad initial conditions

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

    m_Error = ERROR::NONE;

    m_ObjectId    = globalObjectId++;
    m_ObjectType  = OBJECT_TYPE::BASE_BINARY_STAR;
    m_StellarType = STELLAR_TYPE::BINARY_STAR;
    m_RandomSeed  = p_Seed;
    m_Id          = p_Id;

    if (OPTIONS->PopulationDataPrinting()) {                                                            // user wants to see details of binary?
        SAY("Using supplied random seed " << m_RandomSeed << " for Binary Star id = " << m_ObjectId);   // yes - show them
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

    // if CHE enabled, update rotational frequency for constituent stars - assume tidally locked

    if (OPTIONS->CHEMode() != CHE_MODE::NONE) {

        m_Star1->SetOmega(OrbitalAngularVelocity());
        m_Star2->SetOmega(OrbitalAngularVelocity());

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
        else if (m_Star1->MZAMS() <= 0.7) {                                                                                                             // no - MS - initial mass determines actual type  JR: don't use utils::Compare() here
            if (m_Star1->StellarType() != STELLAR_TYPE::MS_LTE_07) (void)m_Star1->SwitchTo(STELLAR_TYPE::MS_LTE_07, true);                              // MS <= 0.7 Msol - switch if necessary
        }
        else {
            if (m_Star1->StellarType() != STELLAR_TYPE::MS_GT_07) (void)m_Star1->SwitchTo(STELLAR_TYPE::MS_GT_07, true);                                // MS > 0.7 Msol - switch if necessary
        }

        // star 2
        if (utils::Compare(m_Star1->Omega(), m_Star2->OmegaCHE()) >= 0) {                                                                               // star 2 CH?
            if (m_Star2->StellarType() != STELLAR_TYPE::CHEMICALLY_HOMOGENEOUS) (void)m_Star2->SwitchTo(STELLAR_TYPE::CHEMICALLY_HOMOGENEOUS, true);    // yes, switch if not already Chemically Homogeneous
        }
        else if (m_Star2->MZAMS() <= 0.7) {                                                                                                             // no - MS - initial mass determines actual type  JR: don't use utils::Compare() here
            if (m_Star2->StellarType() != STELLAR_TYPE::MS_LTE_07) (void)m_Star2->SwitchTo(STELLAR_TYPE::MS_LTE_07, true);                              // MS <= 0.0 Msol - switch if necessary
        }
        else {
            if (m_Star2->StellarType() != STELLAR_TYPE::MS_GT_07) (void)m_Star2->SwitchTo(STELLAR_TYPE::MS_GT_07, true);                                // MS > 0.7 Msol - switch if necessary
        }
    }

    double gyrationRadius1                       = m_Star1->CalculateGyrationRadius();
    double gyrationRadius2                       = m_Star2->CalculateGyrationRadius();

    m_TotalEnergy                                = CalculateTotalEnergy(m_SemiMajorAxis,
                                                                        m_Star1->Mass(),
                                                                        m_Star2->Mass(),
                                                                        m_Star1->RZAMS(),
                                                                        m_Star2->RZAMS(),
                                                                        m_Star1->Omega(),
                                                                        m_Star2->Omega(),
                                                                        gyrationRadius1,
                                                                        gyrationRadius2);

    m_TotalAngularMomentum                      = CalculateAngularMomentum(m_SemiMajorAxis,
                                                                            m_Eccentricity,
                                                                            m_Star1->Mass(),
                                                                            m_Star2->Mass(),
                                                                            m_Star1->RZAMS(),
                                                                            m_Star2->RZAMS(),
                                                                            m_Star1->Omega(),
                                                                            m_Star2->Omega(),
                                                                            gyrationRadius1,
                                                                            gyrationRadius2);

	double totalMass 					         = m_Star1->Mass() + m_Star2->Mass();
	double reducedMass					         = (m_Star1->Mass() * m_Star2->Mass()) / totalMass;
	m_OrbitalEnergy 			                 = CalculateOrbitalEnergy(reducedMass, totalMass, m_SemiMajorAxis);
	m_OrbitalEnergyPrev 			             = m_OrbitalEnergy;

	m_OrbitalAngularMomentum 	                 = CalculateOrbitalAngularMomentum(reducedMass, totalMass, m_SemiMajorAxis);
	m_OrbitalAngularMomentumPrev 	             = m_OrbitalAngularMomentum;

    m_Time                                       = DEFAULT_INITIAL_DOUBLE_VALUE;
	m_Dt                                         = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_TimePrev                                   = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_DCOFormationTime                           = DEFAULT_INITIAL_DOUBLE_VALUE;

    m_aMassLossDiff                              = DEFAULT_INITIAL_DOUBLE_VALUE;

    m_aMassTransferDiff                          = DEFAULT_INITIAL_DOUBLE_VALUE;

	m_MassTransferTrackerHistory                 = MT_TRACKING::NO_MASS_TRANSFER;
    m_MassTransfer                               = false;

    m_JLoss                                      = OPTIONS->MassTransferJloss();

	m_FractionAccreted                           = OPTIONS->MassTransferFractionAccreted();

    // Common Envelope
    m_CEDetails.CEEcount                         = 0;
    m_CEDetails.CEEnow                           = false;
    m_CEDetails.doubleCoreCE                     = false;
	m_CEDetails.optimisticCE                     = false;
	m_CEDetails.postCEE.eccentricity             = DEFAULT_INITIAL_DOUBLE_VALUE;
	m_CEDetails.postCEE.rocheLobe1to2            = DEFAULT_INITIAL_DOUBLE_VALUE;
	m_CEDetails.postCEE.rocheLobe2to1            = DEFAULT_INITIAL_DOUBLE_VALUE;
	m_CEDetails.postCEE.semiMajorAxis            = DEFAULT_INITIAL_DOUBLE_VALUE;
	m_CEDetails.preCEE.eccentricity              = DEFAULT_INITIAL_DOUBLE_VALUE;
	m_CEDetails.preCEE.rocheLobe1to2             = DEFAULT_INITIAL_DOUBLE_VALUE;
	m_CEDetails.preCEE.rocheLobe2to1             = DEFAULT_INITIAL_DOUBLE_VALUE;
	m_CEDetails.preCEE.semiMajorAxis             = DEFAULT_INITIAL_DOUBLE_VALUE;

    m_Flags.stellarMerger                        = false;
    m_Flags.stellarMergerAtBirth                 = false;

	m_Mass1Final                                 = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_Mass2Final                                 = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_MassEnv1                                   = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_MassEnv2                                   = DEFAULT_INITIAL_DOUBLE_VALUE;

    m_ZetaLobe                                   = DEFAULT_INITIAL_DOUBLE_VALUE;
	m_ZetaStar	                                 = DEFAULT_INITIAL_DOUBLE_VALUE;

    // Initialise other parameters to 0
    m_uK                                         = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_CosIPrime                                  = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_IPrime                                     = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_TimeToCoalescence                          = DEFAULT_INITIAL_DOUBLE_VALUE;

    m_SupernovaState                             = SN_STATE::NONE;

    m_Flags.mergesInHubbleTime                   = false;
    m_Unbound                                    = false;

    m_SystemicVelocity                           = Vector3d();
    m_OrbitalAngularMomentumVector               = Vector3d();
	m_ThetaE                                     = DEFAULT_INITIAL_DOUBLE_VALUE;
	m_PhiE                                       = DEFAULT_INITIAL_DOUBLE_VALUE;
	m_PsiE                                       = DEFAULT_INITIAL_DOUBLE_VALUE;

	m_SynchronizationTimescale                   = DEFAULT_INITIAL_DOUBLE_VALUE;
	m_CircularizationTimescale                   = DEFAULT_INITIAL_DOUBLE_VALUE;

	// RLOF details
    m_RLOFDetails.experiencedRLOF                          = false;
    m_RLOFDetails.immediateRLOFPostCEE                     = false;
    m_RLOFDetails.isRLOF                                   = false;
    m_RLOFDetails.simultaneousRLOF                         = false;
    m_RLOFDetails.stableRLOFPostCEE                        = false;

	// RLOF details - properties 1
    m_RLOFDetails.props1.id                                = -1l;

    m_RLOFDetails.props1.stellarType1                      = STELLAR_TYPE::NONE;
    m_RLOFDetails.props1.stellarType2                      = STELLAR_TYPE::NONE;

    m_RLOFDetails.props1.mass1                             = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_RLOFDetails.props1.mass2                             = DEFAULT_INITIAL_DOUBLE_VALUE;

    m_RLOFDetails.props1.radius1                           = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_RLOFDetails.props1.radius2                           = DEFAULT_INITIAL_DOUBLE_VALUE;

    m_RLOFDetails.props1.starToRocheLobeRadiusRatio1       = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_RLOFDetails.props1.starToRocheLobeRadiusRatio2       = DEFAULT_INITIAL_DOUBLE_VALUE;

    m_RLOFDetails.props1.semiMajorAxis                     = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_RLOFDetails.props1.eccentricity                      = DEFAULT_INITIAL_DOUBLE_VALUE;
    
    m_RLOFDetails.props1.eventCounter                      = DEFAULT_INITIAL_ULONGINT_VALUE;

    m_RLOFDetails.props1.time                              = DEFAULT_INITIAL_DOUBLE_VALUE;

    m_RLOFDetails.props1.isRLOF1                           = false;
    m_RLOFDetails.props1.isRLOF2                           = false;

    m_RLOFDetails.props1.isCE                              = false;

	// RLOF details - properties 2
    m_RLOFDetails.props2.id = -1l;

    m_RLOFDetails.props2.stellarType1                      = STELLAR_TYPE::NONE;
    m_RLOFDetails.props2.stellarType2                      = STELLAR_TYPE::NONE;

    m_RLOFDetails.props2.mass1                             = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_RLOFDetails.props2.mass2                             = DEFAULT_INITIAL_DOUBLE_VALUE;

    m_RLOFDetails.props2.radius1                           = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_RLOFDetails.props2.radius2                           = DEFAULT_INITIAL_DOUBLE_VALUE;

    m_RLOFDetails.props2.starToRocheLobeRadiusRatio1       = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_RLOFDetails.props2.starToRocheLobeRadiusRatio2       = DEFAULT_INITIAL_DOUBLE_VALUE;


    m_RLOFDetails.props2.semiMajorAxis                     = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_RLOFDetails.props2.eccentricity                      = DEFAULT_INITIAL_DOUBLE_VALUE;
    
    m_RLOFDetails.props2.eventCounter                      = DEFAULT_INITIAL_ULONGINT_VALUE;

    m_RLOFDetails.props2.time                              = DEFAULT_INITIAL_DOUBLE_VALUE;

    m_RLOFDetails.props2.isRLOF1                           = false;
    m_RLOFDetails.props2.isRLOF2                           = false;

    m_RLOFDetails.props2.isCE                              = false;

    // RLOF details - pre/post-MT props pointers
    m_RLOFDetails.propsPostMT                              = &m_RLOFDetails.props1;
    m_RLOFDetails.propsPreMT                               = &m_RLOFDetails.props2;


    // BeBinary details - properties 1
    m_BeBinaryDetails.props1.id                  = -1l;

    m_BeBinaryDetails.props1.dt                  = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_BeBinaryDetails.props1.totalTime           = DEFAULT_INITIAL_DOUBLE_VALUE;

    m_BeBinaryDetails.props1.massNS              = DEFAULT_INITIAL_DOUBLE_VALUE;

    m_BeBinaryDetails.props1.companionMass       = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_BeBinaryDetails.props1.companionLuminosity = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_BeBinaryDetails.props1.companionTeff       = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_BeBinaryDetails.props1.companionRadius     = DEFAULT_INITIAL_DOUBLE_VALUE;

    m_BeBinaryDetails.props1.semiMajorAxis       = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_BeBinaryDetails.props1.eccentricity        = DEFAULT_INITIAL_DOUBLE_VALUE;

    // BeBinary details - properties 2
    m_BeBinaryDetails.props2.id                  = -1l;

    m_BeBinaryDetails.props2.dt                  = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_BeBinaryDetails.props2.totalTime           = DEFAULT_INITIAL_DOUBLE_VALUE;

    m_BeBinaryDetails.props2.massNS              = DEFAULT_INITIAL_DOUBLE_VALUE;

    m_BeBinaryDetails.props2.companionMass       = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_BeBinaryDetails.props2.companionLuminosity = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_BeBinaryDetails.props2.companionTeff       = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_BeBinaryDetails.props2.companionRadius     = DEFAULT_INITIAL_DOUBLE_VALUE;

    m_BeBinaryDetails.props2.semiMajorAxis       = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_BeBinaryDetails.props2.eccentricity        = DEFAULT_INITIAL_DOUBLE_VALUE;

    // BeBinary details - current/prev props pointers
    m_BeBinaryDetails.currentProps               = &m_BeBinaryDetails.props1;
    m_BeBinaryDetails.previousProps              = &m_BeBinaryDetails.props2;

    // pointers

    m_Donor                                      = nullptr;
    m_Accretor                                   = nullptr;

    m_Supernova                                  = nullptr;
    m_Companion                                  = nullptr;
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
 * The functional return is the value of the property requested.  The type of the
 * functional return is a tuple: std::tuple<bool, COMPAS_VARIABLE_TYPE>.  This type
 * is COMPAS_VARIABLE by typedef.
 *
 * The bool returned indicates whether the property value was retrieved ok: true = yes, fales = no
 * The COMPAS_VARIABLE_TYPE variable returned is a boost variant variable, the value of which is the
 * value of the underlying primitive variable.
 *
 *
 * COMPAS_VARIABLE BinaryPropertyValue(const T_ANY_PROPERTY p_Property) const
 *
 * @param   [IN]    p_Property                  The property for which the value is required
 * @return                                      The value of the requested property
 */
COMPAS_VARIABLE BaseBinaryStar::BinaryPropertyValue(const T_ANY_PROPERTY p_Property) const {

    bool ok = true;                                                                                                     // default is no error

    COMPAS_VARIABLE_TYPE value;                                                                                         // property value

    BINARY_PROPERTY property = boost::get<BINARY_PROPERTY>(p_Property);                                                 // get the id of the property required

    switch (property) {                                                                                                 // which property?

        case BINARY_PROPERTY::BE_BINARY_CURRENT_COMPANION_LUMINOSITY:               value = BeBinaryDetails().currentProps->companionLuminosity;                break;
        case BINARY_PROPERTY::BE_BINARY_CURRENT_COMPANION_MASS:                     value = BeBinaryDetails().currentProps->companionMass;                      break;
        case BINARY_PROPERTY::BE_BINARY_CURRENT_COMPANION_RADIUS:                   value = BeBinaryDetails().currentProps->companionRadius;                    break;
        case BINARY_PROPERTY::BE_BINARY_CURRENT_COMPANION_TEFF:                     value = BeBinaryDetails().currentProps->companionTeff * TSOL;               break;
        case BINARY_PROPERTY::BE_BINARY_CURRENT_DT:                                 value = BeBinaryDetails().currentProps->dt;                                 break;
        case BINARY_PROPERTY::BE_BINARY_CURRENT_ECCENTRICITY:                       value = BeBinaryDetails().currentProps->eccentricity;                       break;
        case BINARY_PROPERTY::BE_BINARY_CURRENT_ID:                                 value = BeBinaryDetails().currentProps->id;                                 break;
        case BINARY_PROPERTY::BE_BINARY_CURRENT_NS_MASS:                            value = BeBinaryDetails().currentProps->massNS;                             break;
        case BINARY_PROPERTY::BE_BINARY_CURRENT_SEMI_MAJOR_AXIS:                    value = BeBinaryDetails().currentProps->semiMajorAxis;                      break;
        case BINARY_PROPERTY::BE_BINARY_CURRENT_TOTAL_TIME:                         value = BeBinaryDetails().currentProps->totalTime;                          break;
        case BINARY_PROPERTY::CIRCULARIZATION_TIMESCALE:                            value = CircularizationTimescale();                                         break;
        case BINARY_PROPERTY::COMMON_ENVELOPE_AT_LEAST_ONCE:                        value = CEAtLeastOnce();                                                    break;
        case BINARY_PROPERTY::COMMON_ENVELOPE_EVENT_COUNT:                          value = CommonEnvelopeEventCount();                                         break;
        case BINARY_PROPERTY::DIMENSIONLESS_KICK_MAGNITUDE:                         value = UK();                                                               break;
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

        default:                                                                                                        // unknown property
            ok    = false;                                                                                              // that's not ok...
            value = "UNKNOWN";                                                                                          // default value
            SHOW_WARN(ERROR::UNKNOWN_BINARY_PROPERTY);                                                                  // show warning
    }

    return std::make_tuple(ok, value);
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
 * This function calls the appropriate helper function to retrieve the value.
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
 * COMPAS_VARIABLE PropertyValue(const T_ANY_PROPERTY p_Property) const
 *
 * @param   [IN]    p_Property                  The property for which the value is required
 * @return                                      The value of the requested property
 */
COMPAS_VARIABLE BaseBinaryStar::PropertyValue(const T_ANY_PROPERTY p_Property) const {

    bool ok = false;                                                                                                    // default is failure

    COMPAS_VARIABLE_TYPE value;                                                                                         // property value

    switch (boost::apply_visitor(VariantPropertyType(), p_Property)) {                                                  // which property type?

        case ANY_PROPERTY_TYPE::T_BINARY_PROPERTY:                                                                      // BSE binary star property
            std::tie(ok, value) = BinaryPropertyValue(p_Property);                                                      // get the value
            break;

        case ANY_PROPERTY_TYPE::T_STAR_1_PROPERTY:                                                                      // star 1 of BSE binary star property
            if (m_Star1) std::tie(ok, value) = m_Star1->StellarPropertyValue(p_Property);                               // if have pointer to primary, get the value
            break;

        case ANY_PROPERTY_TYPE::T_STAR_2_PROPERTY:                                                                      // star 2 of BSE binary star property
            if (m_Star2) std::tie(ok, value) = m_Star2->StellarPropertyValue(p_Property);                               // if have pointer to secondary, get the value
            break;

        case ANY_PROPERTY_TYPE::T_SUPERNOVA_PROPERTY:                                                                   // supernova star of BSE binary star property
            if (m_Supernova) std::tie(ok, value) = m_Supernova->StellarPropertyValue(p_Property);                       // if have pointer to supernova, get the value
            break;

        case ANY_PROPERTY_TYPE::T_COMPANION_PROPERTY:                                                                   // companion star of BSE binary star property
            if (m_Companion) std::tie(ok, value) = m_Companion->StellarPropertyValue(p_Property);                       // if have pointer to companion, get the value
            break;

        case ANY_PROPERTY_TYPE::T_PROGRAM_OPTION:                                                                       // program option
            std::tie(ok, value) = OPTIONS->OptionValue(p_Property);                                                     // get the value
            break;

        default:                                                                                                        // unknown property type
            SHOW_WARN(ERROR::UNKNOWN_PROPERTY_TYPE);                                                                    // show warning
    }

    return std::make_tuple(ok, value);
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

    StashRLOFProperties(MASS_TRANSFER_TIMING::POST_MT);                 // stash properties immediately post-Mass Transfer 

    if (m_Star1->IsRLOF() || m_Star2->IsRLOF()) {                       // print if either star is in RLOF
        m_RLOFDetails.propsPostMT->eventCounter += 1;                   // every time we print a MT event happened, increment counter
        ok = LOGGING->LogRLOFParameters(this, p_RecordType);    // yes - write to log file
    }

    if (OPTIONS->HMXRBinaries()) {
        if (IsHMXRBinary()) {                                           // print if star is HMXRB candidate
            ok = LOGGING->LogRLOFParameters(this, p_RecordType); 
        }
    }

    return ok;
}

/*
 * Write Be binary parameters to logfile if required
 *
 *
 * bool PrintBeBinary(const BE_BINARY_RECORD_TYPE p_RecordType)
 * 
 * @param   [IN]    p_RecordType                Record type to be written
 * @return                                      Boolean status (true = success, false = failure)
 * 
 */
bool BaseBinaryStar::PrintBeBinary(const BE_BINARY_RECORD_TYPE p_RecordType) {
    
    if (!OPTIONS->BeBinaries()) return true;                // do not print if printing option off
    
    StashBeBinaryProperties();                              // stash Be binary properties
    
    return LOGGING->LogBeBinary(this, p_RecordType);        // write to log file
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
 * @param   [IN]    p_StashPostMassTransfer     Boolean - true if post-MT values should be stored (false for pre-MT values)
 */
void BaseBinaryStar::StashRLOFProperties(const MASS_TRANSFER_TIMING p_Which) {

    if (!OPTIONS->RLOFPrinting()) return;                                                                           // nothing to do

    // set whether to update pre-MT or post-MT parameters depending on input argument
    RLOFPropertiesT* rlofPropertiesToReset;
    rlofPropertiesToReset = (p_Which == MASS_TRANSFER_TIMING::PRE_MT) ?
                             m_RLOFDetails.propsPreMT  :
                             m_RLOFDetails.propsPostMT ;

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
    rlofPropertiesToReset->semiMajorAxis               = m_SemiMajorAxis * AU_TO_RSOL;                               // semi-major axis - change units to Rsol
    rlofPropertiesToReset->time                        = m_Time;
    rlofPropertiesToReset->timePrev                    = m_TimePrev;
    rlofPropertiesToReset->isRLOF1                     = m_Star1->IsRLOF();
    rlofPropertiesToReset->isRLOF2                     = m_Star2->IsRLOF();
    rlofPropertiesToReset->isCE                        = m_CEDetails.CEEnow;
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
void BaseBinaryStar::StashBeBinaryProperties() {

    if (!OPTIONS->BeBinaries() || !IsBeBinary()) return;                                                            // nothing to do;

    // switch previous<->current (preserves existing current as (new) previous)
    BeBinaryPropertiesT* tmp;
    tmp                             = m_BeBinaryDetails.previousProps;                                              // save pointer to existing previous props
    m_BeBinaryDetails.previousProps = m_BeBinaryDetails.currentProps;                                               // existing current props become new previous props (values will be preserved)
    m_BeBinaryDetails.currentProps  = tmp;                                                                          // new current props points at existing previous (values will be replaced)

    // now save (new) current
    m_BeBinaryDetails.currentProps->id            = m_ObjectId;                                                      // object id
    m_BeBinaryDetails.currentProps->dt            = m_Dt;                                                            // timestep
    m_BeBinaryDetails.currentProps->totalTime     = m_BeBinaryDetails.previousProps->dt + m_Dt;                      // total time - accumulate, don't just replace
    m_BeBinaryDetails.currentProps->semiMajorAxis = m_SemiMajorAxis * AU_TO_RSOL;                                    // semi-major axis - change units to Rsol
    m_BeBinaryDetails.currentProps->eccentricity  = m_Eccentricity;                                                  // eccentricity

    BinaryConstituentStar* neutronStar   = m_Star1->IsOneOf({ STELLAR_TYPE::NEUTRON_STAR }) ? m_Star1 : m_Star2;    // pointer to neutron star
    BinaryConstituentStar* companionStar = m_Star1->IsOneOf({ STELLAR_TYPE::NEUTRON_STAR }) ? m_Star2 : m_Star1;    // pointer to companion

    m_BeBinaryDetails.currentProps->massNS              = neutronStar->Mass();                                      // neutron star mass
    m_BeBinaryDetails.currentProps->companionMass       = companionStar->Mass();                                    // companion mass
    m_BeBinaryDetails.currentProps->companionLuminosity = companionStar->Luminosity();                              // companion luminosity
    m_BeBinaryDetails.currentProps->companionTeff       = companionStar->Temperature();                             // companion temperature
    m_BeBinaryDetails.currentProps->companionRadius     = companionStar->Radius();                                  // companion radius
}


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
 * Resolves supernova event - one of the stars has gone supernova!
 *
 * Assign a random supernova kick according to the user specified options and then update the orbit and velocities.
 * Vector algebra is directly based on Pfahl, Rappaport, Podsiadlowski 2002, Appendix B:
 * https://arxiv.org/abs/astro-ph/0106141 
 * The change of reference basis angles, ThetaE, PhiE, and PsiE, are the standard Euler angles (see vector3d.h)
 *
 * Note: the systemic speed is only valid for intact binaries, and component speeds are only valid for disrupted binaries.
 * 
 *  /////////////////////////////////
 *  // Logic
 *  // 
 *  // If (Unbound before SN):
 *  //
 *  //         Must be 2nd SN, only need to update starSN component velocity (rotated into previous reference frame).
 *  //
 *  // Else: (Intact before SN)
 *  //
 *  //        Evolve binary according to vector algebra to determine centerofmass velocity, h', e', a', and whether bound or unbound.
 *  //
 *  //        Update binary systemic velocity (even if disrupted, just for consistency) - rotate into previous reference frame if needed.
 *  // 
 *  //        If now unbound:
 *  //
 *  //                Set m_Unbound to True - should be the only place in the code this is done.
 *  //
 *  //                Continue vector algebra to find v1inf and v2inf.
 *  //                Add these values to previous component velocities (rotated if need be) which will be the systemic velocity if this is the 2nd SN. 
 *  //
 *  //                For unbound binary, new Euler Angles should be randomized (see vector3d.cpp).
 *  //
 *  //        If still intact:
 *  //
 *  //                Binary systemic velocity has already been set, so just set the component velocities to the same vector.
 *  //                (this is to make it easier to add just a component velocity later).
 *  //
 *  //                For intact binary, Euler Angles must be calculated according to the vector algebra (see vector3d.h).
 *  //
 *  /////////////////////////////////////////////////////////////////////////////
 *
 *
 * bool ResolveSupernova()
 *
 * @return                                      True if a supernova event occurred, otherwise false
 */
bool BaseBinaryStar::ResolveSupernova() {

    if (!m_Supernova->IsSNevent()) {
        SHOW_WARN(ERROR::RESOLVE_SUPERNOVA_IMPROPERLY_CALLED);
        return false;                                                                                                   // not a supernova event - bail out 
    }

    // Set relevant preSN parameters 
    m_EccentricityPreSN = m_Eccentricity;                                                 
    m_SemiMajorAxisPreSN = m_SemiMajorAxis;                                               

    double totalMassPreSN = m_Supernova->SN_TotalMassAtCOFormation() + m_Companion->Mass();                             // Total Mass preSN
    double reducedMassPreSN = m_Supernova->SN_TotalMassAtCOFormation() * m_Companion->Mass() / totalMassPreSN;          // Reduced Mass preSN
    m_Supernova->SetOrbitalEnergyPreSN(CalculateOrbitalEnergy(reducedMassPreSN, totalMassPreSN, m_SemiMajorAxisPreSN)); // Orbital energy preSN

    // Define the natal kick vector (see above for precise definitions of the angles)
    double theta = m_Supernova->SN_Theta();                                                                             // Angle out of the binary plane
    double phi   = m_Supernova->SN_Phi();                                                                               // Angle in the binary plane
    Vector3d natalKickVector = m_Supernova->SN_KickMagnitude() *Vector3d(cos(theta)*cos(phi), 
                                                                         cos(theta)*sin(phi),
                                                                         sin(theta));
    // Check if the system is already unbound
    if (IsUnbound()) {                                                                                                  // Is system already unbound?

        m_Supernova->UpdateComponentVelocity( natalKickVector.RotateVector(m_ThetaE, m_PhiE, m_PsiE));                  // yes - only need to update the velocity of the star undergoing SN

        // The quantities below are meaningless in this context, so they are set to nan to avoid misuse
        m_OrbitalVelocityPreSN = -nan("");
        m_uK = nan("");                                                                                                 // -- - Dimensionless kick magnitude

    }
    else {                                                                                                              // no - evaluate orbital changes and calculate velocities
        
        //////////////////////////////////////////////////////////////////////////////////////////////////
        // 
        // Evolve SN out of binary
        // 
        //////////////////////////////////////////////////////////////////////////////////////////////////
        
        
        // Functions defined in vector3d.h
        #define cross(x,y)          linalg::cross(x,y)
        #define dot(x,y)            linalg::dot(x,y) 
        #define angleBetween(x,y)   linalg::angleBetween(x,y)
        #define mag                 Magnitude()
        #define hat                 UnitVector()

        // Pre-SN parameters
        double semiMajorAxisPrev_km = m_SemiMajorAxis * AU_TO_KM;                                                       // km  - Semi-Major axis
        double eccentricityPrev = m_Eccentricity;                                                                       // --  - Eccentricity, written with a prev to distinguish from later use
        double sqrt1MinusEccPrevSquared = std::sqrt(1 - eccentricityPrev * eccentricityPrev);                           // useful function of eccentricity

        double m1Prev = m_Supernova->SN_TotalMassAtCOFormation();                                                       // Mo  - SN star pre-SN mass
        double m2Prev = m_Companion->Mass();                                                                            // Mo  - CP star pre-SN mass
        double totalMassPrev = m1Prev + m2Prev;                                                                         // Mo  - Total binary pre-SN mass
        
        // Functions of eccentric anomaly
        m_Supernova->CalculateSNAnomalies(eccentricityPrev);
        double cosEccAnomaly = cos(m_Supernova->SN_EccentricAnomaly());        
        double sinEccAnomaly = sin(m_Supernova->SN_EccentricAnomaly());

        // Derived quantities
        double omega = std::sqrt(G_SN*totalMassPrev / (semiMajorAxisPrev_km * semiMajorAxisPrev_km*semiMajorAxisPrev_km));   // rad/s  - Keplerian orbital frequency

        Vector3d separationVectorPrev = Vector3d( semiMajorAxisPrev_km * (cosEccAnomaly - eccentricityPrev),            
                                                  semiMajorAxisPrev_km * (sinEccAnomaly) * sqrt1MinusEccPrevSquared,
                                                  0.0                    );                                             // km        - Relative position vector, from m1Prev to m2Prev
        double   separationPrev = separationVectorPrev.mag;                                                             // km        - Instantaneous Separation

        Vector3d relativeVelocityVectorPrev = Vector3d(-((semiMajorAxisPrev_km * semiMajorAxisPrev_km) * omega / separationPrev) * sinEccAnomaly,   
                                                        ((semiMajorAxisPrev_km * semiMajorAxisPrev_km) * omega / separationPrev) * cosEccAnomaly * sqrt1MinusEccPrevSquared,  
                                                        0.0                                        );                   // km/s      - Relative velocity vector, in the m1Prev rest frame

        Vector3d orbitalAngularMomentumVectorPrev = cross(separationVectorPrev, relativeVelocityVectorPrev);            // km^2 s^-1 - Specific orbital angular momentum vector 

        Vector3d eccentricityVectorPrev = cross(relativeVelocityVectorPrev, orbitalAngularMomentumVectorPrev) / (G_SN * totalMassPrev) - separationVectorPrev.hat;                                                 // --        - Laplace-Runge-Lenz vector (magnitude = eccentricity)

        m_OrbitalVelocityPreSN = relativeVelocityVectorPrev.mag;                                                        // km/s      - Set the Pre-SN orbital velocity and 
        m_uK = m_Supernova->SN_KickMagnitude() / m_OrbitalVelocityPreSN;                                                // --        - Dimensionless kick magnitude

        /////////////////////////////////////////////////////////////////////////////////////////
        // Note: In the following,
        // orbitalAngularMomentumVectorPrev defines the Z-axis, 
        // eccentricityVectorPrev defines the X-axis, and
        // (orbitalAngularMomentumVectorPrev x eccentricityVectorPrev) defines the Y-axis
        /////////////////////////////////////////////////////////////////////////////////////////
        

        /////////////////////////////////////////////////////////////////////////////////////////
        // Apply supernova natal kick and mass loss  
        //
        // Note: the code allows for mass loss and kick in the companion 
        // (due to ablation), though we currently do not apply these.
        /////////////////////////////////////////////////////////////////////////////////////////
        
        Vector3d companionRecoilVector = Vector3d(0.0, 0.0, 0.0);                                                       // km/s - The recoil of the companion due to ablation
        double m1 = m_Supernova->Mass();                                                                                // Mo   - supernova star postSN mass
        double m2 = m_Companion->Mass();                                                                                // Mo   - companion star postSN mass
        double totalMass = m1 + m2;                                                                                     // Mo   - Total binary postSN mass

        double dm1 = (m1Prev - m1);                                                                                     // Mo   - Mass difference of supernova star
        double dm2 = (m2Prev - m2);                                                                                     // Mo   - Mass difference of companion star

        Vector3d centerOfMassVelocity = (-m2Prev * dm1 / (totalMassPrev*totalMass) + m1Prev * dm2 / (totalMassPrev * totalMass)) * relativeVelocityVectorPrev 
                                         + (m1 / totalMass) * natalKickVector 
                                         + (m2 / totalMass) * companionRecoilVector;                                    // km/s       - PostSN center of mass velocity vector

        Vector3d relativeVelocityVector = relativeVelocityVectorPrev + (natalKickVector - companionRecoilVector);       // km/s       - PostSN relative velocity vector

        Vector3d orbitalAngularMomentumVector = cross(separationVectorPrev, relativeVelocityVector);                    // km^2 s^-1  - PostSN specific orbital angular momentum vector
        double   orbitalAngularMomentum = orbitalAngularMomentumVector.mag;                                             // km^2 s^-1  - PostSN specific orbital angular momentum 
        m_OrbitalAngularMomentumVector = orbitalAngularMomentumVector/orbitalAngularMomentum;                           // set unit vector here to make printing out the inclination vector easier

        Vector3d eccentricityVector = cross(relativeVelocityVector, orbitalAngularMomentumVector) / (G_SN * totalMass) 
                                      - separationVectorPrev / separationPrev;                                          // PostSN Laplace-Runge-Lenz vector
        m_Eccentricity = eccentricityVector.mag;                                                                        // PostSN eccentricity
        double eccSquared = m_Eccentricity * m_Eccentricity;                                                            // useful function of eccentricity

        double semiMajorAxis_km = (orbitalAngularMomentum*orbitalAngularMomentum) / (G_SN * totalMass * (1 - eccSquared));  // km         - PostSN semi-major axis
        m_SemiMajorAxis = semiMajorAxis_km * KM_TO_AU;                                                                  // AU         - PostSN semi-major axis 


        /////////////////////////////////////////////////////////////////////////////////////////
        // Note: similar to above,
        // orbitalAngularMomentumVector defines the Z'-axis, 
        // eccentricityVector defines the X'-axis, and
        // (orbitalAngularMomentumVector x eccentricityVector) defines the Y'-axis
        /////////////////////////////////////////////////////////////////////////////////////////
         
        UpdateSystemicVelocity(centerOfMassVelocity.RotateVector(m_ThetaE, m_PhiE, m_PsiE));                            // Update the system velocity with the new center of mass velocity


        /////////////////////////////////////////////////////////////////////////////////////////
        // Split off and evaluate depending on whether the binary is now bound or unbound
	    if (utils::Compare(m_Eccentricity, 1.0) >= 0) {                                                                     
            
            ////////////////////////////////////////
            // 
            // Binary has become unbound
            // 
            ////////////////////////////////////////

            m_Unbound = true;

            // Calculate the asymptotic Center of Mass velocity 
            double   relativeVelocityAtInfinity = (G_SN*totalMass/orbitalAngularMomentum) * std::sqrt(eccSquared - 1);
            Vector3d relativeVelocityVectorAtInfinity = relativeVelocityAtInfinity 
                                                        * (-1 * (eccentricityVector.hat / m_Eccentricity) 
                                                        + std::sqrt(1 - 1.0 / eccSquared) * cross(orbitalAngularMomentumVector.hat, eccentricityVector.hat));

            // Calculate the asymptotic velocities of Star1 (SN) and Star2 (CP)
            Vector3d component1VelocityVectorAtInfinity =  (m2 / totalMass) * relativeVelocityVectorAtInfinity + centerOfMassVelocity;
            Vector3d component2VelocityVectorAtInfinity = -(m1 / totalMass) * relativeVelocityVectorAtInfinity + centerOfMassVelocity;

            // Update the component velocities 
            m_Supernova->UpdateComponentVelocity(component1VelocityVectorAtInfinity.RotateVector(m_ThetaE, m_PhiE, m_PsiE));
            m_Companion->UpdateComponentVelocity(component2VelocityVectorAtInfinity.RotateVector(m_ThetaE, m_PhiE, m_PsiE));

            // Set Euler Angles 
            m_ThetaE = angleBetween(orbitalAngularMomentumVectorPrev, orbitalAngularMomentumVector);                   // Angle between the angular momentum unit vectors, always well defined
            m_PhiE   = _2_PI * RAND->Random(); 
            m_PsiE   = _2_PI * RAND->Random(); 
        }
        else {                     

            ////////////////////////////////////////
            // 
            // Binary is still bound 
            // 
            ////////////////////////////////////////

            // Set the component velocites to the system velocity. System velocity was already correctly set above.
             
            m_Supernova->UpdateComponentVelocity(centerOfMassVelocity.RotateVector(m_ThetaE, m_PhiE, m_PsiE));
            m_Companion->UpdateComponentVelocity(centerOfMassVelocity.RotateVector(m_ThetaE, m_PhiE, m_PsiE));

            ////////////////////////////////////////////////////////////////////////////////////
            // Calculate Euler angles - see RotateVector() in vector.cpp for details

            m_ThetaE = angleBetween(orbitalAngularMomentumVector, orbitalAngularMomentumVectorPrev);                    // Angle between the angular momentum unit vectors, always well defined

            // If the new orbital A.M. is parallel or anti-parallel to the previous orbital A.M., 
            // then the cross product is not well-defined, and we need to account for degeneracy between eccentricity vectors.
            // Also, if either eccentricity is 0.0, then the eccentricity vector is not well defined.

            if ((utils::Compare(m_ThetaE, 0.0) == 0) &&                                                                 // Is orbitalAngularMomentumVectorPrev parallel to orbitalAngularMomentumVector ...
               ((utils::Compare(eccentricityPrev,  0.0) > 0) &&                                                         // ...
                (utils::Compare(m_Eccentricity, 0.0) > 0))) {                                                           // ...and both eccentricityVectorPrev and eccentricityVector are well defined?

                double psiPlusPhi = angleBetween(eccentricityVector, eccentricityVectorPrev);                           // yes - then psi + phi is constant
                m_PhiE = _2_PI * RAND->Random();    
                m_PsiE = psiPlusPhi - m_PhiE;
            }
            else if ((utils::Compare(m_ThetaE, M_PI) == 0) &&                                                           // Is orbitalAngularMomentumVectorPrev anti-parallel to orbitalAngularMomentumVector ...
                    ((utils::Compare(eccentricityPrev,  0.0) > 0) &&                                                    // ...
                     (utils::Compare(m_Eccentricity, 0.0) > 0))) {                                                      // ...and both eccentricityVectorPrev and eccentricityVector are well defined?

                                                                                                                        // yes - then psi - phi is constant
                double psiMinusPhi = angleBetween(eccentricityVector, eccentricityVectorPrev); 
                m_PhiE = _2_PI * RAND->Random();    
                m_PsiE = psiMinusPhi + m_PhiE;
            }
            else {                                                                                                      // Neither - the cross product of the orbit normals is well-defined

                Vector3d orbitalPivotAxis = cross(orbitalAngularMomentumVectorPrev, orbitalAngularMomentumVector);      // Cross product of the orbit normals

                if ( utils::Compare(eccentricityPrev, 0.0) == 0 ) {                                                     // Is eccentricityVectorPrev well-defined?
                    m_PhiE = _2_PI * RAND->Random();                                                                    // no - set phi random
                }
                else {                                                                                                  // yes - phi is +/- angle between eccentricityVectorPrev and orbitalPivotAxis
                    
                    m_PhiE = utils::Compare( dot(eccentricityVectorPrev, orbitalAngularMomentumVector), 0.0) >= 0 ?     // Are eccentricityVectorPrev and orbitalAngularMomentumVector in the same hemisphere?
                         angleBetween(eccentricityVectorPrev, orbitalPivotAxis):                                        // yes - phi in [0,pi)
                        -angleBetween(eccentricityVectorPrev, orbitalPivotAxis);                                        // no  - phi in [-pi,0)
                }

                if ( utils::Compare(m_Eccentricity, 0.0) == 0 ) {                                                       // Is eccentricityVector well-defined?
                    m_PsiE = _2_PI * RAND->Random();                                                                    // no - set psi random 
                }                                                                                              
                else {                                                                                                  // yes - psi is +/- angle between eccentricityVector and orbitalPivotAxis

                    m_PsiE = utils::Compare( dot(eccentricityVector, orbitalAngularMomentumVectorPrev), 0.0) >= 0 ?     // Are eccentricityVector and orbitalAngularMomentumVectorPrev in the same hemisphere?
                         angleBetween(eccentricityVector, orbitalPivotAxis):                                            // yes - psi in [0,pi)
                        -angleBetween(eccentricityVector, orbitalPivotAxis);                                            // no  - psi in [-pi,0)
                }
            }

            // Note: There is some evidence for evolution of periapsis in mass transferring binaries (see e.g Dosopoulou & Kalogera 2016, 2018). 
            // This should be investigated in more depth, but until then, we assume that the periapsis *may* evolve, 
            // and accordingly randomize the angle of periapsis around the new orbital angular momentum, (i.e, Psi)
            // - RTW 15/05/20
            m_PsiE = _2_PI * RAND->Random();
        }

        // Undefine the pre-processor commands 
        #undef cross
        #undef dot
        #undef angleBetween
        #undef mag        
        #undef hat
    }

    //////////////////////////
    // Do for all systems 

    // Set remaining post-SN values
    double totalMass = m_Supernova->Mass() + m_Companion->Mass();                                                       // Total Mass 
    double reducedMass = m_Supernova->Mass() * m_Companion->Mass() / totalMass;                                         // Reduced Mass
    m_Supernova->SetOrbitalEnergyPostSN(CalculateOrbitalEnergy(reducedMass, totalMass, m_SemiMajorAxis));               // Orbital energy

    m_IPrime    = m_ThetaE;                                                                                             // Inclination angle between preSN and postSN orbital planes 
    m_CosIPrime = cos(m_IPrime);

    (void)PrintSupernovaDetails();                                                                                      // Log record to supernovae logfile
    m_Supernova->ClearCurrentSNEvent();

    return true;
}


/*
 * Update the Center of Mass velocity and speed of the binary system following a Supernova.
 *
 * This simply adds a new CoM vector to the existing one, but note that the new vector
 * must be rotated into the old coordinate frame (see vector3d.h)
 *
 * void UpdateSystemicVelocity(const double p_newVelocity[3] )
 *
 * @param   [IN]    p_newVelocity(3)             3D velocity vector in km/s to add to current velocity vector
 */
void BaseBinaryStar::UpdateSystemicVelocity(Vector3d p_newVelocity) {

    // Update the systemic velocity
    m_SystemicVelocity += p_newVelocity;            
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

    m_SupernovaState = SN_STATE::NONE;                                                                                  // not yet determined
    
    if (m_Star1->IsSNevent()) {                                                                                         // star1 supernova
        m_SupernovaState = SN_STATE::STAR1;                                                                             // star1

        // resolve star1 supernova
        m_Supernova = m_Star1;                                                                                          // supernova
        m_Companion = m_Star2;                                                                                          // companion
        (void)ResolveSupernova();                                                                                       // resolve supernova
    }

    if (m_Star2->IsSNevent()) {                                                                                         // star2 supernova                                                                                                        
        m_SupernovaState = m_SupernovaState == SN_STATE::NONE                                                           // star1 not supernova?
                            ? SN_STATE::STAR2                                                                           // yes - just star2
                            : SN_STATE::BOTH;                                                                           // no - both 

        m_Supernova = m_Star2;                                                                                          // supernova
        m_Companion = m_Star1;                                                                                          // companion
        (void)ResolveSupernova();                                                                                       // resolve supernova
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

	double eccentricity     = Eccentricity();								                                            // current eccentricity (before CEE)
    double semiMajorAxisRsol= SemiMajorAxisRsol();                                                                      // current semi-major axis in default units, Rsol (before CEE)
    double periastronRsol   = PeriastronRsol();                                                                         // periastron, Rsol (before CEE)
    double rRLd1Rsol = periastronRsol * CalculateRocheLobeRadius_Static(m_Star1->Mass(), m_Star2->Mass());              // Roche-lobe radius at periastron in Rsol at the moment where CEE begins, seen by star1
    double rRLd2Rsol = periastronRsol * CalculateRocheLobeRadius_Static(m_Star2->Mass(), m_Star1->Mass());              // Roche-lobe radius at periastron in Rsol at the moment where CEE begins, seen by star2
    
    bool isDonorMS = false;                                                                                             // check for main sequence donor
    if (OPTIONS->AllowMainSequenceStarToSurviveCommonEnvelope()) {                                                      // allow main sequence stars to survive CEE?
        if (m_Star1->IsOneOf(ALL_MAIN_SEQUENCE)) {                                                                      // yes - star1 MS_LTE_07, MS_GT_07 or NAKED_HELIUM_STAR_MS?
            isDonorMS      = isDonorMS || m_Star1->IsRLOF();                                                            // yes - donor MS?
            m_Mass1Final = m_Star1->Mass();                                                                             // set mass
            m_MassEnv1   = 0.0;                                                                                         // no envelope
        }
        else {                                                                                                          // no, star1 not MS_LTE_07, MS_GT_07 or NAKED_HELIUM_STAR_MS
            m_Mass1Final = m_Star1->CoreMass();                                                                         // set mass
            m_MassEnv1   = m_Star1->Mass() - m_Star1->CoreMass();                                                       // and envelope
        }

        if (m_Star2->IsOneOf(ALL_MAIN_SEQUENCE)) {                                                                      // star2 MS_LTE_07, MS_GT_07 or NAKED_HELIUM_STAR_MS?
            isDonorMS      = isDonorMS || m_Star2->IsRLOF();                                                            // yes - donor MS?
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
    
    if( OPTIONS->CommonEnvelopeFormalism() == CE_FORMALISM::ENERGY ) {
        double k1            = m_Star1->IsOneOf(COMPACT_OBJECTS) ? 0.0 : (2.0 / (lambda1 * alphaCE)) * m_Star1->Mass() * m_MassEnv1 / m_Star1->Radius();
        double k2            = m_Star2->IsOneOf(COMPACT_OBJECTS) ? 0.0 : (2.0 / (lambda2 * alphaCE)) * m_Star2->Mass() * m_MassEnv2 / m_Star2->Radius();
        double k3            = m_Star1->Mass() * m_Star2->Mass() / periastronRsol;                                      //assumes immediate circularisation at periastron at start of CE
        double k4            = (m_Mass1Final * m_Mass2Final);
        double aFinalRsol    = k4 / (k1 + k2 + k3);
        m_SemiMajorAxis      = aFinalRsol*RSOL_TO_AU;
    }
    
    // Two-stage common envelope, Hirai & Mandel (2022)
    else if( OPTIONS->CommonEnvelopeFormalism() == CE_FORMALISM::TWO_STAGE ) {
        double convectiveEnvelopeMass1  = m_Star1->CalculateConvectiveEnvelopeMass();
        double radiativeIntershellMass1 = m_MassEnv1 - convectiveEnvelopeMass1;
        double endOfFirstStageMass1     = m_Mass1Final + radiativeIntershellMass1;
        double convectiveEnvelopeMass2  = m_Star2->CalculateConvectiveEnvelopeMass();
        double radiativeIntershellMass2 = m_MassEnv2 - convectiveEnvelopeMass2;
        double endOfFirstStageMass2     = m_Mass2Final + radiativeIntershellMass2;
        
        // Stage 1: convective envelope removal on a dynamical timescale; assumes lambda = 1.5, motivated by bottom panel of Figure 3 of Hirai & Mandel 2022, including internal energy
        double k1            = m_Star1->IsOneOf(COMPACT_OBJECTS) ? 0.0 : (2.0 / (1.5 * alphaCE)) * m_Star1->Mass() * convectiveEnvelopeMass1 / m_Star1->Radius();
        double k2            = m_Star2->IsOneOf(COMPACT_OBJECTS) ? 0.0 : (2.0 / (1.5 * alphaCE)) * m_Star2->Mass() * convectiveEnvelopeMass2 / m_Star2->Radius();
        double k3            = m_Star1->Mass() * m_Star2->Mass() / periastronRsol;                                      //assumes immediate circularisation at periastron at start of CE
        double k4            = (endOfFirstStageMass1 * endOfFirstStageMass2);
        double aFinalRsol    = k4 / (k1 + k2 + k3);
        m_SemiMajorAxis      = aFinalRsol*RSOL_TO_AU;
        
        // Stage 2: radiative envelope removal on a thermal timescale; assumed to be fully non-conservative
        if( utils::Compare(radiativeIntershellMass1, 0.0) > 0 ) {
            m_SemiMajorAxis = CalculateMassTransferOrbit(endOfFirstStageMass1, -radiativeIntershellMass1, *m_Star2, 0.0);
        }
        if( utils::Compare(radiativeIntershellMass2, 0.0) > 0 ) {
            m_SemiMajorAxis = CalculateMassTransferOrbit(endOfFirstStageMass2, -radiativeIntershellMass2, *m_Star1, 0.0);
        }
    }
    
    else {                                                                                                              // Invalid CE formalism
        SHOW_WARN_STATIC(ERROR::UNKNOWN_CE_FORMALISM,                                                                   // show warning
                         "Orbital properties unchanged by CE",
                         OBJECT_TYPE::BASE_BINARY_STAR,
                         STELLAR_TYPE::BINARY_STAR);
    }
    

	double rRLdfin1        = m_SemiMajorAxis * CalculateRocheLobeRadius_Static(m_Mass1Final, m_Mass2Final);             // Roche-lobe radius in AU after CEE, seen by star1
	double rRLdfin2        = m_SemiMajorAxis * CalculateRocheLobeRadius_Static(m_Mass2Final, m_Mass1Final);             // Roche-lobe radius in AU after CEE, seen by star2
    double rRLdfin1Rsol    = rRLdfin1 * AU_TO_RSOL;                                                                     // Roche-lobe radius in Rsol after CEE, seen by star1
    double rRLdfin2Rsol    = rRLdfin2 * AU_TO_RSOL;                                                                     // Roche-lobe radius in Rsol after CEE, seen by star2
    m_Eccentricity         = 0.0;                                                                                       // We assume that a common envelope event (CEE) circularises the binary

    m_Star1->ResolveCommonEnvelopeAccretion(m_Mass1Final);                                                              // update star1's mass after CE accretion
    m_Star2->ResolveCommonEnvelopeAccretion(m_Mass2Final);                                                              // update star2's mass after CE accretion

    // update stellar type after losing its envelope. Star1, Star2 or both if double CEE.

    if (isDonorMS || (!envelopeFlag1 && !envelopeFlag2)) {                                                              // stellar merger
        m_MassTransferTrackerHistory = MT_TRACKING::MERGER; 
        m_Flags.stellarMerger        = true;
    }
    else if ( (m_Star1->DetermineEnvelopeType()==ENVELOPE::RADIATIVE && !m_Star1->IsOneOf(ALL_MAIN_SEQUENCE)) ||
              (m_Star2->DetermineEnvelopeType()==ENVELOPE::RADIATIVE && !m_Star2->IsOneOf(ALL_MAIN_SEQUENCE)) ) {       // check if we have a non-MS radiative-envelope star
        m_CEDetails.optimisticCE = true;
        if(!OPTIONS->AllowRadiativeEnvelopeStarToSurviveCommonEnvelope() ) {                                            // stellar merger
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
    if (OPTIONS->CHEMode() != CHE_MODE::NONE) m_Star1->SetOmega(OrbitalAngularVelocity());
    if (OPTIONS->CHEMode() != CHE_MODE::NONE) m_Star2->SetOmega(OrbitalAngularVelocity());
    
    m_Star1->SetPostCEEValues();                                                                                        // squirrel away post CEE stellar values for star 1
    m_Star2->SetPostCEEValues();                                                                                        // squirrel away post CEE stellar values for star 2
    SetPostCEEValues(m_SemiMajorAxis * AU_TO_RSOL, m_Eccentricity, rRLdfin1Rsol, rRLdfin2Rsol);                         // squirrel away post CEE binary values (checks for post-CE RLOF, so should be done at end)

    if (m_RLOFDetails.immediateRLOFPostCEE == true && !OPTIONS->AllowImmediateRLOFpostCEToSurviveCommonEnvelope()) {    // Is there immediate post-CE RLOF which is not allowed?
            m_MassTransferTrackerHistory = MT_TRACKING::MERGER;
            m_Flags.stellarMerger = true;
    }

    (void)PrintCommonEnvelope();                                                                                        // print (log) common envelope details
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
    double q = p_MassPrimary / p_MassSecondary;
    double qCubeRoot = PPOW(q, 1.0 / 3.0);                                                                              // cube roots are expensive, only compute once
    return 0.49 / (0.6 + log(1.0 + qCubeRoot) / qCubeRoot / qCubeRoot);
}


/*
 * Calculate the fraction of specific angular momentum with which the non-accreted mass leaves the system
 *
 * This is gamma (as in Pols's notes) or jloss (as in Belczynski et al. 2008
 * which is the fraction of specific angular momentum with which the non-accreted mass leaves the system.
 * Macleod_linear comes from Willcox et al. (2022)
 *
 * Updates class member variable m_Error      JR: todo: revisit error handling (this could be a const function)
 * 
 * 
 * Calculation is based on user-specified Angular Momentum Loss prescription
 *
 *
 * double CalculateGammaAngularMomentumLoss(const double p_DonorMass, const double p_AccretorMass)
 *
 * @param   [IN]    p_DonorMass                 The mass of the donor (Msol)
 * @param   [IN]    p_AccretorMass              The mass of the accretor (Msol)
 * @return                                      The fraction of specific angular momentum with which the non-accreted mass leaves the system
 */
double BaseBinaryStar::CalculateGammaAngularMomentumLoss(const double p_DonorMass, const double p_AccretorMass) {

	double gamma;

	switch (OPTIONS->MassTransferAngularMomentumLossPrescription()) {                                                       // which precription?
        case MT_ANGULAR_MOMENTUM_LOSS_PRESCRIPTION::JEANS                : gamma = p_AccretorMass / p_DonorMass; break;     // vicinity of the donor
        case MT_ANGULAR_MOMENTUM_LOSS_PRESCRIPTION::ISOTROPIC_RE_EMISSION: gamma = p_DonorMass / p_AccretorMass; break;     // vicinity of the accretor
        case MT_ANGULAR_MOMENTUM_LOSS_PRESCRIPTION::CIRCUMBINARY_RING    : gamma = (M_SQRT2 * (p_DonorMass + p_AccretorMass) * (p_DonorMass + p_AccretorMass)) / (p_DonorMass * p_AccretorMass); break; // Based on the assumption that a_ring ~= 2*a*, Vinciguerra+, 2020 
        case MT_ANGULAR_MOMENTUM_LOSS_PRESCRIPTION::MACLEOD_LINEAR       : {                                                // Linear interpolation on separation between accretor and L2 point
            double q = p_AccretorMass / p_DonorMass;
            // interpolate in separation between a_acc and a_L2, both normalized to units of separation a
            double aL2 = std::sqrt(M_SQRT2);  // roughly, coincides with CIRCUMBINARY_RING def above
            double aAcc = 1/(1+q);
            double aGamma = aAcc + (aL2 - aAcc)*OPTIONS->MassTransferJlossMacLeodLinearFraction();
            gamma = aGamma*aGamma*(1+q)*(1+q)/q;
            break;
        }
        case MT_ANGULAR_MOMENTUM_LOSS_PRESCRIPTION::ARBITRARY            : gamma = OPTIONS->MassTransferJloss(); break;
        default:                                                                                                            // unknown mass transfer angular momentum loss prescription - shouldn't happen
            gamma = 1.0;                                                                                                    // default value
            m_Error = ERROR::UNKNOWN_MT_ANGULAR_MOMENTUM_LOSS_PRESCRIPTION;                                                 // set error value
            SHOW_WARN(m_Error);                                                                                             // warn that an error occurred
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

    double semiMajorAxis   = m_SemiMajorAxis;                                                                   // new semi-major axis value - default is no change
    double massA           = p_Accretor.Mass();                                                                 // accretor mass
    double massD           = p_DonorMass;                                                                       // donor mass
    double massAtimesMassD = massA * massD;                                                                     // accretor mass * donor mass
    double massAplusMassD  = massA + massD;                                                                     // accretor mass + donor mass
    double jOrb            = (massAtimesMassD / massAplusMassD) * std::sqrt(semiMajorAxis * G1 * massAplusMassD); // orbital angular momentum
    double jLoss;                                                                                               // specific angular momentum carried away by non-conservative mass transfer
    
    if (utils::Compare(p_DeltaMassDonor, 0.0) >= 0) {                                                           // no mass loss from donor, nothing to do here
        return semiMajorAxis;
    }
    int numberIterations   = fmax( floor (fabs(p_DeltaMassDonor/(MAXIMUM_MASS_TRANSFER_FRACTION_PER_STEP*massD))), 1); // number of iterations

    double dM              = p_DeltaMassDonor / numberIterations;                                               // mass change per time step

    for(int i = 0; i < numberIterations ; i++) {
        
        jLoss = CalculateGammaAngularMomentumLoss(massD, massA);
        jOrb = jOrb + ((jLoss * jOrb * (1.0 - p_FractionAccreted) / massAplusMassD) * dM);
        semiMajorAxis = semiMajorAxis + (((-2.0 * dM / massD) * (1.0 - (p_FractionAccreted * (massD / massA)) - ((1.0 - p_FractionAccreted) * (jLoss + 0.5) * (massD / massAplusMassD)))) * semiMajorAxis);

        massD          = massD + dM;
        massA          = massA - (dM * p_FractionAccreted);
        massAplusMassD = massA + massD;
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
 * double CalculateZetaRocheLobe()
 *
 * @param   [IN]    p_jLoss                     Specific angular momentum with which mass is lost during non-conservative mass transfer
 *                                              (Podsiadlowski et al. 1992, Beta: specific angular momentum of matter [2Pia^2/P])
 * @return                                      Roche Lobe response
 */
double BaseBinaryStar::CalculateZetaRocheLobe(const double p_jLoss) const {

    double donorMass    = m_Donor->Mass();                  // donor mass
    double accretorMass = m_Accretor->Mass();               // accretor mass
    double beta         = m_FractionAccreted;               // fraction of mass accreted by accretor
    double gamma        = p_jLoss;

    double q = donorMass / accretorMass;

    double q_1_3 = PPOW(q, 1.0 / 3.0);

    double k1 = -2.0 * (1.0 - (beta * q) - (1.0 - beta) * (gamma + 0.5) * (q / (1.0 + q)));
    double k2 = (2.0 / 3.0) - q_1_3 * (1.2 * q_1_3 + 1.0 / (1.0 + q_1_3)) / (3.0 * (0.6 * q_1_3 * q_1_3 + log(1.0 + q_1_3)));
    double k3 = 1.0 + (beta * q);

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
    
    if (Unbound())
        return;                                                                                                                 // do nothing for unbound binaries
    
    if (!OPTIONS->UseMassTransfer()) return;                                                                                    // mass transfer not enabled - nothing to do
    
    if (!m_Star1->IsRLOF() && !m_Star2->IsRLOF()) return;                                                                       // neither star is overflowing its Roche Lobe - no mass transfer - nothing to do
    
    if (OPTIONS->CHEMode() != CHE_MODE::NONE && HasTwoOf({STELLAR_TYPE::CHEMICALLY_HOMOGENEOUS}) && HasStarsTouching()) {       // CHE enabled and both stars CH?
        m_Flags.stellarMerger = true;
        return;
    }

    if (m_Star1->IsRLOF() && m_Star2->IsRLOF()) {                                                                               // both stars overflowing their Roche Lobe?
        m_CEDetails.CEEnow = true;                                                                                              // yes - common envelope event - no mass transfer
        return;                                                                                                                 // and return - nothing (else) to do
    }

    // one, and only one, star is overflowing its Roche Lobe - resolve mass transfer

    m_Donor    = m_Star2->IsRLOF() ? m_Star2 : m_Star1;                                                                         // donor is primary unless secondary is overflowing its Roche Lobe
    m_Accretor = m_Star2->IsRLOF() ? m_Star1 : m_Star2;                                                                         // accretor is secondary unless secondary is overflowing its Roche Lobe

    // Add event to MT history of the donor
    m_Donor->UpdateMassTransferDonorHistory();

    // Calculate accretion fraction if stable
    // This passes the accretor's Roche lobe radius to m_Accretor->CalculateThermalMassAcceptanceRate()
    // just in case MT_THERMALLY_LIMITED_VARIATION::RADIUS_TO_ROCHELOBE is used; otherwise, the radius input is ignored
    double accretorRLradius = CalculateRocheLobeRadius_Static(m_Accretor->Mass(), m_Donor->Mass()) * AU_TO_RSOL * m_SemiMajorAxis * (1.0 - m_Eccentricity);
    bool donorIsHeRich = m_Donor->IsOneOf(He_RICH_TYPES); 
    std::tie(std::ignore, m_FractionAccreted) = m_Accretor->CalculateMassAcceptanceRate(m_Donor->CalculateThermalMassLossRate(),
                                                                                        m_Accretor->CalculateThermalMassAcceptanceRate(accretorRLradius),
                                                                                        donorIsHeRich);

    double aInitial = m_SemiMajorAxis;                                                                                          // semi-major axis in default units, AU, current timestep
    double aFinal;                                                                                                              // semi-major axis in default units, AU, after next timestep
    double jLoss    = m_JLoss;                            		                                                                // specific angular momentum with which mass is lost during non-conservative mass transfer, current timestep

    if (OPTIONS->MassTransferAngularMomentumLossPrescription() != MT_ANGULAR_MOMENTUM_LOSS_PRESCRIPTION::ARBITRARY) {           // arbitrary angular momentum loss prescription?
        jLoss = CalculateGammaAngularMomentumLoss();                                                                            // no - re-calculate angular momentum
    }

    // Calculate conditions for automatic (in)stability for case BB
    bool caseBBAlwaysStable           = OPTIONS->CaseBBStabilityPrescription() == CASE_BB_STABILITY_PRESCRIPTION::ALWAYS_STABLE;
    bool caseBBAlwaysUnstable         = OPTIONS->CaseBBStabilityPrescription() == CASE_BB_STABILITY_PRESCRIPTION::ALWAYS_UNSTABLE;
    bool caseBBAlwaysUnstableOntoNSBH = OPTIONS->CaseBBStabilityPrescription() == CASE_BB_STABILITY_PRESCRIPTION::ALWAYS_STABLE_ONTO_NSBH;
    bool donorIsHeHGorHeGB            = m_Donor->IsOneOf({ STELLAR_TYPE::NAKED_HELIUM_STAR_HERTZSPRUNG_GAP, STELLAR_TYPE::NAKED_HELIUM_STAR_GIANT_BRANCH });
    bool accretorIsNSorBH             = m_Accretor->IsOneOf({ STELLAR_TYPE::NEUTRON_STAR, STELLAR_TYPE::BLACK_HOLE });
    bool accretorIsWD                 = m_Accretor->IsOneOf(WHITE_DWARFS); 

    // Determine stability
    bool isUnstable;
    if (donorIsHeHGorHeGB && (caseBBAlwaysStable || caseBBAlwaysUnstable || (caseBBAlwaysUnstableOntoNSBH && accretorIsNSorBH))) { // Determine stability based on case BB 
        isUnstable = (caseBBAlwaysUnstable || (caseBBAlwaysUnstableOntoNSBH && accretorIsNSorBH));                              // Already established that donor is HeHG or HeGB - need to check if new case BB prescriptions are added
    } 
    else if (accretorIsWD && (m_Accretor->WhiteDwarfAccretionRegime() == ACCRETION_REGIME::HELIUM_WHITE_DWARF_HYDROGEN_ACCUMULATION)) { 
        isUnstable = true;
        if (!m_Donor->IsOneOf(GIANTS)) m_Flags.stellarMerger = true;
    }
    else if (OPTIONS->QCritPrescription() != QCRIT_PRESCRIPTION::NONE) {                                                           // Determine stability based on critical mass ratios

        // NOTE: Critical mass ratio is defined as mAccretor/mDonor
        double qCrit = m_Donor->CalculateCriticalMassRatio(m_Accretor->IsDegenerate());

        isUnstable = (m_Accretor->Mass()/m_Donor->Mass()) < qCrit;
        m_FractionAccreted = 1.0;                                                                                               // Accretion is assumed fully conservative for qCrit calculations
    }
    else {                                                                                                                      // Determine stability based on zetas

        m_ZetaLobe = CalculateZetaRocheLobe(jLoss);
        m_ZetaStar = m_Donor->CalculateZetaAdiabatic(); 

        isUnstable = (utils::Compare(m_ZetaStar, m_ZetaLobe) < 0);
    }

    // Evaluate separately for stable / unstable MT
    if (isUnstable) {                                                                                                           // Unstable Mass Transfer
         m_CEDetails.CEEnow = true;
    }
    else {                                                                                                                      // Stable MT
            
        m_MassTransferTrackerHistory = m_Donor == m_Star1                                                                       // record what happened - for later printing
            ? MT_TRACKING::STABLE_1_TO_2_SURV
            : MT_TRACKING::STABLE_2_TO_1_SURV; 

        double massDiffDonor;
        double envMassDonor  = m_Donor->Mass() - m_Donor->CoreMass();
        bool isEnvelopeRemoved = false;

        if (utils::Compare(m_Donor->CoreMass(), 0) > 0 && utils::Compare(envMassDonor, 0) > 0) {                                // donor has a core and an envelope
            massDiffDonor = -envMassDonor;                                                                                      // set donor mass loss to (negative of) the envelope mass
            isEnvelopeRemoved = true;
        }
        else{                                                                                                                   // donor has no envelope
            massDiffDonor = -MassLossToFitInsideRocheLobe(this, m_Donor, m_Accretor, m_FractionAccreted);                       // use root solver to determine how much mass should be lost from the donor to allow it to fit within the Roche lobe
            m_Donor->UpdateMinimumCoreMass();                                                                                   // reset the minimum core mass following case A
        } 
        double massGainAccretor = -massDiffDonor * m_FractionAccreted;                                                          // set accretor mass gain to mass loss * conservativeness

        m_Donor->SetMassTransferDiffAndResolveWDShellChange(massDiffDonor);                                                                            // set new mass of donor
        m_Accretor->SetMassTransferDiffAndResolveWDShellChange(massGainAccretor);                                                                      // set new mass of accretor

        aFinal = CalculateMassTransferOrbit(m_Donor->Mass(), massDiffDonor, *m_Accretor, m_FractionAccreted);                   // calculate new orbit
        m_aMassTransferDiff = aFinal - aInitial;                                                                                // set change in orbit (semi-major axis)
                                                                                                                    
        STELLAR_TYPE stellarTypeDonor = m_Donor->StellarType();                                                                 // donor stellar type before resolving envelope loss
        if (isEnvelopeRemoved) m_Donor->ResolveEnvelopeLossAndSwitch();                                                         // if this was an envelope stripping episode, resolve envelope loss
        if (m_Donor->StellarType() != stellarTypeDonor) {                                                                       // stellar type change?
            (void)PrintDetailedOutput(m_Id, BSE_DETAILED_RECORD_TYPE::STELLAR_TYPE_CHANGE_DURING_MT);                           // yes - print (log) detailed output
        }
        
        // Check if this was stable mass transfer after a CEE
        if (m_CEDetails.CEEcount > 0 && !m_RLOFDetails.stableRLOFPostCEE) {
            m_RLOFDetails.stableRLOFPostCEE = m_MassTransferTrackerHistory == MT_TRACKING::STABLE_2_TO_1_SURV ||
                                              m_MassTransferTrackerHistory == MT_TRACKING::STABLE_1_TO_2_SURV;
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
		    m_MassTransfer = true;                                                                                              // ... mass transfer
            m_CEDetails.CEEnow = false;                                                                                         // no common envelope

		    if (OPTIONS->CirculariseBinaryDuringMassTransfer()) {                                                               // circularise binary to the periapsis separation?
                m_SemiMajorAxis *= OPTIONS->AngularMomentumConservationDuringCircularisation()                                  // yes - conserve angular momentum?
                                        ? (1.0 - (m_Eccentricity * m_Eccentricity))                                             // yes - conserve angular momentum
                                        : (1.0 - m_Eccentricity);                                                               // no - angular momentum not conserved, circularise at periapsis

			    m_Eccentricity = 0.0;

                m_Star1->InitialiseMassTransfer(m_CEDetails.CEEnow, m_SemiMajorAxis, m_Eccentricity);                           // re-initialise mass transfer for star1
                m_Star2->InitialiseMassTransfer(m_CEDetails.CEEnow, m_SemiMajorAxis, m_Eccentricity);                           // re-initialise mass transfer for star2
                
			    // ALEJANDRO - 23/11/2016 - Bug fix for systems which enter MT being eccentric.
			    // Previous values have to be the ones for periastron as later orbit is modified according to previous values.
			    // If you don't do this, you end up modifying pre-MT pre-circularisation orbit
			    // JR: todo: check that this is proper functionality, or just a kludge - if kludge, resolve it
			    m_SemiMajorAxisPrev = m_SemiMajorAxis;
			    m_EccentricityPrev = m_Eccentricity;
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
 *                             const double p_Star1Radius,
 *                             const double p_Star2Radius,
 *                             const double p_Star1_SpinAngularVelocity,
 *                             const double p_Star1_SpinAngularVelocity,
 *                             const double p_Star1_GyrationRadius,
 *                             const double p_Star2_GyrationRadius)
 *
 * @param   [IN]    p_SemiMajorAxis             Semi-major axis of the binary
 * @param   [IN]    p_Star1Mass                 Mass of star 1
 * @param   [IN]    p_Star2Mass                 Mass of star 2
 * @param   [IN]    p_Star1Radius               Radius of star 1
 * @param   [IN]    p_Star2Radius               Radius of star 2
 * @param   [IN]    p_Star1_SpinAngularVelocity Spin angular velocity of star 1
 * @param   [IN]    p_Star1_SpinAngularVelocity Spin angular velocity of star 1
 * @param   [IN]    p_Star1_GyrationRadius      Gyration radius of star 1
 * @param   [IN]    p_Star2_GyrationRadius      Gyration radius of star 2
 * @return                                      Total energy of the binary
 */
double BaseBinaryStar::CalculateTotalEnergy(const double p_SemiMajorAxis,
                                            const double p_Star1Mass,
                                            const double p_Star2Mass,
                                            const double p_Star1Radius,
                                            const double p_Star2Radius,
                                            const double p_Star1_SpinAngularVelocity,
                                            const double p_Star2_SpinAngularVelocity,
                                            const double p_Star1_GyrationRadius,
                                            const double p_Star2_GyrationRadius) const {
	double m1  = p_Star1Mass;
	double m2  = p_Star2Mass;

	double R1  = p_Star1Radius;
	double R2  = p_Star2Radius;

	double w1  = p_Star1_SpinAngularVelocity;
	double w2  = p_Star2_SpinAngularVelocity;

	double ks1 = p_Star1_GyrationRadius;
	double ks2 = p_Star2_GyrationRadius;

    constexpr double RSOL_TO_AU_2 = RSOL_TO_AU * RSOL_TO_AU;

	double 	Is1  = ks1 * m1 * R1 * R1 * RSOL_TO_AU_2;
	double 	Is2  = ks2 * m2 * R2 * R2 * RSOL_TO_AU_2;

	return (0.5 * Is1 * w1 * w1) + (0.5 * Is2 * w2 * w2) - (0.5 * G1 * m1 * m2 / p_SemiMajorAxis);
}


/*
 * Calculate the angular momentum of the binary
 *
 * The angular momentum consists of the spin angular momenta of the two stars and the orbital angular momentum of the binary
 *
 *
 * double CalculateAngularMomentum(const double p_SemiMajorAxis,
 *                             const double p_Eccentricity,
 *                             const double p_Star1Mass,
 *                             const double p_Star2Mass,
 *                             const double p_Star1Radius,
 *                             const double p_Star2Radius,
 *                             const double p_Star1_SpinAngularVelocity,
 *                             const double p_Star1_SpinAngularVelocity,
 *                             const double p_Star1_GyrationRadius,
 *                             const double p_Star2_GyrationRadius)
 *
 * @param   [IN]    p_SemiMajorAxis             Semi-major axis of the binary
 * @param   [IN]    p_Eccentricity              Eccentricity of the binary
 * @param   [IN]    p_Star1Mass                 Mass of the primary
 * @param   [IN]    p_Star2Mass                 Mass of the secondary
 * @param   [IN]    p_Star1Radius               Radius of the primary
 * @param   [IN]    p_Star2Radius               Radius of the secondary
 * @param   [IN]    p_Star1_SpinAngularVelocity Orbital frequency of the primary
 * @param   [IN]    p_Star1_SpinAngularVelocity Orbital frequency of the secondary
 * @param   [IN]    p_Star1_GyrationRadius      Gyration radius of the primary
 * @param   [IN]    p_Star2_GyrationRadius      Gyration radius of the secondary
 * @return                                      Angular momentum of the binary
 */
double BaseBinaryStar::CalculateAngularMomentum(const double p_SemiMajorAxis,
                                                const double p_Eccentricity,
                                                const double p_Star1Mass,
                                                const double p_Star2Mass,
                                                const double p_Star1Radius,
                                                const double p_Star2Radius,
                                                const double p_Star1_SpinAngularVelocity,
                                                const double p_Star2_SpinAngularVelocity,
                                                const double p_Star1_GyrationRadius,
                                                const double p_Star2_GyrationRadius) const {
	double m1 = p_Star1Mass;
	double m2 = p_Star2Mass;

	double R1 = p_Star1Radius * RSOL_TO_AU;
	double R2 = p_Star2Radius * RSOL_TO_AU;

	double w1 = p_Star1_SpinAngularVelocity;
	double w2 = p_Star2_SpinAngularVelocity;

	double ks1 = p_Star1_GyrationRadius;
	double ks2 = p_Star2_GyrationRadius;

	double Is1  = ks1 * m1 * R1 * R1;
	double Is2  = ks2 * m2 * R2 * R2;
    double Jorb = ((m1 * m2) / (m1 + m2)) * std::sqrt(G1 * (m1 + m2) * p_SemiMajorAxis * (1.0 - (p_Eccentricity * p_Eccentricity)));

	return (Is1 * w1) + (Is2 * w2) + Jorb;
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
    m_OrbitalEnergyPrev                = m_OrbitalEnergy;
    m_OrbitalAngularMomentumPrev       = m_OrbitalAngularMomentum;

    double totalMass                        = m_Star1->Mass() + m_Star2->Mass();
    double reducedMass                      = (m_Star1->Mass() * m_Star2->Mass()) / totalMass;
    m_OrbitalEnergy                    = CalculateOrbitalEnergy(reducedMass, totalMass, m_SemiMajorAxis);
    m_OrbitalAngularMomentum           = CalculateOrbitalAngularMomentum(reducedMass, totalMass, m_SemiMajorAxis);

    // Calculate total energy and angular momentum using regular conservation of energy, especially useful for checking tides and rotational effects
    m_TotalEnergy                 = CalculateTotalEnergy();
    m_TotalAngularMomentum        = CalculateAngularMomentum();
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
    (void)m_Star1->UpdateAttributes(m_Star1->MassPrev() - m_Star1->Mass() + m_Star1->MassLossDiff() + m_Star1->MassTransferDiff(), 0.0);    // update mass for star1
    m_Star1->UpdateInitialMass();                                                                       // update effective initial mass of star1 (MS, HG & HeMS)
    m_Star1->UpdateAgeAfterMassLoss();                                                                  // update age of star1
    m_Star1->ApplyMassTransferRejuvenationFactor();                                                     // apply age rejuvenation factor for star1
    m_Star1->UpdateAttributes(0.0, 0.0, true);

    // rinse and repeat for star2
    (void)m_Star2->UpdateAttributes(m_Star2->MassPrev() - m_Star2->Mass() + m_Star2->MassLossDiff() + m_Star2->MassTransferDiff(), 0.0);    // update mass for star2
    m_Star2->UpdateInitialMass();                                                                       // update effective initial mass of star 2 (MS, HG & HeMS)
    m_Star2->UpdateAgeAfterMassLoss();                                                                  // update age of star2
    m_Star2->ApplyMassTransferRejuvenationFactor();                                                     // apply age rejuvenation factor for star2
    m_Star2->UpdateAttributes(0.0, 0.0, true);
    
    // update binary
    m_SemiMajorAxis = m_SemiMajorAxisPrev + m_aMassLossDiff + m_aMassTransferDiff;
    
    //Envelope ejection for convective envelope stars exceeding threshold luminosity to mass ratio: assume the entire envelope was lost on timescales long relative to the orbit
    if(m_Star1->EnvelopeJustExpelledByPulsations() || m_Star2->EnvelopeJustExpelledByPulsations()) {
        m_SemiMajorAxis /=  (2.0 - ((m_Star1->MassPrev() + m_Star2->MassPrev()) / (m_Star1->Mass() + m_Star2->Mass())));    // update separation in response to pulsational mass loss
        m_Star1->ResetEnvelopeExpulsationByPulsations();
        m_Star2->ResetEnvelopeExpulsationByPulsations();
    }

    // if CHE enabled, update rotational frequency for constituent stars - assume tidally locked
    if (OPTIONS->CHEMode() != CHE_MODE::NONE) m_Star1->SetOmega(OrbitalAngularVelocity());
    if (OPTIONS->CHEMode() != CHE_MODE::NONE) m_Star2->SetOmega(OrbitalAngularVelocity());

    CalculateEnergyAndAngularMomentum();                                                                // perform energy and angular momentum calculations

    if ((m_Star1->StellarType() != stellarType1) || (m_Star2->StellarType() != stellarType2)) {         // stellar type change?
        (void)PrintDetailedOutput(m_Id, BSE_DETAILED_RECORD_TYPE::STELLAR_TYPE_CHANGE_DURING_MASS_RESOLUTION); // yes - print (log) detailed output
    }
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
 * JR: This is the binary evolution (post stellar evolution).  Maybe we should rename this EvolveBinary()?
 *     The constituent srats have been evolved for a single timestep before entering this function, and here
 *     we update (evolve?) the binary in response to the stellar evolution of the components.
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
        !(OPTIONS->CHEMode() != CHE_MODE::NONE && HasTwoOf({STELLAR_TYPE::CHEMICALLY_HOMOGENEOUS}))) {                  // yes - avoid CEE if CH+CH

        ResolveCommonEnvelopeEvent();                                                                                   // resolve CEE - immediate event
        (void)PrintDetailedOutput(m_Id, BSE_DETAILED_RECORD_TYPE::POST_CEE);                                            // print (log) detailed output
    }
    else if (m_Star1->IsSNevent() || m_Star2->IsSNevent()) {
        EvaluateSupernovae();                                                                                           // evaluate supernovae (both stars) - immediate event
        (void)PrintDetailedOutput(m_Id, BSE_DETAILED_RECORD_TYPE::POST_SN);                                             // print (log) detailed output
        if (HasOneOf({ STELLAR_TYPE::NEUTRON_STAR })) {
            (void)PrintPulsarEvolutionParameters(PULSAR_RECORD_TYPE::DEFAULT);                                                                         // print (log) pulsar evolution parameters 
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
            (void)PrintPulsarEvolutionParameters(PULSAR_RECORD_TYPE::DEFAULT);                                                                         // print (log) pulsar evolution parameters 
        }
    }

    // assign new values to "previous" values, for following timestep
    m_EccentricityPrev  = m_Eccentricity;
    m_SemiMajorAxisPrev = m_SemiMajorAxis;

    CalculateEnergyAndAngularMomentum();                                                                                // perform energy and angular momentum calculations

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
 * Evolve the binary a single timestep - timestep is provided    JR: todo: fix this documentation - this is for SSE version
 *
 * Each individual star is aged for the same timestep
 *
 * See AgeOneTimestep() documentation in BaseStar.cpp for details
 *
 *
 * void EvolveOneTimestep(const double p_Dt, const int p_LogFileId)
 *
 * @param   [IN]    p_Dt                        The suggested timestep to evolve
 */
void BaseBinaryStar::EvolveOneTimestep(const double p_Dt) {

    EvolveOneTimestepPreamble(p_Dt);

    m_Star1->AgeOneTimestep(p_Dt, true);    // Age the primary one timestep and switch to the new stellar type if necessary
    m_Star2->AgeOneTimestep(p_Dt, true);    // Age the secondary one timestep and switch to the new stellar type if necessary
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

    EVOLUTION_STATUS evolutionStatus = EVOLUTION_STATUS::CONTINUE;

    if (m_Error != ERROR::NONE) {                                                                                                           // check for error creating binary
        SHOW_ERROR(m_Error);                                                                                                                // no - show error
        return EVOLUTION_STATUS::ERROR;                                                                                                     // return without evolving
    }

    if (HasStarsTouching()) {                                                                                                               // check if stars are touching
        m_Flags.stellarMerger        = true;
        m_Flags.stellarMergerAtBirth = true;
        evolutionStatus              = EVOLUTION_STATUS::STELLAR_MERGER_AT_BIRTH;                                                           // binary components are touching - merger at birth
    }

    (void)PrintDetailedOutput(m_Id, BSE_DETAILED_RECORD_TYPE::INITIAL_STATE);                                                               // print (log) detailed output: this is the initial state of the binary

    if (OPTIONS->PopulationDataPrinting()) {
        SAY("\nGenerating a new binary - " << m_Id);
        SAY("Binary has masses " << m_Star1->Mass() << " & " << m_Star2->Mass() << " Msol");
        SAY("Binary has initial semiMajorAxis " << m_SemiMajorAxis << " AU");
        SAY("RandomSeed " << m_RandomSeed);
    }

    if (evolutionStatus == EVOLUTION_STATUS::CONTINUE) {                                                                                    // continue evolution
        // evolve the current binary up to the maximum evolution time (and number of steps)
        double dt      = std::min(m_Star1->CalculateTimestep(), m_Star2->CalculateTimestep()) / 1000.0;                                     // initialise the timestep
        int    stepNum = 1;                                                                                                                 // initialise step number
        while (evolutionStatus == EVOLUTION_STATUS::CONTINUE) {                                                                             // perform binary evolution - iterate over timesteps until told to stop

            EvolveOneTimestep(dt);                                                                                                          // evolve the binary system one timestep

            // check for problems
            if (m_Error != ERROR::NONE) {                                                                                                   // SSE error for either constituent star?
                evolutionStatus = EVOLUTION_STATUS::SSE_ERROR;                                                                              // yes - stop evolution
            }
            else if (HasOneOf({ STELLAR_TYPE::MASSLESS_REMNANT })) {                                                                        // at least one massless remnant?
                evolutionStatus = EVOLUTION_STATUS::MASSLESS_REMNANT;                                                                       // yes - stop evolution
            }
            else if (StellarMerger()) {                                                                                                     // have stars merged?
                evolutionStatus = EVOLUTION_STATUS::STELLAR_MERGER;                                                                         // for now, stop evolution
            }
            else if (HasStarsTouching()) {                                                                                                  // binary components touching? (should usually be avoided as MT or CE or merger should happen prior to this)
                evolutionStatus = EVOLUTION_STATUS::STARS_TOUCHING;                                                                         // yes - stop evolution
            }
            else if (IsUnbound() && !OPTIONS->EvolveUnboundSystems()) {                                                                     // binary is unbound and we don't want unbound systems?
                m_Unbound       = true;                                                                                                     // yes - set the unbound flag (should already be set)
                evolutionStatus = EVOLUTION_STATUS::UNBOUND;                                                                                // stop evolution
            }
            else {                                                                                                                          // continue evolution

                (void)PrintDetailedOutput(m_Id, BSE_DETAILED_RECORD_TYPE::POST_STELLAR_TIMESTEP);                                           // print (log) detailed output

                if (OPTIONS->RLOFPrinting()) StashRLOFProperties(MASS_TRANSFER_TIMING::PRE_MT);                                             // stash properties immediately pre-Mass Transfer 

                EvaluateBinary(dt);                                                                                                         // evaluate the binary at this timestep

                (void)PrintDetailedOutput(m_Id, BSE_DETAILED_RECORD_TYPE::POST_BINARY_TIMESTEP);                                            // print (log) detailed output
                
                (void)PrintRLOFParameters();                                                                                                // print (log) RLOF parameters
                
                // check for problems
                if (StellarMerger()) {                                                                                                      // have stars merged?
                    evolutionStatus = EVOLUTION_STATUS::STELLAR_MERGER;                                                                     // for now, stop evolution
                }
                else if (HasStarsTouching()) {                                                                                              // binary components touching? (should usually be avoided as MT or CE or merger should happen prior to this)
                    evolutionStatus = EVOLUTION_STATUS::STARS_TOUCHING;                                                                     // yes - stop evolution
                }
                else if (IsUnbound()) {                                                                                                     // binary is unbound?
                    m_Flags.mergesInHubbleTime = false;                                                                                     // yes - won't merge in a Hubble time

                    if (IsDCO()) {                                                                                                          // DCO (has two COs)?
                        if (m_DCOFormationTime == DEFAULT_INITIAL_DOUBLE_VALUE) {                                                           // DCO not yet evaluated
                            m_DCOFormationTime = m_Time;                                                                                    // set the DCO formation time
                        }
                    }

                    if (!OPTIONS->EvolveUnboundSystems() || IsDCO()) {                                                                      // should we evolve unbound systems?
                        evolutionStatus = EVOLUTION_STATUS::UNBOUND;                                                                        // no - stop evolution
                    }
                }
                
                if (m_Error != ERROR::NONE) evolutionStatus = EVOLUTION_STATUS::BINARY_ERROR;                                               // error in binary evolution

                if (evolutionStatus == EVOLUTION_STATUS::CONTINUE) {                                                                        // continue evolution?

                    if (HasOneOf({ STELLAR_TYPE::NEUTRON_STAR })) {
                        (void)PrintPulsarEvolutionParameters(PULSAR_RECORD_TYPE::POST_BINARY_TIMESTEP);                                                                             // print (log) pulsar evolution parameters 
                    }

                    //(void)PrintBeBinary();                                                                                                  // print (log) BeBinary properties
                        
                    if (IsDCO() && !IsUnbound()) {                                                                                          // bound double compact object?
                        if (m_DCOFormationTime == DEFAULT_INITIAL_DOUBLE_VALUE) {                                                           // DCO not yet evaluated -- to ensure that the coalescence is only resolved once
                            ResolveCoalescence();                                                                                           // yes - resolve coalescence
                            m_DCOFormationTime = m_Time;                                                                                    // set the DCO formation time
                        }

                        if (!(OPTIONS->EvolvePulsars() && HasOneOf({ STELLAR_TYPE::NEUTRON_STAR }))) {
                            if (!OPTIONS->Quiet()) SAY(ERR_MSG(ERROR::BINARY_EVOLUTION_STOPPED) << ": Double compact object");              // announce that we're stopping evolution
                            evolutionStatus = EVOLUTION_STATUS::STOPPED;                                                                    // stop evolving
                        }
                    }

                    // check for problems
                    if (evolutionStatus == EVOLUTION_STATUS::CONTINUE) {                                                                    // continue evolution?
                        if (m_Error != ERROR::NONE) {                                                                                       // error in binary evolution?
                            evolutionStatus = EVOLUTION_STATUS::BINARY_ERROR;                                                               // yes - stop evolution
                        }
                        else if (!OPTIONS->EvolveDoubleWhiteDwarfs() && IsWDandWD()) {                                                      // double WD and their evolution is not enabled?
                            evolutionStatus = EVOLUTION_STATUS::WD_WD;                                                                      // yes - do not evolve double WD systems
                        }
                        else if (IsDCO() && m_Time > (m_DCOFormationTime + m_TimeToCoalescence) && !IsUnbound()) {                          // evolution time exceeds DCO merger time?
                            evolutionStatus = EVOLUTION_STATUS::STOPPED;                                                                    // yes - stop evolution
                        }
                        else if (m_Time > OPTIONS->MaxEvolutionTime()) {                                                                    // evolution time exceeds maximum?
                            evolutionStatus = EVOLUTION_STATUS::TIMES_UP;                                                                   // yes - stop evolution
                        }
                    }
                }
            }

            (void)PrintDetailedOutput(m_Id, BSE_DETAILED_RECORD_TYPE::TIMESTEP_COMPLETED);                                                  // print (log) detailed output: this is after all changes made in the timestep

            if (stepNum >= OPTIONS->MaxNumberOfTimestepIterations()) evolutionStatus = EVOLUTION_STATUS::STEPS_UP;                          // number of timesteps for evolution exceeds maximum

            if (evolutionStatus == EVOLUTION_STATUS::CONTINUE) {                                                                            // continue evolution?

                dt = std::min(m_Star1->CalculateTimestep(), m_Star2->CalculateTimestep()) * OPTIONS->TimestepMultiplier();                  // new timestep
                if ((m_Star1->IsOneOf({ STELLAR_TYPE::MASSLESS_REMNANT }) || m_Star2->IsOneOf({ STELLAR_TYPE::MASSLESS_REMNANT })) || dt < NUCLEAR_MINIMUM_TIMESTEP) {
                    dt = NUCLEAR_MINIMUM_TIMESTEP;                                                                                          // but not less than minimum
		        }
                stepNum++;                                                                                                                  // increment stepNum
            }
        }

        if (evolutionStatus == EVOLUTION_STATUS::STEPS_UP) {                                                                                // stopped because max timesteps reached?
            SHOW_ERROR(ERROR::BINARY_EVOLUTION_STOPPED);                                                                                    // show error
        }
    }

    (void)PrintDetailedOutput(m_Id, BSE_DETAILED_RECORD_TYPE::FINAL_STATE);                                                                 // print (log) detailed output: this is the final state of the binary

    (void)PrintBinarySystemParameters();                                                                                                    // print (log) binary system parameters

    return evolutionStatus;
}

