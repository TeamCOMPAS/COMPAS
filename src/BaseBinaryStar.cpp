#include "BaseBinaryStar.h"

// gsl includes
#include <gsl/gsl_poly.h>


/* Constructor
 *
 * Parameter p_Id is optional, and is only included so that comparison tests can
 * be run against the legacy Compas code.  If a fixed random seed is being used
 * (program option) the legacy code effectivley adds the loop index of the binary
 * (from COMPASBinary() in main.cpp) to the user-specified fixed random seed so
 * that each binary has a repeatable random seed.
 *
 * Notes: the legacy code doesn't actually use the loop index - it uses a generated
 * object id that is the same as the loop index.  The new code also assigns objects
 * object ids, but the ids are assigned to all objects, not just binary stars, so
 * the ids generated in the new code won't match the legacy code ids - hence the
 * need to use the loop index here.  The parameter is optional - if no comparison
 * testing against the legacy code is required, the p_Id parameter can be let default
 * (in which case it is not used to generate the random seed - the generated object
 * id is used instead).
 */


// binary is generated according to distributions specified in program options
BaseBinaryStar::BaseBinaryStar(const AIS &p_AIS, const long int p_Id) {

    SetInitialCommonValues(p_AIS, p_Id);                                                                                                        // start construction of the binary

    m_CEDetails.alpha = OPTIONS->SampleCommonEnvelopeAlpha()
                        ? RAND->Random(OPTIONS->SampleCommonEnvelopeAlphaMin(), OPTIONS->SampleCommonEnvelopeAlphaMax())
                        : OPTIONS->CommonEnvelopeAlpha();

    m_LBVfactor       = OPTIONS->SampleLuminousBlueVariableMultiplier()
                        ? RAND->Random(OPTIONS->SampleLuminousBlueVariableMultiplierMin(), OPTIONS->SampleLuminousBlueVariableMultiplierMax())
                        : OPTIONS->LuminousBlueVariableFactor();

    m_WolfRayetFactor = OPTIONS->SampleWolfRayetMultiplier()
                        ? RAND->Random(OPTIONS->SampleWolfRayetMultiplierMin(), OPTIONS->SampleWolfRayetMultiplierMax())
                        : OPTIONS->WolfRayetFactor();


    // generate initial properties of binary
    // check that the constituent stars are not touching
    // also check m2 > m2min
    // also check that when we are using AIS we are sampling inside the parameter space

    bool merger                                 = false;
    bool rlof                                   = false;
    bool secondarySmallerThanMinimumMass        = false;
    bool initialParametersOutsideParameterSpace = false;

    do {

        if(OPTIONS->AIS_RefinementPhase()) {                                                                                                    // JR: todo: Floor, do we need to do this inside the loop?
            m_AIS.Initialise();                                                                                                                 // run AIS step 2 and sample from importance sampling distribution
        }

        double mass1        = SampleInitialMassDistribution();
        double massRatio    = SampleQDistribution();
        double mass2        = massRatio * mass1;

        double metallicity1 = std::min(std::max(SampleMetallicityDistribution(), 0.0), 1.0);
        double metallicity2 = std::min(std::max(SampleMetallicityDistribution(), 0.0), 1.0);

        m_SemiMajorAxis     = SampleSemiMajorAxisDistribution(mass1, mass2);
        m_Eccentricity      = SampleEccentricityDistribution();

        // binary star contains two instances of star to hold masses, radii and luminosities.
        // star 1 initially more massive
        m_Star1 = new BinaryConstituentStar(m_RandomSeed, mass1, metallicity1, m_LBVfactor, m_WolfRayetFactor);
        m_Star2 = new BinaryConstituentStar(m_RandomSeed, mass2, metallicity2, m_LBVfactor, m_WolfRayetFactor);

        double factor            = m_SemiMajorAxis * (1.0 - m_Eccentricity);
        double rocheLobeTracker1 = (m_Star1->Radius() * RSOL_TO_AU) / (factor * CalculateRocheLobeRadius_Static(mass1, mass2));
        double rocheLobeTracker2 = (m_Star2->Radius() * RSOL_TO_AU) / (factor * CalculateRocheLobeRadius_Static(mass2, mass1));

        m_MassesEquilibrated = false;                                                                                                           // default

        if ((OPTIONS->CHE_Option() != CHE_OPTION::NONE || OPTIONS->AllowRLOFAtBirth()) &&                                                       // CHE enabled or over-contact binaries at birth allowed?
           (utils::Compare(rocheLobeTracker1, 1.0) > 0 || utils::Compare(rocheLobeTracker2, 1.0) > 0)) {                                        // either star overflowing Roche Lobe?
            rlof                 = false;                                                                                                       // over-contact at birth allowed - set this false
            m_MassesEquilibrated = true;                                                                                                        // record that we've equilbrated

            mass1                = (mass1 + mass2) / 2.0;                                                                                       // equilibrate masses
            mass2                = mass1;                                                                                                       // ditto
            m_SemiMajorAxis     *= (1.0 - (m_Eccentricity * m_Eccentricity));                                                                   // circularise; conserve angular momentum
            m_Eccentricity       = 0.0;                                                                                                         // now circular

            // create new stars with equal masses - eveything else is recalculated
            delete m_Star1;
            m_Star1 = new BinaryConstituentStar(m_RandomSeed, mass1, metallicity1, m_LBVfactor, m_WolfRayetFactor);
            delete m_Star2;
            m_Star2 = new BinaryConstituentStar(m_RandomSeed, mass2, metallicity2, m_LBVfactor, m_WolfRayetFactor);
        
            rocheLobeTracker1 = (m_Star1->Radius() * RSOL_TO_AU) / (m_SemiMajorAxis * CalculateRocheLobeRadius_Static(mass1, mass2));
            rocheLobeTracker2 = (m_Star2->Radius() * RSOL_TO_AU) / (m_SemiMajorAxis * CalculateRocheLobeRadius_Static(mass2, mass1));
        }

        rlof = utils::Compare(rocheLobeTracker1, 1.0) > 0 || utils::Compare(rocheLobeTracker2, 1.0) > 0;

        m_Star1->SetCompanion(m_Star2);
        m_Star2->SetCompanion(m_Star1);

        merger                                 = (m_SemiMajorAxis * AU_TO_RSOL) < (m_Star1->Radius() + m_Star2->Radius());
        secondarySmallerThanMinimumMass        = utils::Compare(mass2, OPTIONS->MinimumMassSecondary()) < 0;                                    // JR: todo: original code missed the == case - assume == is ok
        initialParametersOutsideParameterSpace = false;

        if(OPTIONS->AIS_RefinementPhase()) {                                                                                                    // when using Adaptive Importance Sampling (step 2) check if drawns from Gaussians are inside the COMPAS parameter space
            initialParametersOutsideParameterSpace = utils::Compare(mass1,           OPTIONS->InitialMassFunctionMin())       < 0 ||            // mass1 is outside (below) parameter space
                                                     utils::Compare(mass1,           OPTIONS->InitialMassFunctionMax())       > 0 ||            // mass1 is outside (above) parameter space
                                                     utils::Compare(massRatio,       OPTIONS->MassRatioDistributionMin())     < 0 ||            // massRatio is outside (below) parameter space
                                                     utils::Compare(massRatio,       OPTIONS->MassRatioDistributionMax())     > 0 ||            // massRatio is outside (above) parameter space
                                                     utils::Compare(m_SemiMajorAxis, OPTIONS->SemiMajorAxisDistributionMin()) < 0 ||            // semiMajorAxis is outside (below) parameter space
                                                     utils::Compare(m_SemiMajorAxis, OPTIONS->SemiMajorAxisDistributionMax()) > 0;              // semiMajorAxis is outside (above) parameter space
        }
    } while (rlof || (!OPTIONS->AllowTouchingAtBirth() && merger) || secondarySmallerThanMinimumMass || initialParametersOutsideParameterSpace);

    SetRemainingCommonValues(p_Id);                                                                                                             // complete the construction of the binary
}


// binary is generated according to parameters passed
BaseBinaryStar::BaseBinaryStar(const AIS     &p_AIS,
                               const double   p_Mass1,
                               const double   p_Mass2,
                               const double   p_Metallicity1,
                               const double   p_Metallicity2,
                               const double   p_SemiMajorAxis,
                               const double   p_Eccentricity,
                               const long int p_Id) {

    SetInitialCommonValues(p_AIS, p_Id);                                                                                                        // start construction of the binary

    double mass1 = p_Mass1;                                                                                                                     // specified mass of the primary
    double mass2 = p_Mass2;                                                                                                                     // specified mass of the secondary

    double metallicity1 = std::min(std::max(p_Metallicity1, 0.0), 1.0);                                                                         // specified metallicity of the primary
    double metallicity2 = std::min(std::max(p_Metallicity2, 0.0), 1.0);                                                                         // specified metallicity of the secondary

    m_SemiMajorAxis = p_SemiMajorAxis;                                                                                                          // specified separation
    m_Eccentricity  = p_Eccentricity;                                                                                                           // specified eccentricity

    m_CEDetails.alpha = OPTIONS->CommonEnvelopeAlpha();
    m_LBVfactor       = OPTIONS->LuminousBlueVariableFactor();
    m_WolfRayetFactor = OPTIONS->WolfRayetFactor();

    // binary star contains two instances of star to hold masses, radii and luminosities.
    // star 1 initially more massive (JR: todo: this is not guaranteed...)
    m_Star1 = new BinaryConstituentStar(m_RandomSeed, mass1, metallicity1, m_LBVfactor, m_WolfRayetFactor);
    m_Star2 = new BinaryConstituentStar(m_RandomSeed, mass2, metallicity2, m_LBVfactor, m_WolfRayetFactor);

    m_Star1->SetCompanion(m_Star2);
    m_Star2->SetCompanion(m_Star1);

    double factor            = m_SemiMajorAxis * (1.0 - m_Eccentricity);
    double rocheLobeTracker1 = (m_Star1->Radius() * RSOL_TO_AU) / (factor * CalculateRocheLobeRadius_Static(mass1, mass2));
    double rocheLobeTracker2 = (m_Star2->Radius() * RSOL_TO_AU) / (factor * CalculateRocheLobeRadius_Static(mass2, mass1));

    m_MassesEquilibrated = false;                                                                                                               // default

    if ((OPTIONS->CHE_Option() != CHE_OPTION::NONE || OPTIONS->AllowRLOFAtBirth()) &&                                                           // CHE enabled or over-contact binaries at birth allowed?
       (utils::Compare(rocheLobeTracker1, 1.0) > 0 || utils::Compare(rocheLobeTracker2, 1.0) > 0)) {                                            // either star overflowing Roche Lobe?

        m_MassesEquilibrated = true;                                                                                                            // record that we've equilbrated

        double newMass1      = (mass1 + mass2) / 2.0;                                                                                           // equilibrate masses
        double newMass2      = newMass1;                                                                                                        // ditto
        m_SemiMajorAxis     *= (1.0 - (m_Eccentricity * m_Eccentricity));                                                                       // circularise; conserve angular momentum
        m_Eccentricity       = 0.0;                                                                                                             // now circular

        // equilibrate masses - recalculate everything else
        m_Star1->UpdateAttributesAndAgeOneTimestep(newMass1 - mass1, newMass1 - mass1, 0.0, true);
        m_Star2->UpdateAttributesAndAgeOneTimestep(newMass2 - mass2, newMass2 - mass2, 0.0, true);
    }

    SetRemainingCommonValues(p_Id);                                                                                                             // complete the construction of the binary
}


/*
 * Initiate the construction of the binary - initial common values
 *
 *
 * void SetInitialCommonValues(const AIS &p_AIS, const long int p_Id)
 *
 * @param   [IN]    p_AIS                       AIS object passed to the constructor
 * @param   [IN]    p_Id                        Ordinal value of binary - see constructor notes above
 */
void BaseBinaryStar::SetInitialCommonValues(const AIS &p_AIS, const long int p_Id) {

    m_Error = ERROR::NONE;

    m_ObjectId    = globalObjectId++;
    m_ObjectType  = OBJECT_TYPE::BASE_BINARY_STAR;
    m_StellarType = STELLAR_TYPE::BINARY_STAR;


    // binary stars generate their own random seed, and pass that to their constituent stars

    OBJECT_ID id = p_Id < 0 ? m_ObjectId : p_Id;                                                                                                // for legacy testing

    if (OPTIONS->FixedRandomSeed()) {                                                                                                           // user supplied seed for the random number generator?

        m_RandomSeed = RAND->Seed(OPTIONS->RandomSeed() + id);                                                                                  // yes - this allows the user to reproduce results for each binary

        if (OPTIONS->PopulationDataPrinting()) {                                                                                                // JR: todo: what is the aim of PopulationDataPrinting?
            SAY("Using supplied random seed " << m_RandomSeed << " for Binary Star id = " << m_ObjectId);
        }
    }
    else {                                                                                                                                      // no

        m_RandomSeed = RAND->Seed(RAND->DefaultSeed() + id);                                                                                    // use default seed (based on system time) + id

        if (OPTIONS->PopulationDataPrinting()) {                                                                                                // JR: todo: what is the aim of PopulationDataPrinting?
            SAY("Using default random seed " << m_RandomSeed << " for Binary Star id = " << m_ObjectId);
        }
    }

    m_AIS = p_AIS;                                                                                                                              // Adaptive Importance Sampling
}


/*
 * Complete the construction of the binary - remaining common values
 *
 *
 * void SetRemainingCommonValues(const long int p_Id)
 *
 * @param   [IN]    p_Id                        Ordinal value of binary - see constructor notes above
 */
void BaseBinaryStar::SetRemainingCommonValues(const long int p_Id) {

    // Initialise other parameters

    m_SemiMajorAxisPrime          = m_SemiMajorAxis;
    m_SemiMajorAxisPrev           = m_SemiMajorAxis;

    m_EccentricityPrime           = m_Eccentricity;
    m_EccentricityPrev            = m_Eccentricity;

    // initial binary parameters - kept constant as a record of the initial parameters of the binary
    m_SemiMajorAxisInitial        = m_SemiMajorAxis;
    m_EccentricityInitial         = m_Eccentricity;

    // initialise variables to hold parameters prior to 2nd supernova
    m_SemiMajorAxisPre2ndSN       = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_EccentricityPre2ndSN        = DEFAULT_INITIAL_DOUBLE_VALUE;

    // initialise variables to hold parameters at DCO formation
    m_SemiMajorAxisAtDCOFormation = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_EccentricityAtDCOFormation  = DEFAULT_INITIAL_DOUBLE_VALUE;


    m_OrbitalVelocity             = sqrt(G1 * (m_Star1->Mass() + m_Star2->Mass()) / (m_SemiMajorAxis * m_SemiMajorAxis * m_SemiMajorAxis)); // rads/year

    m_OrbitalVelocityPrime        = m_OrbitalVelocity;
    m_OrbitalVelocityPrev         = m_OrbitalVelocity;
    m_OrbitalVelocityPre2ndSN     = DEFAULT_INITIAL_DOUBLE_VALUE;


    // if CHE enabled, update rotational frequency for constituent stars - assume tidally locked
    if (OPTIONS->CHE_Option() != CHE_OPTION::NONE) {

        m_Star1->SetOmega(m_OrbitalVelocity);
        m_Star2->SetOmega(m_OrbitalVelocity);

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
        if (utils::Compare(m_OrbitalVelocity, m_Star1->OmegaCHE()) >= 0) {                                                                      // star 1 CH?
            if (m_Star1->StellarType() != STELLAR_TYPE::CHEMICALLY_HOMOGENEOUS) m_Star1->SwitchTo(STELLAR_TYPE::CHEMICALLY_HOMOGENEOUS, true);  // yes, switch if not alread Chemically Homogeneous
        }
        else if (m_Star1->MZAMS() <= 0.7) {                                                                                                     // no - MS - initial mass determines actual type  JR: don't use utils::Compare() here
            if (m_Star1->StellarType() != STELLAR_TYPE::MS_LTE_07) m_Star1->SwitchTo(STELLAR_TYPE::MS_LTE_07, true);                            // MS <= 0.7 Msol - switch if necessary
        }
        else {
            if (m_Star1->StellarType() != STELLAR_TYPE::MS_GT_07) m_Star1->SwitchTo(STELLAR_TYPE::MS_GT_07, true);                              // MS > 0.7 Msol - switch if necessary
        }

        // star 2
        if (utils::Compare(m_OrbitalVelocity, m_Star2->OmegaCHE()) >= 0) {                                                                      // star 2 CH?
            if (m_Star2->StellarType() != STELLAR_TYPE::CHEMICALLY_HOMOGENEOUS) m_Star2->SwitchTo(STELLAR_TYPE::CHEMICALLY_HOMOGENEOUS, true);  // yes, switch if not alread Chemically Homogeneous
        }
        else if (m_Star2->MZAMS() <= 0.7) {                                                                                                     // no - MS - initial mass determines actual type  JR: don't use utils::Compare() here
            if (m_Star2->StellarType() != STELLAR_TYPE::MS_LTE_07) m_Star2->SwitchTo(STELLAR_TYPE::MS_LTE_07, true);                            // MS <= 0.0 Msol - switch if necessary
        }
        else {
            if (m_Star2->StellarType() != STELLAR_TYPE::MS_GT_07) m_Star2->SwitchTo(STELLAR_TYPE::MS_GT_07, true);                              // MS > 0.7 Msol - switch if necessary
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
                                                                        m_OrbitalVelocityPrime,
                                                                        gyrationRadius1,
                                                                        gyrationRadius2);

    m_TotalEnergyPrime                           = m_TotalEnergy;

    m_TotalAngularMomentumPrime                  = CalculateAngularMomentum(m_SemiMajorAxis,
                                                                            m_Eccentricity,
                                                                            m_Star1->Mass(),
                                                                            m_Star2->Mass(),
                                                                            m_Star1->RZAMS(),
                                                                            m_Star2->RZAMS(),
                                                                            m_Star1->Omega(),
                                                                            m_Star2->Omega(),
                                                                            gyrationRadius1,
                                                                            gyrationRadius2);

    m_TotalAngularMomentumPrev                   = m_TotalAngularMomentumPrime;
	m_TotalMassPrime 					         = m_Star1->Mass() + m_Star2->Mass();
	m_TotalMassPrev						         = m_TotalMassPrime;
	m_ReducedMassPrime					         = (m_Star1->Mass() * m_Star2->Mass()) / m_TotalMassPrime;
	m_ReducedMassPrev					         = m_ReducedMassPrime;
	m_TotalOrbitalEnergyPrime 			         = CalculateOrbitalEnergy(m_ReducedMassPrime, m_TotalMassPrime, m_SemiMajorAxisPrime);
	m_TotalOrbitalEnergyPrev 			         = m_TotalOrbitalEnergyPrime;

	m_TotalOrbitalAngularMomentumPrime 	         = CalculateOrbitalAngularMomentum(m_ReducedMassPrime, m_TotalMassPrime, m_SemiMajorAxisPrime);
	m_TotalOrbitalAngularMomentumPrev 	         = m_TotalOrbitalAngularMomentumPrime;

    m_Time                                       = DEFAULT_INITIAL_DOUBLE_VALUE;
	m_Dt                                         = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_TimePrev                                   = DEFAULT_INITIAL_DOUBLE_VALUE;

    m_SecondaryTooSmallForDCO                    = false;

    // Differential quantities
    m_aTidesDiff                                 = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_OmegaTidesDiff                             = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_OmegaTides                                 = DEFAULT_INITIAL_DOUBLE_VALUE;

    m_aMassLossDiff                              = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_OmegaMassLossDiff                          = DEFAULT_INITIAL_DOUBLE_VALUE;

    m_aMassTransferDiff                          = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_OmegaMassTransferDiff                      = DEFAULT_INITIAL_DOUBLE_VALUE;

	m_MassTransferTrackerHistory                 = MT_TRACKING::NO_MASS_TRANSFER;
    m_MassTransfer                               = false;

    m_JLoss                                      = DEFAULT_INITIAL_DOUBLE_VALUE;

	m_FractionAccreted                           = DEFAULT_INITIAL_DOUBLE_VALUE;

    m_FastPhaseCaseA                             = false;

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

    m_StellarMerger                              = false;
    m_StellarMergerAtBirth                       = false;

	m_Mass1Final                                 = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_Mass2Final                                 = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_MassEnv1                                   = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_MassEnv2                                   = DEFAULT_INITIAL_DOUBLE_VALUE;

    m_ZetaRLOFAnalytic                           = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_ZetaRLOFNumerical                          = DEFAULT_INITIAL_DOUBLE_VALUE;
	m_ZetaStarCompare	                         = DEFAULT_INITIAL_DOUBLE_VALUE;

    // Misalignments
    m_Theta1_i                                   = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_Theta2_i                                   = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_Theta1                                     = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_Theta2                                     = DEFAULT_INITIAL_DOUBLE_VALUE;

    // Initialise other parameters to 0
    m_MSN                                        = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_MSNPrime                                   = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_MC                                         = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_MCPrime                                    = DEFAULT_INITIAL_DOUBLE_VALUE;

    m_VRel                                       = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_uK                                         = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_Radius                                     = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_EPrime                                     = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_CosIPrime                                  = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_IPrime                                     = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_TimeToCoalescence                          = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_Beta                                       = DEFAULT_INITIAL_DOUBLE_VALUE;

    m_SupernovaState                             = SN_STATE::NONE;

    m_Survived                                   = true;
    m_Merged                                     = false;
    m_MergesInHubbleTime                         = false;
    m_Disbound                                   = false;

    m_SystemicVelocity                           = DEFAULT_INITIAL_DOUBLE_VALUE;

	m_SynchronizationTimescale                   = DEFAULT_INITIAL_DOUBLE_VALUE;
	m_CircularizationTimescale                   = DEFAULT_INITIAL_DOUBLE_VALUE;

	// RLOF details
    m_RLOFDetails.experiencedRLOF                = false;
    m_RLOFDetails.immediateRLOFPostCEE           = false;
    m_RLOFDetails.isRLOF                         = false;
    m_RLOFDetails.simultaneousRLOF               = false;
    m_RLOFDetails.stableRLOFPostCEE              = false;

	// RLOF details - current properties
    m_RLOFDetails.currentProps.id                = -1l;
    m_RLOFDetails.currentProps.randomSeed        = DEFAULT_INITIAL_ULONGINT_VALUE;

    m_RLOFDetails.currentProps.stellarType1      = STELLAR_TYPE::NONE;
    m_RLOFDetails.currentProps.stellarType2      = STELLAR_TYPE::NONE;

    m_RLOFDetails.currentProps.mass1             = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_RLOFDetails.currentProps.mass2             = DEFAULT_INITIAL_DOUBLE_VALUE;

    m_RLOFDetails.currentProps.radius1           = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_RLOFDetails.currentProps.radius2           = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_RLOFDetails.currentProps.separation        = DEFAULT_INITIAL_DOUBLE_VALUE;

    m_RLOFDetails.currentProps.eventCounter      = DEFAULT_INITIAL_ULONGINT_VALUE;

    m_RLOFDetails.currentProps.time              = DEFAULT_INITIAL_DOUBLE_VALUE;

    m_RLOFDetails.currentProps.isRLOF1           = false;
    m_RLOFDetails.currentProps.isRLOF2           = false;

    m_RLOFDetails.currentProps.isCE              = false;

    m_RLOFDetails.currentProps.monitorMS1        = false;
    m_RLOFDetails.currentProps.monitorMS2        = false;
    m_RLOFDetails.currentProps.monitorHeMS1      = false;
    m_RLOFDetails.currentProps.monitorHeMS2      = false;

	// RLOF details - previous properties
    m_RLOFDetails.previousProps.id = -1l;
    m_RLOFDetails.previousProps.randomSeed       = DEFAULT_INITIAL_ULONGINT_VALUE;

    m_RLOFDetails.previousProps.stellarType1     = STELLAR_TYPE::NONE;
    m_RLOFDetails.previousProps.stellarType2     = STELLAR_TYPE::NONE;

    m_RLOFDetails.previousProps.mass1            = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_RLOFDetails.previousProps.mass2            = DEFAULT_INITIAL_DOUBLE_VALUE;

    m_RLOFDetails.previousProps.radius1          = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_RLOFDetails.previousProps.radius2          = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_RLOFDetails.previousProps.separation       = DEFAULT_INITIAL_DOUBLE_VALUE;

    m_RLOFDetails.previousProps.eventCounter     = DEFAULT_INITIAL_ULONGINT_VALUE;

    m_RLOFDetails.previousProps.time             = DEFAULT_INITIAL_DOUBLE_VALUE;

    m_RLOFDetails.previousProps.isRLOF1          = false;
    m_RLOFDetails.previousProps.isRLOF2          = false;

    m_RLOFDetails.previousProps.isCE             = false;

    m_RLOFDetails.previousProps.monitorMS1       = false;
    m_RLOFDetails.previousProps.monitorMS2       = false;
    m_RLOFDetails.previousProps.monitorHeMS1     = false;
    m_RLOFDetails.previousProps.monitorHeMS2     = false;

    // BeBinary details - properties 1
    m_BeBinaryDetails.props1.id                  = -1l;
    m_BeBinaryDetails.props1.randomSeed          = DEFAULT_INITIAL_ULONGINT_VALUE;

    m_BeBinaryDetails.props1.dt                  = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_BeBinaryDetails.props1.totalTime           = DEFAULT_INITIAL_DOUBLE_VALUE;

    m_BeBinaryDetails.props1.massNS              = DEFAULT_INITIAL_DOUBLE_VALUE;

    m_BeBinaryDetails.props1.companionMass       = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_BeBinaryDetails.props1.companionLuminosity = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_BeBinaryDetails.props1.companionTeff       = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_BeBinaryDetails.props1.companionRadius     = DEFAULT_INITIAL_DOUBLE_VALUE;

    m_BeBinaryDetails.props1.separation          = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_BeBinaryDetails.props1.eccentricity        = DEFAULT_INITIAL_DOUBLE_VALUE;

    // BeBinary details - properties 2
    m_BeBinaryDetails.props2.id                  = -1l;
    m_BeBinaryDetails.props2.randomSeed          = DEFAULT_INITIAL_ULONGINT_VALUE;

    m_BeBinaryDetails.props2.dt                  = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_BeBinaryDetails.props2.totalTime           = DEFAULT_INITIAL_DOUBLE_VALUE;

    m_BeBinaryDetails.props2.massNS              = DEFAULT_INITIAL_DOUBLE_VALUE;

    m_BeBinaryDetails.props2.companionMass       = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_BeBinaryDetails.props2.companionLuminosity = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_BeBinaryDetails.props2.companionTeff       = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_BeBinaryDetails.props2.companionRadius     = DEFAULT_INITIAL_DOUBLE_VALUE;

    m_BeBinaryDetails.props2.separation          = DEFAULT_INITIAL_DOUBLE_VALUE;
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
        case BINARY_PROPERTY::BE_BINARY_CURRENT_COMPANION_TEFF:                     value = BeBinaryDetails().currentProps->companionTeff;                      break;
        case BINARY_PROPERTY::BE_BINARY_CURRENT_DT:                                 value = BeBinaryDetails().currentProps->dt;                                 break;
        case BINARY_PROPERTY::BE_BINARY_CURRENT_ECCENTRICITY:                       value = BeBinaryDetails().currentProps->eccentricity;                       break;
        case BINARY_PROPERTY::BE_BINARY_CURRENT_ID:                                 value = BeBinaryDetails().currentProps->id;                                 break;
        case BINARY_PROPERTY::BE_BINARY_CURRENT_NS_MASS:                            value = BeBinaryDetails().currentProps->massNS;                             break;
        case BINARY_PROPERTY::BE_BINARY_CURRENT_RANDOM_SEED:                        value = BeBinaryDetails().currentProps->randomSeed;                         break;
        case BINARY_PROPERTY::BE_BINARY_CURRENT_SEPARATION:                         value = BeBinaryDetails().currentProps->separation;                         break;
        case BINARY_PROPERTY::BE_BINARY_CURRENT_TOTAL_TIME:                         value = BeBinaryDetails().currentProps->totalTime;                          break;
        case BINARY_PROPERTY::CIRCULARIZATION_TIMESCALE:                            value = CircularizationTimescale();                                         break;
        case BINARY_PROPERTY::COMMON_ENVELOPE_ALPHA:                                value = CEAlpha();                                                          break;
        case BINARY_PROPERTY::COMMON_ENVELOPE_AT_LEAST_ONCE:                        value = CEAtLeastOnce();                                                    break;
        case BINARY_PROPERTY::COMMON_ENVELOPE_EVENT_COUNT:                          value = CommonEnvelopeEventCount();                                         break;
        case BINARY_PROPERTY::DIMENSIONLESS_KICK_VELOCITY:                          value = UK();                                                               break;
        case BINARY_PROPERTY::DISBOUND:                                             value = Disbound();                                                         break;
        case BINARY_PROPERTY::DOUBLE_CORE_COMMON_ENVELOPE:                          value = DoubleCoreCE();                                                     break;
        case BINARY_PROPERTY::DT:                                                   value = Dt();                                                               break;
        case BINARY_PROPERTY::ECCENTRICITY:                                         value = Eccentricity();                                                     break;
        case BINARY_PROPERTY::ECCENTRICITY_AT_DCO_FORMATION:                        value = EccentricityAtDCOFormation();                                       break;
        case BINARY_PROPERTY::ECCENTRICITY_INITIAL:                                 value = EccentricityInitial();                                              break;
        case BINARY_PROPERTY::ECCENTRICITY_POST_COMMON_ENVELOPE:                    value = EccentricityPostCEE();                                              break;
        case BINARY_PROPERTY::ECCENTRICITY_PRE_2ND_SUPERNOVA:                       value = EccentricityPre2ndSN();                                             break;
        case BINARY_PROPERTY::ECCENTRICITY_PRE_COMMON_ENVELOPE:                     value = EccentricityPreCEE();                                               break;
        case BINARY_PROPERTY::ECCENTRICITY_PRIME:                                   value = EccentricityPrime();                                                break;
        case BINARY_PROPERTY::ERROR:                                                value = Error();                                                            break;
        case BINARY_PROPERTY::ID:                                                   value = ObjectId();                                                         break;
        case BINARY_PROPERTY::IMMEDIATE_RLOF_POST_COMMON_ENVELOPE:                  value = ImmediateRLOFPostCEE();                                             break;
        case BINARY_PROPERTY::LUMINOUS_BLUE_VARIABLE_FACTOR:                        value = LBV_Factor();                                                       break;
        case BINARY_PROPERTY::MASS_1_FINAL:                                         value = Mass1Final();                                                       break;
        case BINARY_PROPERTY::MASS_1_POST_COMMON_ENVELOPE:                          value = Mass1PostCEE();                                                     break;
        case BINARY_PROPERTY::MASS_1_PRE_COMMON_ENVELOPE:                           value = Mass1PreCEE();                                                      break;
        case BINARY_PROPERTY::MASS_2_FINAL:                                         value = Mass2Final();                                                       break;
        case BINARY_PROPERTY::MASS_2_POST_COMMON_ENVELOPE:                          value = Mass2PostCEE();                                                     break;
        case BINARY_PROPERTY::MASS_2_PRE_COMMON_ENVELOPE:                           value = Mass2PreCEE();                                                      break;
        case BINARY_PROPERTY::MASS_ENV_1:                                           value = MassEnv1();                                                         break;
        case BINARY_PROPERTY::MASS_ENV_2:                                           value = MassEnv2();                                                         break;
        case BINARY_PROPERTY::MASSES_EQUILIBRATED:                                  value = MassesEquilibrated();                                               break;
        case BINARY_PROPERTY::MASS_TRANSFER_TRACKER_HISTORY:                        value = MassTransferTrackerHistory();                                       break;
        case BINARY_PROPERTY::MERGES_IN_HUBBLE_TIME:                                value = MergesInHubbleTime();                                               break;
        case BINARY_PROPERTY::OPTIMISTIC_COMMON_ENVELOPE:                           value = OptimisticCommonEnvelope();                                         break;
        case BINARY_PROPERTY::ORBITAL_VELOCITY:                                     value = OrbitalVelocity();                                                  break;
        case BINARY_PROPERTY::ORBITAL_VELOCITY_PRE_2ND_SUPERNOVA:                   value = OrbitalVelocityPre2ndSN();                                          break;
        case BINARY_PROPERTY::RADIUS_1_POST_COMMON_ENVELOPE:                        value = Radius1PostCEE();                                                   break;
        case BINARY_PROPERTY::RADIUS_1_PRE_COMMON_ENVELOPE:                         value = Radius1PreCEE();                                                    break;
        case BINARY_PROPERTY::RADIUS_2_POST_COMMON_ENVELOPE:                        value = Radius2PostCEE();                                                   break;
        case BINARY_PROPERTY::RADIUS_2_PRE_COMMON_ENVELOPE:                         value = Radius2PreCEE();                                                    break;
        case BINARY_PROPERTY::RANDOM_SEED:                                          value = RandomSeed();                                                       break;
        case BINARY_PROPERTY::RLOF_CURRENT_COMMON_ENVELOPE:                         value = RLOFDetails().currentProps.isCE;                                    break;
        case BINARY_PROPERTY::RLOF_CURRENT_EVENT_COUNTER:                           value = RLOFDetails().currentProps.eventCounter;                            break;
        case BINARY_PROPERTY::RLOF_CURRENT_ID:                                      value = RLOFDetails().currentProps.id;                                      break;
        case BINARY_PROPERTY::RLOF_CURRENT_RANDOM_SEED:                             value = RLOFDetails().currentProps.randomSeed;                              break;
        case BINARY_PROPERTY::RLOF_CURRENT_SEPARATION:                              value = RLOFDetails().currentProps.separation;                              break;
        case BINARY_PROPERTY::RLOF_CURRENT_STAR1_MASS:                              value = RLOFDetails().currentProps.mass1;                                   break;
        case BINARY_PROPERTY::RLOF_CURRENT_STAR2_MASS:                              value = RLOFDetails().currentProps.mass2;                                   break;
        case BINARY_PROPERTY::RLOF_CURRENT_STAR1_RADIUS:                            value = RLOFDetails().currentProps.radius1;                                 break;
        case BINARY_PROPERTY::RLOF_CURRENT_STAR2_RADIUS:                            value = RLOFDetails().currentProps.radius2;                                 break;
        case BINARY_PROPERTY::RLOF_CURRENT_STAR1_RLOF:                              value = RLOFDetails().currentProps.isRLOF1;                                 break;
        case BINARY_PROPERTY::RLOF_CURRENT_STAR2_RLOF:                              value = RLOFDetails().currentProps.isRLOF2;                                 break;
        case BINARY_PROPERTY::RLOF_CURRENT_STAR1_STELLAR_TYPE:                      value = RLOFDetails().currentProps.stellarType1;                            break;
        case BINARY_PROPERTY::RLOF_CURRENT_STAR1_STELLAR_TYPE_NAME:                 value = STELLAR_TYPE_LABEL.at(RLOFDetails().currentProps.stellarType1);     break;
        case BINARY_PROPERTY::RLOF_CURRENT_STAR2_STELLAR_TYPE:                      value = RLOFDetails().currentProps.stellarType2;                            break;
        case BINARY_PROPERTY::RLOF_CURRENT_STAR2_STELLAR_TYPE_NAME:                 value = STELLAR_TYPE_LABEL.at(RLOFDetails().currentProps.stellarType2);     break;
        case BINARY_PROPERTY::RLOF_CURRENT_TIME:                                    value = RLOFDetails().currentProps.time;                                    break;
        case BINARY_PROPERTY::RLOF_PREVIOUS_EVENT_COUNTER:                          value = RLOFDetails().previousProps.eventCounter;                           break;
        case BINARY_PROPERTY::RLOF_PREVIOUS_SEPARATION:                             value = RLOFDetails().previousProps.separation;                             break;
        case BINARY_PROPERTY::RLOF_PREVIOUS_STAR1_MASS:                             value = RLOFDetails().previousProps.mass1;                                  break;
        case BINARY_PROPERTY::RLOF_PREVIOUS_STAR2_MASS:                             value = RLOFDetails().previousProps.mass2;                                  break;
        case BINARY_PROPERTY::RLOF_PREVIOUS_STAR1_RADIUS:                           value = RLOFDetails().previousProps.radius1;                                break;
        case BINARY_PROPERTY::RLOF_PREVIOUS_STAR2_RADIUS:                           value = RLOFDetails().previousProps.radius2;                                break;
        case BINARY_PROPERTY::RLOF_PREVIOUS_STAR1_RLOF:                             value = RLOFDetails().previousProps.isRLOF1;                                break;
        case BINARY_PROPERTY::RLOF_PREVIOUS_STAR2_RLOF:                             value = RLOFDetails().previousProps.isRLOF2;                                break;
        case BINARY_PROPERTY::RLOF_PREVIOUS_STAR1_STELLAR_TYPE:                     value = RLOFDetails().previousProps.stellarType1;                           break;
        case BINARY_PROPERTY::RLOF_PREVIOUS_STAR1_STELLAR_TYPE_NAME:                value = STELLAR_TYPE_LABEL.at(RLOFDetails().previousProps.stellarType1);    break;
        case BINARY_PROPERTY::RLOF_PREVIOUS_STAR2_STELLAR_TYPE:                     value = RLOFDetails().previousProps.stellarType2;                           break;
        case BINARY_PROPERTY::RLOF_PREVIOUS_STAR2_STELLAR_TYPE_NAME:                value = STELLAR_TYPE_LABEL.at(RLOFDetails().previousProps.stellarType2);    break;
        case BINARY_PROPERTY::RLOF_PREVIOUS_TIME:                                   value = RLOFDetails().previousProps.time;                                   break;
        case BINARY_PROPERTY::RLOF_SECONDARY_POST_COMMON_ENVELOPE:                  value = RLOFSecondaryPostCEE();                                             break;
        case BINARY_PROPERTY::ROCHE_LOBE_RADIUS_1:                                  value = RocheLobeRadius1();                                                 break;
        case BINARY_PROPERTY::ROCHE_LOBE_RADIUS_1_POST_COMMON_ENVELOPE:             value = RocheLobe1to2PostCEE();                                             break;
        case BINARY_PROPERTY::ROCHE_LOBE_RADIUS_1_PRE_COMMON_ENVELOPE:              value = RocheLobe1to2PreCEE();                                              break;
        case BINARY_PROPERTY::ROCHE_LOBE_RADIUS_2:                                  value = RocheLobeRadius2();                                                 break;
        case BINARY_PROPERTY::ROCHE_LOBE_RADIUS_2_POST_COMMON_ENVELOPE:             value = RocheLobe2to1PostCEE();                                             break;
        case BINARY_PROPERTY::ROCHE_LOBE_RADIUS_2_PRE_COMMON_ENVELOPE:              value = RocheLobe2to1PreCEE();                                              break;
        case BINARY_PROPERTY::ROCHE_LOBE_TRACKER_1:                                 value = RocheLobeTracker1();                                                break;
        case BINARY_PROPERTY::ROCHE_LOBE_TRACKER_2:                                 value = RocheLobeTracker2();                                                break;
        case BINARY_PROPERTY::SECONDARY_TOO_SMALL_FOR_DCO:                          value = SecondaryTooSmallForDCO();                                          break;
        case BINARY_PROPERTY::SEMI_MAJOR_AXIS_AT_DCO_FORMATION:                     value = SemiMajorAxisAtDCOFormation();                                      break;
        case BINARY_PROPERTY::SEMI_MAJOR_AXIS_INITIAL:                              value = SemiMajorAxisInitial();                                             break;
        case BINARY_PROPERTY::SEMI_MAJOR_AXIS_POST_COMMON_ENVELOPE:                 value = SemiMajorAxisPostCEE();                                             break;
        case BINARY_PROPERTY::SEMI_MAJOR_AXIS_PRE_2ND_SUPERNOVA:                    value = SemiMajorAxisPre2ndSN();                                            break;
        case BINARY_PROPERTY::SEMI_MAJOR_AXIS_PRE_2ND_SUPERNOVA_RSOL:               value = SemiMajorAxisPre2ndSN() * AU_TO_RSOL;                               break;
        case BINARY_PROPERTY::SEMI_MAJOR_AXIS_PRE_COMMON_ENVELOPE:                  value = SemiMajorAxisPreCEE();                                              break;
        case BINARY_PROPERTY::SEMI_MAJOR_AXIS_PRIME:                                value = SemiMajorAxisPrime();                                               break;
        case BINARY_PROPERTY::SEMI_MAJOR_AXIS_PRIME_RSOL:                           value = SemiMajorAxisPrime() * AU_TO_RSOL;                                  break;
        case BINARY_PROPERTY::SIMULTANEOUS_RLOF:                                    value = SimultaneousRLOF();                                                 break;
        case BINARY_PROPERTY::STABLE_RLOF_POST_COMMON_ENVELOPE:                     value = StableRLOFPostCEE();                                                break;
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
        case BINARY_PROPERTY::SUPERNOVA_STATE:                                      value = SN_State();                                                         break;
        case BINARY_PROPERTY::SURVIVED_SUPERNOVA_EVENT:                             value = SurvivedSNEvent();                                                  break;
        case BINARY_PROPERTY::SYNCHRONIZATION_TIMESCALE:                            value = SynchronizationTimescale();                                         break;
        case BINARY_PROPERTY::SYSTEMIC_VELOCITY:                                    value = SystemicVelocity();                                                 break;
        case BINARY_PROPERTY::TIME:                                                 value = Time();                                                             break;
        case BINARY_PROPERTY::TIME_TO_COALESCENCE:                                  value = TimeToCoalescence();                                                break;
        case BINARY_PROPERTY::TOTAL_ANGULAR_MOMENTUM_PRIME:                         value = TotalAngularMomentumPrime();                                        break;
        case BINARY_PROPERTY::TOTAL_ENERGY_PRIME:                                   value = TotalEnergyPrime();                                                 break;
        case BINARY_PROPERTY::WOLF_RAYET_FACTOR:                                    value = WolfRayetFactor();                                                  break;
        case BINARY_PROPERTY::ZETA_RLOF_ANALYTIC:                                   value = ZetaRLOFAnalytic();                                                 break;
        case BINARY_PROPERTY::ZETA_RLOF_NUMERICAL:                                  value = ZetaRLOFNumerical();                                                break;
        case BINARY_PROPERTY::ZETA_STAR_COMPARE:                                    value = ZetaStarCompare();                                                  break;

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
 * The bool returned indicates whether the property value was retireved ok: true = yes, fales = no
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
            SHOW_WARN(ERROR::UNKNOWN_PROPERTY_TYPE  );                                                                  // show warning
    }

    return std::make_tuple(ok, value);
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
 * Draw semi-major axis from the distribution specified by the user
 * (SemiMajorAxisDistribution program option; will use AIS distribution if specified (AIS.DrawingFromAISDistributions))
 *
 *
 * double SampleSemiMajorAxisDistribution(const double p_Mass1, const double p_Mass2)
 *
 * @param   [IN]    p_Mass1                     Mass of the primary
 * @param   [IN]    p_Mass1                     Mass of the secondary
 * @return                                      Semi-major axis in AU
 */
double BaseBinaryStar::SampleSemiMajorAxisDistribution(const double p_Mass1, const double p_Mass2) {

    double semiMajorAxis;

    if (!m_AIS.DrawingFromAISDistributions()) {                                                                                                 // draw from priors (not from AIS distributions)

        switch (OPTIONS->SemiMajorAxisDistribution()) {                                                                                         // which distribution?

            case SEMI_MAJOR_AXIS_DISTRIBUTION::FLATINLOG:                                                                                       // FLAT IN LOG

                semiMajorAxis = utils::InverseSampleFromPowerLaw(-1.0, OPTIONS->SemiMajorAxisDistributionMax(), OPTIONS->SemiMajorAxisDistributionMin());
                break;

            case SEMI_MAJOR_AXIS_DISTRIBUTION::DUQUENNOYMAYOR1991:                                                                              // Duquennoy & Mayor (1991) period distribution
                // http://adsabs.harvard.edu/abs/1991A%26A...248..485D
                // See also the period distribution (Figure 1) of M35 in Geller+ 2013 https://arxiv.org/abs/1210.1575
                // See also the period distribution (Figure 13) of local solar type binaries from Raghavan et al 2010 https://arxiv.org/abs/1007.0414
                // They have log-normal distribution with a mean of 5.03 and a standard deviation of 2.28, with a minimum period of around 0.1 days
                // Sampling function taken from binpop.f in NBODY6

                // Make sure that the drawn semi-major axis is in the range specified by the user
                do {                                                                                                                            // JR: todo: catch for non-convergence?
                    double periodInDays = pow(10.0, 2.3 * sqrt(-2.0 * log(RAND->Random())) * cos(2.0 * M_PI * RAND->Random()) + 4.8);
                    semiMajorAxis = utils::ConvertPeriodInDaysToSemiMajorAxisInAU(p_Mass1, p_Mass2, periodInDays);                              // convert period in days to semi-major axis in AU
                } while (semiMajorAxis < OPTIONS->SemiMajorAxisDistributionMin() || semiMajorAxis > OPTIONS->SemiMajorAxisDistributionMax());   // JR: don't use utils::Compare() here
                break;

            case SEMI_MAJOR_AXIS_DISTRIBUTION::CUSTOM:                                                                                          // CUSTOM

                semiMajorAxis = utils::InverseSampleFromPowerLaw(OPTIONS->SemiMajorAxisDistributionPower(), OPTIONS->SemiMajorAxisDistributionMax(), OPTIONS->SemiMajorAxisDistributionMin());
                break;

            case SEMI_MAJOR_AXIS_DISTRIBUTION::SANA2012: {                                                                                      // Sana et al 2012
                // http://science.sciencemag.org/content/sci/337/6093/444.full.pdf
                // distribution of semi-major axes. Sana et al fit for the orbital period, which we sample in here, before returning the semi major axis
                // Taken from table S3 in http://science.sciencemag.org/content/sci/suppl/2012/07/25/337.6093.444.DC1/1223344.Sana.SM.pdf
                // See also de Mink and Belczynski 2015 http://arxiv.org/pdf/1506.03573v2.pdf

                if (OPTIONS->PeriodDistributionMin() <= 1.0 || OPTIONS->PeriodDistributionMax() <= 1.0) {                                       // bounds check  JR: don't use utils::Compare() here
                    SHOW_WARN(ERROR::OUT_OF_BOUNDS, "Period distribution requires period > 1 day")
                }

                double logPeriodMin = OPTIONS->PeriodDistributionMin() > 1.0 ? log(OPTIONS->PeriodDistributionMin()) : 0.0;                     // smallest initial log period  JR: don't use utils::Compare() here
                double logPeriodMax = OPTIONS->PeriodDistributionMax() > 1.0 ? log(OPTIONS->PeriodDistributionMax()) : 0.0;                     // largest initial log period   JR: don't use utils::Compare() here

                double periodInDays = exp(utils::InverseSampleFromPowerLaw(-0.55, logPeriodMax, logPeriodMin));                                 // draw a period in days from their distribution

                semiMajorAxis = utils::ConvertPeriodInDaysToSemiMajorAxisInAU(p_Mass1, p_Mass2, periodInDays);                                  // convert period in days to semi-major axis in AU
                } break;

            default:                                                                                                                            // unknown distribution
                SHOW_WARN(ERROR::UNKNOWN_A_DISTRIBUTION, "Using default");                                                                      // show warning

                semiMajorAxis = utils::InverseSampleFromPowerLaw(-1.0, 100.0, 0.5);                                                             // calculate semiMajorAxis using power law with default values
        }
    }
    else {                                                                                                                                      // draw from AIS distributions
        // Mass ratio distribution from Adaptive Importance Sampling v1 from Broekgaarden et al. (in prep 2018)
        // Function Returns a random semiMajorAxis drawn from one of the random gaussians defined bu vectors mu_loga & cov_loga
        // Notice-> the mu and cov are in log10(a) space so range is e.g. (-1,3) instead of (0.1, 1000).

        // draw randomly from the random Gaussian chosen with RandomGaussianDraw
        // MuLogA()  = m_MuLogA[aisvariables.RandomGaussianDraw]  = mean of the RandomGaussianDraw-th Gaussian
        // CovLogA() = m_CovLogA[aisvariables.RandomGaussianDraw] = cov of the RandomGaussianDraw-th Gaussin

        semiMajorAxis = pow(10, RAND->RandomGaussian(m_AIS.CovLogA()) + m_AIS.MuLogA());                                                        // draw random number from Gaussian
    }

    return semiMajorAxis;
}


/*
 * Draw mass ratio q from the distribution specified by the user
 * (MassRatioDistribution, EccentricityDistribution program options; will use AIS distribution if specified (AIS.DrawingFromAISDistributions))
 *
 *
 * double SampleQDistribution()
 *
 * @return                                      Mass ratio q
 */
double BaseBinaryStar::SampleQDistribution() {

    double q;

    if (!m_AIS.DrawingFromAISDistributions()) {                                                                                         // draw from priors (not from AIS distributions)
        switch (OPTIONS->MassRatioDistribution()) {

            case MASS_RATIO_DISTRIBUTION::FLAT:                                                                                         // FLAT mass ratio distriution
                q = utils::InverseSampleFromPowerLaw(0.0, OPTIONS->MassRatioDistributionMax(), OPTIONS->MassRatioDistributionMin());
                break;

            case MASS_RATIO_DISTRIBUTION::DUQUENNOYMAYOR1991:                                                                           // mass ratio distribution from Duquennoy & Mayor (1991) (http://adsabs.harvard.edu/abs/1991A%26A...248..485D)

                do {                                                                                                                    // JR: todo: catch non-convergence?
                    q = 0.42 * sqrt(-2.0 * log(RAND->Random())) * cos(2.0 * M_PI * RAND->Random()) + 0.23;
                } while (q < 0.0 || q > 1.0);                                                                                           // JR: don't use utils::Compare() here
                break;

            case MASS_RATIO_DISTRIBUTION::SANA2012:                                                                                     // Sana et al 2012 (http://science.sciencemag.org/content/sci/337/6093/444.full.pdf) distribution of eccentricities.
                // Taken from table S3 in http://science.sciencemag.org/content/sci/suppl/2012/07/25/337.6093.444.DC1/1223344.Sana.SM.pdf
                // See also de Mink and Belczynski 2015 http://arxiv.org/pdf/1506.03573v2.pdf

                q = utils::InverseSampleFromPowerLaw(-0.1, OPTIONS->MassRatioDistributionMax(), OPTIONS->MassRatioDistributionMin());   // de Mink and Belczynski use min = 0.1, max = 1.0
                break;

            default:            // unknown q-distribution - reset to default
                SHOW_WARN(ERROR::UNKNOWN_Q_DISTRIBUTION, "Using default");                                                              // show warning
                q = utils::InverseSampleFromPowerLaw(0.0, 1.0, 0.0);                                                                    // calculate q using power law with default values
        }
    }
    else {                                                                                                                              // draw from AIS distributions
        // draw randomly from the random Gaussian chosen with RandomGaussianDraw
        // MuQ()  = m_MuQ[aisvariables.RandomGaussianDraw]  = mean of the RandomGaussianDraw-th Gaussian
        // CovQ() = m_CovQ[aisvariables.RandomGaussianDraw] = cov of the RandomGaussianDraw-th Gaussin

        q = RAND->RandomGaussian(m_AIS.CovQ()) + m_AIS.MuQ();                                                                           // draw random number from Gaussian
    }

    return q;
}


//JR: todo: talk to Floor about using utils::Compare() in this function
/*
 * Calculate the value of the CDF of the Kroupa (2001) IMF at p_Mass
 *
 *
 * double CalculateCDFKroupa(const double p_Mass)
 *
 * @param   [IN]    p_Mass                      Mass value (in Msol) at which to calculate the CDF
 * @return                                      CDF value
 */
double BaseBinaryStar::CalculateCDFKroupa(const double p_Mass) {

    double CDF = 0.0;

    if (OPTIONS->InitialMassFunctionMin() <= KROUPA_BREAK_1 &&
        OPTIONS->InitialMassFunctionMax() >  KROUPA_BREAK_1 &&
        OPTIONS->InitialMassFunctionMax() <= KROUPA_BREAK_2) {

        double term1 = ONE_OVER_KROUPA_POWER_1_PLUS1 * (KROUPA_BREAK_1_PLUS1_1 - pow(OPTIONS->InitialMassFunctionMin(), KROUPA_POWER_PLUS1_1));
        double term2 = ONE_OVER_KROUPA_POWER_2_PLUS1 * KROUPA_BREAK_1_POWER_1_2 * (pow(OPTIONS->InitialMassFunctionMax(), KROUPA_POWER_PLUS1_2) - KROUPA_BREAK_1_PLUS1_2);

        double C1 = 1.0 / (term1 + term2);
        double C2 = C1 * KROUPA_BREAK_1_POWER_1_2;

        if (p_Mass >= OPTIONS->InitialMassFunctionMin() && p_Mass < KROUPA_BREAK_1) {

            CDF = ONE_OVER_KROUPA_POWER_1_PLUS1 * C1 * (pow(p_Mass, KROUPA_POWER_PLUS1_1) - pow(OPTIONS->InitialMassFunctionMin(), KROUPA_POWER_PLUS1_1));
        }
        else if (p_Mass >= KROUPA_BREAK_1 && p_Mass < KROUPA_BREAK_2) {

            CDF = ONE_OVER_KROUPA_POWER_1_PLUS1 * C1 * (KROUPA_BREAK_1_PLUS1_1 - pow(OPTIONS->InitialMassFunctionMin(), KROUPA_POWER_PLUS1_1)) +
                  ONE_OVER_KROUPA_POWER_2_PLUS1 * C2 * (pow(p_Mass, KROUPA_POWER_PLUS1_2) - KROUPA_BREAK_1_PLUS1_2);
        }
        else {
            SHOW_WARN(ERROR::OUT_OF_BOUNDS, "Using CDF = 0.0 (1)");
        }

    }
    else if (OPTIONS->InitialMassFunctionMin() <= KROUPA_BREAK_1 &&
             OPTIONS->InitialMassFunctionMax() >  KROUPA_BREAK_2) {

        double term1 = ONE_OVER_KROUPA_POWER_1_PLUS1 * (KROUPA_BREAK_1_PLUS1_1 - pow(OPTIONS->InitialMassFunctionMin(), KROUPA_POWER_PLUS1_1));
        double term2 = ONE_OVER_KROUPA_POWER_2_PLUS1 * KROUPA_BREAK_1_POWER_1_2 * (KROUPA_BREAK_2_PLUS1_2 - KROUPA_BREAK_1_PLUS1_2);
        double term3 = ONE_OVER_KROUPA_POWER_3_PLUS1 * KROUPA_BREAK_1_POWER_1_2 * KROUPA_BREAK_2_POWER_2_3 * (pow(OPTIONS->InitialMassFunctionMax(), KROUPA_POWER_PLUS1_3) - KROUPA_BREAK_2_PLUS1_3);

        double C1 = 1.0 / (term1 + term2 + term3);
        double C2 = C1 * KROUPA_BREAK_1_POWER_1_2;
        double C3 = C2 * KROUPA_BREAK_2_POWER_2_3;

        if (p_Mass >= OPTIONS->InitialMassFunctionMin() && p_Mass < KROUPA_BREAK_1) {

            CDF = ONE_OVER_KROUPA_POWER_1_PLUS1 * C1 * (pow(p_Mass, KROUPA_POWER_PLUS1_1) - pow(OPTIONS->InitialMassFunctionMin(), KROUPA_POWER_PLUS1_1));
        }
        else if (p_Mass >= KROUPA_BREAK_1 && p_Mass < KROUPA_BREAK_2) {

            CDF = ONE_OVER_KROUPA_POWER_1_PLUS1 * C1 * (KROUPA_BREAK_1_PLUS1_1 - pow(OPTIONS->InitialMassFunctionMin(), KROUPA_POWER_PLUS1_1)) +
                  ONE_OVER_KROUPA_POWER_2_PLUS1 * C2 * (pow(p_Mass, KROUPA_POWER_PLUS1_2) - KROUPA_BREAK_1_PLUS1_2);
        }
        else if (p_Mass >= KROUPA_BREAK_2 && p_Mass < OPTIONS->InitialMassFunctionMax()) {

            CDF = ONE_OVER_KROUPA_POWER_1_PLUS1 * C1 * (KROUPA_BREAK_1_PLUS1_1 - pow(OPTIONS->InitialMassFunctionMin(), KROUPA_POWER_PLUS1_1)) +
                  ONE_OVER_KROUPA_POWER_2_PLUS1 * C2 * (KROUPA_BREAK_2_PLUS1_2 - KROUPA_BREAK_1_PLUS1_2) +
                  ONE_OVER_KROUPA_POWER_3_PLUS1 * C3 * (pow(p_Mass, KROUPA_POWER_PLUS1_3) - KROUPA_BREAK_2_PLUS1_3);
        }
        else {
            SHOW_WARN(ERROR::OUT_OF_BOUNDS, "Using CDF = 0.0 (2)");
        }

    }
    else if (OPTIONS->InitialMassFunctionMin() >  KROUPA_BREAK_1 &&
             OPTIONS->InitialMassFunctionMin() <= KROUPA_BREAK_2 &&
             OPTIONS->InitialMassFunctionMax() >  KROUPA_BREAK_2) {

        double term1 = ONE_OVER_KROUPA_POWER_2_PLUS1 * (KROUPA_BREAK_2_PLUS1_2 - pow(OPTIONS->InitialMassFunctionMin(), KROUPA_POWER_PLUS1_2));
        double term2 = ONE_OVER_KROUPA_POWER_3_PLUS1 * KROUPA_BREAK_2_POWER_2_3 * (pow(OPTIONS->InitialMassFunctionMax(), KROUPA_POWER_PLUS1_3) - KROUPA_BREAK_2_PLUS1_3);

        double C2 = 1.0 / (term1 + term2);
        double C3 = C2 * KROUPA_BREAK_2_POWER_2_3;

        if (p_Mass >= OPTIONS->InitialMassFunctionMin() && p_Mass < KROUPA_BREAK_2) {

            CDF = ONE_OVER_KROUPA_POWER_2_PLUS1 * C2 * (pow(p_Mass, KROUPA_POWER_PLUS1_2) - pow(OPTIONS->InitialMassFunctionMin(), KROUPA_POWER_PLUS1_2));
        }
        else if (p_Mass >= KROUPA_BREAK_2 && p_Mass < OPTIONS->InitialMassFunctionMax()) {

            CDF = ONE_OVER_KROUPA_POWER_2_PLUS1 * C2 * (KROUPA_BREAK_2_PLUS1_2 - pow(OPTIONS->InitialMassFunctionMin(), KROUPA_POWER_PLUS1_2)) +
                  ONE_OVER_KROUPA_POWER_3_PLUS1 * C3 * (pow(p_Mass, KROUPA_POWER_PLUS1_3) - KROUPA_BREAK_2_PLUS1_3);
        }
        else {
            SHOW_WARN(ERROR::OUT_OF_BOUNDS, "Using CDF = 0.0 (3)");
        }
    }

    return CDF;
}


//JR: todo: talk to Floor about using utils::Compare() in this function
/*
 * Draw mass from the distribution specified by the user
 * (InitialMassFunction program option; will use AIS distribution if specified (AIS.DrawingFromAISDistributions))
 *
 *
 * double SampleInitialMassDistribution()
 *
 * @return                                      Mass
 */
double BaseBinaryStar::SampleInitialMassDistribution() {

    double thisMass = 0.0;

    if (!m_AIS.DrawingFromAISDistributions()) {                                                                                                         // draw from priors (not from AIS distributions)

        switch (OPTIONS->InitialMassFunction()) {                                                                                                       // which IMF?

            case INITIAL_MASS_FUNCTION::SALPETER:                                                                                                       // SALPETER

                thisMass = utils::InverseSampleFromPowerLaw(SALPETER_POWER, OPTIONS->InitialMassFunctionMax(), OPTIONS->InitialMassFunctionMin());
                break;

            case INITIAL_MASS_FUNCTION::POWERLAW:                                                                                                       // POWER LAW

                thisMass = utils::InverseSampleFromPowerLaw(OPTIONS->InitialMassFunctionPower(), OPTIONS->InitialMassFunctionMax(), OPTIONS->InitialMassFunctionMin());
                break;

            case INITIAL_MASS_FUNCTION::UNIFORM:                                                                                                        // UNIFORM - convienience function for POWERLAW with slope of 0

                thisMass = RAND->Random(OPTIONS->InitialMassFunctionMin(), OPTIONS->InitialMassFunctionMax());
                break;

            case INITIAL_MASS_FUNCTION::KROUPA:                                                                                                         // KROUPA

                // find out where the user specificed their minimum and maximum masses to generate
                if (utils::Compare(OPTIONS->InitialMassFunctionMin(), KROUPA_BREAK_1) <= 0 && utils::Compare(OPTIONS->InitialMassFunctionMax(), KROUPA_BREAK_1) <= 0) {
                    thisMass = utils::InverseSampleFromPowerLaw(KROUPA_POWER_1, OPTIONS->InitialMassFunctionMax(), OPTIONS->InitialMassFunctionMin());    // draw mass using inverse sampling
                }
                else if (utils::Compare(OPTIONS->InitialMassFunctionMin(), KROUPA_BREAK_1) > 0 && utils::Compare(OPTIONS->InitialMassFunctionMin(), KROUPA_BREAK_2) <= 0 &&
                         utils::Compare(OPTIONS->InitialMassFunctionMax(), KROUPA_BREAK_1) > 0 && utils::Compare(OPTIONS->InitialMassFunctionMax(), KROUPA_BREAK_2) <= 0) {

                    thisMass = utils::InverseSampleFromPowerLaw(KROUPA_POWER_2, OPTIONS->InitialMassFunctionMax(), OPTIONS->InitialMassFunctionMin());    // draw mass using inverse sampling
                }
                else if (utils::Compare(OPTIONS->InitialMassFunctionMin(), KROUPA_BREAK_2) > 0 && utils::Compare(OPTIONS->InitialMassFunctionMax(), KROUPA_BREAK_2) > 0) {

                    thisMass = utils::InverseSampleFromPowerLaw(KROUPA_POWER_3, OPTIONS->InitialMassFunctionMax(), OPTIONS->InitialMassFunctionMin());    // draw mass using inverse sampling
                }
                else if (utils::Compare(OPTIONS->InitialMassFunctionMin(), KROUPA_BREAK_1) <= 0 &&
                         utils::Compare(OPTIONS->InitialMassFunctionMax(), KROUPA_BREAK_1)  > 0 && utils::Compare(OPTIONS->InitialMassFunctionMax(), KROUPA_BREAK_2) <= 0) {

                    double term1 = ONE_OVER_KROUPA_POWER_1_PLUS1 * (KROUPA_BREAK_1_PLUS1_1 - pow(OPTIONS->InitialMassFunctionMin(), KROUPA_POWER_PLUS1_1));
                    double term2 = ONE_OVER_KROUPA_POWER_2_PLUS1 * KROUPA_BREAK_1_POWER_1_2 * (pow(OPTIONS->InitialMassFunctionMax(), KROUPA_POWER_PLUS1_2) - KROUPA_BREAK_1_PLUS1_2);

                    double C1    = 1.0 / (term1 + term2);
                    double C2    = C1 * KROUPA_BREAK_1_POWER_1_2;
                    double A     = ONE_OVER_KROUPA_POWER_1_PLUS1 * C1 * (KROUPA_BREAK_1_PLUS1_1 - pow(OPTIONS->InitialMassFunctionMin(), KROUPA_POWER_PLUS1_1));

                    double rand  = RAND->Random();                                                                                                      // draw a random number between 0 and 1
                    thisMass = utils::Compare(rand, CalculateCDFKroupa(KROUPA_BREAK_1)) < 0
                                ? pow(rand * (KROUPA_POWER_PLUS1_1 / C1) + pow(OPTIONS->InitialMassFunctionMin(), KROUPA_POWER_PLUS1_1), ONE_OVER_KROUPA_POWER_1_PLUS1)
                                : pow((rand - A) * (KROUPA_POWER_PLUS1_2 / C2) + KROUPA_BREAK_1_PLUS1_2, ONE_OVER_KROUPA_POWER_2_PLUS1);
                }
                else if (utils::Compare(OPTIONS->InitialMassFunctionMin(), KROUPA_BREAK_1) <= 0 && utils::Compare(OPTIONS->InitialMassFunctionMax(), KROUPA_BREAK_2_POWER_2_3) > 0) {

                    double term1 = ONE_OVER_KROUPA_POWER_1_PLUS1 * (KROUPA_BREAK_1_PLUS1_1 - pow(OPTIONS->InitialMassFunctionMin(), KROUPA_POWER_PLUS1_1));
                    double term2 = ONE_OVER_KROUPA_POWER_2_PLUS1 * KROUPA_BREAK_1_POWER_1_2 * (KROUPA_BREAK_2_PLUS1_2 - KROUPA_BREAK_1_PLUS1_2);
                    double term3 = ONE_OVER_KROUPA_POWER_3_PLUS1 * KROUPA_BREAK_1_POWER_1_2 * KROUPA_BREAK_2_POWER_2_3 * (pow(OPTIONS->InitialMassFunctionMax(), KROUPA_BREAK_2_POWER_2_3) - KROUPA_BREAK_2_PLUS1_3);

                    double C1    = 1.0 / (term1 + term2 + term3);
                    double C2    = C1 * KROUPA_BREAK_1_POWER_1_2;
                    double C3    = C2 * KROUPA_BREAK_2_POWER_2_3;

                    double A     = ONE_OVER_KROUPA_POWER_1_PLUS1 * C1 * (KROUPA_BREAK_1_PLUS1_1 - pow(OPTIONS->InitialMassFunctionMin(), KROUPA_POWER_PLUS1_1));
                    double B     = ONE_OVER_KROUPA_POWER_2_PLUS1 * C2 * (KROUPA_BREAK_2_PLUS1_2 - KROUPA_BREAK_1_PLUS1_2);

                    double rand  = RAND->Random();                                                                                                      // draw a random number between 0 and 1

                    if (utils::Compare(rand, CalculateCDFKroupa(KROUPA_BREAK_1)) < 0)
                        thisMass = pow(rand * (KROUPA_POWER_PLUS1_1 / C1) + pow(OPTIONS->InitialMassFunctionMin(), KROUPA_POWER_PLUS1_1), ONE_OVER_KROUPA_POWER_1_PLUS1);
                    else if (utils::Compare(rand, CalculateCDFKroupa(KROUPA_BREAK_2)) < 0)
                        thisMass = pow((rand - A) * (KROUPA_POWER_PLUS1_2 / C2) + KROUPA_BREAK_1_PLUS1_2, ONE_OVER_KROUPA_POWER_2_PLUS1);
                    else
                        thisMass = pow((rand - A - B) * (KROUPA_POWER_PLUS1_3 / C3) + KROUPA_BREAK_2_PLUS1_3, ONE_OVER_KROUPA_POWER_3_PLUS1);
                }
                else if (utils::Compare(OPTIONS->InitialMassFunctionMin(), KROUPA_BREAK_1)  > 0 &&
                         utils::Compare(OPTIONS->InitialMassFunctionMin(), KROUPA_BREAK_2) <= 0 && utils::Compare(OPTIONS->InitialMassFunctionMax(), KROUPA_BREAK_2) > 0) {

                    double term1 = ONE_OVER_KROUPA_POWER_2_PLUS1 * (KROUPA_BREAK_2_PLUS1_2 - pow(OPTIONS->InitialMassFunctionMin(), KROUPA_POWER_PLUS1_2));
                    double term2 = ONE_OVER_KROUPA_POWER_3_PLUS1 * KROUPA_BREAK_2_POWER_2_3 * (pow(OPTIONS->InitialMassFunctionMax(), KROUPA_POWER_PLUS1_3) - KROUPA_BREAK_2_PLUS1_3);

                    double C2    = 1.0 / (term1 + term2);
                    double C3    = C2 * KROUPA_BREAK_2_POWER_2_3;
                    double B     = ONE_OVER_KROUPA_POWER_2_PLUS1 * C2 * (KROUPA_BREAK_2_PLUS1_2 - pow(OPTIONS->InitialMassFunctionMin(), KROUPA_POWER_PLUS1_2));

                    double rand  = RAND->Random();                                                                                                      // draw a random number between 0 and 1

                    thisMass = utils::Compare(rand, CalculateCDFKroupa(KROUPA_BREAK_2)) < 0
                                ? pow(rand * (KROUPA_POWER_PLUS1_2 / C2) + pow(OPTIONS->InitialMassFunctionMin(), KROUPA_POWER_PLUS1_2), ONE_OVER_KROUPA_POWER_2_PLUS1)
                                : pow((rand - B) * (KROUPA_POWER_PLUS1_3 / C3) + KROUPA_BREAK_2_PLUS1_3, ONE_OVER_KROUPA_POWER_3_PLUS1);
                }
                // JR: no other case possible - as long as OPTIONS->InitialMassFunctionMin() < OPTIONS->InitialMassFunctionMax() (currently enforced in Options.cpp)
                break;

            default:                                                                                                                                    // unknown IMF
                SHOW_WARN(ERROR::UNKNOWN_INITIAL_MASS_FUNCTION, "Using default");                                                                       // show warning

                thisMass = utils::InverseSampleFromPowerLaw(KROUPA_POWER, KROUPA_MAXIMUM, KROUPA_MINIMUM);                                              // calculate mass using power law with default values
        }
    }
    else {                                                                                                                                              // draw from AIS distributions if  DrawingFromAISDistributions = true
        // draw a random Mass from the random Gaussian chosen with RandomGaussianDraw
        // MuM1()  = m_MuM1[m_RandomGaussianDraw]  = mean of the RandomGaussianDraw-th Gaussian
        // CovM1() = m_CovM1[m_RandomGaussianDraw] = cov of the RandomGaussianDraw-th Gaussin

        thisMass = RAND->RandomGaussian(m_AIS.CovM1()) + m_AIS.MuM1();                                                                                  // draw random number from Gaussian
    }

    return thisMass;
}


/*
 * Draw eccentricity from the distribution specified by the user
 * (EccentricityDistribution program option; will use AIS distribution if specified (AIS.DrawingFromAISDistributions))
 *
 *
 * double SampleEccentricityDistribution()
 *
 * @return                                      Eccentricity
 */
double BaseBinaryStar::SampleEccentricityDistribution() {

    double eccentricity;

    switch (OPTIONS->EccentricityDistribution()) {                                                  // which distribution?

        case ECCENTRICITY_DISTRIBUTION::ZERO:                                                       // ZERO - all systems are initially circular i.e. have zero eccentricity

            eccentricity = 0.0;
            break;

        case ECCENTRICITY_DISTRIBUTION::FIXED:                                                      // FIXED - all systems have same initial eccentricity - not implemented

            SHOW_WARN(ERROR::UNSUPPORTED_ECCENTRICITY_DISTRIBUTION, "Using eccentricity = 0.0");    // show warning
            eccentricity = 0.0;
            break;

        case ECCENTRICITY_DISTRIBUTION::FLAT:                                                       // FLAT

            eccentricity = utils::InverseSampleFromPowerLaw(0.0, OPTIONS->EccentricityDistributionMax(), OPTIONS->EccentricityDistributionMin());
            break;

        case ECCENTRICITY_DISTRIBUTION::THERMALISED:                                                // THERMA eccentricity distribution p(e) = 2e
        case ECCENTRICITY_DISTRIBUTION::THERMAL:

            eccentricity = utils::InverseSampleFromPowerLaw(1.0, OPTIONS->EccentricityDistributionMax(), OPTIONS->EccentricityDistributionMin());
            break;

        case ECCENTRICITY_DISTRIBUTION::GELLER_2013:                                                // M35 eccentricity distribution from Geller, Hurley and Mathieu 2013
            // Gaussian with mean 0.38 and sigma 0.23
            // http://iopscience.iop.org/article/10.1088/0004-6256/145/1/8/pdf
            // Sampling function taken from binpop.f in NBODY6

            do {                                                                                    // JR: todo: catch non-convergence?
                eccentricity = 0.23 * sqrt(-2.0 * log(RAND->Random())) * cos(2.0 * M_PI * RAND->Random()) + 0.38;
            } while(eccentricity < 0.0 || eccentricity > 1.0);                                      // JR: don't use utils::Compare() here
            break;

         case ECCENTRICITY_DISTRIBUTION::DUQUENNOYMAYOR1991:                                        // eccentricity distribution from Duquennoy & Mayor (1991)
            // http://adsabs.harvard.edu/abs/1991A%26A...248..485D
            // Sampling function taken from binpop.f in NBODY6

            do {                                                                                    // JR: todo: catch non-convergence?
                eccentricity = 0.15 * sqrt(-2.0 * log(RAND->Random())) * cos(2.0 * M_PI * RAND->Random()) + 0.3;
            } while(eccentricity < 0.0 or eccentricity > 1.0);                                      // JR: don't use utils::Compare() here
            break;

        case ECCENTRICITY_DISTRIBUTION::SANA2012:                                                   // Sana et al 2012
            // (http://science.sciencemag.org/content/sci/337/6093/444.full.pdf) distribution of eccentricities.
            // Taken from table S3 in http://science.sciencemag.org/content/sci/suppl/2012/07/25/337.6093.444.DC1/1223344.Sana.SM.pdf
            // See also de Mink and Belczynski 2015 http://arxiv.org/pdf/1506.03573v2.pdf

            eccentricity = utils::InverseSampleFromPowerLaw(-0.42, OPTIONS->EccentricityDistributionMax(), OPTIONS->EccentricityDistributionMin());
            break;

        case ECCENTRICITY_DISTRIBUTION::IMPORTANCE:                                                 // IMPORTANCE - not implemented

            SHOW_WARN(ERROR::UNSUPPORTED_ECCENTRICITY_DISTRIBUTION, "Using eccentricity = 0.0");    // show warning
            eccentricity = 0.0;
            break;

        default:                                                                                    // unknown distribution
            SHOW_WARN(ERROR::UNKNOWN_ECCENTRICITY_DISTRIBUTION, "Using eccentricity = 0.0");        // show warning
            eccentricity = 0.0;
    }

    return eccentricity;
}


/*
 * Draw spin from the distribution specified by the user
 * (SpinDistribution program option; will use AIS distribution if specified (AIS.DrawingFromAISDistributions))
 *
 *
 * double SampleSpinDistribution()
 *
 * @return                                      Spin
 */
double BaseBinaryStar::SampleSpinDistribution() {

    double spin;

    switch (OPTIONS->SpinDistribution()) {                                          // which distribution? .

        case SPIN_DISTRIBUTION::ZERO:                                               // ZERO

            spin = 0.0;
            break;

        case SPIN_DISTRIBUTION::FLAT:                                               // FLAT

            spin = RAND->Random(OPTIONS->SpinDistributionMin(), OPTIONS->SpinDistributionMax());
            break;

        case SPIN_DISTRIBUTION::FIXED:                                              // FIXED

            spin = 0.7;                                                             // could have this set by the user - check bounds if so
            break;

        default:                                                                    // unknown distribution
            SHOW_WARN(ERROR::UNKNOWN_SPIN_DISTRIBUTION, "Using spin = 0.0");        // show warning
            spin = 0.0;
    }

    return spin;
}


/*
 * Choose metallicity based on program option (not really drawing from a distribution here...)
 *
 * Chooses metallicity from user-supplied metallicity (program options FixedMetallicity, Metallicity) or ZSOL (constant)
 *
 *
 * double SampleMetallicityDistribution()
 *
 * @return                                      Metallicity
 */
double BaseBinaryStar::SampleMetallicityDistribution() {
    return OPTIONS->FixedMetallicity() ? OPTIONS->Metallicity() : ZSOL;         // user specified value if provided, else solar metallicity by defaulte
}


/*
 * Determine whether to turn monitoring off for CASE A mass transfer
 * If monitoring is not required (VASE A transfer is not/no longer happening)
 * turn monitoring off (set flag)
 *
 *
 * bool RLOFCheckMonitorEndCaseA()
 *
 * @return                                      Boolean indicating if monitoring was disabled (true = disabled, false = not changed)
 */
bool BaseBinaryStar::RLOFCheckMonitorEndCaseA() {

    bool endMonitoring = false;                                                                                 // default is leave monitoring on

    // check need to stop monitoring star 1 for CASE A transfer
    if (!m_Star1->IsOneOf(MAIN_SEQUENCE) && m_RLOFDetails.currentProps.monitorMS1) {                            // if not a (non-HeMS) Main Sequence star and monitoring is on...
        m_RLOFDetails.currentProps.monitorMS1 = false;                                                          // ... turn it off
        endMonitoring                         = true;                                                           // flag that monitoring has ended
    }

    if (!m_Star1->IsOneOf({ STELLAR_TYPE::NAKED_HELIUM_STAR_MS }) && m_RLOFDetails.currentProps.monitorHeMS1) { // if not a HeMS star and monitoring is on...
        m_RLOFDetails.currentProps.monitorHeMS1 = false;                                                        // ... turn it off
        endMonitoring                           = true;                                                         // flag that monitoring has ended
    }

    // rinse and repeat for star 2
    if (!m_Star2->IsOneOf(MAIN_SEQUENCE) && m_RLOFDetails.currentProps.monitorMS2) {                            // if not a (non-HeMS) Main Sequence star and monitoring is on...
        m_RLOFDetails.currentProps.monitorMS2 = false;                                                          // ... turn it off
        endMonitoring                         = true;                                                           // flag that monitoring has ended
    }

    if (!m_Star2->IsOneOf({ STELLAR_TYPE::NAKED_HELIUM_STAR_MS }) && m_RLOFDetails.currentProps.monitorHeMS2) { // if not a HeMS star and monitoring is on...
        m_RLOFDetails.currentProps.monitorHeMS2 = false;                                                        // ... turn it off
        endMonitoring                           = true;                                                         // flag that monitoring has ended
    }

    return endMonitoring;
}


/*
 * Determine whether RLOF parameters should be written to RLOF logfile
 *
 * Determinant is whether RLOF is occurring and the star undergoing RLOF is not a Main Sequence star (including HeMS)
 * If the determinant is true, RLOF monitoring of MS stars is disabled
 *
 *
 * bool ShouldPrintRLOFParameters()
 *
 * @return                                      Boolean indicating if RLOF parameters should be written to logfile (true = yes, false = no)
 */
bool BaseBinaryStar::ShouldPrintRLOFParameters() {

    bool print = false;                                                                                     // default is don't print

    if (!OPTIONS->RLOFPrinting()) return print;                                                             // nothing to do

    // print if not MS (including HeMS) mass transfer
    if (m_Star1->IsRLOF() && !utils::IsOneOf(m_RLOFDetails.previousProps.stellarType1, ALL_MAIN_SEQUENCE)) {
        print = true;                                                                                       // print
        m_RLOFDetails.currentProps.monitorMS1   = false;                                                    // stop monitoring
        m_RLOFDetails.currentProps.monitorHeMS1 = false;                                                    // stop monitoring
    }

    // rinse and repeat for star 2
    if (m_Star2->IsRLOF() && !utils::IsOneOf(m_RLOFDetails.previousProps.stellarType2, ALL_MAIN_SEQUENCE)) {
        print = true;                                                                                       // print
        m_RLOFDetails.currentProps.monitorMS2   = false;                                                    // stop monitoring
        m_RLOFDetails.currentProps.monitorHeMS2 = false;                                                    // stop monitoring
    }

    return print || RLOFCheckMonitorEndCaseA();
}


/*
 * Write RLOF parameters to RLOF logfile if necessary
 *
 * If RLOF printing is enabled, check if the parameters should be written (via ShouldPrintRLOFParameters())
 * and write the parameters to the file if necessary
 *
 *
 * void PrintRLOFParameters()
 */
void BaseBinaryStar::PrintRLOFParameters() {
    if (!OPTIONS->RLOFPrinting()) return;                       // nothing to do

    if (ShouldPrintRLOFParameters()) {                          // need to print?
        LOGGING->LogRLOFParameters(this);                       // yes - write to log file
        m_RLOFDetails.currentProps.eventCounter += 1;           // every time we print a MT event happened
    }
}


/*
 * Squirrel RLOF properties away
 *
 * Various binary property values are stashed into the p_RLOFProperties struct for use/printing later
 * The p_RLOFProperties will point to either the m_RLOFDetails.currentProps or m_RLOFDetails.previousProps struct
 *
 *
 * void StashRLOFProperties(RLOFPropertiesT &p_RLOFProperties)
 */
void BaseBinaryStar::StashRLOFProperties(RLOFPropertiesT &p_RLOFProperties) {

    p_RLOFProperties.id           = m_ObjectId;
    p_RLOFProperties.randomSeed   = m_RandomSeed;
    p_RLOFProperties.mass1        = m_Star1->Mass();
    p_RLOFProperties.mass2        = m_Star2->Mass();
    p_RLOFProperties.radius1      = m_Star1->Radius();
    p_RLOFProperties.radius2      = m_Star2->Radius();
    p_RLOFProperties.stellarType1 = m_Star1->StellarType();
    p_RLOFProperties.stellarType2 = m_Star2->StellarType();
    p_RLOFProperties.separation   = m_SemiMajorAxisPrime * AU_TO_RSOL;
    p_RLOFProperties.time         = m_Time;
    p_RLOFProperties.isRLOF1      = m_Star1->IsRLOF();
    p_RLOFProperties.isRLOF2      = m_Star2->IsRLOF();
    p_RLOFProperties.isCE         = m_CEDetails.CEEnow;
}


/*
 * Squirrel the previous values of RLOF properties away
 *
 * Calls StashRLOFProperties() with the p_RLOFProperties pointing at the m_RLOFDetails.previousProps struct
 *
 *
 * void StashPreviousRLOFProperties()
 */
void BaseBinaryStar::StashPreviousRLOFProperties() {

    if (!OPTIONS->RLOFPrinting()) return;                                                                   // nothing to do

    if (!m_RLOFDetails.currentProps.monitorMS1   && !m_RLOFDetails.currentProps.monitorMS2 &&               // not monitoring non-HeMS Main Sequence Case A transfer, and ...
        !m_RLOFDetails.currentProps.monitorHeMS1 && !m_RLOFDetails.currentProps.monitorHeMS2) {             // not monitoring HeMS Case A transfer, so ...
        StashRLOFProperties(m_RLOFDetails.previousProps);                                                   // stash the properties
    }
}


/*
 * Squirrel the previous values of RLOF properties away
 *
 * Calls StashRLOFProperties() with the p_RLOFProperties pointing at the m_RLOFDetails.currentProps struct
 *
 *
 * void StashCurrentRLOFProperties()
 */
void BaseBinaryStar::StashCurrentRLOFProperties() {

    if (!OPTIONS->RLOFPrinting()) return;                                                                                   // nothing to do

    // check need to start monitoring star 1 for CASE A transfer
    if (m_Star1->IsRLOF()) {                                                                                                // star 1 overflowing Roche Lobe
             if (m_Star1->IsOneOf(MAIN_SEQUENCE))                          m_RLOFDetails.currentProps.monitorMS1   = true;  // if (non-HeMS) Main Sequence star then monitor
        else if (m_Star1->IsOneOf({ STELLAR_TYPE::NAKED_HELIUM_STAR_MS })) m_RLOFDetails.currentProps.monitorHeMS1 = true;  // if HeMS star then monitor
    }

    // rinse and repeat for star 2
    if (m_Star2->IsRLOF()) {                                                                                                // star 2 overflowing Roche Lobe
             if (m_Star2->IsOneOf(MAIN_SEQUENCE))                          m_RLOFDetails.currentProps.monitorMS2   = true;  // if (non-HeMS) Main Sequence star then monitor
        else if (m_Star2->IsOneOf({ STELLAR_TYPE::NAKED_HELIUM_STAR_MS })) m_RLOFDetails.currentProps.monitorHeMS2 = true;  // if HeMS star then monitor
    }

    StashRLOFProperties(m_RLOFDetails.currentProps);                                                                        // stash the properties
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

    if (!OPTIONS->BeBinaries() || !IsBeBinary()) return;                                                            // nothing to do

    // switch previous<->current (preserves existing current as (new) previous)
    BeBinaryPropertiesT* tmp;
    tmp                             = m_BeBinaryDetails.previousProps;                                              // save pointer to existing previous props
    m_BeBinaryDetails.previousProps = m_BeBinaryDetails.currentProps;                                               // existing current props become new previous props (values will be preserved)
    m_BeBinaryDetails.currentProps  = tmp;                                                                          // new current props points at existing prevous (values will be replaced)

    // now save (new) current
    m_BeBinaryDetails.currentProps->id           = m_ObjectId;                                                      // object id
    m_BeBinaryDetails.currentProps->randomSeed   = m_RandomSeed;                                                    // random seed
    m_BeBinaryDetails.currentProps->dt           = m_Dt;                                                            // timestep
    m_BeBinaryDetails.currentProps->totalTime    = m_BeBinaryDetails.previousProps->dt + m_Dt;                      // total time - accumulate, don't just replace
    m_BeBinaryDetails.currentProps->separation   = m_SemiMajorAxisPrime * AU_TO_RSOL;                               // semi-major axis - change units to Rsol
    m_BeBinaryDetails.currentProps->eccentricity = m_Eccentricity;                                                  // eccentricity

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
 *    m_CommonEnvelopeDetails.postCEE.eccentricity
 *    m_CommonEnvelopeDetails.postCEE.semiMajorAxis
 *    m_CommonEnvelopeDetails.postCEE.rocheLobe1to2
 *    m_CommonEnvelopeDetails.postCEE.rocheLobe2to1
 *    m_RLOFDetails.immediateRLOFPostCEE
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

	m_CEDetails.preCEE.semiMajorAxis = p_SemiMajorAxis * AU_TO_RSOL;        // convert to Rsol
	m_CEDetails.preCEE.eccentricity  = p_Eccentricity;
	m_CEDetails.preCEE.rocheLobe1to2 = p_RocheLobe1to2 * AU_TO_RSOL;        // convert to Rsol
	m_CEDetails.preCEE.rocheLobe2to1 = p_RocheLobe2to1 * AU_TO_RSOL;        // convert to Rsol
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

	m_CEDetails.postCEE.semiMajorAxis = p_SemiMajorAxis * AU_TO_RSOL;                               // convert to Rsol
    m_CEDetails.postCEE.eccentricity  = p_Eccentricity;
	m_CEDetails.postCEE.rocheLobe1to2 = p_RocheLobe1to2 * AU_TO_RSOL;                               // convert to Rsol
	m_CEDetails.postCEE.rocheLobe2to1 = p_RocheLobe2to1 * AU_TO_RSOL;                               // convert to Rsol

    if (utils::Compare(m_Star1->RadiusPostCEE(), m_CEDetails.postCEE.rocheLobe1to2) >= 0 ||         // ALEJANDRO - 28/01/2019 - Check for RLOF immediatedly after the CEE.
        utils::Compare(m_Star2->RadiusPostCEE(), m_CEDetails.postCEE.rocheLobe2to1) >= 0) {         // A check for it during the next timestep is done in evaluateBinary funtion.
        m_RLOFDetails.immediateRLOFPostCEE = true;			                                        // JR: todo: is else false (always)?  or do we want to preserve previous true value?
    }
}


/*
 * Calculate the time to coalescence for a binary with arbitrary eccentricity using interpolation
 *
 * Peters 1964 http://journals.aps.org/pr/pdf/10.1103/PhysRev.136.B1224, eq 5.14
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
 * @return                                      Time to coalescence in SI units (s)
 */
double BaseBinaryStar::CalculateTimeToCoalescence(const double p_SemiMajorAxis,
                                                  const double p_Eccentricity,
                                                  const double p_Mass1,
                                                  const double p_Mass2) {

    double beta    = (64.0 / 5.0) * G * G * G * p_Mass1 * p_Mass2 * (p_Mass1 + p_Mass2) / (C * C * C * C * C);  // defined in Equation 5.9 in Peters 1964 http://journals.aps.org/pr/pdf/10.1103/PhysRev.136.B1224
    double _4_beta = 4.0 * beta;

    double tC = p_SemiMajorAxis * p_SemiMajorAxis * p_SemiMajorAxis * p_SemiMajorAxis / _4_beta;                // time for a circular binary to merge

    // calculate t/tc using the interpolated function

    if (utils::Compare(p_Eccentricity, 0) != 0) {

        double e0_2  = p_Eccentricity * p_Eccentricity;
        double c0    = p_SemiMajorAxis * (1.0 - e0_2) * pow(p_Eccentricity, -12.0/19.0) * pow(1.0 + (121.0 * e0_2 / 304.0), -870.0/2299.0);

		double _4_c0 = c0 * c0 * c0 * c0;

        if (utils::Compare(p_Eccentricity, 0.01) < 0) {
            tC = _4_c0 *pow(p_Eccentricity, 48.0/19.0) / _4_beta;
        }
        else if (utils::Compare(p_Eccentricity, 0.99) > 0) {

            double _1_e0_2 = 1.0 - e0_2;
            tC *= (768.0 / 425.0) * (_1_e0_2 * _1_e0_2 * _1_e0_2 * sqrt(_1_e0_2));                              // approximation of eq. 5.14 of Peters 1964, for high eccentricities
        }
        else {

            double sum = 0.0;
            double de  = p_Eccentricity / 10000;

            for (double e = 0.0; utils::Compare(e, p_Eccentricity) < 0; e += de) {
                double _1_e_2 = 1.0 - (e * e);
                sum += de * pow(e, 29.0 / 19.0) * pow((1.0 + (121.0 / 304.0) * e * e), 1181.0 / 2299.0) / ( _1_e_2 * sqrt( _1_e_2));
            }

            tC = (12.0 / 19.0) * (_4_c0 / beta) * sum;
        }
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
    m_SemiMajorAxisAtDCOFormation = m_SemiMajorAxisPrime;
    m_EccentricityAtDCOFormation  = m_Eccentricity;

    double tC           = CalculateTimeToCoalescence(m_SemiMajorAxisPrime * AU, m_Eccentricity, m_Star1->Mass() * MSOL, m_Star2->Mass() * MSOL);
    m_TimeToCoalescence = (tC / SECONDS_IN_YEAR) * YEAR_TO_MYR;                                                                                 // coalescence time in Myrs

    if (utils::Compare(tC, HUBBLE_TIME) < 0) {                                                                                                  // shorter than HubbleTime (will need to worry about time delays eventually and time when born)
        m_Merged = true;                                                                                                                        // merged in hubble time
        m_MergesInHubbleTime = true;                                                                                                            // why do we have 2 flags that do the same thing?       JR: todo: ...why?

//        if (!OPTIONS->Quiet()) {
//            SAY("Binary merges in less than Hubble time, tc = " << tC << "s, seed = " << m_RandomSeed);                                         // JR: todo: do we want to keep this?  Is it really a debug statement?  Use DBG_WARN()?
//        }
    }
    else {
        m_Merged             = false;                                                                                                           // did not merge
        m_MergesInHubbleTime = false;                                                                                                           // why do we have 2 flags that do the same thing?       JR: todo: ...why?

//        if (!OPTIONS->Quiet()) {
//            SAY("Binary doesn't merge in Hubble time, tc = " << tC << "s, tc/t_Hubble = " << tC / HUBBLE_TIME << ", seed: " << m_RandomSeed);   // JR: todo: do we want to keep this?  Is it really a debug statement?  Use DBG_WARN()?
//        }
    }

    PrintDoubleCompactObjects();                                                                                                                // print (log) double compact object details
}


/*
 * Resolve tidal interactions
 *
 * Created by Alejandro Vigna-Gomez on 11/2015.
 *
 * Tides function was initially developed comparing the previous timestep to the new timestep and evolving
 * the system using our simple tides prescription (see COMPAS notes). We later noticed this created some
 * problems when we had events such as CE or MT. The way to implement tides now is to solve them at the end of
 * binary evolution, tidally locking them. This is still work in progress. ALEJANDRO - 04/10/2016
 *
 *
 * void ResolveTides()
 */
void BaseBinaryStar::ResolveTides() {

    m_aTidesDiff     = 0.0;
	m_OmegaTidesDiff = 0.0;

	// assign new values to "previous" values, for following timestep
    m_EccentricityPrev	  = m_EccentricityPrime;
    m_SemiMajorAxisPrev   = m_SemiMajorAxisPrime;
    m_OrbitalVelocityPrev = m_OrbitalVelocityPrime;

	m_Star1->CalculateOmegaTidesIndividualDiff(m_Star1->OmegaPrev());                                                                                       // set star1 m_OmegaTidesIndividualDiff = 0.0
	m_Star2->CalculateOmegaTidesIndividualDiff(m_Star2->OmegaPrev());                                                                                       // set star2 m_OmegaTidesIndividualDiff = 0.0

    if (utils::Compare(m_OrbitalVelocityPrime, 0.0) <= 0) return;                                                                                           // nothing to do
	if (utils::Compare(m_SemiMajorAxisPrime, 0.0)   <= 0) return;                                                                                           // nothing to do
	if (m_StellarMerger)                                  return;                                                                                           // nothing to do

    switch (OPTIONS->TidesPrescription()) {                                                                                                                 // which tides prescription?

        case TIDES_PRESCRIPTION::NONE: break;                                                                                                               // NONE - nothing to do

        case TIDES_PRESCRIPTION::HUT:                                                                                                                       // HUT - not supported              JR: todo: added this - why don't we just remove this option?
        case TIDES_PRESCRIPTION::LOCKED_ENERGY:                                                                                                             // LOCKED_ENERGY - not supported    JR: todo: why don't we just remove this option?

            SHOW_WARN(ERROR::UNSUPPORTED_TIDES_PRESCRIPTION, "Prescription: " + TIDES_PRESCRIPTION_LABEL.at(OPTIONS->TidesPrescription()));                 // show warning
            break;

        case TIDES_PRESCRIPTION::LOCKED_ANG_MOMENTUM: {

            double totalMass = m_Star1->Mass() + m_Star2->Mass();                                                                                           // total mass of binary in Msol
            double mu 		 = (m_Star1->Mass() * m_Star2->Mass()) / totalMass;                                                                             // reduced mass in Msol

            double star1AngularVelocity = m_Star1->Omega();                                                                                                 // angular velocity in yr-1 units of star1
            double star2AngularVelocity	= m_Star2->Omega();                                                                                                 // angular velocity in yr-1 units of star2

            BinaryConstituentStar* star1Copy = new BinaryConstituentStar(*m_Star1);
            BinaryConstituentStar* star2Copy = new BinaryConstituentStar(*m_Star2);
            star1Copy->SetCompanion(star2Copy);
            star2Copy->SetCompanion(star1Copy);

            star1Copy->ResolveRemnantAfterEnvelopeLossAndSwitch();
            star2Copy->ResolveRemnantAfterEnvelopeLossAndSwitch();

            // calculate moments of inertia
            double star1MomentOfInertia = m_Star1->CalculateMomentOfInertiaAU(star1Copy->Radius());
            double star2MomentOfInertia = m_Star2->CalculateMomentOfInertiaAU(star2Copy->Radius());

            delete star1Copy; star1Copy = nullptr;
            delete star2Copy; star2Copy = nullptr;

            // define gyration radius 'k' using fit from de Mink et al. 2013
            double star1GyrationRadius = m_Star1->CalculateGyrationRadius();
            double star2GyrationRadius = m_Star2->CalculateGyrationRadius();

            double L = CalculateAngularMomentumPrime();
            double E = CalculateTotalEnergyPrime();                                                                                                         // JR: todo: note original code passed star1.m_Radius * RSOL_TO_AU and star2.m_Radius * RSOL_TO_AU, but CalculateTotalEnergy() multiplies by RSOL_TO_AU internally

            // solve polynomial of form: ax^4 + bx^3 + cx^2 + dx + e = 0
            // number of coefficients (is order of polynomial + 1)
            constexpr size_t polynomialOrder   = 4;
            constexpr size_t nCoefficients     = polynomialOrder + 1;
            constexpr size_t nSolutions        = 2 * polynomialOrder;                                                                                       // real part + imaginary part for each solution

            double coefficients[nCoefficients] = {0.0, 0.0, 0.0, 0.0, 0.0};                                                                                 // initialise coefficients {e, d, c, b, a}

            // calculate coefficients
            coefficients[0] = (star1MomentOfInertia * sqrt(G1 * totalMass)) + (star2MomentOfInertia * sqrt(G1 * totalMass));                                // e

            double spin1    = star1MomentOfInertia * star1AngularVelocity;                                                                                  // Initial spin angular momenta of object 1
            double spin2    = star2MomentOfInertia * star2AngularVelocity;                                                                                  // Initial spin angular momenta of object 2
            coefficients[3] = (m_OrbitalVelocityPrime * m_SemiMajorAxisPrime * m_SemiMajorAxisPrime * mu) + spin1 + spin2;                                  // b

            coefficients[4] = sqrt(G1 * totalMass) * mu;                                                                                                    // a

            // do GSL rootfinding magic - solutions in 'solutions' array
            double solutions[nSolutions];
            gsl_poly_complex_workspace *workspace = gsl_poly_complex_workspace_alloc (nCoefficients);
            gsl_poly_complex_solve(coefficients, nCoefficients, workspace, solutions);
            gsl_poly_complex_workspace_free(workspace);

            // now solve for values required
            // JR: todo: the following arrays were declared with dimension [3] but the dimension should
            // be [polynomialOrder], otherwise the loops below could walk off the end of the arrays...
            // Setting the dimension at [3] relies on the first loop ending due to the check condition
            // and 'haveSolution' being set true - not guaranteed, especially since there is a catch for
            // that not happening after the first loop...

            double aSolution[polynomialOrder];                                                                                                              // semi-major axes
            double wsync[polynomialOrder];                                                                                                                  // synchronised orbital velocities
            double E1[polynomialOrder];                                                                                                                     // total energies
            double L1[polynomialOrder];                                                                                                                     // angular momenta
            double deltaE1[polynomialOrder];                                                                                                                // energy differences
            double deltaL1[polynomialOrder];                                                                                                                // angular momenta differences

            constexpr double errorPermitted = 0.5;                                                                                                          // how much slop we'll accept...

            bool haveSolution = false;
            unsigned int  i   = 0;
            do {
                aSolution[i] = solutions[2 * i] * solutions[2 * i];
                wsync[i]     = sqrt(totalMass * G1 / (aSolution[i] * aSolution[i] * aSolution[i]));
                E1[i]	     = CalculateTotalEnergy(aSolution[i], wsync[i], wsync[i], wsync[i], star1GyrationRadius, star2GyrationRadius);
                L1[i]	     = CalculateAngularMomentum(aSolution[i], m_EccentricityPrime, wsync[i], wsync[i], star1GyrationRadius, star2GyrationRadius);   // to verify for conservation of angular momentum
                deltaE1[i]   = std::abs(E - E1[i]);
                deltaL1[i]   = std::abs(L - L1[i]);

                if (deltaE1[i] < (errorPermitted * std::abs(E)) && utils::Compare(solutions[(2 * i) + 1], 0.0) == 0) {                                      // JR: todo: don't use utils::Compare() for the error epsilon?  That's not really an epsilon...
                    m_Star1->CalculateOmegaTidesIndividualDiff(wsync[i]);
                    m_Star2->CalculateOmegaTidesIndividualDiff(wsync[i]);
                    m_aTidesDiff     = aSolution[i] - m_SemiMajorAxisPrev;
                    m_OmegaTidesDiff = wsync[i] - m_OrbitalVelocityPrev;
                    m_OmegaTides     = wsync[i];

                    haveSolution     = true;
                }
            } while (++i < polynomialOrder && !haveSolution);

            if (!haveSolution) {                                                                                                                            // didn't converge - check for solutions which maybe not be exact, but close to instant tidal locking
                i = 0;
                do {
                    if (deltaL1[i] < (errorPermitted * L) && utils::Compare(solutions[(2 * i) + 1], 0.0) == 0) {                                            // JR: todo: don't use utils::Compare() for the error epsilon?  That's not really an epsilon...
                        m_Star1->CalculateOmegaTidesIndividualDiff(wsync[i]);
                        m_Star2->CalculateOmegaTidesIndividualDiff(wsync[i]);
                        m_aTidesDiff     = aSolution[i] - m_SemiMajorAxisPrev;
                        m_OmegaTidesDiff = wsync[i] - m_OrbitalVelocityPrev;
                        m_OmegaTides     = wsync[i];

                        haveSolution     = true;
                    }
                } while (++i < polynomialOrder && !haveSolution);

                if (!haveSolution) {                                                                                                                        // still didn't converge
                    m_Error = ERROR::NO_CONVERGENCE;                                                                                                        // set error
                    SHOW_ERROR(m_Error);                                                                                                                    // show error
                }
            }
            } break;

        default:                                                                                                                                            // unknown prescription
            m_Error = ERROR::UNKNOWN_TIDES_PRESCRIPTION;                                                                                                    // set error
            SHOW_ERROR(m_Error);                                                                                                                            // show error
    }

    // update values for star rotation
    m_Star1->CalculateAngularMomentum();                                                                                                                    // calculate angular momentum for star1
    m_Star2->CalculateAngularMomentum();                                                                                                                    // calculate angular momentum for star2

    m_Star1->IncrementOmega(m_Star1->OmegaTidesIndividualDiff());                                                                                           // update star1 omega for tides difference
    m_Star2->IncrementOmega(m_Star2->OmegaTidesIndividualDiff());                                                                                           // update star2 omega for tides difference

    // update binary
    m_OrbitalVelocityPrime += m_OmegaTidesDiff;                                                                                                             // should here be a diff quantity because of MB?     JR: todo: ?
    m_SemiMajorAxisPrime   += m_aTidesDiff;

    // if CHE enabled, update rotational frequency for constituent stars - assume tidally locked
    if (OPTIONS->CHE_Option() != CHE_OPTION::NONE) m_Star1->SetOmega(m_OrbitalVelocityPrime);
    if (OPTIONS->CHE_Option() != CHE_OPTION::NONE) m_Star2->SetOmega(m_OrbitalVelocityPrime);
}


/*
 * Assign misalignments to S1 and S2 based on assumptions for spin study
 *
 * Set spin misalignment angles - will need some model switch here depending on what you are running -
 * either as given by code for both primary and secondary, secondary given by code and primary = 0,
 * or secondary given by code and primary isotropic (uniform in cos(theta))  JR: todo: get clarity on this description from Alejandro
 *
 * DBL_DBL BaseBinaryStar::CalculateMisalignments()
 *
 * @return                                      Tuple containing misalignment angles theta1 & theta2
 */
DBL_DBL BaseBinaryStar::CalculateMisalignments() {

    double theta1;
    double theta2;

    switch (OPTIONS->SpinAssumption()) {                                                // which spin assumption?

        case SPIN_ASSUMPTION::SAME:                                                     // both same
            theta1 = m_IPrime;
            theta2 = m_IPrime;
            break;

        case SPIN_ASSUMPTION::ALIGNED:                                                  // both aligned
            theta1 = 0.0;
            theta2 = 0.0;
            break;

        case SPIN_ASSUMPTION::MISALIGNED:                                               // secondary misaligned
            theta1 = 0.0;
            theta2 = m_IPrime;
            break;

        case SPIN_ASSUMPTION::ISOTROPIC:                                                // both isotropic
            theta1 = acos((RAND->Random() * 2.0) - 1.0);                                // initial misalignment of star 1 uniform in cos(theta)
            theta2 = acos((RAND->Random() * 2.0) - 1.0);                                // initial misalignment of star 2 uniform in cos(theta)
            break;

         case SPIN_ASSUMPTION::GEROSA:                                                  // Gerosa inspired

            // inspired by Gerosa et al. 2013 who assume that after the 1st SN, the secondary is realigned but the primary is not
            // after the second SN, the secondary is therefore slightly misaligned whilst the primary is misaligned as a function of both kicks
            // to model this, we choose the primary misalignment isotropically, whilst the secondary is given from my code          // JR: todo:  whose code?

            theta1 = acos((RAND->Random() * 2.0) - 1.0);
            theta2 = m_IPrime;
            break;

        default:                                                                        // unknown spin assumption

            SHOW_WARN(ERROR::UNKNOWN_SPIN_ASSUMPTION, "Using theta1 = theta2 = 0.0");   // show warning

            theta1  = 0;                                                                // initial misalignment of star 1
            theta2  = 0;                                                                // initial misalignment of star 2
    }

    return std::make_tuple(theta1, theta2);
}


/*
 * Calculate the systemic velocity (centre-of-mass velocity) of the binary after the supernova
 *
 * Brandt & Podsiadlowski 1995 https://arxiv.org/pdf/astro-ph/9412023.pdf, eq Equation 2.10, or
 * Hurley et al 2002 https://arxiv.org/pdf/astro-ph/0201220.pdf, eq A.14
 *
 *
 * double CalculatePostSNSystemicVelocity(const double p_SNMass,
 *                                        const double p_SNDeltaMass,
 *                                        const double p_CompanionMass,
 *                                        const double p_TotalMassPreSN,
 *                                        const double p_TotalMassPostSN,
 *                                        const double p_KickTheta,
 *                                        const double p_KickPhi)
 *
 * @param   [IN]    p_SNMass                    Mass of the supernoa
 * @param   [IN]    p_SNDeltaMass               Change in mass of the supernova from last timestep
 * @param   [IN]    p_CompanionMass             Mass of the companion
 * @param   [IN]    p_TotalMassPreSN            Total mass of binary before supernova event
 * @param   [IN]    p_TotalMassPostSN           Total mass of binary after supernova event
 * @param   [IN]    p_KickTheta                 Kick direction angle out of the plane
 * @param   [IN]    p_KickPhi                   Kick direction angle in the plane
 * @return                                      Post supernova systemic velocity
 */
double BaseBinaryStar::CalculatePostSNSystemicVelocity(const double p_SNMass,
                                                       const double p_SNDeltaMass,
                                                       const double p_CompanionMass,
                                                       const double p_TotalMassPreSN,
                                                       const double p_TotalMassPostSN,
                                                       const double p_KickTheta,
                                                       const double p_KickPhi) {
    // calculate these once for later use
    double cosPhi    = cos(p_KickPhi);
    double term1_1   = p_SNDeltaMass * p_CompanionMass / p_TotalMassPreSN;
    double term3_4_1 = p_SNMass * term1_1;

    // calculate the systemic velocity
    double term1     = term1_1 * term1_1;
    double term2     = p_SNMass * p_SNMass * m_uK * m_uK;
    double term3     = 2.0 * term3_4_1 * m_uK * sin(p_KickTheta) * cosPhi * cos(m_Beta);
    double term4     = 2.0 * term3_4_1 * m_uK * cos(p_KickTheta) * cosPhi * sin(m_Beta);

    return (m_VRel / p_TotalMassPostSN) * sqrt(term1 + term2 + term3 + term4);
}


/*
 * Calculate cos(i), where i = the tilt between the pre and post SN orbital planes (as defined by the angular momentum)
 * Eq (40) in post-SN orbital characteristics 2 notes (Alejandro's?  JR: todo: get proper reference)
 *
 *
 * double CalculateCosFinalPlaneTilt(const double p_KickTheta, const double p_KickPhi)
 *
 * @param   [IN]    p_KickTheta                 Kick direction angle out of the plane
 * @param   [IN]    p_KickPhi                   Kick direction angle in the plane
 * @return                                      cos(i)
 */
double BaseBinaryStar::CalculateCosFinalPlaneTilt(const double p_KickTheta, const double p_KickPhi) {

    // calculate these once for use later
    double sinTheta                = sin(p_KickTheta);
    double cosTheta                = cos(p_KickTheta);
    double sinPhi                  = sin(p_KickPhi);
    double cosPhi                  = cos(p_KickPhi);
    double sinBeta                 = sin(m_Beta);
    double cosBeta                 = cos(m_Beta);
    double ukCosThetaCosPhiPlus1   = m_uK * cosTheta * cosPhi + 1.0;
    double ukCosThetaSinPhiCosBeta = m_uK * cosTheta * sinPhi * cosBeta;

    // calculate cos(tilt)
    double top    = (sinBeta  * ukCosThetaCosPhiPlus1) - ukCosThetaSinPhiCosBeta;
    double bottom = sqrt(
                        (m_uK * m_uK * sinTheta * sinTheta) +
                        (ukCosThetaSinPhiCosBeta * ukCosThetaSinPhiCosBeta) +
                        (sinBeta * sinBeta * ukCosThetaCosPhiPlus1 * ukCosThetaCosPhiPlus1) -
                        (2.0 * cosTheta * sinPhi * cosBeta * sinBeta * ukCosThetaCosPhiPlus1)
                    );

    return top / bottom;
}


/*
 * Calculate the post-supernova orbital eccentricity
 *
 * Post-SN orbital characteristics 2 notes, eq 31  (Alejandro's?  JR: todo: get proper reference)
 *
 * Simplifies to Brandt & Podsiadlowski 1995 (http://arxiv.org/abs/astro-ph/9412023), eq 2.8 when e = 0, beta = pi/2
 * Also given in Hurley et al 2002 (http://arxiv.org/pdf/astro-ph/0201220v1.pdf), eq A.12
 *
 *
 * double CalculateOrbitalEccentricityPostSupernova(const double p_KickVelocity,
 *                                                  const double p_TotalMassPreSN,
 *                                                  const double p_TotalMassPostSN,
 *                                                  const double p_KickTheta,
 *                                                  const double p_KickPhi)
 *
 * @param   [IN]    p_KickVelocity              Dimensionless kick velocity vk/vrel
 * @param   [IN]    p_TotalMassPreSN            Total mass of binary before supernova event
 * @param   [IN]    p_TotalMassPostSN           Total mass of binary after supernova event
 * @param   [IN]    p_KickTheta                 Kick direction angle out of the plane
 * @param   [IN]    p_KickPhi                   Kick direction angle in the plane
 * @return                                      Orbital eccentricity after a supernova
 */
double BaseBinaryStar::CalculateOrbitalEccentricityPostSupernova(const double p_KickVelocity,
                                                                 const double p_TotalMassPreSN,
                                                                 const double p_TotalMassPostSN,
                                                                 const double p_KickTheta,
                                                                 const double p_KickPhi) {
    // calculate these once for use later
    double mOverMprime           = p_TotalMassPreSN / p_TotalMassPostSN;
    double uk_2                  = p_KickVelocity * p_KickVelocity;
    double _2_r_Minus_1_a        = (2.0 / m_Radius) - (1.0 / m_SemiMajorAxisPrime);
    double sinTheta              = sin(p_KickTheta);
    double cosTheta              = cos(p_KickTheta);
    double sinPhi                = sin(p_KickPhi);
    double cosPhi                = cos(p_KickPhi);
    double sinBeta               = sin(m_Beta);
    double cosBeta               = cos(m_Beta);
    double ukCosTheta            = p_KickVelocity * cosTheta;
    double ukCosThetaCosPhi      = p_KickVelocity * cosTheta * cosPhi;
    double ukCosThetaCosPhiPlus1 = ukCosThetaCosPhi + 1.0;

    // calculate orbital eccentricity
    double quadraticTerm         = 1.0 + (2.0 * ukCosThetaCosPhi) + uk_2;
    double firstSquareBrackets   = (uk_2 * sinTheta * sinTheta) + (((ukCosTheta * sinPhi * cosBeta) - (sinBeta * ukCosThetaCosPhiPlus1)) * (ukCosTheta * sinPhi * cosBeta - sinBeta * ukCosThetaCosPhiPlus1));
    double secondSquareBrackets  = (2.0 / m_Radius) - (mOverMprime * _2_r_Minus_1_a * quadraticTerm);
    double oneMinusESquared      = m_Radius * m_Radius * mOverMprime * _2_r_Minus_1_a * firstSquareBrackets*secondSquareBrackets;
    double eSquared              = 1.0 - oneMinusESquared;

    if(eSquared < 1E-8) eSquared = 0.0;     // Deal with small number rounding problems - don't use utils::Compare() here      JR: todo: this should be fixed

    return sqrt(eSquared);

}


/*
 * Calculate the post-supernova semi-major axis
 *
 * Post-SN orbital characteristics 2 document, eq 22        JR: todo get reference to document
 *
 *
 * double CalculateSemiMajorAxisPostSupernova(const double p_KickVelocity,
 *                                            const double p_TotalMassPreSN,
 *                                            const double p_TotalMassPostSN,
 *                                            const double p_KickTheta,
 *                                            const double p_KickPhi)
 *
 * @param   [IN]    p_KickVelocity              Dimensionless kick velocity vk/vrel
 * @param   [IN]    p_TotalMassPreSN            Total mass of binary before supernova event
 * @param   [IN]    p_TotalMassPostSN           Total mass of binary after supernova event
 * @param   [IN]    p_KickTheta                 Kick direction angle out of the plane
 * @param   [IN]    p_KickPhi                   Kick direction angle in the plane
 * @return                                      Semi major axis of the orbit after the supernova
 */
double BaseBinaryStar::CalculateSemiMajorAxisPostSupernova(const double p_KickVelocity,
                                                           const double p_TotalMassPreSN,
                                                           const double p_TotalMassPostSN,
                                                           const double p_KickTheta,
                                                           const double p_KickPhi) {
    double r_2           = 2.0 / m_Radius;
    double quadraticTerm = 1.0 + (2.0 * p_KickVelocity * cos(p_KickTheta) * cos(p_KickPhi)) + (p_KickVelocity * p_KickVelocity);

    return 1.0 / (r_2 - ((p_TotalMassPreSN / p_TotalMassPostSN) * (r_2 - (1.0 / m_SemiMajorAxisPrime)) * quadraticTerm));
}


/*
 * Resolves supernova event - one of the stars has gone supernova!
 *
 * Assign a random supernova kick according to the user specified options and then update the orbit
 *
 * JR: todo: flesh-out this documentation
 *
 *
 * bool ResolveSupernova()
 *
 * @return                                      True if a supernova event occurred, otherwise false
 */
bool BaseBinaryStar::ResolveSupernova() {

    if (!m_Supernova->IsSN() && !m_Supernova->IsECSN() && !m_Supernova->IsUSSN()) return false;                                 // not a supernova event - bail out (or bale out depending whence you hail...) passively

	// Masses should already be correct, mass before SN given by star.m_MassPrev
    // Generate true anomaly - (for e=0, should be a flat distribution) - updates Eccentric anomaly and True anomaly automatically
    // ALEJANDRO - 09/05/2018 - If statement to avoid solving Kepler's equation for an unbound orbit; it may be of interest to have SN of unbound stars in the supernovae.txt file.

    if (IsUnbound()) {      // JR: todo: check this - was just "if (m_SemiMajorAxisPrime > 0.0)"
        // ALEJANDRO - 09/05/2018 - Following 3 lines copied from else statement in the end.                                    // JR: todo: are these going to be executed twice...? (I removed one... not required)
        m_Disbound = true;
        m_Survived = false;
    }
    else {
        m_Supernova->CalculateSNAnomalies(m_Eccentricity);
    }

	m_Radius = (m_SemiMajorAxisPrime * (1.0 - (m_Eccentricity * m_Eccentricity))) / (1.0 + m_Eccentricity * cos(m_Supernova->SN_TrueAnomaly()));   // radius of orbit at current time in AU as a function of the true anomaly psi

	double totalMass        = m_Supernova->MassPrev() + m_Companion->MassPrev();                                                    // total mass of binary before supernova event
	double reducedMass      = (m_Supernova->MassPrev() * m_Companion->MassPrev()) / totalMass;                                      // reduced mass before supernova event
	double totalMassPrime   = m_Supernova->Mass() + m_Companion->Mass();                                                            // total mass of binary after supernova event
	double reducedMassPrime = (m_Supernova->Mass() * m_Companion->Mass()) / totalMassPrime;                                         // reduced mass after supernova event

    #define a m_SemiMajorAxisPrime  // for convenience - undefined below
    #define e m_Eccentricity        // for convenience - undefined below
    #define r m_Radius              // for convenience - undefined below

    m_Beta = utils::Compare(e, 0.0) == 0 ? M_PI_2 : asin(sqrt((a * a * (1.0 - (e * e))) / ((2.0 * r * a) - (r * r))));              // angle between the position and velocity vectors

    #undef r
    #undef e
    #undef a

    double vK = m_Supernova->CalculateSNKickVelocity(m_Supernova->Mass(), m_Supernova->MassPrev() - m_Supernova->Mass());           // draw kick velocoty from distribution

    ///////////////////////////////////////////////////////////////////////////////////
	//          AT THE MOMENT, QUANTITIES BEYOND HERE ARE IN SI (NOT IDEAL)          //                                             // JR: todo: do we need to change this?
	///////////////////////////////////////////////////////////////////////////////////

	// Calculate orbital velocity at some true anomaly psi - default is a circular orbit, V = sqrt(gm/a) = const.
    // Since this equation contains 'G', all other quantities must be in SI to get answer in ms^-1

	vK                       *= KM;                                                                                                 // convert vK to m s^-1.  Would be nice to draw this in nicer units to avoid this secion
	m_VRel                    = sqrt(G * (totalMass * MSOL) * ((2.0 / (m_Radius * AU)) - (1.0 / (m_SemiMajorAxisPrime * AU))));     // orbital velocity
	m_uK                      = OPTIONS->UseFixedUK() ? OPTIONS->FixedUK() : vK / m_VRel;                                           // fix uK to user-defined value if required, otherwise calculate it.  uK is dimensionless
	m_OrbitalVelocityPre2ndSN = m_VRel;                                                                                             // since the kick velocity always occurs in equations as vk/vrel, we need to know vrel

	///////////////////////////////////////////////////////////////////////////////////
	//                       SHOULD BE BACK TO NICE UNITS NOW                        //
	///////////////////////////////////////////////////////////////////////////////////

    m_Supernova->SetPreSNeOrbitalEnergy(CalculateOrbitalEnergy(reducedMass, totalMass, m_SemiMajorAxisPrime));                      // pre-SN orbital energy - should be -ve by construction

	// seemed to be getting into this loop occasionally with E > 0 but E ~ 0 (1e-37 for example) -- what's going on?
    // JR: todo: remove this if we're not seeing the problem...
    // don't use utils::Compare() here - let's see if this turns up as a problem
    DBG_ID_IF(m_Supernova->PreSNeOrbitalEnergy() > 0.0, "orbitalEnergy > 0! totalMass = " << totalMass << ", reducedMass = " << reducedMass << ", m_SemiMajorAxisPrime = " << m_SemiMajorAxisPrime);

	// calculate post-SN orbital properties

    // Record the semi major axis and eccentricity just before each supernova
    m_SemiMajorAxisPre2ndSN = m_SemiMajorAxisPrime;
    m_EccentricityPre2ndSN  = m_Eccentricity;

    m_SemiMajorAxisPrime    = CalculateSemiMajorAxisPostSupernova(m_uK, totalMass, totalMassPrime, m_Supernova->SN_Theta(), m_Supernova->SN_Phi());
    double ePrime           = CalculateOrbitalEccentricityPostSupernova(m_uK, totalMass, totalMassPrime, m_Supernova->SN_Theta(), m_Supernova->SN_Phi());

    m_Supernova->SetPostSNeOrbitalEnergy(CalculateOrbitalEnergy(reducedMassPrime, totalMassPrime, m_SemiMajorAxisPrime));           // post-SN orbital energy, check if still bound
    double epsilon     = -m_Supernova->PostSNeOrbitalEnergy() / m_Supernova->PreSNeOrbitalEnergy();                                 // dimensionless post-SN orbital energy

    m_CosIPrime        = 0.0;
    m_IPrime           = 0.0;
    m_SystemicVelocity = 0.0;

    if (utils::Compare(epsilon, 0.0) < 0) {		                                                                                    // still bound?
                                                                                                                                    // yes
        m_Survived = true;                                                                                                          // it survived

        // Calculate post-SN orbital inclination using the equation for arbitrary eccentricity orbits
        m_CosIPrime   = CalculateCosFinalPlaneTilt(m_Supernova->SN_Theta(), m_Supernova->SN_Phi());
        m_IPrime      = acos(m_CosIPrime);

        std::tie(m_Theta1_i, m_Theta2_i) = CalculateMisalignments();                                                                // assign the spins.  TODO: I think this function is currently broken for two supernovae -- check!!  JR: todo" check this

        // variables to evolve
        m_Theta1 = m_Theta1_i;
        m_Theta2 = m_Theta2_i;

        m_SystemicVelocity = CalculatePostSNSystemicVelocity(m_Supernova->Mass(),                                                   // post-SN systemic (center-of-mass) velocity in ms s^-1
                                                             m_Supernova->MassPrev() - m_Supernova->Mass(),
                                                             m_Companion->Mass(),
                                                             totalMass,
                                                             totalMassPrime,
                                                             m_Supernova->SN_Theta(),
                                                             m_Supernova->SN_Phi());
        m_SystemicVelocity /= KM;                                                                                                   // convert to km s^-1
    }
    else {                                                                                                                          // no longer bound
        m_Disbound = true;
        m_Survived = false;                                                                                                         // it did not survive
    }

    m_Companion->CheckRunaway(m_Disbound, m_Survived);                                                                              // flag companion if runaway


    // update some binary parameters
    m_MSN               = m_Supernova->MassPrev();                                                                                  // exploding star pre-SN mass
    m_MSNPrime          = m_Supernova->Mass();                                                                                      // exploding star post-SN mass
    m_MC                = m_Companion->MassPrev();                                                                                  // companion star pre-SN mass
    m_MCPrime           = m_Companion->Mass();                                                                                      // companion star post-SN mass

    m_EPrime            = m_Supernova->PostSNeOrbitalEnergy();

    m_Eccentricity      = ePrime;
    m_EccentricityPrime = ePrime;

    PrintSupernovaDetails();                                                                                                        // log record to supernovae logfile

    m_Supernova->ClearCurrentSNEvent();

    return true;
}


/*
 * Determine if one or both of the stars are undergoing a supernova event,
 * and if so resolve the event(s) by calling ResolveSupernova() for each of
 * the stars as appropriate
 *
 *
 * void EvaluateSupernovae()
 */
void BaseBinaryStar::EvaluateSupernovae() {

    m_SupernovaState = SN_STATE::NONE;                                                                                                  // not yet determined

    if (m_Star1->IsSNevent()) {                                                                                                         // star1 supernova
        m_SupernovaState = m_Star2->IsSNevent() && utils::Compare(m_SemiMajorAxisPrime, 0.0) > 0 ? SN_STATE::BOTH : SN_STATE::STAR1;    // star1 or both

        // resolve star1 supernova
        m_Supernova = m_Star1;                                                                                                          // supernova
        m_Companion = m_Star2;                                                                                                          // companion
        (void)ResolveSupernova();                                                                                                       // resolve supernova
    }

    if (m_Star2->IsSNevent()) {                                                                                                         // star2 supernova
        m_SupernovaState = m_SupernovaState == SN_STATE::NONE ? SN_STATE::STAR2 : SN_STATE::BOTH;                                       // star2 or both

        // resolve star2 supernova
        m_Supernova = m_Star2;                                                                                                          // supernova
        m_Companion = m_Star1;                                                                                                          // companion
        (void)ResolveSupernova();                                                                                                       // resolve supernova
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

    BinaryConstituentStar* star1Copy = new BinaryConstituentStar(*m_Star1);                                             // clone star1 before CEE
	BinaryConstituentStar* star2Copy = new BinaryConstituentStar(*m_Star2);                                             // clone star2 before CEE
	star1Copy->SetCompanion(star2Copy);                                                                                 // need companion for CalculateSynchronisationTimescale() later
	star2Copy->SetCompanion(star1Copy);                                                                                 // need companion for CalculateSynchronisationTimescale() later

    double alphaCE = m_CEDetails.alpha;                                                                                 // CE efficiency parameter

    double semiMajorAxis = m_SemiMajorAxisPrime;                                                                        // current semi-major axis in default units, AU (before CEE)
	double eccentricity  = m_EccentricityPrime;								                                            // current eccentricity (before CEE)

    bool donorMS = false;                                                                                               // check for main sequence donor
    if (OPTIONS->AllowMainSequenceStarToSurviveCommonEnvelope()) {                                                      // allow main sequence stars to survive CEE?
        if (m_Star1->IsOneOf(ALL_MAIN_SEQUENCE)) {                                                                      // yes - star1 MS_LTE_07, MS_GT_07 or NAKED_HELIUM_STAR_MS?
            donorMS      = donorMS || m_Star1->IsRLOF();                                                                // yes - donor MS?
            m_Mass1Final = m_Star1->Mass();                                                                             // set mass
            m_MassEnv1   = 0.0;                                                                                         // no envelope
        }
        else {                                                                                                          // no, star1 not MS_LTE_07, MS_GT_07 or NAKED_HELIUM_STAR_MS
            m_Mass1Final = m_Star1->CoreMass();                                                                         // set mass
            m_MassEnv1   = m_Star1->Mass() - m_Star1->CoreMass();                                                       // and envelope
        }

        if (m_Star2->IsOneOf(ALL_MAIN_SEQUENCE)) {                                                                      // star2 MS_LTE_07, MS_GT_07 or NAKED_HELIUM_STAR_MS?
            donorMS      = donorMS || m_Star2->IsRLOF();                                                                // yes - donor MS?
            m_Mass2Final = m_Star2->Mass();                                                                             // yes - set mass
            m_MassEnv2   = 0.0;                                                                                         // no envelope
        }
        else {                                                                                                          // no, star2 not MS_LTE_07, MS_GT_07 or NAKED_HELIUM_STAR_MS
            m_Mass2Final = m_Star2->CoreMass();                                                                         // set mass
            m_MassEnv2   = m_Star2->Mass() - m_Star2->CoreMass();                                                       // and envelope
        }
    }
    else {                                                                                                              // no don't allow main sequence stars to survive CEE
        m_Mass1Final = m_Star1->CoreMass();                                                                             // set mass1
        m_MassEnv1   = m_Star1->Mass() - m_Star1->CoreMass();                                                           // and envelope1
        m_Mass2Final = m_Star2->CoreMass();                                                                             // set mass2
        m_MassEnv2   = m_Star2->Mass() - m_Star2->CoreMass();                                                           // and envelope2
    }

    bool envelopeFlag1 = utils::Compare(m_MassEnv1, 0.0) > 0 && utils::Compare(m_Mass1Final, 0.0) > 0;                  // star1 not massless remnant and has envelope?
    bool envelopeFlag2 = utils::Compare(m_MassEnv2, 0.0) > 0 && utils::Compare(m_Mass2Final, 0.0) > 0;                  // star1 not massless remnant and has envelope?

    double rRLd1 = semiMajorAxis * CalculateRocheLobeRadius_Static(m_Star1->Mass(), m_Star2->Mass());                   // Roche Lobe radius in AU at the moment where CEE begins, seen by star1
    double rRLd2 = semiMajorAxis * CalculateRocheLobeRadius_Static(m_Star2->Mass(), m_Star1->Mass());                   // Roche Lobe radius in AU at the moment where CEE begins, seen by star2

    m_CEDetails.CEEcount++;                                                                                             // increment CEE count
    m_RLOFDetails.simultaneousRLOF = m_Star1->IsRLOF() && m_Star2->IsRLOF();                                            // ALEJANDRO - 29/01/2019 - Check for simultaneous RLOF

    // SIMON: Do we really need to calculate all of these?
	// Alejandro: Not really, but I like to have them handy for the common envelopes file.
	// It is only done once or twice per run, so I think it should be fine.
	m_Star1->CalculateLambdas(m_MassEnv1);                                                                              // calculate lambdas for star1
	m_Star2->CalculateLambdas(m_MassEnv2);                                                                              // calculate lambdas for star2

    m_Star1->CalculateBindingEnergies(m_Mass1Final, m_MassEnv1, m_Star1->Radius());                                     // calculate binding energies for star1 (uses lambdas)
    m_Star2->CalculateBindingEnergies(m_Mass2Final, m_MassEnv2, m_Star2->Radius());                                     // calculate binding energies for star2 (uses lambdas)

    m_Star1->CalculateCommonEnvelopeValues();                                                                           // calculate common envelope values for star1
    m_Star2->CalculateCommonEnvelopeValues();                                                                           // calculate common envelope values for star2

    double lambda1 = m_Star1->LambdaAtCEE();                                                                            // measures "the central concentration of the donor"
    double lambda2 = m_Star2->LambdaAtCEE();                                                                            // measures "the central concentration of the accretor"

    if (HasOneOf(ALL_HERTZSPRUNG_GAP)) {                                                                                // check if we have an HG star, and if are allowing such a star to survive the CE       JR: todo: aren't we only checking the first condition?
        m_CEDetails.optimisticCE = true;			                                                                    // JR: todo: is else false (always)?  or do we want to preserve previous true value?
	}

    m_Star1->SetPreCEEValues();                                                                                         // squirrel away pre CEE stellar values for star 1
    m_Star2->SetPreCEEValues();                                                                                         // squirrel away pre CEE stellar values for star 2
  	SetPreCEEValues(semiMajorAxis, eccentricity, rRLd1, rRLd1);                                                         // squirrel away pre CEE binary values

    m_Star1->SetPostCEEValues();                                                                                        // squirrel away (initial) post CEE stellar values for star 1 - default is just pre CEE values
    m_Star2->SetPostCEEValues();                                                                                        // squirrel away (initial) post CEE stellar values for star 2 - default is just pre CEE values

	// double common envelope phase prescription (Brown 1995) to calculate new semi-major axis
	// due to the CEE as described in Belczynsky et al. 2002, eq. (12)
    double k1            = m_Star1->IsOneOf(COMPACT_OBJECTS) ? 0.0 : (2.0 / (lambda1 * alphaCE)) * m_Star1->Mass() * m_MassEnv1 / rRLd1;
    double k2            = m_Star2->IsOneOf(COMPACT_OBJECTS) ? 0.0 : (2.0 / (lambda2 * alphaCE)) * m_Star2->Mass() * m_MassEnv2 / rRLd2;
    double k3            = m_Star1->Mass() * m_Star2->Mass() / semiMajorAxis;
    double k4            = (m_Mass1Final * m_Mass2Final);
    double aFinal        = k4 / (k1 + k2 + k3);                                                                         // semi-major axis after CEE
    m_SemiMajorAxisPrime = aFinal;

    m_CEDetails.doubleCoreCE = utils::Compare(k1, 0.0) > 0 && utils::Compare(k2, 0.0) > 0 && utils::Compare(k3, 0.0) > 0 && utils::Compare(k4, 0.0) > 0;

    double omegafinal      = sqrt((m_Mass1Final + m_Mass2Final) * G1 / (aFinal * aFinal * aFinal));                     // orbital velocity after CEE   SIMON: Should be a function for this equation     JR: todo: is it used elsewhere?
    m_OrbitalVelocityPrime = omegafinal;

    // if CHE enabled, update rotational frequency for constituent stars - assume tidally locked
    if (OPTIONS->CHE_Option() != CHE_OPTION::NONE) m_Star1->SetOmega(m_OrbitalVelocityPrime);
    if (OPTIONS->CHE_Option() != CHE_OPTION::NONE) m_Star2->SetOmega(m_OrbitalVelocityPrime);

	double rRLdfin1 = aFinal * CalculateRocheLobeRadius_Static(m_Mass1Final, m_Mass2Final);                             // Roche Lobe radius in AU after CEE, seen by star1
	double rRLdfin2 = aFinal * CalculateRocheLobeRadius_Static(m_Mass2Final, m_Mass1Final);                             // Roche Lobe radius in AU after CEE, seen by star2

    // Correct for stellar types, stellar mass, separation and period
    // We assume that a common envelope event (CEE) circularises the binary
    m_Eccentricity      = 0.0;
    m_EccentricityPrime = 0.0;

    m_Star1->ResolveCommonEnvelopeAccretion(m_Mass1Final);                                                              // update star's mass after accretion
    m_Star2->ResolveCommonEnvelopeAccretion(m_Mass2Final);                                                              // update star's mass after accretion

    // update stellar type after losing its envelope. Star1, Star2 or both if double CEE.
	// Alejandro - 18/02018 - Calculate tidal timescales

    if (donorMS || (!envelopeFlag1 && !envelopeFlag2)) {                                                                // stellar merger
        m_MassTransferTrackerHistory = HasTwoOf({ STELLAR_TYPE::NAKED_HELIUM_STAR_MS }) ? MT_TRACKING::CE_BOTH_MS : MT_TRACKING::CE_MS_WITH_CO; // Here MS-WD systems are flagged as CE_BOTH_MS
        m_StellarMerger              = true;
    }
	else {
        double periastronRsol = semiMajorAxis * AU_TO_RSOL* (1.0 - eccentricity);                                       // periastron in Rsol

        if (envelopeFlag1) {                                                                                            // star1 donor
            m_Star1->ResolveEnvelopeLossAndSwitch();                                                                    // resolve envelope loss for star1 and switch to new stellar type

            m_SynchronizationTimescale = star1Copy->CalculateSynchronisationTimescale(periastronRsol);
            m_CircularizationTimescale = m_SynchronizationTimescale;

            if (envelopeFlag2) {                                                                                        // correction - double CEE   JR: todo: why do we check envelopeFlags and not value of m_CEDetails.doubleCoreCE calculated above?
                m_Star2->ResolveEnvelopeLossAndSwitch();                                                                // resolve envelope loss for star2 and switch to new stellar type
                m_MassTransferTrackerHistory = MT_TRACKING::CE_DOUBLE_CORE;                                             // record history - double CEE
            }
            else {
                m_MassTransferTrackerHistory = MT_TRACKING::CE_FROM_1_TO_2;                                             // record history - star1 -> star2
            }
        }
        else {                                                                                                          // star2 donor
            m_Star2->ResolveEnvelopeLossAndSwitch();                                                                    // resolve envelope loss for star2 and switch to new stellar type

            m_SynchronizationTimescale   = star2Copy->CalculateSynchronisationTimescale(periastronRsol);
            m_CircularizationTimescale   = m_SynchronizationTimescale;
            m_MassTransferTrackerHistory = MT_TRACKING::CE_FROM_2_TO_1;                                                 // record history - star2 -> star1
        }

        m_Star1->SetPostCEEValues();                                                                                    // squirrel away post CEE stellar values for star 1 - update default values
        m_Star2->SetPostCEEValues();                                                                                    // squirrel away post CEE stellar values for star 2 - update default values
	}

	delete star1Copy; star1Copy = nullptr;
	delete star2Copy; star2Copy = nullptr;

    if (utils::Compare(aFinal, 0.0) <= 0 || utils::Compare(m_Star1->RadiusPostCEE() + m_Star2->RadiusPostCEE(), aFinal * AU_TO_RSOL) > 0) {
        m_StellarMerger = true;			                                                                                // JR: todo: is else false (always)?  or do we want to preserve previous true value?
    }

	SetPostCEEValues(aFinal, m_Eccentricity, rRLdfin1, rRLdfin2);                                                       // squirrel away post CEE binary values.  ALEJANDRO - 06/12/2016 - for populations studies. All separations in Rsol.

    PrintCommonEnvelope();
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
    return 0.49 / (0.6 + pow(q, -2.0 / 3.0) * log(1.0 + pow(q, 1.0 / 3.0)));
}


/*
 * Calculate the amount of mass that a star with radiative envelope needs to lose in order to just fill its Roche Lobe.
 *
 * Based on the function –rlof_method in binary_c.      JR: todo: reference?
 *
 * JR: todo: felsh-out this documentation
 *
 *
 * double CalculateAdaptiveRocheLobeOverFlow(const double p_JLoss)
 * @param   [IN]    p_JLoss                     Specific angular momentum with which mass is lost during non-conservative mass transfer
 * @return                                      Amount of mass lost from the donor in order to barely fill its own Roche Lobe
 */
double BaseBinaryStar::CalculateAdaptiveRocheLobeOverFlow(const double p_JLoss) {

    // JR: todo: should this be a program option?
    bool fixedRL = true;                                                                                                                                            // recalulate Roche Lobe radius after mass loss? (conservative vs non-conservative mass transfer)

    // record properties of the star before any changes
	// JR: todo: check this - should this really be updateDonor instead?
    double donorMass    = m_Donor->Mass();                                                                                                                          // donor mass as passed
    double accretorMass = m_Accretor->Mass();                                                                                                                       // accretor mass as passed

    double RLRadius     = m_SemiMajorAxisPrime * CalculateRocheLobeRadius_Static(donorMass, accretorMass) * AU_TO_RSOL;                                             // Roche Lobe radius - if fixedRL = true
    double jInitial     = (donorMass * accretorMass) * sqrt(G1 * (donorMass + accretorMass) * m_SemiMajorAxisPrime) / (donorMass + accretorMass);                   // initial orbital angular momentum

    // allow star to respond to previous mass loss changes
    BinaryConstituentStar* updatedDonor = new BinaryConstituentStar(*m_Donor);	                                                                                    // copy of donor star - about to be updated for mass loss

    double deltaMass = -donorMass * FAKE_MASS_LOSS_PERCENTAGE;                                                                                                      // mass loss amount
    (void)updatedDonor->UpdateAttributes(deltaMass, 0.0);                                                                                                           // apply fake mass loss: update mass - no change to mass0

    // record properties of the star before fake mass loss
    double radiusBeforeMassLoss = updatedDonor->RadiusPrev();                                                                                                       // radius before fake mass loss  JR: todo: why is this PREV, but not in CalculateMassTransferFastPhaseCaseA?
    double massBeforeMassLoss   = updatedDonor->MassPrev();                                                                                                         // mass before fake mass loss  JR: todo: why is this PREV, but not in CalculateMassTransferFastPhaseCaseA?

    SHOW_ERROR_IF(utils::Compare(radiusBeforeMassLoss, 0.0) <= 0, ERROR::RADIUS_NOT_POSITIVE_ONCE, "Before fake mass loss");                                        // show error if radius <= 0
    SHOW_ERROR_IF(utils::Compare(massBeforeMassLoss, 0.0) <= 0, ERROR::MASS_NOT_POSITIVE_ONCE, "Before fake mass loss");                                            // show error if mass <= 0

    double radiusAfterMassLoss;                                                                                                                                     // for error checking/comparison
    double massAfterMassLoss;                                                                                                                                       // what we need to calculate

    constexpr double ABSOLUTE_ERROR_THRESHOLD = 0.001;                                                                                                              // JR: todo: should this be in constants.h?
    constexpr double ITERATIONS               = 100;                                                                                                                // 100 iterations - JR: todo: should this be in constants.h?

    double absoluteErrorThreshold = std::numeric_limits<double>::max();                                                                                             // guaranteed to be >= any double - set to constant later
    double absoluteError          = std::numeric_limits<double>::max();                                                                                             // guaranteed to be >= any double - set to actual value later;

    double percentageMassLossPerIteration = OPTIONS->MaxPercentageAdaptiveMassTransfer() / ITERATIONS;                                                              // percentage mass loss per iteration
    double percentageMassLoss             = percentageMassLossPerIteration;                                                                                         // cumulative percentage mass loss (initialised for first iteration)

    for (int i = 0; i < ITERATIONS; i++) {

        BinaryConstituentStar* donorCopy = new BinaryConstituentStar(*updatedDonor);                                                                                // new temporary copy of star - set to updated donor

        if (!fixedRL) {                                                                                                                                             // calculate new Roche Lobe radius according to the mass lost
            double thisDonorMass    = (1.0 - percentageMassLoss) * donorMass;                                                                                       // new donor mass
            double thisAccretorMass = accretorMass + (m_FractionAccreted * percentageMassLoss * donorMass);                                                         // new accretor mass
            double dJ               = p_JLoss * ((1.0 - m_FractionAccreted) * (thisDonorMass - donorMass) / (donorMass + accretorMass)) * jInitial;                 // change in orbital angular momentum
            double jInitial_dJ      = jInitial + dJ;
            double aTop             = (thisDonorMass + thisAccretorMass) * (jInitial_dJ * jInitial_dJ);                                                             // new semi-major axis - numerator
            double aBottom          = G1 * thisDonorMass * thisDonorMass * thisAccretorMass * thisAccretorMass;                                                     // new semi-major axis - denominator
            double semiMajorAxis    = aTop / aBottom;                                                                                                               // new simi-major axis
            RLRadius                = semiMajorAxis * CalculateRocheLobeRadius_Static(thisDonorMass, thisAccretorMass) * AU_TO_RSOL;                                // new Roche Lobe radius
        }

        // Reduce mass of star by percentageMassLoss and recalculate the radius
        double deltaMass  = -donorCopy->Mass() * percentageMassLoss;                                                                                                // mass loss amount
        double deltaMass0 = -donorCopy->Mass0() * percentageMassLoss;                                                                                               // mass loss amount
        (void)donorCopy->UpdateAttributes(deltaMass, deltaMass0);                                                                                                   // apply fake mass loss: update mass and mass0

        double thisMassAfterMassLoss = donorCopy->Mass();                                                                                                           // record mass loss for donor

        // Modify donor Mass0 and Age for MS (including HeMS) and HG stars
        donorCopy->UpdateInitialMass();                                                                                                                             // update initial mass (MS, HG & HeMS)  JR: todo: fix this kludge - mass0 is overloaded, and isn't always "initial mass"
        donorCopy->UpdateAgeAfterMassLoss();                                                                                                                        // update age (MS, HG & HeMS)

        (void)donorCopy->AgeOneTimestep(0.0);                                                                                                                       // recalculate radius of star - don't age - just update values

        double thisRadiusAfterMassLoss = donorCopy->Radius();                                                                                                       // record new radius for donor
        double thisAbsoluteError       = std::abs(thisRadiusAfterMassLoss - RLRadius);                                                                              // absolute error in radius

        // record values for lowest relative error
        if (utils::Compare(thisAbsoluteError, absoluteErrorThreshold) <= 0) {                                                                                       // absolute error below threshold?  (will always be on first iteration)
            absoluteErrorThreshold = ABSOLUTE_ERROR_THRESHOLD;                                                                                                      // set threshold to actual constant value
            if (utils::Compare(thisAbsoluteError, absoluteError) < 0) {                                                                                             // this iteration absolute error < lowest absolute error?  (will always be on first iteration)
                absoluteError       = thisAbsoluteError;                                                                                                            // set lowest absolute error
                massAfterMassLoss   = thisMassAfterMassLoss;                                                                                                        // lowest error mass after mass loss
                radiusAfterMassLoss = thisRadiusAfterMassLoss;                                                                                                      // lowest error radius after mass loss
            }
        }

        percentageMassLoss += percentageMassLossPerIteration;                                                                                                       // increment percentage mass loss for next iteration

        SHOW_ERROR_IF(utils::Compare(radiusAfterMassLoss, radiusBeforeMassLoss) > 0, ERROR::INVALID_RADIUS_INCREASE_ONCE, "After fake mass loss");                  // show error if radius increased

        delete donorCopy; donorCopy = nullptr;
    }

    delete updatedDonor; updatedDonor = nullptr;

    SHOW_ERROR_IF(utils::Compare(radiusAfterMassLoss, 0.0) <= 0, ERROR::RADIUS_NOT_POSITIVE_ONCE, "After fake mass loss");                                          // show error if updated radius <= 0
    SHOW_ERROR_IF(utils::Compare(massAfterMassLoss, 0.0) <= 0, ERROR::MASS_NOT_POSITIVE_ONCE, "After fake mass loss");                                              // show error if updated mass <= 0

    return massAfterMassLoss - massBeforeMassLoss;                                                                                                                  // return change in mass
}


/*
 * Calculate the amount of mass that a star with radiative envelope needs to lose in order to just fill its Roche Lobe.
 *
 * Based on the function –rlof_method in binary_c.      JR: todo: reference?
 *
 * For fast phase case A MT, we solve the orbit numerically for the thermal timescale of the donor
 *
 *
 * double CalculateMassTransferFastPhaseCaseA(const double p_JLoss)
 * @param   [IN]    p_JLoss                     Specific angular momentum with which mass is lost during non-conservative mass transfer
 * @return                                      Amount of mass lost from the donor in order to barely fill its own Roche Lobe
 */
double BaseBinaryStar::CalculateMassTransferFastPhaseCaseA(const double p_JLoss) {

    // record properties of the star before any changes
	// JR: todo: check this - should this really be updateDonor instead?
    double donorMass    = m_Donor->Mass();                                                                                                              // donor mass
    double accretorMass = m_Accretor->Mass();                                                                                                           // accretor mass

    double dt = m_Donor->CalculateThermalTimescale();                                                                                                   // dt is donor's thermal timescale

    // allow star to respond to previous mass loss changes
    BinaryConstituentStar* updatedDonor = new BinaryConstituentStar(*m_Donor);	                                                                        // copy of donor star - about to be updated for mass loss

    double deltaMass = -donorMass * FAKE_MASS_LOSS_PERCENTAGE;                                                                                          // mass loss amount
    (void)updatedDonor->UpdateAttributes(deltaMass, 0.0);                                                                                               // apply fake mass loss: update mass - no change to mass0

    // record properties of the star before fake mass loss
    double radiusBeforeMassLoss = updatedDonor->Radius();                                                                                               // radius before fake mass loss (but after previous mass loss)
    double massBeforeMassLoss   = updatedDonor->Mass();                                                                                                 // mass before fake mass loss (but after previous mass loss)

    SHOW_ERROR_IF(utils::Compare(radiusBeforeMassLoss, 0.0) <= 0, ERROR::RADIUS_NOT_POSITIVE_ONCE, "Before fake mass loss");                            // show error if radius <= 0
    SHOW_ERROR_IF(utils::Compare(massBeforeMassLoss, 0.0) <= 0, ERROR::MASS_NOT_POSITIVE_ONCE, "Before fake mass loss");                                // show error if mass <= 0

    double radiusAfterMassLoss;                                                                                                                         // for error checking/comparison
    double massAfterMassLoss;                                                                                                                           // what we need to calculate

    constexpr double RELATIVE_ERROR_THRESHOLD = 0.1;                                                                                                    // JR: todo: should this be in constants.h?
    constexpr double MAX_PERCENTAGE_MASS_LOSS = 0.99;                                                                                                   // In fast phase case A mass transfer, we allow to strip almost all of the star in order to see if it fits in it's Roche lobe         JR: todo: should this be in constants.h?
    constexpr double ITERATIONS               = 100;                                                                                                    // 100 iterations - JR: todo: should this be in constants.h?

    double relativeErrorThreshold = std::numeric_limits<double>::max();                                                                                 // guaranteed to be >= any double - set to constant later
    double relativeError          = std::numeric_limits<double>::max();                                                                                 // guaranteed to be >= any double - set to actual value later;

    double percentageMassLossPerIteration = MAX_PERCENTAGE_MASS_LOSS / ITERATIONS;                                                                      // percentage mass loss per iteration
    double percentageMassLoss             = percentageMassLossPerIteration;                                                                             // cumulative percentage mass loss (initialised for first iteration)

    for (int i = 0; i < ITERATIONS; i++) {

        BinaryConstituentStar* donorCopy = new BinaryConstituentStar(*updatedDonor);                                                                    // new temporary copy of star - set to updated donor

        double dM		     = percentageMassLoss * donorMass;                                                                                          // change in donor mass

        double semiMajorAxis = CalculateMassTransferOrbit(*updatedDonor, *m_Accretor, -dM / dt, dt, p_JLoss);                                           // new semi major axis     JR: todo: check donor values ok (before/after one time step...)
        double RLRadius      = semiMajorAxis * CalculateRocheLobeRadius_Static(donorMass - dM, accretorMass + (m_FractionAccreted * dM)) * AU_TO_RSOL;  // calculate new Roche Lobe Radius according to the mass lost

        // Reduce mass of star by percentageMassLoss and recalculate the radius
        double deltaMass  = -donorCopy->Mass() * percentageMassLoss;                                                                                    // mass loss amount
        double deltaMass0 = -donorCopy->Mass0() * percentageMassLoss;                                                                                   // mass loss amount
        (void)donorCopy->UpdateAttributes(deltaMass, deltaMass0);                                                                                       // apply fake mass loss: update mass and mass0

        double thisMassAfterMassLoss = donorCopy->Mass();                                                                                               // record mass loss for donor
        // Modify donor Mass0 and Age for MS (including HeMS) and HG stars
        donorCopy->UpdateInitialMass();                                                                                                                 // update initial mass (MS, HG & HeMS)  JR: todo: fix this kludge - mass0 is overloaded, and isn't always "initial mass"
        donorCopy->UpdateAgeAfterMassLoss();                                                                                                            // update age (MS, HG & HeMS)

        (void)donorCopy->AgeOneTimestep(0.0);                                                                                                           // recalculate radius of star - don't age - just update values

        double thisRadiusAfterMassLoss = donorCopy->Radius();                                                                                           // record new radius for donor
        double thisRelativeError       = std::abs(thisRadiusAfterMassLoss - RLRadius) / RLRadius;                                                       // relative error in radius

        // record values for lowest relative error
        if (utils::Compare(thisRelativeError, relativeErrorThreshold) <= 0) {                                                                           // relative error below threshold?  (will always be on first iteration)
            relativeErrorThreshold = RELATIVE_ERROR_THRESHOLD;                                                                                          // set threshold to actual constant value
            if (utils::Compare(thisRelativeError, relativeError) < 0) {                                                                                 // this iteration relative error < lowest relative error?  (will always be on first iteration)
                relativeError       = thisRelativeError;                                                                                                // set lowest relative error
                massAfterMassLoss   = thisMassAfterMassLoss;                                                                                            // lowest error mass after mass loss
                radiusAfterMassLoss = thisRadiusAfterMassLoss;                                                                                          // lowest error radius after mass loss
            }
        }

        percentageMassLoss += percentageMassLossPerIteration;                                                                                           // increment percentage mass loss for next iteration

        SHOW_ERROR_IF(utils::Compare(radiusAfterMassLoss, radiusBeforeMassLoss) > 0, ERROR::INVALID_RADIUS_INCREASE_ONCE, "After fake mass loss");      // show error if radius increased

        delete donorCopy; donorCopy = nullptr;
    }

    delete updatedDonor; updatedDonor = nullptr;

    SHOW_ERROR_IF(utils::Compare(radiusAfterMassLoss, 0.0) <= 0, ERROR::RADIUS_NOT_POSITIVE_ONCE, "After fake mass loss");                              // show error if updated radius <= 0
    SHOW_ERROR_IF(utils::Compare(massAfterMassLoss, 0.0) <= 0, ERROR::MASS_NOT_POSITIVE_ONCE, "After fake mass loss");                                  // show error if updated mass <= 0

    // JR: todo: is this still needed?
	// ALEJANDRO - 19/01/2017 - Following error message may not classify as an error per se, but want
	// to check it in the error file to see if there is any correlation with the error in the slow phase.
    SHOW_WARN_IF(utils::Compare(massAfterMassLoss, massBeforeMassLoss * MAX_PERCENTAGE_MASS_LOSS) == 0, ERROR::MAXIMUM_MASS_LOST);

    return massAfterMassLoss - massBeforeMassLoss;                                                                                                      // return change in mass
}


/*
 * Calculate the fraction of specific angular momentum with which the non-accreted mass leaves the system
 *
 * This is gamma (as in Pols's notes) or jloss (as in Belczynski et al. 2008
 * which is the fraction of specific angular momentum with which the non-accreted mass leaves the system.
 *
 * Calculation is based on user-specified Angular Momentum Loss prescription
 *
 *
 * double CalculateGammaAngularMomentumLoss(const double p_DonorMass, const double p_AccretorMass)
 *
 * @param   [IN]    p_DonorMass                 The mass of the donor (Msol)
 * @param   [IN]    p_JLoss                     The mass of the accretor (Msol)
 * @return                                      The fraction of specific angular momentum with which the non-accreted mass leaves the system
 */
double BaseBinaryStar::CalculateGammaAngularMomentumLoss(const double p_DonorMass, const double p_AccretorMass) {

	double gamma;

	switch (OPTIONS->MassTransferAngularMomentumLossPrescription()) {                                                       // which precription?
        case MT_ANGULAR_MOMENTUM_LOSS_PRESCRIPTION::JEANS                : gamma = p_AccretorMass / p_DonorMass; break;
        case MT_ANGULAR_MOMENTUM_LOSS_PRESCRIPTION::ISOTROPIC_RE_EMISSION: gamma = p_DonorMass / p_AccretorMass; break;
        case MT_ANGULAR_MOMENTUM_LOSS_PRESCRIPTION::CIRCUMBINARY_RING    : gamma = (M_SQRT2 * (p_DonorMass + p_AccretorMass) * (p_DonorMass + p_AccretorMass)) / (p_DonorMass * p_AccretorMass); break; // Based on the assumption that a_ring ~= a*sqrt(2) in Evernote Notebook based on talks with deMink, or as tricky people call it, "private communication"
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
 * Belczynski et al. 2008, eq 32, 33, de Mink ????          JR: todo: find deMink reference
 *
 * Calculation based on user-specified Mass Loss prescription
 *
 *
 * double CalculateMassTransferOrbit(const double p_MDotDonor, const double p_Dt, const double p_Jloss)
 *
 * @param   [IN]    p_MDotDonor                 Mass loss rate of the donor
 * @param   [IN]    p_Dt                        Timestep
 * @param   [IN]    p_JLoss                     Specific angular momentum with which mass is lost during non-conservative mass transfer
 * @return                                      Semi-major axis
 */
double BaseBinaryStar::CalculateMassTransferOrbit(BinaryConstituentStar& p_Donor, BinaryConstituentStar& p_Accretor, const double p_MDotDonor, const double p_Dt, const double p_Jloss) {

    double semiMajorAxis   = m_SemiMajorAxisPrime;                                                              // new semi-major axis value - default is no change
    double massA           = p_Accretor.Mass();                                                                 // accretor mass
    double massD           = p_Donor.Mass();                                                                    // donor mass
    double massAtimesMassD = massA * massD;                                                                     // accretor mass * donor mass
    double massAplusMassD  = massA + massD;                                                                     // accretor mass + donor mass
    double jOrb            = (massAtimesMassD / massAplusMassD) * sqrt(semiMajorAxis * G1 * massAplusMassD);    // orbital angular momentum
    double jLoss           = p_Jloss;                                                                           // specific angular momentum (for non-conservative mass transfer)

    // determine which mass loss prescription to use:
    // default is the prescription specified by the program options,
    // but we also choose to use BELCZYNSKI if the donor has a RADIATIVE envelope

    MT_PRESCRIPTION prescription = p_Donor.DetermineEnvelopeType() == ENVELOPE::RADIATIVE ? MT_PRESCRIPTION::BELCZYNSKI : OPTIONS->MassTransferPrescription();

    // calculate new semi-major axis value using the chosen prescription
    switch (prescription) {                                                                                     // which mass transfer prescription?

        case MT_PRESCRIPTION::DEMINK: {                                                                         // using de Mink mass transfer prescription
                                                                                                                // degenerate and non-degenerate accretor solutions same for de Mink
            double jPrime;
            double aPrime;

            double fractionAccreted    = m_FractionAccreted;
            double thermalRateAccretor = p_Accretor.CalculateThermalMassLossRate();
            double thermalRateDonor    = p_Donor.CalculateThermalMassLossRate();
            double dt                  = p_Dt / DEMINK_ORBIT_ITERATIONS;                                        // delta t per iteration

            for(int i = 0; i < DEMINK_ORBIT_ITERATIONS ; i++) {

                jPrime = jOrb + ((jLoss * jOrb * (1.0 - fractionAccreted) * p_MDotDonor / massAplusMassD) * dt);
                aPrime = semiMajorAxis + (((-2.0 * (p_MDotDonor / massD)) * (1.0 - (fractionAccreted * (massD / massA)) - ((1.0 - fractionAccreted) * (jLoss + 0.5) * (massD / massAplusMassD)))) * semiMajorAxis * dt);

                // JR: todo: check convergence and stop early?
                // If we're not going to check for convergence then we don't really need jPrime or aPrime - just assign jOrb and semiMajorAxis directly
                jOrb           = jPrime;
                semiMajorAxis  = aPrime;

                massD          = massD + (p_MDotDonor * dt);
                massA          = massA - (p_MDotDonor * dt * fractionAccreted);
                massAplusMassD = massA + massD;

                jLoss = CalculateGammaAngularMomentumLoss(massD, massA);

                std::tie(std::ignore, fractionAccreted) = p_Accretor.CalculateMassAcceptanceRate(thermalRateDonor, fractionAccreted, thermalRateAccretor);
            }

            double mPrime    = massD + massA;
            double jPrime_mm = jPrime / (massD * massA);

            semiMajorAxis    = jPrime_mm * jPrime_mm * mPrime / G1;                                             // change in the orbit due to angular momentum loss from MT

            } break;

        case MT_PRESCRIPTION::BELCZYNSKI: {                                                                     // using Belczynski mass transfer prescription

            if (p_Accretor.IsDegenerate()) {                                                                    // degenerate accretor - Belczynski et al. 2008, eq 32

                // ALEJANDRO - 29/11/2016 -
                // Disabled this for RADIATIVE ENVELOPE donor as it didn't allow case A MT
                // to solve properly for the orbit, specially in the fast phase.
                // This should be revisited in order to make it efficient.
                //
                // JR: todo: check this

                if (p_Donor.DetermineEnvelopeType() != ENVELOPE::RADIATIVE) {

                    double rCom   = (semiMajorAxis * massD) / massAplusMassD;                                   // distance from the compact object to the centre of mass
                    double jPrime = jOrb + ((rCom * rCom) * m_OrbitalVelocityPrime * (1.0 - m_FractionAccreted) * p_MDotDonor * p_Dt);

                    double massDprime = (p_MDotDonor * p_Dt ) + massD;                                          // mass of the donor after mass transfer
                    double massAprime = (-m_FractionAccreted * p_MDotDonor * p_Dt) + massA;                     // mass of the accretor after mass transfer
                    double mPrime     = massDprime + massAprime;
                    double mAmD_Prime = massDprime * massAprime;

                    semiMajorAxis      = jPrime * jPrime * mPrime / (G1 * mAmD_Prime * mAmD_Prime);             // change in the orbit due to angular momentum loss from MT
                }
            }
            else {                                                                                              // non-degenerate accretor - Belczynski et al. 2008, eq 33

                double massDprime = (p_MDotDonor * p_Dt) + massD;                                               // mass of the donor after mass transfer
                double massAprime = (-m_FractionAccreted * p_MDotDonor * p_Dt) + massA;                         // mass of the accretor after mass transfer
                double mPrime     = massDprime + massAprime;

                double jPrime     = jOrb + ((jLoss * jOrb * (1.0 - m_FractionAccreted) * p_MDotDonor / massAplusMassD) * p_Dt);
                double jPrime_mm  = jPrime / (massDprime * massAprime);

                semiMajorAxis     = jPrime_mm * jPrime_mm * mPrime / G1;                                        // change in the orbit due to angular momentum loss from MT
            }
            } break;

        case MT_PRESCRIPTION::NONE:                                                                             // unknown prescription
        default:                                                                                                // unknown prescription
            SHOW_WARN(ERROR::UNKNOWN_MT_PRESCRIPTION);                                                          // warn that an error occurred
    }

    return semiMajorAxis;
}


/*
 * Calculate the response of the donor Roche Lobe to mass loss during mass transfer per Belczynski et al. 2008
 *
 * Numerical calculation of the Roche Lobe after mass transfer as in StarTrack.
 * Described in Belczynski et al. 2008. Used for a regular star accretor and non-conservative Mass Transfer.
 * Belczynsky et al. (2008), eq 41, using Woods et al. (2012) formula
 *
 *    dJ=Beta*((1.0-Fa)*(Ma2-Ma1)/(Ma1+Mb1))*Jorb1;
 *    a2=((Ma2+Mb2)*pow(Jorb1+dJ,2.0))/(GGG*Ma2*Ma2*Mb2*Mb2);  (non-conservative MT assumption)
 *
 * JR: todo: What does "(Numerical) ZRocheLobe" mean?  Why don't we call this function "CalculateRocheLobResponseToMasslossBelczynski" (or something)?
 *
 *
 * double CalculateNumericalZRocheLobe(const double p_jLoss)
 *
 * @param   [IN]    p_JLoss                     Specific angular momentum with which mass is lost during non-conservative mass transfer
 *                                              (Podsiadlowski et al. 1992, Beta: specific angular momentum of matter [2Pia^2/P])
 * @return                                      Roche Lobe response
 */
double BaseBinaryStar::CalculateNumericalZRocheLobe(const double p_jLoss) {

    double initialDonorMass    = m_Donor->Mass();                                                                                           // donor mass before mass transfer
    double initialAccretorMass = m_Accretor->Mass();                                                                                        // accretor mass before mass transfer

    double initialMassAplusMassD = initialAccretorMass + initialDonorMass;                                                                  // accretor mass + donor mass before mass transfer

    double donorMass    = (1.0 - RL_MASS_LOSS_FRACTION) * initialDonorMass;                                                                 // donor mass after mass transfer
    double accretorMass = initialAccretorMass + m_FractionAccreted * RL_MASS_LOSS_FRACTION * initialDonorMass;                              // accretor mass after mass transfer

    double J  = initialAccretorMass * initialDonorMass * sqrt(G1 * initialMassAplusMassD * m_SemiMajorAxisPrime) / initialMassAplusMassD;   // angular momentum before mass transfer
           J += J * p_jLoss * ((1.0 - m_FractionAccreted) * (donorMass - initialDonorMass) / initialMassAplusMassD);                        // angular momentum after mass transfer

    double semiMajorAxis = (donorMass + accretorMass) * (J * J) / (G1 * donorMass * donorMass * accretorMass * accretorMass);               // semi major axis adjusted for mass transfer

    double RLRadiusBefore = m_SemiMajorAxisPrime * CalculateRocheLobeRadius_Static(initialDonorMass, initialAccretorMass);                  // Roche Lobe radius before mass transfer
    double RLRadiusAfter  = semiMajorAxis * CalculateRocheLobeRadius_Static(donorMass, accretorMass);                                       // Roche Lobe radius after mass transfer

    return (log(RLRadiusAfter) - log(RLRadiusBefore)) / (log(donorMass) - log(initialDonorMass));                                           // Response of the donor Roche Lobe to mass loss during mass transfer
}


/*
 * Calculate the response of the donor Roche Lobe to mass loss during mass transfer per Sluys 2013, Woods et al., 2012
 *
 * Sluys 2013, eq 60, Woods et al., 2012
 * Formula from M. Sluys notes "Binary evolution in a nutshell"
 *
 * JR: todo: What does "ZRocheLobe" mean?  Why don't we call this function "CalculateRocheLobResponseToMasslossSluys" (or something)?
 *
 *
 * double CalculateZRocheLobe()
 */
double BaseBinaryStar::CalculateZRocheLobe() {

    double donorMass    = m_Donor->Mass();                  // donor mass
    double accretorMass = m_Accretor->Mass();               // accretor mass
    double beta         = m_FractionAccreted;               // fraction of mass accreted by accretor

    double q = donorMass / accretorMass;

    double massAplusMassD = accretorMass + donorMass;
    double massAtimeMassD = accretorMass * donorMass;

    double q_1_3 = pow(q, 1.0 / 3.0);

    double k1 = ((2.0 * donorMass * donorMass)-(2.0 * accretorMass * accretorMass) - (massAtimeMassD * (1.0 - beta))) / (accretorMass * massAplusMassD);
    double k2 = (2.0 / 3.0) - ((q_1_3 * ((1.2 * q_1_3) + (1.0 / (1.0 + q_1_3)))) / (3.0 * ((0.6 * pow(q, 2.0 / 3.0)) + (log(1.0 + q_1_3)))));
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
    m_OmegaMassLossDiff = 0.0;                                                                                                  // initially - no change to orbital speed due to winds mass loss

    if (OPTIONS->UseMassTransfer() && m_MassTransfer) {                                                                         // used for halting winds when in mass transfer (first approach).
            m_Star1->SetMassLossDiff(0.0);                                                                                      // JR: todo: find a better way?
            m_Star2->SetMassLossDiff(0.0);                                                                                      // JR: todo: find a better way?
    }
    else {
        if (OPTIONS->UseMassLoss()) {                                                                                           // mass loss enabled?

            double mWinds1 = m_Star1->CalculateMassLossValues(true);                                                            // calculate new values assuming mass loss applied
            double mWinds2 = m_Star2->CalculateMassLossValues(true);                                                            // calculate new values assuming mass loss applied

            double aWinds = m_SemiMajorAxisPrev / (2.0 - ((m_Star1->MassPrev() + m_Star2->MassPrev()) / (mWinds1 + mWinds2)));  // new semi-major axis for circularlised orbit

            m_Star1->SetMassLossDiff(mWinds1 - m_Star1->Mass());                                                                // JR: todo: find a better way?
            m_Star2->SetMassLossDiff(mWinds2 - m_Star2->Mass());                                                                // JR: todo: find a better way?

            m_aMassLossDiff     = aWinds - m_SemiMajorAxisPrev;                                                                 // change to orbit (semi-major axis) due to winds mass loss
            m_OmegaMassLossDiff = sqrt(G1 * (mWinds1 + mWinds2) / (aWinds * aWinds * aWinds)) - m_OrbitalVelocityPrev;          // change to orbital speed due to winds mass loss
        }
    }
}


/*
 * Calculate mass transfer and update individual stars accordingly
 *
 * JR: todo: flesh-out this documentation
 *
 *
 * void CalculateMassTransfer(const double p_Dt)
 *
 * @param   [IN]    p_Dt                        timestep in Myr
 */
void BaseBinaryStar::CalculateMassTransfer(const double p_Dt) {

    if (!OPTIONS->UseMassTransfer()) return;                                                                                                                // mass transfer not enabled - nothing to do

    if (!m_Star1->IsRLOF() && !m_Star2->IsRLOF()) return;                                                                                                   // neither star is overflowing its Roche Lobe - no mass transfer - nothing to do

    if (m_Star1->IsRLOF() && m_Star2->IsRLOF()) {                                                                                                           // both stars overflowing their Roche Lobe?
        m_CEDetails.CEEnow = true;                                                                                                                          // yes - common envelope event - no mass transfer
        return;                                                                                                                                             // and return - nothing (else) to do
    }

    // one, and only one, star is overflowing its Roche Lobe - resolve mass transfer

    m_Donor    = m_Star2->IsRLOF() ? m_Star2 : m_Star1;                                                                                                     // donor is primary unless secondary is overflowing its Roche Lobe
    m_Accretor = m_Star2->IsRLOF() ? m_Star1 : m_Star2;                                                                                                     // accretor is secondary unless secondary is overflowing its Roche Lobe

    m_Donor->BecomePrimary();                                                                                                                               // tell the donor it is the primary
    m_Accretor->BecomeSecondary();                                                                                                                          // tell the accretor it is not the primary

    double aInitial = m_SemiMajorAxisPrime;                                                                                                                 // semi-major axis in default units, AU, current timestep
    double aFinal;                                                                                                                                          // semi-major axis in default units, AU, after next timestep
    double wInitial = m_OrbitalVelocityPrime;                                                                                                               // initial orbital speed of the system in yr-1 units, current timestep
    double wFinal;                                                                                                                                          // initial orbital speed of the system in yr-1 units, after next timestep
    double jLoss    = m_JLoss;                            		                                                                                            // specific angular momentum with which mass is lost during non-conservative mass transfer, current timestep
	bool   isCEE    = false;									                                                                                            // is there a CEE in this MT episode?

    //Addition Coen 18-10-2017
    m_Donor->CalculateZetas();                                                                                                                              // calculate Zetas for donor
    m_Accretor->CalculateZetas();                                                                                                                           // calculate Zetas for accretor

	// Check for stability
	bool qCritFlag = OPTIONS->MassTransferCriticalMassRatioMSLowMass()   || OPTIONS->MassTransferCriticalMassRatioMSHighMass()  ||
	                 OPTIONS->MassTransferCriticalMassRatioHG()          || OPTIONS->MassTransferCriticalMassRatioGiant()       ||
	                 OPTIONS->MassTransferCriticalMassRatioHeliumGiant() || OPTIONS->MassTransferCriticalMassRatioHeliumMS()    ||
                     OPTIONS->MassTransferCriticalMassRatioHeliumHG()    || OPTIONS->MassTransferCriticalMassRatioHeliumGiant() ||
                     OPTIONS->MassTransferCriticalMassRatioWhiteDwarf();

    if (m_Donor->IsMassRatioUnstable(m_Accretor->Mass(), m_Accretor->IsDegenerate()) && qCritFlag) {
        m_CEDetails.CEEnow = true;
    }
    else {

        m_Donor->DetermineInitialMassTransferCase();                                                                                                        // record first mass transfer event type (case A, B or C)

		// Begin Mass Transfer
        switch (OPTIONS->MassTransferPrescription()) {                                                                                                      // which mass transfer prescription?

            case MT_PRESCRIPTION::DEMINK: {                                                                                                                 // de Mink

                switch (m_Donor->DetermineMassTransferCase()) {                                                                                             // which MT case?

                    case MT_CASE::A:                                                                                                                        // A, or
                    case MT_CASE::B: {                                                                                                                       // B

                        double thermalRateDonor    = m_Donor->CalculateThermalMassLossRate();
                        double thermalRateAccretor = OPTIONS->MassTransferThermallyLimitedVariation() == MT_THERMALLY_LIMITED_VARIATION::RADIUS_TO_ROCHELOBE
                                                        ? m_Accretor->Mass() / m_Accretor->CalculateThermalTimescale(m_Accretor->Mass(), m_Accretor->RocheLobeRadius() * AU_TO_RSOL, m_Accretor->Luminosity(), m_Accretor->Mass() - m_Accretor->CoreMass()) // assume Radius = RL
                                                        : m_Accretor->CalculateThermalMassLossRate();



                        std::tie(std::ignore, m_FractionAccreted) = m_Accretor->CalculateMassAcceptanceRate(thermalRateDonor, m_FractionAccreted, thermalRateAccretor);

                        } break;

                    case MT_CASE::C:                                                                                                                        // C
                        std::tie(std::ignore, m_FractionAccreted) = m_Accretor->CalculateMassAcceptanceRate(m_Donor->CalculateDynamicalMassLossRate(),
                                                                                                            m_FractionAccreted,
                                                                                                            m_Accretor->CalculateDynamicalMassLossRate());
                        break;

                    default:                                                                                                                                // unknown MT_CASE
                        SHOW_ERROR(ERROR::UNKNOWN_MT_CASE);                                                                                                 // show error
                }


                if (OPTIONS->MassTransferAngularMomentumLossPrescription() != MT_ANGULAR_MOMENTUM_LOSS_PRESCRIPTION::ARBITRARY) {                           // arbitray angular momentum loss prescription?
                    jLoss = CalculateGammaAngularMomentumLoss();                                                                                            // no - re-calculate angular momentum
                }

                double ZlobAna      = CalculateZRocheLobe();
                m_ZetaRLOFNumerical = CalculateNumericalZRocheLobe(jLoss);
                m_ZetaRLOFAnalytic  = ZlobAna;                                                                                                              // addition by Coen 18-10-2017 for zeta study.  ALEJANDRO - 04/10/2017 - Moved this to calculate it for deMink MT.

                switch (m_Donor->DetermineEnvelopeType()) {                                                                                                 // which enveleope type?

                    case ENVELOPE::RADIATIVE: {                                                                                                             // RADIATIVE: case A

                        // Need to know which is the donor star, make it lose enough mass to stay within its Roche lobe

                        // ALEJANDRO - 20/10/2017 - Code arbitrary zeta cut for case A mass transfer. To use in BNS paper. Should be properly coded.
                        // JR: todo: is this still required?  If so, find out what "properley coded" means...
                        if (m_Donor->IsOneOf({ STELLAR_TYPE::MS_GT_07 })) {
                            if (utils::Compare(OPTIONS->ZetaMainSequence(), ZlobAna) < 0) {
                                m_StellarMerger = true;
                                isCEE           = true;
                            }
                        }

                        m_MassTransferTrackerHistory = m_Donor->IsPrimary() ? MT_TRACKING::STABLE_FROM_1_TO_2 : MT_TRACKING::STABLE_FROM_2_TO_1;            // record what happened - for later printing

                        double dM    = m_Donor->FastPhaseCaseA() ? CalculateAdaptiveRocheLobeOverFlow(jLoss) : CalculateMassTransferFastPhaseCaseA(jLoss);  // amount of mass transferred by donor
                        double MdDot = dM / p_Dt;                                                                                                           // mass transfer rate of donor

                        m_Donor->SetFastPhaseCaseA();                                                                                                       // will be true when we get here   JR: todo: revisit this
                        m_Donor->SetMassTransferDiff(dM);                                                                                                   // mass transferred by donor
                        m_Accretor->SetMassTransferDiff(-dM * m_FractionAccreted);                                                                          // mass accreted by accretor

                        aFinal = CalculateMassTransferOrbit(*m_Donor, *m_Accretor, MdDot, p_Dt, jLoss);                                                     // before solving for orbit
                        wFinal = sqrt(G1 * (m_Donor->Mass() + m_Accretor->Mass()) / (aFinal * aFinal * aFinal));                                            // after solving for orbit

                        m_aMassTransferDiff     = aFinal - aInitial;                                                                                        // change in orbit (semi-major axis)
                        m_OmegaMassTransferDiff = wFinal - wInitial;                                                                                        // change in orbital speed

                        } break;

                    case ENVELOPE::CONVECTIVE: {                                                                                                            // CONVECTIVE: case B or case C

                        BinaryConstituentStar* donorCopy = new BinaryConstituentStar(*m_Donor);
                        BinaryConstituentStar* accretorCopy = new BinaryConstituentStar(*m_Accretor);

                        double zetaCompare = m_Donor->CalculateZeta(OPTIONS->CommonEnvelopeZetaPrescription());

                        if (OPTIONS->ForceCaseBBBCStabilityFlag()) {

                            // ALEJANDRO - 24/08/2017 - Check for case BB or BC mass transfer; particularly for BNS project
                            if (m_Donor->IsOneOf({ STELLAR_TYPE::NAKED_HELIUM_STAR_HERTZSPRUNG_GAP, STELLAR_TYPE::NAKED_HELIUM_STAR_GIANT_BRANCH })) {

                                if (m_Accretor->IsOneOf({ STELLAR_TYPE::NEUTRON_STAR })) {
                                    m_Donor->SetSNCurrentEvent(SN_EVENT::USSN);                                                                             // donor ultra-stripped SN
                                }                                                                                                                           // JR: todo: check "else false"

                                // Hard code stability
                                zetaCompare = OPTIONS->AlwaysStableCaseBBBCFlag() ? 1.0 : 0.0;
                                ZlobAna		= OPTIONS->AlwaysStableCaseBBBCFlag() ? 0.0 : 1.0;
                            }
                        }

                        // ALEJANDRO - 20/10/2017 - Code arbitrary zeta cut for case B mass transfer. To use in BNS paper. Should be properly coded.
                        if (m_Donor->IsOneOf({ STELLAR_TYPE::HERTZSPRUNG_GAP })) {
                            zetaCompare = OPTIONS->ZetaHertzsprungGap();
                        }

                        // ALEJANDRO - 07/08/2018 - Added m_zetaStarCompare for CE study.       JR: todo: is this still required?
                        m_ZetaStarCompare = zetaCompare;                                                                                                    // beware as this variable may have different values, e.g. fixed for MS, fixed for HG, different for Soberman, '1' for case BB.

                        if (zetaCompare > ZlobAna) {                                                                                                        // stable Mass Transfer

                            // Check for conservative or non-conservative MT
                            double envMassDonor  = donorCopy->Mass() - donorCopy->CoreMass();
                            double mdEnvAccreted = envMassDonor * m_FractionAccreted;
                            double MdDot         = -envMassDonor / p_Dt;

                            m_Donor->SetMassTransferDiff(-envMassDonor);
                            m_Accretor->SetMassTransferDiff(mdEnvAccreted);

                            m_Donor->ResolveEnvelopeLossAndSwitch();                                                                                        // only other interaction that adds/removes mass is winds. So, think its safe to update star here.

                            m_MassTransferTrackerHistory = m_Donor->IsPrimary() ? MT_TRACKING::STABLE_FROM_1_TO_2 : MT_TRACKING::STABLE_FROM_2_TO_1;

                            // Update and solve the orbit
                            aFinal                  = CalculateMassTransferOrbit(*donorCopy, *accretorCopy, MdDot, p_Dt, jLoss);
                            wFinal                  = sqrt(G1 * (donorCopy->Mass() + accretorCopy->Mass()) / (aFinal * aFinal * aFinal));
                            m_aMassTransferDiff     = aFinal - aInitial;
                            m_OmegaMassTransferDiff = wFinal - wInitial;

                            // Check for stable mass transfer after any CEE
                            if (m_CEDetails.CEEcount > 0 && !m_RLOFDetails.stableRLOFPostCEE) {
                                m_RLOFDetails.stableRLOFPostCEE = m_MassTransferTrackerHistory == MT_TRACKING::STABLE_FROM_2_TO_1 ||
                                                                  m_MassTransferTrackerHistory == MT_TRACKING::STABLE_FROM_1_TO_2;
                            }
                        }
                        else {                                                                                                                              // Unstable Mass Transfer
                            if (m_Donor->IsOneOf( MAIN_SEQUENCE )) {                                                                                        // How to deal with CEE here? Just worry about donor?
                                m_StellarMerger    = true;
                                isCEE              = true;
                            }
                            else {
                                m_CEDetails.CEEnow = true;
                                isCEE              = true;
                            }
                        }

                        delete donorCopy; donorCopy = nullptr;

                        } break;

                    case ENVELOPE::REMNANT:                                                                                                                 // invalid - donor is NS, BH or Massless remnant
                        m_Error = ERROR::INVALID_MASS_TRANSFER_DONOR;                                                                                       // set error value
                        SHOW_WARN(m_Error);                                                                                                                 // warn that an error occurred
                        break;

                    default:
                        m_Error = ERROR::INVALID_ENVELOPE_TYPE;                                                                                             // set error value
                        SHOW_WARN(m_Error);                                                                                                                 // warn that an error occurred
                }
                } break;

            case MT_PRESCRIPTION::BELCZYNSKI:                                                                                                               // Belczynski - not yet supported
                m_Error = ERROR::UNSUPPORTED_MT_PRESCRIPTION;                                                                                               // set error value
                SHOW_WARN(m_Error);                                                                                                                         // warn that an error occurred
                break;

            default:                                                                                                                                        // unknown mass transfer prescription - shouldn't happen
                m_Error = ERROR::UNKNOWN_MT_PRESCRIPTION;                                                                                                   // set error value
                SHOW_WARN(m_Error);                                                                                                                         // warn that an error occurred
        }
    }

	// Check for recycled pulsars. Not considering CEE as a way of recycling NSs.
	if (!isCEE && m_Accretor->IsOneOf({ STELLAR_TYPE::NEUTRON_STAR })) {                                                                                    // accretor is a neutron star
            m_Donor->SetSNPastEvent(SN_EVENT::RLOF_ONTO_NS);                                                                                                // donor donated mass to a neutron star
            m_Accretor->SetSNPastEvent(SN_EVENT::RECYCLED_NS);                                                                                              // accretor is (was) a recycled NS
	}
}


/*
 * Setup parameters for mass transfer/common envelope event
 *
 *
 * void InitialiseMassTransfer()
 */
void BaseBinaryStar::InitialiseMassTransfer() {

	m_MassTransferTrackerHistory = MT_TRACKING::NO_MASS_TRANSFER;	                                                            // ALEJANDRO - 16/11/2016 - Initiating flag, every timestep, to NO_MASS_TRANSFER. If it undergoes to MT or CEE, it should change.

    m_Star1->InitialiseMassTransfer(m_CEDetails.CEEnow, m_SemiMajorAxisPrime, m_Eccentricity);                                  // initialise mass transfer for star1
    m_Star2->InitialiseMassTransfer(m_CEDetails.CEEnow, m_SemiMajorAxisPrime, m_Eccentricity);                                  // initialise mass transfer for star2

    if (m_Star1->IsRLOF() || m_Star2->IsRLOF()) {                                                                               // either star overflowing its Roche Lobe?
                                                                                                                                // yes - mass transfer if not both CH
        if (OPTIONS->CHE_Option() != CHE_OPTION::NONE && HasTwoOf({STELLAR_TYPE::CHEMICALLY_HOMOGENEOUS})) {                    // CHE enabled and both stars CH?
                                                                                                                                // yes
            // equilibrate masses, circularise, and check for merger

            double mass = (m_Star1->Mass() + m_Star2->Mass()) / 2.0;                                                            // share mass equally
            m_Star1->UpdateAttributes(mass - m_Star1->Mass(), mass - m_Star1->Mass0(), true);                                   // set new mass, mass0 for star 1
            m_Star1->UpdateAttributes(mass - m_Star2->Mass(), mass - m_Star2->Mass0(), true);                                   // set new mass, mass0 for star 2

            m_MassesEquilibrated = true;                                                                                        // record that we've equilbrated

            m_SemiMajorAxis *= (1.0 - (m_Eccentricity * m_Eccentricity));                                                       // circularise; conserve angular momentum
            m_Eccentricity   = 0.0;                                                                                             // now circular

            m_Star1->InitialiseMassTransfer(m_CEDetails.CEEnow, m_SemiMajorAxisPrime, m_Eccentricity);                          // re-initialise mass transfer for star1
            m_Star2->InitialiseMassTransfer(m_CEDetails.CEEnow, m_SemiMajorAxisPrime, m_Eccentricity);                          // re-initialise mass transfer for star2

            m_MassTransfer     = false;                                                                                         // no mass transfer
            m_CEDetails.CEEnow = false;                                                                                         // no common envelope
        }
        else {                                                                                                                  // not both CH, so ...
		    m_MassTransfer = true;                                                                                              // ... mass transfer
            m_CEDetails.CEEnow = false;                                                                                         // no common envelope

		    if (OPTIONS->CirculariseBinaryDuringMassTransfer()) {                                                               // circularise binary to the periapsis separation?
                m_SemiMajorAxisPrime *= OPTIONS->AngularMomentumConservationDuringCircularisation()                             // yes - conserve angular momentum?
                                        ? (1.0 - (m_Eccentricity * m_Eccentricity))                                             // yes - conserve angular momentum
                                        : (1.0 - m_Eccentricity);                                                               // no - angular momentum not coneserved

			    m_Eccentricity        = 0.0;			                                                                        // ALEJANDRO - 22/11/2016 - Think shouldn't use m_Eccentricity but m_EccentrictyPrime. Right now setting both. Check later.     JR: todo: check this
			    m_EccentricityPrime   = 0.0;                                                                                    // JR: todo: check comment above

			    // ALEJANDRO - 23/11/2016 - Bug fix for systems which enter MT being eccentric.
			    // Previous values have to be the ones for periastron as later orbit is modified according to previous values.
			    // If you don't do this, you end up modifying pre-MT pre-circularisation orbit
			    // JR: todo: check that this is proper functionality, or just a kludge - if kludge, resolve it
			    m_SemiMajorAxisPrev   = m_SemiMajorAxisPrime;
			    m_EccentricityPrev    = m_EccentricityPrime;
			    m_OrbitalVelocityPrev = m_OrbitalVelocityPrime;
		    }
        }
    }
    else {
        m_MassTransfer     = false;                                                                                             // no mass transfer
        m_CEDetails.CEEnow = false;                                                                                             // no common envelope
    }

    m_aMassTransferDiff     = 0.0;                                                                                              // iniitially - no changle to orbit (semi-major axis) due to mass transfer
    m_OmegaMassTransferDiff = 0.0;                                                                                              // initially - no change to orbital speed due to mass transfer
}


/*
 * Check if mass transfer should happen (either star, but not both, overflowing Roche Lobe)
 * Perform mass transfer if required
 *
 *
 * void CheckMassTransfer(const double p_Dt)
 *
 * @param   [IN]    p_Dt                        timestep in Myr
 */
void BaseBinaryStar::CheckMassTransfer(const double p_Dt) {

    InitialiseMassTransfer();                                                                                                   // initialise - even if not using mass transfer (sets some flags we might need)

    if (OPTIONS->CHE_Option() != CHE_OPTION::NONE && HasTwoOf({STELLAR_TYPE::CHEMICALLY_HOMOGENEOUS}) && HasStarsTouching()) {  // CHE enabled and both stars CH?
        m_StellarMerger = true;                                                                                                 // just merge
    }
    else {                                                                                                                      // not both CH
        if (OPTIONS->UseMassTransfer()) {                                                                                       // only if using mass transfer (program option)

            if (m_Star1->IsRLOF() && m_Star2->IsRLOF()) {                                                                       // both stars overflowing their Roche Lobe?
			    m_CEDetails.CEEnow = true;                                                                                      // yes - common envelope event - no mass transfer
            }
            else if (m_Star1->IsRLOF() || m_Star2->IsRLOF()) {                                                                  // either star overflowing its Roche Lobe?
                CalculateMassTransfer(p_Dt);                                                                                    // yes - mass transfer
            }
        }
    }
}


/*
 * Calculate the total energy of the binary
 *
 * JR: todo: flesh-out this documentation
 *
 *
 * double CalculateTotalEnergy(const double p_SemiMajorAxis,
 *                             const double p_Star1Mass,
 *                             const double p_Star2Mass,
 *                             const double p_Star1Radius,
 *                             const double p_Star2Radius,
 *                             const double p_Star1_OrbitalFrequency,
 *                             const double p_Star1_OrbitalFrequency,
 *                             const double p_OrbitalVelocity,
 *                             const double p_Star1_GyrationRadius,
 *                             const double p_Star2_GyrationRadius)
 *
 * @param   [IN]    p_SemiMajorAxis             Semi-major axis of the binary
 * @param   [IN]    p_Star1Mass                 Mass of the primary
 * @param   [IN]    p_Star2Mass                 Mass of the secondary
 * @param   [IN]    p_Star1Radius               Radius of the primary
 * @param   [IN]    p_Star2Radius               Radius of the secondary
 * @param   [IN]    p_Star1_OrbitalFrequency    Orbital frequency of the primary
 * @param   [IN]    p_Star1_OrbitalFrequency    Orbital frequency of the secondary
 * @param   [IN]    p_OrbitalVelocity           Orbital velocity of the binary
 * @param   [IN]    p_Star1_GyrationRadius      Gyration radius of the primary
 * @param   [IN]    p_Star2_GyrationRadius      Gyration radius of the secondary
 * @return                                      Total energy of the binary
 */
double BaseBinaryStar::CalculateTotalEnergy(const double p_SemiMajorAxis,
                                            const double p_Star1Mass,
                                            const double p_Star2Mass,
                                            const double p_Star1Radius,
                                            const double p_Star2Radius,
                                            const double p_Star1_OrbitalFrequency,
                                            const double p_Star2_OrbitalFrequency,
                                            const double p_OrbitalVelocity,
                                            const double p_Star1_GyrationRadius,
                                            const double p_Star2_GyrationRadius) {
	double m1  = p_Star1Mass;
	double m2  = p_Star2Mass;

	double R1  = p_Star1Radius;
	double R2  = p_Star2Radius;

	double w1  = p_Star1_OrbitalFrequency;
	double w2  = p_Star2_OrbitalFrequency;

	double w   = p_OrbitalVelocity;

	double ks1 = p_Star1_GyrationRadius;
	double ks2 = p_Star2_GyrationRadius;

    double m1PlusM2 = m1 + m2;
    constexpr double RSOL_TO_AU_2 = RSOL_TO_AU * RSOL_TO_AU;


	double 	Is1  = ks1 * m1 * R1 * R1 * RSOL_TO_AU_2;
	double 	Is2  = ks2 * m2 * R2 * R2 * RSOL_TO_AU_2;
	double	d1   = p_SemiMajorAxis * m2 / m1PlusM2;
	double	d2   = p_SemiMajorAxis * m1 / m1PlusM2;
	double	Itot = (m1 * d1 * d1) + (m2 * d2 * d2);

	return (0.5 * Is1 * w1 * w1) + (0.5 * Is2 * w2 * w2) + (0.5 * Itot * w * w) - (G1 * m1 * m2 / p_SemiMajorAxis);
}


/*
 * Calculate the angular momentum of the binary
 *
 * JR: todo: flesh-out this documentation
 *
 *
 * double CalculateTotalEnergy(const double p_SemiMajorAxis,
 *                             const double p_Eccentricity,
 *                             const double p_Star1Mass,
 *                             const double p_Star2Mass,
 *                             const double p_Star1Radius,
 *                             const double p_Star2Radius,
 *                             const double p_Star1_OrbitalFrequency,
 *                             const double p_Star1_OrbitalFrequency,
 *                             const double p_Star1_GyrationRadius,
 *                             const double p_Star2_GyrationRadius)
 *
 * @param   [IN]    p_SemiMajorAxis             Semi-major axis of the binary
 * @param   [IN]    p_Eccentricity              Eccentricity of the binary
 * @param   [IN]    p_Star1Mass                 Mass of the primary
 * @param   [IN]    p_Star2Mass                 Mass of the secondary
 * @param   [IN]    p_Star1Radius               Radius of the primary
 * @param   [IN]    p_Star2Radius               Radius of the secondary
 * @param   [IN]    p_Star1_OrbitalFrequency    Orbital frequency of the primary
 * @param   [IN]    p_Star1_OrbitalFrequency    Orbital frequency of the secondary
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
                                                const double p_Star1_OrbitalFrequency,
                                                const double p_Star2_OrbitalFrequency,
                                                const double p_Star1_GyrationRadius,
                                                const double p_Star2_GyrationRadius) {
	double m1 = p_Star1Mass;
	double m2 = p_Star2Mass;

	double R1 = p_Star1Radius * RSOL_TO_AU;
	double R2 = p_Star2Radius * RSOL_TO_AU;

	double w1 = p_Star1_OrbitalFrequency;
	double w2 = p_Star2_OrbitalFrequency;

	double ks1 = p_Star1_GyrationRadius;
	double ks2 = p_Star2_GyrationRadius;

	double Is1  = ks1 * m1 * R1 * R1;
	double Is2  = ks2 * m2 * R2 * R2;
    double Jorb = ((m1 * m2) / (m1 + m2)) * sqrt(G1 * (m1 + m2) * p_SemiMajorAxis * (1.0 - (p_Eccentricity * p_Eccentricity)));

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

    // ALEJANDRO - 16/11/2016 - calculate orbital energy and angular momentum
    m_TotalMassPrev                    = m_TotalMassPrime;
    m_ReducedMassPrev                  = m_ReducedMassPrime;
    m_TotalOrbitalEnergyPrev           = m_TotalOrbitalEnergyPrime;
    m_TotalOrbitalAngularMomentumPrev  = m_TotalOrbitalAngularMomentumPrime;

    m_TotalMassPrime                   = m_Star1->Mass() + m_Star2->Mass();
    m_ReducedMassPrime                 = (m_Star1->Mass() * m_Star2->Mass()) / m_TotalMassPrime;
    m_TotalOrbitalEnergyPrime          = CalculateOrbitalEnergy(m_ReducedMassPrime, m_TotalMassPrime, m_SemiMajorAxisPrime);
    m_TotalOrbitalAngularMomentumPrime = CalculateOrbitalAngularMomentum(m_ReducedMassPrime, m_TotalMassPrime, m_SemiMajorAxisPrime);

    // ALEJANDRO - 16/11/2016 - calculate energy and angular momentum using regular conservation of energy, specially useful for checking tides and rotational effects
    m_TotalEnergyPrime                 = CalculateTotalEnergyPrime();
    m_TotalAngularMomentumPrime        = CalculateAngularMomentumPrime();
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

    // update mass of star1 according to mass loss and mass transfer, then update age accordingly
    (void)m_Star1->UpdateAttributes(m_Star1->MassPrev() - m_Star1->Mass() + m_Star1->MassLossDiff() + m_Star1->MassTransferDiff(), 0.0);        // update mass for star1
    m_Star1->UpdateInitialMass();                                                                       // update initial mass of star1 (MS, HG & HeMS)  JR: todo: fix this kludge one day - mass0 is overloaded, and isn't always "initial mass"
    m_Star1->UpdateAgeAfterMassLoss();                                                                  // update age of star1
    m_Star1->ApplyMassTransferRejuvenationFactor();                                                     // apply age rejuvenation factor for star1

    // rinse and repeat for star2
    (void)m_Star2->UpdateAttributes(m_Star2->MassPrev() - m_Star2->Mass() +m_Star2->MassLossDiff() + m_Star2->MassTransferDiff(), 0.0);        // update mass for star2
    m_Star2->UpdateInitialMass();                                                                       // update initial mass of star 2 (MS, HG & HeMS)  JR: todo: fix this kludge one day - mass0 is overloaded, and isn't always "initial mass"
    m_Star2->UpdateAgeAfterMassLoss();                                                                  // update age of star2
    m_Star2->ApplyMassTransferRejuvenationFactor();                                                     // apply age rejuvenation factor for star2

    // update binary
    m_OrbitalVelocityPrime = m_OrbitalVelocityPrev + m_OmegaMassLossDiff + m_OmegaMassTransferDiff;     // should here be a diff quantity because of MB?    JR: todo: ?
    m_SemiMajorAxisPrime   = m_SemiMajorAxisPrev + m_aMassLossDiff + m_aMassTransferDiff;

    // if CHE enabled, update rotational frequency for constituent stars - assume tidally locked
    if (OPTIONS->CHE_Option() != CHE_OPTION::NONE) m_Star1->SetOmega(m_OrbitalVelocityPrime);
    if (OPTIONS->CHE_Option() != CHE_OPTION::NONE) m_Star2->SetOmega(m_OrbitalVelocityPrime);

    CalculateEnergyAndAngularMomentum();                                                                // perform energy and angular momentum calculations
}


/*
 * Perform calculations required before evaluating the binary
 *
 * Calculates:
 *
 *    Lambdas (if necessary)
 *    Zetas (if necessary)
 *    Total angular momentum (previous) - m_TotalAngularMomentumPrev
 *
 * void EvaluateBinaryPreamble()
 */
void BaseBinaryStar::EvaluateBinaryPreamble() {
	// ALEJANDRO - 14/11/2016 - After some profiling done by Jim,
	// seems like calculate the Loveridge lambda at each timestep takes a lot of time.
	// Therefore, I moved it to this function and only calculated it if explicitely indicated.
	//
	// JR: I did some work on calculating Loveridge lambda - it will be less computationally
	// expensive, but still better to do it fewer times if possible

	if (OPTIONS->LambdaCalculationEveryTimeStep()) {
        m_Star1->CalculateLambdas();
        m_Star2->CalculateLambdas();
    }

	if (OPTIONS->ZetaCalculationEveryTimeStep()) {                      // ALEJANDRO - 16/10/2017 - Calculate zetas if specified        JR: todo: check - also called elsewhere
        m_Star1->CalculateZetas();
        m_Star2->CalculateZetas();
    }

    m_TotalAngularMomentumPrev = CalculateAngularMomentumPrev();        // squirrel away previous value for total angular momentum
}


/*
 * Evaluate the binary system
 *
 *    - caclulate any mass transfer
 *    - calculate mass loss due to wonds
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

    EvaluateBinaryPreamble();                                                                                               // get things ready - do some house-keeping

    CheckMassTransfer(p_Dt);                                                                                                // calculate mass transfer if necessary

    if (!(OPTIONS->CHE_Option() != CHE_OPTION::NONE && HasTwoOf({STELLAR_TYPE::CHEMICALLY_HOMOGENEOUS}))) {                 // CHE enabled and both stars CH
                                                                                                                            // no - proceed
        CalculateWindsMassLoss();                                                                                           // calculate mass loss dues to winds

        if ( (m_CEDetails.CEEnow || m_StellarMerger) &&                                                                     // CEE or merger?
            !(OPTIONS->CHE_Option() != CHE_OPTION::NONE && HasTwoOf({STELLAR_TYPE::CHEMICALLY_HOMOGENEOUS}))) {             // yes - avoid CEE if CH+CH
                ResolveCommonEnvelopeEvent();                                                                               // resolve CEE - immediate event
        }
        else if (m_Star1->IsSNevent() || m_Star2->IsSNevent()) {
            EvaluateSupernovae();                                                                                           // evaluate supernovae (both stars) - immediate event
        }
        else {
            ResolveMassChanges();                                                                                           // apply mass loss and mass transfer as necessary
        }

        (void)m_Star1->UpdateAttributes(0.0, 0.0, true);                                                                    // recalculate stellar attributes for star1
        (void)m_Star2->UpdateAttributes(0.0, 0.0, true);                                                                    // recalculate stellar attributes for star2

        EvaluateSupernovae();                                                                                               // evaluate supernovae (both stars)   JR: todo: ?
        ResolveTides();                                                                                                     // resolve tides
        CalculateEnergyAndAngularMomentum();                                                                                // perform energy and angular momentum calculations

        m_Star1->UpdateMagneticFieldAndSpin(m_CEDetails.CEEnow, m_Dt * MYR_TO_YEAR * SECONDS_IN_YEAR, EPSILON_PULSAR);      // update pulsar parameters for star1
        m_Star2->UpdateMagneticFieldAndSpin(m_CEDetails.CEEnow, m_Dt * MYR_TO_YEAR * SECONDS_IN_YEAR, EPSILON_PULSAR);      // update pulsar parameters for star2
    }
}


/*
 * Calculate timestep (dt)
 *
 * Calculates the new timestep
 *
 * JR: todo: flesh-out this documentation
 *
 *
 * double CalculateDt(const double p_Dt, const Star* const p_Primary, const Star* const p_Secondary)
 *
 * @param   [in]        p_Dt                    Existing timestep value (in Myr)
 * @param   [in]        p_Primary               Primary star (for this calculation)
 * @param   [in]        p_Secondary             Secondary star (for this calculation)
 * @return                                      New timestep value (in Myr)
 */
double BaseBinaryStar::CalculateDt(const double p_Dt, const Star* const p_Primary, const Star* const p_Secondary) {

    double newDt = p_Dt;                                                                                                    // default is unchanged

    double tmp     = RSOL_TO_AU / (m_SemiMajorAxisPrime * (1.0 - m_Eccentricity));
    double RLRatio = tmp * p_Primary->Radius() / CalculateRocheLobeRadius_Static(p_Primary->Mass(), p_Secondary->Mass());

    // calculate new dt
    // JR: todo: this code is broken in the original code
    //
    // JR: todo: NOTE: m_FastPhaseCaseA (binary variable, not stellar) is NEVER set TRUE.  cf STELLAR value
    if (utils::Compare(RLRatio, MASS_TRANSFER_THRESHOLD) >= 0) {                                                            // should be done in this or in previous timestep values?
        if (p_Primary->DetermineEnvelopeType() == ENVELOPE::RADIATIVE) {
            newDt = m_FastPhaseCaseA ? p_Primary->CalculateThermalTimescale() : p_Primary->CalculateDynamicalTimescale();   // should only do this for fastPhase aRLOF, just one timestep.
        }
        else {
            newDt = p_Primary->CalculateThermalTimescale();
        }
    }
    else if (utils::Compare(RLRatio, ROCHE_LOBE_UPPER_THRESHOLD) >= 0) {
        newDt /= 10.0;
    }
    else if (utils::Compare(RLRatio, ROCHE_LOBE_LOWER_THRESHOLD) >= 0) {
        newDt /= 5.0;
    }

    return newDt;
}


/*
 * Calculate the timestepbased on:
 *
 *    stellar nuclear timescales,
 *    angular momentum loss timescales, or
 *    the specific case of the shortest timescales with RLOF
 *
 * Based onHurley et al., 2002, section 2.8
 *
 *
 * double ChooseTimestep(const double p_Dt)
 *
 * @param   [IN]    p_Dt                        The suggested timestep to evolve
 * @return                                      Calculated/recommended timestep
 */
double BaseBinaryStar::ChooseTimestep(const double p_Dt) {

    double timestep = 0.0;      // ALEJANDRO - 02/12/2015 - MASSLESS_REMNANT timescale = 0.0

	if (!m_Star1->IsOneOf({ STELLAR_TYPE::MASSLESS_REMNANT }) && !m_Star2->IsOneOf({ STELLAR_TYPE::MASSLESS_REMNANT })) {
        double dt1 = CalculateDt(p_Dt, m_Star1, m_Star2);
        double dt2 = CalculateDt(p_Dt, m_Star2, m_Star1);

        timestep = std::min(dt1, dt2);
	}

    return std::max(timestep, NUCLEAR_MINIMUM_TIMESTEP);
}


/*
 * Calculate next timestep for binary evolution
 *
 * Timesteps are calculated for each individual star, based on stellar type, age, etc.
 * The minimum of the calculated timesteps is returned as the timestep.
 *
 * Rather than evolve and revert here we just create copies of the constituent stars,
 * evolve them for one times step to get the new timestep, then discared them.  Maybe
 * one day we can figure out how to do this without the overhead of evolving and
 * discarding,
 *
 *
 * double CalculateTimestep(const double p_Dt)
 *
 * @param   [IN]    p_Dt                        The suggested timestep to evolve
 * @return                                      Calculated/recommended timestep
 */
double BaseBinaryStar::CalculateTimestep(const double p_Dt) {

    BinaryConstituentStar* star1Copy = new BinaryConstituentStar(*m_Star1);             // copy star1
    double dt1 = star1Copy->EvolveOneTimestep(p_Dt);                                    // evolve the copy one timestep and get suggested timestep

    delete star1Copy; star1Copy = nullptr;                                              // nuke the copy

    // rinse and repeat for star2
    BinaryConstituentStar* star2Copy = new BinaryConstituentStar(*m_Star2);             // copy star2
    double dt2 = star2Copy->EvolveOneTimestep(p_Dt);                                    // evolve the copy one timestep and get suggested timestep

    delete star2Copy; star2Copy = nullptr;                                              // nuke the copy

    return std::min(dt1, dt2);
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

    m_TimePrev = m_Time;    // Remember current simulation time
    m_Time    += p_Dt;      // Advance physical simulation time
    m_Dt       = p_Dt;      // Set timestep
}


/*
 * Evolve the binary a single timestep - timestep is provided    JR: todo: fix this documetation - this is for SSE version
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
 * EVOLUTION_STATUS Evolve(const int p_Index)
 *
 * @param   [IN]    p_Index                     The ordinal number of the binary
 * @return                                      Status of the evolution (EVOLUTION_STATUS)
 */
EVOLUTION_STATUS BaseBinaryStar::Evolve(const int p_Index) {

    EVOLUTION_STATUS evolutionStatus = EVOLUTION_STATUS::CONTINUE;

    if (OPTIONS->OnlyDoubleCompactObjects() && m_Star2->Mass() < MINIMUM_MASS_SECONDARY) {                                                  // check size of secondary
        m_SecondaryTooSmallForDCO = true;
        evolutionStatus           = EVOLUTION_STATUS::SECONDARY_TOO_SMALL_FOR_DCO;                                                          // too small - don't bother - no possibility of forming a double compact object
    }

    if (HasStarsTouching()) {                                                                                                               // check if stars are touching
        m_StellarMerger        = true;
        m_StellarMergerAtBirth = true;
        evolutionStatus        = EVOLUTION_STATUS::STELLAR_MERGER_AT_BIRTH;                                                                 // binary components are touching - merger at birth
    }

    PrintDetailedOutput(p_Index);                                                                                                           // print (log) detailed output for binary

    if (OPTIONS->PopulationDataPrinting()) {                                                                                                // JR: todo: what is the aim of PopulationDataPrinting?
        SAY("\nGenerating a new binary - " << p_Index);
        SAY("Binary has masses " << m_Star1->Mass() << " & " << m_Star2->Mass());
        SAY("Binary has initial separation " << m_SemiMajorAxisPrime);
        SAY("RandomSeed " << m_RandomSeed);
    }

    if (evolutionStatus == EVOLUTION_STATUS::CONTINUE) {                                                                                    // continue evolution

        // evolve the current binary up to the maximum evolution time (and number of steps)
        double dt      = std::min(m_Star1->CalculateTimestep(), m_Star2->CalculateTimestep()) / 1000.0;                                     // initialise the timestep
        int    stepNum = 1;                                                                                                                 // initialise step number
        while (evolutionStatus == EVOLUTION_STATUS::CONTINUE) {                                                                             // perform binary evolution - iterate over timesteps until told to stop

            m_TotalAngularMomentumPrev = m_TotalAngularMomentumPrime;   // Is this line ok here?        JR: todo - this probably should be in evaluateBinary(), except that evaluateBinary() may not be executed at each timestep - maybe this has to stay here

            EvolveOneTimestep(dt);                                                                                                          // evolve the binary system one timestep

            StashPreviousRLOFProperties();                                                                                                  // stash RLOF properties

            // check for problems
                 if (m_Error != ERROR::NONE)                                        evolutionStatus = EVOLUTION_STATUS::SSE_ERROR;          // SSE error for either constituent star
            else if (HasOneOf({ STELLAR_TYPE::MASSLESS_REMNANT }))                  evolutionStatus = EVOLUTION_STATUS::MASSLESS_REMNANT;   // at least one massless remnant has been formed
            else if (HasStarsTouching())                                            evolutionStatus = EVOLUTION_STATUS::STARS_TOUCHING;     // binary components are touching (should usually be avoided as MT or CE should happen prior to this)
            else if (IsUnbound() && !OPTIONS->EvolveUnboundSystems())               evolutionStatus = EVOLUTION_STATUS::UNBOUND;            // binary is unbound and we don't want unbound systems
            else if (!IsGravitationallyBound() && !OPTIONS->EvolveUnboundSystems()) evolutionStatus = EVOLUTION_STATUS::UNBOUND;            // binary is not gravitationally bound and we don't want unbound systems
            else {                                                                                                                          // continue evolution

                EvaluateBinary(dt);                                                                                                         // evaluate the binary at this timestep

                StashBeBinaryProperties();                                                                                                  // stash BeBinary properties
                PrintBeBinary();                                                                                                            // print (log) BeBinary properties

                StashCurrentRLOFProperties();                                                                                               // stash RLOF properties
                PrintRLOFParameters();                                                                                                      // print (log) RLOF parameters

                // check for problems
                     if (HasStarsTouching())                              evolutionStatus = EVOLUTION_STATUS::STARS_TOUCHING;               // binary components are touching (should usually be avoided as MT or CE should happen prior to this)
                else if (!OPTIONS->EvolveUnboundSystems() && IsUnbound()) evolutionStatus = EVOLUTION_STATUS::UNBOUND;                      // we're not evolving unbound systems and binary is unbound

                if (evolutionStatus == EVOLUTION_STATUS::CONTINUE) {                                                                        // continue evolution?

                    PrintDetailedOutput(p_Index);                                                                                           // print (log) detailed output for binary

                    // check for problems
                         if (m_Error != ERROR::NONE) evolutionStatus = EVOLUTION_STATUS::BINARY_ERROR;                                      // error in binary evolution
                    else if (m_StellarMerger)        evolutionStatus = EVOLUTION_STATUS::STELLAR_MERGER;                                    // constituent stars have merged

                    if (evolutionStatus == EVOLUTION_STATUS::CONTINUE) {                                                                    // continue evolution?

                        if (HasOneOf({ STELLAR_TYPE::NEUTRON_STAR })) PrintPulsarEvolutionParameters();                                     // print (log) pulsar evolution parameters    JR: todo: WD?

                        if (IsDCO()) {                                                                                                      // double compact object?
                            ResolveCoalescence();                                                                                           // yes - resolve coalescence

                            if (OPTIONS->AIS_ExploratoryPhase()) (void)m_AIS.CalculateDCOHit(this);                                         // track if we have an AIS DCO hit - internal counter is updated (don't need return value here)

                            if (!OPTIONS->Quiet()) SAY(ERR_MSG(ERROR::BINARY_EVOLUTION_STOPPED) << ": Double compact object");              // announce that we're stopping evolution
                            evolutionStatus = EVOLUTION_STATUS::STOPPED;                                                                    // stop evolving
                        }

                        // check for problems
                        if (evolutionStatus == EVOLUTION_STATUS::CONTINUE) {                                                                // continue evolution?
                                 if (m_Error != ERROR::NONE)               evolutionStatus = EVOLUTION_STATUS::BINARY_ERROR;                // error in binary evolution
                            else if (IsWDandWD())                          evolutionStatus = EVOLUTION_STATUS::WD_WD;                       // ALEJANDRO - 16/03/2017 - Check for double WD systems so we stop evolving them as they take time to process and we are currently not interested in dealing with them.
                            else if (m_Time > OPTIONS->MaxEvolutionTime()) evolutionStatus = EVOLUTION_STATUS::TIMES_UP;                    // evolution time exceeds maximum
                        }
                    }
                }
            }

            if (stepNum >= OPTIONS->MaxNumberOfTimestepIterations()) evolutionStatus = EVOLUTION_STATUS::STEPS_UP;                          // number of timesteps for evolution exceeds maximum

            if (evolutionStatus == EVOLUTION_STATUS::CONTINUE) {                                                                            // continue evolution?

                dt = std::min(m_Star1->CalculateTimestep(), m_Star2->CalculateTimestep());                                                  // new timestep
                                                                                                                                            // yes - prepare for next timestep
                stepNum++;                                                                                                                  // increment stepNum
                dt = CalculateTimestep(ChooseTimestep(dt));                                                                                 // calculate next timestep
            }


        }

        if (evolutionStatus == EVOLUTION_STATUS::STEPS_UP) {                                                                                // stopped because max timesteps reached?
            SHOW_ERROR(ERROR::BINARY_EVOLUTION_STOPPED);                                                                                    // show error
        }
    }

    PrintBinarySystemParameters();                                                                                                          // print (log) binary system parameters

    return evolutionStatus;
}
