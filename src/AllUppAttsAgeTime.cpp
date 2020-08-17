/************************/
/*   BaseStar.h         */
/************************/
STELLAR_TYPE UpdateAttributesAndAgeOneTimestep(
        const double p_DeltaMass,
        const double p_DeltaMass0,
        const double p_DeltaTime,
        const bool   p_ForceRecalculate);
void   UpdateAttributesAndAgeOneTimestepPreamble(
        const double p_DeltaMass, 
        const double p_DeltaMass0, 
        const double p_DeltaTime);
void   AgeOneTimestepPreamble(const double p_DeltaTime);
virtual void UpdateAgeAfterMassLoss() { }                                                                                                                            // Default is NO-OP

double m_Age;    // Current effective age (changes with mass loss/gain)(myrs)
double Age() const  { return m_Age; } // RTW 10/7/20 - should rename to GetAge - it is somewhat ambiguous if this should actually age the star
void   ApplyMassTransferRejuvenationFactor()  { m_Age *= CalculateMassTransferRejuvenationFactor(); }             // Apply age rejuvenation factor
double CalculateMassLossValues(const bool p_UpdateMDot = false, const bool p_UpdateMDt = false);                                                               // JR: todo: better name?
virtual void UpdateInitialMass() { }                                                                                                                                 // Default is NO-OP

/************************/
/*   Star.h             */
/************************/
void UpdateAttributes() { (void)UpdateAttributes(0.0, 0.0, true); }
STELLAR_TYPE UpdateAttributes(
        const double p_DeltaMass,
        const double p_DeltaMass0,
        const bool   p_ForceRecalculate = false);
STELLAR_TYPE UpdateAttributesAndAgeOneTimestep(
        const double p_DeltaMass,
        const double p_DeltaMass0,
        const double p_DeltaTime,
        const bool   p_Switch = true,
        const bool   p_ForceRecalculate = false);
STELLAR_TYPE AgeOneTimestep(const double p_Dt, bool p_Switch = true);

double CalculateMassLossValues(
        const bool p_UpdateMDot = false, 
        const bool p_UpdateMDt = false)  
        { return m_Star->CalculateMassLossValues(p_UpdateMDot, p_UpdateMDt); }

/************************/
/*   BaseStar.cpp       */
/************************/


/*
 * Calculate the radius-time exponent zeta
 * Derived by the nuclear evolution of the star according to SSE, calculated using a fake time step and recomputing the stars radius.
 *
 *
 * double CalculateZetaNuclear(const double p_DeltaTime)
 *
 * @param   [IN]    p_DeltaTime                 Timestep (Myr)
 * @return                                      Radius-Time exponent Zeta for the nuclear timescale
 */
double BaseStar::CalculateZetaNuclear(const double p_DeltaTime) {

    BaseStar* starCopy = new BaseStar(*this);	                                                                                // copy of star - about to be updated for fake mass loss

    double radiusBeforeTimeStep = starCopy->Radius();                                                                           // radius before timestep
    double ageBeforeTimeStep    = starCopy->Age();

    SHOW_ERROR_IF(utils::Compare(radiusBeforeTimeStep, 0.0) <= 0, ERROR::RADIUS_NOT_POSITIVE_ONCE, "Before fake timestep");     // show error if radius <= 0
    SHOW_ERROR_IF(utils::Compare(ageBeforeTimeStep,    0.0) <  0, ERROR::AGE_NEGATIVE_ONCE,        "Before fake timestep");     // show error if age < 0

    double logRadiusBeforeTimeStep = utils::Compare(radiusBeforeTimeStep, 0.0) > 0 ? log(radiusBeforeTimeStep) : 0.0;

    starCopy->UpdateAttributesAndAgeOneTimestep(0.0, 0.0, 0.0, true);                                                           // allow star to respond to previous mass loss changes      JR: todo: is this really necessary?

    // With the current way of doing things, this star has already lost a bunch of mass due to the mass loss caller,
    // need to update radius of star to respond to this, before we evolve it for a short amount of time assuming no mass loss.
    // When Alejandro moves the massLossCaller function, we will need to do things differently, and this will need fixing.
    // JR: todo: what does this actually mean?

    ageBeforeTimeStep = starCopy->Age();                                                                                        // age before timestep      JR: todo: why don't we get the (possibly) new radius here too?
    starCopy->UpdateAttributesAndAgeOneTimestep(0.0, 0.0, p_DeltaTime, false);                                                  // age the star 'p_DeltaTime' Myr

    double radiusAfterTimeStep = starCopy->Radius();
    double ageAfterTimeStep    = starCopy->Age();

    delete starCopy; starCopy = nullptr;

    SHOW_ERROR_IF(utils::Compare(radiusAfterTimeStep, 0.0) <= 0, ERROR::RADIUS_NOT_POSITIVE_ONCE, "After fake timestep");       // show error if radius <= 0
    SHOW_ERROR_IF(utils::Compare(ageAfterTimeStep,    0.0) <  0, ERROR::AGE_NEGATIVE_ONCE,        "After fake timestep");       // show error if age < 0

    double logRadiusAfterTimeStep = utils::Compare(radiusAfterTimeStep, 0.0) > 0 ? log(radiusAfterTimeStep) : 0.0;

    return (logRadiusAfterTimeStep - logRadiusBeforeTimeStep) / (ageAfterTimeStep - ageBeforeTimeStep);

}



/*
 * Calculate the radius-mass exponent Zeta, assuming the star has had time to recover its thermal equilibrium.
 * Zeta is calculated using a fake mass change step and recomputing the star's attributes.
 *
 *
 * double CalculateZetaThermal(double p_PercentageMassChange)
 *
 * @param   [IN]    p_PercentageMassChange      Percentage of mass the star should artificially lose/gain in order to calcluate Zeta (thermal)
 *                                              (sign of p_PercentageMassChange determine if loss or gain - -ve is loss, +ve is gain)
 * @return                                      Radius-Mass exponent Zeta for the thermal timescale
 */
double BaseStar::CalculateZetaThermal(double p_PercentageMassChange) {

    BaseStar* starCopy = new BaseStar(*this);	                                                                                // copy of star - about to be updated for fake mass loss

    starCopy->UpdateAttributesAndAgeOneTimestep(0.0, 0.0, 0.0, true);                                                           // allow star to respond to previous mass loss changes      JR: todo: is this really necessary?
    starCopy->UpdateAttributesAndAgeOneTimestep(-(starCopy->Mass() * FAKE_MASS_LOSS_PERCENTAGE / 100.0), 0.0, 0.0, false);      // apply fake mass loss                                     JR: todo: why do we do this...?

    // record properties of the star before fake mass change
    double radiusBeforeMassLoss = starCopy->Radius();                                                                           // radius before fake mass change       JR: todo: didn't we just do fake mass loss...?
    double massBeforeMassLoss   = starCopy->MassPrev();                                                                         // mass before fake mass change - due to order of updating radius bug     JR: todo: check this

    SHOW_ERROR_IF(utils::Compare(radiusBeforeMassLoss, 0.0) <= 0, ERROR::RADIUS_NOT_POSITIVE_ONCE, "Before fake mass change");  // show error if radius <= 0
    SHOW_ERROR_IF(utils::Compare(massBeforeMassLoss,   0.0) <= 0, ERROR::MASS_NOT_POSITIVE_ONCE,   "Before fake mass change");  // show error if mass <= 0

    double logRadiusBeforeMassLoss = utils::Compare(radiusBeforeMassLoss, 0.0) > 0 ? log(radiusBeforeMassLoss) : 0.0;
    double logMassBeforeMassLoss   = utils::Compare(massBeforeMassLoss,   0.0) > 0 ? log(massBeforeMassLoss  ) : 0.0;

    double deltaMass           = starCopy->Mass() * p_PercentageMassChange / 100.0;                                             // fake mass change - sign of p_PercentageMassChange determines whether loss or gain
    double massAfterMassLoss   = starCopy->Mass() + deltaMass;                                                                  // mass after (just) fake mass change
    starCopy->UpdateAttributesAndAgeOneTimestep(starCopy->Mass() * p_PercentageMassChange / 100.0, 0.0, 0.0, false);            // apply fake mass change and recalculate attributes of star
    double radiusAfterMassLoss = starCopy->m_Radius;                                                                            // radius after fake mass change
    delete starCopy; starCopy  = nullptr;

    SHOW_ERROR_IF(utils::Compare(radiusAfterMassLoss, 0.0) <= 0, ERROR::RADIUS_NOT_POSITIVE_ONCE, "After fake mass change");    // show error if radius <= 0
    SHOW_ERROR_IF(utils::Compare(massAfterMassLoss,   0.0) <= 0, ERROR::MASS_NOT_POSITIVE_ONCE,   "After fake mass change");    // show error if mass <= 0

    double logRadiusAfterMassLoss = utils::Compare(radiusAfterMassLoss, 0.0) > 0 ? log(radiusAfterMassLoss) : 0.0;
    double logMassAfterMassLoss   = utils::Compare(massAfterMassLoss,   0.0) > 0 ? log(massAfterMassLoss  ) : 0.0;

    return (logRadiusAfterMassLoss - logRadiusBeforeMassLoss) / (logMassAfterMassLoss - logMassBeforeMassLoss);                 // zeta thermal
}



/*
 * Calculate values for dt, mDot and mass assuming mass loss is applied
 *
 * Class member variables m_Mdot and m_Dt are updated directly by this function if required (see paramaters)
 * Class member variables m_Mass is not updated directly by this function - the calculated mass is returned as the functional return
 *
 * - calculates (and limits) mass loss
 * - calculate new timestep (dt) and mass loss rate (mDot) to match (possibly limited) mass loss
 * - calculates new mass (mass) based on (possibly limited) mass loss
 *
 * Returns existing value for mass if mass loss not being used (program option)
 *
 *
 * double CalculateMassLossValues()                                                             // JR: todo: pick a better name for this...
 * @return                                      calculated mass (mSol)
 */
double BaseStar::CalculateMassLossValues(const bool p_UpdateMDot, const bool p_UpdateMDt) {

    double dt   = m_Dt;
    double mDot = m_Mdot;
    double mass = m_Mass;

    if (OPTIONS->UseMassLoss()) {                                           // only if using mass loss (program option)

        mDot = CalculateMassLossRate();                                     // calculate mass loss rate
        double massLoss = CalculateMassLoss_Static(mass, mDot, dt);         // calculate mass loss - limited to (mass * MAXIMUM_MASS_LOSS_FRACTION)

        // could do this without the test - we know the mass loss may already
        // have been limited.  This way is probably marginally faster
        if (utils::Compare(massLoss, (mass * MAXIMUM_MASS_LOSS_FRACTION)) < 0) {
            mass -= massLoss;                                               // new mass based on mass loss
        }
        else {
            dt    = massLoss / (mDot * 1.0E6);                              // new timestep to match limited mass loss
            mDot  = massLoss / (dt * 1.0E6);                                // new mass loss rate to match limited mass loss
            mass -= massLoss;                                               // new mass based on limited mass loss

            if (p_UpdateMDt) m_Dt = dt;                                     // update class member variable if necessary
        }

        if (p_UpdateMDot) m_Mdot = mDot;                                    // update class member variable if necessary
    }

    return mass;
}


/*
 * Resolve mass loss
 *
 * - calculates mass loss rate
 * - calculates (and limits) mass loss
 * - resets timestep (m_Dt) and mass loss rate (m_Mdot) to match (possibly limited) mass loss
 * - calculates and sets new mass (m_Mass) based on (possibly limited) mass loss
 * - applies mass rejuvenation factor and calculates new age
 *
 *
 * double ResolveMassLoss()
 */
void BaseStar::ResolveMassLoss() {

    if (OPTIONS->UseMassLoss()) {
        m_Mass = CalculateMassLossValues(true, true);                           // calculate new values assuming mass loss applied

        UpdateInitialMass();                                                    // update initial mass (MS, HG & HeMS)  JR: todo: fix this kludge one day - mass0 is overloaded, and isn't always "initial mass"
        UpdateAgeAfterMassLoss();                                               // update age (MS, HG & HeMS)
        ApplyMassTransferRejuvenationFactor();                                  // apply age rejuvenation factor
    }
}



/*
 * Calculate next timestep for stellar evolution
 *
 * Timestep based on stellar type, age, etc.
 *
 *
 * double CalculateTimestep()
 *
 * @return                                      Timestep
 */
double BaseStar::CalculateTimestep() {

    // the GBParams and Timescale calculations need to be done
    // before the timestep calculation - since the binary code
    // calls this functiom, the GBParams and Timescale functions
    // are called here

    // RTW 05/07/20 - We should have a uniquely single star evolution code which does this separately

    CalculateGBParams();                                                                // calculate giant branch parameters
    CalculateTimescales();                                                              // calculate timescales

    double dt = ChooseTimestep(m_Age);
    if (utils::Compare(TIMESTEP_REDUCTION_FACTOR, 1.0) != 0) {                          // timestep reduction factor == 1.0?
        if (!IsOneOf({ STELLAR_TYPE::MS_LTE_07, STELLAR_TYPE::MS_GT_07 })) {            // no - check stellar type
            // RTW 01/07/20 - Why do MS stars not get reduced?
            dt /= TIMESTEP_REDUCTION_FACTOR;                                            // apply timestep reduction factor
        }
    }

    return LimitTimestep(dt);
}


/*
 * Set parameters required before updating stellar attributes (via evolution) - modify star attributes
 *
 * Will apply mass changes to m_Mass and/or m_Mass0.  Note discussion in documentation for
 * UpdateAttributes() in Star.cpp - changing m_Mass0 is a bit of a kludge and should be fixed.
 *
 * If the timestep (p_DeltaTime) is > 0 then m_Mass and m_Radius will be saved as m_MassPrev and m_RadiusPrev respectively.
 *
 *
 * void UpdateAttributesAndAgeOneTimestepPreamble(const double p_DeltaMass, const double p_DeltaMass0, const double p_DeltaTime)
 *
 * @param   [IN]    p_DeltaMass                 The change in mass to apply in Msol
 * @param   [IN]    p_DeltaMass0                The change in mass0 to apply in Msol
 * @param   [IN]    p_DeltaTime                 Timestep in Myr
 */
void BaseStar::UpdateAttributesAndAgeOneTimestepPreamble(const double p_DeltaMass, const double p_DeltaMass0, const double p_DeltaTime) {

    if (utils::Compare(p_DeltaMass,  0.0) != 0) { m_Mass  = max(0.0, m_Mass  + p_DeltaMass);  }     // update mass as required (only change if delta != 0 with tolerance) and prevent -ve
    if (utils::Compare(p_DeltaMass0, 0.0) != 0) { m_Mass0 = max(0.0, m_Mass0 + p_DeltaMass0); }     // update mass0 as required (only change if delta != 0 with tolerance) and prevent -ve

    // record some current values before they are (possibly) changed by evolution
    if (p_DeltaTime > 0.0) {                                                                        // don't use utils::Compare() here
        m_StellarTypePrev = m_StellarType;
        m_DtPrev          = m_Dt;
        m_MassPrev        = m_Mass;
        m_RadiusPrev      = m_Radius;
    }

    // the GBParams and Timescale calculations need to be done
    // before taking the timestep - since the binary code ultimately
    // calls this UpdateAttributesAndAgeOneTimestep, the GBParams and
    // Timescale functions are called here.
    //
    // JR: todo: we should revisit where and how often we recalulate
    // GBParams and Timescales.  The problem is that there are multiple
    // entry points into the calculate/take timestep code that it isn't
    // always obvious where we need to do this...  A project for another
    // time.

    CalculateGBParams();                                                                            // calculate giant branch parameters
    CalculateTimescales();                                                                          // calculate timescales
}


/*
 * Apply mass changes if required, age the star one timestep, advance the simulation time, and update the
 * attributes of the star.
 *
 * The star's attributes (Age, Radius, Luminosity etc.) are calculated and updated as required.
 *
 * Free parameters in the update process are the star's mass (m_Mass), initial mass (m_Mass0), the star's age
 * (m_Age) and the simulation time attribute (m_Time):
 *
 *    - if required, the star's mass is changed by the amount passed as the p_DeltaMass parameter before other
 *      attributes are updated.  The p_DeltaMass parameter may be zero, in which case no change is made to the
 *      star's mass before the attributes of the star are calculated.
 *
 *    - if required, the star's initial mass is changed by the amount passed as the p_DeltaMass0 parameter before
 *      other attributes are updated.  The p_DeltaMass parameter may be zero, in which case no change is made to
 *      the star's mass before the attributes of the star are calculated.  This should be used infrequently, and
 *      is really a kludge because the Mass0 attribute in Hurley et al. 2000 was overloaded after the introduction
 *      of mass loss (see section 7.1).  We should really separate the different uses of Mass0 in the code and
 *      use a different variable - initial mass shouldn't change (other than to initial mass upon entering a
 *      stellar phase - it doesn't make a lot of sense for initial mass to change during evolution through the
 *      phase).         JR: todo: action this paragraph.
 *
 *    - if required, the star is aged by the amount passed as the p_DeltaTime parameter, and the simulation time is
 *      advanced by the same amount, before other attributes are updated.  The p_deltaTime parameter may be zero,
 *      in which case no change is made to the star's age or the physical time attribute.
 *
 *
 * Checks whether the star:
 *    - is a massless remnant (checked after applying p_DeltaMass and p_DeltaMass0, but before applying p_DeltaTime)
 *    - has become a supernova (checked after applying p_DeltaMass and p_DeltaMass0, but before applying p_DeltaTime)
 *    - should skip this phase for this timestep (checked after applying p_DeltaMass, p_DeltaMass0 and p_DeltaTime)
 *
 * If none of the above are true the star evolves on phase for the specified timestep (which may be 0, in which case
 * the star's attributes other than age are re-calculated), then the need to evolve the star off phase is checked.
 *
 * If p_DeltaMass, p_DeltaMass0 and p_DeltaTime are all passed as zero the checks for massless remnant and supernova
 * are performed (and consequential changes made), but no other changes to the star's attributes are made - unless
 * the p_ForceRecalculate parameter is set true.        JR: todo: I'm not convinced p_ForceRecalculate is necessary - check it
 *
 * The functional return is the stellar type to which the star should evolve.  The returned stellar type is just the
 * stellar type of the star upon entry if it should remain on phase.  The star's stellar type is not changed here.
 *
 *
 * STELLAR_TYPE UpdateAttributesAndAgeOneTimestep(const double p_DeltaMass,
 *                                                const double p_DeltaMass0,
 *                                                const double p_DeltaTime,
 *                                                const bool   p_ForceRecalculate)
 *
 * @param   [IN]    p_DeltaMass                 The change in mass to apply in Msol
 * @param   [IN]    p_DeltaMass0                The change in mass0 to apply in Msol
 * @param   [IN]    p_DeltaTime                 The timestep to evolve in Myr
 * @param   [IN]    p_ForceRecalculate          Specifies whether the star's attributes should be recalculated even if the three deltas are 0.0
 *                                              (optional, default = false)
 * @return                                      Stellar type to which star should evolve
 */
STELLAR_TYPE BaseStar::UpdateAttributesAndAgeOneTimestep(const double p_DeltaMass,
                                                         const double p_DeltaMass0,
                                                         const double p_DeltaTime,
                                                         const bool   p_ForceRecalculate) {
    STELLAR_TYPE stellarType = m_StellarType;                                               // default is no change


    if (ShouldBeMasslessRemnant()) {                                                        // ALEJANDRO - 02/12/2016 - Attempt to fix updating the star if it lost all of its mass
        stellarType = STELLAR_TYPE::MASSLESS_REMNANT;                                       // JR: should also pick up already massless remnant
    }
    else {
        // RTW 30/06/20 - Why is this here?
        stellarType = ResolveSupernova();                                                   // handle supernova     JR: moved this to start of timestep
        
        if (stellarType == m_StellarType) {                                                 // still on phase?
            
            UpdateAttributesAndAgeOneTimestepPreamble(p_DeltaMass, p_DeltaMass0, p_DeltaTime);      // apply mass changes and save current values if required

            if (p_ForceRecalculate                     ||                                   // force recalculate?
                utils::Compare(p_DeltaMass,  0.0) != 0 ||                                   // mass change? or...
                utils::Compare(p_DeltaMass0, 0.0) != 0 ||                                   // mass0 change? or...
                utils::Compare(p_DeltaTime,  0.0)  > 0) {                                   // age/time advance?
                                                                                            // yes - update attributes
                AgeOneTimestepPreamble(p_DeltaTime);                                        // advance dt, age, simulation time

                if (ShouldSkipPhase()) {                                                    // skip phase?
                    stellarType = ResolveSkippedPhase();                                    // phase skipped - do what's required
                }
                else {                                                                      // not skipped - execute phase
                    stellarType = EvolveOnPhase();                                          // evolve on phase

                    // RTW 30/06/20 - why is this the flag for when ResolveEndOfPhase should get called?
                    if (stellarType == m_StellarType) {                                     // still on phase?
                        stellarType = ResolveEndOfPhase();                                  // check for need to move off phase
                    }
                }
            }
        }
    }

    return stellarType;                                                                     // stellar type to which star should evolve
}


/*
 * Set parameters required before ageing one timestep - modify star attributes
 *
 *
 * void AgeOneTimestepPreamble(const double p_DeltaTime)
 *
 * @param   [IN]    p_DeltaTime                 Timestep in Myr
 */
void BaseStar::AgeOneTimestepPreamble(const double p_DeltaTime) {

    if (p_DeltaTime > 0.0) {                        // if dt > 0    (don't use utils::Compare() here)
        m_Time += p_DeltaTime;                      // advance physical simulation time
        m_Age  += p_DeltaTime;                      // advance age of star
        m_Dt    = p_DeltaTime;                      // set timestep
    }

    EvolveOneTimestepPreamble();
}


/*
 * Evolve the star on it's current phase - take one timestep on the current phase
 *
 *
 * STELLAR_TYPE EvolveOnPhase()
 *
 * @return                                      Stellar Type to which star should evolve - unchanged if not moving off current phase
 */
STELLAR_TYPE BaseStar::EvolveOnPhase() {

    STELLAR_TYPE stellarType = m_StellarType;

    if (ShouldEvolveOnPhase()) {                                                    // Evolve timestep on phase

        m_Tau         = CalculateTauOnPhase();

        m_COCoreMass  = CalculateCOCoreMassOnPhase();
        m_CoreMass    = CalculateCoreMassOnPhase();
        m_HeCoreMass  = CalculateHeCoreMassOnPhase();

        m_Luminosity  = CalculateLuminosityOnPhase();

        std::tie(m_Radius, stellarType) = CalculateRadiusAndStellarTypeOnPhase();   // Radius and possibly new stellar type
        // RTW 30/06/20 - New stellar type only comes if its HeMS and needs to be updated
        // this should be made clearer, especially since it says OnPhase 

        ResolveEnvelopeMassOnPhase(m_Tau);

        m_Mu          = CalculatePerturbationMuOnPhase();

        PerturbLuminosityAndRadiusOnPhase();

        m_Temperature = CalculateTemperatureOnPhase();

        STELLAR_TYPE thisStellarType = ResolveEnvelopeLoss();                       // Resolve envelope loss if it occurs - possibly new stellar type
        if (thisStellarType != m_StellarType) {                                     // thisStellarType overrides stellarType (from CalculateRadiusAndStellarTypeOnPhase())
            stellarType = thisStellarType;
        }
    }

    return stellarType;
}


/*
 * Evolve the star onto the next phase if necessary - take one timestep at the end of the current phase
 *  // RTW 30/06/20 - what does it mean to take one timestep at the end?
 *
 *
 * STELLAR_TYPE ResolveEndOfPhase()
 *
 * @return                                      Stellar Type to which star should evolve - unchanged if not moving off current phase
 */
STELLAR_TYPE BaseStar::ResolveEndOfPhase() {

    STELLAR_TYPE stellarType = m_StellarType;

    // RTW test
    //std::cout << "ResolveEndOfPhase() is called" << std::endl;
    if (IsEndOfPhase()) {                                                       // End of phase

        // RTW test
        std::cout << "ResolveEndOfPhase() is used" << std::endl;
        stellarType = ResolveEnvelopeLoss();                                    // Resolve envelope loss if it occurs
        if (stellarType == m_StellarType) {                                     // Staying on phase?

            m_Tau         = CalculateTauAtPhaseEnd();

            m_COCoreMass  = CalculateCOCoreMassAtPhaseEnd();
            m_CoreMass    = CalculateCoreMassAtPhaseEnd();
            m_HeCoreMass  = CalculateHeCoreMassAtPhaseEnd();

            m_Luminosity  = CalculateLuminosityAtPhaseEnd();
            
            m_Radius      = CalculateRadiusAtPhaseEnd();

            ResolveEnvelopeMassAtPhaseEnd(m_Tau);

            m_Mu          = CalculatePerturbationMuAtPhaseEnd();

            PerturbLuminosityAndRadiusAtPhaseEnd();

            m_Temperature = CalculateTemperatureAtPhaseEnd();

            stellarType   = EvolveToNextPhase();                                // determine the stellar type to which the star should evolve
        }
    }

    return stellarType;
}


/************************/
/*   Star.cpp           */
/************************/



/*
 * Apply mass changes if required, age the star one timestep, advance the simulation time, and update the
 * attributes of the star.
 *
 * The star's attributes (Age, Radius, Luminosity etc.) are calculated and updated as required.
 *
 * Free parameters in the update process are the star's mass (m_Mass), initial mass (m_Mass0), the star's age
 * (m_Age) and the simulation time attribute (m_Time):
 *
 *    - if required, the star's mass is changed by the amount passed as the p_DeltaMass parameter before other
 *      attributes are updated.  The p_DeltaMass parameter may be zero, in which case no change is made to the
 *      star's mass before the attributes of the star are calculated.
 *
 *    - if required, the star's initial mass is changed by the amount passed as the p_DeltaMass0 parameter before
 *      other attributes are updated.  The p_DeltaMass parameter may be zero, in which case no change is made to
 *      the star's mass before the attributes of the star are calculated.  This should be used infrequently, and
 *      is really a kludge because the Mass0 attribute in Hurley et al. 2000 was overloaded after the introduction
 *      of mass loss (see section 7.1).  We should really separate the different uses of Mass0 in the code and
 *      use a different variable - initial mass shouldn't change (other than to initial mass upon entering a
 *      stellar phase - it doesn't make a lot of sense for initial mass to change during evolution through the
 *      phase).         JR: todo: action this paragraph.
 *
 *    - if required, the star is aged by the amount passed as the p_DeltaTime parameter, and the simulation time is
 *      advanced by the same amount, before other attributes are updated.  The p_deltaTime parameter may be zero,
 *      in which case no change is made to the star's age or the physical time attribute.
 *
 *
 * Checks whether the star:
 *    - is a massless remnant (checked after applying p_DeltaMass and p_DeltaMass0, but before applying p_DeltaTime)
 *    - has become a supernova (checked after applying p_DeltaMass and p_DeltaMass0, but before applying p_DeltaTime)
 *    - should skip this phase for this timestep (checked after applying p_DeltaMass, p_DeltaMass0 and p_DeltaTime)
 *
 * If none of the above are true the star evolves on phase for the specified timestep (which may be 0, in which case
 * the star's attributes other than age are re-calculated), then the need to evolve the star off phase is checked.
 *
 * If p_DeltaMass, p_DeltaMass0 and p_DeltaTime are all passed as zero the checks for massless remnant and supernova
 * are performed (and consequential changes made), but no other changes to the star's attributes are made - unless
 * the p_ForceRecalculate parameter is set true.        JR: todo: I'm not convinced p_ForceRecalculate is necessary - check it
 *
 * The functional return is the stellar type to which the star should evolve.  The returned stellar type is just the
 * stellar type of the star upon entry if it should remain on phase.
 *
 * If the parameter p_Switch is true the star will switch to the new stellar type before returning.
 *
 *
 * STELLAR_TYPE UpdateAttributesAndAgeOneTimestep(const double p_DeltaMass,
 *                                                const double p_DeltaMass0,
 *                                                const double p_DeltaTime,
 *                                                const bool   p_Switch,
 *                                                const bool   p_ForceRecalculate)
 *
 * @param   [IN]    p_DeltaMass                 The change in mass to apply in Msol
 * @param   [IN]    p_DeltaMass0                The change in mass0 to apply in Msol
 * @param   [IN]    p_DeltaTime                 The timestep to evolve in Myr
 * @param   [IN]    p_Switch                    Specifies whether the star should switch to new stellar type before returning
 *                                              (optional, default = true)
 * @param   [IN]    p_ForceRecalculate          Specifies whether the star's attributes should be recalculated even if the three deltas are 0.0
 *                                              (optional, default = false)
 * @return                                      New stellar type for star
 */
STELLAR_TYPE Star::UpdateAttributesAndAgeOneTimestep(const double p_DeltaMass,
                                                     const double p_DeltaMass0,
                                                     const double p_DeltaTime,
                                                     const bool   p_Switch,
                                                     const bool   p_ForceRecalculate) {

    STELLAR_TYPE stellarType = m_Star->UpdateAttributesAndAgeOneTimestep(p_DeltaMass, p_DeltaMass0, p_DeltaTime, p_ForceRecalculate);

    if (p_Switch && (stellarType != m_Star->StellarType())) {                               // switch to new stellar type if necessary?
        STELLAR_TYPE stellarTypePrev = SwitchTo(stellarType);                               // yes - switch
        m_Star->SetStellarTypePrev(stellarTypePrev);                                        // record previous stellar type
        m_Star->CalculateAllTimescales();                                                   // calculate dynamical, thermal, nuclear and radial expansion timescales
        
        // recalculate stellar attributes after switching if necessary - transition may not be continuous (e.g. CH -> HeMS)
        // (this could get recursive, but shouldn't...)
        if (stellarTypePrev == STELLAR_TYPE::CHEMICALLY_HOMOGENEOUS && 
                stellarType == STELLAR_TYPE::NAKED_HELIUM_STAR_MS) {                        // discontinuous transition?
            UpdateAttributes(0.0, 0.0, true);                                               // yes - recalculate stellar attributes
        }
    }

    return stellarType;                                                                     // return new stellar type
}


/*
 * Update the attributes of the star without ageing the star or advancing the simulation time.
 *
 * Apply mass changes as required and update the attributes of the star.
 *
 * The star's attributes (Radius, Luminosity etc. (not Age or simulation time)) are calculated and updated as
 * required.
 *
 * Free parameters in the update process are the star's mass (m_Mass) and initial mass (m_Mass0):
 *
 *    - if required, the star's mass is changed by the amount passed as the p_DeltaMass parameter before other
 *      attributes are updated.  The p_DeltaMass parameter may be zero, in which case no change is made to the
 *      star's mass before the attributes of the star are calculated.
 *
 *    - if required, the star's initial mass is changed by the amount passed as the p_DeltaMass parameter before
 *      other attributes are updated.  The p_DeltaMass parameter may be zero, in which case no change is made to
 *      the star's mass before the attributes of the star are calculated.  This should be used infrequently, and
 *      is really a kludge because the Mass0 attribute in Hurley et al. 2000 was overloaded after the introduction
 *      of mass loss (see section 7.1).  We should really separate the different uses of Mass0 in the code and
 *      use a different variable - initial mass shouldn't change (other than to initial mass upon entering a
 *      stellar phase - it doesn't make a lot of sense for initial mass to change during evolution through the
 *      phase).         JR: todo: action this paragraph.
 *
 *
 * Checks whether the star:
 *    - is a massless remnant (checked after applying p_DeltaMass and p_DeltaMass0)
 *    - has become a supernova (checked after applying p_DeltaMass and p_DeltaMass0)
 *    - should skip this phase for this timestep (checked after applying p_DeltaMass and p_DeltaMass0)
 *
 * If none of the above are true the star's attributes are updated based on the mass changes required and the star's
 * current phase, then the need to evolve the star off phase is checked.
 *
 * If p_DeltaMass and p_DeltaMass0 are both passed as zero the checks for massless remnant and supernova are performed
 * (and consequential changes made), but no other changes to the star's attributes are made - unless the p_ForceRecalculate
 * parameter is set true.        JR: todo: I'm not convinced p_ForceRecalculate is necessary - check it
 *
 * The functional return is the stellar type to which the star has evolved (if it evolved off the current phase).
 * The returned stellar type is just the stellar type of the star upon entry if it remained on the current phase.
 *
 *
 * STELLAR_TYPE UpdateAttributes(const double p_DeltaMass, const double p_DeltaMass0, const bool p_ForceRecalculate)
 *
 * @param   [IN]    p_DeltaMass                 The change in mass to apply in Msol
 * @param   [IN]    p_DeltaMass0                The change in mass0 to apply in Msol
 * @param   [IN]    p_ForceRecalculate          Specifies whether the star's attributes should be recalculated even if the two deltas are 0.0
 *                                              (optional, default = false)
 * @return                                      New stellar type for star
 */
STELLAR_TYPE Star::UpdateAttributes(const double p_DeltaMass, const double p_DeltaMass0, const bool p_ForceRecalculate) {
   return UpdateAttributesAndAgeOneTimestep(p_DeltaMass, p_DeltaMass0, 0.0, true, p_ForceRecalculate);
}


/*
 * Age the star a single timestep - timestep is provided as parameter
 *
 * Checks whether the star:
 *    - is a massless remnant
 *    - has become a supernova
 *    - should skip this phase for this timestep
 *
 * If none of the above are true, the star evolves on phase for this timestep,
 * then the need to evolve the star off phase is checked.
 *
 * The functional return is the stellar type to which the star should evolve.  The returned stellar
 * type is just the stellar type of the star upon entry if it should remain on phase.
 *
 * If the parameter p_Switch is true the star will switch to the new stellar type before returning
 *
 *
 * STELLAR_TYPESTELLAR_TYPE AgeOneTimestep(const double p_Dt, bool p_Switch)
 *
 * @param   [IN]    p_DeltaTime                 The timestep to evolve in Myr
 * @param   [IN]    p_Switch                    Specifies whether the star should switch to new stellar type before returning
 *                                              (optional, default = true)
 * @return                                      New stellar type for star
 */
STELLAR_TYPE Star::AgeOneTimestep(const double p_DeltaTime, bool p_Switch) {
    return UpdateAttributesAndAgeOneTimestep(0.0, 0.0, p_DeltaTime, p_Switch, false);
}


/*
 * Evolve the star a single timestep - suggested timestep is provided
 *
 *    - save current state
 *    - age star timestep
 *    - if ageing caused too much change, revert state, shorten timestep and try again
 *    - loop until suitable timestep found (and taken)
 *
 * The functional return is the timestep actually taken (in Myr)
 *
 * double EvolveOneTimestep(const double p_Dt)
 *
 * @param   [IN]    p_Dt                        The suggested timestep to evolve
 * @return                                      The timestep actually taken
 */
double Star::EvolveOneTimestep(const double p_Dt) {

    double       dt = p_Dt;

    STELLAR_TYPE stellarType;

    bool         takeTimestep = false;
    int          retryCount   = 0;

    while (!takeTimestep) {                                                                                     // do this until a suitable timestep is found (or the maximum retry count is reached)

        SaveState();                                                                                            // save the state of the star - in case we want to revert

        double minTimestep = std::max(m_Star->CalculateDynamicalTimescale(), ABSOLUTE_MINIMUM_TIMESTEP);        // calculate the minimum timestep - maximum of dynamical timescale for this star and the aboslute minimum timestep

        // evolve the star a single timestep

        stellarType = AgeOneTimestep(dt, false);                                                                // age the star one time step - modify stellar attributes as appropriate, but do not switch stellar type

        // determine if the timestep should be taken
        // don't take the timestep if we stepped too far

        takeTimestep = true;                                                                                    // flag to determine if the timestep should be taken
        if (utils::Compare(m_Star->CalculateRadialChange(), MAXIMUM_RADIAL_CHANGE) >= 0) {                      // too much change?
            if (utils::Compare(dt, minTimestep) <= 0) {                                                         // yes - already at or below minimum timestep?
                takeTimestep = true;                                                                            // yes - just take the last timestep
//                if (!OPTIONS->Quiet()) SAY("WARNING: Timestep below minimum - timestep taken!");                // announce the problem if required and plough on regardless...
            }
            else {                                                                                              // not at or below dynamical - reduce timestep and try again
                retryCount++;                                                                                   // increment retry count
                if (retryCount > MAX_TIMESTEP_RETRIES) {                                                        // too many retries?
                    takeTimestep = true;                                                                        // yes - take the last timestep anyway
//                    if (!OPTIONS->Quiet()) SAY("WARNING: Too many retries finding timestep - timestep taken!"); // announce the problem if required and plough on regardless...
                }
                else {                                                                                          // not too many retries - retry with smaller timestep
                    if (RevertState()) {                                                                        // revert to last state ok?
                        dt = dt / 2.0;                                                                          // yes - halve the timestep (limit to minimum)      JR: probably should be dt = max(dt / 2.0, minTimestep);
                        takeTimestep = false;                                                                   // previous timestep discared - use new one
                    }
                    else {                                                                                      // revert failed
                        takeTimestep = true;                                                                    // take the last timestep anyway
//                        if (!OPTIONS->Quiet()) SAY("Revert of star failed - timestep taken!");                  // announce the problem if required and plough on regardless...
                    }
                }
            }
        }
    }

    // take the timestep

    (void)SwitchTo(stellarType);                                                                                // switch phase if required  JR: whether this goes before or after the log record is a little problematic, but in the end probably doesn't matter too much

    (void)m_Star->ResolveMassLoss();                                                                            // apply wind mass loss if required     JR: should this really be before the call to SwitchTo()?  It isn't in the original code

     m_Star->CalculateAllTimescales();                                                                          // calculate dynamical, thermal, nuclear and radial expansion timescales

    return dt;                                                                                                  // return the timestep actually taken
}


/************************/
/*  BaseBinaryStar.cpp  */
/************************/



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
            // equilibrate masses and circularise (check for merger is done later)

            if (utils::Compare(m_Star1->Mass(), m_Star2->Mass()) != 0) {                                                        // masses already equal?
                                                                                                                                // no - makethem equal
                STELLAR_TYPE stellarType1 = m_Star1->StellarType();                                                             // star 1 stellar type before updating attributes
                STELLAR_TYPE stellarType2 = m_Star2->StellarType();                                                             // star 2 stellar type before updating attributes

                double mass = (m_Star1->Mass() + m_Star2->Mass()) / 2.0;                                                        // share mass equally
                if ((m_Star1->UpdateAttributes(mass - m_Star1->Mass(), mass - m_Star1->Mass0(), true) != stellarType1) ||       // set new mass, mass0 for star 1
                    (m_Star2->UpdateAttributes(mass - m_Star2->Mass(), mass - m_Star2->Mass0(), true) != stellarType2)) {       // set new mass, mass0 for star 2
                    m_PrintExtraDetailedOutput = true;                                                                          // print detailed output record if stellar type changed
                }
                m_MassesEquilibrated = true;                                                                                    // record that we've equilbrated
            }

            // circularise if not already

            if (utils::Compare(m_Eccentricity, 0.0) != 0) {                                                                     // eccentricity = 0.0?
                                                                                                                                // no - circularise
                // conserve angular momentum
                // use J = m1 * m2 * sqrt(G * a * (1 - e^2) / (m1 + m2))

                double M              = m_Star1->Mass() + m_Star2->Mass();
                double m1m2           = m_Star1->Mass() * m_Star2->Mass();
                m_SemiMajorAxisPrime *= 16.0 * m1m2 * m1m2 / (M * M * M * M) * (1.0 - (m_Eccentricity * m_Eccentricity));       // circularise; conserve angular momentum
                m_Eccentricity        = 0.0;                                                                                    // now circular
            }

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

    STELLAR_TYPE stellarType1 = m_Star1->StellarTypePrev();                                                 // star 1 stellar type before updating attributes
    STELLAR_TYPE stellarType2 = m_Star2->StellarTypePrev();                                                 // star 2 stellar type before updating attributes

    // update mass of star1 according to mass loss and mass transfer, then update age accordingly
    (void)m_Star1->UpdateAttributes(m_Star1->MassPrev() - m_Star1->Mass() + m_Star1->MassLossDiff() + m_Star1->MassTransferDiff(), 0.0);    // update mass for star1
    m_Star1->UpdateInitialMass();                                                                       // update initial mass of star1 (MS, HG & HeMS)  JR: todo: fix this kludge one day - mass0 is overloaded, and isn't always "initial mass"
    m_Star1->UpdateAgeAfterMassLoss();                                                                  // update age of star1
    m_Star1->ApplyMassTransferRejuvenationFactor();                                                     // apply age rejuvenation factor for star1

    // rinse and repeat for star2
    (void)m_Star2->UpdateAttributes(m_Star2->MassPrev() - m_Star2->Mass() + m_Star2->MassLossDiff() + m_Star2->MassTransferDiff(), 0.0);    // update mass for star2
    m_Star2->UpdateInitialMass();                                                                       // update initial mass of star 2 (MS, HG & HeMS)  JR: todo: fix this kludge one day - mass0 is overloaded, and isn't always "initial mass"
    m_Star2->UpdateAgeAfterMassLoss();                                                                  // update age of star2
    m_Star2->ApplyMassTransferRejuvenationFactor();                                                     // apply age rejuvenation factor for star2
   
    if ((m_Star1->StellarType() != stellarType1) || (m_Star2->StellarType() != stellarType2)) {         // stellar type change?
        m_PrintExtraDetailedOutput = true;                                                              // yes - print detailed output record
    }

    // update binary
    m_OrbitalVelocityPrime = m_OrbitalVelocityPrev + m_OmegaMassLossDiff + m_OmegaMassTransferDiff;     // should here be a diff quantity because of MB?    JR: todo: ?
    m_SemiMajorAxisPrime   = m_SemiMajorAxisPrev + m_aMassLossDiff + m_aMassTransferDiff;

    // if CHE enabled, update rotational frequency for constituent stars - assume tidally locked
    if (OPTIONS->CHE_Option() != CHE_OPTION::NONE) m_Star1->SetOmega(m_OrbitalVelocityPrime);
    if (OPTIONS->CHE_Option() != CHE_OPTION::NONE) m_Star2->SetOmega(m_OrbitalVelocityPrime);

    CalculateEnergyAndAngularMomentum();                                                                // perform energy and angular momentum calculations
}


/*
 * Evaluate the binary system
 *
 *    - caclulate any mass transfer
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

    // RTW 30/06/20 - Do we want to make it impossible for two type-shifts in one timestep? If it happens, that suggests other problems, maybe we should flag it instead?
    
    EvaluateBinaryPreamble();                                                                                           // get things ready - do some house-keeping

    CheckMassTransfer(p_Dt);                                                                                            // calculate mass transfer if necessary

    // RTW 30/06/20 - why is this in binary evolution?
    CalculateWindsMassLoss();                                                                                           // calculate mass loss dues to winds

    if ((m_CEDetails.CEEnow || m_StellarMerger) &&                                                                      // CEE or merger?
       !(OPTIONS->CHE_Option() != CHE_OPTION::NONE && HasTwoOf({STELLAR_TYPE::CHEMICALLY_HOMOGENEOUS}))) {              // yes - avoid CEE if CH+CH
        ResolveCommonEnvelopeEvent();                                                                                   // resolve CEE - immediate event
    }
    else if (m_Star1->IsSNevent() || m_Star2->IsSNevent()) {
        EvaluateSupernovae(true);                                                                                       // evaluate supernovae (both stars) - immediate event
    }
    else {
        // Dangerous - actually changes the mass of the stars!
        ResolveMassChanges();                                                                                           // apply mass loss and mass transfer as necessary
    }

    STELLAR_TYPE stellarType1 = m_Star1->StellarType();                                                                 // star 1 stellar type before updating attributes
    STELLAR_TYPE stellarType2 = m_Star2->StellarType();                                                                 // star 2 stellar type before updating attributes

    if ((m_Star1->UpdateAttributes(0.0, 0.0, true) != stellarType1) ||                                                  // recalculate stellar attributes for star 1
        (m_Star2->UpdateAttributes(0.0, 0.0, true) != stellarType2)) {                                                  // recalculate stellar attributes for star 2
        m_PrintExtraDetailedOutput = true;                                                                              // print detailed output record if stellar type changed
    }

    if (m_PrintExtraDetailedOutput == true) { PrintDetailedOutput(m_Id); }                                              // print detailed output record if stellar type changed
    m_PrintExtraDetailedOutput = false;                                                                                 // reset detailed output printing flag for the next timestep

    EvaluateSupernovae(false);                                                                                          // evaluate supernovae (both stars)

    // assign new values to "previous" values, for following timestep
    m_EccentricityPrev	  = m_EccentricityPrime;
    m_SemiMajorAxisPrev   = m_SemiMajorAxisPrime;
    m_OrbitalVelocityPrev = m_OrbitalVelocityPrime;    

    CalculateEnergyAndAngularMomentum();                                                                                // perform energy and angular momentum calculations
    
    // RTW 30/06/20 - TODO need to double check that this is in the correct location
    m_TotalAngularMomentumPrev = m_TotalAngularMomentumPrime;   // Is this line ok here?        JR: todo - this probably should be in evaluateBinary(), except that evaluateBinary() may not be executed at each timestep - maybe this has to stay here


    if (!(m_Star1->IsOneOf({ STELLAR_TYPE::MASSLESS_REMNANT })))
        m_Star1->UpdateMagneticFieldAndSpin(m_CEDetails.CEEnow, m_Dt * MYR_TO_YEAR * SECONDS_IN_YEAR, EPSILON_PULSAR);  // update pulsar parameters for star1
    if (!(m_Star2->IsOneOf({ STELLAR_TYPE::MASSLESS_REMNANT })))
        m_Star2->UpdateMagneticFieldAndSpin(m_CEDetails.CEEnow, m_Dt * MYR_TO_YEAR * SECONDS_IN_YEAR, EPSILON_PULSAR);  // update pulsar parameters for star2
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

