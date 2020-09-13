#include "Star.h"
#include <algorithm>
#include <csignal>

// Default constructor
Star::Star() : m_Star(new BaseStar()) {

    m_ObjectId   = globalObjectId++;                                                                // set object id
    m_ObjectType = OBJECT_TYPE::STAR;                                                               // set object type

    m_SaveStar = nullptr;
}


// Regular constructor - with parameters for RandomSeed, MZAMS, Metallicity, LBVFactor and WolfRayetFactor
Star::Star(const unsigned long int p_RandomSeed,
           const double            p_MZAMS,
           const double            p_Metallicity,
           const KickParameters    p_KickParameters,
           const double            p_LBVfactor,
           const double            p_WolfRayetFactor) {

    m_ObjectId   = globalObjectId++;                                                                                // set object id
    m_ObjectType = OBJECT_TYPE::STAR;                                                                               // set object type

    m_Star = new BaseStar(p_RandomSeed, p_MZAMS, p_Metallicity, p_KickParameters, p_LBVfactor, p_WolfRayetFactor);  // create underlying BaseStar object

    // star begins life as a main sequence star, unless it is
    // spinning fast enough for it to be chemically homogeneous

    if (OPTIONS->CHE_Option() != CHE_OPTION::NONE && utils::Compare(m_Star->Omega(), m_Star->OmegaCHE()) >= 0) {    // CHE?
        (void)SwitchTo(STELLAR_TYPE::CHEMICALLY_HOMOGENEOUS, true);                                                 // yes
    }
    else if (p_MZAMS <= 0.7) {                                                                                      // no - MS - initial mass determines actual type  JR: don't use utils::Compare() here
        (void)SwitchTo(STELLAR_TYPE::MS_LTE_07, true);                                                              // MS <= 0.0 Msol
    }
    else {
        (void)SwitchTo(STELLAR_TYPE::MS_GT_07, true);                                                               // MS > 0.7 Msol
    }

    m_SaveStar = nullptr;
}


/*
 * Clone underlying BaseStar
 *
 * Instantiates new object of current underlying star class and initialises
 * it with the star object passed as p_Star
 *
 *
 * BaseStar* Clone()
 */
BaseStar* Star::Clone(const BaseStar& p_Star) {

    BaseStar *ptr = nullptr;

    switch (p_Star.StellarType()) {
        case STELLAR_TYPE::MS_LTE_07                                : {ptr = new MS_lte_07(p_Star, false);} break;
        case STELLAR_TYPE::MS_GT_07                                 : {ptr = new MS_gt_07(p_Star, false);} break;
        case STELLAR_TYPE::CHEMICALLY_HOMOGENEOUS                   : {ptr = new CH(p_Star, false);} break;
        case STELLAR_TYPE::HERTZSPRUNG_GAP                          : {ptr = new HG(p_Star, false);} break;
        case STELLAR_TYPE::FIRST_GIANT_BRANCH                       : {ptr = new FGB(p_Star, false);} break;
        case STELLAR_TYPE::CORE_HELIUM_BURNING                      : {ptr = new CHeB(p_Star, false);} break;
        case STELLAR_TYPE::EARLY_ASYMPTOTIC_GIANT_BRANCH            : {ptr = new EAGB(p_Star, false);} break;
        case STELLAR_TYPE::THERMALLY_PULSING_ASYMPTOTIC_GIANT_BRANCH: {ptr = new TPAGB(p_Star, false);} break;
        case STELLAR_TYPE::NAKED_HELIUM_STAR_MS                     : {ptr = new HeMS(p_Star, false);} break;
        case STELLAR_TYPE::NAKED_HELIUM_STAR_HERTZSPRUNG_GAP        : {ptr = new HeHG(p_Star, false);} break;
        case STELLAR_TYPE::NAKED_HELIUM_STAR_GIANT_BRANCH           : {ptr = new HeGB(p_Star, false);} break;
        case STELLAR_TYPE::HELIUM_WHITE_DWARF                       : {ptr = new HeWD(p_Star, false);} break;
        case STELLAR_TYPE::CARBON_OXYGEN_WHITE_DWARF                : {ptr = new COWD(p_Star, false);} break;
        case STELLAR_TYPE::OXYGEN_NEON_WHITE_DWARF                  : {ptr = new ONeWD(p_Star, false);} break;
        case STELLAR_TYPE::NEUTRON_STAR                             : {ptr = new NS(p_Star, false);} break;
        case STELLAR_TYPE::BLACK_HOLE                               : {ptr = new BH(p_Star, false);} break;
        case STELLAR_TYPE::MASSLESS_REMNANT                         : {ptr = new MR(p_Star, false);} break;
        default: break;                                             // avoids compiler warning - this should never happen
    }

    return ptr;
}

// Copy constructor - deep copy so dynamic variables are also copied
Star::Star(const Star& p_Star) {

    m_ObjectId   = globalObjectId++;                                                    // set object id
    m_ObjectType = OBJECT_TYPE::STAR;                                                   // set object type

    m_Star     = p_Star.m_Star ? Clone(*(p_Star.m_Star)) : nullptr;                     // copy underlying BasStar object
    m_SaveStar = p_Star.m_SaveStar ? Clone(*(p_Star.m_Star)) : nullptr;                 // and the saved copy
}


// Assignment overload - deep copy so dynamic variables are also copied
Star& Star::operator = (const Star& p_Star) {

    if (this != &p_Star) {                                                              // make sure we're not not copying ourselves...

        m_ObjectId   = globalObjectId++;                                                // set object id
        m_ObjectType = OBJECT_TYPE::STAR;                                               // set object type

        delete m_Star;
        m_Star = p_Star.m_Star ? Clone(*(p_Star.m_Star)) : nullptr;                     // copy underlying BasStar object

        delete m_SaveStar;
        m_SaveStar = p_Star.m_SaveStar ? Clone(*(p_Star.m_SaveStar)) : nullptr;         // and the saved copy
    }
    return *this;
}


/*
 * Switch to required star type
 *
 * Instantiates new object of required class, deletes existing pointer to star object and
 * replaces it with pointer to newly instantiated object
 *
 *
 * STELLAR_TYPE SwitchTo(const STELLAR_TYPE p_StellarType, bool p_SetInitialState)
 *
 * @param   [IN]    p_StellarType               StellarType to switch to
 * @param   [IN]    p_SetInitialType            Indicates whether the initial stellar type of the star should be set to p_StellarType
 *                                              (optional, default = false)
 * @return                                      Stellar type of star before switch (previous stellar type)
 */
STELLAR_TYPE Star::SwitchTo(const STELLAR_TYPE p_StellarType, bool p_SetInitialType) {

    STELLAR_TYPE stellarTypePrev = m_Star->StellarType();

    // don't switch if stellarTypePrev == p_StellarType
    // (the call to SwitchTo() in Star::EvolveOneTimestep() doesn't check - it relies on the check here)

    if (p_StellarType != m_Star->StellarType()) {
        BaseStar *ptr = nullptr;

        switch (p_StellarType) {
            case STELLAR_TYPE::MS_LTE_07                                : {ptr = new MS_lte_07(*m_Star);} break;
            case STELLAR_TYPE::MS_GT_07                                 : {ptr = new MS_gt_07(*m_Star);} break;
            case STELLAR_TYPE::CHEMICALLY_HOMOGENEOUS                   : {ptr = new CH(*m_Star);} break;
            case STELLAR_TYPE::HERTZSPRUNG_GAP                          : {ptr = new HG(*m_Star);} break;
            case STELLAR_TYPE::FIRST_GIANT_BRANCH                       : {ptr = new FGB(*m_Star);} break;
            case STELLAR_TYPE::CORE_HELIUM_BURNING                      : {ptr = new CHeB(*m_Star);} break;
            case STELLAR_TYPE::EARLY_ASYMPTOTIC_GIANT_BRANCH            : {ptr = new EAGB(*m_Star);} break;
            case STELLAR_TYPE::THERMALLY_PULSING_ASYMPTOTIC_GIANT_BRANCH: {ptr = new TPAGB(*m_Star);} break;
            case STELLAR_TYPE::NAKED_HELIUM_STAR_MS                     : {ptr = new HeMS(*m_Star);} break;
            case STELLAR_TYPE::NAKED_HELIUM_STAR_HERTZSPRUNG_GAP        : {ptr = new HeHG(*m_Star);} break;
            case STELLAR_TYPE::NAKED_HELIUM_STAR_GIANT_BRANCH           : {ptr = new HeGB(*m_Star);} break;
            case STELLAR_TYPE::HELIUM_WHITE_DWARF                       : {ptr = new HeWD(*m_Star);} break;
            case STELLAR_TYPE::CARBON_OXYGEN_WHITE_DWARF                : {ptr = new COWD(*m_Star);} break;
            case STELLAR_TYPE::OXYGEN_NEON_WHITE_DWARF                  : {ptr = new ONeWD(*m_Star);} break;
            case STELLAR_TYPE::NEUTRON_STAR                             : {ptr = new NS(*m_Star);} break;
            case STELLAR_TYPE::BLACK_HOLE                               : {ptr = new BH(*m_Star);} break;
            case STELLAR_TYPE::MASSLESS_REMNANT                         : {ptr = new MR(*m_Star);} break;
            default: break;                                             // avoids compiler warning - this should never happen
        }

        if (ptr) {
            delete m_Star;
            m_Star = ptr;

            if (p_SetInitialType) m_Star->SetInitialType(p_StellarType);
        }

        // write to switch log file if required

        if (utils::IsOneOf(stellarTypePrev, EVOLVABLE_TYPES)) {                             // star should be evolving from one of the evolvable types (We don't want the initial switch from Star->MS.  Not necessary for BSE (handled differently), but no harm)
        
            if (!OPTIONS->SingleStar() && OPTIONS->BSESwitchLog()) {                        // BSE and BSE Switch Log enabled?
                LOGGING->SetSwitchParameters(m_ObjectId, stellarTypePrev, p_StellarType);   // yes - store switch details to LOGGING service
                raise(SIGUSR1);                                                             // signal to BSE that switch is occurring
            }

            if (OPTIONS->SingleStar() && OPTIONS->SSESwitchLog()) {                         // SSE and SSE Switch Log enabled?
                LOGGING->SetSwitchParameters(m_ObjectId, stellarTypePrev, p_StellarType);   // yes - store switch details to LOGGING service
                m_Star->PrintSwitchLog(m_Id);                                               // no need for the BSE signal shenaningans - just call the function
            }
        }
    }

    return stellarTypePrev;
}


/*
 * Save current state of star
 *
 * Instantiates new object of current star class, deletes existing pointer to saved star
 * object (if it exists) and replaces it with pointer to newly instantiated object
 *
 *
 * void SaveState()
 */
void Star::SaveState() {

    delete m_SaveStar;
    m_SaveStar = Clone(*m_Star);
}


/*
 * Revert to the saved state of the star
 *
 * Changes the current star pointer (m_Star) to point to the save star object and
 * set the saved star pointer (m_SaveStar) to null.  Setting the saved state pointer
 * to null means there will be no saved state after calling this function - so state
 * needs to be saved if necessary (I don't do it here because we may not need to).
 *
 *
 * bool RevertState()
 *
 * @return                                      Boolean flag indicating success/failure (true = success)
 */
bool Star::RevertState() {
    bool result = false;

    if (m_SaveStar) {
        delete m_Star;
        m_Star     = m_SaveStar;
        m_SaveStar = nullptr;
        result     = true;
    }

    return result;
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

    m_Star->PrintStashedSupernovaDetails();                                                                     // print stashed SSE Supernova log record if necessary

    (void)SwitchTo(stellarType);                                                                                // switch phase if required  JR: whether this goes before or after the log record is a little problematic, but in the end probably doesn't matter too much

    (void)m_Star->ResolveMassLoss();                                                                            // apply wind mass loss if required     JR: should this really be before the call to SwitchTo()?  It isn't in the original code

    return dt;                                                                                                  // return the timestep actually taken
}


/*
 * Evolve the star through its entire lifetime
 *
 *
 * void Evolve(const long int p_Id)
 *
 * @param   [IN]    p_Id                        The mass id (e.g. step number) for this star - can be used to name logfiles for this star
 */
void Star::Evolve(const long int p_Id) {

    m_Id = p_Id;                                    // store the id

    // evolve the star

    m_Star->CalculateGBParams();                    // calculate giant branch parameters - in case for some reason star is initially not MS

    double dt = 0.0;

    // JR: todo:
    // m_Error seems to be set ad hoc for SSE, and doesn't actually stop the evolution
    // we should be more rigorous in checking/setting error conditions, and stop the evolution for catastrophic errors

    while (m_Star->IsOneOf({ STELLAR_TYPE::MS_LTE_07, STELLAR_TYPE::MS_GT_07, STELLAR_TYPE::HERTZSPRUNG_GAP,
                             STELLAR_TYPE::FIRST_GIANT_BRANCH, STELLAR_TYPE::CORE_HELIUM_BURNING, STELLAR_TYPE::EARLY_ASYMPTOTIC_GIANT_BRANCH, STELLAR_TYPE::THERMALLY_PULSING_ASYMPTOTIC_GIANT_BRANCH,
                             STELLAR_TYPE::NAKED_HELIUM_STAR_MS, STELLAR_TYPE::NAKED_HELIUM_STAR_HERTZSPRUNG_GAP, STELLAR_TYPE::NAKED_HELIUM_STAR_GIANT_BRANCH })) {

        dt = m_Star->CalculateTimestep();           // calculate new timestep

        dt = EvolveOneTimestep(dt);                 // evolve for timestep

        m_Star->PrintParameters(m_Id);              // log record  JR: this should probably be before the star switches type, but this way matches the original code
    }

}
