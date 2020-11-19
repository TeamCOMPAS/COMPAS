#include "HeHG.h"
#include "HeGB.h"
#include "HeWD.h"


/*
 * Calculate timescales in units of Myr
 *
 * Timescales depend on a star's mass, so this needs to be called at least each timestep
 *
 * Vectors are passed by reference here for performance - preference would be to pass const& and
 * pass modified value back by functional return, but this way is faster - and this function is
 * called many, many times.
 *
 *
 * void CalculateTimescales(const double p_Mass, DBL_VECTOR &p_Timescales)
 *
 * @param   [IN]        p_Mass                  Mass in Msol
 * @param   [IN/OUT]    p_Timescales            Timescales
 */
void HeHG::CalculateTimescales(const double p_Mass, DBL_VECTOR &p_Timescales) {
#define timescales(x) p_Timescales[static_cast<int>(TIMESCALE::x)]  // for convenience and readability - undefined at end of function
#define gbParams(x) m_GBParams[static_cast<int>(GBP::x)]            // for convenience and readability - undefined at end of function

    HeMS::CalculateTimescales(p_Mass, p_Timescales);    // calculate common values

    double p1   = gbParams(p) - 1.0;
    double q1   = gbParams(q) - 1.0;
    double p1_p = p1 / gbParams(p);
    double q1_q = q1 / gbParams(q);

    double LTHe = HeMS::CalculateLuminosityAtPhaseEnd(p_Mass);

    timescales(tinf1_HeGB) = timescales(tHeMS) + (1.0 / ((p1 * gbParams(AHe) * gbParams(D))) * PPOW((gbParams(D) / LTHe), p1_p));
    timescales(tx_HeGB) = timescales(tinf1_HeGB) - (timescales(tinf1_HeGB) - timescales(tHeMS)) * PPOW((LTHe / gbParams(Lx)), p1_p);
    timescales(tinf2_HeGB) = timescales(tx_HeGB) + ((1.0 / (q1 * gbParams(AHe) * gbParams(B))) * PPOW((gbParams(B) / gbParams(Lx)), q1_q));

#undef gbParams
#undef timescales
}


/*
 * Calculate Giant Branch (GB) parameters per Hurley et al. 2000
 *
 * Giant Branch Parameters depend on a star's mass, so this needs to be called at least each timestep
 *
 * Vectors are passed by reference here for performance - preference would be to pass const& and
 * pass modified value back by functional return, but this way is faster - and this function is
 * called many, many times.
 *
 * void CalculateGBParams(const double p_Mass, DBL_VECTOR &p_GBParams)
 *
 * @param   [IN]        p_Mass                  Mass in Msol
 * @param   [IN/OUT]    p_GBParams              Giant Branch Parameters - calculated here
 */
void HeHG::CalculateGBParams(const double p_Mass, DBL_VECTOR &p_GBParams) {
#define gbParams(x) p_GBParams[static_cast<int>(GBP::x)]    // for convenience and readability - undefined at end of function
    GiantBranch::CalculateGBParams(p_Mass, p_GBParams);                         // calculate common values (actually, all)

    // recalculate HeHG specific values

	gbParams(B)      = CalculateCoreMass_Luminosity_B_Static();
	gbParams(D)      = CalculateCoreMass_Luminosity_D_Static(p_Mass);

    gbParams(Mx)     = GiantBranch::CalculateCoreMass_Luminosity_Mx_Static(p_GBParams);      // depends on B, D, p & q - recalculate if any of those are changed
    gbParams(Lx)     = GiantBranch::CalculateCoreMass_Luminosity_Lx_Static(p_GBParams);      // JR: Added this - depends on B, D, p, q & Mx - recalculate if any of those are changed

	gbParams(McBAGB) = CalculateCoreMassAtBAGB();
	gbParams(McBGB)  = CalculateCoreMassAtBGB(p_Mass, p_GBParams);

    gbParams(McSN)   = CalculateCoreMassAtSupernova_Static(gbParams(McBAGB));   // JR: Added this

#undef gbParams
}


/*
 * Calculate Giant Branch (GB) parameters per Hurley et al. 2000
 *
 * Giant Branch Parameters depend on a star's mass, so this needs to be called at least each timestep
 *
 * Vectors are passed by reference here for performance - preference would be to pass const& and
 * pass modified value back by functional return, but this way is faster - and this function is
 * called many, many times.
 *
 *
 * This function exists to facilitate the calculation of gbParams in EAGB::ResolveEnvelopeLoss() for the
 * stellar type to which the star will evolve.  The calculations of some of the stellar attributes there
 * depend on new gbParams.  
 * JR: This really needs to be revisited one day - these calculations should really be performed after
 *     switching to the new stellar type, but other calculations are done (in the legacy code) before the 
 *     switch (see evolveOneTimestep() in star.cpp for EAGB stars in the legacy code).
 *
 *
 * void CalculateGBParams_Static(const double      p_Mass0, 
 *                               const double      p_Mass, 
 *                               const double      p_LogMetallicityXi, 
 *                               const DBL_VECTOR &p_MassCutoffs, 
 *                               const DBL_VECTOR &p_AnCoefficients, 
 *                               const DBL_VECTOR &p_BnCoefficients,
 *                                     DBL_VECTOR &p_GBParams)
 *
 * @param   [IN]        p_Mass0                 Mass0 in Msol
 * @param   [IN]        p_Mass                  Mass in Msol
 * @param   [IN]        p_LogMetallicityXi      log10(Metallicity / Zsol) - called xi in Hurley et al 2000
 * @param   [IN]        p_MassCutoffs           Mass cutoffs
 * @param   [IN]        p_AnCoefficients        a(n) coefficients
 * @param   [IN]        p_BnCoefficients        b(n) coefficients
 * @param   [IN/OUT]    p_GBParams              Giant Branch Parameters - calculated here
 */
void HeHG::CalculateGBParams_Static(const double      p_Mass0, 
                                    const double      p_Mass, 
                                    const double      p_LogMetallicityXi, 
                                    const DBL_VECTOR &p_MassCutoffs, 
                                    const DBL_VECTOR &p_AnCoefficients, 
                                    const DBL_VECTOR &p_BnCoefficients, 
                                          DBL_VECTOR &p_GBParams) {

#define gbParams(x) p_GBParams[static_cast<int>(GBP::x)]    // for convenience and readability - undefined at end of function

    GiantBranch::CalculateGBParams_Static(p_Mass, p_LogMetallicityXi, p_MassCutoffs, p_AnCoefficients, p_BnCoefficients, p_GBParams);                                     // calculate common values (actually, all)

    // recalculate HeHG specific values

	gbParams(B)      = CalculateCoreMass_Luminosity_B_Static();
	gbParams(D)      = CalculateCoreMass_Luminosity_D_Static(p_Mass);

    gbParams(Mx)     = GiantBranch::CalculateCoreMass_Luminosity_Mx_Static(p_GBParams);      // depends on B, D, p & q - recalculate if any of those are changed
    gbParams(Lx)     = GiantBranch::CalculateCoreMass_Luminosity_Lx_Static(p_GBParams);      // JR: Added this - depends on B, D, p, q & Mx - recalculate if any of those are changed

	gbParams(McBAGB) = p_Mass0;
	gbParams(McBGB)  = GiantBranch::CalculateCoreMassAtBGB_Static(p_Mass, p_MassCutoffs, p_AnCoefficients, p_GBParams);

    gbParams(McSN)   = CalculateCoreMassAtSupernova_Static(gbParams(McBAGB));               // JR: Added this

#undef gbParams
}


/*
 * Calculate luminosity for a Helium HertzSprung Gap star
 *
 * Uses Helium Giant Branch luminosity
 *
 *
 * double CalculateLuminosityOnPhase_Static()
 *
 * @return                                      Luminosity for a Helium HertzSprung Gap star
 */
double HeHG::CalculateLuminosityOnPhase() {
#define gbParams(x) m_GBParams[static_cast<int>(GBP::x)]    // for convenience and readability - undefined at end of function

    return HeGB::CalculateLuminosityOnPhase_Static(m_COCoreMass, gbParams(B), gbParams(D));

#undef gbParams
}


/*
 * Calculate radius of a Helium HertzSprung Gap star
 *
 * Uses Helium Giant Branch radius
 *
 *
 * double CalculateRadiusOnPhase()
 *
 * @return                                      Radius of a Helium HertzSprung Gap star
 */
double HeHG::CalculateRadiusOnPhase() {

    double R1, R2;
    std::tie(R1, R2) = HeGB::CalculateRadiusOnPhase_Static(m_Mass, m_Luminosity);

    return std::min(R1, R2);
}


/*
 * Calculate the giant branch radius for a helium star and determine new stellar type
 *
 * Hurley et al. 2000, eqs 85, 86, 87 & 88
 *
 * Calls CalculateRadiusOnPhase_Static() and returns the minimum of R1 and R2.  
 * Returns stellar type to which star should evolve based on radius calculated.
 *
 *
 * std::tuple <double, STELLAR_TYPE> CalculateRadiusAndStellarTypeOnPhase(const double p_Mass, const double p_Luminosity)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @param   [IN]    p_Luminosity                Luminosity in Lsol
 * @return                                      Radius on the helium giant branch / post-HeMs
 */
std::tuple <double, STELLAR_TYPE> HeHG::CalculateRadiusAndStellarTypeOnPhase(const double p_Mass, const double p_Luminosity) {

    double       radius;
    STELLAR_TYPE stellarType = m_StellarType;

    double R1, R2;
    std::tie(R1, R2) = HeGB::CalculateRadiusOnPhase_Static(p_Mass, p_Luminosity);

    radius = std::min(R1, R2);
    
    if (utils::Compare(R1, R2) >= 0) {
        stellarType = STELLAR_TYPE::NAKED_HELIUM_STAR_GIANT_BRANCH;
    }

   return std::make_tuple(radius, stellarType);
}


/*
 * Calculate CO Core Mass for a Helium Hertzsprung Gap star
 *
 *
 * double CalculateCOCoreMassOnPhase()
 *
 * @return                                      HeHG CoCoreMass in Msol
 */
double HeHG::CalculateCOCoreMassOnPhase() {
#define timescales(x) m_Timescales[static_cast<int>(TIMESCALE::x)]  // for convenience and readability - undefined at end of function

    return HeGB::CalculateCoreMassOnPhase_Static(m_Mass0, m_Age, timescales(tHeMS), m_GBParams);

#undef timescales
}


/*
 * Calculate rejuvenation factor for stellar age based on mass lost/gained during mass transfer
 *
 * JR: Description?
 *
 * Always returns 1.0 for HeHG - the rejuvenation factor is unity for convective main sequence stars.
 * The rest of the code is here so sanity checks can be made and an error displayed if a bad prescription
 * was specified in the program options
 *
 *
 * double CalculateMassTransferRejuvenationFactor()
 *
 * @return                                      Rejuvenation factor
 */
double HeHG::CalculateMassTransferRejuvenationFactor() {

    double fRej = 1.0;

    switch (OPTIONS->MassTransferRejuvenationPrescription()) {          // which rejuvenation prescription?

        case MT_REJUVENATION_PRESCRIPTION::NONE:                        // use default Hurley et al. 2000 prescription = 1.0
        case MT_REJUVENATION_PRESCRIPTION::STARTRACK:                   // StarTrack 2008 prescription - section 5.6 of http://arxiv.org/pdf/astro-ph/0511811v3.pdf

            if (utils::Compare(m_Mass, m_MassPrev) <= 0) {              // Rejuvenation factor is unity for mass losing stars
                fRej = 1.0;
            }
            break;

        default:                                                        // unknow prescription - use default Hurley et al. 2000 prescription = 1.0
            SHOW_WARN(ERROR::UNKNOWN_MT_REJUVENATION_PRESCRIPTION);     // show warning
    }

    return fRej;
}


/*
 * Calculate the perturbation parameter mu
 *
 * Hurley et al. 2000, eqs 97 & 98
 *
 *
 * double CalculatePerturbationMu()
 *
 * @return                                      Perturbation parameter mu
 */
double HeHG::CalculatePerturbationMu() {
    double McMax = CalculateMaximumCoreMass(m_Mass);
    return std::max(5.0 * ((McMax - m_CoreMass) / McMax), 0.0);         //return non-negative value to avoid round-off issues
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                LAMBDA CALCULATIONS                                //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * Calculate the common envelope lambda parameter using the "Nanjing" prescription
 * from X.-J. Xu and X.-D. Li arXiv:1004.4957 (v1, 28Apr2010) as implemented in STARTRACK
 *
 * This implementation adapted from the STARTRACK implementation (STARTRACK courtesy Chris Belczynski)
 *
 * This function good for HeHG and HeGB stars (for Helium stars: always use Natasha's fit)
 *
 *
 * double CalculateLambdaNanjing()
 *
 * @return                                      Nanjing lambda for use in common envelope
 */
double HeHG::CalculateLambdaNanjing() {

    double rMin = 0.25;                              // minimum considered radius: Natasha       JR: todo: should this be in constants.h?
	double rMax = 120.0;                             // maximum considered radius: Natasha       JR: todo: should this be in constants.h?

	double rMinLambda = 0.3 * PPOW(rMin, -0.8);       // JR: todo: should this be in constants.h?
	double rMaxLambda = 0.3 * PPOW(rMax, -0.8);       // JR: todo: should this be in constants.h?

	return m_Radius < rMin ? rMinLambda : (m_Radius > rMax ? rMaxLambda : 0.3 * PPOW(m_Radius, -0.8));
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                    MISCELLANEOUS FUNCTIONS / CONTROL FUNCTIONS                    //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

/*
 * Determine the star's envelope type.
 *
 *
 *
 * ENVELOPE DetermineEnvelopeType()
 *
 * @return                                      ENVELOPE::{ RADIATIVE, CONVECTIVE, REMNANT }
 */
ENVELOPE HeHG::DetermineEnvelopeType() {
    
    ENVELOPE envelope = ENVELOPE::CONVECTIVE;                                                        // default envelope type  is CONVECTIVE
    
    switch (OPTIONS->EnvelopeStatePrescription()) {                                         // which envelope prescription?
            
        case ENVELOPE_STATE_PRESCRIPTION::LEGACY:
        case ENVELOPE_STATE_PRESCRIPTION::HURLEY: // Eq. (39,40) of Hurley+ (2002) and end of section 7.2 of Hurley+ (2000) describe gradual growth of convective envelope over HG, but we approximate it as already convective here
            envelope = ENVELOPE::CONVECTIVE;
            break;
            
        case ENVELOPE_STATE_PRESCRIPTION::FIXED_TEMPERATURE:
            envelope =  utils::Compare(Temperature()*TSOL, CONVECTIVE_BOUNDARY_TEMPERATURE) ? ENVELOPE::RADIATIVE : ENVELOPE::CONVECTIVE;  // Envelope is radiative if temperature exceeds fixed threshold, otherwise convective
            break;
            
        default:                                                                                    // unknown prescription - use default envelope type
            SHOW_WARN(ERROR::UNKNOWN_ENVELOPE_STATE_PRESCRIPTION, "Using Envelope = CONVECTIVE");   // show warning
    }
    
    return envelope;
}


/*
 * Determines if mass transfer produces a wet merger
 *
 * According to the mass ratio limit discussed by de Mink et al. 2013 and Claeys et al. 2014
 *
 * Assumes this star is the donor; relevant accretor details are passed as parameters
 *
 *
 * bool IsMassRatioUnstable(const double p_AccretorMass, const bool p_AccretorIsDegenerate)
 *
 * @param   [IN]    p_AccretorMass              Mass of accretor in Msol
 * @param   [IN]    p_AccretorIsDegenerate      Boolean indicating if accretor in degenerate (true = degenerate)
 * @return                                      Boolean indicating stability of mass transfer (true = unstable)
 */
bool HeHG::IsMassRatioUnstable(const double p_AccretorMass, const bool p_AccretorIsDegenerate) {

    bool result = false;                                                                                                    // default is stable

    if (OPTIONS->MassTransferCriticalMassRatioHeliumHG()) {
        result = p_AccretorIsDegenerate
                    ? (p_AccretorMass / m_Mass) < OPTIONS->MassTransferCriticalMassRatioHeliumHGDegenerateAccretor()        // degenerate accretor
                    : (p_AccretorMass / m_Mass) < OPTIONS->MassTransferCriticalMassRatioHeliumHGNonDegenerateAccretor();    // non-degenerate accretor
    }

    return result;
}

/*
 * Choose timestep for evolution
 *
 * Follows the Discussion in Hurley et al. 2000
 *
 *
 * ChooseTimestep(const double p_Time)
 *
 * @param   [IN]    p_Time                      Current age of star in Myr
 * @return                                      Suggested timestep (dt)
 */
double HeHG::ChooseTimestep(const double p_Time) {
#define timescales(x) m_Timescales[static_cast<int>(TIMESCALE::x)]  // for convenience and readability - undefined at end of function

    // Implementation of timestep recommendation from Section 8 of Hurley et al., 2000
    double dt = utils::Compare(p_Time, timescales(tx_HeGB)) > 0
                    ? 0.02 * (timescales(tinf2_HeGB) - p_Time)
                    : 0.02 * (timescales(tinf1_HeGB) - p_Time);

    return std::max(dt, NUCLEAR_MINIMUM_TIMESTEP);

#undef timescales
}


/*
 * Determine if evolution should continue on this phase, or whether evolution
 * on this phase should end (and so evolve to next phase)
 *
 *
 * bool ShouldEvolveOnPhase()
 *
 * @return                                      Boolean flag: true if evolution on this phase should continue, false if not
 */
bool HeHG::ShouldEvolveOnPhase() {
#define gbParams(x) m_GBParams[static_cast<int>(GBP::x)]    // for convenience and readability - undefined at end of function

    double McMax = CalculateMaximumCoreMass(m_Mass);
    return utils::Compare(m_COCoreMass, McMax) <= 0 || utils::Compare(McMax, gbParams(McSN)) >= 0;    // Evolve on HeHG phase if McCO <= McMax or McMax >= McSN

#undef gbParams
}


/*
 * Modify the star after it loses its envelope
 *
 * Hurley et al. 2000, section 6 just before eq 76
 *
 *
 * STELLAR_TYPE ResolveEnvelopeLoss()
 *
 * @return                                      Stellar Type to which star should evolve after losing envelope
 */
STELLAR_TYPE HeHG::ResolveEnvelopeLoss(bool p_NoCheck) {

    STELLAR_TYPE stellarType = m_StellarType;
    
    if (p_NoCheck || utils::Compare(m_COCoreMass, m_Mass) > 0) {        // Envelope lost - determine what type of star to form

        m_CoreMass  = m_COCoreMass;
        m_HeCoreMass= m_COCoreMass;
        m_Mass      = m_CoreMass;
        m_Mass0     = m_Mass;
        
        if (!(IsSupernova())) {
            stellarType = (utils::Compare(m_COCoreMass, OPTIONS->MCBUR1() ) < 0) ? STELLAR_TYPE::CARBON_OXYGEN_WHITE_DWARF : STELLAR_TYPE::OXYGEN_NEON_WHITE_DWARF;
            m_Age       = 0.0;
            m_Radius    = HeWD::CalculateRadiusOnPhase_Static(m_Mass);
        }
    }

    return stellarType;
}

/*
 * Determine if star should continue evolution as a supernova
 *
 *
 * bool IsSupernova()
 *
 * @return                                      Boolean flag: true if star has gone Supernova, false if not
 */
bool HeHG::IsSupernova() {
    if (utils::Compare(m_COCoreMass, m_Mass) == 0) {    // special case of ultra-stripped-star -- go SN immediately if over ECSN limit
        return (utils::Compare(m_Mass, MECS) > 0);
    }
        
    return (utils::Compare(m_COCoreMass, CalculateCoreMassAtSupernova_Static(m_GBParams[static_cast<int>(GBP::McBAGB)]))>= 0); // Go supernova if CO core mass large enough
}

/*
 * Assistant function for determining the supernova explosion type
 *
 *
 * double       CalculateInitialSupernovaMass()
 *
 * @return                                      double: Initial supernova supernova mass variable
 */
double  HeHG::CalculateInitialSupernovaMass() {
    if (utils::Compare(m_COCoreMass, m_Mass) == 0) {    // special case of ultra-stripped-star -- use current mass
        return std::max(m_Mass, m_GBParams[static_cast<int>(GBP::McBAGB)]);
    }
    return GiantBranch::CalculateInitialSupernovaMass();
}

/*
 * Set parameters for evolution to next phase and return Stellar Type for next phase
 *
 *
 * STELLAR_TYPE EvolveToNextPhase()
 *
 * @return                                      Stellar Type for next phase
 */
STELLAR_TYPE HeHG::EvolveToNextPhase() {
    return STELLAR_TYPE::CARBON_OXYGEN_WHITE_DWARF;
}
