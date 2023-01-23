#include "FGB.h"
#include "HeMS.h"
#include "HeWD.h"


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                              LUMINOSITY CALCULATIONS                              //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * Calculate luminosity on the First Giant Branch
 *
 * Uses Core Mass : Luminosity relation
 *
 * Hurley et al. 2000, eq 37
 *
 *
 * double CalculateLuminosityOnPhase(const double p_Time)
 *
 * @param   [IN]    p_Time                      Time in Myr
 * @return                                      Luminosity on the First Giant Branch in Lsol
 */
double FGB::CalculateLuminosityOnPhase(const double p_Time) const {
#define timescales(x) m_Timescales[static_cast<int>(TIMESCALE::x)]  // for convenience and readability - undefined at end of function
#define gbParams(x) m_GBParams[static_cast<int>(GBP::x)]            // for convenience and readability - undefined at end of function

    // The following declarations are for convenience and readability
    // (could be changed to #defines if performence is an issue - but the optimiser should optimise this away)
    double AH  = gbParams(AH);
    double B   = gbParams(B);
    double D   = gbParams(D);
    double p   = gbParams(p);
    double q   = gbParams(q);

    // Calculate the core mass according to Hurley et al. 2000, eq 39, regardless
    // of whether it is the correct expression to use given the star's mass
    double coreMass = utils::Compare(p_Time, timescales(tMx_FGB)) <= 0
                        ? PPOW(((p - 1.0) * AH * D * (timescales(tinf1_FGB) - p_Time)), (1.0 / (1.0 - p)))
                        : PPOW(((q - 1.0) * AH * B * (timescales(tinf2_FGB) - p_Time)), (1.0 / (1.0 - q)));

    return std::min((B * PPOW(coreMass, q)), (D * PPOW(coreMass, p)));

#undef gbParams
#undef timescales
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                 MASS CALCULATIONS                                 //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * Calculate core mass on the First Giant Branch
 *
 * Hurley et al. 2000, eqs 39 & 45
 *
 *
 * double CalculateCoreMassOnPhase(const double p_Mass, const double p_Time)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @param   [IN]    p_Time                      Time after ZAMS in MYRS (tBGB <= time <= tHeI)
 * @return                                      Core mass on the First Giant Branch in Msol
 */
double FGB::CalculateCoreMassOnPhase(const double p_Mass, const double p_Time) const {
#define timescales(x) m_Timescales[static_cast<int>(TIMESCALE::x)]      // for convenience and readability - undefined at end of function
#define gbParams(x) m_GBParams[static_cast<int>(GBP::x)]                // for convenience and readability - undefined at end of function
#define massCutoffs(x) m_MassCutoffs[static_cast<int>(MASS_CUTOFF::x)]  // for convenience and readability - undefined at end of function

    double McGB  = utils::Compare(p_Time, timescales(tMx_FGB)) <= 0
                    ? PPOW(((gbParams(p) - 1.0) * gbParams(AH) * gbParams(D) * (timescales(tinf1_FGB) - p_Time)), (1.0 / (1.0 - gbParams(p))))
                    : PPOW(((gbParams(q) - 1.0) * gbParams(AH) * gbParams(B) * (timescales(tinf2_FGB) - p_Time)), (1.0 / (1.0 - gbParams(q))));

    double tau   = std::max(0.0, std::min(1.0, (p_Time - timescales(tBGB)) / (timescales(tHeI) - timescales(tBGB))));

    return utils::Compare(p_Mass, massCutoffs(MHeF)) < 0 ? McGB : gbParams(McBGB) + ((CalculateCoreMassAtHeIgnition(p_Mass) - gbParams(McBGB)) * tau);

#undef massCutOffs
#undef gbParams
#undef timescales
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                            LIFETIME / AGE CALCULATIONS                            //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * Calculate relative age on the First Giant Branch
 *
 * Hurley et al. 2000, just after eq 45
 * Naturally bounded by [0, 1], but clamp here anyway
 *
 * double CalculateTauOnPhase()
 *
 * @return                                      FGB relative age, clamped to [0, 1]
 */

double FGB::CalculateTauOnPhase() const {
#define timescales(x) m_Timescales[static_cast<int>(TIMESCALE::x)]  // for convenience and readability - undefined at end of function
    return std::max(0.0, std::min(1.0, (m_Age - timescales(tBGB)) / (timescales(tHeI) - timescales(tBGB))));
#undef timescales
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                    MISCELLANEOUS FUNCTIONS / CONTROL FUNCTIONS                    //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * Choose timestep for evolution
 *
 * Can obviously do this your own way
 * Given in the discussion in Hurley et al. 2000
 *
 *
 * ChooseTimestep(const double p_Time)
 *
 * @param   [IN]    p_Time                      Current age of star in Myr
 * @return                                      Suggested timestep (dt)
 */
double FGB::ChooseTimestep(const double p_Time) const {
#define timescales(x) m_Timescales[static_cast<int>(TIMESCALE::x)]      // for convenience and readability - undefined at end of function

    double dtk = utils::Compare(p_Time, timescales(tMx_FGB)) <= 0       // ah because timescales[4,5,6] are not calculated yet   JR: todo: ?but... timescales[4] is used if this is true...? (and 5 if not)
            ? 0.02 * (timescales(tinf1_FGB) - p_Time)
            : 0.02 * (timescales(tinf2_FGB) - p_Time);

    double dte      = timescales(tHeI) - p_Time;
    double timestep = std::max(std::min(dtk, dte), NUCLEAR_MINIMUM_TIMESTEP);

    return timestep;

#undef timescales
}


/*
 * Modify the star after it loses its envelope
 *
 * Hurley et al. 2000, section 6 just before eq 76 (see also after Eq. 105)
 *
 * Where necessary updates attributes of star (depending upon stellar type):
 *
 *     - m_StellarType
 *     - m_Timescales
 *     - m_GBParams
 *     - m_Luminosity
 *     - m_Radius
 *     - m_Mass
 *     - m_Mass0
 *     - m_CoreMass
 *     - m_HeCoreMass
 *     - m_COCoreMass
 *     - m_Age
 *
 * STELLAR_TYPE ResolveEnvelopeLoss()
 *
 * @return                                      Stellar Type to which star shoule evolve after losing envelope
 */
STELLAR_TYPE FGB::ResolveEnvelopeLoss(bool p_NoCheck) {
#define timescales(x) m_Timescales[static_cast<int>(TIMESCALE::x)]                                  // for convenience and readability - undefined at end of function
#define massCutoffs(x) m_MassCutoffs[static_cast<int>(MASS_CUTOFF::x)]                              // for convenience and readability - undefined at end of function

    STELLAR_TYPE stellarType = m_StellarType;

    if (p_NoCheck || utils::Compare(m_CoreMass, m_Mass) >= 0) {                                      // Envelope loss

        m_Mass      = std::min(m_CoreMass, m_Mass);
        m_CoreMass   = m_HeCoreMass;
        m_Mass       = m_CoreMass;
        m_COCoreMass = 0.0;
        
        if (utils::Compare(m_Mass0, massCutoffs(MHeF)) < 0) {                                       // Star evolves to Helium White Dwarf

            stellarType  = STELLAR_TYPE::HELIUM_WHITE_DWARF;

            m_Age        = 0.0;
            m_Radius     = HeWD::CalculateRadiusOnPhase_Static(m_Mass);
        }
        else {                                                                                      // Star evolves to Zero age Naked Helium Main Star

            stellarType  = STELLAR_TYPE::NAKED_HELIUM_STAR_MS;

            m_Mass0      = m_Mass;
            m_Age        = 0.0;
            m_Radius     = HeMS::CalculateRadiusAtZAMS_Static(m_Mass);
        }
    }

    return stellarType;

#undef massCutoffs
#undef timescales
}


/*
 * Modify the star due to (possible) helium flash
 *
 *
 * void ResolveHeliumFlash()
 *
 * Deletermine if Helium Flash occurs, and if so set m_Mass0 equal to current mass as described in Hurley+ (2000), last paragraph before start of 7.1.1.
 */
void FGB::ResolveHeliumFlash() {
#define timescales(x) m_Timescales[static_cast<int>(TIMESCALE::x)]      // for convenience and readability - undefined at end of function
#define massCutoffs(x) m_MassCutoffs[static_cast<int>(MASS_CUTOFF::x)]  // for convenience and readability - undefined at end of function

    if (utils::Compare(m_Mass0, massCutoffs(MHeF)) < 0) {               // Helium flash if initial mass < Helium Flash cutoff
        m_Mass0 = m_Mass;                                               // for LM star at ZAHB (end of GB/begining of CHeB) due to helium flash when doing mass loss
    }

#undef massCutoffs
#undef timescales
}


/*
 * Set parameters for evolution to next phase and return Stellar Type for next phase
 *
 *
 * STELLAR_TYPE EvolveToNextPhase()
 *
 * @return                                      Stellar Type for next phase
 */
STELLAR_TYPE FGB::EvolveToNextPhase() {
    ResolveHeliumFlash();   // ...before we evolve

    return STELLAR_TYPE::CORE_HELIUM_BURNING;
}
