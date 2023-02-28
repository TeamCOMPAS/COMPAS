#include "HeMS.h"
#include "HeWD.h"


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//             PARAMETERS, MISCELLANEOUS CALCULATIONS AND FUNCTIONS ETC.             //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


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
void HeMS::CalculateTimescales(const double p_Mass, DBL_VECTOR &p_Timescales) {
#define timescales(x) p_Timescales[static_cast<int>(TIMESCALE::x)]  // for convenience and readability - undefined at end of function

    TPAGB::CalculateTimescales(p_Mass, p_Timescales);               // calculate common values

    timescales(tHeMS) = CalculateLifetimeOnPhase_Static(p_Mass);    // recalculate tHeMS

#undef timescales
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                              LUMINOSITY CALCULATIONS                              //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * Calculate luminosity at ZAMS for a Helium Main Sequence star
 *
 * Hurley et al. 2000, eq 77
 *
 *
 * double CalculateLuminosityAtZAMS_Static(const double p_Mass)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @return                                      Luminosity at ZAMS for a Helium Main Sequence star in Lsol
 */
double HeMS::CalculateLuminosityAtZAMS_Static(const double p_Mass) {

    // pow() is slow - use multiplication (sqrt() is much faster than pow())
    double m_0_5   = std::sqrt(p_Mass);
    double m_3     = p_Mass * p_Mass * p_Mass;
    double m_6     = m_3 * m_3;
    double m_7_5   = m_6 * p_Mass * m_0_5;
    double m_9     = m_6 * m_3;
    double m_10_25 = m_9 * p_Mass * std::sqrt(m_0_5);

    return (15262.0 * m_10_25) / (m_9 + (29.54 * m_7_5) + (31.18 * m_6) + 0.0469);
}


/*
 * Calculate luminosity at the end of the Helium Main Sequence
 *
 * Hurley et al. 2000, eq 80 at tau = 1
 *
 *
 * double CalculateLuminosityAtPhaseEnd_Static(const double p_Mass)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @return                                      Luminosity at the end of the Helium Main Sequence in Lsol
 */
double HeMS::CalculateLuminosityAtPhaseEnd_Static(const double p_Mass) {
    return CalculateLuminosityAtZAMS_Static(p_Mass) * (1.0 + 0.45 + std::max(0.0, 0.85 - (0.08 * p_Mass)));
}


/*
 * Calculate luminosity for a Helium Main Seqeuence star (during central He burning)
 *
 * Hurley et al. 2000, eqs 80 & 82
 *
 *
 * double CalculateLuminosityOnPhase_Static(const double Mass, const double p_Tau)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @param   [IN]    p_Tau                       tHeMS relative age
 * @return                                      Luminosity for a Helium Main Sequence star in Lsol
 */
double HeMS::CalculateLuminosityOnPhase_Static(const double p_Mass, const double p_Tau) {

    double alpha = std::max(0.0, 0.85 - 0.08 * p_Mass);

    return CalculateLuminosityAtZAMS_Static(p_Mass) * (1.0 + (p_Tau * (0.45 + (alpha * p_Tau))));
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                RADIUS CALCULATIONS                                //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * Calculate radius at ZAMS for a Helium Main Sequence star
 *
 * Hurley et al. 2000, eq 78
 *
 *
 * double CalculateRadiusAtZAMS_Static(const double p_Mass)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @return                                      Radius at ZAMS for a Helium Main Sequence star in Rsol
 */
double HeMS::CalculateRadiusAtZAMS_Static(const double p_Mass) {
    // pow() is slow - use multiplication
    double m_3 = p_Mass * p_Mass * p_Mass;
    double m_4 = m_3 * p_Mass;

    return (0.2391 * PPOW(p_Mass, 4.6)) / (m_4 + (0.162 * m_3) + 0.0065);
}


/*
 * Calculate radius for a Helium Main Sequence star (during central He burning)
 *
 * Hurley et al. 2000, eq 81
 *
 *
 * double CalculateRadiusOnPhase_Static(const double p_Mass, const double p_Tau)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @param   [IN]    p_Tau                       tHeMS relative age
 * @return                                      Radius for a Helium Main Sequence star in Rsol
 */
double HeMS::CalculateRadiusOnPhase_Static(const double p_Mass, const double p_Tau) {

    double tau_6 = p_Tau * p_Tau * p_Tau * p_Tau * p_Tau * p_Tau;   // pow() is slow - use multiplication
    double beta  = std::max(0.0, 0.4 - 0.22 * log10(p_Mass));

    return CalculateRadiusAtZAMS_Static(p_Mass) * (1.0 + (beta * (p_Tau - tau_6)));
}


/*
 * Calculate the radius at the end of the helium main sequence
 *
 * Hurley et al. 2000, eq 81 at tau = 1
 *
 *
 * double CalculateRadiusAtPhaseEnd_Static(const double p_Mass)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @return                                      Radius at the end of the helium main sequence (RTHe)
 */
double HeMS::CalculateRadiusAtPhaseEnd_Static(const double p_Mass) {
    return CalculateRadiusAtZAMS_Static(p_Mass);
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                 MASS CALCULATIONS                                 //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * Calculate rejuvenation factor for stellar age based on mass lost/gained during mass transfer
 *
 * Description?
 *
 *
 * double CalculateMassTransferRejuvenationFactor()
 *
 * @return                                      Rejuvenation factor
 */
double HeMS::CalculateMassTransferRejuvenationFactor() const {

    double fRej = 1.0;                                                                              // default value

    switch (OPTIONS->MassTransferRejuvenationPrescription()) {

        case MT_REJUVENATION_PRESCRIPTION::NONE:                                                    // use default Hurley et al. 2000 prescription = 1.0
            break;

        case MT_REJUVENATION_PRESCRIPTION::STARTRACK:                                               // StarTrack 2008 prescription - section 5.6 of http://arxiv.org/pdf/astro-ph/0511811v3.pdf
            fRej = utils::Compare(m_Mass, m_MassPrev) <= 0 ? 1.0 : m_MassPrev / m_Mass;             // rejuvenation factor is unity for mass losing stars
            break;

        default:                                                                                    // shouldn't get here - use default Hurley et al. 2000 prescription = 1.0
            SHOW_WARN(ERROR::UNKNOWN_MT_REJUVENATION_PRESCRIPTION, "Using default fRej = 1.0");     // show warning
    }

    return fRej;
}


/*
 * Calculate the dominant mass loss mechanism and associated rate for the star
 * at the current evolutionary phase.
 *
 * According to Hurley et al. 2000
 *
 * double CalculateMassLossRateHurley()
 *
 * @return                                      Mass loss rate in Msol per year
 */
double HeMS::CalculateMassLossRateHurley() {
    double rateNJ = CalculateMassLossRateNieuwenhuijzenDeJager();
    double rateKR = CalculateMassLossRateKudritzkiReimers();
    double rateWR = CalculateMassLossRateWolfRayet(0.0); // use mu=0.0 for Helium stars
    double dominantRate;

    if (utils::Compare(rateNJ, rateKR) > 0) {
        dominantRate = rateNJ;
        m_DominantMassLossRate = MASS_LOSS_TYPE::NIEUWENHUIJZEN_DE_JAGER;
    } else {
        dominantRate = rateKR;
        m_DominantMassLossRate = MASS_LOSS_TYPE::KUDRITZKI_REIMERS;
    }
    if (utils::Compare(rateWR, dominantRate) > 0) {
        dominantRate = rateWR;
        m_DominantMassLossRate = MASS_LOSS_TYPE::WOLF_RAYET_LIKE;
    }

    return dominantRate;
}


/*
 * Calculate the dominant mass loss mechanism and associated rate for the star at the current evolutionary phase
 *
 * According to Vink - based on implementation in StarTrack
 *
 * double CalculateMassLossRateVink()
 *
 * @return                                      Mass loss rate in Msol per year
 */
double HeMS::CalculateMassLossRateVink() {
    m_DominantMassLossRate = MASS_LOSS_TYPE::WOLF_RAYET_LIKE;
    return CalculateMassLossRateWolfRayetZDependent(0.0);
}


/*
 * Determines if mass transfer is unstable according to the critical mass ratio.
 *
 * See e.g de Mink et al. 2013, Claeys et al. 2014, and Ge et al. 2010, 2015, 2020 for discussions.
 *
 * Assumes this star is the donor; relevant accretor details are passed as parameters.
 * Critical mass ratio is defined as qCrit = mAccretor/mDonor.
 *
 * double HeMS::CalculateCriticalMassRatioClaeys14(const bool p_AccretorIsDegenerate) 
 *
 * @param   [IN]    p_AccretorIsDegenerate      Boolean indicating if accretor in degenerate (true = degenerate)
 * @return                                      Critical mass ratio for unstable MT 
 */
double HeMS::CalculateCriticalMassRatioClaeys14(const bool p_AccretorIsDegenerate) const {

    double qCrit;
                                                                                                                            
    qCrit = p_AccretorIsDegenerate
                ? OPTIONS->MassTransferCriticalMassRatioHeliumMSDegenerateAccretor()        // degenerate accretor
                : OPTIONS->MassTransferCriticalMassRatioHeliumMSNonDegenerateAccretor();    // non-degenerate accretor
                                                                                                                        
    return qCrit;
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                            LIFETIME / AGE CALCULATIONS                            //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * Calculate lifetime for a Helium Main Sequence star
 *
 * Hurley et al. 2000, eq 79
 *
 *
 * double CalculateLifetimeOnPhase_Static(const double p_Mass)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @return                                      Lifetime for a Helium Main Sequence star in Myr
 */
double HeMS::CalculateLifetimeOnPhase_Static(const double p_Mass) {

    // pow() is slow - use multiplication (sqrt() is much faster than pow())
    double m_4   = p_Mass * p_Mass * p_Mass * p_Mass;
    double m_6   = m_4 * p_Mass * p_Mass;
    double m_6_5 = m_6 * std::sqrt(p_Mass);

    return (0.4129 + (18.81 * m_4) + (1.853 * m_6)) / m_6_5;
}


/*
 * Recalculates the star's age after mass loss
 *
 * Hurley et al. 2000, section 7.1
 *
 * Modifies attribute m_Age
 *
 *
 * UpdateAgeAfterMassLoss()
 *
 */
void HeMS::UpdateAgeAfterMassLoss() {

    double tHeMS      = m_Timescales[static_cast<int>(TIMESCALE::tHeMS)];
    double tHeMSprime = CalculateLifetimeOnPhase_Static(m_Mass);

    m_Age *= tHeMSprime / tHeMS;
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
double HeMS::ChooseTimestep(const double p_Time) const {
#define timescales(x) m_Timescales[static_cast<int>(TIMESCALE::x)]  // for convenience and readability - undefined at end of function

    double dtk = 0.05 * timescales(tHeMS);
    double dte = timescales(tHeMS) - p_Time;

    return std::max(std::min(dtk, dte), NUCLEAR_MINIMUM_TIMESTEP);

#undef timescales
}


/*
 * Modify the star after it loses its envelope
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
 *
 * STELLAR_TYPE ResolveEnvelopeLoss()
 *
 * @return                                      Stellar Type to which star shoule evolve after losing envelope
 */
STELLAR_TYPE HeMS::ResolveEnvelopeLoss(bool p_NoCheck) {

    STELLAR_TYPE stellarType = m_StellarType;

    if (p_NoCheck || utils::Compare(m_Mass, 0.0) <= 0) {
        stellarType = STELLAR_TYPE::MASSLESS_REMNANT;
        m_Radius = 0.0;   // massless remnant
        m_Mass = 0.0;
    }

    return stellarType;
}


/*
 * Set parameters for evolution to next phase and return Stellar Type for next phase
 *
 *
 * STELLAR_TYPE EvolveToNextPhase()
 *
 * @return                                      Stellar Type for next phase
 */
STELLAR_TYPE HeMS::EvolveToNextPhase() {
#define timescales(x) m_Timescales[static_cast<int>(TIMESCALE::x)]  // for convenience and readability - undefined at end of function

    m_Age = timescales(tHeMS);

    return STELLAR_TYPE::NAKED_HELIUM_STAR_HERTZSPRUNG_GAP;

#undef timescales
}
