#include "MS_gt_07.h"


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
double MS_gt_07::CalculateMassLossRateHurley() {
    double rateNJ = CalculateMassLossRateNieuwenhuijzenDeJager();
    if (utils::Compare(rateNJ, 0.0) > 0) {
        m_DominantMassLossRate = MASS_LOSS_TYPE::NIEUWENHUIJZEN_DE_JAGER;
    } else {
        m_DominantMassLossRate = MASS_LOSS_TYPE::NONE;
    }
    return rateNJ;
}


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
double MS_gt_07::CalculateMassTransferRejuvenationFactor() const {

    double fRej = 1.0;                                                                              // default - Hurley et al. 2000

    switch (OPTIONS->MassTransferRejuvenationPrescription()) {                                      // which prescription?

        case MT_REJUVENATION_PRESCRIPTION::NONE:                                                    // NONE: use default Hurley et al. 2000 prescription = 1.0
            break;

        case MT_REJUVENATION_PRESCRIPTION::STARTRACK:                                               // StarTrack 2008 prescription - section 5.6 of http://arxiv.org/pdf/astro-ph/0511811v3.pdf
            fRej = utils::Compare(m_Mass, m_MassPrev) <= 0 ? 1.0 : m_MassPrev / m_Mass;             // rejuvenation factor is unity for mass losing stars
            break;

        default:                                                                                    // unknown prescription - use default Hurley et al. 2000 prescription = 1.0
            SHOW_WARN(ERROR::UNKNOWN_MT_REJUVENATION_PRESCRIPTION, "Using default fRej = 1.0");     // show warning
    }

    return fRej;
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
bool MS_gt_07::IsMassRatioUnstable(const double p_AccretorMass, const bool p_AccretorIsDegenerate) const {

    bool result = false;                                                                                                    // default is stable

    if (OPTIONS->MassTransferCriticalMassRatioMSHighMass()) {
        result = p_AccretorIsDegenerate
                    ? (p_AccretorMass / m_Mass) < OPTIONS->MassTransferCriticalMassRatioMSHighMassDegenerateAccretor()      // degenerate accretor
                    : (p_AccretorMass / m_Mass) < OPTIONS->MassTransferCriticalMassRatioMSHighMassNonDegenerateAccretor();  // non-degenerate accretor
    }

    return result;
}


/*
 * Determine the star's envelope type.
 *
 *
 *
 * ENVELOPE DetermineEnvelopeType()
 *
 * @return                                      ENVELOPE::{ RADIATIVE, CONVECTIVE, REMNANT }
 */
ENVELOPE MS_gt_07::DetermineEnvelopeType() const {
    
    ENVELOPE envelope = ENVELOPE::RADIATIVE;                                                        // default envelope type  is RADIATIVE
    
    switch (OPTIONS->EnvelopeStatePrescription()) {                                                 // which envelope prescription?
            
        case ENVELOPE_STATE_PRESCRIPTION::LEGACY:
            envelope = ENVELOPE::RADIATIVE;
            break;
            
        case ENVELOPE_STATE_PRESCRIPTION::HURLEY:
            // there is some convective envelope for stars below 1.25 solar masses according to Eq. (36) of Hurley+ (2002), but we simplify
            envelope =  utils::Compare(m_Mass, 1.25) < 0 ? ENVELOPE::CONVECTIVE : ENVELOPE::RADIATIVE;
            break;
            
        case ENVELOPE_STATE_PRESCRIPTION::FIXED_TEMPERATURE:
            envelope =  utils::Compare(Temperature() * TSOL, CONVECTIVE_BOUNDARY_TEMPERATURE) > 0 ? ENVELOPE::RADIATIVE : ENVELOPE::CONVECTIVE;  // Envelope is radiative if temperature exceeds fixed threshold, otherwise convective
            break;
            
        default:                                                                                    // unknown prescription - use default envelope type
            SHOW_WARN(ERROR::UNKNOWN_ENVELOPE_STATE_PRESCRIPTION, "Using Envelope = RADIATIVE");    // show warning
    }
    
    return envelope;
}
