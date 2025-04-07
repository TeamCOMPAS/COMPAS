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

    if (utils::Compare(rateNJ, 0.0) > 0) m_DominantMassLossRate = MASS_LOSS_TYPE::GB;
    else                                 m_DominantMassLossRate = MASS_LOSS_TYPE::NONE;

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
double MS_gt_07::CalculateMassTransferRejuvenationFactor() {

    double fRej = 1.0;                                                                          // default - Hurley et al. 2000

    switch (OPTIONS->MassTransferRejuvenationPrescription()) {                                  // which prescription?

        case MT_REJUVENATION_PRESCRIPTION::HURLEY:                                              // use default Hurley et al. 2000 prescription = 1.0
            break;

        case MT_REJUVENATION_PRESCRIPTION::STARTRACK:                                           // StarTrack 2008 prescription - section 5.6 of http://arxiv.org/pdf/astro-ph/0511811v3.pdf
            fRej = utils::Compare(m_Mass, m_MassPrev) <= 0 ? 1.0 : m_MassPrev / m_Mass;         // rejuvenation factor is unity for mass losing stars
            break;

        default:                                                                                // unknown prescription
            // the only way this can happen is if someone added a MT_REJUVENATION_PRESCRIPTION
            // and it isn't accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a prescription this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing prescription or, if the missing
            // prescription is superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNKNOWN_MT_REJUVENATION_PRESCRIPTION);                           // throw error
    }

    return fRej;
}


/*
 * Determines if mass transfer is unstable according to the critical mass ratio.
 *
 * See e.g de Mink et al. 2013, Claeys et al. 2014, and Ge et al. 2010, 2015, 2020 for discussions.
 *
 * Assumes this star is the donor; relevant accretor details are passed as parameters.
 * Critical mass ratio is defined as qCrit = mAccretor/mDonor.
 *
 * double MS_gt_07::CalculateCriticalMassRatioClaeys14(const bool p_AccretorIsDegenerate) 
 *
 * @param   [IN]    p_AccretorIsDegenerate      Boolean indicating if accretor in degenerate (true = degenerate)
 * @return                                      Critical mass ratio for unstable MT 
 */
double MS_gt_07::CalculateCriticalMassRatioClaeys14(const bool p_AccretorIsDegenerate) const {

    double qCrit;
                                                                                                                            
    qCrit = p_AccretorIsDegenerate
                ? OPTIONS->MassTransferCriticalMassRatioMSHighMassDegenerateAccretor()      // degenerate accretor
                : OPTIONS->MassTransferCriticalMassRatioMSHighMassNonDegenerateAccretor();  // non-degenerate accretor
                                                                                                                        
    return qCrit;
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
    
    ENVELOPE envelope = ENVELOPE::RADIATIVE;                                                        // default envelope type
    
    switch (OPTIONS->EnvelopeStatePrescription()) {                                                 // which envelope prescription?
            
        case ENVELOPE_STATE_PRESCRIPTION::LEGACY:
            envelope = ENVELOPE::RADIATIVE;
            break;
            
        case ENVELOPE_STATE_PRESCRIPTION::HURLEY:
            // there is some convective envelope for stars below 1.25 solar masses according to Eq. (36) of Hurley+ (2002), but we simplify
            envelope = utils::Compare(m_Mass, 1.25) < 0 ? ENVELOPE::CONVECTIVE : ENVELOPE::RADIATIVE;
            break;
            
        case ENVELOPE_STATE_PRESCRIPTION::FIXED_TEMPERATURE:
            // envelope is radiative if temperature exceeds fixed threshold, otherwise convective
            envelope = utils::Compare(Temperature() * TSOL, OPTIONS->ConvectiveEnvelopeTemperatureThreshold()) > 0 ? ENVELOPE::RADIATIVE : ENVELOPE::CONVECTIVE;
            break;
            
        case ENVELOPE_STATE_PRESCRIPTION::CONVECTIVE_MASS_FRACTION:
            // envelope is labeled convective when the convective mass exceeds a fixed fraction of the envelope mass
            double convectiveEnvelopeMass, convectiveEnvelopeMassMax;
            std::tie(convectiveEnvelopeMass, convectiveEnvelopeMassMax) = CalculateConvectiveEnvelopeMass();
            envelope = utils::Compare(convectiveEnvelopeMass / (m_Mass - m_CoreMass), OPTIONS->ConvectiveEnvelopeMassThreshold()) > 0 ? ENVELOPE::CONVECTIVE : ENVELOPE::RADIATIVE;
            break;

        default:                                                                                    // unknown prescription
            // the only way this can happen is if someone added an ENVELOPE_STATE_PRESCRIPTION
            // and it isn't accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a prescription this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing prescription or, if the missing
            // prescription is superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNKNOWN_ENVELOPE_STATE_PRESCRIPTION);                                // throw error            
    }
    
    return envelope;
}
