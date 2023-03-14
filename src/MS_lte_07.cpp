#include "MS_lte_07.h"


/*
 * Calculate rejuvenation factor for stellar age based on mass lost/gained during mass transfer
 *
 * Description?
 *
 * Always returns 1.0 for MS_lte_07 - the rejuvenation factor is unity for convective main sequence stars.
 * The rest of the code is here so sanity checks can be made and an error displayed if a bad prescription
 * was specified in the program options
 *
 *
 * double CalculateMassTransferRejuvenationFactor()
 *
 * @return                                      Rejuvenation factor
 */
double MS_lte_07::CalculateMassTransferRejuvenationFactor() const {

    double fRej = 1.0;                                                                              // default - Hurley et al. 2000

    switch (OPTIONS->MassTransferRejuvenationPrescription()) {                                      // which prescription?

        case MT_REJUVENATION_PRESCRIPTION::NONE:                                                    // NONE: use default Hurley et al. 2000 prescription = 1.0
        case MT_REJUVENATION_PRESCRIPTION::STARTRACK:                                               // StarTrack 2008 prescription - section 5.6 of http://arxiv.org/pdf/astro-ph/0511811v3.pdf

            if (utils::Compare(m_Mass, m_MassPrev) <= 0) {                                          // Rejuvenation factor is unity for mass losing stars
                fRej = 1.0;
            }
            break;

        default:                                                                                    // unknown prescription - use default Hurley et al. 2000 prescription = 1.0
            SHOW_WARN(ERROR::UNKNOWN_MT_REJUVENATION_PRESCRIPTION, "Using default fRej = 1.0");     // show warning
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
 * double MS_lte_07::CalculateCriticalMassRatioClaeys14(const bool p_AccretorIsDegenerate) 
 *
 * @param   [IN]    p_AccretorIsDegenerate      Boolean indicating if accretor in degenerate (true = degenerate)
 * @return                                      Critical mass ratio for unstable MT 
 */
double MS_lte_07::CalculateCriticalMassRatioClaeys14(const bool p_AccretorIsDegenerate) const {

    double qCrit;
                                                                                                                            
    qCrit = p_AccretorIsDegenerate
                ? OPTIONS->MassTransferCriticalMassRatioMSLowMassDegenerateAccretor()       // degenerate accretor
                : OPTIONS->MassTransferCriticalMassRatioMSLowMassNonDegenerateAccretor();   // non-degenerate accretor
                                                                                                                        
    return qCrit;
}
