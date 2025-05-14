#include "ONeWD.h"

/* For ONeWDs, calculate:
 *
 *     (a) the maximum mass acceptance rate of this star, as the accretor, during mass transfer, and
 *     (b) the retention efficiency parameter
 *
 * This currently uses the same prescription as for COWDs, but we may consider different 
 * prescriptions in the future.
 *
 * For a given mass transfer rate, this function computes the amount of mass a WD would retain after
 * flashes, as given by appendix B of Claeys+ 2014. 
 * https://ui.adsabs.harvard.edu/abs/2014A%26A...563A..83C/abstract 
 *
 *
 * DBL_DBL CalculateMassAcceptanceRate(const double p_DonorMassRate, const bool p_IsHeRich)
 *
 * @param   [IN]    p_DonorMassRate             Mass transfer rate from the donor (Msun/Myr)
 * @param   [IN]    p_IsHeRich                  Material is He-rich or not
 * @return                                      Tuple containing the Maximum Mass Acceptance Rate (Msun/yr) and Retention Efficiency Parameter
 */
DBL_DBL ONeWD::CalculateMassAcceptanceRate(const double p_DonorMassRate, const bool p_IsHeRich) {

    m_AccretionRegime = DetermineAccretionRegime(p_DonorMassRate, p_IsHeRich); 
                                                                               
    double acceptanceRate   = 0.0;                                                       // acceptance mass rate - default = 0.0
    double fractionAccreted = 0.0;                                                       // accretion fraction - default = 0.0

    acceptanceRate = p_DonorMassRate * CalculateEtaHe(p_DonorMassRate);
    if (!p_IsHeRich) acceptanceRate *= CalculateEtaH(p_DonorMassRate);

    fractionAccreted = acceptanceRate / p_DonorMassRate;

    return std::make_tuple(acceptanceRate, fractionAccreted);
}

/*
 * Allow evolution to a new phase (currently, only SN)
 *
 * bool ShouldEvolveOnPhase()
 *
 * @return                               Whether the WD should evolve on phase or towards a SN.
 */
bool ONeWD::ShouldEvolveOnPhase() const {
    return !IsSupernova();
}


/*
 * List all conditions for SN (AIC) for ONeWD.
 * Each condition should also be a separate clause in EvolveToNextPhase.
 *
 * bool IsSupernova()
 *
 * @return                               Whether WD should undergo AIC
 */
bool ONeWD::IsSupernova() const {
    return IsMassAboveChandrasekhar();
}


/*
 * Specifies next stage, if the star changes its phase.
 *
 * STELLAR_TYPE EvolveToNextPhase()
 *
 * @return                               Stellar type of the upcoming stage.
 */
STELLAR_TYPE ONeWD::EvolveToNextPhase() {
    return ResolveAIC();
}

