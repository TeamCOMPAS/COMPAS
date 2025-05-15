#include "COWD.h"

/* For COWDs, calculate:
 *
 *     (a) the maximum mass acceptance rate of this star, as the accretor, during mass transfer, and
 *     (b) the retention efficiency parameter
 *
 *
 * For a given mass transfer rate, this function computes the amount of mass a WD would retain after
 * flashes, as given by appendix B of Claeys+ 2014. 
 * https://ui.adsabs.harvard.edu/abs/2014A%26A...563A..83C/abstract 
 *
 *
 * DBL_DBL CalculateMassAcceptanceRate(const double p_DonorMassRate, const bool p_IsHeRich)
 *
 * @param   [IN]    p_DonorMassRate             Mass transfer rate from the donor
 * @param   [IN]    p_IsHeRich                  Material is He-rich or not
 * @return                                      Tuple containing the Maximum Mass Acceptance Rate (Msun/yr) and Retention Efficiency Parameter
 */
DBL_DBL COWD::CalculateMassAcceptanceRate(const double p_DonorMassRate, const bool p_IsHeRich) {

    m_AccretionRegime = DetermineAccretionRegime(p_DonorMassRate, p_IsHeRich); 
                                                                               
    double acceptanceRate   = 0.0;                                                       // acceptance mass rate - default = 0.0
    double fractionAccreted = 0.0;                                                       // accretion fraction - default = 0.0

    acceptanceRate = p_DonorMassRate * CalculateEtaHe(p_DonorMassRate);
    if (!p_IsHeRich) acceptanceRate *= CalculateEtaH(p_DonorMassRate);

    fractionAccreted = acceptanceRate / p_DonorMassRate;

    return std::make_tuple(acceptanceRate, fractionAccreted);
}


/*
 * Specifies next stage, if the star changes its phase.
 *
 * STELLAR_TYPE EvolveToNextPhase()
 *
 * @return                               Stellar type of the upcoming stage.
 */

STELLAR_TYPE COWD::EvolveToNextPhase() {

    STELLAR_TYPE stellarType;

    if (m_OffCenterIgnition) {
        stellarType = STELLAR_TYPE::OXYGEN_NEON_WHITE_DWARF;
    }
    else {                                         
        stellarType = ResolveSNIa(); 
    }
    return stellarType;
}
