#include "Remnants.h"


/*
 * Calculate:
 *
 *     (a) the maximum mass acceptance rate of this star, as the accretor, during mass transfer, and
 *     (b) the accretion efficiency parameter
 *
 *
 * The maximum acceptance rate of the accretor star during mass transfer is based on stellar type: this function
 * is for compact remnants (NS, BH).
 *
 * Mass transfer is assumed Eddington limited for BHs and NSs.  The formalism of Nomoto/Claeys is used for WDs.
 *
 * For non compact objects:
 *
 *    1) Kelvin-Helmholtz (thermal) timescale if THERMAL (thermally limited) mass transfer efficiency
 *    2) Choose a fraction of the mass rate that will be effectively accreted for FIXED fraction mass transfer (as in StarTrack)
 *
 *
 * DBL_DBL CalculateMassAcceptanceRate(const double p_DonorMassRate, const double p_AccretorMassRate)
 *
 * @param   [IN]    p_DonorMassRate             Mass transfer rate of the donor
 * @param   [IN]    p_AccretorMassRate          Thermal mass loss rate of the accretor (this star) - ignored here
 * @return                                      Tuple containing the Maximum Mass Acceptance Rate and the Accretion Efficiency Parameter
 */
DBL_DBL Remnants::CalculateMassAcceptanceRate(const double p_DonorMassRate, const double p_AccretorMassRate) {

    double acceptanceRate   = 0.0;                                                          // acceptance mass rate - default = 0.0
    double fractionAccreted = 0.0;                                                          // accretion fraction - default=0.0

    double thisMassRate = CalculateEddingtonCriticalRate(); 

    acceptanceRate   = std::min(thisMassRate, p_DonorMassRate);
    fractionAccreted = acceptanceRate / p_DonorMassRate;

    return std::make_tuple(acceptanceRate, fractionAccreted);
}


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
double Remnants::ChooseTimestep(const double p_Time) const {

    double dtk = std::min(std::max(0.1, 10.0 * std::max(0.1, 10.0 * p_Time)), 5.0E2);
    double dte = dtk;

    return std::max(dte, NUCLEAR_MINIMUM_TIMESTEP);
}


/*
 * Calculate the critical mass ratio for unstable mass transfer
 *
 * double Remnants::CalculateCriticalMassRatio(const bool p_AccretorIsDegenerate) 
 *
 * @param   [IN]    p_AccretorIsDegenerate      Whether or not the accretor is a degenerate star
 * @return                                      Critical mass ratio
 */
double Remnants::CalculateCriticalMassRatio(const bool p_AccretorIsDegenerate) {                                           
    
        double qCrit = 0.0;
        QCRIT_PRESCRIPTION qCritPrescription = OPTIONS->QCritPrescription();

        switch (qCritPrescription) {
            case QCRIT_PRESCRIPTION::GE20: 
            case QCRIT_PRESCRIPTION::GE20_IC:
            case QCRIT_PRESCRIPTION::CLAEYS:
                break;
            case QCRIT_PRESCRIPTION::HURLEY_HJELLMING_WEBBINK:
                qCrit = 1.6; 
                break;
            default:
                m_Error = ERROR::UNKNOWN_QCRIT_PRESCRIPTION;                                     // set error value
                SHOW_ERROR(m_Error);                                                             // warn that an error occurred
        }
        return qCrit;
}