#include "BH.h"

/*
 * Calculate the gravitational mass of the Black Hole in Msol
 *
 *
 * double CalculateNeutrinoMassLoss_Static(const double p_BaryonicMass)
 *
 * @param   [IN]    p_BaryonicMass              Baryonic remnant mass in Msol
 * @return                                      Gravitational mass of remnant in Msol
 */
double BH::CalculateNeutrinoMassLoss_Static(const double p_BaryonicMass) {

    double gravitationalMass = 0.0;

    switch (OPTIONS->NeutrinoMassLossAssumptionBH()) {                                      // which assumption?

        case NEUTRINO_MASS_LOSS_PRESCRIPTION::FIXED_FRACTION:                               // FIXED FRACTION
            gravitationalMass = p_BaryonicMass * (1.0 - OPTIONS->NeutrinoMassLossValueBH());
            break;

        case NEUTRINO_MASS_LOSS_PRESCRIPTION::FIXED_MASS:                                   // FIXED MASS
            gravitationalMass = p_BaryonicMass - OPTIONS->NeutrinoMassLossValueBH();
            break;

        default:                                                                            // unknown assumption
            SHOW_WARN_STATIC(ERROR::UNKNOWN_NEUTRINO_MASS_LOSS_PRESCRIPTION,                // show warning
                             "Using gravitational mass = baryonic mass = " + std::to_string(p_BaryonicMass),
                             OBJECT_TYPE::BASE_STAR,
                             STELLAR_TYPE::BLACK_HOLE);
            gravitationalMass = p_BaryonicMass;
    }

    return gravitationalMass;
}


/*
 * Calculate core collapse Supernova parameters
 *
 *
 * DBL_DBL_DBL CalculateCoreCollapseSNParams_Static(const double p_Mass)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @return                                      Tuple containing Luminosity, Radius and Temperature of Black Hole
 */
DBL_DBL_DBL BH::CalculateCoreCollapseSNParams_Static(const double p_Mass) {
    double luminosity  = BH::CalculateLuminosityOnPhase_Static();                            // Luminosity of BH
    double radius      = BH::CalculateRadiusOnPhase_Static(p_Mass);                          // Schwarzschild radius (not correct for rotating BH)
    double temperature = BaseStar::CalculateTemperatureOnPhase_Static(luminosity, radius);   // Temperature of BH

    return std::make_tuple(luminosity, radius, temperature);
}


/*
 * Calculate the kick given to a black hole based on users chosen assumptions about black hole kicks,
 * fallback and the magnitude of a kick drawn from a distribution
 *
 * Current options are:
 *
 *    FULL    : Black holes receive the same kicks as neutron stars
 *    REDUCED : Black holes receive the same momentum kick as a neutron star, but downweighted by the black hole mass
 *    ZERO    : Black holes receive zero natal kick
 *    FALLBACK: Black holes receive a kick downweighted by the amount of mass falling back onto them
 *
 *
 *  double ReweightSupernovaKickByMass(const double p_vK, const double p_FallbackFraction, const double p_BlackHoleMass)
 *
 * @param   [IN]    p_vK                        Kick magnitude that would otherwise be applied to a neutron star
 * @param   [IN]    p_FallbackFraction          Fraction of mass that falls back onto the proto-compact object
 * @param   [IN]    p_BlackHoleMass             Mass of remnant (in Msol)
 * @return                                      Kick magnitude
 */
 double BH::ReweightSupernovaKickByMass(const double p_vK, const double p_FallbackFraction, const double p_BlackHoleMass) {

    double vK;

    switch (OPTIONS->BlackHoleKicks()) {                            // which BH kicks option specified?

        case BLACK_HOLE_KICKS::FULL:                                // BH receives full kick - no adjustment necessary
            vK = p_vK;
            break;

        case BLACK_HOLE_KICKS::REDUCED:                             // Kick is reduced by the ratio of the black hole mass to neutron star mass i.e. v_bh = ns/bh  *v_ns
            vK = p_vK * NEUTRON_STAR_MASS / p_BlackHoleMass;
            break;

        case BLACK_HOLE_KICKS::ZERO:
            vK = 0.0;                                               // BH Kicks are set to zero regardless of BH mass or kick magnitude drawn.
            break;

        case BLACK_HOLE_KICKS::FALLBACK:                            // Using the so-called 'fallback' prescription for BH kicks
            vK = p_vK * (1.0 - p_FallbackFraction);
            break;

        default:                                                    // unknown BH kick option - shouldn't happen
            vK = p_vK;                                              // return vK unchanged
            m_Error = ERROR::UNKNOWN_BH_KICK_OPTION;                // set error value
            SHOW_WARN(m_Error);                                     // warn that an error occurred
    }

    return vK;
}
