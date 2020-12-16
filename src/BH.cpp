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
                             "Using mass = 0.0",
                             OBJECT_TYPE::BASE_STAR,
                             STELLAR_TYPE::BLACK_HOLE);
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
