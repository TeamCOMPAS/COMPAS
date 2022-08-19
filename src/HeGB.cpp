#include "HeGB.h"


/*
 * Calculate luminosity on the Helium Giant Branch
 *
 * Hurley et al. 2000, eq 84
 *
 *
 * double CalculateLuminosityOnPhase_Static(const double p_CoreMass, const double p_GBPB, const double p_GBPD)
 *
 * @param   [IN]    p_CoreMass                  Core Mass in Msol
 * @param   [IN]    p_GBPB                      Giant Branch parameter B (Hurley et al. 2000, eq 38)
 * @param   [IN]    p_GBPD                      Giant Branch parameter D (Hurley et al. 2000, eq 38)
 * @return                                      Helium Giant Branch (post-HeMS) luminosity in Lsol
 *
 * p_GBPB and p_GBPD passed as parameters so function can be declared static
 */
double HeGB::CalculateLuminosityOnPhase_Static(const double p_CoreMass, const double p_GBPB, const double p_GBPD) {
    double Mc_3 = p_CoreMass * p_CoreMass * p_CoreMass;
    double Mc_5 = Mc_3 * p_CoreMass * p_CoreMass;

    return std::min((p_GBPB * Mc_3), (p_GBPD * Mc_5));
}


/*
 * Calculate the giant branch radius for a helium star
 *
 * Hurley et al. 2000, eqs 85, 86, 87 & 88
 *
 * Calculates and returns R1 and R2 - the caller can then choose the radius and the
 * resultant stellar type based on the radius chosen
 *
 *
 * std::tuple<double, double> CalculateRadiusOnPhase_Static(const double p_Mass, const double p_Luminosity)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @param   [IN]    p_Luminosity                Luminosity in Lsol
 * @return                                      Tuple containing R1 & R2: Radius on the helium giant branch / post-HeMs
 */
std::tuple<double, double> HeGB::CalculateRadiusOnPhase_Static(const double p_Mass, const double p_Luminosity) {

    double RZHe = HeMS::CalculateRadiusAtZAMS_Static(p_Mass);
    double LTHe = CalculateLuminosityAtPhaseEnd_Static(p_Mass);

    // pow() is slow - use multiplication (srqt() is much faster than pow())
    double m_2   = p_Mass * p_Mass;
    double m_2_5 = m_2 * std::sqrt(p_Mass);
    double m_5   = m_2_5 * m_2_5;

    double lamda = 500.0 * (2.0 + m_5) / m_2_5;

    double R1 = RZHe * PPOW((p_Luminosity / LTHe), 0.2) + (0.02 * (exp(p_Luminosity / lamda) - exp(LTHe / lamda)));
    double R2 = 0.08 * PPOW(p_Luminosity, 0.75);     // Mimics the Hayashi track

    // May need to change this according to section 2.1.2 of http://iopscience.iop.org/article/10.1086/340304/pdf which talks about updated helium star evolution
    // or replace with helium star tracks from binary_c
    // Especially important for low mass helium stars, BNS progenitors
    // This paper suggest mass below which envelope is convective is Mconv = 4.5 Msol, they leave it as an uncertain model parameter-- here it is set in constants (as usual, may want to edit from commandLine)
    // Rapid expansion given by the second term in R1

    return std::make_tuple(R1, R2);
}


/*
 * Calculate the giant branch radius for a helium star
 *
 * Hurley et al. 2000, eqs 85, 86, 87 & 88
 *
 * Calls CalculateRadiusOnPhase_Static() and returns the minimum of R1 and R2.
 *
 *
 * double CalculateRadiusOnPhase(const double p_Mass, const double p_Luminosity)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @param   [IN]    p_Luminosity                Luminosity in Lsol
 * @return                                      Radius on the helium giant branch / post-HeMs
 */
double HeGB::CalculateRadiusOnPhase(const double p_Mass, const double p_Luminosity) const {
    double R1, R2;
    std::tie(R1, R2) = CalculateRadiusOnPhase_Static(p_Mass, p_Luminosity);

    return std::min(R1, R2);
}


/*
 * Calculate the giant branch radius for a helium star and determine new stellar type
 *
 * Hurley et al. 2000, eqs 85, 86, 87 & 88
 *
 * Calls CalculateRadiusOnPhase_Static() and returns the minimum of R1 and R2.  
 * Returns stellr type to which star should evolve based on radius calculated.
 *
 *
 * std::tuple <double, STELLAR_TYPE> CalculateRadiusAndStellarTypeOnPhase(const double p_Mass, const double p_Luminosity)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @param   [IN]    p_Luminosity                Luminosity in Lsol
 * @return                                      Radius on the helium giant branch / post-HeMs
 */
std::tuple <double, STELLAR_TYPE> HeGB::CalculateRadiusAndStellarTypeOnPhase(const double p_Mass, const double p_Luminosity) const {

    double       radius;
    STELLAR_TYPE stellarType = m_StellarType;

    double R1, R2;
    std::tie(R1, R2) = CalculateRadiusOnPhase_Static(p_Mass, p_Luminosity);

    radius = std::min(R1, R2);

    if (utils::Compare(R1, R2) < 0) {
        stellarType = STELLAR_TYPE::NAKED_HELIUM_STAR_HERTZSPRUNG_GAP;
    }

   return std::make_tuple(radius, stellarType);
}


/*
 * Calculate the age of an evolved Helium Giant branch star
 *
 * Hurley et al. 2000, eq 39 (core mass - t relation, modified as described after eq 84)
 *
 *
 * double CalculateAgeOnPhase_Static(const double      p_Mass,
 *                                   const double      p_CoreMass,
 *                                   const double      p_tHeMS,
 *                                   const DBL_VECTOR &p_GBParams)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @param   [IN]    p_CoreMass                  Core mass in Msol
 * @param   [IN]    p_tHeMS                     Naked Helium Star central helium burning lifetime
 * @param   [IN]    p_GBParams                  Giant Branch parameters
 * @return                                      Age in Myr
 *
 * p_tHeMS and p_GBParams passed as parameters so function can be declared static
 */
double HeGB::CalculateAgeOnPhase_Static(const double      p_Mass,
                                        const double      p_CoreMass,
                                        const double      p_tHeMS,
                                        const DBL_VECTOR &p_GBParams) {
#define gbParams(x) p_GBParams[static_cast<int>(GBP::x)]    // for convenience and readability - undefined at end of function

    double age;

    double p1    = gbParams(p) - 1.0;
    double p1_p  = p1 / gbParams(p);

    double LtHe  = HeMS::CalculateLuminosityAtPhaseEnd_Static(p_Mass);
    double tinf1 = p_tHeMS + ((1.0 / (p1 * gbParams(AHe) * gbParams(D))) * PPOW(gbParams(D) / LtHe, p1_p));

    if (utils::Compare(p_CoreMass, gbParams(Mx)) > 0) {
        double q1    = gbParams(q) - 1.0;
        double tx    = tinf1 - ((tinf1 - p_tHeMS) * PPOW(LtHe / gbParams(Lx), p1_p));
        double tinf2 = tx + ((1.0 / (q1 * gbParams(AHe) * gbParams(B))) * PPOW(gbParams(B) / gbParams(Lx), q1 / gbParams(q)));

        age = tinf2 - (PPOW(p_CoreMass, 1.0 - gbParams(q)) / (q1 * gbParams(AHe) * gbParams(B)));
    }
    else {
        age = tinf1 - (PPOW(p_CoreMass, 1.0 - gbParams(p)) / (p1 * gbParams(AHe) * gbParams(D)));
    }

    return age;

    // still confused as to why t <= tx translates to Mc >= Mx ?

#undef gbParams
}


/*
 * Calculate the core mass on the Helium Giant Branch
 *
 * Hurley et al. 2000, eq 39 (core mass - t relation, modified as described after eq 84)
 *
 *
 * double CalculateCoreMassOnPhase_Static(const double      p_Mass,
 *                                        const double      p_Time,
 *                                        const double      p_tHeMS,
 *                                        const DBL_VECTOR &p_GBParams)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @param   [IN]    p_Time                      Time in Myr
 * @param   [IN]    p_tHeMS                     Naked Helium Star central helium burning lifetime
 * @param   [IN]    p_GBParams                  Giant Branch parameters
 * @return                                      Age in Myr
 *
 * p_tHeMS and p_GBParams passed as parameters so function can be declared static
 */
double HeGB::CalculateCoreMassOnPhase_Static(const double      p_Mass,
                                             const double      p_Time,
                                             const double      p_tHeMS,
                                             const DBL_VECTOR &p_GBParams) {
#define gbParams(x) p_GBParams[static_cast<int>(GBP::x)]// for convenience and readability - undefined at end of function

    double coreMass;

    double p1   = gbParams(p) - 1.0;
    double p1_p = p1 / gbParams(p);

    double LtHe  = HeMS::CalculateLuminosityAtPhaseEnd_Static(p_Mass);
    double tinf1 = p_tHeMS + ((1.0 / (p1 * gbParams(AHe) * gbParams(D))) * PPOW(gbParams(D) / LtHe, p1_p));
    double tx    = tinf1 - (tinf1 - p_tHeMS) * PPOW((LtHe / gbParams(Lx)), p1_p);
    
    if (utils::Compare(p_Time, tx) > 0) {
        double q1    = gbParams(q) - 1.0;
        double tinf2 = tx + ((1.0 / (q1 * gbParams(AHe) * gbParams(B))) * PPOW(gbParams(B) / gbParams(Lx), q1 / gbParams(q)));
        coreMass = PPOW(q1 * gbParams(AHe) * gbParams(B) * (tinf2 - p_Time), 1.0 / (1.0 - gbParams(q)));
    }
    else {
        coreMass = PPOW(p1 * gbParams(AHe) * gbParams(D) * (tinf1 - p_Time), 1.0 / (1.0 - gbParams(p)));
    }

    return coreMass;

#undef gbParams
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
bool HeGB::IsMassRatioUnstable(const double p_AccretorMass, const bool p_AccretorIsDegenerate) const {

    bool result = false;                                                                                                    // default is stable

    result = p_AccretorIsDegenerate
                ? (p_AccretorMass / m_Mass) < OPTIONS->MassTransferCriticalMassRatioHeliumGiantDegenerateAccretor()     // degenerate accretor
                : (p_AccretorMass / m_Mass) < OPTIONS->MassTransferCriticalMassRatioHeliumGiantNonDegenerateAccretor(); // non-degenerate accretor

    return result;
}
