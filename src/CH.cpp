#include "CH.h"

/*
 * Calculate the log of the ratio of the lifetimes of a CH star and a normal star of the same mass (log(t_CHE/t_MS))
 * 
 * Polynomial fit derived using IZw18 and IZw18CHE BOoST models from Szecsi et al. 2020 (https://arxiv.org/abs/2004.08203)
 *
 * double CalculateLogLifetimeRatio(p_Mass)
 * 
 * @param   [IN]    p_Mass                      Mass in Msol
 *
 * @return                                      Log of the ratio of the lifetimes (log(t_CHE/t_MS)) of a CHE star to a normal star of the same mass
 */
double CH::CalculateLogLifetimeRatio(const double p_Mass) const {

    // Define variables
    double x = log10(p_Mass);

    // Coefficients obtained from np.polynomial.Polynomial.fit
    double term1 = -0.15929168474199387;
    double term2 = 1.050069750483549 * x;
    double term3 = -0.8233601359988406 * x * x;
    double term4 = 0.17772610259473764 * x * x * x;

    double logLifetimeRatio = term1 + term2 + term3 + term4;

    // Make sure we don't decrease the lifetime
    //if (logLifetimeRatio)

    return logLifetimeRatio;
}

/*
 * Calculate the ratio of the MS lifetime of a CH star and a normal star of the same mass/metallicity etc
 *
 * Polynomial fit derived using IZw18 and IZw18CHE BOoST models from Szecsi et al. 2020 (https://arxiv.org/abs/2004.08203)
 *
 * double CalculateLifetimeRatio(const double p_Mass)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 *
 * @return                                      Ratio of MS lifetimes of a CHE star to a normal star of the same mass
 */
double CH::CalculateLifetimeRatio(const double p_Mass) const {

    double lifetimeRatio = 1.0;

    // If user wants to increase CH MS lifetime, calculate the ratio of CH to MS lifetimes
    if (OPTIONS->EnhanceCHELifetimesLuminosities()) {
        double logLifetimeRatio = CalculateLogLifetimeRatio(p_Mass);
        lifetimeRatio = PPOW(10.0, logLifetimeRatio);
    }

    return lifetimeRatio;
}


/*
 * Calculate the ratio of the log luminosities (logL_CH / logL_MS) of a CH star and a normal star of the same mass
 *
 * Polynomial fit derived using IZw18 and IZw18CHE BOoST models from Szecsi et al. 2020 (https://arxiv.org/abs/2004.08203)
 * 
 * double CalculateLogLuminosityRatio(const double p_Mass)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 *
 * @return                                      Ratio of log luminosities of a CHE star to a normal star of the same mass
 */
double CH::CalculateLogLuminosityRatio(const double p_Mass) const {
    
    // Define some variables
    double logLuminosityRatio = 0.0;

    // Polynomial fits to the ratio of logL in models from Szecsi et al., derived using np.polynomial.Polynomial
    double x = log10(p_Mass);

    double term1 =  1.8261540986808193;
    double term2 = -1.1636822407341694 * x;
    double term3 =  0.5876329884434304 * x * x;
    double term4 = -0.10236336828026288 * x * x * x;

    logLuminosityRatio = term1 + term2 + term3 + term4;

    // Make sure we don't reduce the luminosity
    if (logLuminosityRatio < 1.0) {
        logLuminosityRatio = 1.0;
    }

    return logLuminosityRatio;
}

/*
 * Calculate the ratio of the luminosity of a CH star to a normal star of the same mass
 * 
 * double CalculateCHLuminosityRatio()
 * 
 * @param   [IN]    p_Mass                      Mass in Msol
 * 
 * @return                                      Ratio of luminosity of a CHE star to a normal star of the same mass
 */
double CH::CalculateCHLuminosityRatio(const double p_Mass) const {
    
    // Define some variables
    double logLuminosityRatio = 0.0;
    double luminosityRatio    = 1.0;

    // If user wants to increase CH luminosity, then calculate it
    if (OPTIONS->EnhanceCHELifetimesLuminosities()) {        
        logLuminosityRatio = CalculateLogLuminosityRatio(p_Mass); // log(L_CH) / log(L_MS) = x -> log(L_CH) = x * log(L_MS) -> L_CH = L_MS ** x

        std::cout << "logLuminosityRatio = " << logLuminosityRatio << std::endl;

        // luminosityRatio    = PPOW(10.0, logLuminosityRatio);

        // std::cout << "luminosityRatio = " << luminosityRatio << std::endl;
    }

    return luminosityRatio;
}

/*
 * Calculate the luminosity of a CH star on the main sequence
 * 
 * Start with the luminosity of a 'normal' (non-CH) star of the same mass/metallicity etc.
 * Then multiply that luminosity by the ratio of luminosities 
 * 
 * double CalculateLuminosityOnPhase(const double p_Time, const double p_Mass, const double p_LZAMS)
 *
 * @param   [IN]    p_Time                      Time (after ZAMS) in Myr
 * @param   [IN]    p_Mass                      Mass in Msol
 * @param   [IN]    p_LZAMS                     Zero Age Main Sequence (ZAMS) Luminosity
 * @return                                      Luminosity of a CH star as a function of time
 */
double CH::CalculateLuminosityOnPhase(const double p_Time, const double p_Mass, const double p_LZAMS) const {

    // std::cout << "CH::CalculateLuminosityOnPhase" << std::endl;

    // if (OPTIONS->EnhanceCHELifetimesLuminosities()) {     
    //     std::cout << "EnhanceCHELifetimesLuminosities" << std::endl;
    //     std::cout << "CH::CalculateLuminosityOnPhase" << std::endl;
    // }   

    // Firstly, calculate the standard main sequence luminosity 
    double MS_Luminosity = MainSequence::CalculateLuminosityOnPhase(p_Time, p_Mass, p_LZAMS);
    //std::cout << "MS_Luminosity = " << MS_Luminosity << std::endl;

    // Then, calculate the ratio of the standard luminosity to the CHE luminosity [log(L_CH) / log(L_MS)]
    double logLuminosityRatio = CalculateLogLuminosityRatio(p_Mass);
    //std::cout << "logLuminosityRatio = " << logLuminosityRatio << std::endl;

    // logLuminosityRatio is log(L_CH) / log(L_MS). Hence, log(L_CH) = log(L_MS) * logLuminosityRatio;
    double CH_logLuminosity = log10(MS_Luminosity) * logLuminosityRatio;
    //std::cout << "CH_logLuminosity = " << CH_logLuminosity << std::endl;

    // Exponentiate the logLuminosity to get the luminosity L_CH = 10**log(L_CH)
    double CH_Luminosity = PPOW(10.0, CH_logLuminosity);
    // std::cout << "CH_Luminosity" << std::endl;

    // double luminosityRatio = CalculateCHLuminosityRatio(p_Mass);
    // std::cout << "luminosityRatio = " << luminosityRatio << std::endl;

    // Finally, multiply the luminosity by the ratio above to get the CH star luminosity
    //double CH_Luminosity = MS_Luminosity * luminosityRatio;

    std::cout << "CH_Luminosity = " << CH_Luminosity << std::endl;

    return CH_Luminosity;
}

/*
 * Calculate luminosity at the end of the CH main sequence
 *
 * Start with the luminosity of a 'normal' (non-CH) star of the same mass/metallicity etc.
 * Then multiply that luminosity by the ratio of luminosities of a CH star to a normal star
 *
 * double CalculateLuminosityAtPhaseEnd(const double p_Mass)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @return                                      Luminosity at the end of the Main Sequence in Lsol
 */
double CH::CalculateLuminosityAtPhaseEnd(const double p_Mass) const {
    // Firstly, calculate the luminosity at the end of the main sequence for a normal star
    double TAMS_Luminosity = MainSequence::CalculateLuminosityAtPhaseEnd(p_Mass);

    // Then, calculate the ratio of the standard luminosity to the CHE luminosity
    // double luminosityRatio = CalculateCHLuminosityRatio(p_Mass);

    // Then, calculate the ratio of the standard luminosity to the CHE luminosity [log(L_CH) / log(L_MS)]
    double logLuminosityRatio = CalculateLogLuminosityRatio(p_Mass);

    // logLuminosityRatio is log(L_CH) / log(L_MS). Hence, log(L_CH) = log(L_MS) * logLuminosityRatio;
    double CH_logLuminosity = log10(TAMS_Luminosity) * logLuminosityRatio;
    //std::cout << "CH_logLuminosity = " << CH_logLuminosity << std::endl;

    // Exponentiate the logLuminosity to get the luminosity L_CH = 10**log(L_CH)
    double CH_Luminosity = PPOW(10.0, CH_logLuminosity);
    // std::cout << "CH_Luminosity" << std::endl;

    // Finally, multiply the luminosity by the ratio above to get the CH star luminosity
    return CH_Luminosity;
}

// /*
//  * Calculate lifetime of Chemically Homogeneous (CH) star on the Main Sequence
//  * 
//  * The main sequence lifetime is determined Hurley et al. 2000, eq 5
//  * The ratio of the CH star lifetime to a MS star lifetime is determined
//  * using fits to models from Szecsi et al. 2020
//  * 
//  * double CalculateLifetimeOnPhase(const double p_Mass, const double p_TBGB)
//  *
//  * @param   [IN]    p_Mass                      Mass in Msol
//  * @param   [IN]    p_TBGB                      Lifetime to Base of Giant Branch
//  * @return                                      Lifetime of Main Sequence in Myr
//  */
// double CH::CalculateLifetimeOnPhase(const double p_Mass, const double p_TBGB) const {

//     // Firstly, calculate the lifetime of a 'normal' main sequence star
//     double MS_lifetime = MainSequence::CalculateLifetimeOnPhase(p_Mass, p_TBGB);

//     // Now, calculate the ratio of the lifetime of a CH star to a normal MS star
//     double lifetimeRatio = CalculateLifetimeRatio(p_Mass);

//     // Finally, use this to determine the CH star lifetime
//     double CH_lifetime = MS_lifetime * lifetimeRatio;

//     std::cout << "CalculateLifetimeOnPhase" << std::endl;
//     std::cout << "lifetimeRatio = " << lifetimeRatio << std::endl;
//     std::cout << "MS_lifetime = " << MS_lifetime << std::endl;
//     std::cout << "CH_lifetime = " << CH_lifetime << std::endl;

//     return CH_lifetime;
// }

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
void CH::CalculateTimescales(const double p_Mass, DBL_VECTOR &p_Timescales) {
#define timescales(x) p_Timescales[static_cast<int>(TIMESCALE::x)]  // for convenience and readability - undefined at end of function
    
    std::cout << "CH::CalculateTimescales" << std::endl;
    std::cout << "m_Age, m_Tau, m_Time = " << m_Age << " " << m_Tau << " " << m_Time << std::endl;
    
    timescales(tBGB)   = CalculateLifetimeToBGB(p_Mass);

    std::cout << "tBGB = " << timescales(tBGB) << std::endl;

    timescales(tMS)    = CalculateLifetimeOnPhase(p_Mass, timescales(tBGB));

    std::cout << "OG mainSequenceLifetime = " << timescales(tMS) << std::endl;

    // Calculate the ratio of the lifetime of a CH star to a normal MS star
    double lifetimeRatio = CalculateLifetimeRatio(p_Mass);

    // Multiply the main sequence lifetime by this ratio
    timescales(tBGB) *= lifetimeRatio;
    timescales(tMS)  *= lifetimeRatio;

    std::cout << "lifetimeRatio = " << lifetimeRatio << std::endl;
    std::cout << "CH mainSequenceLifetime = " << timescales(tMS) << std::endl;

#undef timescales
}


/*
 * Evolve a CH star to the next stellar phase
 *
 * STELLAR_TYPE EvolveToNextPhase()
 * 
 * @return                                                      Stellar type to evolve to
 */
STELLAR_TYPE CH::EvolveToNextPhase() {

    STELLAR_TYPE stellarType = STELLAR_TYPE::MS_GT_07;

    if (m_Age < m_Timescales[static_cast<int>(TIMESCALE::tMS)]) {           // evolving off because of age?
        stellarType = STELLAR_TYPE::MS_GT_07;                               // no - must have spun down - evolve as MS star now
        m_CHE       = false;                                                // evolved CH->MS
    }
    else {                                                                  // yes
        stellarType = STELLAR_TYPE::NAKED_HELIUM_STAR_MS;                   // evolve as HeMS star now
        m_Age       = 0.0;                                                  // JR: can't use Hurley et al. 2000, eq 76 here - timescales(tHe) not calculated yet
        m_Tau       = 0.0;
        m_CHE       = true;                                                 // stayed on MS as CH
    }

    return stellarType;
}
