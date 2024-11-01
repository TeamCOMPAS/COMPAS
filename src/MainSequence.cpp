#include "MainSequence.h"
#include "MS_gt_07.h"
#include "HG.h"


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//             PARAMETERS, MISCELLANEOUS CALCULATIONS AND FUNCTIONS ETC.             //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * Calculate the helium abundance in the core of the star
 * 
 * Currently just a simple linear model from the initial helium abundance to 
 * the maximum helium abundance (assuming that all hydrogen is converted to
 * helium). 
 * 
 * When tau = 0, heliumAbundanceCore = m_InitialHeliumAbundance
 * When tau = 1, heliumAbundanceCore = heliumAbundanceCoreMax = 1.0 - m_Metallicity
 * 
 * Should be updated to match detailed models.
 *
 * double CalculateHeliumAbundanceCoreOnPhase(const double p_Tau)
 * 
 * @param   [IN]    p_Tau                       Fraction of main sequence lifetime
 *
 * @return                                      Helium abundance in the core (Y_c)
 */
double MainSequence::CalculateHeliumAbundanceCoreOnPhase(const double p_Tau) const {
    double heliumAbundanceCoreMax = 1.0 - m_Metallicity;
    return ((heliumAbundanceCoreMax - m_InitialHeliumAbundance) * p_Tau) + m_InitialHeliumAbundance;
}


/*
 * Calculate the hydrogen abundance in the core of the star
 * 
 * Currently just a simple linear model. Assumes that hydrogen in the core of 
 * the star is burned to helium at a constant rate throughout the lifetime. 
 * 
 * Should be updated to match detailed models.
 *
 * double CalculateHydrogenAbundanceCoreOnPhase(const double p_Tau)
 * 
 * @param   [IN]    p_Tau                       Fraction of main sequence lifetime
 *
 * @return                                      Hydrogen abundance in the core (X_c)
 */
double MainSequence::CalculateHydrogenAbundanceCoreOnPhase(const double p_Tau) const {
    return m_InitialHydrogenAbundance * (1.0 - p_Tau);
}


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
void MainSequence::CalculateTimescales(const double p_Mass, DBL_VECTOR &p_Timescales) {
#define timescales(x) p_Timescales[static_cast<int>(TIMESCALE::x)]  // for convenience and readability - undefined at end of function
    timescales(tBGB) = CalculateLifetimeToBGB(p_Mass);
    timescales(tMS)  = CalculateLifetimeOnPhase(p_Mass, timescales(tBGB));
#undef timescales
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                              LUMINOSITY CALCULATIONS                              //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * Calculate the luminosity perturbation delta_L
 *
 * Hurley et al. 2000, eq 16
 *
 *
 * double CalculateDeltaL(const double p_Mass)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @return                                      Luminosity perturbation (delta_L in Hurley et al. 2000)
 */
double MainSequence::CalculateDeltaL(const double p_Mass) const {
#define a m_AnCoefficients                                              // for convenience and readability - undefined at end of function
#define massCutoffs(x) m_MassCutoffs[static_cast<int>(MASS_CUTOFF::x)]  // for convenience and readability - undefined at end of function

    double deltaL;

    if (utils::Compare(p_Mass, massCutoffs(MHook)) <= 0) {              // per Hurley et al. 2000, eq 16
        deltaL = 0.0;                                                   // this really is supposed to be zero
    }
    else if (utils::Compare(p_Mass, a[33]) < 0) {
        double top    = p_Mass - massCutoffs(MHook);
        double bottom = a[33] - massCutoffs(MHook);
        deltaL        = m_LConstants[static_cast<int>(L_CONSTANTS::B_DELTA_L)] * PPOW((top / bottom), 0.4);
    }
    else {
        deltaL = std::min((a[34] / PPOW(p_Mass, a[35])), (a[36] / PPOW(p_Mass, a[37])));
    }

    return deltaL;

#undef massCutoffs
#undef a
}


/*
 * Calculate the luminosity beta coefficient
 *
 * Hurley et al. 2000, eq 20
 *
 *
 * double CalculateBetaL(const double p_Mass)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @return                                      Luminosity beta coefficient (beta_L in Hurley et al. 2000)
 */
double MainSequence::CalculateBetaL(const double p_Mass) const {
#define a m_AnCoefficients    // for convenience and readability - undefined at end of function

    double betaL  = std::max(0.0, (a[54] - (a[55] * PPOW(p_Mass, a[56]))));
    if ((utils::Compare(p_Mass, a[57]) > 0) && (utils::Compare(betaL, 0.0) > 0)) {
        double bBetaL = m_LConstants[static_cast<int>(L_CONSTANTS::B_BETA_L)];

        betaL = std::max(0.0, (bBetaL - 10.0 * (p_Mass - a[57]) * bBetaL));
    }

    return betaL;

#undef a
}


/*
 * Calculate the luminosity alpha constant alpha_L
 *
 * Hurley et al. 2000, eqs 19a & 19b
 *
 *
 * double CalculateAlphaL(const double p_Mass)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @return                                      Luminosity alpha constant (alpha_L in Hurley et al. 2000)
 */
double MainSequence::CalculateAlphaL(const double p_Mass) const {
#define a m_AnCoefficients    // for convenience and readability - undefined at end of function

    // You might find that functions you give Metallicity to as an argument don't actually need it -- metallicity dependence is in an/MFGB etc.
    // Also, if this is likely to be called in a loop, try to precompute it (only depends on initial values of mass/metallicity right?
    // Of course the only problem with this kind of thing is it makes it less flexible if you have to change one of those)

    double alphaL  = 0.0;

         if (utils::Compare(p_Mass, 0.5)   < 0) alphaL = a[49];
    else if (utils::Compare(p_Mass, 0.7)   < 0) alphaL = a[49] + (5.0 * (0.3 - a[49]) * (p_Mass - 0.5));
    else if (utils::Compare(p_Mass, a[52]) < 0) alphaL = 0.3 + ((a[50] - 0.3) * (p_Mass - 0.7) / (a[52] - 0.7));
    else if (utils::Compare(p_Mass, a[53]) < 0) alphaL = a[50] + ((a[51] - a[50]) * (p_Mass - a[52]) / (a[53] - a[52]));
    else if (utils::Compare(p_Mass, 2.0)   < 0) alphaL = a[51] + ((m_LConstants[static_cast<int>(L_CONSTANTS::B_ALPHA_L)] - a[51]) * (p_Mass - a[53]) / (2.0 - a[53]));
    else                                        alphaL = (a[45] + (a[46] * PPOW(p_Mass, a[48]))) / (PPOW(p_Mass, 0.4) + (a[47] * PPOW(p_Mass, 1.9)));

    return alphaL;

#undef a
}


/*
 * Calculate the exponent eta (for Hurley et al. 2000, eq 12)
 *
 * Hurley et al. 2000, eq 18
 *
 *
 * double CalculateEta(const double p_Mass)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @return                                      The exponent eta (for Hurley et al. 2000, eq 12)
 */
double MainSequence::CalculateEta(const double p_Mass) const {

    double eta = 10.0;

    if (utils::Compare(m_Metallicity, 0.0009) <= 0) {
        if (utils::Compare(p_Mass, 1.1) >= 0) {
            eta = 20.0;
        }
        else if (utils::Compare(p_Mass, 1.0) > 0) {
            eta = (100.0 * p_Mass) - 90.0;  // linear interpolation between end points
        }
    }

    return eta;
}


/*
 * Calculate the gamma coefficient
 *
 * Hurley et al. 2000, eq 23
 *
 *
 * double CalculateGamma(const double p_Mass)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @return                                      The exponent eta (for Hurley et al. 2000, eq 12)
 */
double MainSequence::CalculateGamma(const double p_Mass) const {
#define a m_AnCoefficients                                                      // for convenience and readability - undefined at end of function
#define B_GAMMA m_GammaConstants[static_cast<int>(GAMMA_CONSTANTS::B_GAMMA)]    // for convenience and readability - undefined at end of function
#define C_GAMMA m_GammaConstants[static_cast<int>(GAMMA_CONSTANTS::C_GAMMA)]    // for convenience and readability - undefined at end of function

    double gamma;

         if (utils::Compare(p_Mass,  1.0)          <= 0) gamma = a[76] + (a[77] * PPOW(p_Mass - a[78], a[79]));
    else if (utils::Compare(p_Mass,  a[75])        <= 0) gamma = B_GAMMA + (a[80] - B_GAMMA) * PPOW((p_Mass - 1.0) / (a[75] - 1.0), a[81]);
    else if (utils::Compare(p_Mass, (a[75] + 0.1)) <= 0) gamma = C_GAMMA - (10.0 * (p_Mass - a[75]) * C_GAMMA);                                                             // included = case, missing from Hurley+ 2000
    else                                                 gamma = 0.0;           // this really is zero

    return gamma;

#undef C_GAMMA
#undef B_GAMMA
#undef a
}


/*
 * Calculate luminosity at the end of the Main Sequence
 *
 * Hurley et al. 2000, eq 8
 *
 *
 * double CalculateLuminosityAtPhaseEnd(const double p_Mass)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @return                                      Luminosity at the end of the Main Sequence in Lsol
 */
double MainSequence::CalculateLuminosityAtPhaseEnd(const double p_Mass) const {
#define a m_AnCoefficients    // for convenience and readability - undefined at end of function

    // pow() is slow - use multiplication
    double m_3 = p_Mass * p_Mass * p_Mass;
    double m_4 = m_3 * p_Mass;
    double m_5 = m_4 * p_Mass;

    double top    = (a[11] * m_3) + (a[12] * m_4) + (a[13] * PPOW(p_Mass, (a[16] + 1.8)));
    double bottom = a[14] + (a[15] * m_5) + PPOW(p_Mass, a[16]);

    return top / bottom;

#undef a
}


/*
 * Calculate luminosity on the Main Sequence
 *
 * Hurley et al. 2000, eq 12
 *
 *
 * double CalculateLuminosityOnPhase(const double p_Time, const double p_Mass, const double p_LZAMS)
 *
 * @param   [IN]    p_Time                      Time (after ZAMS) in Myr
 * @param   [IN]    p_Mass                      Mass in Msol
 * @param   [IN]    p_LZAMS0                    Zero Age Main Sequence (ZAMS) Luminosity
 * @return                                      Luminosity on the Main Sequence as a function of time
 */
double MainSequence::CalculateLuminosityOnPhase(const double p_Time, const double p_Mass, const double p_LZAMS) const {
#define a m_AnCoefficients                                          // for convenience and readability - undefined at end of function
#define timescales(x) m_Timescales[static_cast<int>(TIMESCALE::x)]  // for convenience and readability - undefined at end of function

    const double epsilon = 0.01;

    double LTMS   = CalculateLuminosityAtPhaseEnd(p_Mass);
    double alphaL = CalculateAlphaL(p_Mass);
    double betaL  = CalculateBetaL(p_Mass);
    double deltaL = CalculateDeltaL(p_Mass);
    double eta    = CalculateEta(p_Mass);

    double mu     = std::max(0.5, (1.0 - (0.01 * std::max((a[6] / PPOW(p_Mass, a[7])), (a[8] + (a[9] / PPOW(p_Mass, a[10]))))))); // Hurley et al. 2000, eq 7
    double tHook  = mu * timescales(tBGB);                                                                                      // Hurley et al. 2000, just after eq 5
    double tau    = p_Time / timescales(tMS);                                                                                   // Hurley et al. 2000, eq 11
    double tau1   = std::min(1.0, (p_Time / tHook));                                                                            // Hurley et al. 2000, eq 14
    double tau2   = std::max(0.0, std::min(1.0, (p_Time - ((1.0 - epsilon) * tHook)) / (epsilon * tHook)));                     // Hurley et al. 2000, eq 15

    // pow() is slow - use multipliaction where it makes sense
    double logLMS_LZAMS  = alphaL * tau;                                                                                        // Hurley et al. 2000, eq 12, part 1
           logLMS_LZAMS += betaL * PPOW(tau, eta);                                                                              // Hurley et al. 2000, eq 12, part 2
           logLMS_LZAMS += (log10(LTMS / p_LZAMS) - alphaL - betaL) * tau * tau;                                                // Hurley et al. 2000, eq 12, part 3
           logLMS_LZAMS -= deltaL * ((tau1 * tau1) - (tau2 * tau2));                                                            // Hurley et al. 2000, eq 12, part 4

    return p_LZAMS * PPOW(10.0, logLMS_LZAMS);                                                                                  // rewrite Hurley et al. 2000, eq 12 for L(t)

#undef timescales
#undef a
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                RADIUS CALCULATIONS                                //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * Calculate the radius constant alpha_R
 * Hurley et al. 2000, eqs 21a & 21b
 *
 *
 * double CalculateAlphaR(const double p_Mass)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @return                                      Radius constant alpha_R
 */
double MainSequence::CalculateAlphaR(const double p_Mass) const {
#define a m_AnCoefficients    // for convenience and readability - undefined at end of function

    double alphaR = 0.0;

         if (utils::Compare(p_Mass,   0.5) <  0) alphaR = a[62];
    else if (utils::Compare(p_Mass,  0.65) <  0) alphaR = a[62] + (a[63] - a[62]) * (p_Mass - 0.5) / 0.15;
    else if (utils::Compare(p_Mass, a[68]) <  0) alphaR = a[63] + (a[64] - a[63]) * (p_Mass - 0.65) / (a[68] - 0.65);
    else if (utils::Compare(p_Mass, a[66]) <  0) alphaR = a[64] + (m_RConstants[static_cast<int>(R_CONSTANTS::B_ALPHA_R)] - a[64]) * (p_Mass - a[68]) / (a[66] - a[68]);
    else if (utils::Compare(p_Mass, a[67]) <= 0) alphaR = a[58] * PPOW(p_Mass, a[60]) / (a[59] + PPOW(p_Mass, a[61]));
    else                                         alphaR = m_RConstants[static_cast<int>(R_CONSTANTS::C_ALPHA_R)] + a[65] * (p_Mass - a[67]);

    return alphaR;

#undef a
}


/*
 * Calculate the radius constant beta_R
 *
 * Hurley et al. 2000, eqs 22a & 22b
 *
 *
 * double Star::CalculateBetaR(const double p_Mass)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @return                                      Radius constant beta_R
 */
double MainSequence::CalculateBetaR(const double p_Mass) const {
#define a m_AnCoefficients    // for convenience and readability - undefined at end of function

    double betaRPrime = 0.0;

         if (utils::Compare(p_Mass, 1.0)   <= 0) betaRPrime = 1.06;
    else if (utils::Compare(p_Mass, a[74]) <  0) betaRPrime = 1.06 + (a[72] - 1.06) * (p_Mass - 1.0) / (a[74] - 1.06);
    else if (utils::Compare(p_Mass, 2.0)   <  0) betaRPrime = a[72] + (m_RConstants[static_cast<int>(R_CONSTANTS::B_BETA_R)] - a[72]) * (p_Mass - a[74]) / (2.0 - a[74]);
    else if (utils::Compare(p_Mass, 16.0)  <= 0) betaRPrime = (a[69] * p_Mass * p_Mass * p_Mass * std::sqrt(p_Mass)) / (a[70] + PPOW(p_Mass, a[71]));  // pow()is slow - use multiplication (sqrt() is faster than pow())
    else                                         betaRPrime = m_RConstants[static_cast<int>(R_CONSTANTS::C_BETA_R)] + a[73] * (p_Mass - 16.0);

    return betaRPrime - 1.0;

#undef a
}


/*
 * Calculate the value of the radius perturbation DeltaR
 *
 * Hurley et al. 2000, eq 17
 *
 *
 * double CalculateDeltaR(const double p_Mass)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @return                                      The radius perturbation DeltaR
 */
double MainSequence::CalculateDeltaR(const double p_Mass) const {
#define a m_AnCoefficients                                              // for convenience and readability - undefined at end of function
#define massCutoffs(x) m_MassCutoffs[static_cast<int>(MASS_CUTOFF::x)]  // for convenience and readability - undefined at end of function

    double deltaR;

    if (utils::Compare(p_Mass, massCutoffs(MHook)) <= 0) deltaR = 0.0;   // this really is supposed to be 0
    else if (utils::Compare(p_Mass, a[42])         <= 0) deltaR = a[43] * std::sqrt((p_Mass - massCutoffs(MHook)) / (a[42] - massCutoffs(MHook)));
    else if (utils::Compare(p_Mass, 2.0)            < 0) deltaR = a[43] + ((m_RConstants[static_cast<int>(R_CONSTANTS::B_DELTA_R)] - a[43]) * PPOW(((p_Mass - a[42]) / (2.0 - a[42])), a[44]));
    else {
        // pow() is slow - use multiplication (sqrt() is faster than pow())
        double top    = a[38] + (a[39] * p_Mass * p_Mass * p_Mass * std::sqrt(p_Mass));
        double bottom = (a[40] * p_Mass * p_Mass * p_Mass) + PPOW(p_Mass, a[41]);
        deltaR = (top / bottom) - 1.0;
    }

    return deltaR;

#undef massCutoffs
#undef a
}


/*
 * Calculate radius at the end of the Main Sequence
 *
 * Hurley et al. 2000, eqs 9a & 9b
 *
 *
 * double CalculateRadiusAtPhaseEnd(const double p_Mass, const double p_RZAMS)
 *
 * @param   [IN]    p_Mass                      Stellar mass (Msol)
 * @param   [IN]    p_RZAMS                     Zero Age Main Sequence (ZAMS) Radius
 * @return                                      Radius at the end of the Main Sequence in Rsol
 */
double MainSequence::CalculateRadiusAtPhaseEnd(const double p_Mass, const double p_RZAMS) const {
#define a m_AnCoefficients    // for convenience and readability - undefined at end of function

    double RTMS;
    double mAsterisk = a[17] + 0.1;

    if (utils::Compare(p_Mass, a[17]) <= 0) {
        RTMS = (a[18] + (a[19] * PPOW(p_Mass, a[21]))) / (a[20] + PPOW(p_Mass, a[22]));

        if (utils::Compare(p_Mass, 0.5) < 0) {
            RTMS = std::max(RTMS, 1.5 * p_RZAMS);
        }
    }
    else if (utils::Compare(p_Mass, mAsterisk) >= 0) {
        // pow() is slow - use multiplication
        double m_3 = p_Mass * p_Mass * p_Mass;
        double m_5 = m_3 * p_Mass * p_Mass;

        RTMS = ((C_COEFF.at(1) * m_3) + (a[23] * PPOW(p_Mass, a[26])) + (a[24] * PPOW(p_Mass, a[26] + 1.5))) / (a[25] + m_5);
    }
    else {
        // for stars with masses between a17, a17 + 0.1 interpolate between the end points (y = mx + c)

        // pow() is slow - use multiplication
        double mA_3 = mAsterisk * mAsterisk * mAsterisk;
        double mA_5 = mA_3 * mAsterisk * mAsterisk;

        double y2   = ((C_COEFF.at(1) * mA_3) + (a[23] * PPOW(mAsterisk, a[26])) + (a[24] * PPOW(mAsterisk, a[26] + 1.5))) / (a[25] + mA_5);    // RTMS(mAsterisk)
        double y1   = (a[18] + (a[19] * PPOW(a[17], a[21]))) / (a[20] + PPOW(a[17], a[22]));                                                    // RTMS(a17)

        double gradient  = (y2 - y1) / 0.1;
        double intercept = y1 - (gradient * a[17]);

        RTMS = (gradient * p_Mass) + intercept;
    }

    return RTMS;

#undef a
}


/*
 * Calculate radius on the Main Sequence
 *
 * Hurley et al. 2000, eq 13
 *
 *
 * double CalculateRadiusOnPhase(const double p_Mass, const double p_Time, const double p_RZAMS)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @param   [IN]    p_Time                      Time (after ZAMS) in Myr
 * @param   [IN]    p_RZAMS                     Zero Age Main Sequence (ZAMS) Radius
 * @return                                      Radius on the Main Sequence in Rsol
 */
double MainSequence::CalculateRadiusOnPhase(const double p_Mass, const double p_Time, const double p_RZAMS) const {
#define a m_AnCoefficients                                          // for convenience and readability - undefined at end of function
#define timescales(x) m_Timescales[static_cast<int>(TIMESCALE::x)]  // for convenience and readability - undefined at end of function

    const double epsilon = 0.01;

    double RTMS   = CalculateRadiusAtPhaseEnd(p_Mass, p_RZAMS);
    double alphaR = CalculateAlphaR(p_Mass);
    double betaR  = CalculateBetaR(p_Mass);
    double deltaR = CalculateDeltaR(p_Mass);
    double gamma  = CalculateGamma(p_Mass);

    double mu     = std::max(0.5, (1.0 - (0.01 * std::max((a[6] / PPOW(p_Mass, a[7])), (a[8] + (a[9] / PPOW(p_Mass, a[10]))))))); // Hurley et al. 2000, eq 7
    double tHook  = mu * timescales(tBGB);                                                                                      // Hurley et al. 2000, just after eq 5
    double tau    = p_Time / timescales(tMS);                                                                                   // Hurley et al. 2000, eq 11
    double tau1   = std::min(1.0, (p_Time / tHook));                                                                            // Hurley et al. 2000, eq 14
    double tau2   = std::max(0.0, std::min(1.0, (p_Time - ((1.0 - epsilon) * tHook)) / (epsilon * tHook)));                     // Hurley et al. 2000, eq 15

    // pow() is slow - use multipliaction where it makes sense
    double tau_3  = tau * tau * tau;
    double tau_10 = tau_3 * tau_3 * tau_3 * tau;
    double tau_40 = tau_10 * tau_10 * tau_10 * tau_10;
    double tau1_3 = tau1 * tau1 * tau1;
    double tau2_3 = tau2 * tau2 * tau2;

    double logRMS_RZAMS  = alphaR * tau;                                                                                        // Hurley et al. 2000, eq 13, part 1
           logRMS_RZAMS += betaR * tau_10;                                                                                      // Hurley et al. 2000, eq 13, part 2
           logRMS_RZAMS += gamma * tau_40;                                                                                      // Hurley et al. 2000, eq 13, part 3
           logRMS_RZAMS += (log10(RTMS / p_RZAMS) - alphaR - betaR - gamma) * tau_3;                                            // Hurley et al. 2000, eq 13, part 4
           logRMS_RZAMS -= deltaR * (tau1_3 - tau2_3);                                                                          // Hurley et al. 2000, eq 13, part 5

    return p_RZAMS * PPOW(10.0, logRMS_RZAMS);                                                                                   // rewrite Hurley et al. 2000, eq 13 for R(t)

#undef timescales
#undef a
}


/*
 * Calculate radius on the Main Sequence
 *
 * Hurley et al. 2000, eq 13
 *
 *
 * double CalculateRadiusOnPhaseTau(const double p_Mass, const double p_Tau)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @param   [IN]    p_Tau                       Fractional age on Main Sequence
 * @return                                      Radius on the Main Sequence in Rsol
 */
double MainSequence::CalculateRadiusOnPhaseTau(const double p_Mass, const double p_Tau) const {
#define a m_AnCoefficients                                          // for convenience and readability - undefined at end of function

    const double epsilon = 0.01;
    double tBGB = CalculateLifetimeToBGB(p_Mass);
    double tMS  = CalculateLifetimeOnPhase(p_Mass, tBGB);

    double RZAMS  = CalculateRadiusAtZAMS(p_Mass);
    double RTMS   = CalculateRadiusAtPhaseEnd(p_Mass, RZAMS);
    double alphaR = CalculateAlphaR(p_Mass);
    double betaR  = CalculateBetaR(p_Mass);
    double deltaR = CalculateDeltaR(p_Mass);
    double gamma  = CalculateGamma(p_Mass);

    double mu     = std::max(0.5, (1.0 - (0.01 * std::max((a[6] / PPOW(p_Mass, a[7])), (a[8] + (a[9] / PPOW(p_Mass, a[10]))))))); // Hurley et al. 2000, eq 7
    double tHook  = mu * tBGB;                                                                                      // Hurley et al. 2000, just after eq 5
    double time   = tMS * p_Tau;
    double tau1   = std::min(1.0, (time / tHook));                                                                            // Hurley et al. 2000, eq 14
    double tau2   = std::max(0.0, std::min(1.0, (time - ((1.0 - epsilon) * tHook)) / (epsilon * tHook)));                     // Hurley et al. 2000, eq 15

    // pow() is slow - use multipliaction where it makes sense
    double tau_3  = p_Tau * p_Tau * p_Tau;
    double tau_10 = tau_3 * tau_3 * tau_3 * p_Tau;
    double tau_40 = tau_10 * tau_10 * tau_10 * tau_10;
    double tau1_3 = tau1 * tau1 * tau1;
    double tau2_3 = tau2 * tau2 * tau2;

    double logRMS_RZAMS  = alphaR * p_Tau;                                                                                        // Hurley et al. 2000, eq 13, part 1
           logRMS_RZAMS += betaR * tau_10;                                                                                      // Hurley et al. 2000, eq 13, part 2
           logRMS_RZAMS += gamma * tau_40;                                                                                      // Hurley et al. 2000, eq 13, part 3
           logRMS_RZAMS += (log10(RTMS / RZAMS) - alphaR - betaR - gamma) * tau_3;                                            // Hurley et al. 2000, eq 13, part 4
           logRMS_RZAMS -= deltaR * (tau1_3 - tau2_3);                                                                          // Hurley et al. 2000, eq 13, part 5

    return RZAMS * PPOW(10.0, logRMS_RZAMS);                                                                                   // rewrite Hurley et al. 2000, eq 13 for R(t)

#undef a
}




/*
 * Calculate the radial extent of the star's convective envelope (if it has one)
 *
 * Hurley et al. 2002, sec. 2.3, particularly subsec. 2.3.1, eqs 36-38
 *
 *
 * double CalculateRadialExtentConvectiveEnvelope()
 *
 * @return                                      Radial extent of the star's convective envelope in Rsol
 */
double MainSequence::CalculateRadialExtentConvectiveEnvelope() const {
    double radiusEnvelope0 = m_Radius;
    if ( utils::Compare(m_Mass, 1.25) >= 0)
        radiusEnvelope0 = 0.0;
    else if (utils::Compare(m_Mass, 0.35) > 0) {
        double radiusM035 = CalculateRadiusAtZAMS(0.35);        // uses radius of a 0.35 solar mass star at ZAMS rather than at fractional age Tau, but such low-mass stars only grow by a maximum factor of 1.5 [just above Eq. (10) in Hurley, Pols, Tout (2000), so this is a reasonable approximation
        radiusEnvelope0   = radiusM035 * std::sqrt((1.25 - m_Mass) / 0.9);
    }

    return radiusEnvelope0 * std::sqrt(std::sqrt(1.0 - m_Tau));
}

/*
 * Calculate the radial extent of the star's convective core (if it has one)
 *
 * Uses preliminary fit from Minori Shikauchi @ ZAMS, then a smooth interpolation to the HG
 *
 *
 * double CalculateRadialExtentConvectiveEnvelope()
 *
 * @return                                      Radial extent of the star's convective core in Rsol
 */
double MainSequence::CalculateConvectiveCoreRadius() const {
    if(utils::Compare(m_Mass, 1.25) < 0) return 0.0;                                            // low-mass star with a radiative core
       
    double convectiveCoreRadiusZAMS = m_Mass * (0.06 + 0.05 * exp(-m_Mass / 61.57));
    
    // We need TAMSCoreRadius, which is just the core radius at the start of the HG phase.
    // Since we are on the main sequence here, we can clone this object as an HG object
    // and, as long as it is initialised (to correctly set Tau to 0.0 on the HG phase),
    // we can query the cloned object for its core mass.
    //
    // The clone should not evolve, and so should not log anything, but to be sure the
    // clone does not participate in logging, we set its persistence to EPHEMERAL.

    HG *clone = HG::Clone(static_cast<HG&>(const_cast<MainSequence&>(*this)), OBJECT_PERSISTENCE::EPHEMERAL);
    double TAMSCoreRadius = clone->CalculateRemnantRadius();                                    // get core radius from clone
    delete clone; clone = nullptr;                                                              // return the memory allocated for the clone

    return (convectiveCoreRadiusZAMS - m_Tau * (convectiveCoreRadiusZAMS - TAMSCoreRadius));
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                 MASS CALCULATIONS                                 //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * Calculate the mass of the convective core
 *
 * Based on Shikauchi, Hirai, Mandel (2024), core mass shrinks to 60% of initial value over the course of the MS
 *
 *
 * double CalculateConvectiveCoreMass() const
 *
 * @return                                      Mass of convective core in Msol
 */
double MainSequence::CalculateConvectiveCoreMass() const {

    // We need TAMSCoreMass, which is just the core mass at the start of the HG phase.
    // Since we are on the main sequence here, we can clone this object as an HG object
    // and, as long as it is initialised (to correctly set Tau to 0.0 on the HG phase),
    // we can query the cloned object for its core mass.
    //
    // The clone should not evolve, and so should not log anything, but to be sure the
    // clone does not participate in logging, we set its persistence to EPHEMERAL.
      
    HG *clone = HG::Clone(static_cast<HG&>(const_cast<MainSequence&>(*this)), OBJECT_PERSISTENCE::EPHEMERAL);
    double TAMSCoreMass = clone->CoreMass();                                                    // get core mass from clone
    delete clone; clone = nullptr;                                                              // return the memory allocated for the clone

    double finalConvectiveCoreMass   = TAMSCoreMass;
    double initialConvectiveCoreMass = finalConvectiveCoreMass / 0.6;
    return (initialConvectiveCoreMass - m_Tau * (initialConvectiveCoreMass - finalConvectiveCoreMass));
}

/*
 * Calculate the mass of the convective envelope
 *
 * Based on section 7.2 (after Eq. 111) of Hurley, Pols, Tout (2000)
 *
 *
 * double CalculateConvectiveEnvelopeMass() const
 *
 * @return                                      Mass of convective envelope in Msol
 */
DBL_DBL MainSequence::CalculateConvectiveEnvelopeMass() const {
    if (utils::Compare(m_Mass, 1.25) > 0) return std::tuple<double, double> (0.0, 0.0);

    double massEnvelope0 = m_Mass;
    if (utils::Compare(m_Mass, 0.35) > 0) massEnvelope0 = 0.35 * (1.25 - m_Mass) * (1.25 - m_Mass) / 0.81;
    
    double massEnvelope  = massEnvelope0 * sqrt(sqrt(1.0 - m_Tau));
    
    return std::tuple<double, double> (massEnvelope, massEnvelope0);
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                            LIFETIME / AGE CALCULATIONS                            //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * Calculate relative age on the Main Sequence
 *
 * Hurley et al. 2000, eq 11
 * Naturally bounded by [0, 1], but clamp here anyway
 *
 * double CalculateTauOnPhase()
 *
 * @return                                      MS relative age, clamped to [0, 1]
 */
double MainSequence::CalculateTauOnPhase() const {
#define timescales(x) m_Timescales[static_cast<int>(TIMESCALE::x)]  // for convenience and readability - undefined at end of function
    return std::max(0.0, std::min(1.0, m_Age / timescales(tMS)));
#undef timescales
}


/*
 * Calculate lifetime of Main Sequence
 *
 * Hurley et al. 2000, eq 5
 *
 *
 * double CalculatePhaseLifetime(const double p_Mass, const double p_TBGB)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @param   [IN]    p_TBGB                      Lifetime to Base of Giant Branch
 * @return                                      Lifetime of Main Sequence in Myr
 */
double MainSequence::CalculateLifetimeOnPhase(const double p_Mass, const double p_TBGB) const {
#define a m_AnCoefficients    // for convenience and readability - undefined at end of function

    // Calculate time to Hook
    // Hurley et al. 2000, eqs 5, 6 & 7
    double mu    = std::max(0.5, (1.0 - (0.01 * std::max((a[6] / PPOW(p_Mass, a[7])), (a[8] + (a[9] / PPOW(p_Mass, a[10])))))));
    double tHook = mu * p_TBGB;

    // For mass < Mhook, x > mu (i.e. for stars without a hook)
    double x = std::max(0.95, std::min((0.95 - (0.03 * (LogMetallicityXi() + 0.30103))), 0.99));

    return std::max(tHook, (x * p_TBGB));

#undef a
}


/*
 * Recalculates the star's age after mass loss
 *
 * Hurley et al. 2000, section 7.1
 *
 * Modifies attribute m_Age
 *
 *
 * UpdateAgeAfterMassLoss()
 *
 */
void MainSequence::UpdateAgeAfterMassLoss() {

    double tMS       = m_Timescales[static_cast<int>(TIMESCALE::tMS)];
    double tBGBprime = CalculateLifetimeToBGB(m_Mass);
    double tMSprime  = MainSequence::CalculateLifetimeOnPhase(m_Mass, tBGBprime);

    m_Age *= tMSprime / tMS;
}


bool MainSequence::ShouldEvolveOnPhase() const {
    switch (OPTIONS->MainSequenceCoreMassPrescription()) {
        case CORE_MASS_PRESCRIPTION::MANDEL: {
            return (m_Age < m_Timescales[static_cast<int>(TIMESCALE::tMS)]);
        }
        case CORE_MASS_PRESCRIPTION::NONE: {
            return (m_Age < m_Timescales[static_cast<int>(TIMESCALE::tMS)]);
        }
        case CORE_MASS_PRESCRIPTION::SHIKAUCHI: {
            return (m_CentralHeliumFraction <= 1.0 - m_Metallicity & m_Age < m_Timescales[static_cast<int>(TIMESCALE::tMS)]);
        }
    }
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                    MISCELLANEOUS FUNCTIONS / CONTROL FUNCTIONS                    //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * Set parameters required before evolving one timestep - modify star values
 *
 *
 * void EvolveOneTimestepPreamble()
 */
void MainSequence::EvolveOneTimestepPreamble() {
    m_LZAMS0 = CalculateLuminosityAtZAMS(m_Mass0);
    m_RZAMS0 = CalculateRadiusAtZAMS(m_Mass0);
}


/*
 * Choose timestep for evolution
 *
 * Given in the discussion in Hurley et al. 2000
 *
 *
 * ChooseTimestep(const double p_Time)
 *
 * @param   [IN]    p_Time                      Current age of star in Myr
 * @return                                      Suggested timestep (dt)
 */
double MainSequence::ChooseTimestep(const double p_Time) const {
#define timescales(x) m_Timescales[static_cast<int>(TIMESCALE::x)]  // for convenience and readability - undefined at end of function

    double dtk = 1.0E-2 * timescales(tMS);  // 0.01 of MS timescale (sse uses 0.05)
    double dte = timescales(tMS) - p_Time;  // time remaining on MS

    if (utils::Compare(dte, dtk) < 0) {     // short enough to resolve the hook at the end of the MS for HM stars? JAR: why not check for HM star?
        dtk /= 10.0;                        // no - go an order-of-magnitude shorter
    }

    return std::max(std::min(dtk, dte), NUCLEAR_MINIMUM_TIMESTEP);

#undef timescales
}


/*
 * Resolve changes to the remnant after the star loses its envelope
 *
 * Where necessary updates attributes of star (depending upon stellar type):
 *
 *     - m_StellarType
 *     - m_Timescales
 *     - m_GBParams
 *     - m_Luminosity
 *     - m_Radius
 *     - m_Mass
 *     - m_Mass0
 *     - m_CoreMass
 *     - m_HeCoreMass
 *     - m_COCoreMass
 *     - m_Age
 *
 *
 * STELLAR_TYPE ResolveEnvelopeLoss(bool p_Force)
 *
 * @param   [IN]    p_Force                     Boolean to indicate whether the resolution of the loss of the envelope should be performed
 *                                              without checking the precondition(s).
 *                                              Default is false.
 *
 * @return                                      Stellar type to which star should evolve
 */
STELLAR_TYPE MainSequence::ResolveEnvelopeLoss(bool p_Force) {

    STELLAR_TYPE stellarType = m_StellarType;
    
    if (p_Force || utils::Compare(m_Mass, 0.0) <= 0) {      // envelope loss
        stellarType = STELLAR_TYPE::MASSLESS_REMNANT;
        m_Radius    = 0.0;
        m_Mass      = 0.0;
    }
    
    return stellarType;
}


double MainSequence::CalculateMainSequenceCoreMassMandel() {
    // We need TAMSCoreMass, which is just the core mass at the start of the HG phase.
    // Since we are on the main sequence here, we can clone this object as an HG object
    // and, as long as it is initialised (to correctly set Tau to 0.0 on the HG phase),
    // we can query the cloned object for its core mass.
    //
    // The clone should not evolve, and so should not log anything, but to be sure the
    // clone does not participate in logging, we set its persistence to EPHEMERAL.
    
    HG *clone = HG::Clone(*this, OBJECT_PERSISTENCE::EPHEMERAL);
    double TAMSCoreMass = clone->CoreMass();                                                    // get core mass from clone
    delete clone; clone = nullptr;                                                              // return the memory allocated for the clone
    double minimumCoreMass   = std::max(m_MinimumCoreMass, CalculateTauOnPhase() * TAMSCoreMass);    // update minimum core mass
    //double coreMass = std::max(m_MinimumCoreMass, CalculateTauOnPhase() * TAMSCoreMass);
    return minimumCoreMass;
}


double MainSequence::CalculateMainSequenceCoreMassShikauchi() {
    double totalMass = m_Mass;
    double centralHeliumFraction = m_CentralHeliumFraction;
    double mixingCoreMass = m_MinimumCoreMass;
    double lnMixingCoreMass = std::log(mixingCoreMass);
    double logMixingCoreMass = std::log10(mixingCoreMass);
    
    double fmix = SHIKAUCHI_FMIX_COEFFICIENTS[0][0] + SHIKAUCHI_FMIX_COEFFICIENTS[0][1] * std::exp(- m_MZAMS / SHIKAUCHI_FMIX_COEFFICIENTS[0][2]);
    
    double beta = 1.0 - SHIKAUCHI_FMIX_COEFFICIENTS[0][1] * totalMass / (SHIKAUCHI_FMIX_COEFFICIENTS[0][2] * fmix) * std::exp(- totalMass / SHIKAUCHI_FMIX_COEFFICIENTS[0][2]);
    
    double logL = SHIKAUCHI_L_COEFFICIENTS[0][0] * logMixingCoreMass + SHIKAUCHI_L_COEFFICIENTS[0][1] * centralHeliumFraction + SHIKAUCHI_L_COEFFICIENTS[0][2] * logMixingCoreMass * centralHeliumFraction + SHIKAUCHI_L_COEFFICIENTS[0][3] * logMixingCoreMass * logMixingCoreMass + SHIKAUCHI_L_COEFFICIENTS[0][4] * centralHeliumFraction * centralHeliumFraction + SHIKAUCHI_L_COEFFICIENTS[0][5] * logMixingCoreMass * logMixingCoreMass * logMixingCoreMass + SHIKAUCHI_L_COEFFICIENTS[0][6] * centralHeliumFraction * centralHeliumFraction * centralHeliumFraction + SHIKAUCHI_L_COEFFICIENTS[0][7] * logMixingCoreMass * logMixingCoreMass * centralHeliumFraction + SHIKAUCHI_L_COEFFICIENTS[0][8] * logMixingCoreMass * centralHeliumFraction * centralHeliumFraction + SHIKAUCHI_L_COEFFICIENTS[0][9];
    
    double alpha = PPOW(10.0, std::max(-2.0, SHIKAUCHI_ALPHA_COEFFICIENTS[0][1] * mixingCoreMass + SHIKAUCHI_ALPHA_COEFFICIENTS[0][2])) + SHIKAUCHI_ALPHA_COEFFICIENTS[0][0];
    
    double g = -0.0044 * m_MZAMS + 0.27;
    
    double YZAMS = 0.24 + 2.0 * m_Metallicity;
    
    double Yhat = (centralHeliumFraction - YZAMS) / (1.0 - YZAMS - m_Metallicity);
    
    double delta = std::min(PPOW(10.0, -Yhat + g), 1.0);
    
    double currentTimestepInYrs = m_Dt * 1.0E6;
    
    double mDot = m_TotalMassLossRate;  //(m_MassPrev - m_Mass) / previousTimestepInYrs;
    
    controlled_stepper_type controlled_stepper;
    state_type x(3);
    x[0] = totalMass;
    x[1] = centralHeliumFraction;
    x[2] = lnMixingCoreMass;
    
    auto ode = [&](const state_type &x, state_type &dxdt, const double) {
        dxdt[0] = - mDot;
        dxdt[1] = PPOW(10.0, logL) / (Q_CNO * mixingCoreMass);
        dxdt[2] = -alpha/(1.0 - alpha * x[1]) * dxdt[1] + beta * delta * dxdt[0] / x[0];
    };
    
    integrate_adaptive(controlled_stepper, ode, x, 0.0, currentTimestepInYrs, currentTimestepInYrs/100.0);
    
    mixingCoreMass = std::exp(x[2]);
    m_CentralHeliumFraction = x[1];
    
    return mixingCoreMass;
}


/*
 * Update the minimum core mass of a main sequence star that loses mass through Case A mass transfer by
 * setting it equal to the core mass of a TAMS star, scaled by the fractional age.
 *
 * The minimum core mass of the star is updated only if the retain-core-mass-during-caseA-mass-transfer
 * option is specified, otherwise it is left unchanged.
 *
 *
 * void UpdateMinimumCoreMass()
 *
 */
void MainSequence::UpdateMinimumCoreMass() {
    switch (OPTIONS->MainSequenceCoreMassPrescription()) {
        case CORE_MASS_PRESCRIPTION::MANDEL:
            m_MinimumCoreMass = CalculateMainSequenceCoreMassMandel();
            break;
        case CORE_MASS_PRESCRIPTION::NONE:
            m_MinimumCoreMass = 0.0;
            break;
        case CORE_MASS_PRESCRIPTION::SHIKAUCHI:
            m_MinimumCoreMass = CalculateMainSequenceCoreMassShikauchi();
            break;
    }
}


/*
 * Sets the mass and age of a merge product of two main sequence stars
 * (note: treats merger products as main-sequence stars, not CHE, and does not check for MS_lte_07)
 *
 * Uses prescription of Wang et al., 2022, https://www.nature.com/articles/s41550-021-01597-5
 *
 * void UpdateAfterMerger(double p_Mass, double p_HydrogenMass)
 *
 * @param   [IN]    p_Mass                      New value of stellar mass
 * @param   [IN]    p_HydrogenMass              Desired value of hydrogen mass of merger remnant
 *
 */
void MainSequence::UpdateAfterMerger(double p_Mass, double p_HydrogenMass) {
    #define timescales(x) m_Timescales[static_cast<int>(TIMESCALE::x)]  // for convenience and readability - undefined at end of function

    m_Mass            = p_Mass;
    m_Mass0           = m_Mass;
    m_MinimumCoreMass = 0.0;
    
    double initialHydrogenFraction = m_InitialHydrogenAbundance;
    
    CalculateTimescales();
    CalculateGBParams();
            
    m_Tau = (initialHydrogenFraction - p_HydrogenMass / m_Mass) / initialHydrogenFraction;       // assumes uniformly mixed merger product and a uniform rate of H fusion on main sequence
    
    m_Age = m_Tau * timescales(tMS);
    
    UpdateAttributesAndAgeOneTimestep(0.0, 0.0, 0.0, true);
    
    #undef timescales
}
