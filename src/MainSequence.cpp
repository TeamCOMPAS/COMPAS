#include "MainSequence.h"


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//             PARAMETERS, MISCELLANEOUS CALCULATIONS AND FUNCTIONS ETC.             //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


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
    timescales(tBGB)   = CalculateLifetimeToBGB(p_Mass);
    timescales(tMS)    = CalculateLifetimeOnPhase(p_Mass, timescales(tBGB));
    timescales(tMcMax) = 0.0;                                       // JR: todo: this is never used - is it actually needed?

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

    if (utils::Compare(p_Mass, massCutoffs(MHook)) <= 0) {              // JR: todo: added "=" per Hurley et al. 2000, eq 16
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
 * Calcluate the luminosity alpha constant alpha_L
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
           logLMS_LZAMS += betaL * PPOW(tau, eta);                                                                               // Hurley et al. 2000, eq 12, part 2
           logLMS_LZAMS += (log10(LTMS / p_LZAMS) - alphaL - betaL) * tau * tau;                                                // Hurley et al. 2000, eq 12, part 3
           logLMS_LZAMS -= deltaL * ((tau1 * tau1) - (tau2 * tau2));                                                            // Hurley et al. 2000, eq 12, part 4

    return p_LZAMS * PPOW(10.0, logLMS_LZAMS);                                                                                   // rewrite Hurley et al. 2000, eq 12 for L(t)

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
    else if (utils::Compare(p_Mass, 16.0)  <= 0) betaRPrime = (a[69] * p_Mass * p_Mass * p_Mass * sqrt(p_Mass)) / (a[70] + PPOW(p_Mass, a[71]));  // pow()is slow - use multiplication (sqrt() is faster than pow())
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
    else if (utils::Compare(p_Mass, a[42])         <= 0) deltaR = a[43] * PPOW(((p_Mass - massCutoffs(MHook)) / (a[42] - massCutoffs(MHook))), 0.5);
    else if (utils::Compare(p_Mass, 2.0)            < 0) deltaR = a[43] + ((m_RConstants[static_cast<int>(R_CONSTANTS::B_DELTA_R)] - a[43]) * PPOW(((p_Mass - a[42]) / (2.0 - a[42])), a[44]));
    else {
        // pow() is slow - use multiplication (sqrt() is faster than pow())
        double top    = a[38] + (a[39] * p_Mass * p_Mass * p_Mass * sqrt(p_Mass));
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
    else{   // for stars with masses between a17, a17 + 0.1 interpolate between the end points (y = mx + c)

        // pow() is slow - use multiplication
        double mA_3 = mAsterisk * mAsterisk * mAsterisk;
        double mA_5 = mA_3 * mAsterisk * mAsterisk;

        double y2   = ((C_COEFF.at(1) * mA_3) + (a[23] * PPOW(mAsterisk, a[26])) + (a[24] * PPOW(mAsterisk, a[26] + 1.5))) / (a[25] + mA_5); // RTMS(mAsterisk)
        double y1   = (a[18] + (a[19] * PPOW(a[17], a[21]))) / (a[20] + PPOW(a[17], a[22]));                                                  // RTMS(a17)

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
 * double CalculateRadiusOnPhase(const double p_Mass, const double, p_RZAMS, const double p_Time)
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
 * Calculate the radial extent of the star's convective envelope (if it has one)
 *
 * Hurley et al. 2000, sec. 2.3, particularly subsec. 2.3.1, eqs 36-40
 *
 * (Technically not a radius calculation I suppose, but "radial extent" is close enough to put it with the radius calculations...)
 *
 * JR: todo: original code for MS is broken for mass < 1.25 - check this (see calculateRadialExtentConvectiveEnvelope())
 *
 *
 * double CalculateRadialExtentConvectiveEnvelope()
 *
 * @return                                      Radial extent of the star's convective envelope in Rsol
 */
double MainSequence::CalculateRadialExtentConvectiveEnvelope() const {
    return utils::Compare(m_Mass, 0.35) <= 0 ? m_Radius * PPOW(1.0 - m_Tau, 1.0 / 4.0) : 0.0;
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
    double x = std::max(0.95, std::min((0.95 - (0.03 * (m_LogMetallicityXi + 0.30103))), 0.99));

    return std::max(tHook, (x * p_TBGB));

#undef a
}


/*
 * Calculate thermal timescale
 *
 * Kalogera & Webbink 1996, eq 2 [note that (61) of BSE proposes a value a factor of 3 smaller]
 *
 *
 * double CalculateThermalTimescale(const double p_Mass, const double p_Radius, const double p_Luminosity) const
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @param   [IN]    p_Radius                    Radius in Rsol
 * @param   [IN]    p_Luminosity                Luminosity in Lsol
 * @param   [IN]    p_EnvMass                   Envelope mass in Msol (ignored here)
 * @return                                      Thermal timescale in Myr
 */
double MainSequence::CalculateThermalTimescale(const double p_Mass, const double p_Radius, const double p_Luminosity, const double p_EnvMass) const {
    return 30.0 * p_Mass * p_Mass / (p_Radius * p_Luminosity);      // G*Msol^2/(Lsol*Rsol) ~ 30 Myr
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


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                  ROTATIONAL / GYRATION / FREQUENCY CALCULATIONS                   //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

/*
 * Calculate gyration radius
 *
 * Define gyration radius 'k=r_g^2' using fit from de Mink et al. 2013, calling k_definition function
 * Original created by Alejandro Vigna-Gomez on 11/2015.  Rewritten June 2019, JR.
 *
 * The original fits from de Mink+2013 were made for MS stars a Z=0.02.
 *
 * Uses class member variables instaed of passing in parameters
 *
 *
 * double CalculateGyrationRadius()
 *
 * @return                                      Gyration radius in Rsol
 *
 */
double MainSequence::CalculateGyrationRadius() const {

    double log10M = log10(m_Mass);

	double cLower = 0.0;                                                                            // correction factor 'c' (lowercase 'c') in de Mink et al., 2013 eq A1
	if ((utils::Compare(log10M, 1.3) > 0)) {                                                        // log10(M) > 1.3 (de Mink doesn't include '=' - we assume it not to be here))
        double log10M_13 = log10M - 1.3;
        cLower = -0.055 * log10M_13 * log10M_13;
	}

	double CUpper = -2.5;                                                                           // exponent 'C' (uppercase 'C') in de Mink et al., 2013 eq A2
	     if ((utils::Compare(log10M, 0.2) > 0)) CUpper = -1.5;                                      // log10(M) > 0.2
	else if ((utils::Compare(log10M, 0.0) > 0)) CUpper = -2.5 + (5.0 * log10M);                     // 0.2 <= log10(M) > 0.0 (de Mink doesn't include '=' - we assume it to be here (and for log10(M) <= 0.0))

    double k0 = cLower+ std::min(0.21, std::max(0.09 - (0.27 * log10M), 0.037 + (0.033 * log10M))); // gyration radius squared for ZAMS stars

    double radiusRatio = m_Radius / m_RZAMS;

	return ((k0 - 0.025) * PPOW(radiusRatio, CUpper)) + (0.025 * PPOW(radiusRatio, -0.1));          // gyration radius
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
 * Can obviously do this your own way
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

    double dtk = 1.0E-2 * timescales(tMS);
    double dte = timescales(tMS) - p_Time;

    if (utils::Compare(dte, dtk) < 0) {     // short enough to resolve the hook at the end of the MS for HM stars?
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
 * Hurley et al. 2000, just after eq 105
 *
 * JR: todo: why is this different from ResolveEnvelopeLoss()?
 * JR: todo: original code: Star::radiusRemnantStarAfterLosingEnvelope() vs Star::modifyStarAfterLosingEnvelope(int stellarType, double mass)
 * JR: todo: why is stellar type changed for some types, but not others?  CheB and EAGB stars have stellar type changed, but no other types do...
 * JR: todo: probably not a huge issue - only called in TIDES() and ResolveRemnantAfterEnvelopeLoss(), and with a copy of the star - probably ok there that attributes are changed (except maybe TIDES()?)
 *
 *
 * STELLAR_TYPE ResolveRemnantAfterEnvelopeLoss()
 *
 * @return                                      Stellar type to which star should evolve
 */
STELLAR_TYPE MainSequence::ResolveRemnantAfterEnvelopeLoss() {

    if (utils::Compare(m_Mass, 0.0) <= 0) m_Radius = 0.0;   // massless remnant

    return m_StellarType;                                   // no change to stellar type
}
