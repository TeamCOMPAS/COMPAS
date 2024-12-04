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
    
    // If SHIKAUCHI core mass prescription is used, core helium abundance is calculated with the core mass
    if ((OPTIONS->MainSequenceCoreMassPrescription() == CORE_MASS_PRESCRIPTION::SHIKAUCHI) && (m_MZAMS >= 10.0))
        return m_HeliumAbundanceCore;
    
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
    
    // If SHIKAUCHI core mass prescription is used, core helium abundance is calculated with the core mass
    if ((OPTIONS->MainSequenceCoreMassPrescription() == CORE_MASS_PRESCRIPTION::SHIKAUCHI) && (m_MZAMS >= 10.0))
        return 1.0 - m_HeliumAbundanceCore - m_Metallicity;
    
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
    
    // If SHIKAUCHI core prescription is used, return Shikauchi luminosity
    if ((OPTIONS->MainSequenceCoreMassPrescription() == CORE_MASS_PRESCRIPTION::SHIKAUCHI) && (m_MZAMS >= 10.0)) {
            return CalculateLuminosityShikauchi(m_MainSequenceCoreMass, m_HeliumAbundanceCore, p_Time);
    }
    
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


/*
 * Calculate luminosity on the Main Sequence when Shikauchi et al. (2024) core mass prescription is used
 *
 * During core hydrogen burning uses eq (A5) from Shikauchi et al. (2024)
 *
 * When the Main Sequence hook starts (at age 0.99 * tMS) calculates luminosity that smoothly connects the last point
 * of core hydrogen burning with the first point of the HG
 *
 * double CalculateLuminosityShikauchi(const double p_CoreMass, const double p_HeliumAbundanceCore, const double p_Age)
 *
 * @param   [IN]    p_CoreMass                  Main sequence core mass in Msol
 * @param   [IN]    p_HeliumAbundanceCore       Central helium fraction
 * @param   [IN]    p_Age                       Current age in Myr
 * @return                                      Luminosity on the Main Sequence as a function of current core mass and central helium fraction
 */
double MainSequence::CalculateLuminosityShikauchi(const double p_CoreMass, const double p_HeliumAbundanceCore, const double p_Age) const {
    DBL_VECTOR L_COEFFICIENTS    = std::get<2>(SHIKAUCHI_COEFFICIENTS);
    double centralHeliumFraction = p_HeliumAbundanceCore;
    double logMixingCoreMass     = std::log10(p_CoreMass);
    
    double logL = L_COEFFICIENTS[0] * logMixingCoreMass + L_COEFFICIENTS[1] * centralHeliumFraction + L_COEFFICIENTS[2] * logMixingCoreMass * centralHeliumFraction + L_COEFFICIENTS[3] * logMixingCoreMass * logMixingCoreMass + L_COEFFICIENTS[4] * centralHeliumFraction * centralHeliumFraction + L_COEFFICIENTS[5] * logMixingCoreMass * logMixingCoreMass * logMixingCoreMass + L_COEFFICIENTS[6] * centralHeliumFraction * centralHeliumFraction * centralHeliumFraction + L_COEFFICIENTS[7] * logMixingCoreMass * logMixingCoreMass * centralHeliumFraction + L_COEFFICIENTS[8] * logMixingCoreMass * centralHeliumFraction * centralHeliumFraction + L_COEFFICIENTS[9];
    double luminosity = PPOW(10.0, logL);
    
    if ((p_HeliumAbundanceCore == 1.0 - m_Metallicity) && (m_TotalMassLossRate == -m_Mdot)) {                                // star in MS hook and mass transfer not ongoing?
        HG *clone = HG::Clone(static_cast<HG&>(const_cast<MainSequence&>(*this)), OBJECT_PERSISTENCE::EPHEMERAL);
        double luminosityTAMS = clone->Luminosity();                                                                         // Get luminosity from clone (with updated Mass0)
        delete clone; clone = nullptr;                                                                                       // Return the memory allocated for the clone
        
        double tMS                   = m_Timescales[static_cast<int>(TIMESCALE::tMS)];
        double ageAtHookStart        = 0.99 * tMS;
        double luminosityAtHookStart = luminosity;
        
        // Linear interpolation
        luminosity = (luminosityAtHookStart * (tMS - p_Age) + luminosityTAMS * (p_Age - ageAtHookStart)) / (tMS - ageAtHookStart);
    }
    
    return luminosity;
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
    
    // If SHIKAUCHI core prescription is used, return radius that smoothly connects the beginning of MS hook and the beginning of HG
    if ((OPTIONS->MainSequenceCoreMassPrescription() == CORE_MASS_PRESCRIPTION::SHIKAUCHI) && (m_MZAMS >= 10.0)) {
        if ((m_HeliumAbundanceCore == 1.0 - m_Metallicity) && (m_TotalMassLossRate == -m_Mdot))                                 // star in MS hook and mass transfer not ongoing?
            return CalculateRadiusTransitionToHG(p_Mass, p_Time, p_RZAMS);
    }
        
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

    // pow() is slow - use multiplication where it makes sense
    double tau_3  = tau * tau * tau;
    double tau_10 = tau < FLOAT_TOLERANCE_ABSOLUTE ? 0.0: tau_3 * tau_3 * tau_3 * tau;                                          // direct comparison, to avoid underflow
    double tau_40 = tau_10 < FLOAT_TOLERANCE_ABSOLUTE ? 0.0: tau_10 * tau_10 * tau_10 * tau_10;                                 // direct comparison, to avoid underflow
    
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
 * Calculate radius on the transition from the Main Sequence to the HG when Shikauchi et al. (2024) core mass prescription is used
 *
 * SHIKAUCHI core mass prescription cannot be used beyond the MS hook (beyond age 0.99 * tMS), and this function smoothly connects
 * the radius between the beginning of the hook and the beginning of the HG
 *
 *
 * double CalculateRadiusShikauchiTransitionToHG(const double p_Mass, const double p_Age, double const p_RZAMS)
 
 * @param   [IN]    p_Mass                      Mass in Msol
 * @param   [IN]    p_Age                       Age in Myr
 * @param   [IN]    p_RZAMS                     Zero Age Main Sequence (ZAMS) Radius
 * @return                                      Radius on the Main Sequence (for age between tHook and tMS)
 */
double MainSequence::CalculateRadiusTransitionToHG(const double p_Mass, const double p_Age, double const p_RZAMS) const {
    double radiusTAMS        = CalculateRadiusAtPhaseEnd(p_Mass, p_RZAMS);
    double tMS               = m_Timescales[static_cast<int>(TIMESCALE::tMS)];
    double tauAtHookStart    = 0.99;
    double ageAtHookStart    = tauAtHookStart * tMS;
    double radiusAtHookStart = CalculateRadiusOnPhaseTau(p_Mass, tauAtHookStart);
    
    // Linear interpolation
    double radius = (radiusAtHookStart * (tMS - p_Age) + radiusTAMS * (p_Age - ageAtHookStart)) / (tMS - ageAtHookStart);
    return radius;
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
    double finalConvectiveCoreMass   = TAMSCoreMass();                                          // core mass at TAMS
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


/*
 * Calculate the convective core mass of a main sequence star that loses mass either through winds or
 * Case A mass transfer according to Shikauchi et al. (2024)
 *
 * This function also accounts for mass gain by modeling rejuvenation and updates the initial mixing core mass
 * and the helium abundance just outside the core
 *
 * double CalculateMainSequenceCoreMassShikauchi(const double p_Dt)
 *
 * @param   [IN]    p_Dt                        Current time step in Myr
 * @return                                      Tuple containing the mass of the convective core in Msol and core helium abundance
 */
DBL_DBL MainSequence::CalculateMainSequenceCoreMassShikauchi(const double p_Dt) {
    DBL_VECTOR ALPHA_COEFFICIENTS = std::get<0>(SHIKAUCHI_COEFFICIENTS);
    DBL_VECTOR FMIX_COEFFICIENTS  = std::get<1>(SHIKAUCHI_COEFFICIENTS);
    DBL_VECTOR L_COEFFICIENTS     = std::get<2>(SHIKAUCHI_COEFFICIENTS);
    
    double totalMass             = m_Mass;
    double centralHeliumFraction = m_HeliumAbundanceCore;
    double mixingCoreMass        = m_MainSequenceCoreMass;
    double lnMixingCoreMass      = std::log(mixingCoreMass);
    double heliumFractionOut     = m_HeliumAbundanceCoreOut;
        
    // Eq (A3)
    double fmix  = FMIX_COEFFICIENTS[0] + FMIX_COEFFICIENTS[1] * std::exp(-m_MZAMS / FMIX_COEFFICIENTS[2]);
    // Eq (A4)
    double beta  = 1.0 - FMIX_COEFFICIENTS[1] * totalMass / (FMIX_COEFFICIENTS[2] * fmix) * std::exp(-totalMass / FMIX_COEFFICIENTS[2]);
    // Eq (A5)
    double logL  = std::log10(CalculateLuminosityShikauchi(mixingCoreMass, centralHeliumFraction, 0.0));
    // Eq (A2)
    double alpha = PPOW(10.0, std::max(-2.0, ALPHA_COEFFICIENTS[1] * mixingCoreMass + ALPHA_COEFFICIENTS[2])) + ALPHA_COEFFICIENTS[0];
    // Eq (A7)
    double g     = -0.0044 * m_MZAMS + 0.27;
    // Eq (10)
    double Yhat  = (centralHeliumFraction - m_InitialHeliumAbundance) / (1.0 - m_InitialHeliumAbundance - m_Metallicity);
    // Eq (A6)
    double delta = std::min(PPOW(10.0, -Yhat + g), 1.0);
    
    double currentTimestepInYrs = p_Dt * MYR_TO_YEAR;
    double mDot                 = m_TotalMassLossRate;                                          // Positive for mass gain, negative for mass loss
    
    // Use boost adaptive ODE solver
    controlled_stepper_type controlled_stepper;
    state_type x(2);
    x[0] = centralHeliumFraction;
    x[1] = lnMixingCoreMass;
    auto ode = [&](const state_type &x, state_type &dxdt, const double) {
        // Eq (13)
        dxdt[0] = PPOW(10.0, logL) / (Q_CNO * mixingCoreMass);
        // Eq (12)
        dxdt[1] = - alpha / (1.0 - alpha * centralHeliumFraction) * dxdt[0] + beta * delta * mDot / totalMass;
    };
    integrate_adaptive(controlled_stepper, ode, x, 0.0, currentTimestepInYrs, currentTimestepInYrs/100.0);
    
    double newMixingCoreMass        = std::exp(x[1]);                                           // New mixing core mass
    double newCentralHeliumFraction = x[0];                                                     // New central helium fraction
    double deltaCoreMass            = newMixingCoreMass - mixingCoreMass;                       // Difference in core mass

    if (deltaCoreMass > 0.0) {                                                                  // If the core grows, we need to account for rejuvenation
        if (newMixingCoreMass < m_InitialMainSequenceCoreMass) {                                // New core mass less than initial core mass?
            // Use boost adaptive ODE solver
            state_type x(1);
            x[0] = heliumFractionOut;
            auto ode = [&](const state_type &x, state_type &dxdt, const double) {
                // Calculate the change in helium abundance just outside the core, assuming linear profile
                dxdt[0] = (heliumFractionOut - m_InitialHeliumAbundance) / (mixingCoreMass - m_InitialMainSequenceCoreMass) * (deltaCoreMass / currentTimestepInYrs);
            };
            integrate_adaptive(controlled_stepper, ode, x, 0.0, currentTimestepInYrs, currentTimestepInYrs/100.0);
        
            // Calculate the change in helium abundance in the core, assuming linear profile between Yc and Y0, and that the the accreted gas has helium fraction Y0
            double deltaY            = (heliumFractionOut - centralHeliumFraction) / (mixingCoreMass + deltaCoreMass) * deltaCoreMass + 0.5 / (mixingCoreMass + deltaCoreMass) * (heliumFractionOut - m_InitialHeliumAbundance) / (mixingCoreMass - m_InitialMainSequenceCoreMass) * deltaCoreMass * deltaCoreMass;
            newCentralHeliumFraction = centralHeliumFraction + deltaY;
            m_HeliumAbundanceCoreOut = x[0];
        }
        else {                                                                                  // New core mass greater or equal to the initial core mass?
            double firstTerm  = (heliumFractionOut - centralHeliumFraction) / (mixingCoreMass + deltaCoreMass) * deltaCoreMass;
            // Second term is set to zero if the initial core mass had been previously exceeded
            double secondTerm = (m_InitialMainSequenceCoreMass != CalculateInitialMainSequenceCoreMass(m_MZAMS)) ? 0.0 : 0.5 / (mixingCoreMass + deltaCoreMass) * (heliumFractionOut - m_InitialHeliumAbundance) / (mixingCoreMass - m_InitialMainSequenceCoreMass) * deltaCoreMass * deltaCoreMass;
            
            double deltaY                 = firstTerm + secondTerm;                             // Change in helium abundance
            newCentralHeliumFraction      = centralHeliumFraction + deltaY;
            m_HeliumAbundanceCoreOut      = m_InitialHeliumAbundance;
            m_InitialMainSequenceCoreMass = newMixingCoreMass;
        }
    }
    else
        m_HeliumAbundanceCoreOut = newCentralHeliumFraction;                                    // If core did not grow, Y_out = Y_c
    
    return std::tuple<double, double> (newMixingCoreMass, newCentralHeliumFraction);
}


/*
 * Calculate the initial core mass of a main sequence star using Equation (A3) from Shikauchi et al. (2024)
 *
 *
 * double CalculateInitialMainSequenceCoreMass(const double p_MZAMS)
 *
 * @param   [IN]    p_MZAMS                     Mass at ZAMS in Msol
 * @return                                      Mass of the convective core in Msol at ZAMS
 */
double MainSequence::CalculateInitialMainSequenceCoreMass(const double p_MZAMS) const {
    DBL_VECTOR fmixCoefficients = std::get<1>(SHIKAUCHI_COEFFICIENTS);
    double fmix = fmixCoefficients[0] + fmixCoefficients[1] * std::exp(-p_MZAMS / fmixCoefficients[2]);
    return fmix * p_MZAMS;
}


/*
 * Update the core mass of a main sequence star that loses mass through winds or Case A mass transfer
 * When Shikauchi et al. (2024) core prescription is used, also update the core helium abundance and effective age
 *
 *
 * void UpdateMainSequenceCoreMass(const double p_Dt, const double p_TotalMassLossRate)
 *
 * @param   [IN]      p_Dt                      Current timestep
 * @param   [IN]      p_TotalMassLossRate       Mass loss rate either from stellar winds or mass transfer (Msol/yr)
 */
void MainSequence::UpdateMainSequenceCoreMass(const double p_Dt, const double p_TotalMassLossRate) {
    switch (OPTIONS->MainSequenceCoreMassPrescription()) {
        case CORE_MASS_PRESCRIPTION::NONE: {
            m_MainSequenceCoreMass = 0.0;
            break;
        }
        case CORE_MASS_PRESCRIPTION::MANDEL: {                                                                                              // Calculate the minimum core mass of a main sequence star that loses mass through Case A mass transfer as the core mass of a TAMS star, scaled by the fractional age.
            if ((p_Dt != 0.0) && (p_TotalMassLossRate != -m_Mdot) && (p_TotalMassLossRate < 0.0))                                           // Only update core mass if total mass loss rate was updated in binary evolution, not applied to SSE and only applied to donors
                m_MainSequenceCoreMass = std::max(m_MainSequenceCoreMass, CalculateTauOnPhase() * TAMSCoreMass());;
            break;
        }
        case CORE_MASS_PRESCRIPTION::SHIKAUCHI: {                                                                                           // Set core mass following Shikauchi et al. (2024)
            if ((p_Dt == 0.0) || (p_TotalMassLossRate != m_TotalMassLossRate))                                                              // Do not proceed with updating the core mass if time does not advance or calculation was done in binary evolution (total mass loss rate had been updated)
                return;
            double mixingCoreMass;
            double centralHeliumFraction;
            std::tuple<double, double> ShikauchiSolution;
            double tMS = m_Timescales[static_cast<int>(TIMESCALE::tMS)];

            if ((m_MZAMS >= 10.0) && (m_HeliumAbundanceCore < 1.0 - m_Metallicity))  {                                                      // MZAMS >= 10 Msol (SHIKAUCHI prescription) and no calculations past the MS hook (Yc = 1 - Z)
                ShikauchiSolution      = CalculateMainSequenceCoreMassShikauchi(p_Dt);
                mixingCoreMass         = std::get<0>(ShikauchiSolution);
                centralHeliumFraction  = std::get<1>(ShikauchiSolution);
                m_HeliumAbundanceCore  = std::min(centralHeliumFraction, 1.0 - m_Metallicity);                                              // Update the core helium abundance
                m_MainSequenceCoreMass = mixingCoreMass;                                                                                    // Update the core mass
                m_Age                  = (m_HeliumAbundanceCore - m_InitialHeliumAbundance) / m_InitialHydrogenAbundance * 0.99 * tMS;      // Update the effective age based on central helium fraction
            }
            else if ((p_TotalMassLossRate < 0.0) && (m_MZAMS < 10.0) && (p_TotalMassLossRate != -m_Mdot)) {                                 // Mass loss and MZAMS < 10 Msol (Mandel core prescription used, which retains higher core mass after case A mass transfer -- does not take into account mass loss from winds)
                m_MainSequenceCoreMass = CalculateMainSequenceCoreMassMandel();
            }
            else { };                                                                                                                       // Other cases
            break;
        }
    }
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
    CalculateTimescales(m_Mass, m_Timescales);                                      // must update timescales
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


/*
 * Return the expected core mass at terminal age main sequence, i.e., at the start of the HG phase
 *
 * double TAMSCoreMass() const
 *
 *
 * @return                                      TAMS core Mass (Msol)
 *
 */
double MainSequence::TAMSCoreMass() const {
    // Since we are on the main sequence here, we can clone this object as an HG object
    // and, as long as it is initialised (to correctly set Tau to 0.0 on the HG phase),
    // we can query the cloned object for its core mass.
    //
    // The clone should not evolve, and so should not log anything, but to be sure the
    // clone does not participate in logging, we set its persistence to EPHEMERAL.
    
    HG *clone = HG::Clone(static_cast<HG&>(const_cast<MainSequence&>(*this)), OBJECT_PERSISTENCE::EPHEMERAL);
    double TAMSCoreMass = clone->CoreMass();                                                    // get core mass from clone
    delete clone; clone = nullptr;                                                              // return the memory allocated for the clone
    
    return TAMSCoreMass;
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

    m_Mass                 = p_Mass;
    m_Mass0                = m_Mass;
    m_MainSequenceCoreMass = 0.0;
    
    double initialHydrogenFraction = m_InitialHydrogenAbundance;
    
    CalculateTimescales();
    CalculateGBParams();
            
    m_Tau = (initialHydrogenFraction - p_HydrogenMass / m_Mass) / initialHydrogenFraction;       // assumes uniformly mixed merger product and a uniform rate of H fusion on main sequence
    
    m_Age = m_Tau * timescales(tMS);
    
    m_HeliumAbundanceCore   = 1.0 - m_Metallicity - p_HydrogenMass / p_Mass;
    
    m_HydrogenAbundanceCore = 1.0 - m_Metallicity - m_HeliumAbundanceCore;
    
    if ((OPTIONS->MainSequenceCoreMassPrescription() == CORE_MASS_PRESCRIPTION::SHIKAUCHI) && (p_Mass >= 10.0)) {
        m_InitialMainSequenceCoreMass = CalculateInitialMainSequenceCoreMass(p_Mass);                         // update initial mixing core mass
        m_MainSequenceCoreMass        = m_InitialMainSequenceCoreMass;
    }
    
    UpdateAttributesAndAgeOneTimestep(0.0, 0.0, 0.0, true);
    
    #undef timescales
}



/* 
 * Interpolate Ge+ Critical Mass Ratios, for H-rich stars
 * 
 * Function takes input QCRIT_PRESCRIPTION, currently either of the prescriptions for critical mass ratios
 * from Ge et al. (2020), GE or GE_IC. The first is the full adiabatic response, the second assumes
 * artificially isentropic envelopes. From private communication with Ge, we have an updated datatable that
 * includes qCrit for fully conservative and fully non-conservative MT, so we now interpolate on those as well.
 *
 * Interpolation is done linearly in logM, logR, and logZ
 * 
 * double BaseStar::InterpolateGeEtAlQCrit(const QCRIT_PRESCRIPTION p_qCritPrescription, const double p_massTransferEfficiencyBeta) 
 * 
 * @param   [IN]    p_qCritPrescription          Adopted critical mass ratio prescription
 * @param   [IN]    p_massTransferEfficiencyBeta Mass transfer accretion efficiency
 * @return                                       Interpolated value of either the critical mass ratio or zeta for given stellar mass / radius
 */ 
double MainSequence::InterpolateGeEtAlQCrit(const QCRIT_PRESCRIPTION p_qCritPrescription, const double p_massTransferEfficiencyBeta) {

    // Iterate over the two QCRIT_GE tables to get the qcrits at each metallicity
    double qCritPerMetallicity[2];
    std::vector<GE_QCRIT_TABLE> qCritTables = { QCRIT_GE_LOW_Z, QCRIT_GE_HIGH_Z };

    for (int ii=0; ii<2; ii++) { // iterate over the vector of tables to store qCrits per metallicity

        // Get vector of masses from qCritTable
        GE_QCRIT_TABLE &qCritTable = qCritTables[ii];
        DBL_VECTOR massesFromQCritTable = std::get<0>(qCritTable);
        GE_QCRIT_RADII_QCRIT_VECTOR radiiQCritsFromQCritTable = std::get<1>(qCritTable);

        INT_VECTOR indices = utils::BinarySearch(massesFromQCritTable, m_Mass);
        int lowerMassIndex = indices[0];
        int upperMassIndex = indices[1];
    
        if (lowerMassIndex == -1) {                                                   // if masses are out of range, set to endpoints
            lowerMassIndex = 0; 
            upperMassIndex = 1;
        } 
        else if (upperMassIndex == -1) { 
            lowerMassIndex = massesFromQCritTable.size() - 2; 
            upperMassIndex = massesFromQCritTable.size() - 1;
        } 
    
        // Get vector of radii from qCritTable for the lower and upper mass indices
        std::vector<double> logRadiusVectorLowerMass = std::get<0>(radiiQCritsFromQCritTable[lowerMassIndex]);
        std::vector<double> logRadiusVectorUpperMass = std::get<0>(radiiQCritsFromQCritTable[upperMassIndex]);
    
        // Get the qCrit vector for the lower and upper mass bounds 
        std::vector<double> qCritVectorUpperEffLowerMass;
        std::vector<double> qCritVectorUpperEffUpperMass;
        std::vector<double> qCritVectorLowerEffLowerMass;
        std::vector<double> qCritVectorLowerEffUpperMass;
        
        // Set the appropriate qCrit vector, depends on MT eff and whether you use GE STD or IC
        if (p_qCritPrescription == QCRIT_PRESCRIPTION::GE) {
            if (p_massTransferEfficiencyBeta > 0.5) {
                qCritVectorUpperEffLowerMass = std::get<1>(radiiQCritsFromQCritTable[lowerMassIndex]);
                qCritVectorUpperEffUpperMass = std::get<1>(radiiQCritsFromQCritTable[upperMassIndex]);
                qCritVectorLowerEffLowerMass = std::get<2>(radiiQCritsFromQCritTable[lowerMassIndex]);
                qCritVectorLowerEffUpperMass = std::get<2>(radiiQCritsFromQCritTable[upperMassIndex]);
            }
            else {
                qCritVectorUpperEffLowerMass = std::get<2>(radiiQCritsFromQCritTable[lowerMassIndex]);
                qCritVectorUpperEffUpperMass = std::get<2>(radiiQCritsFromQCritTable[upperMassIndex]);
                qCritVectorLowerEffLowerMass = std::get<3>(radiiQCritsFromQCritTable[lowerMassIndex]);
                qCritVectorLowerEffUpperMass = std::get<3>(radiiQCritsFromQCritTable[upperMassIndex]);

            }
        }
        else if (p_qCritPrescription == QCRIT_PRESCRIPTION::GE_IC) {
            if (p_massTransferEfficiencyBeta > 0.5) {
                qCritVectorUpperEffLowerMass = std::get<4>(radiiQCritsFromQCritTable[lowerMassIndex]);
                qCritVectorUpperEffUpperMass = std::get<4>(radiiQCritsFromQCritTable[upperMassIndex]);
                qCritVectorLowerEffLowerMass = std::get<5>(radiiQCritsFromQCritTable[lowerMassIndex]);
                qCritVectorLowerEffUpperMass = std::get<5>(radiiQCritsFromQCritTable[upperMassIndex]);
            }
            else {
                qCritVectorUpperEffLowerMass = std::get<5>(radiiQCritsFromQCritTable[lowerMassIndex]);
                qCritVectorUpperEffUpperMass = std::get<5>(radiiQCritsFromQCritTable[upperMassIndex]);
                qCritVectorLowerEffLowerMass = std::get<6>(radiiQCritsFromQCritTable[lowerMassIndex]);
                qCritVectorLowerEffUpperMass = std::get<6>(radiiQCritsFromQCritTable[upperMassIndex]);

            }
        }
    
        // Get vector of radii from qCritTable for both lower and upper masses
        INT_VECTOR indicesR0          = utils::BinarySearch(logRadiusVectorLowerMass, log10(m_Radius));
        int lowerRadiusLowerMassIndex = indicesR0[0];
        int upperRadiusLowerMassIndex = indicesR0[1];
    
        if (lowerRadiusLowerMassIndex == -1) {                                        // if radii are out of range, set to endpoints
            lowerRadiusLowerMassIndex = 0; 
            upperRadiusLowerMassIndex = 1; 
        }
        else if (upperRadiusLowerMassIndex == -1) {                                                   
            lowerRadiusLowerMassIndex = logRadiusVectorLowerMass.size() - 2; 
            upperRadiusLowerMassIndex = logRadiusVectorLowerMass.size() - 1; 
        }
    
        INT_VECTOR indicesR1          = utils::BinarySearch(logRadiusVectorUpperMass, log10(m_Radius));
        int lowerRadiusUpperMassIndex = indicesR1[0];
        int upperRadiusUpperMassIndex = indicesR1[1];
    
        if (lowerRadiusUpperMassIndex == -1) {                                        // if radii are out of range, set to endpoints
            lowerRadiusUpperMassIndex = 0; 
            upperRadiusUpperMassIndex = 1; 
        }
        else if (upperRadiusUpperMassIndex == -1) {                                                   
            lowerRadiusUpperMassIndex = logRadiusVectorUpperMass.size() - 2; 
            upperRadiusUpperMassIndex = logRadiusVectorUpperMass.size() - 1; 
        }
    
        // Set the 4 boundary points for the 2D interpolation
        double qUppLowLow = qCritVectorUpperEffLowerMass[lowerRadiusLowerMassIndex];
        double qUppLowUpp = qCritVectorUpperEffLowerMass[upperRadiusLowerMassIndex];
        double qUppUppLow = qCritVectorUpperEffUpperMass[lowerRadiusUpperMassIndex];
        double qUppUppUpp = qCritVectorUpperEffUpperMass[upperRadiusUpperMassIndex];
        double qLowLowLow = qCritVectorLowerEffLowerMass[lowerRadiusLowerMassIndex];
        double qLowLowUpp = qCritVectorLowerEffLowerMass[upperRadiusLowerMassIndex];
        double qLowUppLow = qCritVectorLowerEffUpperMass[lowerRadiusUpperMassIndex];
        double qLowUppUpp = qCritVectorLowerEffUpperMass[upperRadiusUpperMassIndex];
    
        double lowerLogRadiusLowerMass = logRadiusVectorLowerMass[lowerRadiusLowerMassIndex];
        double upperLogRadiusLowerMass = logRadiusVectorLowerMass[upperRadiusLowerMassIndex];
        double lowerLogRadiusUpperMass = logRadiusVectorUpperMass[lowerRadiusUpperMassIndex];
        double upperLogRadiusUpperMass = logRadiusVectorUpperMass[upperRadiusUpperMassIndex];

        double logLowerMass   = log10(massesFromQCritTable[lowerMassIndex]);
        double logUpperMass   = log10(massesFromQCritTable[upperMassIndex]);
    
        // Interpolate on logR first, then logM, then on the efficiency, using nearest neighbor for extrapolation
        double logRadius = log10(m_Radius);
        double qCritUpperEffLowerMass = (logRadius < lowerLogRadiusLowerMass) ? qUppLowLow
                                      : (logRadius > upperLogRadiusLowerMass) ? qUppLowUpp
                                      : qUppLowLow + (upperLogRadiusLowerMass - logRadius) / (upperLogRadiusLowerMass - lowerLogRadiusLowerMass) * (qUppLowUpp - qUppLowLow);
        double qCritUpperEffUpperMass = (logRadius < lowerLogRadiusUpperMass) ? qUppUppLow
                                      : (logRadius > upperLogRadiusUpperMass) ? qUppUppUpp
                                      : qUppUppLow + (upperLogRadiusUpperMass - logRadius) / (upperLogRadiusUpperMass - lowerLogRadiusUpperMass) * (qUppUppUpp - qUppUppLow);
        double qCritLowerEffLowerMass = (logRadius < lowerLogRadiusLowerMass) ? qLowLowLow
                                      : (logRadius > upperLogRadiusLowerMass) ? qLowLowUpp
                                      : qLowLowLow + (upperLogRadiusLowerMass - logRadius) / (upperLogRadiusLowerMass - lowerLogRadiusLowerMass) * (qLowLowUpp - qLowLowLow);
        double qCritLowerEffUpperMass = (logRadius < lowerLogRadiusUpperMass) ? qLowUppLow
                                      : (logRadius > upperLogRadiusUpperMass) ? qLowUppUpp
                                      : qLowUppLow + (upperLogRadiusUpperMass - logRadius) / (upperLogRadiusUpperMass - lowerLogRadiusUpperMass) * (qLowUppUpp - qLowUppLow);
    
        double logMass = log10(m_Mass);
        double interpolatedQCritUpperEff = (logMass < logLowerMass) ? qCritUpperEffLowerMass
                                         : (logMass > logUpperMass) ? qCritUpperEffUpperMass
                                         : qCritUpperEffLowerMass + (logUpperMass - logMass) / (logUpperMass - logLowerMass) * (qCritUpperEffUpperMass - qCritUpperEffLowerMass);
        double interpolatedQCritLowerEff = (logMass < logLowerMass) ? qCritLowerEffLowerMass
                                         : (logMass > logUpperMass) ? qCritLowerEffUpperMass
                                         : qCritLowerEffLowerMass + (logUpperMass - logMass) / (logUpperMass - logLowerMass) * (qCritLowerEffUpperMass - qCritLowerEffLowerMass);
    
        double interpolatedQCritForZ = p_massTransferEfficiencyBeta * interpolatedQCritUpperEff + (1.0 - p_massTransferEfficiencyBeta) * interpolatedQCritLowerEff;                 // Don't need to use nearest neighbor for this, beta is always between 0 and 1
        qCritPerMetallicity[ii] = interpolatedQCritForZ;
    }
    double logZlo = -3;         // log10(0.001)
    double logZhi = LOG10_ZSOL; // log10(0.02) 
    
    return qCritPerMetallicity[1] + (m_Log10Metallicity - logZhi)*(qCritPerMetallicity[1] - qCritPerMetallicity[0])/(logZhi - logZlo);
}


/*
 * Linear interpolation/extrapolation for coefficients from Shikauchi et al. (2024), used for core mass calculations
 *
 * std::tuple <DBL_VECTOR, DBL_VECTOR, DBL_VECTOR> MainSequence::InterpolateShikauchiCoefficients(const double p_Metallicity)
 *
 * @param   [IN]     p_Metallicity               Metallicity
 * @return                                       Tuple containing vectors of coefficients for the specified metallicity
 */
std::tuple <DBL_VECTOR, DBL_VECTOR, DBL_VECTOR> MainSequence::InterpolateShikauchiCoefficients(const double p_Metallicity) const {
    double logZ = std::log10(p_Metallicity);
    DBL_VECTOR alphaCoeff(3, 0.0);
    DBL_VECTOR fmixCoeff(3, 0.0);
    DBL_VECTOR lCoeff(10, 0.0);
    
    // Skip calculation if SHIKAUCHI core prescription is not used
    if (OPTIONS->MainSequenceCoreMassPrescription() != CORE_MASS_PRESCRIPTION::SHIKAUCHI)
        return std::tuple<DBL_VECTOR, DBL_VECTOR, DBL_VECTOR> (alphaCoeff, fmixCoeff, lCoeff);
    
    // Coefficients are given for these metallicities
    double low    = std::log10(0.1 * ZSOL_ASPLUND);
    double middle = std::log10(1.0/3.0 * ZSOL_ASPLUND);
    double high   = std::log10(ZSOL_ASPLUND);
    
    if (logZ <= low)                                                      // Linear extrapolation (constant) for metallicity lower than the lowest bound
        return std::tuple<DBL_VECTOR, DBL_VECTOR, DBL_VECTOR> (SHIKAUCHI_ALPHA_COEFFICIENTS[0], SHIKAUCHI_FMIX_COEFFICIENTS[0], SHIKAUCHI_L_COEFFICIENTS[0]);
    
    else if ((logZ > low) && (logZ <= middle)) {                          // Linear interpolation between metallicity low and middle
        for (int i = 0; i < 3; i++)
            alphaCoeff[i] = (SHIKAUCHI_ALPHA_COEFFICIENTS[0][i] * (middle - logZ) + SHIKAUCHI_ALPHA_COEFFICIENTS[1][i] * (logZ - low)) / (middle - low);
        for (int i = 0; i < 3; i++)
            fmixCoeff[i]  = (SHIKAUCHI_FMIX_COEFFICIENTS[0][i] * (middle - logZ) + SHIKAUCHI_FMIX_COEFFICIENTS[1][i] * (logZ - low)) / (middle - low);
        for (int i = 0; i < 10; i++)
            lCoeff[i]     = (SHIKAUCHI_L_COEFFICIENTS[0][i] * (middle - logZ) + SHIKAUCHI_L_COEFFICIENTS[1][i] * (logZ - low)) / (middle - low);
        return std::tuple<DBL_VECTOR, DBL_VECTOR, DBL_VECTOR> (alphaCoeff, fmixCoeff, lCoeff);
    }
    
    else if ((logZ > middle) && (logZ < high)) {                          // Linear interpolation between metallicity middle and high
        for (int i = 0; i < 3; i++)
            alphaCoeff[i] = (SHIKAUCHI_ALPHA_COEFFICIENTS[1][i] * (high - logZ) + SHIKAUCHI_ALPHA_COEFFICIENTS[2][i] * (logZ - middle)) / (high - middle);
        for (int i = 0; i < 3; i++)
            fmixCoeff[i]  = (SHIKAUCHI_FMIX_COEFFICIENTS[1][i] * (high - logZ) + SHIKAUCHI_FMIX_COEFFICIENTS[2][i] * (logZ - middle)) / (high - middle);
        for (int i = 0; i < 10; i++)
            lCoeff[i]     = (SHIKAUCHI_L_COEFFICIENTS[1][i] * (high - logZ) + SHIKAUCHI_L_COEFFICIENTS[2][i] * (logZ - middle)) / (high - middle);
        return std::tuple<DBL_VECTOR, DBL_VECTOR, DBL_VECTOR> (alphaCoeff, fmixCoeff, lCoeff);
    }
    
    else                                                                  // Linear extrapolation (constant) for metallicity equal to solar or higher
        return std::tuple<DBL_VECTOR, DBL_VECTOR, DBL_VECTOR> (SHIKAUCHI_ALPHA_COEFFICIENTS[2], SHIKAUCHI_FMIX_COEFFICIENTS[2], SHIKAUCHI_L_COEFFICIENTS[2]);
}
