#include "GiantBranch.h"
#include "HeMS.h"
#include "BH.h"


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                     COEFFICIENT AND CONSTANT CALCULATIONS ETC.                    //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * Calculate the Hydrogen rate constant AH'
 *
 * Hurley et al. 2000, just after eq 43 (before eq 44)
 *
 *
 * double CalculateHRateConstant_Static(const double p_Mass)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @return                                      Hydrogen rate constant (Msol Lsol^-1 Myr^-1)
 */
double GiantBranch::CalculateHRateConstant_Static(const double p_Mass) {
    return pow(10.0, std::max(-4.8, std::min((-5.7 + (0.8 * p_Mass)), (-4.1 + (0.14 * p_Mass)))));
}


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
void GiantBranch::CalculateTimescales(const double p_Mass, DBL_VECTOR &p_Timescales) {
#define timescales(x) p_Timescales[static_cast<int>(TIMESCALE::x)]      // for convenience and readability - undefined at end of function
#define gbParams(x) m_GBParams[static_cast<int>(GBP::x)]                // for convenience and readability - undefined at end of function

    double p1   = gbParams(p) - 1.0;
    double q1   = gbParams(q) - 1.0;
    double p1_p = p1 / gbParams(p);
    double q1_q = q1 / gbParams(q);

    double LBGB = CalculateLuminosityAtPhaseBase_Static(p_Mass, m_AnCoefficients);

    MainSequence::CalculateTimescales(p_Mass, p_Timescales);   // calculate common values

    timescales(tinf1_FGB) = timescales(tBGB) + ((1.0 / (p1 * gbParams(AH) * gbParams(D))) * pow((gbParams(D) / LBGB), p1_p));
    timescales(tMx_FGB)   = timescales(tinf1_FGB) - ((timescales(tinf1_FGB) - timescales(tBGB)) * pow((LBGB / gbParams(Lx)), p1_p));
    timescales(tinf2_FGB) = timescales(tMx_FGB) + ((1.0 / (q1 * gbParams(AH) * gbParams(B))) * pow((gbParams(B) / gbParams(Lx)), q1_q));

    timescales(tHeI)      = CalculateLifetimeToHeIgnition(p_Mass, timescales(tinf1_FGB), timescales(tinf2_FGB));
    timescales(tHeMS)     = HeMS::CalculateLifetimeOnPhase_Static(p_Mass);

#undef gbParams
#undef timescales
}


/*
 * Calculate the Core mass - Luminosity relation parameter B
 *
 * Hurley et al. 2000, eqs 31 - 38
 *
 *
 * double CalculateCoreMass_Luminosity_B_Static(const double p_Mass)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @return                                      Core mass - Luminosity relation parameter B
 */
double GiantBranch::CalculateCoreMass_Luminosity_B_Static(const double p_Mass) {
    return std::max(3.0E4, (500.0 + (1.75E4 * pow(p_Mass, 0.6))));
}


/*
 * Calculate the Core mass - Luminosity relation parameter D
 *
 * Hurley et al. 2000, eqs 31 - 38
 *
 *
 * double CalculateCoreMass_Luminosity_D_Static(const double p_Mass, const double p_LogMetallicityXi, DBL_VECTOR &p_MassCutoffs)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @param   [IN]    p_LogMetallicityXi          log10(Metallicity / Zsol) - called xi in Hurley et al 2000
 * @param   [IN]    p_MassCutoffs               Mass cutoffs vector
 * @return                                      Core mass - Luminosity relation parameter D
 */
double GiantBranch::CalculateCoreMass_Luminosity_D_Static(const double p_Mass, const double p_LogMetallicityXi, const DBL_VECTOR &p_MassCutoffs) {
#define massCutoffs(x) p_MassCutoffs[static_cast<int>(MASS_CUTOFF::x)]  // for convenience and readability - undefined at end of function

    double D0   = 5.37 + (0.135 * p_LogMetallicityXi);
    double D1   = (0.975 * D0) - (0.18 * p_Mass);
    double logD = D0;

    if (utils::Compare(p_Mass, massCutoffs(MHeF)) > 0) {
        if (utils::Compare(p_Mass, 2.5) >= 0) {
            double D2 = (0.5 * D0) - (0.06 * p_Mass);
            logD      = std::max(std::max(-1.0, D1), D2);
        }
        else {  // Linear interpolation between end points
            double gradient  = (D0 - D1) / (massCutoffs(MHeF) - 2.5);
            double intercept = D0 - (massCutoffs(MHeF) * gradient);
            logD = (gradient * p_Mass) + intercept;
        }
    }

    return pow(10.0, logD);

#undef massCutoffs
}


/*
 * Calculate the Core mass - Luminosity relation parameter p
 *
 * Hurley et al. 2000, eqs 31 - 38
 *
 *
 * double CalculateCoreMass_Luminosity_p_Static(const double p_Mass, DBL_VECTOR &p_MassCutoffs)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @param   [IN]    p_MassCutoffs               Mass cutoffs vector
 * @return                                      Core mass - Luminosity relation parameter p
 */
double GiantBranch::CalculateCoreMass_Luminosity_p_Static(const double p_Mass, const DBL_VECTOR &p_MassCutoffs) {
#define massCutoffs(x) p_MassCutoffs[static_cast<int>(MASS_CUTOFF::x)]  // for convenience and readability - undefined at end of function

    double p = 6.0;

    if (utils::Compare(p_Mass, massCutoffs(MHeF)) > 0) {
        if (utils::Compare(p_Mass, 2.5) >= 0) {
            p = 5.0;
        }
        else {                                                          // linear interpolation between end points
            double gradient  = 1.0 / (massCutoffs(MHeF) - 2.5);         // will be negative
            double intercept = 5.0 - (2.5 * gradient);
            p = (gradient * p_Mass) + intercept;
        }
    }

    return p;

#undef massCutoffs
}


/*
 * Calculate the Core mass - Luminosity relation parameter q
 *
 * Hurley et al. 2000, eqs 31 - 38
 *
 *
 * double CalculateCoreMass_Luminosity_q_Static(const double p_Mass, DBL_VECTOR &p_MassCutoffs)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @param   [IN]    p_MassCutoffs               Mass cutoffs vector
 * @return                                      Core mass - Luminosity relation parameter q
 */
double GiantBranch::CalculateCoreMass_Luminosity_q_Static(const double p_Mass, const DBL_VECTOR &p_MassCutoffs) {
#define massCutoffs(x) p_MassCutoffs[static_cast<int>(MASS_CUTOFF::x)]  // for convenience and readability - undefined at end of function

    double q = 3.0;

    if (utils::Compare(p_Mass, massCutoffs(MHeF)) > 0) {
        if (utils::Compare(p_Mass, 2.5) >= 0) {
            q = 2.0;
        }
        else {                                                          // Linear interpolation between end points
            double gradient  = 1.0 / (massCutoffs(MHeF) - 2.5);
            double intercept = 2.0 - (2.5 * gradient);
            q = (gradient * p_Mass) + intercept;
        }
    }

    return q;

#undef massCutoffs
}


/*
 * Calculate the Core mass - Luminosity relation parameter Mx
 *
 * Mx is the point at which the low- and high-luminosity approximations cross
 *
 * Hurley et al. 2000, eq 38
 *
 *
 * double CalculateCoreMass_Luminosity_Mx_Static(const DBL_VECTOR &p_GBParams)
 *
 * @param   [IN]    p_GBParams                  Giant Branch Parameters
 * @return                                      Core mass - Luminosity relation parameter Mx
 */
double GiantBranch::CalculateCoreMass_Luminosity_Mx_Static(const DBL_VECTOR &p_GBParams) {
#define gbParams(x) p_GBParams[static_cast<int>(GBP::x)]    // for convenience and readability - undefined at end of function

    return pow((gbParams(B) / gbParams(D)), (1.0 / (gbParams(p) - gbParams(q))));

#undef gbParams
}


/*
 * Calculate First Giant Branch (FGB) Core mass - Luminosity relation parameter Lx
 *
 * Hurley et al. 2000, eq 37
 *
 *
 * double CalculateCoreMass_Luminosity_Lx_Static(const DBL_VECTOR &p_GBParams)
 *
 * @param   [IN]    p_GBParams                  Giant Branch Parameters
 * @return                                      Core mass - Luminosity relation parameter Lx
 */
double GiantBranch::CalculateCoreMass_Luminosity_Lx_Static(const DBL_VECTOR &p_GBParams) {
#define gbParams(x) p_GBParams[static_cast<int>(GBP::x)]    // for convenience and readability - undefined at end of function
    // since the mass used here is the mass at crossover (Mx), these
    // should give the same answer - but we'll take the minimum anyway
    return std::min((gbParams(B) * pow(gbParams(Mx), gbParams(q))), (gbParams(D) * pow(gbParams(Mx), gbParams(p))));

#undef gbParams
}


/*
 * Calculate Giant Branch (GB) parameters per Hurley et al. 2000
 *
 * Giant Branch Parameters depend on a star's mass, so this needs to be called at least each timestep
 *
 * Vectors are passed by reference here for performance - preference would be to pass const& and
 * pass modified value back by functional return, but this way is faster - and this function is
 * called many, many times.
 *
 * void CalculateGBParams(const double p_Mass, DBL_VECTOR &p_GBParams)
 *
 * @param   [IN]        p_Mass                  Mass in Msol
 * @param   [IN/OUT]    p_GBParams              Giant Branch Parameters - calculated here
 */
void GiantBranch::CalculateGBParams(const double p_Mass, DBL_VECTOR &p_GBParams) {
#define gbParams(x) p_GBParams[static_cast<int>(GBP::x)]    // for convenience and readability - undefined at end of function

    gbParams(AH)     = CalculateHRateConstant_Static(p_Mass);
    gbParams(AHHe)   = CalculateHHeRateConstant_Static();
    gbParams(AHe)    = CalculateHeRateConstant_Static();

    gbParams(B)      = CalculateCoreMass_Luminosity_B_Static(p_Mass);
    gbParams(D)      = CalculateCoreMass_Luminosity_D_Static(p_Mass, m_LogMetallicityXi, m_MassCutoffs);

    gbParams(p)      = CalculateCoreMass_Luminosity_p_Static(p_Mass, m_MassCutoffs);
    gbParams(q)      = CalculateCoreMass_Luminosity_q_Static(p_Mass, m_MassCutoffs);

    gbParams(Mx)     = CalculateCoreMass_Luminosity_Mx_Static(p_GBParams);      // depends on B, D, p & q - recalculate if any of those are changed
    gbParams(Lx)     = CalculateCoreMass_Luminosity_Lx_Static(p_GBParams);      // JR: Added this - depends on B, D, p, q & Mx - recalculate if any of those are changed

    gbParams(McDU)   = CalculateCoreMassAt2ndDredgeUp_Static(gbParams(McBAGB));
    gbParams(McBAGB) = CalculateCoreMassAtBAGB(p_Mass);
    gbParams(McBGB)  = CalculateCoreMassAtBGB(p_Mass, p_GBParams);

    gbParams(McSN)   = CalculateCoreMassAtSupernova_Static(gbParams(McBAGB));   // JR: Added this

#undef gbParams
}


/*
 * Calculate Giant Branch (GB) parameters per Hurley et al. 2000
 *
 * Giant Branch Parameters depend on a star's mass, so this needs to be called at least each timestep
 *
 * Vectors are passed by reference here for performance - preference would be to pass const& and
 * pass modified value back by functional return, but this way is faster - and this function is
 * called many, many times.
 *
 *
 * This function exists to facilitate the calculation of gbParams in EAGB::ResolveEnvelopeLoss() for the
 * stellar type to which the star will evolve.  The calculations of some of the stellar attributes there
 * depend on new gbParams.  
 * JR: This really needs to be revisited one day - these calculations should really be performed after
 *     switching to the new stellar type, but other calculations are done (in the legacy code) before the 
 *     switch (see evolveOneTimestep() in star.cpp for EAGB stars in the legacy code).
 *
 *
 * void CalculateGBParams_Static(const double      p_Mass, 
 *                               const double      p_LogMetallicityXi, 
 *                               const DBL_VECTOR &p_MassCutoffs, 
 *                               const DBL_VECTOR &p_AnCoefficients, 
 *                               const DBL_VECTOR &p_BnCoefficients,* 
 *                                     DBL_VECTOR &p_GBParams)
 *
 * @param   [IN]        p_Mass                  Mass in Msol
 * @param   [IN]        p_LogMetallicityXi      log10(Metallicity / Zsol) - called xi in Hurley et al 2000
 * @param   [IN]        p_MassCutoffs           Mass cutoffs
 * @param   [IN]        p_AnCoefficients        a(n) coefficients
 * @param   [IN]        p_BnCoefficients        b(n) coefficients
 * @param   [IN/OUT]    p_GBParams              Giant Branch Parameters - calculated here
 */
void GiantBranch::CalculateGBParams_Static(const double      p_Mass, 
                                           const double      p_LogMetallicityXi, 
                                           const DBL_VECTOR &p_MassCutoffs, 
                                           const DBL_VECTOR &p_AnCoefficients, 
                                           const DBL_VECTOR &p_BnCoefficients, 
                                                 DBL_VECTOR &p_GBParams) {

#define gbParams(x) p_GBParams[static_cast<int>(GBP::x)]    // for convenience and readability - undefined at end of function

    gbParams(AH)     = CalculateHRateConstant_Static(p_Mass);
    gbParams(AHHe)   = CalculateHHeRateConstant_Static();
    gbParams(AHe)    = CalculateHeRateConstant_Static();

    gbParams(B)      = CalculateCoreMass_Luminosity_B_Static(p_Mass);
    gbParams(D)      = CalculateCoreMass_Luminosity_D_Static(p_Mass, p_LogMetallicityXi, p_MassCutoffs);

    gbParams(p)      = CalculateCoreMass_Luminosity_p_Static(p_Mass, p_MassCutoffs);
    gbParams(q)      = CalculateCoreMass_Luminosity_q_Static(p_Mass, p_MassCutoffs);

    gbParams(Mx)     = CalculateCoreMass_Luminosity_Mx_Static(p_GBParams);      // depends on B, D, p & q - recalculate if any of those are changed
    gbParams(Lx)     = CalculateCoreMass_Luminosity_Lx_Static(p_GBParams);      // JR: Added this - depends on B, D, p, q & Mx - recalculate if any of those are changed

    gbParams(McDU)   = CalculateCoreMassAt2ndDredgeUp_Static(gbParams(McBAGB));
    gbParams(McBAGB) = CalculateCoreMassAtBAGB_Static(p_Mass, p_BnCoefficients);
    gbParams(McBGB)  = CalculateCoreMassAtBGB_Static(p_Mass, p_MassCutoffs, p_AnCoefficients, p_GBParams);

    gbParams(McSN)   = CalculateCoreMassAtSupernova_Static(gbParams(McBAGB));   // JR: Added this

#undef gbParams
}


/*
 * Calculate the perturbation parameter mu
 *
 * Hurley et al. 2000, eqs 97 & 98
 *
 *
 * double CalculatePerturbationMu()
 *
 * @return                                      Perturbation parameter mu
 */
double GiantBranch::CalculatePerturbationMu() {
    double kappa = -0.5;
    double L0    = 7.0E4;

    return ((m_Mass - m_CoreMass) / m_Mass) * (std::min(5.0, std::max(1.2, pow((m_Luminosity / L0), kappa))));
}


/*
 * Perturb Luminosity and Radius
 *
 * Perturbs Luminosity and Radius per Hurley et al. 2000, section 6.3
 * The attributes of the star are updated.
 *
 * Perturbation is disabled by default when DEBUG is enabled - except when
 * DEBUG_PERTURB is defined (see below).  The stellar class perturbation
 * functions are defined away if DEBUG is defined - so the generic Star
 * function is called (and does nothing). (So far FGB is the only class that
 * defines this function where it actually does anything)
 *
 * If DEBUG_PERTURB is defined then perturbation is not disabled while debbuging.
 * To enable perturbation while DEBUG is enabled, define DEBUG_PERTURB.
 *
 *
 * void PerturbLuminosityAndRadius()
 */
#if !defined(DEBUG) || defined(DEBUG_PERTURB)   // don't perturb if debugging and DEBUG_PERTURB not defined
void GiantBranch::PerturbLuminosityAndRadius() {

    if (utils::Compare(m_Mu, 1.0) < 0) {   // perturb only if mu < 1.0

        double Lc = CalculateRemnantLuminosity();
        double Rc = CalculateRemnantRadius();

        double s = CalculatePerturbationS(m_Mu, m_Mass);
        double r = CalculatePerturbationR(m_Mu, m_Mass, m_Radius, Rc);

        m_Luminosity = Lc * pow((m_Luminosity / Lc), s);        
        m_Radius     = Rc * pow((m_Radius / Rc), r);
    }
}
#else
void GiantBranch::PerturbLuminosityAndRadius() { }
#endif


/*
 * Calculate the mass-radius response exponent Zeta
 *
 * Hurley et al. 2000, eqs 97 & 98
 *
 *
 * double CalculateZeta(CE_ZETA_PRESCRIPTION p_CEZetaPrescription)
 *
 * @return                                      mass-radius response exponent Zeta
 */
double GiantBranch::CalculateZeta(CE_ZETA_PRESCRIPTION p_CEZetaPrescription) {

    double zeta = 0.0;                                              // default value

    switch (p_CEZetaPrescription) {                                 // which prescription?
        case CE_ZETA_PRESCRIPTION::SOBERMAN:                        // SOBERMAN: Soberman, Phinney, and van den Heuvel, 1997, eq 61
            zeta = CalculateZadiabaticSPH(m_HeCoreMass);
            break;

        case CE_ZETA_PRESCRIPTION::HURLEY:                          // HURLEY: Hurley, Tout, and Pols, 2002, eq 56
            zeta = CalculateZadiabaticHurley2002(m_HeCoreMass);
            break;

        case CE_ZETA_PRESCRIPTION::ARBITRARY:                       // ARBITRARY: user program options thermal zeta value
            zeta = OPTIONS->ZetaThermalArbitrary();
            break;

        case CE_ZETA_PRESCRIPTION::STARTRACK:                       // no longer implemented - JR: todo: remove this from program options (then remove from here)?  Currently it is the default option...
            m_Error = ERROR::UNSUPPORTED_CE_ZETA_PRESCRIPTION;      // set error value
            SHOW_ERROR(m_Error);                                    // warn that an error occurred
            break;

        default:                                                    // unknown common envelope prescription - shouldn't happen
            m_Error = ERROR::UNKNOWN_CE_ZETA_PRESCRIPTION;          // set error value
            SHOW_ERROR(m_Error);                                    // warn that an error occurred
    }

    return zeta;
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                              LUMINOSITY CALCULATIONS                              //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * Calculate luminosity at the Base of the Giant Branch
 *
 * Hurley et al. 2000, eq 10
 *
 *
 * double CalculateLuminosityAtPhaseBase_Static(const double p_Mass, const DBL_VECTOR &p_AnCoefficients)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @param   [IN]    p_AnCoefficients            a(n) coefficients
 * @return                                      Luminosity at the Base of the Giant Branch in Lsol
 *
 * p_AnCoefficients passed as parameter so function can be declared static
 */
double GiantBranch::CalculateLuminosityAtPhaseBase_Static(const double p_Mass, const DBL_VECTOR &p_AnCoefficients) {
#define a p_AnCoefficients  // for convenience and readability - undefined at end of function

    double top    = (a[27] * pow(p_Mass, a[31])) + (a[28] * pow(p_Mass, C_COEFF.at(2)));
    double bottom = a[29] + (a[30] * pow(p_Mass, C_COEFF.at(3))) + pow(p_Mass, a[32]);

    return top / bottom;

#undef a
}


/*
 * Calculate luminosity on the Zero Age Horizontal Branch
 * (for Low Mass stars)
 *
 * Hurley et al. 2000, eq 53
 *
 *
 * double CalculateLuminosityOnZAHB_Static(const double      p_Mass,
 *                                         const double      p_CoreMass,
 *                                         const double      p_Alpha1,
 *                                         const double      p_MHeF,
 *                                         const double      p_MFGB,
 *                                         const double      p_MinimumLuminosityOnPhase,
 *                                         const DBL_VECTOR &p_BnCoefficients)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @param   [IN]    p_CoreMass                  Core Mass in Msol
 * @param   [IN]    p_Alpha1                    alpha1 in Hurly et al. 2000 (just after eq 49)
 * @param   [IN]    p_MHeF                      Maximum initial mass for which helium ignites degenerately in a Helium Flash
 * @param   [IN]    p_MFGB                      Maximum initial mass for which helium ignites on the First Giant Branch
 * @param   [IN]    p_BnCoefficients            b(n) coefficients
 * @return                                      Luminosity on the Zero Age Horizontal Branch in Lsol
 *
 * p_Alpha1, p_MHeF, p_MFGB and p_BnCoefficients passed as parameters so function can be declared static
 */
double GiantBranch::CalculateLuminosityOnZAHB_Static(const double      p_Mass,
                                                     const double      p_CoreMass,
                                                     const double      p_Alpha1,
                                                     const double      p_MHeF,
                                                     const double      p_MFGB,
                                                     const double      p_MinimumLuminosityOnPhase,
                                                     const DBL_VECTOR &p_BnCoefficients) {
#define b p_BnCoefficients  // for convenience and readability - undefined at end of function

    double LZHe   = HeMS::CalculateLuminosityAtZAMS_Static(p_CoreMass);
    double LminHe = p_MinimumLuminosityOnPhase;

    double mu     = (p_Mass - p_CoreMass) / (p_MHeF - p_CoreMass);
    double alpha2 = (b[18] + LZHe - LminHe) / (LminHe - LZHe);

    double first  = (1.0 + b[20]) / (1.0 + (b[20] * pow(mu, 1.6479)));
    double second = (b[18] * pow(mu, b[19])) / (1.0 + alpha2 * exp(15.0 * (p_Mass - p_MHeF)));

    return LZHe + (first * second);

#undef b
}


/*
 * Calculate luminosity at Helium Ignition
 *
 * Hurley et al. 2000, eq 49
 *
 *
 * double CalculateLuminosityAtHeIgnition_Static(const double      p_Mass,
 *                                               const double      p_Alpha1,
 *                                               const double      p_MHeF,
 *                                               const DBL_VECTOR &p_BnCoefficients)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @param   [IN]    p_Alpha1                    alpha1 in Hurley et al. 2000 (just after eq 49)
 * @param   [IN]    p_MHeF                      Maximum initial mass for which helium ignites degenerately in a Helium Flash
 * @param   [IN]    p_BnCoefficients            b(n) coefficients
 * @return                                      Luminosity at Helium Ignition in Lsol
 *
 * p_Alpha1, p_MHeF and p_BnCoefficients passed as parameter so function can be declared static
 */
double GiantBranch::CalculateLuminosityAtHeIgnition_Static(const double      p_Mass,
                                                           const double      p_Alpha1,
                                                           const double      p_MHeF,
                                                           const DBL_VECTOR &p_BnCoefficients) {
#define b p_BnCoefficients  // for convenience and readability - undefined at end of function

    return (utils::Compare(p_Mass, p_MHeF) < 0)
            ? (b[9] * pow(p_Mass, b[10])) / (1.0 + (p_Alpha1 * exp(15.0 * (p_Mass - p_MHeF))))
            : (b[11] + (b[12] * pow(p_Mass, 3.8))) / (b[13] + (p_Mass * p_Mass));

#undef b
}


/*
 * Calculate luminosity of the remnant the star would become if it lost all of its
 * envelope immediately (i.e. M = Mc, coreMass)
 *
 * Hurley et al. 2000, just after eq 105
 *
 *
 * double CalculateRemnantLuminosity()
 *
 * @return                                      Luminosity of remnant core in Lsol
 */
double GiantBranch::CalculateRemnantLuminosity() {
#define massCutoffs(x) m_MassCutoffs[static_cast<int>(MASS_CUTOFF::x)]  // for convenience and readability - undefined at end of function

    return (utils::Compare(m_Mass, massCutoffs(MHeF)) < 0)
            ? HeMS::CalculateLuminosityAtZAMS_Static(m_CoreMass)
            : HeWD::CalculateLuminosityOnPhase_Static(m_CoreMass, 0.0, m_Metallicity);

#undef massCutoffs
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                RADIUS CALCULATIONS                                //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * Calculate radius on the Giant Branch
 *
 * Hurley et al. 2000, eq 46
 *
 *
 * double CalculateRadiusOnPhase_Static(const double p_Mass, const double p_Luminosity, const DBL_VECTOR &p_BnCoefficients)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @param   [IN]    p_Luminosity                Luminosity in Lsol
 * @param   [IN]    p_BnCoefficients            b(n) coefficients
 * @return                                      Radius on the First Giant Branch in Rsol
 *
 * p_BnCoefficients passed as parameter so function can be declared static
 */
double GiantBranch::CalculateRadiusOnPhase_Static(const double p_Mass, const double p_Luminosity, const DBL_VECTOR &p_BnCoefficients) {
#define b p_BnCoefficients  // for convenience and readability - undefined at end of function

    double A = std::min((b[4] * pow(p_Mass, -b[5])), (b[6] * pow(p_Mass, -b[7])));  // Hurley et al. 2000, just before eq 47

    return A * (pow(p_Luminosity, b[1]) + (b[2] * pow(p_Luminosity, b[3])));        // Hurley et al. 2000, eq 46

#undef b
}


/*
 * Calculate radius on the Zero Age Horizontal Branch
 *
 * Hurley et al. 2000, eq 54
 *
 *
 * double CalculateRadiusOnZAHB_Static(const double      p_Mass,
 *                                     const double      p_CoreMass,
 *                                     const double      p_Alpha1,
 *                                     const double      p_MHeF,
 *                                     const double      p_MFGB,
 *                                     const double      p_MinimumLuminosityOnPhase,
 *                                     const DBL_VECTOR &p_BnCoefficients)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @param   [IN]    p_CoreMass                  Core Mass in Msol
 * @param   [IN]    p_Alpha1                    alpha1 in Hurly et al. 2000 (just after eq 49)
 * @param   [IN]    p_MHeF                      Maximum initial mass for which helium ignites degenerately in a Helium Flash
 * @param   [IN]    p_MFGB                      Maximum initial mass for which helium ignites don the First Giant Branch
 * @param   [IN]    p_MinimumLuminosityOnPhase  Minimum luminosity on phase (only required for CHeB stars) - calculated once per star
 * @param   [IN]    p_BnCoefficients            b(n) coefficients
 * @return                                      Radius on the Zero Age Horizontal Branch in Rsol
 *
 * p_Alpha1, p_MHeF and p_BnCoefficients passed as parameters so function can be declared static
 */
double GiantBranch::CalculateRadiusOnZAHB_Static(const double      p_Mass,
                                                 const double      p_CoreMass,
                                                 const double      p_Alpha1,
                                                 const double      p_MHeF,
                                                 const double      p_MFGB,
                                                 const double      p_MinimumLuminosityOnPhase,
                                                 const DBL_VECTOR &p_BnCoefficients) {
#define b p_BnCoefficients  // for convenience and readability - undefined at end of function

    double RZHe  = HeMS::CalculateRadiusAtZAMS_Static(p_CoreMass);
    double LZAHB = CalculateLuminosityOnZAHB_Static(p_Mass, p_CoreMass, p_Alpha1, p_MHeF, p_MFGB, p_MinimumLuminosityOnPhase, p_BnCoefficients);
    double RGB   = GiantBranch::CalculateRadiusOnPhase_Static(p_Mass, LZAHB, p_BnCoefficients);

    double mu    = (p_Mass - p_CoreMass) / (p_MHeF - p_CoreMass);
    double f     = ((1.0 + b[21]) * pow(mu, b[22])) / (1.0 + b[21] * pow(mu, b[23]));

    return ((1.0 - f)) * RZHe + (f * RGB);

#undef b
}


/*
 * Calculate radius at Helium Ignition
 *
 * Hurley et al. 2000, eq 50
 *
 *
 * double CalculateRadiusAtHeIgnition(const double p_Mass)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @return                                      Radius at Helium Ignition in Rsol
 */
double GiantBranch::CalculateRadiusAtHeIgnition(const double p_Mass) {
#define massCutoffs(x) m_MassCutoffs[static_cast<int>(MASS_CUTOFF::x)]  // for convenience and readability - undefined at end of function

    double RHeI = 0.0;                                                  // Radius at Helium Ignition

    double LHeI      = CalculateLuminosityAtHeIgnition_Static(p_Mass, m_Alpha1, massCutoffs(MHeF), m_BnCoefficients);
    double RmHe      = CHeB::CalculateMinimumRadiusOnPhase_Static(p_Mass, m_CoreMass, m_Alpha1, massCutoffs(MHeF), massCutoffs(MFGB), m_MinimumLuminosityOnPhase, m_BnCoefficients);
    double RGB_LHeI  = CalculateRadiusOnPhase(p_Mass, LHeI);

    if (utils::Compare(p_Mass, massCutoffs(MFGB)) <= 0) {
        RHeI = RGB_LHeI;
    }
    else if (utils::Compare(p_Mass, std::max(massCutoffs(MFGB), 12.0)) >= 0) {
        double RAGB_LHeI = EAGB::CalculateRadiusOnPhase_Static(p_Mass, LHeI, massCutoffs(MHeF), m_BnCoefficients);
        RHeI             = std::min(RmHe, RAGB_LHeI);                        // Hurley et al. 2000, eq 55
    }
    else {
        double mu = log10(p_Mass / 12.0) / log10(massCutoffs(MFGB) / 12.0);
        RHeI      = RmHe * pow((RGB_LHeI / RmHe), mu);
    }

    return RHeI;

#undef massCutoffs
}


/*
 * Calculate radius of the remnant the star would become if it lost all of its
 * envelope immediately (i.e. M = Mc, coreMass)
 *
 * Hurley et al. 2000, just after eq 105
 *
 *
 * double CalculateRemnantRadius()
 *
 * @return                                      Radius of remnant core in Rsol
 */
double GiantBranch::CalculateRemnantRadius() {
#define massCutoffs(x) m_MassCutoffs[static_cast<int>(MASS_CUTOFF::x)]  // for convenience and readability - undefined at end of function

    return (utils::Compare(m_Mass, massCutoffs(MHeF)) < 0)
            ? HeMS::CalculateRadiusAtZAMS_Static(m_CoreMass)
            : HeWD::CalculateRadiusOnPhase_Static(m_CoreMass);

#undef massCutoffs
}


/*
 * Calculate the radial extent of the star's convective envelope (if it has one)
 *
 * Hurley et al. 2000, sec. 2.3, particularly subsec. 2.3.1, eqs 36-40
 *
 * (Technically not a radius calculation I suppose, but "radial extent" is close enough to put it with the radius calculations...)
 *
 *
 * double CalculateRadialExtentConvectiveEnvelope()
 *
 * @return                                      Radial extent of the star's convective envelope in Rsol
 */
double GiantBranch::CalculateRadialExtentConvectiveEnvelope() {

	BaseStar clone = *this;                         // clone this star so can manipulate without changes persisiting
	clone.ResolveRemnantAfterEnvelopeLoss();        // update clone's attributes after envelope is lost

    return m_Radius - clone.Radius();
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                 MASS CALCULATIONS                                 //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * Calculate core mass at the Base of the Asymptotic Giant Branch
 *
 * Hurley et al. 2000, eq 66
 *
 *
 * double CalculateCoreMassAtBAGB(const double p_Mass)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @return                                      Core mass at the Base of the Asymptotic Giant Branch in Msol
 */
double GiantBranch::CalculateCoreMassAtBAGB(const double p_Mass) {
#define b m_BnCoefficients  // for convenience and readability - undefined at end of function

    return sqrt(sqrt((b[36] * pow(p_Mass, b[37])) + b[38]));   // sqrt() is much faster than pow()

#undef b
}


/*
 * Calculate core mass at the Base of the Asymptotic Giant Branch
 *
 * Hurley et al. 2000, eq 66
 *
 * Static version required by CalculateGBParams_Static()
 * 
 *
 * double CalculateCoreMassAtBAGB_Static(const double p_Mass, const DBL_VECTOR &p_BnCoefficients)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @param   [IN]    p_BnCoefficients            b(n) coefficients
 * @return                                      Core mass at the Base of the Asymptotic Giant Branch in Msol
 */
double GiantBranch::CalculateCoreMassAtBAGB_Static(const double p_Mass, const DBL_VECTOR &p_BnCoefficients) {
#define b p_BnCoefficients  // for convenience and readability - undefined at end of function

    return sqrt(sqrt((b[36] * pow(p_Mass, b[37])) + b[38]));   // sqrt() is much faster than pow()

#undef b
}


/*
 * Calculate core mass at the Base of the Giant Branch
 *
 * Hurley et al. 2000, eq 44
 *
 * For large enough M, we have McBGB ~ 0.098*Mass**(1.35)
 *
 *
 * double CalculateCoreMassAtBGB(const double p_Mass, const DBL_VECTOR &p_GBParams)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @param   [IN]    p_GBParams                  Giant Branch paramaters
 * @return                                      Core mass at the Base of the Giant Branch in Msol
 */
double GiantBranch::CalculateCoreMassAtBGB(const double p_Mass, const DBL_VECTOR &p_GBParams) {
#define gbParams(x) p_GBParams[static_cast<int>(GBP::x)]                // for convenience and readability - undefined at end of function
#define massCutoffs(x) m_MassCutoffs[static_cast<int>(MASS_CUTOFF::x)]  // for convenience and readability - undefined at end of function

    double luminosity = GiantBranch::CalculateLuminosityAtPhaseBase_Static(massCutoffs(MHeF), m_AnCoefficients);
    double Mc_MHeF    = BaseStar::CalculateCoreMassGivenLuminosity_Static(luminosity, p_GBParams);
    double c          = (Mc_MHeF * Mc_MHeF * Mc_MHeF * Mc_MHeF) - (MC_L_C1 * pow(massCutoffs(MHeF), MC_L_C2));  // pow() is slow - use multiplication

    return std::min((0.95 * gbParams(McBAGB)), sqrt(sqrt(c + (MC_L_C1 * pow(p_Mass, MC_L_C2)))));               // sqrt is much faster than pow()

#undef massCutoffs
#undef gbParams
}


/*
 * Calculate core mass at the Base of the Giant Branch
 *
 * Hurley et al. 2000, eq 44
 *
 * For large enough M, we have McBGB ~ 0.098*Mass**(1.35)
 *
 * Static version required by CalculateGBParams_Static()
 *
 *
 * double CalculateCoreMassAtBGB_Static(const double p_Mass, const DBL_VECTOR &p_MassCutoffs, const DBL_VECTOR &p_AnCoefficients, const DBL_VECTOR &p_GBParams)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @param   [IN]    p_MassCutoffs               Mass cutoffs
 * @param   [IN]    p_AnCoefficients            a(n) coefficients
 * @param   [IN]    p_GBParams                  Giant Branch paramaters
 * @return                                      Core mass at the Base of the Giant Branch in Msol
 */
double GiantBranch::CalculateCoreMassAtBGB_Static(const double      p_Mass, 
                                                  const DBL_VECTOR &p_MassCutoffs, 
                                                  const DBL_VECTOR &p_AnCoefficients, 
                                                  const DBL_VECTOR &p_GBParams) {
#define gbParams(x) p_GBParams[static_cast<int>(GBP::x)]                // for convenience and readability - undefined at end of function
#define massCutoffs(x) p_MassCutoffs[static_cast<int>(MASS_CUTOFF::x)]  // for convenience and readability - undefined at end of function

    double luminosity = GiantBranch::CalculateLuminosityAtPhaseBase_Static(massCutoffs(MHeF), p_AnCoefficients);
    double Mc_MHeF    = BaseStar::CalculateCoreMassGivenLuminosity_Static(luminosity, p_GBParams);
    double c          = (Mc_MHeF * Mc_MHeF * Mc_MHeF * Mc_MHeF) - (MC_L_C1 * pow(massCutoffs(MHeF), MC_L_C2));  // pow() is slow - use multiplication

    return std::min((0.95 * gbParams(McBAGB)), sqrt(sqrt(c + (MC_L_C1 * pow(p_Mass, MC_L_C2)))));               // sqrt is much faster than pow()

#undef massCutoffs
#undef gbParams
}


/*
 * Calculate the core mass at which the Asymptotic Giant Branch phase is terminated in a SN/loss of envelope
 *
 * Hurley et al. 2000, eq 75
 *
 *
 * double CalculateCoreMassAtSupernova_Static(const double p_McBAGB)
 *
 * @param   [IN]    p_McBAGB                    Core mass at the Base of the Asymptotic Giant Branch in Msol
 * @return                                      Maximum core mass before supernova on the Asymptotic Giant Branch
 */
double GiantBranch::CalculateCoreMassAtSupernova_Static(const double p_McBAGB) {
    return std::max(MECS, (0.773 * p_McBAGB) - 0.35);
}


/*
 * Calculate core mass at Helium Ignition
 *
 * Hurley et al. 2000, eq 37 for low mass stars (coreMass - Luminosity relation)
 *
 * For M >= X use Hurley et al. 2000, eq 44, replacing Mc(LBGB(MHeF)) with Mc(LHeI(MHeF))
 * (see end of section 5.3 for this note)
 *
 *
 * double CalculateCoreMassAtHeIgnition(const double p_Mass)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @return                                      Core mass at Helium Ignition in Msol
 */
double GiantBranch::CalculateCoreMassAtHeIgnition(const double p_Mass) {
#define massCutoffs(x) m_MassCutoffs[static_cast<int>(MASS_CUTOFF::x)]  // for convenience and readability - undefined at end of function

    double coreMass;

    if (utils::Compare(p_Mass, massCutoffs(MHeF)) < 0) {
        double luminosity = CalculateLuminosityAtHeIgnition_Static(p_Mass, m_Alpha1, massCutoffs(MHeF), m_BnCoefficients);

        coreMass          = BaseStar::CalculateCoreMassGivenLuminosity_Static(luminosity, m_GBParams);
    }
    else {
        double luminosity_MHeF = CalculateLuminosityAtHeIgnition_Static(massCutoffs(MHeF), m_Alpha1, massCutoffs(MHeF), m_BnCoefficients);
        double Mc_MHeF         = BaseStar::CalculateCoreMassGivenLuminosity_Static(luminosity_MHeF, m_GBParams);
        double McBAGB          = CalculateCoreMassAtBAGB(p_Mass);
        double c               = (Mc_MHeF * Mc_MHeF * Mc_MHeF * Mc_MHeF) - (MC_L_C1 * pow(massCutoffs(MHeF), MC_L_C2)); // pow() is slow - use multiplication

        coreMass               = std::min((0.95 * McBAGB), sqrt(sqrt(c + (MC_L_C1 * pow(p_Mass, MC_L_C2)))));           // sqrt() is much faster than pow()
    }

    return coreMass;

#undef massCutoffs
}


/*
 * Calculate the mass of the core during second dredge up
 *
 * Hurley et al. 2000, eq 69 (not numbered in the paper, but between eqs 68 & 70)
 *
 *
 * double CalculateCoreMassAt2ndDredgeUp(const double p_McBAGB)
 *
 * @param       [IN]    p_McBAGB                Core mass at the Base of the Asymptotic Giant Branch
 * @return                                      Core Mass at second dredge up (McDU)
 */
double GiantBranch::CalculateCoreMassAt2ndDredgeUp_Static(const double p_McBAGB) {
    return utils::Compare(p_McBAGB, 0.8) >= 0 ? ((0.44 * p_McBAGB) + 0.448) : p_McBAGB;
}


/*
 * Calculate the dominant mass loss mechanism and associated rate for the star
 * at the current evolutionary phase.
 *
 * According to Hurley et al. 2000
 *
 * double CalculateMassLossRateHurley()
 *
 * @return                                      Mass loss rate in Msol per year
 */
double GiantBranch::CalculateMassLossRateHurley() {

    double dms = (utils::Compare(m_Luminosity, 4.0E3) > 0) ? CalculateMassLossRateNieuwenhuijzenDeJager() : 0.0;
    double dml = CalculateMassLossRateKudritzkiReimers();

    dms = std::max(dml, dms);

    if (utils::Compare(m_Mu, 1.0) < 0) {
        dml = CalculateMassLossRateWolfRayetLike(m_Mu);
        dms = std::max(dml, dms);
    }

    double tmp = m_Radius * sqrt(m_Luminosity) * 1.0E-5;
    if ((utils::Compare(m_Luminosity, 6.0E5) > 0) && (utils::Compare(tmp, 1.0) > 0)) {
        dml = CalculateMassLossRateLBV();
        dms = dms + dml;
    }

    return dms;
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
bool GiantBranch::IsMassRatioUnstable(const double p_AccretorMass, const bool p_AccretorIsDegenerate) {

    bool result = false;                                                                                                    // default is stable

    if (OPTIONS->MassTransferCriticalMassRatioGiant()) {
        result = p_AccretorIsDegenerate
                    ? (p_AccretorMass / m_Mass) < OPTIONS->MassTransferCriticalMassRatioGiantDegenerateAccretor()           // degenerate accretor
                    : (p_AccretorMass / m_Mass) < OPTIONS->MassTransferCriticalMassRatioGiantNonDegenerateAccretor();       // non-degenerate accretor
    }

    return result;
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                            LIFETIME / AGE CALCULATIONS                            //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * Calculate the time to Helium ignition tHeI
 *
 * Hurley et al. 2000, eq 43
 *
 *
 * double CalculateLifetimeToHeIgnition(const double p_Mass, const double p_Tinf1_FGB, const double p_Tinf2_FGB)
 *
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @param   [IN]    p_Tinf1_FGB                 tinf1_FGB (Timescales[TIMESCALE::tinf1_FGB]) - First Giant Branch tinf1
 * @param   [IN]    p_Tinf2_FGB                 tinf2_FGB (Timescales[TIMESCALE::tinf2_FGB]) - First Giant Branch tinf2
 * @return                                      Lifetime to He ignition (tHeI)
 */
double GiantBranch::CalculateLifetimeToHeIgnition(const double p_Mass, const double p_Tinf1_FGB, const double p_Tinf2_FGB) {
#define gbParams(x) m_GBParams[static_cast<int>(GBP::x)]                // for convenience and readability - undefined at end of function
#define massCutoffs(x) m_MassCutoffs[static_cast<int>(MASS_CUTOFF::x)]  // for convenience and readability - undefined at end of function

    double LHeI = CalculateLuminosityAtHeIgnition_Static(p_Mass, m_Alpha1, massCutoffs(MHeF), m_BnCoefficients);
    double p1   = gbParams(p) - 1.0;
    double q1   = gbParams(q) - 1.0;

    return utils::Compare(LHeI, gbParams(Lx)) <= 0
            ? p_Tinf1_FGB - (1.0 / (p1 * gbParams(AH) * gbParams(D))) * pow((gbParams(D) / LHeI), (p1 / gbParams(p)))
            : p_Tinf2_FGB - (1.0 / (q1 * gbParams(AH) * gbParams(B))) * pow((gbParams(B) / LHeI), (q1 / gbParams(q)));

#undef massCutoffs
#undef gbParams
}


/*
 * Calculate thermal timescale
 *
 * Kalogera & Webbink 1996, eq 2
 *
 *
 * double CalculateThermalTimescale(const double p_Mass, const double p_Radius, const double p_Luminosity, const double p_EnvMass)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @param   [IN]    p_Radius                    Radius in Rsol
 * @param   [IN]    p_Luminosity                Luminosity in Lsol
 * @param   [IN]    p_EnvMass                   Envelope mass in Msol
 * @return                                      Thermal timescale in Myr
*/
double GiantBranch::CalculateThermalTimescale(const double p_Mass, const double p_Radius, const double p_Luminosity, const double p_EnvMass) {
    return 30.0 * p_Mass * p_EnvMass / (p_Radius * p_Luminosity);       // G*Msol^2/(Lsol*Rsol) ~ 30
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                              SUPERNOVA CALCULATIONS                               //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * Calculate remnant type given COCoreMass
 *
 * Muller et al. 2016
 *
 *
 * STELLAR_TYPE CalculateRemnantTypeByMuller2016(const double p_COCoreMass)
 *
 * @param   [IN]    p_COCoreMass                COCoreMass in Msol
 * @return                                      Remnant type (stellar type)
 */
STELLAR_TYPE GiantBranch::CalculateRemnantTypeByMuller2016(const double p_COCoreMass) {

    STELLAR_TYPE stellarType;

    if (utils::Compare(p_COCoreMass, MCH) < 0) {
        stellarType = STELLAR_TYPE::OXYGEN_NEON_WHITE_DWARF;
    }

         if (utils::Compare(p_COCoreMass, 3.6 ) < 0) { stellarType = STELLAR_TYPE::NEUTRON_STAR; }
    else if (utils::Compare(p_COCoreMass, 4.05) < 0) { stellarType = STELLAR_TYPE::BLACK_HOLE; }
    else if (utils::Compare(p_COCoreMass, 4.6 ) < 0) { stellarType = STELLAR_TYPE::NEUTRON_STAR; }
    else if (utils::Compare(p_COCoreMass, 5.7 ) < 0) { stellarType = STELLAR_TYPE::BLACK_HOLE; }
    else if (utils::Compare(p_COCoreMass, 6.0 ) < 0) { stellarType = STELLAR_TYPE::NEUTRON_STAR; }
    else                                             { stellarType = STELLAR_TYPE::BLACK_HOLE; }

    return stellarType;
}


/*
 * Calculate remnant mass given COCoreMass and HeCoreMass
 *
 * Mandel & Mueller, 2020
 *
 *
 * double CalculateRemnantMassByMullerMandel (const double p_COCoreMass, const double p_HeCoreMass)
 *
 * @param   [IN]    p_COCoreMass                COCoreMass in Msol
 * @param   [IN]    p_HeCoreMass                HeCoreMass in Msol
 * @return                                      Remnant mass in Msol
 */
double GiantBranch::CalculateRemnantMassByMullerMandel(const double p_COCoreMass, const double p_HeCoreMass){
    double remnantMass=0;   
    double pBH=0;
    double pCompleteCollapse=0;
    

    if (utils::Compare(p_COCoreMass, MULLERMANDEL_M1) < 0) {
	pBH=0;
    }
    else if (utils::Compare(p_COCoreMass, MULLERMANDEL_M3) < 0) {
    	pBH=1.0/(MULLERMANDEL_M3-MULLERMANDEL_M1)*(p_COCoreMass-MULLERMANDEL_M1);
    }
    else {
	pBH=1.0;
    } 
 
    if(utils::Compare(RAND->Random(0,1), pBH) < 0) {  	// this is a BH
        if(utils::Compare(p_COCoreMass, MULLERMANDEL_M4) < 0)
		pCompleteCollapse=1.0/(MULLERMANDEL_M4-MULLERMANDEL_M1)*(p_COCoreMass-MULLERMANDEL_M1);
        else
		pCompleteCollapse=1.0;

	if(utils::Compare(RAND->Random(0,1), pCompleteCollapse) < 0) {
		remnantMass=p_HeCoreMass;
        }
	else {
		while(remnantMass<MULLERMANDEL_MAXNS || remnantMass > (p_COCoreMass+p_HeCoreMass) ){ 
			remnantMass = MULLERMANDEL_MUBH*p_COCoreMass + RAND->RandomGaussian(MULLERMANDEL_SIGMABH);
		}
	}
    }
    else {						// this is an NS
	if (utils::Compare(p_COCoreMass, MULLERMANDEL_M1) < 0) {
		while(remnantMass < MULLERMANDEL_MINNS || remnantMass > MULLERMANDEL_MAXNS || remnantMass > (p_COCoreMass+p_HeCoreMass) ){
			remnantMass = MULLERMANDEL_MU1 + RAND->RandomGaussian(MULLERMANDEL_SIGMA1);
		}
	}
	else if (utils::Compare(p_COCoreMass, MULLERMANDEL_M2) < 0) {
                while(remnantMass < MULLERMANDEL_MINNS || remnantMass > MULLERMANDEL_MAXNS || remnantMass > (p_COCoreMass+p_HeCoreMass) ){
                        remnantMass = MULLERMANDEL_MU2A + 
			MULLERMANDEL_MU2B/(MULLERMANDEL_M2-MULLERMANDEL_M1)*(p_COCoreMass-MULLERMANDEL_M1)+RAND->RandomGaussian(MULLERMANDEL_SIGMA2);
                }
        }
        else {
                while(remnantMass < MULLERMANDEL_MINNS || remnantMass > MULLERMANDEL_MAXNS || remnantMass > (p_COCoreMass+p_HeCoreMass) ){
                        remnantMass = MULLERMANDEL_MU3A + 
			MULLERMANDEL_MU3B/(MULLERMANDEL_M3-MULLERMANDEL_M2)*(p_COCoreMass-MULLERMANDEL_M2)+RAND->RandomGaussian(MULLERMANDEL_SIGMA3);
                }
        }
    }

    return remnantMass;
}


/*
 * Calculate remnant mass given Mass and COCoreMass
 *
 * Muller et al. 2016
 *
 *
 * double CalculateRemnantMassByMuller2016(const double p_Mass, const double p_COCoreMass)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @param   [IN]    p_COCoreMass                COCoreMass in Msol
 * @return                                      Remnant mass in Msol
 */
double GiantBranch::CalculateRemnantMassByMuller2016(const double p_Mass, const double p_COCoreMass) {
	// ALEJANDRO - 14/03/2017 - Updated prescription of Bernhard Muller's models, incorporating lower metallicity models for core masses 1.372 <= m_{C/O} < 1.65
	// Where m_{C/O} = M_{C/O}/M_{\odot}
	// First approach is to model the lower limit using Chnadrasekhar mass, m_{C/O} > Mch = 1.44

    double	remnantMass; 					                                                                        // Limit mass for a White Dwarf units Msun.

    if (utils::Compare(p_COCoreMass, 1.44) < 0) {
        remnantMass = 1.4;
    }
	else if (utils::Compare(p_COCoreMass, 1.49) < 0) { remnantMass = 1.21 - (0.4  * (p_COCoreMass - 1.372)); }
	else if (utils::Compare(p_COCoreMass, 1.65) < 0) { remnantMass = 1.16;                                   }
    else if (utils::Compare(p_COCoreMass, 2.4 ) < 0) { remnantMass = 1.32 + (0.3  * (p_COCoreMass - 1.65));  }
    else if (utils::Compare(p_COCoreMass, 3.2 ) < 0) { remnantMass = 1.42 + (0.7  * (p_COCoreMass - 2.4));   }
    else if (utils::Compare(p_COCoreMass, 3.6 ) < 0) { remnantMass = 1.32 + (0.25 * (p_COCoreMass - 3.2));   }
    else if (utils::Compare(p_COCoreMass, 4.05) < 0) { remnantMass = p_Mass * NEUTRINO_LOSS_FALLBACK_FACTOR; }      // Going to be a Black Hole
    else if (utils::Compare(p_COCoreMass, 4.6 ) < 0) { remnantMass = 1.5;                                    }
    else if (utils::Compare(p_COCoreMass, 5.7 ) < 0) { remnantMass = p_Mass * NEUTRINO_LOSS_FALLBACK_FACTOR; }      // Going to be a Black Hole
    else if (utils::Compare(p_COCoreMass, 6.0 ) < 0) { remnantMass = 1.64 - (0.2  * (p_COCoreMass - 5.7));   }
    else                                             { remnantMass = p_Mass * NEUTRINO_LOSS_FALLBACK_FACTOR; }      // Going to be a Black Hole

    return remnantMass;
}


/*
 * Calculate the baryonic mass of the remnant
 *
 * Fryer et al. 2012, eq 12
 *
 *
 * double CalculateBaryonicRemnantMass(const double p_ProtoMass, double p_FallbackMass)
 *
 * @param   [IN]    p_ProtoMass                 Mass of proto compact object in Msol
 * @param   [IN]    p_FallbackMass              Mass falling back onto proto compact object Mfb = fb*(MpreSN - Mproto)
 * @return                                      Baryonic mass of the remnant in Msol
 */
double GiantBranch::CalculateBaryonicRemnantMass(const double p_ProtoMass, double p_FallbackMass) {
    return p_ProtoMass + p_FallbackMass; //Mfb = fb*(MpreSN - Mproto);
}


/*
 * Calculate the gravitational mass of the remnant
 *
 * Fryer et al. 2012, eqs 13 & 14
 *
 *
 * double CalculateGravitationalRemnantMass(const double p_BaryonicRemnantMass)
 *
 * @param   [IN]    p_BaryonicRemnantMass       Baryonic remnant mass in Msol
 * @return                                      Gravitational mass of the remnant in Msol
 */
double GiantBranch::CalculateGravitationalRemnantMass(const double p_BaryonicRemnantMass) {
    // decide whether to calculate GravitationalRemnantMass from Fryer+2012, Eq.13 for Neutron Star or Black Hole 
    // then calculate GravitationalRemnantMass 
    return (utils::Compare(p_BaryonicRemnantMass, m_BaryonicMassOfMaximumNeutronStarMass) < 0) 
            ? utils::SolveQuadratic(0.075, 1.0, -p_BaryonicRemnantMass)                 // Neutron Star
            : 0.9 * p_BaryonicRemnantMass;                                              // Black Hole
}


/*
 * Calculate the mass falling back onto the proto compact object
 *
 * Fryer et al. 2012, eq 11
 *
 *
 * double CalculateFallbackMass(const double p_PreSNMass, const double p_ProtoMass, const double p_Fallback)
 *
 * @param   [IN]    p_PreSNMass                 Pre supernova stellar mass in Msol
 * @param   [IN]    p_ProtoMass                 Pre supernova Carbon Oxygen (CO) core mass in Msol
 * @param   [IN]    p_Fallback                  Fraction of mass falling back onto proto object
 * @return                                      Mass falling back onto proto object
 */
double GiantBranch::CalculateFallbackMass(const double p_PreSNMass, const double p_ProtoMass, const double p_Fallback) {
    return p_Fallback * (p_PreSNMass - p_ProtoMass);
}


/*
 * Calculate the core mass of the proto compact object using the rapid prescription
 *
 * Fryer et al. 2012, eq 15 (based on Woosley et al. 2002)
 *
 * Sets Mproto = 1.0 Msol regardless of progenitor
 *
 *
 * double CalculateProtoCoreMassRapid()
 *
 * @return                                      Mass of Fe/Ni proto core in Msol (always 1.0)
 */
double GiantBranch::CalculateProtoCoreMassRapid() {
    return 1.0;
}


/*
 * Calculate the fraction of mass falling back onto the proto compact object using the rapid prescription
 *
 * Fryer et al. 2012, eq 15 & 16
 *
 *
 * double CalculateFallbackFractionRapid(const double p_PreSNMass, const double p_ProtoMass, const double p_COCoreMass)
 *
 * @param   [IN]    p_PreSNMass                 Pre supernova stellar mass in Msol
 * @param   [IN]    p_ProtoMass                 Fe/Ni proto compact object mass in Msol
 * @param   [IN]    p_COCoreMass                Pre supernova Carbon Oxygen (CO) core mass in Msol
 * @return                                      Fraction of mass falling back onto proto object
 */
double GiantBranch::CalculateFallbackFractionRapid(const double p_PreSNMass, const double p_ProtoMass, const double p_COCoreMass) {

    double fb;

    double spread = p_PreSNMass - p_ProtoMass;

    if (utils::Compare(p_COCoreMass, 2.5) < 0) {
        fb = 0.2 / spread;
    }
    else if (utils::Compare(p_COCoreMass, 6.0) < 0) {
        fb = ((0.286 * p_COCoreMass) - 0.514) / spread;
    }
    else if (utils::Compare(p_COCoreMass, 7.0) < 0) {
        fb = 1.0;
    }
    else if (utils::Compare(p_COCoreMass, 11.0) < 0) {
        double a1 = 0.25 - (1.275 / spread);
        double b1 = (-11.0 * a1) + 1.0;
        fb        = (a1 * p_COCoreMass) + b1;
    }
    else {
        fb = 1.0;
    }

    return std::max(0.0, std::min(1.0, fb));
}


/*
 * Calculate the mass of the proto core
 *
 * Fryer et al. 2012, eq 18
 *
 *
 * double CalculateProtoCoreMassDelayed(const double p_COCoreMass)
 *
 * @param   [IN]    p_COCoreMass                Pre supernova Carbon Oxygen (CO) core mass in Msol
 * @return                                      Mass of the Fe/Ni proto core in Msol
 */
double GiantBranch::CalculateProtoCoreMassDelayed(const double p_COCoreMass) {

    double protoMass;

         if (utils::Compare(p_COCoreMass,  3.5) < 0) { protoMass = 1.2; }
    else if (utils::Compare(p_COCoreMass,  6.0) < 0) { protoMass = 1.3; }
    else if (utils::Compare(p_COCoreMass, 11.0) < 0) { protoMass = 1.4; }
    else                                             { protoMass = 1.6; }

    return protoMass;
}


/*
 * Calculate the fraction of mass falling back onto the proto compact object using the delayed prescription
 *
 * Fryer et al. 2012, eq 19
 *
 *
 * double CalculateFallbackFractionDelayed(const double p_PreSNMass, const double p_ProtoMass, const double p_COCoreMass)
 *
 * @param   [IN]    p_PreSNMass                 Pre supernova stellar mass in Msol
 * @param   [IN]    p_ProtoMass                 Fe/Ni proto compact object mass in Msol
 * @param   [IN]    p_COCoreMass                Pre supernova Carbon Oxygen (CO) core mass in Msol
 * @return                                      Fraction of mass falling back onto proto object
 */
double GiantBranch::CalculateFallbackFractionDelayed(const double p_PreSNMass, const double p_ProtoMass, const double p_COCoreMass) {

    double fb;

    double spread = p_PreSNMass - p_ProtoMass;

    if (utils::Compare(p_COCoreMass, 2.5) < 0) {
        fb = 0.2 / spread;
    }
    else if (utils::Compare(p_COCoreMass, 3.5) < 0) {
        fb = ((0.5 * p_COCoreMass) - 1.05) / spread;
    }
    else if (utils::Compare(p_COCoreMass, 11.0) < 0) {
        double a2 = 0.133 - (0.093 / spread);
        double b2 = (-11.0 * a2) + 1.0;
        fb        = (a2 * p_COCoreMass) + b2;
    }
    else {
        fb = 1.0;
    }

    return std::max(0.0, std::min(1.0, fb));
}


/*
 * Calculate the remnant mass using one of the Fryer SN Engines
 *
 * Ref?  Fryer et al. 2012?
 *
 *
 * std::tuple<double, double> CalculateRemnantMassByFryer2012(const SNE p_Engine, const double p_Mass, const double p_COCoreMass)
 *
 * @param   [IN]    p_Engine                    The SN Engine to use for the calculation
 * @param   [IN]    p_Mass                      Pre supernova mass in Msol
 * @param   [IN]    p_COCoreMass                Pre supernova Carbon Oxygen (CO) core mass in Msol
 * @return                                      Tuple containing Remnant mass in Msol and updated fraction of mass falling back onto compact object
 */
std::tuple<double, double> GiantBranch::CalculateRemnantMassByFryer2012(const SN_ENGINE p_Engine, const double p_Mass, const double p_COCoreMass) {

    double mProto;
    double fallbackMass;
    double baryonicRemnantMass;

    double fallbackFraction         = 0.0;
    double gravitationalRemnantMass = 0.0;

    switch (p_Engine) {                                                                                     // which SN_ENGINE?

        case SN_ENGINE::DELAYED:                                                                            // DELAYED

            mProto           = CalculateProtoCoreMassDelayed(p_COCoreMass);
            fallbackFraction = CalculateFallbackFractionDelayed(p_Mass, mProto, p_COCoreMass);
            fallbackMass     = CalculateFallbackMass(p_Mass, mProto, fallbackFraction);

            baryonicRemnantMass      = CalculateBaryonicRemnantMass(mProto, fallbackMass);
            gravitationalRemnantMass = CalculateGravitationalRemnantMass(baryonicRemnantMass);
            break;

        case SN_ENGINE::RAPID:                                                                              // RAPID

            mProto           = CalculateProtoCoreMassRapid();
            fallbackFraction = CalculateFallbackFractionRapid(p_Mass, mProto, p_COCoreMass);
            fallbackMass     = CalculateFallbackMass(p_Mass, mProto, fallbackFraction);

            baryonicRemnantMass      = CalculateBaryonicRemnantMass(mProto, fallbackMass);
            gravitationalRemnantMass = CalculateGravitationalRemnantMass(baryonicRemnantMass);
            break;

        default:                                                                                            // unknown SN_ENGINE
            SHOW_WARN(ERROR::UNKNOWN_SN_ENGINE, "Using defaults");                                          // show warning
    }

    return std::make_tuple(gravitationalRemnantMass, fallbackFraction);
}


/*
 * Calculate fallback using linear interpolation
 *
 * Belczynski et al. 2002
 * Description?
 *
 *
 * double CalculateFallbackByBelczynski2002(const double p_COCoreMass)
 *
 * @param   [IN]    p_COCoreMass                Pre supernova Carbon Oxygen (CO) core mass in Msol
 * @return                                      Fraction of mass falling back onto compact object
 */
double GiantBranch::CalculateFallbackByBelczynski2002(const double p_COCoreMass) {
    return utils::Compare(p_COCoreMass, 5.0) <= 0 ? 0.0 : ((utils::Compare(p_COCoreMass, 7.6) < 0) ? (p_COCoreMass - 5.0) / 2.6 : 1.0);
}


/*
 * Calculate remnant mass
 *
 * Formula used by Hurley code not given in Hurley et al 2000 but in Belczynski et al. 2002
 *
 *
 * double CalculateRemnantMassByBelczynski2002(const double p_Mass, const double p_COCoreMass, const double p_FallbackFraction)
 *
 * @param   [IN]    p_COCoreMass                Pre supernova Carbon Oxygen (CO) core mass in Msol
 * @param   [IN]    p_FallbackFraction          Fraction of mass falling back onto compact object
 * @return                                      Remnant mass in Msol
 */
double GiantBranch::CalculateRemnantMassByBelczynski2002(const double p_Mass, const double p_COCoreMass, const double p_FallbackFraction) {
    double McFeNi = (utils::Compare(p_COCoreMass, 2.5) < 0) ? (0.161767 * p_COCoreMass) + 1.067055 : (0.314154 * p_COCoreMass) + 0.686088; // Iron core mass
    return McFeNi + (p_FallbackFraction * (p_Mass - McFeNi));
}


/*
 * Driver function for Core Collapse Supernovas
 *
 * This function determines which prescription is used for the core collapse SN (via program options)
 *
 * The function calls prescription functions that update the following parameters:
 *      Mass, stellarType, drawnKickVelocity, kickVelocity
 *
 * At the end of this function we set the following parameters which are (so far) independent of the
 * ccSN prescriptions (but do depend on the parameters above):
 *      Luminosity, Radius, Temperature, supernova events: current = SN, past = CCSN
 *
 *
 * STELLAR_TYPE IsCoreCollapseSN(const SNE SNEngine)
 *
 * @param   [IN]    SNEngine                    Fryer SN Engine to use (if required)
 * @return                                      The stellar type to which the star should evolve
 */
STELLAR_TYPE GiantBranch::IsCoreCollapseSN(const SN_ENGINE SNEngine) {

    STELLAR_TYPE stellarType = m_StellarType;

    switch (OPTIONS->RemnantMassPrescription()) {                                                           // which prescription?

        case REMNANT_MASS_PRESCRIPTION::POSTITNOTE:                                                         // POSTITNOTE

            m_SupernovaDetails.fallbackFraction = 0.0;                                                      // Not defined
            m_Mass                              = BH::CalculateRemnantMass_Static(m_MZAMS);
            break;

        case REMNANT_MASS_PRESCRIPTION::HURLEY2000:                                                         // Hurley 2000

            m_SupernovaDetails.fallbackFraction = 0.0;                                                      // Not defined
            m_Mass                              = NS::CalculateRemnantMass_Static(m_CoreMass);
            break;

        case REMNANT_MASS_PRESCRIPTION::BELCZYNSKI2002:                                                     // Belczynski 2002

            m_SupernovaDetails.fallbackFraction = CalculateFallbackByBelczynski2002(m_CoreMass);
            m_Mass                              = CalculateRemnantMassByBelczynski2002(m_Mass, m_CoreMass, m_SupernovaDetails.fallbackFraction);
            break;

        case REMNANT_MASS_PRESCRIPTION::FRYER2012:                                                          // Fryer 2012

            std::tie(m_Mass, m_SupernovaDetails.fallbackFraction) = CalculateRemnantMassByFryer2012(SNEngine, m_Mass, m_COCoreMass);
            break;

        case REMNANT_MASS_PRESCRIPTION::MULLER2016:                                                         // Muller 2016
        case REMNANT_MASS_PRESCRIPTION::MULLER2016MAXWELLIAN:                                               // Muller 260016 - Maxwellian

            m_Mass = CalculateRemnantMassByMuller2016(m_Mass, m_COCoreMass);

            // JR: todo: fallback fraction not calculated here?
            break;

	case REMNANT_MASS_PRESCRIPTION::MULLERMANDEL:                                                         // Mandel & Mueller, 2020

            m_Mass = CalculateRemnantMassByMullerMandel(m_COCoreMass, m_HeCoreMass);

            break;

        default:                                                                                            // unknown prescription

            m_Mass                              = 0.0;
            m_SupernovaDetails.fallbackFraction = 0.0;

            m_Error = ERROR::UNKNOWN_REMNANT_MASS_PRESCRIPTION;                                             // set error number
            SHOW_ERROR(ERROR::UNKNOWN_REMNANT_MASS_PRESCRIPTION, "Using default");                          // show error
    }

    // Set the stellar type to which the star should evolve (either use prescription or MAXIMUM_NS_MSS)
    if (OPTIONS->RemnantMassPrescription() == REMNANT_MASS_PRESCRIPTION::MULLER2016) {
        stellarType = CalculateRemnantTypeByMuller2016(m_COCoreMass);
    }
    else if (OPTIONS->RemnantMassPrescription() == REMNANT_MASS_PRESCRIPTION::MULLERMANDEL) {
	if(utils::Compare(m_Mass, MULLERMANDEL_MAXNS ) > 0) 
		stellarType = STELLAR_TYPE::BLACK_HOLE; 
        else
		stellarType = STELLAR_TYPE::NEUTRON_STAR;
    }	

    else if (utils::Compare(m_Mass, OPTIONS->MaximumNeutronStarMass()) > 0) {
        std::tie(m_Luminosity, m_Radius, m_Temperature) = BH::CalculateCoreCollapseSNParams_Static(m_Mass);
        stellarType = STELLAR_TYPE::BLACK_HOLE;
    }
    else {
        std::tie(m_Luminosity, m_Radius, m_Temperature) = NS::CalculateCoreCollapseSNParams_Static(m_Mass);
        stellarType = STELLAR_TYPE::NEUTRON_STAR;
    }

    SetSNCurrentEvent(SN_EVENT::CCSN);                                                                      // flag core-collapse SN happening now
    SetSNPastEvent(SN_EVENT::CCSN);                                                                         // ... and will be a past event

    return stellarType;
}


/*
 * Resolve Electron capture Supernova
 *
 * Calculate the mass of the remnant and set remnant type - always a Neutron Star
 * Updates attributes of star; sets SN flags
 *
 *
 * Short hand wavy story. The core ignites ONeMg and collapses through electron
 * capture. This core is less dense than the heavier Fe Core collapses and thus
 * people expect the neutrinos to escape easier resulting in a zero kick
 * TODO a nice reference (Nomoto 1984 for thorough discussion and a summary in Nomoto 1987)
 *
 *
 * STELLAR_TYPE IsElectronCaptureSN()
 *
 * @return                                      Stellar type of remnant (always STELLAR_TYPE::NEUTRON_STAR)
 */
STELLAR_TYPE GiantBranch::IsElectronCaptureSN() {

    m_Mass       = MECS_REM;                                                        // defined in constant.h
    m_Radius     = NS::CalculateRadiusOnPhase_Static(m_Mass);                       // neutronStarRadius in Rsol
    m_Luminosity = NS::CalculateLuminosityOnPhase_Static(m_Mass, m_Age);

    SetSNCurrentEvent(SN_EVENT::ECSN);                                              // electron capture SN happening now
    SetSNPastEvent(SN_EVENT::ECSN);                                                 // ... and will be a past event

    return STELLAR_TYPE::NEUTRON_STAR;
}


/*
 * Resolve Type IIa Supernova
 *
 * Zero attributes - leaves a Massless remnant
 *
 *
 * STELLAR_TYPE IsTypeIIaSN()
 *
 * @return                                      Stellar type of remnant (always STELLAR_TYPE::MASSLESS_REMNANT)
 */
STELLAR_TYPE GiantBranch::IsTypeIIaSN() {

    m_Mass              = 0.0;
    m_Radius            = 0.0;
    m_Luminosity        = 0.0;
    m_Temperature       = 0.0;

    m_SupernovaDetails.drawnKickVelocity = 0.0;
    m_SupernovaDetails.kickVelocity      = 0.0;

    return STELLAR_TYPE::MASSLESS_REMNANT;
}


/*
 * Resolve Pair-Instability Supernova
 *
 * Calculate the mass of the remnant and set remnant type according to mass
 * Updates attributes of star; sets SN events
 *
 *
 * Short handwavy story.  The core is hot and massive enough that there is a significant
 * amount of high-energetic gamma rays which can create an electron-positron pair.
 * When a significant amounts of pairs are created the photons taken away
 * reduce the radiative pressure. The core contracts becomes hotter creating
 * more pairs. This runaway process will end in a SN that explodes the entire core without
 * leaving a remnant.
 *
 *
 * STELLAR_TYPE IsPulsationalPairInstabilitySN()
 *
 * @return                                      Stellar type of remnant
 */
STELLAR_TYPE GiantBranch::IsPairInstabilitySN() {

    m_Luminosity  = 0.0;
    m_Radius      = 0.0;
    m_Temperature = 0.0;

    m_SupernovaDetails.drawnKickVelocity = 0.0;
    m_SupernovaDetails.kickVelocity      = 0.0;
    m_SupernovaDetails.fallbackFraction  = 0.0;

    SetSNCurrentEvent(SN_EVENT::PISN);                                                                  // pair instability SN happening now
    SetSNPastEvent(SN_EVENT::PISN);                                                                     // ... and will be a past event

    return STELLAR_TYPE::MASSLESS_REMNANT;
}


/*
 * Resolve Pulsational Pair-Instability Supernova
 *
 * Calculate the mass of the remnant and set remnant type according to mass
 * Updates attributes of star; sets SN events
 *
 *
 * STELLAR_TYPE IsPulsationalPairInstabilitySN()
 *
 * @return                                      Stellar type of remnant
 */
STELLAR_TYPE GiantBranch::IsPulsationalPairInstabilitySN() {

    STELLAR_TYPE stellarType = m_StellarType;

    double baryonicMass;
    switch (OPTIONS->PulsationalPairInstabilityPrescription()) {                                        // which prescription?

        case PPI_PRESCRIPTION::COMPAS:                                                                  // Woosley 2017 https://arxiv.org/abs/1608.08939
            baryonicMass = m_HeCoreMass;                                                                // strip off the hydrogen envelope if any was left (factor of 0.9 applied in BH::CalculateNeutrinoMassLoss_Static)
            m_Mass = BH::CalculateNeutrinoMassLoss_Static(baryonicMass);                                // convert to gravitational mass due to neutrino mass loss

            break;

        case PPI_PRESCRIPTION::STARTRACK:                                                               // Belczynski et al. 2016 https://arxiv.org/abs/1607.03116
            baryonicMass = m_HeCoreMass;                                                                // strip off the hydrogen envelope if any was left (factor of 0.9 applied in BH::CalculateNeutrinoMassLoss_Static)
            m_Mass = BH::CalculateNeutrinoMassLoss_Static(baryonicMass);                                // convert to gravitational mass due to neutrino mass loss

            break;

        case PPI_PRESCRIPTION::MARCHANT: {                                                              // Marchant et al. 2018 https://arxiv.org/abs/1810.13412

            // pow() is slow - use multiplication
            double HeCoreMass_2 = m_HeCoreMass * m_HeCoreMass;
            double HeCoreMass_3 = HeCoreMass_2 * m_HeCoreMass;
            double HeCoreMass_4 = HeCoreMass_2 * HeCoreMass_2;
            double HeCoreMass_5 = HeCoreMass_3 * HeCoreMass_2;
            double HeCoreMass_6 = HeCoreMass_3 * HeCoreMass_3;
            double HeCoreMass_7 = HeCoreMass_6 * m_HeCoreMass;

            double ratioOfRemnantToHeCoreMass = std::max(0.0, std::min(1.0, (-1.63057326E-08 * HeCoreMass_7) +
                                                                            ( 5.36316755E-06 * HeCoreMass_6) +
                                                                            (-7.52206933E-04 * HeCoreMass_5) +
                                                                            ( 5.83107626E-02 * HeCoreMass_4) +
                                                                            (-2.69801221E+00 * HeCoreMass_3) +
                                                                            ( 7.45060098E+01 * HeCoreMass_2) +
                                                                            (-1.13694590E+03 * m_HeCoreMass) +
                                                                              7.39643451E+03));

            baryonicMass = ratioOfRemnantToHeCoreMass * m_HeCoreMass;                                   // strip off the hydrogen envelope if any was left (factor of 0.9 applied in BH::CalculateNeutrinoMassLoss_Static)
            m_Mass = BH::CalculateNeutrinoMassLoss_Static(baryonicMass);                                // convert to gravitational mass due to neutrino mass loss

            } break;

        default:                                                                                        // unknown prescription
            SHOW_WARN(ERROR::UNKNOWN_PPI_PRESCRIPTION);                                                 // show warning
            m_Mass = m_HeCoreMass;                                                                      // strip off the hydrogen envelope if any was left -- factor of 0.9 applied later
    }

    if (utils::Compare(m_Mass, 0.0) <= 0) {                                                             // remnant mass <= 0?
        stellarType = IsPairInstabilitySN();                                                            // yes - PISN rather than PPISN
    }
    else {                                                                                              // no - PPISN

        SetSNCurrentEvent(SN_EVENT::PPISN);                                                             // pulsational pair instability SN happening now
        SetSNPastEvent(SN_EVENT::PPISN);                                                                // ... and will be a past event

        stellarType   = STELLAR_TYPE::BLACK_HOLE;                                                       // -> black hole
        
        m_Luminosity  = BH::CalculateLuminosityOnPhase_Static();                                        // black hole luminosity
        m_Radius      = BH::CalculateRadiusOnPhase_Static(m_Mass);                                      // Schwarzschild radius (not correct for rotating BH)
        m_Temperature = CalculateTemperatureOnPhase(m_Luminosity, m_Radius);
        m_SupernovaDetails.fallbackFraction = 1.0;                                                      // fraction of mass that falls back
    }

    return stellarType;
}


/*
 * The main supernova function
 *
 * This function determines the type of the supernova and calls the appropriate functions
 * to calculate attributes correctly, and to determine the type of remnant to which the
 * star should evolve.
 *
 *
 * STELLAR_TYPE ResolveSupernova()
 *
 * @return                                      Stellar type of remnant
 */
STELLAR_TYPE GiantBranch::ResolveSupernova() {

    STELLAR_TYPE stellarType = m_StellarType;

    if (IsSupernova()) {                                                                            // has gone supernova
        // squirrel away some attributes before they get changed...
        m_SupernovaDetails.totalMassAtCOFormation  = m_Mass;
        m_SupernovaDetails.HeCoreMassAtCOFormation = m_HeCoreMass;
        m_SupernovaDetails.COCoreMassAtCOFormation = m_COCoreMass;
        m_SupernovaDetails.coreMassAtCOFormation   = m_CoreMass;

        double snMass = CalculateInitialSupernovaMass();                                            // calculate SN initial mass

        SetSNHydrogenContent();                                                                     // ALEJANDRO - 04/05/2018 - Check if the SN is H-rich or H-poor. For now, classify it for all possible SNe and not only CCSN forming NS.

        if (                             OPTIONS->UsePulsationalPairInstability()              &&
            utils::Compare(m_HeCoreMass, OPTIONS->PulsationalPairInstabilityLowerLimit()) >= 0 &&
            utils::Compare(m_HeCoreMass, OPTIONS->PulsationalPairInstabilityUpperLimit()) <= 0) {   // Pulsational Pair Instability SuperNova

            stellarType = IsPulsationalPairInstabilitySN();
        }
        else if (                        OPTIONS->UsePairInstabilitySupernovae()    &&
            utils::Compare(m_HeCoreMass, OPTIONS->PairInstabilityLowerLimit()) >= 0 &&
            utils::Compare(m_HeCoreMass, OPTIONS->PairInstabilityUpperLimit()) <= 0) {              // Pair Instability SuperNova

            stellarType = IsPairInstabilitySN();
        }
        else if (utils::Compare(snMass, MCBUR1) < 0) {                                                 // Type IIa SuperNova
            stellarType = IsTypeIIaSN();
        }
        else if (utils::Compare(snMass, MCBUR2) < 0) {                                                // Electron Capture SuperNova
            stellarType = IsElectronCaptureSN();
        }
        else {                                                                                      // Core Collapse SuperNova
            stellarType = IsCoreCollapseSN(OPTIONS->FryerSupernovaEngine());
        }

        // SIMON : What is this line for? Why 'reset' things to 0?                                  // JR: todo: check this
        // reset values to zero
        m_COCoreMass  = 0.0;
        m_HeCoreMass  = 0.0;
        m_CoreMass    = 0.0;
        m_Mass0       = 0.0;
        m_EnvMass     = 0.0;
        m_Age         = 0.0;

    	CalculateSNKickVelocity(m_Mass, m_SupernovaDetails.totalMassAtCOFormation - m_Mass, stellarType);
    }

    return stellarType;
}
