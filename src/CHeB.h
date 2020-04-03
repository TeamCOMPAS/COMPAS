#ifndef __CHeB_h__
#define __CHeB_h__

#include "constants.h"
#include "typedefs.h"
#include "utils.h"

#include "FGB.h"


class BaseStar;
class FGB;

class CHeB: virtual public BaseStar, public FGB {

public:

    CHeB(const BaseStar &p_BaseStar, const bool p_Initialise = true) : BaseStar(p_BaseStar), FGB(p_BaseStar, false) {
        if (p_Initialise) Initialise();
    }

    CHeB& operator = (const BaseStar &p_BaseStar) {
        static_cast<BaseStar&>(*this) = p_BaseStar;
        Initialise();
        return *this;
    }


    // member functions - alphabetically

           double CalculateMinimumLuminosityOnPhase(const double      p_Mass,
                                                    const double      p_Alpha1,
                                                    const double      p_MHeF,
                                                    const double      p_MFGB,
                                                    const DBL_VECTOR &p_BnCoefficients);

    static double CalculateMinimumRadiusOnPhase_Static(const double      p_Mass,
                                                       const double      p_CoreMass,
                                                       const double      p_Alpha1,
                                                       const double      p_MHeF,
                                                       const double      p_MFGB,
                                                       const double      p_LuminosityOnPhase,
                                                       const DBL_VECTOR &p_BnCoefficients);


protected:

    void Initialise() {
    #define massCutoffs(x) m_MassCutoffs[static_cast<int>(MASS_CUTOFF::x)]  // for convenience and readability - undefined at end of function

        m_StellarType = STELLAR_TYPE::CORE_HELIUM_BURNING;                                                                                                      // Set stellar type
        CalculateTimescales();                                                                                                                                  // Initialise timescales
        m_Age = m_Timescales[static_cast<int>(TIMESCALE::tHeI)];                                                                                                // Set age appropriately
        m_MinimumLuminosityOnPhase = CalculateMinimumLuminosityOnPhase(massCutoffs(MHeF), m_Alpha1, massCutoffs(MHeF), massCutoffs(MFGB), m_BnCoefficients);    // Calculate once, not many
    #undef massCutoffs
    }


    // member functions - alphabetically

    double          CalculateBluePhaseFBL(const double p_Mass);

    double          CalculateCOCoreMassOnPhase()                                 { return 0.0; }                                                                // McCO(CHeB) = 0.0

    double          CalculateCoreMassAtPhaseEnd()                                { return CalculateCoreMassAtBAGB(m_Mass0); }                                   // Use class member variables
    double          CalculateCoreMassOnPhase(const double p_Mass, const double p_Tau);
    double          CalculateCoreMassOnPhase()                                   { return CalculateCoreMassOnPhase(m_Mass0, m_Tau); }                           // Use class member variables

    double          CalculateEnvelopeMassAtPhaseEnd(const double p_Tau)          { return m_EnvMass; }                                                          // NO-OP
    double          CalculateEnvelopeMassOnPhase(const double p_Tau)             { return BaseStar::CalculateEnvelopeMassOnPhase(p_Tau); }

    double          CalculateHeCoreMassAtPhaseEnd()                              { return m_CoreMass; }

    double          CalculateGyrationRadius()                                    { return 0.21; }                                                               // Hurley et al., 2000, after eq 109 for n=3/2 polytrope or dense convective core. Single number approximation.

    double          CalculateLambdaDewi();
    double          CalculateLambdaNanjing();

    double          CalculateLifetimeOnBluePhase(const double p_Mass);
    double          CalculateLifetimeOnPhase(const double p_Mass);

    double          CalculateLuminosityAtBluePhaseEnd(const double p_Mass);
    double          CalculateLuminosityAtBluePhaseStart(const double p_Mass);

    double          CalculateLuminosityAtPhaseEnd()                              { return CalculateLuminosityAtBAGB(m_Mass0); }
    double          CalculateLuminosityOnPhase(const double p_Mass, const double p_Tau);
    double          CalculateLuminosityOnPhase()                                 { return CalculateLuminosityOnPhase(m_Mass0, m_Tau); }

    double          CalculateRadialExtentConvectiveEnvelope()                    { return BaseStar::CalculateRadialExtentConvectiveEnvelope(); }                // CHeB stars don't have a convective envelope

    double          CalculateRadiusAtBluePhaseEnd(const double p_Mass);
    double          CalculateRadiusAtBluePhaseStart(const double p_Mass);

    double          CalculateRadiusAtPhaseEnd(const double p_Mass, const double p_Luminosity);
    double          CalculateRadiusAtPhaseEnd()                                  { return CalculateRadiusAtPhaseEnd(m_Mass, m_Luminosity); }
    double          CalculateRadiusOnPhase(const double p_Mass, const double p_Luminosity, const double p_Tau);
    double          CalculateRadiusOnPhase()                                     { return CalculateRadiusOnPhase(m_Mass, m_Luminosity, m_Tau); }

    double          CalculateRadiusRho(const double p_Mass, const double p_Tau);

    double          CalculateRemnantLuminosity();
    double          CalculateRemnantRadius();

    double          CalculateTauAtPhaseEnd()                                     { return CalculateTauOnPhase(); }                                              // Same as on phase
    double          CalculateTauOnPhase();

    void            CalculateTimescales(const double p_Mass, DBL_VECTOR &p_Timescales);
    void            CalculateTimescales()                                        { CalculateTimescales(m_Mass0, m_Timescales); }                                // Use class member variables

    double          ChooseTimestep(const double p_Time);

    ENVELOPE        DetermineEnvelopeType()                                      { return ENVELOPE::CONVECTIVE; }                                               // Always CONVECTIVE (JR: should be RADIATIVE?  See http://gitlab.sr.bham.ac.uk/COMPAS/COMPAS/issues/135 and Hurley et al., 2002)
    ENVELOPE        DetermineEnvelopeTypeHurley2002()                            { return ENVELOPE::RADIATIVE; }                                                // Always RADIATIVE

    STELLAR_TYPE    EvolveToNextPhase();

    bool            IsEndOfPhase()                                               { return !ShouldEvolveOnPhase(); }                                             // Phase ends when age at or after He Burning
    bool            IsSupernova()                                                { return false; }                                                              // Not here

    STELLAR_TYPE    ResolveEnvelopeLoss(bool p_NoCheck = false);
    void            ResolveHeliumFlash() {  }                                                                                                                   // NO-OP
    STELLAR_TYPE    ResolveRemnantAfterEnvelopeLoss();

    bool            ShouldEvolveOnPhase()                                        { return (m_Age < (m_Timescales[static_cast<int>(TIMESCALE::tHeI)] + m_Timescales[static_cast<int>(TIMESCALE::tHe)])); }  // Evolve on CHeB phase if age after He Ign and while He Burning
    bool            ShouldSkipPhase()                                            { return false; }                                                              // Never skip CHeB phase

};

#endif // __CHeB_h__
