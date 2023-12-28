#ifndef __CHeB_h__
#define __CHeB_h__

#include "constants.h"
#include "typedefs.h"
#include "profiling.h"
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
                                                    const DBL_VECTOR &p_BnCoefficients) const;

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

    double          CalculateCOCoreMassOnPhase() const                          { return 0.0; }                                                                 // McCO(CHeB) = 0.0

    double          CalculateCoreMassAtPhaseEnd() const                         { return CalculateCoreMassAtBAGB(m_Mass0); }                                    // Use class member variables
    double          CalculateCoreMassOnPhase(const double p_Mass, const double p_Tau) const;
    double          CalculateCoreMassOnPhase() const                            { return CalculateCoreMassOnPhase(m_Mass0, m_Tau); }                            // Use class member variables

    double          CalculateCriticalMassRatioHurleyHjellmingWebbink() const    { return 0.33; }                                                                // As coded in BSE. Using the inverse owing to how qCrit is defined in COMPAS. See Hurley et al. 2002 sect. 2.6.1 for additional details.

    double          CalculateHeCoreMassAtPhaseEnd() const                       { return m_CoreMass; }

    double          CalculateLambdaDewi() const;
    double          CalculateLambdaNanjingStarTrack(const double p_Mass, const double p_Metallicity) const;
    double          CalculateLambdaNanjingEnhanced(const int p_MassInd, const int p_Zind) const;

    double          CalculateLifetimeOnBluePhase(const double p_Mass);
    double          CalculateLifetimeOnPhase(const double p_Mass);

    double          CalculateLuminosityAtBluePhaseEnd(const double p_Mass) const;
    double          CalculateLuminosityAtBluePhaseStart(const double p_Mass) const;

    double          CalculateLuminosityAtPhaseEnd() const                       { return CalculateLuminosityAtBAGB(m_Mass0); }
    double          CalculateLuminosityOnPhase(const double p_Mass, const double p_Tau) const;
    double          CalculateLuminosityOnPhase() const                          { return CalculateLuminosityOnPhase(m_Mass0, m_Tau); }

    double          CalculateRadialExtentConvectiveEnvelope() const             { return BaseStar::CalculateRadialExtentConvectiveEnvelope(); }                 // CHeB stars don't have a convective envelope

    double          CalculateRadiusAtBluePhaseEnd(const double p_Mass) const;
    double          CalculateRadiusAtBluePhaseStart(const double p_Mass) const;

    double          CalculateRadiusAtPhaseEnd(const double p_Mass, const double p_Luminosity) const;
    double          CalculateRadiusAtPhaseEnd() const                           { return CalculateRadiusAtPhaseEnd(m_Mass, m_Luminosity); }
    double          CalculateRadiusOnPhase(const double p_Mass, const double p_Luminosity, const double p_Tau) const;
    double          CalculateRadiusOnPhase() const                              { return CalculateRadiusOnPhase(m_Mass, m_Luminosity, m_Tau); }

    double          CalculateRadiusRho(const double p_Mass, const double p_Tau) const;

    double          CalculateRemnantLuminosity() const;
    double          CalculateRemnantRadius() const;

    double          CalculateTauAtPhaseEnd() const                              { return CalculateTauOnPhase(); }                                               // Same as on phase
    double          CalculateTauOnPhase() const;

    void            CalculateTimescales(const double p_Mass, DBL_VECTOR &p_Timescales);
    void            CalculateTimescales()                                       { CalculateTimescales(m_Mass0, m_Timescales); }                        // Use class member variables

    double          ChooseTimestep(const double p_Time) const;

    ENVELOPE        DetermineEnvelopeType() const;

    STELLAR_TYPE    EvolveToNextPhase();

    bool            IsEndOfPhase() const                                        { return !ShouldEvolveOnPhase(); }                                              // Phase ends when age at or after He Burning
    bool            IsSupernova() const                                         { return false; }                                                               // Not here

    STELLAR_TYPE    ResolveEnvelopeLoss(bool p_NoCheck = false);
    void            ResolveHeliumFlash() {  }                                                                                                                   // NO-OP

    bool            ShouldEnvelopeBeExpelledByPulsations() const { return ( OPTIONS->ExpelConvectiveEnvelopeAboveLuminosityThreshold() && DetermineEnvelopeType() == ENVELOPE::CONVECTIVE && utils::Compare( log10(m_Luminosity/m_Mass), OPTIONS->LuminosityToMassThreshold() ) >= 0 ) ; }                             // Envelope of convective star with luminosity to mass ratio beyond threshold should be expelled
    bool            ShouldEvolveOnPhase() const;
    bool            ShouldSkipPhase() const                                     { return false; }                                                               // Never skip CHeB phase

};

#endif // __CHeB_h__
