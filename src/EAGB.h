#ifndef __EAGB_h__
#define __EAGB_h__

#include "constants.h"
#include "typedefs.h"
#include "profiling.h"
#include "utils.h"

#include "CHeB.h"


class BaseStar;
class CHeB;

class EAGB: virtual public BaseStar, public CHeB {

public:

    EAGB(const BaseStar &p_BaseStar, const bool p_Initialise = true) : BaseStar(p_BaseStar), CHeB(p_BaseStar, false) {
        if (p_Initialise) Initialise();
    }

    EAGB& operator = (const BaseStar &p_BaseStar) {
        static_cast<BaseStar&>(*this) = p_BaseStar;
        Initialise();
        return *this;
    }


    // member functions
    static double CalculateRadiusOnPhase_Static(const double      p_Mass,
                                                const double      p_Luminosity,
                                                const double      p_MHeF,
                                                const DBL_VECTOR &p_BnCoefficients);


protected:

    void Initialise() {
        m_StellarType = STELLAR_TYPE::EARLY_ASYMPTOTIC_GIANT_BRANCH;                                                                                                // Set stellar type
        CalculateTimescales();                                                                                                                                      // Initialise timescales
        m_Age = m_Timescales[static_cast<int>(TIMESCALE::tHeI)] + m_Timescales[static_cast<int>(TIMESCALE::tHe)];                                                   // Set age appropriately
    }


    // member functions - alphabetically
    double          CalculateCOCoreMassAtPhaseEnd()                                             { return m_GBParams[static_cast<int>(GBP::McDU)]; }                 // McCO(EAGB) = McDU at phase end (Hurley et al. 2000, section 5.4)
    double          CalculateCOCoreMassOnPhase(const double p_Time);
    double          CalculateCOCoreMassOnPhase()                                                { return CalculateCOCoreMassOnPhase(m_Age); }                       // Use class member variables

    double          CalculateCoreMassAtPhaseEnd()                                               { return m_GBParams[static_cast<int>(GBP::McDU)]; }                 // Mc(EAGB) = McDU at phase end (Hurley et al. 2000, section 5.4)
    double          CalculateCoreMassOnPhase()                                                  { return m_GBParams[static_cast<int>(GBP::McBAGB)]; }               // Mc(EAGB) = McHe(EAGB) = McBAGB on phase (Hurley et al. 2000, section 5.4)

    double          CalculateGyrationRadius() const                                                   { return 0.1; }                                                     // Hurley et al., 2000, after eq 109 for giants. Single number approximation.   JR: todo: should this be in constants.h?

    double          CalculateHeCoreMassAtPhaseEnd()                                             { return CalculateHeCoreMassOnPhase(); }                            // Same as on phase
    double          CalculateHeCoreMassOnPhase()                                                { return m_HeCoreMass; }                                            // NO-OP

    double          CalculateInitialSupernovaMass()                                             { return m_GBParams[static_cast<int>(GBP::McBAGB)]; }               // For EAGB & TPAGB we use the mass at Base Asymptotic Giant Branch to determine SN type

    double          CalculateLambdaNanjing();

    double          CalculateLifetimeTo2ndDredgeUp(const double p_Tinf1_FAGB, const double p_Tinf2_FAGB);

    double          CalculateLuminosityAtPhaseEnd(const double p_CoreMass)                      { return CalculateLuminosityOnPhase(p_CoreMass); }                  // Same as on phase
    double          CalculateLuminosityAtPhaseEnd()                                             { return CalculateLuminosityAtPhaseEnd(m_COCoreMass); }             // Use class member variables
    double          CalculateLuminosityOnPhase(const double p_CoreMass);
    double          CalculateLuminosityOnPhase()                                                { return CalculateLuminosityOnPhase(m_COCoreMass); }

    double          CalculateMassLossRateHurley();

    double          CalculateRadialExtentConvectiveEnvelope()                                   { return GiantBranch::CalculateRadialExtentConvectiveEnvelope(); }  // Skip CHeB

    double          CalculateRadiusAtPhaseEnd(const double p_Mass, const double p_Luminosity)   { return CalculateRadiusOnPhase(p_Mass, p_Luminosity); }            // Same as on phase
    double          CalculateRadiusAtPhaseEnd()                                                 { return CalculateRadiusAtPhaseEnd(m_Mass, m_Luminosity); }         // Use class member variables
    double          CalculateRadiusOnPhase(const double p_Mass, const double p_Luminosity)      { return CalculateRadiusOnPhase_Static(p_Mass, p_Luminosity, m_MassCutoffs[static_cast<int>(MASS_CUTOFF::MHeF)], m_BnCoefficients); }
    double          CalculateRadiusOnPhase()                                                    { return CalculateRadiusOnPhase(m_Mass, m_Luminosity); }            // Use class member variables

    double          CalculateRemnantLuminosity();
    double          CalculateRemnantRadius();

    double          CalculateTauAtPhaseEnd()                                                    { return m_Tau; }                                                   // NO-OP
    double          CalculateTauOnPhase()                                                       { return m_Tau; }                                                   // NO-OP

    void            CalculateTimescales(const double p_Mass, DBL_VECTOR &p_Timescales);
    void            CalculateTimescales()                                                       { CalculateTimescales(m_Mass0, m_Timescales); }                     // Use class member variables

    double          ChooseTimestep(const double p_Time);

    ENVELOPE        DetermineEnvelopeType()                                                     { return ENVELOPE::CONVECTIVE; }                                    // Always CONVECTIVE

    STELLAR_TYPE    EvolveToNextPhase();

    bool            IsEndOfPhase()                                                              { return !ShouldEvolveOnPhase(); }                                  // Phase ends when age at or after DU timescale, and no TPAGB
    bool            IsSupernova();

    STELLAR_TYPE    ResolveEnvelopeLoss(bool p_NoCheck = false);
    void            ResolveHeliumFlash() {  }                                                                                                                       // NO-OP
    STELLAR_TYPE    ResolveRemnantAfterEnvelopeLoss();
    STELLAR_TYPE    ResolveSkippedPhase()                                                       { return m_StellarType; }                                           // NO-OP

    bool            ShouldEvolveOnPhase();
    bool            ShouldSkipPhase();

};

#endif // __EAGB_h__
