#ifndef __MainSequence_h__
#define __MainSequence_h__

#include "constants.h"
#include "typedefs.h"
#include "utils.h"

#include "BaseStar.h"

class BaseStar;

class MainSequence: virtual public BaseStar {

public:

    MainSequence(const BaseStar& baseStar) : BaseStar(baseStar) {}
    MainSequence& operator = (const BaseStar& baseStar) { static_cast<BaseStar&>(*this) = baseStar; return *this; }


protected:


    // member functions - alphabetically
    double          CalculateAlphaL(const double p_Mass);
    double          CalculateAlphaR(const double p_Mass);

    double          CalculateBetaL(const double p_Mass);
    double          CalculateBetaR(const double p_Mass);

    double          CalculateDeltaL(const double p_Mass);
    double          CalculateDeltaR(const double p_Mass);

    double          CalculateEta(const double p_Mass);

    double          CalculateGamma(const double p_Mass);

    double          CalculateCOCoreMassAtPhaseEnd()                                 { return CalculateCOCoreMassOnPhase(); }                                        // Same as on phase
    double          CalculateCOCoreMassOnPhase()                                    { return 0.0; }                                                                 // McCO(MS) = 0.0

    double          CalculateCoreMassAtPhaseEnd()                                   { return CalculateCoreMassOnPhase(); }                                          // Same as on phase
    double          CalculateCoreMassOnPhase()                                      { return 0.0; }                                                                 // Mc(MS) = 0.0 (Hurley et al. 2000, just before eq 28)

    double          CalculateEnvelopeMassOnPhase(const double p_Tau);

    double          CalculateGyrationRadius();

    double          CalculateHeCoreMassAtPhaseEnd()                                 { return CalculateHeCoreMassOnPhase(); }                                        // Same as on phase
    double          CalculateHeCoreMassOnPhase()                                    { return 0.0; }                                                                 // McHe(MS) = 0.0

    double          CalculateLifetimeOnPhase(const double p_Mass, const double p_TBGB);

    double          CalculateLuminosityAtPhaseEnd(const double p_Mass);
    double          CalculateLuminosityAtPhaseEnd()                                 { return CalculateLuminosityAtPhaseEnd(m_Mass0); }                              // Use class member variables
    double          CalculateLuminosityOnPhase(const double p_Time, const double p_Mass, const double p_LZAMS);
    double          CalculateLuminosityOnPhase()                                    { return CalculateLuminosityOnPhase(m_Age, m_Mass0, m_LZAMS0); }                // Use class member variables

    double          CalculateMomentOfInertia(const double p_RemnantRadius = 0.0)    { return CalculateGyrationRadius() * m_Mass * m_Radius * m_Radius; }
    double          CalculateMomentOfInertiaAU(const double p_RemnantRadius = 0.0)  { return CalculateMomentOfInertia(p_RemnantRadius * RSOL_TO_AU) * RSOL_TO_AU * RSOL_TO_AU; }

    double          CalculatePerturbationMu()                                       { return 5.0; }                                                                 // mu(MS) = 5.0 (Hurley et al. 2000, eqs 97 & 98)

    double          CalculateRadialExtentConvectiveEnvelope();

    double          CalculateRadiusOnPhase(const double p_Mass, const double p_Time, const double p_RZAMS);
    double          CalculateRadiusAtPhaseEnd(const double p_Mass, const double p_RZAMS);
    double          CalculateRadiusAtPhaseEnd()                                     { return CalculateRadiusAtPhaseEnd(m_Mass, m_RZAMS); }                          // Use class member variables
    double          CalculateRadiusOnPhase()                                        { return CalculateRadiusOnPhase(m_Mass, m_Age, m_RZAMS0); }                     // Use class member variables

    double          CalculateTauAtPhaseEnd()                                        { return 1.0; }                                                                 // tau = 1.0 at end of MS
    double          CalculateTauOnPhase();

    double          CalculateThermalTimescale(const double p_Mass, const double p_Radius, const double p_Luminosity, const double p_EnvMass = 1.0);
    double          CalculateThermalTimescale()                                     { return CalculateThermalTimescale(m_Mass, m_Radius, m_Luminosity); }           // Use class member variables

    void            CalculateTimescales(const double p_Mass, DBL_VECTOR &p_Timescales);
    void            CalculateTimescales()                                           { CalculateTimescales(m_Mass0, m_Timescales); }                                 // Use class member variables

    double          CalculateZeta(CE_ZETA_PRESCRIPTION p_CEZetaPrescription)        { m_Error = ERROR::INVALID_TYPE_ZETA_CALCULATION;                               // Set error value
                                                                                      SHOW_WARN(m_Error);                                                           // Warn that an error occurred
                                                                                      return 0.0; }                                                                 // Should never be called...

    double          ChooseTimestep(const double p_Time);

    MT_CASE         DetermineMassTransferCase()                                     { return MT_CASE::A; }                                                          // Mass Transfer Case A for MS stars

    void            EvolveOneTimestepPreamble();
    STELLAR_TYPE    EvolveToNextPhase()                                             { return STELLAR_TYPE::HERTZSPRUNG_GAP; }

    bool            IsEndOfPhase()                                                  { return !ShouldEvolveOnPhase(); }                                              // Phase ends when age at or after MS timescale

    void            PerturbLuminosityAndRadius() { }                                                                                                                // NO-OP

    STELLAR_TYPE    ResolveEnvelopeLoss(bool p_NoCheck = false)                     { return m_StellarType; }                                                       // NO-OP

    STELLAR_TYPE    ResolveRemnantAfterEnvelopeLoss();

    bool            ShouldEvolveOnPhase()                                           { return (m_Age < m_Timescales[static_cast<int>(TIMESCALE::tMS)]); }            // Evolve on MS phase if age in MS timescale

    void            UpdateInitialMass()                                             { m_Mass0 = m_Mass; }                                                           // Per Hurley et al. 2000, section 7.1
    void            UpdateAgeAfterMassLoss();                                                                                                                       // Per Hurley et al. 2000, section 7.1

};

#endif // __MainSequence_h__
