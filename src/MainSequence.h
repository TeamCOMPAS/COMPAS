#ifndef __MainSequence_h__
#define __MainSequence_h__

#include "constants.h"
#include "typedefs.h"
#include "profiling.h"
#include "utils.h"

#include "BaseStar.h"

class BaseStar;

class MainSequence: virtual public BaseStar {

public:

    MainSequence(const BaseStar& baseStar) : BaseStar(baseStar) {}
    MainSequence& operator = (const BaseStar& baseStar) { static_cast<BaseStar&>(*this) = baseStar; return *this; }


protected:


    // member functions - alphabetically
    double          CalculateAlphaL(const double p_Mass) const;
    double          CalculateAlphaR(const double p_Mass) const;

    double          CalculateConvectiveEnvelopeMass() const                                 { return 0.0; }
    double          CalculateBetaL(const double p_Mass) const;
    double          CalculateBetaR(const double p_Mass) const;

    double          CalculateDeltaL(const double p_Mass) const;
    double          CalculateDeltaR(const double p_Mass) const;

    double          CalculateEta(const double p_Mass) const;

    double          CalculateGamma(const double p_Mass) const;

    double          CalculateCOCoreMassAtPhaseEnd() const                                   { return CalculateCOCoreMassOnPhase(); }                                // Same as on phase
    double          CalculateCOCoreMassOnPhase() const                                      { return 0.0; }                                                         // McCO(MS) = 0.0

    double          CalculateCoreMassAtPhaseEnd() const                                     { return OPTIONS->RetainCoreMassDuringCaseAMassTransfer() ? MinimumCoreMass() : 0.0; }                // Accounts for minimal core mass built up prior to mass loss through mass transfer
    double          CalculateCoreMassOnPhase() const                                        { return 0.0; }                                                         // Mc(MS) = 0.0 (Hurley et al. 2000, just before eq 28)

    double          CalculateGyrationRadius() const;

    double          CalculateHeCoreMassAtPhaseEnd() const                                   { return CalculateCoreMassAtPhaseEnd(); }                               // Same as He core mass
    double          CalculateHeCoreMassOnPhase() const                                      { return 0.0; }                                                         // McHe(MS) = 0.0

    double          CalculateLifetimeOnPhase(const double p_Mass, const double p_TBGB) const;

    double          CalculateLuminosityAtPhaseEnd(const double p_Mass) const;
    double          CalculateLuminosityAtPhaseEnd() const                                   { return CalculateLuminosityAtPhaseEnd(m_Mass0); }                      // Use class member variables
    double          CalculateLuminosityOnPhase(const double p_Time, const double p_Mass, const double p_LZAMS) const;
    double          CalculateLuminosityOnPhase() const                                      { return CalculateLuminosityOnPhase(m_Age, m_Mass0, m_LZAMS0); }        // Use class member variables

    double          CalculateMomentOfInertia(const double p_RemnantRadius = 0.0) const      { return CalculateGyrationRadius() * m_Mass * m_Radius * m_Radius; }
    double          CalculateMomentOfInertiaAU(const double p_RemnantRadius = 0.0) const    { return CalculateMomentOfInertia(p_RemnantRadius * RSOL_TO_AU) * RSOL_TO_AU * RSOL_TO_AU; }

    double          CalculatePerturbationMu() const                                         { return 5.0; }                                                         // mu(MS) = 5.0 (Hurley et al. 2000, eqs 97 & 98)

    double          CalculateRadialExtentConvectiveEnvelope() const;

    double          CalculateRadiusOnPhase(const double p_Mass, const double p_Time, const double p_RZAMS) const;
    double          CalculateRadiusAtPhaseEnd(const double p_Mass, const double p_RZAMS) const;
    double          CalculateRadiusAtPhaseEnd() const                                       { return CalculateRadiusAtPhaseEnd(m_Mass, m_RZAMS); }                  // Use class member variables
    double          CalculateRadiusOnPhase() const                                          { return CalculateRadiusOnPhase(m_Mass, m_Age, m_RZAMS0); }             // Use class member variables

    double          CalculateTauAtPhaseEnd() const                                          { return 1.0; }                                                         // tau = 1.0 at end of MS
    double          CalculateTauOnPhase() const;

    void            CalculateTimescales(const double p_Mass, DBL_VECTOR &p_Timescales);
    void            CalculateTimescales()                                                   { CalculateTimescales(m_Mass0, m_Timescales); }                         // Use class member variables

    double          CalculateZetaConstantsByEnvelope(ZETA_PRESCRIPTION p_ZetaPrescription)  { return OPTIONS->ZetaMainSequence(); }

    double          ChooseTimestep(const double p_Time) const;

    void            EvolveOneTimestepPreamble();
    STELLAR_TYPE    EvolveToNextPhase()                                                     { return STELLAR_TYPE::HERTZSPRUNG_GAP; }

    bool            IsEndOfPhase() const                                                    { return !ShouldEvolveOnPhase(); }                                      // Phase ends when age at or after MS timescale

    void            PerturbLuminosityAndRadius() { }                                                                                                                // NO-OP

    STELLAR_TYPE    ResolveEnvelopeLoss(bool p_NoCheck = false);

    bool            ShouldEvolveOnPhase() const                                             { return (m_Age < m_Timescales[static_cast<int>(TIMESCALE::tMS)]); }    // Evolve on MS phase if age in MS timescale

    void            UpdateInitialMass()                                                     { m_Mass0 = m_Mass; }                                                   // Per Hurley et al. 2000, section 7.1
    void            UpdateAgeAfterMassLoss();                                                                                                                       // Per Hurley et al. 2000, section 7.1
    
    void            UpdateMinimumCoreMass();                                                                                                                 // Set minimal core mass following Main Sequence mass transfer to MS age fraction of TAMS core mass

};

#endif // __MainSequence_h__
