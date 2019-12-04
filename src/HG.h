#ifndef __HG_h__
#define __HG_h__

#include "constants.h"
#include "typedefs.h"
#include "utils.h"

#include "GiantBranch.h"


// JR: todo: revisit this one day - sometimes HG works better as GiantBranch, sometimes not...
// Right now it is GiantBranch - figure out which is more appropriate

class BaseStar;
class GiantBranch;

class HG: virtual public BaseStar, public GiantBranch {

public:

    HG(const BaseStar &p_BaseStar, const bool p_Initialise = true) : BaseStar(p_BaseStar), GiantBranch(p_BaseStar) {
        if (p_Initialise) Initialise();
    }

    HG& operator = (const BaseStar &p_BaseStar) {
        static_cast<BaseStar&>(*this) = p_BaseStar;
        Initialise();
        return *this;
    }


protected:

    void Initialise() {
        m_StellarType = STELLAR_TYPE::HERTZSPRUNG_GAP;                                                                                                                      // Set stellar type
        CalculateTimescales();                                                                                                                                              // Initialise timescales
        m_Age = m_Timescales[static_cast<int>(TIMESCALE::tMS)];                                                                                                             // Set age appropriately
    }


    // member functions - alphabetically
    double       CalculateCOCoreMassAtPhaseEnd()                                { return 0.0; }                                                                             // McCO(HG) = 0.0
    double       CalculateCOCoreMassOnPhase()                                   { return 0.0; }                                                                             // McCO(HG) = 0.0

    double       CalculateConvergedMassStepZetaThermal()                        { return BaseStar::CalculateConvergedMassStepZetaThermal(); }                               // Not GiantBranch this time

    double       CalculateCoreMassAt2ndDredgeUp(const DBL_VECTOR &p_GBParams)   { return p_GBParams[static_cast<int>(GBP::McDU)]; }                                         // NO-OP
    double       CalculateCoreMassAtPhaseEnd(const double p_Mass);
    double       CalculateCoreMassAtPhaseEnd()                                  { return CalculateCoreMassAtPhaseEnd(m_Mass0); }                                            // Use class member variables
    double       CalculateCoreMassOnPhase(const double p_Mass, const double p_Time);
    double       CalculateCoreMassOnPhase()                                     { return CalculateCoreMassOnPhase(m_Mass0, m_Age); }                                        // Use class member variables

    double       CalculateEnvelopeMassOnPhase(const double p_Tau);

    double       CalculateGyrationRadius()                                      { return 0.21; }                                                                            // Hurley et al., 2000, after eq 109 for n=3/2 polytrope or dense convective core. Single number approximation.

    double       CalculateHeCoreMassAtPhaseEnd()                                { return m_CoreMass; }                                                                      // McHe(HG) = Core Mass
    double       CalculateHeCoreMassOnPhase()                                   { return m_CoreMass; }                                                                      // McHe(HG) = Core Mass

    double       CalculateLambdaDewi();
    double       CalculateLambdaNanjing();

    double       CalculateLuminosityAtPhaseEnd(const double p_Mass);
    double       CalculateLuminosityAtPhaseEnd()                                { return CalculateLuminosityAtPhaseEnd(m_Mass0);}                                           // Use class member variables
    double       CalculateLuminosityOnPhase(const double p_Age, const double p_Mass);
    double       CalculateLuminosityOnPhase()                                   { return CalculateLuminosityOnPhase(m_Age, m_Mass0); }                                      // Use class member variables

    double       CalculateMassTransferRejuvenationFactor();

    double       CalculateRadialExtentConvectiveEnvelope();

    double       CalculateRadiusAtPhaseEnd(const double p_Mass);
    double       CalculateRadiusAtPhaseEnd ()                                   { return CalculateRadiusAtPhaseEnd(m_Mass); }                                               // Use class member variables
    double       CalculateRadiusOnPhase(const double p_Mass, const double p_Tau, const double p_RZAMS);
    double       CalculateRadiusOnPhase()                                       { return CalculateRadiusOnPhase(m_Mass, m_Tau, m_RZAMS0); }                                 // Use class member variables

    double       CalculateRho(const double p_Mass);

    double       CalculateTauAtPhaseEnd()                                       { return 1.0; }                                                                             // tau = 1.0 at end of HG
    double       CalculateTauOnPhase();

    double       ChooseTimestep(const double p_Time);

    ENVELOPE     DetermineEnvelopeType();
    ENVELOPE     DetermineEnvelopeTypeHurley2002()                              { return ENVELOPE::CONVECTIVE; }                                                            // Always CONVECTIVE

    MT_CASE      DetermineMassTransferCase()                                    { return MT_CASE::B; }                                                                      // Mass Transfer Case B for HG stars

    void         EvolveOneTimestepPreamble()                                    { BaseStar::EvolveOneTimestepPreamble(); }                                                  // Skip MainSequence
    STELLAR_TYPE EvolveToNextPhase();

    bool         IsEndOfPhase()                                                 { return !ShouldEvolveOnPhase(); }                                                          // Phase ends when age at or after Base Giant Branch MS timescale
    bool         IsMassRatioUnstable(const double p_AccretorMass, const bool p_AccretorIsDegenerate);
    bool         IsSupernova()                                                  { return false; }                                                                           // Not here

    STELLAR_TYPE ResolveEnvelopeLoss(bool p_NoCheck = false);
    void         ResolveHeliumFlash() {  }                                                                                                                                  // NO-OP
    STELLAR_TYPE ResolveRemnantAfterEnvelopeLoss();
    STELLAR_TYPE ResolveSkippedPhase()                                          { return m_StellarType; }                                                                   // NO-OP

    bool         ShouldEvolveOnPhase()                                          { return (utils::Compare(m_Age, m_Timescales[static_cast<int>(TIMESCALE::tBGB)]) < 0); }    // Evolve on HG phase if age < Base Giant Branch timescale
    bool         ShouldSkipPhase()                                              { return false; }                                                                           // Never skip HG phase

    void         UpdateAgeAfterMassLoss();                                                                                                                                  // Per Hurley et al. 2000, section 7.1

    void         UpdateInitialMass()                                            { if (m_Mass0 > m_CoreMass) m_Mass0 = m_Mass; }                                             // Per Hurley et al. 2000, section 7.1

};

#endif // __HG_h__
