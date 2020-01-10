#ifndef __FGB_h__
#define __FGB_h__

#include "constants.h"
#include "typedefs.h"
#include "utils.h"

#include "HG.h"


class BaseStar;
class HG;

class FGB: virtual public BaseStar, public HG {

public:

    FGB(const BaseStar &p_BaseStar, const bool p_Initialise = true) : BaseStar(p_BaseStar), HG(p_BaseStar, false) {
        if (p_Initialise) Initialise();
    }

    FGB& operator = (const BaseStar &p_BaseStar) {
        static_cast<BaseStar&>(*this) = p_BaseStar;
        Initialise();
        return *this;
    }


protected:

    void Initialise() {
        m_StellarType = STELLAR_TYPE::FIRST_GIANT_BRANCH;                                                                                                                                       // Set stellar type
        CalculateTimescales();                                                                                                                                                                  // Initialise timescales
        m_Age = m_Timescales[static_cast<int>(TIMESCALE::tBGB)];                                                                                                                                // Set age appropriately
    }


    // member functions - alphabetically
    double          CalculateCOCoreMassAtPhaseEnd()                                             { return CalculateCOCoreMassOnPhase(); }                                                        // Same as on phase
    double          CalculateCOCoreMassOnPhase()                                                { return 0.0; }                                                                                 // McCO(FGB) = 0.0

    double          CalculateConvergedMassStepZetaThermal();

    double          CalculateCoreMassAtPhaseEnd(const double p_Mass, const double p_Time)       { return CalculateCoreMassOnPhase(p_Mass, p_Time); }                                            // Same as on phase
    double          CalculateCoreMassAtPhaseEnd()                                               { return CalculateCoreMassAtPhaseEnd(m_Mass0, m_Age); }                                         // Use class member variables
    double          CalculateCoreMassOnPhase(const double p_Mass, const double p_Time);
    double          CalculateCoreMassOnPhase()                                                  { return CalculateCoreMassOnPhase(m_Mass0, m_Age); }                                            // Use class member variables

    double          CalculateEnvelopeMassOnPhase(const double p_Tau)                            { return m_EnvMass; }                                                                           // NO-OP

    double          CalculateGyrationRadius()                                                   { return 0.1; }                                                                                 // Hurley et al., 2000, after eq 109 for giants. Single number approximation.

    double          CalculateHeCoreMassOnPhase()                                                { return m_CoreMass; }                                                                          // McHe(FGB) = Core Mass
    double          CalculateHeCoreMassAtPhaseEnd()                                             { return CalculateHeCoreMassOnPhase(); }                                                        // Same as on phase

    double          CalculateLuminosityAtPhaseEnd(const double p_Time)                          { return CalculateLuminosityOnPhase(p_Time); }                                                  // Same as on phase
    double          CalculateLuminosityAtPhaseEnd()                                             { return CalculateLuminosityAtPhaseEnd(m_Age); }                                                // Use class member variables
    double          CalculateLuminosityOnPhase(const double p_Time);
    double          CalculateLuminosityOnPhase()                                                { return CalculateLuminosityOnPhase(m_Age); }                                                   // Use class member variables

    double          CalculateRadialExtentConvectiveEnvelope()                                   { return GiantBranch::CalculateRadialExtentConvectiveEnvelope(); }                              // Skip HG

    double          CalculateRadiusAtPhaseEnd(const double p_Mass, const double p_Luminosity)   { return GiantBranch::CalculateRadiusOnPhase(p_Mass, p_Luminosity); }                           // Skip HG - same as on phase
    double          CalculateRadiusAtPhaseEnd()                                                 { return CalculateRadiusAtPhaseEnd(m_Mass, m_Luminosity); }                                     // Use class member variables
    double          CalculateRadiusOnPhase(const double p_Mass, const double p_Luminosity)      { return GiantBranch::CalculateRadiusOnPhase(p_Mass, p_Luminosity); }                           // Skip HG
    double          CalculateRadiusOnPhase()                                                    { return CalculateRadiusOnPhase(m_Mass, m_Luminosity); }                                        // Use class member variables

    double          CalculateTauAtPhaseEnd()                                                    { return m_Tau; }                                                                               // NO-OP
    double          CalculateTauOnPhase();

    double          ChooseTimestep(const double p_Time);

    ENVELOPE        DetermineEnvelopeType()                                                     { return ENVELOPE::CONVECTIVE; }                                                                // Always CONVECTIVE
    ENVELOPE        DetermineEnvelopeTypeHurley2002()                                           { return ENVELOPE::CONVECTIVE; }                                                                // Always CONVECTIVE

    STELLAR_TYPE    EvolveToNextPhase();

    bool            IsEndOfPhase()                                                              { return !ShouldEvolveOnPhase(); }                                                              // Phase ends when age at or after He ignition timescale
    bool            IsMassRatioUnstable(const double p_AccretorMass,
                                        const bool p_AccretorIsDegenerate)                      { return GiantBranch::IsMassRatioUnstable(p_AccretorMass, p_AccretorIsDegenerate); }            // Skip HG
    bool            IsSupernova()                                                               { return false; }                                                                               // Not here

    STELLAR_TYPE    ResolveEnvelopeLoss(bool p_NoCheck = false);
    void            ResolveHeliumFlash();
    STELLAR_TYPE    ResolveSkippedPhase()                                                       { return STELLAR_TYPE::CORE_HELIUM_BURNING; }                                                   // Evolve to CHeB if phase is skipped

    bool            ShouldEvolveOnPhase()                                                       { return (utils::Compare(m_Age, m_Timescales[static_cast<int>(TIMESCALE::tHeI)]) < 0); }        // Evolve on FGB phase if age < He ignition timescale
    bool            ShouldSkipPhase()                                                           { return (utils::Compare(m_Mass0, m_MassCutoffs[static_cast<int>(MASS_CUTOFF::MFGB)]) >= 0); }  // Skip phase if mass >= FGB mass cutoff

    void            UpdateAgeAfterMassLoss()                                                    { GiantBranch::UpdateAgeAfterMassLoss(); }                                                      // Skip HG
    void            UpdateInitialMass()                                                         { GiantBranch::UpdateInitialMass(); }                                                           // Skip HG

};

#endif // __FGB_h__
