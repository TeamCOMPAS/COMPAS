#ifndef __FGB_h__
#define __FGB_h__

#include "constants.h"
#include "typedefs.h"
#include "profiling.h"
#include "utils.h"

#include "HG.h"


class BaseStar;
class HG;

class FGB: virtual public BaseStar, public HG {

public:

    FGB(const BaseStar &p_BaseStar, const bool p_Initialise = true) : BaseStar(p_BaseStar), HG(p_BaseStar, false) {
        m_StellarType = STELLAR_TYPE::FIRST_GIANT_BRANCH;                                                                                                                                           // Set stellar type
        if (p_Initialise) Initialise();                                                                                                                                                             // Initialise if required
    }

    FGB* Clone(const OBJECT_PERSISTENCE p_Persistence, const bool p_Initialise = true) {
        FGB* clone = new FGB(*this, p_Initialise); 
        clone->SetPersistence(p_Persistence); 
        return clone; 
    }

    static FGB* Clone(FGB& p_Star, const OBJECT_PERSISTENCE p_Persistence, const bool p_Initialise = true) {
        FGB* clone = new FGB(p_Star, p_Initialise); 
        clone->SetPersistence(p_Persistence); 
        return clone; 
    }


protected:

    void Initialise() {
        CalculateTimescales();                                                                                                                                                                      // Initialise timescales
        m_Age = m_Timescales[static_cast<int>(TIMESCALE::tBGB)];                                                                                                                                    // Set age appropriately
    }


    // member functions - alphabetically
    double          CalculateCOCoreMassAtPhaseEnd() const                                           { return CalculateCOCoreMassOnPhase(); }                                                        // Same as on phase
    double          CalculateCOCoreMassOnPhase() const                                              { return 0.0; }                                                                                 // McCO(FGB) = 0.0

    double          CalculateCoreMassAtPhaseEnd(const double p_Mass, const double p_Time) const     { return CalculateCoreMassOnPhase(p_Mass, p_Time); }                                            // Same as on phase
    double          CalculateCoreMassAtPhaseEnd() const                                             { return CalculateCoreMassAtPhaseEnd(m_Mass0, m_Age); }                                         // Use class member variables
    double          CalculateCoreMassOnPhase(const double p_Mass, const double p_Time) const;
    double          CalculateCoreMassOnPhase() const                                                { return CalculateCoreMassOnPhase(m_Mass0, m_Age); }                                            // Use class member variables

    double          CalculateCriticalMassRatioClaeys14(const bool p_AccretorIsDegenerate) const     { return GiantBranch::CalculateCriticalMassRatioClaeys14(p_AccretorIsDegenerate); }             // Skip HG 
    double          CalculateCriticalMassRatioHurleyHjellmingWebbink() const                        { return GiantBranch::CalculateCriticalMassRatioHurleyHjellmingWebbink(); }

    double          CalculateHeCoreMassAtPhaseEnd() const                                           { return CalculateHeCoreMassOnPhase(); }                                                        // Same as on phase
    double          CalculateHeCoreMassOnPhase() const                                              { return m_CoreMass; }                                                                          // McHe(FGB) = Core Mass

    double          CalculateLuminosityAtPhaseEnd(const double p_Time) const                        { return CalculateLuminosityOnPhase(p_Time); }                                                  // Same as on phase
    double          CalculateLuminosityAtPhaseEnd() const                                           { return CalculateLuminosityAtPhaseEnd(m_Age); }                                                // Use class member variables
    double          CalculateLuminosityOnPhase(const double p_Time) const;
    double          CalculateLuminosityOnPhase() const                                              { return CalculateLuminosityOnPhase(m_Age); }                                                   // Use class member variables

    double          CalculateRadialExtentConvectiveEnvelope() const                                 { return m_Radius - CalculateConvectiveCoreRadius(); }                                          // Skip HG

    double          CalculateRadiusAtPhaseEnd(const double p_Mass, const double p_Luminosity) const { return GiantBranch::CalculateRadiusOnPhase(p_Mass, p_Luminosity); }                  // Skip HG - same as on phase
    double          CalculateRadiusAtPhaseEnd() const                                               { return CalculateRadiusAtPhaseEnd(m_Mass, m_Luminosity); }                                     // Use class member variables
    double          CalculateRadiusOnPhase(const double p_Mass, const double p_Luminosity) const    { return GiantBranch::CalculateRadiusOnPhase(p_Mass, p_Luminosity); }                           // Skip HG
    double          CalculateRadiusOnPhase() const                                                  { return CalculateRadiusOnPhase(m_Mass, m_Luminosity); }                                        // Use class member variables

    double          CalculateTauAtPhaseEnd() const                                                  { return m_Tau; }                                                                               // NO-OP
    double          CalculateTauOnPhase() const;

    double          ChooseTimestep(const double p_Time) const;

    ENVELOPE        DetermineEnvelopeType() const                                                   { return ENVELOPE::CONVECTIVE; }                                                                // Always CONVECTIVE

    STELLAR_TYPE    EvolveToNextPhase();

    bool            IsEndOfPhase() const                                                            { return !ShouldEvolveOnPhase(); }                                                              // Phase ends when age at or after He ignition timescale
    bool            IsSupernova() const                                                             { return false; }                                                                               // Not here

    STELLAR_TYPE    ResolveEnvelopeLoss(bool p_NoCheck = false);
    void            ResolveHeliumFlash();
    STELLAR_TYPE    ResolveSkippedPhase()                                                           { return STELLAR_TYPE::CORE_HELIUM_BURNING; }                                                   // Evolve to CHeB if phase is skipped

    bool            ShouldEvolveOnPhase() const                                                     { return (utils::Compare(m_Age, m_Timescales[static_cast<int>(TIMESCALE::tHeI)]) < 0); }        // Evolve on FGB phase if age < He ignition timescale
    bool            ShouldSkipPhase() const                                                         { return (utils::Compare(m_Mass0, m_MassCutoffs[static_cast<int>(MASS_CUTOFF::MFGB)]) >= 0); }  // Skip phase if mass >= FGB mass cutoff

    void            UpdateAgeAfterMassLoss()                                                        { GiantBranch::UpdateAgeAfterMassLoss(); }                                                      // Skip HG
    void            UpdateInitialMass()                                                             { GiantBranch::UpdateInitialMass(); }                                                           // Skip HG

};

#endif // __FGB_h__
