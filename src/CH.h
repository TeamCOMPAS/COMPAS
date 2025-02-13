#ifndef __CH_h__
#define __CH_h__

#include "constants.h"
#include "typedefs.h"
#include "profiling.h"
#include "utils.h"

#include "BaseStar.h"
#include "MS_gt_07.h"
#include "HeMS.h"

class BaseStar;
class MS_gt_07;

class CH: virtual public BaseStar, public MS_gt_07 {

public:

    CH() { m_StellarType = STELLAR_TYPE::CHEMICALLY_HOMOGENEOUS; };
    
    CH(const BaseStar &p_BaseStar, const bool p_Initialise = true) : BaseStar(p_BaseStar), MS_gt_07(p_BaseStar) {
        m_StellarType = STELLAR_TYPE::CHEMICALLY_HOMOGENEOUS;                                                                                                                       // Set stellar type
        if (p_Initialise) Initialise();                                                                                                                                             // Initialise if required
    }

    CH* Clone(const OBJECT_PERSISTENCE p_Persistence, const bool p_Initialise = true) {
        CH* clone = new CH(*this, p_Initialise); 
        clone->SetPersistence(p_Persistence); 
        return clone; 
    }

    static CH* Clone(CH& p_Star, const OBJECT_PERSISTENCE p_Persistence, const bool p_Initialise = true) {
        CH* clone = new CH(p_Star, p_Initialise); 
        clone->SetPersistence(p_Persistence); 
        return clone; 
    }


protected:

    void Initialise() {
        CalculateTimescales();                                                                                                                                                      // Initialise timescales
        m_Age = 0.0;                                                                                                                                                                // Set age appropriately
        m_CHE = true;                                                                                                                                                               // initially for CH stars                                                                                                                                                            // Set age appropriately
    }

    // member functions

    // Abundances
    double          CalculateHeliumAbundanceCoreOnPhase(const double p_Tau) const;
    double          CalculateHeliumAbundanceCoreOnPhase() const         { return CalculateHeliumAbundanceCoreOnPhase(m_Tau); };
    double          CalculateHeliumAbundanceSurfaceOnPhase(const double p_Tau) const;
    double          CalculateHeliumAbundanceSurfaceOnPhase() const      { return CalculateHeliumAbundanceSurfaceOnPhase(m_Tau); };

    double          CalculateHydrogenAbundanceCoreOnPhase(const double p_Tau) const;
    double          CalculateHydrogenAbundanceCoreOnPhase() const       { return CalculateHydrogenAbundanceCoreOnPhase(m_Tau); };
    double          CalculateHydrogenAbundanceSurfaceOnPhase(const double p_Tau) const;
    double          CalculateHydrogenAbundanceSurfaceOnPhase() const    { return CalculateHydrogenAbundanceSurfaceOnPhase(m_Tau); };
    
    // Lifetime
    double          CalculateLogLifetimeRatio(const double p_Mass) const;
    double          CalculateLifetimeRatio(const double p_Mass) const;

    // Luminosity
    double          CalculateLogLuminosityRatio(const double p_Mass, const double p_Tau) const;

    double          CalculateLuminosityAtPhaseEnd(const double p_Mass) const;
    double          CalculateLuminosityAtPhaseEnd() const               { return CalculateLuminosityAtPhaseEnd(m_Mass0); }                                                          // Use class member variables

    double          CalculateLuminosityOnPhase(const double p_Time, const double p_Mass, const double p_LZAMS) const;
    double          CalculateLuminosityOnPhase() const                  { return m_Luminosity; }    

    // Mass loss rate
    double          CalculateMassLossRateBelczynski2010();
    double          CalculateMassLossRateMerritt2024();
    double          CalculateMassLossRateWeightOB(const double p_HeliumAbundanceSurface);
    
    // Radius
    double          CalculateRadiusOnPhase() const                      { return m_RZAMS; }                                                                                         // Constant from birth
    double          CalculateRadiusAtPhaseEnd() const                   { return CalculateRadiusOnPhase(); }                                                                        // Same as on phase

    // Timescales
    void            CalculateTimescales(const double p_Mass, DBL_VECTOR &p_Timescales);
    void            CalculateTimescales()                               { return BaseStar::CalculateTimescales(); };  

    STELLAR_TYPE    EvolveToNextPhase();

    bool            ShouldEvolveOnPhase() const                         { return m_Age < m_Timescales[static_cast<int>(TIMESCALE::tMS)] && (OPTIONS->OptimisticCHE() || Omega() >= m_OmegaCHE); } // Evolve on CHE phase if age in MS timescale and spinning at least as fast as CHE threshold

    void            UpdateAgeAfterMassLoss();

};

#endif // __CH_h__
