#ifndef __Remnants_h__
#define __Remnants_h__

#include "constants.h"
#include "typedefs.h"
#include "profiling.h"
#include "utils.h"

#include "HeGB.h"


class BaseStar;
class HeGB;

class Remnants: virtual public BaseStar, public HeGB {

public:

    Remnants(const BaseStar &p_BaseStar, const bool p_Initialise = true) : BaseStar(p_BaseStar), HeGB(p_BaseStar, false) { }

//    Remnants& operator = (const BaseStar &p_BaseStar) { static_cast<BaseStar&>(*this) = p_BaseStar; return *this; }


    // member functions

protected:


    // member functions - alphabetically
    double          CalculateCOCoreMassOnPhase() const                                                          { return m_Mass; }                                                      // Return m_Mass

    double          CalculateConvectiveCoreRadius () const                                                      { return m_Radius; }                                                                                  // All core
    DBL_DBL         CalculateConvectiveEnvelopeMass() const                                                     { return std::tuple<double, double> (0.0, 0.0); }

    double          CalculateCoreMassOnPhase() const                                                            { return m_Mass; }                                                      // Return m_Mass

    double          CalculateCriticalMassRatio(const bool p_AccretorIsDegenerate)                               { return 0.0; }                                                         // Should never be called...

    void            CalculateGBParams(const double p_Mass, DBL_VECTOR &p_GBParams)                              { GiantBranch::CalculateGBParams(p_Mass, p_GBParams); }                 // Default to GiantBranch  
    void            CalculateGBParams()                                                                         { CalculateGBParams(m_Mass0, m_GBParams); }                             // Use class member variables

    double          CalculateHeCoreMassOnPhase() const                                                          { return m_Mass; }                                                      // Return m_Mass

    double          CalculateInitialSupernovaMass() const                                                       { return GiantBranch::CalculateInitialSupernovaMass(); }                // Use GiantBranch

    double          CalculateLambdaDewi() const                                                                 { return BaseStar::CalculateLambdaDewi(); }                             // Not supported - use BaseStar
    double          CalculateLambdaNanjingStarTrack(const double p_Mass, const double p_Metallicity) const      { return BaseStar::CalculateLambdaNanjingStarTrack(0.0, 0.0); }         // Not supported - use BaseStar (0.0 are dummy values)      JR: todo: check this (type 10 not mentioned as not supported in original code)

    DBL_DBL         CalculateMassAcceptanceRate(const double p_DonorMassRate,
                                                const double p_AccretorMassRate);
    DBL_DBL         CalculateMassAcceptanceRate(const double p_DonorMassRate,
                                                const double p_AccretorMassRate,
                                                const bool   p_IsHeRich)                                        { return CalculateMassAcceptanceRate(p_DonorMassRate, p_AccretorMassRate); } // Ignore the He content for non-WDs

    double          CalculateMassLossRateHurley()                                                               { return 0.0; }
    double          CalculateMassLossRateBelczynski2010()                                                                 { return 0.0; }

    double          CalculatePerturbationMuOnPhase() const                                                      { return m_Mu; }                                                        // NO-OP

    double          CalculateRadialExtentConvectiveEnvelope() const                                             { return 0.0; }         // WD stars don't have a convective envelope

    std::tuple <double, STELLAR_TYPE> CalculateRadiusAndStellarTypeOnPhase() const                              { return BaseStar::CalculateRadiusAndStellarTypeOnPhase(); }

    double          CalculateTauOnPhase() const                                                                 { return m_Tau; }                                                       // NO-OP
   
    double          CalculateThermalTimescale(const double p_Radius = 1.0) const                                { return CalculateDynamicalTimescale(); }                               // Parameter is ignored
    double          CalculateThermalTimescale() const                                                           { return CalculateThermalTimescale(m_Radius); }                         // Use inheritance hierarchy

    double          CalculateThermalMassLossRate() const                                                        { return BaseStar::CalculateThermalMassLossRate(); }                    // Set thermal mass gain rate to be effectively infinite, using dynamical timescale (in practice, will be Eddington limited), avoid division by zero

    void            CalculateTimescales(const double p_Mass, DBL_VECTOR &p_Timescales)                          { return TPAGB::CalculateTimescales(p_Mass, p_Timescales); }            // Use TPAGB
    void            CalculateTimescales()                                                                       { CalculateTimescales(m_Mass0, m_Timescales); }                         // Use class member variables

    double          CalculateZetaConstantsByEnvelope(ZETA_PRESCRIPTION p_ZetaPrescription)                      { return 0.0; }                                                         // Should never be called...

    double          ChooseTimestep(const double p_Time) const;
    
    ENVELOPE        DetermineEnvelopeType() const                                                               { return ENVELOPE::REMNANT; }                                           // Always REMNANT

    void            EvolveOneTimestepPreamble()                                                                 { BaseStar::EvolveOneTimestepPreamble(); }                              // Default to BaseStar

    STELLAR_TYPE    EvolveToNextPhase()                                                                         { return BaseStar::EvolveToNextPhase(); }                               // Default to BaseStar

    bool            IsDegenerate() const                                                                        { return true; }                                                        // White Dwarfs, NS and BH are degenerate

    bool            IsEndOfPhase() const                                                                        { return !ShouldEvolveOnPhase(); }                                      // Phase ends when envelope loss or going supernova

    bool            IsSupernova() const                                                                         { return false; }                                                       // Default

    void            PerturbLuminosityAndRadius() { }                                                                                                                                    // NO-OP

    STELLAR_TYPE    ResolveEnvelopeLoss(bool p_NoCheck = false)                                                 { return BaseStar::ResolveEnvelopeLoss(p_NoCheck); }                    // Default to BaseStar

    void            ResolveEnvelopeMassAtPhaseEnd(const double p_Tau) const                                     { ResolveEnvelopeMassOnPhase(p_Tau); }                                  // Same as on phase
    void            ResolveEnvelopeMassOnPhase(const double p_Tau) const { }                                                                                                            // NO-OP

    void            ResolveMassLoss(const bool p_UpdateMDt = true) { }                                                                                                                                               // NO-OP

    STELLAR_TYPE    ResolveSkippedPhase()                                                                       { return BaseStar::ResolveSkippedPhase(); }                             // Default to BaseStar
                                                                                                                                                                                        //
    void            ResolveShellChange(const double p_AccretedMass) { }                                                                                                                 // NO-OP 
                                                                                                                                                                                        //
    STELLAR_TYPE    ResolveSupernova()                                                                          { return BaseStar::ResolveSupernova(); }                                // Default to BaseStar

    void            SetPulsarParameters() const { }                                                                                                                                     // NO-OP

    bool            ShouldEnvelopeBeExpelledByPulsations() const                                                { return false; }                                                       // No envelope to lose by pulsations
    bool            ShouldEvolveOnPhase() const                                                                 { return true; }                                                        // Default
    bool            ShouldSkipPhase() const                                                                     { return false; }                                                       // Don't skip WD phase

};

#endif // __Remnants_h__
