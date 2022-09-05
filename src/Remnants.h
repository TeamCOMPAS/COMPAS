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

    Remnants(const BaseStar &p_BaseStar, const bool p_Initialise = true) : BaseStar(p_BaseStar), HeGB(p_BaseStar, false) {
        if (p_Initialise) Initialise();
    }

    Remnants& operator = (const BaseStar &p_BaseStar) { static_cast<BaseStar&>(*this) = p_BaseStar; return *this; }


    // member functions


protected:


    // member functions - alphabetically
    double          CalculateCOCoreMassOnPhase() const                                                          { return m_Mass; }                                                      // Return m_Mass

    double          CalculateCoreMassOnPhase() const                                                            { return m_Mass; }                                                      // Return m_Mass

    double          CalculateCriticalMassRatio(const bool p_AccretorIsDegenerate) const                         { return 0.0; }                                                         // Should never be called...

    double          CalculateEddingtonCriticalRate() const                                                      { return 2.08E-3 / 1.7 * m_Radius * MYR_TO_YEAR; }                      // Hurley+, 2002, Eq. (67)

    void            CalculateGBParams()                                                                         { GiantBranch::CalculateGBParams(); }                                   // Default to GiantBranch

    double          CalculateGyrationRadius() const                                                             { return 0.21; }                                                        // Hurley et al., 2000, after eq 109 for n=3/2 polytrope or dense convective core. Single number approximation.

    double          CalculateHeCoreMassOnPhase() const                                                          { return m_Mass; }                                                      // Return m_Mass

    double          CalculateInitialSupernovaMass() const                                                       { return GiantBranch::CalculateInitialSupernovaMass(); }                // Use GiantBranch

    double          CalculateLambdaDewi() const                                                                 { return BaseStar::CalculateLambdaDewi(); }                             // Not supported - use BaseStar
    double          CalculateLambdaNanjingStarTrack(const double p_Mass, const double p_Metallicity) const               { return BaseStar::CalculateLambdaNanjingStarTrack(0.0, 0.0); }                  // Not supported - use BaseStar (0.0 are dummy values)      JR: todo: check this (type 10 not mentioned as not supported in original code)

    DBL_DBL         CalculateMassAcceptanceRate(const double p_DonorMassRate,
                                                    const double p_AccretorMassRate = 0.0);

    double          CalculateMassLossRateHurley();
    double          CalculateMassLossRateVink()                                                                 { return 0.0; }

    double          CalculateMomentOfInertia(const double p_RemnantRadius = 0.0) const                          { return GiantBranch::CalculateMomentOfInertia(p_RemnantRadius); }      // Default to GiantBranch
    double          CalculateMomentOfInertiaAU(const double p_RemnantRadius = 0.0) const                        { return GiantBranch::CalculateMomentOfInertiaAU(p_RemnantRadius); }    // Default to GiantBranch

    double          CalculatePerturbationMuOnPhase() const                                                      { return m_Mu; }                                                        // NO-OP

    double          CalculateRadialExtentConvectiveEnvelope() const                                             { return BaseStar::CalculateRadialExtentConvectiveEnvelope(); }         // WD stars don't have a convective envelope

    std::tuple <double, STELLAR_TYPE> CalculateRadiusAndStellarTypeOnPhase() const                              { return BaseStar::CalculateRadiusAndStellarTypeOnPhase(); }

    double          CalculateTauOnPhase() const                                                                 { return m_Tau; }                                                       // NO-OP
   
    double          CalculateThermalTimescale(const double p_Radius = 1.0) const                                { return CalculateDynamicalTimescale(); }                               // Parameter is ignored

    double          CalculateThermalMassLossRate() const                                                        { return BaseStar::CalculateThermalMassLossRate(); }                    // Set thermal mass gain rate to be effectively infinite, using dynamical timescale (in practice, will be Eddington limited), avoid division by zero

    void            CalculateTimescales(const double p_Mass, DBL_VECTOR &p_Timescales)                          { return TPAGB::CalculateTimescales(p_Mass, p_Timescales); }            // Use TPAGB
    void            CalculateTimescales()                                                                       { CalculateTimescales(m_Mass0, m_Timescales); }                         // Use class member variables

    double          CalculateZeta(ZETA_PRESCRIPTION p_ZetaPrescription)                                         { return 0.0; }                                                         // Should never be called...

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

    void            ResolveMassLoss() { }                                                                                                                                         // NO-OP

    STELLAR_TYPE    ResolveSkippedPhase()                                                                       { return BaseStar::ResolveSkippedPhase(); }                             // Default to BaseStar
    STELLAR_TYPE    ResolveSupernova()                                                                          { return BaseStar::ResolveSupernova(); }                                // Default to BaseStar

    void            SetPulsarParameters() const { }                                                                                                                                     // NO-OP

    bool            ShouldEvolveOnPhase() const                                                                 { return true; }                                                        // Default
    bool            ShouldSkipPhase() const                                                                     { return false; }                                                       // Don't skip WD phase

};

#endif // __Remnants_h__
