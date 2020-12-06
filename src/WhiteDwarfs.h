#ifndef __WhiteDwarfs_h__
#define __WhiteDwarfs_h__

#include "constants.h"
#include "typedefs.h"
#include "profiling.h"
#include "utils.h"

#include "Remnants.h"


class BaseStar;
class Remnants;

class WhiteDwarfs: virtual public BaseStar, public Remnants {

public:

    WhiteDwarfs(const BaseStar &p_BaseStar, const bool p_Initialise = true) : BaseStar(p_BaseStar), Remnants(p_BaseStar, false) {
        if (p_Initialise) Initialise();
    }

    WhiteDwarfs& operator = (const BaseStar &p_BaseStar) { static_cast<BaseStar&>(*this) = p_BaseStar; return *this; }


    // member functions
    static  double      CalculateLuminosityOnPhase_Static(const double p_Mass, 
                                                          const double p_Time, 
                                                          const double p_Metallicity, 
                                                          const double p_BaryonNumber);

    static  double      CalculateRadiusOnPhase_Static(const double p_Mass);


protected:


    // member functions - alphabetically
    double          CalculateCOCoreMassOnPhase() const                                                  { return m_COCoreMass; }                                                // NO-OP

    double          CalculateCoreMassOnPhase() const                                                    { return m_Mass; }                                                      // Return m_Mass

    double          CalculateEddingtonCriticalRate() const                                              { return 1.5E-8 * (m_Radius * RSOL_TO_KM / 10.0) * MYR_TO_YEAR; }       // Sluys 2013 ("Binary Evolution in a Nutshell"), eq 70

    void            CalculateGBParams()                                                                 { GiantBranch::CalculateGBParams(); }                                   // Default to GiantBranch

    double          CalculateGyrationRadius() const                                                     { return 0.21; }                                                        // Hurley et al., 2000, after eq 109 for n=3/2 polytrope or dense convective core. Single number approximation.

    double          CalculateHeCoreMassOnPhase() const                                                  { return m_HeCoreMass; }                                                // NO-OP

    double          CalculateInitialSupernovaMass() const                                               { return 0.0; }

    double          CalculateLambdaDewi()                                                               { return BaseStar::CalculateLambdaDewi(); }                             // Not supported - use BaseStar
    double          CalculateLambdaNanjing()                                                            { return BaseStar::CalculateLambdaNanjing(); }                          // Not supported - use BaseStar     JR: todo: check this (type 10 not mentioned as not supported in original code)

    DBL_DBL         CalculateMassAcceptanceRate(const double p_DonorMassRate,
                                                const double p_AccretorMassRate = 0.0 );

    double          CalculateMassLossRateHurley();
    double          CalculateMassLossRateVink() const                                                   { return 0.0; }

    double          CalculateMomentOfInertia(const double p_RemnantRadius = 0.0)                        { return GiantBranch::CalculateMomentOfInertia(p_RemnantRadius); }      // Default to GiantBranch
    double          CalculateMomentOfInertiaAU(const double p_RemnantRadius = 0.0)                      { return GiantBranch::CalculateMomentOfInertiaAU(p_RemnantRadius); }    // Default to GiantBranch

    double          CalculatePerturbationMuOnPhase() const                                              { return m_Mu; }                                                        // NO-OP

    double          CalculateRadialExtentConvectiveEnvelope()                                           { return BaseStar::CalculateRadialExtentConvectiveEnvelope(); }         // WD stars don't have a convective envelope



    double          CalculateRadiusOnPhase(const double p_Mass)                                         { return CalculateRadiusOnPhase_Static(p_Mass); }
    double          CalculateRadiusOnPhase()                                                            { return CalculateRadiusOnPhase(m_Mass); }                              // Use class member variables
    std::tuple <double, STELLAR_TYPE> CalculateRadiusAndStellarTypeOnPhase()                            { return BaseStar::CalculateRadiusAndStellarTypeOnPhase(); }




    double          CalculateTauOnPhase() const                                                         { return m_Tau; }                                                       // NO-OP
   
    double          CalculateThermalTimescale() const                                                   { return CalculateDynamicalTimescale(); }                               // Use dynamical timescale for mass transfer purposes
    double          CalculateThermalTimescale(const double p_Mass,
                                                  const double p_Radius,
                                                  const double p_Luminosity,
                                                  const double p_EnvMass = 1.0) const                   { return CalculateThermalTimescale(); }                                 // Ignore parameters

    double          CalculateThermalMassLossRate()                                                      { return BaseStar::CalculateThermalMassLossRate(); }                    // Set thermal mass gain rate to be effectively infinite, using dynamical timescale (in practice, will be Eddington limited), avoid division by zero

    void            CalculateTimescales(const double p_Mass, DBL_VECTOR &p_Timescales)                  { return TPAGB::CalculateTimescales(p_Mass, p_Timescales); }            // Use TPAGB
    void            CalculateTimescales()                                                               { CalculateTimescales(m_Mass0, m_Timescales); }                         // Use class member variables

    double          CalculateZeta(ZETA_PRESCRIPTION p_ZetaPrescription)                                 { m_Error = ERROR::INVALID_TYPE_ZETA_CALCULATION;                       // Set error value
                                                                                                          SHOW_WARN(m_Error);                                                   // Warn that an error occurred
                                                                                                          return 0.0; }                                                         // Should never be called...

    double          ChooseTimestep(const double p_Time) const;

    ENVELOPE        DetermineEnvelopeType() const                                                       { return ENVELOPE::CONVECTIVE; }                                        // Always CONVECTIVE

    MT_CASE         DetermineMassTransferCase() const                                                   { return MT_CASE::NONE; }                                               // No Mass Transfer Case for WDs/Remnants

    void            EvolveOneTimestepPreamble()                                                         { BaseStar::EvolveOneTimestepPreamble(); }                              // Default to BaseStar

    STELLAR_TYPE    EvolveToNextPhase()                                                                 { return BaseStar::EvolveToNextPhase(); }                               // Default to BaseStar

    bool            IsDegenerate() const                                                                { return true; }                                                        // White Dwarfs, NS and BH are degenerate

    bool            IsSupernova() const                                                                 { return false; }                                                       // Default

    bool            IsMassRatioUnstable(const double p_AccretorMass, const bool p_AccretorIsDegenerate) { m_Error = ERROR::INVALID_TYPE_MT_MASS_RATIO;                          // Set error value
                                                                                                          SHOW_WARN(m_Error);                                                   // Warn that an error occurred
                                                                                                          return false; }                                                       // Should never be called...

    void            PerturbLuminosityAndRadius() const { }                                                                                                                      // NO-OP



    STELLAR_TYPE    ResolveEnvelopeLoss(bool p_NoCheck = false)                                         { return BaseStar::ResolveEnvelopeLoss(p_NoCheck); }                    // Default to BaseStar



    void            ResolveEnvelopeMassAtPhaseEnd(const double p_Tau) const                             { ResolveEnvelopeMassOnPhase(p_Tau); }                                  // Same as on phase
    void            ResolveEnvelopeMassOnPhase(const double p_Tau) const { }                                                                                                    // NO-OP

    void            ResolveMassLoss() const { }                                                                                                                                 // NO-OP

    STELLAR_TYPE    ResolveRemnantAfterEnvelopeLoss()                                                   { return BaseStar::ResolveRemnantAfterEnvelopeLoss(); }                 // Default to BaseStar

    STELLAR_TYPE    ResolveSkippedPhase()                                                               { return BaseStar::ResolveSkippedPhase(); }                             // Default to BaseStar
    STELLAR_TYPE    ResolveSupernova()                                                                  { return BaseStar::ResolveSupernova(); }                                // Default to BaseStar

    void            SetPulsarParameters() const { }                                                                                                                             // NO-OP

    bool            ShouldEvolveOnPhase() const                                                         { return true; }                                                        // Default
    bool            ShouldSkipPhase() const                                                             { return false; }                                                       // Don't skip WD phase

};

#endif // __WhiteDwarfs_h__