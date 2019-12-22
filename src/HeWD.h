#ifndef __HeWD_h__
#define __HeWD_h__

#include "constants.h"
#include "typedefs.h"
#include "utils.h"

#include "HeGB.h"


class BaseStar;
class HeGB;

class HeWD: virtual public BaseStar, public HeGB {

public:

    HeWD(const BaseStar &p_BaseStar, const bool p_Initialise = true) : BaseStar(p_BaseStar), HeGB(p_BaseStar, false) {
        if (p_Initialise) Initialise();
    }

    HeWD& operator = (const BaseStar &p_BaseStar) {
        static_cast<BaseStar&>(*this) = p_BaseStar;
        Initialise();
        return *this;
    }


    // member functions
    static  double      CalculateLuminosityOnPhase_Static(const double p_Mass, const double p_Time, const double p_Metallicity);
    static  double      CalculateRadiusOnPhase_Static(const double p_Mass);


protected:

    void Initialise() {
        m_StellarType = STELLAR_TYPE::HELIUM_WHITE_DWARF;                                                                                                                               // Set stellar type
        CalculateTimescales();                                                                                                                                                          // Initialise timescales
        m_Age = 0.0;                                                                                                                                                                    // Set age appropriately
    }


    // member functions - alphabetically

    // Many of the following prototypes aren't really required - they should never get called because we stop evolution
    // before we get to Helium White Dwarfs - including them just in case someone decides not to stop...
    // Stellar types beyond HeWD (i.e. ONeWD, NS, BH, MR) will default to these by following the inheritance hierarchy
    double          CalculateCOCoreMassOnPhase()                                                        { return m_COCoreMass; }                                                        // NO-OP

    double          CalculateConvergedMassStepZetaThermal()                                             { return BaseStar::CalculateConvergedMassStepZetaThermal(); }                   // Use BaseStar

    double          CalculateCoreMassOnPhase()                                                          { return m_Mass; }                                                              // Return m_Mass

    double          CalculateEddingtonCriticalRate()                                                    { return 1.5E-8 * (m_Radius * RSOL_TO_KM / 10.0) * MYR_TO_YEAR; }               // Sluys 2013 ("Binary Evolution in a Nutshell"), eq 70

    double          CalculateEnvelopeMassOnPhase(const double p_Tau)                                    { return 0.0; }

    void            CalculateGBParams()                                                                 { GiantBranch::CalculateGBParams(); }                                           // Default to GiantBranch

    double          CalculateGyrationRadius()                                                           { return 0.21; }                                                                // Hurley et al., 2000, after eq 109 for n=3/2 polytrope or dense convective core. Single number approximation.


    double          CalculateHeCoreMassOnPhase()                                                        { return m_HeCoreMass; }                                                        // NO-OP

    double          CalculateLambdaDewi()                                                               { return BaseStar::CalculateLambdaDewi(); }                                     // Not supported - use BaseStar
    double          CalculateLambdaNanjing()                                                            { return BaseStar::CalculateLambdaNanjing(); }                                  // Not supported - use BaseStar     JR: todo: check this (type 10 not mentioned as not supported in original code)

    double          CalculateLuminosityOnPhase(const double p_Mass,
                                               const double p_Time,
                                               const double p_Metallicity)                              { return CalculateLuminosityOnPhase_Static(p_Mass, p_Time, p_Metallicity); }
    double          CalculateLuminosityOnPhase()                                                        { return CalculateLuminosityOnPhase_Static(m_Mass, m_Age, m_Metallicity); }     // Use class member variables

    DBL_DBL         CalculateMassAcceptanceRate(const double p_DonorMassRate,
                                                const double p_FractionAccreted,
                                                const double p_AccretorMassRate = 0.0 );

    double          CalculateMassLossRateHurley();
    double          CalculateMassLossRateVink()                                                         { return 0.0; }

    double          CalculatePerturbationMuOnPhase()                                                    { return m_Mu; }                                                                // NO-OP

    double          CalculateRadialExtentConvectiveEnvelope()                                           { return BaseStar::CalculateRadialExtentConvectiveEnvelope(); }                 // HeWD stars don't have a convective envelope

    double          CalculateRadiusOnPhase(const double p_Mass)                                         { return CalculateRadiusOnPhase_Static(p_Mass); }
    double          CalculateRadiusOnPhase()                                                            { return CalculateRadiusOnPhase(m_Mass); }                                      // Use class member variables
    std::tuple <double, STELLAR_TYPE> CalculateRadiusAndStellarTypeOnPhase()                            { return BaseStar::CalculateRadiusAndStellarTypeOnPhase(); }


    double          CalculateTauOnPhase()                                                               { return m_Tau; }                                                               // NO-OP

    double          CalculateThermalTimescale(const double p_Mass,
                                                  const double p_Radius,
                                                  const double p_Luminosity,
                                                  const double p_EnvMass = 1.0)                         { m_Error = ERROR::INVALID_TYPE_MT_THERMAL_TIMESCALE;                           // Set error value
                                                                                                          SHOW_WARN(m_Error);                                                           // Warn that an error occurred
                                                                                                          return 0.0; }                                                                 // Should never be called...

    double          CalculateThermalMassLossRate()                                                      { return BaseStar::CalculateThermalMassLossRate(); }                            // Use BaseStar

    double          CalculateThermalTimescale()                                                         { return CalculateThermalTimescale(m_Mass, m_Radius, m_Luminosity); }           // Should never be called...

    void            CalculateTimescales(const double p_Mass, DBL_VECTOR &p_Timescales);
    void            CalculateTimescales()                                                               { CalculateTimescales(m_Mass0, m_Timescales); }                                 // Use class member variables

    double          CalculateZeta(CE_ZETA_PRESCRIPTION p_CEZetaPrescription)                            { m_Error = ERROR::INVALID_TYPE_ZETA_CALCULATION;                               // Set error value
                                                                                                          SHOW_WARN(m_Error);                                                           // Warn that an error occurred
                                                                                                          return 0.0; }                                                                 // Should never be called...

    void            CheckRunaway(const bool p_Unbound) { }                                                                                                      // NO-OP

    double          ChooseTimestep(const double p_Time);

    ENVELOPE        DetermineEnvelopeType()                                                             { return ENVELOPE::CONVECTIVE; }                                                // Always CONVECTIVE
    ENVELOPE        DetermineEnvelopeTypeHurley2002()                                                   { return ENVELOPE::REMNANT; }                                                   // JR: todo: not convective according to Hurley et al. 2002, but not radiative - so set remnant.  is this right?

    MT_CASE         DetermineMassTransferCase()                                                         { return MT_CASE::NONE; }                                                       // No Mass Transfer Case for WDs/Remnants

    void            EvolveOneTimestepPreamble()                                                         { BaseStar::EvolveOneTimestepPreamble(); }                                      // Default to BaseStar

    STELLAR_TYPE    EvolveToNextPhase()                                                                 { return BaseStar::EvolveToNextPhase(); }                                       // Default to BaseStar

    bool            IsDegenerate() const                                                                { return true; }                                                                // White Dwarfs, NS and BH are degenerate

    bool            IsMassRatioUnstable(const double p_AccretorMass, const bool p_AccretorIsDegenerate) { m_Error = ERROR::INVALID_TYPE_MT_MASS_RATIO;                                  // Set error value
                                                                                                          SHOW_WARN(m_Error);                                                           // Warn that an error occurred
                                                                                                          return false; }                                                               // Should never be called...

    void            PerturbLuminosityAndRadius() { }                                                                                                                                    // NO-OP

    STELLAR_TYPE    ResolveEnvelopeLoss(bool p_NoCheck = false)                                         { return BaseStar::ResolveEnvelopeLoss(p_NoCheck); }                            // Default to BaseStar
    void            ResolveEnvelopeMassAtPhaseEnd(const double p_Tau)                                   { ResolveEnvelopeMassOnPhase(p_Tau); }                                          // Same as on phase
    void            ResolveEnvelopeMassOnPhase(const double p_Tau) { }                                                                                                                  // NO-OP
    void            ResolveMassLoss() { }                                                                                                                                               // NO-OP
    STELLAR_TYPE    ResolveRemnantAfterEnvelopeLoss()                                                   { return BaseStar::ResolveRemnantAfterEnvelopeLoss(); }                         // Default to BaseStar
    STELLAR_TYPE    ResolveSkippedPhase()                                                               { return BaseStar::ResolveSkippedPhase(); }                                     // Default to BaseStar
    STELLAR_TYPE    ResolveSupernova()                                                                  { return BaseStar::ResolveSupernova(); }                                        // Default to BaseStar

    void            SetPulsarParameters() { }                                                                                                                                           // NO-OP

    bool            ShouldEvolveOnPhase()                                                               { return true; }                                                                // Always it seems...  JR: todo: check this
    bool            ShouldSkipPhase()                                                                   { return false; }                                                               // Don't skip

};

#endif // __HeWD_h__
