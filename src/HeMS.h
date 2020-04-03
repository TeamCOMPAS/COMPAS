#ifndef __HeMS_h__
#define __HeMS_h__

#include "constants.h"
#include "typedefs.h"
#include "utils.h"

#include "TPAGB.h"


// JR: todo: revisit this one day - sometimes HeMS is MS, sometimes it is GiantBranch...
// Right now it is GiantBranch - figure out which is more appropriate

class BaseStar;
class TPAGB;

class HeMS: virtual public BaseStar, public TPAGB {

public:

    HeMS(const BaseStar &p_BaseStar, const bool p_Initialise = true) : BaseStar(p_BaseStar), TPAGB(p_BaseStar, false) {
        if (p_Initialise) Initialise();
    }

    HeMS& operator = (const BaseStar &p_BaseStar) {
        static_cast<BaseStar&>(*this) = p_BaseStar;
        Initialise();
        return *this;
    }


    // member functions - alphabetically
    static double   CalculateLifetimeOnPhase_Static(const double p_Mass);

    static double   CalculateLuminosityAtZAMS_Static(const double p_Mass);
    static double   CalculateLuminosityOnPhase_Static(const double p_Mass, const double p_Tau);
    static double   CalculateLuminosityAtPhaseEnd_Static(const double p_Mass);


    static DBL_DBL  CalculateRadiusAtPhaseEnd_Static(const double p_Mass, const double p_Luminosity);
    static double   CalculateRadiusAtZAMS_Static(const double p_Mass);
    static double   CalculateRadiusOnPhase_Static(const double p_Mass, const double p_Tau);


protected:

    void Initialise() {
        m_StellarType = STELLAR_TYPE::NAKED_HELIUM_STAR_MS;                                                                                                                 // Set stellar type
        CalculateTimescales();
        // JR: Age for HeMS is partially calculated before switching -
        // can get here from various places in ResolveEnvelopeLoss(),
        // and Age is calculated differently in those cases
    }


    // member functions - alphabetically
            double          CalculateCOCoreMassAtPhaseEnd()                                      { return CalculateCOCoreMassOnPhase(); }                                // Same as on phase
            double          CalculateCOCoreMassOnPhase()                                         { return 0.0; }                                                         // McCO(HeMS) = 0.0

            double          CalculateConvergedMassStepZetaThermal()                              { return BaseStar::CalculateConvergedMassStepZetaThermal(); }           // Use BaseStar

            double          CalculateCoreMassAtPhaseEnd()                                        { return CalculateHeCoreMassOnPhase(); }                                // Same as on phase
            double          CalculateCoreMassOnPhase()                                           { return m_COCoreMass; }                                                // Mc(HeMS) = McCOMass

            double          CalculateGyrationRadius()                                            { return 0.1; }                                                         // JR: todo: Nobody seems sure about this...

            double          CalculateHeCoreMassOnPhase()                                         { return m_Mass; }                                                      // McHe(HeMS) = Mass
            double          CalculateHeCoreMassAtPhaseEnd()                                      { return CalculateHeCoreMassOnPhase(); }                                // Same as on phase

            double          CalculateInitialSupernovaMass()                                      { return GiantBranch::CalculateInitialSupernovaMass(); }                // Use GiantBranch

            double          CalculateLambdaDewi()                                                { return 0.5; }
            double          CalculateLambdaNanjing()                                             { return BaseStar::CalculateLambdaNanjing(); }                          // Not supported - use BaseStar

            double          CalculateLuminosityAtPhaseEnd(const double p_Mass)                   { return CalculateLuminosityAtPhaseEnd_Static(p_Mass); }
            double          CalculateLuminosityAtPhaseEnd()                                      { return CalculateLuminosityAtPhaseEnd(m_Mass); }                       // Use class member variables
            double          CalculateLuminosityOnPhase(const double p_Mass, const double p_Tau)  { return CalculateLuminosityOnPhase_Static(p_Mass, p_Tau); }
            double          CalculateLuminosityOnPhase()                                         { return CalculateLuminosityOnPhase(m_Mass, m_Tau); }                   // Use class member variables

            double          CalculateMassLossRateHurley();
            double          CalculateMassLossRateVink();

            double          CalculateMassTransferRejuvenationFactor();

            double          CalculateMomentOfInertia(const double p_RemnantRadius = 0.0)         { return MainSequence::CalculateMomentOfInertia(p_RemnantRadius); }     // Use MainSequence
            double          CalculateMomentOfInertiaAU(const double p_RemnantRadius = 0.0)       { return MainSequence::CalculateMomentOfInertiaAU(p_RemnantRadius); }   // Use MainSequence

            double          CalculatePerturbationMu()                                            { return 5.0; }                                                         // Hurley et al. 2000, eqs 97 & 98

            double          CalculateRadialExtentConvectiveEnvelope()                            { return BaseStar::CalculateRadialExtentConvectiveEnvelope(); }         // HeMS stars don't have a convective envelope

            double          CalculateRadiusAtPhaseEnd(const double p_Mass)                       { return CalculateRadiusAtPhaseEnd_Static(p_Mass); }
            double          CalculateRadiusAtPhaseEnd()                                          { return CalculateRadiusAtPhaseEnd(m_Mass); }                           // Use class member variables
    static  double          CalculateRadiusAtPhaseEnd_Static(const double p_Mass);
            double          CalculateRadiusOnPhase(const double p_Mass, const double p_Tau)      { return CalculateRadiusOnPhase_Static(p_Mass, p_Tau); }
            double          CalculateRadiusOnPhase()                                             { return CalculateRadiusOnPhase(m_Mass, m_Tau); }                       // Use class member variables

            double          CalculateTauAtPhaseEnd()                                             { return 1.0; }
            double          CalculateTauOnPhase()                                                { return m_Age / m_Timescales[static_cast<int>(TIMESCALE::tHeMS)]; }

            double          CalculateTemperatureAtPhaseEnd()                                     { return BaseStar::CalculateTemperatureAtPhaseEnd(); }

            double          CalculateThermalMassLossRate()                                       { return BaseStar::CalculateThermalMassLossRate(); }                   // Use BaseStar

            double          CalculateThermalTimescale(const double p_Mass,
                                                      const double p_Radius,
                                                      const double p_Luminosity,
                                                      const double p_EnvMass = 1.0)              { return MainSequence::CalculateThermalTimescale(p_Mass, p_Radius, p_Luminosity); }
            double          CalculateThermalTimescale()                                          { return CalculateThermalTimescale(m_Mass, m_Radius, m_Luminosity); }  // Use class member variables

            void            CalculateTimescales(const double p_Mass, DBL_VECTOR &p_Timescales);
            void            CalculateTimescales()                                                { CalculateTimescales(m_Mass0, m_Timescales); }                         // Use class member variables

            double          ChooseTimestep(const double p_Time);

            MT_CASE         DetermineMassTransferCase()                                          { return MT_CASE::A; }                                                  // Mass Transfer Case A for HeMS stars

            STELLAR_TYPE    EvolveToNextPhase();

            ENVELOPE        DetermineEnvelopeType()                                              { return ENVELOPE::RADIATIVE; }                                         // Always RADIATIVE
            ENVELOPE        DetermineEnvelopeTypeHurley2002()                                    { return ENVELOPE::RADIATIVE; }                                         // Always RADIATIVE

            bool            IsEndOfPhase()                                                       { return !ShouldEvolveOnPhase(); }
            bool            IsMassRatioUnstable(const double p_AccretorMass, const bool p_AccretorIsDegenerate);
            bool            IsSupernova()                                                        { return false; }                                                       // Not here

            void            PerturbLuminosityAndRadius() { }                                                                                                             // NO-OP

            STELLAR_TYPE    ResolveEnvelopeLoss(bool p_NoCheck = false);
            void            ResolveHeliumFlash() { }                                                                                                                     // NO-OP
            STELLAR_TYPE    ResolveSkippedPhase()                                                { return m_StellarType; }                                               // NO-OP

            void            SetSNHydrogenContent()                                               { m_SupernovaDetails.hydrogenContent = HYDROGEN_CONTENT::POOR; }        // Always POOR

            bool            ShouldEvolveOnPhase()                                                { return (utils::Compare(m_Tau, 0.0) >= 0 && utils::Compare(m_Tau, 1.0) < 0); }       // Evolve on HeMS phase if 0 <= tau < 1.0
            bool            ShouldSkipPhase()                                                    { return false; }                                                       // Never skip HeMS phase

            void            UpdateInitialMass()                                                  { m_Mass0 = m_Mass; }                                                   // Per Hurley et al. 2000, section 7.1
            void            UpdateAgeAfterMassLoss();                                                                                                                    // Per Hurley et al. 2000, section 7.1

};

#endif // __HeMS_h__
