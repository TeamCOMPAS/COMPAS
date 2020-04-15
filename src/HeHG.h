#ifndef __HeHG_h__
#define __HeHG_h__

#include "constants.h"
#include "typedefs.h"
#include "utils.h"

#include "HeMS.h"


class BaseStar;
class HeMS;

class HeHG: virtual public BaseStar, public HeMS {

public:

    HeHG(const BaseStar &p_BaseStar, const bool p_Initialise = true) : BaseStar(p_BaseStar), HeMS(p_BaseStar, false) {
        if (p_Initialise) Initialise();
    }

    HeHG& operator = (const BaseStar &p_BaseStar) {
        static_cast<BaseStar&>(*this) = p_BaseStar;
        Initialise();
        return *this;
    }

    static void CalculateGBParams_Static(const double p_Mass0, const double p_Mass, const double p_LogMetallicityXi, const DBL_VECTOR &p_MassCutoffs, const DBL_VECTOR &p_AnCoefficients, const DBL_VECTOR &p_BnCoefficients, DBL_VECTOR &p_GBParams);


protected:

    void Initialise() {
        m_StellarType = STELLAR_TYPE::NAKED_HELIUM_STAR_HERTZSPRUNG_GAP;                                                                                                            // Set stellar type
        CalculateTimescales();                                                                                                                                                      // Initialise timescales
        // JR: Age for HeHG is calculated before switching -
        // can get here via EvolveOneTimestep() and ResolveEnvelopeLoss(),
        // and Age is calculated differently in those cases
    }


    // member functions - aphabetically
            double          CalculateCOCoreMassAtPhaseEnd()                                                  { return m_COCoreMass; }                                               // NO-OP
            double          CalculateCOCoreMassOnPhase();

            double          CalculateConvergedMassStepZetaThermal()                                          { return FGB::CalculateConvergedMassStepZetaThermal(); }               // Skip HeMS

            double          CalculateCoreMassAtBAGB()                                                        { return m_Mass0; }                                                    // McBAGB = M0 (Hurely et al. 2000, discussion just before eq 89)
            double          CalculateCoreMassAtPhaseEnd()                                                    { return m_CoreMass; }                                                 // NO-OP
            double          CalculateCoreMassOnPhase()                                                       { return m_COCoreMass; }                                               // Mc(HeMS) = McCOMass

    static  double          CalculateCoreMass_Luminosity_B_Static()                                          { return 4.1E4; }
    static  double          CalculateCoreMass_Luminosity_D_Static(const double p_Mass)                       { return 5.5E4 / (1.0 + (0.4 * p_Mass * p_Mass * p_Mass * p_Mass)); }  // pow() is slow - use multiplication

            double          CalculateEnvelopeMassOnPhase(const double p_Tau)                                 { return BaseStar::CalculateEnvelopeMassOnPhase(p_Tau); }

            void            CalculateGBParams(const double p_Mass, DBL_VECTOR &p_GBParams);
            void            CalculateGBParams()                                                              { CalculateGBParams(m_Mass0, m_GBParams); }                            // Use class member variables

            double          CalculateGyrationRadius()                                                        { return 0.21; }                                                       // Hurley et al., 2000, after eq 109 for n=3/2 polytrope or dense convective core. Single number approximation.

            double          CalculateHeCoreMassAtPhaseEnd()                                                  { return CalculateHeCoreMassOnPhase(); }                               // Same as on phase
            double          CalculateHeCoreMassOnPhase()                                                     { return m_Mass; }                                               // NO-OP

            double          CalculateLambdaNanjing();

            double          CalculateLuminosityOnPhase();
            double          CalculateLuminosityAtPhaseEnd()                                                  { return m_Luminosity; }                                                                                                // NO-OP

            double          CalculateMassTransferRejuvenationFactor();

   	        double          CalculateMomentOfInertia(const double p_RemnantRadius = 0.0)                     { return GiantBranch::CalculateMomentOfInertia(p_RemnantRadius); }     // Skip HeMS
   	        double          CalculateMomentOfInertiaAU(const double p_RemnantRadius = 0.0)                   { return GiantBranch::CalculateMomentOfInertiaAU(p_RemnantRadius); }   // Skip HeMS

            double          CalculatePerturbationMu();
            double          CalculatePerturbationMuAtPhaseEnd()                                              { return m_Mu; }                                                       // NO-OP

            double          CalculateRadialExtentConvectiveEnvelope()                                        { return GiantBranch::CalculateRadialExtentConvectiveEnvelope(); }     // Skip HeMS

            double          CalculateRadiusAtPhaseEnd()                                                      { return m_Radius; }                                                                                                        // NO-OP
            double          CalculateRadiusOnPhase();

            std::tuple <double, STELLAR_TYPE> CalculateRadiusAndStellarTypeOnPhase(const double p_Mass, const double p_Luminosity);
            std::tuple <double, STELLAR_TYPE> CalculateRadiusAndStellarTypeOnPhase()                         { return CalculateRadiusAndStellarTypeOnPhase(m_Mass, m_Luminosity); }

            double          CalculateTauAtPhaseEnd()                                                         { return m_Tau; }                                                      // NO-OP
            double          CalculateTauOnPhase()                                                            { return 0.0; }

            double          CalculateTemperatureAtPhaseEnd(const double p_Luminosity, const double p_Radius) { return m_Temperature; }                                              // NO-OP
            double          CalculateTemperatureAtPhaseEnd()                                                 { return CalculateTemperatureAtPhaseEnd(m_Luminosity, m_Radius); }     // Use class member variables

            double          CalculateThermalMassLossRate()                                                   { return GiantBranch::CalculateThermalMassLossRate(); }                // Skip HeMS

            void            CalculateTimescales(const double p_Mass, DBL_VECTOR &p_Timescales);
            void            CalculateTimescales()                                                            { CalculateTimescales(m_Mass0, m_Timescales); }                        // Use class member variables

            double          ChooseTimestep(const double p_Time);

            ENVELOPE        DetermineEnvelopeType()                                                          { return ENVELOPE::CONVECTIVE; }                                       // Always CONVECTIVE
            ENVELOPE        DetermineEnvelopeTypeHurley2002()                                                { return ENVELOPE::CONVECTIVE; }                                       // Always CONVECTIVE

            MT_CASE         DetermineMassTransferCase()                                                      { return MT_CASE::B; }                                                 // Mass Transfer Case B for HeHG stars

            STELLAR_TYPE    EvolveToNextPhase();

            bool            IsEndOfPhase()                                                                   { return !ShouldEvolveOnPhase(); }
            bool            IsMassRatioUnstable(const double p_AccretorMass, const bool p_AccretorIsDegenerate);
            bool            IsSupernova()                                                                    { return (utils::Compare(m_COCoreMass, CalculateMaximumCoreMassSN()) > 0); }   // Going supernova if CO core mass large enough

            void            PerturbLuminosityAndRadius()                                                     { GiantBranch::PerturbLuminosityAndRadius(); }                                                                                        // NO-OP

            STELLAR_TYPE    ResolveEnvelopeLoss(bool p_NoCheck = false);
            void            ResolveHeliumFlash() { }                                                                                                                                // NO-OP
            STELLAR_TYPE    ResolveSkippedPhase()                                                            { return m_StellarType; }                                              // NO-OP

            bool            ShouldEvolveOnPhase();
            bool            ShouldSkipPhase()                                                                { return false; }                                                      // Never skip HeMS phase

            void            UpdateAgeAfterMassLoss()                                                         { GiantBranch::UpdateAgeAfterMassLoss(); }                             // Skip HeMS
            void            UpdateInitialMass()                                                              { GiantBranch::UpdateInitialMass(); }                                  // Skip HeMS

};

#endif // __HeHG_h__
