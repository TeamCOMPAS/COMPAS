#ifndef __HeHG_h__
#define __HeHG_h__

#include "constants.h"
#include "typedefs.h"
#include "profiling.h"
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
        m_StellarType = STELLAR_TYPE::NAKED_HELIUM_STAR_HERTZSPRUNG_GAP;                                                                                                                    // Set stellar type
        m_Tau = 0.0;                                                                                      // Start of phase
        CalculateTimescales();                                                                                                                                                              // Initialise timescales
        // JR: Age for HeHG is calculated before switching -
        // can get here via EvolveOneTimestep() and ResolveEnvelopeLoss(),
        // and Age is calculated differently in those cases
        
        //Update stellar properties at start of HeHG phase (since core defintion changes)
        CalculateGBParams();
        m_COCoreMass  = CalculateCOCoreMassOnPhase();
        m_CoreMass    = CalculateCoreMassOnPhase();
        m_HeCoreMass  = CalculateHeCoreMassOnPhase();
        m_Luminosity  = CalculateLuminosityOnPhase();
        std::tie(m_Radius, std::ignore) = CalculateRadiusAndStellarTypeOnPhase();   // Update radius
    }


    // member functions - aphabetically
            double          CalculateCOCoreMassAtPhaseEnd() const                                                   { return m_COCoreMass; }                                                // NO-OP
            double          CalculateCOCoreMassOnPhase() const;

            double          CalculateCoreMassAtBAGB() const                                                         { return m_Mass0; }                                                     // McBAGB = M0 (Hurely et al. 2000, discussion just before eq 89)
            double          CalculateCoreMassAtPhaseEnd() const                                                     { return m_CoreMass; }                                                  // NO-OP
            double          CalculateCoreMassOnPhase() const                                                        { return m_COCoreMass; }                                                // Mc(HeMS) = McCOMass

    static  double          CalculateCoreMass_Luminosity_B_Static()                                                 { return 4.1E4; }
    static  double          CalculateCoreMass_Luminosity_D_Static(const double p_Mass)                              { return 5.5E4 / (1.0 + (0.4 * p_Mass * p_Mass * p_Mass * p_Mass)); }   // pow() is slow - use multiplication

            void            CalculateGBParams(const double p_Mass, DBL_VECTOR &p_GBParams);
            void            CalculateGBParams()                                                                     { CalculateGBParams(m_Mass0, m_GBParams); }                             // Use class member variables

            double          CalculateGyrationRadius() const                                                         { return 0.21; }                                                        // Hurley et al., 2000, after eq 109 for n=3/2 polytrope or dense convective core. Single number approximation.

            double          CalculateHeCoreMassAtPhaseEnd() const                                                   { return CalculateHeCoreMassOnPhase(); }                                // Same as on phase
            double          CalculateHeCoreMassOnPhase() const                                                      { return m_Mass; }                                                      // NO-OP

            double          CalculateLambdaNanjingStarTrack(const double p_Mass, const double p_Metallicity) const;
            double          CalculateLambdaNanjingEnhanced(const int p_MassInd, const int p_Zind) const             { return CalculateLambdaNanjingStarTrack(0.0, 0.0); }                            // 0.0 are dummy values that are not used

            double          CalculateLuminosityOnPhase() const;
            double          CalculateLuminosityAtPhaseEnd() const                                                   { return m_Luminosity; }                                                // NO-OP

            double          CalculateMassTransferRejuvenationFactor() const;

   	        double          CalculateMomentOfInertia(const double p_RemnantRadius = 0.0) const                      { return GiantBranch::CalculateMomentOfInertia(p_RemnantRadius); }      // Skip HeMS
   	        double          CalculateMomentOfInertiaAU(const double p_RemnantRadius = 0.0) const                    { return GiantBranch::CalculateMomentOfInertiaAU(p_RemnantRadius); }    // Skip HeMS

            double          CalculatePerturbationMu() const;
            double          CalculatePerturbationMuAtPhaseEnd() const                                               { return m_Mu; }                                                        // NO-OP

            double          CalculateRadialExtentConvectiveEnvelope() const                                         { return GiantBranch::CalculateRadialExtentConvectiveEnvelope(); }      // Skip HeMS

            double          CalculateRadiusAtPhaseEnd() const                                                       { return m_Radius; }                                                    // NO-OP
            double          CalculateRadiusOnPhase() const;

            std::tuple <double, STELLAR_TYPE> CalculateRadiusAndStellarTypeOnPhase(const double p_Mass, const double p_Luminosity) const;
            std::tuple <double, STELLAR_TYPE> CalculateRadiusAndStellarTypeOnPhase() const                          { return CalculateRadiusAndStellarTypeOnPhase(m_Mass, m_Luminosity); }

            double          CalculateTauAtPhaseEnd() const                                                          { return m_Tau; }                                                       // NO-OP
            double          CalculateTauOnPhase() const                                                             { return 0.0; }

            double          CalculateTemperatureAtPhaseEnd(const double p_Luminosity, const double p_Radius) const  { return m_Temperature; }                                               // NO-OP
            double          CalculateTemperatureAtPhaseEnd() const                                                  { return CalculateTemperatureAtPhaseEnd(m_Luminosity, m_Radius); }      // Use class member variables

            double          CalculateThermalMassLossRate() const                                                    { return GiantBranch::CalculateThermalMassLossRate(); }                 // Skip HeMS

            void            CalculateTimescales(const double p_Mass, DBL_VECTOR &p_Timescales);
            void            CalculateTimescales()                                                                   { CalculateTimescales(m_Mass0, m_Timescales); }                         // Use class member variables
    
            double          CalculateZeta(ZETA_PRESCRIPTION p_ZetaPrescription)                                     { return HG::CalculateZeta(p_ZetaPrescription); }                       // Calculate Zetas as for HG and other giant stars (HeMS stars were an exception)

            double          ChooseTimestep(const double p_Time) const;

            ENVELOPE        DetermineEnvelopeType() const;

            STELLAR_TYPE    EvolveToNextPhase();

            bool            IsEndOfPhase() const                                                                    { return !ShouldEvolveOnPhase(); }
            bool            IsMassRatioUnstable(const double p_AccretorMass, const bool p_AccretorIsDegenerate) const;
            bool            IsSupernova() const;
            double          CalculateInitialSupernovaMass() const;
    
            void            PerturbLuminosityAndRadius()                                                            { GiantBranch::PerturbLuminosityAndRadius(); }                          // NO-OP

            STELLAR_TYPE    ResolveEnvelopeLoss(bool p_NoCheck = false);
            void            ResolveHeliumFlash() { }                                                                                                                                        // NO-OP
            STELLAR_TYPE    ResolveSkippedPhase()                                                                   { return m_StellarType; }                                               // NO-OP

            bool            ShouldEvolveOnPhase() const;
            bool            ShouldSkipPhase() const                                                                 { return false; }                                                       // Never skip HeMS phase

            void            UpdateAgeAfterMassLoss()                                                                { GiantBranch::UpdateAgeAfterMassLoss(); }                              // Skip HeMS
            void            UpdateInitialMass()                                                                     { GiantBranch::UpdateInitialMass(); }                                   // Skip HeMS

};

#endif // __HeHG_h__
