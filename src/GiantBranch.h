#ifndef __GiantBranch_h__
#define __GiantBranch_h__

#include "constants.h"
#include "typedefs.h"
#include "profiling.h"
#include "utils.h"
#include "Rand.h"

#include "MainSequence.h"


class BaseStar;
class MainSequence;

class GiantBranch: virtual public BaseStar, public MainSequence {

public:

    GiantBranch(const BaseStar &p_BaseStar) : BaseStar(p_BaseStar), MainSequence(p_BaseStar) {}


protected:


    // member functions - alphabetically (sort of - some are grouped by functionality)
            double          CalculateConvectiveCoreMass() const { return m_CoreMass; }
            DBL_DBL         CalculateConvectiveEnvelopeMass() const;
    static  double          CalculateCoreMassAt2ndDredgeUp_Static(const double p_McBAGB);
            double          CalculateCoreMassAtBAGB(const double p_Mass) const;
    static  double          CalculateCoreMassAtBAGB_Static(const double p_Mass, const DBL_VECTOR &p_BnCoefficients);
            double          CalculateCoreMassAtBGB(const double p_Mass, const DBL_VECTOR &p_GBParams);
    static  double          CalculateCoreMassAtBGB_Static(const double p_Mass, const DBL_VECTOR &p_MassCutoffs, const DBL_VECTOR &p_AnCoefficients, const DBL_VECTOR &p_GBParams);
            double          CalculateCoreMassAtHeIgnition(const double p_Mass) const;
    static  double          CalculateCoreMassAtSupernova_Static(const double p_McBAGB);

    static  double          CalculateCoreMass_Luminosity_B_Static(const double p_Mass);
    static  double          CalculateCoreMass_Luminosity_D_Static(const double p_Mass, const double p_LogMetallicityXi, const DBL_VECTOR &p_MassCutoffs);
    static  double          CalculateCoreMass_Luminosity_p_Static(const double p_Mass, const DBL_VECTOR &p_MassCutoffs);
    static  double          CalculateCoreMass_Luminosity_q_Static(const double p_Mass, const DBL_VECTOR &p_MassCutoffs);
    static  double          CalculateCoreMass_Luminosity_Lx_Static(const DBL_VECTOR &p_GBParams);
    static  double          CalculateCoreMass_Luminosity_Mx_Static(const DBL_VECTOR &p_GBParams);

            double          CalculateCriticalMassRatioClaeys14(const bool p_AccretorIsDegenerate) const; 
            double          CalculateCriticalMassRatioHurleyHjellmingWebbink() const; 

            void            CalculateGBParams(const double p_Mass, DBL_VECTOR &p_GBParams);
    static  void            CalculateGBParams_Static(const double p_Mass, const double p_LogMetallicityXi, const DBL_VECTOR &p_MassCutoffs, const DBL_VECTOR &p_AnCoefficients, const DBL_VECTOR &p_BnCoefficients, DBL_VECTOR &p_GBParams);
            void            CalculateGBParams()                                                             { CalculateGBParams(m_Mass0, m_GBParams); }                         // Use class member variables

    static  double          CalculateHRateConstant_Static(const double p_Mass);

    virtual double          CalculateInitialSupernovaMass() const                                           { return m_Mass; }                                                  // Use class member variables

            double          CalculateLifetimeToHeIgnition(const double p_Mass, const double p_Tinf1_FGB, const double p_Tinf2_FGB);

    static  double          CalculateLuminosityAtHeIgnition_Static(const double      p_Mass,
                                                                   const double      p_Alpha,
                                                                   const double      p_MHeF,
                                                                   const DBL_VECTOR &p_BnCoefficients);

    static  double          CalculateLuminosityAtPhaseBase_Static(const double p_Mass, const DBL_VECTOR &p_AnCoefficients);
    static  double          CalculateLuminosityOnZAHB_Static(const double      p_Mass,
                                                             const double      p_CoreMass,
                                                             const double      p_Alpha1,
                                                             const double      p_MHeF,
                                                             const double      p_MFGB,
                                                             const double      p_MinimumLuminosityOnPhase,
                                                             const DBL_VECTOR &p_BnCoefficients);

            double          CalculateMassLossRateHurley();

            double          CalculateBaryonicRemnantMass(const double p_ProtoMass, double p_FallbackMass);
            double          CalculateFallbackByBelczynski2002(const double p_COCoreMass);
            double          CalculateFallbackFractionDelayed(const double p_PreSNMass, const double p_ProtoMass, const double p_COCoreMass);
            double          CalculateFallbackFractionRapid(const double p_PreSNMass, const double p_ProtoMass, const double p_COCoreMass);
            double          CalculateFallbackMass(const double p_PreSNMass, const double p_ProtoMass, const double p_Fallback);
            double          CalculateGravitationalRemnantMass(const double p_BaryonicRemnantMass);
            double          CalculateProtoCoreMassDelayed(const double p_COCoreMass);
            double          CalculateProtoCoreMassRapid();
            double          CalculateRemnantMassByBelczynski2002(const double p_Mass, const double p_COCoreMass, const double p_FallbackFraction);
            DBL_DBL         CalculateRemnantMassByFryer2012(const double p_Mass, const double p_COCoreMass);
            DBL_DBL         CalculateRemnantMassByFryer2022(const double p_Mass, const double p_COCoreMass);
            double          CalculateRemnantMassByMuller2016(const double p_Mass, const double p_COCoreMass);
            double          CalculateRemnantMassByMullerMandel(const double p_COCoreMass, const double p_HeCoreMass);
            double          CalculateRemnantMassBySchneider2020(const double p_COCoreMass, const bool p_UseSchneiderAlt = false);
            double          CalculateRemnantMassBySchneider2020Alt(const double p_COCoreMass)               { return CalculateRemnantMassBySchneider2020(p_COCoreMass, true); }

            double          CalculateMomentOfInertia() const;

            double          CalculatePerturbationMu() const;

            double          CalculateRadialExtentConvectiveEnvelope() const                                 { return (m_Radius - CalculateConvectiveCoreRadius()); }    //Hurley et al. 2002, sec. 2.3, particularly subsec. 2.3.1, eqs 36-40

            double          CalculateRadiusAtHeIgnition(const double p_Mass) const;
            double          CalculateRadiusOnPhase(const double p_Mass, const double p_Luminosity) const    { return CalculateRadiusOnPhase_Static(p_Mass, p_Luminosity, m_BnCoefficients); }
            double          CalculateRadiusOnPhase() const                                                  { return CalculateRadiusOnPhase(m_Mass, m_Luminosity); }
    static  double          CalculateRadiusOnPhase_Static(const double p_Mass, const double p_Luminosity, const DBL_VECTOR &p_BnCoefficients);
    static  double          CalculateRadiusOnZAHB_Static(const double      p_Mass,
                                                         const double      p_CoreMass,
                                                         const double      p_Alpha1,
                                                         const double      p_MHeF,
                                                         const double      p_MFGB,
                                                         const double      p_MinimumLuminosityOnPhase,
                                                         const DBL_VECTOR &p_BnCoefficients);

    virtual double          CalculateRemnantLuminosity() const;
            STELLAR_TYPE    CalculateRemnantTypeByMuller2016(const double p_COCoreMass);
	
    virtual double          CalculateRemnantRadius() const;

            double          CalculateThermalMassLossRate() const                                            { return (m_Mass - m_CoreMass) / CalculateThermalTimescale(); }     // Use class member variables

            void            CalculateTimescales(const double p_Mass, DBL_VECTOR &p_Timescales);
            void            CalculateTimescales()                                                           { CalculateTimescales(m_Mass0, m_Timescales); }                     // Use class member variables

            double          CalculateZetaConstantsByEnvelope(ZETA_PRESCRIPTION p_ZetaPrescription);
            double          CalculateZetaConvectiveEnvelopeGiant(ZETA_PRESCRIPTION p_ZetaPrescription);

    virtual void            PerturbLuminosityAndRadius();

            STELLAR_TYPE    ResolveSupernova();
            STELLAR_TYPE    ResolveCoreCollapseSN();
            STELLAR_TYPE    ResolveElectronCaptureSN();
            STELLAR_TYPE    ResolvePairInstabilitySN();
            STELLAR_TYPE    ResolvePulsationalPairInstabilitySN();
            STELLAR_TYPE    ResolveTypeIIaSN();
    
            void            UpdateAgeAfterMassLoss() { }                                                                                                                        // NO-OP for most stellar types

            void            UpdateInitialMass() { }                                                                                                                             // NO-OP for most stellar types
    
            void            UpdateMinimumCoreMass() { }                                                                                                                         // NO-OP for most stellar types

};

#endif // __GiantBranch_h__
