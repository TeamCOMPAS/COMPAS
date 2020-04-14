#ifndef __Star_h__
#define __Star_h__

#include <fstream>

#include "constants.h"
#include "typedefs.h"
#include "utils.h"

#include "BaseStar.h"
#include "MS_lte_07.h"
#include "MS_gt_07.h"
#include "CH.h"
#include "HG.h"
#include "FGB.h"
#include "CHeB.h"
#include "EAGB.h"
#include "TPAGB.h"
#include "HeMS.h"
#include "HeHG.h"
#include "HeGB.h"
#include "HeWD.h"
#include "COWD.h"
#include "ONeWD.h"
#include "NS.h"
#include "BH.h"
#include "MR.h"

class BaseStar;
class MS_lte_07;
class MS_gt_07;
class CH;
class HG;
class FGB;
class CHeB;
class EAGB;
class TPAGB;
class HeMS;
class HeHG;
class HeGB;
class HeWD;
class COWD;
class ONeWD;
class NS;
class BH;
class MR;


class Star {

public:

    Star();
    Star(const unsigned long int p_RandomSeed, 
         const double            p_MZAMS, 
         const double            p_Metallicity, 
         const KickParameters    p_KickParameters = {},
         const double            p_LBVfactor = 0.0, 
         const double            p_WolfRayetFactor = 0.0);
    Star(const Star& p_Star);
    Star& operator = (const Star& p_Star);

    virtual ~Star() { delete m_Star; delete m_SaveStar; }


    // object identifiers - all classes have these
    OBJECT_ID                   ObjectId() const                                                                            { return m_ObjectId; }
    OBJECT_TYPE                 ObjectType() const                                                                          { return m_ObjectType; }
    STELLAR_TYPE                InitialStellarType() const                                                                  { return m_Star->InitialStellarType(); }
    STELLAR_TYPE                StellarType() const                                                                         { return m_Star->StellarType(); }


    // getters - alphabetically
    double                      Age() const                                                                                 { return m_Star->Age(); }
    double                      BindingEnergy_Fixed() const                                                                 { return m_Star->BindingEnergy_Fixed(); }
    double                      BindingEnergy_Loveridge() const                                                             { return m_Star->BindingEnergy_Loveridge(); }
    double                      BindingEnergy_Nanjing() const                                                               { return m_Star->BindingEnergy_Nanjing(); }
    double                      BindingEnergy_Kruckow() const                                                               { return m_Star->BindingEnergy_Kruckow(); }
    double                      COCoreMass() const                                                                          { return m_Star->COCoreMass(); }
    double                      CoreMass() const                                                                            { return m_Star->CoreMass(); }
    double                      DynamicalTimescale() const                                                                  { return m_Star->DynamicalTimescale(); }
    double                      EnvMass() const                                                                             { return m_Star->EnvMass(); }
    bool                        ExperiencedCCSN() const                                                                     { return m_Star->ExperiencedCCSN(); }
    bool                        ExperiencedECSN() const                                                                     { return m_Star->ExperiencedECSN(); }
    bool                        ExperiencedPISN() const                                                                     { return m_Star->ExperiencedPISN() ; }
    bool                        ExperiencedPPISN() const                                                                    { return m_Star->ExperiencedPPISN(); }
    bool                        ExperiencedRunaway() const                                                                  { return m_Star->ExperiencedRunaway(); }
    bool                        ExperiencedUSSN() const                                                                     { return m_Star->ExperiencedUSSN(); }
    double                      HeCoreMass() const                                                                          { return m_Star->HeCoreMass(); }
    bool                        IsCCSN() const                                                                              { return m_Star->IsCCSN(); }
    bool                        IsDegenerate() const                                                                        { return m_Star->IsDegenerate(); }
    bool                        IsECSN() const                                                                              { return m_Star->IsECSN(); }
    bool                        IsMassRatioUnstable(const double p_AccretorMass, const double p_IsAccretorDegenerate) const { return m_Star->IsMassRatioUnstable(p_AccretorMass, p_IsAccretorDegenerate); }
    bool                        IsOneOf(STELLAR_TYPE_LIST p_List) const                                                     { return m_Star->IsOneOf(p_List); }
    bool                        IsPISN() const                                                                              { return m_Star->IsPISN(); }
    bool                        IsPPISN() const                                                                             { return m_Star->IsPPISN(); }
    bool                        IsUSSN() const                                                                              { return m_Star->IsUSSN(); }
    double                      Lambda_Fixed() const                                                                        { return m_Star->Lambda_Fixed(); }
    double                      Lambda_Loveridge() const                                                                    { return m_Star->Lambda_Loveridge(); }
    double                      Lambda_Nanjing() const                                                                      { return m_Star->Lambda_Nanjing(); }
    double                      Lambda_Kruckow() const                                                                      { return m_Star->Lambda_Kruckow(); }
    double                      Luminosity() const                                                                          { return m_Star->Luminosity(); }
    double                      Mass() const                                                                                { return m_Star->Mass(); }
    double                      Mass0() const                                                                               { return m_Star->Mass0(); }
    double                      MassPrev() const                                                                            { return m_Star->MassPrev(); }
    double                      Metallicity() const                                                                         { return m_Star->Metallicity(); }
    double                      MZAMS() const                                                                               { return m_Star->MZAMS(); }
    double                      NuclearTimescale() const                                                                    { return m_Star->NuclearTimescale(); }
    double                      Omega() const                                                                               { return m_Star->Omega(); }
    double                      OmegaCHE() const                                                                            { return m_Star->OmegaCHE(); }
    double                      OmegaPrev() const                                                                           { return m_Star->OmegaPrev(); }
    double                      RadialExpansionTimescale() const                                                            { return m_Star->RadialExpansionTimescale(); }
    double                      Radius() const                                                                              { return m_Star->Radius(); }
    double                      RadiusPrev() const                                                                          { return m_Star->RadiusPrev(); }
    double                      RZAMS() const                                                                               { return m_Star->RZAMS(); }
    SupernovaDetailsT           SN_Details() const                                                                          { return m_Star->SN_Details(); }
    double                      SN_Phi() const                                                                              { return m_Star->SN_Phi(); }
    double                      SN_Theta() const                                                                            { return m_Star->SN_Theta(); }
    double                      SN_TrueAnomaly() const                                                                      { return m_Star->SN_TrueAnomaly(); }
    SN_EVENT                    SN_Type() const                                                                             { return m_Star->SN_Type(); }
    COMPAS_VARIABLE             StellarPropertyValue(const T_ANY_PROPERTY p_Property) const                                 { return m_Star->StellarPropertyValue(p_Property); }
    double                      Temperature() const                                                                         { return m_Star->Temperature(); }
    double                      ThermalTimescale() const                                                                    { return m_Star->ThermalTimescale(); }
    double                      Timescale(TIMESCALE p_Timescale) const                                                      { return m_Star->Timescale(p_Timescale); }
    double                      XExponent() const                                                                           { return m_Star->XExponent(); }


    // setters (JR: I don't really like this, but I think unavoidable - at least for now)
    void                        SetOmega(double p_vRot)                                                                     { m_Star->SetOmega(p_vRot); }


    // member functions - alphabetically
    STELLAR_TYPE    AgeOneTimestep(const double p_Dt, bool p_Switch = true);

    void            ApplyMassTransferRejuvenationFactor()                                                       { m_Star->ApplyMassTransferRejuvenationFactor(); }

    void            CalculateAngularMomentum()                                                                  { m_Star->CalculateAngularMomentum(); }

    void            CalculateBindingEnergies(const double p_CoreMass,
                                             const double p_EnvMass,
                                             const double p_Radius)                                             { m_Star->CalculateBindingEnergies(p_CoreMass, p_EnvMass, p_Radius); }
    
    double          CalculateDynamicalMassLossRate()                                                            { return m_Star->CalculateDynamicalMassLossRate(); }

    double          CalculateDynamicalTimescale() const                                                         { return m_Star->CalculateDynamicalTimescale(); }

    double          CalculateEddyTurnoverTimescale()                                                            { return m_Star->CalculateEddyTurnoverTimescale(); }

    double          CalculateGyrationRadius()                                                                   { return m_Star->CalculateGyrationRadius(); }

    void            CalculateLambdas()                                                                          { m_Star->CalculateLambdas(); }
    void            CalculateLambdas(const double p_EnvMass)                                                    { m_Star->CalculateLambdas(p_EnvMass); }

    DBL_DBL         CalculateMassAcceptanceRate(const double p_DonorMassRate,
                                                const double p_FractionAccreted,
                                                const double p_AccretorMassRate)                                { return m_Star->CalculateMassAcceptanceRate(p_DonorMassRate,
                                                                                                                                                             p_FractionAccreted,
                                                                                                                                                             p_AccretorMassRate); }
    DBL_DBL         CalculateMassAcceptanceRate(const double p_DonorMassRate,
                                                const double p_FractionAccreted)                                { return m_Star->CalculateMassAcceptanceRate(p_DonorMassRate,
                                                                                                                                                             p_FractionAccreted,
                                                                                                                                                             CalculateThermalMassLossRate()); }

    double          CalculateMassLossValues(const bool p_UpdateMDot = false, const bool p_UpdateMDt = false)    { return m_Star->CalculateMassLossValues(p_UpdateMDot, p_UpdateMDt); }

    double          CalculateMomentOfInertia(const double p_RemnantRadius = 0.0)                                { return m_Star->CalculateMomentOfInertia(p_RemnantRadius); }
    double          CalculateMomentOfInertiaAU(const double p_RemnantRadius = 0.0)                              { return m_Star->CalculateMomentOfInertiaAU(p_RemnantRadius); }

    void            CalculateSNAnomalies(const double p_Eccentricity)                                           { m_Star->CalculateSNAnomalies(p_Eccentricity); }
    double          CalculateSNKickVelocity(const double p_RemnantMass, const double p_EjectaMass, 
								const STELLAR_TYPE p_StellarType)               { return m_Star->CalculateSNKickVelocity(p_RemnantMass, 
																p_EjectaMass, p_StellarType); }

    double          CalculateThermalMassLossRate()                                                              { return m_Star->CalculateThermalMassLossRate(); }

    double          CalculateThermalTimescale() const                                                           { return m_Star->CalculateThermalTimescale(); }
    double          CalculateThermalTimescale(const double p_Mass,
                                              const double p_Radius,
                                              const double p_Luminosity,
                                              const double p_EnvMass = 1.0)                                     { return m_Star->CalculateThermalTimescale(p_Mass, p_Radius, p_Luminosity, p_EnvMass); }

    double          CalculateTimestep()                                                                         { return m_Star->CalculateTimestep(); }

    double          CalculateZeta(CE_ZETA_PRESCRIPTION p_CEZetaPrescription)                                    { return m_Star->CalculateZeta(p_CEZetaPrescription); }

    void            CalculateZetas()                                                                            { m_Star->CalculateZetas(); }

    void            CheckRunaway(const bool p_Unbound)                                                          { m_Star->CheckRunaway(p_Unbound); }

    void            ClearCurrentSNEvent()                                                                       { m_Star->ClearCurrentSNEvent(); }

    BaseStar*       Clone(const BaseStar& p_Star);

    ENVELOPE        DetermineEnvelopeType() const                                                               { return m_Star->DetermineEnvelopeType(); }
    ENVELOPE        DetermineEnvelopeTypeHurley2002() const                                                     { return m_Star->DetermineEnvelopeTypeHurley2002(); }

    MT_CASE         DetermineMassTransferCase()                                                                 { return m_Star->DetermineMassTransferCase(); }

    void            Evolve(const int p_StepNum);

    double          EvolveOneTimestep(const double p_Dt);

    void            IncrementOmega(const double p_OmegaDelta)                                                   { m_Star->IncrementOmega(p_OmegaDelta); }

    void            ResolveAccretion(const double p_AccretionMass)                                              { m_Star->ResolveAccretion(p_AccretionMass); }

    void            ResolveEnvelopeLossAndSwitch()                                                              { SwitchTo(m_Star->ResolveEnvelopeLoss(true)); }

    void            ResolveRemnantAfterEnvelopeLossAndSwitch()                                                  { SwitchTo(m_Star->ResolveRemnantAfterEnvelopeLoss()); }

    STELLAR_TYPE    ResolveRemnantAfterEnvelopeLoss()                                                           { return m_Star->ResolveRemnantAfterEnvelopeLoss(); }

    bool            RevertState();

    void            SaveState();

    void            SetSNCurrentEvent(const SN_EVENT p_SNEvent)                                                 { m_Star->SetSNCurrentEvent(p_SNEvent); }
    void            SetSNPastEvent(const SN_EVENT p_SNEvent)                                                    { m_Star->SetSNPastEvent(p_SNEvent); }

    double     	    SN_KickVelocity()       									{ return m_Star->SN_KickVelocity() ; }

    void            SwitchTo(const STELLAR_TYPE p_StellarType, bool p_SetInitialType = false);

    void            UpdateAgeAfterMassLoss()                                                                    { m_Star->UpdateAgeAfterMassLoss(); }

    void            UpdateAttributes()                                                                          { (void)UpdateAttributes(0.0, 0.0, true); }
    STELLAR_TYPE    UpdateAttributes(const double p_DeltaMass,
                                     const double p_DeltaMass0,
                                     const bool   p_ForceRecalculate = false);

    STELLAR_TYPE    UpdateAttributesAndAgeOneTimestep(const double p_DeltaMass,
                                                      const double p_DeltaMass0,
                                                      const double p_DeltaTime,
                                                      const bool   p_Switch = true,
                                                      const bool   p_ForceRecalculate = false);

    void            UpdateInitialMass()                                                                         { m_Star->UpdateInitialMass(); }

    void            UpdateMagneticFieldAndSpin(const bool   p_CommonEnvelope,
                                               const double p_Stepsize,
                                               const double p_MassGainPerTimeStep,
                                               const double p_Epsilon)                                          { m_Star->UpdateMagneticFieldAndSpin(p_CommonEnvelope,
                                                                                                                                                     p_Stepsize,
                                                                                                                                                     p_MassGainPerTimeStep,
                                                                                                                                                     p_Epsilon);}










private:

    OBJECT_ID   m_ObjectId;                                                                                     // instantiated object's unique object id
    OBJECT_TYPE m_ObjectType;                                                                                   // instantiated object's object type

    BaseStar   *m_Star;                                                                                         // pointer to current star
    BaseStar   *m_SaveStar;                                                                                     // pointer to saved star

};

#endif // __Star_h__
