#ifndef __Star_h__
#define __Star_h__

#include <fstream>

#include "constants.h"
#include "typedefs.h"
#include "profiling.h"
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
         const KickParameters    p_KickParameters,
         const double            p_RotationalVelocity = -1.0); 

    Star(const Star& p_Star);

    virtual ~Star() { delete m_Star; delete m_SaveStar; }

/*
 * Switch to required star type
 *
 * Instantiates new object of required class, deletes existing pointer to star object and
 * replaces it with pointer to newly instantiated object
 *
 *
 * STELLAR_TYPE SwitchTo(const STELLAR_TYPE p_StellarType, bool p_SetInitialState)
 *
 * @param   [IN]    p_StellarType               StellarType to switch to
 * @param   [IN]    p_SetInitialType            Indicates whether the initial stellar type of the star should be set to p_StellarType
 *                                              (optional, default = false)
 * @return                                      Stellar type of star before switch (previous stellar type)
 */

BaseStar* NewStar(const STELLAR_TYPE p_StellarType) {
std::cout << "NewStar@1, ST = " << (int)p_StellarType << "\n";

        BaseStar *ptr;

        switch (p_StellarType) {
            case STELLAR_TYPE::MS_LTE_07                                : {ptr = new MS_lte_07();} break;
            case STELLAR_TYPE::MS_GT_07                                 : {ptr = new MS_gt_07();} break;
            case STELLAR_TYPE::CHEMICALLY_HOMOGENEOUS                   : {ptr = new CH();} break;
            case STELLAR_TYPE::HERTZSPRUNG_GAP                          : {ptr = new HG();} break;
            case STELLAR_TYPE::FIRST_GIANT_BRANCH                       : {ptr = new FGB();} break;
            case STELLAR_TYPE::CORE_HELIUM_BURNING                      : {ptr = new CHeB();} break;
            case STELLAR_TYPE::EARLY_ASYMPTOTIC_GIANT_BRANCH            : {ptr = new EAGB();} break;
            case STELLAR_TYPE::THERMALLY_PULSING_ASYMPTOTIC_GIANT_BRANCH: {ptr = new TPAGB();} break;
            case STELLAR_TYPE::NAKED_HELIUM_STAR_MS                     : {ptr = new HeMS();} break;
            case STELLAR_TYPE::NAKED_HELIUM_STAR_HERTZSPRUNG_GAP        : {ptr = new HeHG();} break;
            case STELLAR_TYPE::NAKED_HELIUM_STAR_GIANT_BRANCH           : {ptr = new HeGB();} break;
            case STELLAR_TYPE::HELIUM_WHITE_DWARF                       : {ptr = new HeWD();} break;
            case STELLAR_TYPE::CARBON_OXYGEN_WHITE_DWARF                : {ptr = new COWD();} break;
            case STELLAR_TYPE::OXYGEN_NEON_WHITE_DWARF                  : {ptr = new ONeWD();} break;
            case STELLAR_TYPE::NEUTRON_STAR                             : {ptr = new NS();} break;
            case STELLAR_TYPE::BLACK_HOLE                               : {ptr = new BH();} break;
            case STELLAR_TYPE::MASSLESS_REMNANT                         : {ptr = new MR();} break;

            default:                                                // not ok... (this can only happen if someone added a new stellar type)
                THROW_ERROR_STATIC(ERROR::UNKNOWN_STELLAR_TYPE);    // throw error
        }
std::cout << "NewStar@2, ST = " << (int)(ptr->StellarType()) << "\n";

    return ptr;
}

    // object identifiers - all classes have these
    OBJECT_ID           ObjectId() const                                                                            { return m_ObjectId; }
    OBJECT_TYPE         ObjectType() const                                                                          { return OBJECT_TYPE::STAR; }
    OBJECT_PERSISTENCE  ObjectPersistence() const                                                                   { return m_ObjectPersistence; }
    STELLAR_TYPE        InitialStellarType() const                                                                  { return m_Star->InitialStellarType(); }
    STELLAR_TYPE        StellarType() const                                                                         { return m_Star->StellarType(); }

    // getters - alphabetically
    double              Age() const                                                                                 { return m_Star->Age(); }
    double              BindingEnergy_Fixed() const                                                                 { return m_Star->BindingEnergy_Fixed(); }
    double              BindingEnergy_Loveridge() const                                                             { return m_Star->BindingEnergy_Loveridge(); }
    double              BindingEnergy_Nanjing() const                                                               { return m_Star->BindingEnergy_Nanjing(); }
    double              BindingEnergy_Kruckow() const                                                               { return m_Star->BindingEnergy_Kruckow(); }
    double              BindingEnergy_Dewi() const                                                                  { return m_Star->BindingEnergy_Dewi(); }
    double              CalculateCriticalMassRatio(const bool p_AccretorIsDegenerate, 
                                                   const double p_massTransferEfficiencyBeta) const                 { return m_Star->CalculateCriticalMassRatio(p_AccretorIsDegenerate, p_massTransferEfficiencyBeta); }
    double              CalculateCriticalMassRatioClaeys14(const bool p_AccretorIsDegenerate) const                 { return m_Star->CalculateCriticalMassRatioClaeys14(p_AccretorIsDegenerate); }
    double              CalculateCriticalMassRatioHurleyHjellmingWebbink() const                                    { return m_Star->CalculateCriticalMassRatioHurleyHjellmingWebbink(); }
    double              CalculateDynamicalTimescale() const                                                         { return m_Star->CalculateDynamicalTimescale(); }
    double              CalculateRadialExpansionTimescale() const                                                   { return m_Star->CalculateRadialExpansionTimescale(); }
    double              CalculateThermalTimescale() const                                                           { return m_Star->CalculateThermalTimescale(); }
    double              COCoreMass() const                                                                          { return m_Star->COCoreMass(); }
    double              CoreMass() const                                                                            { return m_Star->CoreMass(); }
    bool                EnvelopeJustExpelledByPulsations() const                                                    { return m_Star->EnvelopeJustExpelledByPulsations(); }
    bool                ExperiencedAIC() const                                                                      { return m_Star->ExperiencedAIC(); }
    bool                ExperiencedCCSN() const                                                                     { return m_Star->ExperiencedCCSN(); }
    bool                ExperiencedECSN() const                                                                     { return m_Star->ExperiencedECSN(); }
    bool                ExperiencedPISN() const                                                                     { return m_Star->ExperiencedPISN() ; }
    bool                ExperiencedPPISN() const                                                                    { return m_Star->ExperiencedPPISN(); }
    bool                ExperiencedUSSN() const                                                                     { return m_Star->ExperiencedUSSN(); }
    double              HeCoreMass() const                                                                          { return m_Star->HeCoreMass(); }
    bool                IsAIC() const                                                                               { return m_Star->IsAIC(); }
    bool                IsCCSN() const                                                                              { return m_Star->IsCCSN(); }
    bool                IsDegenerate() const                                                                        { return m_Star->IsDegenerate(); }
    bool                IsECSN() const                                                                              { return m_Star->IsECSN(); }
    bool                IsHeSD() const                                                                              { return m_Star->IsHeSD(); }
    bool                IsOneOf(STELLAR_TYPE_LIST p_List) const                                                     { return m_Star->IsOneOf(p_List); }
    bool                IsPISN() const                                                                              { return m_Star->IsPISN(); }
    bool                IsPPISN() const                                                                             { return m_Star->IsPPISN(); }
    bool                IsSNIA() const                                                                              { return m_Star->IsSNIA(); }
    bool                IsSNII() const                                                                              { return m_Star->IsSNII(); }
    bool                IsUSSN() const                                                                              { return m_Star->IsUSSN(); }
    double              Lambda_Fixed() const                                                                        { return m_Star->Lambda_Fixed(); }
    double              Lambda_Loveridge() const                                                                    { return m_Star->Lambda_Loveridge(); }
    double              Lambda_Nanjing() const                                                                      { return m_Star->Lambda_Nanjing(); }
    double              Lambda_Kruckow() const                                                                      { return m_Star->Lambda_Kruckow(); }
    double              Lambda_Dewi() const                                                                         { return m_Star->Lambda_Dewi(); }
    double              Luminosity() const                                                                          { return m_Star->Luminosity(); }
    double              Mass() const                                                                                { return m_Star->Mass(); }
    double              Mass0() const                                                                               { return m_Star->Mass0(); }
    double              MassPrev() const                                                                            { return m_Star->MassPrev(); }
    double              Metallicity() const                                                                         { return m_Star->Metallicity(); }
    double              MZAMS() const                                                                               { return m_Star->MZAMS(); }
    double              Omega() const                                                                               { return m_Star->Omega(); }
    double              OmegaCHE() const                                                                            { return m_Star->OmegaCHE(); }
    double              OmegaPrev() const                                                                           { return m_Star->OmegaPrev(); }
    double              Radius() const                                                                              { return m_Star->Radius(); }
    double              RadiusPrev() const                                                                          { return m_Star->RadiusPrev(); }
    unsigned long int   RandomSeed() const                                                                          { return m_Star->RandomSeed(); }
    double              RZAMS() const                                                                               { return m_Star->RZAMS(); }
    SupernovaDetailsT   SN_Details() const                                                                          { return m_Star->SN_Details(); }
    double              SN_Phi() const                                                                              { return m_Star->SN_Phi(); }
    double              SN_Theta() const                                                                            { return m_Star->SN_Theta(); }
    double              SN_TotalMassAtCOFormation() const                                                           { return m_Star->SN_TotalMassAtCOFormation(); }
    double              SN_TrueAnomaly() const                                                                      { return m_Star->SN_TrueAnomaly(); }
    double              SN_EccentricAnomaly() const                                                                 { return m_Star->SN_EccentricAnomaly(); }
    SN_EVENT            SN_Type() const                                                                             { return m_Star->SN_Type(); }
    double              Speed() const                                                                               { return m_Star->Speed(); }
    COMPAS_VARIABLE     StellarPropertyValue(const T_ANY_PROPERTY p_Property) const                                 { return m_Star->StellarPropertyValue(p_Property); }
    STELLAR_TYPE        StellarTypePrev() const                                                                     { return m_Star->StellarTypePrev(); }
    double              Tau() const                                                                                 { return m_Star->Tau(); }
    double              Temperature() const                                                                         { return m_Star->Temperature(); }
    double              Timescale(TIMESCALE p_Timescale) const                                                      { return m_Star->Timescale(p_Timescale); }
    double              XExponent() const                                                                           { return m_Star->XExponent(); }

    
    // setters
    void                SetOmega(double p_vRot)                                                                     { m_Star->SetOmega(p_vRot); }
    void                SetObjectId(const OBJECT_ID p_ObjectId)                                                     { m_ObjectId = p_ObjectId; }
    void                SetPersistence(const OBJECT_PERSISTENCE p_Persistence)                                      { m_ObjectPersistence = p_Persistence; }
    void                UpdateMassTransferDonorHistory()                                                            { m_Star->UpdateMassTransferDonorHistory(); }
    void                ResetEnvelopeExpulsationByPulsations()                                                      { m_Star->ResetEnvelopeExpulsationByPulsations(); }


    // member functions - alphabetically
    STELLAR_TYPE    AgeOneTimestep(const double p_Dt, bool p_Switch = true);

    void            ApplyMassTransferRejuvenationFactor()                                                           { m_Star->ApplyMassTransferRejuvenationFactor(); }

    void            CalculateBindingEnergies(const double p_CoreMass,
                                             const double p_EnvMass,
                                             const double p_Radius)                                                 { m_Star->CalculateBindingEnergies(p_CoreMass, p_EnvMass, p_Radius); }

    double          CalculateConvectiveCoreMass()                                                                   { return m_Star->CalculateConvectiveCoreMass(); }
    double          CalculateConvectiveCoreRadius()                                                                 { return m_Star->CalculateConvectiveCoreRadius(); }

    double          CalculateConvectiveEnvelopeBindingEnergy(const double p_TotalMass,
                                                             const double p_ConvectiveEnvelopeMass,
                                                             const double p_Radius,
                                                             const double p_Lambda)                                 { return m_Star->CalculateConvectiveEnvelopeBindingEnergy(p_TotalMass, p_ConvectiveEnvelopeMass, p_Radius, p_Lambda); }
    double          CalculateConvectiveEnvelopeLambdaPicker(const double p_convectiveEnvelopeMass, const double p_maxConvectiveEnvelopeMass ) const     { return m_Star->CalculateConvectiveEnvelopeLambdaPicker(p_convectiveEnvelopeMass, p_maxConvectiveEnvelopeMass); }
    DBL_DBL         CalculateConvectiveEnvelopeMass()                                                               { return m_Star->CalculateConvectiveEnvelopeMass(); }
    
    double          CalculateEddyTurnoverTimescale()                                                                { return m_Star->CalculateEddyTurnoverTimescale(); }
    
    DBL_DBL_DBL_DBL CalculateImKlmDynamical(const double p_Omega, const double p_SemiMajorAxis, const double p_M2)  { return m_Star->CalculateImKlmDynamical(p_Omega, p_SemiMajorAxis, p_M2); }
    DBL_DBL_DBL_DBL CalculateImKlmEquilibrium(const double p_Omega, const double p_SemiMajorAxis, const double p_M2){ return m_Star->CalculateImKlmEquilibrium(p_Omega, p_SemiMajorAxis, p_M2); }
    DBL_DBL_DBL_DBL CalculateImKlmTidal(const double p_Omega, const double p_SemiMajorAxis, const double p_M2)      { return m_Star->CalculateImKlmTidal(p_Omega, p_SemiMajorAxis, p_M2); }
    
    void            CalculateLambdas()                                                                              { m_Star->CalculateLambdas(); }
    void            CalculateLambdas(const double p_EnvMass)                                                        { m_Star->CalculateLambdas(p_EnvMass); }

    DBL_DBL         CalculateMassAcceptanceRate(const double p_DonorMassRate, 
                                                const double p_AccretorMassRate,
                                                const bool   p_IsHeRich)                                            { return m_Star->CalculateMassAcceptanceRate(p_DonorMassRate, p_AccretorMassRate, p_IsHeRich); }

    double          CalculateMassLossValues(const bool p_UpdateMDot = false, const bool p_UpdateMDt = false)        { return m_Star->CalculateMassLossValues(p_UpdateMDot, p_UpdateMDt); }

    double          CalculateMomentOfInertia() const                                                                { return m_Star->CalculateMomentOfInertia(); }
    double          CalculateMomentOfInertiaAU() const                                                              { return m_Star->CalculateMomentOfInertiaAU(); }
    
    double          CalculateNuclearMassLossRate()                                                                  { return m_Star->CalculateNuclearMassLossRate(); }
    
    double          CalculateRadialExtentConvectiveEnvelope() { return m_Star->CalculateRadialExtentConvectiveEnvelope(); }

    double          CalculateRadiusOnPhaseTau(const double p_Mass, const double p_Tau) const                        { return m_Star->CalculateRadiusOnPhaseTau(p_Mass, p_Tau); }
    
    void            CalculateSNAnomalies(const double p_Eccentricity)                                               { m_Star->CalculateSNAnomalies(p_Eccentricity); }
    
    double          CalculateSNKickMagnitude(const double p_RemnantMass, 
                                             const double p_EjectaMass, 
								             const STELLAR_TYPE p_StellarType)                                      { return m_Star->CalculateSNKickMagnitude(p_RemnantMass, p_EjectaMass, p_StellarType); }


    double          CalculateThermalMassAcceptanceRate(const double p_Radius)                                       { return m_Star->CalculateThermalMassAcceptanceRate(p_Radius); }
    double          CalculateThermalMassAcceptanceRate()                                                            { return m_Star->CalculateThermalMassAcceptanceRate(); }

    double          CalculateThermalMassLossRate() const                                                            { return m_Star->CalculateThermalMassLossRate(); }

    double          CalculateThermalTimescale(const double p_Radius) const                                          { return m_Star->CalculateThermalTimescale(p_Radius); }

    double          CalculateTimestep()                                                                             { return m_Star->CalculateTimestep(); }

    double          CalculateZetaAdiabatic()                                                                        { return m_Star->CalculateZetaAdiabatic(); }
    double          CalculateZetaConstantsByEnvelope(ZETA_PRESCRIPTION p_ZetaPrescription)                          { return m_Star->CalculateZetaConstantsByEnvelope(p_ZetaPrescription); }
    double          CalculateZetaEquilibrium()                                                                      { return m_Star->CalculateZetaEquilibrium(); }

    void            ClearCurrentSNEvent()                                                                           { m_Star->ClearCurrentSNEvent(); }

    ACCRETION_REGIME DetermineAccretionRegime(const bool p_HeRich,
                                              const double p_DonorThermalMassLossRate)                              { return m_Star->DetermineAccretionRegime(p_HeRich, p_DonorThermalMassLossRate); }  // Used in WDs

    ENVELOPE        DetermineEnvelopeType() const                                                                   { return m_Star->DetermineEnvelopeType(); }

    EVOLUTION_STATUS Evolve(const long int p_Id);

    double          EvolveOneTimestep(const double p_Dt, const bool p_Force = false);

    double          InterpolateGe20QCrit(const QCRIT_PRESCRIPTION p_qCritPrescription, 
                                         const double p_massTransferEfficiencyBeta)                                 { return m_Star->InterpolateGe20QCrit(p_qCritPrescription, p_massTransferEfficiencyBeta); }
    void            HaltWinds()                                                                                     { m_Star->HaltWinds(); }

    void            ResolveAccretion(const double p_AccretionMass)                                                  { m_Star->ResolveAccretion(p_AccretionMass); }

    void            ResolveAccretionRegime(const ACCRETION_REGIME p_Regime,
                                           const double p_DonorThermalMassLossRate)                                 { m_Star->ResolveAccretionRegime(p_Regime, p_DonorThermalMassLossRate); }  // Used in WDs

    void            ResolveEnvelopeLossAndSwitch()                                                                  { (void)SwitchTo(m_Star->ResolveEnvelopeLoss(true)); }

    void            ResolveShellChange(const double p_AccretedMass)                                                 { m_Star->ResolveShellChange(p_AccretedMass); }  // Used in WDs

    bool            RevertState();

    void            SaveState();

    void            SetSNCurrentEvent(const SN_EVENT p_SNEvent)                                                     { m_Star->SetSNCurrentEvent(p_SNEvent); }
    void            SetSNPastEvent(const SN_EVENT p_SNEvent)                                                        { m_Star->SetSNPastEvent(p_SNEvent); }

    double     	    SN_KickMagnitude()       									                                    { return m_Star->SN_KickMagnitude() ; }
    double     	    SN_RocketKickMagnitude()       									                                { return m_Star->SN_RocketKickMagnitude(); }
    double     	    SN_RocketKickPhi()       									                                    { return m_Star->SN_RocketKickPhi(); }
    double     	    SN_RocketKickTheta()       									                                    { return m_Star->SN_RocketKickTheta(); }
    

    STELLAR_TYPE    SwitchTo(const STELLAR_TYPE p_StellarType, bool p_SetInitialType = false);

    void            UpdateAgeAfterMassLoss()                                                                        { m_Star->UpdateAgeAfterMassLoss(); }

    void            UpdateAttributes()                                                                              { (void)UpdateAttributes(0.0, 0.0, true); }
    STELLAR_TYPE    UpdateAttributes(const double p_DeltaMass,
                                     const double p_DeltaMass0,
                                     const bool   p_ForceRecalculate = false);

    STELLAR_TYPE    UpdateAttributesAndAgeOneTimestep(const double p_DeltaMass,
                                                      const double p_DeltaMass0,
                                                      const double p_DeltaTime,
                                                      const bool   p_Switch = true,
                                                      const bool   p_ForceRecalculate = false);

    void            UpdateComponentVelocity(const Vector3d p_newVelocity)                                           { m_Star->UpdateComponentVelocity(p_newVelocity); }

    void            UpdateInitialMass()                                                                             { m_Star->UpdateInitialMass(); }

    void            UpdateMagneticFieldAndSpin(const bool   p_CommonEnvelope,
                                               const bool   p_RecycledNS,
                                               const double p_Stepsize,
                                               const double p_MassGainPerTimeStep,
                                               const double p_Epsilon)                                              { m_Star->UpdateMagneticFieldAndSpin(p_CommonEnvelope,
                                                                                                                                                         p_RecycledNS,
                                                                                                                                                         p_Stepsize,
                                                                                                                                                         p_MassGainPerTimeStep,
                                                                                                                                                         p_Epsilon);}

    void            UpdateMinimumCoreMass()                                                                         { m_Star->UpdateMinimumCoreMass(); }
    
    void            UpdatePreviousTimestepDuration()                                                                { m_Star->UpdatePreviousTimestepDuration(); }
    
    ACCRETION_REGIME WhiteDwarfAccretionRegime() const                                                              { return m_Star->WhiteDwarfAccretionRegime(); }

private:

    OBJECT_ID          m_ObjectId;              // instantiated object's unique object id
    OBJECT_PERSISTENCE m_ObjectPersistence;     // instantiated object's persistence
    long int           m_Id;                    // id used to name output files - uses p_Id as passed (usually the step number of multiple single stars being produced)

    ERROR              m_Error;                 // Records most recent error encountered for this star - should be propagated to m_Star

    BaseStar          *m_Star;                  // pointer to current star
    BaseStar          *m_SaveStar;              // pointer to saved star

    std::vector<double> m_Timesteps;            // timesteps vector - for debugging/testing

};

#endif // __Star_h__
