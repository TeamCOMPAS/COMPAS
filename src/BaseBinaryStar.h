#ifndef __BaseBinaryStar_h__
#define __BaseBinaryStar_h__

#include "constants.h"
#include "typedefs.h"
#include "utils.h"

#include "Log.h"
#include "Star.h"
#include "AIS.h"
#include "BinaryConstituentStar.h"


class Log;
class Star;
class AIS;
class BinaryConstituentStar;


class BaseBinaryStar {

public:

    BaseBinaryStar(const AIS &p_AIS, const long int p_Id = -1l);

    BaseBinaryStar(const AIS     &p_AIS,
                   const double   p_Mass1,
                   const double   p_Mass2,
                   const double   p_Metallicity1,
                   const double   p_Metallicity2,
                   const double   p_SemiMajorAxis,
                   const double   p_Eccentricity,
                   const long int p_Id = -1l);


    void CopyMemberVariables(const BaseBinaryStar& p_Star) {

        m_Error                            = p_Star.m_Error;

        m_RandomSeed                       = p_Star.m_RandomSeed;

        m_AIS                              = p_Star.m_AIS;

        m_BeBinaryDetails                  = p_Star.m_BeBinaryDetails;

        m_BeBinaryDetails.currentProps     = p_Star.m_BeBinaryDetails.currentProps  == &(p_Star.m_BeBinaryDetails.props1) ? &(m_BeBinaryDetails.props1) : &(m_BeBinaryDetails.props2);
        m_BeBinaryDetails.previousProps    = p_Star.m_BeBinaryDetails.previousProps == &(p_Star.m_BeBinaryDetails.props1) ? &(m_BeBinaryDetails.props1) : &(m_BeBinaryDetails.props2);

        m_Beta                             = p_Star.m_Beta;

        m_CircularizationTimescale         = p_Star.m_CircularizationTimescale;

        m_CommonEnvelope                   = p_Star.m_CommonEnvelope;
        m_CommonEnvelopeAlpha              = p_Star.m_CommonEnvelopeAlpha;
        m_CommonEnvelopeOccuredAtLeastOnce = p_Star.m_CommonEnvelopeOccuredAtLeastOnce;
        m_CounterCEE                       = p_Star.m_CounterCEE;
        m_DoubleCoreCommonEnvelope         = p_Star.m_DoubleCoreCommonEnvelope;
        m_OptimisticCommonEnvelope         = p_Star.m_OptimisticCommonEnvelope;

        m_Disbound                         = p_Star.m_Disbound;

        m_Dt                               = p_Star.m_Dt;

        m_EPrime                           = p_Star.m_EPrime;

        m_Eccentricity                     = p_Star.m_Eccentricity;
        m_EccentricityAtDCOFormation       = p_Star.m_EccentricityAtDCOFormation;
        m_EccentricityInitial              = p_Star.m_EccentricityInitial;
        m_EccentricityPostCEE              = p_Star.m_EccentricityPostCEE;
        m_EccentricityPre2ndSN             = p_Star.m_EccentricityPre2ndSN;
        m_EccentricityPreCEE               = p_Star.m_EccentricityPreCEE;
        m_EccentricityPrev                 = p_Star.m_EccentricityPrev;
        m_EccentricityPrime                = p_Star.m_EccentricityPrime;

        m_FastPhaseCaseA                   = p_Star.m_FastPhaseCaseA;

        m_FractionAccreted                 = p_Star.m_FractionAccreted;

        m_CosIPrime                        = p_Star.m_CosIPrime;
        m_IPrime                           = p_Star.m_IPrime;

        m_JLoss                            = p_Star.m_JLoss;

        m_LBVfactor                        = p_Star.m_LBVfactor;

        m_MassesEquilibrated               = p_Star.m_MassesEquilibrated;

        m_Mass1Final                       = p_Star.m_Mass1Final;
        m_Mass2Final                       = p_Star.m_Mass2Final;

        m_Mass1PreCEE                      = p_Star.m_Mass1PreCEE;
        m_Mass2PreCEE                      = p_Star.m_Mass2PreCEE;

        m_MassEnv1                         = p_Star.m_MassEnv1;
        m_MassEnv2                         = p_Star.m_MassEnv2;

        m_aMassLossDiff                    = p_Star.m_aMassLossDiff;
        m_OmegaMassLossDiff                = p_Star.m_OmegaMassLossDiff;

        m_MassTransfer                     = p_Star.m_MassTransfer;
        m_aMassTransferDiff                = p_Star.m_aMassTransferDiff;
        m_OmegaMassTransferDiff            = p_Star.m_OmegaMassTransferDiff;

        m_MassTransferTrackerHistory       = p_Star.m_MassTransferTrackerHistory;

        m_ReducedMassPrev                  = p_Star.m_ReducedMassPrev;
        m_ReducedMassPrime                 = p_Star.m_ReducedMassPrime;

        m_TotalMassPrev                    = p_Star.m_TotalMassPrev;
        m_TotalMassPrime                   = p_Star.m_TotalMassPrime;

        m_MC                               = p_Star.m_MC;
        m_MCPrime                          = p_Star.m_MCPrime;

        m_Merged                           = p_Star.m_Merged;
        m_MergesInHubbleTime               = p_Star.m_MergesInHubbleTime;

        m_MSN                              = p_Star.m_MSN;
        m_MSNPrime                         = p_Star.m_MSNPrime;

        m_OrbitalVelocity                  = p_Star.m_OrbitalVelocity;
        m_OrbitalVelocityPre2ndSN          = p_Star.m_OrbitalVelocityPre2ndSN;
        m_OrbitalVelocityPrev              = p_Star.m_OrbitalVelocityPrev;
        m_OrbitalVelocityPrime             = p_Star.m_OrbitalVelocityPrime;

        m_Radius                           = p_Star.m_Radius;

        m_Radius1PostStripping             = p_Star.m_Radius1PostStripping;
        m_Radius2PostStripping             = p_Star.m_Radius2PostStripping;

        m_Radius1PreCEE                    = p_Star.m_Radius1PreCEE;
        m_Radius2PreCEE                    = p_Star.m_Radius2PreCEE;

        m_ImmediateRLOFPostCEE             = p_Star.m_ImmediateRLOFPostCEE;
        m_RLOFDetails                      = p_Star.m_RLOFDetails;
        m_SimultaneousRLOF                 = p_Star.m_SimultaneousRLOF;

        m_RocheLobe1to2PreCEE              = p_Star.m_RocheLobe1to2PreCEE;
        m_RocheLobe1to2PostCEE             = p_Star.m_RocheLobe1to2PostCEE;
        m_RocheLobe2to1PreCEE              = p_Star.m_RocheLobe2to1PreCEE;
        m_RocheLobe2to1PostCEE             = p_Star.m_RocheLobe2to1PostCEE;

        m_SemiMajorAxis                    = p_Star.m_SemiMajorAxis;
        m_SemiMajorAxisAtDCOFormation      = p_Star.m_SemiMajorAxisAtDCOFormation;
        m_SemiMajorAxisInitial             = p_Star.m_SemiMajorAxisInitial;
        m_SemiMajorAxisPostCEE             = p_Star.m_SemiMajorAxisPostCEE;
        m_SemiMajorAxisPre2ndSN            = p_Star.m_SemiMajorAxisPre2ndSN;
        m_SemiMajorAxisPreCEE              = p_Star.m_SemiMajorAxisPreCEE;
        m_SemiMajorAxisPrev                = p_Star.m_SemiMajorAxisPrev;
        m_SemiMajorAxisPrime               = p_Star.m_SemiMajorAxisPrime;

        m_StellarMerger                    = p_Star.m_StellarMerger;

        m_StellarType1PreCEE               = p_Star.m_StellarType1PreCEE;
        m_StellarType2PreCEE               = p_Star.m_StellarType2PreCEE;

        m_SupernovaState                   = p_Star.m_SupernovaState;

        m_Survived                         = p_Star.m_Survived;

        m_SynchronizationTimescale         = p_Star.m_SynchronizationTimescale;

        m_SystemicVelocity                 = p_Star.m_SystemicVelocity;

        m_Time                             = p_Star.m_Time;
        m_TimePrev                         = p_Star.m_TimePrev;
        m_TimeToCoalescence                = p_Star.m_TimeToCoalescence;

        m_TotalAngularMomentumPrev         = p_Star.m_TotalAngularMomentumPrev;
        m_TotalAngularMomentumPrime        = p_Star.m_TotalAngularMomentumPrime;

        m_TotalEnergy                      = p_Star.m_TotalEnergy;
        m_TotalEnergyPrime                 = p_Star.m_TotalEnergyPrime;

        m_TotalOrbitalAngularMomentumPrev  = p_Star.m_TotalOrbitalAngularMomentumPrev;
        m_TotalOrbitalAngularMomentumPrime = p_Star.m_TotalOrbitalAngularMomentumPrime;

        m_TotalOrbitalEnergyPrev           = p_Star.m_TotalOrbitalEnergyPrev;
        m_TotalOrbitalEnergyPrime          = p_Star.m_TotalOrbitalEnergyPrime;

        m_uK                               = p_Star.m_uK;

        m_VRel                             = p_Star.m_VRel;

        m_WolfRayetFactor                  = p_Star.m_WolfRayetFactor;

        m_ZetaRLOFAnalytic                 = p_Star.m_ZetaRLOFAnalytic;
        m_ZetaRLOFNumerical                = p_Star.m_ZetaRLOFNumerical;
        m_ZetaStarCompare                  = p_Star.m_ZetaStarCompare;

        // copy the constituent stars and pointers

        m_Star1     = p_Star.m_Star1 ? new BinaryConstituentStar(*(p_Star.m_Star1)) : nullptr;
        m_Star2     = p_Star.m_Star2 ? new BinaryConstituentStar(*(p_Star.m_Star2)) : nullptr;

        m_Donor     = p_Star.m_Donor    ? (p_Star.m_Donor    == p_Star.m_Star1 ? m_Star1 : m_Star2) : nullptr;
        m_Accretor  = p_Star.m_Accretor ? (p_Star.m_Accretor == p_Star.m_Star1 ? m_Star1 : m_Star2) : nullptr;

        m_Supernova = p_Star.m_Supernova ? (p_Star.m_Supernova == p_Star.m_Star1 ? m_Star1 : m_Star2) : nullptr;
        m_Companion = p_Star.m_Companion ? (p_Star.m_Companion == p_Star.m_Star1 ? m_Star1 : m_Star2) : nullptr;

    }

    // Copy constructor
    BaseBinaryStar(const BaseBinaryStar& p_Star) {

        m_ObjectId    = globalObjectId++;                           // get unique object id (don't copy source)
        m_ObjectType  = OBJECT_TYPE::BASE_BINARY_STAR;              // can only copy from BASE_BINARY_STAR
        m_StellarType = STELLAR_TYPE::BINARY_STAR;                  // always

        CopyMemberVariables(p_Star);                                // copy member variables
    }


    // Assignment overload
    BaseBinaryStar& operator = (const BaseBinaryStar& p_Star) {

        if (this != &p_Star) {                                      // make sure we're not not copying ourselves...

            m_ObjectId    = globalObjectId++;                      // get unique object id (don't copy source)
            m_ObjectType  = OBJECT_TYPE::BASE_BINARY_STAR;         // can only copy from BASE_BINARY_STAR
            m_StellarType = STELLAR_TYPE::BINARY_STAR;             // always

            CopyMemberVariables(p_Star);                            // copy member variables
        }
        return *this;
    }


    virtual ~BaseBinaryStar() { delete m_Star1; delete m_Star2; }


    // object identifiers - all classes have these
    OBJECT_ID           ObjectId() const                            { return m_ObjectId; }
    OBJECT_TYPE         ObjectType() const                          { return m_ObjectType; }
    STELLAR_TYPE        StellarType() const                         { return m_StellarType; }


    // getters - alphabetically
    BeBinaryDetailsT    BeBinaryDetails() const                     { return m_BeBinaryDetails; }
	double              CEAlpha() const                             { return m_CommonEnvelopeAlpha; }
	bool                CEAtLeastOnce() const                       { return m_CommonEnvelopeOccuredAtLeastOnce; }
	double              CircularizationTimescale() const            { return m_CircularizationTimescale; }
    unsigned int        CommonEnvelopeEventCount() const            { return m_CounterCEE; }
    bool                Disbound() const                            { return m_Disbound; }
    bool                DoubleCoreCE() const                        { return m_DoubleCoreCommonEnvelope; }
    double              Dt() const                                  { return m_Dt; }
    double              Eccentricity() const                        { return m_Eccentricity; }
    double              EccentricityAtDCOFormation() const          { return m_EccentricityAtDCOFormation; }
    double              EccentricityInitial() const                 { return m_EccentricityInitial; }
    double              EccentricityPostCEE() const                 { return m_EccentricityPostCEE; }
    double              EccentricityPre2ndSN() const                { return m_EccentricityPre2ndSN; }
    double              EccentricityPreCEE() const                  { return m_EccentricityPreCEE; }
    double              EccentricityPrime() const                   { return m_EccentricityPrime; }
    ERROR               Error() const                               { return m_Error; }
    bool                HasOneOf(STELLAR_TYPE_LIST p_List) const;
    bool                HasStarsTouching() const                    { return (utils::Compare(m_SemiMajorAxisPrime, 0.0) > 0) && (m_SemiMajorAxisPrime <= RSOL_TO_AU * (m_Star1->Radius() + m_Star2->Radius())); }
    bool                HasTwoOf(STELLAR_TYPE_LIST p_List) const;
    bool                ImmediateRLOFPostCEE() const                { return m_ImmediateRLOFPostCEE; }
    STELLAR_TYPE        InitialStellarType1() const                 { return m_Star1->InitialStellarType(); }
    STELLAR_TYPE        InitialStellarType2() const                 { return m_Star2->InitialStellarType(); }
    bool                IsBeBinary() const                          { return HasOneOf({STELLAR_TYPE::NEUTRON_STAR}) && HasOneOf({STELLAR_TYPE::MS_LTE_07, STELLAR_TYPE::MS_GT_07}); }
    bool                IsBHandBH() const                           { return HasTwoOf({STELLAR_TYPE::BLACK_HOLE}); }
    bool                IsDCO() const                               { return HasTwoOf({STELLAR_TYPE::NEUTRON_STAR, STELLAR_TYPE::BLACK_HOLE}); }
    bool                IsGravitationallyBound() const              { return (utils::Compare(-(m_TotalOrbitalEnergyPrime / m_TotalOrbitalEnergyPrev), 0.0) <= 0.0); }                 // epsilon = -EPrime / E;
    bool                IsNSandBH() const                           { return HasOneOf({STELLAR_TYPE::NEUTRON_STAR}) && HasOneOf({STELLAR_TYPE::BLACK_HOLE}); }
    bool                IsNSandNS() const                           { return HasTwoOf({STELLAR_TYPE::NEUTRON_STAR}); }
    bool                IsUnbound() const                           { return (utils::Compare(m_SemiMajorAxisPrime, 0.0) <= 0 || (utils::Compare(m_Eccentricity, 1.0) > 0)); }         // semi major axis <= 0.0 means unbound, presumably by SN)
    bool                IsWDandWD() const                           { return HasTwoOf({STELLAR_TYPE::HELIUM_WHITE_DWARF, STELLAR_TYPE::CARBON_OXYGEN_WHITE_DWARF, STELLAR_TYPE::OXYGEN_NEON_WHITE_DWARF}); }
    double              LBV_Factor() const                          { return m_LBVfactor; }
	double              Mass1Final() const                          { return m_Mass1Final; }
    double              Mass2Final() const                          { return m_Mass2Final; }
    double              Mass1PreCEE() const                         { return m_Mass1PreCEE; }
    double              Mass2PreCEE() const                         { return m_Mass2PreCEE; }
    double              MassEnv1() const                            { return m_MassEnv1; }
    double              MassEnv2() const                            { return m_MassEnv2; }
    double              MassesEquilibrated() const                  { return m_MassesEquilibrated; }
    MT_TRACKING         MassTransferTrackerHistory() const          { return m_MassTransferTrackerHistory; }
    bool                MergesInHubbleTime() const                  { return m_MergesInHubbleTime; }
    bool                OptimisticCommonEnvelope() const            { return m_OptimisticCommonEnvelope; }
    double              OrbitalVelocity() const                     { return m_OrbitalVelocity; }
    double              OrbitalVelocityPre2ndSN() const             { return m_OrbitalVelocityPre2ndSN; }
	double              Radius1PostStripping() const                { return m_Radius1PostStripping; }
	double              Radius2PostStripping() const                { return m_Radius2PostStripping; }
	double              Radius1PreCEE() const                       { return m_Radius1PreCEE; }
	double              Radius2PreCEE() const                       { return m_Radius2PreCEE; }
	unsigned long int   RandomSeed() const                          { return m_RandomSeed; }
	BinaryRLOFDetailsT  RLOFDetails() const                         { return m_RLOFDetails; }
    bool                RLOFSecondaryPostCEE() const                { return m_Star2->RLOFPostCEE(); }
    double              RocheLobe1to2PostCEE() const                { return m_RocheLobe1to2PostCEE; }
    double              RocheLobe1to2PreCEE() const                 { return m_RocheLobe1to2PreCEE; }
    double              RocheLobe2to1PostCEE() const                { return m_RocheLobe2to1PostCEE; }
    double              RocheLobe2to1PreCEE() const                 { return m_RocheLobe2to1PreCEE; }
    double              RocheLobeRadius1() const                    { return m_Star1->RocheLobeRadius(); }
    double              RocheLobeRadius2() const                    { return m_Star2->RocheLobeRadius(); }
    double              RocheLobeTracker1() const                   { return m_Star1->RocheLobeTracker(); }
    double              RocheLobeTracker2() const                   { return m_Star2->RocheLobeTracker(); }
    double              SemiMajorAxisAtDCOFormation() const         { return m_SemiMajorAxisAtDCOFormation; }
    double              SemiMajorAxisInitial() const                { return m_SemiMajorAxisInitial; }
    double              SemiMajorAxisPostCEE() const                { return m_SemiMajorAxisPostCEE; }
    double              SemiMajorAxisPre2ndSN() const               { return m_SemiMajorAxisPre2ndSN; }
    double              SemiMajorAxisPreCEE() const                 { return m_SemiMajorAxisPreCEE; }
    double              SemiMajorAxisPrime() const                  { return m_SemiMajorAxisPrime; }
    bool                SimultaneousRLOF() const                    { return m_SimultaneousRLOF; }
    bool                StableRLOFPostCEE() const                   { return m_RLOFDetails.stableRLOFPostCEE; }
    bool                StellarMerger() const                       { return m_StellarMerger; }
    STELLAR_TYPE        StellarType1() const                        { return m_Star1->StellarType(); }
    STELLAR_TYPE        StellarType1PreCEE() const                  { return m_StellarType1PreCEE; }
    STELLAR_TYPE        StellarType2() const                        { return m_Star2->StellarType(); }
    STELLAR_TYPE        StellarType2PreCEE() const                  { return m_StellarType2PreCEE; }
    SN_STATE            SN_State() const                            { return m_SupernovaState; }
    bool                SurvivedSNEvent() const                     { return m_Survived; }
	double              SynchronizationTimescale() const            { return m_SynchronizationTimescale; }
    double              SystemicVelocity() const                    { return m_SystemicVelocity; }
    double              Time() const                                { return m_Time; }
    double              TimeToCoalescence() const                   { return m_TimeToCoalescence; }
    double              TotalAngularMomentumPrime() const           { return m_TotalAngularMomentumPrime; }
    double              TotalEnergyPrime() const                    { return m_TotalEnergyPrime; }
    double              UK() const                                  { return m_uK; }
    double              WolfRayetFactor() const                     { return m_WolfRayetFactor; }
    double              ZetaRLOFAnalytic() const                    { return m_ZetaRLOFAnalytic; }
    double              ZetaRLOFNumerical() const                   { return m_ZetaRLOFNumerical; }
    double              ZetaStarCompare() const                     { return m_ZetaStarCompare; }


    // member functions - alphabetically
            COMPAS_VARIABLE     BinaryPropertyValue(const T_ANY_PROPERTY p_Property) const;

    static  double              CalculateRocheLobeRadius_Static(const double p_MassPrimary, const double p_MassSecondary);

            EVOLUTION_STATUS    Evolve(const int p_Index);

            COMPAS_VARIABLE     PropertyValue(const T_ANY_PROPERTY p_Property) const;


private:

    BaseBinaryStar() { }

    OBJECT_ID    m_ObjectId;                                                                // Instantiated object's unique object id
    OBJECT_TYPE  m_ObjectType;                                                              // Instantiated object's object type
    STELLAR_TYPE m_StellarType;                                                             // Stellar type defined in Hurley et al. 2000

    ERROR m_Error;                                                                          // Records most recent error encountered for this binary

    // member variables - alphabetical in groups (sort of...)

    unsigned long int   m_RandomSeed;                                                       // Random seed for this binary

    AIS                 m_AIS;

    BeBinaryDetailsT    m_BeBinaryDetails;                                                  // BeBinary details

    double              m_Beta;                                                             // Angle between r and v, related to eccentricity (= pi/2 for circular e = 0)

	double              m_CircularizationTimescale;

    bool                m_CommonEnvelope;
    double              m_CommonEnvelopeAlpha;
    bool                m_CommonEnvelopeOccuredAtLeastOnce;
    unsigned int        m_CounterCEE;                                                       // Common Envelope Event count
    bool                m_DoubleCoreCommonEnvelope;
    bool                m_OptimisticCommonEnvelope;

    bool                m_Disbound;                                                         // Binary disbound?

    double              m_Dt;                                                               // Timestep

    double              m_EPrime;

    double              m_Eccentricity;                                                     // Initial eccentricity
    double              m_EccentricityAtDCOFormation;                                       // Eccentricity at DCO formation
    double              m_EccentricityInitial;                                              // Record initial eccentricity              JR: todo: check necessary
	double              m_EccentricityPostCEE;                                              // Eccentricity post common envelope event
    double              m_EccentricityPre2ndSN;                                             // Eccentricity prior to 2nd supernova
	double              m_EccentricityPreCEE;                                               // Eccentricity pre common envelope event
    double              m_EccentricityPrev;                                                 // Eccentricity at previous timestep
    double              m_EccentricityPrime;                                                // Initial eccentricity

    bool 	            m_FastPhaseCaseA;                                                   // Indicates if the system just entered a case A MT for the first time

    double	            m_FractionAccreted;	                                                // Fraction of mass accreted from the donor during mass transfer

    double              m_CosIPrime;
    double              m_IPrime;

   	double	            m_JLoss;			                                                // Specific angular momentum with which mass is lost during non-conservative mass transfer

    double              m_LBVfactor;

    double              m_MassesEquilibrated;                                               // Indicates whether CH stars had masses equilbrated at birth

	double              m_Mass1Final;                                                       // Star1 mass in Msol after losing its envelope (in this case, we asume it loses all of its envelope)
    double              m_Mass2Final;                                                       // Star2 mass in Msol after losing its envelope (in this case, we asume it loses all of its envelope)

    double              m_Mass1PreCEE;                                                      // Star1 mass in Msol before CEE
    double              m_Mass2PreCEE;                                                      // Star2 mass in Msol before CEE

    double              m_MassEnv1;                                                         // Star1 envelope mass in Msol
    double              m_MassEnv2;                                                         // Star2 envelope mass in Msol

    double              m_aMassLossDiff;
    double              m_OmegaMassLossDiff;

    bool                m_MassTransfer;
    double              m_aMassTransferDiff;
    double              m_OmegaMassTransferDiff;

    MT_TRACKING         m_MassTransferTrackerHistory;

	double              m_ReducedMassPrev;
	double              m_ReducedMassPrime;

	double              m_TotalMassPrev;
	double              m_TotalMassPrime;

    double              m_MC;
    double              m_MCPrime;

    bool                m_Merged;                                                           // Indicates if the stars merged
    bool                m_MergesInHubbleTime;                                               // Indicates if the stars merge in Hubble Time

    double              m_MSN;
    double              m_MSNPrime;

    double              m_OrbitalVelocity;
    double              m_OrbitalVelocityPre2ndSN;
    double              m_OrbitalVelocityPrev;
    double              m_OrbitalVelocityPrime;

    double              m_Radius;

	double              m_Radius1PostStripping;							                    // Star1 radius in Rsol after being stripped (= pre-CEE radius if not stripped)
	double              m_Radius2PostStripping;							                    // Star2 fadius in Rsol after being stripped (= pre-CEE radius if not stripped)

	double              m_Radius1PreCEE;									                // Star1 radius in Rsol before CEE
	double              m_Radius2PreCEE;									                // Star2 radius in Rsol before CEE

    bool                m_ImmediateRLOFPostCEE;
    BinaryRLOFDetailsT  m_RLOFDetails;                                                      // RLOF details
    bool                m_SimultaneousRLOF;

	double              m_RocheLobe1to2PreCEE;
	double              m_RocheLobe1to2PostCEE;
	double              m_RocheLobe2to1PreCEE;
	double              m_RocheLobe2to1PostCEE;

    double              m_SemiMajorAxis;                                                    // Semi-major axis
    double              m_SemiMajorAxisAtDCOFormation;                                      // Semi-major axis at DCO formation
    double              m_SemiMajorAxisInitial;                                             // Record initial semi-major axis              JR: todo: check necessary
	double              m_SemiMajorAxisPostCEE;                                             // Semi-major axis post common envelope event
    double              m_SemiMajorAxisPre2ndSN;                                            // Semi-major axis prior to 2nd supernova
  	double              m_SemiMajorAxisPreCEE;                                              // Semi-major axis pre common envelope event
    double              m_SemiMajorAxisPrev;                                                // Semi-major axis at previous timestep
    double              m_SemiMajorAxisPrime;                                               // Semi-major axis

    bool                m_StellarMerger;

    STELLAR_TYPE        m_StellarType1PreCEE;                                               // Stellar type of star1 prior to the Common Envelope Event
    STELLAR_TYPE        m_StellarType2PreCEE;                                               // Stellar type of star2 prior to the Common Envelope Event

    SN_STATE            m_SupernovaState;                                                   // Indicates which star (or stars) are undergoing / hove undergone a supernova event

    bool                m_Survived;

	double              m_SynchronizationTimescale;

    double              m_SystemicVelocity;                                                 // Post supernova systemic velocity

    double              m_Time;                                                             // Physical simulation time
    double              m_TimePrev;                                                         // Previous simulation time
    double              m_TimeToCoalescence;                                                // Coalescence time

    double              m_TotalAngularMomentumPrev;
    double              m_TotalAngularMomentumPrime;

    double              m_TotalEnergy;
    double              m_TotalEnergyPrime;

	double              m_TotalOrbitalAngularMomentumPrev;
	double              m_TotalOrbitalAngularMomentumPrime;

	double              m_TotalOrbitalEnergyPrev;
	double              m_TotalOrbitalEnergyPrime;

    double              m_uK;

    double              m_VRel;

    double              m_WolfRayetFactor;

    double              m_ZetaRLOFAnalytic;
    double              m_ZetaRLOFNumerical;
    double              m_ZetaStarCompare;


    // Binaries contain two stars
    BinaryConstituentStar *m_Star1;                                                         // Initially more massive star - the primary
    BinaryConstituentStar *m_Star2;                                                         // Initially less massive star - the secondary

    BinaryConstituentStar *m_Donor;                                                         // Pointer to the donor for mass transfer
    BinaryConstituentStar *m_Accretor;                                                      // Pointer to the accretor for mass transfer

    BinaryConstituentStar *m_Supernova;                                                     // Pointer to the star that is undergoing / has undergone a supernova event
    BinaryConstituentStar *m_Companion;                                                     // Pointer to the companion star to the supernova


    // member functions - alphabetical in groups (sort of...)

    // CalculateAngularMomentum - the actual function takes 10 parameters because of the various calling permutations
    //                          - various signatures are defined here - they just assemble the parameters as required
    //                            and call the actual function
    // JR: todo: note in the orginal code the binary orbital velicity was passed in as a parameter but never used - I removed it

    void    SetInitialCommonValues(const AIS &p_AIS, const long int p_Id);
    void    SetRemainingCommonValues(const long int p_Id);


    double  CalculateAngularMomentum(const double p_SemiMajorAxis,
                                     const double p_Eccentricity,
                                     const double p_Star1Mass,
                                     const double p_Star2Mass,
                                     const double p_Star1Radius,
                                     const double p_Star2Radius,
                                     const double p_Star1OrbitalFrequency,
                                     const double p_Star2OrbitalFrequency,
                                     const double p_Star1GyrationRadius,
                                     const double p_Star2GyrationRadius);

    double  CalculateAngularMomentum()                                      { return CalculateAngularMomentum(m_SemiMajorAxisPrime, m_EccentricityPrime, m_Star1->Mass(), m_Star2->Mass(), m_Star1->Radius(), m_Star2->Radius(), m_Star1->Omega(), m_Star2->Omega(), m_Star1->CalculateGyrationRadius(), m_Star2->CalculateGyrationRadius()); }

    double  CalculateAngularMomentum(const double p_SemiMajorAxis,
                                     const double p_Eccentricity,
                                     const double p_Star1_OrbitalFrequency,
                                     const double p_Star2_OrbitalFrequency,
                                     const double p_Star1_GyrationRadius,
                                     const double p_Star2_GyrationRadius)   { return CalculateAngularMomentum(p_SemiMajorAxis, p_Eccentricity, m_Star1->Mass(), m_Star2->Mass(), m_Star1->Radius(), m_Star2->Radius(), p_Star1_OrbitalFrequency, p_Star2_OrbitalFrequency, p_Star1_GyrationRadius, p_Star2_GyrationRadius); }

    double  CalculateAngularMomentumPrev()                                  { return CalculateAngularMomentum(m_SemiMajorAxisPrev, m_EccentricityPrev, m_Star1->MassPrev(), m_Star2->MassPrev(), m_Star1->RadiusPrev(), m_Star2->RadiusPrev(), m_Star1->OmegaPrev(), m_Star2->OmegaPrev(), m_Star1->CalculateGyrationRadius(), m_Star2->CalculateGyrationRadius()); }

    double  CalculateAngularMomentumPrime()                                 { return CalculateAngularMomentum(m_SemiMajorAxisPrime, m_EccentricityPrime, m_Star1->Mass(), m_Star2->Mass(), m_Star1->Radius(), m_Star2->Radius(), m_Star1->Omega(), m_Star2->Omega(), m_Star1->CalculateGyrationRadius(), m_Star2->CalculateGyrationRadius()); }


    void    CalculateCommonEnvelopeValues();

    double  CalculateCosFinalPlaneTilt(const double p_KickTheta, const double p_KickPhi);

    void    CalculateEnergyAndAngularMomentum();

    double  CalculateGammaAngularMomentumLoss(const double p_DonorMass, const double p_AccretorMass);
    double  CalculateGammaAngularMomentumLoss()                             { return CalculateGammaAngularMomentumLoss(m_Donor->Mass(), m_Accretor->Mass()); }


    void    CalculateMassTransfer(const double p_Dt);
    double  CalculateMassTransferFastPhaseCaseA(const double p_JLoss);
    double  CalculateMassTransferOrbit(BinaryConstituentStar& p_Donor, BinaryConstituentStar& p_Accretor, const double p_MDotDonor, const double p_Dt, const double p_Jloss);
    void    CalculateWindsMassLoss();
    void    CheckMassTransfer(const double p_Dt);
    void    InitialiseMassTransfer();

    double  CalculateOrbitalAngularMomentum(const double p_Mu,
                                            const double p_Mass,
                                            const double p_SemiMajorAxis)   { return p_Mu * sqrt(G1 * p_Mass * p_SemiMajorAxis); }

    double  CalculateOrbitalEccentricityPostSupernova(const double p_KickVelocity,
                                                      const double p_TotalMassPreSN,
                                                      const double p_TotalMassPostSN,
                                                      const double p_KickTheta,
                                                      const double p_KickPhi);

    double  CalculateOrbitalEnergy(const double p_Mu,
                                   const double p_Mass,
                                   const double p_SemiMajorAxis)            { return -(G1 * p_Mu * p_Mass) / (2.0 * p_SemiMajorAxis); }

    double  CalculatePostSNSystemicVelocity(const double p_SNMass,
                                            const double p_SNDeltaMass,
                                            const double p_CompanionMass,
                                            const double p_TotalMassPreSN,
                                            const double p_TotalMassPostSN,
                                            const double p_KickTheta,
                                            const double p_KickPhi);

    double  CalculateAdaptiveRocheLobeOverFlow(const double p_JLoss);
    double  CalculateNumericalZRocheLobe(const double p_jLoss);
    double  CalculateZRocheLobe();

    double  CalculateSemiMajorAxisPostSupernova(const double p_KickVelocity,
                                                const double p_TotalMassPreSN,
                                                const double p_TotalMassPostSN,
                                                const double p_KickTheta,
                                                const double p_KickPhi);

    double  CalculateDt(const double p_Dt, const Star* const p_Primary, const Star* const p_Secondary);
    double  CalculateTimestep(const double p_Dt);
    double  ChooseTimestep(const double p_Dt);

    double  CalculateTimeToCoalescence(double a0, double e0, double m1, double m2);

    // CalculateTotalEnergy - the actual function takes 10 parameters because of the various calling permutations
    //                      - various signatures are defined here - they just assemble the parameters as required
    //                        and call the actual function
    double  CalculateTotalEnergy(const double p_SemiMajorAxis,
                                 const double p_Star1Mass,
                                 const double p_Star2Mass,
                                 const double p_Star1Radius,
                                 const double p_Star2Radius,
                                 const double p_Star1_OrbitalFrequency,
                                 const double p_Star2_OrbitalFrequency,
                                 const double p_OrbitalVelocity,
                                 const double p_Star1GyrationRadius,
                                 const double p_Star2GyrationRadius);

    double  CalculateTotalEnergy()                                          { return CalculateTotalEnergy(m_SemiMajorAxis, m_Star1->Mass(), m_Star2->Mass(), m_Star1->Radius(), m_Star2->Radius(), m_Star1->Omega(), m_Star2->Omega(), m_OrbitalVelocity, m_Star1->CalculateGyrationRadius(), m_Star2->CalculateGyrationRadius()); }

    double  CalculateTotalEnergy(const double p_SemiMajorAxis,
                                 const double p_Star1_OrbitalFrequency,
                                 const double p_Star2_OrbitalFrequency,
                                 const double p_OrbitalVelocity,
                                 const double p_Star1_GyrationRadius,
                                 const double p_Star2_GyrationRadius)       { return CalculateTotalEnergy(p_SemiMajorAxis, m_Star1->Mass(), m_Star2->Mass(), m_Star1->Radius(), m_Star2->Radius(), p_Star1_OrbitalFrequency, p_Star2_OrbitalFrequency, p_OrbitalVelocity, p_Star1_GyrationRadius, p_Star2_GyrationRadius); }

    double  CalculateTotalEnergyPrime()                                     { return CalculateTotalEnergy(m_SemiMajorAxisPrime, m_Star1->Mass(), m_Star2->Mass(), m_Star1->Radius(), m_Star2->Radius(), m_Star1->Omega(), m_Star2->Omega(), m_OrbitalVelocityPrime, m_Star1->CalculateGyrationRadius(), m_Star2->CalculateGyrationRadius()); }

    void    EvaluateBinary(const double p_Dt);
    void    EvaluateBinaryPreamble();

    void    EvaluateSupernovae();

    void    EvolveOneTimestep(const double p_Dt);
    void    EvolveOneTimestepPreamble(const double p_Dt);

    void    ResolveCoalescence();
    void    ResolveCommonEnvelopeEvent();
    void    ResolveMassChanges();
    bool    ResolveSupernova();

    bool    RLOFCheckMonitorEndCaseA();

    double  SampleSemiMajorAxisDistribution(const double p_Mass1, const double p_Mass2);
    double  SampleEccentricityDistribution();
    double  SampleInitialMassDistribution();
    double  SampleMetallicityDistribution();
    double  SampleQDistribution();

    bool    ShouldPrintRLOFParameters();

    void    StashBeBinaryProperties();
    void    StashCurrentRLOFProperties();
    void    StashPreviousRLOFProperties();
    void    StashRLOFProperties(RLOFProperties &p_RLOFProperties);


    // printing functions
    void PrintRLOFParameters();
    void PrintBinarySystemParameters()          {                                   LOGGING->LogBinarySystemParameters(this); }
    void PrintDetailedOutput(const int p_Id)    { if (OPTIONS->DetailedOutput())    LOGGING->LogDetailedOutput(this, p_Id); }
    void PrintDoubleCompactObjects()            {                                   LOGGING->LogDoubleCompactObject(this); }
    void PrintCommonEnvelope()                  {                                   LOGGING->LogCommonEnvelope(this); }
    void PrintBeBinary()                        { if (OPTIONS->BeBinaries())        LOGGING->LogBeBinary(this); }
    void PrintPulsarEvolutionParameters()       { if (OPTIONS->EvolvePulsars())     LOGGING->LogPulsarEvolutionParameters(this); }
    void PrintSupernovaDetails()                {                                   LOGGING->LogSupernovaDetails(this); }

};

#endif // __BaseBinaryStar_h__
