#ifndef __BaseBinaryStar_h__
#define __BaseBinaryStar_h__

#include "constants.h"
#include "typedefs.h"
#include "profiling.h"
#include "utils.h"
#include "vector3d.h"

#include "Log.h"
#include "Star.h"
#include "BinaryConstituentStar.h"

#include <boost/math/tools/roots.hpp>

#include <boost/math/special_functions/next.hpp>    // For float_distance.
#include <boost/math/special_functions/cbrt.hpp>    // For boost::math::cbrt.

#include <tuple>                                    // for std::tuple and std::make_tuple.


class Log;
class Star;
class BinaryConstituentStar;


class BaseBinaryStar {

public:

    BaseBinaryStar(const unsigned long int p_Seed, const long int p_Id);

    void CopyMemberVariables(const BaseBinaryStar& p_Star) {

        m_Id                               = p_Star.m_Id;

        m_Error                            = p_Star.m_Error;

        m_RandomSeed                       = p_Star.m_RandomSeed;

        m_BeBinaryDetails                  = p_Star.m_BeBinaryDetails;

        m_BeBinaryDetails.currentProps     = p_Star.m_BeBinaryDetails.currentProps  == &(p_Star.m_BeBinaryDetails.props1) ? &(m_BeBinaryDetails.props1) : &(m_BeBinaryDetails.props2);
        m_BeBinaryDetails.previousProps    = p_Star.m_BeBinaryDetails.previousProps == &(p_Star.m_BeBinaryDetails.props1) ? &(m_BeBinaryDetails.props1) : &(m_BeBinaryDetails.props2);

        m_CircularizationTimescale         = p_Star.m_CircularizationTimescale;

        m_CEDetails                        = p_Star.m_CEDetails;

        m_Unbound                          = p_Star.m_Unbound;

        m_DCOFormationTime                 = p_Star.m_DCOFormationTime;
        
        m_Dt                               = p_Star.m_Dt;

        m_Eccentricity                     = p_Star.m_Eccentricity;
        m_EccentricityAtDCOFormation       = p_Star.m_EccentricityAtDCOFormation;
        m_EccentricityInitial              = p_Star.m_EccentricityInitial;
        m_EccentricityPreSN                = p_Star.m_EccentricityPreSN;
        m_EccentricityPrev                 = p_Star.m_EccentricityPrev;

        m_Flags                            = p_Star.m_Flags;
        
        m_FractionAccreted                 = p_Star.m_FractionAccreted;

        m_CosIPrime                        = p_Star.m_CosIPrime;
        m_IPrime                           = p_Star.m_IPrime;

        m_JLoss                            = p_Star.m_JLoss;

        m_Mass1Final                       = p_Star.m_Mass1Final;
        m_Mass2Final                       = p_Star.m_Mass2Final;

        m_MassEnv1                         = p_Star.m_MassEnv1;
        m_MassEnv2                         = p_Star.m_MassEnv2;

        m_aMassLossDiff                    = p_Star.m_aMassLossDiff;

        m_MassTransfer                     = p_Star.m_MassTransfer;
        m_aMassTransferDiff                = p_Star.m_aMassTransferDiff;

        m_MassTransferTrackerHistory       = p_Star.m_MassTransferTrackerHistory;

        m_OrbitalVelocityPreSN             = p_Star.m_OrbitalVelocityPreSN;

        m_PrintExtraDetailedOutput         = p_Star.m_PrintExtraDetailedOutput;

        m_RLOFDetails                      = p_Star.m_RLOFDetails;
        m_RLOFDetails.propsPreMT           = p_Star.m_RLOFDetails.propsPreMT == &(p_Star.m_RLOFDetails.props1) ? &(m_RLOFDetails.props1) : &(m_RLOFDetails.props2);
        m_RLOFDetails.propsPostMT          = p_Star.m_RLOFDetails.propsPostMT == &(p_Star.m_RLOFDetails.props1) ? &(m_RLOFDetails.props1) : &(m_RLOFDetails.props2);

        m_SemiMajorAxis                    = p_Star.m_SemiMajorAxis;
        m_SemiMajorAxisAtDCOFormation      = p_Star.m_SemiMajorAxisAtDCOFormation;
        m_SemiMajorAxisInitial             = p_Star.m_SemiMajorAxisInitial;
        m_SemiMajorAxisPreSN               = p_Star.m_SemiMajorAxisPreSN;
        m_SemiMajorAxisPrev                = p_Star.m_SemiMajorAxisPrev;

        m_SupernovaState                   = p_Star.m_SupernovaState;

        m_SynchronizationTimescale         = p_Star.m_SynchronizationTimescale;

        m_SystemicVelocity                 = p_Star.m_SystemicVelocity;

        m_ThetaE                           = p_Star.m_ThetaE;
        m_PhiE                             = p_Star.m_PhiE;  
        m_PsiE                             = p_Star.m_PsiE;  

        m_Time                             = p_Star.m_Time;
        m_TimePrev                         = p_Star.m_TimePrev;
        m_TimeToCoalescence                = p_Star.m_TimeToCoalescence;

        m_TotalAngularMomentum             = p_Star.m_TotalAngularMomentum;

        m_TotalEnergy                      = p_Star.m_TotalEnergy;

        m_OrbitalAngularMomentumPrev       = p_Star.m_OrbitalAngularMomentumPrev;
        m_OrbitalAngularMomentum           = p_Star.m_OrbitalAngularMomentum;

        m_OrbitalEnergyPrev                = p_Star.m_OrbitalEnergyPrev;
        m_OrbitalEnergy                    = p_Star.m_OrbitalEnergy;

        m_uK                               = p_Star.m_uK;

        m_ZetaLobe                         = p_Star.m_ZetaLobe;
        m_ZetaStar                         = p_Star.m_ZetaStar;

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

            m_ObjectId    = globalObjectId++;                       // get unique object id (don't copy source)
            m_ObjectType  = OBJECT_TYPE::BASE_BINARY_STAR;          // can only copy from BASE_BINARY_STAR
            m_StellarType = STELLAR_TYPE::BINARY_STAR;              // always

            CopyMemberVariables(p_Star);                            // copy member variables
        }
        return *this;
    }


    virtual ~BaseBinaryStar() { delete m_Star1; delete m_Star2; }


    // object identifiers - all classes have these
    OBJECT_ID           ObjectId() const                            { return m_ObjectId; }
    OBJECT_TYPE         ObjectType() const                          { return m_ObjectType; }
    STELLAR_TYPE        StellarType() const                         { return m_StellarType; }
    long int            Id() const                                  { return m_Id; }


    // getters - alphabetically
    BeBinaryDetailsT    BeBinaryDetails() const                     { return m_BeBinaryDetails; }
    bool                CEAtLeastOnce() const                       { return m_CEDetails.CEEcount > 0; }
    unsigned int        CEEventCount() const                        { return m_CEDetails.CEEcount; }
    double              CircularizationTimescale() const            { return m_CircularizationTimescale; }
    unsigned int        CommonEnvelopeEventCount() const            { return m_CEDetails.CEEcount; }
    bool                Unbound() const                             { return m_Unbound; }
    bool                DoubleCoreCE() const                        { return m_CEDetails.doubleCoreCE; }
    double              Dt() const                                  { return m_Dt; }
    double              Eccentricity() const                        { return m_Eccentricity; }
    double              EccentricityAtDCOFormation() const          { return m_EccentricityAtDCOFormation; }
    double              EccentricityInitial() const                 { return m_EccentricityInitial; }
    double              EccentricityPostCEE() const                 { return m_CEDetails.postCEE.eccentricity; }
    double              EccentricityPreSN() const                   { return m_EccentricityPreSN; }
    double              EccentricityPreCEE() const                  { return m_CEDetails.preCEE.eccentricity; }
    ERROR               Error() const                               { return m_Error; }
    double              FractionAccreted() const                    { return m_FractionAccreted; }
    bool                HasOnlyOneOf(STELLAR_TYPE_LIST p_List) const;
    bool                HasOneOf(STELLAR_TYPE_LIST p_List) const;
    bool                HasStarsTouching() const                    { return (utils::Compare(m_SemiMajorAxis, 0.0) > 0) && (m_SemiMajorAxis <= RSOL_TO_AU * (m_Star1->Radius() + m_Star2->Radius())); }
    bool                HasTwoOf(STELLAR_TYPE_LIST p_List) const;
    bool                ImmediateRLOFPostCEE() const                { return m_RLOFDetails.immediateRLOFPostCEE; }
    STELLAR_TYPE        InitialStellarType1() const                 { return m_Star1->InitialStellarType(); }
    STELLAR_TYPE        InitialStellarType2() const                 { return m_Star2->InitialStellarType(); }
    bool                IsBeBinary() const                          { return HasOneOf({STELLAR_TYPE::NEUTRON_STAR}) && HasOneOf({STELLAR_TYPE::MS_LTE_07, STELLAR_TYPE::MS_GT_07}); }
    bool                IsHMXRBinary() const;
    bool                IsBHandBH() const                           { return HasTwoOf({STELLAR_TYPE::BLACK_HOLE}); }
    bool                IsDCO() const                               { return HasTwoOf({STELLAR_TYPE::NEUTRON_STAR, STELLAR_TYPE::BLACK_HOLE}); }
    bool                IsNSandBH() const                           { return HasOneOf({STELLAR_TYPE::NEUTRON_STAR}) && HasOneOf({STELLAR_TYPE::BLACK_HOLE}); }
    bool                IsNSandNS() const                           { return HasTwoOf({STELLAR_TYPE::NEUTRON_STAR}); }
    bool                IsUnbound() const                           { return (utils::Compare(m_SemiMajorAxis, 0.0) <= 0 || (utils::Compare(m_Eccentricity, 1.0) > 0)); }         // semi major axis <= 0.0 means unbound, presumably by SN)
    bool                IsWDandWD() const                           { return HasTwoOf({STELLAR_TYPE::HELIUM_WHITE_DWARF, STELLAR_TYPE::CARBON_OXYGEN_WHITE_DWARF, STELLAR_TYPE::OXYGEN_NEON_WHITE_DWARF}); }
    double              Mass1PostCEE() const                        { return m_Star1->MassPostCEE(); }
    double              Mass1PreCEE() const                         { return m_Star1->MassPreCEE(); }
    double              Mass2PostCEE() const                        { return m_Star2->MassPostCEE(); }
    double              Mass2PreCEE() const                         { return m_Star2->MassPreCEE(); }
    double              MassEnv1() const                            { return m_MassEnv1; }
    double              MassEnv2() const                            { return m_MassEnv2; }
    bool                MassesEquilibrated() const                  { return m_Flags.massesEquilibrated; }
    bool                MassesEquilibratedAtBirth() const           { return m_Flags.massesEquilibratedAtBirth; }
    MT_TRACKING         MassTransferTrackerHistory() const          { return m_MassTransferTrackerHistory; }
    bool                MergesInHubbleTime() const                  { return m_Flags.mergesInHubbleTime; }
    bool                OptimisticCommonEnvelope() const            { return m_CEDetails.optimisticCE; }
    double              OrbitalAngularVelocity() const              { return std::sqrt(G1 * (m_Star1->Mass() + m_Star2->Mass()) / (m_SemiMajorAxis * m_SemiMajorAxis * m_SemiMajorAxis)); }      // rads/year
    double              OrbitalVelocityPreSN() const                { return m_OrbitalVelocityPreSN; }
    double              Periastron() const                          { return m_SemiMajorAxis * (1.0-m_Eccentricity); }
    double              PeriastronRsol() const                      { return Periastron() * AU_TO_RSOL; }
    double              Radius1PostCEE() const                      { return m_Star1->RadiusPostCEE(); }
    double              Radius2PostCEE() const                      { return m_Star2->RadiusPostCEE(); }
    double              Radius1PreCEE() const                       { return m_Star1->RadiusPreCEE(); }
    double              Radius2PreCEE() const                       { return m_Star2->RadiusPreCEE(); }
    unsigned long int   RandomSeed() const                          { return m_RandomSeed; }
    BinaryRLOFDetailsT  RLOFDetails() const                         { return m_RLOFDetails; }
    bool                RLOFSecondaryPostCEE() const                { return m_Star2->RLOFPostCEE(); }
    double              RocheLobe1to2PostCEE() const                { return m_CEDetails.postCEE.rocheLobe1to2; }
    double              RocheLobe1to2PreCEE() const                 { return m_CEDetails.preCEE.rocheLobe1to2; }
    double              RocheLobe2to1PostCEE() const                { return m_CEDetails.postCEE.rocheLobe2to1; }
    double              RocheLobe2to1PreCEE() const                 { return m_CEDetails.preCEE.rocheLobe2to1; }
    double              RocheLobeRadius1() const                    { return CalculateRocheLobeRadius_Static(m_Star1->Mass(), m_Star2->Mass()); }
    double              RocheLobeRadius2() const                    { return CalculateRocheLobeRadius_Static(m_Star2->Mass(), m_Star1->Mass()); }
    double              StarToRocheLobeRadiusRatio1() const         { return m_Star1->StarToRocheLobeRadiusRatio(m_SemiMajorAxis, m_Eccentricity); }
    double              StarToRocheLobeRadiusRatio2() const         { return m_Star2->StarToRocheLobeRadiusRatio(m_SemiMajorAxis, m_Eccentricity); }
    double              SemiMajorAxisAtDCOFormation() const         { return m_SemiMajorAxisAtDCOFormation; }
    double              SemiMajorAxisInitial() const                { return m_SemiMajorAxisInitial; }
    double              SemiMajorAxisPostCEE() const                { return m_CEDetails.postCEE.semiMajorAxis; }
    double              SemiMajorAxisPreSN() const                  { return m_SemiMajorAxisPreSN; }
    double              SemiMajorAxisPreCEE() const                 { return m_CEDetails.preCEE.semiMajorAxis; }
    double              SemiMajorAxis() const                       { return m_SemiMajorAxis; }
    double              SemiMajorAxisRsol() const                   { return m_SemiMajorAxis * AU_TO_RSOL; }
    bool                SimultaneousRLOF() const                    { return m_RLOFDetails.simultaneousRLOF; }
    bool                StableRLOFPostCEE() const                   { return m_RLOFDetails.stableRLOFPostCEE; }
    bool                StellarMerger() const                       { return m_Flags.stellarMerger; }
    bool                StellarMergerAtBirth() const                { return m_Flags.stellarMergerAtBirth; }
    STELLAR_TYPE        StellarType1() const                        { return m_Star1->StellarType(); }
    STELLAR_TYPE        StellarType1PostCEE() const                 { return m_Star1->StellarTypePostCEE(); }
    STELLAR_TYPE        StellarType1PreCEE() const                  { return m_Star1->StellarTypePreCEE(); }
    STELLAR_TYPE        StellarType2() const                        { return m_Star2->StellarType(); }
    STELLAR_TYPE        StellarType2PostCEE() const                 { return m_Star2->StellarTypePostCEE(); }
    STELLAR_TYPE        StellarType2PreCEE() const                  { return m_Star2->StellarTypePreCEE(); }
    SN_STATE            SN_State() const                            { return m_SupernovaState; }
    double              SynchronizationTimescale() const            { return m_SynchronizationTimescale; }
    double              SystemicSpeed() const                       { return m_SystemicVelocity.Magnitude(); }
    double              Time() const                                { return m_Time; }
    double              TimeToCoalescence() const                   { return m_TimeToCoalescence; }
    double              TotalAngularMomentum() const                { return m_TotalAngularMomentum; }
    double              TotalEnergy() const                         { return m_TotalEnergy; }
    double              UK() const                                  { return m_uK; }
    double              ZetaLobe() const                    	    { return m_ZetaLobe; }
    double              ZetaStar() const                            { return m_ZetaStar; }


    // member functions - alphabetically
            COMPAS_VARIABLE     BinaryPropertyValue(const T_ANY_PROPERTY p_Property) const;

    static  double              CalculateRocheLobeRadius_Static(const double p_MassPrimary, const double p_MassSecondary);

            EVOLUTION_STATUS    Evolve();

            bool                PrintSwitchLog(const bool p_PrimarySwitching) { return OPTIONS->SwitchLog() ? LOGGING->LogBSESwitchLog(this, p_PrimarySwitching) : true; }

            COMPAS_VARIABLE     PropertyValue(const T_ANY_PROPERTY p_Property) const;

            BinaryConstituentStar* Star1() { return m_Star1; }                              // Returns a pointer to the primary - here mainly to support the BSE Switch Log. Be careful!
            BinaryConstituentStar* Star2() { return m_Star2; }                              // Returns a pointer to the secondary - here mainly to support the BSE Switch Log. Be careful!

private:

    BaseBinaryStar() { }

    OBJECT_ID    m_ObjectId;                                                                // Instantiated object's unique object id
    OBJECT_TYPE  m_ObjectType;                                                              // Instantiated object's object type
    STELLAR_TYPE m_StellarType;                                                             // Stellar type defined in Hurley et al. 2000
    long int     m_Id;                                                                      // Id used to name detailed output file - uses p_Id as passed (usually the index number of multiple binaries being produced)

    ERROR m_Error;                                                                          // Records most recent error encountered for this binary

    // member variables - alphabetical in groups (sort of...)

    unsigned long int   m_RandomSeed;                                                       // Random seed for this binary

    BeBinaryDetailsT    m_BeBinaryDetails;                                                  // BeBinary details

    BinaryCEDetailsT    m_CEDetails;                                                        // Common Event details

    double              m_CircularizationTimescale;

    bool                m_Unbound;                                                          // Binary unbound?

    double              m_Dt;                                                               // Timestep

    double              m_Eccentricity;                                                     // Initial eccentricity
    double              m_EccentricityAtDCOFormation;                                       // Eccentricity at DCO formation
    double              m_EccentricityInitial;                                              // Record initial eccentricity              JR: todo: check necessary
    double              m_EccentricityPreSN;                                                // Eccentricity prior to supernova
    double              m_EccentricityPrev;                                                 // Eccentricity at previous timestep

    struct FLAGS {                                                                          // Miscellaneous flags

        bool massesEquilibrated;                                                            // Indicates whether stars had masses equilbrated at some stage after birth
        bool massesEquilibratedAtBirth;                                                     // Indicates whether stars had masses equilbrated at birth

        bool mergesInHubbleTime;                                                            // Indicates if the stars merge in Hubble Time

        bool stellarMerger;                                                                 // Indicates that the constituent stars merged
        bool stellarMergerAtBirth;                                                          // Indicates that the constituent stars were touching at bierth

    }                   m_Flags;

    double	            m_FractionAccreted;	                                                // Fraction of mass accreted from the donor during mass transfer

    double              m_CosIPrime;
    double              m_IPrime;

    double	            m_JLoss;			                                                // Specific angular momentum with which mass is lost during non-conservative mass transfer

    double              m_Mass1Final;                                                       // Star1 mass in Msol after losing its envelope (in this case, we asume it loses all of its envelope)
    double              m_Mass2Final;                                                       // Star2 mass in Msol after losing its envelope (in this case, we asume it loses all of its envelope)

    double              m_MassEnv1;                                                         // Star1 envelope mass in Msol
    double              m_MassEnv2;                                                         // Star2 envelope mass in Msol

    double              m_aMassLossDiff;

    bool                m_MassTransfer;
    double              m_aMassTransferDiff;

    MT_TRACKING         m_MassTransferTrackerHistory;

    double              m_OrbitalVelocityPreSN;

    bool                m_PrintExtraDetailedOutput;                                         // Flag to ensure that detailed output only gets printed once per timestep

    BinaryRLOFDetailsT  m_RLOFDetails;                                                      // RLOF details

    double              m_SemiMajorAxis;                                                    // Semi-major axis
    double              m_SemiMajorAxisAtDCOFormation;                                      // Semi-major axis at DCO formation
    double              m_SemiMajorAxisInitial;                                             // Record initial semi-major axis              JR: todo: check necessary
    double              m_SemiMajorAxisPreSN;                                               // Semi-major axis prior to supernova
    double              m_SemiMajorAxisPrev;                                                // Semi-major axis at previous timestep

    SN_STATE            m_SupernovaState;                                                   // Indicates which star (or stars) are undergoing / have undergone a supernova event

    double              m_SynchronizationTimescale;

    Vector3d            m_SystemicVelocity;                                                 // Systemic velocity vector, relative to ZAMS Center of Mass
    double              m_ThetaE;                                                           // Euler Theta
    double              m_PhiE;                                                             // Euler Phi                
    double              m_PsiE;                                                             // Euler Psi
    
    double              m_Time;                                                             // Physical simulation time
    double              m_TimePrev;                                                         // Previous simulation time
    double              m_TimeToCoalescence;                                                // Coalescence time
    double              m_DCOFormationTime;                                                 // Time of DCO formation

    double              m_TotalAngularMomentum;

    double              m_TotalEnergy;

	double              m_OrbitalAngularMomentumPrev;
	double              m_OrbitalAngularMomentum;

	double              m_OrbitalEnergyPrev;
	double              m_OrbitalEnergy;

    double              m_uK;

    double              m_ZetaLobe;
    double              m_ZetaStar;


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

    void    SetInitialValues(const unsigned long int p_Seed, const long int p_Id);
    void    SetRemainingValues();


    double  CalculateAngularMomentum(const double p_SemiMajorAxis,
                                     const double p_Eccentricity,
                                     const double p_Star1Mass,
                                     const double p_Star2Mass,
                                     const double p_Star1Radius,
                                     const double p_Star2Radius,
                                     const double p_Star1_SpinAngularVelocity,
                                     const double p_Star2_SpinAngularVelocity,
                                     const double p_Star1_GyrationRadius,
                                     const double p_Star2_GyrationRadius) const;

    double  CalculateAngularMomentum() const                                    { return CalculateAngularMomentum(m_SemiMajorAxis, m_Eccentricity, m_Star1->Mass(), m_Star2->Mass(), m_Star1->Radius(), m_Star2->Radius(), m_Star1->Omega(), m_Star2->Omega(), m_Star1->CalculateGyrationRadius(), m_Star2->CalculateGyrationRadius()); }

    double  CalculateAngularMomentum(const double p_SemiMajorAxis,
                                     const double p_Eccentricity,
                                     const double p_Star1_SpinAngularVelocity,
                                     const double p_Star2_SpinAngularVelocity,
                                     const double p_Star1_GyrationRadius,
                                     const double p_Star2_GyrationRadius) const { return CalculateAngularMomentum(p_SemiMajorAxis, p_Eccentricity, m_Star1->Mass(), m_Star2->Mass(), m_Star1->Radius(), m_Star2->Radius(), p_Star1_SpinAngularVelocity, p_Star2_SpinAngularVelocity, p_Star1_GyrationRadius, p_Star2_GyrationRadius); }

    double  CalculateAngularMomentumPrev() const                                { return CalculateAngularMomentum(m_SemiMajorAxisPrev, m_EccentricityPrev, m_Star1->MassPrev(), m_Star2->MassPrev(), m_Star1->RadiusPrev(), m_Star2->RadiusPrev(), m_Star1->OmegaPrev(), m_Star2->OmegaPrev(), m_Star1->CalculateGyrationRadius(), m_Star2->CalculateGyrationRadius()); }

    void    CalculateEnergyAndAngularMomentum();

    double  CalculateGammaAngularMomentumLoss(const double p_DonorMass, const double p_AccretorMass);
    double  CalculateGammaAngularMomentumLoss()                                 { return CalculateGammaAngularMomentumLoss(m_Donor->Mass(), m_Accretor->Mass()); }


    void    CalculateMassTransfer(const double p_Dt);

    double  CalculateMassTransferOrbit(const double                 p_DonorMass, 
                                       const double                 p_DeltaMassDonor, 
                                             BinaryConstituentStar& p_Accretor, 
                                       const double                 p_FractionAccreted);

    void    CalculateWindsMassLoss();
    void    InitialiseMassTransfer();

    double  CalculateOrbitalAngularMomentum(const double p_Mu,
                                            const double p_Mass,
                                            const double p_SemiMajorAxis) const { return p_Mu * std::sqrt(G1 * p_Mass * p_SemiMajorAxis); }

    double  CalculateOrbitalEnergy(const double p_Mu,
                                   const double p_Mass,
                                   const double p_SemiMajorAxis) const          { return -(G1 * p_Mu * p_Mass) / (2.0 * p_SemiMajorAxis); }

    double  CalculateZRocheLobe(const double p_jLoss) const;

    double  CalculateTimeToCoalescence(double a0, double e0, double m1, double m2) const;

    // CalculateTotalEnergy - the actual function takes 9 parameters because of the various calling permutations
    //                      - various signatures are defined here - they just assemble the parameters as required
    //                        and call the actual function
    double  CalculateTotalEnergy(const double p_SemiMajorAxis,
                                 const double p_Star1Mass,
                                 const double p_Star2Mass,
                                 const double p_Star1Radius,
                                 const double p_Star2Radius,
                                 const double p_Star1_SpinAngularVelocity,
                                 const double p_Star2_SpinAngularVelocity,
                                 const double p_Star1GyrationRadius,
                                 const double p_Star2GyrationRadius) const;

    double  CalculateTotalEnergy() const                                    { return CalculateTotalEnergy(m_SemiMajorAxis, m_Star1->Mass(), m_Star2->Mass(), m_Star1->Radius(), m_Star2->Radius(), m_Star1->Omega(), m_Star2->Omega(), m_Star1->CalculateGyrationRadius(), m_Star2->CalculateGyrationRadius()); }

    double  CalculateTotalEnergy(const double p_SemiMajorAxis,
                                 const double p_Star1_SpinAngularVelocity,
                                 const double p_Star2_SpinAngularVelocity,
                                 const double p_Star1_GyrationRadius,
                                 const double p_Star2_GyrationRadius) const { return CalculateTotalEnergy(p_SemiMajorAxis, m_Star1->Mass(), m_Star2->Mass(), m_Star1->Radius(), m_Star2->Radius(), p_Star1_SpinAngularVelocity, p_Star2_SpinAngularVelocity, p_Star1_GyrationRadius, p_Star2_GyrationRadius); }


    void    EvaluateBinary(const double p_Dt);

    void    EvaluateSupernovae();

    void    EvolveOneTimestep(const double p_Dt);
    void    EvolveOneTimestepPreamble(const double p_Dt);

    void    ResolveCoalescence();
    void    ResolveCommonEnvelopeEvent();
    void    ResolveMassChanges();
    bool    ResolveSupernova();

    void    SetPostCEEValues(const double p_SemiMajorAxis,
                             const double p_Eccentricity,
                             const double p_RocheLobe1to2,
                             const double p_RocheLobe2to1);

    void    SetPreCEEValues(const double p_SemiMajorAxis,
                            const double p_Eccentricity,
                            const double p_RocheLobe1to2,
                            const double p_RocheLobe2to1);

    void    StashBeBinaryProperties();
    void    StashRLOFProperties(const MASS_TRANSFER_TIMING p_Which);

    void    UpdateSystemicVelocity(Vector3d p_newVelocity);

    // printing functions
    bool PrintRLOFParameters();
    bool PrintBinarySystemParameters() const            { return LOGGING->LogBSESystemParameters(this); }
    bool PrintDetailedOutput(const long int p_Id) const { return OPTIONS->DetailedOutput() ? LOGGING->LogBSEDetailedOutput(this, p_Id) : true; }
    bool PrintDoubleCompactObjects() const              { return LOGGING->LogDoubleCompactObject(this); }
    bool PrintCommonEnvelope() const                    { return LOGGING->LogCommonEnvelope(this); }
    bool PrintBeBinary();
    bool PrintPulsarEvolutionParameters() const         { return OPTIONS->EvolvePulsars() ? LOGGING->LogBSEPulsarEvolutionParameters(this) : true; }
    bool PrintSupernovaDetails() const                  { return LOGGING->LogBSESupernovaDetails(this); }

    
    //Functor for the boost root finder to determine how much mass needs to be lost from a donor without an envelope in order to fit inside the Roche lobe
    template <class T>
    struct RadiusEqualsRocheLobeFunctor
    {
        RadiusEqualsRocheLobeFunctor(BaseBinaryStar *p_Binary, BinaryConstituentStar *p_Donor, BinaryConstituentStar *p_Accretor, ERROR *p_Error, double p_FractionAccreted)
        {
            m_Binary           = p_Binary;
            m_Donor            = p_Donor;
            m_Accretor         = p_Accretor;
            m_Error            = p_Error;
            m_FractionAccreted = p_FractionAccreted;
        }
        T operator()(double const& dM)
        {
            if (dM >= m_Donor->Mass()) {            // Can't remove more than the donor's mass
                *m_Error = ERROR::TOO_MANY_RLOF_ITERATIONS;
                return m_Donor->Radius();
            }

            double donorMass    = m_Donor->Mass();
            double accretorMass = m_Accretor->Mass();

            BinaryConstituentStar* donorCopy = new BinaryConstituentStar(*m_Donor);
            double semiMajorAxis = m_Binary->CalculateMassTransferOrbit(donorCopy->Mass(), -dM , *m_Accretor, m_FractionAccreted);
            double RLRadius      = semiMajorAxis * (1 - m_Binary->Eccentricity()) * CalculateRocheLobeRadius_Static(donorMass - dM, accretorMass + (m_Binary->FractionAccreted() * dM)) * AU_TO_RSOL;
            
            (void)donorCopy->UpdateAttributes(-dM, -dM*donorCopy->Mass0()/donorCopy->Mass());
            
            // Modify donor Mass0 and Age for MS (including HeMS) and HG stars
            donorCopy->UpdateInitialMass();         // update initial mass (MS, HG & HeMS)  
            donorCopy->UpdateAgeAfterMassLoss();    // update age (MS, HG & HeMS)
            
            (void)donorCopy->AgeOneTimestep(0.0);   // recalculate radius of star - don't age - just update values
            
            double thisRadiusAfterMassLoss = donorCopy->Radius();
            
            delete donorCopy; donorCopy = nullptr;
            
            return (RLRadius-thisRadiusAfterMassLoss);
        }
    private:
        BaseBinaryStar *m_Binary;
        BinaryConstituentStar *m_Donor;
        BinaryConstituentStar *m_Accretor;
        ERROR *m_Error;
        double m_FractionAccreted;
    };
    
  
    //Root solver to determine how much mass needs to be lost from a donor without an envelope in order to fit inside the Roche lobe
    double MassLossToFitInsideRocheLobe(BaseBinaryStar * p_Binary, BinaryConstituentStar * p_Donor, BinaryConstituentStar * p_Accretor, double p_FractionAccreted)
    {
        using namespace std;                                                    // Help ADL of std functions.
        using namespace boost::math::tools;                                     // For bracket_and_solve_root.
        
        double guess  = ADAPTIVE_RLOF_FRACTION_DONOR_GUESS * p_Donor->Mass();   // Rough guess at solution
        double factor = ADAPTIVE_RLOF_SEARCH_FACTOR;                            // Size of search steps
        
        const boost::uintmax_t maxit = ADAPTIVE_RLOF_MAX_ITERATIONS;            // Limit to maximum iterations.
        boost::uintmax_t it = maxit;                                            // Initally our chosen max iterations, but updated with actual.
        bool is_rising = true;                                                  // So if result with guess is too low, then try increasing guess.
        int digits = std::numeric_limits<double>::digits;                       // Maximum possible binary digits accuracy for type T.

        // Some fraction of digits is used to control how accurate to try to make the result.
        int get_digits = digits - 5;                                            // We have to have a non-zero interval at each step, so

        // maximum accuracy is digits - 1.  But we also have to
        // allow for inaccuracy in f(x), otherwise the last few
        // iterations just thrash around.
        eps_tolerance<double> tol(get_digits);                                  // Set the tolerance.
        
        std::pair<double, double> root;
        try {
            ERROR error = ERROR::NONE;
            root = bracket_and_solve_root(RadiusEqualsRocheLobeFunctor<double>(p_Binary, p_Donor, p_Accretor, &error, p_FractionAccreted), guess, factor, is_rising, tol, it);
            if (error != ERROR::NONE) SHOW_WARN(error);
        }
        catch(exception& e) {
            SHOW_ERROR(ERROR::TOO_MANY_RLOF_ITERATIONS, e.what());              // Catch generic boost root finding error
            m_Donor->Radius();
        }
        SHOW_WARN_IF(it>=maxit, ERROR::TOO_MANY_RLOF_ITERATIONS);
        
        return root.first + (root.second - root.first)/2;                       // Midway between brackets is our result, if necessary we could return the result as an interval here.
    }
};

#endif // __BaseBinaryStar_h__
