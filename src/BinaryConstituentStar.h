#ifndef __BinaryConstituentStar_h__
#define __BinaryConstituentStar_h__

#include "constants.h"
#include "typedefs.h"
#include "utils.h"

#include "Star.h"


class Star;


class BinaryConstituentStar: virtual public Star {

public:


    BinaryConstituentStar() : Star() {
        m_ObjectId   = globalObjectId++;
        m_ObjectType = OBJECT_TYPE::BINARY_CONSTITUENT_STAR;
    };


    BinaryConstituentStar(const unsigned long int p_RandomSeed,
                          const double            p_Mass,
                          const double            p_Metallicity,
                          const KickParameters    p_KickParameters,
                          const double            p_LBVfactor,
                          const double            p_WolfRayetFactor) : Star(p_RandomSeed, p_Mass, p_Metallicity, p_KickParameters, p_LBVfactor, p_WolfRayetFactor) {

        m_ObjectId                 = globalObjectId++;
        m_ObjectType               = OBJECT_TYPE::BINARY_CONSTITUENT_STAR;

        m_Companion                = nullptr;


        m_CEDetails.bindingEnergy                    = DEFAULT_INITIAL_DOUBLE_VALUE;
        m_CEDetails.COCoreMass                       = DEFAULT_INITIAL_DOUBLE_VALUE;
        m_CEDetails.CoreMass                         = DEFAULT_INITIAL_DOUBLE_VALUE;
        m_CEDetails.HeCoreMass                       = DEFAULT_INITIAL_DOUBLE_VALUE;
        m_CEDetails.lambda                           = DEFAULT_INITIAL_DOUBLE_VALUE;

        m_CEDetails.preCEE.bindingEnergy             = DEFAULT_INITIAL_DOUBLE_VALUE;
        m_CEDetails.preCEE.luminosity                = DEFAULT_INITIAL_DOUBLE_VALUE;
        m_CEDetails.preCEE.mass                      = DEFAULT_INITIAL_DOUBLE_VALUE;
        m_CEDetails.preCEE.radius                    = DEFAULT_INITIAL_DOUBLE_VALUE;
        m_CEDetails.preCEE.stellarType               = STELLAR_TYPE::NONE;
        m_CEDetails.preCEE.temperature               = DEFAULT_INITIAL_DOUBLE_VALUE;
        m_CEDetails.preCEE.dynamicalTimescale        = DEFAULT_INITIAL_DOUBLE_VALUE;
        m_CEDetails.preCEE.thermalTimescale          = DEFAULT_INITIAL_DOUBLE_VALUE;
        m_CEDetails.preCEE.nuclearTimescale          = DEFAULT_INITIAL_DOUBLE_VALUE;
        m_CEDetails.preCEE.radialExpansionTimescale  = DEFAULT_INITIAL_DOUBLE_VALUE;

        m_CEDetails.postCEE.bindingEnergy            = DEFAULT_INITIAL_DOUBLE_VALUE;
        m_CEDetails.postCEE.luminosity               = DEFAULT_INITIAL_DOUBLE_VALUE;
        m_CEDetails.postCEE.temperature              = DEFAULT_INITIAL_DOUBLE_VALUE;
        m_CEDetails.postCEE.dynamicalTimescale       = DEFAULT_INITIAL_DOUBLE_VALUE;
        m_CEDetails.postCEE.thermalTimescale         = DEFAULT_INITIAL_DOUBLE_VALUE;
        m_CEDetails.postCEE.nuclearTimescale         = DEFAULT_INITIAL_DOUBLE_VALUE;
        m_CEDetails.postCEE.radialExpansionTimescale = DEFAULT_INITIAL_DOUBLE_VALUE;


        m_IsPrimary                                  = false;

        m_RocheLobeRadius                            = DEFAULT_INITIAL_DOUBLE_VALUE;

        m_MassTransferDiff                           = DEFAULT_INITIAL_DOUBLE_VALUE;
        m_MassLossDiff                               = DEFAULT_INITIAL_DOUBLE_VALUE;

        m_PreSNeOrbitalEnergy                        = DEFAULT_INITIAL_DOUBLE_VALUE;
        m_PostSNeOrbitalEnergy                       = DEFAULT_INITIAL_DOUBLE_VALUE;

        m_FastPhaseCaseA                             = false;

        m_FirstMassTransferEpisode                   = false;
        m_MassTransferCaseInitial                    = MT_CASE::NONE;

        m_OmegaTidesIndividualDiff                   = DEFAULT_INITIAL_DOUBLE_VALUE;

        m_RocheLobeTracker                           = DEFAULT_INITIAL_DOUBLE_VALUE;

        m_RLOFDetails.isRLOF                         = false;
        m_RLOFDetails.experiencedRLOF                = false;
        m_RLOFDetails.RLOFPostCEE                    = false;
    }


    // Copy constructor
    BinaryConstituentStar(const BinaryConstituentStar& p_Star) : Star(p_Star) {

        m_ObjectId                 = globalObjectId++;                      // get unique object id (don't copy source)
        m_ObjectType               = OBJECT_TYPE::BINARY_CONSTITUENT_STAR;  // can only copy from BINARY_CONSTITUENT_STAR

        m_Companion                = nullptr;                               // don't point at source companion - this can be updated separately later

        m_CEDetails                = p_Star.m_CEDetails;

        m_IsPrimary                = p_Star.m_IsPrimary;

        m_RocheLobeRadius          = p_Star.m_RocheLobeRadius;

        m_MassTransferDiff         = p_Star.m_MassTransferDiff;
        m_MassLossDiff             = p_Star.m_MassLossDiff;

        m_PreSNeOrbitalEnergy      = p_Star.m_PreSNeOrbitalEnergy;
        m_PostSNeOrbitalEnergy     = p_Star.m_PostSNeOrbitalEnergy;

        m_FastPhaseCaseA           = p_Star.m_FastPhaseCaseA;

        m_FirstMassTransferEpisode = p_Star.m_FirstMassTransferEpisode;
        m_MassTransferCaseInitial  = p_Star.m_MassTransferCaseInitial;

        m_OmegaTidesIndividualDiff = p_Star.m_OmegaTidesIndividualDiff;

        m_RocheLobeTracker         = p_Star.m_RocheLobeTracker;

        m_RLOFDetails              = p_Star.m_RLOFDetails;
    }


    // Assignment overload
    BinaryConstituentStar& operator = (const BinaryConstituentStar& p_Star) {

        if (this != &p_Star) {                                                  // make sure we're not not copying ourselves...

            m_ObjectId                 = globalObjectId++;                      // get unique object id (don't copy source)
            m_ObjectType               = OBJECT_TYPE::BINARY_CONSTITUENT_STAR;  // can only copy from BINARY_CONSTITUENT_STAR

            m_Companion                = nullptr;                               // don't point at source companion - this can be updated separately later

            m_CEDetails                = p_Star.m_CEDetails;

            m_IsPrimary                = p_Star.m_IsPrimary;

            m_RocheLobeRadius          = p_Star.m_RocheLobeRadius;

            m_MassTransferDiff         = p_Star.m_MassTransferDiff;
            m_MassLossDiff             = p_Star.m_MassLossDiff;

            m_PreSNeOrbitalEnergy      = p_Star.m_PreSNeOrbitalEnergy;
            m_PostSNeOrbitalEnergy     = p_Star.m_PostSNeOrbitalEnergy;

            m_FastPhaseCaseA           = p_Star.m_FastPhaseCaseA;

            m_FirstMassTransferEpisode = p_Star.m_FirstMassTransferEpisode;
            m_MassTransferCaseInitial  = p_Star.m_MassTransferCaseInitial;

            m_OmegaTidesIndividualDiff = p_Star.m_OmegaTidesIndividualDiff;

            m_RocheLobeTracker         = p_Star.m_RocheLobeTracker;

            m_RLOFDetails              = p_Star.m_RLOFDetails;
        }
        return *this;
    }

    ~BinaryConstituentStar() { }


    // object identifiers - all classes have these
    OBJECT_ID           ObjectId() const                                                { return m_ObjectId; }
    OBJECT_TYPE         ObjectType() const                                              { return m_ObjectType; }


    // getters - alphabetically
    double          BindingEnergyAtCEE() const                                          { return m_CEDetails.bindingEnergy; }
    double          COCoreMassAtCEE() const                                             { return m_CEDetails.COCoreMass; }
    double          CoreMassAtCEE() const                                               { return m_CEDetails.CoreMass; }
    double          HeCoreMassAtCEE() const                                             { return m_CEDetails.HeCoreMass; }
    double          LambdaAtCEE() const                                                 { return m_CEDetails.lambda; }

    double          BindingEnergyPostCEE() const                                        { return m_CEDetails.postCEE.bindingEnergy; }
    double          BindingEnergyPreCEE() const                                         { return m_CEDetails.preCEE.bindingEnergy; }
    double          DynamicalTimescalePostCEE() const                                   { return m_CEDetails.postCEE.dynamicalTimescale; }
    double          DynamicalTimescalePreCEE() const                                    { return m_CEDetails.preCEE.dynamicalTimescale; }
    double          LuminosityPostCEE() const                                           { return m_CEDetails.postCEE.luminosity; }
    double          LuminosityPreCEE() const                                            { return m_CEDetails.preCEE.luminosity; }
    double          MassPostCEE() const                                                 { return m_CEDetails.postCEE.mass; }
    double          MassPreCEE() const                                                  { return m_CEDetails.preCEE.mass; }
    double          NuclearTimescalePostCEE() const                                     { return m_CEDetails.postCEE.nuclearTimescale; }
    double          NuclearTimescalePreCEE() const                                      { return m_CEDetails.preCEE.nuclearTimescale; }
    double          RadialExpansionTimescalePostCEE() const                             { return m_CEDetails.postCEE.radialExpansionTimescale; }
    double          RadialExpansionTimescalePreCEE() const                              { return m_CEDetails.preCEE.radialExpansionTimescale; }
    double          RadiusPostCEE() const                                               { return m_CEDetails.postCEE.radius; }
    double          RadiusPreCEE() const                                                { return m_CEDetails.preCEE.radius; }
    STELLAR_TYPE    StellarTypePostCEE() const                                          { return m_CEDetails.postCEE.stellarType; }
    STELLAR_TYPE    StellarTypePreCEE() const                                           { return m_CEDetails.preCEE.stellarType; }
    double          TemperaturePostCEE() const                                          { return m_CEDetails.postCEE.temperature; }
    double          TemperaturePreCEE() const                                           { return m_CEDetails.preCEE.temperature; }
    double          ThermalTimescalePostCEE() const                                     { return m_CEDetails.postCEE.thermalTimescale; }
    double          ThermalTimescalePreCEE() const                                      { return m_CEDetails.preCEE.thermalTimescale; }

    bool            ExperiencedRLOF() const                                             { return m_RLOFDetails.experiencedRLOF; }
    bool            FastPhaseCaseA() const                                              { return m_FastPhaseCaseA ; }
    bool            IsPrimary() const                                                   { return m_IsPrimary; }
    bool            IsSNevent() const                                                   { return IsCCSN() || IsECSN() || IsPISN() || IsPPISN(); }
    bool            IsRLOF() const                                                      { return m_RLOFDetails.isRLOF; }
    double          MassLossDiff() const                                                { return m_MassLossDiff; }
    MT_CASE         MassTransferCaseInitial() const                                     { return m_MassTransferCaseInitial; }
    double          MassTransferDiff() const                                            { return m_MassTransferDiff; }
    double          OmegaTidesIndividualDiff() const                                    { return m_OmegaTidesIndividualDiff; }
    double          PostSNeOrbitalEnergy() const                                        { return m_PostSNeOrbitalEnergy; };
    double          PreSNeOrbitalEnergy() const                                         { return m_PreSNeOrbitalEnergy; };
    bool            RLOFPostCEE() const                                                 { return m_RLOFDetails.RLOFPostCEE; }
    double          RocheLobeRadius() const                                             { return m_RocheLobeRadius; }
    double          RocheLobeTracker() const                                            { return m_RocheLobeTracker; }


    // setters
    void            SetCompanion(BinaryConstituentStar* p_Companion)                    { m_Companion = p_Companion; }                              // this star's companion star

    void            SetFastPhaseCaseA()                                                 { m_FastPhaseCaseA = true; }                                // JR: todo: revisit this

    void            SetMassLossDiff(const double p_MassLossDiff)                        { m_MassLossDiff = p_MassLossDiff; }                        // JR: todo: better way?  JR: todo:  sanity check?
    void            SetMassTransferDiff(const double p_MassTransferDiff)                { m_MassTransferDiff = p_MassTransferDiff; }                // JR: todo: better way?  JR: todo:  sanity check?

    void            SetPostSNeOrbitalEnergy(const double p_PostSNeOrbitalEnergy)        { m_PostSNeOrbitalEnergy = p_PostSNeOrbitalEnergy; };
    void            SetPreSNeOrbitalEnergy(const double p_PreSNeOrbitalEnergy)          { m_PreSNeOrbitalEnergy = p_PreSNeOrbitalEnergy; };


    // member functions - alphabetically
    void            BecomePrimary()                                                     { m_IsPrimary = true; }
    void            BecomeSecondary()                                                   { m_IsPrimary = false; }

    void            CalculateCommonEnvelopeValues();

    void            CalculateOmegaTidesIndividualDiff(const double p_OrbitalVelocity)   { m_OmegaTidesIndividualDiff = p_OrbitalVelocity - OmegaPrev(); }

    double          CalculateCircularisationTimescale(const double p_SemiMajorAxis);

    double          CalculateSynchronisationTimescale(const double p_SemiMajorAxis);

    void            DetermineInitialMassTransferCase();

    void            InitialiseMassTransfer(const bool p_CommonEnvelope, const double p_SemiMajorAxis, const double p_Eccentricity);


    void            ResolveCommonEnvelopeAccretion(const double p_FinalMass);

    void            SetPostCEEValues();
    void            SetPreCEEValues();

    COMPAS_VARIABLE StellarPropertyValue(const T_ANY_PROPERTY p_Property) const;

    void            UpdateMagneticFieldAndSpin(const bool   p_CommonEnvelope,
                                               const double p_Stepsize,
                                               const double p_Epsilon)                  { Star::UpdateMagneticFieldAndSpin(p_CommonEnvelope, p_Stepsize, m_MassTransferDiff * MSOL, p_Epsilon); }  // JR: todo: revisit this


private:

    OBJECT_ID               m_ObjectId;                             // Instantiated object's unique object id
    OBJECT_TYPE             m_ObjectType;                           // Instantiated object's object type

    // member variables - alphabetically

    StellarCEDetailsT       m_CEDetails;                            // Common envelope details

    bool                    m_IsPrimary;

    bool 	                m_FastPhaseCaseA;                       // Indicates if the star just entered a case A MT for the first time        JR: todo: cf binary value of same name

    bool                    m_FirstMassTransferEpisode;             // Activated for the initial Mass Transfer Episode

    double                  m_MassLossDiff;
    MT_CASE                 m_MassTransferCaseInitial;              // Indicator of which Mass Transfer occures when first RLOF, if any
    double                  m_MassTransferDiff;

    double                  m_OmegaTidesIndividualDiff;

    double                  m_PostSNeOrbitalEnergy;
    double                  m_PreSNeOrbitalEnergy;

	StellarRLOFDetailsT     m_RLOFDetails;

    double                  m_RocheLobeRadius;

    double                  m_RocheLobeTracker;


    // the companion - set by calling SetCompanion()
    BinaryConstituentStar *m_Companion;


	// member functions - alphabetically
    void                CalculateInitialMassTransferCase();

    double              CalculateMassAccretedForNS(const double p_CompanionMass, const double p_CompanionRadius);

    void                SetRocheLobeFlags(const bool p_CommonEnvelope, const double p_SemiMajorAxis, const double p_Eccentricity);

};

#endif // __BinaryConstituentStar_h__
