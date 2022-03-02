#ifndef __BinaryConstituentStar_h__
#define __BinaryConstituentStar_h__

#include "constants.h"
#include "typedefs.h"
#include "profiling.h"
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
                          const STELLAR_TYPE      p_InitialStellarType,
                          const double            p_Metallicity, 
                          const KickParameters    p_KickParameters,
                          const double            p_RotationalVelocity = -1.0) : Star(p_RandomSeed, p_Mass, p_InitialStellarType, p_Metallicity, p_KickParameters, p_RotationalVelocity) {

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

        m_CEDetails.postCEE.luminosity               = DEFAULT_INITIAL_DOUBLE_VALUE;
        m_CEDetails.postCEE.temperature              = DEFAULT_INITIAL_DOUBLE_VALUE;
        m_CEDetails.postCEE.dynamicalTimescale       = DEFAULT_INITIAL_DOUBLE_VALUE;
        m_CEDetails.postCEE.thermalTimescale         = DEFAULT_INITIAL_DOUBLE_VALUE;
        m_CEDetails.postCEE.nuclearTimescale         = DEFAULT_INITIAL_DOUBLE_VALUE;
        m_CEDetails.postCEE.radialExpansionTimescale = DEFAULT_INITIAL_DOUBLE_VALUE;

        m_Flags.recycledNS                           = false;
        m_Flags.rlofOntoNS                           = false;

        m_MassTransferDiff                           = DEFAULT_INITIAL_DOUBLE_VALUE;
        m_MassLossDiff                               = DEFAULT_INITIAL_DOUBLE_VALUE;

        m_OrbitalEnergyPreSN                         = DEFAULT_INITIAL_DOUBLE_VALUE;
        m_OrbitalEnergyPostSN                        = DEFAULT_INITIAL_DOUBLE_VALUE;

        m_FirstMassTransferEpisode                   = false;

        m_OmegaTidesIndividualDiff                   = DEFAULT_INITIAL_DOUBLE_VALUE;

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

        m_Flags                    = p_Star.m_Flags;

        m_MassTransferDiff         = p_Star.m_MassTransferDiff;
        m_MassLossDiff             = p_Star.m_MassLossDiff;

        m_OrbitalEnergyPreSN       = p_Star.m_OrbitalEnergyPreSN;
        m_OrbitalEnergyPostSN      = p_Star.m_OrbitalEnergyPostSN;

        m_FirstMassTransferEpisode = p_Star.m_FirstMassTransferEpisode;

        m_OmegaTidesIndividualDiff = p_Star.m_OmegaTidesIndividualDiff;

        m_RLOFDetails              = p_Star.m_RLOFDetails;
    }


    // Assignment overload
    BinaryConstituentStar& operator = (const BinaryConstituentStar& p_Star) {

        if (this != &p_Star) {                                                  // make sure we're not not copying ourselves...

            m_ObjectId                 = globalObjectId++;                      // get unique object id (don't copy source)
            m_ObjectType               = OBJECT_TYPE::BINARY_CONSTITUENT_STAR;  // can only copy from BINARY_CONSTITUENT_STAR

            m_Companion                = nullptr;                               // don't point at source companion - this can be updated separately later

            m_CEDetails                = p_Star.m_CEDetails;

            m_Flags                    = p_Star.m_Flags;
            
            m_MassTransferDiff         = p_Star.m_MassTransferDiff;
            m_MassLossDiff             = p_Star.m_MassLossDiff;

            m_OrbitalEnergyPreSN       = p_Star.m_OrbitalEnergyPreSN;
            m_OrbitalEnergyPostSN      = p_Star.m_OrbitalEnergyPostSN;

            m_FirstMassTransferEpisode = p_Star.m_FirstMassTransferEpisode;

            m_OmegaTidesIndividualDiff = p_Star.m_OmegaTidesIndividualDiff;

            m_RLOFDetails              = p_Star.m_RLOFDetails;
        }
        return *this;
    }

    ~BinaryConstituentStar() { }


    // object identifiers - all classes have these
    OBJECT_ID       ObjectId() const                                                    { return m_ObjectId; }
    OBJECT_TYPE     ObjectType() const                                                  { return m_ObjectType; }


    // getters - alphabetically
    double          BindingEnergyAtCEE() const                                          { return m_CEDetails.bindingEnergy; }
    double          BindingEnergyPreCEE() const                                         { return m_CEDetails.preCEE.bindingEnergy; }

    double          COCoreMassAtCEE() const                                             { return m_CEDetails.COCoreMass; }
    double          CoreMassAtCEE() const                                               { return m_CEDetails.CoreMass; }

    double          DynamicalTimescalePostCEE() const                                   { return m_CEDetails.postCEE.dynamicalTimescale; }
    double          DynamicalTimescalePreCEE() const                                    { return m_CEDetails.preCEE.dynamicalTimescale; }

    bool            ExperiencedRecycledNS() const                                       { return m_Flags.recycledNS; }
    bool            ExperiencedRLOF() const                                             { return m_RLOFDetails.experiencedRLOF; }
    bool            ExperiencedRLOFOntoNS() const                                       { return m_Flags.rlofOntoNS; }

    double          HeCoreMassAtCEE() const                                             { return m_CEDetails.HeCoreMass; }

    bool            IsRLOF() const                                                      { return m_RLOFDetails.isRLOF; }
    bool            IsSNevent() const                                                   { return IsCCSN() || IsECSN() || IsPISN() || IsPPISN(); }

    double          LambdaAtCEE() const                                                 { return m_CEDetails.lambda; }
    double          LuminosityPostCEE() const                                           { return m_CEDetails.postCEE.luminosity; }
    double          LuminosityPreCEE() const                                            { return m_CEDetails.preCEE.luminosity; }

    double          MassLossDiff() const                                                { return m_MassLossDiff; }
    double          MassPostCEE() const                                                 { return m_CEDetails.postCEE.mass; }
    double          MassPreCEE() const                                                  { return m_CEDetails.preCEE.mass; }
    double          MassTransferDiff() const                                            { return m_MassTransferDiff; }

    double          NuclearTimescalePostCEE() const                                     { return m_CEDetails.postCEE.nuclearTimescale; }
    double          NuclearTimescalePreCEE() const                                      { return m_CEDetails.preCEE.nuclearTimescale; }

    double          OmegaTidesIndividualDiff() const                                    { return m_OmegaTidesIndividualDiff; }
    double          OrbitalEnergyPostSN() const                                         { return m_OrbitalEnergyPostSN; };
    double          OrbitalEnergyPreSN() const                                          { return m_OrbitalEnergyPreSN; };

    double          RadialExpansionTimescalePostCEE() const                             { return m_CEDetails.postCEE.radialExpansionTimescale; }
    double          RadialExpansionTimescalePreCEE() const                              { return m_CEDetails.preCEE.radialExpansionTimescale; }
    double          RadiusPostCEE() const                                               { return m_CEDetails.postCEE.radius; }
    double          RadiusPreCEE() const                                                { return m_CEDetails.preCEE.radius; }
    bool            RLOFPostCEE() const                                                 { return m_RLOFDetails.RLOFPostCEE; }
    double          StarToRocheLobeRadiusRatio(const double p_SemiMajorAxis, const double p_Eccentricity);

    STELLAR_TYPE    StellarTypePostCEE() const                                          { return m_CEDetails.postCEE.stellarType; }
    STELLAR_TYPE    StellarTypePreCEE() const                                           { return m_CEDetails.preCEE.stellarType; }

    double          TemperaturePostCEE() const                                          { return m_CEDetails.postCEE.temperature; }
    double          TemperaturePreCEE() const                                           { return m_CEDetails.preCEE.temperature; }
    double          ThermalTimescalePostCEE() const                                     { return m_CEDetails.postCEE.thermalTimescale; }
    double          ThermalTimescalePreCEE() const                                      { return m_CEDetails.preCEE.thermalTimescale; }


    // setters
    void            SetCompanion(BinaryConstituentStar* p_Companion)                    { m_Companion = p_Companion; }                              // this star's companion star

    void            SetMassLossDiff(const double p_MassLossDiff)                        { m_MassLossDiff = p_MassLossDiff; }                        // JR: todo: better way?  JR: todo:  sanity check?
    void            SetMassTransferDiff(const double p_MassTransferDiff)                { m_MassTransferDiff = p_MassTransferDiff; }                // JR: todo: better way?  JR: todo:  sanity check?

    void            SetOrbitalEnergyPostSN(const double p_OrbitalEnergyPostSN)          { m_OrbitalEnergyPostSN = p_OrbitalEnergyPostSN; };
    void            SetOrbitalEnergyPreSN(const double p_OrbitalEnergyPreSN)            { m_OrbitalEnergyPreSN = p_OrbitalEnergyPreSN; };

    void            ClearRecycledNS()                                                   { m_Flags.recycledNS = false; }
    void            SetRecycledNS()                                                     { m_Flags.recycledNS = true; }

    void            ClearRLOFOntoNS()                                                   { m_Flags.rlofOntoNS = false; }
    void            SetRLOFOntoNS()                                                     { m_Flags.rlofOntoNS = true; }

    void            CalculateCommonEnvelopeValues();

    void            CalculateOmegaTidesIndividualDiff(const double p_OrbitalAngularVelocity) { m_OmegaTidesIndividualDiff = p_OrbitalAngularVelocity - OmegaPrev(); }

    double          CalculateCircularisationTimescale(const double p_SemiMajorAxis);

    double          CalculateSynchronisationTimescale(const double p_SemiMajorAxis);

    void            InitialiseMassTransfer(const bool p_CommonEnvelope, const double p_SemiMajorAxis, const double p_Eccentricity);

    void            ResolveCommonEnvelopeAccretion(const double p_FinalMass);

    void            SetPostCEEValues();
    void            SetPreCEEValues();

    COMPAS_VARIABLE StellarPropertyValue(const T_ANY_PROPERTY p_Property) const;

    void            UpdateMagneticFieldAndSpin(const bool   p_CommonEnvelope,
                                               const double p_Stepsize,
                                               const double p_Epsilon)                  { Star::UpdateMagneticFieldAndSpin(p_CommonEnvelope, 
                                                                                                                           ExperiencedRecycledNS(), 
                                                                                                                           p_Stepsize, 
                                                                                                                           m_MassTransferDiff * MSOL_TO_KG, p_Epsilon); }  // JR: todo: revisit this


private:

    OBJECT_ID               m_ObjectId;                             // Instantiated object's unique object id
    OBJECT_TYPE             m_ObjectType;                           // Instantiated object's object type

    // member variables - alphabetically

    StellarCEDetailsT       m_CEDetails;                            // Common envelope details

    bool                    m_FirstMassTransferEpisode;             // Activated for the initial Mass Transfer Episode

    struct FLAGS {                                                  // Miscellaneous flags

        bool recycledNS;                                            // Indicate whether the accretor was a recycled neutron star
        bool rlofOntoNS;                                            // Indicates whether the donor donated mass to neutron star through RLOF

    }                       m_Flags;

    double                  m_MassLossDiff;
    double                  m_MassTransferDiff;

    double                  m_OmegaTidesIndividualDiff;

    double                  m_OrbitalEnergyPostSN;
    double                  m_OrbitalEnergyPreSN;

	StellarRLOFDetailsT     m_RLOFDetails;


    // the companion - set by calling SetCompanion()
    BinaryConstituentStar *m_Companion;


	// member functions - alphabetically
    double              CalculateMassAccretedForNS(const double p_CompanionMass, const double p_CompanionRadius);

    void                SetRocheLobeFlags(const bool p_CommonEnvelope, const double p_SemiMajorAxis, const double p_Eccentricity);

};

#endif // __BinaryConstituentStar_h__
