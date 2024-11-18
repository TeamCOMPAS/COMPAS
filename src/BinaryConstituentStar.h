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

    BinaryConstituentStar() : Star() { };

    BinaryConstituentStar(const unsigned long int p_RandomSeed, 
                          const double            p_Mass, 
                          const double            p_Metallicity, 
                          const KickParameters    p_KickParameters,
                          const double            p_RotationalVelocity = -1.0) : Star(p_RandomSeed, p_Mass, p_Metallicity, p_KickParameters, p_RotationalVelocity) {

        m_Companion                                  = nullptr;


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
        m_CEDetails.preCEE.radialExpansionTimescale  = DEFAULT_INITIAL_DOUBLE_VALUE;

        m_CEDetails.postCEE.luminosity               = DEFAULT_INITIAL_DOUBLE_VALUE;
        m_CEDetails.postCEE.temperature              = DEFAULT_INITIAL_DOUBLE_VALUE;
        m_CEDetails.postCEE.dynamicalTimescale       = DEFAULT_INITIAL_DOUBLE_VALUE;
        m_CEDetails.postCEE.thermalTimescale         = DEFAULT_INITIAL_DOUBLE_VALUE;
        m_CEDetails.postCEE.radialExpansionTimescale = DEFAULT_INITIAL_DOUBLE_VALUE;

        m_Flags.recycledNS                           = false;
        m_Flags.rlofOntoNS                           = false;

        m_MassLossDiff                               = DEFAULT_INITIAL_DOUBLE_VALUE;
        m_MassTransferDiff                           = DEFAULT_INITIAL_DOUBLE_VALUE;

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

        m_CEDetails                = p_Star.m_CEDetails;

        m_Flags                    = p_Star.m_Flags;

        m_MassLossDiff             = p_Star.m_MassLossDiff;
        m_MassTransferDiff         = p_Star.m_MassTransferDiff;

        m_OrbitalEnergyPreSN       = p_Star.m_OrbitalEnergyPreSN;
        m_OrbitalEnergyPostSN      = p_Star.m_OrbitalEnergyPostSN;

        m_FirstMassTransferEpisode = p_Star.m_FirstMassTransferEpisode;

        m_OmegaTidesIndividualDiff = p_Star.m_OmegaTidesIndividualDiff;

        m_RLOFDetails              = p_Star.m_RLOFDetails;

        // This BinaryConstituentStar's companion is an instance of the BinaryConstituentStar class
        // Here we copy the pointer to the companion object - note that it could be a nullptr if the
        // companion object has not yet been set.  Also note that when this BinaryConstituentStar
        // star object is deleted (if it is deleted), the object to which m_Companion points will not
        // be deleted as a result.
        m_Companion                = p_Star.m_Companion;
    }


    /*
     * This function should be used to clone a BinaryConstituentStar.
     * 
     * Important:
     * 
     * This function returns a pointer, created by the 'new' operator.  The 'new' operator dynamically allocates
     * memory on the heap, not the stack, which is why it is available to the caller of this function after this
     * function has exited and its stack frame collapsed.  It is the responsibility of the caller of this function
     * to delete the pointer returned when it is no longer required so that the allocated memory is return to the
     * pool of available memory - failing to do so will cause a memory leak and the program will eventually exhaust
     * available memory and fail.  The preferred usage pattern is:
     * 
     *     T* ptr = Clone(obj, persistence)
     *     ...
     *     ...
     *     delete ptr; ptr = nullptr;
     * 
     * 
     * template <class T1>
     * static BinaryConstituentStar* Clone(BinaryConstituentStar& p_Star, const OBJECT_PERSISTENCE p_Persistence)
     * 
     * @param   [IN]    p_Star                      (address of) The star to be cloned
     * @param   [IN]    p_Persistence               Specifies the object persistence to be assigned to the cloned star.
     *                                              If the cloned star is intended to be used temporarily (e.g. for hypothesis testing),
     *                                              persistence should be EPHEMERAL, otherwise PERMANENT.
     * @return                                      (pointer to) The cloned star
     */
    static BinaryConstituentStar* Clone(BinaryConstituentStar& p_Star, const OBJECT_PERSISTENCE p_Persistence) { 
        BinaryConstituentStar* ptr = new BinaryConstituentStar(p_Star); 
        ptr->SetPersistence(p_Persistence); 
        return ptr;
    }
    static BinaryConstituentStar* Clone(BinaryConstituentStar& p_Star) { return Clone(p_Star, p_Star.ObjectPersistence()); }

    ~BinaryConstituentStar() { }


    // object identifiers - all classes have these
    OBJECT_ID           ObjectId() const                                                { return m_ObjectId; }
    OBJECT_TYPE         ObjectType() const                                              { return OBJECT_TYPE::BINARY_CONSTITUENT_STAR; }
    OBJECT_PERSISTENCE  ObjectPersistence() const                                       { return m_ObjectPersistence; }


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
    bool            IsSNevent() const                                                   { return IsCCSN() || IsECSN() || IsPISN() || IsPPISN() || IsAIC() || IsSNIA() || IsHeSD(); }

    double          LambdaAtCEE() const                                                 { return m_CEDetails.lambda; }
    double          LuminosityPostCEE() const                                           { return m_CEDetails.postCEE.luminosity; }
    double          LuminosityPreCEE() const                                            { return m_CEDetails.preCEE.luminosity; }
    double          MassLossDiff() const                                                { return m_MassLossDiff; }
    double          MassPostCEE() const                                                 { return m_CEDetails.postCEE.mass; }
    double          MassPreCEE() const                                                  { return m_CEDetails.preCEE.mass; }
    double          MassTransferDiff() const                                            { return m_MassTransferDiff; }

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
    void            SetMassTransferDiffAndResolveWDShellChange(const double p_MassTransferDiff);

    void            SetOrbitalEnergyPostSN(const double p_OrbitalEnergyPostSN)          { m_OrbitalEnergyPostSN = p_OrbitalEnergyPostSN; };
    void            SetOrbitalEnergyPreSN(const double p_OrbitalEnergyPreSN)            { m_OrbitalEnergyPreSN = p_OrbitalEnergyPreSN; };

    void            ClearRecycledNS()                                                   { m_Flags.recycledNS = false; }
    void            SetRecycledNS()                                                     { m_Flags.recycledNS = true; }

    void            ClearRLOFOntoNS()                                                   { m_Flags.rlofOntoNS = false; }
    void            SetRLOFOntoNS()                                                     { m_Flags.rlofOntoNS = true; }

    void            CalculateCommonEnvelopeValues();

    double          CalculateCircularisationTimescale(const double p_SemiMajorAxis);

    double          CalculateSynchronisationTimescale(const double p_SemiMajorAxis);

    void            InitialiseMassTransfer(const bool p_CommonEnvelope, const double p_SemiMajorAxis, const double p_Eccentricity);

    void            ResolveCommonEnvelopeAccretion(const double p_FinalMass,
                                                   const double p_CompanionMass     = 0.0,
                                                   const double p_CompanionRadius   = 0.0,
                                                   const double p_CompanionEnvelope = 0.0) { m_MassTransferDiff = Star::ResolveCommonEnvelopeAccretion(p_FinalMass,
                                                                                                                                                       m_Companion->Mass(),
                                                                                                                                                       m_Companion->Radius(),
                                                                                                                                                       m_Companion->MassPreCEE() - m_Companion->CoreMassAtCEE());
                                                                                             ResolveAccretion(m_MassTransferDiff);
                                                     
                                                                                           }  

    void            SetPostCEEValues();
    void            SetPreCEEValues();

    void            SetRocheLobeFlags(const bool p_CommonEnvelope, const double p_SemiMajorAxis, const double p_Eccentricity);

    COMPAS_VARIABLE StellarPropertyValue(const T_ANY_PROPERTY p_Property) const;

    void            UpdateMagneticFieldAndSpin(const bool   p_CommonEnvelope,
                                               const double p_Stepsize,
                                               const double p_Epsilon)                  { Star::UpdateMagneticFieldAndSpin(p_CommonEnvelope, 
                                                                                                                           ExperiencedRecycledNS(), 
                                                                                                                           p_Stepsize, 
                                                                                                                           m_MassTransferDiff * MSOL_TO_KG, p_Epsilon); }

    void            SetMassLossDiff(const double p_MassLossDiff)                        { m_MassLossDiff = p_MassLossDiff; }                        // JR: todo: better way?  Sanity check?
    void            SetObjectId(const OBJECT_ID p_ObjectId)                             { m_ObjectId = p_ObjectId; }
    void            SetPersistence(const OBJECT_PERSISTENCE p_Persistence)              { m_ObjectPersistence = p_Persistence; }


private:

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
    BinaryConstituentStar  *m_Companion;


	// member functions - alphabetically 
};

#endif // __BinaryConstituentStar_h__
