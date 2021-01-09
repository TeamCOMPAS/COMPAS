#ifndef __BinaryStar_h__
#define __BinaryStar_h__

#include "constants.h"
#include "typedefs.h"

#include "BaseBinaryStar.h"


class BaseBinaryStar;


class BinaryStar {

public:


    /* Constructors
     *
     * Parameter p_Id is optional, and is only included so that comparison tests can
     * be run against the legacy Compas code.  If a fixed random seed is being used
     * (program option) the legacy code effectivley adds the loop index of the binary
     * (from COMPASBinary() in main.cpp) to the user-specified fixed random seed so
     * that each binary has a repeatable random seed.
     *
     * Notes: the legacy code doesn't actually use the loop index - it uses a generated
     * object id that is the same as the loop index.  The new code also assigns objects
     * object ids, but the ids are assigned to all objects, not just binary stars, so
     * the ids generated in the new code won't match the legacy code ids - hence the
     * need to use the loop index here.  The parameter is optional - if no comparison
     * testing against the legacy code is required, the p_Id parameter can be let default
     * (in which case it is not used to generate the random seed - the generated object
     * id is used instead).
     */

    BinaryStar(const long int p_Id = -1l);

    // Copy constructor
    BinaryStar(const BinaryStar& p_Star) {

        m_ObjectId   = globalObjectId++;                                                        // get unique object id (don't copy source)
        m_ObjectType = OBJECT_TYPE::BINARY_STAR;                                                // can only copy from BINARY_STAR

        m_BinaryStar     = new BaseBinaryStar(*(p_Star.m_BinaryStar));                          // copy underlying BaseBinaryStar
        m_SaveBinaryStar = new BaseBinaryStar(*(p_Star.m_SaveBinaryStar));                      // copy underlying Saved BaseBinaryStar
    }


    // Assignment overload
    BinaryStar& operator = (const BinaryStar& p_Star) {

        if (this != &p_Star) {                                                                  // make sure we're not not copying ourselves...

            m_ObjectId   = globalObjectId++;                                                    // get unique object id (don't copy source)
            m_ObjectType = OBJECT_TYPE::BINARY_CONSTITUENT_STAR;                                // can only copy from BINARY_CONSTITUENT_STAR

            delete m_BinaryStar;                                                                // delete existing
            m_BinaryStar     = new BaseBinaryStar(*(p_Star.m_BinaryStar));                      // copy underlying BaseBinaryStar
            delete m_SaveBinaryStar;                                                            // delete existing
            m_SaveBinaryStar = new BaseBinaryStar(*(p_Star.m_SaveBinaryStar));                  // copy underlying Saved BaseBinaryStar

        }
        return *this;
    }

    virtual ~BinaryStar() { delete m_BinaryStar; delete m_SaveBinaryStar; }


    // object identifiers - all classes have these
    OBJECT_ID           ObjectId() const            { return m_ObjectId; }
    OBJECT_TYPE         ObjectType() const          { return m_ObjectType; }
    STELLAR_TYPE        StellarType() const         { return m_StellarType; }


    // member functions
    long int            Id()                        { return m_BinaryStar->Id(); }
    EVOLUTION_STATUS    Evolve()                    { return m_BinaryStar->Evolve(); }
    bool                RevertState();
    void                SaveState();
    STELLAR_TYPE        Star1InitialType()          { return m_BinaryStar->InitialStellarType1(); }
    STELLAR_TYPE        Star1Type()                 { return m_BinaryStar->StellarType1(); }
    STELLAR_TYPE        Star2InitialType()          { return m_BinaryStar->InitialStellarType2(); }
    STELLAR_TYPE        Star2Type()                 { return m_BinaryStar->StellarType2(); }

    bool                PrintSwitchLog();

private:

    BinaryStar() { }

    OBJECT_ID       m_ObjectId;                                                                 // Instantiated object's unique object id
    OBJECT_TYPE     m_ObjectType;                                                               // Instantiated object's object type
    STELLAR_TYPE    m_StellarType;                                                              // Stellar type defined in Hurley et al. 2000

    BaseBinaryStar *m_BinaryStar;                                                               // Pointer to current binary star
    BaseBinaryStar *m_SaveBinaryStar;                                                           // Pointer to saved binary star

};

#endif // __BinaryStar_h__
