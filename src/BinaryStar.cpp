#include "BinaryStar.h"


// Constructor where binary is generated according to distributions specified in options
BinaryStar::BinaryStar(const AIS &p_AIS, const long int p_Id) : m_BinaryStar(new BaseBinaryStar(p_AIS, p_Id)) {

    m_ObjectId       = globalObjectId++;
    m_ObjectType     = OBJECT_TYPE::BINARY_STAR;
    m_StellarType    = STELLAR_TYPE::BINARY_STAR;

    m_SaveBinaryStar = nullptr;
}


/*
 * Save current state of binary star
 *
 * Deletes existing pointer to saved BaseBinaryStar and instantiates new, copied, BaseBinaryStar
 *
 *
 * void SaveState()
 */
void BinaryStar::SaveState() {

    delete m_SaveBinaryStar;            // delete existing saved BaseBinaryStar
    m_SaveBinaryStar = m_BinaryStar;    // copy the underlying BaseBinaryStar
}


/*
 * Revert to the saved state of the binary star
 *
 * Changes the current binary star pointer (m_BinaryStar) to point to the save binary star
 * object and set the saved binary star pointer (m_SaveBinaryStar) to null.  Setting the
 * saved state pointer to null means there will be no saved state after calling this function -
 * so state needs to be saved if necessary (I don't do it here because we may not need to).
 *
 *
 * bool RevertState()
 *
 * @return                                      Boolean flag indicating success/failure (true = success)
 */
bool BinaryStar::RevertState() {

    bool result = false;

    if (m_SaveBinaryStar) {
        delete m_BinaryStar;
        m_BinaryStar     = m_SaveBinaryStar;
        m_SaveBinaryStar = nullptr;
        result           = true;
    }

    return result;
}
