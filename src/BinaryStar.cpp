#include "BinaryStar.h"


// binary is generated according to distributions specified in program options
BinaryStar::BinaryStar(const unsigned long int p_Seed, const long int p_Id) : m_BinaryStar(new BaseBinaryStar(p_Seed, p_Id)) {

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


/*
 * Print BSE Switch Log record
 * 
 * Called from main() when SIGUSR1 received - raised by Star::SwitchTo() to indicate a stellar type switch.
 * Here we use the switch parameters stored in the LOGGING service singleton by Star:SwitchTo() to determine
 * whether it is the primary or the secondary switching, then call BinaryStar::PrintSwitchLog() with the
 * appropriate parameters to print the log file record.
 * 
 * BinaryStar::PrintSwitchLog()
 *
 * @return                                      Boolean flag indicating success/failure (true = success)
 */
bool BinaryStar::PrintSwitchLog() { 
    
    bool result = true;

    OBJECT_ID primaryObjectId   = m_BinaryStar->Star1()->StarObjectId();
    OBJECT_ID secondaryObjectId = m_BinaryStar->Star2()->StarObjectId();
    OBJECT_ID objectIdSwitching = LOGGING->ObjectIdSwitching();

         if (objectIdSwitching == primaryObjectId  ) result = m_BinaryStar->PrintSwitchLog(true);   // primary
    else if (objectIdSwitching == secondaryObjectId) result = m_BinaryStar->PrintSwitchLog(false);  // secondary
    else {                                                                                          // otherwise...
        SHOW_ERROR(ERROR::OUT_OF_BOUNDS, "Expected primary or secondary for BSE Switch Log");       // announce error
        result = false;
    }

    return result;
}
