#include "COWD.h"

/*
 * Specifies next stage, if the star changes its phase.
 *
 * STELLAR_TYPE EvolveToNextPhase()
 *
 * @return                               Stellar type of the upcoming stage.
 */

STELLAR_TYPE COWD::EvolveToNextPhase() {

    STELLAR_TYPE stellarType;

    if (m_OffCenterIgnition) {
        stellarType = STELLAR_TYPE::OXYGEN_NEON_WHITE_DWARF;
    }
    else {                                         
        stellarType = ResolveSNIa(); 
    }
    return stellarType;
}
