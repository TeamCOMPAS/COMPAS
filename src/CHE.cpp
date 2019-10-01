#include "CHE.h"


STELLAR_TYPE CHE::EvolveToNextPhase() {

    STELLAR_TYPE stellarType = STELLAR_TYPE::MS_GT_07;

    if (m_Age < m_Timescales[static_cast<int>(TIMESCALE::tMS)]) {           // evolving off because of age?
        stellarType = STELLAR_TYPE::MS_GT_07;                               // no - must have spun down - evolve as MS star now
    }
    else {                                                                  // yes
        stellarType = STELLAR_TYPE::NAKED_HELIUM_STAR_MS;                   // evolve as HeMS star now
        m_Age       = 0.0;                                                  // JR: can't use Hurley et al. 2000, eq 76 here - timescales(tHe) not calculated yet
        m_Tau       = 0.0;
    }

    return stellarType;
}
