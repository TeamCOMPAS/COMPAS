#include "ONeWD.h"

/*
 * Allow evolution to a new phase (currently, only SN)
 *
 * bool ShouldEvolveOnPhase()
 *
 * @return                               Whether the WD should evolve on phase or towards a SN.
 */
bool ONeWD::ShouldEvolveOnPhase() const {
    return !IsSupernova();
}


/*
 * List all conditions for SN (AIC) for ONeWD.
 * Each condition should also be a separate clause in EvolveToNextPhase.
 *
 * bool IsSupernova()
 *
 * @return                               Whether WD should undergo AIC
 */
bool ONeWD::IsSupernova() const {
    return IsMassAboveChandrasekhar();
}


/*
 * Specifies next stage, if the star changes its phase.
 *
 * STELLAR_TYPE EvolveToNextPhase()
 *
 * @return                               Stellar type of the upcoming stage.
 */
STELLAR_TYPE ONeWD::EvolveToNextPhase() {
    return ResolveAIC();
}

