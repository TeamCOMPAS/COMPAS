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
 * List all conditions for SN (AIC or SN Ia) for ONeWD. 
 * Each condition should also be a separate clause in EvolveToNextPhase.
 *
 * bool IsSupernova()
 *
 * @return                               Whether WD should undergo AIC or SN Ia
 */
bool ONeWD::IsSupernova() const {
    return IsMassAboveEcsnThreshold();      
}


/*
 * Specifies next stage, if the star changes its phase.
 *
 * STELLAR_TYPE EvolveToNextPhase()
 *
 * @return                               Stellar type of the upcoming stage.
 */
STELLAR_TYPE ONeWD::EvolveToNextPhase() {
    STELLAR_TYPE stellarType;
    stellarType = ResolveAIC();
    return stellarType;
}

