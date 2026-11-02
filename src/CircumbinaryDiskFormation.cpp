#include "CircumbinaryDisk.h"

#include "constants.h"
#include "Options.h"
#include "utils.h"

#include <cmath>

namespace CircumbinaryDisk {


/*
 * Copy the current global COMPAS options needed by the CBD formation model
 *
 * FormationOptions MakeFormationOptionsFromCompasOptions()
 *
 * @return                                      CBD formation options
 */
FormationOptions MakeFormationOptionsFromCompasOptions() {
    return {
        OPTIONS->CircumbinaryDiskBeta(),
        OPTIONS->CircumbinaryDiskEnvelopeMassFractionSupplied(),
        OPTIONS->CircumbinaryDiskStructureFactor(),
        OPTIONS->CircumbinaryDiskInnerRadiusOverSeparation()
    };
}

/*
 * Copy the current global COMPAS options needed by the CBD evolution model
 *
 * EvolutionOptions MakeEvolutionOptionsFromCompasOptions()
 *
 * @return                                      CBD evolution options
 */
EvolutionOptions MakeEvolutionOptionsFromCompasOptions() {
    return {
        OPTIONS->CircumbinaryDiskEddingtonMode(),
        OPTIONS->CircumbinaryDiskEddingtonCapFactor(),
        OPTIONS->CircumbinaryDiskEddingtonGe23Factor()
    };
}

/*
 * Determine whether a set of supernova-like events destroys an active CBD
 *
 * AIC is treated as a change to the live COMPAS stellar state, so a bound
 * system continues evolving through the same physical disk.
 *
 * bool SupernovaEventsDisruptDisk(const SN_EVENT p_Events)
 *
 * @param   [IN]    p_Events                   Supernova-event bitmask generated during the current event
 * @return                                      true for any non-AIC supernova event; false otherwise
 */
bool SupernovaEventsDisruptDisk(const SN_EVENT p_Events) {
    return p_Events != SN_EVENT::NONE && p_Events != SN_EVENT::AIC;
}

/*
 * Determine whether a CE has the supported compact-object plus envelope-donor geometry
 *
 * bool IsSupportedCommonEnvelopeChannel(const STELLAR_TYPE p_StellarType1PreCE,
 *                                       const STELLAR_TYPE p_StellarType2PreCE,
 *                                       const double       p_Mass1PreCE,
 *                                       const double       p_Mass2PreCE,
 *                                       const bool         p_EnvelopeFlag1,
 *                                       const bool         p_EnvelopeFlag2)
 *
 * @param   [IN]    p_StellarType1PreCE       Stellar type of star 1 immediately before CE
 * @param   [IN]    p_StellarType2PreCE       Stellar type of star 2 immediately before CE
 * @param   [IN]    p_Mass1PreCE              Mass of star 1 immediately before CE (Msol)
 * @param   [IN]    p_Mass2PreCE              Mass of star 2 immediately before CE (Msol)
 * @param   [IN]    p_EnvelopeFlag1           Whether star 1 has an envelope immediately before CE
 * @param   [IN]    p_EnvelopeFlag2           Whether star 2 has an envelope immediately before CE
 * @return                                      Whether the CE channel is supported
 */
bool IsSupportedCommonEnvelopeChannel(const STELLAR_TYPE p_StellarType1PreCE,
                                      const STELLAR_TYPE p_StellarType2PreCE,
                                      const double       p_Mass1PreCE,
                                      const double       p_Mass2PreCE,
                                      const bool         p_EnvelopeFlag1,
                                      const bool         p_EnvelopeFlag2) {

    const STELLAR_TYPE_LIST supportedCompactObjects = {
        STELLAR_TYPE::HELIUM_WHITE_DWARF,
        STELLAR_TYPE::CARBON_OXYGEN_WHITE_DWARF,
        STELLAR_TYPE::OXYGEN_NEON_WHITE_DWARF,
        STELLAR_TYPE::NEUTRON_STAR,
        STELLAR_TYPE::BLACK_HOLE
    };
    const bool star1IsCompactObject = utils::IsOneOf(p_StellarType1PreCE, supportedCompactObjects);
    const bool star2IsCompactObject = utils::IsOneOf(p_StellarType2PreCE, supportedCompactObjects);

    if (!((star1IsCompactObject && p_EnvelopeFlag2 && !p_EnvelopeFlag1) ||
          (star2IsCompactObject && p_EnvelopeFlag1 && !p_EnvelopeFlag2))) {
        return false;
    }

    const double binaryMassPreCEMsol = p_Mass1PreCE + p_Mass2PreCE;
    if (utils::Compare(binaryMassPreCEMsol, 0.0) <= 0) {
        return false;
    }

    return true;
}

/*
 * Evaluate CBD formation criterion
 *
 * The mass supplied to the potential CBD is a fraction f_M of the decrease in total binary mass
 * across CE: (M1 + M2)_pre-CE - (M1 + M2)_post-CE. The resulting initial CBD
 * reservoir is assumed to satisfy m_CBD << M_binary during the evolution; formation does not
 * check this assumption.
 *
 * The criterion compares the total angular momentum assigned to the
 * CBD reservoir against the total angular momentum required for
 * supporting the CBD, given by the Keplerian AM at the adopted CBD inner radius times a disk structure factor.
 *
 * FormationResult EvaluateFormation(const FormationInput& p_Input)
 *
 * @param   [IN]    p_Input                    Physical state and model options for the formation criterion
 * @return                                      Formation outcome and supplied mass
 */
FormationResult EvaluateFormation(const FormationInput& p_Input) {
    FormationResult result{};

    const double massLostDuringCEMsol = p_Input.binaryMassPreCEMsol - p_Input.binaryMassPostCEMsol;
    if (utils::Compare(massLostDuringCEMsol, 0.0) <= 0 || utils::Compare(p_Input.semiMajorAxisPreRLOFRsol, 0.0) <= 0) {
        return result;
    }

    result.initialReservoirMassMsol = p_Input.options.envelopeMassFractionSupplied * massLostDuringCEMsol;
    if (utils::Compare(result.initialReservoirMassMsol, 0.0) <= 0) {
        return result;
    }

    const double totalAngularMomentumPreRLOF = p_Input.totalAngularMomentumPreRLOF;
    // beta = fJ/fM, so the angular momentum assigned to the supplied CBD
    // reservoir is J_CBD = beta * fM * J_available.
    const double suppliedDiskAngularMomentum = p_Input.options.circumbinaryDiskBeta *
                                               p_Input.options.envelopeMassFractionSupplied *
                                               totalAngularMomentumPreRLOF;

    // The adopted inner CBD radius scales with the post-CE separation because
    // the disk is assumed to orbit the post-CE binary, not the
    // larger pre-CE orbit.  Formation is allowed when the assigned angular
    // momentum meets or exceeds the threshold.
    const double diskInnerRadiusAU = p_Input.options.diskInnerRadiusOverSeparation * p_Input.semiMajorAxisPostCERsol * RSOL_TO_AU;
    const double minimumDiskAngularMomentum = result.initialReservoirMassMsol *
                                              p_Input.options.diskStructureFactor *
                                              std::sqrt(G_AU_Msol_yr *
                                                        p_Input.binaryMassPostCEMsol *
                                                        diskInnerRadiusAU);

    result.forms = utils::Compare(suppliedDiskAngularMomentum, minimumDiskAngularMomentum) >= 0;
    return result;
}

}  // namespace CircumbinaryDisk
