#ifndef __stellarUtils_h__
#define __stellarUtils_h__

#include "constants.h"
#include "typedefs.h"
#include "ErrorCatalog.h"

#include "BaseStar.h"
#include "MS_lte_07.h"
#include "MS_gt_07.h"
#include "CH.h"
#include "HG.h"
#include "FGB.h"
#include "CHeB.h"
#include "EAGB.h"
#include "TPAGB.h"
#include "HeMS.h"
#include "HeHG.h"
#include "HeGB.h"
#include "HeWD.h"
#include "COWD.h"
#include "ONeWD.h"
#include "NS.h"
#include "BH.h"
#include "MR.h"


namespace stellarUtils {


/*
 * Instantiates a star of the requested stellar type, cast to BaseStar*.
 *
 * This is not a cloning function - a new, empty, star of the requested stellar type is
 * instantiated.  This is primarily so functions pertaining to the stellar type can be
 * called when only the stellar type is known, and there is no object of that stellar type
 * available.
 * 
 * Currently the new star is not initialised via the Initialise() function of the relevant
 * class - some ov the Initialise() functions depend on class member variables carried over
 * from a prior stellar type, so until we resolve that, this star cannot be initialised in
 * that way.
 * 
 * This function returns a pointer.  It is the responsibility of the caller to return the
 * memory allocation for the new star after use.  e.g.:
 * 
 *     BaseStar* CHeBstar       = stellarUtils::NewStar(STELLAR_TYPE::CORE_HELIUM_BURNING);    // create new CHeB star
 *     MT_CASE massTransferCase = CHeBstar->DetermineMassTransferTypeAsDonor();                // get MT type as donor of CHeB star
 *     delete CHeBstar; CHeBstar = nullptr;                                                    // free memory allocated for CHeBstar
 *
 *
 * BaseStar* NewStar(const STELLAR_TYPE p_StellarType)
 *
 * @param   [IN]    p_StellarType               Stellar type of the new star
 * @return                                      Newly instantiated star, cast to BaseStar*
 */
    BaseStar* NewStar(const STELLAR_TYPE p_StellarType) {

        BaseStar *ptr;

        switch (p_StellarType) {
            case STELLAR_TYPE::MS_LTE_07                                : {ptr = new MS_lte_07();} break;
            case STELLAR_TYPE::MS_GT_07                                 : {ptr = new MS_gt_07();} break;
            case STELLAR_TYPE::CHEMICALLY_HOMOGENEOUS                   : {ptr = new CH();} break;
            case STELLAR_TYPE::HERTZSPRUNG_GAP                          : {ptr = new HG();} break;
            case STELLAR_TYPE::FIRST_GIANT_BRANCH                       : {ptr = new FGB();} break;
            case STELLAR_TYPE::CORE_HELIUM_BURNING                      : {ptr = new CHeB();} break;
            case STELLAR_TYPE::EARLY_ASYMPTOTIC_GIANT_BRANCH            : {ptr = new EAGB();} break;
            case STELLAR_TYPE::THERMALLY_PULSING_ASYMPTOTIC_GIANT_BRANCH: {ptr = new TPAGB();} break;
            case STELLAR_TYPE::NAKED_HELIUM_STAR_MS                     : {ptr = new HeMS();} break;
            case STELLAR_TYPE::NAKED_HELIUM_STAR_HERTZSPRUNG_GAP        : {ptr = new HeHG();} break;
            case STELLAR_TYPE::NAKED_HELIUM_STAR_GIANT_BRANCH           : {ptr = new HeGB();} break;
            case STELLAR_TYPE::HELIUM_WHITE_DWARF                       : {ptr = new HeWD();} break;
            case STELLAR_TYPE::CARBON_OXYGEN_WHITE_DWARF                : {ptr = new COWD();} break;
            case STELLAR_TYPE::OXYGEN_NEON_WHITE_DWARF                  : {ptr = new ONeWD();} break;
            case STELLAR_TYPE::NEUTRON_STAR                             : {ptr = new NS();} break;
            case STELLAR_TYPE::BLACK_HOLE                               : {ptr = new BH();} break;
            case STELLAR_TYPE::MASSLESS_REMNANT                         : {ptr = new MR();} break;

            default:                                                                    // unknown stellar type
                // the only ways this can happen are if someone added a STELLAR_TYPE
                // and it isn't accounted for in this code, or if there is a defect in the code that causes
                // this function to be called with a bad parameter.  We should not default here, with or without
                // a warning.
                // We are here because the user chose a prescription this code doesn't account for, and that should
                // be flagged as an error and result in termination of the evolution of the star or binary.
                // The correct fix for this is to add code for the missing prescription or, if the missing
                // prescription is superfluous, remove it from the option.
                THROW_ERROR_STATIC(ERROR::UNKNOWN_STELLAR_TYPE);                        // throw error
        }

        return ptr;
    }
}

#endif // __stellarUtils_h__
