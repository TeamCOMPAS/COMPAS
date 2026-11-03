#include "CircumbinaryDisk.h"

// Stub implementation used when COMPAS is built without spindler-c support.

#include <stdexcept>

namespace CircumbinaryDisk {

/*
 * Reject CBD evolution when COMPAS was built without the spindler-c backend
 *
 * EvolutionTrack EvolveTrack(const EvolutionInput&)
 *
 * @return                                      This function does not return; it throws a runtime error
 */
EvolutionTrack EvolveTrack(const EvolutionInput&) {
    throw std::runtime_error(
        "Circumbinary-disk evolution was requested, but this COMPAS executable "
        "was built without spindler-c support. Rebuild with SPINDLER_C_ROOT "
        "set to the spindler-c build tree."
    );
}

}  // namespace CircumbinaryDisk
