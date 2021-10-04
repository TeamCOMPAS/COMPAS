#ifndef __WhiteDwarfs_h__
#define __WhiteDwarfs_h__

#include "constants.h"
#include "typedefs.h"
#include "profiling.h"
#include "utils.h"

#include "Remnants.h"


class BaseStar;
class Remnants;

class WhiteDwarfs: virtual public BaseStar, public Remnants {

public:

    WhiteDwarfs(const BaseStar &p_BaseStar, const bool p_Initialise = true) : BaseStar(p_BaseStar), Remnants(p_BaseStar, false) {
        if (p_Initialise) Initialise();
    }

    WhiteDwarfs& operator = (const BaseStar &p_BaseStar) { static_cast<BaseStar&>(*this) = p_BaseStar; return *this; }


    // member functions
    static  double      CalculateLuminosityOnPhase_Static(const double p_Mass, 
                                                         const double p_Time, 
                                                         const double p_Metallicity, 
                                                         const double p_BaryonNumber);

    static  double      CalculateRadiusOnPhase_Static(const double p_Mass);


protected:


    // member functions - alphabetically
            double      CalculateCOCoreMassOnPhase() const                  { return m_COCoreMass; }                                                // NO-OP

            double      CalculateHeCoreMassOnPhase() const                  { return m_HeCoreMass; }                                                // NO-OP

            double      CalculateInitialSupernovaMass() const               { return 0.0; }

            double      CalculateRadiusOnPhase(const double p_Mass) const   { return CalculateRadiusOnPhase_Static(p_Mass); }
            double      CalculateRadiusOnPhase() const                      { return CalculateRadiusOnPhase(m_Mass); }                              // Use class member variables

            ENVELOPE    DetermineEnvelopeType() const                       { return ENVELOPE::CONVECTIVE; }                                        // Always CONVECTIVE

};

#endif // __WhiteDwarfs_h__
