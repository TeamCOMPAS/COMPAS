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


    void                IncrementShell(const double p_AccretedMass,
                                       bool p_HeRich);

protected:


    // member functions - alphabetically
            double      CalculateCOCoreMassOnPhase() const                  { return m_COCoreMass; }                                                // NO-OP

            double      CalculateHeCoreMassOnPhase() const                  { return m_HeCoreMass; }                                                // NO-OP

            double      CalculateetaH(const double p_LogMassRate);

            double      CalculateetaHe(const double p_LogMassRate);

            double      CalculateetaPTY(const double p_MassRate);

            double      Calculatel0() const                                  {return (m_Metallicity > 0.01) ? 1995262.3 : 31622.8; }

            double      CalculateX() const                                   {return (m_Metallicity > 0.01) ? 0.7 : 0.8 ; }

            double      Calculatelambda() const                              {return (m_Metallicity > 0.01) ? 8 : 5 ;  }

            double      CalculateInitialSupernovaMass() const               { return 0.0; }

            DBL_DBL     CalculateWDMassAcceptanceRate(const double p_DonorMassRate,
                                                      const bool   p_IsHeRich);

            double      CalculateRadiusOnPhase(const double p_Mass) const   { return CalculateRadiusOnPhase_Static(p_Mass); }
            double      CalculateRadiusOnPhase() const                      { return CalculateRadiusOnPhase(m_Mass); }                              // Use class member variables

            ENVELOPE    DetermineEnvelopeType() const                       { return ENVELOPE::CONVECTIVE; }                                        // Always CONVECTIVE

};

#endif // __WhiteDwarfs_h__
