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
                                       const bool p_HeRich);

protected:
    // member variables

            bool                    m_DoubleDetonation;                         // Flag
            double                  m_HeShell;                                  // Current WD He-shell size (Msol). Increases through accretion.
            double                  m_HShell;                                   // Current WD H-shell size (Msol). Increases through accretion.
            double                  m_l0Ritter;                                 // Parameter from numerical calculations, see Ritter 1999, section 3. Eqs 10 and 12, as well as table 2. Corresponds to L0.
            double                  m_lambdaRitter;                             // Parameter from numerical calculations, see Ritter 1999, section 3. Eqs 10 and 12, as well as table 2.
            bool                    m_OffCenterIgnition;                        // Flag for CO WD evolution into ONe WD
            bool                    m_Rejuvenate;                               // Flag for evolution of HeWD back into HeMS
            bool                    m_SubChandrasekhar;                         // Flag for SubCh SN of HeWD
            double                  m_XRitter;                                  // Assumed hydrogen-mass fraction of material being accreted by He WD, as in Ritter 1999, table 2.

    // member functions - alphabetically
            double      CalculateCOCoreMassOnPhase() const                  { return m_COCoreMass; }                                                // NO-OP

            double      CalculateHeCoreMassOnPhase() const                  { return m_HeCoreMass; }                                                // NO-OP

            double      CalculateetaH(const double p_LogMassRate);

            double      CalculateetaHe(const double p_LogMassRate);

            double      CalculateetaPTY(const double p_MassRate);

            double      Calculatel0() const                                  {return (m_Metallicity > 0.01) ? 1995262.3 : 31622.8; } // Luminosity constant which depends on metallicity in Ritter 1999, eq 10

            double      CalculateX() const                                   {return (m_Metallicity > 0.01) ? 0.7 : 0.8 ; } // Assumed Hydrogen-mass fraction

            double      Calculatelambda() const                              {return (m_Metallicity > 0.01) ? 8 : 5 ;  } // Exponent for the assumed core-mass and luminosity relationship in Ritter 1999

            double      CalculateInitialSupernovaMass() const               { return 0.0; }

            DBL_DBL     CalculateWDMassAcceptanceRate(const double p_DonorMassRate,
                                                      const bool   p_IsHeRich);

            double      CalculateRadiusOnPhase(const double p_Mass) const   { return CalculateRadiusOnPhase_Static(p_Mass); }
            double      CalculateRadiusOnPhase() const                      { return CalculateRadiusOnPhase(m_Mass); }                              // Use class member variables

            ENVELOPE    DetermineEnvelopeType() const                       { return ENVELOPE::CONVECTIVE; }                                        // Always CONVECTIVE

};

#endif // __WhiteDwarfs_h__
