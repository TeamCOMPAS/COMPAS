#ifndef __ONeWD_h__
#define __ONeWD_h__

#include "constants.h"
#include "typedefs.h"

#include "WhiteDwarfs.h"


class BaseStar;
class WhiteDwarfs;

class ONeWD: virtual public BaseStar, public WhiteDwarfs {

public:

    ONeWD(const BaseStar &p_BaseStar, const bool p_Initialise = true) : BaseStar(p_BaseStar), WhiteDwarfs(p_BaseStar, false) {
        if (p_Initialise) Initialise();
    }

    ONeWD& operator = (const BaseStar &p_BaseStar) {
        static_cast<BaseStar&>(*this) = p_BaseStar;
        Initialise();
        return *this;
    }

    // member functions
    static double       CalculateLuminosityOnPhase_Static(const double p_Mass, 
                                                          const double p_Time, 
                                                          const double p_Metallicity)   { return WhiteDwarfs::CalculateLuminosityOnPhase_Static(p_Mass, 
                                                                                                                                        p_Time, 
                                                                                                                                        p_Metallicity, 
                                                                                                                                        WD_Baryon_Number.at(STELLAR_TYPE::OXYGEN_NEON_WHITE_DWARF)); }

protected:

    void Initialise() {
        m_StellarType = STELLAR_TYPE::OXYGEN_NEON_WHITE_DWARF;                                                                                                      // Set stellar type
        CalculateTimescales();                                                                                                                                      // Initialise timescales
        m_Age = 0.0;                                                                                                                                                // Set age appropriately
    }

    void FastForward() {                                                                                                                                                        // Set stellar attributes for stars initialized to this stellar type

        m_Radius                                   = CalculateRadiusOnPhase();
        m_Luminosity                               = CalculateLuminosityOnPhase();
    
        m_InitialLuminosity                        = m_Luminosity;
        m_InitialRadius                            = m_Radius;
        m_InitialStellarType                       = m_StellarType;
        m_StellarTypePrev                          = m_StellarType;
    }


    // member functions

            double       CalculateLuminosityOnPhase(const double p_Mass,
                                                   const double p_Time,
                                                   const double p_Metallicity) const    { return CalculateLuminosityOnPhase_Static(p_Mass, p_Time, p_Metallicity); }

            double       CalculateLuminosityOnPhase() const                              { return CalculateLuminosityOnPhase(m_Mass, m_Age, m_Metallicity); }        // Use class member variables
            bool         IsSupernova() const                                             { return (utils::Compare(m_Mass, MECS) > 0); }                              // Going supernova if mass large enough

            STELLAR_TYPE ResolveAIC();  

            STELLAR_TYPE ResolveSupernova()                                              { return ResolveAIC(); }                                                   // WDs should not SN, but they can collapse due to accretion

            bool         ShouldEvolveOnPhase() const                                     { return (utils::Compare(m_Mass, MECS) <= 0); }                             // Evolve on phase unless mass > ECSN threshold mass

};

#endif // __ONeWD_h__
