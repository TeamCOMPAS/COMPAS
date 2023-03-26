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
        m_HShell = 0.0; // Initialize hydrogen shell
        m_HeShell = 0.0; // Initialize helium shell
        m_AccretionRegime = ACCRETION_REGIME::NONE;
    }


    // member functions

            double          CalculateLuminosityOnPhase(const double p_Mass,
                                                       const double p_Time,
                                                       const double p_Metallicity) const    { return CalculateLuminosityOnPhase_Static(p_Mass, p_Time, p_Metallicity); }

            double          CalculateLuminosityOnPhase() const                              { return CalculateLuminosityOnPhase(m_Mass, m_Age, m_Metallicity); }        // Use class member variables

            STELLAR_TYPE    EvolveToNextPhase();
            bool            IsSupernova() const;                                             
            bool            ShouldEvolveOnPhase();                                                  
};

#endif // __ONeWD_h__
