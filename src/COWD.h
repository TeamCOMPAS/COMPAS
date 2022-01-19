#ifndef __COWD_h__
#define __COWD_h__

#include "constants.h"
#include "typedefs.h"

#include "WhiteDwarfs.h"


class BaseStar;
class WhiteDwarfs;

class COWD: virtual public BaseStar, public WhiteDwarfs {

public:

    COWD(const BaseStar &p_BaseStar, const bool p_Initialise = true) : BaseStar(p_BaseStar), WhiteDwarfs(p_BaseStar, false) {
        if (p_Initialise) Initialise();
    }

    COWD& operator = (const BaseStar &p_BaseStar) {
        static_cast<BaseStar&>(*this) = p_BaseStar;
        Initialise();
        return *this;
    }


    // member functions
           void   CalculateAngularMomentum() const { }                                                                                                          // NO-OP

    static double CalculateLuminosityOnPhase_Static(const double p_Mass, 
                                                    const double p_Time, 
                                                    const double p_Metallicity)     { return WhiteDwarfs::CalculateLuminosityOnPhase_Static(p_Mass, 
                                                                                                                                            p_Time, 
                                                                                                                                            p_Metallicity, 
                                                                                                                                            WD_Baryon_Number.at(STELLAR_TYPE::CARBON_OXYGEN_WHITE_DWARF)); }


protected:

    void Initialise() {
        m_StellarType = STELLAR_TYPE::CARBON_OXYGEN_WHITE_DWARF;                                                                                                // Set stellar type
        CalculateTimescales();                                                                                                                                  // Initialise timescales
    }

    void FastForward() {                                                                                                                                                        // Set stellar attributes for stars initialized to this stellar type

        m_Radius                                   = CalculateRadiusOnPhase();
        m_Luminosity                               = CalculateLuminosityOnPhase();
    
        m_InitialLuminosity                        = m_Luminosity;
        m_InitialRadius                            = m_Radius;
        m_InitialStellarType                       = m_StellarType;
        m_StellarTypePrev                          = m_StellarType;
    }


    double          CalculateLuminosityOnPhase(const double p_Mass,
                                               const double p_Time,
                                               const double p_Metallicity) const    { return CalculateLuminosityOnPhase_Static(p_Mass, p_Time, p_Metallicity); }
    double          CalculateLuminosityOnPhase() const                              { return CalculateLuminosityOnPhase(m_Mass, m_Age, m_Metallicity); }        // Use class member variables

    STELLAR_TYPE    EvolveToNextPhase()                                             { m_Mass = m_Radius = m_Luminosity = m_Age = 0.0; return STELLAR_TYPE::MASSLESS_REMNANT; }

    bool            ShouldEvolveOnPhase() const                                     { return (m_Mass <= MCH); }                                                 // Evolve on phase unless mass > Chandrasekhar mass

};

#endif // __COWD_h__
