#ifndef __COWD_h__
#define __COWD_h__

#include "constants.h"
#include "typedefs.h"

#include "HeWD.h"


class BaseStar;
class HeWD;

class COWD: virtual public BaseStar, public HeWD {

public:

    COWD(const BaseStar &p_BaseStar, const bool p_Initialise = true) : BaseStar(p_BaseStar), HeWD(p_BaseStar, false) {
        if (p_Initialise) Initialise();
    }

    COWD& operator = (const BaseStar &p_BaseStar) {
        static_cast<BaseStar&>(*this) = p_BaseStar;
        Initialise();
        return *this;
    }


    // member functions
           void   CalculateAngularMomentum() { }                        // NO-OP
    static double CalculateLuminosityOnPhase_Static(const double p_Mass, const double p_Time, const double p_Metallicity);


protected:

    void Initialise() {
        m_StellarType = STELLAR_TYPE::CARBON_OXYGEN_WHITE_DWARF;                                                                                            // Set stellar type
        CalculateTimescales();                                                                                                                              // Initialise timescales
    }


    double          CalculateLuminosityOnPhase(const double p_Mass,
                                               const double p_Time,
                                               const double p_Metallicity)  { return CalculateLuminosityOnPhase_Static(p_Mass, p_Time, p_Metallicity); }
    double          CalculateLuminosityOnPhase()                            { return CalculateLuminosityOnPhase_Static(m_Mass, m_Age, m_Metallicity); }     // Use class member variables

    STELLAR_TYPE    EvolveToNextPhase();

    bool            ShouldEvolveOnPhase()                                   { return (m_Mass <= MCH); }                                                     // Evolve on phase unless mass > Chandrasekhar mass
    bool            ShouldSkipPhase()                                       { return false; }                                                               // Never skip HeMS phase

};

#endif // __COWD_h__
