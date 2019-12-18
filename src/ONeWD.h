#ifndef __ONeWD_h__
#define __ONeWD_h__

#include "constants.h"
#include "typedefs.h"

#include "COWD.h"


class BaseStar;
class COWD;

class ONeWD: virtual public BaseStar, public COWD {

public:

    ONeWD(const BaseStar &p_BaseStar, const bool p_Initialise = true) : BaseStar(p_BaseStar), COWD(p_BaseStar, false) {
        if (p_Initialise) Initialise();
    }

    ONeWD& operator = (const BaseStar &p_BaseStar) {
        static_cast<BaseStar&>(*this) = p_BaseStar;
        Initialise();
        return *this;
    }


protected:

    void Initialise() {
        m_StellarType = STELLAR_TYPE::OXYGEN_NEON_WHITE_DWARF;                                                                                                  // Set stellar type
        CalculateTimescales();                                                                                                                                  // Initialise timescales
        m_Age = 0.0;                                                                                                                                            // Set age appropriately
    }


    // member functions

           double       CalculateInitialSupernovaMass()                         { return 5.0; }                                                                 // Force ONeWD to ccSN, 5.0 doesn't change a physical parameter of the star.

           double       CalculateLuminosityOnPhase(const double p_Mass,
                                                   const double p_Time,
                                                   const double p_Metallicity)  { return CalculateLuminosityOnPhase_Static(p_Mass, p_Time, p_Metallicity); }
           double       CalculateLuminosityOnPhase()                            { return CalculateLuminosityOnPhase_Static(m_Mass, m_Age, m_Metallicity); }     // Use class member variables
    static double       CalculateLuminosityOnPhase_Static(const double p_Mass, const double p_Time, const double p_Metallicity);

           STELLAR_TYPE EvolveToNextPhase()                                     { return BaseStar::EvolveToNextPhase(); }                                       // Default to BaseStar

           bool         IsSupernova()                                           { return (utils::Compare(m_Mass, MCH) > 0); }                                   // Going supernova if mass large enough

           STELLAR_TYPE ResolveSupernova()                                      { return GiantBranch::ResolveSupernova(); }                                     // Use GiantBranch

           bool         ShouldEvolveOnPhase()                                   { return (utils::Compare(m_Mass, MCH) <= 0); }                                  // Evolve on phase unless mass > Chandrasekhar mass
           bool         ShouldSkipPhase()                                       { return false; }                                                               // Never skip HeMS phase

};

#endif // __ONeWD_h__
