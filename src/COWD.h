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

    std::tuple<double,int> DetermineAccretionRegime(bool p_HeRich,
        const double p_AccretedMass,
        const double p_Dt);                                                                                                                                      // NRS: To check the current accretion regime and mass retention. Also activates flags for type change in some situations.
    void ResolveAccretionRegime(const int p_Regime, const double p_AccretedMass, const double p_Dt);

protected:

    void Initialise() {
        m_StellarType = STELLAR_TYPE::CARBON_OXYGEN_WHITE_DWARF;                                                                                                // Set stellar type
        CalculateTimescales();                                                                                                                                  // Initialise timescales
    m_HShell = 0.0; // Initialize hydrogen shell
    m_HeShell = 0.0; // Initialize helium shell
    m_DoubleDetonation = false;
    m_OffCenterIgnition = false;
    }


    double          CalculateLuminosityOnPhase(const double p_Mass,
                                               const double p_Time,
                                               const double p_Metallicity) const    { return CalculateLuminosityOnPhase_Static(p_Mass, p_Time, p_Metallicity); }
    double          CalculateLuminosityOnPhase() const                              { return CalculateLuminosityOnPhase(m_Mass, m_Age, m_Metallicity); }        // Use class member variables

    STELLAR_TYPE    EvolveToNextPhase(); // Modified to include evolution to ONe WD

    bool            IsSupernova() const                                             { return m_DoubleDetonation; }     // Going supernova if mass and He shell are large enough

    STELLAR_TYPE    ResolveSupernova()                                              { return GiantBranch::ResolveSupernova(); }                                 // Use GiantBranch, for now

    bool            ShouldEvolveOnPhase();                                                  // Evolve on phase unless mass > Chandrasekhar mass. Modified to include evolution to ONe WD

};

#endif // __COWD_h__
