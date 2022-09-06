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

    ACCRETION_REGIME DetermineAccretionRegime(const bool p_HeRich,
                                              const double p_DonorThermalMassLossRate);                                                                                                                                      // To check the current accretion regime and mass retention. Also activates flags for type change in some situations.

    //void ResolveAccretionRegime(const ACCRETION_REGIME p_Regime, const double p_DonorThermalMassLossRate);

protected:

    void Initialise() {
        m_StellarType = STELLAR_TYPE::CARBON_OXYGEN_WHITE_DWARF;                                                                                                // Set stellar type
        CalculateTimescales();                                                                                                                                  // Initialise timescales
        // RTW is this the right place for these? - Do we want to reset them if switch from another type?
        m_HShell = 0.0; // Initialize hydrogen shell
        m_HeShell = 0.0; // Initialize helium shell
        m_DoubleDetonation = false;
        m_OffCenterIgnition = false;
    }


    double          CalculateLuminosityOnPhase(const double p_Mass,
                                               const double p_Time,
                                               const double p_Metallicity) const    { return CalculateLuminosityOnPhase_Static(p_Mass, p_Time, p_Metallicity); }
    double          CalculateLuminosityOnPhase() const                              { return CalculateLuminosityOnPhase(m_Mass, m_Age, m_Metallicity); }        // Use class member variables

    STELLAR_TYPE    EvolveToNextPhase();

    bool            IsSupernova() const                                             { return m_DoubleDetonation; }     // Going supernova if mass and He shell are large enough

    // RTW why are we using the giantbranch version?
    STELLAR_TYPE    ResolveSupernova()                                              { return GiantBranch::ResolveSupernova(); }                                 // Use GiantBranch, for now

    bool            ShouldEvolveOnPhase();                                                  // Evolve on phase unless mass > Chandrasekhar mass.

};

#endif // __COWD_h__
