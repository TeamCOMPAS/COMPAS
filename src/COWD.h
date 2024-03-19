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

    static double CalculateLuminosityOnPhase_Static(const double p_Mass, 
                                                    const double p_Time, 
                                                    const double p_Metallicity)     { return WhiteDwarfs::CalculateLuminosityOnPhase_Static(p_Mass, 
                                                                                                                                            p_Time, 
                                                                                                                                            p_Metallicity, 
                                                                                                                                            WD_Baryon_Number.at(STELLAR_TYPE::CARBON_OXYGEN_WHITE_DWARF)); }

    ACCRETION_REGIME DetermineAccretionRegime(const bool   p_HeRich,
                                              const double p_DonorThermalMassLossRate);                                                                         // To check the current accretion regime and mass retention. Can also change flags related to SN events.

protected:

    void Initialise() {
        m_StellarType       = STELLAR_TYPE::CARBON_OXYGEN_WHITE_DWARF;                                                                                          // Set stellar type
        CalculateTimescales();                                                                                                                                  // Initialise timescales
        m_HShell            = 0.0;                                                                                                                              // Initialize hydrogen shell
        m_HeShell           = 0.0;                                                                                                                              // Initialize helium shell
        m_HeShellDetonation = false;
        m_OffCenterIgnition = false;
        m_AccretionRegime   = ACCRETION_REGIME::NONE;
    }


    double          CalculateLuminosityOnPhase(const double p_Mass,
                                               const double p_Time,
                                               const double p_Metallicity) const    { return CalculateLuminosityOnPhase_Static(p_Mass, p_Time, p_Metallicity); }
    double          CalculateLuminosityOnPhase() const                              { return CalculateLuminosityOnPhase(m_Mass, m_Age, m_Metallicity); }        // Use class member variables

    DBL_DBL         CalculateMassAcceptanceRate(const double p_DonorMassRate,
                                                const bool   p_IsHeRich);  
    DBL_DBL         CalculateMassAcceptanceRate(const double p_DonorMassRate,
                                                const double p_AccretorMassRate,
                                                const bool   p_IsHeRich)            { return CalculateMassAcceptanceRate(p_DonorMassRate, p_IsHeRich); }        // Ignore the input accretion rate for WDs

    STELLAR_TYPE    EvolveToNextPhase();
    bool            IsSupernova() const                                             { return m_HeShellDetonation || IsMassAboveChandrasekhar(); };                                             
    bool            ShouldEvolveOnPhase() const                                     { return m_OffCenterIgnition ? false : !IsSupernova(); };                   // From https://ui.adsabs.harvard.edu/abs/2017MNRAS.472.1593W/abstract around the end of section 3.2. Also, allows SN.                                               

};

#endif // __COWD_h__
