#ifndef __HeWD_h__
#define __HeWD_h__

#include "constants.h"
#include "typedefs.h"
#include "profiling.h"
#include "utils.h"

#include "WhiteDwarfs.h"


class BaseStar;
class WhiteDwarfs;

class HeWD: virtual public BaseStar, public WhiteDwarfs {

public:

    HeWD(const BaseStar &p_BaseStar, const bool p_Initialise = true) : BaseStar(p_BaseStar), WhiteDwarfs(p_BaseStar, false) {
        if (p_Initialise) Initialise();
    }

    HeWD& operator = (const BaseStar &p_BaseStar) {
        static_cast<BaseStar&>(*this) = p_BaseStar;
        Initialise();
        return *this;
    }


    // member functions
    static  double      CalculateLuminosityOnPhase_Static(const double p_Mass, 
                                                          const double p_Time, 
                                                          const double p_Metallicity)   { return WhiteDwarfs::CalculateLuminosityOnPhase_Static(p_Mass, 
                                                                                                                                                p_Time, 
                                                                                                                                                p_Metallicity, 
                                                                                                                                                WD_Baryon_Number.at(STELLAR_TYPE::HELIUM_WHITE_DWARF)); }
    ACCRETION_REGIME DetermineAccretionRegime(const bool p_HeRich,
                                              const double p_DonorThermalMassLossRate);

protected:

    void Initialise() {
        m_StellarType              = STELLAR_TYPE::HELIUM_WHITE_DWARF;  // Set stellar type
        CalculateTimescales();                                          // Initialise timescales
        m_Age                      = 0.0;                               // Set age appropriately
        m_HShell                   = 0.0;                               // Initialize Hydrogen Shell
        m_HeShell                  = 0.0;                               // Initialize Helium Shell
        m_l0Ritter                 = Calculatel0Ritter();
        m_XRitter                  = CalculateXRitter();
        m_lambdaRitter             = CalculatelambdaRitter();
        m_IsSubChandrasekharTypeIa = false; 
        m_ShouldRejuvenate         = false;
        m_AccretionRegime          = ACCRETION_REGIME::NONE;
    }


    // member functions - alphabetically
    double          CalculateLambdaDewi() const                                                             { return BaseStar::CalculateLambdaDewi(); }                                     // Not supported - use BaseStar
    double          CalculateLambdaNanjingStarTrack(const double p_Mass, const double p_Metallicity) const  { return BaseStar::CalculateLambdaNanjingStarTrack(0.0, 0.0); }                 // Not supported - use BaseStar (0.0 are dummy values)    JR: todo: check this (type 10 not mentioned as not supported in original code)
    double          CalculateLambdaNanjingEnhanced(const int p_MassInd, const int p_Zind) const             { return CalculateLambdaNanjingStarTrack(0.0, 0.0); }                           // 0.0 are dummy values that are not used
    double          CalculateLuminosityOnPhase(const double p_Mass,
                                               const double p_Time,
                                               const double p_Metallicity) const                            { return CalculateLuminosityOnPhase_Static(p_Mass, p_Time, p_Metallicity); }
    double          CalculateLuminosityOnPhase() const                                                      { return CalculateLuminosityOnPhase_Static(m_Mass, m_Age, m_Metallicity); }     // Use class member variables

    DBL_DBL         CalculateMassAcceptanceRate(const double p_DonorMassRate,
                                                const bool   p_IsHeRich);  
    DBL_DBL         CalculateMassAcceptanceRate(const double p_DonorMassRate,
                                                const double p_AccretorMassRate,
                                                const bool   p_IsHeRich)                                    { return CalculateMassAcceptanceRate(p_DonorMassRate, p_IsHeRich); }            // Ignore the input accretion rate for WDs

    double          CalculateRadiusOnPhase(const double p_Mass) const                                       { return CalculateRadiusOnPhase_Static(p_Mass); }
    double          CalculateRadiusOnPhase() const                                                          { return CalculateRadiusOnPhase(m_Mass); }                                      // Use class member variables
    std::tuple <double, STELLAR_TYPE> CalculateRadiusAndStellarTypeOnPhase() const                          { return BaseStar::CalculateRadiusAndStellarTypeOnPhase(); }

    STELLAR_TYPE    EvolveToNextPhase();                                                                                                                                                    // Allow evolution, either SN or Rejuvenation
    bool            IsSupernova() const                                                                     { return m_IsSubChandrasekharTypeIa; }                                          // Go supernova if mass and He shell are large enough;                                                           
    bool            ShouldEvolveOnPhase() const                                                             { return m_ShouldRejuvenate ? false : !IsSupernova(); }              

};

#endif // __HeWD_h__
