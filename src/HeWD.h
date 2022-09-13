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

    //void ResolveAccretionRegime(const ACCRETION_REGIME p_Regime, const double p_DonorThermalMassLossRate);

protected:

    void Initialise() {
        m_StellarType = STELLAR_TYPE::HELIUM_WHITE_DWARF;                                                                                                           // Set stellar type
        CalculateTimescales();                                                                                                                                      // Initialise timescales
        // RTW: Is this the right place for these? Do we want them to reset if they switch from another type? (currently the case)
        // NRS: Probably the Ritter values can be initialized somewhere else, as they only depend on the initial choice of metallicity.
        // NRS: for the shells, I would argue that there is no option to come back directly to this stage after evolving.
        m_Age = 0.0;                                                                                                                                                // Set age appropriately
        m_HShell = 0.0; // Initialize Hydrogen Shell
        m_HeShell = 0.0; // Initialize Helium Shell
        m_l0Ritter = Calculatel0Ritter();
        m_XRitter = CalculateXRitter();
        m_lambdaRitter = CalculatelambdaRitter();
        m_SubChandrasekhar = false; // RTW: Should this be true at the start? NRS: Probably a confusing name, as this only flags whether there should be a SN Ia-like event at sub-Chandrasekhar mass value.
        m_Rejuvenate = false;
        m_AccretionRegime = ACCRETION_REGIME::NONE;
    }


    // member functions - alphabetically
    double  CalculateLambdaDewi() const                                                   { return BaseStar::CalculateLambdaDewi(); }                                     // Not supported - use BaseStar
    double  CalculateLambdaNanjingStarTrack(const double p_Mass, const double p_Metallicity) const { return BaseStar::CalculateLambdaNanjingStarTrack(0.0, 0.0); }                          // Not supported - use BaseStar (0.0 are dummy values)    JR: todo: check this (type 10 not mentioned as not supported in original code)
    double  CalculateLambdaNanjingEnhanced(const int p_MassInd, const int p_Zind) const   { return CalculateLambdaNanjingStarTrack(0.0, 0.0); }                            // 0.0 are dummy values that are not used
    double  CalculateLuminosityOnPhase(const double p_Mass,
                                       const double p_Time,
                                       const double p_Metallicity) const                  { return CalculateLuminosityOnPhase_Static(p_Mass, p_Time, p_Metallicity); }
    double  CalculateLuminosityOnPhase() const                                            { return CalculateLuminosityOnPhase_Static(m_Mass, m_Age, m_Metallicity); }     // Use class member variables

    DBL_DBL CalculateMassAcceptanceRate(const double p_DonorMassRate,
                                        const bool   p_IsHeRich);  

    double  CalculateRadiusOnPhase(const double p_Mass) const                             { return CalculateRadiusOnPhase_Static(p_Mass); }
    double  CalculateRadiusOnPhase() const                                                { return CalculateRadiusOnPhase(m_Mass); }                                      // Use class member variables
    std::tuple <double, STELLAR_TYPE> CalculateRadiusAndStellarTypeOnPhase() const        { return BaseStar::CalculateRadiusAndStellarTypeOnPhase(); }

    STELLAR_TYPE    EvolveToNextPhase(); // Allow evolution, either SN or Rejuvenation

    // RTW: Might want to rename to m_IsSubChandrasekhar for clarity. Also, should this be the opposite, as in, go SN if !m_SubChandrasekhar ?
    // NRS: I agree with renaming, that is what causes confusion.
    // RTW: If the SN is AIC, you can copy the code in ONeWD.h
    bool    IsSupernova() const                                                           { return m_SubChandrasekhar; }     // Going supernova if mass and He shell are large enough

    bool    ShouldEvolveOnPhase(); // Modified so now it includes rejuvenation and Sub-Ch SN Ia
};

#endif // __HeWD_h__
