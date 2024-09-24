#ifndef __ONeWD_h__
#define __ONeWD_h__

#include "constants.h"
#include "typedefs.h"

#include "WhiteDwarfs.h"


class BaseStar;
class WhiteDwarfs;

class ONeWD: virtual public BaseStar, public WhiteDwarfs {

public:

    ONeWD() { m_StellarType = STELLAR_TYPE::OXYGEN_NEON_WHITE_DWARF; };
    
    ONeWD(const BaseStar &p_BaseStar, const bool p_Initialise = true) : BaseStar(p_BaseStar), WhiteDwarfs(p_BaseStar) {
        m_StellarType = STELLAR_TYPE::OXYGEN_NEON_WHITE_DWARF;                                                                                                      // Set stellar type
        if (p_Initialise) Initialise();                                                                                                                             // Initialise if required
    }

    ONeWD* Clone(const OBJECT_PERSISTENCE p_Persistence, const bool p_Initialise = true) {
        ONeWD* clone = new ONeWD(*this, p_Initialise); 
        clone->SetPersistence(p_Persistence); 
        return clone; 
    }

    static ONeWD* Clone(ONeWD& p_Star, const OBJECT_PERSISTENCE p_Persistence, const bool p_Initialise = true) {
        ONeWD* clone = new ONeWD(p_Star, p_Initialise); 
        clone->SetPersistence(p_Persistence); 
        return clone; 
    }


    // member functions
    static double CalculateLuminosityOnPhase_Static(const double p_Mass, 
                                                    const double p_Time, 
                                                    const double p_Metallicity)     { return WhiteDwarfs::CalculateLuminosityOnPhase_Static(p_Mass, 
                                                                                                                                            p_Time, 
                                                                                                                                            p_Metallicity, 
                                                                                                                                            WD_Baryon_Number.at(STELLAR_TYPE::OXYGEN_NEON_WHITE_DWARF)); }
    

protected:

    void Initialise() {
        m_StellarType = STELLAR_TYPE::OXYGEN_NEON_WHITE_DWARF;                                                                                                      // Set stellar type
        CalculateTimescales();                                                                                                                                      // Initialise timescales
        m_Age             = 0.0;                                                                                                                                    // Set age appropriately
        m_HShell          = 0.0;                                                                                                                                    // Initialize hydrogen shell
        m_HeShell         = 0.0;                                                                                                                                    // Initialize helium shell
        m_AccretionRegime = ACCRETION_REGIME::ZERO;
    }


    // member functions

    double          CalculateHeliumAbundanceCoreOnPhase() const                     { return 0.0; };
    double          CalculateHeliumAbundanceSurfaceOnPhase() const                  { return 0.0; };
    
    double          CalculateHydrogenAbundanceCoreOnPhase() const                   { return 0.0; };
    double          CalculateHydrogenAbundanceSurfaceOnPhase() const                { return 0.0; };
    
    double          CalculateLuminosityOnPhase(const double p_Mass,
                                                const double p_Time,
                                                const double p_Metallicity) const   { return CalculateLuminosityOnPhase_Static(p_Mass, p_Time, p_Metallicity); }

    double          CalculateLuminosityOnPhase() const                              { return CalculateLuminosityOnPhase(m_Mass, m_Age, m_Metallicity); }    // Use class member variables

    STELLAR_TYPE    EvolveToNextPhase();
    bool            IsSupernova() const;                                             
    bool            ShouldEvolveOnPhase() const;   
                                                           
};

#endif // __ONeWD_h__
