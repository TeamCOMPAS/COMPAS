#ifndef __HG_h__
#define __HG_h__

#include "constants.h"
#include "typedefs.h"
#include "profiling.h"
#include "utils.h"
#include <boost/math/tools/roots.hpp>

#include "GiantBranch.h"


// JR: todo: revisit this one day - sometimes HG works better as GiantBranch, sometimes not...
// Right now it is GiantBranch - figure out which is more appropriate

class BaseStar;
class GiantBranch;

class HG: virtual public BaseStar, public GiantBranch {

public:

    HG(const BaseStar &p_BaseStar, const bool p_Initialise = true) : BaseStar(p_BaseStar), GiantBranch(p_BaseStar) {
        if (p_Initialise) Initialise();
    }

    HG& operator = (const BaseStar &p_BaseStar) {
        static_cast<BaseStar&>(*this) = p_BaseStar;
        Initialise();
        return *this;
    }


protected:

    void Initialise() {
        m_StellarType = STELLAR_TYPE::HERTZSPRUNG_GAP;                                                                                                                          // Set stellar type
        m_Tau = 0.0;                                                                                                                                                            // Start of phase
        CalculateTimescales();                                                                                                                                                  // Initialise timescales
        m_Age = m_Timescales[static_cast<int>(TIMESCALE::tMS)];                                                                                                                 // Set age appropriately
        
        //Update stellar properties at start of HG phase (since core defintion changes)
        CalculateGBParams();
        
        //update effective "initial" mass m_Mass0 so that the core mass is at least equal to the minimum core mass (only relevant if RetainCoreMassDuringCaseAMassTransfer() ) but no more than total mass
        if(utils::Compare(CalculateCoreMassOnPhase(m_Mass0, m_Age), std::min(m_Mass, MinimumCoreMass())) < 0) {
            m_Mass0 = Mass0ToMatchDesiredCoreMass(this, std::min(m_Mass,MinimumCoreMass()));
            CalculateTimescales();
            m_Age = m_Timescales[static_cast<int>(TIMESCALE::tMS)];
            CalculateGBParams();
        }
        m_CoreMass    = CalculateCoreMassOnPhase();
        m_COCoreMass  = CalculateCOCoreMassOnPhase();
        m_HeCoreMass  = CalculateHeCoreMassOnPhase();
        m_Luminosity  = CalculateLuminosityOnPhase();
        std::tie(m_Radius, std::ignore) = CalculateRadiusAndStellarTypeOnPhase();   // Update radius
    }


    // member functions - alphabetically
    double          CalculateCOCoreMassAtPhaseEnd() const                           { return 0.0; }                                                                             // McCO(HG) = 0.0
    double          CalculateCOCoreMassOnPhase() const                              { return 0.0; }                                                                             // McCO(HG) = 0.0

    double          CalculateCoreMassAt2ndDredgeUp(const DBL_VECTOR &p_GBParams)    { return p_GBParams[static_cast<int>(GBP::McDU)]; }                                         // NO-OP
    double          CalculateCoreMassAtPhaseEnd(const double p_Mass) const;
    double          CalculateCoreMassAtPhaseEnd() const                             { return CalculateCoreMassAtPhaseEnd(m_Mass0); }                                            // Use class member variables
    double          CalculateCoreMassOnPhase(const double p_Mass, const double p_Time) const;
    double          CalculateCoreMassOnPhase() const                                { return CalculateCoreMassOnPhase(m_Mass0, m_Age); }                                        // Use class member variables

    double          CalculateCriticalMassRatio(const bool p_AccretorIsDegenerate) const;

    double          CalculateGyrationRadius() const                                 { return 0.21; }                                                                            // Hurley et al., 2000, after eq 109 for n=3/2 polytrope or dense convective core. Single number approximation.

    double          CalculateHeCoreMassAtPhaseEnd() const                           { return m_CoreMass; }                                                                      // McHe(HG) = Core Mass
    double          CalculateHeCoreMassOnPhase() const                              { return m_CoreMass; }                                                                      // McHe(HG) = Core Mass

    double          CalculateLambdaDewi() const;
    double          CalculateLambdaNanjingStarTrack(const double p_Mass, const double p_Metallicity) const;
    double          CalculateLambdaNanjingEnhanced(const int p_MassInd, const int p_Zind) const;

    double          CalculateLuminosityAtPhaseEnd(const double p_Mass) const;
    double          CalculateLuminosityAtPhaseEnd() const                           { return CalculateLuminosityAtPhaseEnd(m_Mass0);}                                           // Use class member variables
    double          CalculateLuminosityOnPhase(const double p_Age, const double p_Mass) const;
    double          CalculateLuminosityOnPhase() const                              { return CalculateLuminosityOnPhase(m_Age, m_Mass0); }                                      // Use class member variables

    double          CalculateMassTransferRejuvenationFactor() const;

    double          CalculateRadialExtentConvectiveEnvelope() const;

    double          CalculateRadiusAtPhaseEnd(const double p_Mass) const;
    double          CalculateRadiusAtPhaseEnd() const                               { return CalculateRadiusAtPhaseEnd(m_Mass); }                                               // Use class member variables
    double          CalculateRadiusOnPhase(const double p_Mass, const double p_Tau, const double p_RZAMS) const;
    double          CalculateRadiusOnPhase() const                                  { return CalculateRadiusOnPhase(m_Mass, m_Tau, m_RZAMS0); }                                 // Use class member variables

    double          CalculateRho(const double p_Mass) const;

    double          CalculateTauAtPhaseEnd() const                                  { return 1.0; }                                                                             // tau = 1.0 at end of HG
    double          CalculateTauOnPhase() const;

    double          ChooseTimestep(const double p_Time) const;

    ENVELOPE        DetermineEnvelopeType() const;

    void            EvolveOneTimestepPreamble()                                     { BaseStar::EvolveOneTimestepPreamble(); }                                                  // Skip MainSequence
    STELLAR_TYPE    EvolveToNextPhase();

    bool            IsEndOfPhase() const                                            { return !ShouldEvolveOnPhase(); }                                                          // Phase ends when age at or after Base Giant Branch MS timescale
    bool            IsSupernova() const                                             { return false; }                                                                           // Not here

    STELLAR_TYPE    ResolveEnvelopeLoss(bool p_NoCheck = false);
    void            ResolveHeliumFlash() {  }                                                                                                                                   // NO-OP
    STELLAR_TYPE    ResolveSkippedPhase()                                           { return m_StellarType; }                                                                   // NO-OP

    bool            ShouldEvolveOnPhase() const                                     { return (utils::Compare(m_Age, m_Timescales[static_cast<int>(TIMESCALE::tBGB)]) < 0); }    // Evolve on HG phase if age < Base Giant Branch timescale
    bool            ShouldSkipPhase() const                                         { return false; }                                                                           // Never skip HG phase

    void            UpdateAgeAfterMassLoss();                                                                                                                                   // Per Hurley et al. 2000, section 7.1

    void            UpdateInitialMass();                                                   // Per Hurley et al. 2000, section 7.1

    
    //Functor for the boost root finder to determine the "initial mass" m_Mass0 based on desired core mass
    template <class T>
    struct Mass0YieldsDesiredCoreMassFunctor
    {
        Mass0YieldsDesiredCoreMassFunctor(HG *p_Star, double p_DesiredCoreMass, ERROR *p_Error)
        {
            m_Star             = p_Star;
            m_DesiredCoreMass  = p_DesiredCoreMass;
            m_Error            = p_Error;
        }
        T operator()(double const& guessMass0)
        {
            HG * copy = new HG(*m_Star, false);
            copy->UpdateAttributesAndAgeOneTimestep(0.0, guessMass0 - copy->Mass0(), 0.0, true);
            double coreMassEstimate=copy->CalculateCoreMassOnPhase(guessMass0, copy->Age());
            delete copy; copy = nullptr;
            return (coreMassEstimate - m_DesiredCoreMass);
        }
    private:
        HG *m_Star;
        double m_DesiredCoreMass;
        ERROR *m_Error;
    };
    
    
    //Root solver to determine "initial mass" m_Mass0 based on desired core mass
    double Mass0ToMatchDesiredCoreMass(HG * p_Star, double p_DesiredCoreMass)
    {
        using namespace std;                                                    // Help ADL of std functions.
        using namespace boost::math::tools;                                     // For bracket_and_solve_root.
        
        double guess  = p_Star->Mass();                                         // Rough guess at solution
        double factor = ADAPTIVE_MASS0_SEARCH_FACTOR;                           // Size of search steps
        
        const boost::uintmax_t maxit = ADAPTIVE_MASS0_MAX_ITERATIONS;            // Limit to maximum iterations.
        boost::uintmax_t it = maxit;                                            // Initally our chosen max iterations, but updated with actual.
        bool is_rising = true;                                                  // So if result with guess is too low, then try increasing guess.
        int digits = std::numeric_limits<double>::digits;                       // Maximum possible binary digits accuracy for type T.
        
        // Some fraction of digits is used to control how accurate to try to make the result.
        int get_digits = digits - 5;                                            // We have to have a non-zero interval at each step, so
        
        // maximum accuracy is digits - 1.  But we also have to
        // allow for inaccuracy in f(x), otherwise the last few
        // iterations just thrash around.
        eps_tolerance<double> tol(get_digits);                                  // Set the tolerance.
        
        std::pair<double, double> root;
        try {
            ERROR error = ERROR::NONE;
            root = bracket_and_solve_root(Mass0YieldsDesiredCoreMassFunctor<double>(p_Star, p_DesiredCoreMass, &error), guess, factor, is_rising, tol, it);
            if (error != ERROR::NONE) SHOW_WARN(error);
        }
        catch(exception& e) {
            SHOW_ERROR(ERROR::TOO_MANY_MASS0_ITERATIONS, e.what());              // Catch generic boost root finding error
        }
        SHOW_WARN_IF(it>=maxit, ERROR::TOO_MANY_MASS0_ITERATIONS);
        
        return root.first + (root.second - root.first)/2;                       // Midway between brackets is our result, if necessary we could return the result as an interval here.
    }
};

#endif // __HG_h__
