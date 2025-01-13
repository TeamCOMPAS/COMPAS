#ifndef __HG_h__
#define __HG_h__

#include "constants.h"
#include "typedefs.h"
#include "profiling.h"
#include "utils.h"

#include <boost/math/tools/roots.hpp>

#include "GiantBranch.h"


class BaseStar;
class GiantBranch;

class HG: virtual public BaseStar, public GiantBranch {
    
public:
    
    HG() { m_StellarType = STELLAR_TYPE::HERTZSPRUNG_GAP; };
    
    HG(const BaseStar &p_BaseStar, const bool p_Initialise = true) : BaseStar(p_BaseStar), GiantBranch(p_BaseStar) {
        m_StellarType = STELLAR_TYPE::HERTZSPRUNG_GAP;                                                                                                                          // Set stellar type
        if (p_Initialise) Initialise();                                                                                                                                         // Initialise if required
    }
    
    HG* Clone(const OBJECT_PERSISTENCE p_Persistence, const bool p_Initialise = true) {
        HG* clone = new HG(*this, p_Initialise);
        clone->SetPersistence(p_Persistence);
        return clone;
    }
    
    static HG* Clone(HG& p_Star, const OBJECT_PERSISTENCE p_Persistence, const bool p_Initialise = true) {
        HG* clone = new HG(p_Star, p_Initialise);
        clone->SetPersistence(p_Persistence);
        return clone;
    }
    
    MT_CASE DetermineMassTransferTypeAsDonor() const { return MT_CASE::B; }                                                                                                     // Always case B
    
    
protected:
    
    void Initialise() {
        
        m_Tau = 0.0;                                                                                                                                                            // Start of phase
        
        // update stellar properties at start of HG phase (since core definition changes)
        CalculateGBParams();
        CalculateTimescales();
        // Initialise timescales
        m_Age = m_Timescales[static_cast<int>(TIMESCALE::tMS)];                                                                                                                 // Set age appropriately
        
        // update effective "initial" mass (m_Mass0) so that the core mass is at least equal to the minimum core mass but no more than total mass
        // (only relevant if MANDEL or SHIKAUCHI main sequence core mass prescription is used)
        if (utils::Compare(CalculateCoreMassOnPhase(m_Mass0, m_Age), std::min(m_Mass, MainSequenceCoreMass())) < 0) {
            double desiredCoreMass = std::min(m_Mass, MainSequenceCoreMass());                                                                                                  // desired core mass
            m_Mass0                = Mass0ToMatchDesiredCoreMass(this, desiredCoreMass);                                                                                        // use root finder to find new core mass estimate
            if (m_Mass0 <= 0.0) {                                                                                                                                               // no root found - no solution for estimated core mass
                m_Mass0 = m_Mass;                                                                                                                                               // if no root found we keep m_Mass0 equal to the total mass
            }
            CalculateGBParams();
            CalculateTimescales();
            m_Age = m_Timescales[static_cast<int>(TIMESCALE::tMS)];
        }
        EvolveOnPhase(0.0);
    }
    
    
    // member functions - alphabetically
    double          CalculateCOCoreMassAtPhaseEnd() const                           { return 0.0; }                                                                             // McCO(HG) = 0.0
    double          CalculateCOCoreMassOnPhase() const                              { return 0.0; }                                                                             // McCO(HG) = 0.0
    
    double          CalculateCoreMassAt2ndDredgeUp(const DBL_VECTOR &p_GBParams)    { return p_GBParams[static_cast<int>(GBP::McDU)]; }                                         // NO-OP
    double          CalculateCoreMassAtPhaseEnd(const double p_Mass) const;
    double          CalculateCoreMassAtPhaseEnd() const                             { return CalculateCoreMassAtPhaseEnd(m_Mass0); }                                            // Use class member variables
    double          CalculateCoreMassOnPhase(const double p_Mass, const double p_Time) const;
    double          CalculateCoreMassOnPhase() const                                { return CalculateCoreMassOnPhase(m_Mass0, m_Age); }                                        // Use class member variables
    double          CalculateCoreMassOnPhaseIgnoringPreviousCoreMass(const double p_Mass, const double p_Time) const;                                                           //  Ignore previous core mass constraint when computing expected core mass
    
    double          CalculateCriticalMassRatioClaeys14(const bool p_AccretorIsDegenerate) const;
    double          CalculateCriticalMassRatioHurleyHjellmingWebbink() const        { return 0.25; }                                                                            // As coded in BSE. Using the inverse owing to how qCrit is defined in COMPAS. See Hurley et al. 2002 sect. 2.6.1 for additional details.
    
    double          CalculateHeCoreMassAtPhaseEnd() const                           { return m_CoreMass; }                                                                      // McHe(HG) = Core Mass
    double          CalculateHeCoreMassOnPhase() const                              { return m_CoreMass; }                                                                      // McHe(HG) = Core Mass
    
    double          CalculateHeliumAbundanceCoreAtPhaseEnd() const                  { return 1.0 - m_Metallicity; }
    double          CalculateHeliumAbundanceCoreOnPhase() const                     { return 1.0 - m_Metallicity; }                                                             // Use class member variables
    
    double          CalculateHeliumAbundanceSurfaceAtPhaseEnd() const               { return CalculateHeliumAbundanceSurfaceOnPhase(); }
    double          CalculateHeliumAbundanceSurfaceOnPhase() const                  { return m_InitialHeliumAbundance; }                                                        // Use class member variables
    
    double          CalculateHydrogenAbundanceCoreAtPhaseEnd() const                { return CalculateHydrogenAbundanceCoreOnPhase(); }
    double          CalculateHydrogenAbundanceCoreOnPhase(const double p_Tau) const;
    double          CalculateHydrogenAbundanceCoreOnPhase() const                   { return 0.0; }                                                                             // Star has exhausted hydrogen in its core
    
    double          CalculateHydrogenAbundanceSurfaceAtPhaseEnd() const             { return CalculateHydrogenAbundanceSurfaceOnPhase(); }
    double          CalculateHydrogenAbundanceSurfaceOnPhase() const                { return m_InitialHydrogenAbundance; }                                                      // Use class member variables
    
    
    
    double          CalculateLambdaDewi() const;
    double          CalculateLambdaNanjingStarTrack(const double p_Mass, const double p_Metallicity) const;
    double          CalculateLambdaNanjingEnhanced(const int p_MassIndex, const STELLAR_POPULATION p_StellarPop) const;
    
    double          CalculateLuminosityAtPhaseEnd(const double p_Mass) const;
    double          CalculateLuminosityAtPhaseEnd() const                           { return CalculateLuminosityAtPhaseEnd(m_Mass0);}                                           // Use class member variables
    double          CalculateLuminosityOnPhase(const double p_Age, const double p_Mass) const;
    double          CalculateLuminosityOnPhase() const                              { return CalculateLuminosityOnPhase(m_Age, m_Mass0); }                                      // Use class member variables
    
    double          CalculateMassTransferRejuvenationFactor()                       { return 1.0; }
    
    double          CalculateRadiusAtPhaseEnd(const double p_Mass) const;
    double          CalculateRadiusAtPhaseEnd() const                               { return CalculateRadiusAtPhaseEnd(m_Mass); }                                               // Use class member variables
    double          CalculateRadiusOnPhase(const double p_Mass, const double p_Tau, const double p_RZAMS) const;
    double          CalculateRadiusOnPhase() const                                  { return CalculateRadiusOnPhase(m_Mass0, m_Tau, m_RZAMS0); }                                // Use class member variables
    
    double          CalculateRho(const double p_Mass) const;

    double          CalculateTauAtPhaseEnd() const                                  { return 1.0; }                                                                             // tau = 1.0 at end of HG
    double          CalculateTauOnPhase() const;

    double          ChooseTimestep(const double p_Time) const;

    ENVELOPE        DetermineEnvelopeType() const;

    void            EvolveOneTimestepPreamble()                                     { BaseStar::EvolveOneTimestepPreamble(); }                                                  // Skip MainSequence
    STELLAR_TYPE    EvolveToNextPhase();

    bool            IsEndOfPhase() const                                            { return !ShouldEvolveOnPhase(); }                                                          // Phase ends when age at or after Base Giant Branch MS timescale
    bool            IsSupernova() const                                             { return false; }                                                                           // Not here

    STELLAR_TYPE    ResolveEnvelopeLoss(bool p_Force = false);
    void            ResolveHeliumFlash() {  }                                                                                                                                   // NO-OP
    STELLAR_TYPE    ResolveSkippedPhase()                                           { return m_StellarType; }                                                                   // NO-OP

    bool            ShouldEvolveOnPhase() const                                     { return (utils::Compare(m_Age, m_Timescales[static_cast<int>(TIMESCALE::tBGB)]) < 0); }    // Evolve on HG phase if age < Base Giant Branch timescale
    bool            ShouldSkipPhase() const                                         { return false; }                                                                           // Never skip HG phase

    double          TAMSCoreMass() const                                            { return 0.0; }
    void            UpdateAfterMerger(double p_Mass, double p_HydrogenMass) { }                                                                                                 // Nothing to do for stars beyond the Main Sequence for now
    void            UpdateAgeAfterMassLoss();                                                                                                                                   // Per Hurley et al. 2000, section 7.1

    void            UpdateInitialMass();                                                                                                                                        // Per Hurley et al. 2000, section 7.1

       
    /*
     * Functor for Mass0ToMatchDesiredCoreMass()
     *
     *
     * Constructor: initialise the class
     * template <class T> Mass0YieldsDesiredCoreMassFunctor(HG *p_Star, double p_DesiredCoreMass)
     *
     * @param   [IN]    p_Star                      (Pointer to) The star under examination
     * @param   [IN]    p_DesiredCoreMass           The desired core mass
     * 
     * Function: core mass difference after mass loss
     * T RadiusEqualsRocheLobeFunctor(double const& p_dM)
     * 
     * @param   [IN]    p_GuessMass0                Guess for Mass0
     * @return                                      Difference between estimated core mass (based on p_GuessMass0) and desired core mass (p_DesiredCoreMass)
     */
    template <class T>
    struct Mass0YieldsDesiredCoreMassFunctor {
        Mass0YieldsDesiredCoreMassFunctor(HG *p_Star, double p_DesiredCoreMass) {
            m_Star            = p_Star;
            m_DesiredCoreMass = p_DesiredCoreMass;
        }
        T operator()(double const& p_GuessMass0) {
        
            // We need an estimate of the core mass of the star so we clone the star without
            // initialisation (i.e. we leave it where it is on the phase) so we can calculate
            // and query its core mass.
            //
            // To ensure the clone does not participate in logging, we set its persistence to EPHEMERAL.

            HG *clone = m_Star->Clone(OBJECT_PERSISTENCE::EPHEMERAL, false);
            clone->UpdateAttributesAndAgeOneTimestep(0.0, p_GuessMass0 - clone->Mass0(), 0.0, true);        // update clone's mass and age it one timestep 
            double coreMassEstimate = clone->CalculateCoreMassOnPhase(p_GuessMass0, clone->Age());          // calculate clone's core mass
            delete clone; clone = nullptr;                                                                  // return the memory allocated for the clone

            return (coreMassEstimate - m_DesiredCoreMass);
        }
    private:
        HG *m_Star;
        double m_DesiredCoreMass;
    };
    
    
    /*
     * Root solver to determine "initial" mass (m_Mass0) based on desired core mass
     *
     * Uses boost::math::tools::bracket_and_solve_root()
     *
     *
     * double Mass0ToMatchDesiredCoreMass(HG *p_Star, double p_DesiredCoreMass)
     *
     * @param   [IN]    p_Star                      (Pointer to) The star under examination
     * @param   [IN]    p_DesiredCoreMass           The desired core mass
     * @return                                      Root found: will be -1.0 if no acceptable real root found
     */
    double Mass0ToMatchDesiredCoreMass(HG *p_Star, double p_DesiredCoreMass) {

        const boost::uintmax_t maxit = ADAPTIVE_MASS0_MAX_ITERATIONS;                                       // limit to maximum iterations.
        boost::uintmax_t it          = maxit;                                                               // initially our chosen max iterations, but updated with actual.

        // find root
        // we use an iterative algorithm to find the root here:
        //    - if the root finder throws an exception, we stop and return a negative value for the root (indicating no root found)
        //    - if the root finder reaches the maximum number of (internal) iterations, we stop and return a negative value for the root (indicating no root found)
        //    - if the root finder returns a solution, we check that func(solution) = 0.0 +/ ROOT_ABS_TOLERANCE
        //       - if the solution is acceptable, we stop and return the solution
        //       - if the solution is not acceptable, we reduce the search step size and try again
        //       - if we reach the maximum number of search step reduction iterations, or the search step factor reduces to 1.0 (so search step size = 0.0),
        //         we stop and return a negative value for the root (indicating no root found)

        double guess      = p_Star->Mass();                                                                 // rough guess at solution
        
        double factorFrac = ADAPTIVE_MASS0_SEARCH_FACTOR_FRAC;                                              // search step size factor fractional part
        double factor     = 1.0 + factorFrac;                                                               // factor to determine search step size (size = guess * factor)

        std::pair<double, double> root(p_DesiredCoreMass, 0.0);                                             // initialise root - default return
        std::size_t tries = 0;                                                                              // number of tries
        bool done         = false;                                                                          // finished (found root or exceed maximum tries)?
        Mass0YieldsDesiredCoreMassFunctor<double> func = Mass0YieldsDesiredCoreMassFunctor<double>(p_Star, p_DesiredCoreMass);
        while (!done) {                                                                                     // while no acceptable root found
        
            bool isRising = func((const double)guess) >= func((const double)guess * factor) ? false : true; // gradient direction from guess to upper search increment

            // run the root finder
            // regardless of any exceptions or errors, display any problems as a warning, then
            // check if the root returned is within tolerance - so even if the root finder
            // bumped up against the maximum iterations, or couldn't bracket the root, use
            // whatever value it ended with and check if it's good enough for us - not finding
            // an acceptable root should be the exception rather than the rule, so this strategy
            // shouldn't cause undue performance issues.
            try {
                root = boost::math::tools::bracket_and_solve_root(func, guess, factor, isRising, utils::BracketTolerance, it); // find root
                // root finder returned without raising an exception
                if (it >= maxit) { SHOW_WARN(ERROR::TOO_MANY_MASS0_ITERATIONS); }                           // too many root finder iterations
            }
            catch(std::exception& e) {                                                                      // catch generic boost root finding error
                // root finder exception
                // could be too many iterations, or unable to bracket root - it may not
                // be a hard error - so no matter what the reason is that we are here,
                // we'll just emit a warning and keep trying
                if (it >= maxit) { SHOW_WARN(ERROR::TOO_MANY_MASS0_ITERATIONS); }                           // too many root finder iterations
                else             { SHOW_WARN(ERROR::ROOT_FINDER_FAILED, e.what()); }                        // some other problem - show it as a warning
            }

            // we have a solution from the root finder - it may not be an acceptable solution
            // so we check if it is within our preferred tolerance
            if (fabs(func(root.first + (root.second - root.first) / 2.0)) <= ROOT_ABS_TOLERANCE) {          // solution within tolerance?
                done = true;                                                                                // yes - we're done
            }
            else if (fabs(func(root.first)) <= ROOT_ABS_TOLERANCE) {                                        // solution within tolerance at endpoint 1?
                root.second=root.first;
                done = true;                                                                                // yes - we're done
            }
            else if (fabs(func(root.second)) <= ROOT_ABS_TOLERANCE) {                                       // solution within tolerance at endpoint 2?
                root.first=root.second;
                done = true;                                                                                // yes - we're done
            }
            else {                                                                                          // no - try again
                // we don't have an acceptable solution - reduce search step size and try again
                factorFrac /= 2.0;                                                                          // reduce fractional part of factor
                factor      = 1.0 + factorFrac;                                                             // new search step size
                tries++;                                                                                    // increment number of tries
                if (tries > ADAPTIVE_MASS0_MAX_TRIES || fabs(factor - 1.0) <= ROOT_ABS_TOLERANCE) {         // too many tries, or step size 0.0?
                    // we've tried as much as we can - fail here with -ve return value
                    root.first  = -1.0;                                                                     // yes - set error return
                    root.second = -1.0;
                    SHOW_WARN(ERROR::TOO_MANY_MASS0_TRIES);                                                 // show warning
                    done = true;                                                                            // we're done
                }
            }
        }
        
        return root.first + (root.second - root.first) / 2.0;                                               // midway between brackets is our result, if necessary we could return the result as an interval here.
    }
};

#endif // __HG_h__
