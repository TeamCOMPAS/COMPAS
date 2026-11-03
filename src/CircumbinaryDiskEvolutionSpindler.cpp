#include "CircumbinaryDisk.h"

#include "constants.h"
#include "utils.h"

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <unistd.h>

#include <boost/numeric/odeint.hpp>

// spindler-c is a C library. Its public header currently does not wrap the
// declarations in extern "C", so we do it here for C++/COMPAS builds.
extern "C" {
#include "spindler.h"
}

namespace odeint = boost::numeric::odeint;

// Physics scope of this module:
//   * Disk formation is handled in CircumbinaryDiskFormation.cpp.  This file
//     assumes a CBD reservoir has already formed and evolves the binary using
//     the Siwek/Spindler tabulated derivatives for circumbinary-disk interaction.
//   * The supplied reservoir is delivered at a constant rate over the chosen
//     disk lifetime.  Each component accretes the table-implied fraction unless
//     limited by the selected Eddington prescription; unretained mass removes
//     angular momentum from the binary using the accretor's specific angular momentum.

namespace {

// Backend details and ODE solver settings are local to the Spindler implementation.
constexpr char CBD_SPINDLER_MODEL_NAME[] = "Siwek23";
constexpr double CBD_SPINDLER_ODE_RELATIVE_TOLERANCE = 1.0E-4;
constexpr double CBD_SPINDLER_ODE_ABSOLUTE_TOLERANCE = 1.0E-5;

// Physical and dimensionless inputs used by the internal Spindler integration.
// suppliedMassFraction is the total mass supplied by the CBD in units of the
// initial post-CE binary mass and is assumed to be much smaller than unity; a
// small integration increment does not by itself make an initially large
// reservoir consistent with that assumption. initialMassRatio follows the ordered
// Spindler convention q = M2/M1 <= 1, not necessarily the original COMPAS
// star-2/star-1 label ratio. initialBinaryMassMsol converts between physical Msol
// units and the dimensionless ODE variables.
struct SpindlerEvolutionParams {
    double suppliedMassFraction;
    double initialMassRatio;
    double initialEccentricity;
    double initialBinaryMassMsol;
    double durationYears;

    STELLAR_TYPE stellarType1 = STELLAR_TYPE::NONE;
    STELLAR_TYPE stellarType2 = STELLAR_TYPE::NONE;
    double radius1Rsol = 0.0;
    double radius2Rsol = 0.0;

    CircumbinaryDisk::EvolutionOptions options{};
};

// Final dimensionless state produced by the internal Spindler solver before
// component ordering is converted back to COMPAS labels.
struct SpindlerEvolutionResult {
    double semiMajorAxisFactor = 1.0;
    double mass1Fraction = 0.5;
    double mass2Fraction = 0.5;
    double eccentricity = 0.0;
};

// Owns one initialized spindler-c table set and evaluates CBD-driven binary
// evolution with COMPAS Eddington-limit and angular-momentum-loss extensions.
class SpindlerBackend {
public:
    explicit SpindlerBackend(const std::string& p_ModelName = CBD_SPINDLER_MODEL_NAME);
    ~SpindlerBackend();

    // The backend uniquely owns spindler-c table data, so copying and moving are disabled.
    SpindlerBackend(const SpindlerBackend&) = delete;
    SpindlerBackend& operator=(const SpindlerBackend&) = delete;
    SpindlerBackend(SpindlerBackend&&) = delete;
    SpindlerBackend& operator=(SpindlerBackend&&) = delete;
    // Evolve with Eddington-limited, non-conservative accretion.
    SpindlerEvolutionResult Evolve(const SpindlerEvolutionParams& p_Params);

private:
    // Spindler logarithmic orbital and mass-ratio derivatives returned for one
    // labelled binary state.
    struct Derivatives {
        double da;
        double dq;
        double de;
    };

    Derivatives DerivativesForLabeledQ(double p_QLabeled, double p_Eccentricity) const;
    std::pair<double, double> ComponentMassFractions(double p_QLabeled, double p_Eccentricity) const;
    double ComputeEddingtonRate(
        const STELLAR_TYPE p_StellarType,
        double p_MassFraction,
        double p_InitialBinaryMassMsol,
        double p_RadiusRsol
    ) const;

    static double ApplyEddingtonLimit(
        double p_RequestedRate,
        double p_EddingtonRate,
        CIRCUMBINARY_DISK_EDDINGTON_MODE p_Mode,
        double p_CapFactor,
        double p_Ge23Factor
    );

    // Dimensionless Spindler-internal orbital angular momentum: G = 1, a0 = 1, and Mtot0 = 1.
    static double SpindlerOrbitalAngularMomentum(double p_Mass1Fraction, double p_Mass2Fraction, double p_Eccentricity, double p_SemiMajorAxisFactor);

    struct spindler_data_t m_SpindlerData{};
};

#ifndef SPINDLER_C_ROOT
#define SPINDLER_C_ROOT ""
#endif

/*
 * Validate the ordered initial mass ratio and eccentricity before constructing the ODE state
 *
 * void ValidateCommonInputs(const double p_InitialMassRatio, const double p_InitialEccentricity)
 *
 * @param   [IN]    p_InitialMassRatio        Initial ordered mass ratio, M2/M1 <= 1
 * @param   [IN]    p_InitialEccentricity     Initial eccentricity
 */
void ValidateCommonInputs(const double p_InitialMassRatio, const double p_InitialEccentricity) {
    if (utils::Compare(p_InitialMassRatio, 0.0) <= 0 || !std::isfinite(p_InitialMassRatio)) {
        throw std::invalid_argument("q0 must be positive and finite");
    }
    if (!std::isfinite(p_InitialEccentricity) || utils::Compare(p_InitialEccentricity, 0.0) < 0 || utils::Compare(p_InitialEccentricity, 1.0) >= 0) {
        throw std::invalid_argument("e0 must satisfy 0 <= e0 < 1");
    }
}

/*
 * Convert the dimensionless ODE state into the internal Spindler result
 *
 * SpindlerEvolutionResult MakeResultFromState(const std::vector<double>& p_State)
 *
 * @param   [IN]    p_State                    ODE state [a/a0, M1/M0, M2/M0, e]
 * @return                                      Internal evolution result for this state
 */
SpindlerEvolutionResult MakeResultFromState(const std::vector<double>& p_State) {
    SpindlerEvolutionResult result{};
    result.semiMajorAxisFactor = p_State[0];
    result.mass1Fraction = p_State[1];
    result.mass2Fraction = p_State[2];
    result.eccentricity = p_State[3];
    return result;
}

std::mutex spindlerInitMutex;

// spindler_init() resolves its table paths relative to the process working
// directory. This RAII guard temporarily switches to the spindler-c root and
// guarantees that COMPAS's original working directory is restored on exit.
class ScopedWorkingDirectory {
public:
    /*
     * Temporarily change the process working directory
     *
     * explicit ScopedWorkingDirectory(const std::string& p_NewDirectory)
     *
     * @param   [IN]    p_NewDirectory             Directory required while spindler_init() resolves table paths
     */
    explicit ScopedWorkingDirectory(const std::string& p_NewDirectory) {
        if (p_NewDirectory.empty()) return;

        char* cwd = getcwd(nullptr, 0);
        if (cwd == nullptr) {
            std::ostringstream msg;
            msg << "getcwd failed before spindler_init: " << std::strerror(errno);
            throw std::runtime_error(msg.str());
        }
        m_OldDirectory = cwd;
        std::free(cwd);

        if (chdir(p_NewDirectory.c_str()) != 0) {
            std::ostringstream msg;
            msg << "Could not chdir to SPINDLER_C_ROOT='" << p_NewDirectory
                << "' before spindler_init: " << std::strerror(errno);
            throw std::runtime_error(msg.str());
        }
        m_Changed = true;
    }

    /*
     * Restore the process working directory saved by the constructor
     *
     * ~ScopedWorkingDirectory()
     */
    ~ScopedWorkingDirectory() {
        if (m_Changed) {
            (void)chdir(m_OldDirectory.c_str());
        }
    }

    // Keep directory-restoration ownership unique to a single scope guard.
    ScopedWorkingDirectory(const ScopedWorkingDirectory&) = delete;
    ScopedWorkingDirectory& operator=(const ScopedWorkingDirectory&) = delete;

private:
    std::string m_OldDirectory;
    bool m_Changed = false;
};

/*
 * Initialise one spindler-c table backend
 *
 * SpindlerBackend(const std::string& p_ModelName)
 *
 * @param   [IN]    p_ModelName                Name of the spindler-c table model to load
 */
SpindlerBackend::SpindlerBackend(const std::string& p_ModelName) {
    std::vector<char> modelNameCString(p_ModelName.begin(), p_ModelName.end());
    modelNameCString.push_back('\0');

    // If enabled, the Makefile defines SPINDLER_C_ROOT to the user-supplied
    // spindler-c root.  spindler-c currently looks for tables/<model>/...
    // relative to the current working directory during spindler_init(), so we
    // temporarily initialize from that directory and then restore COMPAS's
    // original working directory.

    std::lock_guard<std::mutex> initLock(spindlerInitMutex);
    ScopedWorkingDirectory cwdGuard{SPINDLER_C_ROOT};

    const int err = spindler_init(modelNameCString.data(), &m_SpindlerData);
    if (err != SPINDLER_NO_ERROR) {
        std::ostringstream msg;
        msg << "spindler_init failed for model '" << p_ModelName
            << "' with error code " << err
            << ". Check SPINDLER_C_ROOT='" << SPINDLER_C_ROOT
            << "' and verify that tables/" << p_ModelName << " exists there.";
        throw std::runtime_error(msg.str());
    }
}

/*
 * Release the table data allocated by spindler_init()
 *
 * ~SpindlerBackend()
 */
SpindlerBackend::~SpindlerBackend() {
    spindler_free_data(&m_SpindlerData);
}

/*
 * Query the tabulated logarithmic orbital derivatives
 *
 * Derivatives DerivativesForLabeledQ(double p_QLabeled, double p_Eccentricity) const
 *
 * @param   [IN]    p_QLabeled                 Mass ratio in the caller's current component labels
 * @param   [IN]    p_Eccentricity             Current eccentricity
 * @return                                      Da, Dq, and De mapped back to the caller's labels
 */
SpindlerBackend::Derivatives SpindlerBackend::DerivativesForLabeledQ(double p_QLabeled, double p_Eccentricity) const {
    if (utils::Compare(p_QLabeled, 0.0) <= 0 || !std::isfinite(p_QLabeled)) {
        throw std::runtime_error("Encountered non-positive or non-finite mass ratio during Spindler integration");
    }

    const bool swapped = utils::Compare(p_QLabeled, 1.0) > 0;
    const double qOrderedRaw = swapped ? 1.0 / p_QLabeled : p_QLabeled;
    const double eccentricityRaw = std::max(p_Eccentricity, 0.0);

    // The table coefficients are logarithmic responses of the binary to CBD
    // accretion/torques.  Da multiplies a, De multiplies e, and Dq describes the
    // evolution of the ordered mass ratio.  If COMPAS's labelled mass ratio is
    // inverted relative to q <= 1, the sign of the q derivative must be inverted
    // before returning to the caller's labels.
    const double qLookup = std::clamp(qOrderedRaw, CBD_SIWEK23_Q_MIN, CBD_SIWEK23_Q_MAX);
    const double eccentricityLookup = std::clamp(eccentricityRaw, CBD_SIWEK23_E_MIN, CBD_SIWEK23_E_MAX);

    Derivatives derivatives{};
    derivatives.da = spindler_get_Da(qLookup, eccentricityLookup, const_cast<struct spindler_data_t*>(&m_SpindlerData));

    const double deLookup = spindler_get_De(qLookup, eccentricityLookup, const_cast<struct spindler_data_t*>(&m_SpindlerData));
    derivatives.de = utils::Compare(eccentricityRaw, 0.0) > 0 ? deLookup * eccentricityLookup / eccentricityRaw : 0.0;

    const double dqLookup = spindler_get_Dq(qLookup, eccentricityLookup, const_cast<struct spindler_data_t*>(&m_SpindlerData));
    const double dqOrdered = dqLookup * qLookup / qOrderedRaw;
    derivatives.dq = swapped ? -dqOrdered : dqOrdered;
    return derivatives;
}

/*
 * Infer the component shares of the supplied CBD mass
 *
 * std::pair<double, double> ComponentMassFractions(double p_QLabeled, double p_Eccentricity) const
 *
 * @param   [IN]    p_QLabeled                 Mass ratio in the caller's current component labels
 * @param   [IN]    p_Eccentricity             Current eccentricity
 * @return                                      Fractions of supplied mass assigned to components 1 and 2
 */
std::pair<double, double> SpindlerBackend::ComponentMassFractions(double p_QLabeled, double p_Eccentricity) const {
    if (utils::Compare(p_QLabeled, 0.0) <= 0 || !std::isfinite(p_QLabeled)) {
        throw std::runtime_error("Encountered non-positive or non-finite mass ratio while computing accretion split");
    }

    const bool swapped = utils::Compare(p_QLabeled, 1.0) > 0;
    const double qOrderedRaw = swapped ? 1.0 / p_QLabeled : p_QLabeled;
    const double eccentricityRaw = std::max(p_Eccentricity, 0.0);
    const double qLookup = std::clamp(qOrderedRaw, CBD_SIWEK23_Q_MIN, CBD_SIWEK23_Q_MAX);
    const double eccentricityLookup = std::clamp(eccentricityRaw, CBD_SIWEK23_E_MIN, CBD_SIWEK23_E_MAX);

    const double dqLookup = spindler_get_Dq(qLookup, eccentricityLookup, const_cast<struct spindler_data_t*>(&m_SpindlerData));
    const double dqOrdered = dqLookup * qLookup / qOrderedRaw;

    // Infer how the supplied CBD mass is split between the two components from
    // the tabulated q evolution.  This keeps the mass growth consistent with the
    // same table derivative that drives the orbital response instead of imposing
    // an independent accretion-ratio prescription.
    double secondaryAccretionFraction = qOrderedRaw * (dqOrdered + 1.0 + qOrderedRaw) / std::pow(1.0 + qOrderedRaw, 2.0);
    secondaryAccretionFraction = std::clamp(secondaryAccretionFraction, 0.0, 1.0);
    const double primaryAccretionFraction = 1.0 - secondaryAccretionFraction;

    if (swapped) {
        return {secondaryAccretionFraction, primaryAccretionFraction};
    }
    return {primaryAccretionFraction, secondaryAccretionFraction};
}

/*
 * Calculate the Eddington accretion rate for one component
 *
 * double ComputeEddingtonRate(const STELLAR_TYPE p_StellarType,
 *                             double             p_MassFraction,
 *                             double             p_InitialBinaryMassMsol,
 *                             double             p_RadiusRsol) const
 *
 * @param   [IN]    p_StellarType             Current stellar type of the component
 * @param   [IN]    p_MassFraction             Component mass divided by the initial binary mass
 * @param   [IN]    p_InitialBinaryMassMsol   Initial binary mass used for ODE normalisation (Msol)
 * @param   [IN]    p_RadiusRsol               Current component radius (Rsol)
 * @return                                      Eddington rate in dimensionless Spindler mass units per year
 */
double SpindlerBackend::ComputeEddingtonRate(
    const STELLAR_TYPE p_StellarType,
    double p_MassFraction,
    double p_InitialBinaryMassMsol,
    double p_RadiusRsol
) const {
    if (p_StellarType == STELLAR_TYPE::BLACK_HOLE) {
        //  Black-hole Eddington mass accretion rate = 4.42e-7 * (M_phys / 10) Msol/yr.  Because the ODE masses
        // are normalised by the initial binary mass, this is returned in
        // dimensionless Spindler mass units per year.
        return CBD_BH_EDDINGTON_RATE_AT_REFERENCE_MASS * (p_MassFraction / CBD_BH_EDDINGTON_REFERENCE_MASS_MSOL);
    }

    const bool isWhiteDwarf = p_StellarType == STELLAR_TYPE::HELIUM_WHITE_DWARF ||
                              p_StellarType == STELLAR_TYPE::CARBON_OXYGEN_WHITE_DWARF ||
                              p_StellarType == STELLAR_TYPE::OXYGEN_NEON_WHITE_DWARF;

    if (p_StellarType == STELLAR_TYPE::NEUTRON_STAR || isWhiteDwarf) {
        if (utils::Compare(p_RadiusRsol, 0.0) <= 0 || !std::isfinite(p_RadiusRsol)) {
            throw std::invalid_argument("radius must be positive and finite for radius-scaled Eddington limits");
        }
        if (utils::Compare(p_InitialBinaryMassMsol, 0.0) <= 0 || !std::isfinite(p_InitialBinaryMassMsol)) {
            throw std::invalid_argument("initialBinaryMassMsol must be positive for radius-scaled Eddington limits");
        }

        // COMPAS uses a radius-scaled Eddington rate for white dwarfs and other
        // non-BH accretors. Preserve the existing CBD neutron-star normalisation
        // while applying the standard COMPAS/Hurley coefficient to white dwarfs.
        const double eddingtonRateMsolPerYear = p_StellarType == STELLAR_TYPE::NEUTRON_STAR
                                     ? CBD_NS_EDDINGTON_RATE_AT_REFERENCE_RADIUS * (p_RadiusRsol / CBD_NS_EDDINGTON_REFERENCE_RADIUS_RSOL)
                                     : CBD_WD_EDDINGTON_RATE_PER_RADIUS * p_RadiusRsol;
        return eddingtonRateMsolPerYear / p_InitialBinaryMassMsol;
    }

    // The Eddington limit is assumed to be irrelevant for non-compact accretors.
    return std::numeric_limits<double>::infinity();
}

/*
 * Apply the selected Eddington-limited retention prescription
 *
 * CAP imposes a hard multiple of the Eddington rate, TM23 reduces the
 * super-Eddington rate following Tuna & Metzger (2023), and GE23 retains a
 * fixed fraction of the requested super-Eddington rate.
 *
 * double ApplyEddingtonLimit(double                                p_RequestedRate,
 *                            double                                p_EddingtonRate,
 *                            CIRCUMBINARY_DISK_EDDINGTON_MODE      p_Mode,
 *                            double                                p_CapFactor,
 *                            double                                p_Ge23Factor)
 *
 * @param   [IN]    p_RequestedRate            Accretion rate requested by the Spindler mass split
 * @param   [IN]    p_EddingtonRate            Eddington accretion rate in matching units
 * @param   [IN]    p_Mode                      Eddington-limiting prescription
 * @param   [IN]    p_CapFactor                 Multiple of the Eddington rate allowed in CAP mode
 * @param   [IN]    p_Ge23Factor                Retained fraction above Eddington in GE23 mode
 * @return                                      Retained accretion rate
 */
double SpindlerBackend::ApplyEddingtonLimit(
    double p_RequestedRate,
    double p_EddingtonRate,
    CIRCUMBINARY_DISK_EDDINGTON_MODE p_Mode,
    double p_CapFactor,
    double p_Ge23Factor
) {
    switch (p_Mode) {
        case CIRCUMBINARY_DISK_EDDINGTON_MODE::CAP: {
            return std::min(p_RequestedRate, p_CapFactor * p_EddingtonRate);
        }

        case CIRCUMBINARY_DISK_EDDINGTON_MODE::TM23: {
            if (std::isfinite(p_EddingtonRate) && utils::Compare(p_EddingtonRate, 0.0) > 0 && utils::Compare(p_RequestedRate, p_EddingtonRate) > 0) {
                const double ratio = p_RequestedRate / p_EddingtonRate;
                const double tm23RetentionFraction = std::min(1.0, std::pow(CBD_TM23_RETENTION_COEFFICIENT / ratio, CBD_TM23_RETENTION_EXPONENT));
                return p_RequestedRate * tm23RetentionFraction;
            }
            return p_RequestedRate;
        }

        case CIRCUMBINARY_DISK_EDDINGTON_MODE::GE23: {
            return utils::Compare(p_RequestedRate, p_EddingtonRate) > 0 ? p_RequestedRate * p_Ge23Factor : p_RequestedRate;
        }
    }

    throw std::invalid_argument("eddingtonMode must be 'cap', 'ge23', or 'tm23'");
}

/*
 * Calculate dimensionless orbital angular momentum in the Spindler normalisation
 *
 * double SpindlerOrbitalAngularMomentum(double p_Mass1Fraction,
 *                                       double p_Mass2Fraction,
 *                                       double p_Eccentricity,
 *                                       double p_SemiMajorAxisFactor)
 *
 * @param   [IN]    p_Mass1Fraction            Ordered primary mass in Spindler units
 * @param   [IN]    p_Mass2Fraction            Ordered secondary mass in Spindler units
 * @param   [IN]    p_Eccentricity             Current eccentricity
 * @param   [IN]    p_SemiMajorAxisFactor      Current semi-major axis divided by its initial value
 * @return                                      Orbital angular momentum for G = Mtot,0 = a0 = 1
 */
double SpindlerBackend::SpindlerOrbitalAngularMomentum(double p_Mass1Fraction, double p_Mass2Fraction, double p_Eccentricity, double p_SemiMajorAxisFactor) {
    const double totalMassFraction = p_Mass1Fraction + p_Mass2Fraction;
    const double reducedMassFraction = p_Mass1Fraction * p_Mass2Fraction / totalMassFraction;
    const double oneMinusE2 = std::max(1.0 - p_Eccentricity * p_Eccentricity, std::numeric_limits<double>::min());
    return reducedMassFraction * std::sqrt(totalMassFraction * p_SemiMajorAxisFactor * oneMinusE2);
}

/*
 * Evolve the Spindler system with Eddington-limited accretion
 *
 * State vector: [a/a0, M1/M0, M2/M0, e] in the ordered Spindler labels.  The
 * independent variable is time in years.  The supplied CBD mass is spread
 * uniformly over p_Params.durationYears, then split between the two stars using
 * the Spindler qdot-derived accretion fractions. The encompassing CBD model assumes the initial
 * reservoir, Initial_CBD_Mass, satisfies m_CBD << M_binary.
 *
 * The orbital terms have two pieces:
 *   - CBD torque/accretion response from the Spindler Da and De tables.
 *   - Extra non-conservative widening/shrinking from angular momentum carried
 *     away by material that was supplied but ejected due to Eddington limits.
 *
 * The table torques scale with the supplied mass flux.  Eccentricity is evolved
 * only through the tabulated CBD response; the simple wind-loss term affects a
 * through Jdot but does not add an independent de/dt prescription.
 *
 * SpindlerEvolutionResult Evolve(const SpindlerEvolutionParams& p_Params)
 *
 * @param   [IN]    p_Params                   Ordered, dimensionless initial state and model parameters
 * @return                                      State reached at the requested duration
 */
SpindlerEvolutionResult SpindlerBackend::Evolve(const SpindlerEvolutionParams& p_Params) {
    ValidateCommonInputs(p_Params.initialMassRatio, p_Params.initialEccentricity);

    if (utils::Compare(p_Params.suppliedMassFraction, 0.0) < 0 || !std::isfinite(p_Params.suppliedMassFraction)) {
        throw std::invalid_argument("suppliedMassFraction must be non-negative and finite");
    }
    if (utils::Compare(p_Params.initialBinaryMassMsol, 0.0) <= 0 || !std::isfinite(p_Params.initialBinaryMassMsol)) {
        throw std::invalid_argument("initialBinaryMassMsol must be positive and finite");
    }
    if (utils::Compare(p_Params.durationYears, 0.0) <= 0 || !std::isfinite(p_Params.durationYears)) {
        throw std::invalid_argument("duration must be positive and finite");
    }
    if (utils::Compare(p_Params.options.eddingtonCapFactor, 0.0) <= 0 || !std::isfinite(p_Params.options.eddingtonCapFactor)) {
        throw std::invalid_argument("eddingtonCapFactor must be positive and finite");
    }
    if (!(utils::Compare(p_Params.options.eddingtonGe23Factor, 0.0) > 0 && utils::Compare(p_Params.options.eddingtonGe23Factor, 1.0) <= 0) ||
        !std::isfinite(p_Params.options.eddingtonGe23Factor)) {
        throw std::invalid_argument("eddingtonGe23Factor must satisfy 0 < f <= 1");
    }

    using StateType = std::vector<double>;
    StateType state = {1.0, 1.0 / (1.0 + p_Params.initialMassRatio), p_Params.initialMassRatio / (1.0 + p_Params.initialMassRatio), p_Params.initialEccentricity};

    if (utils::Compare(p_Params.suppliedMassFraction, 0.0) == 0) {
        return MakeResultFromState(state);
    }

    const double supplyRateFractionPerYear = p_Params.suppliedMassFraction / p_Params.durationYears;

    // Define the ODE derivatives for orbital size, component masses, and
    // eccentricity under CBD torques plus non-conservative mass loss.
    auto system = [this, &p_Params, supplyRateFractionPerYear](const StateType& p_XState, StateType& p_DxDt, const double /*p_Time*/) {
        p_DxDt.resize(4);
        const double semiMajorAxisFactor = p_XState[0];
        const double mass1Fraction = p_XState[1];
        const double mass2Fraction = p_XState[2];
        const double eccentricity = p_XState[3];

        if (utils::Compare(semiMajorAxisFactor, 0.0) <= 0 || utils::Compare(mass1Fraction, 0.0) <= 0 || utils::Compare(mass2Fraction, 0.0) <= 0) {
            throw std::runtime_error("Non-positive a, M1, or M2 encountered during Spindler integration");
        }

        const double totalMassFraction = mass1Fraction + mass2Fraction;
        const double massRatio = mass2Fraction / mass1Fraction;
        const double orbitalAngularMomentum = SpindlerOrbitalAngularMomentum(mass1Fraction, mass2Fraction, eccentricity, semiMajorAxisFactor);

        const double mdotEdd1 = ComputeEddingtonRate(
            p_Params.stellarType1, mass1Fraction, p_Params.initialBinaryMassMsol, p_Params.radius1Rsol
        );
        const double mdotEdd2 = ComputeEddingtonRate(
            p_Params.stellarType2, mass2Fraction, p_Params.initialBinaryMassMsol, p_Params.radius2Rsol
        );

        const auto componentMassFractions = ComponentMassFractions(massRatio, eccentricity);
        const double mdotReq1 = componentMassFractions.first * supplyRateFractionPerYear;
        const double mdotReq2 = componentMassFractions.second * supplyRateFractionPerYear;

        const double mdotAcc1 = ApplyEddingtonLimit(
            mdotReq1,
            mdotEdd1,
            p_Params.options.eddingtonMode,
            p_Params.options.eddingtonCapFactor,
            p_Params.options.eddingtonGe23Factor
        );
        const double mdotAcc2 = ApplyEddingtonLimit(
            mdotReq2,
            mdotEdd2,
            p_Params.options.eddingtonMode,
            p_Params.options.eddingtonCapFactor,
            p_Params.options.eddingtonGe23Factor
        );
        const double mdotLost1 = std::max(mdotReq1 - mdotAcc1, 0.0);
        const double mdotLost2 = std::max(mdotReq2 - mdotAcc2, 0.0);

        // The Spindler tables give the CBD-driven response per unit supplied mass
        // through the disk/binary interaction.
        const auto derivatives = DerivativesForLabeledQ(massRatio, eccentricity);
        const double dadtDisk = derivatives.da * semiMajorAxisFactor * supplyRateFractionPerYear / totalMassFraction;
        const double dedt = derivatives.de * eccentricity * supplyRateFractionPerYear / totalMassFraction;

        const double jBinary = orbitalAngularMomentum / totalMassFraction;
        const double j1 = massRatio * jBinary;
        const double j2 = jBinary / massRatio;
        const double dJdtWind = -(j1 * mdotLost1 + j2 * mdotLost2);
        const double dadtWind = utils::Compare(orbitalAngularMomentum, 0.0) > 0 ? 2.0 * semiMajorAxisFactor * dJdtWind / orbitalAngularMomentum : 0.0;

        p_DxDt[0] = dadtDisk + dadtWind;
        p_DxDt[1] = mdotAcc1;
        p_DxDt[2] = mdotAcc2;
        p_DxDt[3] = dedt;
    };

    auto stepper = odeint::make_controlled(
        CBD_SPINDLER_ODE_ABSOLUTE_TOLERANCE,
        CBD_SPINDLER_ODE_RELATIVE_TOLERANCE,
        odeint::runge_kutta_cash_karp54<StateType>()
    );

    double time = 0.0;
    double dt = std::max(p_Params.durationYears / 1000.0, 1.0e-8);

    // Deliberately use strict comparisons here rather than utils::Compare():
    // the ODE controller must advance to, and clip exactly at, the requested
    // endpoint instead of treating a nearby time as equivalent.
    while (time < p_Params.durationYears) {
        if (time + dt > p_Params.durationYears) {
            dt = p_Params.durationYears - time;
        }

        boost::numeric::odeint::controlled_step_result stepResult;
        do {
            stepResult = stepper.try_step(system, state, time, dt);
        } while (stepResult == boost::numeric::odeint::fail);
    }

    return MakeResultFromState(state);
}

/*
 * Translate COMPAS-labelled input to Spindler variables, evolve it, and translate the result back
 *
 * CircumbinaryDisk::EvolutionTrack EvolveTrackWithBackend(const CircumbinaryDisk::EvolutionInput& p_Input,
 *                                                         SpindlerBackend&                         p_Spindler)
 *
 * @param   [IN]    p_Input                    COMPAS-labelled physical input for one binary-CBD interaction segment
 * @param   [IN]    p_Spindler                 Initialised Spindler table backend
 * @return                                      Final COMPAS-labelled state after the requested increment
 */
CircumbinaryDisk::EvolutionTrack EvolveTrackWithBackend(
    const CircumbinaryDisk::EvolutionInput& p_Input,
    SpindlerBackend& p_Spindler
) {
    const double binaryMassMsol = p_Input.mass1Msol + p_Input.mass2Msol;
    const bool componentsSwapped = utils::Compare(p_Input.mass1Msol, p_Input.mass2Msol) < 0;
    const double primaryMassMsol = componentsSwapped ? p_Input.mass2Msol : p_Input.mass1Msol;
    const double secondaryMassMsol = componentsSwapped ? p_Input.mass1Msol : p_Input.mass2Msol;

    // Convert the public COMPAS-labelled, physical-unit input into the ordered,
    // dimensionless variables expected by the Spindler ODE. componentsSwapped
    // records whether this translation exchanged the COMPAS labels, so the
    // output can be mapped back to the original star ordering.
    SpindlerEvolutionParams params{};
    params.suppliedMassFraction = p_Input.suppliedMassMsol / binaryMassMsol;
    params.initialMassRatio = secondaryMassMsol / primaryMassMsol;
    params.initialEccentricity = p_Input.eccentricity;
    params.initialBinaryMassMsol = binaryMassMsol;
    params.durationYears = p_Input.durationYears;
    params.stellarType1 = componentsSwapped ? p_Input.stellarType2 : p_Input.stellarType1;
    params.stellarType2 = componentsSwapped ? p_Input.stellarType1 : p_Input.stellarType2;
    params.radius1Rsol = componentsSwapped ? p_Input.radius2Rsol : p_Input.radius1Rsol;
    params.radius2Rsol = componentsSwapped ? p_Input.radius1Rsol : p_Input.radius2Rsol;
    params.options = p_Input.options;

    const SpindlerEvolutionResult spindlerResult{p_Spindler.Evolve(params)};

    CircumbinaryDisk::EvolutionTrack track{};
    track.semiMajorAxisFactor = spindlerResult.semiMajorAxisFactor;
    track.eccentricity = spindlerResult.eccentricity;
    track.mass1Msol = (componentsSwapped ? spindlerResult.mass2Fraction : spindlerResult.mass1Fraction) * binaryMassMsol;
    track.mass2Msol = (componentsSwapped ? spindlerResult.mass1Fraction : spindlerResult.mass2Fraction) * binaryMassMsol;

    return track;
}

}  // namespace

namespace CircumbinaryDisk {

/*
 * Evolve a formed CBD and return COMPAS-labelled final state
 *
 * spindler-c is initialized lazily on the first CBD call in a COMPAS process
 * and then reused for all subsequent binaries in that process. Spindler evolves
 * a/a0, so the returned semi-major-axis value is a multiplicative factor that
 * the caller applies to its current binary separation.
 *
 * EvolutionTrack EvolveTrack(const EvolutionInput& p_Input)
 *
 * @param   [IN]    p_Input                    COMPAS-labelled physical input; per-call supplied mass is assumed to be << binary mass
 * @return                                      Final COMPAS-labelled state after the requested increment
 */
EvolutionTrack EvolveTrack(const EvolutionInput& p_Input) {
    static SpindlerBackend spindler{CBD_SPINDLER_MODEL_NAME};
    return EvolveTrackWithBackend(p_Input, spindler);
}


}  // namespace CircumbinaryDisk
