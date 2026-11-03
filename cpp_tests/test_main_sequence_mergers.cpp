#include <cmath>
#include <iostream>
#include <string>

#include "CH.h"
#include "MS_gt_07.h"
#include "MS_lte_07.h"
#include "Options.h"
#include "Rand.h"

// Normally supplied by main.cpp, which is not linked into this regression executable.
OBJECT_ID globalObjectId = 1;
OBJECT_ID m_ObjectId = 0;

namespace {
int checks = 0;
int failures = 0;

void Check(bool condition, const std::string& label) {
    checks++;
    if (!condition) {
        failures++;
        std::cerr << "FAIL: " << label << '\n';
    }
}

bool Close(double actual, double expected, double tolerance = 1.0E-11) {
    return std::isfinite(actual) && std::abs(actual - expected) <= tolerance * std::max(1.0, std::abs(expected));
}

template <typename Phase>
class MergerProbe : public Phase {
public:
    explicit MergerProbe(double mass)
        : BaseStar(42, mass, 0.0142, KickParameters{}), Phase(static_cast<const BaseStar&>(*this)) {
        this->SetPersistence(OBJECT_PERSISTENCE::EPHEMERAL);
        this->m_SurfaceMagneticFlux = 1.0E5;  // Exercise flux continuity with a nonzero field.
        // Keep a CHE probe above its rotation threshold during the forced attribute update.
        if (this->StellarType() == STELLAR_TYPE::CHEMICALLY_HOMOGENEOUS) this->SetOmega(100.0 * std::max(1.0, this->OmegaCHE()));
    }

    void Merge(double mass, double hydrogenMass) { MainSequence::UpdateAfterMerger(mass, hydrogenMass); }
    void Refresh() { this->UpdateAttributesAndAgeOneTimestep(0.0, 0.0, 0.0, true); }
    double Lifetime() const { return this->m_Timescales[static_cast<int>(TIMESCALE::tMS)]; }
    double Flux() const { return this->m_SurfaceMagneticFlux; }
    double ProcessedFraction(double tau, double mass) const { return MainSequence::CalculateMinimumCoreMassFraction(tau, mass); }
    double CoreHydrogenFraction(double tau, double mass) const { return MainSequence::CalculateHydrogenAbundanceCoreFraction(tau, mass); }
};

template <typename Phase>
void CheckMixedRemnants(double initialMass) {
    for (double fraction : {1.0, 0.95, 0.8, 0.2, 0.0}) {
        MergerProbe<Phase> star(initialMass);
        double finalMass = 1.2 * initialMass;
        double mixedHydrogen = star.InitialHydrogenAbundance() * fraction;
        double hydrogenMass = finalMass * mixedHydrogen;
        star.Merge(finalMass, hydrogenMass);

        std::string label = "M=" + std::to_string(initialMass) + ", X/X0=" + std::to_string(fraction);
        Check(Close(star.Mass(), finalMass), label + ": remnant mass");
        Check(Close(star.HydrogenAbundanceCore(), mixedHydrogen), label + ": mixed core hydrogen");
        Check(Close(star.HydrogenAbundanceSurface(), mixedHydrogen), label + ": mixed surface hydrogen");
        Check(Close(star.HydrogenAbundanceCore() + star.HeliumAbundanceCore() + star.Metallicity(), 1.0), label + ": core abundance sum");
        Check(Close(star.Mass() * star.HydrogenAbundanceCore(), hydrogenMass), label + ": homogeneous hydrogen budget");
        Check(std::isfinite(star.Age()) && star.Tau() >= 0.0 && star.Tau() <= 1.0, label + ": bounded age");

        double age = star.Age();
        double flux = star.Flux();
        star.Refresh();
        Check(Close(star.HydrogenAbundanceCore(), mixedHydrogen), label + ": zero-time refresh preserves hydrogen");
        Check(Close(star.Age(), age), label + ": zero-time refresh preserves age");
        Check(Close(star.Flux(), flux), label + ": zero-time refresh preserves flux");

        bool brcekCore = OPTIONS->MainSequenceCoreMassPrescription() == CORE_MASS_PRESCRIPTION::BRCEK &&
                         initialMass >= OPTIONS->MainSequenceCoreMassBrcekMinimumMass();
        if (fraction == 0.8 && !brcekCore) {
            double before = star.HydrogenAbundanceCore();
            star.UpdateAttributesAndAgeOneTimestep(0.0, 0.0, 1.0E-6 * star.Lifetime(), true);
            Check(star.HydrogenAbundanceCore() <= before + 1.0E-12, label + ": burning does not restore hydrogen");
            if (OPTIONS->MainSequenceMagneticFluxDecay() == MAIN_SEQUENCE_MAGNETIC_FLUX_DECAY::EXPONENTIAL)
                Check(star.Flux() > 0.0 && star.Flux() <= flux, label + ": aging does not restore magnetic flux");
        }
    }
}
} // namespace

int main(int argc, char* argv[]) {
    RAND->Initialise();
    RAND->Seed(0);
    if (!OPTIONS->Initialise(argc, argv)) return 2;

    CheckMixedRemnants<MS_lte_07>(0.7);
    for (double mass : {1.0, 1.15, 1.2, 1.3, 5.0, 8.0, 20.0, 60.0}) CheckMixedRemnants<MS_gt_07>(mass);

    MergerProbe<MS_gt_07> example(12.0);
    example.Merge(20.0, 20.0 * example.InitialHydrogenAbundance() * 0.8);
    bool mist = OPTIONS->MainSequenceCoreHydrogenPrescription() == CORE_HYDROGEN_PRESCRIPTION::MIST;
    bool explicitCore = OPTIONS->MainSequenceCoreMassPrescription() == CORE_MASS_PRESCRIPTION::BRCEK &&
                        example.MZAMS() >= OPTIONS->MainSequenceCoreMassBrcekMinimumMass();
    // 0.3503778302887089 with the MIST fit of v03.55.02 (fitted to 1.3 - 20 Msol); v03.55.03 adds the high-mass terms above 8 Msol
    Check(Close(example.Tau(), mist && !explicitCore ? 0.3248657515089267 : 0.2), "20-Msol review example: fractional age");

    // MIST convective-core fit of X_c / X_0, independently: the v03.55.00 coefficients below 8 Msol (L = 0, unchanged) and the
    // v03.55.03 high-mass terms c2 L, d2 L above it
    for (double mass : {1.3, 5.0, 8.0, 20.0, 40.0, 100.0}) {
        double s = std::pow(1.3 / mass, 3.0);
        double L = std::max(0.0, std::log10(mass / 8.0));
        double p = 1.2467 - 0.1960 * s - 0.0926 * L;
        double q = 0.7074 + 0.1585 * s + 0.1133 * L;
        for (double tau : {0.1, 0.5, 0.9}) {
            double expected = mist ? std::pow(std::max(0.0, 1.0 - std::pow(tau, p)), q) : 1.0 - tau;
            Check(Close(example.CoreHydrogenFraction(tau, mass), expected), "core hydrogen fit at M=" + std::to_string(mass) + ", tau=" + std::to_string(tau));
        }
    }

    for (double mass : {1.0, 8.0, 20.0}) {
        double tau = 0.5;
        // Independent normalisation check of the existing MIST hydrogen-processed-mass fit.
        double exponent = 1.15 + 0.72 * std::max(0.0, 1.4 - mass);
        Check(Close(example.ProcessedFraction(tau, mass), mist ? std::pow(tau, exponent) : tau), "processed hydrogen fraction at M=" + std::to_string(mass));
        Check(Close(example.ProcessedFraction(0.0, mass), 0.0), "ZAMS hydrogen budget");
        Check(Close(example.ProcessedFraction(1.0, mass), 1.0), "TAMS hydrogen budget");
    }

    // CHE overrides ordinary-MS depletion with its own linear composition law even when MIST is selected.
    MergerProbe<CH> homogeneous(12.0);
    homogeneous.Merge(20.0, 20.0 * homogeneous.InitialHydrogenAbundance() * 0.8);
    Check(Close(homogeneous.Tau(), 0.2), "CHE retains linear age mapping");
    Check(Close(homogeneous.HydrogenAbundanceCore(), 0.8 * homogeneous.InitialHydrogenAbundance()), "CHE mixed composition");

    std::cout << checks - failures << '/' << checks << " merger checks passed; " << failures << " failures\n";
    return failures ? 1 : 0;
}
