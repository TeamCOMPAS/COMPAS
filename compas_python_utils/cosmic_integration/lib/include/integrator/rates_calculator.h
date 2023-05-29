#include <iostream>
#include <vector>
#include <cmath>

#ifndef _INTEGRATOR_RATES_CALCULATOR_H_
#define _INTEGRATOR_RATES_CALCULATOR_H_


struct FormationAndMergerRates {
    std::vector<std::vector<double>> formation_rate;
    std::vector<std::vector<double>> merger_rate;
};

FormationAndMergerRates find_formation_and_merger_rates(
    int n_binaries, std::vector<double> redshifts,
    std::vector<double> times,
    double time_first_SF,
    double n_formed,
    std::vector<std::vector<double>> dPdlogZ,
    std::vector<double> metallicities,
    double p_draw_metallicity,
    std::vector<double> COMPAS_metallicites,
    std::vector<double> COMPAS_delay_times,
    std::vector<double> COMPAS_weights
    ) {
    // check if weights were provided, if not use uniform weights
    if (COMPAS_weights.empty()) {
        COMPAS_weights.resize(n_binaries, 1.0);
    }

    // initialize rates to zero
    int n_redshifts = redshifts.size();
    double redshift_step = redshifts[1] - redshifts[0];
    std::vector<std::vector<double>> formation_rate(n_binaries, std::vector<double>(n_redshifts, 0.0));
    std::vector<std::vector<double>> merger_rate(n_binaries, std::vector<double>(n_redshifts, 0.0));

    // interpolate times to redshifts
    auto times_to_redshifts = [&times, &redshifts](double t) {
        // TODO: use an interpolation library
    };

    // make note of the first time at which star formation occurred
    double t0 = time_first_SF;

    // go through each binary in the COMPAS data
    for (int i = 0; i < n_binaries; i++) {
        // if metallicities array is empty, assume all SFR happened at one fixed metallicity
        if (metallicities.empty()) {
            for (int j = 0; j < n_redshifts; j++) {
                formation_rate[i][j] = n_formed * COMPAS_weights[i];
            }
        } else {
            int metallicity_index = 0;
            for (int j = 0; j < n_redshifts; j++) {
                metallicity_index = std::min(static_cast<int>(metallicities.size() - 1), static_cast<int>(std::ceil(COMPAS_metallicites[i])));
                formation_rate[i][j] = n_formed * dPdlogZ[j][metallicity_index] / p_draw_metallicity * COMPAS_weights[i];
            }
        }

        // calculate the time at which the binary formed if it merges at this redshift
        std::vector<double> tform(times.size());
        for (size_t j = 0; j < times.size(); j++) {
            tform[j] = times[j] - COMPAS_delay_times[i];
        }

        // find the index above which the binary would have formed before z = max(redshifts)
        int t0_to_tform_idx = -1;
        for (size_t j = 0; j < tform.size(); j++) {
            if (t0 <= tform[j]) {
                t0_to_tform_idx = j;
                break;
            }
        }

        // include the whole array if t0 is beyond the range
        if (t0_to_tform_idx == -1) {
            t0_to_tform_idx = tform.size();
        }

        // calculate merger rates for formation times at z < max(redshifts)
        if (t0_to_tform_idx > 0) {
            std::vector<double> z_of_formation(t0_to_tform_idx - 1);
            std::vector<int> z_of_formation_index(t0_to_tform_idx - 1);
            for (int j = 0; j < t0_to_tform_idx - 1; j++) {
                z_of_formation[j] = times_to_redshifts(tform[j]);
                z_of_formation_index[j] = std::ceil(z_of_formation[j] / redshift_step);
            }
            for (int j = 0; j < t0_to_tform_idx - 1; j++) {
                merger_rate[i][j] = formation_rate[i][z_of_formation_index[j]];
            }
        }
    }

    FormationAndMergerRates result;
    result.formation_rate = formation_rate;
    result.merger_rate = merger_rate;
    return result;
}

#endif