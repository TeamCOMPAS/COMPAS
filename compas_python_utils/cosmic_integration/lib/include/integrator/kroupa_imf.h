/**
 * Kroupa IMF
 *
 * Sample masses from the Kroupa 2001 IMF
 * [http://adsabs.harvard.edu/abs/2001MNRAS.322..231K, eqn 2]
 *
 */

#ifndef _INTEGRATOR_KROUPA_IMF_H_
#define _INTEGRATOR_KROUPA_IMF_H_

#include <random>
#include <vector>
#include <cmath>

namespace integrator {
  namespace kroupa_imf {

    double * sample_n_masses(int n);
    double powerlaw_cdf(double x, double x0, double b, double a);
    double * get_normalisation_constants(double * bounds, double * slopes);
    double generate_mass_from_inv_cdf(double a, double b, double U, double F, double m);
    double * inverse_CDF_IMF(double * U, int n);
    double CDF_IMF(double m, double * bounds, double * slopes, double * norms);






    // Sample n masses from the Kroupa IMF using the inverse CDF method
    double * sample_n_masses(int n) {
      double * U = new double[n];
      double * masses;

      for (int i = 0; i < n; i++) {
        U[i] = (double) rand() / RAND_MAX;
      }
      masses = inverse_CDF_IMF(U, n);
      return masses;
    }

    // Inverse Cumulative Distribution Function for the Kroupa IMF
    double * inverse_CDF_IMF(double * U, int n) {
      // Mass Bounds and slopes for the Kroupa IMF
      double bounds[4] = {
        0.01,
        0.08,
        0.5,
        200
      };
      double slopes[3] = {
        0.3,
        1.3,
        2.3
      };

      double * norms = get_normalisation_constants(bounds, slopes);
      double * F = new double[4];
      for (int j = 0; j < 4; j++) {
        F[j] = CDF_IMF(bounds[j], bounds, slopes, norms);
      }

      double * masses = new double[n];
      for (int i = 0; i < (int) n; i++) {
        masses[i] = 0;
        for (int j = 0; j < 3; j++) {
          bool in_range = (F[j] < U[i] && U[i] <= F[j + 1]);
          if (in_range) {
            masses[i] = generate_mass_from_inv_cdf(slopes[j], norms[j], U[i], F[j], bounds[j]);
          }
        }
      }
      return masses;
    }

    double CDF_IMF(double m, double * bounds, double * slopes, double * norms) {
      if (m <= bounds[0]) {
        return 0;
      } else if (m <= bounds[1]) {
        return powerlaw_cdf(m, bounds[0], norms[0], slopes[0]);
      } else if (m <= bounds[2]) {
        double previous_cdf = powerlaw_cdf(bounds[1], bounds[0], norms[0], slopes[0]);
        double current_cdf = powerlaw_cdf(m, bounds[1], norms[1], slopes[1]);
        return previous_cdf + current_cdf;
      } else if (m <= bounds[3]) {
        return 1;
      }
      return 0;
    }

    double powerlaw_cdf(double x, double x0, double b, double a) {
      return b / (1 - a) * (pow(x, 1 - a) - pow(x0, 1 - a));
    }

    double * get_normalisation_constants(double * bounds, double * slopes) {
      double m1 = bounds[0], m2 = bounds[1], m3 = bounds[2], m4 = bounds[3];
      double a1 = slopes[0], a2 = slopes[1], a3 = slopes[2];
      double b1 = 1 / ((pow(m2, 1 - a1) - pow(m1, 1 - a1)) / (1 - a1) + pow(m2, -(a1 - a2)) * (pow(m3, 1 - a2) - pow(m2, 1 - a2)) / (1 - a2) + pow(m2, -(a1 - a2)) * pow(m3, -(a2 - a3)) * (pow(m4, 1 - a3) - pow(m3, 1 - a3)) / (1 - a3));
      double b2 = b1 * pow(m2, -(a1 - a2));
      double b3 = b2 * pow(m3, -(a2 - a3));

      double * norms = new double[3];
      norms[0] = b1;
      norms[1] = b2;
      norms[2] = b3;
      return norms;
    }

    double generate_mass_from_inv_cdf(double a, double b, double U, double F, double m) {
      double v = pow((1 - a) / b * (U - F) + pow(m, 1 - a), 1 / (1 - a));
      return v;
    }



    double compute_star_forming_mass_per_binary(double binaryFraction, double Mlower, double Mupper, double m2_min, int n) {
        double* primary_mass = sample_n_masses(n);
        double total_mass = 0.0;
        double total_mass_in_cut = 0.0;
        int count_in_cut = 0;

        // Seed the random number generator
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<double> dis(0.0, 1.0);

        for (int i = 0; i < n; i++) {
            double mass_ratio = dis(gen);
            bool is_binary = dis(gen) < binaryFraction;
            double secondary_mass = 0.0;
            if (is_binary) {
                secondary_mass = primary_mass[i] * mass_ratio;
            }

            total_mass += primary_mass[i] + secondary_mass;

            bool primary_in_cut = (primary_mass[i] >= Mlower) && (primary_mass[i] <= Mupper);
            bool secondary_in_cut = secondary_mass > m2_min;
            bool binary_in_cut = primary_in_cut && secondary_in_cut;
            if (binary_in_cut) {
                total_mass_in_cut += primary_mass[i] + secondary_mass;
                count_in_cut++;
            }
        }

        double f_mass_sampled = total_mass_in_cut / total_mass;
        double average_mass_in_cut = total_mass_in_cut / count_in_cut;

        // find the average star forming mass evolved per binary in the Universe
        double average_mass_evolved_per_binary = average_mass_in_cut / f_mass_sampled;
        return average_mass_evolved_per_binary;
    }


  }
}
#endif