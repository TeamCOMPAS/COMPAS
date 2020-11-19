#include <algorithm>
#include <time.h>

#include "Rand.h"

Rand* Rand::m_Instance = nullptr;


Rand* Rand::Instance() {

    if (!m_Instance) {
      m_Instance = new Rand();
    }
    return m_Instance;
}


/*
 * Initialise the random number generator
 *
 *
 * void Initialise()
 */
void Rand::Initialise() {

    // Set up the gsl random number generator
    gsl_rng_env_setup();

    // Seed the random number generator
    // Preferably use environment GSL_RNG_SEED, otherwise time(0)
    if (!getenv("GSL_RNG_SEED")) {
        gsl_rng_default_seed = time(NULL);
    }

    if (!m_Rng) {
        m_Rng = gsl_rng_alloc(gsl_rng_default);
    }
}


/*
 * Free the dynamically allocated memory
 *
 *
 * void Free()
 */
void Rand::Free() {
    gsl_rng_free(m_Rng);
}


/*
 * Return a random floating point number uniformly distributed in the range [0.0, 1.0)
 *
 *
 * double Random()
 *
 * @return                                      Random floating point number uniformly distributed in the range [0.0, 1.0)
 */
double Rand::Random() {
    return gsl_rng_uniform(m_Rng);
}


/*
 * Return a random floating point number uniformly distributed in the range [p_Lower, p_Upper), where p_Lower <= p_Upper
 * (p_Lower and p_Upper will be swapped if p_Lower > p_Upper as passed)
 *
 *
 * double Random(const double p_Lower, const double p_Upper)
 *
 * @param   [IN]    p_Lower                     Inclusive lower bound of range of distribution
 * @param   [IN]    p_Upper                     Exclusive upper bound of range of distribution
 * @return                                      Random floating point number uniformly distributed in the range [p_Lower, p_Upper)
 */
double Rand::Random(const double p_Lower, const double p_Upper) {

    double lower = std::min(p_Lower, p_Upper);
    double upper = std::max(p_Lower, p_Upper);

    return (gsl_rng_uniform(m_Rng) * (upper - lower)) + lower;
}


/*
 * Return a random integer number uniformly distributed in the range [p_Lower, p_Upper), where p_Lower <= p_Upper
 * (p_Lower and p_Upper will be swapped if p_Lower > p_Upper as passed)
 *
 *
 * int RandomInt(const int p_Lower, const int p_Upper)
 *
 * @param   [IN]    p_Lower                     Inclusive lower bound of range of distribution
 * @param   [IN]    p_Upper                     Exclusive upper bound of range of distribution
 * @return                                      Random integer number uniformly distributed in the range [p_Lower, p_Upper)
 */
int Rand::RandomInt(const int p_Lower, const int p_Upper) {

    int lower = std::min(p_Lower, p_Upper);
    int upper = std::max(p_Lower, p_Upper);

    return gsl_rng_uniform_int(m_Rng, (upper - lower)) + lower;
}


/*
 * Return a Gaussian random variate, with mean zero and standard deviation p_Sigma
 *
 *
 * double RandomGaussian(const double p_Sigma)
 *
 * @param   [IN]    p_Sigma                     Standard deviation of distribution
 * @return                                      Gaussian random variate, with mean zero and standard deviation p_Sigma
 */
double Rand::RandomGaussian(const double p_Sigma) {
    return gsl_ran_gaussian(m_Rng, p_Sigma);
}
