#ifndef __Rand_H__
#define __Rand_H__

#define RAND Rand::Instance()

#include <gsl/gsl_rng.h>                                    // GSL random number generator
#include <gsl/gsl_randist.h>


/*
 * Rand Singleton - interface to GSL rng (random number generator)
 *
 * Singletons and global variables are sometimes frowned-upon, but doing it this
 * way means the objects don't need to be passed around to all and sundry.
 * I think convenience and clarity sometimes trump dogma.
 */

class Rand {

private:

    Rand() { m_Rng = NULL; };
    Rand(Rand const&) = delete;
    Rand& operator = (Rand const&) = delete;

    static Rand*      m_Instance;

    gsl_rng*          m_Rng;                                                                           // GSL random number generator

    unsigned long int m_Seed;


public:

    static Rand*  Instance();

    void          Initialise();
    void          Free();

    unsigned long int CurrentSeed()                    { return m_Seed; }
    unsigned long int DefaultSeed()                    { return gsl_rng_default_seed; }
    unsigned long int Seed(const unsigned long p_Seed) { gsl_rng_set(m_Rng, p_Seed); m_Seed = p_Seed; return p_Seed; }

    double        Random();
    double        Random(const double p_Lower, const double p_Upper);
    int           RandomInt(const int p_Lower, const int p_Upper);
    int           RandomInt(const int p_Upper) { return p_Upper < 0 ? 0 : RandomInt(0, p_Upper); }
    double        RandomGaussian(const double p_Sigma);
};


#endif // __Rand_H__
