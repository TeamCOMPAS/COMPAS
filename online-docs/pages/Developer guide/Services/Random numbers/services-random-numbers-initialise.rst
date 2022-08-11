Rand::Initialise()
==================

::

    VOID Rand::Initialise()

Initialises the ``gsl`` random number generator. If the environment variable ``GSL_RNG_SEED`` exists, the ``gsl`` random number generator
is seeded with the value of the environment variable, otherwise it is seeded with the current time.