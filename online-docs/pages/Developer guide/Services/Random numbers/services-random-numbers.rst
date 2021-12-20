Random numbers service
======================

A Random Number service is provided, with the ``gsl`` Random Number Generator encapsulated in a singleton object (an instantiation of
the ``Rand`` class).

The ``Rand`` class member variables are private, and public functions have been created for random number functionality required by 
the code.

The Rand service can be accessed by referring to the ``Rand::Instance()`` object. For example, to generate a uniform random floating 
point number in the range [0, 1), call the ``Rand::Random()`` function::

     double u = Rand::Instance()→Random();

Since that could become unwieldy, there is a convenience macro to access the Rand service. The macro just defines "RAND" as 
"Rand::Instance()", so calling the ``Rand::Random()`` function can be written as::

    double u = RAND→Random();

The Rand service must be initialised before use. Initialise the Rand service by calling the ``Rand::Initialise()`` function::

    RAND→Initialise();

Dynamically allocated memory associated with the ``gsl`` random number generator should be returned to the system by calling the 
``Rand::Free()`` function::

    RAND→Free();

before exiting the program.

The Rand service provides the following public member functions:


.. toctree::
   :maxdepth: 1

   ./services-random-numbers-initialise
   ./services-random-numbers-free
   ./services-random-numbers-seed
   ./services-random-numbers-defaultSeed
   ./services-random-numbers-random1
   ./services-random-numbers-random2
   ./services-random-numbers-randomGaussian
   ./services-random-numbers-randomInt1
   ./services-random-numbers-randomInt2

