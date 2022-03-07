Services
========

A number of services have been provided to help simplify the code. The code for each service is encapsulated in a singleton object
(an instantiation of the relevant class). The singleton design pattern allows the definition of a class that can only be instantiated
once, and that instance effectively exists as a global object available to all the code without having to be passed around as a 
parameter. Singletons are a little anti-OO, but provided they are used judiciously are not necessarily a bad thing, and can be very 
useful in certain circumstances.

Services provided are:


.. toctree::
   :maxdepth: 1

   ./Program options/services-program-options
   ./Random numbers/services-random-numbers
   ./Logging and debugging/services-logging-debugging
   ./services-error-handling