Configuration
=============

Run-time Configuration
----------------------

COMPAS is configured at run-time via command-line options, referred to as "program options" in this documentation, and, optionally, 
grid files.

Configuring COMPAS via command-line options and grid files gives users the flexibility to specify both the initial conditions and the 
evolutionary parameters that define single and binary stars at birth, and the conditions under which the stars and systems evolve over 
their lifetime.

COMPAS has a rich set of command-line options that can be set by the user, and these, coupled with the flexibility afforded by the 
COMPAS grid file implementation, allow users to configure COMPAS over a broad range of initial conditions and evolutionary parameters.
This provides users with the means to explore a broad range of physics and physical assumptions.


Compile-time Configuration
--------------------------

The values of some physical constants, bounds for some initial conditions, evolutionary parameters, and physical processes, etc., are 
specified in the COMPAS source file `constants.h`.  While it is unlikely that these constants would need to be changed in most ordinary 
COMPAS runs, the possibility exists that users may want to change some of them.  Should that be the case, the user should change the 
value(s) required in ``constants.h`` and rebuild the COMPAS executable. A makefile is provided in the source directory.

See :doc:`../Getting started/building-COMPAS` for details of how to build the COMPAS executable.
