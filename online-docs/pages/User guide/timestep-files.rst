Timestep files
==============

A timestep file allows users to specify the timesteps to be taken during Single Star Evolution (SSE) or Binary Star Evolution (BSE).

If a timestep file is provided (via the ``timesteps-filename`` program option), the COMPAS code uses the user-provided timeteps to
evolve stars/systems instead of calculating a new timestep at each iteration. The COMPAS code will read the first non-blank and
non-comment line of the timesteps file and set the first timestep to that value, then for each iteration during evolution, each
subsequent non-blank and non-comment line of the timesteps file is read and the timestep for that iteration set to the value read.

Note that COMPAS will use the timestep exactly as read - the timesteps read from the timesteps file are not quantised by COMPAS, and 
neither will they be truncated to any limit (e.g. ``ABSOLUTE_MINIMUM_TIMESTEP`` from ``constants.h``) or multiplied by any timestep
multiplier set by the ``timestep-multiplier`` option.  Furthermore, no check will be made after the timestep is taken to limit the
radial expansion of the star being evolved.

If the timesteps provided in the timesteps file are completely consumed (i.e. each timestep in the timesteps file has been taken
during evolution), but evolution is not complete (i.e. no stopping condition for evolution has been met), evolution is stopped and
a warning issued, and the status of the star/system being evolved is set to ``"User-provided timesteps exhausted"``.

If the timesteps provided in the timesteps file are not completely consumed when evolution is complete (i.e. a stopping condition for
evolution has been met), a warning issued, and the status of the star/system evolved is set to ``"User-provided timesteps not consumed"``.

Timesteps files can have blank lines and/or comment lines. Comments begin with a hash/pound character ('#') - the hash/pound character
and text beyond it are ignored by COMPAS.

The ``timesteps-filename`` option is valid both on the command line and in a grid file.  If the ``timesteps-filename`` option is specified
in conjunction with the ``--number-of-systems`` (``-n``) option, or with option ranges and/or sets, the same timeteps file specified will
be used for each star or system being evolved. On the other hand, if the ``timesteps-filename`` option is specified in a grid file, the
timesteps file specified on a grid file line will be read and used for the star or system evolved by that grid file line.

Each line of a timesteps file which is not a blank line or comment line must contain a single, non-negative, floating-point number (in
plain text - COMPAS will read the plain text and convert it to a floating-point number). An excerpt from a timesteps file is shown below::

    # timesteps files can have blank lines and/or comment lines
    # comment lines begin with the '#' character

    # timesteps file for my experiment...
    # description of file here
    
    0.201713464000000
    0.201935435000000

    0.020216053100000
    0.020218312700000
    0.020220604600000
    0.002430922190000
    0.002430922190000

    # The following timestep is the NUCLEAR_MINIMUM_TIMESTEP (10^6 * TIMESTEP_QUANTUM = 10^6 * 10^-12 = 10^-6 Myr)
    0.000001000000000

    0.000171981412000
    0.000168541784000
    0.000165170948000
    0.000161867529000
    0.000099674841500
    0.000097681344700
    0.000095727717800
    0.000093813163400


COMPAS imposes a hard limit of ``1,000,000`` timesteps in a timesteps file.

