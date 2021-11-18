Grid files
==========

A grid file allows users to specify initial values for multiple systems for both Single Star Evolution (SSE) and Binary Star Evolution 
(BSE).  Each line of a grid file is used by COMPAS to set the initial conditions and evolutionary parameters for an individual single 
star (SSE) or binary star (BSE), and each single star or binary star defined by a grid file line is evolved using those values.

Each line of a grid file is a set of program option specifications, with the specifications being exactly as they would appear on the 
command line if running COMPAS from the command line.

For example, a grid file could contain the following two lines:

    --metallicity 0.001 --eccentricity 0.0 --remnant-mass-prescription fryer2012 |br|
    --remnant-mass-prescription mullermandel --metallicity 0.02 --semi-major-axis 45.678

in which case COMPAS would evolve two binaries, with the option values set per the grid file lines.

Grid files can have blank lines and comments. Comments begin with a hash/pound character ('#') - the hash/pound character and text 
beyond it are ignored by COMPAS.

Not all program options can be specified in a grid file. Options that should remain constant for a single execution of COMPAS, such as
options that specify the mode of evolution (e.g. ``--mode``), or the name or path of output files (e.g. ``--output-path``, 
``--logfile-detailed-output`` etc.) can only be specified on the command line.  COMPAS will issue an error message if an option that is
not supported in a grid file is specified on a grid file line.

COMPAS imposes no limit on the number of grid file lines: the size of a grid file is limited only by the filesystem of the host system.


Specifying a Subset of the Grid File to be Processed
----------------------------------------------------

Users can instruct COMPAS to process only a subset of a specified grid file. This is achieved via the program options ``--grid-start-line``,
and ``--grid-lines-to-process``:

    The ``--grid-start-line`` program option takes a single parameter: an integer specifying the zero-based line number of the first line of 
    the grid file to be processed. The default value is 0 - the first line of the grid file.  Specifying a start line beyond the end of the 
    grid file will result in an UNEXPECTED-END-OF-FILE error.

    The ``--grid-lines-to-process`` program option takes a single parameter: an integer specifying the number of lines of the grid file to be 
    processed. The default is to process all lines in the grid file from the start line (which may have been specified by the 
    ``--grid-start-line option``) through to the end of the grid file. Specifying a number of lines to be processed that, when coupled with 
    the start line, would result in attempting to process lines beyond the end of the grid file will result in an INCOMPLETE-GRID error.
    
    Note that blank lines and comments count towards the number of grid lines processed by COMPAS when deciding if the number specified by the
    user has been reached.

Both option ``--grid-start-line`` and ``--grid-lines-to-process`` are ignored if no grid file is specified via the ``--grid`` program option.
