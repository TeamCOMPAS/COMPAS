Program options
===============

COMPAS provides a rich set of configuration parameters via program options, allowing users to vary many
parameters that affect the evolution of single and binary stars, and the composition of the population of
stars being evolved. Furthermore, COMPAS allows some parameters to be specified as ranges, or sets, of
values via the program options, allowing users to specify a grid of parameter values on the command line.

Combining command-line program options ranges and sets with a :doc:`grid file <../grid-files>` allows users
more flexibility and the ability to specify more complex combinations of parameter values.

Not all program options can be specified as ranges or sets of values. Options for which mixing different
values in a single execution of COMPAS would either not be meaningful, or might cause undesirable results,
such as options that specify the mode of evolution (e.g. ``--mode``), or the name or path of output files
(e.g. ``--output-path``, ``--logfile-detailed-output`` etc.), cannot be specified as a range or set of values.
COMPAS will issue an error message if ranges or sets are specified for options for which they are not supported.

.. toctree::
   :maxdepth: 1

   program-options-ranges
   program-options-sets
   program-options-mixing-ranges-sets
   program-options-vector-options
   program-options-list-defaults

See :doc:`./program-options-list-defaults` for a full list of available program options and their default valaues.

