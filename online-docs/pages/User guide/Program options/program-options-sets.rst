Program option sets
===================

A set of values can be specified for options of any data type that are not excluded from set specifications (see
note above).

Option value sets are specified by:

    --option-name set-specifier

where `set-specifier` is defined as:

    set-identifier[value\ :sub:`1` ,value\ :sub:`2` ,value\ :sub:`3` , ... ,value\ :sub:`n`]

    and

    .. list-table::
       :widths: 24 76 
       :header-rows: 0
       :class: aligned-text

       * - `set-identifier`
         - is one of {’s’, ’set’} (case is not significant)
       * - `value`:sub:`i`
         - is a value for the option

    Note that:

        `set-identifier` is mandatory for `set-specifier`. |br|
        `value`:sub:`i` must be the same data type as `option-name`.
        
There should be no spaces inside the brackets ([]). Spaces on the command line are interpreted as argument delimiters
by the shell parser before passing the command-line arguments to the COMPAS executable, so if spaces are present inside
the brackets the shell parser breaks the set specification into multiple command-line arguments.

Valid values for boolean options are {1|0, TRUE|FALSE, YES|NO, ON|OFF}, and all set values must be of
the same type (i.e. all 1|0, or all YES|NO etc.).

There is no limit to the number of values specified for a set, values can be repeated, and neither order nor
case is significant.

To specify a set of values for the ``--eccentricity-distribution`` option, a user, if running COMPAS from the
command line and with no grid file, would type any of the following::

    ./COMPAS --eccentricity-distribution s[THERMALISED,FIXED,FLAT]

    ./COMPAS --eccentricity-distribution set[THERMALISED,FIXED,FLAT]

In each of the examples above the user has specified, by the use of the `set-specifier`, that three binary stars
should be evolved, using the eccentricity distributions ’THERMALISED’, ’FIXED’, and ’FLAT’.

Note that when a set is, or sets are, specified on the command line, the ``--number-of-systems`` command-line option is ignored.
This is to avoid multiple systems with identical initial values being evolved.  Ranges and sets can be mixed with grid files, and
in that case ranges and sets specified on the command line will be played out for each grid file line.

   