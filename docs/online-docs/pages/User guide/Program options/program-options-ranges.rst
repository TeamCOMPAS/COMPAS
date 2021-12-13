Program option ranges
=====================

A range of values can be specified for any numeric options (i.e. integer (or integer variant), and floating
point (or floating point variant) data types) that are not excluded from range specifications.

Option value ranges are specified by:

    --option-name range-specifier

where `range-specifier` is defined as:

    range-identifier[start,count,increment]

    .. list-table::
       :widths: 19 81 
       :header-rows: 0
       :class: aligned-text

       * - `range-identifier`
         - is one of {’r’, ’range’} (case is not significant)
       * - `start`
         - is the starting value of the range
       * - `count`
         - is the number of values in the range (must be an unsigned long int)
       * - `increment`
         - is the amount by which the value increments for each value in the range

    Note that:

        `range-identifier` is optional for range-specifier. |br|
        `start` and `increment` must be the same data type as option-name. |br|
        `count` must be a positive integer value.

There should be no spaces inside the brackets ([]). Spaces on the command line are interpreted as argument delimiters
by the shell parser before passing the command-line arguments to the COMPAS executable, so if spaces are present inside
the brackets the shell parser breaks the range specification into multiple command-line arguments.

To specify a range of values for the ``--metallicity`` option, a user, if running COMPAS from the command line
and with no grid file, would type any of the following::

    ./COMPAS --metallicity [0.0001,5,0.0013]

    ./COMPAS --metallicity r[0.0001,5,0.0013]
    
    ./COMPAS --metallicity range[0.0001,5,0.0013]

In each of the examples above the user has specified, by the use of the `range-specifier`, that five binary stars
should be evolved, with constituent star metallicities = 0.0001, 0.0014, 0.0027, 0.0040, and 0.0053.

To evolve a grid of binaries with ten different metallicities, starting at 0.0001 and incrementing by 0.0002,
and five different common envelope alpha values, starting at 0.1 and incrementing by 0.2, the user would
type::

    ./COMPAS --metallicity [0.0001,10,0.0013] --common-envelope-alpha [0.1,5,0.2]

and COMPAS would evolve a grid of 50 binaries using the 10 metallicity values and 5 common envelope alpha values.

Note that when a range is, or ranges are, specified on the command line, the ``--number-of-systems`` command-line option is ignored.
This is to avoid multiple systems with identical initial values being evolved.  Ranges and sets can be mixed with grid files, and
in that case ranges and sets specified on the command line will be played out for each grid file line.

   