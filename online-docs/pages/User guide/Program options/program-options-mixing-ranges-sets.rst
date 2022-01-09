Mixing ranges and sets
======================

Ranges and sets can be specified together, and there is no limit to the number of ranges or sets that can be
specified on the command line, or in the :doc:`grid file <../grid-files>`.

Running COMPAS with the command::

    ./COMPAS --metallicity r[0.0001,5,0.0013] --common-envelope-alpha s[0.1,0.2,0.6,0.9]

would result in 20 binaries being evolved: 5 for the range of metallicities, times 4 for the set of CE alpha values.


Consider the following grid file, named `gridfile.txt`::

    --metallicity r[0.0001,5,0.0013] --common-envelope-alpha s[0.1,0.2,0.6,0.9]
    --fryer-supernova-engine s[rapid,delayed] --eccentricity r[0.0001,3,0.0003]

Running COMPAS with the command::

    ./COMPAS --grid gridfile.txt

would result in 26 binaries being evolved:

- 20 for the first grid line (5 for the range of metallicities, times 4 for the set of CE alpha values), and |br|
- 6 for the second grid line (2 for the set of Fryer SN engine values, and 3 for the range of eccentricities)


Running COMPAS with the command::

    ./COMPAS --remnant-mass-prescription s[mullermandel,fryer2012,hurley2000,muller2016] --grid gridfile.txt

would result in 104 binaries being evolved: the grid file would be ‘executed’ for each of the four remnant
mass prescriptions specified on the command line.

Multiple ranges and/or sets can be specified on the command line, and on each line of the grid file – so very
large numbers of stars/binaries can be evolved with just a few range/set specifications.
