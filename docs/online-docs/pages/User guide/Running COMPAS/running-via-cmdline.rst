Running COMPAS from the command line
====================================

COMPAS is a command-line application.  Interaction with COMPAS is entirely through the terminal and shell - there is no
visual or graphical user interface (GUI).  COMPAS interacts with the user by accepting input via the keyboard, and providing
ouput by writing plain text to the terminal. COMPAS reads input files where necessary: a ``grid`` file (see :doc:`../grid-files`),
and a log file definitions file (see :doc:`../COMPAS output/standard-logfiles-record-specification`), and produces output files
(see :doc:`../COMPAS output/output`), but these are not interactive.

Command-line applications accept interactive input from the user in a number of ways: one of those is via command-line switches and 
arguments, or, more generally, command-line options. This is the method COMPAS uses to interact with the user.  

A few example COMPAS runs, using a small sample of available program option, are shown below. For detailed information regarding program 
option use, and a full list of program options available including their default values see :doc:`../Program options/program-options`.


To run COMPAS from the command line, and have COMPAS output its version string, type either of::

    ./compas -v
    ./compas --version

This should produce output that looks something like::

    COMPAS v02.22.00
    Compact Object Mergers: Population Astrophysics and Statistics
    by Team COMPAS (http://compas.science/index.html)
    A binary star simulator


To run COMPAS from the command line, and have COMPAS output a list of the command-line options available, type::

    ./compas -h

Note the single dash before the option name.  This should produce output that looks something like::

    COMPAS v02.22.00
    Compact Object Mergers: Population Astrophysics and Statistics
    by Team COMPAS (http://compas.science/index.html)
    A binary star simulator

    Options:
    --PISN-lower-limit
    --PISN-upper-limit
    --PPI-lower-limit
    --PPI-upper-limit
    --add-options-to-sysparms
    --allow-rlof-at-birth
    --allow-touching-at-birth
    --angular-momentum-conservation-during-circularisation
    --black-hole-kicks
    --case-BB-stability-prescription
    --check-photon-tiring-limit
    --chemically-homogeneous-evolution

    ...

The options are listed in `ASCIIbetical` order (digits come before uppercase alpha characters, which come before lowercase alpha characters).


To run COMPAS from the command line, and have COMPAS output a list of the command-line options available, together with a brief description
of  the option, and its default value in the COMPAS code, type::

    ./compas --help

Note two dashes before the option name.  This should produce output that looks something like::

    COMPAS v02.22.00
    Compact Object Mergers: Population Astrophysics and Statistics
    by Team COMPAS (http://compas.science/index.html)
    A binary star simulator

    Options:
    --PISN-lower-limit
      Minimum core mass for PISN (default = 60.000000)
    --PISN-upper-limit
      Maximum core mass for PISN (default = 135.000000)
    --PPI-lower-limit
      Minimum core mass for PPI (default = 35.000000)
    --PPI-upper-limit
      Maximum core mass for PPI (default = 60.000000)
    --add-options-to-sysparms
      Add program options columns to BSE/SSE SysParms file (options: [ALWAYS, GRID, NEVER], default = GRID)
    --allow-rlof-at-birth
      Allow binaries that have one or both stars in RLOF at birth to evolve (default = TRUE)
    --allow-touching-at-birth
      Allow binaries that are touching at birth to evolve (default = FALSE)
    --angular-momentum-conservation-during-circularisation
      Conserve angular momentum when binary is circularised when entering a Mass Transfer episode (default = FALSE)
    --black-hole-kicks
      Black hole kicks relative to NS kicks (options: [FULL, REDUCED, ZERO, FALLBACK], default = FALLBACK)
    --case-BB-stability-prescription
      Case BB/BC mass transfer stability prescription (options: [ALWAYS_STABLE, ALWAYS_STABLE_ONTO_NSBH, TREAT_AS_OTHER_MT, ALWAYS_UNSTABLE], default = ALWAYS_STABLE)
    --check-photon-tiring-limit
      Check the photon tiring limit hasn't been exceeded by wind mass loss (default = FALSE)
    --chemically-homogeneous-evolution
      Chemically Homogeneous Evolution (options: [NONE, OPTIMISTIC, PESSIMISTIC], default = PESSIMISTIC)

    ...

Again, the options are listed in `ASCIIbetical` order.


To run a default COMPAS run of 10 binary systems with default initial conditions and evolutionary parameters, type::

    ./compas

This should produce an output put similar to::

    COMPAS v02.22.00
    Compact Object Mergers: Population Astrophysics and Statistics
    by Team COMPAS (http://compas.science/index.html)
    A binary star simulator

    Start generating binaries at Tue Sep  7 17:34:49 2021

    0: Stars merged: (Main_Sequence_>_0.7 -> Main_Sequence_>_0.7) + (Main_Sequence_>_0.7 -> Main_Sequence_>_0.7)
    1: Stars merged: (Main_Sequence_>_0.7 -> Main_Sequence_>_0.7) + (Main_Sequence_>_0.7 -> Main_Sequence_>_0.7)
    2: Double White Dwarf: (Main_Sequence_>_0.7 -> Carbon-Oxygen_White_Dwarf) + (Main_Sequence_>_0.7 -> Helium_White_Dwarf)
    3: Stars merged: (Main_Sequence_>_0.7 -> Naked_Helium_Star_MS) + (Main_Sequence_<=_0.7 -> Main_Sequence_<=_0.7)
    4: Unbound binary: (Main_Sequence_>_0.7 -> Neutron_Star) + (Main_Sequence_>_0.7 -> Main_Sequence_>_0.7)
    5: Unbound binary: (Main_Sequence_>_0.7 -> Naked_Helium_Star_MS) + (Main_Sequence_>_0.7 -> Neutron_Star)
    6: Unbound binary: (Main_Sequence_>_0.7 -> Neutron_Star) + (Main_Sequence_>_0.7 -> Main_Sequence_>_0.7)
    7: Unbound binary: (Main_Sequence_>_0.7 -> Neutron_Star) + (Main_Sequence_>_0.7 -> Core_Helium_Burning)
    8: Allowed time exceeded: (Main_Sequence_>_0.7 -> Carbon-Oxygen_White_Dwarf) + (Main_Sequence_>_0.7 -> Neutron_Star)
    9: Unbound binary: (Main_Sequence_>_0.7 -> Black_Hole) + (Main_Sequence_>_0.7 -> Main_Sequence_>_0.7)

    Generated 10 of 10 binaries requested

    Simulation completed

    End generating binaries at Tue Sep  7 17:34:49 2021

    Clock time = 0.078125 CPU seconds
    Wall time  = 0000:00:00 (hhhh:mm:ss)


To run COMPAS and evolve five binary system with somespecific initial conditions and evolutionary parameters (default values for the remainder), type::

    ./compas --number-of-systems 5 --initial-mass-1 8.5 --initial-mass-2 13.7 --metallicity 0.015 --mass-loss-prescription VINK --common-envelope-alpha 0.8 --common-envelope-lambda 0.2

This should produce an output put similar to::

    COMPAS v02.22.00
    Compact Object Mergers: Population Astrophysics and Statistics
    by Team COMPAS (http://compas.science/index.html)
    A binary star simulator

    Start generating binaries at Tue Sep  7 17:49:40 2021

    0: Unbound binary: (Main_Sequence_>_0.7 -> Main_Sequence_>_0.7) + (Main_Sequence_>_0.7 -> Neutron_Star)
    1: Unbound binary: (Main_Sequence_>_0.7 -> Main_Sequence_>_0.7) + (Main_Sequence_>_0.7 -> Neutron_Star)
    2: Stars merged: (Main_Sequence_>_0.7 -> Main_Sequence_>_0.7) + (Main_Sequence_>_0.7 -> Naked_Helium_Star_MS)
    3: Unbound binary: (Main_Sequence_>_0.7 -> Main_Sequence_>_0.7) + (Main_Sequence_>_0.7 -> Neutron_Star)
    4: Unbound binary: (Main_Sequence_>_0.7 -> Main_Sequence_>_0.7) + (Main_Sequence_>_0.7 -> Neutron_Star)

    Generated 5 of 5 binaries requested

    Simulation completed

    End generating binaries at Tue Sep  7 17:49:40 2021

    Clock time = 0.0625 CPU seconds
    Wall time  = 0000:00:00 (hhhh:mm:ss)

