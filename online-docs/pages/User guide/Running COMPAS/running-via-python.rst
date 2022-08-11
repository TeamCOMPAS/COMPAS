Running COMPAS via Python
=========================

A convenient method of managing the many program options provided by COMPAS is to run COMPAS via Python, using a script to manage and 
specify the values of the program options.

An example Python script is provided in the COMPAS suite on github: ``runSubmit.py``. Additionally, the default COMPAS options are specified on ``compasConfigDefault.yaml``. Users should copy the ``runSubmit.py`` and ``runSubmit.py`` scripts and modify the ``compasConfigDefault.yaml`` copy to match their experimental requirements. Refer to the :doc:`Getting started guide <../../Getting started/getting-started>` for more details.

To run COMPAS via Python using the ``runSubmit.py`` script provided, set the shell environment variable ``COMPAS-ROOT-DIR``
to the parent directory of the directory in which the COMPAS executable resides, then type `python /path-to-runSubmit/runSubmit.py`. 
For example, for Ubuntu Linux, type::

    export COMPAS_ROOT_DIR=/path-to-dir

    python /path-to-runSubmit/runSubmit.py

This should produce an output similar to::

    python_version = 3
    ../src/COMPAS --enable-warnings False --use-mass-loss  --mass-transfer  --detailed-output False --evolve-unbound-systems False --population-data-printing False --rlof-printing  --circularise-binary-during-mass-transfer  --angular-momentum-conservation-during-circularisation False --pair-instability-supernovae  --pulsational-pair-instability  --quiet False --common-envelope-allow-main-sequence-survive  --evolve-pulsars False --debug-to-file False --errors-to-file False --allow-rlof-at-birth  --allow-touching-at-birth False --store-input-files  --switch-log False --check-photon-tiring-limit False --number-of-systems 10 --metallicity 0.0142 --common-envelope-alpha 1.0 --common-envelope-lambda 0.1 --common-envelope-slope-kruckow -0.8333333333333334 --common-envelope-alpha-thermal 1.0 --common-envelope-lambda-multiplier 1.0 --luminous-blue-variable-multiplier 1.5 --overall-wind-mass-loss-multiplier 1.0 --wolf-rayet-multiplier 1.0 --cool-wind-mass-loss-multiplier 1.0 --mass-transfer-fa 0.5 --mass-transfer-jloss 1.0 --maximum-evolution-time 13700.0 --maximum-number-timestep-iterations 99999 --timestep-multiplier 1.0 --initial-mass-min 5.0 --initial-mass-max 150.0 --initial-mass-power 0.0 --semi-major-axis-min 0.01 --semi-major-axis-max 1000.0 --mass-ratio-min 0.01 --mass-ratio-max 1.0 --minimum-secondary-mass 0.1 --eccentricity-min 0.0 --eccentricity-max 1.0 --metallicity-min 0.0001 --metallicity-max 0.03 --pulsar-birth-magnetic-field-distribution-min 11.0 --pulsar-birth-magnetic-field-distribution-max 13.0 --pulsar-birth-spin-period-distribution-min 10.0 --pulsar-birth-spin-period-distribution-max 100.0 --pulsar-magnetic-field-decay-timescale 1000.0 --pulsar-magnetic-field-decay-massscale 0.025 --pulsar-minimum-magnetic-field 8.0 --orbital-period-min 1.1 --orbital-period-max 1000 --kick-magnitude-sigma-CCSN-NS 265.0 --kick-magnitude-sigma-CCSN-BH 265.0 --fix-dimensionless-kick-magnitude -1 --kick-direction-power 0.0 --random-seed 0 --mass-transfer-thermal-limit-C 10.0 --eddington-accretion-factor 1 --pisn-lower-limit 60.0 --pisn-upper-limit 135.0 --ppi-lower-limit 35.0 --ppi-upper-limit 60.0 --maximum-neutron-star-mass 2.5 --kick-magnitude-sigma-ECSN 30.0 --kick-magnitude-sigma-USSN 30.0 --kick-scaling-factor 1.0 --maximum-mass-donor-nandez-ivanova 2.0 --common-envelope-recombination-energy-density 15000000000000.0 --common-envelope-mass-accretion-max 0.1 --common-envelope-mass-accretion-min 0.04 --zeta-main-sequence 2.0 --zeta-radiative-envelope-giant 6.5 --kick-magnitude-max -1.0 --muller-mandel-kick-multiplier-BH 200.0 --muller-mandel-kick-multiplier-NS 400.0 --log-level 0 --debug-level 0 --hdf5-chunk-size 100000 --hdf5-buffer-size 1 --neutrino-mass-loss-BH-formation-value 0.1 --mode BSE --case-BB-stability-prescription ALWAYS_STABLE --chemically-homogeneous-evolution PESSIMISTIC --luminous-blue-variable-prescription HURLEY_ADD --mass-loss-prescription VINK --mass-transfer-angular-momentum-loss-prescription ISOTROPIC --mass-transfer-accretion-efficiency-prescription THERMAL --mass-transfer-rejuvenation-prescription STARTRACK --initial-mass-function KROUPA --semi-major-axis-distribution FLATINLOG --orbital-period-distribution FLATINLOG --mass-ratio-distribution FLAT --eccentricity-distribution ZERO --metallicity-distribution ZSOLAR --rotational-velocity-distribution ZERO --remnant-mass-prescription FRYER2012 --fryer-supernova-engine DELAYED --black-hole-kicks FALLBACK --kick-magnitude-distribution MAXWELLIAN --kick-direction ISOTROPIC --output-path /d/Jeff/User_Files/compas/dev/my_fork/compas/src --common-envelope-lambda-prescription LAMBDA_NANJING --stellar-zeta-prescription SOBERMAN --mass-transfer-thermal-limit-accretor CFACTOR --pulsational-pair-instability-prescription MARCHANT --neutron-star-equation-of-state SSE --pulsar-birth-magnetic-field-distribution ZERO --pulsar-birth-spin-period-distribution ZERO --common-envelope-mass-accretion-prescription ZERO --envelope-state-prescription LEGACY --logfile-type HDF5 --neutrino-mass-loss-BH-formation FIXED_MASS

    COMPAS v02.22.00
    Compact Object Mergers: Population Astrophysics and Statistics
    by Team COMPAS (http://compas.science/index.html)
    A binary star simulator

    Start generating binaries at Tue Sep  7 18:14:40 2021

    0: Stars merged: (Main_Sequence_>_0.7 -> Main_Sequence_>_0.7) + (Main_Sequence_>_0.7 -> Main_Sequence_>_0.7)
    1: Unbound binary: (Main_Sequence_>_0.7 -> Neutron_Star) + (Main_Sequence_>_0.7 -> Thermally_Pulsing_Asymptotic_Giant_Branch)
    2: Unbound binary: (Main_Sequence_>_0.7 -> Black_Hole) + (Main_Sequence_>_0.7 -> Main_Sequence_>_0.7)
    3: Unbound binary: (Main_Sequence_>_0.7 -> Neutron_Star) + (Main_Sequence_>_0.7 -> Main_Sequence_>_0.7)
    4: Unbound binary: (Main_Sequence_>_0.7 -> Black_Hole) + (Main_Sequence_>_0.7 -> Main_Sequence_>_0.7)
    5: Allowed time exceeded: (Main_Sequence_>_0.7 -> Carbon-Oxygen_White_Dwarf) + (Main_Sequence_<=_0.7 -> Main_Sequence_<=_0.7)
    6: Unbound binary: (Main_Sequence_>_0.7 -> Neutron_Star) + (Main_Sequence_>_0.7 -> Main_Sequence_>_0.7)
    7: Double White Dwarf: (Main_Sequence_>_0.7 -> Carbon-Oxygen_White_Dwarf) + (Main_Sequence_>_0.7 -> Carbon-Oxygen_White_Dwarf)
    8: Unbound binary: (Main_Sequence_>_0.7 -> Black_Hole) + (Main_Sequence_<=_0.7 -> Main_Sequence_<=_0.7)
    9: Stars merged: (Main_Sequence_>_0.7 -> Naked_Helium_Star_Giant_Branch) + (Main_Sequence_>_0.7 -> Main_Sequence_>_0.7)

    Generated 10 of 10 binaries requested

    Simulation completed

    End generating binaries at Tue Sep  7 18:14:40 2021

    Clock time = 0.109375 CPU seconds
    Wall time  = 0000:00:00 (hhhh:mm:ss)

Note that Python prints the Python version, the executes the command to run COMPAS.  The command exceuted is echoed to the stdout.  COMPAS
then runs and produces its usual output.

When using Python and a script file (such as `runSubmit.py`) to run COMPAS, care must be taken to specify program option values correctly in the ``compasConfigDefault.yaml`` file.
For example, ranges and sets can be specified for options in the ``compasConfigDefault.yaml`` file, but the range or set parameter must be enclosed in quotes – 
otherwise python tries to parse the construct. For example, to specify a set of metallicity values in the Python script file, use::

    metallicity = 's[0.001,0.002,0.003,0.007,0.01,0.015,0.02]'

If the set parameter is not enclosed in quotes, Python will attempt to parse it, and will fail.
