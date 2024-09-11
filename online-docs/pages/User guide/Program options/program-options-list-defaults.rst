Program option list and default values
======================================

Any program options that are not specified take default values.

- On the command line, program options that are not explicitly specified default to the COMPAS default value for the option (as specified in the COMPAS code - may be sampled from a distribution).

- On a :doc:`grid file <../grid-files>` line, program options that are not explicitly specified default to the value specified for that option on the command line. If the program option was not explicitly specified on the command line, it will default to the COMPAS default value for the option, as described above. That is, the value for any option not specified on a grid file line option falls back to the value specified on the command line, which falls back to the COMPAS default if it was not specified on the command line.


.. _options-props-top:

The full list of program options with brief explanations and their default values is shown below.  We also include a listing of options (this time, by name only) grouped by category.

**Alphabetical listing: jump to**
:ref:`A <options-props-A>` :ref:`B <options-props-B>` :ref:`C <options-props-C>` :ref:`D <options-props-D>`
:ref:`E <options-props-E>` :ref:`F <options-props-F>` :ref:`G <options-props-G>` :ref:`H <options-props-H>`
:ref:`I <options-props-I>` :ref:`J <options-props-J>` :ref:`K <options-props-K>` :ref:`L <options-props-L>`
:ref:`M <options-props-M>` :ref:`N <options-props-N>` :ref:`O <options-props-O>` :ref:`P <options-props-P>`
:ref:`Q <options-props-Q>` :ref:`R <options-props-R>` :ref:`S <options-props-S>` :ref:`T <options-props-T>`
:ref:`U <options-props-U>` :ref:`V <options-props-V>` :ref:`W <options-props-W>` :ref:`X <options-props-X>`
:ref:`Y <options-props-Y>` :ref:`Z <options-props-Z>`

**Category listing: jump to**
:ref:`Initial conditions <options-initial-conditions>`
:ref:`Stellar evolution and winds <options-stellar-evolution>`
:ref:`Mass transfer physics <options-mass-transfer>`
:ref:`Supernovae <options-supernovae>`
:ref:`Administrative <options-admin>`

COMPAS information
------------------

**--help [ -h ]** |br|
Prints COMPAS help.

**--version [ -v ]** |br|
Prints COMPAS version string.


Alphabetical listing
--------------------

.. _options-props-A:

**--add-options-to-sysparms** |br|
Add columns for program options to SSE System Parameters/BSE System Parameters file (mode dependent). |br|
Options: { ALWAYS, GRID, NEVER } |br|
Default = GRID

.. list-table::
   :widths: 11 80 
   :header-rows: 0
   :class: aligned-text

   * - ALWAYS
     - indicates that the program options should be added to the sysparms file
   * - GRID
     - indicates that the program options should be added to the sysparms file `only if`
   * -  
     - a GRID file is specified, or RANGEs or SETs are specified for options
   * - NEVER
     - indicates that the program options should `not` be added to the sysparms file

**--allow-non-stripped-ECSN** |br|
Allow ECSNe in effectively single progenitors. |br|
Default = FALSE

**--allow-rlof-at-birth** |br|
Allow binaries that have one or both stars in RLOF at birth to evolve as over-contact systems. |br|
Default = TRUE

**--allow-touching-at-birth** |br|
Allow binaries that are touching at birth to be included in the sampling. |br|
Default = FALSE

**--angular-momentum-conservation-during-circularisation** |br|
Conserve angular momentum when binary is circularised when entering a Mass Transfer episode. |br|
Default = FALSE

.. _options-props-B:

:ref:`Back to Top <options-props-top>`

**--black-hole-kicks** |br|
Black hole kicks relative to NS kicks (not relevant for `MANDELMUELLER` ``--remnant-mass-prescription``). |br|
Options: { FULL, REDUCED, ZERO, FALLBACK } |br|
Default = FALLBACK |br|
DEPRECATION NOTICE: this option has been deprecated and will soon be removed. Please use ``--black-hole-kicks-mode`` in future.

**--black-hole-kicks-mode** |br|
Black hole kicks relative to NS kicks (not relevant for `MANDELMUELLER` ``--remnant-mass-prescription``). |br|
Options: { FULL, REDUCED, ZERO, FALLBACK } |br|
Default = FALLBACK

.. _options-props-C:

:ref:`Back to Top <options-props-top>`

**--case-bb-stability-prescription** |br|
Prescription for the stability of case BB/BC mass transfer. |br|
Options: { ALWAYS_STABLE, ALWAYS_STABLE_ONTO_NSBH, TREAT_AS_OTHER_MT, ALWAYS_UNSTABLE } |br|
Case BB mass transfer is treated as always stable, always stable only for mass transfer onto neutron stars or black holes, with stability as determined for all other mass transfer, or always unstable, respectively |br|
Default = ALWAYS_STABLE

**--check-photon-tiring-limit** |br|
Check the photon tiring limit is not exceeded during mass loss. |br|
Default = FALSE

**--chemically-homogeneous-evolution** |br|
Chemically Homogeneous Evolution mode. See :cite:`Riley2021` for details of the implementation
of Chemically Homogeneous Evolution in COMPAS |br|
Options: { NONE, OPTIMISTIC, PESSIMISTIC } |br|
Default = PESSIMISTIC |br|
DEPRECATION NOTICE: this option has been deprecated and will soon be removed. Please use ``--chemically-homogeneous-evolution-mode`` in future.

**--chemically-homogeneous-evolution-mode** |br|
Chemically Homogeneous Evolution mode. See :cite:`Riley2021` for details of the implementation
of Chemically Homogeneous Evolution in COMPAS |br|
Options: { NONE, OPTIMISTIC, PESSIMISTIC } |br|
Default = PESSIMISTIC |br|

**--circulariseBinaryDuringMassTransfer** |br|
Circularise binary when it enters a Mass Transfer episode. |br|
Default = TRUE

**--common-envelope-allow-immediate-RLOF-post-CE-survive** |br|
Allow binaries that experience Roche lobe overflow immediately at the end of the CE phase to survive. |br|
Default = FALSE

**--common-envelope-allow-main-sequence-survive** |br|
Allow main sequence accretors to survive common envelope evolution if other criteria point to survival. |br|
Default = TRUE

**--common-envelope-allow-radiative-envelope-survive** |br| 
Allow binaries with an evolved component with a radiative envelope to survive the common envelope phase (they always survive in the 
`--common-envelope-formalism TWO_STAGE` option). |br|
Default = FALSE

**--common-envelope-alpha** |br|
Common Envelope efficiency alpha. |br|
Default = 1.0

**--common-envelope-alpha-thermal** |br|
Thermal energy contribution to the total envelope binding energy. |br|
Defined such that :math:`\lambda = \alpha_{th} \times \lambda_{b} + (1.0 - \alpha_{th}) \times \lambda_{g}`. |br|
Default = 1.0

**--common-envelope-formalism** |br|
CE formalism prescription. |br|
Options: { ENERGY, TWO_STAGE } |br|
``ENERGY`` is the standard alpha-lambda formalism; ``TWO_STAGE`` is the formalism of Hirai & Mandel (2022) -- the latter always allows radiative-envelope 
donors to survive CE, so `--common-envelope-allow-radiative-envelope-survive` option is ignored |br| 
Default = ENERGY

**--common-envelope-lambda** |br|
Common Envelope lambda. |br|
Only used when ``--common-envelope-lambda-prescription = LAMBDA_FIXED``. |br|
Default = 0.1

**--common-envelope-lambda-multiplier** |br|
Multiplicative constant to be applied to the common envelope lambda parameter for any prescription. |br|
Default = 1.0

**--common-envelope-lambda-nanjing-enhanced** |br|
Continuous extrapolation beyond maximum radius range in Nanjing lambda's as implemented in StarTrack. Only used when ``--common-envelope-lambda-prescription = LAMBDA_NANJING``. |br|
Default = FALSE

**--common-envelope-lambda-nanjing-interpolate-in-mass** |br|
Interpolate Nanjing lambda parameters across different mass models. Only used when ``--common-envelope-lambda-prescription = LAMBDA_NANJING``. |br|
Default = FALSE

**--common-envelope-lambda-nanjing-interpolate-in-metallicity** |br|
Interpolate Nanjing lambda parameters across population I and population II metallicity models. Only used when ``--common-envelope-lambda-prescription = LAMBDA_NANJING``. |br|
Default = FALSE

**--common-envelope-lambda-nanjing-use_rejuvenated-mass** |br|
Use rejuvenated or effective ZAMS mass instead of true birth mass when computing Nanjing lambda parameters. Only used when ``--common-envelope-lambda-prescription = LAMBDA_NANJING``. |br|
Default = FALSE

**--common-envelope-lambda-prescription** |br|
CE lambda (envelope binding energy) prescription. |br|
Options: { LAMBDA_FIXED, LAMBDA_LOVERIDGE, LAMBDA_NANJING, LAMBDA_KRUCKOW, LAMBDA_DEWI } |br|
``LAMBDA_FIXED`` is a constant; ``LAMBDA_LOVERIDGE`` is the prescription from Loveridge et al., 2011; ``LAMBDA_NANJING`` is from Xu & Li, 2010; ``LAMBDA_KRUCKOW`` is from Kruckow et al., 2016; and ``LAMBDA_DEWI`` is the fit from Appendix A in Claeys et al. 2014, based on Dewi & Tauris 2000 |br|
Default = LAMBDA_NANJING

**--common-envelope-mass-accretion-constant** |br|
Value of mass accreted by NS/BH during common envelope evolution if assuming all NS/BH accrete same amount of mass. |br|
Used when ``--common-envelope-mass-accretion-prescription = CONSTANT``, ignored otherwise. |br|
Default = 0.0

**--common-envelope-mass-accretion-max** |br|
Maximum amount of mass accreted by NS/BHs during common envelope evolution (:math:`M_\odot`). |br|
Default = 0.1

**--common-envelope-mass-accretion-min** |br|
Minimum amount of mass accreted by NS/BHs during common envelope evolution (:math:`M_\odot`). |br|
Default = 0.04

**--common-envelope-mass-accretion-prescription** |br|
Assumption about whether NS/BHs can accrete mass during common envelope evolution. |br|
``ZERO`` is no accretion; ``CONSTANT`` means a fixed amount of accretion determined by ``--common-envelope-mass-accretion-constant``; ``UNIFORM`` means a uniform random draw between ``--common-envelope-mass-accretion-min`` and ``--common-envelope-mass-accretion-max`` (Oslowski et al., 2011);, ``MACLEOD`` follows the prescription of MacLeod et al., 2015, and ``CHEVALIER`` follows the accretion assumptions in Chevalier et al. 1993 as in Model 2 from van Son et al. 2020 |br|
Options: { ZERO, CONSTANT, UNIFORM, MACLEOD, CHEVALIER } |br|
Default = ZERO

**--common-envelope-recombination-energy-density** |br|
Recombination energy density (erg g−1). |br|
Default = :math:`1.5 \times 10^{13}`

**--common-envelope-slope-kruckow** |br|
Slope for the Kruckow lambda (see Kruckow et al. 2016 as implemented by Vigna-Gomez et al. 2018). |br|
Default = −0.833333

**--convective-envelope-temperature-threshold** |br|
Temperature [K] threshold, below which the envelopes of giants are convective. 
Only used for --envelope-state-prescription = FIXED_TEMPERATURE, ignored otherwise. |br|
Default = 5370

**--cool-wind-mass-loss-multiplier** |br|
Multiplicative constant for wind mass loss of cool stars, i.e. those with temperatures below the
VINK_MASS_LOSS_MINIMUM_TEMP (default 12500K). |br|
Only applicable when ``--mass-loss-prescription = VINK``. |br|
Default = 1.0

**--create-YAML-file** |br|
Creates new YAML file.  Argument is filename for new YAML file. |br|
Default = None - name must be supplied if option is present.

**--critical-mass-ratio-HG-degenerate-accretor** |br|
Critical mass ratio for MT from a HG star to a degenerate accretor. |br|
0 is always stable, < 0 is disabled. |br|
Only used for ``--critical-mass-ratio-prescription CLAEYS``, ignored otherwise. |br|
Default = 0.21

**--critical-mass-ratio-HG-non-degenerate-accretor** |br|
Critical mass ratio for MT from a HG star to a non-degenerate accretor. |br|
0 is always stable, < 0 is disabled. |br|
Only used for ``--critical-mass-ratio-prescription CLAEYS``, ignored otherwise. |br|
Default = 0.25

**--critical-mass-ratio-MS-high-mass-degenerate-accretor** |br|
Critical mass ratio for MT from a MS star to a degenerate accretor. |br|
0 is always stable, < 0 is disabled. |br|
Only used for ``--critical-mass-ratio-prescription CLAEYS``, ignored otherwise. |br|
Default = 0.0

**--critical-mass-ratio-MS-high-mass-non-degenerate-accretor** |br|
Critical mass ratio for MT from a MS star to a non-degenerate accretor. |br|
0 is always stable, < 0 is disabled. |br|
Only used for ``--critical-mass-ratio-prescription CLAEYS``, ignored otherwise. |br|
Default = 0.625

**--critical-mass-ratio-MS-low-mass-degenerate-accretor** |br|
Critical mass ratio for MT from a MS star to a degenerate accretor. |br|
0 is always stable, < 0 is disabled. |br|
Only used for ``--critical-mass-ratio-prescription CLAEYS``, ignored otherwise. |br|
Default = 1.0

**--critical-mass-ratio-MS-low-mass-non-degenerate-accretor** |br|
Critical mass ratio for MT from a MS star to a non-degenerate accretor. |br|
0 is always stable, < 0 is disabled. |br|
Only used for ``--critical-mass-ratio-prescription CLAEYS``, ignored otherwise. |br|
Default = 1.44

**--critical-mass-ratio-giant-degenerate-accretor** |br|
Critical mass ratio for MT from a giant star to a degenerate accretor. |br|
0 is always stable, < 0 is disabled. |br|
Only used for ``--critical-mass-ratio-prescription CLAEYS``, ignored otherwise. |br|
Default = 0.87

**--critical-mass-ratio-giant-non-degenerate-accretor** |br|
Critical mass ratio for MT from a giant star to a non-degenerate accretor. |br|
0 is always stable, < 0 is disabled. |br|
Only used for ``--critical-mass-ratio-prescription CLAEYS``, ignored otherwise. |br|
Default shows -1, but this translates to a function of the core mass ratio, as described in Claeys+ 2014. 

**--critical-mass-ratio-helium-HG-degenerate-accretor** |br|
Critical mass ratio for MT from a helium HG star to a degenerate accretor. |br|
0 is always stable, < 0 is disabled. |br|
Only used for ``--critical-mass-ratio-prescription CLAEYS``, ignored otherwise. |br|
Default = 0.21

**--critical-mass-ratio-helium-HG-non-degenerate-accretor** |br|
Critical mass ratio for MT from a helium HG star to a non-degenerate accretor. |br|
0 is always stable, < 0 is disabled. |br|
Only used for ``--critical-mass-ratio-prescription CLAEYS``, ignored otherwise. |br|
Default = 0.25

**--critical-mass-ratio-helium-MS-degenerate-accretor** |br|
Critical mass ratio for MT from a helium MS star to a degenerate accretor. |br|
0 is always stable, < 0 is disabled. |br|
Only used for ``--critical-mass-ratio-prescription CLAEYS``, ignored otherwise. |br|
Default = 0.0

**--critical-mass-ratio-helium-MS-non-degenerate-accretor** |br|
Critical mass ratio for MT from a helium MS star to a non-degenerate accretor. |br|
0 is always stable, < 0 is disabled. |br|
Only used for ``--critical-mass-ratio-prescription CLAEYS``, ignored otherwise. |br|
Default = 0.0

**--critical-mass-ratio-helium-giant-degenerate-accretor** |br|
Critical mass ratio for MT from a helium giant star to a degenerate accretor. |br|
0 is always stable, < 0 is disabled. |br|
Only used for ``--critical-mass-ratio-prescription CLAEYS``, ignored otherwise. |br|
Default = 0.87

**--critical-mass-ratio-helium-giant-non-degenerate-accretor** |br|
Critical mass ratio for MT from a helium giant star to a non-degenerate accretor. |br|
0 is always stable, < 0 is disabled. |br|
Only used for ``--critical-mass-ratio-prescription CLAEYS``, ignored otherwise. |br|
Default = 1.28

**--critical-mass-ratio-prescription** |br|
Critical mass ratio stability prescription (if any). |br|
Options: { NONE, ZERO, CLAEYS, GE20, GE20_IC, HURLEY_HJELLMING_WEBBINK } |br|
``NONE``    defaults to the zeta prescription for stability. |br|
``CLAEYS``  uses qCrit values from Claeys et al. 2014. |br|
``GE20``    uses qCrit values from Ge et al. 2020 (adiabatic assumption). |br|
``GE20_IC`` uses qCrit values from Ge et al. 2020 (isentropic envelope assumption). |br|
``HURLEY_HJELLMING_WEBBINK`` uses qCrit values from Hurley et al. 2002 (Hjellming & Webbink 1987 for mass transfer from a giant primary). |br|
Warning: if running with ``--critical-mass-ratio-prescription``, zetas will not be computed, so should not be trusted in the outputs. |br|
Default = NONE |br|

**--critical-mass-ratio-white-dwarf-degenerate-accretor** |br|
Critical mass ratio for MT from a white dwarf to a degenerate accretor. |br|
0 is always stable, < 0 is disabled. |br|
Only used for ``--critical-mass-ratio-prescription CLAEYS``, ignored otherwise. |br|
Default = 1.6

**--critical-mass-ratio-white-dwarf-non-degenerate-accretor** |br|
Critical mass ratio for MT from a white dwarf to a non-degenerate accretor. |br|
0 is always stable, < 0 is disabled. |br|
Only used for ``--critical-mass-ratio-prescription CLAEYS``, ignored otherwise. |br|
Default = 0.0

.. _options-props-D:

:ref:`Back to Top <options-props-top>`

**--debug-classes** |br|
Developer-defined debug classes to enable (vector). |br|
Default = `All debug classes enabled (e.g. no filtering)`

**--debug-level** |br|
Determines which print statements are displayed for debugging. |br|
Default = 0

**--debug-to-file** |br|
Write debug statements to file. |br|
Default = FALSE

**--detailed-output** |br|
Print BSE detailed information to file. |br|
Default = FALSE

.. _options-props-E:

:ref:`Back to Top <options-props-top>`

**--eccentricity [ -e ]** |br|
Initial eccentricity for a binary star when evolving in BSE mode.
Default = 0.0 |br|

**--eccentricity-distribution** |br|
Initial eccentricity distribution. |br|
Options: { ZERO, FLAT, GELLER+2013, THERMAL, DUQUENNOYMAYOR1991, SANA2012 } |br|
``ZERO`` always circular, ``FLAT`` is uniform on [``--eccentricity-min``,``--eccentricity-max``], ``THERMAL`` is p(e) proportional to e, and the other options refer to the distributions of Geller et al. 2013, Duqennoy & Mayor 1991, and Sana et al. 2012. |br|
Default = ZERO

**--eccentricity-max** |br|
Maximum eccentricity to generate. |br|
Default = 1.0

**--eccentricity-min** |br|
Minimum eccentricity to generate. |br|
Default = 0.0

**--eddington-accretion-factor** |br|
Multiplication factor for Eddington accretion for NS & BH (i.e. > 1 is super-eddington and 0 is no accretion). |br|
Default = 1.0

**--emit-gravitational-radiation**  |br|
Emit gravitational radiation at each timestep of binary evolution according to Peters 1964. |br|
Default = FALSE

**--enable-warnings** |br|
Display warning messages to stdout. |br|
Default = FALSE

**--enable-rotationally-enhanced-mass-loss** |br|
Enable rotationally enhanced mass loss for rapidly rotating stars following Langer (1998)
Default = FALSE

**--enhance-CHE-lifetimes-luminosities** |br|
Enhance lifetimes and luminosities of CH stars using a fit to detailed models from Szecsi et al. (2015)
Default = FALSE

**--envelope-state-prescription** |br|
Prescription for determining whether the envelope of the star is convective or radiative. |br|
Options: { LEGACY, HURLEY, FIXED_TEMPERATURE } |br|
``LEGACY`` refers to the model used in Stevenson et al., 2017 |br|
``HURLEY`` refers to the model of Hurley, Pols, Tout, 2002 |br|
``FIXED_TEMPERATURE`` assumes that a deep convective envelope developes only when the temperature drops below ``CONVECTIVE_BOUNDARY_TEMPERATURE`` (Klencki et al., 2020) |br|
Default = LEGACY

**--errors-to-file** |br|
Write error messages to file. |br|
Default = FALSE

**--expel-convective-envelope-above-luminosity-threshold** |br|
Expel convective envelope in a pulsation if the luminosity to mass ratio exceeds the threshold given by ``--luminosity-to-mass-threshold`` |br|
Default = FALSE

**--evolve-double-white-dwarfs** |br|
Continue evolving double white dwarf systems after their formation. |br|
Default = FALSE

**--evolve-main-sequence-mergers** |br|
Continue evolving the remnant after a main sequence merger. |br|
Default = FALSE

**--evolve-pulsars** |br|
Evolve pulsar properties of Neutron Stars. |br|
Default = FALSE

**--evolve-unbound-systems** |br|
Continue evolving stars even if the binary is disrupted. |br|
Default = TRUE

.. _options-props-F:

:ref:`Back to Top <options-props-top>`

**--fix-dimensionless-kick-magnitude** |br|
Fix dimensionless kick magnitude to this value. |br|
Default = n/a (not used if option not present)

**--fp-error-mode** |br|
Specifies the floating-point error handling mode. |br|
Options: { OFF, ON, DEBUG } |br|
Default = OFF |br|
Refer to :doc:`../Handling errors/handling-errors` for a discussion of the option values.

**--fryer-supernova-engine** |br|
Supernova engine type if using the remnant mass prescription from :cite:`Fryer2012`. |br|
Options: { DELAYED, RAPID } |br|
Default = DELAYED

**--fryer-22-fmix** |br|
Parameter describing the mixing growth time when using the 'FRYER2022' remnant mass distribution  :cite:`Fryer2022`. |br|
Default = 0.5, which is closest to the 'DELAYED' remnant mass prescription from :cite:`Fryer2012`. |br|
A value of 4.0 is closest to  the 'RAPID' remnant mass prescription from :cite:`Fryer2012`. |br|
If the FALLBACK option is used for the kicks, then the proto core masses will be determined by the fryer-supernova-engine option.

**--fryer-22-mcrit** |br|
Critical CO core mass for black hole formation when using the 'FRYER2022' remnant mass distribution :cite:`Fryer2022`. |br|
Default = 5.75



.. _options-props-G:

:ref:`Back to Top <options-props-top>`

**--grid** |br|
Grid filename. (See :doc:`Grid files <../grid-files>`) |br|
Default = ’’ (None)

**--grid-lines-to-process** |br|
The number of grid file lines to be processed. |br|
Default = Process to EOF

**--grid-start-line** |br|
The first line of the grid file to be processed. |br|
Default = 0

.. _options-props-H:

:ref:`Back to Top <options-props-top>`

**--hdf5-buffer-size** |br|
The ``HDF5`` IO buffer size for writing to ``HDF5`` logfiles (number of ``HDF5`` chunks). |br|
Default = 1

**--hdf5-chunk-size** |br|
The ``HDF5`` dataset chunk size to be used when creating ``HDF5`` logfiles (number of logfile entries). |br|
Default = 100000

**--help [ -h ]** |br|
Prints COMPAS help (-h is short form, --help includes more information).

.. _options-props-I:

:ref:`Back to Top <options-props-top>`

**--include-WD-binaries-as-DCO** |br|
When enabled, changes the definition of "Double Compact Object" from a binary comprised of any two of {Neutron Star, Black Hole} |br|
to a binary star comprised of any two of {Helium White Dwarf, Carbon-Oxygen White Dwarf, Oxygen-Neon White Dwarf, Neutron Star, Black Hole}. |br|
Default = FALSE

**--initial-mass** |br|
Initial mass for a single star when evolving in SSE mode (:math:`M_\odot`). |br|
Default = Sampled from IMF

**--initial-mass-1** |br|
Initial mass for the primary star when evolving in BSE mode (:math:`M_\odot`). |br|
Default = Sampled from IMF

**--initial-mass-2** |br|
Initial mass for the secondary star when evolving in BSE mode (:math:`M_\odot`). |br|
Default = Sampled from IMF

**--initial-mass-function [ -i ]** |br|
Initial mass function. |br|
Options: { SALPETER, POWERLAW, UNIFORM, KROUPA } |br|
``SALPETER`` and ``KROUPA`` use the IMFs of Salpeter 1955 and Kroupa 2001 |br|
``POWERLAW`` samples from a single power law with slope ``--initial-mass-power`` |br|
``UNIFORM`` samples uniformly between ``--initial-mass-min`` and ``--initial-mass-min`` |br|
Default = KROUPA

**--initial-mass-max** |br|
Maximum mass to generate using given IMF (:math:`M_\odot`). |br|
Default = 150.0

**--initial-mass-min** |br|
Minimum mass to generate using given IMF (:math:`M_\odot`). |br|
Default = 5.0

**--initial-mass-power** |br|
Single power law power to generate primary mass using ``POWERLAW`` IMF. |br|
Default = 0.0

.. _options-props-J:

.. _options-props-K:

:ref:`Back to Top <options-props-top>`

**--kick-direction** |br|
Natal kick direction distribution. |br|
Options: { ISOTROPIC, INPLANE, PERPENDICULAR, POWERLAW, WEDGE, POLES } |br|
Kick angles are defined relative to the spin axis. |br|
``INPLANE`` and ``PERPENDICULAR`` are strictly in the equatorial plane or in polar directions. |br|
``WEDGE`` and ``POLES`` are preferentially but exactly in the equatorial plane or in polar directions with 1 degree scales, respectively. |br|
``POWERLAW`` quantifies the preference for polar vs planar kicks with the ``--kick-direction-power`` parameter. |br|
Default = ISOTROPIC |br|
DEPRECATION NOTICE: this option has been deprecated and will soon be removed. Please use ``--kick-direction-distribution`` in future.

**--kick-direction-distribution** |br|
Natal kick direction distribution. |br|
Options: { ISOTROPIC, INPLANE, PERPENDICULAR, POWERLAW, WEDGE, POLES } |br|
Kick angles are defined relative to the spin axis. |br|
``INPLANE`` and ``PERPENDICULAR`` are strictly in the equatorial plane or in polar directions. |br|
``WEDGE`` and ``POLES`` are preferentially but exactly in the equatorial plane or in polar directions with 1 degree scales, respectively. |br|
``POWERLAW`` quantifies the preference for polar vs planar kicks with the ``--kick-direction-power`` parameter. |br|
Default = ISOTROPIC |br|



**--kick-direction-power** |br|
Power for power law kick direction distribution, where 0.0 = isotropic, +ve = polar, -ve = in plane. |br|
Default = 0.0 (isotropic)

**--kick-magnitude** |br|
Value to be used as the (drawn) kick magnitude for a single star when evolving in SSE mode, should the star
undergo a supernova event (:math:`km s^{−1}`). |br|
If a value for option ``--kick-magnitude-random`` is specified, it will be used in preference to ``--kick-magnitude``. |br|
Default = 0.0

**--kick-magnitude-1** |br|
Value to be used as the (drawn) kick magnitude for the primary star of a binary system when evolving in
BSE mode, should the star undergo a supernova event (:math:`km s^{−1}`). |br|
If a value for option ``--kick-magnitude-random-1`` is specified, it will be used in preference to ``--kick-magnitude-1``. |br|
Default = 0.0

**--kick-magnitude-2** |br|
Value to be used as the (drawn) kick magnitude for the secondary star of a binary system when evolving in
BSE mode, should the star undergo a supernova event (:math:`km s^{−1}`). |br|
If a value for option ``--kick-magnitude-random-2`` is specified, it will be used in preference to ``--kick-magnitude-2``. |br|
Default = 0.0

**--kick-magnitude-distribution** |br|
Natal kick magnitude distribution. |br|
Options: { ZERO, FIXED, FLAT, MAXWELLIAN, BRAYELDRIDGE, MULLER2016, MULLER2016MAXWELLIAN, MULLERMANDEL } |br|
``ZERO`` assigns kick magnitudes of 0.0. |br|
``FIXED`` always sets the magnitude to a fixed value based on supernova type. |br|
``FLAT`` and ``MAXWELLIAN`` draw kicks from uniform or Maxwellian (e.g., Hobbs et al., 2005) distributions, respectively. |br|
``BRAYELDRIDGE`` and ``MULLERMANDEL`` use momentum-preserving kicks from Bray & Eldrigde 2018 and Mandel & Mueller 2020, respectively. |br|
``MULLER2016`` and ``MULLER2016MAXWELLIAN`` use kicks from Mueller 2016 as implemented in Vigna-Gomez et al., 2018 (reduced by a factor of sqrt(3) in the latter case). |br|
Note that this is independent from ``--remnant-mass-prescription`` to provide flexibility; 
however, the ``MULLERMANDEL`` kick prescription is intended to be consistently used with the ``MULLERMANDEL`` remnant mass prescription. |br|
Default = MULLERMANDEL

**--kick-magnitude-max** |br|
Maximum drawn kick magnitude (:math:`km s^{−1}`). |br|
Must be > 0 if using ``--kick-magnitude-distribution = FLAT``. |br|
Default = −1.0

**--kick-magnitude-random** |br|
CDF value to be used to draw the kick magnitude for a single star when evolving in SSE mode, should the star undergo a
supernova event and should the chosen distribution sample from a cumulative distribution function. |br|
Must be a floating-point number in the range :math:`[0.0, 1.0)`. |br|
The specified value for this option will be used in preference to any specified value for ``--kick-magnitude``. |br|
Default = Random number drawn uniformly from :math:`[0.0, 1.0)`

**--kick-magnitude-random-1** |br|
CDF value to be used to draw the kick magnitude for the primary star of a binary system when evolving in BSE mode,
should the star undergo a supernova event and should the chosen distribution sample from a cumulative distribution function. |br|
Must be a floating-point number in the range :math:`[0.0, 1.0)`. |br|
The specified value for this option will be used in preference to any specified value for ``--kick-magnitude-1``. |br|
Default = Random number drawn uniformly from :math:`[0.0, 1.0)`

**--kick-magnitude-random-2** |br|
CDF value to be used to draw the kick magnitude for the secondary star of a binary system when evolving in BSE mode,
should the star undergo a supernova event and should the chosen distribution sample from a cumulative distribution function. |br|
Must be a floating-point number in the range :math:`[0.0, 1.0)`. |br|
The specified value for this option will be used in preference to any specified value for ``--kick-magnitude-2``. |br|
Default = Random number drawn uniformly from :math:`[0.0, 1.0)`

**--kick-magnitude-sigma-CCSN-BH** |br|
Sigma for chosen kick magnitude distribution for black holes (:math:`km s^{−1}`). |br|
Ignored if not needed for the chosen kick magnitude distribution. |br|
Default = 265.0

**--kick-magnitude-sigma-CCSN-NS** |br|
Sigma for chosen kick magnitude distribution for neutron stars (:math:`km s^{−1}`). |br|
Ignored if not needed for the chosen kick magnitude distribution. |br|
Default = 265.0

**--kick-magnitude-sigma-ECSN** |br|
Sigma for chosen kick magnitude distribution for ECSN (:math:`km s^{−1}`). |br|
Ignored if not needed for the chosen kick magnitude distribution. |br|
Default = 30.0

**--kick-magnitude-sigma-USSN** |br|
Sigma for chosen kick magnitude distribution for USSN (:math:`km s^{−1}`). |br|
Ignored if not needed for the chosen kick magnitude distribution. |br|
Default = 30.0

**--kick-mean-anomaly-1** |br|
The mean anomaly at the instant of the supernova for the primary star of a binary system when evolving in
BSE mode, should it undergo a supernova event. |br|
Must be a floating-point number in the range :math:`[0.0, 2\pi)`. |br|
Default = Random number drawn uniformly from :math:`[0.0, 2\pi)`

**--kick-mean-anomaly-2** |br|
The mean anomaly at the instant of the supernova for the secondary star of a binary system when evolving
in BSE mode, should it undergo a supernova event. |br|
Must be a floating-point number in the range :math:`[0.0, 2\pi)`. |br|
Default = Random number drawn uniformly from :math:`[0.0, 2\pi)`

**--kick-phi-1** |br|
The angle between ’x’ and ’y’, both in the orbital plane of the supernova vector, for the primary star of
a binary system when evolving in BSE mode, should it undergo a supernova event (radians). |br|
Default = Drawn according to specified ``--kick-direction`` distribution

**--kick-phi-2** |br|
The angle between ’x’ and ’y’, both in the orbital plane of the supernova vector, for the secondary
star of a binary system when evolving in BSE mode, should it undergo a supernova event (radians). |br|
Default = Drawn according to specified ``--kick-direction`` distribution

**--kick-scaling-factor** |br|
Arbitrary factor used to scale kicks. |br|
Default = 1.0

**--kick-theta-1** |br|
The angle between the orbital plane and the ’z’ axis of the supernova vector for the primary star of a
binary system when evolving in BSE mode, should it undergo a supernova event (radians). |br|
Default = Drawn according to specified ``--kick-direction`` distribution

**--kick-theta-2** |br|
The angle between the orbital plane and the ’z’ axis of the supernova vector for the secondary star of a
binary system when evolving in BSE mode, should it undergo a supernova event (radians). |br|
Default = Drawn according to specified ``--kick-direction`` distribution

.. _options-props-L:

:ref:`Back to Top <options-props-top>`

**--LBV-mass-loss-prescription** |br|
Luminous blue variable mass loss prescription. |br|
Options: { NONE, ZERO, HURLEY, HURLEY_ADD, BELCZYNSKI } |br|
``NONE``       : No LBV winds |br|
``ZERO``       : No LBV winds |br|
``HURLEY``     : Hurley, Pols, Tout (2000) |br|
``HURLEY_ADD`` : ``HURLEY`` in addition to other winds |br|
``BELCZYNSKI`` : Belzcynski et al. 2010 |br|
Default = HURLEY_ADD |br|
DEPRECATION NOTICE: the value of ``NONE`` for this option has been deprecated and will soon be removed. Please use ``ZERO`` in future.

**--log-classes** |br|
Logging classes to be enabled (vector). |br|
Default = `All debug classes enabled (e.g. no filtering)`

**--logfile-common-envelopes** |br|
Filename for Common Envelopes logfile (BSE mode). |br|
Default = ’BSE_Common_Envelopes’

**--logfile-common-envelopes-record-types** |br|
Enabled record types for Common Envelopes logfile (BSE mode). |br|
Default = -1 (all record types) |br|

**--logfile-definitions** |br|
Filename for logfile record definitions file. |br|
Default = ’’ (None)

**--logfile-detailed-output** |br|
Filename for the Detailed Output logfile. |br|
Default = ’SSE_Detailed_Output’ for SSE mode; ’BSE_Detailed_Output’ for BSE mode |br|

**--logfile-detailed-output-record-types** |br|
Enabled record types for the Detailed Output logfile. |br|
Default = -1 (all record types) |br|

**--logfile-double-compact-objects** |br|
Filename for the Double Compact Objects logfile (BSE mode). |br|
Default = ’BSE_Double_Compact_Objects’

**--logfile-double-compact-objects-record-types** |br|
Enabled record types for the Double Compact Objects logfile (BSE mode). |br|
Default = -1 (all record types) |br|

**--logfile-name-prefix** |br|
Prefix for logfile names. |br|
Default = ’’ (None)

**--logfile-pulsar-evolution** |br|
Filename for the Pulsar Evolution logfile (BSE mode). |br|
Default = ’BSE_Pulsar_Evolution’

**--logfile-pulsar-evolution-record-types** |br|
Enabled record types for the Pulsar Evolution logfile (BSE mode). |br|
Default = -1 (all record types) |br|

**--logfile-rlof-parameters** |br|
Filename for the RLOF Printing logfile (BSE mode). |br|
Default = ’BSE_RLOF’

**--logfile-rlof-parameters-record-types** |br|
Enabled record types for the RLOF Printing logfile (BSE mode). |br|
Default = -1 (all record types) |br|

**--logfile-supernovae** |br|
Filename for the Supernovae logfile. |br|
Default = ’SSE_Supernovae’ for SSE mode; ’BSE_Supernovae’ for BSE mode |br|

**--logfile-supernovae-record-types** |br|
Enabled record types for the Supernovae logfile. |br|
Default = -1 (all record types) |br|

**--logfile-switch-log** |br|
Filename for the Switch Log logfile. |br|
Default = ’SSE_Switch_Log’ for SSE mode; ’BSE_Switch_Log’ for BSE mode |br|

**--logfile-system-parameters** |br|
Filename for the System Parameters logfile (BSE mode). |br|
Default = ’SSE_System_Parameters’ for SSE mode; ’BSE_System_Parameters’ for BSE mode |br|

**--logfile-system-parameters-record-types** |br|
Enabled record types for the System Parameters logfile (BSE mode). |br|
Default = -1 (all record types) |br|

**--logfile-type** |br|
The type of logfile to be produced by COMPAS. Options are: HDF5, CSV, TSV, TXT. |br|
Default = ’HDF5’

**--log-level** |br|
Determines which print statements are included in the logfile. |br|
Default = 0

**--luminous-blue-variable-multiplier** |br|
Multiplicative constant for LBV mass loss. (Use 10 for Mennekens & Vanbeveren (2014)). |br|
Note that wind mass loss will also be multiplied by the ``--overall-wind-mass-loss-multiplier``. |br|
Default = 1.5

**--luminous-blue-variable-prescription** |br|
Luminous blue variable mass loss prescription. |br|
Options: { NONE, ZERO, HURLEY, HURLEY_ADD, BELCZYNSKI } |br|
``NONE``       : No LBV winds |br|
``ZERO``       : No LBV winds |br|
``HURLEY``     : Hurley, Pols, Tout (2000) LBV winds only for LBV stars |br|
``HURLEY_ADD`` : ``HURLEY`` in addition to other winds |br|
``BELCZYNSKI`` : Belzcynski et al. 2010 |br|
Default = HURLEY_ADD |br|
DEPRECATION NOTICE: this option has been deprecated and will soon be removed. Please use ``--LBV-mass-loss-prescription`` in future. |br|
DEPRECATION NOTICE: the value of ``NONE`` for this option has been deprecated and will soon be removed. Please use ``ZERO`` in future.

**--luminosity-to-mass-threshold** |br|
Threshold :math:`\log_{10}(Luminosity/Mass)` (in solar units) above which, if the option
``expel-convective-envelope-above-luminosity-threshold`` is set to TRUE, pulsations eject the convective envelope |br|
Default = 4.2

.. _options-props-M:

:ref:`Back to Top <options-props-top>`

**--mass-change-fraction** |br|
Approximate desired fractional change in stellar mass on phase when setting SSE and BSE timesteps (applied before ``--timestep--multiplier``). |br|
Recommended value is 0.005. |br|
A value of 0.0 means that this choice is ignored, and timestep estimates will be calculated by COMPAS. |br|
Default = 0.0

**--mass-loss-prescription** |br|
Mass loss prescription. |br|
Options: { ZERO, NONE, HURLEY, BELCZYNSKI2010, MERRITT2024 } |br|
``ZERO/NONE``     : No wind mass loss |br|
``HURLEY``        : Mass loss according to Hurley (2000) |br|
``BELCZYNSKI2010``: Mass loss as per Belczynski (2010), and the default prior to 2023 |br|
``MERRITT2024``   : Flexible mass loss with phase specific options: (OB, RSG, WR, VMS) |br|
Default = MERRITT2024 |br|
DEPRECATION NOTICE: the value of ``NONE`` for this option has been deprecated and will soon be removed. Please use ``ZERO`` in future.


**--mass-ratio [ -q ]** |br|
Mass ratio :math:`\frac{m2}{m1}` used to determine secondary mass if not specified via ``--initial-mass-2``. |br|
Default: value is sampled if option not specified.

**--mass-ratio-distribution** |br|
Initial mass ratio distribution for :math:`q = \frac{m2}{m1}`. |br|
Options: { FLAT, DUQUENNOYMAYOR1991, SANA2012 } |br|
``FLAT`` is uniform in the mass ratio between ``--mass-ratio-min`` and ``--mass-ratio-max``. |br|
Other prescriptions follow Duquennoy & Mayor 1991 and Sana et al. 2012 |br|
Default = FLAT

**--mass-ratio-max** |br|
Maximum mass ratio :math:`\frac{m2}{m1}` to generate. |br|
Default = 1.0

**--mass-ratio-min** |br|
Minimum mass ratio :math:`\frac{m2}{m1}` to generate. |br|
Default = 0.01

**--mass-transfer** |br|
Enable mass transfer. |br|
Default = TRUE |br|
DEPRECATION NOTICE: this option has been deprecated and will soon be removed. Please use ``--use-mass-transfer`` in future.

**--mass-transfer-accretion-efficiency-prescription** |br|
Mass transfer accretion efficiency prescription. |br|
Options: { THERMAL, FIXED, CENTRIFUGAL } |br|
Default = THERMAL

**--mass-transfer-angular-momentum-loss-prescription** |br|
Mass Transfer Angular Momentum Loss prescription. |br|
Options: { JEANS, ISOTROPIC, CIRCUMBINARY, MACLEOD_LINEAR, ARBITRARY } |br|
Default = ISOTROPIC

**--mass-transfer-fa** |br|
Mass Transfer fraction accreted. |br|
Used when ``--mass-transfer-accretion-efficiency-prescription = FIXED_FRACTION``. |br|
Default = 0.5

**--mass-transfer-jloss** |br|
Specific angular momentum with which the non-accreted system leaves the system. |br|
Used when ``--mass-transfer-angular-momentum-loss-prescription = ARBITRARY``, ignored otherwise. |br|
Default = 1.0

**--mass-transfer-jloss-macleod-linear-fraction-degen** |br|
Specific angular momentum interpolation fraction for degenerate accretors, linear between 0 and 1 corresponding to the accretor and L2 point. |br|
Used when ``--mass-transfer-angular-momentum-loss-prescription = MACLEOD_LINEAR``, ignored otherwise. |br|
Default = 0.5

**--mass-transfer-jloss-macleod-linear-fraction-non-degen** |br|
Specific angular momentum interpolation fraction for non-degenerate accretors, linear between 0 and 1 corresponding to the accretor and L2 point. |br|
Used when ``--mass-transfer-angular-momentum-loss-prescription = MACLEOD_LINEAR``, ignored otherwise. |br|
Default = 0.5

**--mass-transfer-rejuvenation-prescription** |br|
Mass Transfer Rejuvenation prescription. |br|
Options: { HURLEY, STARTRACK } |br|
``HURLEY`` uses the Hurley, Pols, Tout (2000) model. |br|
``STARTRACK`` uses the model from Belczynski et al. 2008 |br|
Default = STARTRACK

**--mass-transfer-thermal-limit-accretor** |br|
Mass Transfer Thermal Accretion limit multiplier. |br|
Options: { CFACTOR, ROCHELOBE } |br|
Default = CFACTOR |br|
DEPRECATION NOTICE: this option has been deprecated and will soon be removed. Please use ``--mass-transfer-thermal-limit-accretor-multiplier`` in future.

**--mass-transfer-thermal-limit-accretor-multiplier** |br|
Mass Transfer Thermal Accretion limit multiplier. |br|
Options: { CFACTOR, ROCHELOBE } |br|
Default = CFACTOR

**--mass-transfer-thermal-limit-C** |br|
Mass Transfer Thermal rate factor for the accretor. |br|
Default = 10.0

**--maximum-evolution-time** |br|
Maximum time to evolve binaries (Myr). Evolution of the binary will stop if this number is reached. |br|
Default = 13700.0

**--maximum-mass-donor-nandez-ivanova** |br|
Maximum donor mass allowed for the revised common envelope formalism of Nandez & Ivanova (:math:`M_\odot`). |br|
Default = 2.0

**--maximum-neutron-star-mass** |br|
Maximum mass of a neutron star (:math:`M_\odot`). |br|
Default = 2.5

**--maximum-number-timestep-iterations** |br|
Maximum number of timesteps to evolve binary. Evolution of the binary will stop if this number is reached. |br|
Default = 99999

**--mcbur1** |br|
Minimum core mass at base of AGB to avoid fully degenerate CO core formation (:math:`M_\odot`). |br|
e.g. 1.6 in :cite:`Hurley2000` presciption; 1.83 in :cite:`Fryer2012` and :doc:`Belczynski et al. (2008) <../../references>` models. |br|
Default = 1.6

**--metallicity [ -z ]** |br|
Metallicity. |br|
The value specified for metallicity is applied to both stars for BSE mode. |br|
Default = 0.0142

**--metallicity-distribution** |br|
Metallicity distribution. |br|
Options: { ZSOLAR, LOGUNIFORM } |br|
``ZSOLAR`` uses ``ZSOL_ASPLUND`` for all initial metallicities. |br|
``LOGUNIFORM`` draws the metallicity uniformly in the log between ``--metallicity-min`` and ``--metallicity-max`` |br|
Default = ZSOLAR

**--metallicity-max** |br|
Maximum metallicity to generate. |br|
Default = 0.03

**--metallicity-min** |br|
Minimum metallicity to generate. |br|
Default = 0.0001

**--minimum-secondary-mass** |br|
Minimum mass of secondary to generate (:math:`M_\odot`). |br|
Defaults to 0.1 if ``--initial-mass-2`` is specified, otherwise value of ``--initial-mass-min``.

**--mode** |br|
The mode of evolution. |br|
Options: { SSE, BSE } |br|
Default = BSE

**--muller-mandel-kick-multiplier-BH** |br|
Scaling prefactor for BH kicks when using the `MULLERMANDEL` kick magnitude distribution |br|
Default = 200.0

**--muller-mandel-kick-multiplier-NS** |br|
Scaling prefactor for NS kicks when using the `MULLERMANDEL` kick magnitude distribution |br|
Default = 520.0

**--muller-mandel-sigma-kick** |br|
Scatter width for NS and BH kicks when using the `MULLERMANDEL` kick magnitude distribution |br|
Default = 0.3

.. _options-props-N:

:ref:`Back to Top <options-props-top>`

**--natal-kick-for-PPISN**
TRUE indicates PPISN remnants will receive natal kicks via the same prescription as CCSN remnants. |br|
FALSE indicates PPISN remnants will receive no natal kicks. |br|
Default = TRUE  

**--neutrino-mass-loss-BH-formation** |br|
Assumption about neutrino mass loss during BH formation (works with `FRYER2012` or `FRYER2022` ``--remnant-mass-prescription``, but not `MANDELMUELLER`). |br|
Options: { FIXED_FRACTION, FIXED_MASS } |br|
Default = FIXED_MASS

**--neutrino-mass-loss-BH-formation-value** |br|
Amount of mass lost in neutrinos during BH formation (either as fraction or in solar masses, depending on the value of ``--neutrino-mass-loss-bh-formation``). |br|
Default = 0.1

**--neutron-star-equation-of-state** |br|
Neutron star equation of state. |br|
Options: { SSE, ARP3 } |br|
Default = SSE

**--notes** |br|
Annotation strings (vector). |br|
Default = ""

**--notes-hdrs** |br|
Annotations header strings (vector). |br|
Default = `No annotations`

**--number-of-systems [ -n ]** |br|
The number of systems to simulate. |br|
Single stars for SSE mode; binary stars for BSE mode. |br|
This option is ignored if either of the following is true: |br|

    - the user specified a grid file |br|
    - the user specified a range or set for any options - this implies a grid |br|

In both cases the number of objects evolved will be the number specified by the grid. |br|
Default = 10

.. _options-props-O:

:ref:`Back to Top <options-props-top>`

**--OB-mass-loss** |br|
Main sequence mass loss prescription. |br|
Options: { NONE, ZERO, VINK2001, VINK2021, BJORKLUND2022, KRTICKA2018 } |br|
NONE/ZERO turns off mass loss for main sequence stars.|br|
Default = VINK2021 |br|
DEPRECATION NOTICE: this option has been deprecated and will soon be removed. Please use ``--OB-mass-loss-prescription`` in future. |br|
DEPRECATION NOTICE: the value of ``NONE`` for this option has been deprecated and will soon be removed. Please use ``ZERO`` in future.

**--OB-mass-loss-prescription** |br|
Main sequence mass loss prescription. |br|
Options: { NONE, ZERO, VINK2001, VINK2021, BJORKLUND2022, KRTICKA2018 } |br|
NONE/ZERO turns off mass loss for main sequence stars.|br|
Default = VINK2021 |br|
DEPRECATION NOTICE: the value of ``NONE`` for this option has been deprecated and will soon be removed. Please use ``ZERO`` in future.

**--orbital-period** |br|
Initial orbital period for a binary star when evolving in BSE mode (days). |br|
Used only if the semi-major axis is not specified via ``--semi-major-axis``. |br|
Default: value is sampled if option not specified.

**--orbital-period-distribution** |br|
Initial orbital period distribution. |br|
Options: { FLATINLOG } |br|
Default = FLATINLOG

**--orbital-period-max** |br|
Maximum period to generate (days). |br|
Default = 1000.0

**--orbital-period-min** |br|
Minimum period to generate (days). |br|
Default = 1.1

**--output-container [ -c ]** |br|
Container (directory) name for output files. |br|
Default = ’COMPAS_Output’

**--output-path [ -o ]** |br|
Path to which output is saved (i.e. directory in which the output container is created). |br|
Default = Current working directory (CWD)

**--overall-wind-mass-loss-multiplier** |br|
Multiplicative constant for overall wind mass loss. |br|
Note that this multiplication factor is applied after the ``luminous-blue-variable-multiplier``,
the ``wolf-rayet-multiplier``, and the ``cool-wind-mass-loss-multiplier``. |br|
Default = 1.0

.. _options-props-P:

:ref:`Back to Top <options-props-top>`

**--pair-instability-supernovae** |br|
Enable pair instability supernovae (PISN). |br|
Default = TRUE

**--PISN-lower-limit** |br|
Minimum core mass for PISN (:math:`M_\odot`). |br|
Default = 60.0

**--PISN-upper-limit** |br|
Maximum core mass for PISN (:math:`M_\odot`). |br|
Default = 135.0

**--population-data-printing** |br|
Print details of population. |br|
Default = FALSE

**--PPI-CO-Core-Shift-Hendriks** |br|
Shift in CO core mass for PPI (in Msol) for the Hendriks+23 PPI prescription
Default = 0.0

**--PPI-lower-limit** |br|
Minimum core mass for PPI (:math:`M_\odot`). |br|
Default = 35.0

**--PPI-upper-limit** |br|
Maximum core mass for PPI (:math:`M_\odot`). |br|
Default = 60.0

**--print-bool-as-string** |br|
Print boolean properties as ’TRUE’ or ’FALSE’. |br|
Default = FALSE

**--pulsar-birth-magnetic-field-distribution** |br|
Pulsar birth magnetic field distribution. |br|
Options: { ZERO, FLATINLOG, UNIFORM, LOGNORMAL } |br|
Default = ZERO

**--pulsar-birth-magnetic-field-distribution-max** |br|
Maximum (:math:`log_{10}`) pulsar birth magnetic field. |br|
Default = 13.0

**--pulsar-birth-magnetic-field-distribution-min** |br|
Minimum (:math:`log_{10}`) pulsar birth magnetic field. |br|
Default = 11.0

**--pulsar-birth-spin-period-distribution** |br|
Pulsar birth spin period distribution. |br|
Options: { ZERO, UNIFORM, NORMAL } |br|
Default = ZERO

**--pulsar-birth-spin-period-distribution-max** |br|
Maximum pulsar birth spin period (ms). |br|
Default = 100.0

**--pulsar-birth-spin-period-distribution-min** |br|
Minimum pulsar birth spin period (ms). |br|
Default = 10.0

**--pulsar-magnetic-field-decay-massscale** |br|
Mass scale on which magnetic field decays during accretion (:math:`M_\odot`). |br|
Default = 0.025

**--pulsar-magnetic-field-decay-timescale** |br|
Timescale on which magnetic field decays (Myr). |br|
Default = 1000.0

**--pulsar-minimum-magnetic-field** |br|
:math:`log_{10}` of the minimum pulsar magnetic field (Gauss). |br|
Default = 8.0

**--pulsational-pair-instability** |br|
Enable mass loss due to pulsational-pair-instability (PPI). |br|
Default = TRUE

**--pulsational-pair-instability-prescription** |br|
Pulsational pair instability prescription. |br|
Options: { HENDRIKS, COMPAS, STARTRACK, MARCHANT, FARMER } |br|
``HENDRIKS`` implements the prescription from Hendriks et al. 2023 |br|
``COMPAS``, ``STARTRACK`` and ``MARCHANT`` follow Woosley 2017, Belczynski et al. 2016, and Marchant et al. 2018, all as implemented in Stevenson et al. 2019. |br|
``FARMER`` follows Farmer et al. 2019 |br|
Default = MARCHANT

.. _options-props-Q:

:ref:`Back to Top <options-props-top>`

**--quiet** |br|
Suppress printing to stdout. |br|
Default = FALSE

.. _options-props-R:

:ref:`Back to Top <options-props-top>`

**--radial-change-fraction** |br|
Approximate desired fractional change in stellar radius on phase when setting SSE and BSE timesteps (applied before ``--timestep--multiplier``). |br|
Recommended value is 0.005. |br|
A value of 0.0 means that this choice is ignored and timestep estimates will be calculated by COMPAS. |br|
Default = 0

**--random-seed** |br|
Value to use as the seed for the random number generator. |br|
Default = 0

**--remnant-mass-prescription** |br|
Remnant mass prescription. |br|
Options: { HURLEY2000, BELCZYNSKI2002, FRYER2012, FRYER2022, MULLER2016, MULLERMANDEL, SCHNEIDER2020, SCHNEIDER2020ALT } |br|
Default = MULLERMANDEL

**--retain-core-mass-during-caseA-mass-transfer** |br|
If TRUE, preserve a larger donor core mass following case A mass transfer. |br|
The core is set equal to the expected core mass of a newly formed HG star with mass equal to that of the donor,
scaled by the fraction of the donor's MS lifetime at mass transfer. |br|
Default = TRUE

**--revised-energy-formalism-nandez-ivanova** |br|
Enable revised energy formalism of Nandez & Ivanova. |br|
Default = FALSE

**--rlof-printing** |br|
Print RLOF events to logfile. |br|
Default = TRUE

**--rocket-kick-magnitude-1** |br|
Magnitude of post-SN pulsar rocket kick for the primary, in km/s. |br|
Default = 0.0 

**--rocket-kick-magnitude-2** |br|
Magnitude of post-SN pulsar rocket kick for the secondary, in km/s. |br|
Default = 0.0 

**--rocket-kick-phi-1** |br|
The in-plane angle :math:`[0.0, 2pi)` of the rocket kick velocity that primary neutron star receives following the supernova. |br|
Default = 0.0 

**--rocket-kick-phi-2** |br|
The in-plane angle :math:`[0.0, 2pi)` of the rocket kick velocity that secondary neutron star receives following the supernova. |br|
Default = 0.0 

**--rocket-kick-theta-1** |br|
The polar angle [0, pi] of the rocket kick velocity that primary neutron star receives following the supernova. 0 is aligned with orbital AM. |br|
Default = 0.0 

**--rocket-kick-theta-2** |br|
The polar angle :math:`[0, pi]`` of the rocket kick velocity that secondary neutron star receives following the supernova. 0 is aligned with orbital AM. |br|
Default = 0.0 

**--rotational-frequency** |br|
Initial rotational frequency of the star for SSE (Hz). |br|
Default = 0.0 (``--rotational-velocity-distribution`` used if ``--rotational-frequency`` not specified)

**--rotational-frequency-1** |br|
Initial rotational frequency of the primary star for BSE (Hz). |br|
Default = 0.0 (``--rotational-velocity-distribution`` used if ``--rotational-frequency-1`` not specified)

**--rotational-frequency-2** |br|
Initial rotational frequency of the secondary star for BSE (Hz). |br|
Default = 0.0 (``--rotational-velocity-distribution`` used if ``--rotational-frequency-2`` not specified)

**--rotational-velocity-distribution** |br|
Initial rotational velocity distribution. |br|
Options: { ZERO, HURLEY, VLTFLAMES } |br|
``ZERO`` sets all initial rotational velocities to 0.0. |br|
``HURLEY`` and ``VLTFLAMES`` sample initial rotational velocities from the Hurley, Pols, Tout (2000) and Ramirez-Agudelo et al. (2013,2015), respectively |br|
Default = ZERO

**--RSG-mass-loss** |br|
Red supergiant mass loss prescription. |br|
Options: { NONE, ZERO, VINKSABHAHIT2023, BEASOR2020, DECIN2023, YANG2023, KEE2021, NJ90 } |br|
NONE/ZERO turns off mass loss for giant stars (stellar types CHeB, FGB, EAGB, TPAGB) below the ``RSG_MAXIMUM_TEMP``. |br|
Default = DECIN2023 |br|
DEPRECATION NOTICE: this option has been deprecated and will soon be removed. Please use ``--RSG-mass-loss-prescription`` in future. |br|
DEPRECATION NOTICE: the value of ``NONE`` for this option has been deprecated and will soon be removed. Please use ``ZERO`` in future.

**--RSG-mass-loss-prescription** |br|
Red supergiant mass loss prescription. |br|
Options: { NONE, ZERO, VINKSABHAHIT2023, BEASOR2020, DECIN2023, YANG2023, KEE2021, NJ90 } |br|
NONE/ZERO turns off mass loss for giant stars (stellar types CHeB, FGB, EAGB, TPAGB) below the ``RSG_MAXIMUM_TEMP``. |br|
Default = DECIN2023 |br|
DEPRECATION NOTICE: the value of ``NONE`` for this option has been deprecated and will soon be removed. Please use ``ZERO`` in future.

.. _options-props-S:

:ref:`Back to Top <options-props-top>`

**--scale-CHE-mass-loss-with-surface-helium-abundance** |br|
Scale mass loss for chemically homogeneously evolving (CHE) stars with the surface helium abundance. 
Transition from OB to WR mass loss towards the end of the main sequence.
Default = False

**--scale-terminal-wind-velocity-with-metallicity-power** |br|
Scale terminal wind velocity with metallicity to this power
Default = 0.0

**--semi-major-axis** |br|
Initial semi-major axis for a binary star when evolving in BSE mode (AU). |br|
Default = 0.1

**--semi-major-axis-distribution [ -a ]** |br|
Initial semi-major axis distribution. |br|
Options: { FLATINLOG, DUQUENNOYMAYOR1991, SANA2012 } |br|
Default = FLATINLOG

**--semi-major-axis-max** |br|
Maximum semi-major axis to generate (AU). |br|
Default = 1000.0

**--semi-major-axis-min** |br|
Minimum semi-major axis to generate (AU). |br|
Default = 0.01

**--stellar-zeta-prescription** |br|
Prescription for convective donor radial response zeta. |br|
Options: { SOBERMAN, HURLEY, ARBITRARY } |br|
The prescription only applies to stars with convective envelopes. |br|
Stars with radiative envelopes take the values from ``--zeta-main-sequence`` or ``--zeta-radiative-giant-star``. |br|
``SOBERMAN`` uses zeta from Soberman, Phinney, and van den Heuvel (1997). |br|
``HURLEY`` uses zeta from Hurley, Pols, Tout (2002). |br|
``ARBITRARY`` uses fixed value set in ``--zeta-adiabatic-arbitrary``. |br|
Default = SOBERMAN

**--store-input-files** |br|
Enables copying of any specified grid file and/or logfile-definitios file to the COMPAS output container. |br|
Default = TRUE

**--switch-log** |br|
Enables printing of the Switch Log logfile. |br|
Default = FALSE

.. _options-props-T:

:ref:`Back to Top <options-props-top>`

**--tides-prescription** |br|
Prescription for tidal evolution of the binary. |br|
Options: { NONE, PERFECT, KAPIL2024 } |br|
``NONE`` disables tidal interactions. |br|
``PERFECT`` evolves the binary assuming instantaneous synchronization and circularization. |br|
``KAPIL2024`` uses the prescription from Kapil+ (2024). |br|
Default = NONE

**--timestep-filename** |br|
User-defined timesteps filename. (See :doc:`Timestep files <../timestep-files>`) |br|
Default = ’’ (None)

**--timestep-multiplier** |br|
Multiplicative factor for timestep duration. |br|
Default = 1.0

.. _options-props-U:

:ref:`Back to Top <options-props-top>`

**--use-mass-loss** |br|
Enable mass loss through winds. |br|
Default = TRUE
Note that setting this option to FALSE can have unexpected consequences, e.g., TPAGB stars that are prevented from losing mass 
cannot become white dwarfs, so will become massless remnants.  This is a useful option for testing, but this setting is not recommended
for production. It is better to use specific wind prescription controls, such as: |br|
``--cool-wind-mass-loss-multiplier`` |br|
``--LBV-mass-loss-prescription`` |br|
``--luminous-blue-variable-multiplier`` |br|
``--mass-loss-prescription`` |br|
``--overall-wind-mass-loss-multiplier`` |br|
``--wolf-rayet-multiplier`` |br|

**--use-mass-transfer** |br|
Enable mass transfer. |br|
Default = TRUE

.. _options-props-V:

:ref:`Back to Top <options-props-top>`

**--version [ -v ]** |br|
Prints COMPAS version string.

**--VMS-mass-loss** |br|
Very massive main sequence mass loss prescription. |br|
Options: { NONE, ZERO, VINK2011, SABHAHIT2023, BESTENLEHNER2020 } |br|
Applied above the VMS_MASS_THRESHOLD (:math:`100 M_\odot`) by default. |br|
NONE/ZERO turns off VMS mass loss. |br|
Default = SABHAHIT2023 |br|
DEPRECATION NOTICE: this option has been deprecated and will soon be removed. Please use ``--VMS-mass-loss-prescription`` in future. |br|
DEPRECATION NOTICE: the value of ``NONE`` for this option has been deprecated and will soon be removed. Please use ``ZERO`` in future.

**--VMS-mass-loss-prescription** |br|
Very massive main sequence mass loss prescription. |br|
Options: { NONE, ZERO, VINK2011, SABHAHIT2023, BESTENLEHNER2020 } |br|
Applied above the VMS_MASS_THRESHOLD (:math:`100 M_\odot`) by default. |br|
NONE/ZERO turns off VMS mass loss. |br|
Default = SABHAHIT2023 |br|
DEPRECATION NOTICE: the value of ``NONE`` for this option has been deprecated and will soon be removed. Please use ``ZERO`` in future.

.. _options-props-W:

:ref:`Back to Top <options-props-top>`

**--wolf-rayet-multiplier** |br|
Multiplicative constant for Wolf Rayet winds. Note that wind mass loss will also be multiplied by the
``overall-wind-mass-loss-multiplier``. |br|
Default = 1.0

**--WR-mass-loss** |br|
Wolf-Rayet mass loss prescription. |br|
Options: { BELCZYNSKI2010, SANDERVINK2023, SHENAR2019 } |br|
Default = SANDERVINK2023 |br|
DEPRECATION NOTICE: this option has been deprecated and will soon be removed. Please use ``--WR-mass-loss-prescription`` in future.

**--WR-mass-loss-prescription** |br|
Wolf-Rayet mass loss prescription. |br|
Options: { BELCZYNSKI2010, SANDERVINK2023, SHENAR2019 } |br|
Default = SANDERVINK2023

.. _options-props-X:
.. _options-props-Y:

:ref:`Back to Top <options-props-top>`

**--YAML-template** |br|
Template filename for creation of YAML file (see also ``--create-YAML-file``). |br|
Default = "" (No template file)

.. _options-props-Z:

:ref:`Back to Top <options-props-top>`

**--zeta-adiabatic-arbitrary** |br|
Value of logarithmic derivative of radius with respect to mass, :math:`\zeta` adiabatic. |br|
Default = :math:`1.0 \times 10^4`

**--zeta-main-sequence** |br|
Value of logarithmic derivative of radius with respect to mass, :math:`\zeta` on the main sequence. |br|
Default = 2.0

**--zeta-radiative-giant-star** |br|
Value of logarithmic derivative of radius with respect to mass, :math:`\zeta` for radiative-envelope giant-like stars
(including Hertzsprung Gap (HG) stars). |br|
Default = 6.5


Category listing
----------------

Go to :ref:`the top of this page <options-props-top>` for the full alphabetical list of options with explanations and default values

.. _options-initial-conditions:

**Initial conditions**

--initial-mass-function, --initial-mass, --initial-mass-1, --initial-mass-2, --initial-mass-min, --initial-mass-max, --initial-mass-power

--mass-ratio-distribution, --mass-ratio, --mass-ratio-min, --mass-ratio-max, --minimum-secondary-mass

--eccentricity-distribution, --eccentricity, --eccentricity-min, --eccentricity-max

--metallicity-distribution, --metallicity, --metallicity-min, --metallicity-max

--orbital-period-distribution, --orbital-period, --orbital-period-min, --orbital-period-max, --semi-major-axis-distribution, --semi-major-axis, 
--semi-major-axis-min, --semi-major-axis-max, --allow-rlof-at-birth, --allow-touching-at-birth

--rotational-velocity-distribution, --rotational-frequency, --rotational-frequency-1, --rotational-frequency-2

:ref:`Back to Top <options-props-top>`

.. _options-stellar-evolution:

**Stellar evolution and winds**

--use-mass-loss, --check-photon-tiring-limit, --cool-wind-mass-loss-multiplier, --luminous-blue-variable-prescription, --LBV-mass-loss-prescription
--luminous-blue-variable-multiplier, --mass-loss-prescription, --overall-wind-mass-loss-multiplier, --wolf-rayet-multiplier, 
--expel-convective-envelope-above-luminosity-threshold, --luminosity-to-mass-threshold,
--OB-mass-loss, --OB-mass-loss-prescription, --RSG-mass-loss, --RSG-mass-loss-prescription, --VMS-mass-loss, --vms-mass-loss-prescription, --WR-mass-loss, --WR-mass-loss-prescription

--chemically-homogeneous-evolution, --chemically-homogeneous-evolution-mode

:ref:`Back to Top <options-props-top>`

.. _options-mass-transfer:

**Mass transfer physics**

--case-bb-stability-prescription, --convective-envelope-temperature-threshold, --critical-mass-ratio-prescription,
--critical-mass-ratio-HG-degenerate-accretor, --critical-mass-ratio-HG-non-degenerate-accretor, --critical-mass-ratio-MS-high-mass-degenerate-accretor,
--critical-mass-ratio-MS-high-mass-non-degenerate-accretor, --critical-mass-ratio-MS-low-mass-degenerate-accretor, --critical-mass-ratio-MS-low-mass-non-degenerate-accretor,
--critical-mass-ratio-giant-degenerate-accretor, --critical-mass-ratio-giant-non-degenerate-accretor, --critical-mass-ratio-helium-HG-degenerate-accretor,
--critical-mass-ratio-helium-HG-non-degenerate-accretor, --critical-mass-ratio-helium-MS-degenerate-accretor, --critical-mass-ratio-helium-MS-non-degenerate-accretor, 
--critical-mass-ratio-helium-giant-degenerate-accretor, --critical-mass-ratio-helium-giant-non-degenerate-accretor, --critical-mass-ratio-white-dwarf-degenerate-accretor, 
--critical-mass-ratio-white-dwarf-non-degenerate-accretor, --eddington-accretion-factor, --mass-transfer, --use-mass-transfer, --mass-transfer-accretion-efficiency-prescription, 
--mass-transfer-angular-momentum-loss-prescription, --mass-transfer-fa, --mass-transfer-jloss, --mass-transfer-jloss-macleod-linear-fraction-degen, --mass-transfer-jloss-macleod-linear-fraction-non-degen, 
--mass-transfer-rejuvenation-prescription, --mass-transfer-thermal-limit-accretor, --mass-transfer-thermal-limit-accretor-multiplier, --mass-transfer-thermal-limit-C, --retain-core-mass-during-caseA-mass-transfer, 
--stellar-zeta-prescription, --zeta-adiabatic-arbitrary, --zeta-main-sequence, --zeta-radiative-giant-star 

--circulariseBinaryDuringMassTransfer, --angular-momentum-conservation-during-circularisation, --tides-prescription

--envelope-state-prescription, --common-envelope-alpha, --common-envelope-alpha-thermal, --common-envelope-formalism,
--common-envelope-lambda-prescription, --common-envelope-lambda, 
--common-envelope-slope-kruckow, --common-envelope-lambda-multiplier, --common-envelope-lambda-nanjing-enhanced, 
--common-envelope-lambda-nanjing-interpolate-in-mass, --common-envelope-lambda-nanjing-interpolate-in-metallicity, 
--common-envelope-lambda-nanjing-use_rejuvenated-mass, --common-envelope-allow-main-sequence-survive, --common-envelope-allow-radiative-envelope-survive, 
--common-envelope-allow-immediate-RLOF-post-CE-survive, --common-envelope-mass-accretion-prescription, --common-envelope-mass-accretion-constant, 
--common-envelope-mass-accretion-min, --common-envelope-mass-accretion-max, --common-envelope-recombination-energy-density, --maximum-mass-donor-nandez-ivanova, 
--revised-energy-formalism-nandez-ivanova

:ref:`Back to Top <options-props-top>`

.. _options-supernovae:

**Supernovae**

--remnant-mass-prescription, --fryer-supernova-engine, --fryer-22-fmix, --fryer-22-mcrit, --maximum-neutron-star-mass, --mcbur1, --allow-non-stripped-ECSN, 
--neutrino-mass-loss-BH-formation, --neutrino-mass-loss-BH-formation-value, --neutron-star-equation-of-state, --pair-instability-supernovae, --PISN-lower-limit, 
--PISN-upper-limit, --PPI-CO-Core-Shift-Hendriks, --PPI-lower-limit, --PPI-upper-limit, --pulsational-pair-instability, --pulsational-pair-instability-prescription

--pulsar-birth-magnetic-field-distribution, --pulsar-birth-magnetic-field-distribution-min, --pulsar-birth-magnetic-field-distribution-max, 
--pulsar-birth-spin-period-distribution, --pulsar-birth-spin-period-distribution-min, --pulsar-birth-spin-period-distribution-max, 
--pulsar-magnetic-field-decay-massscale, --pulsar-magnetic-field-decay-timescale, --pulsar-minimum-magnetic-field

--kick-magnitude-distribution, --kick-magnitude-sigma-CCSN-BH, --kick-magnitude-sigma-CCSN-NS, --kick-magnitude-sigma-ECSN, --kick-magnitude-sigma-USSN, 
--black-hole-kicks, --black-hole-kicks-mode, --fix-dimensionless-kick-magnitude, --kick-magnitude, --kick-magnitude-1, --kick-magnitude-2, --kick-magnitude-min, --kick-magnitude-max, 
--kick-magnitude-random, --kick-magnitude-random-1, --kick-magnitude-random-2, --kick-scaling-factor, -muller-mandel-kick-multiplier-BH, 
--muller-mandel-kick-multiplier-NS, --muller-mandel-sigma-kick

--kick-direction, --kick-direction-distribution, --kick-direction-power, --kick-mean-anomaly-1, --kick-mean-anomaly-2, --kick-phi-1, --kick-phi-2, --kick-theta-1, --kick-theta-2

:ref:`Back to Top <options-props-top>`

.. _options-admin:

**Administrative**

--mode, --number-of-systems, 
--emit-gravitational-radiation, --evolve-double-white-dwarfs, --evolve-main-sequence-mergers, --evolve-pulsars, --evolve-unbound-systems, 
--include-WD-binaries-as-DCO,
--mass-change-fraction, --maximum-evolution-time, --maximum-number-timestep-iterations,
--radial-change-fraction, --random-seed, --timestep-multiplier, --timestep-filename

--fp-error-mode

--grid, --grid-start-line, --grid-lines-to-process

--add-options-to-sysparms, --debug-classes, --debug-level, --debug-to-file, --detailed-output, --detailed-output, --enable-warnings, --errors-to-file, 
--help, --notes, --notes-hdrs, --population-data-printing, --print-bool-as-string, --quiet, --version

--log-classes, --logfile-definitions, --logfile-name-prefix, --logfile-type, --log-level, --logfile-common-envelopes, --logfile-common-envelopes-record-types, 
--logfile-detailed-output, --logfile-detailed-output-record-types, --logfile-double-compact-objects, --logfile-double-compact-objects-record-types, 
--logfile-pulsar-evolution, --logfile-pulsar-evolution-record-type, --logfile-rlof-parameters, --logfile-rlof-parameters-record-types, --logfile-supernovae, 
--logfile-supernovae-record-types, --logfile-switch-log, --logfile-system-parameters, --logfile-system-parameters-record-types, --output-container, 
--output-path, --rlof-printing, --store-input-files, --switch-log, --hdf5-buffer-size, --hdf5-chunk-size

--create-YAML-file, YAML-template

:ref:`Back to Top <options-props-top>`


