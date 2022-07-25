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
Black hole kicks relative to NS kicks. |br|
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
Default = PESSIMISTIC

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
Allow binaries with an evolved component with a radiative envelope to survive the common envelope phase. |br|
Deafult = FALSE

**--common-envelope-alpha** |br|
Common Envelope efficiency alpha. |br|
Default = 1.0

**--common-envelope-alpha-thermal** |br|
Thermal energy contribution to the total envelope binding energy. |br|
Defined such that :math:`\lambda = \alpha_{th} \times \lambda_{b} + (1.0 - \alpha_{th}) \times \lambda_{g}`. |br|
Default = 1.0

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
``ZERO`` is no accretion; ``CONSTANT`` means a fixed amount of accretion determined by ``--common-envelope-mass-accretion-constant``; ``UNIFORM`` means a uniform random draw between ``--common-envelope-mass-accretion-min`` and ``--common-envelope-mass-accretion-max`` (Oslowski et al., 2011); and ``MACLEOD`` follows the prescription of MacLeod et al., 2015 |br|
Options: { ZERO, CONSTANT, UNIFORM, MACLEOD } |br|
Default = ZERO

**--common-envelope-recombination-energy-density** |br|
Recombination energy density (erg g−1). |br|
Default = :math:`1.5 \times 10^{13}`

**--common-envelope-slope-kruckow** |br|
Slope for the Kruckow lambda (see Kruckow et al. 2016 as implemented by Vigna-Gomez et al. 2018). |br|
Default = −0.833333

**--cool-wind-mass-loss-multiplier** |br|
Multiplicative constant for wind mass loss of cool stars, i.e. those with temperatures below the
VINK_MASS_LOSS_MINIMUM_TEMP (default 12500K). |br|
Only applicable when ``--mass-loss-prescription = VINK``. |br|
Default = 1.0

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

**--enable-warnings** |br|
Display warning messages to stdout. |br|
Default = FALSE

**--envelope-state-prescription** |br|
Prescription for determining whether the envelope of the star is convective or radiative. |br|
Options: { LEGACY, HURLEY, FIXED_TEMPERATURE } |br|
``LEGACY`` refers to the model used in Stevenson et al., 2017; ``HURLEY`` refers to the model of Hurley, Pols, Tout, 2002; and ``FIXED_TEMPERATURE`` assumes that a deep convective envelope developes only when the temperature drops below ``CONVECTIVE_BOUNDARY_TEMPERATURE`` (Klencki et al., 2020) |br|
Default = LEGACY


**--errors-to-file** |br|
Write error messages to file. |br|
Default = FALSE

**--evolve-pulsars** |br|
Evolve pulsar properties of Neutron Stars. |br|
Default = FALSE

**--evolve-unbound-systems** |br|
Continue evolving stars even if the binary is disrupted. |br|
Default = FALSE

.. _options-props-F:

:ref:`Back to Top <options-props-top>`

**--fix-dimensionless-kick-magnitude** |br|
Fix dimensionless kick magnitude to this value. |br|
Default = n/a (not used if option not present)

**--fryer-supernova-engine** |br|
Supernova engine type if using the remnant mass prescription from :cite:`Fryer2012`. |br|
Options: { DELAYED, RAPID }
Default = DELAYED

**--fryer-22-fmix** |br|
Paramter describing the mixing growth time when using the 'FRYER2022' remnant mass distribution  :cite:`Fryer2022`. |br|
Default = 0.5, which is closest to the 'DELAYED' remnant mass prescription from :cite:`Fryer2012`. A value of 4.0 is closest to  the 'RAPID' remnant mass prescription from :cite:`Fryer2012`. |br|
If the FALLBACK option is used for the kicks, then the proto core masses will be determined by the fryer-supernova-engine option. |br| 


**--fryer-22-mcrit** |br|
Critical CO core mass for black hole formation when using the 'FRYER2022' remnant mass distribution :cite:`Fryer2022`. |br|
Default = 5.75



.. _options-props-G:

:ref:`Back to Top <options-props-top>`

**--grid** |br|
Grid filename. |br|
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
``SALPETER`` and ``KROUPA`` use the IMFs of Salpeter 1955 and Kroupa 2001, ``POWERLAW`` samples from a single power law with slope ``--initial-mass-power``, and ``UNIFORM`` samples uniformly between ``--initial-mass-min`` and ``--initial-mass-min`` |br|
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
Kick angles are defined relative to the spin axis.  ``INPLANE`` and ``PERPENDICULAR`` are strictly in the equatorial plane or in polar directions, while ``WEDGE`` and ``POLES`` are preferentially but exactly in the equatorial plane or in polar directions with 1 degree scales, respectively;  ``POWERLAW`` quantifies the preference for polar vs planar kicks with the ``--kick-direction-power`` parameter. |br|
Default = ISOTROPIC

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
``ZERO`` assigns kick magnitudes of 0, ``FIXED`` always sets the magnitude to a fixed value based on supernova type, ``FLAT`` and ``MAXWELLIAN`` draw kicks from uniform or Maxwellian (e.g., Hobbs et al., 2005) distributions, respectively, ``BRAYELDRIDGE`` and ``MULLERMANDEL`` use momenum-preserving kicks from Bray & Eldrigde 2018 and Mandel & Mueller 2020, respectively, and ``MULLER2016`` and ``MULLER2016MAXWELLIAN`` use kicks from Mueller 2016 as implemented in Vigna-Gomez et al., 2018 (reduced by a factor of sqrt(3) in the latter case). |br|
Default = MAXWELLIAN

**--kick-magnitude-max** |br|
Maximum drawn kick magnitude (:math:`km s^{−1}`). |br|
Must be > 0 if using ``--kick-magnitude-distribution = FLAT``. |br|
Default = −1.0

**--kick-magnitude-random** |br|
CDF value to be used to draw the kick magnitude for a single star when evolving in SSE mode, should the star undergo a supernova event and should the chosen distribution sample from a cumulative distribution function. |br|
Must be a floating-point number in the range :math:`[0.0, 1.0)`. |br|
The specified value for this option will be used in preference to any specified value for ``--kick-magnitude``. |br|
Default = Random number drawn uniformly from :math:`[0.0, 1.0)`

**--kick-magnitude-random-1** |br|
CDF value to be used to draw the kick magnitude for the primary star of a binary system when evolving in BSE mode, should the star undergo a supernova event and should the chosen distribution sample from a cumulative distribution function. |br|
Must be a floating-point number in the range :math:`[0.0, 1.0)`. |br|
The specified value for this option will be used in preference to any specified value for ``--kick-magnitude-1``. |br|
Default = Random number drawn uniformly from :math:`[0.0, 1.0)`

**--kick-magnitude-random-2** |br|
CDF value to be used to draw the kick magnitude for the secondary star of a binary system when evolving in BSE mode, should the star undergo a supernova event and should the chosen distribution sample from a cumulative distribution function. |br|
Must be a floating-point number in the range :math:`[0.0, 1.0)`. |br|
The specified value for this option will be used in preference to any specified value for ``--kick-magnitude-2``. |br|
Default = Random number drawn uniformly from :math:`[0.0, 1.0)`

**--kick-magnitude-sigma-CCSN-BH** |br|
Sigma for chosen kick magnitude distribution for black holes (:math:`km s^{−1}`); ignored if not needed for the chosen kick magnitude distribution. |br|
Default = 265.0

**--kick-magnitude-sigma-CCSN-NS** |br|
Sigma for chosen kick magnitude distribution for neutron stars (:math:`km s^{−1}`); ignored if not needed for the chosen kick magnitude distribution. |br|
Default = 265.0

**--kick-magnitude-sigma-ECSN** |br|
Sigma for chosen kick magnitude distribution for ECSN (:math:`km s^{−1}`); ignored if not needed for the chosen kick magnitude distribution. |br|
Default = 30.0

**--kick-magnitude-sigma-USSN** |br|
Sigma for chosen kick magnitude distribution for USSN (:math:`km s^{−1}`); ignored if not needed for the chosen kick magnitude distribution. |br|
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
The angle between ’x’ and ’y’, both in the orbital plane of the supernova vector, for the primary star of a binary system when evolving in BSE mode, should it undergo a supernova event (radians). |br|
Default = Drawn according to specified ``--kick-direction`` distribution

**--kick-phi-2** |br|
The angle between ’x’ and ’y’, both in the orbital plane of the supernova vector, for the secondary
star of a binary system when evolving in BSE mode, should it undergo a supernova event (radians). |br|
Default = Drawn according to specified ``--kick-direction`` distribution

**--kick-scaling-factor** |br|
Arbitrary factor used to scale kicks. |br|
Default = 1.0

**--kick-theta-1** |br|
The angle between the orbital plane and the ’z’ axis of the supernova vector for the primary star of a binary system when evolving in BSE mode, should it undergo a supernova event (radians). |br|
Default = Drawn according to specified ``--kick-direction`` distribution

**--kick-theta-2** |br|
The angle between the orbital plane and the ’z’ axis of the supernova vector for the secondary star of a binary system when evolving in BSE mode, should it undergo a supernova event (radians). |br|
Default = Drawn according to specified ``--kick-direction`` distribution

.. _options-props-L:

:ref:`Back to Top <options-props-top>`

**--log-classes** |br|
Logging classes to be enabled (vector). |br|
Default = `All debug classes enabled (e.g. no filtering)`

**--logfile-common-envelopes** |br|
Filename for Common Envelopes logfile (BSE mode). |br|
Default = ’BSE_Common_Envelopes’

**--logfile-definitions** |br|
Filename for logfile record definitions file. |br|
Default = ’’ (None)

**--logfile-detailed-output** |br|
Filename for the Detailed_Output logfile. |br|
Default = ’SSE_Detailed_Output’ for SSE mode; ’BSE_Detailed_Output’ for BSE mode |br|

**--logfile-double-compact-objects** |br|
Filename for the Double Compact Objects logfile (BSE mode). |br|
Default = ’BSE_Double_Compact_Objects’

**--logfile-name-prefix** |br|
Prefix for logfile names. |br|
Default = ’’ (None)

**--logfile-pulsar-evolution** |br|
Filename for the Pulsar Evolution logfile (BSE mode). |br|
Default = ’BSE_Pulsar_Evolution’

**--logfile-rlof-parameters** |br|
Filename for the RLOF Printing logfile (BSE mode). |br|
Default = ’BSE_RLOF’

**--logfile-supernovae** |br|
Filename for the Supernovae logfile. |br|
Default = ’SSE_Supernovae’ for SSE mode; ’BSE_Supernovae’ for BSE mode |br|

**--logfile-switch-log** |br|
Filename for the Switch Log logfile. |br|
Default = ’SSE_Switch_Log’ for SSE mode; ’BSE_Switch_Log’ for BSE mode |br|

**--logfile-system-parameters** |br|
Filename for the System Parameters logfile (BSE mode). |br|
Default = ’SSE_System_Parameters’ for SSE mode; ’BSE_System_Parameters’ for BSE mode |br|

**--logfile-type** |br|
The type of logfile to be produced by COMPAS. |br|
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
Options: { NONE, HURLEY, HURLEY_ADD, BELCZYNSKI } |br|
No LBV winds for ``NONE``,  Hurley, Pols, Tout (2000) LBV winds only for ``HURLEY`` LBV stars (or in addition to other winds for ``HURLEY_ADD``, Belzcynski et al. 2010 winds for ``BELCZYNSKI`` |br|
Default = HURLEY_ADD

.. _options-props-M:

:ref:`Back to Top <options-props-top>`

**--mass-loss-prescription** |br|
Mass loss prescription. |br|
Options: { NONE, HURLEY, VINK } |br|
Default = VINK

**--mass-ratio [ -q ]** |br|
Mass ratio :math:`\frac{m2}{m1}` used to determine secondary mass if not specified via ``--initial-mass-2``. |br|
Default = Value is sampled if option not specified.

**--mass-ratio-distribution** |br|
Initial mass ratio distribution for :math:`q = \frac{m2}{m1}`. |br|
Options: { FLAT, DUQUENNOYMAYOR1991, SANA2012 } |br|
``FLAT`` is uniform in the mass ratio between ``--mass-ratio-min`` and ``--mass-ratio-max``, the other prescriptions follow Duquennoy & Mayor 1991 and Sana et al. 2012 |br|
Default = FLAT

**--mass-ratio-max** |br|
Maximum mass ratio :math:`\frac{m2}{m1}` to generate. |br|
Default = 1.0

**--mass-ratio-min** |br|
Minimum mass ratio :math:`\frac{m2}{m1}` to generate. |br|
Default = 0.01

**--mass-transfer** |br|
Enable mass transfer. |br|
Default = TRUE

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

**--mass-transfer-jloss-macleod-linear-fraction** |br|
Specific angular momentum interpolation fraction, linear between 0 and 1 corresponding to the accretor and L2 point. |br|
Used when ``--mass-transfer-angular-momentum-loss-prescription = MACLEOD_LINEAR``, ignored otherwise. |br|
Default = 0.5

**--mass-transfer-rejuvenation-prescription** |br|
Mass Transfer Rejuvenation prescription. |br|
Options: { NONE, STARTRACK } |br|
``NONE`` uses the Hurley, Pols, Tout (2000) model, ``STARTRACK`` uses the model from Belczynski et al. 2008 |br|
Default = STARTRACK

**--mass-transfer-thermal-limit-accretor** |br|
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
Minimum core mass at base of AGB to avoid fully degnerate CO core formation (:math:`M_\odot`). |br|
e.g. 1.6 in :cite:`Hurley2000` presciption; 1.83 in :cite:`Fryer2012` and :doc:`Belczynski et al. (2008) <../../references>` models. |br|
Default = 1.6

**--metallicity [ -z ]** |br|
Metallicity. |br|
The value specified for metallicity is applied to both stars for BSE mode. |br|
Default = 0.0142

**--metallicity-distribution** |br|
Metallicity distribution. |br|
Options: { ZSOLAR, LOGUNIFORM } |br|
``ZSOLAR`` uses ``ZSOL_ASPLUND`` for all initial metallicities, ``LOGUNIFORM`` draws the metallicity uniformly in the log between ``--metallicity-min`` and ``--metallicity-max`` |br|
Default = ZSOLAR

**--metallicity-max** |br|
Maximum metallicity to generate. |br|
Default = 0.03

**--metallicity-min** |br|
Minimum metallicity to generate. |br|
Default = 0.0001

**--minimum-secondary-mass** |br|
Minimum mass of secondary to generate (:math:`M_\odot`). |br|
Default = 0.1 if ``--initial-mass-2`` specified; value of ``--initial-mass-min`` if ``--initial-mass-2`` not specified.

**--mode** |br|
The mode of evolution. |br|
Options: { SSE, BSE } |br|
Default = BSE

**--muller-mandel-kick-multiplier-BH** |br|
Scaling prefactor for BH kicks when using the `MULLERMANDEL` kick magnitude distribution |br|
Default = 200.0

**--muller-mandel-kick-multiplier-NS** |br|
Scaling prefactor for NS kicks when using the `MULLERMANDEL` kick magnitude distribution |br|
Default = 400.0

**--muller-mandel-sigma-kick** |br|
Scatter width for NS and BH kicks when using the `MULLERMANDEL` kick magnitude distribution |br|
Default = 0.3

.. _options-props-N:

:ref:`Back to Top <options-props-top>`

**--neutrino-mass-loss-BH-formation** |br|
Assumption about neutrino mass loss during BH formation. |br|
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

**--orbital-period** |br|
Initial orbital period for a binary star when evolving in BSE mode (days). |br|
Used only if the semi-major axis is not specified via ``--semi-major-axis``. |br|
Default = Value is sampled if option not specified.

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
the ``wolf-rayet-multiplier`` and the ``cool-wind-mass-loss-multiplier``. |br|
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
Options: { ZERO, FIXED, FLATINLOG, UNIFORM, LOGNORMAL } |br|
Default = ZERO

**--pulsar-birth-magnetic-field-distribution-max** |br|
Maximum (:math:`log_{10}`) pulsar birth magnetic field. |br|
Default = 13.0

**--pulsar-birth-magnetic-field-distribution-min** |br|
Minimum (:math:`log_{10}`) pulsar birth magnetic field. |br|
Default = 11.0

**--pulsar-birth-spin-period-distribution** |br|
Pulsar birth spin period distribution. |br|
Options: { ZERO, FIXED, UNIFORM, NORMAL } |br|
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
Options: { COMPAS, STARTRACK, MARCHANT, FARMER } |br|
``COMPAS``, ``STARTRACK`` and ``MARCHANT`` follow Woosley 2017, Belczynski et al. 2016, and Marchant et al. 2018, all as implemented in Stevenson et al. 2019, ``FARMER`` follows Farmer et al. 2019 |br|
Default = MARCHANT

.. _options-props-Q:

:ref:`Back to Top <options-props-top>`

**--quiet** |br|
Suppress printing to stdout. |br|
Default = FALSE

.. _options-props-R:

:ref:`Back to Top <options-props-top>`

**--random-seed** |br|
Value to use as the seed for the random number generator. |br|
Default = 0

**--remnant-mass-prescription** |br|
Remnant mass prescription. |br|
Options: { HURLEY2000, BELCZYNSKI2002, FRYER2012, FRYER2022, MULLER2016, MULLERMANDEL, SCHNEIDER2020, SCHNEIDER2020ALT } |br|
Remnant mass recipes from Hurley, Pols, Tout (2000) for ``HURLEY2000``, Belczynski et al. 2002, Fryer et al. 2012,  Fryer et al. 2022, Mueller 2016, Mandel & Mueller 2020, and Schneider et al. 2020 (with the alternative prescription for effectively single stars from the same paper in the ``SCHNEIDER2020ALT`` case) |br|
Default = FRYER2012

**--retain-core-mass-during-caseA-mass-transfer |br|
If set to true, preserve a larger donor core mass following case A mass transfer.  The core is set equal to the expected core mass of a newly formed HG star with mass equal to that of the donor, scaled by the fraction of the donor's MS lifetime at mass transfer. |br|
Default = FALSE

**--revised-energy-formalism-nandez-ivanova** |br|
Enable revised energy formalism of Nandez & Ivanova. |br|
Default = FALSE

**--rlof-printing** |br|
Print RLOF events to logfile. |br|
Default = TRUE

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
``ZERO`` sets all initial rotational velocities to 0, while ``HURLEY`` and ``VLTFLAMES`` sample them from the Hurley, Pols, Tout (2000) and Ramirez-Agudelo et al. (2013,2015), respectively |br|
Default = ZERO

.. _options-props-S:

:ref:`Back to Top <options-props-top>`

**--semi-major-axis** |br|
Initial semi-major axis for a binary star when evolving in BSE mode (AU). |br|
Default = 0.1

**--semi-major-axis-distribution [ -a ]** |br|
Initial semi-major axis distribution. |br|
Options: { FLATINLOG, DUQUENNOYMAYOR1991, SANA2012 } |br|
Flat-in-log (Opik 1924), Duquennoy & Mayor (1991) or Sana et al. (2012) distributions |BR|
Default = FLATINLOG

**--semi-major-axis-max** |br|
Maximum semi-major axis to generate (AU). |br|
Default = 1000.0

**--semi-major-axis-min** |br|
Minimum semi-major axis to generate (AU). |br|
Default = 0.01

**--stellar-zeta-prescription** |br|
Prescription for stellar zeta. |br|
Options: { SOBERMAN, HURLEY, ARBITRARY } |br|
Use Soberman, Phinney, and van den Heuvel (1997) or Hurley, Pols, Tout (2002) or the fixed value specified via ``--zeta-adiabatic-arbitrary`` for the stellar radial response to mass loss for convective-envelope giant-like stars |br|
Default = SOBERMAN

**--store-input-files** |br|
Enables copying of any specified grid file and/or logfile-definitios file to the COMPAS output container. |br|
Default = TRUE

**--switch-log** |br|
Enables printing of the Switch Log logfile. |br|
Default = FALSE

.. _options-props-T:

:ref:`Back to Top <options-props-top>`

**--timestep-multiplier** |br|
Multiplicative factor for timestep duration. |br|
Default = 1.0

.. _options-props-U:

:ref:`Back to Top <options-props-top>`

**--use-mass-loss** |br|
Enable mass loss. |br|
Default = TRUE

.. _options-props-V:

:ref:`Back to Top <options-props-top>`

**--version [ -v ]** |br|
Prints COMPAS version string.

.. _options-props-W:

:ref:`Back to Top <options-props-top>`

**--wolf-rayet-multiplier** |br|
Multiplicative constant for Wolf Rayet winds. Note that wind mass loss will also be multiplied by the
``overall-wind-mass-loss-multiplier``. |br|
Default = 1.0

.. _options-props-X:
.. _options-props-Y:
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
--------------------

Go to :ref:`the top of this page <options-props-top>` for the full alphabetical list of options with explanations and default values

.. _options-initial-conditions:

**Initial conditions**

--initial-mass-function, --initial-mass, --initial-mass-1, --initial-mass-2, --initial-mass-min, --initial-mass-max, --initial-mass-power

--mass-ratio-distribution, --mass-ratio, --mass-ratio-min, --mass-ratio-max, --minimum-secondary-mass

--eccentricity-distribution, --eccentricity, --eccentricity-min, --eccentricity-max

--metallicity-distribution, --metallicity, --metallicity-min, --metallicity-max

--orbital-period-distribution, --orbital-period, --orbital-period-min, --orbital-period-max, --semi-major-axis-distribution, --semi-major-axis, --semi-major-axis-min, --semi-major-axis-max, --allow-rlof-at-birth, --allow-touching-at-birth

--rotational-velocity-distribution, --rotational-frequency, --rotational-frequency-1, --rotational-frequency-2

:ref:`Back to Top <options-props-top>`

.. _options-stellar-evolution:

**Stellar evolution and winds**

--use-mass-loss, --check-photon-tiring-limit, --cool-wind-mass-loss-multiplier, --luminous-blue-variable-prescription, --luminous-blue-variable-multiplier, --mass-loss-prescription, --overall-wind-mass-loss-multiplier, --wolf-rayet-multiplier

--chemically-homogeneous-evolution

:ref:`Back to Top <options-props-top>`

.. _options-mass-transfer:

**Mass transfer physics**

--mass-transfer, --mass-transfer-accretion-efficiency-prescription, --mass-transfer-angular-momentum-loss-prescription, --mass-transfer-fa, --mass-transfer-jloss, --mass-transfer-jloss-macleod-linear-fraction, --mass-transfer-rejuvenation-prescription, --mass-transfer-thermal-limit-accretor, --mass-transfer-thermal-limit-C, --stellar-zeta-prescription, --zeta-adiabatic-arbitrary, --zeta-main-sequence, --zeta-radiative-giant-star, --case-bb-stability-prescription, --eddington-accretion-factor, --retain-core-mass-during-caseA-mass-transfer

--circulariseBinaryDuringMassTransfer, --angular-momentum-conservation-during-circularisation

--envelope-state-prescription, --common-envelope-alpha, --common-envelope-alpha-thermal, --common-envelope-lambda-prescription, --common-envelope-lambda, --common-envelope-slope-kruckow, --common-envelope-lambda-multiplier, --common-envelope-lambda-nanjing-enhanced, --common-envelope-lambda-nanjing-interpolate-in-mass, --common-envelope-lambda-nanjing-interpolate-in-metallicity, --common-envelope-lambda-nanjing-use_rejuvenated-mass, --common-envelope-allow-main-sequence-survive, --common-envelope-allow-radiative-envelope-survive*, --common-envelope-allow-immediate-RLOF-post-CE-survive, --common-envelope-mass-accretion-prescription, --common-envelope-mass-accretion-constant, --common-envelope-mass-accretion-min, --common-envelope-mass-accretion-max, --common-envelope-recombination-energy-density, --maximum-mass-donor-nandez-ivanova, --revised-energy-formalism-nandez-ivanova

:ref:`Back to Top <options-props-top>`

.. _options-supernovae:

**Supernovae**

--remnant-mass-prescription, --fryer-supernova-engine, --fryer-22-fmix, --fryer-22-mcrit, --maximum-neutron-star-mass, --mcbur1, --allow-non-stripped-ECSN, --neutrino-mass-loss-BH-formation, --neutrino-mass-loss-BH-formation-value, --neutron-star-equation-of-state

--pair-instability-supernovae, --PISN-lower-limit, --PISN-upper-limit, --PPI-lower-limit, --PPI-upper-limit, --pulsational-pair-instability, --pulsational-pair-instability-prescription

--pulsar-birth-magnetic-field-distribution, --pulsar-birth-magnetic-field-distribution-min, --pulsar-birth-magnetic-field-distribution-max, --pulsar-birth-spin-period-distribution, --pulsar-birth-spin-period-distribution-min, --pulsar-birth-spin-period-distribution-max, --pulsar-magnetic-field-decay-massscale, --pulsar-magnetic-field-decay-timescale, --pulsar-minimum-magnetic-field

--kick-magnitude-distribution, --kick-magnitude-sigma-CCSN-BH, --kick-magnitude-sigma-CCSN-NS, --kick-magnitude-sigma-ECSN, --kick-magnitude-sigma-USSN, --black-hole-kicks, --fix-dimensionless-kick-magnitude, --kick-magnitude, --kick-magnitude-1, --kick-magnitude-2, --kick-magnitude-min, --kick-magnitude-max, --kick-magnitude-random, --kick-magnitude-random-1, --kick-magnitude-random-2, --kick-scaling-factor, -muller-mandel-kick-multiplier-BH, --muller-mandel-kick-multiplier-NS, --muller-mandel-sigma-kick

--kick-direction, --kick-direction-power, --kick-mean-anomaly-1, --kick-mean-anomaly-2, --kick-phi-1, --kick-phi-2, --kick-theta-1, --kick-theta-2

:ref:`Back to Top <options-props-top>`

.. _options-admin:

**Administrative**

--mode, --number-of-systems, --evolve-pulsars, --evolve-unbound-systems, --maximum-evolution-time, --maximum-number-timestep-iterations, --random-seed, --timestep-multiplier

--grid, --grid-start-line, --grid-lines-to-process

--add-options-to-sysparms, --debug-classes, --debug-level, --debug-to-file, --detailed-output, --enable-warnings, --errors-to-file, --help, --notes, --notes-hdrs, --population-data-printing, --print-bool-as-string, --quiet, --version

--log-classes, --logfile-definitions, --logfile-name-prefix, --logfile-type, --log-level, --logfile-common-envelopes, --logfile-detailed-output, --logfile-double-compact-objects, --logfile-pulsar-evolution, --logfile-rlof-parameters, --logfile-supernovae, --logfile-switch-log, --logfile-system-parameters, --output-container, --output-path, --rlof-printing, --store-input-files, --switch-log, --hdf5-buffer-size, --hdf5-chunk-size

:ref:`Back to Top <options-props-top>`


