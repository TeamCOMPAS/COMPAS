What's new
==========

Following is a brief list of important updates to the COMPAS code.  A complete record of changes can be found in the file ``changelog.h``.


**LATEST RELEASE** |br|

**03.01.00 Aug 24, 2024**

* New option to emit gravitational radiation at each timestep of binary evolution: ``--emit-gravitational-radiation``. The effects of radiation are approximated by the change in semimajor axis and eccentricity from Peters 1946 equations 5.6 and 5.7.  Reduce timestep if required to keep orbital separation change per step due to GW radiation within ~ 1%.

**03.00.00 Jul 26, 2024**

This is a major release of COMPAS. There are some significant changes in COMPAS operation and functionality in this release. The major change, and the impetus for
the release, is the implementation of a new, coherent, robust, error handling strategy.

Early versions of COMPAS (versions prior to v03.00.00) did not have a coherent, robust, error-handling strategy. In those versions, errors were typically displayed
as either errors or warnings (depending on the severity) as they occurred, and evolution of the star (SSE mode) or binary (BSE mode) continued - users were expected
to check errors or warnings displayed and use results with appropriate caution.  This was not ideal.

In addition, earlier versions of COMPAS did not handle floating-point errors (divide-by-zero, invalid-parameters (e.g. :math:`\sqrt{-2.0}`), etc.) well - in fact the
code effectively ignored these errors (for a detailed explanation of why this was so, refer to the Error Handling documentation in the Developer Guide
(See :ref:`developer-guide-fp-errors`).

In COMPAS version 03.00.00 the error handling philosophy has changed, and more coherent and robust error-handling code implemented. The new error-handling philosophy
is to stop evolution of a star or binary if an error occurs (including, optionally by a program option, floating-point errors), and record in the (SSE/BSE) system
parameters file the fact that an error occurred, and an error number identifying the error that occurred. This way users can check the system parameters file 
at the completion of a run for the disposition of a star or binary and, if the evolution of that star or binary was stopped because an error occurred, the 
actual error that occurred.

Users should refer to the Error Handling documentation in the User Guide (See :doc:`./User guide/Handling errors/handling-errors`).
Developers should refer to the Error Handling documentation in the Developer Guide (See :doc:`./Developer guide/Services/services-error-handling`).

In addition to the new error handling strategy, a number of program options have been deprecated and replaced by new program options that are more consistent with the
naming convention we are trying to maintain. The program options deprecated, and their replacements, are:

1. `--mass-transfer`                       , replaced by `--use-mass-transfer` (to be consistent with `--use-mass-loss`)
#. `--luminous-blue-variable-prescription` , replaced by `--LBV-mass-loss-prescription` (to be consistent with other 'prescription'-type options)
#. `--OB-mass-loss`                        , replaced by `--OB-mass-loss-prescription` (to be consistent with other 'prescription'-type options)
#. `--RSG-mass-loss`                       , replaced by `--RSG-mass-loss-prescription` (to be consistent with other 'prescription'- type options)
#. `--VMS-mass-loss`                       , replaced by `--VMS-mass-loss-prescription` (to be consistent with other 'prescription'- type options)
#. `--WR-mass-loss`                        , replaced by `--WR-mass-loss-prescription` (to be consistent with other 'prescription'- type options)
#. `--kick-direction`                      , replaced by `--kick-direction-distribution` (to be consistent with other 'distribution'-type options)
#. `--mass-transfer-thermal-limit-accretor`, replaced by `--mass-transfer-thermal-limit-accretor-multiplier` (for consistency and to better describe the option)
#. `--black-hole-kicks`                    , replaced by `--black-hole-kicks-mode` (for consistency and to better describe the option) 
#. `--chemically-homogeneous-evolution`    , replaced by `--chemically-homogeneous-evolution-mode` (for consistency and to better describe the option)

Deprecated program options will still be available, in tandem with their replacements, for some time (at least six months from the release date of v03.00.00),
after which time the deprecated options will be removed and only their replacements will be valid options. Appropriate warning messages will be displayed in
the period of deprecation if deprecated program options are used.

As well as deprecating some program options, some program option values have been deprecated and replaced by new values that are more consistent with the
naming convention we are trying to maintain. The program option value deprecated, and their replacements, are:

1. `--LBV-mass-loss-prescription`         , value `NONE` replaced by `ZERO`
#. `--luminous-blue-variable-prescription`, value `NONE` replaced by `ZERO`
#. `--OB-mass-loss`                       , value `NONE` replaced by `ZERO`
#. `--OB-mass-loss-prescription`          , value `NONE` replaced by `ZERO`
#. `--RSG-mass-loss`                      , value `NONE` replaced by `ZERO`
#. `--RSG-mass-loss-prescription`         , value `NONE` replaced by `ZERO`
#. `--VMS-mass-loss`                      , value `NONE` replaced by `ZERO`
#. `--VMS-mass-loss-prescription`         , value `NONE` replaced by `ZERO`
#. `--WR-mass-loss`                       , value `NONE` replaced by `ZERO`
#. `--WR-mass-loss-prescription`          , value `NONE` replaced by `ZERO`

Deprecated program option values will still be available, in tandem with their replacements, for some time (at least six months from the release date of v03.00.00),
after which time the deprecated option values will be removed and only their replacements will be valid option values. Appropriate warning messages will be displayed
in the period of deprecation if deprecated program option values are used.

Finally, for program option `--mt-rejuvenation-prescription`, the value `NONE` was replaced by `HURLEY`


**02.48.01 May 24, 2024**

* changed functionality of ``--output-path`` option so that missing directories in the specified path are created.
* Added "Quick Links" to online documentation.

**02.48.00 May 22, 2024**

* added options ``--mass-transfer-jloss-macleod-linear-fraction-degen`` and ``--mass-transfer-jloss-macleod-linear-fraction-non-degen`` to allow for different accretor AM response for degenerate and non-degenerate companions.
* removed option ``--mass-transfer-jloss-macleod-linear-fraction`` (no longer required - see above).

**02.46.00 May 13, 2024**

* added options ``--radial-change-fraction`` and ``--mass-change-fraction``, as approximate desired fractional changes in stellar radius and mass on phase when setting SSE and BSE timesteps
* the recommended values for both parameters are 0.005, but the default remains 0, which reproduces previous timestep choices
* mass transfer from main sequence donors (including HeMS) can now proceed on nuclear timescales -- approximated as the radial expansion timescales -- if equilibrium zetas are greater than Roche lobe zetas

**02.45.00 Apr 09, 2024**

* Changed compiler standard in Makefile from ``c++11`` to ``c++17``.  This is required for ``boost v1.82`` and above. ``c++11`` can still be used if boost version is below ``v1.82``, but moving to ``c++17`` and boost ``v1.8x`` is preferred (and will eventually be mandatory). Tested with ``Ubuntu v20.04, g++ v11.04, and boost v1.74``; and ``macOS v14.1.1, clang v15.0.0, and boost v1.85``.

**02.44.00 Apr 04, 2024**

* Added 'realistic' tides option, which implements dynamical and equilibrium tides using the formalism described in Kapil et al. (2024). 
* Functionality enabled with new option ``--tides-prescription KAPIL2024`` (default is ``NONE``).
* Removed old option ``--enable-tides``, which can now be enabled by setting ``--tides-prescription PERFECT``.


**02.43.00 Mar 29, 2024**

* Implementation of the neutrino rocket kick.

**02.42.00 Jan 04, 2023**

* Timesteps are now quantised to an integral multiple of 1e-12Myr.
* New option provided to allow user-defined timesteps: ``--timesteps-filename`` (See :doc:`./User guide/timestep-files`).
* Code changes to make SSE and BSE evolution more consistent (See `PR 1052 <https://github.com/TeamCOMPAS/COMPAS/pull/1052>`_).

**02.41.03 Dec 28, 2023**

* The functions ``BaseBinaryStar::CalculateAngularMomentum()``, ``BaseBinaryStar::CalculateTotalEnergy()``, and ``BaseStar::AngularMomentum()`` changed to use moment of inertia instead of gyration radius.
* Changed CalculateMomentOfInertia() to properly implement Hurley et al., 2000 eq 109.
* This change may change DCO yields slightly when compared to previous versions of the code.

**02.41.00 Nov 02, 2023**

* Added a naive tides implementation.
* Added program option ``enable-tides`` to enable the tides implementation (default is ``false``).

**02.40.00 Oct 20, 2023**

* Added ``FLEXIBLE2023`` as a new default, and ``BELCZYNSKI2010`` as a replacement for the previous ``VINK`` mass loss prescription. The following new sub-wrappers are overridden when selecting ``BELCZYNSKI2010``:
* Added ``--OB-mass-loss`` program option, applying to main sequence stars, with default ``VINK2021``, and options ``NONE``, ``VINK2001`` (previous default), ``BJORKLUND2022``, and ``KRTICKA2018``.
* Added ``--RSG-mass-loss`` program option, applying to stars below 8kK in giant branch stellar types, with default ``DECIN2023``, and options ``NONE``, ``VINISABHAHIT2023``, ``BEASOR2020``, ``YANG2023``, ``KEE2021``, ``NJ90`` (previous default).
* Added ``--VMS-mass-loss`` program option, applying to stars over 100 Msol, with default ``SABHAHIT2023``, and options ``NONE``, ``VINK2011``, and ``BESTENLEHNER2020``.
* Added ``--WR-mass-loss`` program option, with default ``SANDERVINK2023``, and options ``BELCZYNSKI2010``, and ``SHENAR2019``.
* Changed default value for option ``--wolf-rayet-multiplier`` from 0.1 to 1.0

**02.39.00 Jul 4, 2023**

* Added 'Evolution_Status' columns to both SSE and BSE default system parameters records - records final status of evolution (reason evolution stopped).

**02.38.03 Apr 20, 2023**

* Changed some of the default options, see issue # 957 and PR # 961 for explanations

**02.37.00 Mar 26, 2023**

* Added functionality for WDs to accrete in different regimes. 
* New supernova types: SNIA (Type Ia), and HeSD (Helium shell detonation). 

**02.36.00 Mar 15, 2023**

* Added functionality to automatically create COMPAS YAML file - adds two new options: ``--create-YAML-file`` and ``YAML-template``. See documentation for details.  

  **Note:** From this release, the default COMPAS YAML file (``compasConfigDefault.yaml``), as distributed, has all COMPAS option entries commented so that the COMPAS default value for the option is used by default. To use a value other than the COMPAS default value, users must uncomment the entry and change the option value to the desired value.

**02.35.03 Feb 27, 2023**

Added mass accretion prescription during CE ``CHEVALIER`` for option ``--common-envelope-mass-accretion-prescription``, following model 2 from van Son + 2020

**02.35.02 Feb 19, 2023**

* Changed ``BINARY_PROPERTY::ROCHE_LOBE_RADIUS_1`` and ``BINARY_PROPERTY::ROCHE_LOBE_RADIUS_2`` to be the Roche lobe radius as computed at periapsis, in units of :math:`R_\odot`.
* Changed header string for ``BINARY_PROPERTY::ROCHE_LOBE_RADIUS_1`` from ``'RocheLobe(1)|a'`` to ``'RocheLobe(1)'`` - same change made for ``BINARY_PROPERTY::ROCHE_LOBE_RADIUS_2``.
* Removed ``BINARY_PROPERTY::STAR_TO_ROCHE_LOBE_RADIUS_RATIO_1`` (header string ``'Radius(1)|RL'``) and ``BINARY_PROPERTY::STAR_TO_ROCHE_LOBE_RADIUS_RATIO_2`` (header string ``'Radius(2)|RL'``) from ``BSE_DETAILED_OUTPUT_REC`` (BSE detailed output file default record).  Note that both variables are still selectable for output via the logfile-definitions file.

  **Note:** These changes will affect post-processing code that consumes the affected variables - users should check their post-processing code. 

**02.35.00 Dec 8, 2022**

* Added critical mass ratios from Ge+ 2020 for determining if MT is unstable.

**02.34.01 Dec 7, 2022**

* Fixed Time<MT in BSE_RLOF, which previously was identical with Time>MT.

**02.33.00 Aug 28, 2022**

* Added simplified (constant per stellar type) critical mass ratios from Claeys+ 2014 for determining if MT is unstable

**02.32.00 Aug 27, 2022**

* Added 'record type' functionality to all standard log files.  **Note:** This changes default behaviour: only Detailed Output log files affected in this release
* Added/rationalised Detailed Output records printed for binary systems
* Added new program option for each standard log file to allow specification of which record types to print. See e.g. ``--logfile-detailed-output-record-types``
* Changed case on column header strings for switch log files (SSE and BSE. ``SWITCHING_FROM``, ``SWITCHING_TO``, and ``STAR_SWITCHING`` are now ``Switching_From``, ``Switching_To``, and ``Star_Switching`` respectively).   **Note:** This could affect post-processig code that consumes the switch log files - users should check that their code will recognise the new header strings.
* Added new section to online documentation: 'What's new'

**02.31.10 Aug 12, 2022**

* Added option to set the Temperature boundary between convective/radiative giant envelopes

**02.31.09 Aug 9, 2022**

* Max evolution time and max number of timesteps now read in from gridline as well as commandline

**02.31.08 Aug 3, 2022**

* Added Accretion Induced Collapse (AIC) of ONeWD as another type of SN

**02.31.07 Aug 1, 2022**

* Added print to DetailedOutput after merger, addresses https://github.com/TeamCOMPAS/COMPAS/issues/825
* Ensure no ONeWDs are formed with masses above Chandrasekhar mass

**02.31.06 Aug 2, 2022**

* Added stellar merger to default BSE_RLOF output

**02.31.05 July 25, 2022**

* Renamed program option ``--allow-H-rich-ECSN`` to ``allow-non-stripped-ECSN``
* Fixed check for non-interacting ECSN progenitors to consider MT history instead of H-richness

**02.31.04 Jun 10, 2022**

* Changed MT_TRACKER values to be clearer and complementary to each other
* Updated the relevant section in the detailed plotter that uses MT_TRACKER values
* Removed end states from detailed plotter (Merger, DCO, Unbound) so that they don't over compress the rest

**02.31.03 May 20, 2022**

* Fixed MS+MS unstable MT not getting flagged as a CEE

**02.31.00 May 14, 2022**

* Added new program option ``--retain-core-mass-during-caseA-mass-transfer`` to preserve a larger donor core mass following case A MT, set equal to the expected core mass of a newly formed HG star with mass equal to that of the donor, scaled by the fraction of its MS lifetime

**02.30.00 May 8, 2022**

* Added MACLEOD_LINEAR specific angular momentum gamma loss prescription for stable mass transfer (see ``--mass-transfer-angular-momentum-loss-prescription``)

**02.29.00 May 5, 2022**

* Added new program option to allow for H-rich ECSN (``--allow-H-rich-ECSN``, defaults to FALSE). When the option is TRUE, non-interacting ECSN progenitors do not contribute to the single pulsar population.  Addresses issue https://github.com/TeamCOMPAS/COMPAS/issues/596

**02.28.00 May 11, 2022**

* Added new remnant mass prescription: Fryer+ 2022
* Added new program options ``--fryer-22-fmix`` and ``--fryer-22-mcrit``

**02.27.09 Apr 25, 2022**

* Added new program option ``--muller-mandel-sigma-kick``

**02.27.08 Apr 12, 2022**

* Fix for issue https://github.com/TeamCOMPAS/COMPAS/issues/783

**02.27.07 Apr 5, 2022**

* Fix for issue https://github.com/TeamCOMPAS/COMPAS/issues/773

**02.27.06 Apr 5, 2022**

* Fixed StarTrack PPISN prescription: previously it was doing the same thing as the COMPAS PPISN prescription

**02.27.05 Feb 17, 2022**

* Added new program option ``--hmxr-binaries``, which tells COMPAS to store high-mass x-ray binaries in BSE_RLOF output file
* Added columns for pre- and post-timestep ratio of stars to Roche Lobe radius to BSE_RLOF output file (addressing issue https://github.com/TeamCOMPAS/COMPAS/issues/746)

**02.27.04 Feb 15, 2022**

* Fix for issue https://github.com/TeamCOMPAS/COMPAS/issues/761

**02.27.03 Feb 8, 2022**

* Fix for issue https://github.com/TeamCOMPAS/COMPAS/issues/745

**v02.27.02 Feb 3, 2022**

* Fixed mass change on forced envelope loss in response to issue https://github.com/TeamCOMPAS/COMPAS/issues/743

**v02.27.01 Feb 3, 2022**

* Fixed condition for envelope type when using ENVELOPE_STATE_PRESCRIPTION::FIXED_TEMPERATURE (previously, almost all envelopes were incorrectly declared radiative)

**v02.27.00 Jan 12, 2022**

* Added enhanced Nanjing lambda option that continuously extrapolates beyond radial range
* Added Nanjing lambda option to switch between calculation using rejuvenated mass and true birth mass
* Added Nanjing lambda mass and metallicity interpolation options
* No change in default behaviour

