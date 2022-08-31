What's new
==========

Following is an brief list of important updates to the COMPAS code.  A complete record of changes can be found in the file ``changelog.h``.


**LATEST RELEASE** |br|

**02.33.00 Aug 28, 2022**

* Added simplified (constant per stellar type) critical mass ratios from Claeys+ 2014 for determining if MT is unstable

**02.32.00 Aug 27, 2022**

* Added 'record type' functionality to all standard log files.  **Note:** This changes default behaviour: only Detailed Output log files affected in this release
* Added/rationalised Detailed Output records printed for binary systems
* Added new program option for each standard log file to allow specification of which record types to print. See e.g. '--logfile-detailed-output-record-types'
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

