# ifndef __changelog_h__
# define __changelog_h__

// =====================================================================
// 
// COMPAS Changelog
// 
// =====================================================================
// 
// 02.00.00      JR - Sep 17, 2019 - Initial commit of new version
// 02.00.01      JR - Sep 20, 2019 - Fix compiler warnings. Powwow fixes
// 02.00.02      JR - Sep 21, 2019 - Make code clang-compliant
// 02.00.03      IM - Sep 23, 2019 - Added fstream include
// 02.01.00      JR - Oct 01, 2019 - Support for Chemically Homogeneous Evolution
// 02.02.00      JR - Oct 01, 2019 - Support for Grids - both SSE and BSE
// 02.02.01      JR - Oct 01, 2019 - Changed BaseBinaryStar code to assume tidal locking only if CHE is enabled
// 02.02.02      JR - Oct 07, 2019 - Defect repairs:
//                                       SSE iteration (increment index - Grids worked, range of values wasn't incrementing)
//                                       Errors service (FIRST_IN_FUNCTION errors sometimes printed every time)
//                                       Added code for both SSE and BSE so that specified metallicities be clamped to [0.0, 1.0].  What are reasonable limits?
//                                   Errors service performance enhancement (clean deleted stellar objects from catalog)
//                                   Changed way binary constituent stars masses equilibrated (they now retain their ZAMS mass, but (initial) mass and mass0 changes)
//                                   Add initial stellar type variable - and to some record definitions
//                                   Added change history and version number to constants.h
// 02.02.03      JR - Oct 09, 2019 - Defect repairs:
//                                       Initialised all BaseStar.m_Supernova elements (some had not been initialised)
//                                       Fixed regression in BaseStar.cpp (INITIAL_STELLAR_TYPE & INITIAL_STELLAR_TYPE_NAME in StellarPropertyValue())
//                                   Added max iteration check to Newton-Raphson method in SolveKeplersEquation (see constant MAX_KEPLER_ITERATIONS)
// 02.02.04      JR - Oct 09, 2019 - Defect repairs:
//                                       SN kick direction calculation corrected
//                                       Boolean value output corrected
//                                       Typos fixed
// 02.02.05      JR - Oct 10, 2019 - Defect repairs:
//                                       Determination of Chemically Homogeneous star fixed (threshold calculation)
//                                       Removed checks for RLOF to/from CH stars
//                                       Typos fixed
// 02.02.06      JR - Oct 11, 2019 - Renamed class "CHE" - now class "CH"
//                                   Updated CHE documentation
//                                   Added m_MassesEquilibrated variable to BaseBinaryStar
// 02.02.07      JR - Oct 20, 2019 - Defect repairs:
//                                       CEE printing systems post-stripping - github issue - reworked CE details/pre/post CE - partial fix (BindingEnergy remaining)
//                                       Added RANDOM_SEED to Options::OptionValue() (omitted erroneously)
//                                   Added m_SecondaryTooSmallForDCO variable to BaseBinaryStar - and to some record definitions
//                                   Added m_StellarMergerAtBirth variable to BaseBinaryStar - and to some record definitions
//                                   Added allow-rlof-at-birth program option
//                                       If CHE enabled, or allow-rlof-at-birth option is true, binaries that have one or both stars
//                                       in RLOF at birth will have masses equilibrated, radii recalculated, orbit circularised, and
//                                       semi-major axis recalculated, while conserving angular momentum - then allowed to evolve
//                                   Added allow-touching-at-birth program option
//                                       Binaries that have stars touching at birth (check is done after any equilibration and
//                                       recalculation of radius and separation is done) are allowed to evolve.  Evolve() function
//                                       immediately checks for merger at birth, flags status as such and stops evolution.
//                                   Documentation updated (see updated doc for detailed explanation of new program options)
// 02.03.00      JR - Oct 25, 2019 - Defect repairs:
//                                       removed extraneous delimiter at end of log file records
//                                   Added '--version' option
//                                   Changed minor version number - should have been done at last release - we'll grant the '--version' option minor release status...
// 02.03.01      JR - Nov 04, 2019 - Defect repair:
//                                       removed erroneous initialisation of m_CEDetails.alpha from BaseBinaryStar::SetRemainingCommonValues()
//                                       (CE Alpha was always being initialised to 0.0 regardless of program options)
// 02.03.02      JR - Nov 25, 2019 - Defect repairs:
//                                       added check for active log file before closing in Log::Stop()
//                                       added CH stars to MAIN_SEQUENCE and ALL_MAIN_SEQUENCE initializer_lists defined in constants.h
//                                       moved InitialiseMassTransfer() outside 'if' - now called even if not using mass transfer - sets some flags we might need
//                                       added code to recalculate rlof if CH stars are equilibrated in BaseBinaryStar constructor
//                                   Enhancements:
//                                       moved KROUPA constants from AIS class to constants.h
//                                       moved CalculateCDFKroupa() function from AIS class to BaseBinaryStar class
//                                       added m_CHE variable to BaseStar class - also selectable for printing
//                                       added explicit check to ResolveCommonEnvelope() to merge binary if the donor is a main sequence star
//                                   Chemically Homogeneous Evolution changes:
//                                       added check to CheckMassTransfer() in BaseBinaryStar.cpp to merge if CH+CH and touching - avoid CEE
//                                       added code to InitialiseMassTransfer() in BaseBinaryStar.cpp to equilibrate and possibly merge if both CH stars in RLOF
// (Unchanged)   IM - Nov 29, 2019 - Defect repairs:
//                                       changed Disbound -> Unbounded in header strings in constants.h
//                                       left one line in default/example grid file (Grid.txt)
//                                       fix default PPISN mass limit in python submit: 65 Msol -> 60 Msol
// 02.03.03      JR - Dec 04, 2019 - Defect repairs:
//                                       added code to UpdateAttributesAndAgeOneTimestep() in Star.cpp to recalculate stellar attributes after switching to new stellar type
//                                       (addresses discontinuous transitions e.g. CH -> HeMS)
//                                       changed IsPulsationalPairInstabilitySN() in GiantBranch.cpp to call IsPairInstabilitySN() instead of set MASSLESS_REMNANT if remnant mass <= 0.0
//                                       changed CalculateSNKickVelocity() in BaseStar.cpp to set m_SupernovaDetails.kickVelocity correctly after adjusting for fallback
// 02.03.04      FSB - Dec 04, 2019 - Defect repairs:
//                                       fixed bug in Fryer+2012 CalculateGravitationalRemnantMassadded() function to compare baryon mass of star remnant with
//  									                   baryon mass of MaximumNeutronStarMass instead of just MaximumNeutronStarMass. 
//                                       added m_BaryonicMassOfMaximumNeutronStarMass to BaseStar.h and BaseStar.cpp
// 02.03.05      JR - Dec 05, 2019 - Defect repairs:
//                                       fixed EvolveSingleStars() in main.cpp to print correct initial mass
//                                       fixed TPAGB::CalculateCOCoreMassAtPhaseEnd() - added conditional
// 02.04.00      JR - Dec 18, 2019 - New functionality:
//                                       added columns to BSE grid functionality: Kick_Velocity_1(&2), Kick_Theta_1(&2), Kick_Phi_1(&2), Kick_Mean_Anomaly_1(&2).  Updated documentation.
//                                   Changed functionality:
//                                       removed compiler version checks from Makefile - they seemed to only work for native Ubuntu and were more of a nuisance than anything...  (old version exists as Makefile-checks)
//                                   Defect repairs:
//                                       added recalculation of gbParams Mx & Lx in HeHG calculateGbParams()
//                                       created HeHG::CalculateGBParams_Static() and GiantBranch::CalculateGBParams_Static(), called from EAGB::ResolveEnvelopeLoss() to facilitate calculation of attributes for new stellar type before actually switching.  Needed to rewrite some other functions as static.  Note: this needs to be revisited and a more elegant solution implemented.
//                                       added CalculateRadiusAndStellarTypeOnPhase() for HeHG and HeGBstars, and changed call to calculateRadiusOnPhase() to CalculateRadiusAndStellarTypeOnPhase() in BaseStar::EvolveOnPhase().  This allows for HeHG and HeGB stars to change stellar type based on radius (previously missed).
//                                       set M = McBAGB for EAGB & TPAGB only (was being set for all types >= TPAGB)
//                                       added extra print detailed in BaseBinaryStar:Evolve() - sometimes missing a switch type in detailed output if only 1 timestep
//                                       swapped heading strings for ANY_STAR_PROPERTY::IS_ECSN and ANY_STAR_PROPERTY::IS_USSN (now correct)
//                                       removed condition in BaseBinaryStar::EvaluateSupernovae().  ResolveSupernova() is now called for all stellar types (not sure what I was thinking originally. I'm sure I had a good reason - or maybe I was just tired...)
//                                       changed name of GiantBranch::CalculateProtoCoreMass() to GiantBranch::CalculateProtoCoreMassDelayed() and changed calls to the function
//                                       swapped order of calculations of ePrime (CalculateOrbitalEccentricityPostSupernova()) and m_SemiMajorAxisPrime (CalculateSemiMajorAxisPostSupernova()) in BaseBinaryStar::ResolveSupernova().  Improper order was causing wrong value of m_SeminMajorAxisPrime to be used in calculation of ePrime
//                                       set m_Disbound = true appropriately in BaseBinaryStar::Evolve() (note: m_Disbound will change name to m_Unbound soon...)
//                                       changed return value of CHeB::DetermineEnvelopeType() to CONVECTIVE.  Left CHeB DetermineEnvelopeTypeHurley2002() as RADIATIVE (used in BinaryConstituentStar::CalculateSynchronisationTimescale())
//                                       changed BINARY_PROPERTY::ORBITAL_VELOCITY to BINARY_PROPERTY::ORBITAL_VELOCITY_PRE_2ND_SUPERNOVA in BSE_SUPERNOVAE_REC (6th value printed)
//                                       added p_Erase parameter to Log::CloseStandardFile(); changed Log::CloseAllStandardFiles() to call Log::CloseStandardFile() with p_Erase=false and erase entire map after all files closed (prevent coredump when closing all files)
//                                       added ResolveSupernova() to ONeWD.h - ONeWD stars were previously not checking for SN
//                                       fixed BaseBinaryStar::InitialiseMassTransfer() - star1 was being updated instead of star2 for CH + CH stars when CHE enabled
// 02.04.01      JR - Dec 23, 2019 - Defect repairs:
//                                       Removed SN_EVENT::SN - all occurrences of SN_EVENT::SN replaced by SN_EVENT::CCSN.
//                                           The current SN event ("Is"), and past SN event ("Experienced") are now bit maps (implemented as Enum Classes).  Each can have any of the values: CCSN, ECSN, PISN, PPSIN, USSN, RUNAWAY, RECYCLED_NS, and RLOF_ONTO_NS.  See definition of SN_EVENT Enum Class in constants.h for implementation and explanation.  
//                                       Updated variables selectable for printing:
//                                           Added ANY_STAR_PROPERTY::SN_TYPE (STAR_PROPERTY, SUPERNOVA_PROPERTY, COMPANION_PROPERTY (should always be SN_EVENT::NONE for companion star))
//                                           Added ANY_STAR_PROPERTY::EXPERIENCED_SN_TYPE (STAR_PROPERTY, SUPERNOVA_PROPERTY, COMPANION_PROPERTY)
//                                           All of ANY_STAR_PROPERTY::{CCSN, ECSN, PISN, PPISN, USSN} now selectable
//                                           Removed ANY_STAR_PROPERTY::SN - no longer selectable for printing (replaced by CCSN)
//                                           Updated documentation
//                                       Changed default record specifications for logfiles BSE_DOUBLE_COMPACT_OBJECTS_REC and BSE_SUPERNOVAE_REC
//                                           Removed the individual SN_EVENT columns for both "Is" and "Experienced" conditions (e.g. CCSN, ECSN etc)
//                                           "Is*" and "Experienced*" columns replaced with SN_TYPE & Experienced_SN_TYPE columns that record the SN event type (e.g. CCSN, ECSN, PPSN, PPSIN, USSN).  
//                                           RUNAWAY, RECYCLED_NS, and RLOF_ONTO_NS are still reported in separate, individual columns.
//                                       Added workaround for non-existent CHeB blue loop.  See description in CHeB::CalculateTimescales()
//                                       Removed binary star "survived" flag - it is always the NOT of the "unbound" flag
//                                       Changed initialisation function for HeGB stars (HeGB::Initialise() in HeGB.h) to NOT recalculate m_Age if evolving from HeHG -> HeGB 
//                                       Removed initialisation of m_Age (to 0.0) from COWD::Initialise() in COWD.h
//                                   Changed behaviour:  
//                                       Changed binary star "disbound" flag to "unbound" flag.  Changed all occurrences of "disbound" to "unbound".  Changed "unbound" header flag to "Unbound"
// 02.04.02      JR - Jan 06, 2020 - Defect repairs:
//                                       Added IsPISN() & IsPPISN() to IsSNEvent()
//                                       Fixed check for SN event at top of BaseBinaryStar::ResolveSupenova()
//                                       Changed BaseBinaryStar::EvaluateSupernovae() to more closely match legacy code behaviour (see notes in function description):
//                                          Added p_Calculate2ndSN parameter to determine if 2nd supernova needs to be resolved
//                                          Clear star2 current SN event if necessary
//                                          Check m_SemiMajorAxisPrime value prior to SN events (viz. new aPrime variable)
//                                       Fixed timestep initialisation in BaseStar::CalculateConvergedTimestepZetaNuclear()  (was negative)
//                                       Fixed m_Age calculation in FGB::ResolveEnvelopeLoss()
//                                       Added CalculateInitialSupernovaMass() to NS.h - was setting M = 5.0 for >= ONeWD, should be ONeWD only (introduced in fix in v04.02.00)
//                                       Changed NS functions to return Radius in Rsol instead of km:
//                                          Added function NS:CalculateRadiusOnPhaseInKM_Static() (returns radius in km)
//                                          Changed NS:CalculateRadiusOnPhase_Static() to return Rsol
//                                          Added CalculateRadiusOnPhase() for NS (ns.h) - returns Rsol 
//                                   Changed behaviour:  
//                                       Print detailed output record whenever stellartype changes (after star 2 if both change)
// (Unchanged)   LK - Jan 10, 2020 - Defect repairs:
//                                       Added missing includes to Star.cpp, utils.h and utils.cpp (required for some compiler versions)
// 02.05.00      JR - Jan 23, 2020 - New functionality:
//                                       Grid files:
//                                          Added kick velocity magnitude random number to BSE grid file - see docs re Grids
//                                          Added range check for Kick_Mean_Anomaly_1 and Kick_Mean_Anomaly_2 ([0.0, 2pi)) in BSE grid file
//                                          Cleaned up SSE & BSE grid file code
//                                       Added m_LBVphaseFlag variable to BaseStar class; also added ANY_STAR_PROPERTY::LBV_PHASE_FLAG print variable.
//                                   Deleted functionality:  
//                                       Removed IndividualSystem option and related options - this can now be achieved via a grid file
//                                          Update pythonSubmitDefault.py to remove individual system related parameters
//                                   Changed behaviour:
//                                       Removed check for Options->Quiet() around simulation ended and cpu/wall time displays at end of EvolveSingleStars() and EvolveBinaryStars() in main.cpp
//                                   Defect repairs:
//                                       Removed erroneous check for CH stars in BaseBinaryStar::EvaluateBinary()
//                                       Fix for issue #46 (lower the minimum value of McSN in star.cpp from Mch to 1.38)
//                                          Changed 'MCH' to 'MECS' in 
//                                             BaseStar::CalculateMaximumCoreMassSN()
//                                             GiantBranch::CalculateCoreMassAtSupernova_Static
// 02.05.01      FSB - Jan 27, 2020 -Enhancement:
//                                       Cleaned up default printed headers and parameters constants.h:
//                                           - removed double parameters that were printed in multiple output files 
//                                           - changed some of the header names to more clear / consistent names
//                                           - added some comments in the default printing below for headers that we might want to remove in the near future
// 02.05.02      JR - Feb 21, 2020 - Defect repairs:
//                                       - fixed issue #31: zRocheLobe function does not use angular momentum loss
//                                       - fixed default logfile path (defaulted to '/' instead of './')
//                                       - changed default CE_ZETA_PRESCRIPTION to SOBERMAN (was STARTRACK which is no longer supported)
// 02.05.03      JR - Feb 21, 2020 - Defect repairs:
//                                       - removed extraneous debug print statement from Log.cpp
// 02.05.04      JR - Feb 23, 2020 - Defect repairs:
//                                       - fixed regression introduced in v02.05.00 that incread DNS rate ten-fold
//                                           - changed parameter from m_SupernovaDetails.initialKickParameters.velocityRandom to m_SupernovaDetails.kickVelocityRandom in call to DrawSNKickVelocity() in BaseStar::CalculateSNKickVelocity()
//                                       - reinstated STAR_1_PROPERTY::STELLAR_TYPE and STAR_2_PROPERTY::STELLAR_TYPE in BSE_SYSTEM_PARAMETERS_REC
// 02.05.05      JR - Feb 27, 2020 - Defect repair:
//                                       - fixed age resetting to 0.0 for MS_GT_07 stars after CH star spins down and switches to MS_GT_07
//                                           - ensure m_Age = 0.0 in constructor for BaseStar
//                                           - remove m_Age = 0.0 from Initialise() in MS_gt.07.h 
// 02.05.06      JR - Mar 02, 2020 - Defect repair:
//                                       - fixed m_MassesEquilibrated and associated functions - was erroneously typed as DOUBLE - now BOOL
//                                   Added/changed functionality:
//                                       - added m_MassesEquilibratedAtBirth variable to class BaseBinaryStar and associated property BINARY_PROPERTY::MASSES_EQUILIBRATED_AT_BIRTH
//                                       - tidied up pythonSubmitDefault.py a little:
//                                             - set grid_filename = None (was '' which worked, but None is correct)
//                                             - set logfile_definitions = None (was '' which worked, but None is correct)
//                                             - added logfile names - set to None (COMPAS commandline arguments already exist for these - introduced in v02.00.00)
// 02.05.07      JR - Mar 08, 2020 - Defect repair:
//                                       - fixed circularisation equation in BaseBinaryStar::InitialiseMassTransfer() - now takes new mass values into account
// 02.06.00      JR - Mar 10, 2020 - Changed functionality:
//                                       - removed RLOF printing code & associated pythonSubmitDefault.py options
// 02.06.01      JR - Mar 11, 2020 - Defect repair:
//                                       - removed extraneous debug print statement from Log.cpp (was previously removed in v02.05.03 but we backed-out the change...)
// 02.06.02      JR - Mar 15, 2020 - Defect repairs:
//                                       - removed commented RLOF printing lines in constant.h (somehow that was lost in some out of sync git merges...)
//                                       - removed commented options no longer used from Options.h and Options.cpp
//                                       - fixed units headers in constants.h - there are now no blank units headers, so SPACE delimited files now parse ok (multiple spaces should be treated as a single space)
//                                       - changed file extension for TAB delimited files to 'tsv'
//                                       - removed "useImportanceSampling" option - not used in code
//                                       - fixed typo in zeta-calculation-every-timestep option in Options.cpp
//                                       - removed redundant OPTIONS->MassTransferCriticalMassRatioHeliumGiant() from qcritflag if statement in BaseBinaryStar::CalculateMassTransfer()
//                                       - fixed OPTIONS->FixedMetallicity() - always returned true, now returns actual value
//                                       - fixed OPTIONS->OutputPathString() - was always returning raw option instead of fully qualified path
//                                       - changed the following in BaseBinaryStar::SetRemainingCommonValues() - erroneously not ported from legacy code:
//                                           (a) m_JLoss = OPTIONS->MassTransferJloss();
//                                           (b) m_FractionAccreted = OPTIONS->MassTransferFractionAccreted();
//                                           (both were being set to default value of 0.0)
//                                       - added OPTIONS->ZetaAdiabaticArbitrary() - option existed, but Options code had no function to retrieve value
//                                       - added OPTIONS->MassTransferFractionAccreted() to options - erroneously not ported from legacy code
//                                   Changed functionality:
//                                       - all options now have default values, and those values will be displayed in the help text (rather than string constants which may be incorrect)
//                                       - boolean options can now be provided with an argument (e.g. --massTransfer false)
//                                       - added ProgramOptionDetails() to Options.cpp and OPTIONS->OptionsDetails() in preparation for change in output functionality
// 02.07.00      JR - Mar 16, 2020 - New/changed functionality:
//                                       - COMPAS Logfiles are created in a (newly created) directory - this way they are all kept together
//                                       - new command line option 'output-container' implemented (also in pythonSubmitDefault.py) - this option allows the user to specify the name of the log files container directory (default is 'COMPAS_Output')
//                                       - if detailed log files are created they will be created in a directory named 'Detailed_Output' within the container directory
//                                       - a run details file named 'Run_details' is created in the container directory.  The file records the run details:
//                                             - COMPAS version
//                                             - date & time of run
//                                             - timing details (wall time, CPU seconds)
//                                             - the command line options and parameters used:
//                                                   - the value of options and an indication of whether the option was supplied by the user or the default value was used
//                                                   - other parameters - calculated/determined - are recorded
//                                   Defect repair:
//                                       - changed "--outut" option name to "--outpuPath" in stringCommands in pythonSubmitDefault.py
// 02.08.00		  AVG - Mar 17, 2020 - Changed functionality:
//  									                   - removed post-newtonian spin evolution	code & associated pythonSubmitDefault.py options
//  									                   - removed only_double_compact_objects code & associated pythonSubmitDefault.py options
//  									                   - removed tides code & associated pythonSubmitDefault.py options
//  									                   - removed deprecated options from pythonSubmitDefault.py options
//  									                   - renamed options: mass transfer, iterations -> timestep-iterations
//  									                   - commented AIS Options until fully implemented
// 02.08.01      JR - Mar 18, 2020 - Defect repairs:
//                                      - restored initialisation of AIS options in Options.cpp (AIS now defaults off instead of on)
//                                      - fixed retrieval of values for:
//                                            - ANY_STAR_PROPERTY::LAMBDA_KRUCKOW_BOTTOM, 
//                                            - ANY_STAR_PROPERTY::LAMBDA_KRUCKOW_MIDDLE, and 
//                                            - ANY_STAR_PROPERTY::LAMBDA_KRUCKOW_TOP 
//                                         in BaseStar::StellarPropertyValue().  Were all previously retrieving same value as ANY_STAR_PROPERTY::LAMBDA_KRUCKOW
//                                      - fixed some comments in BAseBinaryStar.cpp (lines 2222 and 2468, "de Mink" -> "HURLEY")
//                                      - fixed description (in comments) of BinaryConstituentStar::SetPostCEEValues() (erroneously had "pre" instead of "post" - in comments only, not code)
//                                      - fixed description of BaseStar::DrawKickDirection()
// 02.08.02      JR - Mar 27, 2020 - Defect repairs:
//                                      - fixed issue #158 RocheLobe_1<CE == RocheLobe_2<CE always
//                                      - fixed issue #160 Circularisation timescale incorrectly calculated
//                                      - fixed issue #161 Splashscreen printed twice - now only prints once
//                                      - fixed issue #162 OPTIONS->UseFixedUK() always returns FALSE.  Now returns TRUE if user supplies a fixed kick velocity via --fix-dimensionless-kick-velocity command line option
// 02.08.03      JR - Mar 28, 2020 - Defect repairs:
//                                      - fixed typo in BaseBinaryStar::ResolveCommonEnvelopeEvent() when calculating circularisation timescale in the case where star2 is the donor: star1Copy was erroneously used instead of star2Copy; changed to star2Copy
//                                      - changed circularisation timescale of binary to be minimum of constituent stars circularisation timescales, clamped to (0.0, infinity)
// 02.09.00      JR - Mar 30, 2020 - Minor enhancements:
//                                      - tightened the conditions under which we allow over-contact binaries - enabling CHE is no longer a sufficient condition after this change: the allow-rlof-at-birth option must also be specified (ussue #164)
//                                      - added printing of number of stars (for SSE) or binaries (for BSE) created to both stdout and Run_Details (issue #165)
//                                      - enhanced grid processing code in main.cpp to better handle TAB characters
// 02.09.01      JR - Mar 30, 2020 - Defect repair:
//                                      - OPTIONS->UseFixedUK() returns TRUE when user supplies -ve value via --fix-dimensionless-kick-velocity.  Now return TRUE iff the user supplies a value >=0 via --fix-dimensionless-kick-velocity
// 02.09.02      DC - Mar 30, 2020 - Defect repairs:
//                                      - Pulsar code fixed by correcting unit of NS radius in NS.cpp (added KM_TO_M constant in constants.h as a part of this),
//                                      correcting initialisation of pulsar birth parameters from GiantBranch.cpp to NS.cpp, adding an extra condition for isolated evolution when the companion loses mass but the NS does not accrete 
//                                      - option MACLEOD was printing wrongly as MACLEOD+2014 for user options, hence corrected it to MACLEOD in Options.cpp
// 02.09.03      JR - Apr 01, 2020 - Defect repairs:
//                                      - reinstated assignment of "prev" values in BaseBinaryStar::EvaluateBinary() (where call to ResolveTides() was removed).  Fixes low DNS count introduced in v02.08.00 caused by removal of ResolveTides() function (and call)
//                                      - commented option --logfile-BSE-be-binaries to match Be-Binary options commented by AVG in v02.08.00
// 02.09.04      JR - Apr 03, 2020 - Defect repair:
//                                      - removed IsUSSN() from IsSNEvent() definition in BinaryConstituentStar.cpp (USSN flag indicates just US, not USSN. Needs to be tidied-up properly)
// 02.09.05	     IM - Apr 03, 2020 - Defect repair:
//  		                            - fixed timescale calculation issue for newly created HeHG stars (from stripped EAGB stars); fixes drop in CO core mass
// 02.09.06      JR - Apr 07, 2020 - Defect repair:
//                                      - corrected calculation in return statement for Rand::Random(const double p_Lower, const double p_Upper) (issue #201)
//                                      - corrected calculation in return statement for Rand::RandomInt(const double p_Lower, const double p_Upper) (issue #201)
// 02.09.07      SS - Apr 07, 2020 - Change eccentricity, semi major axis and orbital velocity pre-2nd supernove to just pre-supernova everywhere in the code
// 02.09.08      SS - Apr 07, 2020 - Update zetaMainSequence=2.0 and zetaHertzsprungGap=6.5 in Options::SetToFiducialValues
// 02.09.09      JR - Apr 11, 2020 - Defect repair:
//                                      - restored property names in COMPASUnorderedMap<STAR_PROPERTY, std::string> STAR_PROPERTY_LABEL in constants.h (issue #218) (was causing logfile definitions files to be parsed incorrectly)
// 02.09.10	     IM - Apr 12, 2020 - Minor enhancement: added Mueller & Mandel 2020 remnant mass and kick prescription, MULLERMANDEL
//  			                     Defect repair: corrected spelling of output help string for MULLER2016 and MULLER2016MAXWELLIAN
// 02.10.01	     IM - Apr 14, 2020 - Minor enhancement: 
//  				                            - moved code so that SSE will also sample SN kicks, following same code branch as BSE 
// 02.10.02      SS - Apr 16, 2020 - Bug Fix for issue #105 ; core and envelope masses for HeHG and TPAGB stars
// 02.10.03      JR - Apr 17, 2020 - Defect repair:
//                                      - added LBV and WR winds to SSE (issue #223)
// 02.10.04	     IM - Apr 25, 2020 - Minor enhancement: moved Mueller & Mandel prescription constants to constants.h, other cleaning of this option
// 02.10.05      JR - Apr 26, 2020 - Enhancements:
//                                      - Issue #239 - added actual random seed to Run_Details
//                                      - Issue #246 - changed Options.cpp to ignore --single-star-mass-max if --single-star-mass-steps = 1.  Already does in main.cpp.
// 02.10.06      JR - Apr 26, 2020 - Defect repair:
//                                      - Issue #233 - corrected cicularisation formalae used in both BaseBinartStar constructors
// 02.11.00      JR - Apr 27, 2020 - Enhancement:
//                                      - Issue #238 - add supernova kick functionality to SSE grid file (+ updated docs)
//                                   Defect repairs:
//                                      - fixed typo in Options.h: changed '#include "rand.h" to '#include "Rand.h"
//                                      - fixed printing of actual random seed in Run_Details file (moved to Log.cpp from Options.cpp: initial random seed is set after options are set)
// 02.11.01	     IM - May 20, 2020 - Defect repair: 
//                                      - changed max NS mass for MULLERMANDEL prescription to a self-consistent value
// 02.11.02      IM - Jun 15, 2020 - Defect repair:
//                                      - added constants CBUR1 and CBUR2 to avoid hardcoded limits for He core masses leading to partially degenerate CO cores
// 02.11.03     RTW - Jun 20, 2020 - Enhancement:
//                                      - Issue #264 - fixed mass transfer printing bug 
// 02.11.04      JR - Jun 25, 2020 - Defect repairs:
//                                      - Issue #260 - Corrected recalculation of ZAMS values after eqilibration and cicularisation at birth when using grid files
//                                      - Issue #266 - Corrected calculation in BaseBinaryStar::SampleInitialMassDistribution() for KROUPA IMF distribution
//                                      - Issue #275 - Previous stellar type not set when stellar type is switched mid-timestep - now fixed
// 02.11.05      IM - Jun 26, 2020 - Defect repair:
//  				                    - Issue #280 - Stars undergoing RLOF at ZAMS after masses are equalised were removed from run even if AllowRLOFatZAMS set
// 02.12.00      IM - Jun 29, 2020 - Defect repair:
//                                      - Issue 277 - move UpdateAttributesAndAgeOneTimestepPreamble() to after ResolveSupernova() to avoid inconsistency
// 02.12.01      IM - Jul 18, 2020 - Enhancement:
//                                      - Starting to clean up mass transfer functionality
// 02.12.02      IM - Jul 23, 2020 - Enhancement:
//                                      - Change to thermal timescale MT for both donor and accretor to determine MT stability
// 02.12.03      IM - Jul 23, 2020 - Enhancement:
//                                      - Introduced a new ENVELOPE_STATE_PRESCRIPTION to deal with different prescriptions for convective vs. radiative envelopes (no actual behaviour changes yet for ENVELOPE_STATE_PRESCRIPTION::LEGACY);
//                                      - Removed unused COMMON_ENVELOPE_PRESCRIPTION
// 02.12.04      IM - Jul 24, 2020 - Enhancement:
//                                      - Changed temperatures to be written in Kelvin (see issue #278)
// 02.12.05      IM - Jul 25, 2020 - Enhancement:
//                                      - Added definition of FIXED_TEMPERATURE prescription to DetermineEnvelopeType()
//                                      - Removed unnecessary (and inaccurate) numerical zeta Roche lobe calculation
// 02.12.06      IM - Jul 26, 2020 - Enhancement:
//                                      - Extended use of zetaRadiativeEnvelopeGiant (formerley zetaHertzsprungGap) for all radiative envelope giant-like stars
// 02.12.07      IM - Jul 26, 2020 - Defect repair:
//                                      - Issue 295: do not engage in mass transfer if the binary is unbound
// 02.12.08   	AVG - Jul 26, 2020 - Defect repair:
//                                      - Issue #269: legacy bug in eccentric RLOF leading to a CEE
// 02.12.09      IM - Jul 30, 2020 - Enhancement:
//                                      - Cleaning of BaseBinaryStar::CalculateMassTransferOrbit(); dispensed with mass-transfer-prescription option
// 02.13.00      IM - Aug 2, 2020  - Enhancements and defect repairs:
//                                      - Simplified timescale calculations in BaseBinaryStar
//                                      - Replaced Fast Phase Case A MT and regular RLOF MT from non-envelope stars with a single function based on a root solver rather than random guesses (significantly improves accuracy)
//                                      - Removed all references to fast phase case A MT
//                                      - Corrected failure to update stars in InitialiseMassTransfer if orbit circularised on mass transfer
//                                      - Corrected incorrect timestep calculation for HeHG stars
// 02.13.01     AVG - Aug 6, 2020  - Defect repair:
//  									- Issue #267: Use radius of the star instead of Roche-lobe radius throughout ResolveCommonEnvelopeEvent()
// 02.13.02      IM - Aug 8, 2020  - Enhancements and defect repairs:
//                                      - Simplified random draw from Maxwellian distribution to use gsl libraries
//                                      - Fixed mass transfer with fixed accretion rate
//                                      - Cleaned up code and removed unused code
//                                      - Updated documentation
// 02.13.03       IM - Aug 9, 2020  - Enhancements and defect repairs:
//                                      - Use total core mass rather than He core mass in calls to CalculateZAdiabtic (see Issue #300)
//                                      - Set He core mass to equal the CO core mass when the He shell is stripped (see issue #277)
//                                      - Ultra-stripped SNe are set at core collapse (do not confusingly refer to stripped stars as previously, see issue #189)
// 02.13.04       IM - Aug 14, 2020 - Enhancements and defect repairs:
//                                      - Catch exception in boost root finder for mass transfer (resolve issue #317)
//                                      - Update core masses during Initialisation of HG and HeHG stars to be consistent with Hurley models
//                                      - Avoid division by zero in mass transfer rates of WDs
//                                      - Remove POSTITNOTE remnant mass prescription
// 02.13.05       IM - Aug 16, 2020 - Enhancements and defect repairs:
//                                      - General code cleaning
//                                      - Removed some redundant variables (e.g., m_EnvMass, which can be computed from m_Mass and m_CoreMass)
//                                      - Removed calculations of ZetaThermal and ZetaNuclear (these were previously incorrect because they relied on the evolution of a stellar copy which reverted to BaseStar and therefore didn't have the correct behaviour)
//                                      - Fixed CalculateZadiabatic to use ZetaAdiabaticArbitrary rather than ZetaThermalArbitrary; removed the latter
//                                      - Capped He core mass gain during shell H burning for CHeB and TPAGB stars, whose on-phase evolution now ends promptly when this limit is reached; this change also resolves issue #315 (higher mass SN remnants than total stellar mass)
// 02.13.06     AVG - Aug 20, 2020  - Defect repair:
//  									- Issue #229: Corrected fitting parameters in Muller 16 SN kick function
// 02.13.07      IM - Aug 20, 2020  - Enhancements:
//                                      - ONeWDs can now undergo ECSN if their mass rises above MECS=1.38 solar masses (previously, they could only undergo CCSN on rising above 1.38 solar masses).  ONeWD::CalculateInitialSupernovaMass now returns MCBUR1 rather than 5.0 to ensure this happens
//                                      - BaseStar::CalculateMaximumCoreMassSN() has been removed - it is superfluous since  GiantBranch::CalculateCoreMassAtSupernova_Static does the same thing
//                                      - Some misleading comments in TPAGB dealing with SNe have been clarified
//                                      - Option to set MCBUR1 [minimum core mass at base of the AGB to avoid fully degenerate CO core formation] to a value different from the Hurley default of 1.6 solar masses added, Issue #65 resolved
//                                      - Removed unused Options::SetToFiducialValues()
//                                      - Documentation updated
// 02.13.08       JR - Aug 20, 2020 - Code cleanup:
//                                      - moved BaseStar::SolveKeplersEquation() to utils
//                                      - changed call to (now) utils::SolveKeplersEquation() in BaseStar::CalculateSNAnomalies() to accept tuple with error and show error/warning as necessary
//                                      - removed call to std::cerr from utils::SolveQuadratic() - now returns error if equation has no real roots
//                                      - changed call to utils::SolveQuadratic() in GiantBranch::CalculateGravitationalRemnantMass() to accept tuple with error and show warning as necessary
//                                      - changed RadiusEqualsRocheLobeFunctor() in BinaryBaseStar.h to not use the SHOW_WARN macro (can't uset ObjectId() function inside a templated function - no object)
//                                      - changed COMMANDLINE_STATUS to PROGRAM_STATUS (better description)
//                                      - moved ERROR:NONE to top of enum in constants.h (so ERROR = 0 = NONE - makes more sense...)
//                                      - added new program option '--enable-warnings' to enable warning messages (via SHOW_WARN macros).  Default is false.  SHOW_WARN macros were previously #undefined
// 02.13.09     RTW - Aug 21, 2020  - Code cleanup:
// 									    - Created changelog.txt and moved content over from constants.h
// 									    - Changed OrbitalVelocity to OrbitalAngularVelocity where that parameter was misnamed
// 									    - Changed Pre/PostSNeOrbitalVelocity to OrbitalVelocityPre/PostSN for consistency
// 									    - Added and updated physical conversion constants for clarity (e.g MSOL to MSOL_TO_KG)
// 									    - Removed ID from output files, it is confusing and superseded by SEED
// 									    - Removed 'Total' from TotalOrbital(Energy/AngularMomentum)
// 									    - Typos
// 02.13.10     IM - Aug 21, 2020   - Enhancement:
//                                      - Added caseBBStabilityPrescription in lieu of forceCaseBBBCStabilityFlag and alwaysStableCaseBBBCFlag to give more options for case BB/BC MT stability (issue #32)
// 02.13.11     IM - Aug 22, 2020   - Enhancement:
//                                      - Removed several stored options (e.g., m_OrbitalAngularVelocity, m_StarToRocheLobeRadiusRatio, etc.) to recompute them on an as-needed basis
//                                      - Removed some inf values in detailed outputs
//                                      - Slight speed-ups where feasible
//                                      - Shift various calculations to only be performed when needed, at printing, and give consistent values there (e.g., OmegaBreak, which was never updated previously)
//                                      - Remove a number of internal variables
//                                      - Declare functions constant where feasible
//                                      - Remove options to calculate Zetas and Lambdas at every timestep; variables that only appear in detailed outputs should not be computed at every timestep in a standard run
//                                      - Update documentation
//                                      - Remove postCEE binding energy (meaningless and wasn't re-computed, anyway)
// 02.13.12     IM - Aug 23, 2020   - Enhancement:
//                                      - More cleaning, removed some of the unnecessary prime quantities like m_SemiMajorAxisPrime, m_EccentricityPrime, etc.
//                                      - Thermal timescales are now correctly computed after the CE phase
//                                      - Detailed output passes a set of self-consistency checks (issue #288)
// 02.13.13     JR - Aug 23, 2020   - Defect repairs:
//                                      - Fixed debugging and logging macros in LogMacros.h
// 02.13.14     IM - Aug 29, 2020   - Defect repairs:
//                                      - Address issue #306 by removing detailed printing of merged binaries
//                                      - Address issue #70 by stopping evolution if the binary is touching
//                                      - Check for merged binaries rather than just touching binaries in Evaluate
//                                      - Minor cleaning (e.g., removed unnecessary CheckMassTransfer, which just repeated the work of CalculateMassTransfer but with a confusing name)
// 02.13.15     IM - Aug 30, 2020   - Defect repairs:
//                                      - Fixed issue #347: CalculateMassTransferOrbit was not correctly accounting for the MT_THERMALLY_LIMITED_VARIATION::RADIUS_TO_ROCHELOBE option
//                                      - Assorted very minor cleaning, including comments
// 02.14.00     IM - Aug 30, 2020   - Enhancement:
//                                      - Recreate RLOF printing (resolve issue #212)
// 02.14.01     ML - Sep 05, 2020   - Code cleanup:
//                                      - Issue #354 - Combine HYDROGEN_RICH and HYDROGEN_POOR supernova output variables into a single boolean variable IS_HYDROGEN_POOR 
// 02.15.00     JR - Sep 09, 2020   - Enhancements and related code cleanup:
//                                      - implemented "DETAILED_OUTPUT" folder inside "COMPAS_Output" container for SSE output
//                                      - SSE Parameters files moved to "DETAILED_OUTPUT" folder (they are analogous to BSE_Detailed_Output files)
//                                      - implemented SSE Switch Log and BSE Switch Log files (record written at the time of stellar type switch - see documentation)
//                                      - implemented SSE Supernova log file - see documentation (issue #253)
//                                      - added TIMESCALE_MS as a valid property in BaseStar::StellarPropertyValue().  The TIMESCALE_MS value in the SSE_Parameters file was being printed as "ERROR!" and nobody noticed :-)  It now prints correctly.
// 02.15.01     RS - Sep 10, 2020   - Enhancement
//                                       - added profiling option to keep track of repeated pow() calls
// 02.15.02     IM - Sep 11, 2020   - Defect repair
//                                       - changed ultra-stripped HeHG and HeGB stars to immediately check for supernovae before collapsing into WDs; this resolves issue #367
// 02.15.03     RTW - Sep 11, 2020   - Code cleanup:
//                                      - Set all references to kick "velocity" to magnitude. This is more correct, and will help distinguish from system and component vector velocities later
// 02.15.04     JR - Sep 11, 2020   - Enhancement
//                                       - refactored profiling code
//                                          - profiling code can now be #defined away for production build
//                                          - added options (via #defines) to profiling code: counts only (no CPU spinning), and print calling function name
//                                       - removed profiling program option
// 02.15.05     JR - Sep 12, 2020   - Code cleanup
//                                       - removed superfluous (and broken) #define guard around profiling.cpp
//                                       - minor change to profiling output (moved header and trailer to better place)
// 02.15.06     IM - Sep 12, 2020   - Defect repair
//                                       - Changed BaseBinaryStar::ResolveSupernova to account only for mass lost by the exploding binary during the SN when correcting the orbit
//                                       - Delayed supernova of ultra-stripped stars so that the orbit is adjusted in response to mass transfer first, before the SN happens
// 02.15.07     RTW - Sep 13, 2020   - Enhancement:
//                                      - Issue #12 - Move enhancement STROOPWAFEL from Legacy COMPAS to new COMPAS
//                                      - Issue #18 - double check STROOPWAFEL works in newCOMPAS
//                                      - Issue #154 - Test compatibility of CompasHPC and BSE_Grid.txt
//                                      - Added in combined functionality of Stroopwafel and pythonSubmit, with support for HPC runs
// 02.15.08     IM - Sep 14, 2020   - Defect repair:
//                                      - Issue #375 Error in Hurley remnant mass calculation
// 02.15.09     RTW - Oct 1, 2020   - Code cleanup:
//                                      - Rewrote ResolveSupernova to match Pfahl, Rappaport, Podsiadlowski 2002, and to allow for vector addition of system and component velocities
//                                      - Changed meaning of Supernova_State (see Docs)
//                                      - PostSN parameters have been removed
//                                      - SN phi has been redefined
// 02.15.10     IM - Oct 3, 2020    - Code cleanup:
//                                      - Removed some unnecessary internal variables and functions (m_TotalMass, m_TotalMassPrev, m_ReducedMass, m_ReducedMassPrev, m_TotalAngularMomentumPrev, CalculateAngularMomentumPrev(), EvaluateBinaryPreamble(),...
//                                      - Cleaned up some unclear comments
//                                      - ResolveCoreCollapseSN() no longer takes the Fryer engine as an argument (Fryer is just one of many possible prescriptions)
// 02.15.11     IM - Oct 3, 2020    - Defect repair and code cleanup:
//                                      - Fixed a number of defects in single stellar evolution (Github issues #381, 382, 383, 384, 385)
//                                      - The Fryer SN engine (delayed vs rapid) is no longer passed around, but read in directly in CalculateRemnantMassByFryer2012()
// 02.15.12     IM - Oct 5, 2020    - Enhancement
//                                      - Added timestep-multiplier option to adjust SSE and BSE timesteps relative to default
//                                      - Added eccentricity printing to RLOF logging
//                                      - Adjusted pythonSubmitDefault.py to include PESSIMISTIC CHE
//                                      - Updated documentation
// 02.15.13     JR - Oct 8, 2020    - Defect repair:
//                                      - Added checks for maximum time and timesteps to SSE code- issue #394
// 02.15.14     IM - Oct 8, 2020    - Defect repair:
//                                      - Added checks for dividing by zero when calculating fractional change in radius
// 02.15.15     IM - Oct 8, 2020    - Defect repair:
//                                      - Added safeguards for R<R_core in radius perturbation for small-envelope stars, complete addressing issue #394
// 02.15.16     RTW - Oct 14, 2020  - Code cleanup
//                                      - Changed separation to semiMajorAxis in RLOF and BeBinary properties
// 02.15.17     IM - Oct 16, 2020   - Defect repair and code cleanup:
//                                      - Issue 236 fixed: SN printing correctly enabled for all SNe
//                                      - Minor code cleaning: Cleaned up EvaluateSupernovae(), removed unnecessary m_Merged variable
// 02.15.18     RTW - Oct 22, 2020  - Code cleanup
//                                      - Removed redundant 'default' extension from files in the "defaults/" folder, and fixed references in the documentation.
//                                      - Added in '0' buffers to the Wall Times output to match the HH:MM:SS format
// 02.15.19     IM - Oct 23, 2020   - Enhancements
//                                      - Continue evolving DCOs until merger if EvolvePulsars is on (Issue #167)
//                                      - Removed m_SecondaryTooSmallForDCO (Issue #337)
// 02.15.20     RTW - Nov 03, 2020  - Code cleanup
//                                      - Removed unnecessary supernova phi rotation - it was added to agree with Simon's original definition, and to allow for seeds to reproduce the same SN final orbit. 
//                                      -   Removing it means seeds won't reproduce the same systems before and after, but populations are unaffected.
// 02.16.00     JR - Nov 03, 2020   - Enhancements
//                                      - Implemented new grid file functionality (see discussion in issue #412); updated docs - see docs (doc v2.3 has new documentation)
//
//                                      - Added all options to printing functionality: all options can now be selected for printing, 
//                                        either in the default log record specifications, or at runtime via the logfile-definitions option
//
//                                      - 'CHE_Option' header string changed to 'CHE_Mode'.  A few typos fixed in header strings.
//
//                                      - Added options
//                                          - initial-mass                          initial mass for single star (SSE)
//                                          - initial-mass-1                        initial mass for primary (BSE)
//                                          - initial-mass-2                        initial mass for secondary (BSE)
//                                          - semi-major-axis, a                    initial semi-major axis (BSE)
//                                          - orbital-period                        initial orbital period – only used if ‘semi-major-axis’ not specified
//                                          - eccentricity, e                       initial eccentricity (BSE)
//                                          - mode                                  mode of evolution: SSE or BSE (default is BSE)
//                                          - number-of-systems                     number of systems (single stars/binary stars) to evolve
//                                          - kick-magnitude-random                 kick magnitude random number for the star (SSE): used to draw the kick magnitude
//                                          - kick-magnitude                        the (drawn) kick magnitude for the star (SSE)
//                                          - kick-magnitude-random-1               kick magnitude random number for the primary star (BSE): used to draw the kick magnitude
//                                          - kick-magnitude-1                      the (drawn) kick magnitude for the primary star (BSE)
//                                          - kick-theta-1                          the angle between the orbital plane and the ’z’ axis of the supernova vector for the primary star (BSE)
//                                          - kick-phi-1                            the angle between ’x’ and ’y’, both in the orbital plane of the supernova vector, for the primary star (BSE)
//                                          - kick-mean-anomaly-1                   the mean anomaly at the instant of the supernova for the primary star (BSE)
//                                          - kick-magnitude-random-2               kick magnitude random number for the secondary star (BSE): used to draw the kick magnitude
//                                          - kick-magnitude-2                      the (drawn) kick magnitude for the secondary star (BSE)
//                                          - kick-theta-2                          the angle between the orbital plane and the ’z’ axis of the supernova vector for the secondary star (BSE)
//                                          - kick-phi-2                            the angle between ’x’ and ’y’, both in the orbital plane of the supernova vector, for the secondary star (BSE)
//                                          - kick-mean-anomaly-2                   the mean anomaly at the instant of the supernova for the secondary star (BSE)
//                                          - muller-mandel-kick-multiplier-BH      scaling prefactor for BH kicks when using 'MULLERMANDEL'
//                                          - muller-mandel-kick-multiplier-NS      scaling prefactor for NS kicks when using 'MULLERMANDEL'
//                                          - switchlog                             replaces ‘BSEswitchLog’ and ‘SSEswitchLog’
//                                          - logfile-rlof-parameters               replaces ‘logfile-BSE-rlof-parameters’
//                                          - logfile-common-envelopes              replaces ‘logfile-BSE-common-envelopes’
//                                          - logfile-detailed-output               replaces ‘logfile-BSE-detailed-output’, and now also used for SSE
//                                          - logfile-double-compact-objects		replaces ‘logfile-BSE-double-compact-objects’
//                                          - logfile-pulsar-evolution              replaces ‘logfile-BSE-pulsar-evolution’
//                                          - logfile-supernovae                    replaces ‘logfile-BSE-supernovae’ and ‘logfile-SSE-supernova’
//                                          - logfile-switch-log                    replaces ‘logfile-BSE-switch-log’ and ‘logfile-SSE-switch-log’
//                                          - logfile-system-parameters             replaces ‘logfile-BSE-system-parameters’
//
//                                      - Removed options
//                                          - number-of-binaries                    replaced by ‘number-of-systems’ for both SSE and BSE
//                                          - single-star-min                       replaced by ‘initial-mass’ and ‘number-of-stars’
//                                          - single-star-max                       replaced by ‘initial-mass’ and ‘number-of-stars’
//                                          - single-star-mass-steps                replaced by ‘initial-mass’ and ‘number-of-stars’
//                                          - BSEswitchLog                          replaced by ‘switchlog’
//                                          - SSEswitchLog                          replaced by ‘switchlog’
//                                          - logfile-BSE-rlof-parameters           replaced by ‘logfile-rlof-parameters’
//                                          - logfile-BSE-common-envelopes          replaced by ‘logfile-common-envelopes’
//                                          - logfile-BSE-detailed-output           replaced by ‘logfile-detailed-output’
//                                          - logfile-BSE-double-compact-objects    replaced by ‘logfile-double-compact-objects’
//                                          - logfile-BSE-pulsar-evolution          replaced by ‘logfile-pulsar-evolution’
//                                          - logfile-BSE-supernovae                replaced by ‘logfile-supernovae’
//                                          - logfile-SSE-supernova                 replaced by ‘logfile-supernovae’
//                                          - logfile-BSE-switch-log                replaced by ‘logfile-switch-log’
//                                          - logfile-SSE-switch-log                replaced by ‘logfile-switch-log’
//                                          - logfile-BSE-system-parameters         replaced by ‘logfile-system-parameters’
//
//                                      - Overloaded Options – these options are context-aware and are used for both SSE and BSE:
//                                          - number-of-systems                     specifies the number of systems (single stars/binary stars) to evolve
//                                          - detailed-output                       switches detailed output on/off for SSE or BSE
//                                          - switchlog                             enables the switch log for SSE or BSE
//                                          - logfile-detailed-ouput                defines filename for SSE or BSE detailed output file
//                                          - logfile-supernovae                    defines filename for SSE or BSE supernovae file
//                                          - logfile-switch-log                    defines filename for SSE or BSE switch log file
// 02.16.01     JR - Nov 04, 2020   - Enhancement
//                                      - changed switchlog implementation so that a single switchlog file is created per run
//                                        (see Issue #387 - note: single '--switch-log' option (shared SSE/BSE) implemented in v02.16.00)
// 02.16.02     IM - Nov 05, 2020   - Enhancements, Defect repairs
//                                      - Updated MT stability criteria for HeMS stars (Issue #425) to use MS zeta value
//                                      - Corrected baryon number for HeWD to match Hurley prescription (Issue #416)
//                                      - Corrected calculation of core mass after 2nd dredge-up (Issue #419)
//                                      - Corrected calculation of minimum radius on CHeB (Issue #420)
// 02.16.03     JR - Nov 08, 2020   - Defect repairs, Enhancements
//                                      - Issue #308
//                                          - added constant for minimum initial mass, maximum initial mass, minim metallicity and maximum metallicity to constants.h
//                                          - added checks to options code (specifically Options::OptionValues::CheckAndSetOptions()) to check option values for
//                                            initial mass and metallicity against constraints in constants.h
//                                      - Issue #342
//                                          - replaced header string suffixes '_1', '_2', '_SN', and '_CP' with '(1)', '(2)', '(SN)', and '(CP)' respectively
//                                          - now header strings ending in '(1)' indicate the value is for Star_1, '(2) for Star_2, '(SN)' for the supernova, and '(CP)' the companion
//                                      - Issue #351
//                                          - moved flags RECYCLED_NS and RLOF_ONTO_NS from SN_EVENT enum - now flags in BinaryConstiuentStar class
//                                          - removed RUNAWAY flag from SN_EVENT enum - removed entirely from code (not required)
//                                      - Issue #362
//                                          - changed header strings for RZAMS (radius at ZAMS) to 'Radius@ZAMS' - now consistent with MZAMS (mass at ZAMS - 'Mass@ZAMS')
//                                      - Issue #363
//                                          - made header strings for Lambdas uniform (all now start with 'Lambda_')
//                                      - Issue #409
//                                          - removed SN_THETA and SN_PHI from default SSE_SUPERNOVAE_REC (don't apply to SSE)
//                                      - Fixed defect that caused semi-major axis to be drawn from distribution rather than calculated from supplied orbital period
//                                        (moved check and calculation from options.cpp to BaseBinaryStar.cpp)
// 02.17.00     JR - Nov 10, 2020   - Enhancement, defect repairs, code cleanup
//                                      - Added SSE System Parameters file
//                                          - records initial parameters and result (final stellar type) 
//                                          - useful when detailed output is not required
//                                      - Fix for Issue #439
//                                      - Fixed typo in LogfileSwitchLog() in Options.h - only affected situation where user specified switchlog filename (overriding default filename)
//                                      - Removed m_LBVfactor variable from BaseBinaryStar - never used in BSE code
//                                      - Removed m_LBVfactor variable from BaseStar - use OPTIONS->LuminousBlueVariableFactor()
//                                      - Removed m_WolfRayetFactor variable from BaseBinaryStar - never used in BSE code
//                                      - Removed m_LBVfactor variable from BaseStar - use OPTIONS->WolfRayetFactor()
// 02.17.01     RTW - Nov 10, 2020  - Enhancement:
//                                      - Added in Schneider 2020 remnant mass prescriptions (standard and alternative)
//                                      - Added parameter MassTransferDonorHistory, as required for above prescription, which tracks the MT donor type (from which the MT Case can be established)
// 02.17.02     RTW - Nov 13, 2020  - Enhancement:
//                                      - Cleaned up the demo plotting routine so that the plot produced is the plot we use in the methods paper
// 02.17.03     JR - Nov 13, 2020   - Enhancements, code cleanup
//                                      - Added metallicity-distribution option: available distributions are ZSOLAR and LOGUNIFORM (see documentation)
//                                          - Added metallicity-min and metallicity-max options (for metallicity-distribution option)
//                                          - Metallicity is sampled if not explicitly specified via the --metallicity option - this was existing functionality, but
//                                            no distribution was implemented: sampling always returned ZSOLAR.  This change adds the LOGUNIFORM distribution, and 'formalises' the ZSOLAR 'distribution'.
//                                      - Added MASS to default SSE_SYSTEM_PARAMETERS_REC
//                                      - Removed AIS code
//                                      - Removed variable 'alpha' from BinaryCEDetails struct - use OPTIONS->CommonEnvelopeAlpha()
//                                          - Removed BINARY_PROPERTY::COMMON_ENVELOPE_ALPHA - use PROGRAM_OPTION::COMMON_ENVELOPE_ALPHA
//                                      - Issue #443: removed eccentricity distribution options FIXED, IMPORTANCE & THERMALISE (THERMALISE = THERMAL, which remains) 
// 02.17.04     JR - Nov 14, 2020   - Defect repairs
//                                      - Added CalculateRadiusOnPhase() and CalculateLuminosityOnPhase() to class BH (increases DNS yield)
//                                      - Added metallicity to sampling conditions in BaseBinaryStar constructor (should have been done when LOGUNIFORM metallicity distribution added)
// 02.17.05     TW - Nov 16, 2020   - Defect repairs
//                                      - Issue #444
//                                          - Fixed typo in synchronisation timescale
// 02.17.06     RTW - Nov 17, 2020  - Bug fix:
//                                      - Fixed Schneider remnant mass inversion from logRemnantMass^10 to 10^logRemnantMass, added some comments in the same section
// 02.17.07     TW - Nov 17, 2020   - Enhancements, code cleanup
//                                      - Issue #431
//                                          - Added option to change LBV wind prescription: choices are NONE, HURLEY_ADD, HURLEY and BELCYZNSKI
//                                      - Replaced numbers with constants for luminosity and temperature limits in mass loss
//                                      - Consolidated checks of luminosity for NJ winds within function
//                                      - NOTE: the above makes sure luminosity is checked before applying NJ winds for MS stars, this was not previously the case but I think it should be
// 02.17.08     JR - Nov 19, 2020   - Enhancements, code cleanup
//                                      - Added orbital-period-distribution option (see note in Options.cpp re orbital period option)
//                                      - Added mass-ratio option
//                                      - Updated default pythonSubmit to reflect new options, plus some previous omissions (by me...)
//                                      - Minor typo/formatting changes throughout
//                                      - Updated docs for new options, plus some typos/fixes/previous omissions
// 02.17.09     RTW - Nov 20, 2020  - Bug fix:
//                                      - Removed corner case for MT_hist=8 stars in the Schneider prescription (these should be considered Ultra-stripped)
// 02.17.10     RTW - Nov 25, 2020  - Enhancement:
//                                      - Cleaned up Schneider remnant mass function (now uses PPOW), and set the HeCore mass as an upper limit to the remnant mass
// 02.17.11     LVS - Nov 27, 2020  - Enhancements:
//                                      - Added option to vary all winds with OverallWindMassLossMultiplier
// 02.17.12     TW - Dec 9, 2020    - Enhancement, code cleanup, bug fix
//                                      - Issue #463
//                                          - Changed variable names from dml, dms etc. to rate_XX where XX is the mass loss recipe
//                                          - No longer overwrite variables with next mass loss recipe for clarity
//                                      - Added a new option to check the photon tiring limit during mass loss (default false for now)
//                                      - Added a new class variable to track the dominant mass loss rate at each timestep
// 02.17.13     JR - Dec 11, 2020   - Defect repair
//                                      - uncomment initialisations of mass transfer critical mass ratios in Options.cpp (erroneously commented in v02.16.00)
// 02.17.14     TW - Dec 16, 2020   - Bug fix
//                                      - fix behaviour at fLBV=0 (had been including other winds but should just ignore them)
// 02.17.15     JR - Dec 17, 2020   - Code and architecture cleanup
//                                      - Architecture changes:
//                                          - Added Remnants class    - inherits from HeGB class
//                                          - Added WhiteDwarfs class - inherits from Remnants class; most of the WD code moved from HeWD, COWD and ONeWD to WhiteDwarfs class
//                                          - Changed HeWD class      - inherits from WhiteDwarfs class (COWD still inherits from HeWD; ONeWD from COWD)
//                                          - Change NS class         - inherits from Remnants class; code added/moved as necessary
//                                          - Change BH class         - inherits from Remnants class; code added/moved as necessary
//                                          - Change MR class         - inherits from Remnants class; code added/moved as necessary
//                                      - Code cleanup:
//                                          - added "const" to many functions (mostly SSE code) that dont modify class variables ("this") (still much to do, but this is a start)
//                                          - added "virtual" to GiantBranch::CalculateCoreMassAtBAGB() and BaseStar::CalculateTemperatureAtPhaseEnd()
//                                              - will have no impact given where they are called, but the keyword should be there (in case of future changes)
//                                          - changed hard-coded header suffixes from _1 -> (1), _2 -> (2)
//                                      - Added call to main() to seed random number generator with seed = 0 before options are processed (and user specified seed is know).  Ensures repeatability.
//                                      - Changed "timestep below minimum" warnings in Star.cpp to be displayed only if --enable-warnings is specified
// 02.17.16     JR - Dec 17, 2020   - Code cleanup
//                                      - Removed "virtual" from GiantBranch::CalculateCoreMassAtBAGB() (incorrectly added in v02.17.15 - I was right the first time)
//                                      - Removed "const" from Remnants::ResolveMassLoss() (inadvertently added in v02.17.15)
//                                      - Removed declarations of variables m_ReducedMass, m_ReducedMassPrev, m_TotalMass, and m_TotalMassPrevfrom BaseBinaryStar.h (cleanup begun in v02.15.10 - these declarations were missed)
// 02.17.17     RTW - Dec 17, 2020  - Code cleanup
//                                      - Removed MassTransferCase related variables in favor of MassTransferDonorHist
// 02.17.18     JR - Dec 18, 2020   - Defect repair
//                                      - Typo in options code for option --switch-log: "switchlog" was incorrectly used instead of "switch-log"
// 02.17.19     LVS - Dec 19, 2020  - Enhancements:
//                                      - Added option to vary winds of cool stars (with T < VINK_MASS_LOSS_MINIMUM_TEMP) via a CoolWindMassLossMultiplier
// 02.18.00     JR - Jan 08, 2021   - Enhancement:
//                                      - Added support for HDF5 logfiles (see notes at top of log.h)
//                                      - Added 'logfile-type' option; allowed values are HDF5, CSV, TSV, TXT; default is HDF5
//                                      - Added 'hdf5-chunk-size' option - specifies the HDF5 chunk size (number of dataset entries)
//                                      - Added 'hdf5-buffer-size' option - specifies the HDF5 IO buffer size (number of chunks)
//                                      - Removed 'logfile-delimiter' option - delimiter now set by logfile type (--logfile-type option described above)
//                                      - Changed header strings containing '/' character: '/' replaced by '|' (header strings become dataset names in HDF5 files, and '/' is a path delimiter...)
// 02.18.01     SS - Jan 11, 2021   - Defect repair
//                                      - Added check if binary is bound when evolving unbound binaries
// 02.18.02     JR - Jan 12, 2021   - Defect repair:
//                                      - Changed "hdf5_chunk_size = 5000" to "hdf5_chunk_size = 100000" in default pythonSubmit (inadvertently left at 5000 after some tests...)
// 02.18.03     SS - Jan 19, 2021   - Enhancement:
// 									    - Added check for neutron star mass against maximum neutron star mass. 
//									      If a neutron star exceeds this mass it should collapse to a black hole. 
//                                        This can be relevant for neutron stars accreting, e.g. during common envelope evolution
// 02.18.04     IM - Jan 28, 2021   - Enhancement:
//                                      - NS to BH collapse preserves mass (see discussion in #514)
//                                      - Fixed comment typo
// 02.18.05     JR - Jan 29, 2021   - Defect repair:
//                                      - Honour '--evolve-unbound-systems' option when specified in a grid file (see issue #519)
//                                      - Honour '--evolve-pulsars' option when specified in a grid file (same as issue #519)
//                                      - Added "maximum-evolution-time", "maximum-number-timestep-iterations", and "timestep-multiplier" to m_GridLineExcluded vector in Options.h (previous oversight)
// 02.18.06     SS - Feb 1, 2021    - Defect repair:
//                                      - Make COMPAS use --neutrino-mass-loss-BH-formation options (resolves issue #453)
// 02.18.07     JR - Feb 18, 2021   - Enhancement:
//                                      - Added 'rotational-frequency' option so users can specify initial rotational frequency of SSE star
//                                      - Added 'rotational-frequency-1' and 'rotational-frequency-2' options so users can specify initial rotational frequency of both BSE stars
//                                      - Changed units of rotational frequencies written to logfiles (omega, omega_break, omega_ZAMS) from rotations per year to Hz
//                                      - Changed program option header strings containing '_1' and '_2' to '(1)' and '(2)' for consistency
// 02.18.08     JR - Feb 26, 2021   - Defect repairs:
//                                      - Remove stray diagnostic print from BaseStar constructor in BaseStar.cpp
//                                      - Fix for issue #530 - some commandline options ignored when a grid file is used
//                                          - the issue here was case-sensitive vs case-insensitive matches (asking BOOST to do case-insensitive matches for option names doesn't propagate to all matches BOOST does...)
//                                          - the options affected were all options that have mixed-case names:
//
//                                              - case-BB-stability-prescription
//                                              - kick-magnitude-sigma-CCSN-BH
//                                              - kick-magnitude-sigma-CCSN-NS
//                                              - kick-magnitude-sigma-ECSN
//                                              - kick-magnitude-sigma-USSN
//                                              - mass-transfer-thermal-limit-C
//                                              - muller-mandel-kick-multiplier-BH
//                                              - muller-mandel-kick-multiplier-NS
//                                              - neutrino-mass-loss-BH-formation
//                                              - neutrino-mass-loss-BH-formation-value
//                                              - PISN-lower-limit
//                                              - PISN-upper-limit
//                                              - PPI-lower-limit
//                                              - PPI-upper-limit
// 02.18.09     ML - Mar 22, 2021   - Defect repair:
//                                      - Correct polynomial evaluation of Nanjing lambda's for EAGB and TPAGB stellar types.
// 02.18.10     LVS - Apr 06, 2021   - Enhancement:
//                                      - Added PPISN prescription option - Farmer 2019
// 02.19.00     JR - Apr 20, 2021   - Enhancements and Defect Repairs:
//                                      - Enhancements:
//                                          - Added option to enable users to add program options values to BSE/SSE system parameters files
//                                              - option is '--add-options-to-sysparms', allowed values are {ALWAYS, GRID, NEVER}.  See docs for details.
//                                          - Included "Run_Details" file in HDF5 output file if logfile type = HDF5.  The text Run_Details file still exists
//                                            so users can still easily look at the contents of the Run_Details file - this enhancements adds a copy of the
//                                            Run_Details file to the HDF5 output file.
//
//                                      - Defect Repairs:
//                                          - fixed a few previously unnoticed typos in PROGRAM_OPTION map in constamts.h, and in Options::OptionValue() function.
//                                            Fairly benign since they had't been noticed, but needed to be fixed.
//
//                                      Modified h5copy.py (in postProcessing/Folders/H5/PythonScripts) so that groups (COMPAS files) will not be copied
//                                      if the group exists in the destination file but has a different number of datasets (columns) from the group in
//                                      the source file.
//
//                                      Also provided h5view.py - an HDF5 file viewer for COMPAS HDF5 files (in postProcessing/Folders/H5/PythonScripts).  See
//                                      documentation as top of source file for details.
// 02.19.01     JR - Apr 30, 2021   - Enhancements and Defect Repairs:
//                                      - Enhancements:
//                                          - changed chunk size for HDF5 files to HDF5_MINIMUM_CHUNK_SIZE for Run_Details group in COMPAS_Output and for detailed output files.
//                                              - Run_Details is a small file, and detailed output files are generally a few thousand records rather than hundreds of thousands, 
//                                                so a smaller chunk size wastes less space and doesn't impact performance significantly
//
//                                      - Defect Repairs:
//                                          - fixed issue #548 - HDF5 detailed output files not created when random-seed specified in a grid file
//                                          - fixed defect where records in HDF5 output files would be duplicated if the number of systems exceeded the HDF5 chunk size
//                                            being used (the default chunk size is 100000 - that might explain why this problem hasn't been reported)
//
//                                      Modified h5view.py (in postProcessing/Folders/H5/PythonScripts) to handle detailed output files
// 02.19.02     LVS - May 04, 2021   - Defect Repairs:
//                                      - Avoid possibility of exceeding total mass in Farmer PPISN prescription
// 02.19.03     TW - May 18, 2021    - Enhancement:
//                                      - Change default LBV wind prescription to HURLEY_ADD
// 02.19.04     JR - May 24, 2021    - Defect Repair:
//                                      - Fixed incrementing of random seed and binary id when grid file contains sets/ranges
//
//                                      Modified h5view.py (in postProcessing/Folders/H5/PythonScripts) to print number of unique seeds (where relevant) in summary output
// 02.20.00     IM - June 14, 2021  - Enhancement:
//                                      - Port defaults from preProcessing/pythonSubmit.py to options.cpp
//                                      - Minor fixes (e.g., documentation)
// 02.20.01     JR - June 21, 2021  - Defect repair:
//                                      - Fix for issue #585: add formatted value and delimiter to logrecord string in Log.h (defect introduced in v02.18.00; only affected SSE_Supernovae logfile)
// 02.20.02     JR - July 26, 2021  - Defect repair:
//                                      - Add HDF5 support to logging code for SSE/BSE switch log files.  Support for HDF5 switch files was inadvertently not added when HDF5 file support as added in v02.18.00 for all standard log files.  Switch log files are 'special' (they have extra columns, not part of the 'standard' log file functionality), and that was missed.
//                                      - Also removed '-lsz' from Makefile and Makefile.docker - library not required
// 02.21.00     JR - July 28, 2021  - Enhancement and Defect Repairs:
//                                      - Added code to copy any grid file and/or logfile-definitions file specified to output container.
//                                      - Copying a large grid file could take time, and take up much space, so added new program option '--store-input-files' which is TRUE by default.  If FALSE, neither the grid file (if specified) nor the logfile-definitions file (if specified) will be copied to the output container (if TRUE, both will be copied (if specified)).
//                                      - Fixed issue #600: changed pythonSubmit.py to treat fully-qualified grid filenames and fully-qualified logfile-definitions filenames correctly (i.e. don't add CWD if the filename is already fully-qualified).
//                                      - Fixed issue #601: changed pythonSubmit.py to put all boolean parameters on the commandline, with "True" or "False" value.
// 02.21.01     RTW - Aug 21, 2021  - Defect Repair:
//                                      - PrintRLOFProperties now gets called immediately before and after the call to EvaluateBinary so that the changes reflect only BSE changes.
//                                      - The function call has also been tidied up to take an argument specifying whether the call was made before or after the MT took place.
// 02.22.00     JR - Aug 26, 2021   - Enhancement:
//                                      - Added functionality to allow users to select a range of lines from the grid file (if specified) to process.  Added program options --grid-start-line and --grid-lines-to-process - see documentation for details.
// 02.22.01     JR - Sep 11, 2021   - Defect repair:
//                                      - Fix for issue #615: defaults for calculated/drawn program options now calculated after random seed is set for the system being evolved.
// 02.22.02     IM - Oct 4, 2021    - Defecr repair:
//                                      - Removed unnecessary IsPrimary() / BecomePrimary() functionality, fixed incorrect MassTransferTrackerHistory (see issue #605)
// 02.22.03     IM - Oct 4, 2022    - Defect repair:
//                                      - Corrected Eddington mass accretion limits, issue #612 (very minor change for WDs and NSs, factor of a few increase for BHs)
// 02.23.00 FSB/JR - Oct 11, 2021   - Enhancement:
//                                      - updated kelvin-helmholtz (thermal) timescale calculation with more accurate pre-factor and updated documentation.
//                                      - rationalised parameters of, and calls to, CalculateThermalTimescale()
// 02.23.01     JR - Oct 11, 2021   - Code cleanup:
//                                      - Typo fixed in version for changes made on October 11, 2021
//                                      - Changed KROUPA_POWER to SALPETER_POWER in utils:SampleInitialMass(); Removed KROUPA_POWER from constants.h
//                                      - Removed p_Id parameter from SSE/BSE switchlog functions - leftover from debugging
//                                      - Added CHEMICALLY_HOMOGENEOUS_MAIN_SEQUENCE property to SSE_SYSTEM_PARAMETERS_REC and BSE_SYSTEM_PARAMETERS_REC (both stars)
//                                      - Tidied up some parameters etc. to better comply with COMPAS coding guidelines
//                                      - Typo fixed in preProcessing/COMPAS_Output_Definitions.txt
// 02.24.00     JR - Oct 12, 2021   - Minor enhancements/optimisations:
//                                      - Added BaseStar::CalculateThermalMassAcceptanceRate() as a first-pass to address issue #595 - can be changed/expanded as required
//                                      - Changed BaseBinaryStar::CalculateTimeToCoalescence() to use Mandel 2021 https://iopscience.iop.org/article/10.3847/2515-5172/ac2d35, eq 5 to address issue #538
// 02.24.01     RTW - Oct 13, 2021  - Enhancements:
//                                      - Added units uniformly to the --help input descriptions
//                                      - Removed the BeBinary- and RLOF-specific random seeds (which were attributes of the events and were printed with e.g <MT) and replaced with system random seed
//                                      - In CE output, changed MASS_2_FINAL (which was sort of a wrapper for core mass) for MASS_2_POST_COMMON_ENVELOPE
//                                      - Removed SN kick angles from SystemParameters output (they are duplicated in SN output) and changed true_anomaly to mean_anomaly in BSE SN output
//                                      - Cosmetic typo fixes and added consistency, in the Event_Counter parameters and some function definitions
//                                      - Added *.eps, *.png to gitignore
// 02.24.02     JR - Oct 13, 2021   - Minor fixes:
//                                      - Fixed a few typos in header strings
//                                      - Changed true_anomaly to mean_anomaly in SSE SN output
// 02.25.00     JR - Oct 30, 2021   - Enhancements and minor fixes:
//                                      - Added ability for users to annotate log files via new program options '--notes-hdrs' and '--notes'.  See docs for details. 
//                                      - Added a shorthand notation for vector program options (e.g. annotations, log-classes, debug-classes).  See docs for details.
//                                      - Added '--notes-hdrs' and '--notes' to pythonSubmit.py (default = None for both)
//                                      - Added HDF5 support to Log::GetLogStandardRecord() (return value) and Log::LogStandardRecord() (input parameter).  This only matters
//                                        to SSE Supernovae file - for delayed writes.  The original implementation may have resulted in minor discrepanicies in SSE Supernovae
//                                        log records, (because of when the values were sampled (i.e. mid-timestep, or end of timestep)), which would only have been evident if
//                                        HDF5 files were compared to e.g. CSV files for the same binary - CSV, TSV, and TXT files had values sampled mid-timestep, HDF5 files 
//                                        at end of timestep).
//                                      - Added Log::Write() and Log::Put() for HDF5 files (better implementation - worked around in original implementation)
//                                      - Added additional checks for bad string -> number conversions throughout (for stoi(), stod(), etc.)
//                                      - Performance enhancement to BaseBinaryStar::CalculateTimeToCoalescence() (return early if e = 0.0)
//                                      - Fixed a few typos in comments
// 02.25.01     IM - Nov 1, 2021    -  Enhancements:
//                                      - Introduced common-envelope-allow-radiative-envelope-survive and common-envelope-allow-immediate-rlof-post-ce-survive options
//                                      - Addresses issue # 637
// 02.25.02     JR - Nov 1 , 2021    - Minor fixes:
//                                      - reinstated "_n" suffix for BSE detailed filenames (inadvertently removed in v02.25.00)
//                                      - updated pythonSubmit files:
//                                          preProcessing/pythonSubmit.py
//                                          examples/methods_paper_plots/detailed_evolution/pythonSubmitDemo.py
//                                          examples/methods_paper_plots/chirpmass_distribution/pythonSubmit.py
//                                          examples/methods_paper_plots/fig_5_HR_diagram/pythonSubmit.py
//                                          examples/methods_paper_plots/fig_6_max_R/pythonSubmit.py
//                                          examples/methods_paper_plots/fig_8_initial_core_final_mass_relations/pythonSubmitDefaults.py
//                                          examples/methods_paper_plots/fig_8_initial_core_final_mass_relations/pythonSubmitFryerRapid.py
//                                          examples/methods_paper_plots/fig_8_initial_core_final_mass_relations/pythonSubmitMandelMueller.py
// 02.25.03     JR - Nov 1 , 2021    - Minor fixes:
//                                      - fixed typo in Options.cpp for option --common-envelope-allow-immediate-RLOF-post-CE-survive (was typed common-envelope-allow-immediate-RLOF-post-CE_survive)
//                                      - fixed typo in Options.cpp for option --common-envelope-allow-radiative-envelope-survive (was typed common-envelope-allow-radiative-envelope-surive)
//                                        (neither of these caused problems because Boost matches only as many characters as necessary to determine the option name - would have if the names were not unique up to the typos)
// 02.25.04     IM - Nov 4, 2021     - Minor fixes
//                                      - More surive->survive typo fixes in python files to address issue #660
//                                      - Documentation edits to reflect new options common-envelope-allow-radiative-envelope-survive and common-envelope-allow-immediate-rlof-post-ce-survive options
// 02.25.05     IM - Nov 4, 2021     - Defect repair:
//                                      - Changed GiantBranch::CalculateRemnantMassByMullerMandel() to ensure that the remnant mass is no greater than the He core mass
// 02.25.06     IM - Nov 7, 2021     - Enhancements:
//                                      - Clarified program option documentation
//                                      - Removed unused CUSTOM semi-major axis initial distribution
//                                      - Removed unused STARTRACK zeta prescription
// 02.25.07     IM - Nov 12, 2021    - Defect repair:
//                                      - Changed EAGB::CalculateLuminosityOnPhase() and EAGB::CalculateLuminosityAtPhaseEnd() to use the helium core mass rather than the CO core mass (see Eq. in second paragraph of section 5.4 of Hurley+, 2000); this fixes a downward step in luminosity and radius on transition to EAGB
// 02.25.08     JR - Nov 15, 2021    - Defect repair:
//                                      - Fixed error introduced in v02.25.00: Added HDF5 support to GetLogStandardRecord().
//                                        Defect introduced was omission of code for HDF5 file support if a specified property is supplied to GetLogStandardRecord(), causing a boost::bad_get error.
//                                        The defect only affected HDF5 SSE_Supernovae files.  This fix adds the omitted code.
//                                      - Changed Options::PrintOptionHelp() to print help (-h/--h) to stdout instead of stderr.
// 02.25.09     IM - Nov 16, 2021    - Defect repair:
//                                      -Revert EAGB treatment to 02.25.06 until a proper fix is introduced
// 02.25.10     JR - Nov 19, 2021    - Defect repairs:
//                                      - clamp timestep returned in BaseStar::CalculateTimestep() to NUCLEAR_MINIMUM_TIMESTEP
//                                      - change NUCLEAR_MINIMUM_TIMESTEP to 1 year (from 100 years) in constants.h
// 02.26.00     IM - Nov 30, 2021    - Defect repairs:
//                                      - only decrease effective initial mass for HG and HeHG stars on mass loss when this decrease would not drive an unphysical decrease in the core mass
//                                      - change mass comparisons (e.g., mass vs. He flash mass threshold) to compare effective initial mass rather than current mass
//                                      - minor code and comment cleanup
// 02.26.01     IM - Dec 5, 2021     - Defect repair, Code cleanup:
//                                      - Removed redundant function ResolveRemnantAfterEnvelopeLoss (ResolveEnvelopeLoss is sufficient)
//                                      - Cleaned / updated ResolveEnvelopeLoss
//                                      - Fixed issue with masses and types of remnants formed from stripped HG stars
// 02.26.02     RTW - Dec 17, 2021   - Defect repair, Code cleanup:
//                                      - Changed all occurrences of PPOW(base, 1.0/3.0) with std::cbrt, as the former could not handle negative bases
//                                      - Changed all occurrences of sqrt with std::sqrt for consistency with the above change
// 02.26.03     IM - Jan 10, 2022    - Defect repair, code cleanup:
//                                      - Cleaned up treatment of HG donors having CONVECTIVE envelopes in LEGACY; fixed an issues with CEs from HG donors introduced in 02.25.01 
// 02.27.00     ML - Jan 12, 2022    - Enhancements:
//                                      - Add enhanced Nanjing lambda option that continuously extrapolates beyond radial range
//                                      - Add Nanjing lambda option to switch between calculation using rejuvenated mass and true birth mass
//                                      - Add Nanjing lambda mass and metallicity interpolation options
//                                      - No change in default behaviour
// 02.27.01     IM - Feb 3, 2022     - Defect repair:
//                                      - Fixed condition for envelope type when using ENVELOPE_STATE_PRESCRIPTION::FIXED_TEMPERATURE (previously, almost all envelopes were incorrecctly declared radiative)
// 02.27.02     IM - Feb 3, 2022     - Defect repair:
//                                      - Fixed mass change on forced envelope loss in response to issue # 743
// 02.27.03     JR - Feb 8, 2022     - Defect repair:
//                                      - Fix for issue # 745 - logfile definition records not updated correctly when using logfile-definitions file (see issue for details)
// 02.27.04     RTW - Feb 15, 2022   - Defect repair:
//                                      - Fix for issue # 761 - USSNe not occurring. See issue for details.
// 02.27.05     IRS - Feb 17, 2022   - Enhancements:
//                                      - Add function HasOnlyOneOf, which returns true if a binary has only one component in the list of stellar types passed, and false if neither or both are in the list
//                                      - Add function IsHMXRBinary, which returns true if HasOnlyOneOf(Black hole, Neutron star) and the companion radius is > 80% of the Roche Lobe radius
//                                      - Add flag --hmxr-binaries, which tells COMPAS to store binaries in BSE_RLOF output file if IsHMXRBinary
//                                      - Add columns for pre- and post-timestep ratio of stars to Roche Lobe radius to BSE_RLOF output file (addressing issue #746)
//                                      - Changed variables named rocheLobeTracker, roche_lobe_tracker etc. to starToRocheLobeRadiusRatio, star_to_roche_lobe_radius_ratio, etc. for clarity
// 02.27.06     SS - Apr 5, 2022     -  Defect repair:
//                                      - Fixed StarTrack PPISN prescription, previously it was doing the same thing as the COMPAS PPISN prescription.
// 02.27.07     RTW - Apr 5, 2022    - Defect repair:
//                                      - Fix for issue # 773 - ONeWD not forming due to incorrect mass comparison in TPAGB. 
// 02.27.08     RTW - Apr 12, 2022   - Defect repair:
//                                      - Fix for issue # 783 - Some mergers involving a massive star were not logged properly in BSE_RLOF, whenever a jump in radius due to changing stellar type within ResolveMassChanges was much greater than the separation.
// 02.27.09     VK - Apr 25, 2022    - Minor Enhancement:
//                                      - Converted constant: MULLERMANDEL_SIGMAKICK into an option: --muller-mandel-sigma-kick
// 02.28.00     Lvs - May 11, 2022   - Enhancements:
//                                      - Introduced new remnant mass prescription: Fryer+ 2022, adding new options --fryer-22-fmix and --fryer-22-mcrit
// 02.29.00     RTW - May 5, 2022    - Enhancement:
//                                      - Fix for issue # 596 - New option to allow for H rich ECSN (defaults to false). This removes non-interacting ECSN progenitors from contributing to the single pulsar population.
// 02.30.00     RTW - May 8, 2022    - Enhancement
//                                      - Added MACLEOD_LINEAR specific angular momentum gamma loss prescription for stable mass transfer
// 02.31.00     IM - May 14, 2022    - Enhancement
//                                      - Added option retain-core-mass-during-caseA-mass-transfer to preserve a larger donor core mass following case A MT, set equal to the expected core mass of a newly formed HG star with mass equal to that of the donor, scaled by the fraction of its MS lifetime
//                                      - Code and comment cleaning
// 02.31.01     RTW - May 16, 2022   - Defect repair:
//                                      - Fixed help string for H rich ECSN option implemented in v2.29.99
// 02.31.02     JR - May 18, 2022    - Defect repairs:
//                                      - Fixed STAR_PROPERTY_LABEL entries in contsants.h for INITIAL_STELLAR_TYPE and INITIAL_STELLAR_TYPE_NAME - both missing the prefix "INITIAL_".
//                                        Only caused a problem if a user wanted to add either of those to the logfile-definitions file - but since they are in the system parameters files (SSE and BSE)
//                                        by default encountering the problem would probably be unlikely.
//                                      - Fixed error identifier in Log::UpdateAllLogfileRecordSpecs() - was (incorrectly) ERROR::UNKNOWN_BINARY_PROPERTY, now (correctly) ERROR::UNKNOWN_STELLAR_PROPERTY 
// 02.31.03     RTW - May 20, 2022   - Defect repair:
//                                      - Fixed MS+MS unstable MT not getting flagged as a CEE
// 02.31.04     RTW - June 10, 2022  - Enhancements
//                                      - Fixed MT_TRACKER values to be clearer and complementary to each other
//                                      - Updated the relevant section in the detailed plotter that uses MT_TRACKER values
//                                      - Removed end states from detailed plotter (Merger, DCO, Unbound) so that they don't over compress the rest
// 02.31.05     RTW - July 25, 2022  - Defect repair:
//                                      - Renamed option '--allow-H-rich-ECSN' to 'allow-non-stripped-ECSN'
//                                      - Fixed check for non-interacting ECSN progenitors to consider MT history instead of H-richness
// 02.31.06     RTW - Aug 2, 2022    - Enhancement:
//                                      - Added stellar merger to default BSE_RLOF output
// 02.31.07     IM - August 1, 2022  - Defect repair:
//                                      - Print to DetailedOutput after merger, addresses https://github.com/TeamCOMPAS/COMPAS/issues/825
//                                      - Ensures no ONeWDs are formed with masses above Chandrasekhar mass
//                                      - Minor comment tweaks and a bit of defensive programming
// 02.31.08     RTW - Aug 3, 2022    - Enhancement:
//                                      - Added Accretion Induced Collapse (AIC) of ONeWD as another type of SN
// 02.31.09     RTW - Aug 9, 2022    - Enhancement:
//                                      - Max evolution time and max number of timesteps now read in from gridline as well as commandline
// 02.31.10     RTW - Aug 12, 2022   - Enhancement:
//                                      - Added option to set the Temperature boundary between convective/radiative giant envelopes
// 02.32.00     JR - Aug 27, 2022    - Enhancement & minor cleanup:
//                                      - Add 'record type' functionality to all standard log files
//                                      - Add/rationalise calls to PrintDetailedOutput() for binary systems
//                                          - remove m_PrintExtraDetailedOutput variable (and associated code) from BaseBinaryStar class
//                                      - Add new option for each standard log file to allow specification of which record types to print
//                                          - see e.g. '--logfile-detailed-output-record-types'
//                                      - Online documentation updated for record types and new options
//                                      - Detailed ploter changed to work with record type column (thanks RTW)
//                                      - Added new section to online documentation: 'What's new'
//                                          - documented record types changes in this new section
//                                      - Minor cleanup:
//                                          - minor formatting and typo fixes (src + docs)
//                                          - removed IncrementOmega() function from the BaseStar and Star classes (anti-patterm and no longer used - if it ever was)
//                                          - tidied up description of MainSequence::UpdateMinimumCoreMass()
// 02.33.00     RTW - Aug 13, 2022   - Enhancement:
//                                      - Added critical mass ratios from Claeys+ 2014 for determining if MT is unstable
//                                      - Cleaned up stability check functions in BaseBinaryStar.cpp for clarity, and to allow for critical mass ratios to be checked correctly
// 02.33.01     RTW - Sep 26, 2022   - Defect repair:
//                                      - Fixed interpolation of MACLEOD_LINEAR gamma for specific angular momentum. Previously interpolated on the gamma value, now interpolates in orbital separation
// 02.33.02      IM - Nov 27, 2022   - Defect repair:
//                                      - Fixed ignored value of input radius when computing the thermal timescale, relevant if using Roche lobe radius instead (issue #853)
//                                      - Cleaned code and comments around the use of MT_THERMALLY_LIMITED_VARIATION::RADIUS_TO_ROCHELOBE vs. C_FACTOR (issue #850)
// 02.34.00      IM - Nov 28, 2022   - Enhancement:
//                                      - Adding framework for Hirai & Mandel 2-stage common envelope formalism
//                                          (placeholders for now -- will have identical results to default version)
//                                      - Placed Dewi CE prescription on parity with others
// 02.34.01     RTW - Nov 30, 2022   - Defect repair:
//                                      - Fixed Time<MT and Time>MT calls in BSE_RLOF. Previously, they were identical. Now, Time<MT correctly reflects the previous time.
// 02.34.02     JR - Nov 30, 2022    - Defect repair:
//                                      - Fixed problem with no content in switchlog files (issue #870 - introduced in v2.32.00).
//                                      - Changed conditional statement in HG::ResolveEnvelopeLoss() and FGB::ResolveEnvelopeLoss() to be consistent with other stellar types ('>' -> '>=').
// 02.34.03     NRS - Jan 9, 2023    - Defect repair:
//                                      - Fixed units for post-CEE semi-major axis in CEE logs (issue #876).
// 02.34.04     RTW - Jan 31, 2023   - Enhancement:
//                                      - Added SN orbit inclination angle to BSE_SUPERNOVAE output
// 02.34.05     JR - Jan 29, 2023    - Code cleanup:
//                                      - Addressed issue #888 - replaced class variables m_LogMetallicityXi, m_LogMetallicitySigma, and m_LogMetallicityRho in BaseStar with getter functions.
// 02.34.06     IM - Feb 1, 2023     - Bug fixes:
//                                      - Re-enabled ResolveMassLoss() for Remnants so that Mdot is correctly reset
//                                      - Set Mdot to 0 in BaseBinaryStar::CalculateWindsMassLoss() when winds are turned off while the binary is in mass trensfer
//                                      - Removed Dutch winds for Remnants
//                                      - Fixed typo in NS::CalculateLuminosityOnPhase_Static()
//                                      - Minor code cleaning
// 02.35.00     RTW - Dec 8, 2022    - Enhancement:
//                                      - Added critical mass ratios from Ge+ 2020 for determining if MT is unstable
// 02.35.01     RTW - Feb 12, 2022   - Enhancement:
//                                      - Added post-SN orbital inclination vector to the output-able BINARY_PROPERTIES (not included in output, by default). 
// 02.35.02     JR - Feb 19, 2023    - Minor change and defect repair:
//                                      - Changed units of ROCHE_LOBE_RADIUS_1 and ROCHE_LOBE_RADIUS_2 from orbital separation to RSOL
//                                      - Changed header string for ROCHE_LOBE_RADIUS_1 from "RocheLobe(1)|a" to "RocheLobe(1)" - ditto for ROCHE_LOBE_RADIUS_2
//                                      - removed STAR_TO_ROCHE_LOBE_RADIUS_RATIO_1 ("Radius(1)|RL")and STAR_TO_ROCHE_LOBE_RADIUS_RATIO_2 ("Radius(2)|RL") from
//                                        the default output for BSE_DETAILED_OUTPUT_REC (can be calculated from other values in the default output)
//                                      - changed plot_detailed_evolution.py to accommodate the removal of STAR_TO_ROCHE_LOBE_RADIUS_RATIO_1 and 
//                                        STAR_TO_ROCHE_LOBE_RADIUS_RATIO_2 from the default output
//                                      - changed online documentation to reflect:
//                                           (a) removal of STAR_TO_ROCHE_LOBE_RADIUS_RATIO_1 and STAR_TO_ROCHE_LOBE_RADIUS_RATIO_2 from the default output
//                                           (b) change of header strings for ROCHE_LOBE_RADIUS_1 and ROCHE_LOBE_RADIUS_2 (units already (accidentally...) correct)
//                                      - fixed minor defect in call to m_Accretor->CalculateMassAcceptanceRate() in BaseBinaryStar::CalculateMassTransfer()
//                                        (only affected runs with mass-transfer-thermal-limit-accretor = RADIUS_TO_ROCHELOBE)
// 02.35.03     LvS - Feb 27, 2023   - Enhancement:
//                                      - Added mass accretion prescription during CE following model 2 from van Son + 2020
// 02.36.00     JR - Mar 15, 2023    - Enhancement, minor defect repairs:
//                                      - Addressed issue #797 - implemented functionality to create YAML file.  Two new options (--create-YAML-file and --YAML-template).  See documentation for details.
//                                      - Modified runSubmit.py to work with new yaml file format (i.e. all options could be commented...)
//                                      - Minor defect repairs in options code
//                                      - Minor fixes to online documentation; also clarified make arguments
// 02.36.01     JR - Mar 20, 2023    - Documentation:
//                                      - Updated documentation for YAML files.
//                                      - Modified YAML template to include notice regarding commented lines in default YAML file.
// 02.37.00     NR,RTW - Mar 26, 2023 - Enhancement:
//                                      - Added functionality for WDs to accrete in different regimes. 
//                                          - This applies to each WD subtype individually, though there is some overlap between COWDs and ONeWDs.
//                                          - Also involves tracking the WD shell mass, to account for shell burning that later increases the WD mass.
//                                          - Includes possible instability, and merger if the donor is a giant, as well as new SN types, 
//                                          - AIC (accretion induced collapse), SNIA (Type Ia), and HeSD (Helium shell detonation). 
//                                      - Tangential but related changes:
//                                          - Cleaned up the call to EddingtonCriticalRate, puttting it in BaseStar along with the optional prefactor.
//                                          - Moved NS radius and luminosity calls into NS.h from elsewhere in the code.
// 02.37.01     JR - Mar 27, 2023    - Defect repair:
//                                      - Updated changelog.h and whats-new.rst to correctly reflect changes to the code and version numbers after a bad fix for merge conflicts
//                                      - Changed "DD" to "HeSD" as appropriate
//                                      - A couple of code-cleanups
// 02.37.02     JR - Mar 27, 2023    - Defect repair:
//                                      - Changed yaml.h to include <algorithm> and <chrono> - not including them causes docker build to fail.
// 02.37.03     IM - Apr 8, 2023     - Defect repair:
//                                      - Resolved issue #855 by using Mass0 rather than Mass to determine ages and timescales
// 02.38.01     IM - Apr 16, 2023    - Enhancement:
//                                      - Added option to eject the convective envelope by pulsations (ExpelConvectiveEnvelopeAboveLuminosityThreshold)
//                                          if log10(m_Luminosity/m_Mass) exceeds LuminosityToMassThreshold
// 02.38.02     NR - Apr 20, 2023    - Defect repair:
//                                      - Added missing const in WD files which was generating warnings when compiling.
//                                   - Enhancement:
//                                      - Added QCRIT_PRESCRIPTION::HURLEY_HJELLMING_WEBBINK based on Hurley+ 2002 and its corresponding documentation.
// 02.38.03     IM - Apr 20, 2023    - Enhancement:
//                                      - Updated defaults following #957
// 02.38.04     IM - Apr 20, 2023    - Enhancement:
//                                      - Included Picker et al. (2023, in prep.) fits for the convective envelope mass in the TWO_STAGE common envelope treatment
// 02.38.05     YS - May 10, 2023    - Updates and changes to NS.cpp:
//                                      - Added NS::ChooseTimeStep(). Detailed time step description and reasoning can be found in NS.cpp
//                                      - Added output options (not default): PULSAR_BIRTH_PERIOD and PULSAR_BIRTH_SPIN_DOWN_RATE, which output the birth spin period and period derivative of a pulsar
//                                      - Updated codes on pulsar evolution, solving the problem of pulsars not evolving properly. This is written in cgs. 
//                                      - Added NS::SpinDownIsolatedPulsar(), describes single pulsar spinning down with magnetic braking. 
//                                          This is later used in NS::UpdateMagneticFieldAndSpin()
//                                      - m_PulsarDetails.spinDownRate was described as Pdot (s s^-1), when it is in fact f-dot(rad s^-2). This is now corrected. 
//                                      - In BSE_Pulsar_Evolution file, the pulsar parameters at birth were not recorded. 
//                                          Pulsar was also evolved an additional time step here with unspecified size.
//                                          Fix to this problem is done by setting the PULSAR_RECORD_TYPE:
//                                           (a) if record_type = 1 (DEFAULT), these are the initial values of the pulsar set at birth
//                                           (b) if record_type = 3 (POST_BINARY_TIMESTEP), these describe normal pulsar evolution
//                                      - Another caveat:
//                                         pulsar recycling mechanisms are not yet fully implemented, so COMPAS cannot produce MSPs for the time being; more updates to come.
// 02.38.06     JR - Jun 04, 2023    - Defect repair:
//                                      - Fixed "hides overloaded virtual function" warnings.
//                                      - Added "-Woverloaded-virtual" to compiler flags to enable warnings for g++ on linux systems.
// 02.38.07     JR - Jun 04, 2023    - Defect repair:
//                                      - Fix for issue #958 - evolving unbound systems that contain two compact objects.  Also added BINARY_PROPERTY::UNBOUND 
//                                        to BSE Detailed Output file default record.
//                                      - Changed makefile to be POSIX compliant for .o suffix rule.  No need to change docker Makefile - it is already POSIX compliant.
//                                         - since GNU Make 4.3 a warning is issued for suffix rules that include prerequisites - in our case the .o rule on line 125:
//                                           "Makefile:125: warning: ignoring prerequisites on suffix rule definition" - and in future releases it will become an error.

const std::string VERSION_STRING = "02.38.07";


# endif // __changelog_h__
