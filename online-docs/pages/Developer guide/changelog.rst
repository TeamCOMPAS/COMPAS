Change log
==========

The COMPAS code includes a change log header file - ``changelog.h``.  The change log records changes to the COMPAS ``C++`` code, as well
as the COMPAS version string.

Following is a fragment of the change log::

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

Recorded in each entry of the  changelog is:

    - the new COMPAS version number
    - initials of the developer making the changes
    - date of the changes
    - a brief description of the changes - including github issue number where appropriate


COMPAS version number
---------------------

Currently the COMPAS version number is set manually whenever changes are made to the code.  A planned enhancement is to have the v ersion number
increment automaticall whenever a github pull request is mmerged.

The version number is formatted as: 

    `major`\.\ `minor`\ .\ `defect`

where each of the components is a 2-digit integer.

The `major` number will increment whenever significant new functionality is added to COMPAS. |br|
The `minor` number will increment whenever minor new functionality is added to COMPAS. |br|
The `defect` number will increment whenever defect repairs are made to COMPAS.

`major` and `minor` (and `significant`) are somewhat subjective terms. Some guidance on what constitures a `major` release vs a `minor` release is
given below. `defect` releases are releases that repair known defects.

The COMPAS version string is recorded towards the end of the change log file, and should be incremented whenver a change is made to the code::

    const std::string VERSION_STRING = "02.22.00";


Major releases
--------------

A major release typically includes substantial changes to the way the application functions, and introduces key improvements in functionality. A major
release might include:

    - refactoring the code base `(as we did when we moved from COMPAS v1 to COMPAS v2)`
    - significantly improved performance
    - significant changes to the user interface `(the new grid file format probably should have been a major release)`
    - significant new features
    - removing deprecated features
    - integration with other applications
    
Major releases typically occur somewhat infrequently.


Minor releases
--------------

Minor releases introduce new features to the application. Minor releases are smaller than major 
releases - they can be regarded as `edits` to the current version of the application. Minor releases are not a total overhaul - they enhance and 
improve existing functionality. A minor release might include:

    - limited new features and functionality
    - small updates to existing features

Minor releases typically occur fairly frequently.
