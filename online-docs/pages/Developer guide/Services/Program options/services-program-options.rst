Program options service
=======================

A Program Options service is provided, encapsulated in a singleton object (an instantiation of the ``Options`` class).

The ``Options`` class member variables are private, and public getter functions have been created for the program options currently
used in the code.

The Options service can be accessed by referring to the ``Options::Instance()`` object. For example, to retrieve the value of 
the ``--quiet`` program option, call the ``Options::Quiet()`` getter function::

    bool quiet = Options::Instance()→Quiet();

Since that could become unwieldy, there is a convenience macro to access the Options service. The macro just defines "OPTIONS" as
"Options::Instance()", so retrieving the value of the ``--quiet`` program option can be written as::

    bool quiet = OPTIONS→Quiet();

The Options service must be initialised before use. Initialise the Options service by calling the ``Options::Initialise()`` function::

    COMMANDLINE_STATUS programStatus = OPTIONS→Initialise(argc, argv);

(see ``constants.h`` for details of the ``COMMANDLINE_STATUS`` type)

See :doc:`../../../User guide/Program options/program-options-list-defaults` for a full list of available program options and their default valaues.


Adding new program options
--------------------------

To add a new program option, take the following steps (these are duplicated at the top of the ``Options.cpp`` source file):

1. Decide on a string for the option - this is the string the user will use on the commandline or in the grid file (e.g. "random-seed"). The convention is hyphenated lower case. Try to be consistent with existing option names.  For example, if you are adding a new prescription for something, make sure the option name ends with "-prescription", or if you are adding a new distribution for something, make sure the option name ends with "-distribution".

2. Decide on a class member variable name for the option (e.g. m_RandomSeed).

3. Decide on a default value for the option.

4. Add the class member variable to the PRIVATE area of the OptionValues class in ``Options.h``.

5. Add a getter for the class member variable to the PUBLIC area of the Options class in ``Options.h``.

   Decide if the getter should always retrieve the value specified on the commandline, or whether it should retrieve the grid line value if one was specified by the user, and only retrieve the commandline value if the user did not specify a grid line value - see the ``OPT_VALUE`` macro defined in ``Options.h``.

6. Add the class member variable initialisation (to the default value) to the ``Options::OptionValues::Initialise()`` function in ``Options.cpp``.

7. Add the option to the Options::AddOptions() function in Options.cpp. This is where we tell Boost about the option - the option string, the variable in which the user-specified value should be stored, the default value to use if the user does not specify a value, and the (text) description of the option.

   ``Options::AddOptions()`` is a bit of a beast - there's no easy, short way to specify the details of as many options as we have.  I have tried to make it semi-readable, but it is, and always will be, just long...  The best we can do is keep it neat so it doesn't become too hard to read.

   When adding options to ``Options::AddOptions()``, the convention I have used is that the option string (e.g. "random-seed") is predominantly in lower case - that's not strictly required, but I think it's easier for users to remember the option names if they don't have to remember a mixture of upper and lower case.  I have configured ``Boost`` to perform a case-insensitive match, so option strings can have mixed case if necessary (e.g. ``--muller-mandel-kick-multiplier-BH``).

   Where appropriate, use ``AllowedOptionValuesFormatted()`` to format the allowed values for an option - doing so ensures all allowed values are included, and will keep the format consistent.

8. Add any sanity checks: constraint/range/dependency checks etc. for the new option, and any affected existing options, to ``Options::OptionValues::CheckAndSetOptions()`` in ``Options.cpp``.  It is also here you can set any final values that, perhaps due to dependencies on options that had not yet been parsed, could not be set directly by ``Boost`` when the options were parsed (also see ``SetCalculatedOptionDefaults()``; viz. ``m_KickPhi1`` etc.).

9. If the option is a string option with multiple-choices - in that the user can select from a list of possible values recorded in ``typedefs.h`` or ``LogTypedefs.h`` in an ``ENUM CLASS`` and corresponding ``COMPASUnorderedMap`` labels map - then add the option to the function ``AllowedOptionValues()`` here so that we can easily extract the allowed values for that option.

10. Add the new option to one or more of the following vectors in ``Options.h``, as required:

        ``m_ShorthandAllowed``: options for which shorthand notation is allowed

        ``m_GridLineExcluded``: option strings that may not be specified on a grid line

        ``m_SSEOnly``         : option strings that apply to SSE only
        ``m_BSEOnly``         : option strings that apply to BSE only

        ``m_RangeExcluded``   : option strings for which a range may not be specified
        ``m_SetExcluded``     : option strings for which a set may not be specified

    Read the explanations for each of the vectors in ``Options.h`` to get a better idea of what they are for and where the new option should go.

11. Add the new option to the following structures in ``LogTypedefs.h`` (only required if the option is required to be available for printing in the logfiles):

       - enum class PROGRAM_OPTION
       - const COMPASUnorderedMap<PROGRAM_OPTION, std::string> PROGRAM_OPTION_LABEL
       - const std::map<PROGRAM_OPTION, PROPERTY_DETAILS> PROGRAM_OPTION_DETAIL

12. Add the new option to ``Options::OptionValue()`` - this enables selection of the option value for printing in the output (log) files.  Only required if the option is required to be available for printing in the logfiles.

13. Add the new option to the default YAML file template in ``yaml.h`` if required. While this is not mandatory, it does help to keep the YAML file ordered and consequently relatively easy to read and search for options.
