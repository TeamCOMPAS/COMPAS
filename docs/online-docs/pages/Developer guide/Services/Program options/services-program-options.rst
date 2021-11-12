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
