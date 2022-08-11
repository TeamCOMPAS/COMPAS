Error handling service
======================

An error handling service is provided encapsulated in a singleton object (an instantiation of the ``Errors`` class).

The Errors service provides global error handling functionality. Following is a brief description of the Errors service
(full documentation coming soon...):

Errors are defined in the error catalog in ``constants.h`` (see ``ERROR_CATALOG``). |br|

Errors defined in the error catalog have a scope and message text. The scope is used to determine if and when an error should
be printed.

The current values for scope are:

    .. list-table::
       :widths: 30 70 
       :header-rows: 0
       :class: aligned-text

       * - **NEVER**
         - the error will not be printed.
       * - **ALWAYS**
         -  the error will always be printed.
       * - **FIRST**
         -  the error will be printed only on the first time it is encountered anywhere in the program in the current execution of COMPAS.
       * - **FIRST_IN_OBJECT_TYPE**
         - the error will be printed only on the first time it is encountered anywhere in objects of the same type (e.g. Binary Star objects) in the current execution of COMPAS.
       * - **FIRST_IN_STELLAR_TYPE**
         - the error will be printed only on the first time it is encountered anywhere in objects of the same stellar type (e.g. HeWD Star objects) in the current execution of COMPAS.
       * - **FIRST_IN_OBJECT_ID**
         - the error will be printed only on the first time it is encountered anywhere in an object instance in the current execution of COMPAS (useful for debugging).
       * - **FIRST_IN_FUNCTION**
         - the error will be printed only on the first time it is encountered anywhere in the same function of an object instance in the current execution of COMPAS (i.e. will print more than once if encountered in the same function name in different objects; useful for debugging).

The Errors service provides methods to print both warnings and errors – essentially the same thing, but warning messages are prepended with 
':boldtext:`WARNING:` ', whereas error messages are prepended with ':boldtext:`ERROR:` '.

Errors and warnings are printed by using the macros defined in ``ErrorsMacros.h``. They are:

Error macros
------------

::
    
    SHOW_ERROR(error_number)                         // prints the error message associated with error
                                                     // number (from the error catalog) prepended by
                                                     // 'ERROR: '

    SHOW_ERROR(error_number, error_string)           // prints the error message associated with error
                                                     // number (from the error catalog) prepended by
                                                     // 'ERROR: ', and appends 'error_string'

    SHOW_ERROR_IF(cond, error_number)                // if 'cond' is TRUE, prints the error message
                                                     // associated with error number (from the error
                                                     // catalog) prepended by 'ERROR: '

    SHOW_ERROR_IF(cond, error_number, error_string)  // if 'cond' is TRUE, prints the error message
                                                     // associated with error number (from the error
                                                     // catalog) prepended by 'ERROR: ', and appends 
                                                     // 'error_string'


Warning macros
--------------

The ``WARNING`` macros function in the same way as the ``ERROR`` macros, with the exception that instead of prepending the
message with ':boldtext:`ERROR:` ', the ``WARNING`` macros prepend the message with ':boldtext:`WARNING:` '.

The ``WARNING`` macros are:

::

    SHOW_WARN(error_number)
    SHOW_WARN(error_number, error_string)
    SHOW_WARN_IF(cond, error_number)
    SHOW_WARN_IF(cond, error_number, error_string)

|br|
Error and warning messages always contain:

- The object id of the calling object.
- The object type of the calling object.
- The stellar type of the calling object (will be ”NONE” if the calling object is not a star-type object).
- The function name of the calling function.

Any object that uses the Errors service (i.e. any of the ``ERROR`` and ``WARNING`` macros) must expose the following functions:

::

    OBJECT_ID ObjectId()** const { return m ObjectId; }
    OBJECT_TYPE ObjectType()** const { return m ObjectType; }
    STELLAR_TYPE StellarType()** const { return m StellarType; }

These functions are called by the ``ERROR`` and ``WARNING`` macros. If any of the functions are not applicable to the object, 
then they must return ':italictext:`type`::NONE' (all objects should implement ObjectId() correctly).

The filename to which error records are written when ``Log::Start()`` parameter ``p_ErrorsToFile`` is TRUE is declared in 
``constants.h`` – see the enum class ``LOGFILE`` and associated descriptor map ``LOGFILE_DESCRIPTOR``. Currently the name is 'Error_Log'.
