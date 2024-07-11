Error catalog (ErrorCatalog.h)
==============================

``ErrorCatalog.h`` is the COMPAS error catalogue.  The error catalogue defines symbolic names for all COMPAS errors, and
corresponding error strings for those errors.

To add a new error, add the symbolic name to the ``ERROR`` enum class, and the corresponding error string to the ``ERROR_CATALOG`` map.
 
The key to the ``ERROR_CATALOG`` map is the symbolic name in the ``ERROR`` enum class.  The map entry is a tuple containing the
``ERROR_SCOPE`` associated with the error (see below), and the error string.

The ``ERROR_SCOPE`` enum class provides a mechanism to allow developers to specify when, if at all, a particular error/warning
should be displayed by the ``SHOW_WARN*`` and ``SHOW_ERROR*`` macros (see :doc:`../../Developer guide/Services/services-error-handling`).

The values for ``ERROR_SCOPE`` and their meanings are:

::

      NEVER                : the error/warning should never be displayed
      ALWAYS               : the error/warning should always be displayed
      FIRST                : the error/warning should only be displayed the first time it is encountered
      FIRST_IN_FUNCTION    : the error/warning should only be displayed the first time it is encountered for a function
      FIRST_IN_STELLAR_TYPE: the error/warning should only be displayed the first time it is encountered for a stellar type (see enum class STELLAR_TYPE in typedefs.h)
      FIRST_IN_OBJECT_TYPE : the error/warning should only be displayed the first time it is encountered for an object type (see enum class OBJECT_TYPE in typedefs.h)
      FIRST_IN_OBJECT_ID   : the error/warning should only be displayed the first time it is encountered for an object id (each object is assigned a unique object id - e.g. a star or binary, each constituent star of a binary)

The ``THROW_ERROR*`` macros (see :doc:`../../Developer guide/Services/services-error-handling`) are not affected by ``ERROR_SCOPE``.

A convenience function for retrieving the error text is #defined here:

    #define ERR_MSG(x) std::get<1>(ERROR_CATALOG.at(x))

