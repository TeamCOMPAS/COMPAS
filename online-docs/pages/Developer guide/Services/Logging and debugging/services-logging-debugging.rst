Logging & debugging service
===========================

A logging and debugging service is provided encapsulated in a singleton object (an instantiation of the ``Log`` class).


.. toctree::
   :maxdepth: 1

   ./Base-level logging/services-base-level-logging
   ./services-extended-logging
   ./Logging macros/services-logging-debugging-macros-logging
   ./services-logging-debugging-macros-debugging


:bolditalictext:`Historical note:`

The logging functionality was first implemented when the Single Star Evolution code was refactored, and the base-level of logging
was sufficient for the needs of the ``SSE`` code. Refactoring the Binary Star Evolution code highlighted the need for expanded logging
functionality. To provide for the logging needs of the ``BSE`` code, new functionality was added almost as a wrapper around the original,
base-level logging functionality. Some of the original base-level logging functionality has almost been rendered redundant by the new
functionality implemented for ``BSE`` code, but it remains (almost) in its entirety because it may still be useful in some circumstances.

When the base-level logging functionality was created, debugging functionality was also provided, as well as a set of macros to make
debugging and the issuing of warning messages easier. A set of logging macros was also provided to make logging easier. The debug
macros are still useful, and their use is encouraged (rather than inserting print statements using ``std::cout`` or ``std::cerr``).

When the ``BSE`` code was refactored, some rudimentary error handling functionality was also provided in the form of the 
:doc:`Errors service <../services-error-handling>` an attempt at making error handling easier. Some of the functionality provided by the
:doc:`Errors service <../services-error-handling>` supersedes the ``DBG_WARN*`` macros provided as part of the Log class, but the ``DBG_WARN*``
macros are still useful in some circumstances (and in fact are still used in various places in the code). The ``LOG*`` macros are somewhat less
useful, but remain in case the original base-level logging functionality (that which underlies the expanded logging functionality) is used in 
the future (as mentioned above, it could still be useful in some circumstances).

The expanded logging functionality introduces Standard Log Files - described in :doc:`services-extended-logging`.
