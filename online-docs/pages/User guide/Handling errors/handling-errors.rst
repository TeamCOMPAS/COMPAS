Handling COMPAS errors
======================

COMPAS output files record detailed information about the evolution of single stars (SSE mode) and binary stars (BSE mode).
In addition to the attributes of the systems evolved, COMPAS output files also record the disposition (evolution status) of
the evolved systems, and any errors that may have occurred during evolution.


Evolution status
----------------

The evolution status is the final disposition of the system evolved. Users should check the final disposition of any system
before using the results recorded in any of the COMPAS output files for that system.

The evolution status is recorded, by default, in the System Parameters file (SSE or BSE) - follow link below for a full list
of possible values for evolution status. If the evolution status recorded for a system indicates that evolution was stopped
because an error occurred, the results recorded in the COMPAS output files for that system are not complete, and not trustable.


Error value
-----------

If the evolution status recorded for a system indicates that evolution was stopped because an error occurred, the error can be
checked by looking at the error number recorded in the System Parameters file - follow link below for a full list of possible
values the error.


Floating-point errors
---------------------

Floating-point errors (divide-by-zero, invalid-argument, overflow, and underflow) will be indicated by setting the error in the
COMPAS output file to an appropriate value (see link below for COMPAS error values).

When floating-point error checking is ON (``--fp-error-mode ON``), if a floating-point error is encountered during evolution of
a system, evolution of that system will be stopped, and the evoltion status and error values set appropriately.

When floating-point error checking is OFF (``--fp-error-mode OFF``), floating-point errors that occur during the evolution of a 
system will not cause the evolution of that system to stop. However, if a floating-point error does occur during evolution while
floating-point error checking is off, the error value in the COMPAS output files will still indicate that a floating-point error
occurred. If that happens, users should use the results of the evolution of an affected system with caution.

When floating-point error checking is in debug mode (``--fp-error-mode DEBUG``), any floating-point errors encountered will cause
COMPAS to display a stack trace and halt execution. The stack trace displayed will allow developers to determine the location of
the code (to function level) that caused the floatig-point error and implement a code repair.

Note that ``--fp-error-mode OFF`` is the default.

.. toctree::
   :maxdepth: 1

   evolution-status-table
   error-table
