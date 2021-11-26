COMPAS output
=============

Summary and status information during the evolution of stars is written to stdout; how much is written depends upon the value of the 
``--quiet`` program option.

Detailed information is written to log files (described below). All COMPAS output files are created inside a container directory, specified 
by the ``--output-container`` program option.

If detailed output log files are created (see the ``--detailed-output`` program option), they will be created inside a containing directory 
named ``Detailed_Output`` within the COMPAS output container directory.

Also created in the COMPAS container directory is a file named ``Run_Details`` in which COMPAS records some details of the run (COMPAS 
version, start time, command line program option values etc.). Note that the option values recorded in the Run details file are the values
specified on the command line, not the values specified in a ``grid`` file (if used).


.. toctree::
   :maxdepth: 1

   standard-logfiles
   standard-logfiles-record-specifiers
   standard-logfiles-record-specification
   standard-logfiles-format
   standard-logfiles-annotations
   standard-logfiles-example-definitions-file
   standard-logfiles-default-record-specifications
