h5view.py
=========

This program displays summary, header, and content information for specified COMPAS ``HDF5`` file(s). It's fairly rudimetary - the 
``HDF5`` package provides ``h5dump`` and ``h5ls`` which have far more options than this program - but this program is somewhat 
COMPAS aware.


h5view usage
------------

::

    h5view.py [-h] [-f FILENAME_FILTER] [-r [RECURSION_DEPTH]] [-S] [-H]
                   [-C [CONTENTS]] [-s] [-x EXCLUDE_GROUP [EXCLUDE_GROUP ...]]
                   [-V SEED_LIST [SEED_LIST ...]]
                   input [input ...]

    HDF5 file content viewer.

    positional arguments:
      input
        input directory and/or file name(s)

    optional arguments:
      -h, --help
        show this help message and exit
      -f FILENAME_FILTER, --filter FILENAME_FILTER
        input filename filter (default = *)
      -r [RECURSION_DEPTH], --recursive [RECURSION_DEPTH]
        recursion depth (default is no recursion)
      -S, --summary
        display summary output for HDF5 file (default is not to displat summary)
      -H, --headers
        display file headers for HDF5 file (default is not to display headers)
      -C [CONTENTS], --contents [CONTENTS]
        display file contents for HDF5 file: argument is number of entries (+ve from top, -ve
        from bottom) (default is not to display contents)
      -s, --stop-on-error
        stop all copying if an error occurs (default is skip to next file and continue)
      -x EXCLUDE_GROUP [EXCLUDE_GROUP ...], --exclude EXCLUDE_GROUP [EXCLUDE_GROUP ...]
        list of input groups to be excluded (default is all groups will be copied)
      -V SEED_LIST [SEED_LIST ...], --seeds SEED_LIST [SEED_LIST ...]
        list of seeds to be printed (for content printing) (default is print all seeds)


Example
-------

Typing::

    python3 h5view.py compas-output-file.h5
    
will result in summary output of the ``HDF5`` file `compas-output-file.h5` that looks something like this::

    Summary of HDF5 file /d/compas/h5out.h5
    =======================================

    File size    : 2.1520 GB
    Last modified: 2021-07-26 16:25:12.928401

    COMPAS Filename              Columns   Entries   Unique SEEDs
    --------------------------   -------   -------   ------------
    Run_Details                      346        30
    BSE_Common_Envelopes              73    582485         476514
    BSE_Double_Compact_Objects        12      8725           8725
    BSE_RLOF                          34   2997481         536332
    BSE_Supernovae                    32    103000          87162
    BSE_Switch_Log                    13   7472234         956623
    BSE_System_Parameters             33   1050000        1050000


Other ``h5view.py`` options (listed above) display headers and file contents.



h5view functionality overview
-----------------------------

``h5view.py`` displays summary, header and content information for specified COMPAS ``HDF5`` file(s). If none of the command-line
options ``--summary [-S]``, ``--headers [-H]``, or ``--contents [-C]`` are specified, ``--summary [-s]`` is assumed. If any of 
``--summary [-S]``, ``--headers [-H]``, or ``--contents [-C]`` are specified, then only the option(s) specified are actioned.

Displaying summary information
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Summary information displays, for each COMPAS file in the ``HDF5`` file:
   - the name of the COMPAS file
   - the number of columns in the COMPAS file
   - the number of entries in the COMPAS file (actually, the maximum number of entries in any column in the COMPAS file)


Displaying header information
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Header information displays, for each COMPAS file in the ``HDF5`` file:
   - the name of each column in the COMPAS file
   - the number of entries in each column of the COMPAS file
   - the data type of each column of the COMPAS file
   - the units associated with each column of the COMPAS file
     (with the exception of the ``Run_Details`` file - there are no units associated with columns in the ``Run_Details`` file)


Displaying contents
~~~~~~~~~~~~~~~~~~~

Contents information displays, for each COMPAS file in the ``HDF5`` file:
   - a header showing the column names in the COMPAS file
   - a row for each entry in the COMPAS file, showing the column values for that row (comma delimited)

   The contents display can be limited in two ways:

      (a) The ``--contents [-C]`` option takes and optional argument: an integer number of rows to display. The argument to 
          ``--contents [-C]`` can be positive or negative: a positive value indicates that the number of rows specified by the 
          argument should be displayed from the start of the file; a negative value indicates that the number of rows specified
          by the (absolute value of the) argument should be displayed from the end of the file.  The +ve and -ve arguments to 
          the ``--contents [-C]`` option are akin the the Unix ``head`` and ``tail`` commands.

      (b) The ``--seeds [-V]`` option allows the user to specify a list of SEED values that should be printed. If the 
          ``--seeds [-V]`` option is specified, only rows containing the seeds specified by the user will be printed - and only 
          if they are in the entries printed if limited by the ``--contents [-C]`` argument  described in (a).

          Note that printing only seeds specified in a list of seeds could be slow - we effectively have to look through the 
          entire dataset looking for the seeds required.

