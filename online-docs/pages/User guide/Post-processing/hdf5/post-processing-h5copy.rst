h5copy.py
=========

This program copies ``COMPAS_Output.h5`` ``HDF5`` file(s) [but not ``Detailed_Ouput`` files] to a designated output ``HDF5`` file. 
If the output file is an existing ``HDF5`` file, the user can specify whether the existing content should be erased before copying 
begins, or whether the copied data should be appended to the existing data. If multiple files are given as input files, the 
resultant ``HDF5`` file is the concatenation of the input files.

Use this program either with:
```bash
python /path/to/h5copy.py [options]
```

or, if you have installed the COMPAS python utilities:
```bash
compas_h5copy [options]
```



Some nomenclature
-----------------

Data in ``HDF5`` files are arranged in ``groups`` and ``datasets``:

    A COMPAS output file (e.g. `BSE_System_Parameters`, `BSE_RLOF`, etc.) maps to an ``HDF5 group`` where the group name is
    the name of the COMPAS output file.

    A column in a COMPAS output file (e.g. `SEED`, `Mass(1)`, `Radius(2)`, etc.) maps to an ``HDF5 dataset``, where the
    dataset name is the column heading string.

    COMPAS column datatype strings are encoded in the ``HDF5 dataset`` meta-details (``dataset.dtype``).

    COMPAS column units strings are attached to ``HDF5 datasets`` as ``attributes``.


h5copy usage
------------

::

    h5copy.py [-h] [-b BUFFER_SIZE] [-c CHUNK_SIZE] [-e] [-f FILENAME_FILTER]
                   [-o OUTPUT_FILENAME] [-r [RECURSION_DEPTH]] [-s]
                   [-x EXCLUDE_GROUP [EXCLUDE_GROUP ...]]
                   input [input ...]

    HDF5 file copier.

    positional arguments:
      input
        input directory and/or file name(s)

    optional arguments:
      -h, --help
        show this help message and exit
      -b BUFFER_SIZE, --buffer-size BUFFER_SIZE
        IO buffer size (number of HDF5 chunks, default = 10)
      -c CHUNK_SIZE, --chunk-size CHUNK_SIZE
        HDF5 output file dataset chunk size (default = 100000)
      -e, --erase-output
        erase existing output file before copying input files (default = False)
      -f FILENAME_FILTER, --filter FILENAME_FILTER
        input filename filter (default = *)
      -o OUTPUT_FILENAME, --output OUTPUT_FILENAME
        output file name (default = h5out.h5)
      -r [RECURSION_DEPTH], --recursive [RECURSION_DEPTH]
        recursion depth (default is no recursion)
      -s, --stop-on-error
        stop all copying if an error occurs (default is skip to next file and continue)
      -x EXCLUDE_GROUP [EXCLUDE_GROUP ...], --exclude EXCLUDE_GROUP [EXCLUDE_GROUP ...]
        list of input groups to be excluded (default is all groups will be copied)


    Note: if the -x option is specified, it should be specified at the end of the options 
          list (i.e. the list of input files can't follow the -x option or the list of input 
          files will be subsumed by the list of groups to be excluded)


h5copy functionality overview
-----------------------------

Output file
~~~~~~~~~~~

- If the specified ``HDF5`` output file does not exist, it will be created, and data from the input file(s) will be appended 
  to the new output file.

- If the specified ``HDF5`` output file does exist, and the ``--erase-ouput [-e]`` command-line option is not specified, the 
  existing output file will be preserved and data from the input file(s) will be appended to the existing output file if the 
  output file was created with chunking enabled (see :ref:`HDF5-chunking` - only files that were created with chunking enabled 
  can be extended).

  Attempting to append data to an existing ``HDF5`` file that was not created with chunking enabled will result in the following 
  message being displayed for each dataset::

      Only chunked datasets can be resized

  If that happens, using ``h5copy.py`` to copy the output file to a new file will copy the existing data and create the new file 
  with chunking enabled - the newly created file can then be used as a base file to which other files can be appended.

- If the specified ``HDF5`` output file does exist, and the ``--erase-ouput [-e]`` command-line option is specified, the existing
  output file will be deleted, and new, empty, file created, and data from the input file(s) will be appended to the new output file.


##############
Appending data
##############

   (a) if a group in an input file already exists in the output file, then the group data (the datasets within the group) will only 
       be copied if the number of datasets in the input file group matches the number of datasets in the output file group. If there 
       is a mismatch a warning will be issued and the group will not be copied (but the file copy will continue, just as though the 
       group had been excluded) 
   
   (b) if a dataset in an input file does not exist in the output file, the dataset will be created, otherwise the data copied will 
       be appended to the existing dataset.


Input files
~~~~~~~~~~~

A list of input filenames and or directory names must be supplied. The list can be a single name. Each name in the list is processed 
in order. If a name is the name of a file with the file extension `.h5`, the contents of the file will be copied to the output file, 
otherwise it will be ignored. If a name is the name of a directory and the specified recursion level requires that the directory be 
processed (see command-line option ``--recursive [-r]``), the program will descend into the directory and process all files and 
directories there (directories will be processed depending upon the value of the ``--recursive [-r]`` option), otherwise it will be 
ignored.

The command-line option ``--recursive [-r]`` specifies whether recursion is enabled for directory processing, and if it is, to what 
depth:

    - If the ``--recursive [-r]`` option is not specified, recursion is not enabled and only files in the specified working directory 
      will be candidates for copying.

    - if ``--recursive [-r]`` is specified with no ``depth`` value, recursion is enabled and the depth is not limited - that is, all 
      files in the specified working directory, and all files in all directories below the specified working directory, will be 
      candidates for copying.
    
    - If ``--recursive [-r]`` is specified with a specified ``depth`` value, recursion is enabled and the depth is limited to the
      depth specified - that is, all files in the specified working directory, and all files in all directories `depth` levels below
      the specified working directory, will be candidates for copying.
    

#####################
Input filename filter
#####################

If the ``--filter [-f]`` command-line option is specified, the names of all candidate files will be checked against the specified filter,
and only files whose names match the filter will be copied.  The specified filter is a filename-only filter - the file's path (i.e. its
location) will not be matched to the filter. The specified filter should not include a file extension, but the program adds the extension
`.h5` to the specified filter - only files that have the file extension `.h5` will match the filter.

If ``--filter [-f]`` is not specified, the program uses a default filter value of `*`, then adds the `.h5` file extension - so all candidate
files with the `.h5` extension will be copied.


Excluding HDF5 groups
~~~~~~~~~~~~~~~~~~~~~

If the ``--exclude [-x]`` command-line option is specified, the specified list of groupnames will be excluded from data copied from all input
files.  If ``--exclude [-x]`` is not specified, all groups in all candidate files will be copied.


Erase output
~~~~~~~~~~~~

If the ``--erase-ouput [-e]`` command-line option is specified and the output file (specified or default) exists, it will be erased before 
copying begins.  The ``--erase-ouput [-e]`` command-line option is ignored if the output file does not exist.

If ``--erase-ouput [-e]`` is not specified and the output file (specified or default) exists, the existing content will be preserved and any
data copied to the file will be appended to the existing data.

