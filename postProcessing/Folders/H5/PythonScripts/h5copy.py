'''

h5copy.py 

This program copies HDF5 file(s) to a designated output HDF5 file.  If the output file is
an existing HDF5 file, the user can specify whether the existing content should be erased
before copying begins, or whether the copied data should be appended to the existing data.


Some nomenclature
=================

A COMPAS output file (e.g. BSE_System_Parameters, BSE_RLOF, etc.) maps to an HDF5 GROUP,
where the group name is the name of the COMPAS output file.

A column in a COMPAS output file (e.g. SEED, Mass(1), Radius(2), etc.) maps to an HDF5 DATASET,
where the dataset name is the column heading string.

COMPAS column datatype strings are encoded in the dataset meta-details (dataset.dtype).
COMPAS column units strings are attached to HDF5 datasets as attributes.


Usage
=====

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




Functionality overview
======================

Output file
-----------

- If the specified HDF5 output file does not exist, it will be created, and data from the input
  file(s) will be appended to the new output file.

- If the specified HDF5 output file does exist, and the -e commandline option is not present, the
  existing output file will be preserved and data from the input file(s) will be appended to the
  existing output file if the output file was created with chunking enabled (see description of
  HDF5 chunking below - only files that were created with chunking enabled can be extended).

  Attempting to append data to an existing HDF5 file that was not created with chunking enabled
  will result in the following message being displayed for each dataset:

      "Only chunked datasets can be resized"

  If that happens, using this program to copy the output file to a new file will copy the existing
  data and create the newf file with chunking enabled - the newly created file can then be used as
  a base file to which other files can be appended.

- If the specified HDF5 output file does exist, and the -e commandline option is present, the
  existing output file will be deleted, and new, empty, file created, and data from the input
  file(s) will be appended to the new output file.

When appending data, if a dataset in an input file does not exist in the output file, the 
dataset will be created, otherwise the data copied will be appended to the existing dataset.


Input files
-----------

A list of input filenames and or directory names must be supplied.  The list can be a single name.
Each name in the list is processed in order.  If a name is the name of a file with the file extension
'.h5', the contents of the file will be copied to the output file, otherwise it will be ignored.
If a name is the name of a directory and the specified recursion level requires that the directory
be processed (see commandline option -r), the program will descend into the directory and process all
files and directories there (directories will be processed depending upon the value of the -r option),
otherwise it will be ignored.

The commandline option '-r' specifies whether recursion is enabled for directory processing, and if
it is, to what depth:

    if the -r option is not present, recursion is not enabled and only files in the specified working 
    directory will be candidates for copying

    if -r is present with no 'depth' value, recursion is enabled and the depth is not limited - that
    is, all files in the specified working directory, and all files in all directories below the 
    specified working directory, will be candidates for copying
    
    if -r is present with a specified 'depth' value, recursion is enabled and the depth is limited to
    the depth specified - that is, all files in the specified working directory, and all files in all 
    directories 'depth' levels below the specified working directory, will be candidates for copying
    

Input filename filter
---------------------

If the -f commandline option is present, the names of all candidate files will be checked against the
specified filter, and only files whose names match the filter will be copied.  The specified filter
is a filename-only filter - the file's path (i.e. it's location) will not be matched to the filter.
The specified filter should not include a file extension, but the program adds the extension '.h5' to
the specified filter - only files that have the file extension '.h5' will match the filter.

If the -f option is not present the program uses a default filter value of '*', then adds the '.h5'
file extension - so all candidate files with the '.h5' extension will be copied.


Exclude group list
------------------

If the -x commandline option is present, the specified list of groupnames will be excluded from data
copied from all input files.  If -x is not present, all groups in all candidate files will be copied.


Erase output
------------

If the -e commandline option is present and the output file (specified or default) exists, it will be
erased before copying begins.  The -e option is ignored if the output file does not exist.

If -e is not present and the output file (specified or default) exists, the existing content will be
preserved and any data copied to the file will be appended to the existing data.


HDF5 file chunking and IO
=========================

(mostly copied from the intro to COMPAS source file log.h, so some bits might be specific to the
COMPAS C++ code)

A brief description of HDF5 files and chunking, in the COMPAS context:
(With the caveat that all I know about HDF5 I learned so that I could write this code - it is
very likely that there are better ways to do things, so if anyone knows a better way, please
either change the code or tell me how to improve it and I'll change it.  Most of what follows
is just a brain dump of my reading/research over the past week or so, and it could very well
be based on my misunderstanding of what I have read - so uf anyone notices something I've
misunderstood please let me know so I can improve the code.)
 

Data in HDF5 files are arranged in groups and datasets.

    A COMPAS output file (e.g. BSE_System_Parameters, BSE_RLOF, etc.) maps to an HDF5 group,
    where the group name is the name of the COMPAS output file.

    A column in a COMPAS output file (e.g. SEED, Mass(1), Radius(2), etc.) maps to an HDF5 dataset,
    where the dataset name is the column heading string.

    COMPAS column datatype strings are encoded in the dataset meta-details (dataset.dtype).
    COMPAS column units strings are attached to HDF5 datasets as attributes.

Each dataset in an HDF5 files is broken into "chunks", where a chunk is defined as a number of dataset
entries.  In COMPAS, all datasets are 1-d arrays (columns), so a chunk is defined as a number of values
in the 1-d array (or column).  Chunking can be enabled or not, but if chunking is not* enabled a dataset
cannot be resized - so if chunking is not enabled the size of the dataset must be known at the time of
creation, and the entire datset created in one go.  That doesn't work for COMPAS - even though we know
the number of systems being evolved, we don't know the number of entries we'll have in each of the output
log files (and therefore the HDF5 datasests if we're logging to HDF5 files).  So, we enable chunking.
 
Chunking can improve, or degrade, performance depending upon how it is implemented - mostly related to
the chunk size chosen.
 
Datasets are stored inside an HDF5 file as a number of chunks - the chunks are not guaranteed (not even
likely) to be contiguous in the file or on the storage media (HDD, SSD etc.).  Chunks are mapped/indexed
in the HDF5 file using a B-tree, and the size of the B-tree, and therefore the traversal time, depends 
directly upon the number of chunks allocated for a dataset - so the access time for a chunk increases as
the number of chunks in the dataset increases.  So many small chunks will degrade performance.
 
Chunks are the unit of IO for HDF5 files - all IO to HDF5 is performed on the basis of chunks.  This means
that whenever dataset values are accessed (read or written (i.e. changed)), if the value is not already in
memory, the entire chunck containing the value must be read from, or written to, the storage media - even
if the dataset value being accessed is the only value in the chunk.  So few large chunks could cause 
empty, "wasted", space in the HDF5 files (at the end of datasets) - but they could also adversely affect
performance by causing unecessary IO traffic (although probably not much in the way we access data in COMPAS
files).
 
HDF5 files implement a chunk cache on a per-dataset basis.  The default size of the chunk cache is 1MB, and
its maximum size is 32MB.  The purpose of the chunk cache is to reduce storage media IO - even with SSDs,
memory access is way faster than storage media access, so the more of the file data that can be kept in
memory and maipulated there, the better.  Assuming the datatype of a particular dataset is DOUBLE, and
therefore consumes 8 bytes of storage space, at its maximum size the chunk cache for that dataset could hold
4,000,000 values - so a single chuck with 4,000,000 values, two chunks with 2,000,000 values, four with 
1,000,000, and so on.  Caching a single chunk defeats the purpose of the cache, so chunk sizes somewhat less
that 4,000,000 would be most appropriate if the chunk cache is to be utilised.  Chunks too big to fit in the
cache simply bypass the cache and are read from, or written to, the storage media directly.
 
However, the chunk cache is really only useful for random access of the dataset.  Most, if not all, of the
access in the COMPAS context (including post-creation analyses) is serial - the COMPAS code writes the
datasets from top to bottom, and later analyses (generally) read the datasets the same way.  Caching the
chunks for serial access just introduces overhead that costs memory (not much, to be sure: up to 32MB per 
open dataset), and degrades performace (albeit it a tiny bit).  For that reason I disable the chunk cache
in COMPAS - so all IO to/from an HDF5 file in COMPAS is directly to/from the storage media.  (To be clear,
post-creation analysis software can disable the cache or not when accessing HDF5 files created by COMPAS - 
disabling the cache here does not affect how other software accesses the files post-creation).
 
So many small chunks is not so good, and neither is just a few very large chunks.  So what's the optimum
chunk size?  It depends upon several things, and probably the most important of those are the final size
of the dataset and the access pattern.
 
As mentioned above, we tend to access datasets serially, and generally from to to bottom, so larger chunks 
would seem appropriate, but not so large that we generate HDF5 files with lots of unused space.  However, 
disk space, even SSD space, is cheap, so trading space against performance is probably a good trade.
 
Also as mentioned above, we don't know the final size of (most of) the datasets when creating the HDF5 in
COMPAS - though we do know the number of systems being generated, which allows us to determine an upper
bound for at least some of the datasets (though not for groups such as BSE_RLOF).
 
One thing we need to keep in mind is that when we create the HDF5 file we write each dataset of a group
in the same iteration - this is analagous to writing a single record in (e.g.) a CSV log file (the HDF5
group corresponds to the CSV file, and the HDF5 datasets in the group correspond to the columns in the
CSV file).  So for each iteration - typically each system evolved (though each timestep for detailed
output files) we do as many IOs to the HDF5 file as there are datasets in the group (columns in the file).
We are not bound to reading or writing a single chunk at a time - but we are bound to reading or writing
an integral multiple of whole chunks at a time.
 
We want to reduce the number of storage media accesses when writing (or later reading) the HDF5 files, so
larger chunk sizes are appropriate, but not so large that we create excessively large HDF5 files that have
lots of unused space (bearing in mind the trade-off mentioned above), especially when were evolving just
a few systems (rather than millions).
 
To really optimise IO performance for HDF5 files we'd choose chunk sizes that are close to multiples of
storage media block sizes, but I chose not to go down that rathole...
 
Based on everything written above, and some tests, I've chosen a default chunk size of 100,000 (dataset
entries) for all datasets (HDF5_DEFAULT_CHUNK_SIZE in constants.h) for the COMPAS C++ code.  This clearly
trades performance against storage space.  For the (current) default logfile record specifications, per-binary
logfile space is about 1K bytes, so in the very worst case we will waste some space at the end of a COMPAS
HDF5 output file, but the performance gain, especially for post-creation analyses, is significant.  Ensuring
the number of systems evolved is an integral multiple of this fixed chunk size will minimise storage space waste.

I have chosen a minimum chunk size of 1000 (HDF5_MINIMUM_CHUNK_SIZE in constants.h) for the COMPAS C++ code.  
If the number of systems being evolved is >= HDF5_MINIMUM_CHUNK_SIZE the chunk size used will be the value of 
the hdf5-chunk-size program option (either HDF5_DEFAULT_CHUNK_SIZE or a value specified by the user), but if 
the number of systems being evolved is < HDF5_MINIMUM_CHUNK_SIZE the chunk size used will be HDF5_MINIMUM_CHUNK_SIZE.
This is just so we don't waste too much storage space when running small tests - and if they are that small
performance is probably not going to be much of an issue, so no real trade-off against storage space.  

The chunk size chosen for the COMPAS C++ code determines the chunk size of the logfiles produced by the COMPAS
C++ code.  If those files are only ever given to this program as input files, their chunk size only matters
in that it affects the read performance of the files (the more chunks, and the more smaller chunks, in a
dataset of an input file means locating the chunks and reading them takes longer).  That may not be a huge
problem depending upon h ow many input files there are and how big they are.  If a COMPAS logfile is used as
a base file and other files are being appended to it, then the chunk size of the base output file will be the
chunk size used for writing to the file - that could affect performance if it is too small.  That's why the
chunk size is an option, both for this program and for the COMPAS C++ code.

Writing to the output HDF5 file is buffered here - we buffer a number of chunks for each open dataset and write
the buffer to the file when the buffer fills (or a partial buffer upon file close if the buffer is not full).
This IO buffering is not HDF5 or filesystem buffering - this is a an internal implementation here to improve
performance.  The IO buffer size is also an option (both here and in t he COMPAS C++ code).
 
Users should bear in mind that the combination of HDF5 chunk size and HDF5 IO buffer size affect performance,
storage space, and memory usage - so they may need to experiment to find a balance that suits their needs.


A note on string values
=======================

COMPAS writes string data to its HDF5 output files as C-type strings.  Python interprets C-type
strings in HDF5 files as byte arrays - regardless of the specified datatype when written (COMPAS
writes the strings as ASCII data (as can be seen with h5dump), but Python ignores that).  Note that
this affects the values in datasets (and attributes) only, not the dataset names (or group names,
attribute names, etc.).

The only real impact of this is that if the byte array is printed by Python, it will be displayed
as (e.g.) "b'abcde'" rather than just "abcde".  All operations on the data work as expected - it
is just the output that is impacted.  If that's an issue, use .decode('utf-8') to decode the byte
array as a Python string variable.  E.g.

    str = h5File[Group][Dataset][0]
    str is a byte array and print(str) will display (e.g.) b'abcde'

    but

    str = h5File[Group][Dataset][0].decode('utf-8')
    str is a Python string and print(str) will display (e.g.) abcde

Note: HDF5 files not created by COMPAS will not exhibit this behaviour, so for HDF5 files created by
the existing post-processing python scripts the use of .decode() is not only not necessary, it will
fail (because the strings in HDF5 files created by python are already python strings, not byte arrays).


JR, January 2021
'''

#!/usr/bin/env python3
import sys
import os
import numpy as np
import h5py as h5
import contextlib
import argparse
import gettext
from fnmatch import fnmatch


CHUNK_SIZE     = 100000         # HDF5 dataset chunk size (number of entries), minimum 1
IO_BUFFER_SIZE = 10             # number of HDF5 chunks, minimum 1


# copyHDF5File()
#
# Copies the contents of the file passed in 'path' parameter
# to the open HDF5 file specified by the 'outFile' parameter

def copyHDF5File(path, outFile, chunkSize = CHUNK_SIZE, bufferSize = IO_BUFFER_SIZE, excludeList = ''):

    ok = False                                                                                                          # result

    outFname = os.path.abspath(outFile.filename)                                                                        # fully-qualified output filename
    
    try:
        with h5.File(path, 'r') as srcFile:                                                                             # open the input HDF5 file

            srcFname = os.path.abspath(srcFile.filename)                                                                # fully-qualified source filename
            if srcFname == outFname:                                                                                    # is this source file the output file?
                ok = True                                                                                               # don't copy the output file...
            else:

                print('Copying file', srcFname)                                                                         # announce file being copied

                srcGroupNames = list(srcFile.keys())                                                                    # list of groups in input file

                for srcGroupName in srcGroupNames:                                                                      # for each source group

                    if srcGroupName in excludeList: continue                                                            # skip it if required

                    srcGroup = srcFile[srcGroupName]                                                                    # source group

                    try:
                        destGroup = outFile.require_group(srcGroupName)                                                 # open group in output file - create it if necessary

                        srcDatasetNames = list(srcGroup.keys())                                                         # list of datasets in srcGroup in input file

                        for srcDatasetName in srcDatasetNames:                                                          # for each source dataset in srcGroup

                            thisChunkSize  = chunkSize                                                                  # default dataset chunk size
                            thisBufferSize = bufferSize                                                                 # default IO block size

                            srcDataset       = srcGroup[srcDatasetName]                                                 # source dataset
                            srcDataset_dtype = srcDataset.dtype                                                         # source dataset data type
                            srcDatasetLen    = srcDataset.size                                                          # source dataset length

                            try:
                                datasetOpen = False                                                                     # no dataset open yet
                                try:                                                                                    # get destination dataset details if it exists
                                    destDataset    = destGroup[srcDatasetName]                                          # destination dataset
                                    destDatasetLen = destDataset.size                                                   # destination dataset length
                                    # use destination dataset chunk size if it is a chunked dataset
                                    # if it is not a chunked dataset the resize (below) will fail
                                    if destDataset.chunks is not None: thisChunkSize = destDataset.chunks[0]
                                    datasetOpen = True                                                                  # existing dataset is now open
                                except Exception as e:                                                                  # dataset does not exist - create it
                                    if chunkSize < 1: chunkSize = 1                                                     # clamp chunk size to minimum

                                    #create the dataset (with chunking enabled)
                                    try:
                                        destDataset    = destGroup.create_dataset(srcDatasetName, (0,), maxshape=(None,), chunks = (thisChunkSize,), dtype = srcDataset_dtype)
                                        destDatasetLen = destDataset.size                                               # destination dataset length
                                        datasetOpen    = True                                                           # newly created dataset is now open
                                    except:
                                        print('Error creating dataset', srcDatasetName, 'in group', srcGroupName, 'in file', outFname, ':', str(e))

                                if datasetOpen:                                                                         # dataset open ok?
                                    if thisBufferSize < 1: thisBufferSize = 1                                           # yes - clamp io minimum block size to 1 chunk
                                    thisBufferSize *= thisChunkSize                                                     # convert to number of entries

                                    srcDataset_attrs = list(srcDataset.attrs.items())                                   # list of dataset attributes

                                    for srcAttr in srcDataset_attrs:
                                        try:
                                            destDataset.attrs[srcAttr[0]] = srcAttr[1]                                  # set dataset attributes in destDataset - overwrites existing

                                            try:

                                                srcStart      = 0                                                       # source start position for copy
                                                srcEnd        = srcStart + thisBufferSize                               # source end position for copy
                                                destStart     = destDatasetLen                                          # destination start position for copy
                                                destEnd       = destStart + thisBufferSize                              # destination end position for copy

                                                while srcEnd <= srcDatasetLen:                                          # while source copy end position is inside source dataset
                                                    destDataset.resize((destStart + thisBufferSize,))                   # resize the destination dataset appropriately
                                                    destDataset[destStart : destEnd] = srcDataset[srcStart : srcEnd]    # copy source block to destination block

                                                    srcStart  = srcEnd                                                  # advance source start position for copy
                                                    srcEnd   += thisBufferSize                                          # advance source end position for copy
                                                    destStart = destEnd                                                 # advance destination start position for copy
                                                    destEnd  += thisBufferSize                                          # advance destination end position for copy

                                                if srcEnd > srcDatasetLen:                                              # source copy end position at or beyond end of source dataset?
                                                                                                                        # yes - last chunk (partial)
                                                    srcEnd  = srcDatasetLen                                             # set source end position for copy to end of dataset
                                                    destEnd = destStart + srcEnd - srcStart                             # set destination end position for copy appropriately
                                                    destDataset.resize((destEnd,))                                      # resize the destination dataset appropriately
                                                    destDataset[destStart : destEnd] = srcDataset[srcStart : srcEnd]    # copy source chunk to destination chunk

                                                ok = True                                                               # all good

                                            except Exception as e:                                                      # error occurred while writing to dataset
                                                print('Error writing to dataset', srcDatasetName, 'in group', srcGroupName, 'in file', outFname, ':', str(e))

                                        except Exception as e:                                                          # error occurred while accessing the dataset attributes
                                            print('Error accessing attribute', srcAttr[0], 'in dataset', srcDatasetName, 'in group', srcGroupName, 'in file', outFname, ':', str(e))

                            except Exception as e:                                                                      # error occurred while accessing the dataset
                                print('Error accessing dataset', srcDatasetName, 'in group', srcGroupName, 'in file', outFname, ':', str(e))

                    except Exception as e:
                        print('Error accessing group', srcGroupName, 'in file', outFname, ':', str(e))                  # error occurred while accessing the group

    except Exception as e:                                                                                              # error occurrd accessing the input file
        print('Error accessing HDF5 file', path, ':', str(e))

    return ok


# processDirectory()
#
# Processes file and directories in directory specified by 'path' parameter
# Recursion is controlled by 'recursive' and 'depth' parameters

def processDirectory(path, 
                     outFile, 
                     recursive   = 0, 
                     fileFilter  = '*.h5', 
                     stopOnError = False, 
                     depth       = 0, 
                     chunkSize   = CHUNK_SIZE, 
                     bufferSize  = IO_BUFFER_SIZE, 
                     excludeList = ''):

    ok = True                                                                                                       # result

    thisPath = os.path.abspath(path)                                                                                # absolute path

    try:
        for dirpath, dirnames, filenames in os.walk(thisPath):                                                      # walk directory
            absDirpath = os.path.abspath(dirpath)                                                                   # absolute path
            print('Processing directory', absDirpath)                                                               # announce directory being processed

            for filename in filenames:                                                                              # for each filename
                if fnmatch(filename, fileFilter):                                                                   # filename matches filter?
                    ok = copyHDF5File(absDirpath + '/' + filename, 
                                      outFile, 
                                      chunkSize   = chunkSize, 
                                      bufferSize  = bufferSize,
                                      excludeList = excludeList)                                                    # yes - copy the file

                if stopOnError and not ok: break                                                                    # check error - stop if required

            if stopOnError and not ok: break                                                                        # check error - stop if required

            if depth < recursive:                                                                                   # recurse?
                depth += 1                                                                                          # increment recursion depth

                for dirname in dirnames:                                                                            # for each directory
                    processDirectory(absDirpath + '/' + dirname, 
                                     outFile, 
                                     recursive   = recursive, 
                                     fileFilter  = fileFilter, 
                                     stopOnError = stopOnError, 
                                     depth       = depth, 
                                     chunkSize   = chunkSize, 
                                     bufferSize  = bufferSize,
                                     excludeList = excludeList)                                                     # process the directory

            break                                                                                                   # control recursion

    except Exception as e:                                                                                          # error occurred accessing directory
        print('Error accessing directory', thisPath, ':', str(e))                                                   # announce error
        ok = False                                                                                                  # we're done

    return ok


def main():

    ok = True

    # setup argument parser
    formatter = lambda prog: argparse.HelpFormatter(prog, max_help_position = 4, width = 90)
    parser = argparse.ArgumentParser(description = 'HDF5 file copier.', formatter_class = formatter)

    # define arguments
    parser.add_argument('inputPaths', metavar = 'input', type = str, nargs = '+', help = 'input directory and/or file name(s)')
    parser.add_argument('-b', '--buffer-size', dest = 'buffer_size', type = int, action = 'store',  default = IO_BUFFER_SIZE, help = 'IO buffer size (number of HDF5 chunks, default = ' + str(IO_BUFFER_SIZE) + ')')
    parser.add_argument('-c', '--chunk-size', dest = 'chunk_size', type = int, action = 'store',  default = CHUNK_SIZE, help = 'HDF5 output file dataset chunk size (default = ' + str(CHUNK_SIZE) + ')')
    parser.add_argument('-e', '--erase-output', dest = 'erase_ouput', action = 'store_true',  default = False, help = 'erase existing output file before copying input files (default = False)')
    parser.add_argument('-f', '--filter', dest = 'filename_filter', type = str, action = 'store',  default = '*', help = 'input filename filter (default = *)')
    parser.add_argument('-o', '--output', dest = 'output_fileName', type = str, action = 'store',  default = 'h5out.h5', help = 'output file name (default = h5out.h5)')
    parser.add_argument('-r', '--recursive', dest = 'recursion_depth', type = int, nargs = '?', action = 'store', default = 0, const = sys.maxsize,  help = 'recursion depth (default is no recursion)')
    parser.add_argument('-s', '--stop-on-error', dest = 'stop_on_error', action = 'store_true',  default = False, help = 'stop all copying if an error occurs (default is skip to next file and continue)')
    parser.add_argument('-x', '--exclude', dest = 'exclude_group', type = str, nargs = '+', action = 'store', default = '', help = 'list of input groups to be excluded (default is all groups will be copied)')

    # parse arguments
    args = parser.parse_args()

    # output filename
    outFname = args.output_fileName
    if outFname[-3:] != '.h5' and outFname[-3:] != '.H5': outFname += '.h5'                                         # add output file extension if necessary

    fileFilter = args.filename_filter + '.h5'                                                                       # add file extension to filer

    excludeList = ' '.join(args.exclude_group)                                                                      # construct exclude list for groups

    # check whether output file already exists
    # if it does exist, check whether it is an existing file or existing directory
    existingOutputFile = False
    fullOutFpath       = os.path.abspath(outFname)
    if os.path.exists(fullOutFpath):
        if os.path.isfile(fullOutFpath):
            if not args.erase_ouput: existingOutputFile = True
        elif os.path.isdir(fullOutFpath):
            print('Error encountered:', fullOutFpath, 'is the name of an existing directory - choose a different output filename')
            ok = False
        else:
            print('Error encountered:', fullOutFpath, 'is the name of an existing filesystem object - choose a different output filename')
            ok = False

    if ok:

        # disable HDF5 dataset cache
        # default cache size is 1MB - maximum is 32MB
        # we can change it here by multiplying the stored (default) value by some number (max 32...)
        # we don't need the cache, so just multiply by 0 to disable it - we save a little bit of memory
        # (the cache is per open dataset), and gain a (tiny) bit in performance
        try:
            h5FileAccessPropertyList = h5.h5p.create(h5.h5p.FILE_ACCESS)
            h5CacheSettings = list(h5FileAccessPropertyList.get_cache())
            h5CacheSettings[2] *= 0
            h5FileAccessPropertyList.set_cache(*h5CacheSettings)
        except Exception as e:                                                                                                  # error disabling cache
            print('Warning: unable to disable the HDF5 dataset cache:', str(e))                                                 # show warning

        # open the output file - create it if necessary
        # using low-level functions here so we can provide the propertly list
        if existingOutputFile:                                                                                                  # output file exists?
            try:                                                                                                                # yes
                outHDF5file = h5.h5f.open(fullOutFpath.encode('utf-8'), fapl = h5FileAccessPropertyList)                        # open it
            except Exception as e:                                                                                              # error opening file
                print('Error occurred while disabling HDF5 dataset cache:', str(e))                                             # announce error
                ok = false                                                                                                      # fail
        else:                                                                                                                   # output file does not exist
            try:
                outHDF5file = h5.h5f.create(fullOutFpath.encode('utf-8'), fapl = h5FileAccessPropertyList)                      # creat it
            except Exception as e:                                                                                              # error creating file
                print('Error occurred while disabling HDF5 dataset cache:', str(e))                                             # announce error
                ok = false                                                                                                      # fail

        if ok:
            try:
                with contextlib.closing(outHDF5file) as h5OutFid:                                                               # close the file when done...
                    # process input files and directories
                    with h5.File(h5OutFid) as outFile:                                                                          # get file id for high-level functions
                        for thisPath in args.inputPaths:                                                                        # for each input path
                            thisFullPath = os.path.abspath(thisPath)                                                            # fully-qualified filename
                            if os.path.exists(thisFullPath):                                                                    # path exists?
                                if os.path.isfile(thisFullPath):                                                                # yes - is it a file?
                                    if fnmatch(thisPath, fileFilter):                                                           # yes - filename matches filter?
                                        ok = copyHDF5File(thisFullPath, 
                                                          outFile, 
                                                          chunkSize   = args.chunk_size, 
                                                          bufferSize  = args.buffer_size, 
                                                          excludeList = excludeList)                                            # yes - copy it
                                    else:                                                                                       # no - does not match filter
                                        print('Warning:', thisPath, 'does not match file filter (', fileFilter, '): ignored')   # show warning
                                elif os.path.isdir(thisFullPath):                                                               # not a file - directory?
                                                                                                                                # yes - process directory
                                    ok = processDirectory(thisFullPath, 
                                                          outFile, 
                                                          recursive   = args.recursion_depth, 
                                                          fileFilter  = fileFilter, 
                                                          stopOnError = args.stop_on_error, 
                                                          chunkSize   = args.chunk_size, 
                                                          bufferSize  = args.buffer_size, 
                                                          excludeList = excludeList)
                                else:                                                                                           # not a file or directory
                                    print('Warning:', thisFullPath, 'is not a file or a directory: ignored')                    # show warning
                            else:                                                                                               # path does not exist
                                print('Warning:', thisFullPath, 'does not exist: ignored')                                      # show warning

                            if args.stop_on_error and not ok:                                                                     # error occurred, and stop-on-error specified?
                                print('Error encountered: copy stopped')                                                        # yes - announce error
                                break                                                                                           # and stop

                if not existingOutputFile and os.path.getsize(fullOutFpath) <= 0:                                               # is output file a new file and empty?
                    os.remove(fullOutFpath)                                                                                     # yes - didn't write to it - delete it

            except Exception as e:                                                                                              # error occurred
                print('Error:', str(e))                                                                                         # announce error

if __name__ == "__main__":
    main()   
