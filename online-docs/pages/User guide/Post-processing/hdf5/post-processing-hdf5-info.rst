HDF5 basics
===========

Here we provide some basic information regarding the COMPAS ``HDF5`` log files, and how COMPAS interacts with ``HDF5``.
Interested readers can learn more about the ``HDF5`` file format at `The HDF Group <https://www.hdfgroup.org/>`__.


.. _HDF5-chunking:

HDF5 file chunking and IO
-------------------------

Following is a brief description of ``HDF5`` files and ``chunking``, in the COMPAS context.

Data in ``HDF5`` files are arranged in ``groups`` and ``datasets``:

    A COMPAS output file (e.g. `BSE_System_Parameters`, `BSE_RLOF`, etc.) maps to an ``HDF5 group`` where the group name is
    the name of the COMPAS output file.

    A column in a COMPAS output file (e.g. `SEED`, `Mass(1)`, `Radius(2)`, etc.) maps to an ``HDF5 dataset``, where the
    dataset name is the column heading string.

    COMPAS column datatype strings are encoded in the ``HDF5 dataset`` meta-details (``dataset.dtype``).

    COMPAS column units strings are attached to ``HDF5 datasets`` as ``attributes``.

Each dataset in an ``HDF5`` file is broken into `chunks`, where a chunk is defined as a number of dataset entries. In COMPAS,
all datasets are 1-d arrays (columns), so a chunk is defined as a number of values in the 1-d array (or column). Chunking can 
be enabled or not, but if chunking is not enabled a dataset cannot be resized - so if chunking is not enabled the size of the 
dataset must be known at the time of creation, and the entire dataset created in one go. That doesn't work for COMPAS - even 
though we know the number of systems being evolved, we don't know the number of entries we'll have in each of the output log 
files (and therefore the ``HDF5`` datasests if we're logging to ``HDF5`` files).  So, we enable chunking.
 
Chunking can improve, or degrade, performance depending upon how it is implemented - mostly related to the chunk size chosen.
 
``Datasets`` are stored inside an ``HDF5`` file as a number of chunks - the chunks are not guaranteed (not even likely) to be
contiguous in the file or on the storage media (HDD, SSD etc.). Chunks are mapped/indexed in the ``HDF5`` file using a B-tree,
and the size of the B-tree, and therefore the traversal time, depends directly upon the number of chunks allocated for a 
dataset - so the access time for a chunk increases as the number of chunks in the dataset increases. So many small chunks will
degrade performance.
 
``Chunks`` are the unit of IO for ``HDF5`` files - all IO to ``HDF5`` is performed on the basis of chunks. This means that 
whenever dataset values are accessed (read or written (i.e. changed)), if the value is not already in memory, the entire chunk
containing the value must be read from, or written to, the storage media - even if the dataset value being accessed is the only
value in the chunk. So few large chunks could cause empty, "wasted", space in the ``HDF5`` files (at the end of datasets) - but
they could also adversely affect performance by causing unnecessary IO traffic (although probably not much in the way we access 
data in COMPAS files).
 
``HDF5`` files implement a chunk cache on a per-dataset basis. The default size of the chunk cache is $\small{1MB}$, and its 
maximum size is $\small{32MB}$. The purpose of the chunk cache is to reduce storage media IO - even with SSDs, memory access is 
much faster than storage media access, so the more of the file data that can be kept in memory and maipulated there, the better.  
Assuming the datatype of a particular dataset is ``DOUBLE``, and therefore consumes $\small{8}$ bytes of storage space, at its 
maximum size the chunk cache for that dataset could hold $\small{4,000,000}$ values - so a single chunk with $\small{4,000,000}$ 
values, two chunks with $\small{2,000,000}$ values, four with $\small{1,000,000}$, and so on. Caching a single chunk defeats the 
purpose of the cache, so chunk sizes somewhat less that $\small{4,000,000}$ would be most appropriate if the chunk cache is to be 
utilised. Chunks too big to fit in the cache simply bypass the cache and are read from, or written to, the storage media directly.
 
However, the chunk cache is really only useful for random access of the dataset. Most, if not all, of the access in the COMPAS 
context (including post-creation analyses) is serial - the COMPAS code writes the datasets from top to bottom, and later analyses
(generally) read the datasets the same way. Caching the chunks for serial access just introduces overhead that costs memory (not 
much, to be sure: up to $\small{32MB}$ per open dataset), and degrades performance (albeit it a tiny bit). For that reason we disable
the chunk cache in COMPAS - so all IO to/from an ``HDF5`` file in COMPAS is directly to/from the storage media. (To be clear, 
post-creation analysis software can disable the cache or not when accessing ``HDF5`` files created by COMPAS - disabling the cache
here does not affect how other software accesses the files post-creation).
 
So many small chunks is not so good, and neither is just a few very large chunks. So what's the optimum chunk size? That depends 
upon several things, and probably the most important of those are the final size of the dataset and the access pattern.
 
As mentioned above, we tend to access datasets serially, and generally from top to bottom, so larger chunks would seem appropriate, 
but not so large that we generate ``HDF5`` files with lots of unused space. However, disk space, even SSD space, is cheap, so 
trading space against performance is probably a good trade.
 
Also as mentioned above, we don't know the final size of (most of) the datasets when creating the ``HDF5`` in COMPAS - though we do
know the number of systems being generated, which allows us to determine an upper bound for at least some of the datasets (though 
not for groups such as `BSE_RLOF`).
 
One thing we need to keep in mind is that when we create the ``HDF5`` file we write each dataset of a group in the same 
iteration - this is analogous to writing a single record in (e.g.) a ``CSV`` log file (the ``HDF5 group`` corresponds to the ``CSV``
file, and the ``HDF5 datasets`` in the group correspond to the columns in the ``CSV`` file). So for each iteration - typically each
system evolved (though each timestep for detailed output files) we do as many IOs to the ``HDF5`` file as there are datasets in the
group (columns in the file). We are not bound to reading or writing a single chunk at a time - but we are bound to reading or writing
an integral multiple of whole chunks at a time.
 
We want to reduce the number of storage media accesses when writing (or later reading) the ``HDF5`` files, so larger chunk sizes are 
appropriate, but not so large that we create excessively large ``HDF5`` files that have lots of unused space (bearing in mind the 
trade-off mentioned above), especially when we're evolving just a few systems (rather than millions).
 
To really optimise IO performance for ``HDF5`` files we'd choose chunk sizes that are close to multiples of storage media block sizes,
but that would be too problematic given the number of disparate systems COMPAS could be run on...
 
Based on everything written above, and some tests, we've chosen a default chunk size of $\small{100,000}$ (dataset entries) for all 
datasets (``HDF5_DEFAULT_CHUNK_SIZE`` in ``constants.h``) for the COMPAS C++ code. This clearly trades performance against storage space.
For the (current) default logfile record specifications, per-binary logfile space is about $\small{1K}$ bytes, so in the very worst case
we will waste some space at the end of a COMPAS ``HDF5`` output file, but the performance gain, especially for post-creation analyses, is 
significant. Ensuring the number of systems evolved is an integral multiple of this fixed chunk size will minimise storage space waste.

We have chosen a minimum chunk size of $\small{1,000}$ (``HDF5_MINIMUM_CHUNK_SIZE`` in ``constants.h``) for the COMPAS C++ code. If the
number of systems being evolved is not less than ``HDF5_MINIMUM_CHUNK_SIZE`` the chunk size used will be the value of the ``hdf5-chunk-size``
program option (either ``HDF5_DEFAULT_CHUNK_SIZE`` or a value specified by the user), but if the number of systems being evolved is less
than ``HDF5_MINIMUM_CHUNK_SIZE`` the chunk size used will be ``HDF5_MINIMUM_CHUNK_SIZE``. This is just so we don't waste too much storage 
space when running small tests - and if they are that small performance is probably not going to be much of an issue, so no real trade-off 
against storage space.  


Copying and concatenating HDF5 files with h5copy.py
---------------------------------------------------

The chunk size chosen for the COMPAS C++ code determines the chunk size of the logfiles produced by the COMPAS C++ code.  If those files
are only ever given to programs such as ``h5copy`` as input files, their chunk size only matters in that it affects the read performance 
of the files (the more chunks, and the more smaller chunks, in a dataset of an input file means locating the chunks and reading them takes 
longer).  That may not be a huge problem depending upon how many input files there are and how big they are.  If a COMPAS logfile is used 
as a base file and other files are being appended to it via ``h5copy``, then the chunk size of the base output file will be the chunk size 
used for writing to the file - that could affect performance if it is too small.  We provide command-line options to specify the chunk size 
in both ``h5copy`` and the COMPAS C++ code so that users have some control over chunksize and performance.

Writing to the output HDF5 file is buffered in both ``h5copy`` and the COMPAS C++ code - we buffer a number of chunks for each open dataset 
and write the buffer to the file when the buffer fills (or a partial buffer upon file close if the buffer is not full). This IO buffering is 
not ``HDF5`` or filesystem buffering - this is a an internal implementation of ``h5copy`` and the COMPAS C++ code to improve performance.  
The IO buffer size can be changed via command-line options in both ``h5copy`` and the COMPAS C++ code.
 
Users should bear in mind that the combination of ``HDF5`` chunk size and ``HDF5`` IO buffer size affect performance, storage space, and 
memory usage - so they may need to experiment to find a balance that suits their needs.


A note on string values
-----------------------

COMPAS writes string data to its ``HDF5`` output files as C-type strings.  Python interprets C-type strings in ``HDF5`` files as byte 
arrays - regardless of the specified datatype when written (COMPAS writes the strings as ASCII data (as can be seen with ``h5dump``), but 
Python ignores that).  Note that this affects the values in datasets (and attributes) only, not the dataset names (or group names,
attribute names, etc.).

The only real impact of this is that if the byte array is printed by Python, it will be displayed as (e.g.) "b'abcde'" rather than just 
"abcde".  All operations on the data work as expected - it is just the output that is impacted.  If that's an issue, use ``.decode('utf-8')``
to decode the byte array as a Python string variable.

For example::

    str = h5File[Group][Dataset][0]

Here str is a byte array and ``print(str)`` will display (e.g.) ``b'abcde'``, but::

    str = h5File[Group][Dataset][0].decode('utf-8')

Here str is a Python string and ``print(str)`` will display (e.g.) ``abcde``

Note: ``HDF5`` files not created by COMPAS will not (necessarily) exhibit this behaviour, so for ``HDF5`` files created by the existing 
post-processing Python scripts the use of ``.decode()`` is not only not necessary, it will fail (because the strings in ``HDF5`` files 
created by Python are already Python strings, not byte arrays).

