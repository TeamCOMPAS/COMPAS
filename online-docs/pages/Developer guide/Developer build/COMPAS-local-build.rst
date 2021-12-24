Building COMPAS locally
=======================

A makefile is provided to build COMPAS locally (``Mkaefile``).  The Makefile provided defines a number of variables that can be 
specified on the command line when ``make`` is run, including variables that allow the user to specify the compiler, the `include` 
directory (for source header files), the `lib` directory (for shared libraries) for each external library required by COMPAS, and 
the COMPAS executable name:

    **CPP** |br|
    Specifies the compiler to be used.  The default value is 'g++'.

    **GSLINCDIR** |br|
    Specifies the `include` directory for the GNU scientific library, ``gsl``.  The default value is '/include'.

    **GSLLIBDIR** |br|
    Specifies the `lib` directory for the GNU scientific library, ``gsl``.  The default value is '/lib'.

    **BOOSTINCDIR** |br|
    Specifies the `include` directory for the BOOST library.  The default value is '/include'
    
    **BOOSTLIBDIR** |br|
    Specifies the `lib` directory for the BOOST library.  The default value is '/lib'.
    
    **HDF5INCDIR** |br|
    Specifies the `include` directory for the HDF5 library.  The default value is '/usr/include/hdf5/serial'

    **HDF5LIBDIR** |br|
    Specifies the `lib` directory for the HDF5 library.  The default value is '/usr/lib/x86_64-linux-gnu/hdf5/serial'
    
    **EXE** |br|
    Specifies the name of the COMPAS exectuable to be built.  The default is 'COMPAS'


For example, typing::

    make GCC=c++ EXE=mycompas -j$(nproc)
    
will cause the `c++` compiler to be used to create the executable file 'mycompas', using all available cores.


The makefile provided also defines several entry points:

    **clean** |br|
    Instructs ``make`` to remove all existing object files (.o), and the COMPAS executable.  A subsequent ``make`` is then forced to 
    compile all source files and link the resultant object files (and external libraries) into a new executable.

    **static** |br|
    Specifies that functions in the external libraries referenced by COMPAS should be statically linked - that is, they are copied into
    the COMPAS executable.  The default executable name for the *static* entry point is `COMPAS_STATIC`.

    **fast** |br| 
    Adds `-march=native` and `-O3` to the compiler flags.

        - specifying `-march=native` enables all instruction subsets supported by the compiling machine, thus producing an executable
          file that will perform better than it would if the full native instruction set was not available.  Care should be taken when
          using the executable created from this entry point - some of the native instructions may not be available for use on machines
          of different architectures, so the resultant executable file is not necessarily portable.
        - specifying `-O3` causes the compiler to perform many optimisations to produce an executable that is optimised for performance
          (at the expense of compile time).\ [#f1]_

    **staticfast** |br|
    The functionality of the **static** and **fast** entry points combined.



.. rubric:: Footnotes

.. [#f1] https://gcc.gnu.org/onlinedocs/gcc/Optimize-Options.html

