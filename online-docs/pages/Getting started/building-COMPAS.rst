Building COMPAS
===============

We first need to define an environment variable for the root directory of COMPAS in your shell start-up file for COMPAS to run properly. For example, 
if you use bash as your shell, open `~/.bashrc` with a text editor and put in the following::

    export COMPAS_ROOT_DIR=~/codes/COMPAS

where `~/codes` should be replaced with the path to the directory where you cloned the COMPAS repository. For this to take effect, either restart your 
bash session or run::

    source ~/.bashrc

If your shell is ``zsh`` (which is the default of macOS 10.15), set the environment variable as above in `~/.zshrc` instead of `~/.bashrc`. If your shell
is ``csh``, set the environment variable in `~/.cshrc` using::

    setenv COMPAS_ROOT_DIR ~/codes/COMPAS
    
Now go to the COMPAS source code directory::

    cd $COMPAS_ROOT_DIR/src

In this directory you will find the file ``Makefile``, which you need to edit to point to your ``gsl``, ``boost``, and ``hdf5`` include files and libraries. 

If you installed the packages with Homebrew, the package files are likely to be found in /usr/local/opt (in directories gsl, boost, and hdf5 respectively),
but if they are not found there you will need to use Homebrew to locate the files::

    $ brew info boost
    boost: stable 1.72.0 (bottled), HEAD
    Collection of portable C++ source libraries
    https://www.boost.org/
    /usr/local/Cellar/boost/1.72.0 (14,466 files, 648.5MB) *
    ...

Copy the path, which in this case is `/usr/local/Cellar/boost/1.72.0`, and add it to the appropriate lines of the Makefile::

    BOOSTINCDIR = /usr/local/Cellar/boost/1.72.0/include
    BOOSTLIBDIR = /usr/local/Cellar/boost/1.72.0/lib
 
To build the COMPAS executable (compile and link) type::

    make -f Makefile

The build process will run much faster if multiple processors/cores are available. To build the COMPAS executable using (e.g.) 4 cores, type::

    make -j 4 -f Makefile

Note that both ``make`` commands shown above will conduct incremental builds: they will only compile source files that have changed. To ensure a clean build
in which all source files are compiled, type::

    make clean
    make -j 4 -f Makefile

The `clean` option instructs ``make`` to remove all existing object files (.o), and the COMPAS executable.  A subsequent ``make`` is then forced to compile
all source files and link the resultant object files (and external libraries) into a new executable.

The executable can be tested with, e.g.,

    ./COMPAS -v

which will display the code version.

See :doc:`../Developer guide/Developer build/COMPAS-local-build` for a detailed description of ``Makefile`` functionality.


:bolditalictext:`A note for Mac users:`

If you are using MacOS and running into linking issues with the boost libraries, try::

    make clean
    make CPP=clang++ -j$(sysctl -n hw.ncpu)

In some Mac installations, the GNU C++ compiler is not installed how we might expect, so trying to compile and link with ``clang++`` might help.

