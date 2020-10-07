# COMPAS_HPC

------------

Here are some basic instructions for setting up and running COMPAS_HPC. COMPAS is our binary stellar evolution code. COMPAS_HPC is a suite of tools written in python to enable running COMPAS on High Performance Computers (HPC), also known as supercomputers (see supported below). 

## Supported HPC Facilities

-------------------------

[OzSTAR](https://supercomputing.swin.edu.au) (Swinburne) - Slurm

## Requirements

--------------

A basic knowledge of linux commands is assumed, as well as access to one of the above HPC facilities. You will also need a working COMPAS installation, meaning the source code has been compiled and you can run `./COMPAS` in the `src/` directory without error.
In addition to the standard COMPAS prerequisites, COMPAS_HPC requires:

- python 3.6 
- numpy 1.14 
- h5py 2.7 
- pandas 0.22 
- astropy 3.1 
- scipy 1.0 
- openblas 0.2 
- fftw 3.3 
- scalapack 2.0 
- sqlite 3.21 
- szip 2.1 
- hdf5 1.10 
- stroopwafel 1.0

*Note:* Other versions of these packages may work but have not necessarily been tested.

*Note:* Some HPCs require that the necessary packages be sourced at login, e.g in the `~/.bash_profile` for bash users. If you run into version issues, ensure that your .bash_profile (or equivalent) is pointed to the correct packages.
For users on OzStar, a COMPAS_SDK file will source all of the required prerequisites at runtime if it is called in ~/.bash_profile, with the exception of stroopwafel 1.0. For that, you will need to run `pip install stroopwafel --user`.

## Instructions

---------------

COMPAS is [hosted on github.](https://github.com/TeamCOMPAS/COMPAS/) 
 
CompasHPC is a suite of python tools written to make running COMPAS on High Performance Computers (HPC) very easy. 
 
The HPC code is in the folder **`$COMPAS_ROOT_DIR/CompasHPC/`** which is abbreviated to **`$CHPC`** below for compactness and readability. 
 
User modifiable input files are provided as templates in the folder **`$COMPAS_ROOT_DIR/defaults/`** abbreviated as **`$Cdefs`**

#### Compiling an optimised version of COMPAS

If you are interested in CompasHPC, you are likely considering simulating many binaries. This can be computationally expensive. The Makefile included with COMPAS includes a mode which will optimise the compilation of COMPAS for the machine it is to be run on; in our benchmarking tests, this can lead to speed ups of a factor of 2-3. 

To utilise this option, recompile COMPAS by doing

	cd $COMPAS_ROOT_DIR/src
	make clean
	make fast

This compilation takes longer than the standard compilation, which is why it is off by default. To speed this up, you can utilise all available cores on your machine by doing (on linux systems)

	make fast -j $(nproc)

On MacOS, `nproc` does not exist by default. You can overcome this either by specifying the number by hand if you know it e.g.

	make fast -j 4

or by adding

	alias nproc='sysctl -n hw.ncpu'

to your `~/.bash_profile` or equivalent.

#### Run COMPAS on HPC
   
  Running COMPAS on an HPC is very simple. 

  Make a copy of the `stroopwafel_interface.py` file in the `defaults/` folder into your current directory, and set the `run_on_helios` parameter to `True`.
  Set any other stroopwafel parameters as you see fit.

  If you have many non-default COMPAS arguments, you are encouraged to set them in a `pythonSubmit.py` file in the same directory, 
  and set the `usePythonSubmit` parameter to `True`.

  See sampling.md for details.

Currently, we only support the Slurm Workload Manager. If you are running on an HPC that does not use Slurm and would like assitance configuring COMPAS for this purporse, 
please contact us at [COMPAS User Google Group](https://groups.google.com/forum/#!members/compas-user).
