# COMPAS_HPC

------------

Here are some basic instructions for setting up and running COMPAS_HPC. COMPAS is our binary stellar evolution code. COMPAS_HPC is a suite of tools written in python to enable running COMPAS on High Performance Computers (HPC), also known as supercomputers (see supported below). 

## Supported HPC Facilities

-------------------------

[OzSTAR](https://supercomputing.swin.edu.au) (Swinburne) - Slurm

Tsunami (Birmingham) - Condor

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

*Note:* Other versions of these packages may work but have not necessarily been tested.

*Note:* Some HPCs require that the necessary packages be sourced at login, e.g in the `~/.bash_profile` for bash users. If you run into version issues, ensure that your .bash_profile (or equivalent) is pointed to the correct packages.

## Instructions

---------------

COMPAS is [hosted on github.](https://github.com/TeamCOMPAS/COMPAS/) 
 
CompasHPC is a suite of python tools written to make running COMPAS on High Performance Computers (HPC) very easy. 
 
The HPC code is in the folder **`$COMPAS_ROOT_DIR/CompasHPC/`** which is abbreviated to **`$CHPC`** below for compactness and readability. 
 
User modifiable input files are provided as templates in the folder **`$COMPAS_ROOT_DIR/defaults/`** abbreviated as **`$Cdefs`**

#### 0) Compiling an optimised version of COMPAS

If you are interested in CompasHPC, you are likely considering simulating many binaries. This can be computationally expensive. The Makefile included with COMPAS includes a mode which will optimise the compilation of COMPAS for the machine it is to be run on; in our benchmarking tests, this can lead to speed ups of a factor of 2--3. 

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

#### 1) CompasHPC
   
  There are 3 input files required for CompasHPC, 

     $CHPC/compas_hpc_input.py

     $CHPC/masterFolder/pythonSubmit.py

	 $CHPC/postProcessing.py
   
  By default, these files do not exist. We provide the following templates
   
     $Cdefs/compas_hpc_input_default.py

     $Cdefs/pythonSubmitDefault.py
   
     $Cdefs/postProcessingDefault.py

  Which you should not edit. Create copies of these files by doing
   
     cp $Cdefs/compas_hpc_input_default.py $CHPC/compas_hpc_input.py

     cp $Cdefs/pythonSubmitDefault.py $CHPC/masterFolder/pythonSubmit.py
    
     cp $Cdefs/postProcessingDefault.py $CHPC/postProcessing.py

  Now edit these with your desired settings (see next sections).

#### 2) `compas_hpc_input.py`

   This contains 'meta settings' such as where you want to send the output of your jobs and how many jobs to run, along with other settings. Make sure to set `cluster = 'ozstar' # or tsunami`.
    
   On OzSTAR, you need to send your output to fred e.g. `rootOutputDir = '/fred/oz101/sstevens/popsynth/test'` and make sure nBatches is small to start with (say, 3). Because COMPAS is embarassingly parallel, we can gain a factor nBatches speed up by splitting any job in to nBatches pieces, which can be very powerful if nBatches is say 100. It can also cause lots of problems very quickly, so be careful!
 
#### 3) `pythonSubmit.py`

   This contains the 'physics settings' that are sent to COMPAS when you run a job. You can choose how many binaries to evolve, what metallicity to use, what mass distribution and so on. There's a lot there, I won't try to explain it all here. Try 
   
     $COMPAS_ROOT_DIR/src/COMPAS --help
   
   for more information. For the default run using self-generated random seeds, note that the total number of binaries evolved will be the product of the variables `number_of_binaries` in pythonSubmit.py and `nBatches` in compas_hpc_input.py.
 
#### 4) `postProcessing.py`

   This contains settings related to the output data files from COMPAS. 
The script converts the plaintext output into a more compact HDF5 file, 
so if you have a small run and plan to work directly with the plaintext files, you may not need to worry about this.
You should modify this only if you changed the output prefix, delimiter, or extension in the pythonSubmit file.

#### 5) Run `compas_hpc.py`

   Finally, run
    
     python $CHPC/compas_hpc.py
    
   to run COMPAS. This will automatically submit the job to the relevant job scheduler for your HPC facility. On OzSTAR that is Slurm, on G2 that is PBS and on tsunami that is Condor. There should be some pretty verbose output showing you the commands the code is running.


