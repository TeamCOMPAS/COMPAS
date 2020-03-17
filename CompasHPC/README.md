COMPAS_HPC
------------
Here are some basic instructions for setting up and running COMPAS_HPC. COMPAS is our binary stellar evolution code. COMPAS_HPC is a suite of tools written in python to enable running COMPAS on High Performance Computers (HPC), also known as supercomputers (see supported below). 

Supported HPC Facilities
-------------------------
[OzSTAR](https://supercomputing.swin.edu.au) (Swinburne) - Slurm

Tsunami (Birmingham) - Condor

Requirements
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

Instructions
---------------
COMPAS is [hosted on github.](https://github.com/TeamCOMPAS/COMPAS/) 

CompasHPC is a suite of python tools written to make running COMPAS on High Performance Computers (HPC) very easy. 

The HPC code is in the folder `$COMPAS_ROOT_DIR/CompasHPC/` which is abbreviated to `$CHPC` below for compactness and readability.
 
1) CompasHPC
   
   There are 2 input files required for CompasHPC, 

     `$CHPC/compas_hpc_input.py` and

     `$CHPC/masterFolder/pythonSubmit.py`
   
   By default, these files do not exist. We provide the following templates
   
     `$CHPC/compas_hpc_input_default.py` 

     `$CHPC/masterFolder/pythonSubmitDefault.py`
   
   Which you should not edit. Create copies of these files by doing
   
     `cp $CHPC/compas_hpc_input_default.py $CHPC/compas_hpc_input.py` and

     `cp $CHPC/masterFolder/pythonSubmitDefault.py $CHPC/masterFolder/pythonSubmit.py`
    
   Now edit these with your desired settings (see next sections).

2) `compas_hpc_input.py`

   This contains 'meta settings' such as where you want to send the output of your jobs and how many jobs to run, along with other settings. Make sure to set `cluster = 'ozstar' # or tsunami`.
    
   On OzSTAR, you need to send your output to fred e.g. `rootOutputDir = '/fred/oz101/sstevens/popsynth/test'` and make sure nBatches is small to start with (say, 3). Because COMPAS is embarassingly parallel, we can gain a factor nBatches speed up by splitting any job in to nBatches pieces, which can be very powerful if nBatches is say 100. It can also cause lots of problems very quickly, so be careful!
 
3) `pythonSubmit.py`

   This contains the 'physics settings' that are sent to COMPAS when you run a job. You can choose how many binaries to evolve, what metallicity to use, what mass distribution and so on. There's a lot there, I won't try to explain it all here. Try 
   
     `$CRD/src/COMPAS --help`
   
   for more information. For the default run using self-generated random seeds, note that the total number of binaries evolved will be the product of the variables `number_of_binaries` in pythonSubmit.py and `nBatches` in compas_hpc_input.py.
 
4) Run `compas_hpc.py`

   Finally, run
    
     `python $CHPC/compas_hpc.py`
    
   to run COMPAS. This will automatically submit the job to the relevant job scheduler for your HPC facility. On OzSTAR that is Slurm, on G2 that is PBS and on tsunami that is Condor. There should be some pretty verbose output showing you the commands the code is running.


