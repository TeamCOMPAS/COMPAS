COMPAS_HPC
------------
Here are some basic instructions for setting up and running COMPAS_HPC. COMPAS is our binary stellar evolution code. COMPAS_HPC is a suite of tools written in python to enable running COMPAS on High Performance Computers (HPC), also known as supercomputers (see supported below). 

Supported HPC Facilities
-------------------------
`OzSTAR <https://supercomputing.swin.edu.au>`_ (Swinburne) - Slurm

Tsunami (Birmingham) - Condor

Requirements
--------------
I will assume here a basic knowledge of linux commands, a working COMPAS installation and access to one of the above HPC facilities.

In addition to the standard COMPAS prerequisites, COMPAS_HPC requires:
python version 3.x
numpy version ..


Instructions
---------------
COMPAS is hosted on github here:
https://github.com/TeamCOMPAS/COMPAS/
 
1) CompasHPC
 
CompasHPC is a suite of python tools written to make running COMPAS on High Performance Computers (HPC) very easy. The code is in:
 
`$COMPAS_ROOT_DIR/CompasHPC/`

There are 2 input files required for CompasHPC, `$COMPAS_ROOT_DIR/CompasHPC/compas_hpc_input.py` and `$COMPAS_ROOT_DIR/CompasHPC/masterFolder/pythonSubmit.py`

By default, these files do not exist. We provide the following templates

`$COMPAS_ROOT_DIR/CompasHPC/compas_hpc_input_default.py` 
`$COMPAS_ROOT_DIR/CompasHPC/masterFolder/pythonSubmitDefailt.py`

Which you should not edit. Create copies of these files by doing

`cp $COMPAS_ROOT_DIR/CompasHPC/compas_hpc_input_default.py COMPAS_ROOT_DIR/CompasHPC/compas_hpc_input.py`
`cp $COMPAS_ROOT_DIR/CompasHPC/masterFolder/pythonSubmitDefault.py $COMPAS_ROOT_DIR/CompasHPC/masterFolder/pythonSubmit.py`
 
Now edit these with your desired settings (see next sections).

2) `compas_hpc_input.py`

This contains 'meta settings' such as where you want to send the output of your jobs and how many jobs to run, along with other settings. Make sure to set
 
`cluster = 'ozstar' # or tsunami`
 
On OzSTAR, you need to send your output to fred e.g.
 
`rootOutputDir = '/fred/oz101/sstevens/popsynth/test'`
 
and make sure nBatches is small to start with (say, 3). Because COMPAS is embarassingly parallel, we can gain a factor nBatches speed up by splitting any job in to nBatches pieces, which can be very powerful if nBatches is say 100. It can also cause lots of problems very quickly, so be careful!
 
3) `pythonSubmit.py`

This contains the 'physics settings' that are sent to COMPAS when you run a job. You can choose how many binaries to evolve, what metallicity to use, what mass distribution and so on. There's a lot there, I won't try to explain it all here. Try 

`$COMPAS_ROOT_DIR/COMPAS/COMPAS --help`

For more information.
 
4) Run `compas_hpc.py`

Finally, run
 
`python $COMPAS_ROOT_DIR/CompasHPC/compas_hpc.py`
 
to run COMPAS. This will automatically submit the job to the relevant job scheduler for your HPC facility. On OzSTAR that is Slurm, on G2 that is PBS and on tsunami that is Condor. There should be some pretty verbose output showing you the commands the code is running.
