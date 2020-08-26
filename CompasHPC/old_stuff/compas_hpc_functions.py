#####################################################
#
# compas_hpc_functions.py
#
# This file contains functions used by compas_hpc.py
# 
#####################################################

import numpy as np
import os
import time
from subprocess import Popen, PIPE
import sys
import pickle

#####################################################
#-- Templates

#-- old slurm string template without using array
# SlurmJobStringTemplate="""#!/bin/bash
# #SBATCH --job-name=%s
# #SBATCH -N %s
# #SBATCH -D %s
# #SBATCH --output=%s
# #SBATCH --error=%s
# #SBATCH --time=%s
# #SBATCH --mail-user=%s
# #SBATCH --mail-type=ALL
# echo $SLURM_JOB_ID
# echo $SLURM_JOB_NAME
# echo $SLURM_SUBMIT_DIR
# cd $SLURM_SUBMIT_DIR
# %s""" 

SlurmJobStringTemplateWithEmail="""#!/bin/bash
#SBATCH --job-name=%s
#SBATCH --nodes=%s
#SBATCH --ntasks=%s
#SBATCH --output=%s
#SBATCH --error=%s
#SBATCH --time=%s
#SBATCH --mem=%s
#SBATCH --mail-user=%s
#SBATCH --mail-type=ALL
#
echo $SLURM_JOB_ID
echo $SLURM_JOB_NAME
echo $SLURM_ARRAY_TASK_ID
cd %s
#
%s"""

SlurmJobStringTemplateWithoutEmail="""#!/bin/bash
#SBATCH --job-name=%s
#SBATCH --nodes=%s
#SBATCH --ntasks=%s
#SBATCH --output=%s
#SBATCH --error=%s
#SBATCH --time=%s
#SBATCH --mem=%s
#
echo $SLURM_JOB_ID
echo $SLURM_JOB_NAME
echo $SLURM_ARRAY_TASK_ID
cd %s
#
%s"""

PBSJobStringTemplateWithEmail="""#!/bin/bash
#PBS -q sstar
#PBS -A p003_swin
#PBS -N %s
#PBS -l walltime=%s
#PBS -l %s
#PBS -o %s
#PBS -e %s
#PBS -W depend=%s
echo $PBS_O_WORKDIR
cd $PBS_O_WORKDIR
%s"""

PBSJobStringTemplateWithoutEmail="""#!/bin/bash
#PBS -q sstar
#PBS -A p003_swin
#PBS -N %s
#PBS -l walltime=%s
#PBS -l %s
#PBS -o %s
#PBS -e %s
#PBS -W depend=%s
echo $PBS_O_WORKDIR
cd $PBS_O_WORKDIR
%s"""

#####################################################

# Check if we are using python 3
python_version = sys.version_info[0]
print("python_version =", python_version)

#####################################################

def runBashCommand(bashCommand, verbose=True):
	"""
	Run bash command

	Parameters
	-----------
	bashCommand : str
		Bash command to be run
	verbose : boolean
		Whether to echo command being run

	Returns
	--------

	"""
	if(verbose):
		print(bashCommand)
	os.system(bashCommand)
	return

def createBashExecutableFile(directory, command, executable_name, venvActivatePath=None, run_command=False, verbose=False):
	"""
	Create an executable Bash file to be run
	
	Parameters
	------------
	directory : path
		Path to where to create the executable file
	command : str
		The command to be run
	executable_name : str
		The name of executable to create
	venvActivatePath : boolean or None
		Whether to source a virtual environment or not
	verbose : boolean
		Whether or not to run bash command

	Returns
	---------
	
	"""
	this_bash_executable = os.path.join(directory, executable_name + '.bash')
	exeFile = open(this_bash_executable, 'w')
	exeFile.write('#!/bin/bash\n') 
	#exeFile.write('#!/usr/bin/env bash\n') ##!/usr/bin/env bash

	if not directory == None:
		exeFile.write('cd ' + directory + '\n')

	if not venvActivatePath == None:
		exeFile.write('source ' + venvActivatePath + '\n')

	if not command == None:
		exeFile.write(command + '\n')

	exeFile.close()

	#-- Make post-processing executable
	bashCommand = 'chmod 744 '+ this_bash_executable
	print(bashCommand)
	os.system(bashCommand)
	
	if run_command:
		bashCommand = 'bash ' + this_bash_executable
		runBashCommand(bashCommand, verbose)

	return

def setupGridListRunDirectories(gridDictionary=None, nGridPoints=1, isGridRun=False, isListRun=False, dirNames=None, pklName="", gridListSubDir="", verbose=True):
	"""
	Set up structure of folders for grid/list runs

	Parameters
	-----------
	gridDictionary : dictionary
		Dictionary containing options for grid/list runs
	nGridPoints : int
		Number of grid/list runs
	isGridRun : boolean
		Whether this is a grid run
	isListRun : boolean
		Whether this is a list run
	dirNames : list
		List of names of run directories
	pklName : str
		Name of pickled object
	gridListSubDir : str
		Name of subdirectory containing grid/list output folders

	Returns
	--------
	"""
	if(gridDictionary is not None):
		if(isGridRun or isListRun):
			for dirName in dirNames:

				if(verbose):
					print(dirName)

				with open(os.path.join(dirName,pklName),'wb') as pg:
					pickle.dump(gridDictionary,pg)
				
				for j in range(nGridPoints):
					path = dirName + gridListSubDir + 'output-' + str(j)

					if(verbose):
						print(path)

					os.makedirs(path)

	return

def writeCondorSubmitFile(condor_submit_path="", send_email=False, user_email=""):
	"""
	Write a Condor submit file 
	
	Parameters
	------------
	condor_submit_path : str
		Path to condor submit file with no extension
    send_email : bool
        Whether to send an email to the user
	user_email : str
		Users email
	
	Returns
	--------

	"""
	submitFile = open(condor_submit_path + '.submit','w')
	if send_email:
		submitFile.write('notify_user=' + user_email + '\n')
		submitFile.write('notification=Always\n')
	submitFile.write('executable='+condor_submit_path+'.bash\n')
	submitFile.write('getenv=True\n')
	submitFile.write('universe=vanilla\n')
	submitFile.write('log='+condor_submit_path+'-log.txt\n')
	submitFile.write('error='+condor_submit_path+'-error.txt\n')
	submitFile.write('output='+condor_submit_path+'-output.txt\n')
	#submitFile.write('requirements = machine != "cyclone.sr.bham.ac.uk"\n')
	#submitFile.write('request_cpus = 2' + '\n')
	submitFile.write('queue')
	submitFile.close()
	return

def SubmitPBSJob(dirName, job_name, bashScript, dependencies, walltime, user_email):
	"""
	Submit a COMPAS job to the job scheduler PBS

	Parameters
	----------
	dirName : str
		Name of directory where run is performed
	job_name : str
		Name of this job. By default is COMPAS_i where i is in [0, nBatches]
	bashScript : str
		Name of the bash script you want PBS to run
	dependencies : str
		Dependencies of this job 
	walltime : str
		As a string, in the format HH:MM:SS, the requested walltime
	user email : str
		Users email address

	Returns
	-------
	job_id : str
		Job id as a string
	"""
	# Open a pipe to the qsub command.
	proc = Popen('qsub', shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE, close_fds=True)

	# Customize your options here
	processors = "nodes=1:ppn=1" # hardcoded for now

	logfile = os.path.join(dirName, job_name) + ".log"
	outfile = os.path.join(dirName, job_name) + ".out"
	errfile = os.path.join(dirName, job_name) + ".err"

	command = os.path.join(dirName, bashScript) + " > " + logfile

	job_string = PBSJobStringTemplate % (job_name, walltime, processors, outfile, errfile, dependencies, command)

	# Send job_string to qsub
	if (sys.version_info > (3, 0)):
		proc.stdin.write(job_string.encode('utf-8'))
	else:
		proc.stdin.write(job_string)
	
	out, err = proc.communicate()

	job_id = out.split()[-1]

	if(verbose):
	    # Print your job and the system response to the screen as it's submitted
		print(job_string)
		print("out = ", out)

	time.sleep(0.1)

	return job_id


def GenerateSlurmJobString(job_name, number_of_nodes, number_of_cores, outfile, errfile, walltime, memory, send_email=False, user_email="", run_directory="./", command="./"):
	"""
	Parameters
	-----------
	job_name: str
		Name of the job
	number_of_nodes : int
		Number of nodes to request
	number_of_cores : int
		Number of cores to request
	outfile : str
		Path to file to print output
	errfile : str
		Path to file to print errors
	walltime : str
		Walltime (HH:MM:SS)
	memory : str
		Memory to request (in MB)
	send_email : bool
		Whether to send email or not
	user_email : str
		User's email address
	run_directory : str
		Where the job will run
	command : str
		Job to be run

	Returns
	--------
	job_string : str
		The string which will be submitted to sbatch
	"""
	if send_email:
		job_string = SlurmJobStringTemplateWithEmail % (job_name, number_of_nodes, number_of_cores, outfile, errfile, walltime, memory, user_email, run_directory, command)
	else:
		job_string = SlurmJobStringTemplateWithoutEmail % (job_name, number_of_nodes, number_of_cores, outfile, errfile, walltime, memory, run_directory, command)
	
	return job_string

def write_dag_retry_line(dag, job_name, nRetries):
	"""
	Parameters
	-----------
	dag : file
		DAG file to write to
	job_name : str
		Name of job to retry
	nRetries : int
		Number of times to retry job
		
	Returns
	--------

	"""
	dag.write('RETRY ' + job_name + ' ' + str(nRetries) + '\n')
	dag.write('\n')
	return 

def sbatchCommandDependency(dependencyString='afterok', dependencyID=''):
	"""
	Write an sbatch command for a job with a dependency

	Parameters
	-----------
	dependencyString : str
		Type of job dependency. Default = 'afterok'
	dependencyID : str 
		ID of job to be dependent on
	
	Returns
	---------
	"""
	# Check if we are using python 3
	if python_version >= 3:
		dependencyID = dependencyID.decode("utf-8")
		print("Dependency ID =", dependencyID)

	sbatchCommand = 'sbatch --dependency=afterok:' + str(dependencyID)

	return sbatchCommand
