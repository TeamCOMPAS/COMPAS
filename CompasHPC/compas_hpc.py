##############################################
#
# compas_hpc.py
#
# Pipeline to run COMPAS on a High Performance Computer (HPC)
# 
##############################################

import numpy as np
import os
import time
from subprocess import Popen, PIPE
import sys
import pickle
import compas_hpc_functions as chpc

###############################################
#-- Get environment variable
compasRootDir = os.environ.get('COMPAS_ROOT_DIR')

compasHPCDir = os.path.join(compasRootDir, 'CompasHPC') #compasRootDir + 'CompasHPC'  #-- Path to compasHPC

AISRootDir = os.path.join(compasRootDir, "AdaptiveImportanceSampling/compasHPC")

#-- Set up commands
masterFolderDir = os.path.join(compasHPCDir, "masterFolder") #-- Path to the folder masterFolder which this script will copy multiple times

#-- get a copy of the pythonSubmit
sys.path.append(masterFolderDir)
import pythonSubmit as ps
programOptions = ps.pythonProgramOptions()

###############################################
#-- Get command line arguments
import argparse
parser = argparse.ArgumentParser(description='''
								Pipeline to run COMPAS on a High Performance Computer (HPC)
								''',
								epilog='''
								''') #prog='compas_hpc',
# parser.add_argument("--masterFolderDir", help="Directory to masterFolder") # example of adding an option
args = parser.parse_args()

###############################################

#-- Load the user specified settings
#   Change this if you change the name of the input
#   Edit this to be a class?
from compas_hpc_input import *

###############################################

# Check if we are using python 3
python_version = sys.version_info[0]
print("python_version =", python_version)

###############################################
#
# Check number of batches is sensible
#
###############################################

if(nBatches <= 0):
	raise ValueError("We don't think you really meant to request nBatches = " + str(nBatches) + ". Please use nBatches >= 1")

if(np.logical_not(type(nBatches) == int)):
	raise TypeError("You must provide an integer number of batches. You requested " + str(nBatches) + ". Maybe try int(number).")

###############################################

#-- Make the output directory
if(np.logical_not(os.path.isdir(rootOutputDir))):
	os.makedirs(rootOutputDir) #os.mkdir(rootOutputDir)

###############################################

# What options has the used specified for this run?
isGridRun = programOptions.hyperparameterGrid
isListRun = programOptions.hyperparameterList
isAISRun  = programOptions.AIS_exploratory_phase

if(isAISRun and (isGridRun or isListRun)):
	raise ValueError("CompasHPC is not tested using both grid/list run and AIS. Contact Simon if you want to do this.")

if(isGridRun and isListRun):
	raise ValueError("You can't have both a list and a grid!")

if(isGridRun):
	nGridPoints = 1
	for key in gridDictionary.keys():
	    nGridPoints *= len(gridDictionary[key])
	pklName='pickledGrid.pkl'
	gridListSubDir='/gridOutputs/'

if(isListRun):
	lengths = []
	for key in gridDictionary.keys():
	    lengths.append(len(gridDictionary[key]))

	nGridPoints = lengths[0]

	if np.any(np.diff(lengths)):
		raise ValueError("all hyperparameters need the same number of points to make a list!")
	
	pklName='pickledList.pkl'
	gridListSubDir='/listOutputs/'

if(not isListRun and not isGridRun):
	nGridPoints = 0
	pklName='pickledList.pkl'
	gridListSubDir='/listOutputs/'

################################################

################################################

#-- Generate random seeds
if(isAISRun):
	if generate_random_seeds:
		seeds = np.random.randint(2,maxRandomSeed,2*nBatches)
	else:
		seeds = np.loadtxt(seedsFile)
	assert len(seeds) == 2*nBatches

else:
	if generate_random_seeds:
		seeds = np.random.randint(2,maxRandomSeed,nBatches) 
	else:
		seeds = np.loadtxt(seedsFile)
	assert len(seeds) == nBatches

ppFileHead   = 'postProcessing'
ppPythonFile = ppFileHead + '.py'
ppName       = "COMPAS_PP"
ppWalltime   = "10:00:00"
ppMemory     = "8000"

#-- Names of AIS executables
AISstep2Name         = "AISstep2"
AISStep2PythonName   = "ImportanceSampling-step2.py"
AISStep2PythonFile   = os.path.join(AISRootDir, AISStep2PythonName)

AISstep3Name         = "AISstep3"
AISStep3PythonName   = "ImportanceSampling-step3.py"
AISStep3PythonFile   = os.path.join(AISRootDir, AISStep3PythonName)

AISCombineFileHead   = 'AISCombine'
AISCombinePythonName = 'combineHDF5FilesScript.py'
AISCombinePythonFile = os.path.join(AISRootDir, AISCombinePythonName)
AISCombineName       = 'AISCombine'
AISCombineWalltime   = "10:00:00"
AISCombineMemory     = "8000"

AISPostProcessingName       = "AISpostProcessing"
AISPostProcessingPythonName = "ImportanceSampling-postProcessing.py"
AISPostProcessingPythonFile = os.path.join(AISRootDir, AISPostProcessingPythonName)

if isAISRun:
	#-- Set up folders for Exploratory and Sampling phases, and directory for combined results
	AISExploratoryDir = os.path.join(rootOutputDir, "AIS_exploratory")
	AISSamplingDir    = os.path.join(rootOutputDir, "AIS_sampling")
	AISCombinedDir    = os.path.join(rootOutputDir, "AIS_combined")

	os.makedirs(AISExploratoryDir)
	os.makedirs(AISSamplingDir)
	os.makedirs(AISCombinedDir)

	bashCommand = 'cp ' + os.path.join(compasHPCDir, ppPythonFile) + " " + os.path.join(AISExploratoryDir, '.')
	chpc.runBashCommand(bashCommand, verbose=True)

	bashCommand = 'cp ' + os.path.join(compasHPCDir, ppPythonFile) + " " + os.path.join(AISSamplingDir, '.')
	chpc.runBashCommand(bashCommand, verbose=True)

else:
	bashCommand = 'cp ' + os.path.join(compasHPCDir, ppPythonFile) + " " + os.path.join(rootOutputDir, '.')
	chpc.runBashCommand(bashCommand, verbose=True)

bashCommand = 'cp ' + os.path.join(compasHPCDir, "plottingRoutines.py") + " " + os.path.join(rootOutputDir, '.')
chpc.runBashCommand(bashCommand, verbose=True)

bashCommand = 'cp ' + os.path.join(compasHPCDir, "template.html") + " " + os.path.join(rootOutputDir, '.')
chpc.runBashCommand(bashCommand, verbose=True)

#-- Make a copy of this script in the output directory
bashCommand = 'cp ' + os.path.join(compasHPCDir, "compas_hpc.py") + " " + os.path.join(rootOutputDir, '.')
chpc.runBashCommand(bashCommand, verbose=True)

#-- Make a copy of the input file in the output directory
bashCommand = 'cp ' + os.path.join(compasHPCDir, "compas_hpc_input.py") + " " + os.path.join(rootOutputDir, '.')
chpc.runBashCommand(bashCommand, verbose=True)

if(isAISRun):
	ppDir  = AISExploratoryDir
	COMPASOutuptFileName = "COMPASOutput.h5"
	COMPASOutput1 = os.path.join(AISExploratoryDir, COMPASOutuptFileName)
	COMPASOutput2 = os.path.join(AISSamplingDir, COMPASOutuptFileName)
	COMPASOutputCombined = os.path.join(AISCombinedDir, COMPASOutuptFileName)

	AISCombineCommand  = 'python ' + AISCombinePythonFile + " --input-file1 " + COMPASOutput1 + " --input-file2 " + COMPASOutput2 + " --output-file " + COMPASOutputCombined

else:
	ppDir  = rootOutputDir

ppPath = os.path.join(ppDir, ppFileHead)

newMasterFolderDir = os.path.join(rootOutputDir,'masterFolder')

dirNames = []
exploratory_dirNames = []
sampling_dirNames = []

#-- Number of times to retry failed jobs (default = 3)
number_of_times_to_retry_failed_jobs = 3

slurm_cluster_list  = ['ozstar', 'slurm', 'helios']
condor_cluster_list = ['tsunami', 'condor']
pbs_cluster_list    = ['G2', 'PBS']

#####################################################
# 	Run based on chosen cluster
#####################################################

if cluster in condor_cluster_list:

	#-- Set up bash file for post-processing
	ppLogFile = os.path.join(ppDir, ppFileHead + '.log')
	ppCommand = 'python ' + os.path.join(ppDir, ppPythonFile) + " --masterFolderDir " + newMasterFolderDir + ' > ' + ppLogFile + '\n'
	chpc.createBashExecutableFile(directory=ppDir, command=ppCommand, executable_name=ppFileHead, venvActivatePath=venvActivatePath, run_command=False, verbose=True)

	#-- Set up post-processing condor submit file
	chpc.writeCondorSubmitFile(condor_submit_path=ppPath, send_email=send_email, user_email=user_email)

	if(isAISRun):

		ppStep3Command = 'python ' + os.path.join(AISSamplingDir, ppPythonFile) + " --masterFolderDir " + newMasterFolderDir + '\n'
		chpc.createBashExecutableFile(directory=AISSamplingDir, command=ppStep3Command, executable_name=ppFileHead, venvActivatePath=venvActivatePath, run_command=False, verbose=True)

		#-- Set up post-processing condor submit file
		ppStep3Path = os.path.join(AISSamplingDir, ppFileHead)
		chpc.writeCondorSubmitFile(condor_submit_path=ppStep3Path, send_email=send_email, user_email=user_email)

		#-- Set up condor submit file to combine exploratory and sampling results
		AISCombinePath = os.path.join(AISCombinedDir, AISCombineName)
		chpc.writeCondorSubmitFile(condor_submit_path=AISCombinePath, send_email=send_email, user_email=user_email)

		chpc.createBashExecutableFile(directory=AISCombinedDir, command=AISCombineCommand, executable_name=AISCombineName, venvActivatePath=venvActivatePath, run_command=False, verbose=True)
	
	#-- Copy master folder
	bashCommand = 'cp -r ' + masterFolderDir + ' ' + newMasterFolderDir
	chpc.runBashCommand(bashCommand, verbose=True)

	#-- Write condor submit file
	# haven't replaced this because of the log/error/output files
	submitFile = open(os.path.join(newMasterFolderDir, 'compasHPC.submit'), 'w')
	submitFile.write('notify_user=' + user_email + '\n')
	submitFile.write('notification=Always\n')
	submitFile.write('executable=$(directory)/compasHPC.py' + '\n')
	submitFile.write('getenv=True\n')
	submitFile.write('universe=vanilla\n')
	submitFile.write('log=$(directory)/log.txt' + '\n')
	submitFile.write('error=$(directory)/error.txt' + '\n')
	submitFile.write('output=$(directory)/output.txt' + '\n')
	#submitFile.write('requirements = machine != "cyclone.sr.bham.ac.uk"\n')
	submitFile.write('request_cpus = 2' + '\n')
	submitFile.write('queue')
	submitFile.close()

	#-- Make compasHPC executable
	bashCommand = 'chmod +x ' + os.path.join(newMasterFolderDir, 'compasHPC.py')
	chpc.runBashCommand(bashCommand, verbose=True)

	if(isAISRun):

		for i in range(nBatches):
			#-- Directories for exploratory phase
			exploratory_dirName = os.path.join(AISExploratoryDir, 'output' + str(i))
			exploratory_dirNames.append(str(exploratory_dirName))

			bashCommand = 'cp -r ' + newMasterFolderDir + ' ' + exploratory_dirName
			chpc.runBashCommand(bashCommand, verbose=True)

			seedFile = open(exploratory_dirName + '/randomSeed.txt','w')
			seedFile.write(str(seeds[i]))
			seedFile.close()
            
			#-- Directories for sampling phase
			sampling_dirName = os.path.join(AISSamplingDir, 'output' + str(i))
			sampling_dirNames.append(str(sampling_dirName))

			os.makedirs(sampling_dirName)

			seedFile = open(sampling_dirName + '/randomSeed.txt','w')
			seedFile.write(str(seeds[i+nBatches]))
			seedFile.close()

			#-- AIS Step 3 -- sample from AIS distributions
			AISstep3path = os.path.join(sampling_dirName, AISstep3Name)
			AISstep3Python = os.path.join(AISRootDir, AISStep3PythonName)

			#-- Make a copy of the AIS step 3 script in each AIS_sampling output directory
			bashCommand = 'cp ' + AISstep3Python + " " + sampling_dirName + '/.'
			chpc.runBashCommand(bashCommand, verbose=True)

			#-- Make a copy of the pythonSubmit in each AIS_sampling output directory
			bashCommand = 'cp ' + os.path.join(newMasterFolderDir, "pythonSubmit.py") + " " + os.path.join(sampling_dirName, '.')
			chpc.runBashCommand(bashCommand, verbose=True)

			#-- Set up bash file for AIS step 3
			AISstep3Command = 'python ' + os.path.join(sampling_dirName, AISStep3PythonName) + " --masterFolderDir " + newMasterFolderDir + '\n'
			chpc.createBashExecutableFile(directory=sampling_dirName, command=AISstep3Command, executable_name=AISstep3Name, venvActivatePath=venvActivatePath, run_command=False, verbose=True)

			#-- Write condor submit file for AIS step 3
			chpc.writeCondorSubmitFile(condor_submit_path=AISstep3path, user_email=user_email)

	else:
		#-- loop over number of batches
		for i in range(nBatches):

			#-- output directory name
			dirName = os.path.join(rootOutputDir, 'output' + str(i))
			dirNames.append(str(dirName))

			bashCommand = 'cp -r ' + newMasterFolderDir + ' ' + dirName
			chpc.runBashCommand(bashCommand, verbose=True)

			seedFile = open(dirName + '/randomSeed.txt','w')
			seedFile.write(str(seeds[i]))
			seedFile.close()

	#-- set up subdirectories for grid/list runs
	chpc.setupGridListRunDirectories(gridDictionary, nGridPoints, isGridRun, isListRun, dirNames, pklName, gridListSubDir, verbose=True)

	if(isAISRun):

		#-- AIS Step 2 -- generate gaussians
		AISstep2path = os.path.join(newMasterFolderDir, AISstep2Name)
		AISstep2Python = os.path.join(AISRootDir, AISStep2PythonName)

		#-- Make a copy of the AIS step 2 script in the new master folder directory
		bashCommand = 'cp ' + AISstep2Python + " " + newMasterFolderDir + '/.'
		chpc.runBashCommand(bashCommand, verbose=True)

		#-- Set up bash file for AIS step 2
		AISstep2Command = 'python ' + os.path.join(newMasterFolderDir, AISStep2PythonName) + '\n'
		chpc.createBashExecutableFile(directory=newMasterFolderDir, command=AISstep2Command, executable_name=AISstep2Name, venvActivatePath=venvActivatePath, run_command=False, verbose=True)
        
		#-- Write condor submit file for AIS step 2
		chpc.writeCondorSubmitFile(condor_submit_path=AISstep2path, user_email=user_email)

	#-- Set up dag file to run job
	dagFilePath = os.path.join(rootOutputDir, 'compas.dag')

	dag = open(dagFilePath, 'w')

	parentsList = 'PARENT '

	#-- Line in the dag file for post processing
	this_job_name = 'JPP'
	dag.write('JOB ' + this_job_name + ' ' + ppPath + '.submit\n')
	chpc.write_dag_retry_line(dag, this_job_name, number_of_times_to_retry_failed_jobs)

	if(isAISRun):

		#-- Line in the dag file for running AIS step 2 script to generate sampling distributions
		this_job_name = 'JAIS2'
		dag.write('JOB ' + this_job_name + ' ' + os.path.join(newMasterFolderDir, AISstep2Name + ".submit\n"))
		chpc.write_dag_retry_line(dag, this_job_name, number_of_times_to_retry_failed_jobs)

		AISstep3Joblist = ''

		for i in range(nBatches):

			this_job_name = 'J' + str(i)

			#-- Add AIS step 1 exploratory runs
			dag.write('JOB ' + this_job_name + ' ' + exploratory_dirNames[i] + '/compasHPC.submit\n')
			dag.write('VARS ' + this_job_name + ' directory="' + exploratory_dirNames[i] +'"\n')
			chpc.write_dag_retry_line(dag, this_job_name, number_of_times_to_retry_failed_jobs)
			parentsList += this_job_name + ' '

			this_job_name = 'JAIS3' + str(i)

			#-- Add AIS step 3 sampling runs
			dag.write('JOB ' + this_job_name + ' ' + sampling_dirNames[i] + '/AISstep3.submit\n')
			dag.write('VARS ' + this_job_name + ' directory="' + sampling_dirNames[i] + '"\n')
			chpc.write_dag_retry_line(dag, this_job_name, number_of_times_to_retry_failed_jobs)
			AISstep3Joblist += this_job_name + ' '

		#-- Add post processing after step 3 (AIS sampling runs)
		this_job_name = 'JPP2'
		dag.write('JOB ' + this_job_name + ' ' + ppStep3Path + '.submit\n')
		chpc.write_dag_retry_line(dag, this_job_name, number_of_times_to_retry_failed_jobs)

		#-- AIS post-processing
		AISPostProcessingPath   = os.path.join(newMasterFolderDir, AISPostProcessingName)
		AISPostProcessingPython = os.path.join(AISRootDir, AISPostProcessingPythonName)

		#-- Make a copy of the AIS post processing script in the new master folder directory
		bashCommand = 'cp ' + AISPostProcessingPython + " " + newMasterFolderDir + '/.'
		chpc.runBashCommand(bashCommand, verbose=True)
		
		#-- Set up bash file for AIS post processing script
		AISPostProcessingCommand = 'python ' + os.path.join(newMasterFolderDir, AISPostProcessingPythonName) + '\n'
		chpc.createBashExecutableFile(directory=newMasterFolderDir, command=AISPostProcessingCommand, executable_name=AISPostProcessingName, venvActivatePath=venvActivatePath, run_command=False, verbose=True)

		#-- Write condor submit file for AIS post processing
		chpc.writeCondorSubmitFile(condor_submit_path=AISPostProcessingPath, user_email=user_email)
		
		#-- Line in the dag file for updating AIS weights
		this_job_name = 'JAISPP'
		dag.write('JOB ' + this_job_name + ' ' + os.path.join(newMasterFolderDir, AISPostProcessingName + ".submit\n"))
		chpc.write_dag_retry_line(dag, this_job_name, number_of_times_to_retry_failed_jobs)
		
		#-- Line in the dag file for Combining AIS exploratory and sampling output files
		this_job_name = 'JAISPPC'
		dag.write('JOB ' + this_job_name + ' ' + os.path.join(AISCombinedDir, AISCombineName + ".submit\n"))
		chpc.write_dag_retry_line(dag, this_job_name, number_of_times_to_retry_failed_jobs)

	else:
		for i in range(nBatches):
			this_job_name = 'J' + str(i)
			dag.write('JOB ' + this_job_name + ' ' + dirNames[i] + '/compasHPC.submit\n')
			dag.write('VARS ' + this_job_name + ' directory="' + dirNames[i] +'"\n')
			chpc.write_dag_retry_line(dag, this_job_name, number_of_times_to_retry_failed_jobs)
			parentsList += this_job_name + ' '

	#-- Add post processing as a child of main jobs
	dag.write(parentsList + ' CHILD JPP\n')
	
	if(isAISRun):
		#-- Add AIS step 2 as a child of post processing of step 1
		dag.write('PARENT JPP CHILD JAIS2\n')

		#-- Add AIS sampling runs as a child of 
		dag.write('PARENT JAIS2 CHILD ' + AISstep3Joblist + '\n')

		#-- Add post-processing as a child of step 3
		dag.write('PARENT ' + AISstep3Joblist + ' CHILD JPP2' + '\n')

		#-- Add AIS post processing as a child of step 3
		dag.write('PARENT JPP2 CHILD JAISPP\n')

		#-- Add AIS combining exploratory and sampling as child of step 3
		dag.write('PARENT JAISPP CHILD JAISPPC\n')

	dag.close()

	bashCommand = 'condor_submit_dag ' + dagFilePath
	chpc.runBashCommand(bashCommand, verbose=True)

elif cluster in pbs_cluster_list:

	bashCommand = 'cp -r ' + masterFolderDir + ' ' + newMasterFolderDir
	chpc.runBashCommand(bashCommand, verbose=True)

	job_ids = []

	for i in range(nBatches):
		
		#-- Set up a run directory for this run
		dirName = os.path.join(rootOutputDir, 'output' + str(i))
		dirNames.append(str(dirName))
		print(dirName)

		#-- copy the contents of the master folder to the run directory
		bashCommand = 'cp -r ' + newMasterFolderDir + ' ' + dirName
		chpc.runBashCommand(bashCommand, verbose=True)

		#-- Write the seed for this run in a txt file		
		seedFile = open(dirName+'/randomSeed.txt','w')
		seedFile.write(str(seeds[i]))
		seedFile.close()

		if(isGridRun or isListRun):

			print(dirName)

			with open(os.path.join(dirName,pklName),'wb') as pg:
				pickle.dump(gridDictionary,pg)
			
			for j in range(nGridPoints):
				path = dirName + gridListSubDir + 'output-' + str(j)
				print(path)
				os.makedirs(path)
		
		compasHPCcommand = 'python ' + os.path.join(dirName, 'compasHPC.py')
		chpc.createBashExecutableFile(directory=dirName, command=compasHPCcommand, executable_name='compasHPC', venvActivatePath=venvActivatePath, run_command=False, verbose=False)
		
		#-- Go to the run directory
		os.chdir(dirName)
		print(os.getcwd())

		# Customize your options here
		job_name = "COMPAS_%d" % i

		#-- Write the submit script to file
		chpc.createBashExecutableFile(directory=dirName, command=job_string, executable_name='PBSScript', venvActivatePath=venvActivatePath, run_command=False, verbose=False)

		#-- Submit COMPAS job to PBS
		job_id = chpc.SubmitPBSJob(dirName=dirName, job_name=job_name, bashScript='compasHPC.bash', dependencies='', walltime=walltime, user_email=user_email)

		#-- save the job id to a list
		job_ids.append(job_id)

	#-- Run post processing once the rest of the jobs have finished
	if(isAISRun):
		print("Not yet implemented for G2")
	else:
		os.chdir(rootOutputDir)	

	ppCommand = 'python ' + os.path.join(ppDir, ppPythonFile) + " --masterFolderDir " + newMasterFolderDir

	chpc.createBashExecutableFile(directory=rootOutputDir, command=ppCommand, executable_name=ppName, venvActivatePath=venvActivatePath, run_command=False, verbose=False)

	deps = 'afterok:' + ":".join(job_ids)

	#-- Submit COMPAS job to PBS
	ppJobID = chpc.SubmitPBSJob(dirName=ppDir, job_name=ppName, bashScript=ppName + '.bash', dependencies=deps, walltime=walltime, user_email=user_email)

	print("Post-processing job ID = ", ppJobID)

elif cluster in slurm_cluster_list:

	bashCommand = 'cp -r ' + masterFolderDir + ' ' + newMasterFolderDir
	chpc.runBashCommand(bashCommand, verbose=True)

	# Customize your options here
	job_name = "COMPAS"
	number_of_nodes='1'
	number_of_cores = '1'
	#walltime is defined in compas_hpc_input.py
	#memory is defined in compas_hpc_input.py
	#user_email is defined in compas_hpc_input.py

	if(isAISRun):

		run_directory = os.path.join(AISExploratoryDir, 'output${SLURM_ARRAY_TASK_ID}')

		logfile = os.path.join(run_directory, 'COMPAS_${SLURM_ARRAY_TASK_ID}.log')
		errfile = os.path.join(AISExploratoryDir, 'output%a/COMPAS_%a.out')
		outfile = os.path.join(AISExploratoryDir, 'output%a/COMPAS_%a.err')

		command = os.path.join(run_directory, 'compasHPC.bash') + ' > ' + logfile

		job_string = chpc.GenerateSlurmJobString(job_name, number_of_nodes, number_of_cores, outfile, errfile, walltime, memory, send_email=send_email, user_email=user_email, run_directory=run_directory, command=command)

		print(job_string)

		# Save to a file
		sbatchFile = open(AISExploratoryDir+'/compas_hpc.sbatch','w')
		sbatchFile.write(job_string)
		sbatchFile.close()

	else:
		
		run_directory = os.path.join(rootOutputDir, 'output${SLURM_ARRAY_TASK_ID}')
		
		logfile = os.path.join(run_directory, 'COMPAS_${SLURM_ARRAY_TASK_ID}.log')
		errfile = os.path.join(rootOutputDir, 'output%a/COMPAS_%a.out')
		outfile = os.path.join(rootOutputDir, 'output%a/COMPAS_%a.err')

		command = os.path.join(run_directory, 'compasHPC.bash') + ' > ' + logfile

		job_string = chpc.GenerateSlurmJobString(job_name, number_of_nodes, number_of_cores, outfile, errfile, walltime, memory, send_email=send_email, user_email=user_email, run_directory=run_directory, command=command)
		
		print("run directory = ", run_directory)
		print("command = ", command)
		print("Job string:")
		print(job_string)

		# Save to a file
		sbatchFile = open(rootOutputDir+'/compas_hpc.sbatch','w')
		sbatchFile.write(job_string)
		sbatchFile.close()

	if(isAISRun):

		for i in range(nBatches):
			#-- Directories for exploratory phase
			exploratory_dirName = os.path.join(AISExploratoryDir, 'output' + str(i))
			exploratory_dirNames.append(str(exploratory_dirName))

			bashCommand = 'cp -r ' + newMasterFolderDir + ' ' + exploratory_dirName
			chpc.runBashCommand(bashCommand, verbose=True)

			seedFile = open(exploratory_dirName + '/randomSeed.txt','w')
			seedFile.write(str(seeds[i]))
			seedFile.close()

			compasHPCcommand = 'python ' + os.path.join(exploratory_dirName, 'compasHPC.py')

			chpc.createBashExecutableFile(directory=exploratory_dirName, command=compasHPCcommand, executable_name='compasHPC', venvActivatePath=venvActivatePath, run_command=False, verbose=False)

			#-- Directories for sampling phase
			sampling_dirName = os.path.join(AISSamplingDir, 'output' + str(i))
			sampling_dirNames.append(str(sampling_dirName))

			os.makedirs(sampling_dirName)

			seedFile = open(sampling_dirName + '/randomSeed.txt','w')
			seedFile.write(str(seeds[i+nBatches]))
			seedFile.close()

			compasHPCcommand = 'python ' + os.path.join(sampling_dirName, 'compasHPC.py')

			chpc.createBashExecutableFile(directory=sampling_dirName, command=compasHPCcommand, executable_name='compasHPC', venvActivatePath=venvActivatePath, run_command=False, verbose=False)

			#-- AIS Step 3 -- sample from AIS distributions
			AISstep3path = os.path.join(sampling_dirName, AISstep3Name)
			AISstep3Python = os.path.join(AISRootDir, AISStep3PythonName)

			#-- Make a copy of the AIS step 3 script in each AIS_sampling output directory
			bashCommand = 'cp ' + AISstep3Python + " " + os.path.join(sampling_dirName, '.')
			chpc.runBashCommand(bashCommand, verbose=True)

			#-- Make a copy of the pythonSubmit in each AIS_sampling output directory
			bashCommand = 'cp ' + os.path.join(newMasterFolderDir, "pythonSubmit.py") + " " + os.path.join(sampling_dirName, '.')
			chpc.runBashCommand(bashCommand, verbose=True)

			#-- Set up bash file for AIS step 3
			AISstep3Command = 'python ' + os.path.join(sampling_dirName, AISStep3PythonName) + " --masterFolderDir " + newMasterFolderDir + '\n'
			chpc.createBashExecutableFile(directory=sampling_dirName, command=AISstep3Command, executable_name=AISstep3Name, venvActivatePath=venvActivatePath, run_command=False, verbose=True)

		#-- run sbatch array with
		sbatchArrayCommand = 'sbatch --array=0-' + str(nBatches-1) + '%' + str(maxNumJobsRun) + ' ' + os.path.join(AISExploratoryDir, 'compas_hpc.sbatch')
		print(sbatchArrayCommand)

	else: 
		for i in range(nBatches):

			#-- Set up a run directory for this run
			dirName = os.path.join(rootOutputDir, 'output' + str(i))
			dirNames.append(str(dirName))
			print(dirName)

			#-- copy the contents of the master folder to the run directory
			bashCommand = 'cp -r ' + newMasterFolderDir + ' ' + dirName
			chpc.runBashCommand(bashCommand, verbose=True)

			#-- Write the seed for this run in a txt file
			seedFile = open(dirName+'/randomSeed.txt','w')
			seedFile.write(str(seeds[i]))
			seedFile.close()

			compasHPCcommand = 'python ' + os.path.join(dirName, 'compasHPC.py')

			chpc.createBashExecutableFile(directory=dirName, command=compasHPCcommand, executable_name='compasHPC', venvActivatePath=venvActivatePath, run_command=False, verbose=False)

		#setup output folders for grid/list run
		chpc.setupGridListRunDirectories(gridDictionary, nGridPoints, isGridRun, isListRun, dirNames, pklName, gridListSubDir, verbose=True)

		#-- run sbatch array with
		sbatchArrayCommand = 'sbatch --array=0-' + str(nBatches-1) + '%' + str(maxNumJobsRun) + ' ' + os.path.join(rootOutputDir, 'compas_hpc.sbatch')
		print(sbatchArrayCommand)

	proc = Popen(sbatchArrayCommand, shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE, close_fds=True)

	## Send job_string to sbatch
	if (sys.version_info > (3, 0)):
	    proc.stdin.write(sbatchArrayCommand.encode('utf-8'))
	else:
	    proc.stdin.write(sbatchArrayCommand)        

	out, err = proc.communicate()

	time.sleep(0.5)

	print("out = ", out)

	main_job_id = out.split()[-1]

	print("Main job ID = ", main_job_id)

	time.sleep(1.0)
	
	if(isAISRun):

		#-- Run post processing once the exploratory run jobs have finished
		os.chdir(AISExploratoryDir)	

		ppOutputFile = os.path.join(AISExploratoryDir, ppName + '.out')
		ppErrorFile = os.path.join(AISExploratoryDir, ppName + '.err')
		ppLogFile = os.path.join(AISExploratoryDir, ppName + '.log')

		ppCommand = 'python ' + os.path.join(ppDir, ppPythonFile) + " --masterFolderDir " + newMasterFolderDir
		ppCommand += " > " + ppLogFile

		ppString = chpc.GenerateSlurmJobString(ppName, number_of_nodes, number_of_cores, ppOutputFile, ppErrorFile, ppWalltime, ppMemory, user_email=user_email, run_directory=AISExploratoryDir, command=ppCommand)

		sbatchPPCommand = chpc.sbatchCommandDependency(dependencyString='afterok', dependencyID=main_job_id)
		#sbatchPPCommand = 'sbatch --dependency=afterok:' + str(main_job_id)

	else:

		#-- Run post processing once the rest of the jobs have finished
		os.chdir(rootOutputDir)	

		ppOutputFile = os.path.join(rootOutputDir, ppName + '.out')
		ppErrorFile = os.path.join(rootOutputDir, ppName + '.err')
		ppLogFile = os.path.join(rootOutputDir, ppName + '.log')

		ppCommand = 'python ' + os.path.join(ppDir, ppPythonFile) + " --masterFolderDir " + newMasterFolderDir
		ppCommand += " > " + ppLogFile

		ppString = chpc.GenerateSlurmJobString(ppName, number_of_nodes, number_of_cores, ppOutputFile, ppErrorFile, ppWalltime, ppMemory, user_email=user_email, run_directory=rootOutputDir, command=ppCommand)

		sbatchPPCommand = chpc.sbatchCommandDependency(dependencyString='afterok', dependencyID=main_job_id)
		#sbatchPPCommand = 'sbatch --dependency=afterok:' + str(main_job_id)

	print(ppString)

	# Save post-processing job to a file
	sbatchPPFile = open(rootOutputDir+'/postProcessing.sbatch','w')
	sbatchPPFile.write(ppString)
	sbatchPPFile.close()

	print(sbatchPPCommand)

	# Open a pipe to the sbatch command.
	proc = Popen(sbatchPPCommand, shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE, close_fds=True)
	
	# Send job_string to sbatch
	if (sys.version_info > (3, 0)):
		proc.stdin.write(ppString.encode('utf-8'))
	else:
		proc.stdin.write(ppString)
	
	ppOut, err = proc.communicate()
	
	print("ppOut = ", ppOut)
	print("err = ", err)
	
	ppJobID = ppOut.split()[-1]

	print("Post-processing job ID = ", ppJobID)

	if(isAISRun):
		#-- Run AIS step 2 - generate Gaussians
		AISstep2path = os.path.join(newMasterFolderDir, AISstep2Name)
		AISstep2Python = os.path.join(AISRootDir, AISStep2PythonName)
		
		#-- Make a copy of the AIS step 2 script in the new master folder directory
		bashCommand = 'cp ' + AISstep2Python + " " + os.path.join(newMasterFolderDir, '.')
		chpc.runBashCommand(bashCommand, verbose=True)

		AISStep2OutputFile = os.path.join(newMasterFolderDir, AISstep2Name + '.out')
		AISStep2ErrorFile = os.path.join(newMasterFolderDir, AISstep2Name + '.err')
		AISStep2LogFile = os.path.join(newMasterFolderDir, AISstep2Name + '.log')

		#-- Set up bash file for AIS step 2
		AISstep2Command = 'python ' + os.path.join(newMasterFolderDir, AISStep2PythonName) + ' > ' + AISStep2LogFile + '\n'
		chpc.createBashExecutableFile(directory=newMasterFolderDir, command=AISstep2Command, executable_name=AISstep2Name, venvActivatePath=venvActivatePath, run_command=False, verbose=True)

		AISStep2SlurmString = chpc.GenerateSlurmJobString(AISstep2Name, number_of_nodes, number_of_cores, AISStep2OutputFile, AISStep2ErrorFile, ppWalltime, ppMemory, user_email=user_email, run_directory=newMasterFolderDir, command=AISstep2Command)

		sbatchCommand = chpc.sbatchCommandDependency(dependencyString='afterok', dependencyID=ppJobID)
		#sbatchCommand = 'sbatch --dependency=afterok:' + str(ppJobID)

		# Open a pipe to the sbatch command.
		proc = Popen(sbatchCommand, shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE, close_fds=True)
		
		# Send job_string to sbatch
		if (sys.version_info > (3, 0)):
			proc.stdin.write(AISStep2SlurmString.encode('utf-8'))
		else:
			proc.stdin.write(AISStep2SlurmString)
		
		AISStep2Out, AISStep2Err = proc.communicate()

		AISStep2JobID = AISStep2Out.split()[-1]

		print("AIS step 2 job ID = ", AISStep2JobID)

		# #-- run sbatch array with
		run_directory = os.path.join(AISSamplingDir, 'output${SLURM_ARRAY_TASK_ID}')
		
		logfile = os.path.join(run_directory, 'COMPAS_${SLURM_ARRAY_TASK_ID}.log')
		errfile = os.path.join(AISSamplingDir, 'output%a/COMPAS_%a.out')
		outfile = os.path.join(AISSamplingDir, 'output%a/COMPAS_%a.err')

		command = os.path.join(run_directory, 'AISstep3.bash') + ' > ' + logfile

		job_string = chpc.GenerateSlurmJobString(job_name, number_of_nodes, number_of_cores, outfile, errfile, walltime, memory, user_email=user_email, run_directory=run_directory, command=command)
	
		print(job_string)

		# Save to a file
		sbatchFile = open(AISSamplingDir+'/AISstep3.sbatch','w')
		sbatchFile.write(job_string)
		sbatchFile.close()

		if (sys.version_info > (3, 0)):
                    AISStep2JobID = AISStep2JobID.decode('utf-8')
		
		sbatchArrayCommand = 'sbatch --array=0-' + str(nBatches-1) + '%' + str(maxNumJobsRun) + ' ' + '--dependency=afterok:' + str(AISStep2JobID) + ' ' + os.path.join(AISSamplingDir, 'AISstep3.sbatch')
		print(sbatchArrayCommand)

		proc = Popen(sbatchArrayCommand, shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE, close_fds=True)

		## Send job_string to sbatch
		if (sys.version_info > (3, 0)):
		    proc.stdin.write(sbatchArrayCommand.encode('utf-8'))
		else:
		    proc.stdin.write(sbatchArrayCommand)        

		out, err = proc.communicate()

		time.sleep(0.5)

		print("out = ", out)
		print("err = ", err)
		
		AISStep3JobID = out.split()[-1]

		print("AIS step 3 job ID = ", AISStep3JobID)

		time.sleep(1.0)

		#- Run post processing again
		# -- Run post processing once the sampling run jobs have finished
		os.chdir(AISSamplingDir)	

		ppOutputFile = os.path.join(AISSamplingDir, ppName + '.out')
		ppErrorFile  = os.path.join(AISSamplingDir, ppName + '.err')
		ppLogFile    = os.path.join(AISSamplingDir, ppName + '.log')
		ppDir        = AISSamplingDir

		ppCommand = 'python ' + os.path.join(ppDir, ppPythonFile) + " --masterFolderDir " + newMasterFolderDir
		ppCommand += " > " + ppLogFile

		ppString = chpc.GenerateSlurmJobString(ppName, number_of_nodes, number_of_cores, ppOutputFile, ppErrorFile, ppWalltime, ppMemory, user_email=user_email, run_directory=AISSamplingDir, command=ppCommand)

		print(ppString)
		
		sbatchPPCommand = chpc.sbatchCommandDependency(dependencyString='afterok', dependencyID=jobID)

		#sbatchPPCommand = 'sbatch --dependency=afterok:' + str(AISStep3JobID)

		print(sbatchPPCommand)

		# Open a pipe to the sbatch command.
		proc = Popen(sbatchPPCommand, shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE, close_fds=True)
		
		# Send job_string to sbatch
		if (sys.version_info > (3, 0)):
			proc.stdin.write(ppString.encode('utf-8'))
		else:
			proc.stdin.write(ppString)
		
		pp2Out, err = proc.communicate()

		print("pp2Out = ", pp2Out)
		print("err = ", err)

		pp2JobID = pp2Out.split()[-1]

		print("Post-processing 2 job ID = ", pp2JobID)

		#-- Run AIS post processing to update weights
		#-- Make a copy of the AIS post processing script in the new master folder directory
		bashCommand = 'cp ' + AISPostProcessingPythonFile + " " + newMasterFolderDir + '/.'
		chpc.runBashCommand(bashCommand, verbose=True)

		AISUpdateWeightsOutputFile = os.path.join(AISCombinedDir, AISPostProcessingName + '.out')
		AISUpdateWeightsErrorFile  = os.path.join(AISCombinedDir, AISPostProcessingName + '.err')
		AISUpdateWeightsLogFile    = os.path.join(AISCombinedDir, AISPostProcessingName + '.log')
		AISUpdateWeightsCommand    = 'python ' + os.path.join(newMasterFolderDir, AISPostProcessingPythonName) + " > " + AISUpdateWeightsLogFile		

		chpc.createBashExecutableFile(directory=newMasterFolderDir, command=AISUpdateWeightsCommand, executable_name=AISPostProcessingName, venvActivatePath=venvActivatePath, run_command=False, verbose=True)
		
		#-- Run bash file using Slurm
		AISUpdateWeightsSlurmString = chpc.GenerateSlurmJobString(AISPostProcessingName, number_of_nodes, number_of_cores, AISUpdateWeightsOutputFile, AISUpdateWeightsErrorFile, ppWalltime, ppMemory, user_email=user_email, run_directory=newMasterFolderDir, command=AISUpdateWeightsCommand)

		
		sbatchCommand = chpc.sbatchCommandDependency(dependencyString='afterok', dependencyID=pp2JobID)
		#sbatchCommand = 'sbatch --dependency=afterok:' + str(pp2JobID)

		# Open a pipe to the sbatch command.
		proc = Popen(sbatchCommand, shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE, close_fds=True)
		
		# Send job_string to sbatch
		if (sys.version_info > (3, 0)):
			proc.stdin.write(AISUpdateWeightsSlurmString.encode('utf-8'))
		else:
			proc.stdin.write(AISUpdateWeightsSlurmString)
		
		AISUpdateWeightsOut, AISUpdateWeightsErr = proc.communicate()

		AISUpdateWeightsJobID = AISUpdateWeightsOut.split()[-1]

		print("AIS Update Weights job ID = ", AISUpdateWeightsJobID)

		#-- Combine exploratory and sampling phases after AIS postprocessing
		os.chdir(AISCombinedDir)

		AISCombineOutputFile = os.path.join(AISCombinedDir, AISCombineName + '.out')
		AISCombineErrorFile  = os.path.join(AISCombinedDir, AISCombineName + '.err')
		AISCombineLogFile    = os.path.join(AISCombinedDir, AISCombineName + '.log')

		AISCombineCommand += " > " + AISCombineLogFile

		AISCombineString = chpc.GenerateSlurmJobString(AISCombineName, number_of_nodes, number_of_cores, AISCombineOutputFile, AISCombineErrorFile, AISCombineWalltime, AISCombineMemory, send_email, user_email=user_email, run_directory=AISCombinedDir, command=AISCombineCommand)

		sbatchAISCommand = chpc.sbatchCommandDependency(dependencyString='afterok', dependencyID=AISUpdateWeightsJobID)
		#sbatchAISCommand = 'sbatch --dependency=afterok:' + str(AISUpdateWeightsJobID)

		# Open a pipe to the sbatch command.
		proc = Popen(sbatchAISCommand, shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE, close_fds=True)
		
		# Send job_string to sbatch
		if (sys.version_info > (3, 0)):
			proc.stdin.write(AISCombineString.encode('utf-8'))
		else:
			proc.stdin.write(AISCombineString)
		
		AISCombineOut, AISCombineErr = proc.communicate()

		AISCombineJobID = AISCombineOut.split()[-1]

		print("AIS Combine job ID = ", AISCombineJobID)

elif cluster is 'local':

	import multiprocessing as mp

	number_of_local_cores = mp.cpu_count()

	if(nBatches != 1):
		raise ValueError("You have " + str(number_of_local_cores) + " cores locally. However, currently only nBatches = 1 is supported when running locally. Sorry. Maybe Email Simon.")

	#-- Copy master folder
	bashCommand = 'cp -r ' + masterFolderDir + ' ' + newMasterFolderDir
	chpc.runBashCommand(bashCommand, verbose=True)

	if(isAISRun):

		for i in range(nBatches):
			#-- Directories for exploratory phase
			exploratory_dirName = os.path.join(AISExploratoryDir, 'output' + str(i))
			exploratory_dirNames.append(str(exploratory_dirName))

			bashCommand = 'cp -r ' + newMasterFolderDir + ' ' + exploratory_dirName
			chpc.runBashCommand(bashCommand, verbose=True)

			seedFile = open(exploratory_dirName + '/randomSeed.txt','w')
			seedFile.write(str(seeds[i]))
			seedFile.close()

			#-- Directories for sampling phase
			sampling_dirName = os.path.join(AISSamplingDir, 'output' + str(i))
			sampling_dirNames.append(str(sampling_dirName))

			os.makedirs(sampling_dirName)

			seedFile = open(sampling_dirName + '/randomSeed.txt','w')
			seedFile.write(str(seeds[i+nBatches]))
			seedFile.close()

			#-- AIS Step 3 -- sample from AIS distributions
			AISstep3path = os.path.join(sampling_dirName, AISstep3Name)
			AISstep3Python = os.path.join(AISRootDir, AISStep3PythonName)

			#-- Make a copy of the AIS step 3 script in each AIS_sampling output directory
			bashCommand = 'cp ' + AISstep3Python + " " + os.path.join(sampling_dirName, '.')
			chpc.runBashCommand(bashCommand, verbose=True)

			#-- Make a copy of the pythonSubmit in each AIS_sampling output directory
			bashCommand = 'cp ' + os.path.join(newMasterFolderDir, "pythonSubmit.py") + " " + os.path.join(sampling_dirName, '.')
			chpc.runBashCommand(bashCommand, verbose=True)

			AISstep3Command = 'python ' + os.path.join(sampling_dirName, AISStep3PythonName) + " --masterFolderDir " + newMasterFolderDir + '\n'
			chpc.createBashExecutableFile(directory=sampling_dirName, command=AISstep3Command, executable_name=AISstep3Name, venvActivatePath=venvActivatePath, run_command=False, verbose=False)

			#-- change to run directory
			os.chdir(exploratory_dirName)

			#-- run COMPAS
			runCommand = "python " + os.path.join(exploratory_dirName, "pythonSubmit.py")
			chpc.runBashCommand(runCommand, verbose=True)

	else:
		#-- loop over number of batches
		for i in range(nBatches):

			#-- output directory name
			dirName = os.path.join(rootOutputDir, 'output' + str(i))
			dirNames.append(str(dirName))

			bashCommand = 'cp -r ' + newMasterFolderDir + ' ' + dirName
			chpc.runBashCommand(bashCommand, verbose=True)

			seedFile = open(dirName + '/randomSeed.txt','w')
			seedFile.write(str(seeds[i]))
			seedFile.close()

			os.chdir(dirName)

			runCommand = "python " + os.path.join(dirName, "pythonSubmit.py")
			chpc.runBashCommand(runCommand, verbose=True)

	chpc.setupGridListRunDirectories(gridDictionary, nGridPoints, isGridRun, isListRun, dirNames, pklName, gridListSubDir, verbose=True)

	#-- run postProcessing
	os.chdir(ppDir)

	ppCommand = 'python ' + os.path.join(ppDir, ppPythonFile) + " --masterFolderDir " + newMasterFolderDir
	chpc.runBashCommand(ppCommand, verbose=True)

	if(isAISRun):

		#-- AIS Step 2 -- generate gaussians
		AISstep2path = os.path.join(newMasterFolderDir, AISstep2Name)
		AISstep2Python = os.path.join(AISRootDir, AISStep2PythonName)

		#-- Make a copy of the AIS step 2 script in the new master folder directory
		bashCommand = 'cp ' + AISstep2Python + " " + os.path.join(newMasterFolderDir, '.')
		chpc.runBashCommand(bashCommand, verbose=True)

		AISstep2Command = 'python ' + os.path.join(newMasterFolderDir, AISStep2PythonName) + '\n'
		chpc.createBashExecutableFile(directory=newMasterFolderDir, command=AISstep2Command, executable_name=AISstep2Name, venvActivatePath=venvActivatePath, run_command=True, verbose=True)

		#-- Run AIS step 3
		for i in range(nBatches):
			AISstep3path = os.path.join(sampling_dirNames[i], AISstep3Name)
			bashCommand = 'bash ' + AISstep3path + '.bash'
			chpc.runBashCommand(bashCommand, verbose=True)
		
		#-- run postProcessing on Sampling phase
		ppDir = AISSamplingDir
		os.chdir(ppDir)

		ppCommand = 'python ' + os.path.join(ppDir, ppPythonFile) + " --masterFolderDir " + newMasterFolderDir
		chpc.runBashCommand(ppCommand, verbose=True)

		#-- Run AIS post processing to update weights
		AISPostProcessingPath   = os.path.join(newMasterFolderDir, AISPostProcessingName)
		
		#-- Make a copy of the AIS post processing script in the new master folder directory
		bashCommand = 'cp ' + AISPostProcessingPythonFile + " " + newMasterFolderDir + '/.'
		chpc.runBashCommand(bashCommand, verbose=True)
		
		#-- Set up bash file for AIS post processing script
		AISPostProcessingCommand = 'python ' + os.path.join(newMasterFolderDir, AISPostProcessingPythonName) + '\n'
		chpc.createBashExecutableFile(directory=newMasterFolderDir, command=AISPostProcessingCommand, executable_name=AISPostProcessingName, venvActivatePath=venvActivatePath, run_command=True, verbose=True)

		#-- Combine exploratory and sampling phases
		os.chdir(AISCombinedDir)

		AISCombineLogFile = os.path.join(AISCombinedDir, AISCombineName + '.log')

		AISCombineCommand += " > " + AISCombineLogFile
		
		chpc.runBashCommand(AISCombineCommand, verbose=True)
		
		print("Combined AIS exploratory and sampling phases.")

elif cluster is 'bluebear':
	print("Not currently supported")
	exit()

else:
	print("No valid cluster chosen")
	exit()

print("compas_hpc.py completed successfully")
