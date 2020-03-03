#
#  compas_hpc_input.py
#
#  This is the script that should be modified by the user to give your inputs 
#  to compas_hpc.py
#  
#  The names of the variables in this script should not be changed as they
#  are expected by compas_hpc.py
#  
#  Once you have edited this script, run
#  
#  python compas_hpc.py
#

import numpy as np
import os
import time
from subprocess import Popen, PIPE
import sys
import pickle

################################################
#
#	The user should modify the following options
#
#################################################

nBatches = 3 #-- Choose how many batches to split the run into
maxRandomSeed = int(2**32-1)
venvActivatePath = None # path to source or None

send_email = False # Whether you want to recieve an email when your jobs are done
user_email = 'your.email.address@your.domain'

#-- Path where the output will go
# Does not have to exist, compas_hpc will attempt to create it
rootOutputDir = '/fred/oz101/user_name/folder_name'

#-- on G2 this should be something like
# /lustre/projects/p027/user_name/foler_name

#-- on OzSTAR this should be something like
# /fred/oz003/user_name/folder_name

#-- on tsunami this should be something like
# /home/user_name/folder_name

#-- Whether to use your own list of nBatches random seeds or generate a new list
generate_random_seeds = True
seedsFile = "/path/to/seeds/file"

#-- Which cluster to use. Current possible clusters are tsunami, g2, ozstar, helios and local
cluster = 'ozstar'

#-- Request walltime in HH:MM:SS -- used on ozstar and g2
walltime = "1:00:00"

#-- Request memory in MB -- used on ozstar and g2
memory='4000'

#-- Set maximum number of jobs to run at any one time
maxNumJobsRun = '100'

#
#  Example of how to set up the grid dictionary. 
#  Modify to what you want
#
#  To see available options to change do 
#  $COMPAS_ROOT_DIR/COMPAS/COMPAS --help
#
#  Example of a grid of common envelope alphas
#  gridDictionary = {}
#  gridDictionary['--common-envelope-alpha'] = np.linspace(0.,2.,10)
#  
#  Example of a metallicity grid
#  gridDictionary = {}
#  n_metallicities = 50
#  gridDictionary['--metallicity'] = np.logspace(log_metallicity_lower_limit,log_metallicity_upper_limit,n_metallicities)
#
#  These are the SSE (Hurley et al 2000) limits for the 
#  range of metallicities
#
metallicity_lower_limit = 1E-4
metallicity_upper_limit = 3E-2

log_metallicity_lower_limit = np.log10(metallicity_lower_limit)
log_metallicity_upper_limit = np.log10(metallicity_upper_limit)

common_envelope_alpha_lower_limit = 0.01
common_envelope_alpha_upper_limit = 2.0

sigma_kick_black_hole_lower_limit = 0.0
sigma_kick_black_hole_upper_limit = 400.0

flbv_lower_limit = 0.0
flbv_upper_limit = 10.0

ranges = np.array([[log_metallicity_lower_limit, log_metallicity_upper_limit],
                   [common_envelope_alpha_lower_limit, common_envelope_alpha_upper_limit],
                   [sigma_kick_black_hole_lower_limit, sigma_kick_black_hole_upper_limit],
                   [flbv_lower_limit, flbv_upper_limit]])

gridDictionary = None
