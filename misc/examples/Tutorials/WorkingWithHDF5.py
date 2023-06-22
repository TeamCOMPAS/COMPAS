# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.14.1
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# # Introduction
#
# COMPAS simulations produce all output by default in the form of an [HDF5 file](https://www.hdfgroup.org/solutions/hdf5/), which is a compact and memory efficient, but non-human-readable data file format. In order to interact with these output files, we use the python module `h5py`. 
#
# Within each file are a variety of HDF5 Groups, representing a specific common event in binary evolution, e.g Roche-Lobe Overflow or Supernovae. These are described throughout the post-processing jupyter notebooks.

# ## Material
#
# ### [1) Producing HDF5 output ](#1.-Producing-HDF5-output)
# How to run default COMPAS and produce a simple output file.
#         
# ### [2) Reading HDF5 files ](#2.-Reading-HDF5-files)
# The basics and syntax of loading an HDF5 file.
#         
# ### [3) Rewriting HDF5 files ](#3.-Rewriting-HDF5-files)
# How to rewrite/reduce the HDF5 data.



# ### For the following sections, you will need to have the following python packages installed. 
# ### `numpy, h5py, pandas`
#

# +
import os, sys    # for handling paths
import h5py as h5  #for handling data format

# Import COMPAS root and python script directories
compasRootDir = os.environ['COMPAS_ROOT_DIR'] 
tutorialDir = compasRootDir + '/postProcessing/Tutorial/'

# Import COMPAS specific scripts
sys.path.append(compasRootDir + 'postProcessing/PythonScripts')
import h5rewrite


# -

# # 1. Producing HDF5 output
#
# ## Here we show an example of how to run a basic COMPAS simulation in order to produce output in an HDF5 format. To run COMPAS, you need to set the environment variable $COMPAS_ROOT_DIR to the top level COMPAS directory. Ensure that you have compiled COMPAS, and that there is a COMPAS executable. We will store the output in the Tutorial directory for easy access later on.
#
# ## For simplicity, we will run COMPAS on all default settings except for the number of binaries produced. Here, we will run 1k binaries, which is relatively small. Later, we will look at production of double compact objects (DCOs), but these objects are sufficiently rare that in order to study them, we will need to run a simulation of ~1M binaries, but this can take some time depending on your hardware resources.
#
# ### *Note:* It is required to have some COMPAS output to complete the rest of the post-processing notebooks, but if you already have an output file, you can skip this section.

# To run terminal commands in a jupyter notebook, prepend the command with an exclamation mark
# The curly brackets here are used to input python variables into the bash command
# !COMPAS -n 1000 -o {tutorialDir}

# ### *Note:* This will produce output in the COMPAS_Output directory only the first time it is run. If you run this multiple times, new directories will be created as COMPAS_Output_X. Ensure that you are using the data in the desired directory.

# ## Non-default settings can be seen by running the following in the command line

# !COMPAS --help
print()



# # 2. Reading HDF5 files
#
# ## Here we show the basic h5 file syntax for how to load and close the file, and extract/inspect the data.
#
# While working with this notebook, feel free to change the data file used to see how it affects the output. For the purpose of consistency and ensuring that our examples work properly, the default data file `'COMPAS_Output_tutorial.h5'` has been pre-run using precisely chosen input parameters, specified in `'tutorial_grid.txt'`. 
#
# If the default data file is accidentally deleted, it can be reconstructed by running:
#
# `COMPAS --grid tutorial_grid.txt`

# ## Load the h5 file

# +
# Set the appropriate path to the data file

pathToH5 = tutorialDir + 'COMPAS_Tutorial_Output.h5' 
Data  = h5.File(pathToH5)
# -

# ## Inspect the data

list(Data.keys())

# The list of output here depends on your simulation. Here you may see 'BSE_System_Parameters', which collects all of the information of the systems at the beginning of the evolution, as well as other major events such as 'BSE_RLOF' for all mass transfer events, 'BSE_Common_Envelopes' for all unstable mass transfer events, 'BSE_Supernovae' if you have sufficiently massive stars to undergo a supernova event, and 'BSE_Double_Compact_Objects' if you happen to form any intact binaries composed of either neutron stars or black holes (though this is quite rare and may require a large simulation).
#
# Running Single Stellar Evolution (SSE) will produce generally different outputs to Binary Stellar Evolution (BSE), but both will create a Run_Details file, which captures all of the input settings.

# ### To show all the parameters in a given file

# It is a pain to write the entire group each time so we define shorthands
SPs = Data['BSE_System_Parameters']
list(SPs.keys())

# ### To find the unit of a single parameter

print(SPs['Mass@ZAMS(1)'].attrs['units']) # attrs refers to attributes

# ### To access the values of a column

#Giving me the actual array
mZams1 = SPs['Mass@ZAMS(1)'][()]
print(mZams1.shape)                   # number of systems in this file
print(mZams1[:3])                    # the values of the first 3 entries

# ## Closing the Data
#
# Accessing a single h5file from multiple scripts is not always possible.
# With notebooks, sometimes closing the notebook is not enough to have it 
# close the h5data. Therefore we recommend to close the h5file explicitly after
# the calculations are done
#

Data.close()




# # 3. Rewriting HDF5 files
#
# The COMPAS simulations might be very large in data size while the actual data you need to reproduce your results could be small. Hence it might make sense to reduced the number of files and columns based on some criteria.
#
# ## Here we show how you can reduce your data. You need to decide:
#
# ### - Which data categories (or HDF5 groups) you want to include
# ### - Which seeds and parameters from each data category. 

# ## Load the Data

# +
# Set the appropriate paths to the input and output data files
pathToDataInput = pathToH5 # use the output from the run above
pathToDataOutput = tutorialDir + '/COMPAS_Output/COMPAS_Output_reduced.h5' 

Data  = h5.File(pathToDataInput)
print("The main files I have at my disposal are:\n",list(Data.keys()))
# -

# ## Specify the data categories and parameters and columns
#
# We use dictionaries to specifically link all the entries.
#
# The columnsOfInterest dictionary maps the data categories you wish to include to their appropriate parameters of interest. The seedsOfInterest dictionary maps the data categories to a list of desired seeds. If you wish to exclude some of the systems, that can be accomplished by applying a mask onto the seeds array.

# ### Hypothetical Example
#
# Suppose you are studying Double Neutron Star systems, and you want to know the initial parameters of both components. Suppose you are separately curious about the eccentricity of systems following a Supernova that leaves the binary intact, and you want to use the same COMPAS run to save on CPU*hours. 
#
# To be safe, you should probably keep the entire BSE_System_Parameters file, which contains all of the initial system settings. 
#
# To get information about only Double Neutron Stars, you will need to create a mask for them from the BSE_Double_Compact_Objects file.
#
# Information on post-SN eccentricity and whether or not the system disrupted is found in the BSE_Supernovae file. 
#
# You will not need any other files. You will also want to grab the system 'SEED's column from any file, since that is the unique identifier of the binaries. 

# +
### For each data category, give a list of parameters you want to include
columnsOfInterest = {'BSE_System_Parameters':      ['All'],
                     'BSE_Double_Compact_Objects': ['All'],
                     'BSE_Supernovae':             ['SEED', 'Eccentricity']
                    }

# The seedsOfInterest are a little more involved
# -

# ## Select the desired seeds

# +
### BSE_System_Parameters - keep all seeds
SPs = Data['BSE_System_Parameters']
seedsSP = SPs['SEED'][()]


### BSE_Double_Compact_Objects - keep only Double Neutron Stars
DCs = Data['BSE_Double_Compact_Objects']
seedsDC =  DCs['SEED'][()]
stellarType1   =  DCs['Stellar_Type(1)'][()]
stellarType2   =  DCs['Stellar_Type(2)'][()]

# Stellar type 13 corresponds to Neutron Stars
maskDNS        =  (stellarType1 == 13) & (stellarType2 == 13)
seedsDNS       =  seedsDC[maskDNS]


### BSE_Supernovae - keep only binaries which stay intact post-SN
SNe = Data['BSE_Supernovae']
seedsSN = SNe['SEED']
isUnbound    = SNe['Unbound'][()] == 1
isIntact     = ~isUnbound

seedsIntact  = seedsSN[isIntact]



### Create seedsOfInterest dictionary 
seedsOfInterest   = {'BSE_System_Parameters':      seedsSP,
                     'BSE_Double_Compact_Objects': seedsDNS,
                     'BSE_Supernovae':             seedsIntact
                    }


### Don't forget to close the original h5 data file
Data.close()
# -

# ## Call the function which creates the h5 file

# +
h5rewrite.reduceH5(pathToOld = pathToDataInput, pathToNew = pathToDataOutput,\
                     dictColumns=columnsOfInterest, dictSeeds=seedsOfInterest)

h5rewrite.printAllColumnsInH5(pathToDataOutput)
