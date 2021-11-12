# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.12.0
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# # Introduction
#
# ## The COMPAS simulations might be very large in data size while the actual data you need to reproduce your results could be small. Hence it might make sense to reduced the number of files and columns based on some criteria.
#
# ## Here we show how you can reduce your data. The main things you need are:
#     
# ### 1. The seeds you want to have in your data
# ### 2. The files you want in your data
# ### 3. The columns (parameters) you want for each file
#
# ## The python script to do this is found in:
# ### `$COMPAS_ROOT_DIR/postProcessing/PythonScripts/rewrite_H5.py`
# ## Here we just show an example of how to call the script in order to reduce the data.

# +
import os, sys    # for handling paths
import h5py as h5  #for handling data format

# Import local directories and scripts
compasRootDir = os.environ['COMPAS_ROOT_DIR'] 
pythonScriptDir = compasRootDir + 'postProcessing/PythonScripts'
sys.path.append(pythonScriptDir)

import rewrite_H5
# -

# # 1  Load the Data

# +
# Set the appropriate paths to the input and output data files
pathToDataInput = compasRootDir + '/postProcessing/Tutorial/COMPAS_Output/COMPAS_Output.h5' 
pathToDataOutput = compasRootDir + '/postProcessing/Tutorial/COMPAS_Output/COMPAS_Output_reduced.h5' 

Data  = h5.File(pathToDataInput)
print("The main files I have at my disposal are:\n",list(Data.keys()))

# +
# To see the parameter choices in each file, use, e.g:
#print(list(Data['SystemParameters']))
# -

# # 2 Specify which files and columns you want
#
# We use dictionaries to specifically link all the entries.
#
# The filesOfInterest dictionary should contain all files which hold any relevant data. The columnsOfInterest dictionary specifies the parameters in each file that you want to be included in the new output h5. Any filters or masks should be used to determine the seedsOfInterest (on a per file basis), and so do not need to be included in the columnsOfInterest.

# ### Hypothetical Example
#
# Suppose you are studying Double Neutron Star systems, and you want to know the initial parameters of both components. Suppose you are separately curious about the eccentricity of systems following a Supernova that leaves the binary intact, and you want to use the same COMPAS run to save on CPU*hours. 
#
# To be safe, you should probably keep the entire SystemParameters file, which contains all of the initial system settings. 
#
# To get information about only Double Neutron Stars, you will need to create a mask for them from the DoubleCompactObjects file.
#
# Information on post-SN eccentricity and whether or not the system disrupted is found in the Supernovae file. 
#
# You will not need any other files. You will also want to grab the system 'SEED's column from any file, since that is the unique identifier of the binaries. 

# +
# Which files do you want?

# For the files of interest, create 2 dictionary mappings, 
# One mapping file to columns of interest in that file
# And other other mapping file to seeds of interest in that file

#filesOfInterest   = {1:'SystemParameters',\
#                     2:'DoubleCompactObjects',\
#                     3:'Supernovae'}

print(Data.keys())
# Give a list of columns you want, if you want all, say ['All']
columnsOfInterest = {'BSE_System_Parameters':      ['All'],
                     'BSE_Double_Compact_Objects': ['All'],
                     'BSE_Supernovae':             ['SEED', 'Eccentricity']
                    }

# The seedsOfInterest are a little more involved
# -

# # 3 Which seeds do I want per file?

# +
### Do not filter out any systems/seeds from SystemParameters

SPs = Data['BSE_System_Parameters']
seedsSP = SPs['SEED'][()]



### Of all the double compact objects, keep only the DNSs

DCs = Data['BSE_Double_Compact_Objects']
seedsDC =  DCs['SEED'][()]

stellarType1   =  DCs['Stellar_Type(1)'][()]
stellarType2   =  DCs['Stellar_Type(2)'][()]

# Stellar type 13 corresponds to Neutron Stars
maskDNS        =  (stellarType1 == 13) & (stellarType2 == 13)
seedsDNS       =  seedsDC[maskDNS]



### From Supernovae, keep only binaries which stay intact post-SN

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



# Don't forget to close the original h5 data file
Data.close()
# -

# # 4 Call the function which creates the h5 file

rewrite_H5.reduceH5(pathToOld = pathToDataInput, pathToNew = pathToDataOutput,\
                     dictColumns=columnsOfInterest, dictSeeds=seedsOfInterest)

rewrite_H5.printAllColumnsInH5(pathToDataOutput)


