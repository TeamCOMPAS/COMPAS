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
# ## Here we show the basic h5 file syntax for how to:
#     
# ### 1 - Load the file.
# ### 2 - Inspect the data.
# ### 3 - Close the file.
#
# ## For this and the following sections, you will need to have the python package h5py installed. 
# ## Running `pip install h5py` on the command line is one way to do this.

# +
import os, sys    # for handling paths
import h5py as h5  #for handling data format

# Import local directories and scripts
compasRootDir = os.environ['COMPAS_ROOT_DIR'] 
pythonScriptDir = compasRootDir + 'postProcessing/PythonScripts'
sys.path.append(pythonScriptDir)

from printCompasDetails import printCompasDetails
# -

# # 1 - Load the h5 file

# +
# Set the appropriate path to the data file

outputDir = compasRootDir + '/postProcessing/Tutorial/'
pathToH5 = outputDir + 'COMPAS_Output/COMPAS_Output.h5' 
Data  = h5.File(pathToH5)
# -

# # 2 - Inspect the data

list(Data.keys())

# ### The list of output here depends on your simulation. Here you may see 'BSE_System_Parameters', which collects all of the information of the systems at the beginning of the evolution, as well as other major events such as 'BSE_RLOF' for all mass transfer events, 'BSE_Common_Envelopes' for all unstable mass transfer events, 'BSE_Supernovae' if you have sufficiently massive stars to undergo a supernova event, and 'BSE_Double_Compact_Objects' if you happen to form any intact binaries composed of either neutron stars or black holes (though this is quite rare and may require a large simulation).
#
# ### Running Single Stellar Evolution (SSE) will produce generally different outputs to Binary Stellar Evolution (BSE), but both will create a Run_Details file, which captures all of the input settings.

# ### To show all the parameters in a given file

# It is a pain to write the entire group each time so we define shorthands
SPs = Data['BSE_System_Parameters']
list(SPs.keys())

# ### To find the unit of a single parameter

print(SPs['Mass@ZAMS(1)'].attrs['units']) # attrs refers to attributes

# ### To view all of the contents of a given HDF5 group, use the printCompasDetails function

printCompasDetails(SPs) # Note - the output of this is a pandas dataframe

# ### To access the values of a column

#Giving me the actual array
mZams1 = SPs['Mass@ZAMS(1)'][()]
print(mZams1.shape)                   # number of systems in this file
print(mZams1[:10])                    # the values of the first 10 entries

# # 3 - Closing the Data
#
# Accessing a single h5file from multiple scripts is not always possible.
# With notebooks, sometimes closing the notebook is not enough to have it 
# close the h5data. Therefore we recommend to close the h5file explicitly after
# the calculations are done
#

Data.close()


