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
# In order to calculate the formation rate of a type of double compact object at a given redshift, we need to know:
#
# 1 - The grid of metallicities we assume for the integral
#
# 2 - The amount of solar mass evolved per metallicity per system $\frac{d}{dMsun}(Z)$
#
# 3 - The type of double compact object (DCO) we are interested in 
#
# 4 - The metallicity  $Z$ at which each DCO formed and the
#     delay time (time from formation till merger) $t_c$ for each DCO
#
#
# Given a time at which it merges we can then calculate the time at which it formed to recover the MSSFR ($\frac{dMsun}{dt})$ 
#
#
# In addition we need to know the component masses of the system in order to calculate any selection effects
#
#
# The ClassCOMPAS.py is there to store the information that we need such that we can access it quickly when calculating the rates.

# # Paths

import os
pathNoteBook     = os.getcwd()
pathClassCOMPAS  = pathNoteBook + '/PythonScripts/'
pathData        = '/Users/lieke/surfdrive/Documents/test_CI/COMPAS_Output/'#"/home/cneijssel/Desktop/Test/"

# # Imports

import numpy as np
import sys
sys.path.append(pathClassCOMPAS)
import ClassCOMPAS

# # Creating an instance of the COMPAS data class

# To create an instance of the COMPAS output class, we need to specify the following attributes:
#
#     path      = None
#     fileName  = 'COMPAS_output.h5'
#     
#     lazyData  = True
#
#     Mlower    = None
#     Mupper    = None
#     binaryFraction=None

# Path: 
#     
#     path to the h5-data. 
#     
# fileName:
#
#     name of the h5 data
#
# lazyData:
#
#     means we store additional info like the mass ratios, chirpmasses of 
#     each DCO system. In principle this could be done externally but this
#     is slightly easier when plotting/combingin info later on, 
#     but it does use more memory
#     
#     
# Mlower:
#
#     lower limit used for M1 in the pythonSubmit of the simulation. 
#     Needed to recover `true' amount of Msun evolved (see step 4)
#     
# Mupper:
#
#     upper limit used for M1 in the pythonSubmit of the simulation. 
#     Needed to recover `true' amount of Msun evolved (see step 4)
#
# binaryFraction:
#
#
#     assumed fraction of stars in binaries. 
#     Needed to recover `true' amount of Msun evolved (see step 4)
#

#I assume all the defaults and just set the path
COMPASData = ClassCOMPAS.COMPASData(path=pathData)

# The output are  reminders which will be explained in next steps

# # 1 Total mass Evolved
#
# In the COMPAS simulation we often only evolve massive stars.
# This means that the total mass in our simulation does not represent the total mass evolved
# in all stars. Here we recover an estimate of what that total mass is using the
# lower and upper mass for the primary from your python submit and assuming a binary fraction.
#
# The code will then check your data and test per metallicity how much mass is evolved.
# It assumes the metallicities are subject to the same pytonSubmit, but maybe due to sampling had a different number of systems. It also recovers from the data what metallicity grid is used

# +
COMPASData.Mlower = 15
COMPASData.Mupper = 150
COMPASData.binaryFraction =0.7

COMPASData.setGridAndMassEvolved()
# -

# # 1.1 The grid of metallicities we assume for the integral

# By default the ClassCOMPAS will automatically try to recover the metallicity grid from the data/
# It assumes that metallicities of all the systems in the h5-data represent the assumed metallicity grid
# for the calculation. 
#
#         metallicities =Data['SystemParameters']['Metallicity@ZAMS_1'][()]
#         self.metallicityGrid     = np.unique(metallicities)

# In principle you could instead overwrite this with your own metallicity grid. However
# remember to reassign the metallicities of each DCO and the amount of solar mass evolved per metallicity. However, we leave it at that for now. You can acces the grid by printing

print(COMPASData.metallicityGrid)

# # 2 The amount of solar mass evolved per system per Z

# Again by default the ClassCOMPAS will automatically recover the amount of
# `true'  amount of solar mass evolved per system using the totalMassEvolvedPerZ script and
# by reading the total mass per system in the simulation. This recovers an amount of solar mass per metallicity of the metallicity grid (units Msun).

print(COMPASData.totalMassEvolvedPerZ)

# # 3 The select type of DCO to calculate the merger rate for

# To recover the metallicities delaytimes and other parameters of the DCOs you are interested in we use a boolean mask. The boolean mask, which has the same length as the DCO h5 group, selects the systems we want to include in the calculation. 
#
# You could set your own mask using any combination you want by
#
#     maskDCO = some criteria you like on the h5 data
#     COMPASData.DCOmask = maskDCO
#     
# However, usually we are interested in a specific group of merging DCOs assuming a type of physics. The setCOMPASDCOmask()  allows you to quickly set the mask without doing the slicing yourself and takes the following arguments (default then all options)
#
#     argument         = default   / options
#     types            =  'BBH'    / 'BBH', 'BNS', 'BHNS', 'All' (BBH, BNS, or BHNS)
#     withinHubbleTime =  True     / True, False 
#     pessimistic      =  True     / True, False
#     noRLOFafterCEE   =  True     / True, False

# type: Type of double compact object (DCO) to mask for. Can also take the argument 'All' to mask for BBHs, BNSs, and BHNSs.
#
# withinHubbleTime: If True, only use DCOs that merge within a Hubble time.
#
# pessimistic: If True, mask out DCOs that have formed through a common-envelope event involving a Hertzsprung-gap donor. 
#
# noRLOFafterCEE: If True, mask out DCOs that have at some point experienced RLOF immediately after a common-envelope event. 

COMPASData.setCOMPASDCOmask(types='BBH', pessimistic=True)

#Check if we have any system meeting the criteria
print('nr systems =%s ' %(np.sum(COMPASData.DCOmask)))

# # 4 - Get the metallicities and delay times
#
# using the DCO mask defined in step 3 the class can now get the parameters of interest
# for each mergingg DCO
#

COMPASData.setCOMPASData()

# Now the data is set and you are ready to go

#
# # For different Data

# If you have your own simulation which is different then the COMPAS data, or you want
# to test a toy model, then you can still use the set of pipelines for the cosmic integration.
#
# The only thing you need to do is construct your own arrays.
#
# Create an instance of the clasCOMPAS without a path

MockData = ClassCOMPAS.COMPASData(path=None)

# Then manually set each array for
#
#     #grid for integral
#     MockData.metallicityGrid
#     MockData.totalMassEvolvedPerZ #same length array as grid
#     
#     #Metallicity of each system corresponding to a grid-point
#     MockData.metallicitySystems
#     MockData.delayTimes  #Myr
#     MockData.mass1       #Msun
#     MockData.mass2       #Msun    
#     #All four arrays are same length since they correspond to number of systems

# All other pipelines just read arrays and are independent of the data

# So you could sample from the IMF for M1 and M2, 
# create a grid in metallicities and uniformly sample from the grid and then sample
# from a $t^{-1}$ delay time distribution, just note the assumed units :)
#
# If you then set the totalMassEvolvedperZ to an array of ones, then you can at least create
# predictions for the shape of the merger rate distributions. Since this array is only a normalization affecting the absolute rates.


