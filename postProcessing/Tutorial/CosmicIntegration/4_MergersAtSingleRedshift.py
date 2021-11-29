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

# Magic function to set the backend of matplotlib to the 'inline' backend
# %matplotlib inline

# # Introduction
#
# Here we show how to calculate the merger rate density of systems
# merging at a single redshift z. By now we assume you understand the indiviual pipelines of;
#
# ClassCOMPAS:
#
#     -- (handling the (mock) data needed)
#     
# ClassMSSFR:
#
#     -- defining the model for metallicity specific SFR
#     
# selection_effects:
#
#     -- module to calculate probility detecting a system
#
#
#
# Here we show which additional steps are needed for the calculation.
# All these steps are done by default in the ClassCosmicIntegrator.
# We highlight the steps here outside the function for clarity
# since the ClassCosmicIntegrator merely acts as a giant for loop over multiple redshifts and a way to conveniently store the results

# # Paths

# +
import os

pathNoteBook    = os.getcwd()
pathScripts     = pathNoteBook + '/PythonScripts/'


pathData        = '/Users/lieke/surfdrive/Documents/test_CI/COMPAS_Output/'#"/home/cneijssel/Desktop/Test/"

# -

# # Imports

# +
import numpy as np
import sys
import matplotlib.pyplot as plt

from   astropy.cosmology import WMAP9 as cosmology
from   scipy.optimize import newton
#custom scripts
sys.path.append(pathScripts)
import ClassCOMPAS 
import ClassMSSFR  



# -

# # 1- Set up data and MSSFR model

# +
# Create instance COMPAS data class

COMPAS = ClassCOMPAS.COMPASData(path=pathData, fileName='COMPAS_Output.h5')


# +
#Set the type of DCO of interest and recover their parameters
COMPAS.Mlower = 15
COMPAS.Mupper = 150
COMPAS.binaryFraction =0.7

COMPAS.setGridAndMassEvolved()

COMPAS.setCOMPASDCOmask() #Pessimistic BBHs
COMPAS.setCOMPASData()

# +
# The MSSFR model

#use the metallicityGrid of the Data
metallicityGrid = COMPAS.metallicityGrid

#Create instance Class
MSSFR = ClassMSSFR.MSSFR(metallicityGrid=metallicityGrid)

# -

#Set the MSSFR model
MSSFR.SFRprescription = 'Neijssel et al. (2019)'
MSSFR.Zprescription = 'logNormal'
MSSFR.logNormalPrescription ='Neijssel Phenomenological'

# # 2 - Define the redshifts
#
# The entire calculation depends on defining a redshift
# at which the DCOs merge. Then using the delaytimes and astropy
# you need to recover when the systems were born
#
# First define our cosmology

# +
# see top notebook

# +
mergerRedshift = 0.2

#Define an age when the first stars formed based on redshift
firstSFR = cosmology.age(10).value
#single value in units of Gyr
ageUniverseAtMergerGyr = cosmology.age(mergerRedshift)
#Recover the array delaytimes in Units of Gyr
delayTimeGyr = np.divide(COMPAS.delayTimes, 1000)
#Calculate the age at which they were born
ageBirth   = ageUniverseAtMergerGyr.value - delayTimeGyr
#If they are born before first SFR mask them
maskUnreal = ageBirth<firstSFR
#I set those to minus 1 to label them
#This way we can use this as a mask everywhere instead of having
#to slice the data in different ways. Useful for larger calculations
ageBirth[maskUnreal] = -1

# -

# Note that the above might further reduce the number of DCOs in the data, despite we use the flag merger in a Hubble time. This because that flag assumes redshift is zero. When we change our reference frame to higher redshifts, then the universe is younger and therefore fewer systems will be able to merge in  time
#
# We set the unphysical systems to -1. This because later when we loop over redshifts the number of possible systems can vary. However I want to fill in the rates in a predefined 2D array of fixed shape (nr systems, nr of redshift bins). Hence I assume the largest array (size of nr of DCOs) and set the rate to zero in case.  Note that the MSSFR script also depends on the mask of systems being unreal of -1. (see returnMSSFR())
#
# To recover the redshift it is a bit tricky. Astropy can quicly
# calculate the age from redshift, but the inverse is more difficult.
# In the code we use look up the nearest value in a dense precalculated table.
# Here we use a older method to calculate (for credits see source code classCosmicintegrator) which is considerably slower. 

redshiftsBirth = np.zeros(len(ageBirth))
for nr, age in enumerate(ageBirth):
    if age != -1:
        redshiftsBirth[nr] = newton(lambda x: cosmology.age(x).value-age, 0)
        
    else:
        redshiftsBirth[nr] = -1

print("nr of DCOs %s, nr DCOs merging %s"\
     %(len(COMPAS.delayTimes), np.sum(ageBirth!=-1)))


# # Calculate the rate of systems per metallicity
#
# The code is structured to do the calculation per subpopulation of DCOs of a single metallicity. Note that if the system was not physically possible (age == -1) then the rate is set to zero.

#create an array for rate per system merging at redshift z
ratePerSystem = np.zeros(len(COMPAS.delayTimes))

for nrZ, Z in enumerate(metallicityGrid):
    maskZ    = COMPAS.metallicitySystems == Z
    #give MSSFR per system which has metallicity Z [Msun/dGpc3/dyr]
    mssfr    = MSSFR.returnMSSFR(metallicity=Z,\
                                 agesBirth  =ageBirth[maskZ],
                                 redshiftBirth=redshiftsBirth[maskZ])
    #Calculate rate using amount of Msun evolved [dN/dGpc3/dyr]
    RatesZ   = np.divide(mssfr, COMPAS.totalMassEvolvedPerZ[nrZ])
    #Fill the rates in the defined array according to mask
    ratePerSystem[maskZ] = RatesZ

print(metallicityGrid)

print(ratePerSystem)

print(np.sum(ratePerSystem))

# +
# Using the rates in a histogram

chirpMasses = COMPAS.mChirp

binsM = np.linspace(0,30,100)
dM    = np.diff(binsM)
center= (binsM[1:]+binsM[:-1])/2.

#remember that the reight is essentially a weight per system

y , _     = np.histogram(chirpMasses, bins=binsM, \
                     weights=ratePerSystem)

dydMchirp = np.divide(y, dM)

fig, axes = plt.subplots(1,1, figsize=(9,8))
axes.plot(center, dydMchirp)
axes.set_xlabel('chirp mass [Msun]', fontsize=20)
axes.set_ylabel('rate [yr-1 Gpc-3]', fontsize=20)
axes.set_title('merger rate density at z=%s'\
               %(mergerRedshift), fontsize=20)
plt.tight_layout()
plt.show()

# -

# Now here we only have the chirp mass distribution at a single redshift
# per unit volume. The next step is to do it over a range of redshifts
# and get the absolute rates.


