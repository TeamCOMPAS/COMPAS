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
#
# The most important number in the COMPAS data is the seed. The seed represents the unique identifier to a specific system in a simulation. Therefore the properties of a single system can be recovered by looking at seeds in different types of files. 
#
# Here we introduce the basics of manipulating the data using the seeds. We provide an example on how we get the initial parameters of systems that ended up forming double compact objects.
#
# Naively, we might try to use For Loops with Conditions to extract systems of interest to a list. However, this can potentially be computationally expensive.
#
# Here we present a method to more efficiently 'slice' the data using boolean masks. These are slightly more involved but are computationally quick and use intuitive logic.
#
# If you do not already have a COMPAS_Output.h5 ready, see Section 1 on how to create the h5 file using the csv data of your simulation, or download some data from [compas.science.](https://compas.science/)

# ** *Note:* These cells may take a long time if you test them on large datasets.**

# # Path to be set by user
#

pathToData = '/home/cneijssel/Desktop/Test/COMPAS_output.h5'

# # Imports

#python libraries
import numpy as np               # for handling arrays
import h5py as h5                # for reading the COMPAS data
import time                      # for finding computation time

Data  = h5.File(pathToData)
print(list(Data.keys()))


# The print statement shows the different types of files that are combined in your h5file.
#
# The system seed links, e.g, information about the Supernovae to information about the initial SystemParameters.

# # Question: What were the initial total masses of the double compact objects?

def calculateTotalMassesNaive(pathData=None):
    Data  = h5.File(pathToData)
    
    totalMasses = []
    
    #for syntax see section 1 
    seedsDCOs     = Data['DoubleCompactObjects']['SEED'][()]
    
    #get info from ZAMS
    seedsSystems  = Data['SystemParameters']['SEED'][()]
    M1ZAMSs       = Data['SystemParameters']['Mass@ZAMS_1'][()]
    M2ZAMSs       = Data['SystemParameters']['Mass@ZAMS_2'][()]

    for seedDCO in seedsDCOs:
        for nrseed in range(len(seedsSystems)):
            seedSystem = seedsSystems[nrseed]
            if seedSystem == seedDCO:
                M1 = M1ZAMSs[nrseed]
                M2 = M2ZAMSs[nrseed]
                Mtot = M1 + M2
                totalMasses.append(Mtot)

    Data.close()
    return totalMasses


# +
# calculate function run time
start   = time.time()
MtotOld = calculateTotalMassesNaive(pathData=pathToData)
end     = time.time()
timeDiffNaive = end-start

print('%s seconds, using for loops.' %(timeDiffNaive)) 
# -

# # Optimizing the above loop
#
# ## 0 - Use built-in numpy routines

# Numpy arrays can make use of a powerful library of optimization tools which allow the user to bypass computationally heavy for loops. 
#
# For example, we can speed up the calculation of the element-wise sum of two arrays with:

# +
M1ZAMS  = Data['SystemParameters']['Mass@ZAMS_1'][()]
M2ZAMS  = Data['SystemParameters']['Mass@ZAMS_2'][()]
    
mTotalAllSystems  = np.add(M1ZAMS, M2ZAMS)
# -

# ## 1 - Use boolean masks in a single file

# There is a useful trick for when you want only those elements which satisfy a specific condition. 
#
# Where previously we put the condition in an if statement nested within a for loop, now we will use an array of booleans to mask out the undesired elements. 
#
# The boolean array will have the same length as the input array, with 

# +
# Create a boolean array from the total mass array which is True
# if the total mass of the corrresponding system is less than 40. 

maskMtot = (mTotalAllSystems <= 40)
# -

# **Crucially, you can apply this mask to all other columns in the same file because, by construction, they all have the same length.**

# seeds of systems with total mass below 40
seeds  = Data['SystemParameters']['SEED'][()]
seedsMtotBelow40 = seeds[maskMtot]

# Note that this works because the order of the two columns (seeds and total masses) are the same. 
#
# For example, the total mass of the third system entry corresponds to the seed at the third system entry.

# ## 2 - Use seeds as masks between files

# ### Example 1
#
# Before we continue it is useful to understand how the COMPAS-popsynth printing works.
#
# Each simulated system will be initialized only once and so will have only one line in the SystemParameters file. However, lines in CommonEnvelopes are created whenever a system goes through CE, which might happen multiple times for a single system, or potentially not at all. Similarly, in the Supernovae file, you will find at most two lines per system, but possibly none. DoubleCompactObject lines are printed only when both remnants are either Neutron Stars or Black Holes (but also includes disrupted systems), which happens at most once per system. 
#
# For this reason, it is in general not the case that the system on line $n$ of one file corresponds will match the system on line $n$ of another file.
#
# In order to match systems across files, we need to extract the seeds of desired systems from one file, and apply them as a mask in the other file. 

# +
# example mock data from two files
SystemSeeds = np.array([1,  2,  3,  4 ])
SystemMass1 = np.array([1, 20,  5, 45 ])
DCOSeeds    = np.array([    2,      4 ])

# Calculate mask for which elements of SystemSeeds are found in DCOSeeds - see numpy.in1d documentation for details
mask = np.in1d(SystemSeeds, DCOSeeds)

print(mask)
print(SystemSeeds[mask])
print(SystemMass1[mask])


# -

# # Optimized loop

def calculateTotalMassesOptimized(pathData=None):
    Data  = h5.File(pathToData)
    
    totalMasses = []
    
    #for syntax see section 1 with basic syntax
    seedsDCOs     = Data['DoubleCompactObjects']['SEED'][()]
    #get info from ZAMS
    seedsSystems  = Data['SystemParameters']['SEED'][()]
    M1ZAMSs       = Data['SystemParameters']['Mass@ZAMS_1'][()]
    M2ZAMSs       = Data['SystemParameters']['Mass@ZAMS_2'][()]
    
    MZAMStotal    = np.add(M1ZAMS, M2ZAMS)
    
    maskSeedsBecameDCO  = np.in1d(seedsSystems, seedsDCOs)
    totalMassZAMSDCO    = MZAMStotal[maskSeedsBecameDCO]
    
    Data.close()
    return totalMassZAMSDCO


# +
# calculate function run time
start   = time.time()
MtotNew = calculateTotalMassesNaive(pathData=pathToData)
end     = time.time()
timeDiffOptimized = end-start

# calculate number of Double Compact Objects
nrDCOs = len(Data['DoubleCompactObjects']['SEED'][()])

print('Compare')
print('%s seconds, using Optimizations.' %(timeDiffOptimized)) 
print('%s seconds, using For Loops.'     %(timeDiffNaive)) 
print('Using %s DCO systems'             %(nrDCOs))
# -

# *Note:* The time difference will depend on the number of systems under investigation, as well as the number of bypassed For Loops.

# test that the two arrays are in fact identical
print(np.array_equal(MtotOld, MtotNew))


# Note that the above loop can easily be expanded with more conditions.
#
# If you do not want all the DCO initial total masses but only of the double neutron stars, then you just need to apply another mask to the seedsDCOs.

def calculateTotalMassesDNS(pathToData=None):
    Data  = h5.File(pathToData)
    
    totalMasses = []
    
    #for syntax see section 1 with basic syntax
    seedsDCOs     = Data['DoubleCompactObjects']['SEED'][()]
    type1         = Data['DoubleCompactObjects']['Stellar_Type_1'][()]
    type2         = Data['DoubleCompactObjects']['Stellar_Type_2'][()]
    maskDNS       = (type1 == 13) & (type2 == 13)
    seedsDNS      = seedsDCOs[maskDNS]
    
    #get info from ZAMS
    seedsSystems  = Data['SystemParameters']['SEED'][()]
    M1ZAMSs       = Data['SystemParameters']['Mass@ZAMS_1'][()]
    M2ZAMSs       = Data['SystemParameters']['Mass@ZAMS_2'][()]
    
    MZAMStotal    = np.add(M1ZAMS, M2ZAMS)
    
    
    maskSeedsBecameDNS  = np.in1d(seedsSystems, seedsDNS)
    totalMassZAMSDNS    = MZAMStotal[maskSeedsBecameDNS]
    
    Data.close()
    return totalMassZAMSDNS


# +
# calculate function run time
start   = time.time()
MtotDNS = calculateTotalMassesDNS(pathToData=pathToData)
end     = time.time()
timeDiffDNS = end-start

# calculate number of DNS systems
nrDNSs = len(MtotDNS)
    
print('%s seconds for all %s DNS systems.' %(timeDiffDNS, nrDNSs)) 
# -

# ### Example 2
#
# The previous example uses the fact that both SystemParameters and DoubleCompactObjects only print at most one line per system. However, as mentioned above, events such as supernovae or common envelopes might happen multiple times to a given system, and as a result there would be multiple occurences of a given seed in the relevant file. 
#
# To account for this, we will need to modify the previous method. Consider again the 4 seeds of the previous example. Both 2 and 4 formed a DCO and hence both stars in these binaries went SN. Seeds 1 and 3 are low mass stars hence they did not go SN. (Note that we do not specify the companion masses for any of these systems, but for simplicity we assume that the companions to 1 and 3 are also sufficiently low mass to not produce a supernova). The SN file prints one line per SN and therefore seeds 2 and 4 appear twice each.
#
# Imagine you want the primary masses of systems that experienced at any point a core collapse supernova (CCSN). We'll reuse our mock data, with additional information about the types of SN which occured in each star. Here, PPISN refers to Pulsational Pair Instability Supernovae.

# +
# example mock data from above
SystemSeeds = np.array([1,  2,  3,  4 ])
SystemMass1 = np.array([1, 20,  5, 45 ])
DCOSeeds    = np.array([    2,      4 ])

SNSeeds     = np.array([     2,      2,      4,       4 ])  
SNTypes     = np.array(['CCSN', 'CCSN', 'CCSN', 'PPISN' ])

# get seeds which had a CCSN
maskCCSN  = SNTypes == 'CCSN'
seedsCCSN = SNSeeds[maskCCSN]
print('CCSN seeds =%s' %(seedsCCSN))

#compare which element of 1-d array are in other
#this because in 

seedsCCSN = np.unique(seedsCCSN)
# in this particular case, it is not necessary to reduce seedsCCSN to it's unique entries.
# the numpy.in1d function will work with duplicate seeds, but we include it explicitly here
# as other more complicated scenarios might rely on unique sets of seeds

mask = np.in1d(SystemSeeds, seedsCCSN)
print(SystemMass1[mask])
# -

# Always remember to close your data file
Data.close()
