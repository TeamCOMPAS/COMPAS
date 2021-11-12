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
# Most of the post-processing material here uses python code.
# Here we show two important basics of python code in the context
# of investigating your COMPAS simulation.


# ## Material
#
# ### [1. Slicing the data](#1.-Slicing-the-data)
# How to select specific systems and their parameters using seeds.
#
# ### [2. Visualizing the data](#2.-Visualizing-the-data)
# The basics of binning and visualising your data.




# +
#python libraries
import numpy as np               # for handling arrays
import h5py as h5                # for reading the COMPAS data
import time                      # for finding computation time
import matplotlib.pyplot as plt  #for plotting

# %matplotlib inline

# Choose an output hdf5 file to work with
compasRootDir = os.environ['COMPAS_ROOT_DIR'] 
pathToData = compasRootDir + 'postProcessing/Tutorial/COMPAS_Output/COMPAS_Output.h5'

# -

# # 1. Slicing the data
#
# One of the most important numbers in the COMPAS output is the system seed. The seed represents the unique identifier to a specific system in a simulation. It is also used as the seed value in random number generation, which is useful when trying to reproduce an identical system. Therefore the properties of a single system can be recovered by looking at its seed across different output categories. 
#
# Here we introduce the basics of manipulating the data using the seeds. We provide an example on how we get the initial parameters of systems that ended up forming double compact objects.
#
# Naively, we might try to use For Loops with Conditions to extract systems of interest to a list. However, this can potentially be computationally expensive.
#
# Here we present a method to more efficiently 'slice' the data using numpy and boolean masks. These are slightly more involved but are computationally quick and use intuitive logic.
#
# --
#
# If you do not already have a COMPAS_Output.h5 ready, see Section 1 [Working With HDF5](./WorkingWithHDF5.py) on how to create your own output file, or download some data from our [Zenodo database](https://zenodo.org/communities/compas/?page=1&size=20).

# *Note:* These cells may take a long time if you test them on large datasets.

Data  = h5.File(pathToData)
print(list(Data.keys()))


# The print statement shows the different categories that are combined in your h5file.
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

# # 2. Visualizing the data
#
# Although math is the fundamental basis of physics and astrophysics, we cannot always easily convert numbers and equations into a coherent picture. Plotting is therefore a vital tool in bridging the gap between raw data and a deeper scientific understanding. 
#
# *Disclaimer:*
#
# There are many ways to make the same plot in matplotlib and there are many ways to bin your data. Often, there is no "best" way to display data in a plot, and the message conveyed can be heavily dependent on the context of the data as well as asthetic plotting decisions.
#
# For example, in histograms, as we discuss below, the relatively subjective choice of bin size can significantly affect the interpretation of the results. It is important to be aware of when and how we make these choices and to try to reduce any unintended bias.
#
# ---
#
# ** Example: inspect the component masses of Double Compact Objects**
#
# In the example below, we use the following conventions:
#
# 1 - We deliberately choose to use the matplotlib.pyplot.subplots routine even when creating a single figure (as opposed to using pyplot.plot). This is because many online forums (e.g Stackoverflow) use this syntax. Furthermore, this means you do not have to learn two different types of syntax when creating either a single or multiple panel figure.
#
# 2 - We choose to do the binning within the numpy/array environment instead of with inbuilt functions such as plt.hist / axes.hist. The reason is that you have more control over what you do, such as custom normalization (using rates, weights, pdf, etc.). It also forces you to have a deeper understanding of what you are calculating, and allows you to check intermediate steps with print statements.  Once you know how to bin your data this way you can also easily expand these routines for more complicated plots (2D binning).

# # Path to be set by user
#

pathToData = '/home/cneijssel/Desktop/Test/COMPAS_output.h5'

# # Get some data to plot

# +
Data  = h5.File(pathToData)

print(list(Data.keys()))
# DCOs = double compact objects


DCOs = Data['DoubleCompactObjects']

M1   = DCOs['Mass_1'][()]
M2   = DCOs['Mass_2'][()]
Mtot = np.add(M1, M2)
# -

Data.close()

# # Histogram

# +
# You can use numpy to create an array with specific min, max and interval values
minMtot = 0
maxMtot = max(Mtot)
nBins   = 50

# Number of bin edges is one more than number of bins
binEdges = np.linspace(minMtot, maxMtot, nBins+1)

# What is the value at the center of the bin?
# add each edge of the side of the bin and divide by 2
xvaluesHist  = (binEdges[:-1] + binEdges[1:])/2.

# What is the width of each bin? (an array in general, if the spacing is non-uniform)
binWidths = np.diff(binEdges)


### Set yvalues to the height of the bins

# Create an array of y-values for each x-value
yvalues = np.zeros(len(xvaluesHist))

# Iterate over the bins to calcuate the number of data points per bin
for iBin in range(nBins):
    mask = (Mtot >= binEdges[iBin]) & (Mtot < binEdges[iBin+1])
    yvalues[iBin] = np.sum(mask)

# You can of course apply any mask you like to get the desired histogram    

## Generally, you can calculate the rate per unit x (dy/dx) using
dYdXHist = np.divide(yvalues, binWidths)

# To convert your distribution to a PDF, normalize in y-values:
PDF = np.divide(yvalues, np.sum(yvalues))

# You can then multiply by, e.g, rates/weights to scale the distribution
# -

# # CDF
#
# Sometimes we want to know what fraction of the data lies below a given value. To find this, we calculate a Cumulative Distribution Function, or CDF.

# +
# Question: How many points have a value less than X? 

# Sort the values of interest
MtotSorted = np.sort(Mtot)   

# These values are your xvalues 
xvaluesCDF = MtotSorted

# The CDF is a non-strictly increasing function from 0 to 1 across the range of x values.
# It should increment by 1/len(xvaluesCDF) at each x in the array, and remain constant otherwise.

# Numpy provides several functions that make this very straightforward
nDataPoints = len(xvaluesCDF)
yvalues = np.cumsum(np.ones(nDataPoints))
CDF = yvalues / nDataPoints
# -

# # A two panel plot

# +
# For two panels side by side:
fig, axes = plt.subplots(1,2, figsize=(18,8))

# axes is an array relating to each panel
# panel1 = axes[0]
# panel2 = axes[1]

largefontsize = 18
smallfontsize = 13




### In the left panel, we want to plot the histogram and CDF overlayed
### with the same x-axis, but different y-axes

# Plot the Histogram first
histAxes = axes[0]
histAxes.plot(xvaluesHist, dYdXHist)

histAxes.set_xlabel(r'Mtot [$M_\odot$]', fontsize=smallfontsize)
histAxes.set_ylabel(r'dN/dMtot [$M_\odot^{-1}$]', fontsize=smallfontsize)

# Overlay the CDF with the same x-axis but different y-axis
cdfAxes =  axes[0].twinx()
cdfAxes.plot(xvaluesCDF, CDF, c='r')

# Dont have to do xlabel since they are the same
cdfAxes.set_ylabel('CDF', fontsize=smallfontsize, labelpad=-40)
cdfAxes.tick_params(axis='y', direction='in', pad=-20) # Adjust the CDF axis for clarity in the plot

axes[0].set_title('Total Mass Histogram and CDF', fontsize=largefontsize)




### In the right panel, we want to display a scatterplot of M1 & M2 

axes[1].scatter(M1, M2)
axes[1].set_xlabel(r'M1 [$M_\odot$]', fontsize=smallfontsize)
axes[1].set_ylabel(r'M2 [$M_\odot$]', fontsize=smallfontsize)

axes[1].set_title('Component Masses', fontsize=largefontsize)




### Clean up and display the plot

# You can force plt to pad enough between plots
# such that the labels fit
plt.tight_layout()

# If you want to save the figure, use:
#plt.savefig(pathToSave)

# To produce the plot, always remember to:
plt.show()
# -


