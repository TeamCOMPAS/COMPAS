# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.5
#   kernelspec:
#     display_name: base
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Basic COMPAS data analysis
#
# Most of the post-processing material presented here is written in python, and makes extensive use of the numpy package for rapid computation on large data arrays.
# Here we show two important basics of python/numpy in the context of investigating your COMPAS simulation.


# %% [markdown]
# ## Material
#
# ### [1) Inspecting the data ](#1.-Inspecting-the-data)
# Look at the data to see which parameters are available and check that it matches expectations.
#
# ### [2) Slicing the data ](#2.-Slicing-the-data)
# Select specific systems and their parameters using seeds.
#
# ### [3) Visualizing the data ](#3.-Visualizing-the-data)
# Binning and visualising your data.

# %%


# %% [markdown]
# ### For the following sections, you will need to have the following packages installed.
# ### `numpy, h5py, time, matplotlib`

# %%
#python libraries
import os, sys
import numpy as np               # for handling arrays
import h5py as h5                # for reading the COMPAS data
import time                      # for finding computation time
import matplotlib.pyplot as plt  #for plotting

# Import COMPAS specific scripts
compas_root_dir = os.environ['COMPAS_ROOT_DIR'] 
sys.path.append(compas_root_dir + 'compas_python_utils')
from debugging_utils import print_compas_details_dataframe, get_event_history, get_event_strings

# Choose an output hdf5 file to work with
path_to_data = 'COMPAS_Tutorial_Output.h5'

# This is known as an ipython magic command, and allows plots to be produced within the notebook
# %matplotlib inline
# %%


# %% [markdown]
# ## 1. Inspecting the data
#
# Often the first thing you want to do with new data is simply to look at it! Getting familiar with the data, including available parameters, size of the data file, etc. will help to inform how best to proceed with the analysis. We provide several useful functions for inspecting the data, `print_compas_details_dataframe`, `get_event_history`, and `get_event_strings`.
#
# --
#
# If you want to create an alternative COMPAS_Output.h5, see Section 1 [Working With HDF5](./WorkingWithHDF5.ipynb), or download some data from our [Zenodo database](https://zenodo.org/communities/compas/?page=1&size=20).

# %% [markdown]
# *Note:* These cells may take a long time if you test them on large datasets.

# %%
data  = h5.File(path_to_data)
print(list(data.keys()))

# %% [markdown]
# The output above represents the event categories available from the particular run. If you used the output produced in the previous tutorial, you should see `['BSE_Common_Envelopes', 'BSE_Double_Compact_Objects', 'BSE_RLOF', 'BSE_Supernovae', 'BSE_System_Parameters', 'Run_Details']`. Note that for smaller runs which do not produce any of a particular type of output, the output category will not be created. 
#
# Brief description of the categories:
# - 'BSE_System_Parameters': Initial state of the binary
# - 'BSE_RLOF': Any mass transfer events that occured within the binary
# - 'BSE_Common_Envelopes': If any of the mass transfer events were unstable, details will be included here.
# - 'BSE_Supernovae': Parameters and outcome of any supernovae that occured in the binary
# - 'BSE_Double_Compact_Objects': Includes key information of all binaries which end their lives as an intact pair of compact obects (either neutron stars or black holes)
# - 'Run_Details': Information on the input settings supplied to the Compas run
#
# To extract the data from these categories, we use the following syntax

# %%
SPs = data['BSE_System_Parameters']
MTs = data['BSE_RLOF']
CEs = data['BSE_Common_Envelopes']
SNe = data['BSE_Supernovae']
DCs = data['BSE_Double_Compact_Objects']

# %% [markdown]
# Each of these is a dictionary mapping parameter names (keys) to an array of values

# %%
print(SPs.keys())

# %% [markdown]
# One of the most important parameters in the COMPAS output is the system seed. The seed represents the unique identifier to a specific system in a simulation. It is also used as the seed value in random number generation, which is useful when trying to reproduce a given system identically. 
#
# If we want to view the random seeds in the system parameters file, we run

# %%
seeds_SP = SPs['SEED'][()]
print(seeds_SP)

# %% [markdown]
# ### print_compas_details_dataframe
#
# This is useful for extracting the arrays of single parameters, but for a more convenient view of the whole system parameters file, we can use the `print_compas_details_dataframe` function.

# %%
print_compas_details_dataframe(SPs) # Note - the output of this is a pandas dataframe

# %% [markdown]
# `print_compas_details_dataframe` optionally also takes seeds as arguments, to focus on specific systems.

# %%
seeds_MT = MTs['SEED'][()]
first_three_unique_seeds = np.unique(seeds_MT)[0:3]
print("Look at these seeds: ", first_three_unique_seeds)
print_compas_details_dataframe(MTs, first_three_unique_seeds)

# %% [markdown]
# ### get_event_history
#
# Often, it is useful to quickly retrieve an overview of the event history of a binary, including all mass transfer, common envelope, and supernova events. For this, we use `get_event_history`

# %%
seeds, events = get_event_history(data)

for ii, seed in enumerate(seeds):
    print(seed, events[ii])

# %% [markdown]
# `get_event_history` takes the h5file as input, and returns an array of the seeds processed as well as the major events for that seed. The format for events depends on the event type. Currently, we only include supernova and mass transfer events, but mass transfer events include a flag for whether the system underwent CEE.
#
# - For MT events, it is ('RL|CE|MG', time, stellar_type_primary, stellar_type_secondary, is_rlof1, is_rlof2, is_cee, is_merger)
# - For SN events, it is ('SN', time, stellar_type_progenitor, stellar_type_remnant, which_star_is_progenitor, is_binary_unbound)
#
#
# There is also an optional argument `exclude_null` which defaults to False. If True, it will skip systems which undergo no events of interest (which may speed up large runs).

# %% [markdown]
# A useful function that builds off of `get_event_history` is `get_event_strings`, which collects the event information into a succint string, which may be easier to read (once you get used to the syntax).
#
# The syntax for the event strings takes the following convention:
#     
# - For MT events:
#     - P>S, P<S, P=S, or P&S
#     - where P is primary type, S is secondary type, and `>`, `<` is RLOF (1->2 or 1<-2), `=` is a (successful) CEE, and `&` is a merger.
#
# - For SN events:
#     - P\*SR for star1 the SN progenitor, or 
#     - R\*SP for star2 the SN progenitor,
#     - where P is progenitor type, R is remnant type, 
#       S is state (`i` for intact, `u` for unbound)
#     - Note that the companion type is not reported in this truncated syntax.
#
# Event strings for the same seed are ordered chronologically and separated by the undesrcore character `_`
#

# %%
event_strings = get_event_strings(all_events=events)
for ii in range(5):
    print(seeds[ii], event_strings[ii])


# %%

# %% [markdown]
# ## 2. Slicing the data
#
# Since the random seed is unique and constant for a given binary, the properties and events of the binary system can be recovered by looking at its seed across different output categories. 
#
# Here we introduce the basics of manipulating the data using the seeds. We provide an example on how we get the initial parameters of systems that ended up forming double compact objects.
#
# Naively, we might try to use For Loops with Conditions to extract systems of interest to a list. However, this can potentially be computationally expensive.
#
# Here we present a method to more efficiently 'slice' the data using numpy and boolean masks. These are slightly more involved but are computationally quick and use intuitive logic.

# %% [markdown]
# ### Question: What were the initial total masses of the double compact objects?

# %%
def calculate_total_masses_naive(path_data=None):
    data  = h5.File(path_to_data)
    
    total_masses = []
    
    # Retrive the categories
    SPs = data['BSE_System_Parameters']
    DCs = data['BSE_Double_Compact_Objects']
    
    # For syntax see section 1 
    
    # Extract parameters of interest
    seeds_DC       = DCs['SEED'][()]
    seeds_SP       = SPs['SEED'][()]
    m1_zams        = SPs['Mass@ZAMS(1)'][()]
    m2_zams        = SPs['Mass@ZAMS(2)'][()]

    for dc_seed in seeds_DC:
        for seed_index in range(len(seeds_SP)):
            sp_seed = seeds_SP[seed_index]
            if sp_seed == dc_seed:
                m1 = m1_zams[seed_index]
                m2 = m2_zams[seed_index]
                m_tot = m1 + m2
                total_masses.append(m_tot)

    data.close()
    return total_masses


# %%
# calculate function run time
start   = time.time()
m_tot_old = calculate_total_masses_naive(path_data=path_to_data)
end     = time.time()
time_diff_naive = end-start

print('%s seconds, using for loops.' %(time_diff_naive)) 

# %% [markdown]
# ### I) Optimizing the above loop
#
# #### a - Use built-in numpy routines

# %% [markdown]
# Numpy arrays can make use of a powerful library of optimization tools which allow the user to bypass computationally heavy for-loops. 
#
# For example, we can speed up the calculation of the element-wise sum of two arrays with:

# %%
SPs = data['BSE_System_Parameters']

m1_zams  = SPs['Mass@ZAMS(1)'][()]
m2_zams  = SPs['Mass@ZAMS(2)'][()]
    
m_total_all_systems  = np.add(m1_zams, m2_zams)

# %% [markdown]
# #### b - Use boolean masks in a single file

# %% [markdown]
# Where previously we put the condition in an if statement nested within a for loop, now we again make use of boolean masks to filter out the undesired elements. 
#
# The boolean array must have the same length as the input array.

# %%
# Create a boolean array from the total mass array which is True
# if the total mass of the corrresponding system is less than 40. 

mask_m_tot_less_than_40 = (m_total_all_systems <= 40)

# %% [markdown]
# **Crucially, you can apply this mask to all other columns in the same file because, by construction, they all have the same length.**

# %%
# seeds of systems with total mass below 40
seeds_m_tot_below_40 = seeds_SP[mask_m_tot_less_than_40]

# %% [markdown]
# Note that this works because the order of the two columns (seeds and total masses) are the same. 
#
# For example, the total mass of the system at index 2 corresponds to the seed at index 2.

# %% [markdown]
# ### II) Use seeds as masks between files

# %% [markdown]
# #### Example 1
#
# Before we continue it is useful to understand how the COMPAS-popsynth printing works.
#
# Each simulated system will be initialized only once and so will have only one line in the `BSE_System_Parameters` file. However, lines in `BSE_RLOF` are created whenever a system goes through a mass transfer event, which might happen multiple times for a single system, or potentially not at all. Similarly, in the `BSE_Supernovae` file, you will find at most two lines per system, but possibly none. `BSE_Double_Compact_Objects` lines are printed only when the final system is intact and composed of either Neutron Stars or Black Holes, which is a rare event that happens at most once per system. 
#
# For this reason, it is generally not the case that the system on line $n$ of one file corresponds to the system on line $n$ of another file.
#
# In order to match systems across files, we need to extract the seeds of desired systems from one file, and apply them as a mask in the other file. 

# %%
# Example: calculate the primary ZAMS mass of systems which become DCOs (Double Compact Objects)
seeds_SP = SPs['SEED'][()]
seeds_DC = DCs['SEED'][()]
m1_zams  = SPs['Mass@ZAMS(1)'][()]

# Calculate mask for which elements of seeds_SP are found in seeds_DC
# - see numpy.isin documentation for details
mask = np.isin(seeds_SP, seeds_DC)

print("The occurence rate of DCOs is {}/{}".format(sum(mask), len(mask)))

# %%
seeds_DC = DCs['SEED'][()]
print_compas_details_dataframe(DCs, seeds_DC[:3])


# %% [markdown]
# #### Optimized loop
#

# %%
def calculate_total_masses_optimized(path_data=None):
    data  = h5.File(path_to_data)
    
    total_masses = []
        
    # Retrive the categories
    SPs = data['BSE_System_Parameters']
    DCs = data['BSE_Double_Compact_Objects']
    
    # For syntax see section 1 
    
    # Extract parameters of interest
    seeds_DC       = DCs['SEED'][()]
    seeds_SP       = SPs['SEED'][()]
    m1_zams        = SPs['Mass@ZAMS(1)'][()]
    m2_zams        = SPs['Mass@ZAMS(2)'][()]
    
    m_zams_tot             = np.add(m1_zams, m2_zams)
    mask_seeds_became_DCO  = np.isin(seeds_SP, seeds_DC)
    m_zams_tot_of_DCOs     = m_zams_tot[mask_seeds_became_DCO]
    
    data.close()
    return m_zams_tot_of_DCOs


# %%
# calculate function run time
start   = time.time()
m_tot_new = calculate_total_masses_optimized(path_data=path_to_data)
end     = time.time()
time_diff_optimized = end-start

# calculate number of Double Compact Objects
n_DCos = len(seeds_DC)

print('Compare')
print('%s seconds, using For Loops.'     %(time_diff_naive)) 
print('%s seconds, using Optimizations.' %(time_diff_optimized)) 
print('Using %s DCO systems'             %(n_DCos))

# %% [markdown]
# *Note:* The time difference will depend heavily on the number of systems under investigation, as well as the number of bypassed For Loops. If you used the path to the pre-generated tutorial data set, you should see very little improvement. 

# %%
# Test that the two arrays are in fact identical
print(np.array_equal(m_tot_old, m_tot_new))


# %% [markdown]
# *Note:* the above loop can easily be expanded with more conditions.
#
# If you do not want all the DCO initial total masses but only of the binary black holes, then you just need to apply another mask to the seeds_DC.

# %%
def calculate_total_masses_bbh(path_to_data=None):
    data  = h5.File(path_to_data)
    
    total_masses = []
    
    SPs = data['BSE_System_Parameters']
    DCs = data['BSE_Double_Compact_Objects']

    seeds_DC = DCs['SEED'][()]
    stype1  = DCs['Stellar_Type(1)'][()]
    stype2  = DCs['Stellar_Type(2)'][()]

    dc_mask_bbh     = (stype1 == 14) & (stype2 == 14)
    seeds_bbh      = seeds_DC[dc_mask_bbh]
    
    # Get info from ZAMS
    seeds_SP  = SPs['SEED'][()]
    m1_zams   = SPs['Mass@ZAMS(1)'][()]
    m2_zams   = SPs['Mass@ZAMS(2)'][()]
    
    m_zams_tot = np.add(m1_zams, m2_zams)    
    
    sp_mask_bbh    = np.isin(seeds_SP, seeds_bbh)
    m_zams_tot_bbh = m_zams_tot[sp_mask_bbh]
    
    data.close()
    return m_zams_tot_bbh


# %%
# calculate function run time
start   = time.time()
m_tot_bbh = calculate_total_masses_bbh(path_to_data=path_to_data)
end     = time.time()
time_diff_bbh = end-start

# calculate number of BBH systems
n_bbh = len(m_tot_bbh)
    
print('%s seconds for all %s BBH systems.' %(time_diff_bbh, n_bbh)) 

# %% [markdown]
# Note that the `print_compas_details_dataframe` function can also optionally take a mask as argument. The mask array must have the same length as the data arrays for the given category.

# %%
mask_merges_hubble_time = DCs['Merges_Hubble_Time'][()] == 1 # Booleans are stored as 0 or 1, so need the final == 1 to get a boolean array
print_compas_details_dataframe(DCs, mask=mask_merges_hubble_time)


# %% [markdown]
# ### Example 2
#
# The previous example uses the fact that both `BSE_System_Parameters` and `BSE_Double_Compact_Objects` only print at most one line per system. However, as mentioned above, events such as supernovae or mass transfer might happen multiple times to a given system, and as a result there would be multiple occurences of a given seed in the relevant file. 

# %%
# Example: Want to investigate CEE events for a given system. 
# 
# To illustrate the point, we find the seed of the system with the most mass transfer events, one of which is a CEE.

seeds_CE = CEs['SEED'][()]
uniq_seeds, uniq_counts = np.unique(seeds_MT[np.isin(seeds_MT, seeds_CE)], return_counts=True)
best_seed = uniq_seeds[np.argmax(uniq_counts)]  
print_compas_details_dataframe(MTs, best_seed)

# %% [markdown]
# There is too much information here to be useful in this case, so we additionally apply a mask to filter out events we're not interested in

# %%
mask_cee = MTs['CEE>MT'][()] == 1 # true for mass transfer events which undergo Common Envelope Evolution
print_compas_details_dataframe(MTs, best_seed, mask=mask_cee)

# %% [markdown]
#
# ### Example 3
#
# Combining masks on the seeds and other data can provide a lot of flexibility to help explore your science case. Imagine you want the primary masses of systems that experienced two core collapse supernovae (CCSNe) and resulted in a double compact object that will merge in a Hubble time. We'll reuse our mock data, with additional information about the types of SN which occured in each star. 

# %%
# Example: get the primary ZAMS masses of systems which experience 2 CCSNe before becoming a DCO

seeds_SP = SPs['SEED'][()]
seeds_SN = SNe['SEED'][()]
seeds_DC = DCs['SEED'][()]

sn_type  = SNe['SN_Type(SN)'][()]
m1_zams  = SPs['Mass@ZAMS(1)'][()]
merges_hubble_time = DCs['Merges_Hubble_Time'][()] == 1

# Note: the SN_Type(SN) parameter maps integers to SN types - see documentation for details
sn_type_dict = {
    1: 'CCSN',
    2: 'ECSN',
    4: 'PISN',
    8: 'PPISN', 
    16: 'USSN',
} # this dictionary is illustrative, but not explicitly used here


# Determine which seeds experienced 2 CCSNe
mask_ccsn = sn_type == 1
seeds_ccsn, counts_ccsn = np.unique(seeds_SN[mask_ccsn], return_counts=True) 
seeds_double_ccsn = seeds_ccsn[counts_ccsn == 2]                              # Seeds with 2 CCSNe will have a counts_ccsn value of 2
mask_SP_double_ccsn = np.isin(seeds_SP, seeds_double_ccsn)

# Determine which systems end up as merging DCOs
mask_merges_hubble_time = DCs['Merges_Hubble_Time'][()] == 1
mask_SP_hubble_time_mergers = np.isin(seeds_SP, seeds_DC[mask_merges_hubble_time]) 

# Combine the 2 masks to get only the correct subset of systems 
combined_mask = mask_SP_double_ccsn & mask_SP_hubble_time_mergers
m1_zams_masked = m1_zams[combined_mask]

print(f"Primary ZAMS masses for the {np.sum(combined_mask)} system(s) which experience 2 CCSNe before becoming a merging DCO:\n {m1_zams_masked}")

# %%
# Always remember to close your data file
data.close()

# %% [markdown]
# ## 3. Visualizing the data
#
# Although math is the fundamental basis of physics and astrophysics, we cannot always easily convert numbers and equations into a coherent picture. Plotting is therefore a vital tool in bridging the gap between raw data and a deeper scientific understanding. 
#
# *Disclaimer:*
#
# There are many ways to make the same plot in matplotlib and there are many ways to bin your data. Often, there is no "best" way to display data in a plot, and the message conveyed can be heavily dependent on the context of the data as well as asthetic plotting decisions.
#
# For example, in histograms, as we discuss below, the relatively subjective choice of bin size can significantly affect the interpretation of the results. It is important to be aware of when and how we make these choices and to try to reduce any unintended bias.
#

# %% [markdown]
#
#
# ### Example: inspect the chirp mass and component masses of Double Compact Objects
#
# In the example below, we use the following conventions:
#
# 1 - We deliberately choose to use the matplotlib.pyplot.subplots routine even when creating a single figure (as opposed to using pyplot.plot). This is because many online forums (e.g Stackoverflow) use this syntax. Furthermore, this means you do not have to learn two different types of syntax when creating either a single or multiple panel figure.
#
# 2 - We choose to do the binning within the numpy/array environment instead of with inbuilt functions such as plt.hist / axes.hist. The reason is that you have more control over what you do, such as custom normalization (using rates, weights, pdf, etc.). It also forces you to have a deeper understanding of what you are calculating, and allows you to check intermediate steps with print statements.  Once you know how to bin your data this way you can also easily expand these routines for more complicated plots (2D binning).

# %% [markdown]
# **Note:** for this exercise, we recommend running your own simulation of at least 100,000 binaries in order to have a sufficient number of DCOs to have an interesting plot. We use the default tutorial data here for illustrative purposes, and because such a file is too large to store on github.

# %% [markdown]
# ### Get some data to plot

# %%
path_to_data = 'COMPAS_Tutorial_Output.h5'

data  = h5.File(path_to_data)
print(list(data.keys()))

DCs = data['BSE_Double_Compact_Objects']

m1 = DCs['Mass(1)'][()]
m2 = DCs['Mass(2)'][()]
def calculate_chirp_mass(m1, m2):
    m_chirp = (m1 * m2)**(3/5) / (m1 + m2)**(1/5)
    return m_chirp
m_chirp = calculate_chirp_mass(m1, m2)

data.close()

# %% [markdown]
# ### Plot histogram and CDF of data on left, and component mass scatter plot on the right

# %%
fig, axes = plt.subplots(ncols=2, figsize=(12,5))
fs_title = 30
fs_label= 20
fs_tick = 15
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.sans-serif": "Times New Roman",
})

# Histogram
ax = axes[0]
bins = np.linspace(0, max(m_chirp), 21) # use 20 bins, up to the maximum chirp mass
ax.hist(m_chirp, bins=bins, color='black', alpha=0.3)
ax.set_ylabel('Histogram counts', fontsize=fs_label)
ax.set_xlabel(r'$\mathcal{M} \; [M_\odot$]', fontsize=fs_label)
ax.set_title('Chirp mass distribution', fontsize=fs_title)
ax.tick_params(axis='both', which='major', labelsize=fs_tick)

# CDF
ax = ax.twinx()
cdf_x_values = np.sort(m_chirp)
cdf_y_values = np.cumsum(cdf_x_values)
cdf_y_values /= cdf_x_values[-1]
np.insert(cdf_x_values, 0, 0) # insert a 0 at the front of the array
cdf_y_values = np.linspace(0, ax.get_ylim()[1], len(cdf_x_values))
ax.plot(cdf_x_values, cdf_y_values, 'r')
ax.set_ylabel('CDF values', rotation=270, labelpad=25, fontsize=fs_label, color='red')
ax.set_ylim(0, 1)
ax.tick_params(axis='y', which='major', labelsize=fs_tick, labelcolor='red')
ax.grid(False)

# Scatter plot 
ax = axes[1]
ax.scatter(m1, m2, color='black')
ax.set_title('Component Masses', fontsize=fs_title)
ax.set_xlabel(r'$M_1 \; [M_\odot$]', fontsize=fs_label)
ax.set_ylabel(r'$M_2 \; [M_\odot$]', fontsize=fs_label)
ax.tick_params(axis='both', which='major', labelsize=fs_tick)

fig.tight_layout()

# %%

# %%
