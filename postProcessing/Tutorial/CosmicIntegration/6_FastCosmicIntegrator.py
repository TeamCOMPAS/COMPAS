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
# This last notebook is meant to introduce the "FastCosmicIntegration.py"
#
# This script works in a very similar way as the other cosmic integrators, but it is optimised to be run from the terminal, and thus can be used to easily run on HPC and as part of grid runs. 
#
# For now this only includes one variation of dP/dZ, namely a skewed log-normal distribution for which the default values result in the log-normal distribution from Neijssel et al. 2019.
# Following Neijssel et al. 2019, the SFR(z) for this function follows the functional form from Madau & Dickinson 2014, but allows the user to manually adjust the parameters of this function
#
#

# # Path definitions

# +
import os

pathNoteBook    = os.getcwd()
pathScripts     = pathNoteBook + '/PythonScripts/'

print(pathScripts)

pathData        = '/Users/lieke/surfdrive/Documents/test_CI/COMPAS_Output/'
# -

# # Imports

# +
import numpy as np
import sys
import matplotlib.pyplot as plt
import astropy.units as u

#custom scripts
sys.path.append(pathScripts)
import FastCosmicIntegration as CI #this imports other routines by itself
import ClassCOMPAS
import selection_effects
import warnings

# To make all your plots look nice
from matplotlib import rc
import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

plt.rc('font', family='serif')

# -

# # The magic happens in find_detection_rate()
#
# The most important function in FasCosmicIntegration.py is find_detection_rate
# We will spend most of this notebook breaking down this function.
#
#
# When you run FasCosmicIntegration.py from your terminal, it is run with a lot of default parameters. 
# These parameters can be changed from your terminal run. 
#
# E.g. if you would like to know the rate for BHNS in stead of for BBHs (which is default) you would run
#
#  ``` python  FasCosmicIntegration.py --dco_type 'BHNS' ```

# +
# First define the parameters
path            = pathData
filename        ="COMPAS_Output.h5"

# For what DCO would you like the rate?  options: ALL, BHBH, BHNS NSNS
dco_type        ="BBH"
weight_column   =None
                        
merges_hubble_time     = True
pessimistic_CEE        = True
no_RLOF_after_CEE      = True

# Options for the redshift evolution 
max_redshift           = 10.0
max_redshift_detection = 2.0
redshift_step          = 0.001
z_first_SF             = 10

# Metallicity of the Universe
min_logZ               =-12.0 
max_logZ               =0.0 
step_logZ              =0.01

#and detector sensitivity
sensitivity            ="O1" 
snr_threshold          =8 

Mc_max                 =300.0 
Mc_step                =0.1 
eta_max                =0.25 
eta_step               =0.01
snr_max                =1000.0 
snr_step               =0.1


# Parameters to calculate the representing SF mass (make sure these match YOUR simulation!)
m1_min          = 15 * u.Msun 
m1_max          = 150 * u.Msun
m2_min          = 0.1 * u.Msun
fbin            = 0.7

# Parameters determining dP/dZ and SFR(z), default options from Neijssel 2019
aSF             = 0.01 
bSF             = 2.77 
cSF             = 2.90 
dSF             = 4.70
#
mu0             =0.035 
muz             =-0.23 
sigma0          =0.39
sigmaz          =0. 
alpha           =0.0 





# -

# # This means you can do all the heavy lifting at once by running find_detection_rate()!
#

# +
detection_rate, formation_rate, merger_rate, redshifts, COMPAS = CI.find_detection_rate(pathData, filename="COMPAS_Output.h5", dco_type="BBH", weight_column=None,
                        merges_hubble_time=True, pessimistic_CEE=True, no_RLOF_after_CEE=True,
                        max_redshift=10.0, max_redshift_detection=1.0, redshift_step=0.001, z_first_SF = 10,
                        m1_min=5 * u.Msun, m1_max=150 * u.Msun, m2_min=0.1 * u.Msun, fbin=0.7,
                        aSF = 0.01, bSF = 2.77, cSF = 2.90, dSF = 4.70,
                        mu0=0.035, muz=-0.23, sigma0=0.39,sigmaz=0., alpha=0.0, 
                        min_logZ=-12.0, max_logZ=0.0, step_logZ=0.01,
                        sensitivity="O1", snr_threshold=8, 
                        Mc_max=300.0, Mc_step=0.1, eta_max=0.25, eta_step=0.01,
                        snr_max=1000.0, snr_step=0.1)


# -

print(detection_rate)



# # We will now spend the rest of this notebook breaking this function down piece by piece, to actually understand what is going on :) 
#
# ## The function starts by checking the values you have supplied

# +
# assert that input will not produce errors
assert max_redshift_detection <= max_redshift, "Maximum detection redshift cannot be below maximum redshift"
assert m1_min <= m1_max, "Minimum sampled primary mass cannot be above maximum sampled primary mass"
assert np.logical_and(fbin >= 0.0, fbin <= 1.0), "Binary fraction must be between 0 and 1"
assert Mc_step < Mc_max, "Chirp mass step size must be less than maximum chirp mass"
assert eta_step < eta_max, "Symmetric mass ratio step size must be less than maximum symmetric mass ratio"
assert snr_step < snr_max, "SNR step size must be less than maximum SNR"

nonnegative_args = [(max_redshift, "max_redshift"), (max_redshift_detection, "max_redshift_detection"), (m1_min.value, "m1_min"), (m1_max.value, "m1_max"),
                    (m2_min.value, "m2_min"), (mu0, "mu0"), (sigma0, "sigma0"),  
                    (step_logZ, "step_logZ"), (snr_threshold, "snr_threshold"), (Mc_max, "Mc_max"),
                    (Mc_step, "Mc_step"), (eta_max, "eta_max"), (eta_step, "eta_step"), (snr_max, "snr_max"), (snr_step, "snr_step")]


for arg, arg_str in nonnegative_args:
    assert arg >= 0.0, "{} must be nonnegative".format(arg_str)

# warn if input is not advisable
if redshift_step > max_redshift_detection:
    warnings.warn("Redshift step is greater than maximum detection redshift", stacklevel=2)
if Mc_step > 1.0:
    warnings.warn("Chirp mass step is greater than 1.0, large step sizes can produce unpredictable results", stacklevel=2)
if eta_step > 0.1:
    warnings.warn("Symmetric mass ratio step is greater than 0.1, large step sizes can produce unpredictable results", stacklevel=2)
if snr_step > 1.0:
    warnings.warn("SNR step is greater than 1.0, large step sizes can produce unpredictable results", stacklevel=2)


# -

# # Use functions in ClassCOMPAS to read your data

# start by getting the necessary data from the COMPAS file
COMPAS = ClassCOMPAS.COMPASData(path, fileName=filename, Mlower=m1_min, Mupper=m1_max, m2_min=m2_min, binaryFraction=fbin, suppress_reminder=True)
COMPAS.setCOMPASDCOmask(types=dco_type, withinHubbleTime=merges_hubble_time, pessimistic=pessimistic_CEE, noRLOFafterCEE=no_RLOF_after_CEE)
COMPAS.setCOMPASData()
COMPAS.set_sw_weights(weight_column)
COMPAS.find_star_forming_mass_per_binary_sampling()


# +
# COMPAS now contains the data from your hdf5 file:
print('DCO mask', COMPAS.DCOmask)

print('Primary masses %s, \nSecondary masses %s '%(COMPAS.mass1, COMPAS.mass2) )

# Hint, you can see all the options in COMPAS by typing COMPAS. and hitting tab to use tab-complete :)


# -

# ## Compute some more useful values

# compute the chirp masses and symmetric mass ratios only for systems of interest
chirp_masses = (COMPAS.mass1*COMPAS.mass2)**(3/5) / (COMPAS.mass1 + COMPAS.mass2)**(1/5)
etas = COMPAS.mass1 * COMPAS.mass2 / (COMPAS.mass1 + COMPAS.mass2)**2
n_binaries = len(chirp_masses)
# another warning on poor input
if max(chirp_masses)*(1+max_redshift_detection) < Mc_max:
    warnings.warn("Maximum chirp mass used for detectability calculation is below maximum binary chirp mass * (1+maximum redshift for detectability calculation)", stacklevel=2)


# ## Now compute the redshift parameters that you will use for the cosmic integration
# mostly a list of redsifts and their corresponding cosmological times

# calculate the redshifts array and its equivalents
redshifts, n_redshifts_detection, times, time_first_SF, distances, shell_volumes = CI.calculate_redshift_related_params(max_redshift, max_redshift_detection, redshift_step, z_first_SF)


print('redshifts', redshifts, '\ntime_first_SF', time_first_SF, '\nshell_volumes', shell_volumes)

# ## compute the SFR(z)
# Following the functional form of Madau & Dickinson (defined by the parameters a,b,c,d), we compute the amount of stellar mass formed at each of the redshifts we supplied
#
# ## And convert this to the number of stars that we needed to form
# By dividing the ```SFR(z) [Msun/Gpc^-3]``` by the star forming mass needed to get your simulation ```(COMPAS.mass_evolved_per_binary.value * COMPAS.n_systems)```, we basically rescale our simulation to represent the number of stars formed at each redshift.
#
#

# find the star forming mass per year per Gpc^3 and convert to total number formed per year per Gpc^3
sfr = CI.find_sfr(redshifts, a = aSF, b = bSF, c = cSF, d = dSF) # functional form from Madau & Dickinson 2014
n_formed = sfr / (COMPAS.mass_evolved_per_binary.value * COMPAS.n_systems) # Divide the star formation rate density by the representative SF mass


print('sfr', sfr, '[$\mathrm{M_{\odot} yr^{-1} Gpc^{-3}}$]')
print('Star forming mass needed to get your simulation:', (COMPAS.mass_evolved_per_binary.value * COMPAS.n_systems), '[$\mathrm{M_{\odot}$]')
print('Number formation rate', n_formed, '[$\mathrm{yr^{-1} Gpc^{-3}}$]')


# ## Get your metallicity density distribution (dP/dZ) at each redshift
#
# This assumes a skewed-log-normal distribution for metallicities at each redshift. The shape of this distribution is controlled by the parameters:
# ```mu0, muz, sigma_0, sigma_z, alpha```, which you can set from the terminal flags 

# +
# work out the metallicity distribution at each redshift and probability of drawing each metallicity in COMPAS
dPdlogZ, metallicities, p_draw_metallicity = CI.find_metallicity_distribution(redshifts, min_logZ_COMPAS = np.log(np.min(COMPAS.initialZ)),
                                                                            max_logZ_COMPAS = np.log(np.max(COMPAS.initialZ)),
                                                                            mu0=mu0, muz=muz, sigma_0=sigma0, sigma_z=sigmaz, alpha = alpha,
                                                                            min_logZ=min_logZ, max_logZ=max_logZ, step_logZ = step_logZ)


# -

print('dPdlogZ=%s, \nmetallicities=%s, \np_draw_metallicity = %s)'%(dPdlogZ, metallicities, p_draw_metallicity) )
# shape dPdlogZ is 
print(np.shape(dPdlogZ))


# ## Do the actual integration
# We are now going to place our DCO systems at each of the chosen redshifts, and check if they have enough time to merge. Together with it's metallicity, the SFR(z) and dP/dZ and the representative star forming mass (n_formed) we can calculate a merger rate!
#

# calculate the formation and merger rates using what we computed above
formation_rate, merger_rate = CI.find_formation_and_merger_rates(n_binaries, redshifts, times, time_first_SF, n_formed, dPdlogZ,
                                                                metallicities, p_draw_metallicity, COMPAS.metallicitySystems,
                                                                COMPAS.delayTimes, COMPAS.sw_weights)


# ## Calculate detection probability
#
# Gravitational wave detectors are not perfect. And since heavy objects make louder gravitational waves, we will see them from farther away. With this function we will compute the probability of detecting each system, which will give us a 'detection_rate'

# +

# create lookup tables for the SNR at 1Mpc as a function of the masses and the probability of detection as a function of SNR
snr_grid_at_1Mpc, detection_probability_from_snr = CI.compute_snr_and_detection_grids(sensitivity, snr_threshold, Mc_max, Mc_step,
                                                                                eta_max, eta_step, snr_max, snr_step)

# use lookup tables to find the probability of detecting each binary at each redshift
detection_probability = CI.find_detection_probability(chirp_masses, etas, redshifts, distances, n_redshifts_detection, n_binaries,
                                                    snr_grid_at_1Mpc, detection_probability_from_snr, Mc_step, eta_step, snr_step)

# finally, compute the detection rate using Neijssel+19 Eq. 2
detection_rate = np.zeros(shape=(n_binaries, n_redshifts_detection))
detection_rate = merger_rate[:, :n_redshifts_detection] * detection_probability \
                * shell_volumes[:n_redshifts_detection] / (1 + redshifts[:n_redshifts_detection])

# -

# # You are done! :D 
#
# The next step in ```FastCosmicItegration.py``` will append your newly calculated rates to the COMPAS_output.hdf5 file. This happens in ```append_rates()```. Because appending rates could lead to a data heavy file, it might be useful to only append your rates binned by redshit. For this purpose you can set ```append_binned_by_z = True``` in  ```append_rates()```
#

# # Now let's go plot your results!
#
# The function ```CI.plot_rates()``` will do the same as what we are going to do in the cells below. 

# +
chirp_masses = (COMPAS.mass1*COMPAS.mass2)**(3./5.) / (COMPAS.mass1 + COMPAS.mass2)**(1./5.)


# sum things up across binaries
total_formation_rate = np.sum(formation_rate, axis=0)
total_merger_rate = np.sum(merger_rate, axis=0)
total_detection_rate = np.sum(detection_rate, axis=0)

# and across redshifts
cumulative_detection_rate = np.cumsum(total_detection_rate)
detection_rate_by_binary = np.sum(detection_rate, axis=1)





# +
###########################
#Start plotting

# set some constants for the plots
plt.rc('font', family='serif')
fs = 20
lw = 3

fig, axes = plt.subplots(2, 2, figsize=(20, 20))

axes[0,0].plot(redshifts, total_formation_rate, lw=lw)
axes[0,0].set_xlabel('Redshift', fontsize=fs)
axes[0,0].set_ylabel(r'Formation rate $[\rm \frac{\mathrm{d}N}{\mathrm{d}Gpc^3 \mathrm{d}yr}]$', fontsize=fs)

axes[0,1].plot(redshifts, total_merger_rate, lw=lw)
axes[0,1].set_xlabel('Redshift', fontsize=fs)
axes[0,1].set_ylabel(r'Merger rate $[\rm \frac{\mathrm{d}N}{\mathrm{d}Gpc^3 \mathrm{d}yr}]$', fontsize=fs)

axes[1,0].plot(redshifts[:len(cumulative_detection_rate)], cumulative_detection_rate, lw=lw)
axes[1,0].set_xlabel('Redshift', fontsize=fs)
axes[1,0].set_ylabel(r'Cumulative detection rate $[\rm \frac{\mathrm{d}N}{\mathrm{d}yr}]$', fontsize=fs)

axes[1,1].hist(chirp_masses, weights=detection_rate_by_binary, bins=25, range=(0, 50))
axes[1,1].set_xlabel(r'Chirp mass, $\mathcal{M}_c$', fontsize=fs)
axes[1,1].set_ylabel(r'Mass distrbution of detections $[\rm \frac{\mathrm{d}N}{\mathrm{d}\mathcal{M}_c \mathrm{d}yr}]$', fontsize=fs)

#########################
#Plotvalues

# Add text upper left corner
axes[0,0].text(0.05,0.8, "mu0=%s \nmuz=%s \nsigma0=%s \nsigmaz=%s \nalpha=%s"%(mu0,muz,sigma0,sigmaz,alpha), transform=axes[0,0].transAxes, size = fs) 

for ax in axes.flatten():
    ax.tick_params(labelsize=0.9*fs)

# Save and show :)
plt.savefig(pathData +'Rate_Info'+"mu0%s_muz%s_alpha%s_sigma0%s_sigmaz%s"%(mu0,muz,alpha,sigma0, sigmaz)+'.png', bbox_inches='tight') 
plt.show()

# -

# # One more plot for the road
#
# Because the test sample used in this notebook is quite small, the detected chirp mass distribution does not look great. 
# Below is an example plotting the same distribution, but adding a KDE to look at a more smooth distribution
#
#

# +
from scipy import stats

#########################
# Get the Hist    
bins = np.arange(0,50, 2)
hist, bin_edge = np.histogram(chirp_masses, weights = detection_rate_by_binary, bins=bins)
center_bins = (bin_edge[:-1] + bin_edge[1:])/2.
# And the KDE
kernel = stats.gaussian_kde(chirp_masses, bw_method='scott', weights=detection_rate_by_binary)
binwidth = np.diff(bin_edge) 

colors = ['blue']
###########################
#Start plotting
fig, ax = plt.subplots(figsize=(12,8))

########################
# Plot the Hist    
norm = 1.
ax.bar(center_bins, hist/norm, width= np.diff(bins), 
       alpha=1.0, fill=False, edgecolor=colors[0],lw = 1., zorder = 0) 

########################
# Add KDE
x_KDE = np.arange(bins[0],bins[-1],0.1)
KDEy_vals = kernel(x_KDE)*sum(hist)*np.diff(bins)[0]/norm #re-normalize the KDE
ax.plot(x_KDE, KDEy_vals, lw=5, color=colors[0], zorder =1,label = 'KDE of Mchirp')
ax.fill_between(x_KDE, y1=0, y2=KDEy_vals, color=colors[0], alpha = 0.05, zorder = 1)

########################
#Plotvalues
########################
# Add text upper left corner
ax.text(0.05,0.8, "mu0=%s \nmuz=%s \nsigma0=%s \nsigmaz=%s \nalpha=%s"%(mu0,muz,sigma0,sigmaz,alpha), 
        transform=ax.transAxes, size = 15) 

# ax.hist(chirp_masses, weights=detection_rate_by_binary, bins=25, range=(0, 50))
ax.set_xlabel(r'Chirp mass, $\mathcal{M}_c$ $\mathrm{[M_{\odot}]}$', fontsize=30)
ax.set_ylabel(r'Mass distrbution of detections $[\rm \frac{\mathrm{d}N}{\mathrm{d}\mathcal{M}_c \mathrm{d}yr}]$', fontsize=30)
plt.rc('xtick', labelsize=25)    # fontsize of the tick labels
plt.rc('ytick', labelsize=25)    # fontsize of the tick labels

# Save and show :)
plt.savefig(pathData +'Rate_Info'+"mu0%s_muz%s_alpha%s_sigma0%s_sigmaz%s"%(mu0,muz,alpha,sigma0, sigmaz)+'.png', bbox_inches='tight') 
plt.show()

# -


