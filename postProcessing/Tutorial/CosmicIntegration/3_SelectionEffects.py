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
# The selection effects assign a probability of observing a system.
#
# We determine if we can confidently select a system by choosing a signal-to-noise ratio (SNR), for which we often use 8. 
#
# The SNR depends on the individual component masses, the distance and the orientation of the binary compared to the gravitational wave detector.  
#
# By sampling uniformly over the possible orientations of the system
# we can assign the fraction of the time the system at that distance can be observed. This fraction is the probability.
#
#
# The SNR also depends on the sensitivity of the detector
#
# combined we need to pass to the function
#
# detection_probability(m1, m2, redshift, distance, snr_threshold,sensitivity='design')
#
# If you use this pipeline we would appreciate it if you cite
# Selection effects   ; https://arxiv.org/pdf/1711.06287
#

# # Paths

# +
import os

pathNoteBook     = os.getcwd()
pathClassCOMPAS  = pathNoteBook + '/PythonScripts/'


# -

# # Imports

import numpy as np
import sys
#custom scripts
sys.path.append(pathClassCOMPAS)
import selection_effects

# # Quick example

m1 = 40 #Msun
m2 = 40 #Msun
redshift = 0.1
distance = 463.4  # Mpc quick estimate for illustration purpuses
                  # code uses astropy toconvert
                  # redshift to luminosity distance
snr_threshold = 8
sensitivity = 'O1'

P = selection_effects.detection_probability(m1, m2, redshift, distance, snr_threshold,sensitivity=sensitivity)

print(P)


