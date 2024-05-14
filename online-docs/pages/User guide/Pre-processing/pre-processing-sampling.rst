Sampling in COMPAS
==================

COMPAS has built-in functionality for sampling using simple distributions for initial masses,
binary mass ratios, semi-major axes, and eccentricities, and metallicities.  

More complex distributions can be sampled from with external scripts which can then generate
a suitable gridfile.  An example is the Moe & DiStefano (2017) distribution for 
the initial parameters of a binary that couples the two masses, orbital period, and eccentricity.

The sampler script for the Moe & DiStefano (2017) distribution
can be found in ``compas_python_utils/preprocessing/sampleMoeDiStefano.py``. As described in the 
header, only the number of systems and upper and lower mass bounds can be
set (the parameter correlations break if you try to set the other bounds).
These values are set at the bottom of the script.

