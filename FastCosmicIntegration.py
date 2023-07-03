#! /usr/bin/env python3
import os
import sys
import time
import psutil
import math
import io
import numpy as np
import h5py
import scipy.interpolate
from scipy.interpolate import interp1d
from scipy.stats import norm as NormDist
import astropy.units as units
from astropy.cosmology import WMAP9 as cosmology
from datetime import date


"""

OBJECT-ORIENTED CODE
--------------------

The code here was written primarily to enable the use of the (original) FastCosmicIntegration code
(written by Ilya Mandel), with a few changes. The SelectionEffects and COMPAS classes were (re)written
here just enough to allow the (modified) FastCosmicIntegration code to run (now encapsulated in its own
class - the CosmicIntegration class) - some of the original code from those the SelectionEffects and
COMPAS classes has been omitted here.  The omitted functionality can mostly be provided in another 
script that would use this one to generate the data (plotting, writing of rates back to the COMPAS HDF5
file, etc.).  The one piece of functionality not included here that would have to be included in this
file is the Stroopwaffel weighting.  We can add that in here when Stroopwaffel releases and we know what
we need to do.

The code in this file is as object-oriented as I could make it (or, perhaps better, as object-oriented
as I wanted to make it...), and OO principles should be used when using the code.  Some functions and
variables in the classes defined in this file are declared private - though "private" is really just a
convention in Python, and "private" functions and variables can still be accessed outside the class if
you're savvy enough...  In keeping with OO principles, the functions and variables declared as private
should not be accessed outside the class, and if they are results are not guaranteed to be consistent.

While, as noted above, this code was written primarily to facilitate a modified version of the
FastCosmicIntegration code, the other classes (SelectionEffects, COMPAS) can be instantiated and used
independently.

The way I've implemented each class is that the class is initialised when constructed, with the
parameters passed to the constructor initialising class member variables.  The parameters to each of
the class constrcutors will default to values set as constants at the top of the file.  Each class
exposes a number of public functions, including getters.  The philosophy I used is:

    . all class member functions and getters can be called using either keyword or positional parameters.
      
    . when calling functions and getters using keyword parameters, any function parameters omitted from
      the call will be set to the current value of the corresponding class member variable. 
      
    . when calling functions and getters using positional parameters, any parameters passed with a value
      of 'None' will be set to the current value of the corresponding class member variable. 
     
    . a successful call to any class member function or getter that takes parameters that were also
      parameters to the constructor will cause the class member variables corresponding to those
      parameters to be updated to the values passed to the public function - and any class member
      variables calculated using those parameters will reatain their values.

    . if a call to a function fails (i.e. returns an error), the class member variables are not updated
      as described above.
      
Implementing this philosophy means that callers can instantiate and initialise a class, then call a public
function a number of times, with each successive call to the public function (or any other function) changing
corresponding class member variables - and so changing the default values to be used in succesive calls
to functions if parameters are not specified (or specified as 'None').  There are a couple of minor
variations - namely 'p_McBins', and 'p_CreatorDeatils' - these are documented in the code.

An example run with multiple calls to CosmicIntegration::FindRates() might look like:

    import FastCosmicIntegration as fCI
    
    McBins, errStr = fCI.Utils.FixedWidthBins(5.0, 50.0, 0.1)
    if errStr is None:

        CI = fCI.CosmicIntegration(p_NoiseFilepath        = './noise_files',
                                   p_NoiseFilename        = 'SNR_Grid_IMRPhenomPv2_FD_all_noise.hdf5',
                                   p_Sensitivity          = 'O3', 
                                   p_SNRthreshold         = 8.0, 
                                   p_COMPASfilepath       = '/e/compas/InferenceProject/newCOMPASruns', 
                                   p_COMPASfilename       = 'h5out.h5', 
                                   p_MaxRedshift          = 10.0, 
                                   p_MaxRedshiftDetection = 1.0, 
                                   p_RedshiftStep         = 0.001,
                                   p_DCOtype              = 'BBH',
                                   p_BinaryFraction       = 0.7,
                                   p_RandomSeed           = 0)

        rates, files, errStr = CI.FindRates(p_MU_Z = -0.23, p_SIGMA_0 = 0.39, p_SFR_a = 0.01, p_SFR_d = 4.7, \
                                            p_BinaryFraction = 0.7, p_McBins = None, p_PerBinaryFiles = False, \
                                            p_PerBinaryFilesPath = None, p_CreatorDetails = 'Fred Flintstone,fred@bedrock.com', \
                                            p_UseSampledMassRanges = False)
        if errStr is None:

            rates, files, errStr = CI.FindRates(p_BinaryFraction = 1.0)
            if errStr is None:

                rates, files, errStr = CI.FindRates(p_BinaryFraction = 1.0, p_RedshiftStep = 0.01)
                if errStr is None:

                    rates, files, errStr = CI.FindRates(p_BinaryFraction = 0.7, p_RedshiftStep = 0.1, p_PerBinaryFiles = True)
                    if errStr is None:

                        rates, files, errStr = CI.FindRates(p_BinaryFraction = 0.7, p_RedshiftStep = 0.1, p_PerBinaryFiles = False, p_Verbose = False)

                        ...
                        ...


The CosmicIntegration class exposes the SelectionEffects and COMPAS objects instantiated, so users can
call the public functions and getters from those classes if required.  e.g. in the code snippet above,
after instantiating the CosmicIntegration class, the user might call COMPAS getters:

    massEvolvedPerBinary = CI.COMPAS().MassEvolvedPerBinary()

    or, for a getter that takes an index parameter:

        chirpMass, errStr = CI.COMPAS().ChirpMasses()(idx), where idx is an interger index into the chirp masses array

    or, alternatively, if the entire chirp mass array is required:

        chirpMasses = CI.COMPAS().ChirpMasses()[0]

Or the SelectionEffects getters:

    sensitivity = CI.SelectionEffects().Sensitivity()

Note that some, but not all, public functions and getters return an error string in addition to the value
requested.


ROBUSTNESS
----------

This code is as robust as I could make it.  The aim is for the code to run to completetion and return
a result to the caller rather than halt execution if an error occurs.  If an error occurs that prevents
the requested data being returned to the caller, an error message will be returned that the caller can
display to the user to help identify the problem.  All functions implemented in this code, except for
the class constructors, return an error string that will be 'None' if no error occurred, otherwise the
content of the string will describe the error.  Class constructors do not return values to the caller,
so can't implement this functionality: instead, the class constructors will raise an exception if an
error occurs and execution will be halted.  Since the class constructors are called early in the
execution path, halting execution if they fail shoukdn't cause too much of an issue.


MEMORY CONSUMPTION
------------------

The main use of this code is to calculate formation, merger, and detection rates for a poulation of
binaries, and return those rates, in arrays, to the caller.  Depending upon the number of binaries in
the population, the rates arrays can consume a great deal of memory - indeed, the rates arrays could
consume more memory than there is physical RAM installed in the system upon which this code is being
executed.  The system may have available swap the OS can use as virtual memory, but even if swap is
configured on an SSD (rather than an HDD) it is still very much slower than RAM access (possibly a couple
or orders of magnitude).  For performance reasosn we choose not to use virtual memory in this code - we
instead check that there is sufficient physivall memory available at execution to to allocate the rates
arrays required.  Because numpy does not allocate the entire array when the array is created, even if
there is sufficient memory when the check is made, that doesn't guarantee that we won't run out - as the
code progresses more of the arrays will be allocated, and if other users/processes have consumed available
memory before we get to allocate it, we could run out of memory - but checking first is useful, especially
if the user has requested an array size that exceeds installed memory.

We could add a commandline option to allow users to choose to use swap, but for now we assume we only want
to use physical memory.


PERFORMANCE
-----------

I really didn't do much in the way of enhancing performance over the original code - just a few minor tweaks.
If performance becomes an issue as we create bigger and bigger data files, we should probably do some code
profiling to see where we might make improvements.  Or maybe we could look at implenting this code in C++.


ASSERT STATEMENTS
-----------------

The python ASSERT statement is used in the original FastCosmicIntegration.py code to check for valid
parameters passed to some functions etc.  The ASSERT statement is a debugging tool for developers, 
primarily used while code is in development, not production.  It really shouldn't be used to check 
for run-time errors that may occur during ordinary production runs.  That doesn't mean it can't be 
used during production of course - and we see that it works ok in the original FastCosmicIntegration.py.
However, if we were to run the original FastCosmicIntegration.py with the -O or -OO commandline options,
the ASSERT statements would be optimised away, and we would expose ourselves to failures that aren't then
caught by ASSERT statements (and if we then tried to duplicate those errors to fix reported problems, we
would presumably run the code without the -O/-OO options and would fail to duplicate the problems...).
The bottom line is that unless we're using ASSERT as instrumentation to find a specific problem, we really
shouldn't have ASSERT statements in production code (and even then they should only be there until the 
problem is found).

Here we set an error string and bubble the error up to the caller where possible.  Where that is not 
possible we use the construct:

   if cond:
       raise exception()

instead of using ASSERT statements.  It's safer (won't be stripped out if -O/-OO is used), and is more 
aligned with good programming practices.

I chose to use the 'RuntimeError' exception throughout just to be consistent - the error string printed
with the exception should explain the problem clearly.

"""


"""
Assumptions made in this code

- We use the KROUPA IMF (Kroup 2021)
- Minimum primary mass must be > 0.5 (KROUPA_BREAK_2)
- Mass ratio in binaries is uniform over [0, 1)
- Metallicities in a COMPAS run are drawn from a flat-in-log distribution

"""



# constants

# miscellaneous
VERBOSE_DEFAULT                = True                                                       # print diagnostics etc. by default

KB                             = int(1024)                                                  # kilobyte
MB                             = KB * KB                                                    # megabyte
GB                             = MB * KB                                                    # gigabyte

RNGType                        = np.random._generator.Generator                             # type for random number generator
RNGType_asSTR                  = 'np.random._generator.Generator'                           # ... as STR

SigFigs4                       = '{:01.4G}'                                                 # floating point format specifier for 4 significant figures
SigFigs                        = SigFigs4                                                   # floating point format specifier used  in writing to per-binary files


# Kroupa IMF is a broken power law with three slopes
# There are two breaks in the Kroupa power law -- they occur here (in solar masses)
# We define upper and lower bounds (KROUPA_BREAK_3 and KROUPA_BREAK_0 respectively)
KROUPA_BREAK_0                 = 0.01
KROUPA_BREAK_1                 = 0.08
KROUPA_BREAK_2                 = 0.5
KROUPA_BREAK_3                 = 200.0

# There are three line segments, each with a specific slope
KROUPA_SLOPE_1                 = 0.3
KROUPA_SLOPE_2                 = 1.3
KROUPA_SLOPE_3                 = 2.3


# Neijssel+19 Eq.7-9
NEIJSSEL_MU_0                  = 0.035                                                      # location (mean in normal) at redshift 0
NEIJSSEL_MU_Z                  = -0.23                                                      # redshift scaling/evolution of the location
NEIJSSEL_SIGMA_0               = 0.39                                                       # scale (variance in normal) at redshift 0
NEIJSSEL_SIGMA_Z               = 0.0                                                        # redshift scaling of the scale (variance in normal)
NEIJSSEL_ALPHA                 = 0.0                                                        # shape (skewness, alpha = 0 retrieves normal dist)
NEIJSSEL_SFR_a                 = 0.01                                                       # scale factor for SFR (see Neijssel+2019)
NEIJSSEL_SFR_b                 = 2.77                                                       # scale factor for SFR (see Neijssel+2019)
NEIJSSEL_SFR_c                 = 2.9                                                        # scale factor for SFR (see Neijssel+2019)
NEIJSSEL_SFR_d                 = 4.7                                                        # scale factor for SFR (see Neijssel+2019)


# selection effects defaults
#
SE_NoiseFilepath               = './'                                                       # noise file path
SE_NoiseFilename               = 'SNR_Grid_IMRPhenomPv2_FD_all_noise.hdf5'                  # noise file name

SE_Dataset_design              = 'SimNoisePSDaLIGODesignSensitivityP1200087'                # noise HDF5 file datset name for LIGO 'design' sensitivity
SE_Dataset_O1                  = 'P1500238_GW150914_H1-GDS-CALIB_STRAIN.txt'                # noise HDF5 file datset name for LIGO 'O1' sensitivity
SE_Dataset_O3                  = 'SimNoisePSDaLIGOMidHighSensitivityP1200087'               # noise HDF5 file datset name for LIGO 'O3' sensitivity
SE_Dataset_ET                  = 'ET_D.txt'                                                 # noise HDF5 file datset name for ET 'ET' sensitivity

SE_Sensitivities               = [ 'design', 'O1', 'O3', 'ET' ]                             # all available detector sensitivities  

SE_sensitivity                 = 'O1'                                                       # detector sensitivity
SE_SNRthreshold                = 8.0                                                        # SNR threshold required for detection

SE_numThetas                   = int(1E6)                                                   # number of thetas to generate
SE_minThetas                   = int(1E4)                                                   # minimum allowable number of thetas

SE_McMax                       = 300.0                                                      # maximum chirp mass in SNR grid
SE_McStep                      = 0.1                                                        # step in chirp mass to use in SNR grid
SE_ETAmax                      = 0.25                                                       # maximum symmetric mass ratio in SNR grid
SE_ETAstep                     = 0.01                                                       # step in symmetric mass ratio to use in SNR grid
SE_SNRmax                      = 1000.0                                                     # maximum SNR in SNR grid
SE_SNRstep                     = 0.1                                                        # step in SNR to use in SNR grid

SE_SNRinterpolatorType         = scipy.interpolate._fitpack2.RectBivariateSpline            # type for SNR interpolator
SE_SNRinterpolatorType_asSTR   = 'scipy.interpolate._fitpack2.RectBivariateSpline'          # ... as STR


# COMPAS defaults
#
COMPAS_MS0                     = 0                                                          # Main Sequence star, initial mass <= 0.7 Msun
COMPAS_MS1                     = 1                                                          # Main Sequence star, initial mass > 0.7 Msun
COMPAS_HG                      = 2                                                          # Hertzsprung Gap star
COMPAS_FGB                     = 3                                                          # First Giant Branch star
COMPAS_CHeB                    = 4                                                          # Core Helium Burning star
COMPAS_EAGB                    = 5                                                          # Early Asymptotic Giant Branch star
COMPAS_TPAG                    = 6                                                          # Thermally Pulsing Asymptotic Giant Branch star
COMPAS_HeMS                    = 7                                                          # Main Sequence Naked Helium star
COMPAS_HeHG                    = 8                                                          # Hertzsprung Gap Naked Helium star
COMPAS_HeGB                    = 9                                                          # Giant Branch Naked Helium star
COMPAS_HeWD                    = 10                                                         # Helium White Dwarf 
COMPAS_COWD                    = 11                                                         # Carbon-Oxygen White Dwarf
COMPAS_ONeWD                   = 12                                                         # Oxygen-Neon White Dwarf
COMPAS_NS                      = 13                                                         # Neutron Star
COMPAS_BH                      = 14                                                         # Black Hole
COMPAS_MR                      = 15                                                         # Massless Remnant
COMPAS_CH                      = 16                                                         # Chemically Homeogeneous star

COMPAS_filepath                = '.'                                                        # data file path
COMPAS_filename                = 'COMPAS_Output.h5'                                         # data file name (expects HDF5 file)
COMPAS_m1Minimum               = 5.0                                                        # minimum value for m1
COMPAS_m1Maximum               = 150.0                                                      # maximum value for m1
COMPAS_m2Minimum               = 0.1                                                        # minimum value for m2
COMPAS_UseSampledMassRanges    = False                                                      # flag to indicate that the mass ranges sampled by COMPAS should be used instead of hard cuts at COMPAS_m1Minimum, COMPAS_m1Maximum, and COMPAS_m2Minimum
COMPAS_binaryFraction          = 1.0                                                        # binary fraction
COMPAS_SFM_popSize             = 20000000                                                   # population size for calculationg start forming mass per binary

COMPAS_flag_withinHubbleTime   = True                                                       # default is only include binaries that merge within a Hubble time 
COMPAS_flag_pessimisticCE      = True                                                       # default is pessimistic CE
COMPAS_flag_noRLOFafterCEE     = True                                                       # default is no RLOF after CE event

# HDF5 file group names
COMPAS_RUN_DETAILS             = 'Run_Details'                                              # Run Details group (file)
COMPAS_BSE_SYSPARMS            = 'BSE_System_Parameters'                                    # BSE System Parameters group (file)
COMPAS_BSE_DCO                 = 'BSE_Double_Compact_Objects'                               # BSE Double Compact Objects group (file)
COMPAS_BSE_CE                  = 'BSE_Common_Envelopes'                                     # BSE Common Envelopes group (file)

# HDF5 file dataset names

COMPAS_seed                    = 'SEED'                                                     # random seed (common dataset name for all groups (files) (except Run_Details))

# Run_Details datasets                                                                      # Note that the COMPAS HDF5 file being used may actually be a file created by
                                                                                            # concatenating many smaller HDF55 files (using e.g. h5copy.py), and if that is
                                                                                            # the case the Run_Details group in the HDF5 file will have an entry for each of
                                                                                            # the smaller files concatenated.  We only read the first entry here to get the
                                                                                            # COMPAS version and run start time and date.  It is possible, though probably
                                                                                            # unlikely, that the user concatenated COMPAS HDF5 files created using different
                                                                                            # versions of COMPAS (maybe h5copy.h5 should check for that...).

COMPAS_version                 = 'COMPAS-Version'                                           # COMPAS version used to create COMPAS HDF5 file
COMPAS_runStart                = 'Run-Start'                                                # start time and date of COMPAS run that created COMPAS HDF5 file

# BSE_SYSPARMS datasets
COMPAS_BSE_SYS_zamsST1         = 'Stellar_Type@ZAMS(1)'                                     # stellar type at ZAMS for star 1
COMPAS_BSE_SYS_zamsST2         = 'Stellar_Type@ZAMS(2)'                                     # stellar type at ZAMS for star 2
COMPAS_BSE_SYS_zamsMass1       = 'Mass@ZAMS(1)'                                             # mass at ZAMS for star 1
COMPAS_BSE_SYS_zamsMass2       = 'Mass@ZAMS(2)'                                             # mass at ZAMS for star 2
COMPAS_BSE_SYS_zamsZ           = 'Metallicity@ZAMS(1)'                                      # metallicity at ZAMS for star 1 (assumed same for star 2)
COMPAS_BSE_SYS_CHE1            = 'CH_on_MS(1)'                                              # star 1 evolved CH on MS flag
COMPAS_BSE_SYS_CHE2            = 'CH_on_MS(2)'                                              # star 2 evolved CH on MS flag

# BSE_DCO datasets
COMPAS_BSE_DCO_ST1             = 'Stellar_Type(1)'                                          # star 1 stellar type
COMPAS_BSE_DCO_ST2             = 'Stellar_Type(2)'                                          # star 2 stellar type
COMPAS_BSE_DCO_Mass1           = 'Mass(1)'                                                  # star 1 mass
COMPAS_BSE_DCO_Mass2           = 'Mass(2)'                                                  # star 2 mass
COMPAS_BSE_DCO_formationTime   = 'Time'                                                     # DCO formation time
COMPAS_BSE_DCO_coalescenceTime = 'Coalescence_Time'                                         # DCO coalescence time
COMPAS_BSE_DCO_mergesHubble    = 'Merges_Hubble_Time'                                       # merges in hubble time flag

# DCO types

COMPAS_DCO_ALL                 = 'ALL'                                                      # all DCO types
COMPAS_DCO_BBH                 = 'BBH'                                                      # Binary Black Hole
COMPAS_DCO_BHNS                = 'BHNS'                                                     # Black Hole + Neutron Star
COMPAS_DCO_BNS                 = 'BNS'                                                      # Binary Neutron Star
COMPAS_DCO_CHE_BBH             = 'CHE_BBH'                                                  # Binary Black Hole formed via CHE
COMPAS_DCO_NON_CHE_BBH         = 'NON_CHE_BBH'                                              # Binary Black Hole formed via NON CHE process

DCO_Types                      = [ COMPAS_DCO_ALL, COMPAS_DCO_BBH, COMPAS_DCO_BHNS, COMPAS_DCO_BNS, COMPAS_DCO_CHE_BBH, COMPAS_DCO_NON_CHE_BBH ] # possible DCO type  

# BSE_CE datasets
COMPAS_BSE_CE_immediateRLOF    = 'Immediate_RLOF>CE'                                        # immedidate RLOF following CE event
COMPAS_BSE_CE_optimisticCE     = 'Optimistic_CE'                                            # optimistice CE

# Metallicity bounds
COMPAS_MINIMUM_METALLICITY     = 0.0001                                                     # Minimum metallicity allowed in COMPAS
COMPAS_MAXIMUM_METALLICITY     = 0.03                                                       # Maximum metallicity allowed in COMPAS
COMPAS_LOG_MINIMUM_METALLICITY = np.log(COMPAS_MINIMUM_METALLICITY)                         # Log(Minimum metallicity allowed in COMPAS)
COMPAS_LOG_MAXIMUM_METALLICITY = np.log(COMPAS_MAXIMUM_METALLICITY)                         # Log(Maximum metallicity allowed in COMPAS)


# cosmic integration defaults
#

CI_PerBinaryFiles              = False                                                      # don't produce per-binaryt mergers and yields files by default
CI_PerBinaryFilesPath          = './'                                                       # default path is current directory

CI_Mergers_FilenamePrefix      = 'COMPAS_Mergers_'                                          # prefix for per-binary mergers file (if required)
CI_Yields_FilenamePrefix       = 'COMPAS_Yields_'                                           # prefix for per-binary yields file (if required)
CI_Meta_FilenamePrefix         = 'COMPAS_Meta_'                                             # prefix for per-binary meta-data file (if required)

CI_maxRedshift                 = 10.0                                                       # maximum formation redshift
CI_maxRedshiftDetection        = 1.0                                                        # maximum detection redshift
CI_redshiftStep                = 0.001                                                      # redshift step for detection matrix

CI_MinLogZ                     = -12.0                                                      # Minimum logZ at which to calculate dPdlogZ (influences normalization)
CI_MaxLogZ                     = 0.0                                                        # Maximum logZ at which to calculate dPdlogZ (influences normalization)
CI_StepLogZ                    = 0.01                                                       # Size of logZ steps to take in finding a Z range


# utils defaults
#

# binning defaults
# ... fixed width chirp mass bins
FIXED_Mc_BINS_MINIMUM_Mc       = 0.0                                                        # Fixed-width chirpmass bins minimum value
FIXED_Mc_BINS_MAXIMUM_Mc       = 125.00                                                     # Fixed-width chirpmass bins maximum value
FIXED_Mc_BINS_WIDTH            = math.exp(0.1)                                              # Fixed-width chirpmass bins bin width
# ... variable width chirp mass bins
VAR_Mc_BINS_MINIMUM_Mc         = 0.0                                                        # Variable-width chirpmass bins minimum value
VAR_Mc_BINS_MAXIMUM_Mc         = 125.00                                                     # Variable-width chirpmass bins maximum value
VAR_Mc_BINS_FIRST_BIN_WIDTH    = 0.5                                                        # Variable-width chirpmass bins first bin width
VAR_Mc_BINS_WIDTH_PERCENT      = 5.0                                                        # Variable-width chirpmass bins subsequent bins width percentage



"""
Class Utils

Just a container for miscellaneous utility functions.

Public interface:

    Utils::FixedWidthBins() : constructs fixed-width bins list according to parameters passed
    Utils::VariableMcBins() : constructs variable-width bins list according to parameters passed
    Utils::CheckBins()      : check validity of bins list passed

"""
class Utils:

    # Binning utilities


    """
    Utils::FixedWidthBins()

    Returns a LIST of fixed-width bin widths and right edges, based on parameters passed.
    Bins are arranged such that right-edge values increase monotonically from first bin to last bin (left to right).

    The LIST returned will contain exactly two elements, each of them also a LIST:

        1. LIST of bin widths (for fixed-width bins, all elements in this list will be the same)
        2. LIST of values at the right edges of the bins

    The right edge of the first (left-most) bin will be p_MinimumValue + p_BinWidth, and subsequent
    bins will have right edges at the previous bin right edge + p_BinWidth.  The last (right-most)
    bin will be first bin whose right edge is >= p_MaximumValue.

    Despite the parameter defaults being default values for chirp mass bins, this is a general purpose
    function that may be used for purposes other than creating chirp mass bins, so bin right edge values
    may be negative.

    Args:
        p_MinimumValue : FLOAT : Value of the left (lower) edge of the first (lowest) bin
                                 Default = FIXED_Mc_BINS_MINIMUM_Mc
        p_MaximumValue : FLOAT : Highest value that must be emcompassed by the last (rightmost) bin
                                 Must be > p_MinimumValue, default = FIXED_Mc_BINS_MAXIMUM_Mc
        p_BinWidth     : FLOAT : Width of each bin
                                 Must be > 0.0, default = FIXED_Mc_BINS_WIDTH    

    Returns:
        1 : LIST of LIST : LIST containing bin width LIST and right-edge LIST
                           'None' if an error occurred.
        2 : STRING       : Error string - set if an error occurred.
                           'None' if no error occurred.    
    """
    @staticmethod
    def FixedWidthBins(p_MinimumValue = FIXED_Mc_BINS_MINIMUM_Mc,
                       p_MaximumValue = FIXED_Mc_BINS_MAXIMUM_Mc,
                       p_BinWidth     = FIXED_Mc_BINS_WIDTH):

        func = 'Utils::FixedWidthBins(): '

        # default return values
        binWidths     = None
        binRightEdges = None
        errStr        = None

        # set local variables from parameter values
        # Doing it this way ensures that we have reasonable values for all parameters and allows users to call the function
        # with positional parameters and just set any parameters they want to default to 'None' rather than specify the value
        minimumValue = FIXED_Mc_BINS_MINIMUM_Mc if p_MinimumValue is None else p_MinimumValue
        maximumValue = FIXED_Mc_BINS_MAXIMUM_Mc if p_MaximumValue is None else p_MaximumValue
        binWidth     = FIXED_Mc_BINS_WIDTH      if p_BinWidth     is None else p_BinWidth

        # check parameter values
        tmpStr = ' must be a floating point number'

        if   not isinstance(minimumValue, float):                                   errStr = func + 'Minimmum value' + tmpStr
        elif not (isinstance(maximumValue, float) and maximumValue > minimumValue): errStr = func + 'Maximum value' + tmpStr + ' > minimmum value ({:.6f})'.format(minimumValue)
        elif not (isinstance(binWidth, float) and binWidth > 0.0):                  errStr = func + 'Bin width' + tmpStr + ' > 0.0'
        else:                                                                                                       # all parameters ok
        
            # first bin left edge is minimumValue, right edge is minimumValue + binWidth
            binRightEdges = [minimumValue + binWidth]
            binWidths     = [binWidth]

            # remaining bins are each binWidth units
            # just keep appending new bins - shouldn't be too onerous here (we don't expect that many bins...)
            thisRightEdge = binRightEdges[0]
            while thisRightEdge < maximumValue:
                thisRightEdge += binWidth
                binRightEdges.append(thisRightEdge)
                binWidths.append(binWidth)
    
        return [binWidths, binRightEdges], errStr


    """
    Utils::VariableWidthBins()

    Returns a LIST of variable-width bin widths and right edges, based on parameters passed.
    Bins are arranged such that right-edge values increase monotonically from first bin to last bin (left to right).

    The LIST returned will contain exactly two elements, each of them also a LIST:

        1. LIST of bin widths
        2. LIST of values at the right edges of the bins

    The right edge of the first (left-most) bin will be p_MinimumValue + p_FirstBinWidth, and subsequent
    bins will have right edges at the previous bin right edge + the calculate bin width.  The last (right-most)
    bin will be first bin whose right edge is >= p_MaximumValue.

    Despite the parameter defaults being default values for chirp mass bins, this is a general purpose
    function that may be used for purposes other than creating chirp mass bins, so bin right edge values
    may be negative.

    Args:
        p_MinimumValue    : FLOAT : Value of the left (lower) edge of the first (lowest) bin
                                    Default = VAR_Mc_BINS_MINIMUM_Mc
        p_MaximumValue    : FLOAT : Highest value that must be emcompassed by the last (rightmost) bin
                                    Must be > p_MinimumValue, default = VAR_Mc_BINS_MAXIMUM_Mc
        p_FirstBinWidth   : FLOAT : Width of first bin
                                    Must be > 0.0, default = VAR_Mc_BINS_FIRST_BIN_WIDTH
        p_BinWidthPercent : FLOAT : Percentage of the bin median value used to calculate bin width
                                    Must be > 0.0 and <= 100.0, default = VAR_Mc_BINS_WIDTH_PERCENT

    Returns:
        1 : LIST of LIST : LIST containing bin width LIST and right-edge LIST
                           'None' if an error occurred.
        2 : STRING       : Error string - set if an error occurred.
                           'None' if no error occurred.    
    """
    @staticmethod
    def VariableWidthBins(p_MinimumValue    = VAR_Mc_BINS_MINIMUM_Mc,
                          p_MaximumValue    = VAR_Mc_BINS_MAXIMUM_Mc,
                          p_FirstBinWidth   = VAR_Mc_BINS_FIRST_BIN_WIDTH,
                          p_BinWidthPercent = VAR_Mc_BINS_WIDTH_PERCENT):

        func = 'Utils::VariableWidthBins(): '

        # default return values
        minimumValue  = None
        binWidths     = None
        binRightEdges = None
        errStr        = None

        # set local variables from parameter values
        # Doing it this way ensures that we have reasonable values for all parameters and allows users to call the function
        # with positional parameters and just set any parameters they want to default to 'None' rather than specify the value
        minimumValue    = VAR_Mc_BINS_MINIMUM_Mc      if p_MinimumValue    is None else p_MinimumValue
        maximumValue    = VAR_Mc_BINS_MAXIMUM_Mc      if p_MaximumValue    is None else p_MaximumValue
        firstBinWidth   = VAR_Mc_BINS_FIRST_BIN_WIDTH if p_FirstBinWidth   is None else p_FirstBinWidth
        binWidthPercent = VAR_Mc_BINS_WIDTH_PERCENT   if p_BinWidthPercent is None else p_BinWidthPercent

        # check parameter values
        tmpStr = ' must be a floating point number'

        if   not isinstance(minimumValue, float):                                   errStr = func + 'Minimmum value' + tmpStr
        elif not (isinstance(maximumValue, float) and maximumValue > minimumValue): errStr = func + 'Maximum value' + tmpStr + ' > minimmum value ({:.6f})'.format(minimumValue)
        elif not (isinstance(firstBinWidth, float) and firstBinWidth > 0.0):        errStr = func + 'First bin width' + tmpStr + ' > 0.0'
        elif not (isinstance(binWidthPercent, float) and \
                  binWidthPercent > 0.0 and binWidthPercent <= 100.0):              errStr = func + 'Bin width percent' + tmpStr + ' > 0.0 and <= 100.0'
        else:                                                                                                       # all parameters ok
        
            # first bin left edge is minimumValue, right edge is minimumValue + firstBinWidth
            binRightEdges = [minimumValue + firstBinWidth]
            binWidths     = [firstBinWidth]

            # remaining bins are each binWidthPercent around a centre value
            # just keep appending new bins - shouldn't be too onerous here (we don't expect that many bins...)
            thisRightEdge = binRightEdges[0]
            while thisRightEdge < maximumValue:
                binLeftEdge   = binRightEdges[len(binRightEdges) - 1]
                thisChirpMass = 100.0 * binLeftEdge / (100.0 - (binWidthPercent / 2.0))
                binHalfWidth  = thisChirpMass - binLeftEdge
                thisRightEdge = thisChirpMass + binHalfWidth
                binRightEdges.append(thisRightEdge)
                binWidths.append(binHalfWidth + binHalfWidth)
    
        return [binWidths, binRightEdges], errStr


    """
    Utils::CheckBins()

    Checks bins passed in for valid bin values.  If p_Bins is not 'None', and is not an empty LIST: 

        1. Must be a LIST of LIST
        2. Must contain exactly 2 elements, both non-empty lists of floating point numbers
        3. If p_AllowNegative = False, values for bin right edgees must be >= 0.0
        4. Values for bin right edges must increase monotonically from lowest to highest (left to right)

    Args:
        p_Bins          : LIST of LIST : LIST containing bin width LIST and right-edge LIST
                                         Default = None
        p_AllowNegative : BOOL         : Flag to indicate whether negative values for bin right edges are allowed
                                         Default = False (bin right edge values must be >= 0.0)

    Returns:
        1 : STRING : Error string - set if an error occurred.
                     'None' if no error occurred.    
    """
    @staticmethod
    def CheckBins(p_Bins = None, p_AllowNegative = False):

        errStr = None

        if   p_Bins is None:                                               errStr = 'No bins specified'
        elif not (isinstance(p_Bins, list) and len(p_Bins) == 2):          errStr = 'Bins must be a LIST containing exactly 2 elements'
        elif not (isinstance(p_Bins[0], list) and len(p_Bins[0]) > 0):     errStr = 'Bin widths must be a non-empty LIST'
        elif not (all(isinstance(elem, float) for elem in p_Bins[0]) and \
                  all(elem > 0.0 for elem in p_Bins[0])):                  errStr = 'Bin widths must be floating point numbers > 0.0'
        elif not (isinstance(p_Bins[1], list) and len(p_Bins[1]) > 0):     errStr = 'Bin right edges must be a non-empty LIST'
        elif not all(isinstance(elem, float) for elem in p_Bins[1]):       errStr = 'Bin right edges must be floating point numbers'
        # if p_AllowNegative is False, check for bin right edges < 0.0
        elif not p_AllowNegative and \
             not (all(elem >= 0.0 for elem in p_Bins[1])):                 errStr = 'Bin right edges must be floating point numbers >= 0.0'
        elif len(p_Bins[0]) != len(p_Bins[1]):                             errStr = 'Bin widths and right edges lists must be same length'
        # check for bin widths monotonically increasing
        # there are many ways to do this, with varying performance impacts
        # this way is easy, and not too onerous for not too many bins
        # we don't expect very many bins, and this will generally only be done once per run
        elif p_Bins[1] != sorted(p_Bins[1]):                               errStr = 'Bin right edges must be monotonically increasing' 

        return errStr


    """
    Utils::IntToLetters()

    Convert integer to string of letters for versioning, thus:

        1  : A
        2  : B
        ...
        26 : Z
        27 : AA
        28 : AB
        ...
        52 : AZ
        53 : BA
        54 : BB
        ...
        702: ZZ
        703: AAA
        704: AAB
        ...
        728: AAZ
        729: ABA
        &30: ABB
        ...
        etc.

    Args:
        p_Value : INT : Integer value to be converted to string of letters
                        Must be > 0, default = 1

    Returns:
        1 : STRING : String of letters corresponding to integer parameter passed
                     'None' if no error occurred.    
        2 : STRING : Error string - set if an error occurred.
                     'None' if no error occurred.    
    """
    @staticmethod
    def IntToLetters(p_Value = 1):

        func = 'Utils::IntToLetters(): '

        # default return values
        resStr = None
        errStr = None

        # check parameter
        value = 1 if p_Value is None else p_Value
        if not (isinstance(value, int) and value > 0): errStr = func + 'Value to be converted must be an integer > 0'
        else:
            value = value - 1

            resStr = ''
            while (value >= 0):   
                resStr = chr(ord('A') + value % 26 ) + resStr
                value  = int(value / 26) - 1
    
        return resStr, errStr


"""
Class SelectionEffects

Provides functionality for calculating SNR values, detection probabilities etc.

Public interface (public functions declared at end of class):

    SelectionEffects() : Class constructor
    Initialise()       : Reinitialises class with new parameters

    Getters:
        DetectionProbabilityFromSNR() : Returns the probability of detecting a CBC with given SNR and threshold (single value at index, or 1D array)
        DetectionProbabilityLen()     : Returns length of detection probability array
        ETAmax()                      : Returns maximum symmetric mass ratio in SNR grid
        ETAstep()                     : Returns step in symmetric mass ratio in SNR grid
        McMax()                       : Returns maximum chirp mass in SNR grid
        McStep()                      : Returns step in chirp mass in SNR grid
        Sensitivity()                 : Returns detector sensitivity being used
        SNRgridAt1Mpc()               : Returns SNR values if binary was at 1Mpc (2D array)
        SNRmax()                      : Returns maximum SNR in detection probability array
        SNRstep()                     : Returns step in SNR in detection probability array
        SNRthreshold()                : Returns SNR threshold being used

Adapted from code written by Sebastian M. Gaebel, last known email: sgaebel@star.sr.bham.ac.uk
"""
class SelectionEffects:

    """
    SelectionEffects::SelectionEffects

    Class constructor

    Initialises class member variables

    Args:
        p_NoiseFilepath : STRING    : Path to SNR noise file
                                      Default = SE_NoiseFilepath
        p_NoiseFilename : STRING    : Name of SNR noise file
                                      Default = SE_NoiseFilename
        p_Sensitivity   : STRING    : Detector sensitivity
                                      Must be one of ['design', 'O1', 'O3', 'ET'], default = SE_sensitivity
        p_MinThetas     : INT       : Minimum number of thetas to generate (thetas are used to calculate the average probability of detecting an event)
                                      Must be positive, default = SE_minThetas
        p_NumThetas     : INT       : Number of thetas to generate
                                      Must be >= p_MinThetas, default = SE_numThetas
        p_SNRthreshold  : FLOAT     : SNR threshold required for detection
                                      Must be >= 0.0, default = SE_SNRthreshold
        p_McMax         : FLOAT     : Maximum chirp mass in SNR grid
                                      Must be positive, default = SE_McMax
        p_McStep        : FLOAT     : Step in chirp mass to use in SNR grid
                                      Must be positive, default = SE_McStep
        p_ETAmax        : FLOAT     : Maximum symmetric mass ratio in SNR grid
                                      Must be positive, default = SE_ETAmax
        p_ETAstep       : FLOAT     : Step in symmetric mass ratio to use in SNR grid
                                      Must be positive, default = SE_ETAstep
        p_SNRmax        : FLOAT     : Maximum SNR in detection probability array
                                      Must be positive, default = SE_SNRmax
        p_SNRstep       : FLOAT     : Step in SNR to use in detection probability array
                                      Must be positive, default = SE_SNRstep
        p_RNG           : Numpy RNG : Numpy random number generator object (type = numpy.random._generator.Generator)
                                      Default = None
        p_Verbose       : BOOL      : Flag to indicate if diagnostics/stats should be printed to console
                                      Default = VERBOSE_DEFAULT

    Errors here are handled by raising an exception - the constructor can't return an error to the caller.
    """
    def __init__(self, 
                 p_NoiseFilepath = SE_NoiseFilepath, 
                 p_NoiseFilename = SE_NoiseFilename, 
                 p_Sensitivity   = SE_sensitivity, 
                 p_MinThetas     = SE_minThetas,
                 p_NumThetas     = SE_numThetas,
                 p_SNRthreshold  = SE_SNRthreshold,
                 p_McMax         = SE_McMax,
                 p_McStep        = SE_McStep,
                 p_ETAmax        = SE_ETAmax,
                 p_ETAstep       = SE_ETAstep,
                 p_SNRmax        = SE_SNRmax,
                 p_SNRstep       = SE_SNRstep,
                 p_RNG           = None,
                 p_Verbose       = VERBOSE_DEFAULT):

        self.__className = 'SelectionEffects'
        func = self.__className + '::SelectionEffects(): '

        self.__initialised = False

        self.__verbose = VERBOSE_DEFAULT if p_Verbose is None else p_Verbose

        if self.__verbose:
            print('...Initialising SelectionEffects')
            mark = time.time()

        # set local variables from parameter values
        # Doing it this way ensures that we have reasonable values for all parameters and allows users to call the function
        # with positional parameters and just set any parameters they want to default to 'None' rather than specify the value
        noiseFilepath = SE_NoiseFilepath if p_NoiseFilepath is None else p_NoiseFilepath
        noiseFilename = SE_NoiseFilename if p_NoiseFilename is None else p_NoiseFilename
        sensitivity   = SE_sensitivity   if p_Sensitivity   is None else p_Sensitivity
        minThetas     = SE_minThetas     if p_MinThetas     is None else p_MinThetas
        numThetas     = SE_numThetas     if p_NumThetas     is None else p_NumThetas
        SNRthreshold  = SE_SNRthreshold  if p_SNRthreshold  is None else p_SNRthreshold
        McMax         = SE_McMax         if p_McMax         is None else p_McMax
        McStep        = SE_McStep        if p_McStep        is None else p_McStep
        ETAmax        = SE_ETAmax        if p_ETAmax        is None else p_ETAmax
        ETAstep       = SE_ETAstep       if p_ETAstep       is None else p_ETAstep
        SNRmax        = SE_SNRmax        if p_SNRmax        is None else p_SNRmax
        SNRstep       = SE_SNRstep       if p_SNRstep       is None else p_SNRstep
        rng           = p_RNG

        # check parameter types and values
        # (only check those required here - other parameters are checked in functions called from here)
        if not (isinstance(minThetas, int) and minThetas > 0):          raise RuntimeError(func + 'Minimum number of thetas must be an integer > 0')
        if not (isinstance(numThetas, int) and numThetas >= minThetas): raise RuntimeError(func + 'Number of thetas must be an integer >= {:d}'.format(minThetas))
        if rng is not None and type(rng) != RNGType:                    raise RuntimeError(func + 'Random number generator must be of type \'' + RNGType_asSTR + '\'')

        # construct a random number generator if necessary - no initial seed
        self.__rng = np.random.default_rng() if rng is None else rng

        # generate thetas - used for probability calculations
        # only need to do this once
        self.__numThetas = numThetas
        self.__GenerateThetas()

        # set remaining class member variables
        self.__noiseFilepath               = None
        self.__noiseFilename               = None
        self.__SNRinterpolator             = None
        self.__SNRgridAt1Mpc               = None
        self.__detectionProbabilityFromSNR = None
        self.__sensitivity                 = None
        self.__SNRthreshold                = None
        self.__McMax                       = None
        self.__McStep                      = None
        self.__ETAmax                      = None
        self.__ETAstep                     = None
        self.__SNRmax                      = None
        self.__SNRstep                     = None

        # initialise class based on parameters passed
        # class can be reinitialised by subsequent calls to __Initialise()
        errStr = self.__Initialise(noiseFilepath, noiseFilename, sensitivity, SNRthreshold, McMax, McStep, ETAmax, ETAstep, SNRmax, SNRstep)
        if errStr is not None: raise RuntimeError(errStr)

        if self.__verbose: 
            timeTaken = time.time() - mark
            memUsage  = psutil.Process().memory_info().rss / GB
            print('   SelectionEffects initialised in {:.2f} seconds.  Process total memory usage = {:.6f} GB'.format(timeTaken, memUsage))


    """
    SelectionEffects::__Initialise()

    Initialises/reinitialises class member variables according to parameters passed.

    Updates the following class member variables if no error occurs:
        __detectionProbabilityFromSNR
        __SNRgridAt1Mpc
        __SNRinterpolator

    Updates the following class member variables to parameter values if no error occurs:
        __ETAmax
        __ETAstep
        __initialised
        __McMax
        __McStep
        __sensitivity
        __SNRmax
        __SNRstep
        __SNRthreshold
        __verbose

    Args:
        p_NoiseFilepath : STRING    : Path to SNR noise file
                                      Default = SE_NoiseFilepath
        p_NoiseFilename : STRING    : Name of SNR noise file
                                      Default = SE_NoiseFilename
        p_Sensitivity   : STRING    : Detector sensitivity
                                      Must be one of ['design', 'O1', 'O3', 'ET'], default = SE_sensitivity
        p_SNRthreshold  : FLOAT     : SNR threshold required for detection
                                      Must be >= 0.0, default = SE_SNRthreshold
        p_McMax         : FLOAT     : Maximum chirp mass in SNR grid
                                      Must be positive, default = SE_McMax
        p_McStep        : FLOAT     : Step in chirp mass to use in SNR grid
                                      Must be positive, default = SE_McStep
        p_ETAmax        : FLOAT     : Maximum symmetric mass ratio in SNR grid
                                      Must be positive, default = SE_ETAmax
        p_ETAstep       : FLOAT     : Step in symmetric mass ratio to use in SNR grid
                                      Must be positive, default = SE_ETAstep
        p_SNRmax        : FLOAT     : Maximum SNR in detection probability array
                                      Must be positive, default = SE_SNRmax
        p_SNRstep       : FLOAT     : Step in SNR to use in detection probability array
                                      Must be positive, default = SE_SNRstep
        p_Verbose       : BOOL      : Flag to indicate if diagnostics/stats should be printed to console
                                      Default = VERBOSE_DEFAULT
    Returns:
        1 : STRING : Error string - set if an error occurred.
                     'None' if no error occurred.    
    """
    def __Initialise(self,
                     p_NoiseFilepath = None, 
                     p_NoiseFilename = None, 
                     p_Sensitivity   = None,
                     p_SNRthreshold  = None,
                     p_McMax         = None,
                     p_McStep        = None,
                     p_ETAmax        = None,
                     p_ETAstep       = None,
                     p_SNRmax        = None,
                     p_SNRstep       = None,
                     p_Verbose       = None): 

        func = self.__className + '::__Initialise(): '

        saveVerbose = self.__verbose                                                                                # in case of error
        self.__verbose = self.__verbose if p_Verbose is None else p_Verbose                                         # reset verbose setting

        # default return values
        errStr = None

        # set local variables from parameter values
        # Doing it this way ensures that we have reasonable values for all parameters and allows users to call the function
        # with positional parameters and just set any parameters they want to default to 'None' rather than specify the value
        noiseFilepath = self.__noiseFilepath if p_NoiseFilepath is None else p_NoiseFilepath
        noiseFilename = self.__noiseFilename if p_NoiseFilename is None else p_NoiseFilename
        sensitivity   = self.__sensitivity   if p_Sensitivity   is None else p_Sensitivity
        SNRthreshold  = self.__SNRthreshold  if p_SNRthreshold  is None else p_SNRthreshold
        McMax         = self.__McMax         if p_McMax         is None else p_McMax
        McStep        = self.__McStep        if p_McStep        is None else p_McStep
        ETAmax        = self.__ETAmax        if p_ETAmax        is None else p_ETAmax
        ETAstep       = self.__ETAstep       if p_ETAstep       is None else p_ETAstep
        SNRmax        = self.__SNRmax        if p_SNRmax        is None else p_SNRmax
        SNRstep       = self.__SNRstep       if p_SNRstep       is None else p_SNRstep
      
        # parameters are checked in called functions...

        SNRinterpolator = self.__SNRinterpolator

        # read noise file and set interpolator
        # only need to read the noise file if this is a reinitialisation and sensitivity has changed
        # (the noice file cannot be changed via reinitialisation - create a new instance if that's required)
        if not self.__initialised or (self.__initialised and sensitivity != self.__sensitivity):
            SNRmassAxis, SNRgrid, errStr = self.__ReadNoiseData(noiseFilepath, noiseFilename, sensitivity)
            if errStr is None:
                # set the interpolator for later use - used for probability calculations
                SNRinterpolator = scipy.interpolate.RectBivariateSpline(np.log(SNRmassAxis), np.log(SNRmassAxis), SNRgrid)
        
        if errStr is None:

            # compute SNR and detection probability grids based on the SNR threshold and grid specifications provided
            # we could check parameters here to see if we really need to do this, but the check would probably take longer
            # than the recalculation...  (well, not quite, but the point is that this is not onerous)
            SNRgridAt1Mpc, detectionProbabilityFromSNR, errStr = self.__ComputeSNRandDetectionGrids(SNRinterpolator, SNRthreshold, McMax, McStep, ETAmax, ETAstep, SNRmax, SNRstep)
            if errStr is None:

                # update class member variables
                self.__SNRinterpolator             = SNRinterpolator
                self.__SNRgridAt1Mpc               = SNRgridAt1Mpc
                self.__detectionProbabilityFromSNR = detectionProbabilityFromSNR

                # update class member variables to match parameters
                self.__sensitivity  = sensitivity
                self.__SNRthreshold = SNRthreshold
                self.__McMax        = McMax
                self.__McStep       = McStep
                self.__ETAmax       = ETAmax
                self.__ETAstep      = ETAstep
                self.__SNRmax       = SNRmax
                self.__SNRstep      = SNRstep

                # set initialised
                self.__initialised  = True

        if errStr is not None: self.__verbose = saveVerbose                                                         # error - restore verbose setting

        return errStr


    """
    SelectionEffects::__GenerateThetas()

    Generate an array of random realisations of inclination and sky position (thetas), used to calculate the average
    probability of detecting an event.  The size of the array (number of thetas) is specified by the class member
    variable 'numThetas'.

    Updates the following class member variables if no error occurs:
        __SNRthetas

    Args: None

    Returns:
        1 : 1D ARRAY of FLOAT : Thetas - random realisations of inclination and sky position.
    """
    def __GenerateThetas(self):
       
        # calculate thetas
        cosThetas = self.__rng.uniform(low = -1.0, high = 1.0, size = self.__numThetas)
        cosIncs   = self.__rng.uniform(low = -1.0, high = 1.0, size = self.__numThetas)
        phis      = self.__rng.uniform(low =  0.0, high = 2.0 * np.pi, size = self.__numThetas)
        zetas     = self.__rng.uniform(low =  0.0, high = 2.0 * np.pi, size = self.__numThetas)

        Fps = (0.5 * np.cos(2 * zetas) * (1.0 + cosThetas**2) * np.cos(2.0 * phis) - np.sin(2.0 * zetas) * cosThetas * np.sin(2.0 * phis))
        Fxs = (0.5 * np.sin(2 * zetas) * (1.0 + cosThetas**2) * np.cos(2.0 * phis) + np.cos(2.0 * zetas) * cosThetas * np.sin(2.0 * phis))

        self.__SNRthetas = np.sort(np.sqrt(0.25 * Fps**2 * (1.0 + cosIncs**2)**2 + Fxs**2 * cosIncs**2))


    """
    SelectionEffects::__ReadNoiseData()

    Reads SNR noise data from the dataset specified by the class member variable 'hdf5DatasetName' (corresponds to a
    specific sensitivity) in the noise data file specified by the class member variables 'NoiseFilepath' and 'NoiseFilename'

    Updates the following class member variables if no error occurs:
        __fqNoiseFilename
        __hdf5DatasetName

    Updates the following class member variables to parameter values if no error occurs:
        __noiseFilename
        __noiseFilepath
        __sensitivity

    Args:
        p_NoiseFilepath : STRING    : Path to SNR noise file
                                      Default = None, self.__NoiseFilepath will be used
        p_NoiseFilename : STRING    : Name of SNR noise file
                                      Default = None, self.__NoiseFilename will be used
        p_Sensitivity   : STRING    : Detector sensitivity
                                      Must be one of ['design', 'O1', 'O3', 'ET']
                                      Default = None, self.__sensitivity will be used

    Returns:
        1 : 1D ARRAY of FLOAT : Mass axis defining the (symmetric) grid.
                                'None' if an error occurred.
        2 : 2D ARRAY of FLOAT : Symmetric SNR grid as refernce points for interpolator - each dimension is the length of the mass axis.
                                'None' if an error occurred.
        3 : STRING            : Error string - set if an error occurred.
                                'None' if no error occurred.
    """
    def __ReadNoiseData(self,
                        p_NoiseFilepath = None, 
                        p_NoiseFilename = None, 
                        p_Sensitivity   = None): 

        func = self.__className + '::__ReadNoiseData(): '

        # default return values
        SNRmassAxis = None
        SNRgrid     = None
        errStr      = None

        # set local variables from parameter values
        # Doing it this way ensures that we have reasonable values for all parameters and allows users to call the function
        # with positional parameters and just set any parameters they want to default to 'None' rather than specify the value
        noiseFilepath = self.__NoiseFilepath if p_NoiseFilepath is None else p_NoiseFilepath
        noiseFilename = self.__NoiseFilename if p_NoiseFilename is None else p_NoiseFilename
        sensitivity   = self.__sensitivity   if p_Sensitivity   is None else p_Sensitivity

        # check parameter types and values
        tmpStr = ' must be a string'

        if   not isinstance(noiseFilepath, str):                              errStr = func + 'Path to noise file' + tmpStr
        elif not (isinstance(noiseFilename, str) and len(noiseFilename) > 0): errStr = func + 'Noise file name must be a non-empty string'
        elif not isinstance(sensitivity, str):                                errStr = func + 'Detector sensitivity' + tmpStr
        elif sensitivity not in SE_Sensitivities:                             errStr = func + 'Unknown detector sensitivity \'{:s}\': must be one of '.format(sensitivity) + SE_Sensitivities
        else:                                                                                                       # all parameters ok

            # open noise file
            if len(noiseFilepath) < 1: noiseFilepath = '.'
            fqNoiseFilename = noiseFilepath + '/' + noiseFilename
            if not os.path.isfile(fqNoiseFilename):                                                                 # opened ok?
                errStr = func + 'SNR noise file \'{:s}\' not found'.format(fqNoiseFilename)
            else:

                # process noise file
                with h5py.File(fqNoiseFilename, 'r') as SNRfile:

                    grouplist = list(SNRfile.keys())                                                                # get HDF5 file keys

                    # get mass axis if present
                    if 'mass_axis' not in grouplist:                                                                # mass axis present?
                        errStr = func + 'Group \'mass_axis\' not found in SNR file \'{:s}\''.format(fqNoiseFilename)
                    else:
                        SNRmassAxis = SNRfile['mass_axis'][...]                                                     # set mass axis

                        # get noise dataset names if present
                        if 'snr_values' not in grouplist:                                                           # noise dataset names present?
                            errStr = func + 'Group \'snr_values\' not found in SNR file \'{:s}\''.format(fqNoiseFilename)
                            SNRmassAxis = None                                                                      # error return
                        else:
                            SNRdatasets = SNRfile['snr_values']                                                     # set dataset names

                            # set SNR grid from noise file if dataset present
                            if   sensitivity == 'design': hdf5DatasetName = SE_Dataset_design
                            elif sensitivity == 'O1'    : hdf5DatasetName = SE_Dataset_O1
                            elif sensitivity == 'O3'    : hdf5DatasetName = SE_Dataset_O3
                            elif sensitivity == 'ET'    : hdf5DatasetName = SE_Dataset_ET
                        
                            if hdf5DatasetName not in SNRdatasets:                                                  # noise dataset present?
                                errStr = func + ' Group \'{:s}\' for sensitivity \'{:s}\' not found in SNR file \'{:s}\''.format(hdf5DatasetName, sensitivity, fqNoiseFilename)
                                SNRmassAxis = None                                                                  # error return
                            else:
                                SNRgrid = SNRdatasets[hdf5DatasetName][...]                                         # set grid

                                # update class member variables
                                self.__fqNoiseFilename = fqNoiseFilename
                                self.__hdf5DatasetName = hdf5DatasetName

                                # update class member variables to match parameters
                                self.__noiseFilepath = noiseFilepath
                                self.__noiseFilename = noiseFilename
                                self.__sensitivity   = sensitivity

        return SNRmassAxis, SNRgrid, errStr


    """
    SelectionEffects::__GetDetectionProbabilityFromSNR()

    Compute the probability of detecting a compact binary coalescence (CBC) with given SNR and threshold,
    averaging over all orientations and sky positions.

    Based on Finn & Chernoff 1993 (https://arxiv.org/abs/gr-qc/9301003).

    Args:
        p_SNRthreshold  : FLOAT  : SNR threshold required for detection
                                   Must be >= 0.0, default = None: self.__SNRthreshold will be used
        p_SNRmax        : FLOAT  : Maximum SNR in detection probability array
                                   Must be positive, default = None: self.__SNRmax will be used
        p_SNRstep       : FLOAT  : Step in SNR to use in detection probability array
                                   Must be positive, default = None: self.__SNRstep will be used

    Returns:
        1 : 1D ARRAY of FLOAT : Probability of the received signal to be above the threshold.
                                'None' if an error occurred.
        2 : STRING            : Error string - set if an error occurred.
                                'None' if no error occurred.
    """
    def __GetDetectionProbabilityFromSNR(self,
                                         p_SNRthreshold = None,
                                         p_SNRmax       = None,
                                         p_SNRstep      = None):

        func = self.__className + '::__GetDetectionProbabilityFromSNR(): '

        # default return values
        detectionProbability = None
        errStr               = None

        # set local variables from parameter values
        # Doing it this way ensures that we have reasonable values for all parameters and allows users to call the function
        # with positional parameters and just set any parameters they want to default to 'None' rather than specify the value
        SNRthreshold = self.__SNRthreshold if p_SNRthreshold is None else p_SNRthreshold
        SNRmax       = self.__SNRmax       if p_SNRmax       is None else p_SNRmax
        SNRstep      = self.__SNRstep      if p_SNRstep      is None else p_SNRstep

        # check parameter types and values
        tmpStr = ' must be a floating point number'                                                                 # common error string

        if   not (isinstance(SNRthreshold, float) and SNRthreshold >= 0.0): errStr = func + 'SNR threshold' + tmpStr + ' >= 0.0'
        elif not (isinstance(SNRmax, float) and SNRmax > 0.0):              errStr = func + 'Maximum SNR' + tmpStr + ' > 0.0'
        elif not (isinstance(SNRstep, float) and SNRstep > 0.0):            errStr = func + 'SNR step size' + tmpStr + ' > 0.0'
        else:

            # calculate detection probabilities
            thetaMin                              = SNRthreshold / np.arange(SNRstep, SNRmax + SNRstep, SNRstep)
            detectionProbability                  = np.zeros_like(thetaMin)
            idx                                   = np.searchsorted(self.__SNRthetas, thetaMin[thetaMin <= 1.0], side = 'right')
            detectionProbability[thetaMin <= 1.0] = 1.0 - ((idx - 1.0) / float(self.__SNRthetas.shape[0]))
 
        return detectionProbability, errStr


    """
    SelectionEffects::__ComputeSNRandDetectionGrids()

    Compute a grid of SNRs and detection probabilities for a range of masses and SNRs.  The grids allow for
    interpolating the values of the SNR and detection probability.

    Updates the following class member variables to parameter values if no error occurs:
        __ETAmax
        __ETAstep
        __McMax
        __McStep
        __SNRmax 
        __SNRstep
        __SNRthreshold

    Args:
        p_SNRthreshold  : FLOAT  : SNR threshold required for detection
                                   Must be >= 0.0, default = None: self.__SNRthreshold will be used
        p_McMax         : FLOAT  : Maximum chirp mass in SNR grid
                                   Must be positive, default = None: self.__McMax will be used
        p_McStep        : FLOAT  : Step in chirp mass to use in SNR grid
                                   Must be positive, default = None: self.__McStep will be used
        p_ETAmax        : FLOAT  : Maximum symmetric mass ratio in SNR grid
                                   Must be positive, default = None: self.__ETAmax will be used
        p_ETAstep       : FLOAT  : Step in symmetric mass ratio to use in SNR grid
                                   Must be positive, default = None: self.__ETAstep will be used
        p_SNRmax        : FLOAT  : Maximum SNR in detection probability array
                                   Must be positive, default = None: self.__SNRmax will be used
        p_SNRstep       : FLOAT  : Step in SNR to use in detection probability array
                                   Must be positive, default = None: self.__SNRstep will be used

    Returns:
        1 : 2D ARRAY of FLOAT : The SNR of a binary with masses (Mc, eta) at a distance of 1 Mpc.
                                'None' if an error occurred.
        2 : 1D ARRAY of FLOAT : Detection probabilities for different SNRs.
                                'None' if an error occurred.
        3 : STRING            : Error string - set if an error occurred.
                                'None' if no error occurred.
    """
    def __ComputeSNRandDetectionGrids(self,
                                      p_SNRinterpolator = None,
                                      p_SNRthreshold    = None, 
                                      p_McMax           = None, 
                                      p_McStep          = None,
                                      p_ETAmax          = None, 
                                      p_ETAstep         = None, 
                                      p_SNRmax          = None, 
                                      p_SNRstep         = None):

        func = self.__className + '::__ComputeSNRandDetectionGrids(): '

        # default return values
        SNRgridAt1Mpc               = None
        detectionProbabilityFromSNR = None
        errStr                      = None
        
        # set local variables from parameter values
        # Doing it this way ensures that we have reasonable values for all parameters and allows users to call the function
        # with positional parameters and just set any parameters they want to default to 'None' rather than specify the value
        SNRinterpolator = self.__SNRinterpolator if p_SNRinterpolator is None else p_SNRinterpolator
        SNRthreshold    = self.__SNRthreshold    if p_SNRthreshold    is None else p_SNRthreshold
        McMax           = self.__McMax           if p_McMax           is None else p_McMax
        McStep          = self.__McStep          if p_McStep          is None else p_McStep
        ETAmax          = self.__ETAmax          if p_ETAmax          is None else p_ETAmax
        ETAstep         = self.__ETAstep         if p_ETAstep         is None else p_ETAstep
        SNRmax          = self.__SNRmax          if p_SNRmax          is None else p_SNRmax
        SNRstep         = self.__SNRstep         if p_SNRstep         is None else p_SNRstep

        # check parameter types and values
        tmpStr = ' must be a floating point number'                                                                 # common error string
        
        if   type(SNRinterpolator) != SE_SNRinterpolatorType:               errStr = func + 'SNR interpolator must be of type \'' + SE_SNRinterpolatorType_asSTR + '\''
        elif not (isinstance(SNRthreshold, float) and SNRthreshold >= 0.0): errStr = func + 'SNR threshold' + tmpStr + ' >= 0.0'
        elif not (isinstance(McMax, float) and McMax > 0.0):                errStr = func + 'Maximum chirp mass' + tmpStr + ' > 0.0'
        elif not (isinstance(McStep, float) and McStep > 0.0):              errStr = func + 'Chirp mass step size' + tmpStr + ' > 0.0'
        elif not (isinstance(ETAmax, float) and ETAmax > 0.0):              errStr = func + 'Maximum symmetric mass ratio' + tmpStr + ' > 0.0'
        elif not (isinstance(ETAstep, float) and ETAstep > 0.0):            errStr = func + 'Symmetric mass ratio step size' + tmpStr + ' > 0.0'
        elif not (isinstance(SNRmax, float) and SNRmax > 0.0):              errStr = func + 'Maximum SNR' + tmpStr + ' > 0.0'
        elif not (isinstance(SNRstep, float) and SNRstep > 0.0):            errStr = func + 'SNR step size' + tmpStr + ' > 0.0'
        else:                                                                                                       # all parameters ok

            # create chirp mass and eta arrays
            Mc  = np.arange(McStep, McMax + McStep, McStep)
            eta = np.arange(ETAstep, ETAmax + ETAstep, ETAstep)

            # convert to total, primary and secondary mass arrays
            Mt = Mc / eta[:, np.newaxis]**0.6
            M1 = Mt * 0.5 * (1.0 + np.sqrt(1.0 - 4.0 * eta[:, np.newaxis]))
            M2 = Mt - M1

            # interpolate to get snr values if binary was at 1Mpc
            SNRgridAt1Mpc = SNRinterpolator(np.log(M1), np.log(M2), grid = False)

            # precompute detection probabilities as a function of snr
            detectionProbabilityFromSNR, errStr = self.__GetDetectionProbabilityFromSNR(SNRthreshold, SNRmax, SNRstep)
            if errStr is not None:
                SNRgridAt1Mpc = None                                                                                # error return

        return SNRgridAt1Mpc, detectionProbabilityFromSNR, errStr
    
    
    # public interface

    def Initialise(self,
                   p_Sensitivity   = None,
                   p_SNRthreshold  = None,
                   p_McMax         = None,
                   p_McStep        = None,
                   p_ETAmax        = None,
                   p_ETAstep       = None,
                   p_SNRmax        = None,
                   p_SNRstep       = None,
                   p_Verbose       = None): 
        return self.__Initialise(None, None, p_Sensitivity, p_SNRthreshold, p_McMax, p_McStep, p_ETAmax, p_ETAstep, p_SNRmax, p_SNRstep, p_Verbose)


    # getters

    """
    DetectionProbabilityFromSNR()

    Returns value of self.__detectionProbabilityFromSNR at index if a valid index is specified
    Returns complete self.__detectionProbabilityFromSNR array if no index specified
    Returns an error if an invalid index is specified
    """
    def DetectionProbabilityFromSNR(self, p_Index = None):

        func = self.__className + '::DetectionProbabilityFromSNR(): '

        # default return values
        prob   = None
        errStr = None

        if p_Index is None:
            prob = self.__detectionProbabilityFromSNR
        elif not (isinstance(p_Index, int) and p_Index >= 0 and p_Index < self.__detectionProbabilityFromSNR.shape[0]):
            errStr = func + 'Index value must be an integer >= 0 and < {:d}'.format(self.__detectionProbabilityFromSNR.shape[0])
        else:
            prob = self.__detectionProbabilityFromSNR[p_Index]

        return prob, errStr


    # the following getters take no parameters

    def DetectionProbabilityLen(self): return self.__detectionProbabilityFromSNR.shape[0]
    def ETAmax(self)                 : return self.__ETAmax
    def ETAstep(self)                : return self.__ETAstep
    def McMax(self)                  : return self.__McMax
    def McStep(self)                 : return self.__McStep
    def Sensitivity(self)            : return self.__sensitivity
    def SNRgridAt1Mpc(self)          : return self.__SNRgridAt1Mpc
    def SNRmax(self)                 : return self.__SNRmax
    def SNRstep(self)                : return self.__SNRstep
    def SNRthreshold(self)           : return self.__SNRthreshold



"""
Class COMPAS

Provides functionality to access and process COMPAS data from HDF5 file.

Public interface (public functions declared at end of class):

    COMPAS::COMPAS()     : Class constructor
    COMPAS::Initialise() : Reinitialises class with new parameters

    Getters:
        ChirpMasses()          : Returns chirpmass for a specified DCO, or entire array
        DCOtype()              : Returns the DCO type(s) requested
        DelayTime()            : Returns delay time for a specified DCO, or entire array
        ETAs()                 : Returns symmetric mass ratio for a specified DCO, or entire array
        MassEvolvedPerBinary() : Returns the star forming mass evolved per binary 
        MaxLogZ()              : Returns the natural log if the maximum metallicity in the COMPAS data
        MinLogZ()              : Returns the natural log if the minimum metallicity in the COMPAS data
        nAllDCOs()             : returns the number of all DCO types in the COMPAS data
        nBinaries()            : Returns the number of binaries in the COMPAS data
        nDCOs()                : Returns the number of DCOs in the COMPAS data that match the DCO type(s) requested
        PrimaryMass()          : Returns primary mass for a specified DCO, or entire array
        SecondaryMass()        : Returns secondary mass for a specified DCO, or entire array
        Zsystems()             : Returns metallicity for a specified DCO, or entire array

Adapted from ClassCOMPAS.py
"""
class COMPAS:       
    
    """
    COMPAS::COMPAS
    
    Class constructor

    Initialises class member variables

    Args:
        p_COMPASfilepath       : STRING    : Path to COMPAS HDF5 data file
                                             Default = COMPAS_filepath
        p_COMPASfilename       : STRING    : Name of COMPAS HDF5 data file
                                             Default = COMPAS_filename
        p_m1Minimum            : FLOAT     : COMPAS minimum value for primary star mass
                                             Must be >= 0.0, default = COMPAS_m1Minimum
        p_m1Maximum            : FLOAT     : COMPAS maximum value for primary star mass
                                             Must be >= p_m1Minimum, default = COMPAS_m1Maximum
        p_m2Minimum            : FLOAT     : COMPAS minimum value for secondary star mass
                                             Must be >= 0.0, default = COMPAS_m2Minimum
        p_UseSampledMassRanges : BOOL      : Flag to indicate that the mass ranges sampled by COMPAS should be used instead of
                                             hard cuts at p_m1Minimum, p_m1Maximum, and p_m2Minimum.
                                             Default = COMPAS_UseSampledMassRanges
        p_BinaryFraction       : FLOAT     : Fraction of stars that form binaries
                                             Must be >= 0.0 and <= 1.0, default = COMPAS_binaryFraction
        p_SFM_popSize          : INT       : Number of samples to draw when creating a mock universe
                                             Must be > 0, default = COMPAS_SFM_popSize
        p_DCOtype              : STRING    : DCO type for which masks should be set
                                             Must be of ['ALL', 'BBH', 'BHNS', 'BNS', 'CHE_BBH', 'NON_CHE_BBH'], default = COMPAS_DCO_ALL ('ALL')
        p_WithinHubbleTime     : BOOL      : Flag to only include DCOs that will merge within a Hubble time
                                             Default = COMPAS_flag_withinHubbleTime
        p_PessimisticCE        : BOOL      : Flag to set pessimistic/optimistic CE mask
                                             Default = COMPAS_flag_pessimisticCE
        p_NoRLOFafterCEE       : BOOL      : Flag to exclude DCOs that underwent RLOF after a CE event
                                             Default = COMPAS_flag_noRLOFafterCEE
        p_RNG                  : Numpy RNG : Numpy random number generator object (type = numpy.random._generator.Generator)
                                             Default = None
        p_Verbose              : BOOL      : Flag to indicate if diagnostics/stats should be printed to console
                                             Default = VERBOSE_DEFAULT

    Errors here are handled by raising an exception - the constructor can't return an error to the caller.
    """
    def __init__(self,
                 p_COMPASfilepath       = COMPAS_filepath, 
                 p_COMPASfilename       = COMPAS_filename, 
                 p_m1Minimum            = COMPAS_m1Minimum, 
                 p_m1Maximum            = COMPAS_m1Maximum, 
                 p_m2Minimum            = COMPAS_m2Minimum,
                 p_UseSampledMassRanges = COMPAS_UseSampledMassRanges,
                 p_BinaryFraction       = COMPAS_binaryFraction,
                 p_SFM_popSize          = COMPAS_SFM_popSize,
                 p_DCOtype              = COMPAS_DCO_ALL,
                 p_WithinHubbleTime     = COMPAS_flag_withinHubbleTime,
                 p_PessimisticCE        = COMPAS_flag_pessimisticCE,
                 p_NoRLOFafterCEE       = COMPAS_flag_noRLOFafterCEE,
                 p_RNG                  = None,
                 p_Verbose              = VERBOSE_DEFAULT):

        self.__className = 'COMPAS'
        func = self.__className + '::COMPAS(): '

        self.__verbose = VERBOSE_DEFAULT if p_Verbose is None else p_Verbose
        
        if self.__verbose:
            print('...Initialising COMPAS')
            mark = time.time()
       
        # set local variables from parameter values
        # Doing it this way ensures that we have reasonable values for all parameters and allows users to call the function
        # with positional parameters and just set any parameters they want to default to 'None' rather than specify the value
        COMPASfilepath       = COMPAS_filepath              if p_COMPASfilepath       is None else p_COMPASfilepath
        COMPASfilename       = COMPAS_filename              if p_COMPASfilename       is None else p_COMPASfilename
        m1Minimum            = COMPAS_m1Minimum             if p_m1Minimum            is None else p_m1Minimum
        m1Maximum            = COMPAS_m1Maximum             if p_m1Maximum            is None else p_m1Maximum
        m2Minimum            = COMPAS_m2Minimum             if p_m2Minimum            is None else p_m2Minimum
        useSampledMassRanges = COMPAS_useSampledMassRanges  if p_UseSampledMassRanges is None else p_UseSampledMassRanges
        binaryFraction       = COMPAS_binaryFraction        if p_BinaryFraction       is None else p_BinaryFraction
        SFM_popSize          = COMPAS_SFM_popSize           if p_SFM_popSize          is None else p_SFM_popSize
        DCOtype              = COMPAS_DCO_BBH               if p_DCOtype              is None else p_DCOtype
        withinHubbleTime     = COMPAS_flag_withinHubbleTime if p_WithinHubbleTime     is None else p_WithinHubbleTime
        pessimisticCE        = COMPAS_flag_pessimisticCE    if p_PessimisticCE        is None else p_PessimisticCE
        noRLOFafterCEE       = COMPAS_flag_noRLOFafterCEE   if p_NoRLOFafterCEE       is None else p_NoRLOFafterCEE
        rng                  = p_RNG

        # check parameter types and values
        if not isinstance(COMPASfilepath, str):                                                                     # all we can do for now - check existence in ReadData()
            raise RuntimeError(func + 'Path to COMPAS data file must be a string')
        if len(COMPASfilepath) < 1: COMPASfilepath = '.'
        
        if not (isinstance(COMPASfilename, str) and len(COMPASfilename) > 0):                                       # all we can do for now - check existence in ReadData()
            raise RuntimeError(func + 'COMPAS data file name must be a non-empty string')
        
        tmpStr1 = ' must be a floating point number >= 0.0'                                                         # common error string               
        tmpStr2 = ' flag must be True or False'                                                                     # common error string      

        if not (isinstance(m1Minimum, float) and m1Minimum >= 0.0):       raise RuntimeError(func + 'Minimum primary star mass' + tmpStr1)
        if not (isinstance(m1Maximum, float) and m1Maximum >= m1Minimum): raise RuntimeError(func + 'Maximum primary star mass must >= minimum primary star mass ({:.6f})'.format(m1Minimum))
        if not (isinstance(m2Minimum, float) and m2Minimum >= 0.0):       raise RuntimeError(func + 'Minimum secondary star mass' + tmpStr1)
        if not isinstance(useSampledMassRanges, bool):                    raise RuntimeError(func + '\'Use sampled mass ranges\'' + tmpStr2)
        if not (isinstance(binaryFraction, float) and \
               binaryFraction >= 0.0 and binaryFraction <= 1.0):          raise RuntimeError(func + 'Binary fraction must' + tmpStr1 + ' and <= 1.0')
        if not (isinstance(SFM_popSize, int) and SFM_popSize > 0):        raise RuntimeError(func + 'SFM population size must be a positive  integer')
        if not (isinstance(DCOtype, str) and DCOtype in DCO_Types):       raise RuntimeError(func + 'DCO type must be one of' + DCO_Types)
        if not isinstance(withinHubbleTime, bool):                        raise RuntimeError(func + '\'Merges within a Hubble time\'' + tmpStr2)
        if not isinstance(pessimisticCE, bool):                           raise RuntimeError(func + '\'Pessimistic CE\'' + tmpStr2)
        if not isinstance(noRLOFafterCEE, bool):                          raise RuntimeError(func + '\'Immediate RLOF after CE\'' + tmpStr2)
        if rng is not None and type(rng) != RNGType:                      raise RuntimeError(func + 'Random number generator must be of type \'' + RNGType_asSTR + '\'')

        # construct a random number generator if necessary - no initial seed
        self.__rng = np.random.default_rng() if rng is None else rng

        # set class member variables
        self.__COMPASfilepath       = COMPASfilepath
        self.__COMPASfilename       = COMPASfilename
        self.__fqCOMPASfilename     = COMPASfilepath + '/' + COMPASfilename
        self.__m1Minimum            = m1Minimum
        self.__m1Maximum            = m1Maximum
        self.__m2Minimum            = m2Minimum
        self.__useSampledMassRanges = useSampledMassRanges
        self.__binaryFraction       = binaryFraction
        self.__SFM_popSize          = SFM_popSize

        # these will be set by self__.ReadData() - directly or indirectly
        self.__version              = None
        self.__runStart             = None
        self.__sysSeeds             = None
        self.__nBinaries            = None
        self.__zamsST1              = None
        self.__zamsST2              = None
        self.__zamsMass1            = None
        self.__zamsMass2            = None
        self.__zamsZ                = None
        self.__CHEonMS1             = None
        self.__CHEonMS2             = None
        self.__allDCOseeds          = None
        self.__st1                  = None
        self.__st2                  = None
        self.__allMass1             = None
        self.__allMmass2            = None
        self.__formationTime        = None
        self.__coalescenceTime      = None
        self.__mergesInHubbleTime   = None
        self.__ceSeeds              = None
        self.__immediateRLOF        = None
        self.__optimisticCE         = None


        # these will be set by self.__Initialise() - directly or indirectly
        self.__DCOtype              = None
        self.__withinHubbleTimeFlag = None
        self.__pessimisticCEFlag    = None
        self.__noRLOFafterCEEFlag   = None
        self.__massEvolvedPerBinary = None
        self.__maskedDCOseeds       = None
        self.__nAllDCOs             = None
        self.__nDCOs                = None
        self.__mass1                = None
        self.__mass2                = None
        self.__formationTime        = None
        self.__coalescenceTime      = None
        self.__mergesInHubbleTime   = None
        self.__Zsystems             = None
        self.__delayTime            = None
        self.__chirpMasses          = None
        self.__ETAs                 = None

        self.__DCOmask              = None
        self.__BBHmask              = None
        self.__BHNSmask             = None
        self.__DNSmask              = None
        self.__CHE_BBHmask          = None
        self.__NON_CHE_BBHmask      = None
        self.__ALL_TYPESmask        = None
        self.__OPTIMISTICmask       = None

        
        # read data file - fail if an error occurs
        errStr = self.__ReadData()
        if errStr is not None: raise RuntimeError(errStr)

        # reset primary max and min, and secondary min, values if caller wants
        # to use mass ranges sampled by COMPAS
        if self.__useSampledMassRanges:                                                                             # use mass ranges sampled by COMPAS instead of hard cuts?
            self.__m1Minimum = np.min(self.__mass1[self.__mass1 != self.__mass2])                                   # yes - don't include masses set equal through RLOF at ZAMS
            self.__m1Maximum = np.max(self.__mass1)
            self.__m2Minimum = np.min(self.__mass2)

        # set min and max logZ values
        self.__minLogZ = np.log(np.min(self.__zamsZ))
        self.__maxLogZ = np.log(np.max(self.__zamsZ))

        # set initial DCO masks - can be changed later
        # calculate population values based on the DCOmask requested
        errStr = self.__Initialise(self.__m1Minimum, self.__m1Maximum, self.__m2Minimum, useSampledMassRanges, binaryFraction, DCOtype, withinHubbleTime, pessimisticCE, noRLOFafterCEE)

        if errStr is not None: raise RuntimeError(errStr)

        if self.__verbose: 
            timeTaken = time.time() - mark
            memUsage  = psutil.Process().memory_info().rss / GB
            print('   COMPAS initialised in {:.2f} seconds.  Process total memory usage = {:.6f} GB'.format(timeTaken, memUsage))


    """
    COMPAS::__ReadRunDetailsData()

    Reads required Run Details data from COMPAS HDF5 file.
    Checks data read as much as possible here.

    Updates the following class member variables if no error occurs:
        __Version
        __RunDate

    Args:
        p_COMPASfile : HDF5 FILE : Open COMPAS HDF5 data file
                                   Default = None

    Returns:
        1 : STRING : Error string - set if an error occurred.
                     'None' if no error occurred.
    """
    def __ReadRunDetailsData(self, p_COMPASfile = None):

        func = self.__className + '::__ReadRunDetailsData(): '

        if self.__verbose: 
            print('         ...Reading Run Details data')
            mark = time.time()

        # default return values
        errStr = None
        
        COMPASfile = p_COMPASfile                                                                                   # COMPAS data file

        # check parameter types and values
        if COMPASfile is None or type(COMPASfile) != h5py._hl.files.File:                                           # valid COMPAS file specified?
            errStr = func + 'No valid COMPAS data file specified'
        else:
            groupList = list(COMPASfile.keys())                                                                     # list of groups in COMPAS data file
            if len(groupList) < 1:                                                                                  # empty list?
                errStr = func + 'COMPAS data file contains no groups'
            else:
                # check that run details group exists
                if COMPAS_RUN_DETAILS not in groupList:                                                             # run details group in data file?
                    errStr = func + 'Group \'{:s}\' not found in COMPAS data file \'{:s}\''.format(COMPAS_RUN_DETAILS, self.__fqCOMPASfilename)
                else:
                    # extract run details data from data file
                    runGroup    = COMPASfile[COMPAS_RUN_DETAILS]                                                    # run details data
                    datasetList = list(runGroup.keys())                                                             # list of keys (datasets)

                    # check that all required datasets exist in group
                    tmpStr = func + 'Dataset \'{:s}\' not found in \'{:s}\' group in COMPAS file \'{:s}\''          # common error string prefix

                    if   COMPAS_version  not in datasetList: errStr = tmpStr.format(COMPAS_version,  COMPAS_RUN_DETAILS, self.__fqCOMPASfilename)
                    elif COMPAS_runStart not in datasetList: errStr = tmpStr.format(COMPAS_runStart, COMPAS_RUN_DETAILS, self.__fqCOMPASfilename)
                    else:                                                                                           # all required datasets exist - extract data
                        self.__version  = runGroup[COMPAS_version][0].decode('utf-8')                               # COMPAS version
                        self.__runStart = runGroup[COMPAS_runStart][0].decode('utf-8')                              # COMPAS run start date
                        if self.__runStart[-1] == '\n': self.__runStart = self.__runStart[:-1]                      # remove newline if present

        if self.__verbose: 
            timeTaken = time.time() - mark
            memUsage  = psutil.Process().memory_info().rss / GB
            print('            Run details data read in {:.2f} seconds.  Process total memory usage = {:.6f} GB'.format(timeTaken, memUsage))

        return errStr


    """
    COMPAS::__ReadSYSPARMSData()

    Reads required System Parameters data from COMPAS HDF5 file.
    Checks data read as much as possible here.

    Updates the following class member variables if no error occurs:
        __CHEonMS1
        __CHEonMS2
        __nBinaries
        __sysSeeds
        __zamsMass1
        __zamsMass2
        __zamsST1
        __zamsST2
        __zamsZ

    Args:
        p_COMPASfile : HDF5 FILE : Open COMPAS HDF5 data file
                                   Default = None

    Returns:
        1 : STRING : Error string - set if an error occurred.
                     'None' if no error occurred.
    """
    def __ReadSYSPARMSData(self, p_COMPASfile = None):

        func = self.__className + '::__ReadSYSPARMSData(): '

        if self.__verbose: 
            print('         ...Reading SYSPARMS data')
            mark = time.time()

        # default return values
        errStr = None
        
        COMPASfile = p_COMPASfile                                                                                   # COMPAS data file

        # check parameter types and values
        if COMPASfile is None or type(COMPASfile) != h5py._hl.files.File:                                           # valid COMPAS file specified?
            errStr = func + 'No valid COMPAS data file specified'
        else:
            groupList = list(COMPASfile.keys())                                                                     # list of groups in COMPAS data file
            if len(groupList) < 1:                                                                                  # empty list?
                errStr = func + 'COMPAS data file contains no groups'
            else:
                # check that sysparms group exists
                if COMPAS_BSE_SYSPARMS not in groupList:                                                            # sys parms group in data file?
                    errStr = func + 'Group \'{:s}\' not found in COMPAS data file \'{:s}\''.format(COMPAS_BSE_SYSPARMS, self.__fqCOMPASfilename)
                else:
                    # extract sys parms data from data file
                    sysGroup    = COMPASfile[COMPAS_BSE_SYSPARMS]                                                   # sys parms data
                    datasetList = list(sysGroup.keys())                                                             # list of keys (datasets)

                    # check that all required datasets exist in group
                    tmpStr = func + 'Dataset \'{:s}\' not found in \'{:s}\' group in COMPAS file \'{:s}\''          # common error string prefix

                    if   COMPAS_seed              not in datasetList: errStr = tmpStr.format(COMPAS_seed, COMPAS_BSE_SYSPARMS, self.__fqCOMPASfilename)
                    elif COMPAS_BSE_SYS_zamsST1   not in datasetList: errStr = tmpStr.format(COMPAS_BSE_SYS_zamsST1, COMPAS_BSE_SYSPARMS, self.__fqCOMPASfilename)
                    elif COMPAS_BSE_SYS_zamsST2   not in datasetList: errStr = tmpStr.format(COMPAS_BSE_SYS_zamsST2, COMPAS_BSE_SYSPARMS, self.__fqCOMPASfilename)
                    elif COMPAS_BSE_SYS_zamsMass1 not in datasetList: errStr = tmpStr.format(COMPAS_BSE_SYS_zamsMass1, COMPAS_BSE_SYSPARMS, self.__fqCOMPASfilename)
                    elif COMPAS_BSE_SYS_zamsMass2 not in datasetList: errStr = tmpStr.format(COMPAS_BSE_SYS_zamsMass2, COMPAS_BSE_SYSPARMS, self.__fqCOMPASfilename)
                    elif COMPAS_BSE_SYS_zamsZ     not in datasetList: errStr = tmpStr.format(COMPAS_BSE_SYS_zamsZ, COMPAS_BSE_SYSPARMS, self.__fqCOMPASfilename)
                    else:                                                                                           # all required datasets exist - extract data
                        sysSeeds = sysGroup[COMPAS_seed][...]                                                       # SEEDS
                        nBinaries = len(sysSeeds)                                                                   # number of binaries (= #SEEDS)
                        if nBinaries < 1:                                                                           # at least 1 system?
                            errStr = func + 'No systems found in \'{:s}\' group in COMPAS data file \'{:s}\''.format(COMPAS_BSE_SYSPARMS, self.__fqCOMPASfilename)
                        else:

                            self.__sysSeeds = sysSeeds
                            self.__nBinaries = nBinaries

                            # we assume that if len(SEEDS) > 0 there are commensurate values in the other datasets
                            # (not necessarily true, but will fail (hopefully gracefully) later if not true)
                            self.__zamsST1   = sysGroup[COMPAS_BSE_SYS_zamsST1][...]                                # primary star stellar type at ZAMS
                            self.__zamsST2   = sysGroup[COMPAS_BSE_SYS_zamsST2][...]                                # secondary star stellar type at ZAMS
                            self.__zamsMass1 = sysGroup[COMPAS_BSE_SYS_zamsMass1][...]                              # primary star mass at ZAMS
                            self.__zamsMass2 = sysGroup[COMPAS_BSE_SYS_zamsMass2][...]                              # secondary star mass at ZAMS
                            self.__zamsZ     = sysGroup[COMPAS_BSE_SYS_zamsZ][...]                                  # primary star metallicity at ZAMS

                            # CHE info may not be available - just disable CHE DCO types if not
                            self.__CHEonMS1 = None                                                                  # default is disabled
                            self.__CHEonMS2 = None                                                                  # default is disabled
                            if COMPAS_BSE_SYS_CHE1 in datasetList:                                                  # have CHE info for primary star?
                                self.__CHEonMS1 = sysGroup[COMPAS_BSE_SYS_CHE1][...].astype(bool)                   # yes - extract the values
            
                            if COMPAS_BSE_SYS_CHE2 in datasetList:                                                  # have CHE info for secondary star?    
                                self.__CHEonMS2 = sysGroup[COMPAS_BSE_SYS_CHE2][...].astype(bool)                   # yes - extract the values

        if self.__verbose: 
            timeTaken = time.time() - mark
            memUsage  = psutil.Process().memory_info().rss / GB
            print('            SYSPARMS data read in {:.2f} seconds.  Process total memory usage = {:.6f} GB'.format(timeTaken, memUsage))

        return errStr


    """
    COMPAS::__ReadDCOData()

    Reads required Double Compact Objects data from COMPAS HDF5 file.
    Checks data read as much as possible here.

    Updates the following class member variables if no error occurs:
        __coalescenceTime
        __dcoSeeds
        __formationTime
        __mass1
        __mass2
        __mergesInHubbleTime
        __st1
        __st2

    Args:
        p_COMPASfile : HDF5 FILE : Open COMPAS HDF5 data file
                                   Default = None

    Returns:
        1 : STRING : Error string - set if an error occurred.
                     'None' if no error occurred.
    """
    def __ReadDCOData(self, p_COMPASfile = None):

        func = self.__className + '::__ReadDCOData(): '

        if self.__verbose: 
            print('         ...Reading DCO data')
            mark = time.time()

        # default return values
        errStr = None
        
        COMPASfile = p_COMPASfile                                                                                   # COMPAS data file

        # check parameter types and values
        if COMPASfile is None or type(COMPASfile) != h5py._hl.files.File:                                           # valid COMPAS file specified?
            errStr = func + 'No valid COMPAS data file specified'
        else:                                                                                                       # check group list
            groupList = list(COMPASfile.keys())                                                                     # list of groups in COMPAS data file
            if len(groupList) < 1:                                                                                  # empty list?
                errStr = func + 'COMPAS data file contains no groups'
            else:                                                                                                   # non-empty group list
                # check that DCO group exists
                if COMPAS_BSE_DCO not in groupList:                                                                 # DCO group in data file?
                    errStr = func + 'Group \'{:s}\' not found in COMPAS data file \'{:s}\''.format(COMPAS_BSE_DCO, self.__fqCOMPASfilename)
                else:
                    # extract sys parms data from data file
                    dcoGroup    = COMPASfile[COMPAS_BSE_DCO]                                                        # DCO data
                    datasetList = list(dcoGroup.keys())                                                             # list of keys (datasets)

                    # check that all required datasets exist in group
                    tmpStr = func + 'Dataset \'{:s}\' not found in \'{:s}\' group in COMPAS file \'{:s}\''          # common error string prefix

                    if   COMPAS_seed                    not in datasetList: errStr = tmpStr.format(COMPAS_seed, COMPAS_BSE_DCO, self.__fqCOMPASfilename)
                    elif COMPAS_BSE_DCO_ST1             not in datasetList: errStr = tmpStr.format(COMPAS_BSE_DCO_ST1, COMPAS_BSE_DCO, self.__fqCOMPASfilename)
                    elif COMPAS_BSE_DCO_ST2             not in datasetList: errStr = tmpStr.format(COMPAS_BSE_DCO_ST2, COMPAS_BSE_DCO, self.__fqCOMPASfilename)
                    elif COMPAS_BSE_DCO_Mass1           not in datasetList: errStr = tmpStr.format(COMPAS_BSE_DCO_Mass1, COMPAS_BSE_DCO, self.__fqCOMPASfilename)
                    elif COMPAS_BSE_DCO_Mass2           not in datasetList: errStr = tmpStr.format(COMPAS_BSE_DCO_Mass2, COMPAS_BSE_DCO, self.__fqCOMPASfilename)
                    elif COMPAS_BSE_DCO_formationTime   not in datasetList: errStr = tmpStr.format(COMPAS_BSE_DCO_formationTime, COMPAS_BSE_DCO, self.__fqCOMPASfilename)
                    elif COMPAS_BSE_DCO_coalescenceTime not in datasetList: errStr = tmpStr.format(COMPAS_BSE_DCO_coalescence, COMPAS_BSE_DCO, self.__fqCOMPASfilename)
                    elif COMPAS_BSE_DCO_mergesHubble    not in datasetList: errStr = tmpStr.format(COMPAS_BSE_DCO_mergesHubble, COMPAS_BSE_DCO, self.__fqCOMPASfilename)
                    else:                                                                                           # all required datasets exist - extract data
                        dcoSeeds = dcoGroup[COMPAS_seed][...]                                                       # SEEDS
                        if len(dcoSeeds) < 1:                                                                       # at least 1 system?
                            errStr = func + 'No DCOs found in \'{:s}\' group in COMPAS data file \'{:s}\''.format(COMPAS_BSE_DCO, self.__fqCOMPASfilename)
                        else:

                            self.__allDCOseeds = dcoSeeds

                            # we assume that if len(SEEDS) > 1 there are commensurate values in the other datasets
                            # (not necessarily true, but will fail (hopefully gracefully) later if not true)
                            self.__st1                = dcoGroup[COMPAS_BSE_DCO_ST1][...]                           # primary star stellar type
                            self.__st2                = dcoGroup[COMPAS_BSE_DCO_ST2][...]                           # secondary star stellar type
                            self.__allMass1           = dcoGroup[COMPAS_BSE_DCO_Mass1][...]                         # primary star mass
                            self.__allMass2           = dcoGroup[COMPAS_BSE_DCO_Mass2][...]                         # secondary star mass
                            self.__formationTime      = dcoGroup[COMPAS_BSE_DCO_formationTime][...]                 # DCO formation time
                            self.__coalescenceTime    = dcoGroup[COMPAS_BSE_DCO_coalescenceTime][...]               # DCO coalescence time
                            self.__mergesInHubbleTime = dcoGroup[COMPAS_BSE_DCO_mergesHubble][...].astype(bool)     # DCO merges in a Hubble time flag

        if self.__verbose: 
            timeTaken = time.time() - mark
            memUsage  = psutil.Process().memory_info().rss / GB
            print('            DCO data read in {:.2f} seconds.  Process total memory usage = {:.6f} GB'.format(timeTaken, memUsage))

        return errStr


    """
    COMPAS::__ReadCEData()

    Reads required Common Envelopes data from COMPAS HDF5 file.
    Checks data read as much as possible here.

    Updates the following class member variables if no error occurs:
        __ceSeeds
        __immediateRLOF
        __optimisticCE

    Args:
        p_COMPASfile : HDF5 FILE : Open COMPAS HDF5 data file
                                   Default = None

    Returns:
        1 : STRING : Error string - set if an error occurred.
                     'None' if no error occurred.
    """
    def __ReadCEData(self, p_COMPASfile = None):

        func = self.__className + '::__ReadCEData(): '

        if self.__verbose: 
            print('         ...Reading CE data')
            mark = time.time()

        # default return values
        errStr = None
        
        COMPASfile = p_COMPASfile                                                                                   # COMPAS data file

        # check parameter types and values
        if COMPASfile is None or type(COMPASfile) != h5py._hl.files.File:                                           # valid COMPAS file specified?
            errStr = func + 'No valid COMPAS data file specified'
        else:                                                                                                       # check group list
            groupList = list(COMPASfile.keys())                                                                     # list of groups in COMPAS data file
            if len(groupList) < 1:                                                                                  # empty list?
                errStr = func + 'COMPAS data file contains no groups'
            else:
                # check that CE group exists
                if COMPAS_BSE_CE not in groupList:                                                                  # CE group in data file?
                    errStr = func + 'Group \'{:s}\' not found in COMPAS data file \'{:s}\''.format(COMPAS_BSE_CE, self.__fqCOMPASfilename)
                else:
                    # extract sys parms data from data file
                    ceGroup    = COMPASfile[COMPAS_BSE_CE]                                                          # CE data
                    datasetList = list(ceGroup.keys())                                                              # list of keys (datasets)

                    # check that all required datasets exist in group
                    tmpStr = func + 'Dataset \'{:s}\' not found in \'{:s}\' group in COMPAS file \'{:s}\''          # common error string prefix

                    if   COMPAS_seed                 not in datasetList: errStr = tmpStr.format(COMPAS_seed, COMPAS_BSE_CE, self.__fqCOMPASfilename)
                    elif COMPAS_BSE_CE_immediateRLOF not in datasetList: errStr = tmpStr.format(COMPAS_BSE_CE_immediateRLOF, COMPAS_BSE_CE, self.__fqCOMPASfilename)
                    elif COMPAS_BSE_CE_optimisticCE  not in datasetList: errStr = tmpStr.format(COMPAS_BSE_CE_optimisticCE, COMPAS_BSE_CE, self.__fqCOMPASfilename)
                    else:                                                                                           # all required datasets exist - extract data
                        ceSeeds = ceGroup[COMPAS_seed][...]                                                         # SEEDS
                        if len(ceSeeds) < 1:                                                                        # at least 1 event?
                            errStr = func + 'No CE evnts found in \'{:s}\' group in COMPAS data file \'{:s}\''.format(COMPAS_BSE_CE, self.__fqCOMPASfilename)
                        else:

                            self.__ceSeeds = ceSeeds

                            # we assume that if len(SEEDS) > 1 there are commensurate values in the other datasets
                            # (not necessarily true, but will fail (hopefully gracefully) later if not true)
                            self.__immediateRLOF = ceGroup[COMPAS_BSE_CE_immediateRLOF][...].astype(bool)           # immediate RLOF>CE flag
                            self.__optimisticCE  = ceGroup[COMPAS_BSE_CE_optimisticCE][...].astype(bool)            # optimistic CE flag

        if self.__verbose: 
            timeTaken = time.time() - mark
            memUsage  = psutil.Process().memory_info().rss / GB
            print('            CE data read in {:.2f} seconds.  Process total memory usage = {:.6f} GB'.format(timeTaken, memUsage))

        return errStr


    """
    COMPAS::__ReadData()

    Reads required COMPAS data from HDF5 file.

    Returns:
        1 : STRING : Error string - set if an error occurred.
                     'None' if no error occurred.
    """
    def __ReadData(self):

        func = self.__className + '::__ReadData(): '

        if self.__verbose: 
            print('      ...Reading data from HDF5 file')
            mark = time.time()

        # default return values
        errStr = None

        # open COMPAS HDF5 file
        if not os.path.isfile(self.__fqCOMPASfilename):
            errStr = func + 'COMPAS HDF5 file \'{:s}\' not found'.format(self.__fqCOMPASfilename)
        else:        
            # extract data
            with h5py.File(self.__fqCOMPASfilename, 'r') as COMPASfile:

                errStr = self.__ReadRunDetailsData(COMPASfile)                                                      # get run details info
                if errStr is None:                                                                                  # ok to proceed?
                    errStr = self.__ReadSYSPARMSData(COMPASfile)                                                    # yes - get system parameters info
                    if errStr is None:                                                                              # ok to proceed?
                        errStr = self.__ReadDCOData(COMPASfile)                                                     # yes - get double compact objects info
                        if errStr is None:                                                                          # ok to proceed?
                            errStr = self.__ReadCEData(COMPASfile)                                                  # yes - get common envelopes info

        if self.__verbose: 
            timeTaken = time.time() - mark
            memUsage  = psutil.Process().memory_info().rss / GB
            print('         Data read from HDF5 file in {:.2f} seconds.  Process total memory usage = {:.6f} GB'.format(timeTaken, memUsage))

        return errStr


    """
    COMPAS::__Initialise()

    Sets DCO masks appropriately for use throughout.
    Calculates population values based on DCO masks.

    Note: does not re-read data from HDF5 file.

    Updates the following class member variables if no error occurs:
        __ALL_TYPESmask
        __BBHmask
        __BHNSmask
        __CHE_BBHmask
        __DCOmask
        __DNSmask
        __NON_CHE_BBHmask
        __OPTIMISTICmask

    Updates the following class member variables to parameter values if no error occurs:
        __DCOtype
        __noRLOFafterCEEFlag
        __pessimisticCEFlag
        __useSampledMassRanges
        __withinHubbleTimeFlag

    Args:
        p_m1Minimum            : FLOAT  : COMPAS minimum value for primary star mass
                                          Must be >= 0.0, default = COMPAS_m1Minimum
        p_m1Maximum            : FLOAT  : COMPAS maximum value for primary star mass
                                          Must be >= p_m1Minimum, default = COMPAS_m1Maximum
        p_m2Minimum            : FLOAT  : COMPAS minimum value for secondary star mass
                                          Must be >= 0.0, default = COMPAS_m2Minimum
        p_UseSampledMassRanges : BOOL   : Flag to indicate that the mass ranges sampled by COMPAS should be used instead of
                                          hard cuts at p_m1Minimum, p_m1Maximum, and p_m2Minimum.
                                          Default = COMPAS_UseSampledMassRanges
        p_BinaryFraction       : FLOAT  : Fraction of stars that form binaries
                                          Must be >= 0.0 and <= 1.0, default = COMPAS_binaryFraction
        p_DCOtype              : STRING : DCO type for which masks should be set
                                          Must be one of ['ALL', 'BBH', 'BHNS', 'BNS', 'CHE_BBH', 'NON_CHE_BBH']
                                          Default = None: self.__DCOtype will be used
        p_WithinHubbleTime     : BOOL   : Flag to only include DCOs that will merge within a Hubble time.
                                          Default = None: self.__withinHubbleTimeFlag will be used
        p_PessimisticCE        : BOOL   : Flag to set pessimistic/optimistic CE mask
                                          Default = None: self.__pessimisticCEFlag will be used
        p_NoRLOFafterCEE       : BOOL   : Flag to exclude DCOs that underwent RLOF after a CE event
                                          Default = None: self.__noRLOFafterCEEFlag will be used
        p_Verbose              : BOOL   : Flag to indicate if diagnostics/stats should be printed to console
                                          Default = VERBOSE_DEFAULT

    Returns:
        1 : STRING : Error string - set if an error occurred.
                     'None' if no error occurred.
    """
    def __Initialise(self,
                     p_m1Minimum            = None,
                     p_m1Maximum            = None,
                     p_m2Minimum            = None,
                     p_UseSampledMassRanges = None,
                     p_BinaryFraction       = None,
                     p_DCOtype              = None,
                     p_WithinHubbleTime     = None,
                     p_PessimisticCE        = None,
                     p_NoRLOFafterCEE       = None,
                     p_Verbose              = None):

        func = self.__className + '::__Initialise(): '

        saveVerbose = self.__verbose                                                                                            # in case of error
        self.__verbose = self.__verbose if p_Verbose is None else p_Verbose                                                     # rset verbose setting

        # default return value
        errStr = None

        # set local variables from parameter values
        # Doing it this way ensures that we have reasonable values for all parameters and allows users to call the function
        # with positional parameters and just set any parameters they want to default to 'None' rather than specify the value
        DCOtype              = self.__DCOtype              if p_DCOtype              is None else p_DCOtype
        withinHubbleTime     = self.__withinHubbleTimeFlag if p_WithinHubbleTime     is None else p_WithinHubbleTime
        pessimisticCE        = self.__pessimisticCEFlag    if p_PessimisticCE        is None else p_PessimisticCE
        noRLOFafterCEE       = self.__noRLOFafterCEEFlag   if p_NoRLOFafterCEE       is None else p_NoRLOFafterCEE
        useSampledMassRanges = self.__useSampledMassRanges if p_UseSampledMassRanges is None else p_UseSampledMassRanges

        # check parameter types and values
        tmpStr = ' flag must be True or False'                                                                              # common error string  

        if   not (isinstance(DCOtype, str) and DCOtype in DCO_Types): errStr = func + 'DCO type must be one of' + DCO_Types
        elif not isinstance(withinHubbleTime, bool):                  errStr = func + '\'Merges within a Hubble time\'' + tmpStr
        elif not isinstance(pessimisticCE, bool):                     errStr = func + '\'Pessimistic CE\'' + tmpStr
        elif not isinstance(noRLOFafterCEE, bool):                    errStr = func + '\'Immediate RLOF after CE\'' + tmpStr
        else:                                                                                                               # all parameters ok

            DCOmask = self.__DCOmask

            # set all DCO type masks

            if self.__verbose:
                print('      ...Setting DCO masks')
                mark = time.time()

            typeMasks = {
                COMPAS_DCO_ALL : np.repeat(True, len(self.__allDCOseeds)),
                COMPAS_DCO_BBH : np.logical_and(self.__st1 == COMPAS_BH, self.__st2 == COMPAS_BH),
                COMPAS_DCO_BHNS: np.logical_or(np.logical_and(self.__st1  == COMPAS_BH, self.__st2 == COMPAS_NS), np.logical_and(self.__st1  == COMPAS_NS, self.__st2 == COMPAS_BH)),
                COMPAS_DCO_BNS : np.logical_and(self.__st1  == COMPAS_NS, self.__st2 == COMPAS_NS),
            }
            typeMasks[COMPAS_DCO_CHE_BBH]     = np.repeat(False, len(self.__allDCOseeds))                               # for now - updated below
            typeMasks[COMPAS_DCO_NON_CHE_BBH] = np.repeat(True, len(self.__allDCOseeds))                                # for now - updated below

            # set CHE type masks if required and if able
            if DCOtype == COMPAS_DCO_CHE_BBH or DCOtype == COMPAS_DCO_NON_CHE_BBH:                                      # if required
                if self.__CHonMS1 is not None and self.__CHonMS2 is not None:                                           # ... and if able
                    mask     = np.logical_and.reduce((self.__zamsST1 == COMPAS_CH, self.__zamsST2 == COMPAS_CH, self.__CHonMS1 == True, self.__CHonMS2 == True))
                    cheSeeds = self.__sysSeeds[...][mask]
                    mask     = np.in1d(self.__allDCOseeds, cheSeeds)

                    if DCOtype == COMPAS_DCO_CHE_BBH    : typeMasks[COMPAS_DCO_CHE_BBH]     = np.logical_and(mask, type_masks[COMPAS_DCO_BBH])
                    if DCOtype == COMPAS_DCO_NON_CHE_BBH: typeMasks[COMPAS_DCO_NON_CHE_BBH] = np.logical_and(np.logical_not(mask), type_masks[COMPAS_DCO_BBH])

            hubbleMask = self.__mergesInHubbleTime if withinHubbleTime else np.repeat(True, len(self.__allDCOseeds))

            rlofMask = np.repeat(True, len(self.__allDCOseeds))
            pessimisticMask = np.repeat(True, len(self.__allDCOseeds))
            if noRLOFafterCEE or pessimisticCE:                                                                         # if required
                if self.__ceSeeds is not None and self.__immediateRLOF is not None and self.__optimisticCE is not None: # ... and if able

                    dcoFromCE  = np.in1d(self.__ceSeeds, self.__allDCOseeds)
                    dcoCEseeds = self.__ceSeeds[dcoFromCE]

                    if noRLOFafterCEE:                                                                                  # if required
                        rlofSeeds = np.unique(dcoCEseeds[self.__immediateRLOF[dcoFromCE]])
                        rlofMask  = np.logical_not(np.in1d(self.__allDCOseeds, rlofSeeds))
                    else:
                        rlofMask = np.repeat(True, len(self.__allDCOseeds))

                    if pessimisticCE:                                                                                   # if required
                        pessimisticSeeds = np.unique(dcoCEseeds[self.__optimisticCE[dcoFromCE]])
                        pessimisticMask  = np.logical_not(np.in1d(self.__allDCOseeds, pessimisticSeeds))
                    else:
                        pessimisticMask = np.repeat(True, len(self.__allDCOseeds))

            # set DCO mask
            DCOmask = typeMasks[DCOtype] * hubbleMask * rlofMask * pessimisticMask

            if self.__verbose: 
                timeTaken = time.time() - mark
                memUsage  = psutil.Process().memory_info().rss / GB
                print('         DCO masks set in {:.2f} seconds.  Process total memory usage = {:.6f} GB'.format(timeTaken, memUsage))

            # calculate population values based on DCO mask
            errStr = self.__CalculatePopulationValues(DCOmask, p_m1Minimum, p_m1Maximum, p_m2Minimum, useSampledMassRanges, p_BinaryFraction)
            if errStr is None:

                # set class member variables for all DCO masks
                self.__DCOmask         = DCOmask
                self.__BBHmask         = typeMasks[COMPAS_DCO_BBH]         * hubbleMask * rlofMask * pessimisticMask
                self.__BHNSmask        = typeMasks[COMPAS_DCO_BHNS]        * hubbleMask * rlofMask * pessimisticMask
                self.__DNSmask         = typeMasks[COMPAS_DCO_BNS]         * hubbleMask * rlofMask * pessimisticMask
                self.__CHE_BBHmask     = typeMasks[COMPAS_DCO_CHE_BBH]     * hubbleMask * rlofMask * pessimisticMask
                self.__NON_CHE_BBHmask = typeMasks[COMPAS_DCO_NON_CHE_BBH] * hubbleMask * rlofMask * pessimisticMask
                self.__ALL_TYPESmask   = typeMasks[COMPAS_DCO_ALL]         * hubbleMask * rlofMask * pessimisticMask
                self.__OPTIMISTICmask  = pessimisticMask

                # update class member variables to match parameters
                self.__DCOtype              = DCOtype
                self.__withinHubbleTimeFlag = withinHubbleTime
                self.__pessimisticCEFlag    = pessimisticCE
                self.__noRLOFafterCEEFlag   = noRLOFafterCEE
                self.__useSampledMassRanges = useSampledMassRanges

        if errStr is not None: self.__verbose = saveVerbose                                                                     # error - restore verbose setting

        return errStr


    """
    COMPAS::__CalculatePopulationValues()

    Calculate population values for use throughout.
    
    Updates the following class member variables if no error occurs:
        __chirpMasses
        __delayTime
        __ETAs
        __mass1
        __mass2
        __massEvolvedPerBinary
        __maskedDCOSeeds
        __nAllDCOs
        __nDCOs
        __Zsystems

    Args:
        p_DCOmask              : 1D ARRAY of BOOL : Mask for DCO data - masks DCO types required
                                                    Default = None - an error
        p_m1Minimum            : FLOAT            : COMPAS minimum value for primary star mass
                                                    Must be >= 0.0, default = None: self.__m1Minimum will be used
        p_m1Maximum            : FLOAT            : COMPAS maximum value for primary star mass
                                                    Must be >= p_m1Minimum, default = None: self.__m1Maximum will be used
        p_m2Minimum            : FLOAT            : COMPAS minimum value for secondary star mass
                                                    Must be >= 0.0, default = None: self.__m2Minimum will be used
        p_UseSampledMassRanges : BOOL             : Flag to indicate that the mass ranges sampled by COMPAS should be used instead of
                                                    hard cuts at p_m1Minimum, p_m1Maximum, and p_m2Minimum.
                                                    Default = None, self.__useSampledMassRanges will be used
        p_BinaryFraction       : FLOAT            : Fraction of stars that form binaries
                                                    Must be >= 0.0 and <= 1.0, default = None: self.__binaryFraction will be used

    Returns:
        1 : STRING : Error string - set if an error occurred.
                     'None' if no error occurred.
    """
    def __CalculatePopulationValues(self, 
                                    p_DCOmask              = None, 
                                    p_m1Minimum            = None,
                                    p_m1Maximum            = None,
                                    p_m2Minimum            = None,
                                    p_UseSampledMassRanges = None,
                                    p_BinaryFraction       = None):

        func = self.__className + '::__CalculatePopulationValues() '

        if self.__verbose:
            print('         ...Calculating population values')
            mark = time.time()

        # default return value
        errStr = None

        # set local variables from parameter values
        # Doing it this way ensures that we have reasonable values for all parameters and allows users to call the function
        # with positional parameters and just set any parameters they want to default to 'None' rather than specify the value
        DCOmask              = p_DCOmask
        useSampledMassRanges = self.__useSampledMassRanges if p_UseSampledMassRanges is None else p_UseSampledMassRanges

        # check parameter types and values
        if p_DCOmask is None:                                    errStr = func + 'No DCO mask specified'
        elif not (isinstance(p_DCOmask, np.ndarray) and \
             p_DCOmask.ndim == 1 and p_DCOmask.dtype == 'bool'): errStr = func + 'DCO mask must be specified as a 1D ARRAY of BOOL'
        elif p_DCOmask.shape[0] < 1:                             errStr = func + 'No DCO mask specified'
        elif not isinstance(useSampledMassRanges, bool):         errStr = func + '\'Use sampled mass ranges\'' + tmpStr
        else:                                                                                                       # all parameters ok

            # use mass ranges sampled by COMPAS if requested
            if useSampledMassRanges:                                                                                # use mass ranges sampled by COMPAS instead of hard cuts?
                m1Minimum = np.min(self.__mass1[self.__mass1 != self.__mass2])                                      # yes - don't include masses set equal through RLOF at ZAMS
                m1Maximum = np.max(self.__mass1)
                m2Minimum = np.min(self.__mass2)
            else:                                                                                                   # no - use parameters passed
                m1Minimum = p_m1Minimum
                m1Maximum = p_m1Maximum
                m2Minimum = p_m2Minimum

            # calculate star forming mass evolved per binary
            massEvolvedPerBinary, errStr = self.__CalculateSFMassPerBinaryKroupa(m1Minimum, m1Maximum, m2Minimum, p_BinaryFraction)
            if errStr is None:

                # set population values

                self.__massEvolvedPerBinary = massEvolvedPerBinary

                self.__nAllDCOs             = self.__allDCOseeds.shape[0]
                self.__maskedDCOSeeds       = self.__allDCOseeds[p_DCOmask]
                self.__nDCOs                = self.__maskedDCOSeeds.shape[0]

                self.__mass1                = self.__allMass1[p_DCOmask]
                self.__mass2                = self.__allMass2[p_DCOmask]

                self.__Zsystems             = self.__zamsZ[np.in1d(self.__sysSeeds, self.__maskedDCOSeeds)]
                self.__delayTime            = np.add(self.__formationTime[p_DCOmask], self.__coalescenceTime[p_DCOmask])

                m1xm2                       = self.__mass1 * self.__mass2
                m1_m2                       = self.__mass1 + self.__mass2
                self.__chirpMasses          = m1xm2**(3.0 / 5.0) / m1_m2**(1.0 / 5.0)
                self.__ETAs                 = m1xm2 / m1_m2**2

        if self.__verbose: 
            timeTaken = time.time() - mark
            memUsage  = psutil.Process().memory_info().rss / GB
            print('            Population values calculated in {:.2f} seconds.  Process total memory usage = {:.6f} GB'.format(timeTaken, memUsage))

        return errStr


    """
    COMPAS::__CDF_IMF()

    Calculate the fraction of stellar mass between 0 and masses specifed for a three-part broken power law.
    Defaults to Kroupa IMF (Kroupa 2001) if power law parameters are not specified.

    This function is implemented as a recursive function in ClassCOMPAS.py - here it is not.

    Args:
        p_Masses : 1D ARRAY of FLOAT : Masses at which to calculate fractions
                                       Must be positive, default = None
        p_M1     : FLOAT             : Lower bound of lower line segment of three-part power law
                                       Must be >= 0.0, default = KROUPA_BREAK_0
        p_M2     : FLOAT             : First break point: lower bound of middle line segment of three-part power law
                                       Must be > p_M1, default = KROUPA_BREAK_1
        p_M3     : FLOAT             : Second break point: lower bound of upper line segment of three-part power law
                                       Must be > p_M2, default = KROUPA_BREAK_2
        p_M4     : FLOAT             : Upper bound of upper line segment of three-part power law
                                       Must be > p_M3, default = KROUPA_BREAK_3
        p_S12    : FLOAT             : Slope of lower line segment of three-part power law
                                       Default = KROUPA_SLOPE_1
        p_S23    : FLOAT             : Slope of middle line segment of three-part power law.
                                       Default = KROUPA_SLOPE_2
        p_S34    : FLOAT             : Slope of upper line segment of three-part power law.
                                       Default = KROUPA_SLOPE_3

    Returns:
        1 : 1D ARRAY of FLOAT : Values of the IMF at p_Masses.
                                'None' if an error occurred
        2 : STRING            : Error string - set if an error occurred.
                                'None' if no error occurred.
    """
    def __CDF_IMF(self,
                  p_Masses = None, 
                  p_M1     = KROUPA_BREAK_0, 
                  p_M2     = KROUPA_BREAK_1, 
                  p_M3     = KROUPA_BREAK_2, 
                  p_M4     = KROUPA_BREAK_3,
                  p_S12    = KROUPA_SLOPE_1, 
                  p_S23    = KROUPA_SLOPE_2, 
                  p_S34    = KROUPA_SLOPE_3):

        func = self.__className + '::__CDF_IMF(): '

        # default return values
        CDF    = None
        errStr = None

        # set local variables from parameter values
        # Doing it this way ensures that we have reasonable values for all parameters and allows users to call the function
        # with positional parameters and just set any parameters they want to default to 'None' rather than specify the value
        masses = p_Masses
        m1     = KROUPA_BREAK_0 if p_M1  is None else p_M1
        m2     = KROUPA_BREAK_1 if p_M2  is None else p_M2
        m3     = KROUPA_BREAK_2 if p_M3  is None else p_M3
        m4     = KROUPA_BREAK_3 if p_M4  is None else p_M4
        s12    = KROUPA_SLOPE_1 if p_S12 is None else p_S12
        s23    = KROUPA_SLOPE_2 if p_S23 is None else p_S23
        s34    = KROUPA_SLOPE_3 if p_S34 is None else p_S34

        # check parameter types and values
        tmpStr1 = ' must be a floating point number'                                                                # common error string
        tmpStr2 = ' line segment of three-part power law'                                                           # common error string
        tmpStr3 = 'Lower bound of '                                                                                 # common error string

        if masses is None:                                              errStr = func + 'No masses specified'
        elif not (isinstance(masses, np.ndarray) and masses.ndim == 1): errStr = func + 'Masses must be specified as a 1D array'
        elif masses.shape[0] < 1:                                       errStr = func + 'No masses specified'
        elif not (masses.dtype == 'float64' and np.all(masses > 0.0)):  errStr = func + 'Masses must be positive floating point numbers'
        elif not (isinstance(m1, float) and m1 >= 0.0):                 errStr = func + tmpStr3 + 'lower' + tmpStr2 + tmpStr1 + ' >= 0.0'
        elif not (isinstance(m2, float) and m2 > m1):                   errStr = func + tmpStr3 + 'middle' + tmpStr2 + tmpStr1 + ' >= ' + tmpStr3 + 'lower line segment'
        elif not (isinstance(m3, float) and m3 > m2):                   errStr = func + tmpStr3 + 'upper' + tmpStr2 + tmpStr1 + ' >= ' + tmpStr3 + 'middle line segment'
        elif not (isinstance(m4, float) and m4 > m3):                   errStr = func + 'Upper bound of upper' + tmpStr2 + tmpStr1 + ' >= ' + tmpStr3 + 'upper line segment'
        elif not isinstance(s12, float):                                errStr = func + 'Slope of lower' + tmpStr2 + tmpStr1
        elif not isinstance(s23, float):                                errStr = func + 'Slope of middle' + tmpStr2 + tmpStr1
        elif not isinstance(s34, float):                                errStr = func + 'Slope of upper' + tmpStr2 + tmpStr1
        else:

            # calculate common variables
            one_s12    = 1.0 - s12
            one_s23    = 1.0 - s23
            one_s34    = 1.0 - s34
            s23_s12    = s23 - s12
            s34_s23    = s34 - s23
            m1_1_s12   = m1**one_s12
            m2_1_s12   = m2**one_s12
            m2_1_s23   = m2**one_s23
            m2_s23_s12 = m2**s23_s12
            m3_1_s23   = m3**one_s23
            m3_1_s34   = m3**one_s34
            m3_s34_s23 = m3**s34_s23
            m4_1_s34   = m4**one_s34

            # calculate normalisation constants that ensure the IMF is continuous
            b1 = 1.0 / ((m2_1_s12 - m1_1_s12) / one_s12 + m2_s23_s12 * (m3_1_s23 - m2_1_s23) / one_s23 + m2_s23_s12 * m3_s34_s23 * (m4_1_s34 - m3_1_s34) / one_s34)
            b2 = b1 * m2_s23_s12
            b3 = b2 * m3_s34_s23

            # calculate IMF values              
            CDF                                            = np.zeros(len(masses))

            CDF[np.logical_and(masses >= m1, masses < m2)] = b1 / one_s12 * (masses[np.logical_and(masses >= m1, masses < m2)]**one_s12 - m1_1_s12)

            CDF[np.logical_and(masses >= m2, masses < m3)] = b1 / one_s12 * (m2_1_s12 - m1_1_s12) + \
                                                             b2 / one_s23 * (masses[np.logical_and(masses >= m2, masses < m3)]**one_s23 - m2_1_s23)

            CDF[np.logical_and(masses >= m3, masses < m4)] = b1 / one_s12 * (m2_1_s12 - m1_1_s12) + \
                                                             b2 / one_s23 * (m3_1_s23 - m2_1_s23) + \
                                                             b3 / one_s34 * (masses[np.logical_and(masses >= m3, masses < m4)]**one_s34 - m3_1_s34)

            CDF[masses >= m4]                              = np.ones(len(masses[masses >= m4]))

        return CDF, errStr


    """
    COMPAS::__SampleInitialMass()

    Samples initial masses by calculating the inverse CDF for a three-part broken power law.
    Defaults to Kroupa IMF (Kroupa 2001) if power law parameters are not specified.

    Args:   
        p_Samples : 1D ARRAY of FLOAT : Samples for which masses should be calculated.  Must be >= 0.0 and < 1.0.
                                        Must be >= 0.0 and < 1.0, default = None
        p_M1     : FLOAT              : Lower bound of lower line segment of three-part power law
                                        Must be >= 0.0, default = KROUPA_BREAK_0
        p_M2     : FLOAT              : First break point: lower bound of middle line segment of three-part power law
                                        Must be > p_M1, default = KROUPA_BREAK_1
        p_M3     : FLOAT              : Second break point: lower bound of upper line segment of three-part power law
                                        Must be > p_M2, default = KROUPA_BREAK_2
        p_M4     : FLOAT              : Upper bound of upper line segment of three-part power law
                                        Must be > p_M3, default = KROUPA_BREAK_3
        p_S12    : FLOAT              : Slope of lower line segment of three-part power law
                                        Default = KROUPA_SLOPE_1
        p_S23    : FLOAT              : Slope of middle line segment of three-part power law.
                                        Default = KROUPA_SLOPE_2
        p_S34    : FLOAT              : Slope of upper line segment of three-part power law.
                                        Default = KROUPA_SLOPE_3

    Returns:
        1 : 1D ARRAY of FLOAT : Mass values (inverse of IMF) at p_Samples (type matches type of p_Samples).
                                'None' if an error occurred.
        2 : STRING            : Error string - set if an error occurred.
                                'None' if no error occurred.
    """
    def __SampleInitialMass(self,
                            p_Samples = None,
                            p_M1      = KROUPA_BREAK_0, 
                            p_M2      = KROUPA_BREAK_1, 
                            p_M3      = KROUPA_BREAK_2, 
                            p_M4      = KROUPA_BREAK_3,
                            p_S12     = KROUPA_SLOPE_1, 
                            p_S23     = KROUPA_SLOPE_2, 
                            p_S34     = KROUPA_SLOPE_3):

        func = self.__className + '::__SampleInitialMass(): '

        # default return values
        masses = None
        errStr = None

        # set local variables from parameter values
        # Doing it this way ensures that we have reasonable values for all parameters and allows users to call the function
        # with positional parameters and just set any parameters they want to default to 'None' rather than specify the value
        samples = p_Samples
        m1      = KROUPA_BREAK_0 if p_M1  is None else p_M1
        m2      = KROUPA_BREAK_1 if p_M2  is None else p_M2
        m3      = KROUPA_BREAK_2 if p_M3  is None else p_M3
        m4      = KROUPA_BREAK_3 if p_M4  is None else p_M4
        s12     = KROUPA_SLOPE_1 if p_S12 is None else p_S12
        s23     = KROUPA_SLOPE_2 if p_S23 is None else p_S23
        s34     = KROUPA_SLOPE_3 if p_S34 is None else p_S34

        # check parameter types and values
        tmpStr1 = ' must be a floating point number'                                                                # common error string
        tmpStr2 = ' line segment of three-part power law'                                                           # common error string
        tmpStr3 = 'Lower bound of '                                                                                 # common error string

        if samples is None:                                               errStr = func + 'No samples specified'
        elif not (isinstance(samples, np.ndarray) and samples.ndim == 1): errStr = func + 'Samples must be specified as a 1D array'
        elif samples.shape[0] < 1:                                        errStr = func + 'No samples specified'
        elif not (samples.dtype == 'float64' and \
             np.all(samples >= 0.0) and np.all(samples < 1.0)):           errStr = func + 'Samples must be floating point numbers >= 0.0 and < 1.0'
        elif not (isinstance(m1, float) and m1 >= 0.0):                   errStr = func + tmpStr3 + 'lower' + tmpStr2 + tmpStr1 + ' >= 0.0'
        elif not (isinstance(m2, float) and m2 > m1):                     errStr = func + tmpStr3 + 'middle' + tmpStr2 + tmpStr1 + ' >= ' + tmpStr3 + 'lower line segment'
        elif not (isinstance(m3, float) and m3 > m2):                     errStr = func + tmpStr3 + 'upper' + tmpStr2 + tmpStr1 + ' >= ' + tmpStr3 + 'middle line segment'
        elif not (isinstance(m4, float) and m4 > m3):                     errStr = func + 'Upper bound of upper' + tmpStr2 + tmpStr1 + ' >= ' + tmpStr3 + 'upper line segment'
        elif not isinstance(s12, float):                                  errStr = func + 'Slope of lower' + tmpStr2 + tmpStr1
        elif not isinstance(s23, float):                                  errStr = func + 'Slope of middle' + tmpStr2 + tmpStr1
        elif not isinstance(s34, float):                                  errStr = func + 'Slope of upper' + tmpStr2 + tmpStr1
        else:
                # calculate common variables
                one_s12 = 1.0 - s12
                one_s23 = 1.0 - s23
                one_s34 = 1.0 - s34
                s23_s12 = s23 - s12
                s34_s23 = s34 - s23
                m1_one_s12 = m1**one_s12
                m2_one_s23 = m2**one_s23
                m2_s23_s12 = m2**s23_s12
                m3_one_s34 = m3**one_s34
                m3_s34_s23 = m3**s34_s23

                # calculate normalisation constants that ensure the IMF is continuous
                b1 = 1.0 / ((m2**one_s12 - m1_one_s12) / one_s12 + m2_s23_s12 * (m3**one_s23 - m2_one_s23) / one_s23 + m2_s23_s12 * m3_s34_s23 * (m4**one_s34 - m3_one_s34) / one_s34)
                b2 = b1 * m2_s23_s12
                b3 = b2 * m3_s34_s23

                # find the probabilities at which the gradient changes
                CDF, errStr = self.__CDF_IMF(np.array([m1, m2, m3, m4]))
                if errStr is None:
                    # calculate masses
                    masses                                                      = np.zeros(len(samples))
                    masses[np.logical_and(samples > CDF[0], samples <= CDF[1])] = np.power(one_s12 / b1 * (samples[np.logical_and(samples > CDF[0], samples <= CDF[1])] - CDF[0]) + m1_one_s12, 1.0 / one_s12)
                    masses[np.logical_and(samples > CDF[1], samples <= CDF[2])] = np.power(one_s23 / b2 * (samples[np.logical_and(samples > CDF[1], samples <= CDF[2])] - CDF[1]) + m2_one_s23, 1.0 / one_s23)
                    masses[np.logical_and(samples > CDF[2], samples <= CDF[3])] = np.power(one_s34 / b3 * (samples[np.logical_and(samples > CDF[2], samples <= CDF[3])] - CDF[2]) + m3_one_s34, 1.0 / one_s34)
    
        return masses, errStr


    """
    COMPAS::__CalculateSFMassPerBinaryKroupa()

    Compute the star forming mass evolved per binary analytically, using Ilya Mandel's derivation.
    Assumes Kroupa IMF (Kroupa 2001):

        p(M) \propto M^-0.3 for M between m1 and m2
        p(M) \propto M^-1.3 for M between m2 and m3
        p(M) = alpha * M^-2.3 for M between m3 and m4

    Updates the following class member variables to parameter values if no error occurs:
        __m1Minimum
        __m1Maximum
        __m2Minimum
        __binaryFraction

    Args:
        p_m1Minimum      : FLOAT  : COMPAS minimum value for primary star mass
                                    Must be >= 0.0, default = None: self.__m1Minimum will be used
        p_m1Maximum      : FLOAT  : COMPAS maximum value for primary star mass
                                    Must be >= p_m1Minimum, default = None: self.__m1Maximum will be used
        p_m2Minimum      : FLOAT  : COMPAS minimum value for secondary star mass
                                    Must be >= 0.0, default = None: self.__m2Minimum will be used
        p_BinaryFraction : FLOAT  : Fraction of stars that form binaries
                                    Must be >= 0.0 and <= 1.0, default = None: self.__binaryFraction will be used

    Returns:
        1 : FLOAT  : Average star forming mass evolved per binary.
                     'None' if an error occurred.
        2 : STRING : Error string - set if an error occurred.
                     'None' if no error occurred.
    """
    def __CalculateSFMassPerBinaryKroupa(self, p_m1Minimum = None, p_m1Maximum = None, p_m2Minimum = None, p_BinaryFraction = None):

        func = self.__className + '::__CalculateSFMassPerBinaryKroupa(): '

        # default return values
        SFmass = None
        errStr = None

        # set local variables from parameter values
        # Doing it this way ensures that we have reasonable values for all parameters and allows users to call the function
        # with positional parameters and just set any parameters they want to default to 'None' rather than specify the value
        m1Minimum      = self.__m1Minimum      if p_m1Minimum      is None else p_m1Minimum 
        m1Maximum      = self.__m1Maximum      if p_m1Maximum      is None else p_m1Maximum 
        m2Minimum      = self.__m2Minimum      if p_m2Minimum      is None else p_m2Minimum 
        binaryFraction = self.__binaryFraction if p_BinaryFraction is None else p_BinaryFraction 

        # check parameter types and values
        tmpStr = ' must be a floating point number'                                                                 # common error string

        if   not (isinstance(m1Minimum, float) and m1Minimum >= KROUPA_BREAK_2): errStr = func + 'Minimum primary star mass' + tmpStr + ' >= {:.3f}'.format(KROUPA_BREAK_2)
        elif not (isinstance(m1Maximum, float) and m1Maximum >= m1Minimum):      errStr = func + 'Maximum primary star mass must be >= minimum primary star mass'
        elif not (isinstance(m2Minimum, float) and m2Minimum >= 0.0):            errStr = func + 'Minimum secondary star mass' + tmpStr + ' >= 0.0'
        elif not (isinstance(binaryFraction, float) and \
                 binaryFraction >= 0.0 and binaryFraction <= 1.0):               errStr = func + 'Binary fraction' + tmpStr + ' >= 0.0 and <= 1.0'
        else:

            # set IMF mass bounds to Kroupa values
            m1 = KROUPA_BREAK_0
            m2 = KROUPA_BREAK_1
            m3 = KROUPA_BREAK_2
            m4 = KROUPA_BREAK_3

            # constants used below are:
            #
            #  0.3 = KROUPA_SLOPE_1
            #  0.7 = 1 - KROUPA_SLOPE_1
            #  1.3 = KROUPA_SLOPE_2
            # -0.3 = 1 - KROUPA_SLOPE_2
            #  2.3 = KROUPA_SLOPE_3
            # -1.3 = 1 - KROUPA_SLOPE_3
            # -2.3 = -KROUPA_SLOPE_3
            #
            #  1.7 comes from integrating Mp(M) dM when p(M) \propto M^{-0.3}
            #  1.5 comes from the avergae mass of all binaries being 1.5 times larger than the average mass of stars
            #
            # ordinarily I might set variables for the constants so they can be easily changed, and set variables equal
            # to common calculated values (e.g. m2_7 = m2**0.7) just so we only compute them once, but that might affect
            # readability here...
   
            alpha = (-(m4**(-1.3) - m3**(-1.3)) / 1.3 - (m3**(-0.3) - m2**(-0.3)) / (m3 * 0.3) + (m2**0.7 - m1**0.7) / (m2 * m3 * 0.7))**(-1) # 1.0/(...) might be faster...

            # average mass of stars
            mAvg = alpha * (-(m4**(-0.3) - m3**(-0.3)) / 0.3 + (m3**0.7 - m2**0.7) / (m3 * 0.7) + (m2**1.7 - m1**1.7)/(m2 * m3 * 1.7))

            # fraction of binaries that COMPAS simulates
            fCOMPAS = -alpha / 1.3 * (m1Maximum**(-1.3) - m1Minimum**(-1.3)) + alpha * m2Minimum / 2.3 * (m1Maximum**(-2.3) - m1Minimum**(-2.3))

            # mass represented by each binary simulated by COMPAS
            SFmass = (1.0 / fCOMPAS) * mAvg * (1.5 + (1.0 - binaryFraction) / binaryFraction)

            # update class member variables to match parameters
            self.__m1Minimum      = m1Minimum
            self.__m1Maximum      = m1Maximum
            self.__m2Minimum      = m2Minimum
            self.__binaryFraction = binaryFraction

        return SFmass, errStr


    # public interface
    
    def Initialise(self, 
                   p_m1Minimum            = None,
                   p_m1Maximum            = None,
                   p_m2Minimum            = None,
                   p_UseSampledMassRanges = None,
                   p_BinaryFraction       = None,
                   p_DCOtype              = None,
                   p_WithinHubbleTime     = None,
                   p_PessimisticCE        = None,
                   p_NoRLOFafterCEE       = None,
                   p_Verbose              = None):
        return self.__Initialise(p_m1Minimum, p_m1Maximum, p_m2Minimum, p_UseSampledMassRanges, p_BinaryFraction, p_DCOtype, p_WithinHubbleTime, p_PessimisticCE, p_NoRLOFafterCEE, p_Verbose)


    # getters

    """
    The following getters return:
    
        - an array value at the index specified, if index parameter is specified
        - the complete array if the index parameter is not specified

    Each will return an error if an invalid index is specified
    """

    def ChirpMasses(self, p_Index = None):

        func = self.__className + "::ChirpMasses(): "

        # default return values
        Mc     = None
        errStr = None

        if p_Index is None:
            Mc = self.__chirpMasses
        elif not (isinstance(p_Index, int) and p_Index >= 0 and p_Index < self.__chirpMasses.shape[0]):
            errStr = func + 'Index value must be an integer >= 0 and < {:d}'.format(self.__chirpMasses.shape[0])
        else:
            Mc = self.__chirpMasses[p_Index]

        return Mc, errStr


    def DelayTime(self, p_Index = None):

        func = self.__className + "::DelayTime(): "

        # default return values
        delay  = None
        errStr = None

        if p_Index is None:
            delay = self.__delayTime
        elif not (isinstance(p_Index, int) and p_Index >= 0 and p_Index < self.__delayTime.shape[0]):
            errStr = func + 'Index value must be an integer >= 0 and < {:d}'.format(self.__delayTime.shape[0])
        else:
            delay = self.__delayTime[p_Index]

        return delay, errStr


    def ETAs(self, p_Index = None):

        func = self.__className + "::ETAs(): "

        # default return values
        eta    = None
        errStr = None

        if p_Index is None:
            eta = self.__ETAs
        elif not (isinstance(p_Index, int) and p_Index >= 0 and p_Index < self.__ETAs.shape[0]):
            errStr = func + 'Index value must be an integer >= 0 and < {:d}'.format(self.__ETAs.shape[0])
        else:
            eta = self.__ETAs[p_Index]

        return eta, errStr


    def PrimaryMass(self, p_Index = None):

        func = self.__className + "::PrimaryMass(): "

        # default return values
        mass   = None
        errStr = None

        if p_Index is None:
            mass = self.__mass1
        elif not (isinstance(p_Index, int) and p_Index >= 0 and p_Index < self.__mass1.shape[0]):
            errStr = func + 'Index value must be an integer >= 0 and < {:d}'.format(self.__mass1.shape[0])
        else:
            mass = self.__mass1[p_Index]

        return mass, errStr


    def SecondaryMass(self, p_Index = None):

        func = self.__className + "::SecondaryMass(): "

        # default return values
        mass   = None
        errStr = None

        if p_Index is None:
            mass = self.__mass2
        elif not (isinstance(p_Index, int) and p_Index >= 0 and p_Index < self.__mass2.shape[0]):
            errStr = func + 'Index value must be an integer >= 0 and < {:d}'.format(self.__mass2.shape[0])
        else:
            mass = self.__mass2[p_Index]

        return mass, errStr


    def Zsystems(self, p_Index = None):

        func = self.__className + "::Zsystems(): "

        # default return values
        Z      = None
        errStr = None

        if p_Index is None:
            Z = self.__Zsystems
        elif not (isinstance(p_Index, int) and p_Index >= 0 and p_Index < self.__Zsystems.shape[0]):
            errStr = func + 'Index value must be an integer >= 0 and < {:d}'.format(self.__Zsystems.shape[0])
        else:
            Z = self.__Zsystems[p_Index]

        return Z, errStr


    # the following getters take no parameters

    def DCOtype(self)             : return self.__DCOtype
    def MassEvolvedPerBinary(self): return self.__massEvolvedPerBinary
    def MaxLogZ(self)             : return self.__maxLogZ
    def MinLogZ(self)             : return self.__minLogZ
    def nAllDCOs(self)            : return self.__nAllDCOs
    def nBinaries(self)           : return self.__nBinaries
    def nDCOs(self)               : return self.__nDCOs
    def Version(self)             : return self.__version
    def RunStart(self)            : return self.__runStart
    def RunEnd(self)              : return self.__runEnd


"""
Class CosmicIntegration

Provides cosmoslogical integration support functions.

Public interface (public functions declared at end of class):

    CosmicIntegration::CosmicIntegration() : Class constructor
    CosmicIntegration::FindRates()         : Finds formation, merger, and detection rates

        Getters:
            COMPAS()           : Returns instantiated COMPAS object
            Redshifts()        : Returns redshifts at which various values were calculated
            SelectionEffects() : Returns instantiated SelectionEffects object
            SFmass()           : Returns representative star formation mass

Adapted from FastCosmicIntegration.py
"""
class CosmicIntegration:

    """
    CosmicIntegration::CosmicIntegration

    Class constructor

    Args:

        Class SelectionEffects related parameters: see class SelectionEffects constructor
        Class COMPAS related parameters: see class COMPAS constructor

        Class CosmicIntegration related parameters:

        p_MaxRedshift          : FLOAT         : Maximum formation redshift
                                                 Must be >= 0.0, default = CI_maxRedshift
        p_MaxRedshiftDetection : FLOAT         : Maximum redshift to calculate detection rates
                                                 Must be >= 0.0 and <= p_MaxRedshift, default = CI_maxRedshiftDetection
        p_RedshiftStep         : FLOAT         : Size of step to take in redshift
                                                 Must be > 0.0, default = CI_redshiftStep
        p_MinLogZ              : FLOAT         : Minimum logZ for dPdlogZ calculation (influences normalization)
                                                 Must be >= CI_MinLogZ and <= CI_MaxLogZ, default = CI_MinLogZ
        p_MaxLogZ              : FLOAT         : Maximum logZ for dPdlogZ calculation (influences normalization)
                                                 Must be > p_MinLogZ and <= CI_MaxLogZ, default = CI_MaxLogZ
        p_StepLogZ             : FLOAT         : Size of logZ steps to take in finding a Z range
                                                 Must be positive, default = CI_StepLogZ
        p_MU_0                 : FLOAT         : Location (mean in normal) at redshift 0
                                                 Must be >= 0.0 and <= 1.0, default = NEIJSSEL_MU_0
        p_MU_Z                 : FLOAT         : Redshift scaling/evolution of the mean metallicity
                                                 Default = NEIJSSEL_MU_Z
        p_SIGMA_0              : FLOAT         : Scale (variance in normal) at redshift 0
                                                 Must be >= 0.0, default = NEIJSSEL_SIGMA_0
        p_SIGMA_Z              : FLOAT         : Redshift scaling of the scale (variance in normal)
                                                 Default = NEIJSSEL_SIGMA_Z
        p_ALPHA                : FLOAT         : Shape (skewness) of metallicity density distribution, p_ALPHA = 0 retrieves log-normal dist)
                                                 Default = NEIJSSEL_ALPHA
        p_SFR_a                : FLOAT         : Scale factor for SFR (see Neijssel+2019)
                                                 Default = NEIJSSEL_SFR_a
        p_SFR_b                : FLOAT         : Scale factor for SFR (see Neijssel+2019)
                                                 Default = NEIJSSEL_SFR_b
        p_SFR_c                : FLOAT         : Scale factor for SFR (see Neijssel+2019)
                                                 Default = NEIJSSEL_SFR_c
        p_SFR_d                : FLOAT         : Scale factor for SFR (see Neijssel+2019)
                                                 Default = NEIJSSEL_SFR_d
        p_McBins               : LIST of LIST  : Chirp mass bins to be used for binning rates - fixed or variable width
                                                 Default = None

                                                 A value of 'None', or an empty list ([]), indicates no binning - FindRates() will be calculated per binary

                                                 If not 'None' or an empty LIST, p_McBins must be a LIST and must contain (in order):

                                                     binWidths     : LIST of FLOAT : Width of each chirp mas bin (Msun)
                                                     binRightEdges : LIST of FLOAT : Chirp mass at right (upper) edge of each chirp mass bin

                                                     (see class UTILS for helper functions to create and check chirp mass bins)
        p_PerBinaryFiles       : BOOL          : Flag to specify whether per-binary files (mergers and yields) are required
                                                 Default = CI_PerBinaryFiles
        p_PerBinaryFilesPath   : STRING        : Path at which per-binary files should be created (if necessary)
                                                 Default = None: current directory will be used
        p_CreatorDetails       : STRING        : Creator details for per-binary file meta-data (if necessary)
                                                 (usually name and email address - can contain commas for CSV file)
                                                 Default = None
                                                 A value of 'None', or an empty string (''), indicates that no creator details should be recorded in per-binary files       
        p_RandomSeed           : INT           : Random seed to initialize numpy randum number generator
                                                 Must be >= 0, default = None: use numpy default seed (note: means results will (very likely) not be reproducible)
        p_Verbose              : BOOL          : Flag to indicate if diagnostics/stats should be printed to console
                                                 Default = VERBOSE_DEFAULT

    Errors here are handled by raising an exception - the constructor can't return an error to the caller.
    """
    def __init__(self,
                 # CosmicIntegration related parameters
                 p_MaxRedshift          = CI_maxRedshift,
                 p_MaxRedshiftDetection = CI_maxRedshiftDetection,
                 p_RedshiftStep         = CI_redshiftStep,
                 p_MinLogZ              = CI_MinLogZ,
                 p_MaxLogZ              = CI_MaxLogZ,
                 p_StepLogZ             = CI_StepLogZ,
                 p_MU_0                 = NEIJSSEL_MU_0, 
                 p_MU_Z                 = NEIJSSEL_MU_Z, 
                 p_SIGMA_0              = NEIJSSEL_SIGMA_0,
                 p_SIGMA_Z              = NEIJSSEL_SIGMA_Z,
                 p_ALPHA                = NEIJSSEL_ALPHA,
                 p_SFR_a                = NEIJSSEL_SFR_a,
                 p_SFR_b                = NEIJSSEL_SFR_b,
                 p_SFR_c                = NEIJSSEL_SFR_c,
                 p_SFR_d                = NEIJSSEL_SFR_d,
                 p_McBins               = None,
                 p_PerBinaryFiles       = CI_PerBinaryFiles,
                 p_PerBinaryFilesPath   = None,
                 p_CreatorDetails       = None,

                 # SelectionEffects related parameters
                 p_NoiseFilepath        = SE_NoiseFilepath, 
                 p_NoiseFilename        = SE_NoiseFilename, 
                 p_Sensitivity          = SE_sensitivity, 
                 p_MinThetas            = SE_minThetas,
                 p_NumThetas            = SE_numThetas,
                 p_SNRthreshold         = SE_SNRthreshold,
                 p_McMax                = SE_McMax,
                 p_McStep               = SE_McStep,
                 p_ETAmax               = SE_ETAmax,
                 p_ETAstep              = SE_ETAstep,
                 p_SNRmax               = SE_SNRmax,
                 p_SNRstep              = SE_SNRstep,

                 # COMPAS related parameters
                 p_COMPASfilepath       = COMPAS_filepath,
                 p_COMPASfilename       = COMPAS_filename,
                 p_m1Minimum            = COMPAS_m1Minimum, 
                 p_m1Maximum            = COMPAS_m1Maximum, 
                 p_m2Minimum            = COMPAS_m2Minimum, 
                 p_UseSampledMassRanges = COMPAS_UseSampledMassRanges,
                 p_BinaryFraction       = COMPAS_binaryFraction,
                 p_DCOtype              = COMPAS_DCO_ALL,
                 p_WithinHubbleTime     = COMPAS_flag_withinHubbleTime,
                 p_PessimisticCE        = COMPAS_flag_pessimisticCE,
                 p_NoRLOFafterCEE       = COMPAS_flag_noRLOFafterCEE,

                 # miscellaneos parameters
                 p_RandomSeed           = None,
                 p_Verbose              = VERBOSE_DEFAULT):

        self.__className = 'CosmicIntegration'
        func = self.__className + '::CosmicIntegration(): '

        self.__verbose = VERBOSE_DEFAULT if p_Verbose is None else p_Verbose

        if self.__verbose:
            print('\n...Initialising CosmicIntegration')
            mark = time.time()

        # set local variables from parameter values
        # only need cosmis integration parameters - others are checked by called functions
        # Doing it this way ensures that we have reasonable values for all parameters and allows users to call the function
        # with positional parameters and just set any parameters they want to default to 'None' rather than specify the value
        maxRedshift          = CI_maxRedshift          if p_MaxRedshift          is None else p_MaxRedshift
        maxRedshiftDetection = CI_maxRedshiftDetection if p_MaxRedshiftDetection is None else p_MaxRedshiftDetection
        redshiftStep         = CI_redshiftStep         if p_RedshiftStep         is None else p_RedshiftStep
        minLogZ              = CI_MinLogZ              if p_MinLogZ              is None else p_MinLogZ
        maxLogZ              = CI_MaxLogZ              if p_MaxLogZ              is None else p_MaxLogZ
        stepLogZ             = CI_StepLogZ             if p_StepLogZ             is None else p_StepLogZ
        mu_0                 = NEIJSSEL_MU_0           if p_MU_0                 is None else p_MU_0 
        mu_Z                 = NEIJSSEL_MU_Z           if p_MU_Z                 is None else p_MU_Z 
        sigma_0              = NEIJSSEL_SIGMA_0        if p_SIGMA_0              is None else p_SIGMA_0
        sigma_Z              = NEIJSSEL_SIGMA_Z        if p_SIGMA_Z              is None else p_SIGMA_Z
        alpha                = NEIJSSEL_ALPHA          if p_ALPHA                is None else p_ALPHA 
        sfr_a                = NEIJSSEL_SFR_a          if p_SFR_a                is None else p_SFR_a
        sfr_b                = NEIJSSEL_SFR_b          if p_SFR_b                is None else p_SFR_b
        sfr_c                = NEIJSSEL_SFR_c          if p_SFR_c                is None else p_SFR_c
        sfr_d                = NEIJSSEL_SFR_d          if p_SFR_d                is None else p_SFR_d
        McBins               = p_McBins                if p_McBins               is None else (None if (isinstance(p_McBins, list) and len(p_McBins) == 0) else p_McBins)
        perBinaryFiles       = CI_PerBinaryFiles       if p_PerBinaryFiles       is None else p_PerBinaryFiles
        perBinaryFilesPath   = CI_PerBinaryFilesPath   if p_PerBinaryFilesPath   is None else p_PerBinaryFilesPath
        creatorDetails       = p_CreatorDetails
        randomSeed           = p_RandomSeed


        # check parameter types and values
        tmpStr = ' must be a floating point number'                                                                 # common error string

        if   not (isinstance(maxRedshift, float) and maxRedshift >= 0.0):              raise RuntimeError(func + 'Maximum formation redshift' + tmpStr + ' >= 0.0')
        elif not (isinstance(maxRedshiftDetection, float) and \
                 maxRedshiftDetection >= 0.0 and maxRedshiftDetection <= maxRedshift): raise RuntimeError(func + 'Maximum detection redshift' + tmpStr + ' >= 0.0 and <= maximum formation redshift ({:.6f})'.format(maxRedshift))
        elif not (isinstance(redshiftStep, float) and redshiftStep > 0.0):             raise RuntimeError(func + 'Redshift step size' + tmpStr + ' > 0.0')
        elif not (isinstance(minLogZ, float) and \
                 minLogZ >= CI_MinLogZ and minLogZ <= CI_MaxLogZ):                     raise RuntimeError(func + 'Minimum logZ for dPdlogZ calculation' + tmpStr + ' >= {:.6f} and <= {:.6f}'.format(CI_MinLogZ, CI_MaxLogZ))
        elif not (isinstance(maxLogZ, float) and \
                 maxLogZ > minLogZ and maxLogZ <= CI_MaxLogZ):                         raise RuntimeError(func + 'Maximum logZ for dPdlogZ calculation' + tmpStr + ' > {:.6f} and <= {:.6f}'.format(minLogZ, CI_MaxLogZ))
        elif not (isinstance(stepLogZ, float) and stepLogZ > 0.0):                     raise RuntimeError(func + 'logZ step size' + tmpStr + ' > 0.0')
        elif not (isinstance(mu_0, float) and mu_0 >= 0.0 and mu_0 <= 1.0):            raise RuntimeError(func + 'Redshift mean at z=0 (mu_0)' + tmpStr + ' >= 0.0 and <= 1.0')
        elif not isinstance(mu_Z, float):                                              raise RuntimeError(func + 'Redshift mean scale factor (mu_z)' + tmpStr)
        elif not (isinstance(sigma_0, float) and sigma_0 >= 0.0):                      raise RuntimeError(func + 'Redshift variance at z=0 (sigma_0)' + tmpStr + ' >= 0.0')
        elif not isinstance(sigma_Z, float):                                           raise RuntimeError(func + 'Redshift variance scale factor (sigma_Z)' + tmpStr)
        elif not isinstance(alpha, float):                                             raise RuntimeError(func + 'Distribution skew (alpha)' + tmpStr)
        elif not isinstance(sfr_a, float):                                             raise RuntimeError(func + 'SFR scale factor SFR_a' + tmpStr)
        elif not isinstance(sfr_b, float):                                             raise RuntimeError(func + 'SFR scale factor SFR_b' + tmpStr)
        elif not isinstance(sfr_c, float):                                             raise RuntimeError(func + 'SFR scale factor SFR_c' + tmpStr)
        elif not isinstance(sfr_d, float):                                             raise RuntimeError(func + 'SFR scale factor SFR_d' + tmpStr)
        elif McBins is not None and (thisErr := Utils.CheckBins(McBins) is not None):  raise RuntimeError(func + thisErr)
        elif perBinaryFiles is not None and not isinstance(perBinaryFiles, bool):      raise RuntimeError(func + 'When specified, flag for per-binary files must be a boolean')
        elif not isinstance(perBinaryFilesPath, str):                                  raise RuntimeError(func + 'When specified, per-binary files path must be a string')
        elif creatorDetails is not None and not isinstance(creatorDetails, str):       raise RuntimeError(func + 'When specified, creator details must be a string')
        elif randomSeed is not None and \
             not (isinstance(randomSeed, int) and randomSeed >= 0):                    raise RuntimeError(func + 'When specified, random seed must be an integer >= 0')
        else:                                                                                                       # all parameters ok

            # fix up per-binary files path if necessary
            if   len(perBinaryFilesPath) < 1:   perBinaryFilesPath = CI_PerBinaryFilesPath                          # handle empty path
            elif perBinaryFilesPath[-1] != '/': perBinaryFilesPath += '/'                                           # add path delimiter if necessary

            # construct a random number generator - initialise with random seed (may be 'None')
            self.__rng = np.random.default_rng(randomSeed)

            # set class member variables
            self.__maxRedshift          = maxRedshift 
            self.__maxRedshiftDetection = maxRedshiftDetection
            self.__redshiftStep         = redshiftStep
            self.__minLogZ              = minLogZ
            self.__maxLogZ              = maxLogZ
            self.__stepLogZ             = stepLogZ
            self.__mu_0                 = mu_0
            self.__mu_Z                 = mu_Z
            self.__sigma_0              = sigma_0
            self.__sigma_Z              = sigma_Z
            self.__alpha                = alpha
            self.__sfr_a                = sfr_a
            self.__sfr_b                = sfr_b
            self.__sfr_c                = sfr_c
            self.__sfr_d                = sfr_d
            
            self.__McBins               = McBins
            self.__perBinaryFiles       = perBinaryFiles
            self.__perBinaryFilesPath   = perBinaryFilesPath
            self.__creatorDetails       = creatorDetails
            self.__randomSeed           = randomSeed

            self.__SFmass               = None                                                                      # calculated later

            self.__today = date.today().strftime("%Y-%m-%d")                                                        # get current date - so filenames match meta-data content   

            # initialise SelectionEffects class
            # no need to check parameters here - the class constructor will do that for us
            self.__SE = SelectionEffects(p_NoiseFilepath = p_NoiseFilepath,
                                         p_NoiseFilename = p_NoiseFilename,
                                         p_Sensitivity   = p_Sensitivity,
                                         p_MinThetas     = p_MinThetas,
                                         p_NumThetas     = p_NumThetas,
                                         p_SNRthreshold  = p_SNRthreshold,
                                         p_McMax         = p_McMax,
                                         p_McStep        = p_McStep,
                                         p_ETAmax        = p_ETAmax,
                                         p_ETAstep       = p_ETAstep,
                                         p_SNRmax        = p_SNRmax,
                                         p_SNRstep       = p_SNRstep,
                                         p_RNG           = self.__rng,
                                         p_Verbose       = p_Verbose)

            # initialise COMPAS class
            # no need to check parameters here - the class constructor will do that for us
            self.__COMPAS = COMPAS(p_COMPASfilepath       = p_COMPASfilepath,
                                   p_COMPASfilename       = p_COMPASfilename,
                                   p_m1Minimum            = p_m1Minimum, 
                                   p_m1Maximum            = p_m1Maximum, 
                                   p_m2Minimum            = p_m2Minimum,
                                   p_UseSampledMassRanges = p_UseSampledMassRanges, 
                                   p_BinaryFraction       = p_BinaryFraction,
                                   p_DCOtype              = p_DCOtype,
                                   p_WithinHubbleTime     = p_WithinHubbleTime,
                                   p_PessimisticCE        = p_PessimisticCE,
                                   p_NoRLOFafterCEE       = p_NoRLOFafterCEE,
                                   p_RNG                  = self.__rng,
                                   p_Verbose              = p_Verbose)

            # save current state for possible future restoration
            self.__SaveState()

        if self.__verbose: 
            timeTaken = time.time() - mark
            memUsage  = psutil.Process().memory_info().rss / GB
            print('   CosmicIntegration initialised in {:.2f} seconds.  Process total memory usage = {:.6f} GB\n'.format(timeTaken, memUsage))


    """
    CosmicIntegration::__SaveState()

    Saves the current state of the class
       - stores the current value of each class member variable
       - stores the SelectionEffects and COMPAS class instances 

    Timed at ~ 1.5E-5s - performance shouldn't be an issue, especially since
    this is only called at the creation of the CosmicIntegration class, and at
    the start of FindRates().  Restoring state is similar overhead, and is only
    called (once) in FindRates() if an error occurred.     
    """
    def __SaveState(self):

        # save values of class member Cosmic Integration related variables

        self.__saved_alpha                = self.__alpha
        self.__saved_creatorDetails       = self.__creatorDetails
        self.__saved_maxLogZ              = self.__maxLogZ
        self.__saved_maxRedshift          = self.__maxRedshift 
        self.__saved_maxRedshiftDetection = self.__maxRedshiftDetection
        self.__saved_McBins               = self.__McBins
        self.__saved_minLogZ              = self.__minLogZ
        self.__saved_mu_0                 = self.__mu_0
        self.__saved_mu_Z                 = self.__mu_Z
        self.__saved_perBinaryFiles       = self.__perBinaryFiles
        self.__saved_perBinaryFilesPath   = self.__perBinaryFilesPath
        self.__saved_randomSeed           = self.__randomSeed
        self.__saved_redshiftStep         = self.__redshiftStep
        self.__saved_rng                  = self.__rng
        self.__saved_SFmass               = self.__SFmass
        self.__saved_sfr_a                = self.__sfr_a
        self.__saved_sfr_b                = self.__sfr_b
        self.__saved_sfr_c                = self.__sfr_c
        self.__saved_sfr_d                = self.__sfr_d
        self.__saved_sigma_0              = self.__sigma_0
        self.__saved_sigma_Z              = self.__sigma_Z
        self.__saved_stepLogZ             = self.__stepLogZ
        self.__saved_today                = self.__today  
        self.__saved_verbose              = self.__verbose

        # save SelectionEffects and COMPAS class instances

        self.__saved_SE     = self.__SE
        self.__saved_COMPAS = self.__COMPAS


    """
    CosmicIntegration::__RestoreState()

    Restores the state of the class to the previously save state
       - restores the value of each class member variable
       - restores the SelectionEffects and COMPAS class instances      
    """
    def __RestoreState(self):

        # restore class member Cosmic Integration related variables to save values

        self.__alpha                = self.__saved_alpha
        self.__creatorDetails       = self.__saved_creatorDetails
        self.__maxLogZ              = self.__saved_maxLogZ
        self.__maxRedshift          = self.__saved_maxRedshift 
        self.__maxRedshiftDetection = self.__saved_maxRedshiftDetection
        self.__McBins               = self.__saved_McBins
        self.__minLogZ              = self.__saved_minLogZ
        self.__mu_0                 = self.__saved_mu_0
        self.__mu_Z                 = self.__saved_mu_Z
        self.__perBinaryFiles       = self.__saved_perBinaryFiles
        self.__perBinaryFilesPath   = self.__saved_perBinaryFilesPath
        self.__randomSeed           = self.__saved_randomSeed
        self.__redshiftStep         = self.__saved_redshiftStep
        self.__rng                  = self.__saved_rng
        self.__SFmass               = self.__saved_SFmass
        self.__sfr_a                = self.__saved_sfr_a
        self.__sfr_b                = self.__saved_sfr_b
        self.__sfr_c                = self.__saved_sfr_c
        self.__sfr_d                = self.__saved_sfr_d
        self.__sigma_0              = self.__saved_sigma_0
        self.__sigma_Z              = self.__saved_sigma_Z
        self.__stepLogZ             = self.__saved_stepLogZ
        self.__today                = self.__saved_today  
        self.__verbose              = self.__saved_verbose

        # restore SelectionEffects and COMPAS class instances

        self.__SE     = self.__saved_SE
        self.__COMPAS = self.__saved_COMPAS


    """ 
    CosmicIntegration::__CalculateRedshiftRelatedParams()

    Given limits on the redshift, create an array of redshifts, times, distances and volumes.
    
    Updates the following class member variables if no error occurs:
        __distances
        __nRedshifts
        __nRedshiftsDetection
        __nTimes
        __redshifts
        __shellVolumes
        __times

    Updates the following class member variables to parameter values if no error occurs:
        __maxRedshift
        __maxRedshiftDetection
        __redshiftStep

    Args:
        p_MaxRedshift          : FLOAT : Maximum redshift
                                         Must be >= 0.0, default = None: self.__maxRedshift will be used
        p_MaxRedshiftDetection : FLOAT : Maximum detection redshift
                                         Must be >= 0.0 and <= p_MaxRedshift), default = None: self.__maxRedshiftDetection will be used
        p_RedshiftStep         : FLOAT : Redshift step size
                                         Must be positive, default = None: self.__redshiftStep will be used

    Returns:
        1 : STRING : Error string - set if an error occurred.
                     'None' if no error occurred.
    """
    def __CalculateRedshiftRelatedParams(self,
                                         p_MaxRedshift          = None, 
                                         p_MaxRedshiftDetection = None, 
                                         p_RedshiftStep         = None):

        func = self.__className + '::__CalculateRedshiftRelatedParams(): '

        # default return value
        errStr = None

        # set local variables from parameter values
        # Doing it this way ensures that we have reasonable values for all parameters and allows users to call the function
        # with positional parameters and just set any parameters they want to default to 'None' rather than specify the value
        maxRedshift          = self.__maxRedshift          if p_MaxRedshift          is None else p_MaxRedshift
        maxRedshiftDetection = self.__maxRedshiftDetection if p_MaxRedshiftDetection is None else p_MaxRedshiftDetection
        redshiftStep         = self.__redshiftStep         if p_RedshiftStep         is None else p_RedshiftStep

        # check parameter types and values
        tmpStr = ' must be a floating point number'                                                                 # common error string

        if   not (isinstance(maxRedshift, float) and maxRedshift >= 0.0):  errStr = func + 'Maximum formation redshift' + tmpStr + ' >= 0.0'
        elif not (isinstance(maxRedshiftDetection, float) and \
                  maxRedshiftDetection >= 0.0             and \
                  maxRedshiftDetection <= maxRedshift):                    errStr = func + 'Maximum detection redshift' + tmpStr + ' >= 0.0 and <= maximum formation redshift ({:.6f})'.format(maxRedshift)
        elif not (isinstance(redshiftStep, float) and redshiftStep > 0.0): errStr = func + 'Redshift step size' + tmpStr
        else:                                                                                                       # all parameters ok

            # create a list of redshifts and record lengths
            self.__redshifts           = np.arange(0.0, maxRedshift + redshiftStep, redshiftStep)
            self.__nRedshifts          = self.__redshifts.shape[0]
            self.__nRedshiftsDetection = int(maxRedshiftDetection / redshiftStep)

            # convert redshifts to times and ensure all times are in Myr
            self.__times  = cosmology.age(self.__redshifts).to(units.Myr).value
            self.__nTimes = self.__times.shape[0]

            # convert redshifts to distances and ensure all distances are in Mpc (also avoid D=0 because division by 0)
            self.__distances    = cosmology.luminosity_distance(self.__redshifts).to(units.Mpc).value
            self.__distances[0] = 0.001

            # convert redshifts to volumes and ensure all volumes are in Gpc^3
            volumes = cosmology.comoving_volume(self.__redshifts).to(units.Gpc**3).value

            # split volumes into shells and duplicate last shell to keep same length
            self.__shellVolumes = np.diff(volumes)
            self.__shellVolumes = np.append(self.__shellVolumes, self.__shellVolumes[-1])

            # update class member variables to match parameters
            self.__maxRedshift          = maxRedshift
            self.__maxRedshiftDetection = maxRedshiftDetection
            self.__redshiftStep         = redshiftStep

        return errStr


    """
    CosmicIntegration::__CalculateStarFormingMass()

    Calculate the star forming mass per unit volume per year using Neijssel+19 Eq. 6

    Updates the following class member variables to parameter values if no error occurs:
        __redshifts
        __sfr_a
        __sfr_b
        __sfr_c
        __sfr_d

    Args:
        p_Redshifts : 1D ARRAY of FLOAT : Redshifts at which star forming mass should be calculated
                                          Must be >= 0.0, default = None: self.__redshifts will be used
        p_SFR_a     : FLOAT             : Scale factor for SFR (see Neijssel+2019)
                                          Default = None: self.__sfr_a will be used
        p_SFR_b     : FLOAT             : Scale factor for SFR (see Neijssel+2019)
                                          Default = None: self.__sfr_b will be used
        p_SFR_c     : FLOAT             : Scale factor for SFR (see Neijssel+2019)
                                          Default = None: self.__sfr_c will be used
        p_SFR_d     : FLOAT             : Scale factor for SFR (see Neijssel+2019)
                                          Default = None: self.__sfr_d will be used

    Returns:
        1 : 1D ARRAY of FLOAT : Star forming mass per unit volume per year for each redshift.
                                'None' if an error occurred.
        2 : STRING            : Error string - set if an error occurred.
                                'None' if no error occurred.          
    """
    def __CalculateStarFormingMass(self,
                                   p_Redshifts = None, 
                                   p_SFR_a     = None, 
                                   p_SFR_b     = None, 
                                   p_SFR_c     = None, 
                                   p_SFR_d     = None):

        func = self.__className + '::__CalculateStarFormingMass(): '

        # default return values
        sfr    = None
        errStr = None

        # set local variables from parameter values
        # Doing it this way ensures that we have reasonable values for all parameters and allows users to call the function
        # with positional parameters and just set any parameters they want to default to 'None' rather than specify the value
        redshifts = self.__redshifts if p_Redshifts is None else p_Redshifts
        sfr_a     = self.__sfr_a     if p_SFR_a     is None else p_SFR_a
        sfr_b     = self.__sfr_b     if p_SFR_b     is None else p_SFR_b
        sfr_c     = self.__sfr_c     if p_SFR_c     is None else p_SFR_c
        sfr_d     = self.__sfr_d     if p_SFR_d     is None else p_SFR_d

        # check parameter types and values
        tmpStr = ' must be a floating point number'                                                                 # common error string

        if redshifts is None:                                                 errStr = func + 'No redshifts specified'
        elif not (isinstance(redshifts, np.ndarray) and redshifts.ndim == 1): errStr = func + 'Redshifts must be specified as a 1D array'
        elif redshifts.shape[0] < 1:                                          errStr = func + 'No redshifts specified'
        elif not (redshifts.dtype == 'float64' and np.all(redshifts >= 0.0)): errStr = func + 'Redshifts must be floating point numbers >= 0.0'
        elif not isinstance(sfr_a, float):                                    errStr = func + 'SFR scale factor SFR_a' + tmpStr
        elif not isinstance(sfr_b, float):                                    errStr = func + 'SFR scale factor SFR_b' + tmpStr
        elif not isinstance(sfr_c, float):                                    errStr = func + 'SFR scale factor SFR_c' + tmpStr
        elif not isinstance(sfr_d, float):                                    errStr = func + 'SFR scale factor SFR_d' + tmpStr
        else:                                                                                                       # all parameters ok

            # use Neijssel+19 to get value in mass per year per cubic Mpc and convert to per cubic Gpc then return
            sfr = sfr_a * ((1.0 + redshifts)**sfr_b) / (1.0 + ((1.0 + redshifts) / sfr_c)**sfr_d) * units.Msun / units.yr / units.Mpc**3
            sfr = sfr.to(units.Msun / units.yr / units.Gpc**3).value

            # update class member variables to match parameters
            self.__redshifts = redshifts
            self.__sfr_a     = sfr_a
            self.__sfr_b     = sfr_b
            self.__sfr_c     = sfr_c
            self.__sfr_d     = sfr_d

        return sfr, errStr


    """
    CosmicIntegration::__CalculateZdistribution()

    Calculate the distribution of metallicities at different redshifts using a log-skew-normal distribution.
    The log-normal distribution is a special case of the log-skew-normal distribution distribution, and is
    retrieved by setting the skewness to zero (p_ALPHA = 0).

    Based on the method in Neijssel+19.  Setting:

        p_MU_0 = 0.035, p_MU_Z = -0.23, p_SIGMA_0 = 0.39, p_SIGMA_Z = 0.0, p_ALPHA =0.0

    will retrieve the dP/dZ distribution used in Neijssel+19.  See van Son+2022 for skewed-log-normal distribution.

    This function assumes that metallicities in COMPAS are drawn from a flat-in-log distribution.

    Updates the following class member variables to parameter values if no error occurs:
        __alpha
        __maxLogZ
        __minLogZ
        __mu_0
        __mu_Z
        __redshifts
        __sigma_0
        __sigma_Z
        __stepLogZ

    Args:
        p_Redshifts : 1D ARRAY of FLOAT : Redshifts at which various values should be calculated
                                          Must be >= 0.0, default = None: self.__redshifts will be used
        p_MinLogZ   : FLOAT             : Minimum logZ for dPdlogZ calculation (influences normalization)
                                          Must be >= CI_MinLogZ and <= CI_MaxLogZ, default = None: self__.minLogZ will be used
        p_MaxLogZ   : FLOAT             : Maximum logZ for dPdlogZ calculation (influences normalization)
                                          Must be > p_MinLogZ and <= CI_MaxLogZ, default = None: self__.maxLogZ will be used
        p_StepLogZ  : FLOAT             : Size of logZ steps to take in finding a Z range
                                          Must be positive, default = default = None: self__.stepLogZ will be used
        p_MU_0      : FLOAT             : Location (mean in normal) at redshift 0
                                          Must be >= 0.0 and <= 1.0, default = None: self__.mu_0 will be used
        p_MU_Z      : FLOAT             : Redshift scaling/evolution of the mean metallicity
                                          Default = None: self__.mu_Z will be used
        p_SIGMA_0   : FLOAT             : Scale (variance in normal) at redshift 0
                                          Must be >= 0.0, default = None: self__.sigma_0 will be used
        p_SIGMA_Z   : FLOAT             : Redshift scaling of the scale (variance in normal)
                                          Default = None: self__.sigma_Z will be used
        p_ALPHA     : FLOAT             : Shape (skewness) of metallicity density distribution, p_ALPHA = 0 retrieves log-normal dist)
                                          Default = None: self__.alpha will be used

    Returns:
        1 : 2D ARRAY of FLOAT : Probability of getting a particular logZ at a certain redshift
                                'None' if an error occurred.
        2 : LIST of FLOAT     : Metallicities at which dPdlogZ is evaluated
                                'None' if an error occurred.
        3 : FLOAT             : Probability of drawing a certain metallicity in COMPAS (float because assuming uniform)
                                'None' if an error occurred.
        4 : STRING            : Error string - set if an error occurred.
                                'None' if no error occurred.          
    """   
    def __CalculateZdistribution(self,
                                 p_Redshifts = None,
                                 p_MinLogZ   = None,
                                 p_MaxLogZ   = None,
                                 p_StepLogZ  = None,
                                 p_MU_0      = None, 
                                 p_MU_Z      = None, 
                                 p_SIGMA_0   = None,
                                 p_SIGMA_Z   = None,
                                 p_ALPHA     = None): 

        func = self.__className + '::__CalculateZdistribution(): '
 
        # default return values
        dPdlogZ = None
        Zs      = None
        pDrawZ  = None
        errStr  = None

        # set local variables from parameter values
        # Doing it this way ensures that we have reasonable values for all parameters and allows users to call the function
        # with positional parameters and just set any parameters they want to default to 'None' rather than specify the value
        redshifts = self.__redshifts if p_Redshifts is None else p_Redshifts
        minLogZ   = self.__minLogZ   if p_MinLogZ   is None else p_MinLogZ
        maxLogZ   = self.__maxLogZ   if p_MaxLogZ   is None else p_MaxLogZ
        stepLogZ  = self.__stepLogZ  if p_StepLogZ  is None else p_StepLogZ
        mu_0      = self.__mu_0      if p_MU_0      is None else p_MU_0 
        mu_Z      = self.__mu_Z      if p_MU_Z      is None else p_MU_Z 
        sigma_0   = self.__sigma_0   if p_SIGMA_0   is None else p_SIGMA_0
        sigma_Z   = self.__sigma_Z   if p_SIGMA_Z   is None else p_SIGMA_Z
        alpha     = self.__alpha     if p_ALPHA     is None else p_ALPHA 

        # check parameter types and values
        tmpStr = ' must be a floating point number'                                                                 # common error string

        if   redshifts is None:                                               errStr = func + 'No redshifts specified'
        elif not (isinstance(redshifts, np.ndarray) and redshifts.ndim == 1): errStr = func + 'Redshifts must be specified as a 1D array'
        elif redshifts.shape[0] < 1:                                          errStr = func + 'No redshifts specified'
        elif not (redshifts.dtype == 'float64' and np.all(redshifts >= 0.0)): errStr = func + 'Redshifts must be floating point numbers >= 0.0'
        elif not (isinstance(minLogZ, float) and \
                 minLogZ >= CI_MinLogZ and minLogZ <= CI_MaxLogZ):            errStr = func + 'Minimum logZ for dPdlogZ calculation' + tmpStr + ' >= {:.6f} and <= {:.6f}'.format(CI_MinLogZ, CI_MaxLogZ)
        elif not (isinstance(maxLogZ, float) and \
                 maxLogZ > minLogZ and maxLogZ <= CI_MaxLogZ):                errStr = func + 'Maximum logZ for dPdlogZ calculation' + tmpStr + ' > {:.6f} and <= {:.6f}'.format(minLogZ, CI_MaxLogZ)
        elif not (isinstance(stepLogZ, float) and stepLogZ > 0.0):            errStr = func + 'logZ step size' + tmpStr + ' > 0.0'
        elif not (isinstance(mu_0, float) and mu_0 >= 0.0 and mu_0 <= 1.0):   errStr = func + 'Redshift mean at z=0 (mu_0)' + tmpStr + ' >= 0.0 and <= 1.0'
        elif not isinstance(mu_Z, float):                                     errStr = func + 'Redshift mean scale factor (mu_z)' + tmpStr
        elif not (isinstance(sigma_0, float) and sigma_0 >= 0.0):             errStr = func + 'Redshift variance at z=0 (sigma_0)' + tmpStr + ' >= 0.0'
        elif not isinstance(sigma_Z, float):                                  errStr = func + 'Redshift variance scale factor (sigma_Z)' + tmpStr
        elif not isinstance(alpha, float):                                    errStr = func + 'Distribution skew (alpha)' + tmpStr
        else:                                                                                                       # all parameters ok

            sigma = sigma_0 * 10**(sigma_Z * redshifts)                                                             # Log-Linear redshift dependence of sigma
            meanZ = mu_0 * 10**(mu_Z * redshifts)                                                                   # Follow Langer & Norman 2007? in assuming that mean metallicities evolve in z
        
            # Now we re-write the expected value of our log-skew-normal to retrieve mu
            beta  = alpha / (np.sqrt(1.0 + (alpha * alpha)))
            phi   = NormDist.cdf(beta * sigma) 
            muZ   = np.log(meanZ / 2.0 * 1.0 / (np.exp(0.5 * (sigma * sigma)) * phi)) 

            # create a range of metallicities (the x-values, or random variables)
            logZ  = np.arange(minLogZ, maxLogZ + stepLogZ, stepLogZ)
            Zs    = np.exp(logZ)

            # probabilities of log-skew-normal (without the factor of 1/Z since this is dp/dlogZ not dp/dZ)
            normPDF = NormDist.pdf((logZ - muZ[:, np.newaxis]) / sigma[:, np.newaxis])
            normCDF = NormDist.cdf(alpha * (logZ - muZ[:,np.newaxis]) / sigma[:, np.newaxis])
            dPdlogZ = 2.0 / (sigma[:, np.newaxis]) * normPDF * normCDF

            # normalise the distribution over all metallicities
            norm    = dPdlogZ.sum(axis = -1) * stepLogZ
            dPdlogZ = dPdlogZ / norm[:, np.newaxis]

            # assume a flat in log distribution in metallicity to find probability of drawing Z in COMPAS
            pDrawZ = 1.0 / (self.__COMPAS.MaxLogZ() - self.__COMPAS.MinLogZ())

            # update class member variables to match parameters
    
            self.__redshifts = redshifts
            self.__minLogZ   = minLogZ
            self.__maxLogZ   = maxLogZ
            self.__stepLogZ  = stepLogZ
            self.__mu_0      = mu_0 
            self.__mu_Z      = mu_Z 
            self.__sigma_0   = sigma_0
            self.__sigma_Z   = sigma_Z
            self.__alpha     = alpha 

        return dPdlogZ, Zs, pDrawZ, errStr


    """
    CosmicIntegration::__GetChirpMassBins()

    Creates chirp mass bins for binning of rates if necessary.

    The caller-defined bins are passed in - we may need to add bins to the extremities of the caller-defined
    bins depending upon what the caller passed in - we want the bins to encompass the minimum and maximum chirp
    masses of all DCOs.

    We also want to record bin counts - number of binaries placed in each bin - so we will add an element to
    the bins list we return.
                    
    The new (returned) bins list has three elements:

        0: list of bin widths (Msun)
        1: list of chirp masses at bin right edges (Msun)
        2: list of bin counts

    Bins in the new bins list will encompass the minimum and maximum chirp masses of all DCOs. The caller-defined
    bins are preserved: we might have just added a lower and/or upper bin to the caller-defined list of bins.

    Note that if upper and/or lower bins are added to the caller-defined list of bins, the widths of those new
    upper and/or lower bins may not be the same, either as each other, or the same as and fixed bin widths the
    caller specified.
    """
    def __GetChirpMassBins(self, p_McBins = None):

        func = self.__className + '::__GetChirpMassBins(): '

        # default return values
        newMcBins = None
        errStr    = None

        if p_McBins is not None:                                                                                    # binning? (should be, but just in case)
                      
            errStr = Utils.CheckBins(p_McBins)                                                                      # check caller-defined bins for validity
            if errStr is None:                                                                                      # ok?
                                                                                                                    # yes
                newMcBins = [p_McBins[0], p_McBins[1], []]                                                          # new bins list

                minMc = np.min(self.__COMPAS.ChirpMasses()[0])                                                      # minimum chirp mass of all DCOs
                maxMc = np.max(self.__COMPAS.ChirpMasses()[0])                                                      # maximum chirp mass of all DCOs

                firstBinLeftEdge = newMcBins[1][0] - (newMcBins[0][0] / 2.0)                                        # first bin left edge
                if firstBinLeftEdge >= minMc:                                                                       # first bin left edge includes minMc?
                                                                                                                    # no - insert new first bin
                    newMcBins[0].insert(0, firstBinLeftEdge - minMc)                                                # insert width of new first bin
                    newMcBins[1].insert(0, firstBinLeftEdge)                                                        # insert right edge of new first bin

                lastBinRightEdge = newMcBins[1][len(newMcBins[1]) - 1]                                              # last bin right edge
                if lastBinRightEdge < maxMc:                                                                        # last bin right edge includes maxMc?
                                                                                                                    # no - append new last bin
                    newMcBins[0].append(maxMc - lastBinRightEdge)                                                   # append width of new last bin
                    newMcBins[1].append(maxMc)                                                                      # append right edge of new last bin

                newMcBins[2] = [0] * len(newMcBins[0])                                                              # initialise bin counts to 0

        return newMcBins, errStr

        
    """
    CosmicIntegration::__ClosePerBinaryFile()

    Close a per-binary .csv file - parameters determine which file.

    Args:
        p_PerBinaryFile     : FILE   : File handle of file for which close is requested.  File type must be 'io.TextIOWrapper'
                                       Default = None (an error)
        p_PerBinaryFilename : STRING : The name of the file for which close is requested
                                       Default = None (an error)
    Returns:
        1 : FILE   : File handle of file for which close was requested.
                     'None' if closed ok.
        4 : STRING : Error string - set if an error occurred.
                     'None' if no error occurred.
    """
    def __ClosePerBinaryFile(self, p_PerBinaryFile = None, p_PerBinaryFilename = None):

        func = self.__className + '::__ClosePerBinaryFile(): '

        # default return values
        thisFile = None
        errStr   = None

        # check parameters
        if   p_PerBinaryFile is None or not isinstance(p_PerBinaryFile, io.TextIOWrapper):
            errStr = func + 'Per-binary file must be of type \'io.TextIOWrapper\''
        elif p_PerBinaryFilename is None or not (isinstance(p_PerBinaryFilename, str) and len(p_PerBinaryFilename) > 0):
            errStr = func + 'Per-binary filename must be a non-empty string'
        else:                                                                                                       # all parameters ok
            thisFile = p_PerBinaryFile
            try:
                thisFile.close()                                                                                    # close file
                thisFile = None
            except OSError as OSerror :                                                                             # close failed
                errStr = func + 'Could not close per-binary file \'' + p_PerBinaryFilename + '\'\n' + OSerror

        return thisFile, errStr


    """
    CosmicIntegration::__WritePerBinaryRecord()

    Write a record to one of the per-binary .csv files

    Prints an error message if write fails

    Args:
        p_PerBinaryFile     : FILE   : File handle of file to which record is to be written.  File type must be 'io.TextIOWrapper'
                                       Default is None
        p_PerBinaryFilename : STRING : The name of the file to which the record is to be written
                                       Default is None
        p_Record            : STRING : The record to be written to the file
                                       Default is None
    Returns:
        1 : STRING : Error string - set if an error occurred.
                     'None' if no error occurred.
    """
    def __WritePerBinaryRecord(self, p_PerBinaryFile = None, p_PerBinaryFilename = None, p_Record = None):

        func = self.__className + '::__WritePerBinaryRecord(): '

        # default return value
        errStr = None

        # check parameters
        if p_PerBinaryFile is None:                                                       errStr = func + 'No per-binary file specified'
        elif not isinstance(p_PerBinaryFile, io.TextIOWrapper):                           errStr = func + 'Per-binary file must be of type \'io.TextIOWrapper\''
        elif p_PerBinaryFilename is None:                                                 errStr = func + 'No per-binary filename specified'
        elif not (isinstance(p_PerBinaryFilename, str) and len(p_PerBinaryFilename) > 0): errStr = func + 'Per-binary filename must be a non-empty string'
        elif p_Record is None:                                                            errStr = func + 'No record to be written specified'
        elif not (isinstance(p_Record, str) and len(p_Record) > 0):                       errStr = func + 'Record to be written must be a non-empty string'
        else:                                                                                                       # all parameters ok

            try:
                p_PerBinaryFile.write(p_Record + '\n')                                                              # write record
            except OSError as OSerror :                                                                             # write failed
                errStr = func + 'Could not write to per-binary file \'' + p_PerBinaryFilename + '\'\n' + OSerror

        return errStr

        
    """
    CosmicIntegration::__OpenPerBinaryFile()

    Open a per-binary .csv file - parameters determine which file - and write preamble to file.

    The filename is determined by the p_FilenameStem parameter.  No existing files will be overwritten - the file 
    will be versioned if a file already exists with the same name.  A version string will be constructed and appended
    to the p_FilenameStem parameter if required.

    See Utils::IntToLetters() for the method of construction of the version string.

    The file-specific preamble, and formatting, can be found in Ilya Mandel's "Binary mergers: a common format"
    document.

    Args:
        p_FilenameStem   : STRING : The filename stem (or prefix) to which the version will be appended (if necessary)
                                    Default = None (an error)
        p_CreatorDetails : STRING : The details of the file creator, to be written to the 'creator' meta-data row
                                    Default = None: self.__creatorDetails will be used
    Returns:
        1 : FILE   : File handle of file opened.  File type will be 'io.TextIOWrapper'
                     'None' if an error occurred.
        2 : STRING : The name of the file to which the record is to be written
                     'None' if an error occurred.
        4 : STRING : Error string - set if an error occurred.
                     'None' if no error occurred.
    """
    def __OpenPerBinaryFile(self, p_FilenameStem = None, p_CreatorDetails = None):

        func = self.__className + '::__OpenPerBinaryFile(): '

        # default return values
        thisFile     = None
        thisFilename = None
        errStr       = None

        # set local variables from parameter values (where necessary - some are checked in called functions)
        # Doing it this way ensures that we have reasonable values for all parameters and allows users to call the function
        # with positional parameters and just set any parameters they want to default to 'None' rather than specify the value
        filenameStem   = p_FilenameStem
        creatorDetails = self.__creatorDetails if p_CreatorDetails is None else p_CreatorDetails

        # check parameters
        if filenameStem is None:                                                 errStr = func + 'No filename stem specified'
        elif not (isinstance(filenameStem, str) and len(filenameStem) > 0):      errStr = func + 'Filename stem must be a non-empty string'
        elif creatorDetails is not None and not isinstance(creatorDetails, str): errStr = func + 'When specified, creator details must be a non-empty string'
        else:                                                                                                       # all parameters ok

            # open the file if possible
            versionNum    = 0                                                                                       # file version number
            versionStr    = ''                                                                                      # file version string
            thisFilename  = filenameStem + '.csv'                                                                   # fully-qualified filename
            while os.path.isfile(thisFilename) and errStr is None:                                                  # already exists?
                                                                                                                    # yes
                versionNum += 1                                                                                     # increment version number
                versionStr, errStr = Utils.IntToLetters(versionNum)                                                 # get new version string
                if errStr is None:
                    thisFilename = filenameStem + '_' + versionStr + '.csv'                                         # fully-qualified filename
                else:
                    thisFilename = None                                                                             # error return

            if errStr is None:
                try:
                    thisFile = open(thisFilename, 'xt')                                                             # create and open file

                    # file is open - write meta-data, header, and units rows

                    try:
                        thisFile.write('M,model-type,isolated-binary\n')                                            # first meta-data row
                        thisFile.write('M,code-name,COMPAS\n')                                                      # second meta-data row
                        thisFile.write('M,code-version,' + self.__COMPAS.Version() + '\n')                          # third meta-data row
                        thisFile.write('M,date,' + self.__today + '\n')                                             # fouth meta-data row

                        if creatorDetails is not None and len(creatorDetails) > 0:                                  # have creator details?
                            thisFile.write('M,creator,' + creatorDetails + '\n')                                    # fifth meta-data row

                    except OSError as OSerror :                                                                     # one of the writes failed
                        errStr = func + 'Could not write meta-data to per-binary file \'' + thisFilename + '\'\n' + OSerror

                except OSError as OSerror :                                                                         # create/open of file failed
                    errStr = func + 'Could not create per-binary file \'' + thisFilename + '\'\n' + OSError
                                                                                                                    # no
            # we could remove files if we had errors here - e.g. if we couldn't write the header row to the mergers
            # file for some reason, we could cleanup - close the file and remove it, but for now I'll leave it there
            # in case it useful as a debugging tool to find out  what went wrong.  If it becomes an issue, or too 
            # tedious to clean up the files manually in case of error, we can revist removing then here.

            if errStr is not None:
                thisFile     = None                                                                                 # flag file not open
                thisFilename = None                                                                                 # error return

        return thisFile, thisFilename, errStr


    """
    CosmicIntegration::__OpenPerBinaryFiles()

    Opens per-binary .csv files for mergers and yields.

    Fully qualified filename stems are constructed by concatenating the path passed as p_PerBinaryFilesPath with
    the defined constants for the per-binary files and the current date - a version string will be added later 
    if necessary.  The file naming convention can be found in Ilya Mandel's "Binary mergers: a common format"
    document.

    Args:
        p_PerBinaryFilesPath : STRING : Path at which the per-binary files should be created
                                        Default = None: self.__perBinaryFilesPath will be used
        p_CreatorDetails     : STRING : The details of the file creator, to be written to the 'creator' meta-data row
                                        Default = None: self.__creatorDetails will be used
    Returns:
        1 : FILE   : File handle for the open mergers file.  File type will be 'io.TextIOWrapper'
                     'None' if an error occurred.
        2 : STRING : The name of the mergers file
                     'None' if an error occurred.
        3 : FILE   : File handle for the open yields file.  File type will be 'io.TextIOWrapper'
                     'None' if an error occurred.
        4 : STRING : The name of the yields file
                     'None' if an error occurred.
        5 : STRING : Error string - set if an error occurred.
                     'None' if no error occurred.
    """
    def __OpenPerBinaryFiles(self, p_PerBinaryFilesPath = None, p_CreatorDetails = None):

        func = self.__className + '::__OpenPerBinaryFiles(): '

        # default return values
        mergersFile     = None
        mergersFilename = None
        yieldsFile      = None
        yieldsFilename  = None
        errStr          = None

        # set local variables from parameter values (where necessary - some are checked in called functions)
        # Doing it this way ensures that we have reasonable values for all parameters and allows users to call the function
        # with positional parameters and just set any parameters they want to default to 'None' rather than specify the value
        perBinaryFilesPath = self.__perBinaryFilesPath if p_PerBinaryFilesPath is None else p_PerBinaryFilesPath

        # check parameters where necessary
        if not isinstance(perBinaryFilesPath, str): errStr = func + 'When specified, per-binary files path must be a string'
        else:                                                                                                           # all parameters ok

            # fix up per-binary files path if necessary
            if   len(perBinaryFilesPath) < 1:   perBinaryFilesPath = CI_PerBinaryFilesPath                              # handle empty path
            elif perBinaryFilesPath[-1] != '/': perBinaryFilesPath += '/'                                               # add path delimiter if necessary

            # construct name stems for the per-binary output files

            stemMergers = perBinaryFilesPath + CI_Mergers_FilenamePrefix + self.__today                                 # file name stem for mergers file
            stemYields  = perBinaryFilesPath + CI_Yields_FilenamePrefix + self.__today                                  # file name stem for yields file
            stemMeta    = perBinaryFilesPath + CI_Meta_FilenamePrefix + self.__today                                    # file name stem for meta-datas file

            mergersFile, mergersFilename, errStr = self.__OpenPerBinaryFile(stemMergers, p_CreatorDetails)              # open mergers file

            # the mergers file is mandatory - if we can't create that we don't bother creating the rest
            if errStr is None and mergersFile is not None:                                                              # have mergers file?
                                                                                                                        # yes - write header and units rows
                errStr = self.__WritePerBinaryRecord(mergersFile, mergersFilename, 'H,Mass1,Mass2,Redshift,MergerRate')
                if errStr is None:
                    errStr = self.__WritePerBinaryRecord(mergersFile, mergersFilename, 'U,Msun,Msun,,Gpc-3yr-1')
                    if errStr is None:                                                                                  # ok?
                                                                                                                        # yes
                        # open yields file
                        # as much as we'd like to produce the yields file, it is not mandatory - if we can't
                        # create it here and write headers etc., we won't return an error - we'll just issue
                        # a warning here and continue.
                        yieldsFile, yieldsFilename, tmpErr = self.__OpenPerBinaryFile(stemYields, p_CreatorDetails)     # open yields file
                        if tmpErr is None and yieldsFile is not None:                                                   # have yields file?
                                                                                                                        # yes - write header and units rows
                            tmpErr = self.__WritePerBinaryRecord(yieldsFile, yieldsFilename, 'H,Mass1,Mass2,DelayTime,Metallicity,Yield')
                            if tmpErr is None:
                                tmpErr = self.__WritePerBinaryRecord(yieldsFile, yieldsFilename, 'U,Msun,Msun,Myr,,MMsun-1')

                        if tmpErr is not None:                                                                          # yields file ok?
                            print('\n!! WARNING !!')                                                                    # no - show warning
                            print(tmpErr)

        if errStr is not None:                                                                                          # ok?
                                                                                                                        # no
            # we could remove files if we had errors here - e.g. if we couldn't write the header row to the mergers
            # file for some reason, we could cleanup - close the file and remove it, but for now I'll leave it there
            # in case it useful as a debugging tool to find out  what went wrong.  If it becomes an issue, or too 
            # tedious to clean up the files manually in case of error, we can revist removing then here.

            mergersFile     = None                                                                                      # error return
            mergersFilename = None                                                                                      # error return
            yieldsFile      = None                                                                                      # error return
            yieldsFilename  = None                                                                                      # error return

        return mergersFile, mergersFilename, yieldsFile, yieldsFilename, errStr


    # Public interface


    """
    CosmicIntegration::FindRates()

    Find the formation, merger, and detection rates for each binary at each redshift

    Note: this is a long function, and could easily indent off the right-side of the page...  To avoid that,
    and (hopefully) to make the function a little more readable, every now and then I bring the indentation
    back a bit to the left.  Because an error condition is always flagged by setting the error string 'errStr',
    every now and then I have a checkpoint, denoted by the text:

        '# checkpoint: no errors - ok to proceed?'

    followed by the conditional statement:

        'if errStr is None:'

    That construct lets me effectively end a block of code, bring the indentation back a bit to the left, and
    only execute the rest of the function if everything is ok and there have been no errors at that point.  I
    think it keeps the code a little more readable for what is a very long function.  I could break it out a
    bit and farm-out some functionality to other functions, but really not that much more than I've already
    done, and as noted in the code below, for performance reasons I wanted to do as much as possible inside the
    per-DCO loop as inline code rather than add fucntion call overhead.

    Args:
    
        Class SelectionEffects related parameters: see class SelectionEffects constructor
        Class COMPAS related parameters: see class COMPAS constructor

        Class CosmicIntegration related parameters:

        p_MaxRedshift          : FLOAT         : Maximum formation redshift
                                                 Must be >= 0.0, default = None: self.__maxRedshift will be used
        p_MaxRedshiftDetection : FLOAT         : Maximum redshift to calculate detection rates
                                                 Must be >= 0.0 and <= p_MaxRedshift, default = None: self.__maxRedshiftDetection will be used
        p_RedshiftStep         : FLOAT         : Size of step to take in redshift
                                                 Must be > 0.0, default = None: self.__redshiftStep will be used
        p_MinLogZ              : FLOAT         : Minimum logZ for dPdlogZ calculation (influences normalization)
                                                 Must be >= CI_MinLogZ and <= CI_MaxLogZ, default = None: self__.minLogZ will be used
        p_MaxLogZ              : FLOAT         : Maximum logZ for dPdlogZ calculation (influences normalization)
                                                 Must be > p_MinLogZ and <= CI_MaxLogZ, default = None: self__.maxLogZ will be used
        p_StepLogZ             : FLOAT         : Size of logZ steps to take in finding a Z range
                                                 Must be positive, default = default = None: self__.stepLogZ will be used
        p_MU_0                 : FLOAT         : Location (mean in normal) at redshift 0
                                                 Must be >= 0.0 and <= 1.0, default = None: self__.mu_0 will be used
        p_MU_Z                 : FLOAT         : Redshift scaling/evolution of the mean metallicity
                                                 Default = None: self__.mu_Z will be used
        p_SIGMA_0              : FLOAT         : Scale (variance in normal) at redshift 0
                                                 Must be >= 0.0, default = None: self__.sigma_0 will be used
        p_SIGMA_Z              : FLOAT         : Redshift scaling of the scale (variance in normal)
                                                 Default = None: self__.sigma_Z will be used
        p_ALPHA                : FLOAT         : Shape (skewness) of metallicity density distribution, p_ALPHA = 0 retrieves log-normal dist)
                                                 Default = None: self__.alpha will be used
        p_SFR_a                : FLOAT         : Scale factor for SFR (see Neijssel+2019)
                                                 Default = None: self.__sfr_a will be used
        p_SFR_b                : FLOAT         : Scale factor for SFR (see Neijssel+2019)
                                                 Default = None: self.__sfr_b will be used
        p_SFR_c                : FLOAT         : Scale factor for SFR (see Neijssel+2019)
                                                 Default = None: self.__sfr_c will be used
        p_SFR_d                : FLOAT         : Scale factor for SFR (see Neijssel+2019)
                                                 Default = None: self.__sfr_d will be used
        p_McBins               : LIST of LIST  : Chirp mass bins to be used for binning rates - fixed or variable width
                                                 Default = None: self.__McBins will be used

                                                 A final value of 'None', or an empty list ([]), indicates no binning - FindRates() will be calculated per binary.
                                                 'final value' here means after assignment to self.__McBins if 'None' was passed as the parameter here - this allows
                                                 the caller to either use previously specified bins, or turn of binning if bins had been previously specified.

                                                 If not 'None', or an empty LIST, p_McBins must be a LIST and must contain (in order):

                                                     binWidths     : LIST of FLOAT : Width of each chirp mas bin (Msun)
                                                     binRightEdges : LIST of FLOAT : Chirp mass at right (upper) edge of each chirp mass bin

                                                     (see class UTILS for helper functions to create and check chirp mass bins)

        p_PerBinaryFiles       : BOOL          : Flag to specify whether per-binary files (mergers and yields) are required
                                                 Default = None: self.__perBinaryFiles will be used
        p_PerBinaryFilesPath   : STRING        : Path at which per-binary files should be created (if necessary)
                                                 Default = None: self.__perBinaryFilesPath will be used
        p_CreatorDetails       : STRING        : Creator details for per-binary file meta-data (if necessary)
                                                 (usually name and email address - can contain commas for CSV file)
                                                 Default = None: no creator details will be recorded in per-binary files       
        p_RandomSeed           : INT           : Random seed to initialize numpy randum number generator
                                                 Must be >= 0, default = None: use numpy default seed (note: means results will (very likely) not be reproducible)
        p_Verbose              : BOOL          : Flag to indicate if diagnostics/stats should be printed to console
                                                 Default = VERBOSE_DEFAULT


        p_CreatorDetails       : STRING        : Creator details for per-binary file meta-data (if necessary)
                                                 (usually name and email address - can contain commas for CSV file)
                                                 Default = None
                                                 A value of 'None', or an empty string (''), indicates that no creator details should be recorded in per-binary files       

    Returns:
        1 : LIST : LIST containing the following LISTs:

            1 : LIST of numpy ARRAYs : LIST containg the following numpy ARRAYs:

                1 : 2D ARRAY of FLOAT : Formation rate for each binary at each redshift
                                'None' if an error occurred.
                2 : 2D ARRAY of FLOAT : Merger rate for each binary at each redshift
                                'None' if an error occurred.
                3 : 2D ARRAY of FLOAT : Detection rate for each binary at each detection redshift
                                'None' if an error occurred.

            2 : LIST or None: Chirp mass bins used for binning rates.  
                              'None' if no binning was performed.
                               If not 'None', will be a LIST and will contain:

                                    1 : LIST of FLOAT : Width of each chirp mas bin (Msun)
                                    2 : LIST of FLOAT : Chirp mass at right (upper) edge of each chirp mass bin

        2 : LIST : List containing the following LISTs:

            1 : LIST : List containing the following:
                1 : INT    : Count of data records written to per-binary mergers file
                             (does not include comment, meta-data, header, or units records)
                             -1 if an error occurred
                2 : STRING : Fully qualified name of per-binary mergers file
                             (will be valid if any data written to file, even in an error occurred)

            2 : LIST : List containing the following:

                1 : INT    : Count of data records written to per-binary yields file
                             (does not include comment, meta-data, header, or units records)
                             -1 if an error occurred
                2 : STRING : Fully qualified name of per-binary yields file
                             (will be valid if any data written to file, even in an error occurred)

        3 : STRING : Error string - set if an error occurred.
                     'None' if no error occurred.          
    """
    def FindRates(self,

                  # CosmicIntegration related parameters
                  p_MaxRedshift          = None,
                  p_MaxRedshiftDetection = None,
                  p_RedshiftStep         = None,
                  p_MinLogZ              = None,
                  p_MaxLogZ              = None,
                  p_StepLogZ             = None,
                  p_MU_0                 = None,
                  p_MU_Z                 = None,
                  p_SIGMA_0              = None,
                  p_SIGMA_Z              = None,
                  p_ALPHA                = None,
                  p_SFR_a                = None,
                  p_SFR_b                = None,
                  p_SFR_c                = None,
                  p_SFR_d                = None,
                  p_McBins               = None,
                  p_PerBinaryFiles       = None,
                  p_PerBinaryFilesPath   = None,
                  p_CreatorDetails       = None,

                  # SelectionEffects related parameters
                  p_Sensitivity          = None,
                  p_SNRthreshold         = None,
                  p_McMax                = None,
                  p_McStep               = None,
                  p_ETAmax               = None,
                  p_ETAstep              = None,
                  p_SNRmax               = None,
                  p_SNRstep              = None,

                  # COMPAS related parameters
                  p_m1Minimum            = None,
                  p_m1Maximum            = None,
                  p_m2Minimum            = None,
                  p_UseSampledMassRanges = None,
                  p_BinaryFraction       = None,    
                  p_DCOtype              = None,
                  p_WithinHubbleTime     = None,
                  p_PessimisticCE        = None,
                  p_NoRLOFafterCEE       = None,

                  # miscellaneos parameters
                  p_Verbose              = None):

        func = self.__className + '::Findrates(): '

        self.__verbose = VERBOSE_DEFAULT if p_Verbose is None else p_Verbose

        # default return values
        formationRate = None
        mergerRate    = None
        detectionRate = None

        if self.__verbose:
            memUsage  = psutil.Process().memory_info().rss / GB
            print('...Beginning CosmicIntegration::FindRates().  Process total memory usage = {:.6f} GB\n'.format(memUsage))

            mark = time.time()

        # save state - this is cleaner here than trying not to change state if an error occurs
        self.__SaveState()

        # set local variables from parameter values (where necessary - some are checked in called functions)
        # Doing it this way ensures that we have reasonable values for all parameters and allows users to call the function
        # with positional parameters and just set any parameters they want to default to 'None' rather than specify the value
        McStep         = self.__SE.McStep()    if p_McStep         is None else p_McStep
        ETAstep        = self.__SE.ETAstep()   if p_ETAstep        is None else p_ETAstep
        SNRstep        = self.__SE.SNRstep()   if p_SNRstep        is None else p_SNRstep
        McBins         = self.__McBins         if p_McBins         is None else (None if (isinstance(p_McBins, list) and len(p_McBins) == 0) else p_McBins)
        perBinaryFiles = self.__perBinaryFiles if p_PerBinaryFiles is None else p_PerBinaryFiles

        # check parameter types and values
        # only some here - most are checked in called functions
        tmpStr = ' must be a floating point number'                                                                     # common error string

        if   not (isinstance(McStep,  float) and McStep  > 0.0): errStr = func + 'Chirp mass step size' + tmpStr + ' > 0.0'
        elif not (isinstance(ETAstep, float) and ETAstep > 0.0): errStr = func + 'Symmetric mass ratio step size' + tmpStr + ' > 0.0'
        elif not (isinstance(SNRstep, float) and SNRstep > 0.0): errStr = func + 'SNR step size' + tmpStr + ' > 0.0'
        elif perBinaryFiles is not None and \
             not isinstance(perBinaryFiles, bool):               errStr = func + 'When specified, flag for per-binary files must be a boolean'
        else:                                                                                                           # all checked parameters ok

            # print preamble if necessary
            if self.__verbose:
                print('   Using data created by COMPAS v' + self.__COMPAS.Version() + ' on ' + self.__COMPAS.RunStart())
                print('   Number of binaries     :', self.__COMPAS.nBinaries())
                print('   Number of DCOs         :', self.__COMPAS.nAllDCOs())
                print('   DCO type requested     :', self.__COMPAS.DCOtype())
                print('   Number of matching DCOs:', self.__COMPAS.nDCOs(), '\n')


            # reinitialise SelectionEffects object with parameters specified
            # SelectionEffects will return an error if parms are bad
            errStr = self.__SE.Initialise(p_Sensitivity, p_SNRthreshold, p_McMax, p_McStep, p_ETAmax, p_ETAstep, p_SNRmax, p_SNRstep, p_Verbose)                
            if errStr is None:

                # reinitialise COMPAS object with parameters specified
                # COMPAS will return an error if parms are bad
                errStr = self.__COMPAS.Initialise(p_m1Minimum, p_m1Maximum, p_m2Minimum, p_UseSampledMassRanges, p_BinaryFraction, p_DCOtype, p_WithinHubbleTime, p_PessimisticCE, p_NoRLOFafterCEE, p_Verbose)


        # checkpoint: no errors - ok to proceed?
        if errStr is None:
            # yes - have both SelectionEffects and COMPAS class instances
            # get some initial COMPAS data

            # set the redshifts array and its equivalents
            # will return an error if parms are bad
            errStr = self.__CalculateRedshiftRelatedParams(p_MaxRedshift, p_MaxRedshiftDetection, p_RedshiftStep)
            if errStr is None:

                # get the star forming mass per unit volume per year
                # will return an error if parms are bad
                sfr, errStr = self.__CalculateStarFormingMass(None, p_SFR_a, p_SFR_b, p_SFR_c, p_SFR_d)
                if errStr is None:

                    # calculate representative star formation mass
                    SFmass = self.__COMPAS.MassEvolvedPerBinary() * float(self.__COMPAS.nBinaries())

                    # calculate the number of binaries formed per year per cubic Gpc
                    nFormed = sfr / SFmass

                    # calculate the metallicity distribution at each redshift, and probability of drawing each metallicity in COMPAS
                    Zs      = None
                    dPdlogZ = 1
                    pDrawZ  = 1
                    if self.__COMPAS.MinLogZ() != self.__COMPAS.MaxLogZ():                                              # metallicity varies between binaries?
                                                                                                                        # yes - perform integral over metallicities   
                        # will return an error if parms are bad
                        dPdlogZ, Zs, pDrawZ, errStr = self.__CalculateZdistribution(None, p_MinLogZ, p_MaxLogZ, p_StepLogZ, p_MU_0, p_MU_Z, p_SIGMA_0, p_SIGMA_Z, p_ALPHA)
                        if errStr is None:

                            # have initial COMPAS data
                            # check for binning - create chirp mass bins if required

                            nDCOs = self.__COMPAS.nDCOs()                                                               # number of DCOs in COMPAS data
                            if McBins is None:                                                                          # binning?
                                McBins = None                                                                           # no - not binning
                                nBins  = nDCOs                                                                          # number of bins is number of DCOs
                            else:                                                                                       # yes - binning
                                McBins, errStr = self.__GetChirpMassBins(McBins)                                        # set up chirp mass bins
                                if errStr is None:
                                    nBins = len(McBins[0])                                                              # number of bins

                    else:
                        Zs      = None                                                                                  # no - single metallicity for all binaries
                        dPdlogZ = 1.0
                        pDrawZ  = 1.0



        # checkpoint: no errors - ok to proceed?
        if errStr is None:
            
            # Check for sufficient available physical memory.
            #    
            # Each binary requires 8 bytes for each entry in each of the formationRate, mergerRate, and detectionRate arrays.
            #
            # The number of entries in each of the formationRate and mergerRate arrays is nBins * nRedshifts, where nRedshifts
            # is equal to the maximum redshift / redshiftStep.
            #
            # The number of entries in the detectionRate array is nBins * nRedshiftsDetection, where nRedshiftsDetection is
            # equal to the maximum detection redshift / redshiftStep.
            #
            # Check here that we have sufficient memory available to fully allocate all of the arrays.  Even if there is sufficient
            # memory available now, that doesn't guarantee that we won't run out.  numpy does not allocate the entire array, so as
            # the loop below progresses more of the array will be allocated, and if other users/processes have consumed available
            # memory before we get to allocate it, we could run out of memory - but checking first is useful, especially if the user
            # has requested an array size that exceeds installed memory
            #
            # This check checks for available physical memory.  The system may have available swap the OS can use as virtual memory,
            # but even if swap is configured on an SSD (rather than an HDD) it is still very much slower than RAM access.  We could
            # add a commandline option to allow users to choose to use swap, but for now we assume we only want to use physical memory.

            createRates  = True                                                                                         # unless we don't have sufficient memory
            memRequired  = 8.0 * nBins * ((2.0 * self.__nRedshifts) + self.__nRedshiftsDetection) / GB                  # memory required in GB
            memAvailable = psutil.virtual_memory().available / GB                                                       # physical memory available in GB - despite the function name, this returns details of physical memeory...
            if memRequired > memAvailable:                                                                              # sufficient physical memory?
                createRates = False                                                                                     # no - don't create rates arrays

                errStr = func + 'Memory required ({:.6f} GB) exceeds physical memory available ({:.6f} GB)'.format(memRequired, memAvailable)

                # if per-binary files were requested we can continue and write the per-binary data there as
                # requested, but we can't produce and return and rates arrays (all will be returned as None)
                if perBinaryFiles:                                                                                      # per-binary files requested?
                    # show a warning that rates arrays won't be created
                    print('\n!! WARNING !!')
                    print(errStr)
                    print('Per-binary data will be written to per-binary files, but no rates arrays will be created\n')

                    errStr = None                                                                                       # ok to proceed
            else:
                # sufficient memory available - initalise rates to zero
                formationRate = np.zeros(shape = (nBins, self.__nRedshifts))
                mergerRate    = np.zeros(shape = (nBins, self.__nRedshifts))
                detectionRate = np.zeros(shape = (nBins, self.__nRedshiftsDetection))


        # checkpoint: no errors - ok to proceed?
        if errStr is None:

            # yes - we know if we're binning, and we either have sufficient physical memory or we've
            # disabled creation of the rates arrays

            mergersWriteCount = -1
            yieldsWriteCount  = -1
            # open per-binary mergers and yields output files as .csv files if necessary
            mergersFile     = None
            mergersFilename = None
            yieldsFile      = None
            yieldsFilename  = None
            if perBinaryFiles:
                mergersFile, mergersFilename, yieldsFile, yieldsFilename, errStr = self.__OpenPerBinaryFiles(p_PerBinaryFilesPath, p_CreatorDetails)

            # if we don't have sufficient memory to create rates arrays, and we couldn't create at least the
            # per-binary merger file, we might as well stop now - nothing to do...
            if errStr is None and (createRates or (mergersFile is not None and yieldsFile is not None)):

                # yes - we have sufficient memory, and/or both of the per-binary output files

                # determine redshift step
                redshiftStep = self.__redshifts[1] - self.__redshifts[0] if self.__nRedshifts > 1 else 1.0

                # interpolate times and redshifts for conversion
                timesToRedshifts = interp1d(self.__times, self.__redshifts)

                # make note of the first time at which star formation occured
                ageFirstSFR = np.min(self.__times)

                # iterate over all DCOs in the COMPAS data
                # we want to do as much as we can inline in this loop - we could have tens, or even
                # hundreds, of thousands (or millions in some cases...) of DCOs, so removing overhead
                # here is a good thing for performance (there is overhead every time a function is
                # called).  Of course we need to balance that with readability and maintainability...
                
                if self.__verbose: print()

                for i in range(nDCOs):

                    # show progress if required
                    if self.__verbose:
                        if i == 0                                             or \
                          (i > 0       and i <= 1000     and i % 100 == 0)    or \
                          (i > 1000    and i <= 10000    and i % 1000 == 0)   or \
                          (i > 10000   and i <= 100000   and i % 10000 == 0)  or \
                          (i > 100000  and i <= 500000   and i % 50000 == 0)  or \
                          (i > 500000  and i <= 1000000  and i % 100000 == 0) or \
                          (i > 1000000 and i <= 10000000 and i % 200000 == 0) or \
                          (i > 10000000                  and i % 1000000 == 0): print('   CosmicIntegration::FindRates(): Processing DCO', i, 'of', nDCOs)
                        
                    # we know 'i' is not out of bounds for any of the arrays here, so no need to bounds check
                    # or check called functions for an error return

                    # if we're binning, we calculate the bin for this binary and update that
                    # if we're not binning, the bin is just the index of the individual binary
                    McBin = i if McBins is None else np.searchsorted(McBins[1], self.__COMPAS.ChirpMasses(i)[0])
                    if McBins is not None: McBins[2][McBin] += 1                                                        # increment bin count if we're binning

                    # calculate formation rate
                    thisFormationRate = np.zeros(shape = (self.__nRedshifts))

                    # if single metallicity assume all SFR happened at that fixed metallicity
                    if Zs is None :
                        thisFormationRate = nFormed
                    else:
                        # (see Neijssel+19 Section 4) - note this uses p_dPdlogZ for *closest* metallicity
                        thisFormationRate = nFormed * dPdlogZ[:, np.searchsorted(Zs, self.__COMPAS.Zsystems(i)[0], side = 'right')] / pDrawZ

                    # update formation rate array if necessary
                    if formationRate is not None: formationRate[McBin, :] += thisFormationRate
                            
                    # initialise merger rate for this DCO
                    thisMergerRate = np.zeros(shape = (self.__nRedshifts))

                    # time at which the binary formed if it merges at this redshift 
                    formationTime = self.__times - self.__COMPAS.DelayTime(i)[0]

                    # we have only calculated formation rate up to z = max(p_Redshifts), so we need to only find
                    # merger rates for formation times at z < max(redshifts)
                    # first locate the index above which the binary would have formed before z = max(redshifts)
                    firstTooEarlyIndex = np.searchsorted(-formationTime, -ageFirstSFR)

                    # include the whole array if searchsorted returns end of array and subtract one so we don't
                    # include the time past the limit
                    firstTooEarlyIndex = firstTooEarlyIndex + 1 if firstTooEarlyIndex == self.__nRedshifts else firstTooEarlyIndex
                                        
                    # as long as that doesn't preclude the whole range
                    if firstTooEarlyIndex > 0:
                        # redshift at time of formation
                        zOfFormation = timesToRedshifts(formationTime[:firstTooEarlyIndex - 1])

                        # index in the redshift array to which these redshifts correspond
                        zOfFormationIndex = np.ceil(zOfFormation / redshiftStep).astype(int) if self.__nRedshifts > 1 else 0

                        # merger rate for this DCO at z (with z < CI_maxRedshift) = formation rate at z_form
                        thisMergerRate[:firstTooEarlyIndex - 1] = thisFormationRate[zOfFormationIndex]

                        # update merger rate array if necessary
                        if mergerRate is not None: mergerRate[McBin, :firstTooEarlyIndex - 1] += thisMergerRate[:firstTooEarlyIndex - 1]

                        # write to per-binary files if necessary
                        # common string for primary and seconday masses
                        tmpStr = 'D,' + SigFigs.format(self.__COMPAS.PrimaryMass(i)[0]) + ',' + SigFigs.format(self.__COMPAS.SecondaryMass(i)[0]) + ','
                                                
                        # we could remove files if we had errors here - e.g. if we couldn't write the header row to the mergers
                        # file for some reason, we could cleanup - close the file and remove it, but for now I'll leave it there
                        # in case it useful as a debugging tool to find out  what went wrong.  If it becomes an issue, or too 
                        # tedious to clean up the files manually in case of error, we can revist removing then here.

                        # per-binary mergers file
                        if mergersFile is not None:
                            for z in range(self.__nRedshifts):                                                          # for each redshift
                                # construct record
                                rec = tmpStr + SigFigs.format(self.__redshifts[z]) + ',' + SigFigs.format(thisMergerRate[z])
                                errStr = self.__WritePerBinaryRecord(mergersFile, mergersFilename, rec)                 # write mergers record
                                if errStr is None:                                                                      # write ok?
                                    mergersWriteCount += 1                                                              # yes - increment write count
                                else:                                                                                   # no
                                    print('\n!! WARNING !!')                                                            # show warning
                                    print(errStr)
                                    print('Per-binary mergers data will not be written\n')
                                    mergersFile, errStr1 = self.__ClosePerBinaryFile(mergersFile, mergersFilename)      # close the mergers file
                                    mergersWriteCount = -1                                                              # indicate error - reset write count
                                    # we close the yields file here too since only the mergers file is mandatory
                                    # if we can't give the user the mergers data, we don't bother giving them the mergers data
                                    # we can change this behaviour later if we decide the yields data is important on its own
                                    print('Per-binary yields data will not be written\n')
                                    yieldsFile, errStr2 = self.__ClosePerBinaryFile(yieldsFile, yieldsFilename)         # close the yields file
                                    yieldsWriteCount = -1                                                               # indicate error - reset write count
                                    errStr = errStr1 if errStr1 is not None else errStr2                                # set the error return

                        # per-binary yields file
                        if yieldsFile is not None:
                            # construct record
                            rec = tmpStr + SigFigs.format(self.__COMPAS.DelayTime(i)[0]) + ',' + SigFigs.format(self.__COMPAS.Zsystems(i)[0]) + ',' + SigFigs.format(1.0E6 / SFmass)
                            errStr = self.__WritePerBinaryRecord(yieldsFile, yieldsFilename, rec)                       # write yields record
                            if errStr is None:                                                                          # write ok?
                                yieldsWriteCount += 1                                                                   # yes - increment write count
                            else:                                                                                       # no
                                print('\n!! WARNING !!')                                                                # show warning
                                print(errStr)
                                print('Per-binary yields data will not be written\n')                                   # no
                                yieldsFile, errStr = self.__ClosePerBinaryFile(yieldsFile, yieldsFilename)              # close the yields file
                                yieldsWriteCount = -1                                                                   # indicate error - reset write count
                                # we don't close the mergers file here because only the mergers file is mandatory
                                # the yields data is a bonus, so no drama if we can't produce it

                    # initialise detection rate for this DCO
                    thisDetectionRate = np.zeros(shape = (self.__nRedshiftsDetection))

                    # calculate detection probability
                    # use lookup tables to find the probability of detecting the binary at each redshift
                    detectionProbability = np.ones(shape = (self.__nRedshiftsDetection))                                # default detection probabilities

                    McShifted = self.__COMPAS.ChirpMasses(i)[0] * (1.0 + self.__redshifts[:self.__nRedshiftsDetection]) # shift frames for the chirp mass

                    etaIndex  = np.round(self.__COMPAS.ETAs(i)[0] / ETAstep).astype(int) - 1                            # the closest index to the given value of eta
                    McIndex   = np.round(McShifted / McStep).astype(int) - 1                                            # the closest index to the given value of Mc

                    # lookup values for SNR
                    SNRs             = np.ones(self.__nRedshiftsDetection) * 0.00001                                    # default values
                    McBelowMax       = McIndex < self.__SE.SNRgridAt1Mpc().shape[1]                                     # mask below max SNR
                    SNRs[McBelowMax] = self.__SE.SNRgridAt1Mpc()[etaIndex, McIndex[McBelowMax]]                         # SNR lookup values           
                    SNRs             = SNRs / self.__distances[:self.__nRedshiftsDetection]                             # correct for distances

                    # lookup values for detection probability
                    detectionListIndex = np.round(SNRs / SNRstep).astype(int) - 1                                       # index values
                    SNRbelowMax        = detectionListIndex < self.__SE.DetectionProbabilityLen()                       # mask below max SNR
                    SNRbelowMin        = detectionListIndex < 0                                                         # mask below min SNR

                    # remember we set probability = 1 by default? Because if we don't set it here, we have
                    # snr > max snr, which is 1000 by default, meaning very detectable
                    detectionProbability[SNRbelowMax] = self.__SE.DetectionProbabilityFromSNR()[0][detectionListIndex[SNRbelowMax]]

                    # on the other hand, if SNR is too low, the detection probability is effectively zero
                    detectionProbability[SNRbelowMin] = 0
                    # detection rate
                    thisDetectionRate = thisMergerRate[:self.__nRedshiftsDetection] * detectionProbability * \
                                        self.__shellVolumes[:self.__nRedshiftsDetection] / (1.0 + self.__redshifts[:self.__nRedshiftsDetection])

                    # update detection rate array if necessary
                    if detectionRate is not None: detectionRate[McBin] += thisDetectionRate

                    if not createRates and mergersFile is None: break                                                   # something went wrong - bail out (yields file is not mandatory)

                # close per-binary mergers file if necessary
                errStr1 = None
                if mergersFile is not None:
                    mergersFile, errStr1 = self.__ClosePerBinaryFile(mergersFile, mergersFilename)                      # close the mergers file

                # close per-binary yields file if necessary
                errStr2 = None
                if yieldsFile is not None:
                    yieldsFile, errStr2 = self.__ClosePerBinaryFile(yieldsFile, yieldsFilename)                         # close the yields file

                if errStr is None and not (errStr1 is None and errStr2 is None):
                    errStr = errStr1 if errStr1 is not None else errStr2                                                # set the error return
                                        
                # set class member variable if necessary
                if errStr is None: self.__SFmass = SFmass                                                               # for getter

            if errStr is not None: self.__RestoreState()                                                                # restore previously saved state if necessary

            if mergersWriteCount >= 0: mergersWriteCount += 1                                                           # actual records written if no error occurred
            if yieldsWriteCount >= 0 : yieldsWriteCount += 1                                                            # actual records written if no error occurred

            # print postmble if necessary
            if self.__verbose:
                timeTaken = time.time() - mark
                memUsage  = psutil.Process().memory_info().rss / GB
                if mergersWriteCount >= 0: print('\n   Per-binary mergers file \'{:s}\' created with {:d} data records'.format(mergersFilename, mergersWriteCount))
                if yieldsWriteCount >= 0: print('   Per-binary yields file \'{:s}\' created with {:d} data records'.format(yieldsFilename, yieldsWriteCount))
                print('\n   CosmicIntegration::FindRates() completed in {:.2f} seconds.  Process total memory usage = {:.6f} GB\n'.format(timeTaken, memUsage))


        return [[formationRate, mergerRate, detectionRate], McBins], [[mergersWriteCount, mergersFilename], [yieldsWriteCount, yieldsFilename]], errStr


    # Getters

    """
    The following getters return:
    
        - an array value at the index specified, if index parameter is specified
        - the complerte array if the index parameter is not specified

    Each will return an error if an invalid index is specified
    """

    def Redshifts(self, p_Index = None):

        func = self.__className + "::Redshifts(): "

        # default return values
        z      = None
        errStr = None

        if p_Index is None:
            z = self.__redshifts
        elif not (isinstance(p_Index, int) and p_Index >= 0 and p_Index < self.__redshifts.shape[0]):
            errStr = func + 'Index value must be an integer >= 0 and < {:d}'.format(self.__redshifts.shape[0])
        else:
            z = self.__redshifts[p_Index]

        return z, errStr


    # the following getters take no parameters

    def COMPAS(self): return self.__COMPAS
    def SE(self)    : return self.__SE
    def SFmass(self): return self.__SFmass
