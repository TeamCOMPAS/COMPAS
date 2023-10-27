# ! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created for Python 3
@author: Sebastian M. Gaebel
@email: sgaebel@star.sr.bham.ac.uk
"""

from __future__ import division, print_function
import h5py
import numpy as np
import os
import scipy.interpolate


# Global variables to reuse objects for multiple calls
_random_thetas = None
_interpolator = None
_sens = None

def detection_probability(m1, m2, redshift, distance, snr_threshold,sensitivity='design'):
    """
    Returns the detection probability of a CBC event with given
    masses and distance.
    This function is a convenience function to generate the
    interpolator with 'SNR_Grid_IMRPhenomPv2_FD_all_noise.hdf5'
    and 'SimNoisePSDaLIGODesignSensitivityP1200087', redshift the
    masses, rescale the SNR using the distance, and then call
    'detection_probability_from_snr' to return the probability.
    Parameters
    ----------
    m1, m2 : float > 0
        Primary and secondary mass in source frame.
    redshift : float >= 0
        Redshift of the waveform through cosmology.
    distance : float > 0
        Luminosity Distance in Mpc to the source.
    snr_threshold : float > 0
        Threshold above which an event is considered to be detected.
    sensitivity : str
        Which detector sensitivity PSD to use. Options are 'design' and 'O1'
    Returns
    -------
    out : float
        Estimate of the detection probability.
    Notes
    -----
    The interpolator is only initialized once and then stored in a
    module level global variable to be reused.
    """

    global _interpolator
    global _sens
    if (_interpolator is None) or (_sens is not sensitivity):
        _sens = sensitivity
        _interpolator = SNRinterpolator(_sens)
    interpolated_snr = _interpolator(m1*(1+redshift), m2*(1+redshift))
    # SNR scales as 1/distance
    interpolated_snr /= distance

    return detection_probability_from_snr(snr_value=interpolated_snr,
                                          snr_threshold=snr_threshold)

def SNRinterpolator(sensitivity='design'):
    """
    Returns an Interpolator class instance for a given sensitivity
    
    This function is a convenience function to generate the
    interpolator with 'SNR_Grid_IMRPhenomPv2_FD_all_noise.hdf5'
    and 'SimNoisePSDaLIGODesignSensitivityP1200087'.
    
    Parameters
    ----------
    sensitivity : str
    Which detector sensitivity PSD to use. Options are 'design' and 'O1'
    
    Returns
    -------
    out : Interpolator
    Interpolator class instance
    
    Notes
    -----
    The interpolator is only initialized once and then stored in a
    module level global variable to be reused.
    """
    path = os.path.dirname(os.path.abspath(__file__))
    
    if sensitivity == 'design':
        hdfDatasetName = 'SimNoisePSDaLIGODesignSensitivityP1200087'
    elif sensitivity == 'O1':
        hdfDatasetName = 'P1500238_GW150914_H1-GDS-CALIB_STRAIN.txt'
    elif sensitivity == 'O3':
        hdfDatasetName = 'SimNoisePSDaLIGOMidHighSensitivityP1200087'

    global _interpolator
    global _sens
    if (_interpolator is None) or (_sens is not sensitivity):
        _sens = sensitivity
        _interpolator = Interpolator(
                                     path+'/SNR_Grid_IMRPhenomPv2_FD_all_noise.hdf5',
                                     hdfDatasetName, mode='scipy')
    return _interpolator

class Interpolator:
    """
    Interpolation class to estimate the values of arbitrary points
    from a given grid. The interpolator class is initialized with
    a mass axis defining the (symmetric) grid, a grid of values
    which are used as reference points, and a operation mode.
    The operation mode determines the interpolation method used and
    may be either 'scipy', which uses 'RectBivariateSpline' from
    'scipy.interpolate', or 'custom' where all non-NaN points adjacent
    to the call are averaged while being weighed by the inverse of
    the distance between the given corner point and queried position.
    The grid is generally assumed to be spaced uniformly in log,
    therefore mass_axis and evaluation masses are transformed to
    log space before interpolating.
    """
    def __init__(self, first_arg, second_arg, mode='scipy'):
        """
        Initialize the interpolator by definig the used reference
        values and interpolation mode.
        The interpolator may be initialized by providing the mass
        axis and SNR grid, or via a path to a stored precomputed
        grid and a string specifying the noise spectrum to use.
        Parameters
        ----------
        first_arg : mass axis as 1D ndarray or path as string
            If given as mass axis, must be a one dimensional
            ndarray. If given as string, must be a valid path
            to a HDF5 file containing mass axis and SNR grid
            as produced by 'generate_SNR_files.py'.
        second_arg : SNR grid as 2D ndarray or noise spectrum as string
            If given as SNR grid, must be a two dimensional
            ndarray with each dimension being the length of the mass
            axis. If given as noise curve, must be the string specifying
            a valid group of the HDF5 file given above.
        mode : str, optional
            Valid values are 'scipy' and 'custom', all other will
            raise a ValueError.
        Returns
        -------
        out : None
        """
        if isinstance(first_arg, str) and isinstance(second_arg, str):
            if not os.path.isfile(first_arg):
                raise FileNotFoundError('HDF5 file expected: %r' % first_arg)
            with h5py.File(first_arg, 'r') as hdf:
                mass_axis = hdf['mass_axis'][...]
                if second_arg not in hdf['snr_values']:
                    raise ValueError('Group %r not found.' % second_arg)
                snr_grid = hdf['snr_values'][second_arg][...]
        else:
            mass_axis = first_arg
            snr_grid = second_arg

        if mode == 'scipy':
            self.use_scipy = True
            self.interpolator = scipy.interpolate.RectBivariateSpline(
                np.log(mass_axis), np.log(mass_axis), snr_grid)
        elif mode == 'custom':
            self.use_scipy = False
            self.mass_axis = np.log(mass_axis)
            self.snr_grid = snr_grid
        else:
            raise ValueError('Invalid mode: %r' % mode)

    def __call__(self, m1, m2):
        """
        Obtain an interpolated value for the given position.
        Parameters
        ----------
        m1, m2 : float
            Positive values for which the interpolated SNR is calculated.
        Returns
        -------
        out : float
            Estimate for the SNR at (m1, m2).
        """
        if self.use_scipy:
            return self.interpolator(np.log(m1), np.log(m2), grid=False)
        else:
            return self._custom(np.log(m1), np.log(m2))

    def _custom(self, m1, m2):
        """
        Custom function to approximate the value of any given point
        within the limits of the grid. The value is calculated as the
        average value of the 4 corners of the grid, weighed by the
        inverse distance to (m1, m2). Points for which the SNR is
        NaN are ignored.
        Parameters
        ----------
        m1, m2 : float
            Positive values for which the interpolated SNR is calculated.
        Returns
        -------
        out : float
            Estimate for the SNR at (m1, m2).
        """
        # The point is in the lower triangle by definition
        if isinstance(m1, np.ndarray):
            ret_values = np.empty_like(m1)
            for i, (m1_, m2_) in enumerate(zip(m1, m2)):
                ret_values[i] = self._custom(m1_, m2_)
            return ret_values
        elif m1 < m2:
            m1, m2 = m2, m1
        assert np.min(self.mass_axis) < m1 < np.max(self.mass_axis)
        assert np.min(self.mass_axis) < m2 < np.max(self.mass_axis)
        # Find the indices of the next higher element in the axis
        # this ensures mass_axis[i-1] < m1 < mass_axis[i]
        # Edge cases are covered by the assertions, which ensure no pair are on
        # the outer boundary.
        i, j = np.searchsorted(self.mass_axis, [m1, m2])
        m1_selection = self.mass_axis[[i-1, i-1, i, i]]
        m2_selection = self.mass_axis[[j-1, j, j-1, j]]
        snr_selection = self.snr_grid[[i-1, i-1, i, i], [j-1, j, j-1, j]]
        distance = np.sqrt((m1_selection - m1)**2 + (m2_selection - m2)**2)
        # If the distance is zero, return the value at that position
        if np.any(distance == 0):
            return snr_selection[distance == 0]
        # Ensure there are no NaN in the SNR values,which may happen for the
        # upper triangle.
        distance = distance[~np.isnan(snr_selection)]
        snr_selection = snr_selection[~np.isnan(snr_selection)]
        if snr_selection.size == 0:
            raise ValueError('No non-NaN values surrounding (m1, m2).')
        return np.average(snr_selection, weights=1./distance)


def detection_probability_from_snr(snr_value, snr_threshold, n_thetas=1e6):
    """
    Compute the probability of detecting an CBC with given SNR and
    threshold, averaging over all orientations and sky positions.
    Based of Finn & Chernoff 1993 (https://arxiv.org/abs/gr-qc/9301003).
    Parameters
    ----------
    snr_value : float or array of floats
        SNR value of the event under ideal conditions. Must be positive.
    snr_threshold : float
        SNR threshold for detection. Must be positive.
    n_thetas : unsigned int
        Number of random realizations of inclination and sky position
        used to calculate the average probability of detecting the
        event. Default: 1e6.
    Returns
    -------
    out : float or array of floats (same as 'snr_value')
        Probability of the received signal to be above the threshold.
    """
    # Check input values
    assert np.all(snr_value > 0), repr(snr_value)
    assert snr_threshold > 0, repr(snr_threshold)
    assert n_thetas > 0, repr(n_thetas)
    n_thetas = int(n_thetas)

    # Generate thetas if necessary
    global _random_thetas
    if _random_thetas is None or _random_thetas.size != n_thetas:
        cos_thetas = np.random.uniform(low=-1, high=1, size=n_thetas)
        cos_incs = np.random.uniform(low=-1, high=1, size=n_thetas)
        phis = np.random.uniform(low=0, high=2*np.pi, size=n_thetas)
        zetas = np.random.uniform(low=0, high=2*np.pi, size=n_thetas)

        Fps = (0.5 * np.cos(2*zetas) * (1+cos_thetas**2)*np.cos(2*phis) -
               np.sin(2*zetas) * cos_thetas*np.sin(2*phis))  # 3.29a
        Fxs = (0.5 * np.sin(2*zetas) * (1+cos_thetas**2)*np.cos(2*phis) +
               np.cos(2*zetas) * cos_thetas*np.sin(2*phis))  # 3.29b
        _random_thetas = np.sqrt(0.25 * Fps**2 * (1 + cos_incs**2)**2 +
                                 Fxs**2 * cos_incs**2)  # 3.31
        _random_thetas = np.sort(_random_thetas)

    # Calculate and return the probability of detection
    # From Finn & Chernoff 1993, we have SNR ~ theta*integrand, assuming
    # that the polarisations are orthogonal.
    # Remember: snr_value is the maximum possible SNR
    theta_min = snr_threshold / snr_value
    # theta_min is converted to a 1-dimensional array so this function
    # can be used with array and scalar arguments.


    #I have replaced this easy to read for loop with an unreadable (but much faster) one liner. Leaving it commented so
    #that future people can understand what the hell it's doing

    # theta_min = np.atleast_1d(theta_min).flatten()
    # detection_prob = np.zeros_like(theta_min)
    # for i, theta_min_value in enumerate(theta_min):
    #     if theta_min_value <= 1:
    #         detection_prob[i] = np.mean(_random_thetas > theta_min_value)

    theta_min = np.atleast_1d(theta_min).flatten()
    detection_prob = np.zeros_like(theta_min)
    detection_prob[theta_min <= 1] = 1.-((np.digitize(theta_min[theta_min <= 1], _random_thetas)-1.)/float(n_thetas))

    if len(detection_prob) == 1 and np.array(snr_value).ndim == 0:
        # If the original arg was 0-dimensional (scalar), return a scalar.
        return detection_prob[0]
    return detection_prob.reshape(np.array(snr_value).shape)

