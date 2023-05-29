from compas_python_utils.cosmic_integration import FastCosmicIntegration
from compas_python_utils.cosmic_integration.ClassCOMPAS import COMPASData, CDF_IMF, inverse_CDF_IMF
from compas_python_utils import cdriver
import matplotlib.pyplot as plt
import numpy as np
import time
import os
import pytest
from tqdm.auto import trange
from astropy import units as u

from compas_python_utils.cosmic_integration.ClassCOMPAS import find_star_forming_mass_per_binary_sampling
from compas_python_utils.cosmic_integration.imf import plot_kroupa_imf_pdf, get_star_forming_mass_per_binary

def test_fast_cosmic_integration(example_compas_output_path, capsys, test_archive_dir):
    """Test that fast cosmic integration works"""
    # example_compas_output_path = "/home/avaj040/Documents/projects/data/COMPAS_data/jeff_data/h5out_5M.h5"

    t0 = time.time()
    (
        detection_rate,
        formation_rate,
        merger_rate,
        redshifts,
        COMPAS,
    ) = FastCosmicIntegration.find_detection_rate(
        path=example_compas_output_path,
    )
    sout = capsys.readouterr().out
    runtime = time.time() - t0

    logfname = os.path.join(test_archive_dir, "test_fast_cosmic_integration.log")
    with open(logfname, "w") as f:
        f.write(sout)

    assert runtime < 10

    # check that the shape of the detection rate, formation rate and merger rate are the same
    num_rows = [detection_rate.shape[0], formation_rate.shape[0], merger_rate.shape[0]]
    assert len(set(num_rows)) == 1

    # check that the len of the redshifts is the same as the number of cols in formation and merger rate
    assert len(redshifts) == formation_rate.shape[1]
    assert len(redshifts) == merger_rate.shape[1]

    # check that the COMPAS object is a COMPASData object
    assert isinstance(COMPAS, COMPASData)


def test_imf(test_archive_dir):
    """Test that the imf is correctly calculated"""
    np.random.seed(0)
    N, bins = 300000, 30
    U = np.random.uniform(size=N)

    masses_py = inverse_CDF_IMF(U)
    hist_py, bins = np.histogram(masses_py, bins=np.geomspace(min(masses_py), max(masses_py), bins))
    hist_py = hist_py / np.sum(hist_py)

    t0 = time.time()
    masses_c = cdriver.sample_from_imf(N)
    time_c = time.time() - t0
    hist_c, _ = np.histogram(masses_c, bins=bins)
    hist_c = hist_c / np.sum(hist_c)

    plt.plot(bins[:-1], hist_py, color="red", ls="--", label=f'python ({time_py:.3E}s)')
    plt.plot(bins[:-1], hist_c, color="blue", ls="-", alpha=0.5, label=f'c ({time_c:.3E}s)')
    plt.xscale("log")
    plt.xlabel("Mass")
    plt.ylabel("PDF")
    plt.legend()
    plt.savefig(f"{test_archive_dir}/test_imf.png")

    assert np.allclose(hist_py, hist_c, rtol=1e-2, atol=1e-2)


@pytest.mark.slow
def test_find_star_forming_mass_per_binary_sampling(test_archive_dir):
    """Test that the star forming mass per binary is correctly calculated"""
    kwgs = dict(binaryFraction=0.7, Mlower=5 * u.Msun, Mupper=150 * u.Msun, m2_min=0.1 * u.Msun)

    n_samples = np.geomspace(1e5, 2e7, 10).astype(int)
    repetitions = 1
    avg_mass = np.zeros((repetitions, len(n_samples)))
    for i in trange(repetitions):
        avg_mass[i, :] = np.array([
            find_star_forming_mass_per_binary_sampling(**kwgs, SAMPLES=n).value for n in n_samples])

    ci = np.quantile(avg_mass, [0.025, 0.5, 0.975], axis=0)
    plt.plot(n_samples, ci[1])
    plt.fill_between(n_samples, ci[0], ci[2], alpha=0.5)
    error_percent = (ci[2] - ci[0]) / ci[1] * 100


    kwgs = dict(binaryFraction=0.7, Mlower=5, Mupper=150, m2_min=0.1)
    avg_mass_c = np.array([cdriver.compute_star_forming_mass_per_binary(**kwgs, n=n) for n in n_samples])*u.Msun

    plt.plot(n_samples, avg_mass_c, color="red", ls="--")
    plt.xscale("log")
    plt.xlabel("Number of samples")
    plt.ylabel("Average mass/Binary [Msun]")
    plt.savefig(f"{test_archive_dir}/test_find_star_forming_mass_per_binary_sampling.png")





from scipy.interpolate import interp1d

def find_formation_and_merger_rates(n_binaries, redshifts, times, time_first_SF, n_formed, dPdlogZ, metallicities, p_draw_metallicity,
                                    COMPAS_metallicites, COMPAS_delay_times, COMPAS_weights=None):
    """
        Find both the formation and merger rates for each binary at each redshift

        Args:
            n_binaries          --> [int]            Number of DCO binaries in the arrays
            redshifts           --> [list of floats] Redshifts at which to evaluate the rates
            times               --> [list of floats] Equivalent of the redshifts in terms of age of the Universe
            n_formed            --> [float]          Binary formation rate (number of binaries formed per year per cubic Gpc) represented by each simulated COMPAS binary
            dPdlogZ             --> [2D float array] Probability of getting a particular logZ at a certain redshift
            metallicities       --> [list of floats] Metallicities at which dPdlogZ is evaluated; if this is None, assume that metallicity weighting is 1 (corresponds to all SFR happening at one fixed metallicity)
            p_draw_metallicity  --> [float]          Probability of drawing a certain metallicity in COMPAS (float because assuming uniform)
            COMPAS_metallicites --> [list of floats] Metallicity of each binary in COMPAS data
            COMPAS_delay_times  --> [list of floats] Delay time of each binary in COMPAS data
            COMPAS_weights      --> [list of floats] Adaptive sampling weights for each binary in COMPAS data (defaults to all 1s for unweighted samples)

        Returns:
            formation_rate      --> [2D float array] Formation rate for each binary at each redshift
            merger_rate         --> [2D float array] Merger rate for each binary at each redshift
    """
    # check if weights were provided, if not use uniform weights
    if COMPAS_weights is None:
        COMPAS_weights = np.ones(n_binaries)

    # initalise rates to zero
    n_redshifts = len(redshifts)
    redshift_step = redshifts[1] - redshifts[0]
    formation_rate = np.zeros(shape=(n_binaries, n_redshifts))
    merger_rate = np.zeros(shape=(n_binaries, n_redshifts))

    times_to_redshifts = interp1d(times, redshifts)

    # make note of the first time at which star formation occured
    t0 = time_first_SF

    # go through each binary in the COMPAS data
    for i in range(n_binaries):
        # if metallicities array is None, assume all SFR happened at one fixed metallicity
        if metallicities is None :
            formation_rate[i, :] = n_formed * COMPAS_weights[i]
        # calculate formation rate (see Neijssel+19 Section 4) - note this uses dPdlogZ for *closest* metallicity
        else:
            formation_rate[i, :] = n_formed * dPdlogZ[:, np.digitize(COMPAS_metallicites[i], metallicities)] / p_draw_metallicity * COMPAS_weights[i]

        # calculate the time at which the binary formed if it merges at this redshift
        tform = times - COMPAS_delay_times[i]

        # we have only calculated formation rate up to z=max(redshifts),
        # so we need to only find merger rates for formation times at z<max(redshifts)
        # first locate the index above which the binary would have formed before z=max(redshifts)
        t0_to_tform_idx = np.digitize(t0, tform)

        # include the whole array if digitize returns end of array and subtract one so we don't include the time past the limit
        t0_to_tform_idx = t0_to_tform_idx + 1 if t0_to_tform_idx == n_redshifts else t0_to_tform_idx

        # as long as that doesn't preclude the whole range
        if t0_to_tform_idx > 0:
            z_of_formation = times_to_redshifts(tform[:t0_to_tform_idx - 1])
            z_of_formation_index = np.ceil(z_of_formation / redshift_step).astype(int)
            merger_rate[i, :t0_to_tform_idx - 1] = formation_rate[i, z_of_formation_index]
    return formation_rate, merger_rate