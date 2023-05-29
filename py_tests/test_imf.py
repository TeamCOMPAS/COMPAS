from compas_python_utils.cosmic_integration.ClassCOMPAS import COMPASData, CDF_IMF, inverse_CDF_IMF
from compas_python_utils import cdriver
import matplotlib.pyplot as plt
import numpy as np
import pytest
from tqdm.auto import trange, tqdm
from astropy import units as u
import time

from compas_python_utils.cosmic_integration.ClassCOMPAS import find_star_forming_mass_per_binary_sampling
from compas_python_utils.cosmic_integration.imf import plot_kroupa_imf_pdf


def test_imf(test_archive_dir):
    """Test that the imf is correctly calculated"""
    np.random.seed(0)
    n, bins = 300000, 30
    masses_py = inverse_CDF_IMF(np.random.uniform(size=n))
    hist_py, bins = np.histogram(masses_py, bins=np.geomspace(min(masses_py), max(masses_py), bins))
    hist_py = hist_py / np.sum(hist_py)

    masses_c = cdriver.sample_from_imf(n)
    hist_c, bins = np.histogram(masses_c, bins=bins)
    hist_c = hist_c / np.sum(hist_c)

    plt.plot(bins[:-1], hist_py, color="tab:red", ls=":", label='python', zorder=10)
    plot_kroupa_imf_pdf(
        bins, fig=plt.gcf(), n=n,
        plot_kwgs=dict(color="tab:blue", alpha=0.5, label='cpp',lw=3))
    plt.legend()
    plt.savefig(f"{test_archive_dir}/test_imf.png")

    assert np.allclose(hist_py, hist_c, rtol=1e-2, atol=1e-2)


def test_find_star_forming_mass_per_binary_sampling(test_archive_dir):
    """Test that the star forming mass per binary is correctly calculated"""
    py_kwgs = dict(binaryFraction=0.7, Mlower=5 * u.Msun, Mupper=150 * u.Msun, m2_min=0.1 * u.Msun)
    c_kwgs = dict(binaryFraction=0.7, Mlower=5, Mupper=150, m2_min=0.1)

    n_samples = np.geomspace(1e5, 2e7, 5).astype(int)
    repetitions = 5
    avg_mass = np.zeros((repetitions, len(n_samples)))
    avg_mass_c = np.zeros((repetitions, len(n_samples)))
    for i in trange(repetitions, desc="repetitions"):
        for j, n in tqdm(enumerate(n_samples), desc="n_samples", total=len(n_samples)):
            avg_mass[i, j] = find_star_forming_mass_per_binary_sampling(**py_kwgs, SAMPLES=n).value
            avg_mass_c[i, j] = cdriver.compute_star_forming_mass_per_binary(**c_kwgs, n=n)

    assert avg_mass_c[0, -1] - avg_mass_c[0, -1] < 5 * u.Msun

    def plot_ci(avg_mass, color, label):
        ci = np.quantile(avg_mass, [0.025, 0.5, 0.975], axis=0)
        plt.plot(n_samples, ci[1], color=color, label=label)
        plt.fill_between(n_samples, ci[0], ci[2], alpha=0.5, color=color)

    plot_ci(avg_mass, color="tab:red", label="python")
    plot_ci(avg_mass_c, color="tab:blue", label="cpp")
    plt.xscale("log")
    plt.xlabel("Number of IMF samples")
    plt.ylabel("Average mass/Binary [Msun]")
    plt.legend()
    plt.savefig(f"{test_archive_dir}/test_find_star_forming_mass_per_binary_sampling.png")
