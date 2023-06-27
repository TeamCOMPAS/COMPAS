from compas_python_utils.cosmic_integration.totalMassEvolvedPerZ import (
    IMF, get_COMPAS_fraction, analytical_star_forming_mass_per_binary_using_kroupa_imf,
    star_forming_mass_per_binary, draw_samples_from_kroupa_imf, inverse_sample_IMF,
)
import numpy as np
import matplotlib.pyplot as plt
import h5py as h5

import pytest

MAKE_PLOTS = False


M1_MIN = 5
M1_MAX = 150
M2_MIN = 0.1

def test_imf(test_archive_dir):
    mmin, mmax = 0.01, 200
    ms = np.geomspace(mmin, mmax, 100)
    m_breaks = [0.01, 0.08, 0.5, 200]
    powers = [0.3, 1.3, 2.3]
    imf_vals = IMF(ms, *m_breaks, *powers)

    assert IMF(5e-3, *m_breaks, *powers) == 0, "IMF(m<min-m) ==0"
    assert np.sum(imf_vals < 0) == 0, "IMF must be > 0"

    imf_samples = inverse_sample_IMF(n_samples=int(1e5), m_min=mmin, m_max=mmax)
    hist, edges = np.histogram(imf_samples, bins=10)
    assert np.sum(hist[0]) > 0, "sample_IMF must return some samples in [10,100]"

    if MAKE_PLOTS:
        fig = plot_imf(m_breaks, imf_vals, imf_samples, mmin, mmax, ms)
        fig.savefig(f"{test_archive_dir}/IMF.png")


def test_compas_fraction():
    assert get_COMPAS_fraction(
        m1_low=10, m1_upp=100, m2_low=10, f_bin=0
    ) == 0
    assert get_COMPAS_fraction(
        m1_low=0.01, m1_upp=200, m2_low=0, f_bin=1,
        m1=0.01, m2=0.08, m3=0.5, m4=200
    ) == 1
    assert 0 < (get_COMPAS_fraction(
        m1_low=0.01, m1_upp=100, m2_low=0, f_bin=0.5,
        m1=0.01, m2=0.08, m3=0.5, m4=200
    )) < 1


def test_analytical_function():
    default_case = analytical_star_forming_mass_per_binary_using_kroupa_imf(
        m1_max=150,
        m1_min=5,
        m2_min=0.1,
        fbin=1
    )
    assert 79.0 < default_case < 79.2


def test_analytical_vs_numerical_star_forming_mass_per_binary(fake_compas_output, tmpdir, test_archive_dir):
    np.random.seed(42)
    m1_max = M1_MAX
    m1_min = M1_MIN
    m2_min = M2_MIN
    fbin = 1

    numerical = star_forming_mass_per_binary(fake_compas_output, m1_min, m1_max, m2_min, fbin)
    analytical = analytical_star_forming_mass_per_binary_using_kroupa_imf(m1_min, m1_max, m2_min, fbin)

    assert numerical > 0
    assert analytical > 0

    assert np.isclose(numerical, analytical, rtol=1)
    if MAKE_PLOTS:
        fig = plot_star_forming_mass_per_binary_comparison(tmpdir, analytical, m1_min, m1_max, m2_min, fbin)
        fig.savefig(f"{test_archive_dir}/analytical_vs_numerical.png")


@pytest.fixture
def fake_compas_output(tmpdir)->str:
    return generate_fake_result(tmpdir, n_samples=int(1e4))


def generate_fake_result(tmpdir, n_samples):
    """
    Create a fake COMPAS output file with a group 'BSE_System_Parameters' containing
    the following 1D datasets:
    - Metallicity@ZAMS(1)
    - Mass@ZAMS(1)
    - Mass@ZAMS(2)

    The values of these datasets should be from the IMF function.
    """
    compas_path = f"{tmpdir}/COMPAS.h5"
    m1, m2 = draw_samples_from_kroupa_imf(
        n_samples=n_samples, Mlower=M1_MIN, Mupper=M1_MAX, m2_low=M2_MIN)
    with h5.File(compas_path, "w") as f:
        f.create_group("BSE_System_Parameters")
        f["BSE_System_Parameters"].create_dataset(
            "Metallicity@ZAMS(1)", data=np.linspace(0, 1, len(m1))
        )
        f["BSE_System_Parameters"].create_dataset(
            "Mass@ZAMS(1)", data=m1
        )
        f["BSE_System_Parameters"].create_dataset(
            "Mass@ZAMS(2)", data=m2
        )
    return compas_path



def plot_star_forming_mass_per_binary_comparison(
        tmpdir, analytical, m1_min, m1_max, m2_min, fbin,
        nreps = 5, nsamps = 5
):
    plt.axhline(analytical, color='tab:blue', label="analytical", ls='--')
    n_samps = np.geomspace(1e3, 5e4, nsamps)
    numerical_vals = []
    for _ in range(nreps):
        vals = np.zeros(len(n_samps))
        for i, n in enumerate(n_samps):
            res = generate_fake_result(tmpdir, n_samples=int(n))
            vals[i] = (star_forming_mass_per_binary(res, m1_min, m1_max, m2_min, fbin))
        numerical_vals.append(vals)

    # plot the upper and lower bounds of the numerical values
    numerical_vals = np.array(numerical_vals)
    lower = np.percentile(numerical_vals, 5, axis=0)
    upper = np.percentile(numerical_vals, 95, axis=0)
    plt.fill_between(
        n_samps, lower, upper, color='tab:orange', alpha=0.3,
        linewidth=0
    )
    plt.plot(n_samps, np.median(numerical_vals, axis=0), color='tab:orange', label="numerical")
    plt.xscale("log")
    plt.ylabel("Star forming mass per binary [M$_{\odot}$]")
    plt.xlabel("Number of samples")
    plt.xlim(min(n_samps), max(n_samps))
    plt.legend()
    return plt.gcf()

def plot_imf(m_breaks, imf_vals, imf_samples, mmin, mmax, ms):
    for mi in m_breaks:
        plt.axvline(mi, zorder=-10, color='gray', alpha=0.2)
    plt.plot(ms, imf_vals, label="IMF function")
    plt.hist(imf_samples, bins=np.geomspace(mmin, mmax, 50), histtype='step',
             density=True, alpha=0.5, label="sampled")
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel(r"Mass [M$_{\odot}]$")
    plt.ylabel("IMF")
    plt.legend()
    return plt.gcf()