from compas_python_utils.cosmic_integration import FastCosmicIntegration
from compas_python_utils.cosmic_integration.ClassCOMPAS import COMPASData, CDF_IMF, inverse_CDF_IMF
from compas_python_utils import cdriver
import matplotlib.pyplot as plt
import numpy as np
import time
import os


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

    t0 = time.time()
    masses_py = inverse_CDF_IMF(U)
    time_py = time.time() - t0
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
