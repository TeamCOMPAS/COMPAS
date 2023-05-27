from compas_python_utils.cosmic_integration import FastCosmicIntegration
from compas_python_utils.cosmic_integration.ClassCOMPAS import COMPASData, CDF_IMF, inverse_CDF_IMF

import imf

import time
import os


def test_fast_cosmic_integration(example_compas_output_path, capsys, test_archive_dir):
    """Test that fast cosmic integration works"""
    example_compas_output_path = "/home/avaj040/Documents/projects/data/COMPAS_data/jeff_data/h5out_5M.h5"

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

    # write logs from run to file in OUTDIR
    assert False


def test_imf():
    """Test that the imf is correctly calculated"""
    import matplotlib.pyplot as plt
    import numpy as np
    import time

    np.random.seed(0)
    U = np.random.uniform(size=200000)

    t0 = time.time()
    masses =    imf.inverse_imf(U)
    print(time.time() - t0)

    hist, bins = np.histogram(masses, bins=np.geomspace(min(masses), max(masses), 100))
    hist = hist / np.sum(hist)
    plt.plot(bins[:-1], hist)
    plt.xscale("log")

    mmax = 200

    # twin axis
    t0 = time.time()
    masses = inverse_CDF_IMF(U)
    print(time.time() - t0)

    hist, bins = np.histogram(masses, bins=np.geomspace(min(masses), max(masses), 100))


    hist = hist / np.sum(hist)
    ax2 = plt.gca().twinx()
    ax2.plot(bins[:-1], hist, color="red", ls="--")
    ax2.set_xscale("log")
    plt.show()
    assert False




from compas_python_utils.cosmic_integration.ClassCOMPAS import u


def test_c():
    from compas_python_utils import cdriver
    vals = cdriver.sample_from_imf(100) * u.Msun
    print(vals)