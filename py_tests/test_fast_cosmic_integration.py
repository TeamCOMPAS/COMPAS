from compas_python_utils.cosmic_integration import FastCosmicIntegration
from compas_python_utils.cosmic_integration.ClassCOMPAS import COMPASData

import time
import os


def test_fast_cosmic_integration(example_compas_output_path, capsys, test_archive_dir):
    """Test that fast cosmic integration works"""
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
