from compas_python_utils.cosmic_integration.FastCosmicIntegration import find_detection_rate
from compas_python_utils.cosmic_integration.ClassCOMPAS import COMPASData

import io
import time
import h5py
from pytest_utils import get_compas_output_path, time_func_and_capture_sout

import unittest
import os


OUTDIR = os.path.join(os.path.dirname(__file__), "output_test")

class TestFastIntegrator(unittest.TestCase):
    def test_fast_cosmic_integration(self):
        """Test that fast cosmic integration works"""
        func_out, runtime, sout = time_func_and_capture_sout(
            find_detection_rate,
            get_compas_output_path()
        )
        (
            detection_rate,
            formation_rate,
            merger_rate,
            redshifts,
            COMPAS,
        ) = func_out

        self.assertLess(runtime, 10)

        # check that the shape of the detection rate, formation rate and merger rate are the same
        num_rows = [detection_rate.shape[0], formation_rate.shape[0], merger_rate.shape[0]]
        self.assertTrue(len(set(num_rows)) == 1)

        # check that the len of the redshifts is the same as the number of cols in formation and merger rate
        self.assertTrue(len(redshifts) == formation_rate.shape[1])
        self.assertTrue(len(redshifts) == merger_rate.shape[1])

        # check that the COMPAS object is a COMPASData object
        self.assertIsInstance(COMPAS, COMPASData)

        # write logs from run to file in OUTDIR
        if not os.path.exists(OUTDIR):
            os.mkdir(OUTDIR)
        with open(os.path.join(OUTDIR, "test_fast_cosmic_integration.log"), "w") as f:
            f.write(sout)


