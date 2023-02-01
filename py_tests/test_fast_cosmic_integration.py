from compas_python_utils.cosmic_integration.FastCosmicIntegration import find_detection_rate
from compas_python_utils.cosmic_integration.ClassCOMPAS import COMPASData

import io
import time
from contextlib import redirect_stdout

import unittest
import os


class TestFastIntegrator(unittest.TestCase):
    def setUp(self) -> None:
        self.res_dir = os.path.join(
            os.path.dirname(__file__), "test_data", "COMPAS_Output")
        self.res_h5_file = os.path.join(self.res_dir, "COMPAS_Output.h5")
        self.res_h5_file = "/Users/avaj0001/Documents/projects/compas_dev/quasir_compass_blocks/data/COMPAS_Output.h5"

    def test_fast_cosmic_integration(self):
        trap = io.StringIO()
        with redirect_stdout(trap):
            start_time = time.time()
            (
                detection_rate,
                formation_rate,
                merger_rate,
                redshifts,
                COMPAS,
            ) = find_detection_rate(path=self.res_h5_file)
            runtime = time.time() - start_time

        # check that the shape of the detection rate, formation rate and merger rate are the same
        num_rows = [detection_rate.shape[0], formation_rate.shape[0], merger_rate.shape[0]]
        self.assertTrue(len(set(num_rows)) == 1)

        # check that the len of the redshifts is the same as the number of cols in formation and merger rate
        self.assertTrue(len(redshifts) == formation_rate.shape[1])
        self.assertTrue(len(redshifts) == merger_rate.shape[1])

        # check that the COMPAS object is a COMPASData object
        self.assertIsInstance(COMPAS, COMPASData)
