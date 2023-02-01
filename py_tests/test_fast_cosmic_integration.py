from compas_python_utils.cosmic_integration.FastCosmicIntegration import find_detection_rate

import io
import time
from contextlib import redirect_stdout

from astropy import units

import unittest
import os
import shutil

from compas_python_utils import h5copy
import h5py


class TestFastIntegrator(unittest.TestCase):
    def setUp(self) -> None:
        self.res_dir = os.path.join(
            os.path.dirname(__file__), "test_data", "COMPAS_Output")
        self.res_h5_file = os.path.join(self.res_dir, "COMPAS_Output.h5")


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
