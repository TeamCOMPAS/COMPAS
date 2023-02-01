from compas_python_utils.cosmic_integration.FastCosmicIntegration import find_detection_rate
from compas_python_utils.cosmic_integration.ClassCOMPAS import COMPASData

import io
import time
import h5py
from data_path import get_data_path
from contextlib import redirect_stdout

import unittest
import os


class TestFastIntegrator(unittest.TestCase):
    def test_fast_cosmic_integration(self):
        self.res_h5_file = h5py.File(get_data_path(), 'r')
        trap = io.StringIO()
        with redirect_stdout(trap):
            start_time = time.time()
            (
                detection_rate,
                formation_rate,
                merger_rate,
                redshifts,
                COMPAS,
            ) = find_detection_rate(path=get_data_path())
            runtime = time.time() - start_time

        # check that the shape of the detection rate, formation rate and merger rate are the same
        num_rows = [detection_rate.shape[0], formation_rate.shape[0], merger_rate.shape[0]]
        self.assertTrue(len(set(num_rows)) == 1)

        # check that the len of the redshifts is the same as the number of cols in formation and merger rate
        self.assertTrue(len(redshifts) == formation_rate.shape[1])
        self.assertTrue(len(redshifts) == merger_rate.shape[1])

        # check that the COMPAS object is a COMPASData object
        self.assertIsInstance(COMPAS, COMPASData)
