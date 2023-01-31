"""
Tests that the utils/h5view.py module works as expected.
"""

import unittest
import os
import h5py
from compas_python_utils import h5view


class TestH5view(unittest.TestCase):
    def setUp(self) -> None:
        self.res_dir = os.path.join(
            os.path.dirname(__file__), "test_data", "COMPAS_Output")

    def test_h5view_print_contents_doest_fail(self):
        """Test that h5view print_contents doesn't fail"""
        h5name = os.path.join(self.res_dir, "COMPAS_Output.h5")
        h5file = h5py.File(h5name, 'r')
        h5view.printContents(h5name=h5name, h5file=h5file)
        ## NOTE: this is a bit of a weak test -- as everything goes to stdout
        ## ideally the functions should be a bit more modularised...
