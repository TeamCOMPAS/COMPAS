"""
Tests that the utils/h5view.py module works as expected.
"""

import unittest
import os
import h5py
from pytest_utils import get_compas_output_path, time_func_and_capture_sout

from compas_python_utils import h5view
from contextlib import redirect_stdout
from io import StringIO

class TestH5view(unittest.TestCase):
    def test_h5view_print_contents_doest_fail(self):
        """Test that h5view print_contents doesn't fail"""
        h5file = h5py.File(get_compas_output_path(), 'r')
        _, dur, sout = time_func_and_capture_sout(
            h5view.printContents,
            h5name=get_compas_output_path(), h5file=h5file
        )
        self.assertTrue(sout != "")
        # This is a pretty weak test atm -- just checks that something is printed,
        # this test is not actually checking the contents of the output

