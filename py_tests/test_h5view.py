"""
Tests that the utils/h5view.py module works as expected.
"""

import unittest
import os
import h5py
from data_path import DATA_PATH

from compas_python_utils import h5view
from contextlib import redirect_stdout
from io import StringIO

class TestH5view(unittest.TestCase):
    def test_h5view_print_contents_doest_fail(self):
        """Test that h5view print_contents doesn't fail"""
        h5file = h5py.File(DATA_PATH, 'r')
        with redirect_stdout(StringIO()) as sout:
            h5view.printContents(h5name=DATA_PATH, h5file=h5file)
        sout = sout.getvalue()
        assert sout != ""
        print(sout)
        # This is a pretty weak test atm -- just checks that something is printed,
        # this test is not actually checking the contents of the output

