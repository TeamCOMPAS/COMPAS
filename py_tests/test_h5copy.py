import unittest
import os
import shutil

from pytest_utils import get_compas_output_path, time_func_and_capture_sout
from compas_python_utils import h5copy
import h5py

OUTDIR = os.path.join(os.path.dirname(__file__), "output_test")

class TestH5Copy(unittest.TestCase):
    def setUp(self) -> None:
        self.init_file = get_compas_output_path()
        self.new_file = f"{OUTDIR}/copied_compas_out.h5"
        if not os.path.exists(OUTDIR):
            os.mkdir(OUTDIR)
        shutil.copy(self.init_file, self.new_file)
        # for some reason the copyHDF5File function doesnt work if the file doesnt exist...?

    def test_h5copy_copyHDF5File(self):
        """Test that h5copy Copy file can copy a file"""
        new_h5file = h5py.File(self.new_file, 'r')
        h5copy.copyHDF5File(path=self.init_file, outFile=new_h5file)
        self.assertTrue(os.path.exists(self.new_file))
        # TODO: read in files and do a DeepDiff to verify they are the same?
        # ATM its not clear is there is a COMPASOutput reader (the h5view doenst return an COMPASOutput object)
