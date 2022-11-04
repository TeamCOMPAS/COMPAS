import unittest
import os
import shutil

from compas import h5copy
import h5py


class TestH5Copy(unittest.TestCase):
    def setUp(self) -> None:
        self.res_dir = os.path.join(
            os.path.dirname(__file__), "test_data", "COMPAS_Output")
        self.out_dir = "copy_outdir"
        self.new_file = os.path.join(self.out_dir, "new_COMPAS_Output.h5")
        os.makedirs(self.out_dir, exist_ok=True)

    def tearDown(self) -> None:
        if os.path.exists(self.out_dir):
            shutil.rmtree(self.out_dir)

    def test_h5copy_copyHDF5File(self):
        """Test that h5copy Copy file can copy a file"""
        h5name = os.path.join(self.res_dir, "COMPAS_Output.h5")
        new_h5file = h5py.File(self.new_file, 'r')
        h5copy.copyHDF5File(path=h5name, outFile=new_h5file)
        self.assertTrue(os.path.exists(self.new_file))
        # TODO: read in files and do a DeepDiff to verify they are the same?
        # ATM its not clear is there is a COMPASOutput reader (the h5view doenst return an COMPASOutput object)
