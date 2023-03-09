import os
from compas_python_utils import h5copy
import h5py


def test_h5copy_copyHDF5File(tmp_path, example_compas_output_path):
    """Test that h5copy Copy file can copy a file"""
    init_file = example_compas_output_path

    # creating new h5 file to copy contents to
    new_file = f"{tmp_path}/copied_compas_out.h5"
    new_h5file = h5py.File(new_file, 'w')

    h5copy.copyHDF5File(path=init_file, outFile=new_h5file)
    assert os.path.exists(new_file)
    # TODO: read in files and do a DeepDiff to verify they are the same?
    # ATM its not clear is there is a COMPASOutput reader (the h5view doenst return an COMPASOutput object)
