import os
from typing import Any, Dict

import h5py
from conftest import get_compas_data
from deepdiff import DeepDiff

from compas_python_utils import h5copy


def test_h5copy_copyHDF5File(tmp_path, example_compas_output_path):
    """Test that h5copy Copy file can copy a file"""
    init_file = example_compas_output_path

    # creating new h5 file to copy contents to
    new_file = f"{tmp_path}/copied_compas_out.h5"
    with h5py.File(new_file, "w") as new_h5file:
        h5copy.copyHDF5File(path=init_file, outFile=new_h5file)
    assert os.path.exists(new_file), f"File {new_file} does not exist"
    init_data = get_compas_data(init_file)
    new_data = get_compas_data(new_file)
    diff = DeepDiff(init_data, new_data)
    assert len(diff) == 0, f"The copied file is not the same as the original: {diff}"
