"""
Tests that the utils/h5view.py module works as expected.
"""

import h5py

from compas_python_utils import h5view


def test_h5view_print_contents_doest_fail(example_compas_output_path, capsys):
    """Test that h5view print_contents doesn't fail"""
    h5file = h5py.File(example_compas_output_path, 'r')
    h5view.printContents(
        h5name=example_compas_output_path, h5file=h5file
    )
    sout = capsys.readouterr().out
    assert f"Contents of HDF5 file {example_compas_output_path}" in sout
    assert "COMPAS file: Run_Details" in sout
