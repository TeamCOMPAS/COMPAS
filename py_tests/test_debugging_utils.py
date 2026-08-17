import h5py
import numpy as np
import pandas as pd

from compas_python_utils.debugging_utils import (
    build_event_string,
    convert_bytes_array_to_strings,
    get_event_history,
    get_event_strings,
    get_mt_data_tuple,
    get_sn_data_tuple,
    print_compas_details_dataframe,
)


def _make_mt_group(file_handle, name="BSE_RLOF"):
    group = file_handle.create_group(name)
    group.create_dataset("SEED", data=np.array([1, 1, 2], dtype=np.int64))
    group.create_dataset("Time<MT", data=np.array([1.0, 2.0, 3.0]))
    group.create_dataset("Stellar_Type(1)<MT", data=np.array([1, 2, 3], dtype=np.int64))
    group.create_dataset("Stellar_Type(2)<MT", data=np.array([2, 3, 4], dtype=np.int64))
    group.create_dataset("RLOF(1)>MT", data=np.array([1, 0, 0], dtype=np.int64))
    group.create_dataset("RLOF(2)>MT", data=np.array([0, 1, 0], dtype=np.int64))
    group.create_dataset("CEE>MT", data=np.array([0, 1, 0], dtype=np.int64))
    group.create_dataset("Merger", data=np.array([0, 0, 1], dtype=np.int64))
    for key in group.keys():
        group[key].attrs["units"] = "dummy"
    return group


def _make_sn_group(file_handle, name="BSE_Supernovae"):
    group = file_handle.create_group(name)
    group.create_dataset("SEED", data=np.array([1, 2], dtype=np.int64))
    group.create_dataset("Time", data=np.array([5.0, 6.0]))
    group.create_dataset("Stellar_Type_Prev(SN)", data=np.array([1, 3], dtype=np.int64))
    group.create_dataset("Stellar_Type(SN)", data=np.array([13, 14], dtype=np.int64))
    group.create_dataset("Supernova_State", data=np.array([1, 2], dtype=np.int64))
    group.create_dataset("Unbound", data=np.array([0, 1], dtype=np.int64))
    for key in group.keys():
        group[key].attrs["units"] = "dummy"
    return group


def _make_system_group(file_handle, name="BSE_System_Parameters"):
    group = file_handle.create_group(name)
    group.create_dataset("SEED", data=np.array([1, 2], dtype=np.int64))
    for key in group.keys():
        group[key].attrs["units"] = "dummy"
    return group


def test_convert_bytes_array_to_strings():
    """Test that the function correctly converts a numpy array of bytes to strings."""
    values = np.array([b"a", b"b"], dtype="S1")
    assert np.array_equal(convert_bytes_array_to_strings(values), np.array(["a", "b"]))

    text_values = np.array(["c", "d"], dtype=str)
    assert np.array_equal(convert_bytes_array_to_strings(text_values), text_values)


def test_print_compas_details_dataframe_accepts_seed_mask(tmp_path):
    """Test that the function correctly filters the dataframe based on a seed mask."""
    path = tmp_path / "COMPAS_debugging.h5"
    with h5py.File(path, "w") as f:
        mt_group = _make_mt_group(f)
        mask = np.array([True, False, True])
        df = print_compas_details_dataframe(mt_group, mask=mask)

    assert isinstance(df, pd.DataFrame)
    expected_keys = {"Time<MT", "Stellar_Type(1)<MT", "Stellar_Type(2)<MT", "RLOF(1)>MT", "RLOF(2)>MT", "CEE>MT", "Merger"}
    assert set(df.index) == expected_keys
    assert df.columns.tolist()[:2] == ["(units)", 1]
    assert df.loc["Time<MT", 1] == 1.0


def test_print_compas_details_dataframe_run_details_without_seed(tmp_path):
    """Test that the function correctly handles a run details group without a seed dataset."""
    path = tmp_path / "COMPAS_run_details.h5"
    with h5py.File(path, "w") as f:
        run_group = f.create_group("Run_Details")
        run_group.create_dataset("Mass_1", data=np.array([10.0]))
        run_group["Mass_1"].attrs["units"] = "Msun"
        run_group.create_dataset("Mass_1-Derivation", data=np.array([b"input"]))

        df = print_compas_details_dataframe(run_group)

    assert isinstance(df, pd.DataFrame)
    assert "Parameter" in df.columns
    assert "Derivation" in df.columns
    assert df.loc["Mass_1", "Derivation"] == "input"


def test_get_mt_data_tuple_deduplicates_and_sorts(tmp_path):
    """Test that the function correctly deduplicates and sorts the data based on the SEED dataset."""
    path = tmp_path / "COMPAS_mt.h5"
    with h5py.File(path, "w") as f:
        mt_group = _make_mt_group(f)
        seeds, events, times = get_mt_data_tuple(mt_group)

    assert seeds == [1, 2]
    assert len(events[0]) == 2
    assert times[0] == [1.0, 2.0]
    assert events[1] == [(3, 4, False, False, False, True)]


def test_get_event_history_and_event_strings(tmp_path):
    """Test that the function correctly retrieves the event history and generates event strings."""
    path = tmp_path / "COMPAS_events.h5"
    with h5py.File(path, "w") as f:
        _make_system_group(f)
        mt_group = _make_mt_group(f)
        sn_group = _make_sn_group(f)
        all_seeds, all_events = get_event_history({
            "BSE_System_Parameters": f["BSE_System_Parameters"],
            "BSE_RLOF": mt_group,
            "BSE_Supernovae": sn_group,
        })
        strings = get_event_strings(all_events=all_events)

    assert all_seeds == [1, 2]
    assert all_events[0][0][0] == "RL"
    assert isinstance(strings[0], str)
    assert "MS" in strings[0] or "M" in strings[0]
    assert len(strings) == 2


def test_build_event_string_for_sn_and_mt_events():
    """Test that the function correctly builds event strings for both MT and SN events."""
    mt_event = ("RL", 1.5, 1, 2, True, False, False, False)
    sn_event = ("SN", 10.0, 1, 13, 1, False)

    mt_string = build_event_string([mt_event]).item()
    sn_string = build_event_string([sn_event]).item()

    assert "MS" in mt_string and ">" in mt_string
    assert "MS" in sn_string and "*" in sn_string and "NS" in sn_string
