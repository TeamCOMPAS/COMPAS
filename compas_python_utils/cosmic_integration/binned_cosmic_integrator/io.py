import h5py
from typing import Dict
import numpy as np


def recursively_load_dict_contents_from_group(h5file: h5py.File, group: str):
    output = dict()
    for key, item in h5file[group].items():
        if isinstance(item, h5py.Dataset):
            output[key] = decode_from_hdf5(item[()])
        elif isinstance(item, h5py.Group):
            output[key] = recursively_load_dict_contents_from_group(
                h5file, group + key + "/"
            )
    return output


def recursively_save_dict_contents_to_group(h5file: h5py.File, group: str, dic: Dict):
    for key, item in dic.items():
        item = encode_for_hdf5(key, item)
        if isinstance(item, dict):
            recursively_save_dict_contents_to_group(h5file, group + key + "/", item)
        else:
            h5file[group + key] = item


def encode_for_hdf5(key, item):
    if isinstance(item, (np.generic, int, float, complex)):
        if isinstance(item, np.integer):
            item = int(item)
        elif isinstance(item, np.floating):
            item = float(item)
        elif isinstance(item, np.complexfloating):
            item = complex(item)
    if isinstance(item, np.ndarray):
        if item.dtype.kind == 'U':
            item = np.array(item, dtype='S')
    if isinstance(item, (np.ndarray, int, float, complex, str, bytes)):
        output = item
    elif item is None:
        output = "__none__"
    elif isinstance(item, list):
        if len(item) == 0:
            output = item
        elif isinstance(item[0], (str, bytes)) or item[0] is None:
            output = []
            for value in item:
                if isinstance(value, str):
                    output.append(value.encode("utf-8"))
                elif isinstance(value, bytes):
                    output.append(value)
                else:
                    output.append(b"__none__")
        elif isinstance(item[0], (int, float, complex)):
            output = np.array(item)
        else:
            raise ValueError(f'Cannot save {key}: {type(item)} type')
    elif isinstance(item, dict):
        output = item.copy()
    else:
        raise ValueError(f'Cannot save {key}: {type(item)} type')
    return output


def decode_from_hdf5(item):
    if isinstance(item, str) and item == "__none__":
        output = None
    elif isinstance(item, bytes) and item == b"__none__":
        output = None
    elif isinstance(item, (bytes, bytearray)):
        output = item.decode("utf-8")
    elif isinstance(item, np.ndarray):
        if item.size == 0:
            output = item
        elif "|S" in str(item.dtype) or (item.dtype.kind == 'S') or \
                (item.size > 0 and isinstance(item.flat[0], bytes)):
            output = [it.decode("utf-8") for it in item]
        else:
            output = item
    elif isinstance(item, (np.bool_, bool)):
        output = bool(item)
    else:
        output = item
    return output