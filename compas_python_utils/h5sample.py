"""Module to sample an COMPAS h5 file.

This allows users to sample a COMPAS h5 file to contain a smaller set of systems,
or upsample a COMPAS h5 file (sample with replacements) to contain a larger set of systems.
"""
import argparse
from typing import Optional, List
import sys

import h5py
import numpy as np

from compas_python_utils.h5view import printSummary
from compas_python_utils.h5copy import copyHDF5File


def sample_h5(
        compas_h5_filepath: str,
        output_filepath: Optional[str] = None,
        n: Optional[int] = None,
        frac: Optional[float] = None,
        replace: Optional[bool] = False,
        seed_group: Optional[str] = "BSE_System_Parameters",
        seed_key: Optional[str] = "SEED",
):
    """Sample an COMPAS h5 file.

    Parameters
    ----------
    compas_h5_filepath : str
        Path to the COMPAS h5 file.
    output_filepath : Optional[str], optional
        Path to the output COMPAS h5 file, by default None.
        If None, the output filepath will be the input filepath with "_sampled" appended.
    n : Optional[int], optional
        Number of binaries to sample, by default None.
    frac : Optional[float], optional
        Fraction of binaries to sample, by default None.
    replace : Optional[bool], optional
        Sample with or without replacement, by default False.
    seed_group : Optional[str], optional
        Group to get binary seed list to sample from, by default "BSE_System_Parameters".
    seed_key : Optional[str], optional
        Key to get binary seed list to sample from, by default "SEED".

    Raises
    ------
    ValueError
        Must specify either n or frac.
        n must be greater than 0.
        Cannot sample without replacement more than the number of binaries. Set replace=True.

    Returns
    -------
    None
        A new COMPAS h5 file will be created at the output filepath.
    """
    if n is None and frac is None:
        raise ValueError("Must specify either n or frac.")
    if n is not None and frac is not None:
        raise ValueError("Must specify either n or frac, not both.")

    if output_filepath is None:
        output_filepath = compas_h5_filepath.replace(".h5", "_sampled.h5")

    with h5py.File(compas_h5_filepath, "r") as compas_h5_file:
        binary_seeds = compas_h5_file[seed_group][seed_key][:]

    if frac is not None:
        n = int(frac * len(binary_seeds))
    else:
        frac = n / len(binary_seeds)

    if n <= 0:
        raise ValueError("n must be greater than 0.")

    if n > len(binary_seeds) and not replace:
        raise ValueError(
            "Cannot sample without replacement more than the number of binaries. "
            "Set replace=True."
        )

    print(f"Sampling {n} ({frac * 100:.2f}%) binaries from {compas_h5_filepath}.")

    print("Original file summary:")
    with h5py.File(compas_h5_filepath, "r") as compas_h5_file:
        printSummary(compas_h5_filepath, compas_h5_file)

    with h5py.File(output_filepath, 'w') as out_h5_file:
        copyHDF5File(compas_h5_filepath, out_h5_file)
        sampled_binary_seeds = np.random.choice(binary_seeds, size=n, replace=replace)
        _sample(out_h5_file, seed_key, sampled_binary_seeds)

    print("Sampled file summary:")
    with h5py.File(output_filepath, "r") as output_h5_file:
        printSummary(output_filepath, output_h5_file)


def _sample(h5_file: h5py.File, sample_key: str, sample_values: np.ndarray):
    for group_name, group in h5_file.items():

        if sample_key not in group:
            continue  # if the group doesn't have the axis key, skip this group

        group_values = group[sample_key][:]

        # get group_seed index for values in sampled_binary_seeds
        sample_idx = []
        for i, s in enumerate(sample_values):
            sample_idx.append(np.where(group_values == s))
        sample_idx = np.concatenate(sample_idx).ravel()
        new_shape = (len(sample_idx),)

        for dataset_name, dataset in group.items():
            # get dataset data at the sampled binary seed index
            resampled_data = dataset[:][sample_idx]
            dataset.resize(new_shape)
            dataset[:] = resampled_data
    return h5_file


def create_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Sample an COMPAS h5 file.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "compas_h5_filepath",
        type=str,
        help="Path to the COMPAS h5 file.",
    )
    parser.add_argument(
        "--output_filepath",
        type=str,
        help="Path to the output COMPAS h5 file. "
             "If None, the output filepath will be the input filepath with '_sampled' appended.",
        default=None,
    )
    parser.add_argument(
        "--n",
        type=int,
        help="Number of binaries to sample.",
        default=None,
    )
    parser.add_argument(
        "--frac",
        type=float,
        help="Fraction of binaries to sample.",
        default=None,
    )
    parser.add_argument(
        "--replace",
        type=bool,
        help="Sample with or without replacement.",
        default=False,
    )
    parser.add_argument(
        "--seed_group",
        type=str,
        help="Group to get binary seed list to sample from.",
        default="Run Details",
    )
    parser.add_argument(
        "--seed_key",
        type=str,
        help="Key to get binary seed list to sample from.",
        default="SEED",
    )
    return parser


def parse_args(args: List[str]) -> argparse.Namespace:
    return create_parser().parse_args(args)


def main():  # pragma: no cover
    args = parse_args(sys.argv[1:])
    sample_h5(
        args.compas_h5_filepath,
        output_filepath=args.output_filepath,
        n=args.n,
        frac=args.frac,
        replace=args.replace,
        seed_group=args.seed_group,
        seed_key=args.seed_key,
    )


if __name__ == "__main__":  # pragma: no cover
    main()
