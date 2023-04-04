"""Sample an COMPAS h5 file."""
import argparse
from typing import List, Optional, Tuple, Union

import h5py
import numpy as np

from compas_python_utils.h5view import printSummary


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

    if n > len(binary_seeds) and not replace:
        raise ValueError(
            "Cannot sample without replacement more than the number of binaries. "
            "Set replace=True."
        )

    sampled_binary_seeds = np.random.choice(binary_seeds, size=n, replace=replace)

    print("Sampling {} binaries from {}.".format(n, compas_h5_filepath))

    with h5py.File(output_filepath, "w") as output_h5_file:
        with h5py.File(compas_h5_filepath, "r") as compas_h5_file:
            for group_name, group in compas_h5_file.items():
                output_h5_file.create_group(group_name)
                indicies = None
                if seed_key in group:
                    group_seeds = group[seed_key][:]
                    indicies = np.isin(group_seeds, sampled_binary_seeds)

                for dataset_name, dataset in group.items():

                    if indicies is None:  # direct copy
                        output_h5_file[group_name].create_dataset(
                            dataset_name, data=dataset[:]
                        )
                    else:  # copy the values at the sampled indicies
                        if dataset_name == seed_key:
                            output_h5_file[group_name].create_dataset(
                                dataset_name, data=dataset[:][indicies]
                            )

    print("Original file summary:")
    with h5py.File(compas_h5_filepath, "r") as compas_h5_file:
        printSummary(compas_h5_filepath, compas_h5_file)

    print("Sampled file summary:")
    with h5py.File(output_filepath, "r") as output_h5_file:
        printSummary(output_filepath, output_h5_file)


def main():
    parser = argparse.ArgumentParser(description="Sample an COMPAS h5 file.")
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
    args = parser.parse_args()
    sample_h5(
        args.compas_h5_filepath,
        output_filepath=args.output_filepath,
        n=args.n,
        frac=args.frac,
        replace=args.replace,
        seed_group=args.seed_group,
        seed_key=args.seed_key,
    )


if __name__ == "__main__":
    main()
