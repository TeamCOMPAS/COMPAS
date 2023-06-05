import numpy as np
import h5py as h5
from typing import List, Optional

from ..totalMassEvolvedPerZ import analytical_star_forming_mass_per_binary_using_kroupa_imf


class BBHPopulation(object):
    def __init__(
            self,
            path=None,
            m1_min=None,
            m1_max=None,
            m2_min=None,
            binary_fraction=None,
    ):
        self.path = path
        self.m1_min = m1_min
        self.m1_max = m1_max
        self.m2_min = m2_min
        self.binaryFraction = binary_fraction

        self.mass_evolved_per_binary = analytical_star_forming_mass_per_binary_using_kroupa_imf(
            m1_max=self.m1_max,
            m1_min=self.m1_min,
            m2_min=self.m2_min,
            fbin=self.binaryFraction,
        )

        mask = self.generate_mask(self.path)
        m1, m2, t_form, t_merge, dco_seeds = _load_data(
            self.path, "BSE_Double_Compact_Objects",
            ["Mass(1)", "Mass(2)", "Time", "Coalescence_Time", "SEED"], mask=mask
        )
        self.n_binaries = len(m1)
        self.t_delay = t_form + t_merge
        self.m1 = m1
        self.m2 = m2

        all_seeds, z_zams = _load_data(self.path, "BSE_System_Parameters", ["SEED", "Metallicity@ZAMS(1)"])
        dco_mask = np.in1d(all_seeds, dco_seeds)
        self.n_systems = len(all_seeds)
        self.z_zams = z_zams[dco_mask]

    @property
    def chirp_mass(self):
        if not hasattr(self, "_chirp_mass"):
            m1, m2 = self.m1, self.m2
            self._chirp_mass = (m1 * m2) ** (3 / 5) / (m1 + m2) ** (1 / 5)
        return self._chirp_mass

    @property
    def eta(self):
        if not hasattr(self, "_eta"):
            m1, m2 = self.m1, self.m2
            self._eta = m1 * m2 / (m1 + m2) ** 2
        return self._eta

    def __repr__(self):
        return f"<BBHPopulation {self.path} ({self.n_binaries:,} BBH /{self.n_systems:,} systems)>"

    def __str__(self):
        return self.__repr__()

    @staticmethod
    def generate_mask(path):

        type_1, type_2, hubble_mask, dco_seeds = _load_data(
            path, "BSE_Double_Compact_Objects",
            ["Stellar_Type(1)", "Stellar_Type(2)", "Merges_Hubble_Time", "SEED"]
        )
        hubble_mask = hubble_mask.astype(bool)
        bbh_mask = type_1 == type_2 == 14
        del type_1, type_2

        # get the flags and unique seeds from the Common Envelopes file
        ce_seeds, rlof_flag, optimistic_ce = _load_data(
            path, "BSE_Common_Envelopes", ["SEED", "Immediate_RLOF>CE", "Optimistic_CE"])
        dco_from_ce = np.in1d(ce_seeds, dco_seeds)
        dco_ce_seeds = ce_seeds[dco_from_ce]
        del ce_seeds

        # mask out all DCOs that have RLOF after CE
        rlof_flag = rlof_flag[dco_from_ce].astype(bool)
        rlof_seeds = np.unique(dco_ce_seeds[rlof_flag])
        mask_out_with_rlof_seeds = np.logical_not(np.in1d(dco_seeds, rlof_seeds))
        del rlof_flag, rlof_seeds

        # mask out all DCOs that have an "optimistic CE"
        optimistic_ce_flag = optimistic_ce[dco_from_ce].astype(bool)
        optimistic_ce_seeds = np.unique(dco_ce_seeds[optimistic_ce_flag])
        mask_out_optimistic_ce_seeds = np.logical_not(np.in1d(dco_seeds, optimistic_ce_seeds))
        del optimistic_ce_flag, optimistic_ce_seeds

        return bbh_mask * hubble_mask * mask_out_with_rlof_seeds * mask_out_optimistic_ce_seeds


def _load_data(path: str, group: str, var_names: List[str], mask: Optional[np.ndarray] = None):
    """Load variables from a hdf5 file"""
    with h5.File(path, "r") as compas_file:
        data = [
            compas_file[group][var_name][...].squeeze().flatten() for var_name in var_names
        ]
    if mask is not None:  # apply mask if given
        for i in range(len(data)):  # not using list comprehension to avoid copying data
            data[i] = data[i][mask]
    return data
