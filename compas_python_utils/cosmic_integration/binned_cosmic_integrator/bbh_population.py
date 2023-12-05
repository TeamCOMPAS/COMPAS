import numpy as np
import h5py

from .gpu_utils import xp
import h5py as h5
from typing import List, Optional

from ..totalMassEvolvedPerZ import (
    analytical_star_forming_mass_per_binary_using_kroupa_imf,
    draw_samples_from_kroupa_imf)
from .conversions import m1_m2_to_chirp_mass, m1_m2_to_eta
from .plotting import plot_bbh_population

# Default values for BBH population
M1_MIN = 5
M1_MAX = 150
M2_MIN = 0.1


class BBHPopulation(object):

    def __init__(
            self,
            m1,
            m2,
            t_delay,
            z_zams,
            n_systems,
            m1_min=M1_MIN,
            m1_max=M1_MAX,
            m2_min=M2_MIN,
            binary_fraction=0.7,
    ):
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
        self.m1 = m1
        self.m2 = m2
        self.t_delay = t_delay
        self.z_zams = z_zams
        self.n_systems = n_systems



    @classmethod
    def from_compas_h5(
            cls,
            path=None,
            m1_min=M1_MIN,
            m1_max=M1_MAX,
            m2_min=M2_MIN,
            binary_fraction=0.7,
    ):
        """
        Loads the BBH population from COMPAS output.
        Requires the COMPAS output to contain the following datasets:
        Size N (number of systems)
        - BSE_System_Parameters/SEED
        - BSE_System_Parameters/Metallicity@ZAMS(1)

        Size N_CE (number of CE events) < N
        - BSE_Common_Envelopes/SEED
        - BSE_Common_Envelopes/Immediate_RLOF>CE
        - BSE_Common_Envelopes/Optimistic_CE

        Size N_BBH (number of BBHs) < N
        - BSE_Double_Compact_Objects/SEED
        - BSE_Double_Compact_Objects/Mass(1)
        - BSE_Double_Compact_Objects/Mass(2)
        - BSE_Double_Compact_Objects/Time
        - BSE_Double_Compact_Objects/Coalescence_Time
        - BSE_Double_Compact_Objects/Stellar_Type(1)
        - BSE_Double_Compact_Objects/Stellar_Type(2)
        - BSE_Double_Compact_Objects/Merges_Hubble_Time

        """
        mask = BBHPopulation.__generate_mask(path)
        m1, m2, t_form, t_merge, dco_seeds = _load_data(
            path, "BSE_Double_Compact_Objects",
            ["Mass(1)", "Mass(2)", "Time", "Coalescence_Time", "SEED"],
            mask=mask
        )
        all_seeds, z_zams = _load_data(path, "BSE_System_Parameters", ["SEED", "Metallicity@ZAMS(1)"])
        dco_mask = xp.in1d(all_seeds, dco_seeds)
        return cls(
            m1_min = m1_min,
            m1_max = m1_max,
            m2_min = m2_min,
            binary_fraction = binary_fraction,
            m1 = m1,
            m2 = m2,
            t_delay = t_form + t_merge,
            z_zams = z_zams[dco_mask],
            n_systems = len(all_seeds),
        )

    @property
    def avg_sf_mass_needed(self):
        return self.mass_evolved_per_binary * self.n_systems

    @property
    def chirp_mass(self):
        if not hasattr(self, "_chirp_mass"):
            self._chirp_mass = m1_m2_to_chirp_mass(self.m1, self.m2)
        return self._chirp_mass

    @property
    def eta(self):
        if not hasattr(self, "_eta"):
            self._eta = m1_m2_to_eta(self.m1, self.m2)
        return self._eta

    def __repr__(self):
        return f"<BBHPopulation ({self.n_bbh:,} BBH /{self.n_systems:,} systems)>"

    def __str__(self):
        return self.__repr__()

    @staticmethod
    def __generate_mask(path):
        type_1, type_2, hubble_mask, dco_seeds = _load_data(
            path, "BSE_Double_Compact_Objects",
            ["Stellar_Type(1)", "Stellar_Type(2)", "Merges_Hubble_Time", "SEED"]
        )
        hubble_mask = hubble_mask.astype(bool)
        bbh_mask = (type_1 == 14) & (type_2 == 14)
        del type_1, type_2

        # get the flags and unique seeds from the Common Envelopes file
        ce_seeds, rlof_flag, optimistic_ce = _load_data(
            path, "BSE_Common_Envelopes", ["SEED", "Immediate_RLOF>CE", "Optimistic_CE"])
        dco_from_ce = xp.in1d(ce_seeds, dco_seeds)
        dco_ce_seeds = ce_seeds[dco_from_ce]
        del ce_seeds

        # mask out all DCOs that have RLOF after CE
        rlof_flag = rlof_flag[dco_from_ce].astype(bool)
        rlof_seeds = xp.unique(dco_ce_seeds[rlof_flag])
        mask_out_with_rlof_seeds = xp.logical_not(xp.in1d(dco_seeds, rlof_seeds))
        del rlof_flag, rlof_seeds

        # mask out all DCOs that have an "optimistic CE"
        optimistic_ce_flag = optimistic_ce[dco_from_ce].astype(bool)
        optimistic_ce_seeds = xp.unique(dco_ce_seeds[optimistic_ce_flag])
        mask_out_optimistic_ce_seeds = xp.logical_not(xp.in1d(dco_seeds, optimistic_ce_seeds))
        del optimistic_ce_flag, optimistic_ce_seeds

        return bbh_mask * hubble_mask * mask_out_with_rlof_seeds * mask_out_optimistic_ce_seeds

    @property
    def label(self):
        return f"n{self.n_bbh}_bbh_population"

    @property
    def n_bbh(self):
        return len(self.m1)

    def plot(self):
        return plot_bbh_population(
            data=xp.asarray([self.m1, self.m2, self.chirp_mass, np.log(self.z_zams), np.log(self.t_delay)]).T,
            params=[
                r"$m_1\ [M_{\odot}]$", r"$m_2\ [M_{\odot}]$", r"$\mathcal{M}_{\rm chirp}\ [M_{\odot}]$",
                r"$\ln z_{\rm ZAMS}$", r"$\ln t_{\rm delay}\ [\ln {\rm Myr}]$",
            ])


    def bootstrap_population(self):
        """Artificially generate a new population by drawing from the original population with replacement"""
        n_bbh = np.random.poisson(self.n_bbh)
        idx = np.random.choice(np.arange(self.n_bbh), size=n_bbh, replace=True)
        return BBHPopulation(
            m1=self.m1[idx],
            m2=self.m2[idx],
            t_delay=self.t_delay[idx],
            z_zams=self.z_zams[idx],
            n_systems=self.n_systems,
            m1_min=self.m1_min,
            m1_max=self.m1_max,
            m2_min=self.m2_min,
            binary_fraction=self.binaryFraction,
        )


def _load_data(path: str, group: str, var_names: List[str], mask: Optional[xp.ndarray] = None):
    """Load variables from a hdf5 file"""
    with h5.File(path, "r") as compas_file:
        data = [
            compas_file[group][var_name][...].squeeze().flatten() for var_name in var_names
        ]
    if mask is not None:  # apply mask if given
        for i in range(len(data)):  # not using list comprehension to avoid copying data
            data[i] = data[i][mask]
    return data


def generate_mock_bbh_population_file(
        filename: str = "", n_systems: int = 2000, frac_bbh: float = 0.7, m1_min: float = M1_MIN,
        m1_max: float = M1_MAX, m2_min: float = M2_MIN,
):
    """Generate a mock BBH population file with the given number of systems and BBHs.

    Args:
        filename: The filename of the mock population file.
        n_systems: The number of systems in the population.
        frac_bbh: The fraction of systems that are BBHs.
        m1_min: The minimum mass of the primary star in solar masses.
        m1_max: The maximum mass of the primary star in solar masses.
        m2_min: The minimum mass of the secondary star in solar masses.

    Returns: filename
    """
    if filename == "":
        filename = f"bbh_mock_population.h5"

    m1, m2 = draw_samples_from_kroupa_imf(
        n_samples=n_systems, Mlower=m1_min, Mupper=m1_max, m2_low=m2_min)
    n_systems = len(m1)

    # draw BBH masses
    mask = np.random.uniform(size=n_systems) < frac_bbh
    n_bbh = np.sum(mask)

    SYS_PARAM = "BSE_System_Parameters"
    DCO = "BSE_Double_Compact_Objects"
    CE = "BSE_Common_Envelopes"
    with h5py.File(filename, "w") as f:
        f.create_group(SYS_PARAM)
        f.create_group(DCO)
        f.create_group(CE)
        f[SYS_PARAM].create_dataset("SEED", data=np.arange(n_systems))
        f[SYS_PARAM].create_dataset("Metallicity@ZAMS(1)", data=np.random.uniform(1e-4, 1e-2, size=n_systems))
        f[SYS_PARAM].create_dataset("Mass@ZAMS(1)", data=m1)
        f[SYS_PARAM].create_dataset("Mass@ZAMS(2)", data=m2)
        f[CE].create_dataset("SEED", data=np.arange(n_systems))
        f[CE].create_dataset("Immediate_RLOF>CE", data=np.array([False] * n_systems))
        f[CE].create_dataset("Optimistic_CE", data=np.array([False] * n_systems))
        f[DCO].create_dataset("SEED", data=np.arange(n_bbh))
        f[DCO].create_dataset("Mass(1)", data=m1[mask])
        f[DCO].create_dataset("Mass(2)", data=m2[mask])
        f[DCO].create_dataset("Time", data=np.random.uniform(4, 13.8, size=n_bbh))
        f[DCO].create_dataset("Coalescence_Time", data=np.random.uniform(0, 14000, size=n_bbh))
        f[DCO].create_dataset("Stellar_Type(1)", data=np.array([14] * n_bbh))
        f[DCO].create_dataset("Stellar_Type(2)", data=np.array([14] * n_bbh))
        f[DCO].create_dataset("Merges_Hubble_Time", data=np.array([True] * n_bbh))
    return filename
