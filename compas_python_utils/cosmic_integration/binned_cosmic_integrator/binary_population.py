import numpy as np
import h5py
from .gpu_utils import xp
from typing import List, Optional

from ..totalMassEvolvedPerZ import (
    analytical_star_forming_mass_per_binary_using_kroupa_imf,
    draw_samples_from_kroupa_imf,
)
from .conversions import m1_m2_to_chirp_mass, m1_m2_to_eta
from .plotting import plot_binary_population
from .stellar_type import BH, NS, WD

# Default IMF limits
M1_MIN = 5
M1_MAX = 150
M2_MIN = 0.1

DCO_GROUPS = dict(
    BBH=[BH, BH],
    BNS=[NS, NS],
    BWD=[WD, WD],
    NSBH=[NS, BH],
    WDNS=[WD, NS],
    WDBH=[WD, BH],
)

VALID_DCO_TYPES = list(DCO_GROUPS.keys())


class BinaryPopulation:
    """
    General DCO population class supporting BBH, BNS, NSBH.

    Loads the DCOs population from COMPAS output.
    Requires the COMPAS output to contain the following datasets:
    Size N (number of systems)
    - BSE_System_Parameters/SEED
    - BSE_System_Parameters/Metallicity@ZAMS(1)

    Size N_CE (number of CE events) <= N
    - BSE_Common_Envelopes/SEED
    - BSE_Common_Envelopes/Immediate_RLOF>CE
    - BSE_Common_Envelopes/Optimistic_CE

    Size N_DCOs (number of N_DCOs) <= N
    - BSE_Double_Compact_Objects/SEED
    - BSE_Double_Compact_Objects/Mass(1)
    - BSE_Double_Compact_Objects/Mass(2)
    - BSE_Double_Compact_Objects/Time
    - BSE_Double_Compact_Objects/Coalescence_Time
    - BSE_Double_Compact_Objects/Stellar_Type(1)
    - BSE_Double_Compact_Objects/Stellar_Type(2)
    - BSE_Double_Compact_Objects/Merges_Hubble_Time

    """

    def __init__(
            self,
            m1: np.ndarray,
            m2: np.ndarray,
            t_delay: np.ndarray,
            z_zams: np.ndarray,
            n_systems: int,
            dcos_included: List[str],
            m1_min: float = M1_MIN,
            m1_max: float = M1_MAX,
            m2_min: float = M2_MIN,
            binary_fraction: float = 0.7,
    ):
        # Population selection
        self.dcos_included = dcos_included
        if not any([x in VALID_DCO_TYPES for x in dcos_included]):
            raise ValueError(
                f"Invalid DCO types: {dcos_included}. "
                f"Valid types are: {VALID_DCO_TYPES}"
            )

        # IMF parameters
        self.m1_min = m1_min
        self.m1_max = m1_max
        self.m2_min = m2_min
        self.binary_fraction = binary_fraction
        self.mass_evolved_per_binary = analytical_star_forming_mass_per_binary_using_kroupa_imf(
            m1_max=self.m1_max,
            m1_min=self.m1_min,
            m2_min=self.m2_min,
            fbin=self.binary_fraction,
        )

        # Data arrays
        self.m1 = m1
        self.m2 = m2
        self.t_delay = t_delay
        self.z_zams = z_zams
        self.n_systems = n_systems

    @classmethod
    def from_compas_h5(
            cls,
            path: str,
            dcos_included: List[str] = ["BBH"],
            m1_min: float = M1_MIN,
            m1_max: float = M1_MAX,
            m2_min: float = M2_MIN,
            binary_fraction: float = 0.7,
    ) -> "BinaryPopulation":
        mask = cls._generate_mask(path, dcos_included)

        m1, m2, t0, tC, seeds = _load_data(
            path,
            "BSE_Double_Compact_Objects",
            ["Mass(1)", "Mass(2)", "Time", "Coalescence_Time", "SEED"],
            mask=mask,
        )
        seeds = seeds.flatten()

        all_seeds, z_zams = _load_data(
            path,
            "BSE_System_Parameters",
            ["SEED", "Metallicity@ZAMS(1)"],
        )
        dco_mask = xp.in1d(all_seeds, seeds)

        return cls(
            m1=m1,
            m2=m2,
            t_delay=t0 + tC,
            z_zams=z_zams[dco_mask],
            n_systems=len(all_seeds),
            dcos_included=dcos_included,
            m1_min=m1_min,
            m1_max=m1_max,
            m2_min=m2_min,
            binary_fraction=binary_fraction,
        )

    @staticmethod
    def _generate_mask(
            path: str,
            dcos_included: List[str],
    ) -> xp.ndarray:
        type_mask = _generate_dco_mask(path, dcos_included)

        # Load the Hubble time flag + BSE seeds
        hubble_flag, dco_seeds = _load_data(
            path,
            "BSE_Double_Compact_Objects",
            ["Merges_Hubble_Time", "SEED"],
        )

        # Hubble time filter
        hubble_mask = hubble_flag.astype(bool)

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

        lens = dict(
            type=len(type_mask),
            hubble=len(hubble_mask),
            rlof=len(mask_out_with_rlof_seeds),
            ce=len(mask_out_optimistic_ce_seeds)
        )
        # assert all lens are equal
        assert len(set(lens.values())) == 1, f"Length mismatch in masks: {lens}"

        return type_mask * hubble_mask * mask_out_with_rlof_seeds * mask_out_optimistic_ce_seeds

    @property
    def n_dcos(self) -> int:
        return len(self.m1)

    @property
    def chirp_mass(self) -> np.ndarray:
        if not hasattr(self, "_chirp_mass"):
            self._chirp_mass = m1_m2_to_chirp_mass(self.m1, self.m2)
        return self._chirp_mass

    @property
    def eta(self) -> np.ndarray:
        if not hasattr(self, "_eta"):
            self._eta = m1_m2_to_eta(self.m1, self.m2)
        return self._eta

    @property
    def avg_sf_mass_needed(self) -> float:
        return self.mass_evolved_per_binary * self.n_systems

    def plot(self):
        arr = xp.asarray([
            self.m1,
            self.m2,
            self.chirp_mass,
            np.log(self.z_zams),
            np.log(self.t_delay),
        ]).T
        labels = [
            r"$m_1\ [M_{\odot}]$",
            r"$m_2\ [M_{\odot}]$",
            r"$\mathcal{M}_{\rm chirp}\ [M_{\odot}]$",
            r"$\ln z_{\rm ZAMS}$",
            r"$\ln t_{\rm delay}\ [\ln {\rm Myr}]$",
        ]
        return plot_binary_population(data=arr, params=labels)

    def bootstrap_population(self) -> "BinaryPopulation":
        n = np.random.poisson(self.n_dcos)
        idx = np.random.choice(self.n_dcos, size=n, replace=True)
        return BinaryPopulation(
            m1=self.m1[idx],
            m2=self.m2[idx],
            t_delay=self.t_delay[idx],
            z_zams=self.z_zams[idx],
            n_systems=self.n_systems,
            dcos_included=self.dcos_included,
            m1_min=self.m1_min,
            m1_max=self.m1_max,
            m2_min=self.m2_min,
            binary_fraction=self.binary_fraction,
        )

    @property
    def label(self):
        return f"n{self.n_dcos}_dco_population"

    def __repr__(self):
        return f"<BinaryPopulation ({self.n_dcos:,} DCOs /{self.n_systems:,} systems)>"

    def __str__(self):
        return self.__repr__()


def _generate_dco_mask(
        compas_path: str,
        dcos_included: List[str]
) -> xp.ndarray:
    # Load fundamental DCO variables
    t1, t2 = _load_data(
        compas_path,
        "BSE_Double_Compact_Objects",
        ["Stellar_Type(1)", "Stellar_Type(2)"],
    )

    masks = []
    for dco in dcos_included:
        if dco in DCO_GROUPS:
            a = [type.value for type in DCO_GROUPS[dco][0]]
            b = [type.value for type in DCO_GROUPS[dco][1]]
            # check if t1 in set 'a' of stellar types and t2 in set 'b' of stellar types
            masks.append((np.isin(t1, a) & np.isin(t2, b)) | (np.isin(t1, b) & np.isin(t2, a)))

    if not masks:
        raise ValueError("At least one DCO type must be included.")

    # Combine all masks using logical OR
    type_mask = np.logical_or.reduce(masks)
    return type_mask


def _load_data(path: str, group: str, var_names: List[str], mask: Optional[xp.ndarray] = None):
    with h5py.File(path, "r") as f:
        data = [f[group][v][...].squeeze().flatten() for v in var_names]
    if mask is not None:
        data = [d[mask] for d in data]
    return data


# Mock generation utility

def generate_mock_population(
        filename: str = "",
        n_systems: int = 2000,
        frac_bbh: float = 0.7,
        frac_bns: float = 0.2,
        frac_bhns: float = 0.1,
        m1_min: float = M1_MIN,
        m1_max: float = M1_MAX,
        m2_min: float = M2_MIN,
):
    if filename == "":
        filename = "dco_mock_population.h5"

    # sample masses and assign types
    m1, m2 = draw_samples_from_kroupa_imf(n_samples=n_systems, Mlower=m1_min, Mupper=m1_max, m2_low=m2_min)
    n_systems = len(m1)
    n_dcos = n_systems // 2
    n_ce = n_systems * 2
    types = np.random.choice(["BBH", "BNS", "NSBH"], size=n_dcos,
                             p=[frac_bbh, frac_bns, frac_bhns])

    # Define the type-to-mass mapping
    type_to_pair = {
        "BBH": [14, 14],
        "BNS": [13, 13],
        "NSBH": [13, 14]
    }

    # Create a 2D array by mapping each type to its corresponding mass pair
    mass_pairs = np.array([type_to_pair[t] for t in types]).T

    # create file structure
    with h5py.File(filename, "w") as f:
        f.create_group("BSE_System_Parameters")
        f.create_group("BSE_Common_Envelopes")
        f.create_group("BSE_Double_Compact_Objects")
        seeds = np.arange(n_systems)
        f["BSE_System_Parameters"].create_dataset("SEED", data=seeds)
        f["BSE_System_Parameters"].create_dataset("Metallicity@ZAMS(1)", data=np.random.uniform(1e-4, 1e-2, n_systems))
        f["BSE_System_Parameters"].create_dataset("Mass@ZAMS(1)", data=m1)
        f["BSE_System_Parameters"].create_dataset("Mass@ZAMS(2)", data=m2)
        # CE
        ce_seeds = np.arange(n_ce)
        f["BSE_Common_Envelopes"].create_dataset("SEED", data=ce_seeds)
        f["BSE_Common_Envelopes"].create_dataset("Immediate_RLOF>CE", data=np.zeros(n_ce, dtype=bool))  # no RLOF after CE
        f["BSE_Common_Envelopes"].create_dataset("Optimistic_CE", data=np.zeros(n_ce, dtype=bool))  # no optimistic CE
        # DCOs
        dco_seeds = np.arange(n_dcos)
        f["BSE_Double_Compact_Objects"].create_dataset("Stellar_Type(1)", data=mass_pairs[0, :])
        f["BSE_Double_Compact_Objects"].create_dataset("Stellar_Type(2)", data=mass_pairs[1, :])
        f["BSE_Double_Compact_Objects"].create_dataset("SEED", data=dco_seeds)
        f["BSE_Double_Compact_Objects"].create_dataset("Mass(1)", data=m1[:n_dcos])
        f["BSE_Double_Compact_Objects"].create_dataset("Mass(2)", data=m2[:n_dcos])
        f["BSE_Double_Compact_Objects"].create_dataset("Time", data=np.random.uniform(4, 13.8, n_dcos))
        f["BSE_Double_Compact_Objects"].create_dataset("Coalescence_Time", data=np.random.uniform(0, 14000, n_dcos))
        f["BSE_Double_Compact_Objects"].create_dataset("Merges_Hubble_Time", data=np.ones(n_dcos, dtype=bool))
    return filename
