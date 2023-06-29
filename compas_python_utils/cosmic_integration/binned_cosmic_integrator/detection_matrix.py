import numpy as np
import os
from typing import Dict
import h5py as h5

from .bbh_population import BBHPopulation
from .cosmological_model import CosmologicalModel
from .snr_grid import SNRGrid
from .gpu_utils import xp
from .detection_rate_computer import compute_binned_detection_rates
from .io import recursively_load_dict_contents_from_group, recursively_save_dict_contents_to_group
from .plotting import plot_detection_rate_matrix
from .bin_2d_data import bin_2d_data


class DetectionMatrix:
    def __init__(
            self,
            compas_path: str,
            cosmological_parameters: Dict,
            rate_matrix: np.ndarray,
            chirp_mass_bins: np.array,
            redshift_bins: np.array,
            outdir: str
    ):
        self.compas_path = compas_path
        self.cosmological_parameters = cosmological_parameters
        self.rate_matrix = rate_matrix
        self.chirp_mass_bins = chirp_mass_bins
        self.redshift_bins = redshift_bins
        self.outdir = outdir
        os.makedirs(self.outdir, exist_ok=True)

    @classmethod
    def from_compas_output(
            cls,
            compas_path: str,
            cosmological_parameters: Dict = dict(aSF=0.01, dSF=4.70, mu_z=-.23, sigma_z=0),
            max_detectable_redshift: float = 1.0,
            chirp_mass_bins: int = None,
            redshift_bins: int = None,
            outdir: str = None,
            save_plots: bool = False,
    ) -> "DetectionMatrix":
        outdir = 'detection_matrix' if outdir is None else outdir
        bbh_population = BBHPopulation(compas_path)
        cosmological_model = CosmologicalModel(**cosmological_parameters)
        snr_grid = SNRGrid()

        sorted_idx = xp.argsort(bbh_population.chirp_mass)
        redshift = cosmological_model.redshift

        if chirp_mass_bins is None:
            chirp_mass_bins = bbh_population.chirp_mass[sorted_idx]
        elif isinstance(chirp_mass_bins, int):
            chirp_mass_bins = xp.linspace(3, 40, chirp_mass_bins)

        if redshift_bins is None:
            redshift_bins = redshift[redshift < max_detectable_redshift]
        elif isinstance(redshift_bins, int):
            redshift_bins = xp.linspace(0, max_detectable_redshift, redshift_bins)

        rate_matrix = compute_binned_detection_rates(
            bbh_population, cosmological_model, snr_grid,
            max_detectable_redshift=max_detectable_redshift,
            chirp_mass_bins=chirp_mass_bins,
            redshift_bins=redshift_bins,
        )

        detection_matrix = cls(
            compas_path=os.path.abspath(compas_path),
            cosmological_parameters=cosmological_parameters,
            rate_matrix=rate_matrix,
            chirp_mass_bins=chirp_mass_bins,
            redshift_bins=redshift_bins,
            outdir=outdir
        )

        if save_plots:
            detection_matrix.plot().savefig(f"{outdir}/plot_{detection_matrix.label}.png")
            cosmological_model.plot().savefig(f"{outdir}/plot_{cosmological_model.label}.png")
            snr_grid.plot().savefig(f"{outdir}/plot_{snr_grid.label}.png")
            bbh_population.plot().savefig(f"{outdir}/plot_{bbh_population.label}.png")
        return detection_matrix

    @classmethod
    def from_h5(cls, path) -> "DetectionMatrix":
        with h5.File(path, "r") as f:
            data = recursively_load_dict_contents_from_group(f, '/')
        return cls(**data, outdir=os.path.dirname(path))

    def save(self):
        with h5.File(f"{self.outdir}/{self.label}.h5", "w") as f:
            recursively_save_dict_contents_to_group(f, '/', self.to_dict())

    def to_dict(self) -> Dict:
        return dict(
            compas_path=self.compas_path,
            cosmological_parameters=self.cosmological_parameters,
            rate_matrix=self.rate_matrix,
            chirp_mass_bins=self.chirp_mass_bins,
            redshift_bins=self.redshift_bins,
        )

    @property
    def label(self):
        compas_fname = os.path.basename(self.compas_path).split(".")[0]
        return f"detmatrix_{compas_fname}_{self.param_str}"

    @property
    def param_str(self):
        return "_".join([f"{k}_{v:.4f}" for k, v in self.cosmological_parameters.items()])

    def plot(self):
        return plot_detection_rate_matrix(self.rate_matrix, self.chirp_mass_bins, self.redshift_bins)

    def bin_data(self, mc_bins=50, z_bins=100):
        """Allows users to bin data in post-processing"""
        binned_data = self.rate_matrix.copy()
        mc = self.chirp_mass_bins
        z = self.redshift_bins

        # new bins
        if isinstance(mc_bins, int):
            mc_bins = np.linspace(3, 40, mc_bins)
        if isinstance(z_bins, int):
            z_bins = np.linspace(0, 1, z_bins)

        # bin rate data
        binned_data = bin_2d_data(binned_data, mc, mc_bins, axis=0)
        binned_data = bin_2d_data(binned_data, z, z_bins, axis=1)

        # get bin centers
        # mc_bins = 0.5 * (mc_bins[1:] + mc_bins[:-1])
        # z_bins = 0.5 * (z_bins[1:] + z_bins[:-1])

        self.rate_matrix = binned_data
        self.chirp_mass_bins = mc_bins
        self.redshift_bins = z_bins
