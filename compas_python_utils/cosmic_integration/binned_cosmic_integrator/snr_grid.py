import matplotlib.pyplot as plt
from .conversions import chirp_mass_eta_to_m1_m2, m1_m2_to_eta_chirp_mass
from .gpu_utils import xp
from ..selection_effects import SNRinterpolator, detection_probability_from_snr
from .plotting import plot_snr_grid

MC_MAX = 300.0
MC_STEP = 0.1

ETA_MAX = 0.25
ETA_STEP = 0.01

SNR_THRESHOLD = 8.0
SNR_MAX = 1000.0
SNR_STEP = 0.1


class SNRGrid:
    def __init__(
            self,
            sensitivity="O1",
    ):
        self.sensitivity = sensitivity
        self.snr_threshold = SNR_THRESHOLD
        self.chirp_mass = xp.arange(MC_STEP, MC_MAX + MC_STEP, MC_STEP)
        self.eta = xp.arange(ETA_STEP, ETA_MAX + ETA_STEP, ETA_STEP)
        self.m1, self.m2 = chirp_mass_eta_to_m1_m2(*xp.meshgrid(self.chirp_mass, self.eta))
        # TODO: SNRInterpolator can be GPU-ized
        self.snr_grid_at_1Mpc = SNRinterpolator(sensitivity)(self.m1, self.m2)
        self.snr = xp.arange(SNR_STEP, SNR_MAX + SNR_STEP, SNR_STEP)
        self.pdetection = detection_probability_from_snr(self.snr, SNR_THRESHOLD)

        self.Mc_step = MC_STEP
        self.eta_step = ETA_STEP
        self.snr_step = SNR_STEP

    def plot(self):
        return plot_snr_grid(
            snr_grid_at_1Mpc=self.snr_grid_at_1Mpc,
            m1=self.m1, m2=self.m2,
            snr=self.snr, pdetection=self.pdetection,
            snr_threshold=SNR_THRESHOLD,
        )

    @property
    def label(self):
        return f"snr_{self.sensitivity}"
