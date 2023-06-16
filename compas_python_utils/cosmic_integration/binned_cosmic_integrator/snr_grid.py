import matplotlib.pyplot as plt
from .conversions import chirp_mass_eta_to_m1_m2, m1_m2_to_eta_chirp_mass
from .gpu_utils import xp
from ..selection_effects import SNRinterpolator, detection_probability_from_snr

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
        self.snr_grid_at_1Mpc = SNRinterpolator(sensitivity)(self.m1, self.m2)
        self.snr = xp.arange(SNR_STEP, SNR_MAX + SNR_STEP, SNR_STEP)
        self.pdetection = detection_probability_from_snr(self.snr, SNR_THRESHOLD)

        self.Mc_step = MC_STEP
        self.eta_step = ETA_STEP
        self.snr_step = SNR_STEP

    def plot(self):
        fig, axes = plt.subplots(3, 1, figsize=(5, 9))
        ax = axes[0]
        im = ax.imshow(
            self.snr_grid_at_1Mpc,
            extent=[self.m1.min(), self.m1.max(), self.m2.min(), self.m2.max()],
            aspect="auto",
            origin="lower",
            cmap="viridis",
            vmin=SNR_THRESHOLD,
        )
        ax.set_xlabel(r"$m_1$ [$M_{\odot}$]")
        ax.set_ylabel(r"$m_2$ [$M_{\odot}$]")
        cbar = fig.colorbar(im, ax=ax)
        cbar.set_label(r"SNR @ 1 Mpc")

        ax = axes[1]
        eta, mc = m1_m2_to_eta_chirp_mass(self.m1, self.m2)
        im = ax.imshow(
            self.snr_grid_at_1Mpc,
            extent=[eta.min(), eta.max(), mc.min(), mc.max()],
            aspect="auto",
            origin="lower",
            cmap="viridis",
            vmin=SNR_THRESHOLD,
        )
        ax.set_xlabel(r"$\eta$")
        ax.set_ylabel(r"$\mathcal{M}_{\rm chirp}$ [$M_{\odot}$]")
        cbar = fig.colorbar(im, ax=ax)
        cbar.set_label(r"SNR @ 1 Mpc")

        ax = axes[2]
        ax.plot(self.snr, self.pdetection)
        ax.set_xlabel(r"SNR")
        ax.set_ylabel(r"$P({\rm detection})$")
        # determine pdet at snr threshold
        pdet_at_thresh = self.pdetection[self.snr == SNR_THRESHOLD]
        snr_at_pdet_max = self.snr[self.pdetection >= 0.985][0]
        ax.set_xlim(left=0, right=snr_at_pdet_max)
        ax.set_ylim(bottom=pdet_at_thresh, top=1)
        fig.tight_layout()
        return plt.gcf()
