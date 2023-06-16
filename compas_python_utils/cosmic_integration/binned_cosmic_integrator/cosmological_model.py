import matplotlib.pyplot as plt
from .gpu_utils import xp

from ..FastCosmicIntegration import find_sfr, find_metallicity_distribution, calculate_redshift_related_params

MAX_REDSHIFT = 10.0
REDSHIFT_STEP = 0.001
REDSHIFT_FIRST_STAR_FORMATION = 10.0

MIN_LOG_METALLICITY = -12.0
MAX_LOG_METALLICITY = 0.0
LOG_METALLICITY_STEP = 0.01


class CosmologicalModel:
    def __init__(
            self,
            aSF=0.01, bSF=2.77, cSF=2.90, dSF=4.70,
            mu_0=0.035, sigma_0=0.39, mu_z=-.23, sigma_z=0, alpha=0,
            min_observed_log_metallicity=-4,
            max_observed_log_metallicity=0,
    ):
        self.aSF = aSF
        self.bSF = bSF
        self.cSF = cSF
        self.dSF = dSF
        self.mu_0 = mu_0
        self.sigma_0 = sigma_0
        self.mu_z = mu_z
        self.sigma_z = sigma_z
        self.alpha = alpha

        redshifts, _, times, time_first_SF, distances, shell_volumes = calculate_redshift_related_params(
            max_redshift=MAX_REDSHIFT,
            redshift_step=REDSHIFT_STEP,
            z_first_SF=REDSHIFT_FIRST_STAR_FORMATION,
        )

        self.redshift = redshifts
        self.time = times
        self.time_first_star_formation = time_first_SF
        self.distance = distances
        self.shell_volume = shell_volumes
        self.comoving_volume = self.shell_volume / (1 + self.redshift)

        self.sfr = find_sfr(self.redshift, a=self.aSF, b=self.bSF, c=self.cSF, d=self.dSF)
        dPdlogZ, metallicities, p_draw_metallicity = find_metallicity_distribution(
            redshifts=self.redshift,
            min_logZ_COMPAS=min_observed_log_metallicity,
            max_logZ_COMPAS=max_observed_log_metallicity,
            mu0=self.mu_0,
            muz=self.mu_z,
            sigma_0=self.sigma_0,
            sigma_z=self.sigma_z,
            alpha=self.alpha,
            min_logZ=MIN_LOG_METALLICITY,
            max_logZ=MAX_LOG_METALLICITY,
            step_logZ=LOG_METALLICITY_STEP,
        )
        self.dPdlogZ = dPdlogZ
        self.metallicities = metallicities
        self.p_draw_metallicity = p_draw_metallicity

    def plot(self):
        fig, axes = plt.subplots(3, 1, figsize=(5, 8))
        ax = axes[0]
        ax.plot(self.redshift, self.sfr, label="SFR")
        ax.set_xlabel("Redshift")
        ax.set_ylabel(r"SFR [$M_{\odot}/\rm{yr}/\rm{Mpc}^3$]")
        ax.set_yscale("log")
        ax.set_xlim(0, MAX_REDSHIFT)
        textstr = self.sf_string.replace(" ", "\n")
        ax.text(0.05, 0.4, textstr, transform=ax.transAxes, fontsize=14,
                verticalalignment='top', horizontalalignment='left')

        ax = axes[1]
        dPdlogZ = self.dPdlogZ / xp.max(self.dPdlogZ)
        im = ax.imshow(dPdlogZ, extent=[0, MAX_REDSHIFT, MIN_LOG_METALLICITY, MAX_LOG_METALLICITY], aspect="auto")
        ax.set_xlabel("Redshift")
        ax.set_ylabel("logZ")
        cbar = fig.colorbar(im, ax=ax)
        cbar.set_label(r"$\frac{dP}{d\log Z}$", rotation=0, labelpad=-10)
        # set the color bar ticks to be 0 and 1 and label them min and max
        cbar.set_ticks([0, 1])
        cbar.set_ticklabels([r"min", r"max"])
        textstr = self.metallicity_string.replace(";", "\n")
        textstr = textstr.replace("N", "\mathcal{N}")
        textstr = f"${textstr}$"
        props = dict(boxstyle='round', facecolor='white', alpha=0.5)
        ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=14,
                verticalalignment='top', horizontalalignment='left', bbox=props)

        ax = axes[2]
        # plot the probability of drawing a metallicity
        ps = [self.p_draw_metallicity] * len(self.metallicities)
        ax.fill_between(self.metallicities, ps, 0, alpha=0.5)
        ax.set_ylim(0, 1)
        ax.set_xlabel("logZ")
        ax.set_ylabel("p(Z)")
        # fig.tight_layout(rect=[0, 0, 0.9, 1])
        fig.tight_layout()
        return fig

    @property
    def sf_string(self):
        return f"a={self.aSF:.3f} b={self.bSF:.3f} c={self.cSF:.3f} d={self.dSF:.3f}"

    @property
    def metallicity_string(self):
        n0 = f"N_0({self.mu_0:.3f}, {self.sigma_0:.3f})"
        nz = f"N_z({self.mu_z:.3f}, {self.sigma_z:.3f})"
        return f"{n0};{nz}"

    def __repr__(self):
        return f"CosmologicalModel({self.sf_string}, {self.metallicity_string})"
