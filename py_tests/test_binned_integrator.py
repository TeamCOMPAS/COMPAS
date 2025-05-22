from compas_python_utils.cosmic_integration.binned_cosmic_integrator.cosmological_model import CosmologicalModel
from compas_python_utils.cosmic_integration.binned_cosmic_integrator.binary_population import BinaryPopulation
from compas_python_utils.cosmic_integration.binned_cosmic_integrator.snr_grid import SNRGrid
from compas_python_utils.cosmic_integration.binned_cosmic_integrator.conversions import *
from compas_python_utils.cosmic_integration.binned_cosmic_integrator import DetectionMatrix
import numpy as np
import os
import glob



def test_cosmological_models(test_archive_dir):
    model = CosmologicalModel()
    assert model.dPdlogZ.shape == (len(model.redshift), len(model.metallicities))
    fig = model.plot()
    fn = os.path.join(test_archive_dir, "cosmological_model_plot.png")
    fig.savefig(fn)
    assert os.path.exists(fn)


def test_bbh_population(fake_compas_output):
    population = BinaryPopulation.from_compas_h5(fake_compas_output)
    assert population.n_systems > 2
    assert population.mass_evolved_per_binary > 0
    assert population.n_dcos > 0

def test_SNR_grid(test_archive_dir):
    snr_grid = SNRGrid()
    fig = snr_grid.plot()
    fig.show()
    fn = os.path.join(test_archive_dir, "snr_grid_plot.png")
    fig.savefig(fn)
    assert os.path.exists(fn)

def test_conversions():
    m1, m2 = 10, 3
    chirp_mass = m1_m2_to_chirp_mass(m1, m2)
    eta = m1_m2_to_eta(m1, m2)
    mtot = chirp_mass_eta_to_total_mass(chirp_mass, eta)
    assert np.isclose(mtot, m1 + m2)
    m1_new, m2_new = chirp_mass_eta_to_m1_m2(chirp_mass, eta)
    assert np.isclose(m1_new, m1)
    assert np.isclose(m2_new, m2)

def test_binned_cosmic_integration(fake_compas_output, test_archive_dir,):
    detection_matrix = DetectionMatrix.from_compas_output(
        fake_compas_output, outdir=test_archive_dir, save_plots=True,
        chirp_mass_bins=None, redshift_bins=None, n_bootstrapped_matrices=1
    )
    assert detection_matrix.rate_matrix.shape == (len(detection_matrix.chirp_mass_bins), len(detection_matrix.redshift_bins))
    detection_matrix.save()
    det_matrix_fn = glob.glob(f'{test_archive_dir}/{detection_matrix.label}.h5')[0]
    loaded_det_matrix = DetectionMatrix.from_h5(det_matrix_fn)
    assert np.allclose(detection_matrix.rate_matrix, loaded_det_matrix.rate_matrix)
    loaded_det_matrix.bin_data(mc_bins=50, z_bins=100)
    fig = loaded_det_matrix.plot()
    fig.suptitle("Binning after FastCosmicIntegrator")
    fig.savefig(os.path.join(test_archive_dir, "binned_detection_matrix_plot.png"))

    detection_matrix = DetectionMatrix.from_compas_output(
        fake_compas_output, outdir=test_archive_dir, save_plots=False,
        chirp_mass_bins=50, redshift_bins=100,
    )
    assert np.allclose(detection_matrix.rate_matrix, loaded_det_matrix.rate_matrix)
    fig = detection_matrix.plot()
    fig.suptitle("Binning during FastCosmicIntegrator")
    fig.savefig(os.path.join(test_archive_dir, "binned_detection_matrix_plot_v2.png"))








