from compas_python_utils.cosmic_integration import cosmology
from astropy import cosmology as cosmo

from astropy.cosmology import WMAP9, Planck18


def test_set_cosmology():
    assert isinstance(cosmology.DEFAULT_COSMOLOGY, cosmo.FLRW)
    assert cosmology.DEFAULT_COSMOLOGY.name == "Planck18"


def test_setting_cosmology_with_string():
    cosmology.get_cosmology("WMAP9")
    assert cosmology.COSMOLOGY[1] == "WMAP9"
    cosmology.get_cosmology("Planck18")


def test_setting_cosmology_with_astropy_object():
    cosmology.get_cosmology(WMAP9)
    assert cosmology.COSMOLOGY[1] == "WMAP9"
    cosmology.get_cosmology(Planck18)


def test_setting_cosmology_with_default():
    cosmology.get_cosmology()
    assert cosmology.COSMOLOGY[1] == cosmology.DEFAULT_COSMOLOGY.name


def test_setting_cosmology_with_flat_lambda_cdm_dict():
    cosmo_dict = dict(H0=67.7, Om0=0.3)
    cosmology.get_cosmology(cosmo_dict)
    assert cosmology.COSMOLOGY[1][:13] == "FlatLambdaCDM"


def test_setting_cosmology_with_lambda_cdm_dict():
    cosmo_dict = dict(H0=67.7, Om0=0.3, Ode0=0.7)
    cosmology.get_cosmology(cosmo_dict)
    assert cosmology.COSMOLOGY[1][:9] == "LambdaCDM"


def test_setting_cosmology_with_w_cdm_dict():
    cosmo_dict = dict(H0=67.7, Om0=0.3, Ode0=0.7, w0=-1.0)
    cosmology.get_cosmology(cosmo_dict)
    assert cosmology.COSMOLOGY[1][:4] == "wCDM"
