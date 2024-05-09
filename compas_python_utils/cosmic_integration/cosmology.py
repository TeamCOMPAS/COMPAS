from typing import Union
from astropy import cosmology as cosmo

DEFAULT_COSMOLOGY = cosmo.Planck18
COSMOLOGY = [
    DEFAULT_COSMOLOGY,
    DEFAULT_COSMOLOGY.name
]
COSMO_TYPE = Union[cosmo.FLRW, str, dict]

def get_cosmology(cosmology: COSMO_TYPE = None) -> cosmo.FLRW:
    """
    Get an instance of a astropy.cosmology.FLRW subclass.

    Eg:
        cosmology.get_cosmology()
        cosmology.get_cosmology(astropy.cosmology.WMAP9)
        cosmology.get_cosmology("WMAP9")
        cosmology.get_cosmology(dict(H0=67.7, Om0=0.3, Ode0=0.7, w0=-1.0) --> wCDM
        cosmology.get_cosmology(dict(H0=67.7, Om0=0.3, Ode0=0.7) --> LambdaCDM



    Parameters
    ==========
    cosmology: astropy.cosmology.FLRW, str, dict
        Description of cosmology, one of:
            None - Use DEFAULT_COSMOLOGY
            Instance of astropy.cosmology.FLRW subclass
            String with name of known Astropy cosmology, e.g., "Planck18"
            Dictionary with arguments required to instantiate the cosmology
            class.
    """
    global COSMOLOGY


    if cosmology is None:
        cosmology = DEFAULT_COSMOLOGY

    # if already a 'astropy.cosmology' object
    elif isinstance(cosmology, cosmo.FLRW):
        cosmology = cosmology

    elif isinstance(cosmology, str):
        cosmology = getattr(cosmo, cosmology)

    elif isinstance(cosmology, dict):
        if 'Ode0' in cosmology.keys():
            if 'w0' in cosmology.keys():
                cosmology = cosmo.wCDM(**cosmology)
            else:
                cosmology = cosmo.LambdaCDM(**cosmology)
        else:
            cosmology = cosmo.FlatLambdaCDM(**cosmology)

    # cache the cosmology
    COSMOLOGY[0] = cosmology
    COSMOLOGY[1] = repr(cosmology) if not cosmology.name else cosmology.name

    return cosmology


def set_cosmology(cosmology: COSMO_TYPE = None):
    """

    """
    _set_default_cosmology()
    if cosmology is None:
        cosmology = DEFAULT_COSMOLOGY

