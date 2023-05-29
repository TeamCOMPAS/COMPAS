from compas_python_utils import cdriver
import matplotlib.pyplot as plt
import os
import numpy as np
import pandas as pd
from astropy import units as u

DATA_DIR = os.path.join(os.path.dirname(__file__), "imf_data")
AVG_STAR_FORMING_MASS_PER_BINARY_CACHE = os.path.join(DATA_DIR, "star_forming_mass_per_binary_cache.csv")


def plot_kroupa_imf_pdf(bins=30, fig=None, label="Kroupa IMF", n=int(1e6)):
    """Plot the imf pdf"""
    kroupa_masses = cdriver.sample_from_imf(n)

    if isinstance(bins, int):
        bins = np.geomspace(min(kroupa_masses), max(kroupa_masses), bins)

    kroupa_counts, mass_bins = np.histogram(kroupa_masses, bins=bins)
    kroupa_pdf = kroupa_counts / np.sum(kroupa_counts)

    if fig is None:
        fig = plt.figure()

    plt.plot(mass_bins[:-1], kroupa_pdf, label=label)
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel(r"Mass [M$_\odot$]")
    plt.ylabel("PDF")
    return fig


def get_star_forming_mass_per_binary(binaryFraction, Mupper, Mlower, m2_min):
    """Get the star forming mass per binary"""
    _round_m = lambda m: np.round((m / u.M_sun).value, 2)
    query = dict(binaryFraction=binaryFraction, Mupper=_round_m(Mupper), Mlower=_round_m(Mlower),
                 m2_min=_round_m(m2_min))
    avg_mass = _StarAvgMassPerBinaryCache.get_item(query)
    return avg_mass * u.M_sun


class _StarAvgMassPerBinaryCache:
    def __init__(self):
        self._file = AVG_STAR_FORMING_MASS_PER_BINARY_CACHE

    @property
    def exists(self):
        return os.path.exists(self._file)

    @property
    def cache(self):
        return pd.read_csv(self._file)

    def __update(self, data_dict):
        df = pd.DataFrame(data_dict, index=[0])
        if self.exists:
            df = self.cache.append(df, ignore_index=True)
        df.to_csv(self._file, index=False)

    def __query(self, query):
        """query cache, return None if not in cache, otherwise return the value"""
        if self.exists:
            result = self.cache.query(" & ".join([f"{key} == {val}" for key, val in query.items()]))
            if not result.empty:
                return result["avg_mass"].values[0]
        return None

    @classmethod
    def get_item(cls, query):
        cache = cls()
        avg_mass = cache.__query(query)
        if avg_mass is None:
            print("Computing star forming mass per binary")
            avg_mass = cdriver.compute_star_forming_mass_per_binary(**query, n=int(1e7))
            cache.__update(dict(query, avg_mass=avg_mass))
        return avg_mass
