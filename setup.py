import codecs
import os
import re
import sys

from setuptools import find_packages, setup

python_version = sys.version_info
if python_version < (3, 8):
    sys.exit("Python < 3.8 is not supported, aborting setup")

NAME = "compas_python_utils"
PACKAGES = find_packages()
HERE = os.path.dirname(os.path.realpath(__file__))
META_PATH = os.path.join(NAME, "__init__.py")
CPP_VERSION_FILE = os.path.join("src", "changelog.h")
CLASSIFIERS = [
    "Development Status :: 5 - Production/Stable",
    "Intended Audience :: Developers",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: MIT License",
    "Operating System :: OS Independent",
    "Programming Language :: Python",
    "Programming Language :: Python :: 3",
]
INSTALL_REQUIRES = [
    "numpy",
    "h5py",
    "argparse",
    "stroopwafel",
    "pytest>=3.6",
    "pre-commit",
    "flake8",
    "black==22.10.0",
    "isort",
    "matplotlib",
    "pandas",
    "astropy",
    "scipy",
    "latex",
    "PyYAML",
    "tqdm",
    "corner"
]
EXTRA_REQUIRE = dict(
    docs=[
        "sphinx",
        "numpydoc",
        "nbsphinx",
        "sphinx_rtd_theme",
        "sphinx-tabs",
        "sphinx-argparse",
        "sphinx-math-dollar",
        "sphinxcontrib.bibtex",
        "linuxdoc",
        "ipython"
    ],
    dev=[
        "pytest-cov",
        "pre-commit",
        "flake8",
        "black==22.10.0",
        "isort",
        "coverage-badge",
        "deepdiff",
        "jupytext",
        "jupyter-autotime",
        "memory_profiler",
        "nbconvert",
        "ipykernel",
    ],
    gpu=["cupy"],
)


def read(*parts):
    with codecs.open(os.path.join(HERE, *parts), "rb", "utf-8") as f:
        return f.read()


def find_meta(meta, meta_file=read(META_PATH)):
    meta_match = re.search(
        r"^__{meta}__ = ['\"]([^'\"]*)['\"]".format(meta=meta), meta_file, re.M
    )
    if meta_match:
        return meta_match.group(1)
    raise RuntimeError("Unable to find __{meta}__ string.".format(meta=meta))


def find_version(version_file=read(CPP_VERSION_FILE)):
    version_match = re.search(
        r"VERSION_STRING = ['\"]([^'\"]*)['\"]", version_file, re.M
    )
    if version_match:
        return version_match.group(1)
    raise RuntimeError("Unable to find version string.")


if __name__ == "__main__":
    setup(
        name=NAME,
        version=find_meta("version"),
        author=find_meta("author"),
        author_email=find_meta("email"),
        maintainer=find_meta("author"),
        maintainer_email=find_meta("email"),
        url=find_meta("uri"),
        license=find_meta("license"),
        description=find_meta("description"),
        long_description=read("README.md"),
        long_description_content_type="text/markdown",
        packages=PACKAGES,
        package_data={
            f"{NAME}.preprocessing": ["*.txt", "*.yaml"],
            f"{NAME}.detailed_evolution_plotter": ["van_den_heuvel_figures/*"],
            f"{NAME}.cosmic_integration": ["SNR_Grid*"],
        },
        include_package_data=True,
        install_requires=INSTALL_REQUIRES,
        extras_require=EXTRA_REQUIRE,
        classifiers=CLASSIFIERS,
        zip_safe=True,
        entry_points={
            "console_scripts": [
                f"compas_h5view= {NAME}.h5view:main",
                f"compas_h5copy= {NAME}.h5copy:main",
                f"compas_h5sample= {NAME}.h5sample:main",
                f"compas_plot_detailed_evolution={NAME}.detailed_evolution_plotter.plot_detailed_evolution:main",
                f"compas_run_submit={NAME}.preprocessing.runSubmit:main",
                f"compas_sample_stroopwafel={NAME}.preprocessing.stroopwafelInterface:main",
                f"compas_sample_moe_di_stefano={NAME}.preprocessing.sampleMoeDiStefano:main",
                f"compas_fast_cosmic_integrator={NAME}.cosmic_integration.FastCosmicIntegration:main",
            ]
        },
    )
