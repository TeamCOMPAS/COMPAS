import codecs
import os
import re
import sys

from setuptools import find_packages, setup

def find_packages_in_src(src_dir, package_name):
    source_package_regex = re.compile(rf'^{src_dir}')
    source_packages = find_packages(include=[src_dir, f'{src_dir}.*'])
    proj_packages = [source_package_regex.sub(package_name, name) for name in source_packages]
    return proj_packages


python_version = sys.version_info
if python_version < (3, 8):
    sys.exit("Python < 3.8 is not supported, aborting setup")

NAME = "compas_python_utils"
PY_SRC = "python_src"
PACKAGES = find_packages(NAME)
HERE = os.path.dirname(os.path.realpath(__file__))
META_PATH = os.path.join(PY_SRC, "__init__.py")
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
]
EXTRA_REQUIRE = dict(dev=[
    "pytest>=3.6",
    "pre-commit",
    "flake8",
    "black==22.10.0",
    "isort"
])


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
    print("find_packages_in_src: ", find_packages_in_src(PY_SRC, NAME))
    print(f"find_packages {PY_SRC}: ", find_packages(PY_SRC))
    print(f"find_packages {NAME}: ", find_packages(NAME))

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
        package_dir={
            NAME: PY_SRC,
        },
        package_data={},
        include_package_data=True,
        install_requires=INSTALL_REQUIRES,
        extras_require=EXTRA_REQUIRE,
        classifiers=CLASSIFIERS,
        zip_safe=True,
        entry_points={
            "console_scripts": [
                f"compas_h5view= {NAME}.h5view:main",
            ]
        },
    )