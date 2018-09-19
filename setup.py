import os
from setuptools import setup

HERE = os.path.abspath(os.path.dirname(__file__))

about = {}
with open(os.path.join(HERE, "fluxpart", "__version__.py")) as f:
    exec(f.read(), about)

LONG = (
    "Python 3 module implementing the Scanlon and Sahu (2008) procedure "
    "for partitioning eddy covariance measurements of water vapor and "
    "carbon dioxide fluxes into stomatal (transpiration, photosynthesis) "
    "and nonstomatal (evaporation, respiration) components."
)

SHORT = "Module for partitioning eddy covariance flux measurements."

setup(
    name="fluxpart",
    version=about["__version__"],
    description=SHORT,
    long_description=LONG,
    url="https://github.com/usda-ars-ussl/fluxpart",
    author="Todd Skaggs",
    author_email="todd.skaggs@ars.usda.gov",
    license="Public Domain CC0",
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering",
        "License :: CC0 1.0 Universal (CC0 1.0) Public Domain Dedication",
        "Natural Language :: English",
        "Programming Language :: Python :: 3.7",
    ],
    packages=["fluxpart"],
    install_requires=["numpy", "matplotlib", "pandas", "pywavelets", "attrs"],
    zip_safe=False,
    tests_require=["numpy", "pandas"],
    test_suite="tests",
    extras_require={"tests": ["pytest"], "docs": ["sphinx >= 1.4"]},
)
