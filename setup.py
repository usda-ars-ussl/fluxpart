import os
from setuptools import setup


VERSION = '0.1.0-beta-2'

LONG = ("Python 3 module implementing the Scanlon and Sahu (2008) procedure "
        "for partitioning eddy covariance measurements of water vapor and "
        "carbon dioxide fluxes into stomatal (transpiration, photosynthesis) "
        "and nonstomatal (evaporation, respiration) components.")

SHORT = ("Module for partitioning eddy covariance flux measurements.")

setup(
    name='fluxpart',
    version=VERSION,
    description=SHORT,
    long_description=LONG,
    url='https://github.com/usda-ars-ussl/fluxpart',
    author='Todd Skaggs',
    author_email='todd.skaggs@ars.usda.gov',
    license='Public Domain CC0',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering',
        'License :: CC0 1.0 Universal (CC0 1.0) Public Domain Dedication',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.5',
    ],
    packages=['fluxpart'],
    install_requires=['numpy', 'matplotlib', 'pywavelets'],
    zip_safe=False,
    tests_require=['numpy'],
    test_suite='tests',
    extras_require={
        'tests': ['pytest'],
        'docs': ['sphinx >= 1.4']}
)
