#!/usr/bin/env python
# -*- coding: utf-8 -*-

# AROSICS - Automated and Robust Open-Source Image Co-Registration Software
#
# Copyright (C) 2017-2022
# - Daniel Scheffler (GFZ Potsdam, daniel.scheffler@gfz-potsdam.de)
# - Helmholtz Centre Potsdam - GFZ German Research Centre for Geosciences Potsdam,
#   Germany (https://www.gfz-potsdam.de/)
#
# This software was developed within the context of the GeoMultiSens project funded
# by the German Federal Ministry of Education and Research
# (project grant code: 01 IS 14 010 A-C).
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""The setup script."""

from setuptools import setup, find_packages


with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

version = {}
with open("arosics/version.py") as version_file:
    exec(version_file.read(), version)

req = [
    'cartopy',
    'cmocean',
    'dill',
    'folium>=0.6.0,!=0.12.0',
    'gdal',
    'geojson',
    'geoarray>=0.15.0',
    'geopandas',
    'matplotlib',
    'numpy',
    'packaging',
    'pandas',
    'plotly',
    'pyfftw',
    'pykrige',
    'pyproj>2.2.0',
    'py_tools_ds>=0.18.0',
    'scikit-image>=0.16.0',
    'shapely',
]

req_setup = [
    'setuptools',
    'setuptools-git'
]

req_intplot = ['holoviews', 'ipython']

req_test = ['pytest', 'pytest-cov', 'pytest-reporter-html1', 'urlchecker'] + req_intplot

req_doc = ['sphinx-argparse', 'sphinx_rtd_theme', 'sphinx-autodoc-typehints']

req_lint = ['flake8', 'pycodestyle', 'pydocstyle']

req_dev = req_setup + req_test + req_doc + req_lint

setup(
    author="Daniel Scheffler",
    author_email='daniel.scheffler@gfz-potsdam.de',
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering',
        'License :: OSI Approved :: Apache Software License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11'
    ],
    description="An Automated and Robust Open-Source Image Co-Registration Software for Multi-Sensor Satellite Data",
    entry_points={
        'console_scripts': [
            'arosics=arosics.arosics_cli:main',
        ],
    },
    extras_require={
        "interactive_plotting": req_intplot,
        "doc": req_doc,
        "test": req_test,
        "lint": req_lint,
        "dev": req_dev
    },
    include_package_data=True,
    install_requires=req,
    keywords=['arosics', 'image co-registration', 'geometric pre-processing', 'remote sensing', 'sensor fusion'],
    license="Apache-2.0",
    long_description=readme,
    name='arosics',
    packages=find_packages(exclude=['tests*']),
    project_urls={
        "Source code": "https://git.gfz-potsdam.de/danschef/arosics",
        "Issue Tracker": "https://git.gfz-potsdam.de/danschef/arosics/-/issues",
        "Documentation": "https://danschef.git-pages.gfz-potsdam.de/arosics/doc/",
        "Change log": "https://git.gfz-potsdam.de/danschef/arosics/-/blob/master/HISTORY.rst",
        "Algorithm paper": "https://www.mdpi.com/2072-4292/9/7/676",
        "Zenodo": "https://zenodo.org/record/5093940"
    },
    python_requires='>3.8',
    setup_requires=req_setup,
    test_suite='tests',
    tests_require=req + req_test,
    url='https://git.gfz-potsdam.de/danschef/arosics',
    version=version['__version__'],
    zip_safe=False,
)
