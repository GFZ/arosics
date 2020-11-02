#!/usr/bin/env python
# -*- coding: utf-8 -*-

# AROSICS - Automated and Robust Open-Source Image Co-Registration Software
#
# Copyright (C) 2017-2020  Daniel Scheffler (GFZ Potsdam, daniel.scheffler@gfz-potsdam.de)
#
# This software was developed within the context of the GeoMultiSens project funded
# by the German Federal Ministry of Education and Research
# (project grant code: 01 IS 14 010 A-C).
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License along
# with this program.  If not, see <http://www.gnu.org/licenses/>.

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
    'folium>=0.6.0',
    'gdal',
    'geojson',
    'geoarray>=0.9.0',
    'geopandas',
    'matplotlib',
    'numpy',
    'pandas',
    'plotly',
    'pyfftw',
    'pykrige',
    'py_tools_ds>=0.15.10',
    'scikit-image>=0.16.0',
    'shapely',
    'six',
]

req_setup = [
    'setuptools',
    'setuptools-git'
]

req_intplot = ['holoviews', 'ipython']

req_test = ['coverage', 'nose', 'nose2', 'nose-htmloutput', 'rednose'] + req_intplot

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
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
    ],
    description="An Automated and Robust Open-Source Image Co-Registration Software for Multi-Sensor Satellite Data",
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
    license="GPL-3.0-or-later",
    long_description=readme + '\n\n' + history,
    name='arosics',
    packages=find_packages(exclude=['tests*']),
    scripts=["bin/arosics_cli.py"],
    setup_requires=req_setup,
    test_suite='tests',
    tests_require=req + req_test,
    url='https://gitext.gfz-potsdam.de/danschef/arosics',
    version=version['__version__'],
    zip_safe=False,
)
