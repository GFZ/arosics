#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = [
    'pyfftw', 'numpy', 'gdal', 'shapely', 'scikit-image', 'matplotlib', 'geopandas', 'spectral', 'pykrige',
]

setup_requirements = [
    # TODO(danschef): put setup requirements (distutils extensions, etc.) here
]

test_requirements = ['coverage']

setup(
    name='arosics',
    version='0.1.0',
    description="An Automated and Robust Open-Source Image Co-Registration Software for Multi-Sensor Satellite Data",
    long_description=readme + '\n\n' + history,
    author="Daniel Scheffler",
    author_email='daniel.scheffler@gfz-potsdam.de',
    url='https://github.com/danschef/arosics',
    packages=find_packages(include=['arosics']),
    include_package_data=True,
    install_requires=requirements,
    license="Apache Software License 2.0",
    zip_safe=False,
    keywords='arosics',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: Apache Software License',
        'Natural Language :: English',
        "Programming Language :: Python :: 2",
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
    ],
    test_suite='tests',
    tests_require=test_requirements,
    setup_requires=setup_requirements,
)
