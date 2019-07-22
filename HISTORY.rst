=======
History
=======

0.5.0 (2017-09-19)
------------------

New features:

* Added two test cases for local co-registration and the respective test data.

* Added test cases for global co-registration

* Added test of output writer and tie point grid visualization.

* Added nosetests. Resolved some setup requirements by conda during test_arosics_install.

* PEP8 code style now checked with automatic style checkers

Fixes and improvements:

* Coverage now also working in multiprocessing.

* Replaced test data of test case INTER1 with LZW compressed GeoTIFFs to speed up testing.

* Revised docker container builder.

* Bugfix for unexpected FFTW return value that caused the matching to fail

* Added some docstrings.

* Refactored command line interface 'arosics.py' to 'arosics_cli.py' to fix import issues.

* Added usage documentation for command line interface.

* Removed pykrige from automatically installed libraries during setup. It is now optional (Fixes issue #12)

* Bugfix in connection with optional library pyfftw.

* Revised installation guidelines within README.rst, README.md and installation.rst. Added link for nosetests HTML report.

* Fixed exception in case no arguments are provided to command line interface.

* Revised error handling and added additional check for projection.

* GDAL_DATA environment variable is now handled within py_tools_ds. Updated minimal version of py_tools_ds in setup.py.

* Fixed pickling error when running COREG_LOCAL in multiprocessing under a Windows environment.

* Replaced all occurrences of "quality grid" with "tie point grid". Updated version info.


0.4.0 (2017-07-07)
------------------

New features:

* added a logo

* added auto-deploy to PyPI

* added test cases for local co-registration


Fixes and improvements:

* fixed warping issues in case only very few tie points could be identified


0.2.1 (2017-07-03)
------------------

* First release on PyPI.


0.1.0 (2017-06-15)
------------------

* Package creation.
