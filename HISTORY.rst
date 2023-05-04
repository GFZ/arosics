=======
History
=======

1.9.0 (2023-05-04)
------------------

* Re-implemented interpolation of tie point attributes into space (!23).
  Three techniques are now supported: RBF, GPR, and Ordinary Kriging. Added scikit-learn to dependencies.
* Fixed GDAL warning related to PROJ_DATA/PROJ_LIB environment variables (!34).
* Added parameter 'rs_random_state' (RANSAC random state) and set it to 0 by default (!35).


1.8.1 (2023-03-10)
------------------

* Fixed https://github.com/GFZ/arosics/issues/20 (out_crea_options has no effect) (!32).
* Updated copyright (!33).


1.8.0 (2022-12-02)
------------------

* Replaced CentOS CI container with an Ubuntu-based one (!31).
* Removed version pinning of pyfftw. AROSICS now uses numpy instead of pyfftw>=0.13.0 (!30).
* Added Python 3.11 classifiers and dropped support for Python 3.7 due to near EOL status (!30).
* Fixed DeprecationWarnings (2x matplotlib, 1x pandas) (!30).
* Dropped deprecated 'arosics_cli.py' entry point.


1.7.9 (2022-11-21)
------------------

* Allow pyfftw=0.13.0=*0 (!28).
* Replaced deprecated URL (!28).
* Updated copyright (!29).


1.7.8 (2022-08-19)
------------------

* Fixed missing grid and frame lines of Tie_Point_Grid.plot_shift_distribution() (!27).
* Remove deprecated sphinx_autodoc_typehints option.


1.7.7 (2022-07-21)
------------------

* Added a validation to reject images that only contain nodata values.


1.7.6 (2022-01-31)
------------------

* Fixed typo in requirements.


1.7.5 (2022-01-28)
------------------

* Implemented a workaround for #71 (pyFFTW RuntimeError: Undefined plan with nthreads. This is a bug).
  Pin pyFFTW to <0.13.0. (!25)


1.7.4 (2021-12-15)
------------------

* Migrated test calls from nosetests to pytest and implemented new test report (!24).
* Removed folium 0.12.1 exclusion due to fix on PyPI and conda-forge.
* Fixed dead link.


1.7.3 (2021-12-02)
------------------

* Avoid folium 0.12.1 as requirement due to https://github.com/python-visualization/folium/issues/1513.


1.7.2 (2021-11-09)
------------------

* Listed dill in package requirements.
* Resolved inconsistency in documenting boolean default values.
* Improved error message when trying to compute statistics and all tie points are flagged as false-positives.


1.7.1 (2021-10-13)
------------------

* Added more statistics to Tie_Point_Grid.calc_overall_stats().


1.7.0 (2021-09-30)
------------------

* Added method Tie_Point_Grid.calc_overall_stats() + test to compute overall statistics from all tie points found.


1.6.2 (2021-09-29)
------------------

* Fixed 'too many values to unpack' exception in COREG_LOCAL.view_CoRegPoints().
* Added new parameters to Tie_Point_Grid.plot_shift_distribution().
* Added documentation to COREG_LOCAL.view_CoRegPoints().


1.6.1 (2021-09-29)
------------------

* The output map of COREG_LOCAL.view_CoRegPoints() is now cropped nicely. Added parameter 'figsize_multiplier'.


1.6.0 (2021-09-27)
------------------

* CI now uses Mambaforge.
* Switched to Apache 2.0 license.


1.5.1 (2021-08-11)
------------------

* Added project URLs to setup.py.
* Revised test_arosics_install CI job (now uses mamba).
* Updated minimal version of py_tools_ds which fixes #63 (Computing footprint error).


1.5.0 (2021-07-13)
------------------

* COREG.show_matchWin() now also works when no valid shift was found yet.
* Updated minimal version of geoarray to fix a sporadic TypeError when writing the coregistered result.
* Added basic compatibility with images that have a rotation in the map info (fixes #60).
* Fixed incorrect tolerance in COREG.equalize_pixGrids() which speeds up COREG_LOCAL a lot.


1.4.8 (2021-07-06)
------------------

* Updated HISTORY.rst (closes #46).
* Updated minimal version of py_tools_ds which fixes #61 (Multiprocessing issue).


1.4.7 (2021-07-03)
------------------

* Set minimal Python version to 3.6 in setup.py.
* Updated minimal version of py_tools_ds which fixes #59 (argmax error when not specifying poly WKT).
* Removed requirements_dev.txt and requirements_pip.txt as they are not needed anymore.


1.4.6 (2021-06-22)
------------------

* Updated minimal version of py_tools_ds which fixes #49 (Deadlock in Multiprocessing on Linux,
  in case GDAL 3.2.1 is installed.).


1.4.5 (2021-06-21)
------------------

* Fixed issue #56 (Setting align_grids to False in case of global co-registration does not avoid resampling.).


1.4.4 (2021-06-18)
------------------

* Removed non-operational tests.
* Apply 'make lint' to the whole package.
* Moved arosics_cli.py from bin to arosics folder to avoid entry point issues at conda-forge.
* Replaced some http links with https.
* The tie point grid now also includes the last x any y grid points (bugfix).


1.4.3 (2021-05-08)
------------------

* Clearer error message in case of unequal input projections that have equal names.
* Updated minimal version of geoarray.
* Reverted .gitlab-ci.yml.


1.4.2 (2021-05-08)
------------------

* Temporarily disabled installation check.
* Remove updates of py_tools_ds and geoarray for now.
* Updated minimal versions of geoarray and py_tools_ds. Updated installation.rst.
* Caught a warning in case the footprint polygon has disjunct polygons.
* Increased tolerance when snapping clip extent to output grid in DESHIFTER
  (fixes https://github.com/GFZ/arosics/issues/4 (SSIM fails with window shape mismatch error)).


1.4.1 (2021-04-30)
------------------

* Deprecated the 'arosics_cli.py' argument parser call and added a new entrypoint called 'arosics'.
* Temporarily increase the URLChecker timeout to 60 secs due to Zenodo latency.


1.4.0 (2021-04-23)
------------------

* Temporarily set URLChecker timeout to 30 secs due to Zenodo latency.
* Updated projection related documentation.
* Added test which covers an input projection other than UTM and LonLat and which has no corresponding EPSG code.
* Updated minimal version of geoarray to v0.11.0 which brings no-EPSG compatibility.
* Replaced all UTM specific code and refactored tie point grid table columns 'X_UTM' and 'Y_UTM' to 'X_MAP' and 'Y_MAP'.
  This allows to run arosics on input images with projections other than UTM and Lon/Lat.
* Replaces EPSG related parameters by WKT string to get rid of missing EPSG issue.
* Updated copyright.
* Improved documentation.
* Added tolerance when comparing coordinate grids with float values.
* Improved docstring for COREG parameters wp and ws.
* Fixed inaccurate type hint.
* Removed unnecessary import.
* Fix typo.


1.3.0 (2021-03-12)
------------------

* Replaced Python 2 compatible type hints by PEP 484 type hints.
* Revised docstrings.
* Fixed a lot of Sphinx build warnings. Some code style improvements.
* Dropped support for Python versions below 3.6.
* Replaced deprecated np.object type.
* Make lint now directly prints its output in case of exceptions.
* Removed deprecated coreg_cmd.py. Fixed typo.


1.2.6 (2021-02-16)
------------------

* Fixed CI job.
* Fixed pyproj DeprecationWarning related to proj4 string. Added pyproj to dependencies
  (which was already used under the hood).
* Fixed DeprecationWarning related to deprecated numpy data types that are only aliases for builtin types.
* Added type hints for COREG_LOCAL.tiepoint_grid and COREG_LOCAL.CoRegPoints_table.


1.2.5 (2021-02-02)
------------------

* Excluded folium 0.12.0 from requirements due to https://github.com/python-visualization/folium/issues/1438.
* Fixed incompletely deleted coverage artifacts after running 'make clean'.
* Fixed wrong dependency name.
* Updated URLs due to changes on the server side.
* Removed tests for issue 70 and 47.
* Commented rever CI job out.
* Added URL checker CI job and fixed all dead URLs. Removed travis related file
* Fixed issue of remaining coverage artifacts after running 'make clean-test'.


1.2.4 (2021-02-02)
------------------

* Caught the no-tie-points-found-case in some methods of Tie_Point_Grid.


1.2.3 (2020-11-13)
------------------

* Fixed KeyError 'ABS_SHIFT' in  Tie_Point_Grid.plot_shift_distribution() in case no tie points have been found at all.


1.2.2 (2020-11-13)
------------------

* Fixed issue #47 (COREG_LOCAL.view_CoRegPoints() raises KeyError: 'X_SHIFT_M' error when there are too many clouds).
* Increased default figsize of COREG_LOCAL.view_CoRegPoints().


1.2.1 (2020-11-11)
------------------

* Added 'coverage erase' to clean-test.
* Fixed issue #45 (CoReg gives ValueError: `min_samples` must be in range (0, <number-of-samples>)`).
* Replaced deprecated osgeo imports.


1.2.0 (2020-11-02)
------------------

* Fixed issue 44
  (SSIM filtering flags too much tie points in case of completely different data ranges of the input images).


1.1.1 (2020-11-02)
------------------

* Replaced deprecated osgeo imports.


1.1.0 (2020-10-30)
------------------

* Added a warning in case the input image consists of multiple patches and AROSICS processes only the largest one.
* Added a warning in case the reliability filtering filters more than 70% of the tie points.
* Fixed issue #43 (AttributeError in case COREG_LOCAL.tieP_filter_level = 0).


1.0.6 (2020-10-27)
------------------

* Updated minimal version of py_tools_ds (fixes issue #41 (Sporadic AssertionErrors in case the matching window
  crosses the image edge)).
* Revised requirements and environment_arosics.yml.
* Replaced deprecated 'source activate' by 'conda activate'. Updated installation instructions.
* Unittests are now also executable on Windows.


1.0.5 (2020-10-21)
------------------

* Added shebang to bin files to ensure they Python executable (fixes issue #16).


1.0.4 (2020-10-21)
------------------

* Fix for not passing the quiet mode parameter to Tie_Point_Refiner class when using CORE_LOCAL.


1.0.3 (2020-10-19)
------------------

* Fixed linting.
* Fixed an unhelpful error message in case no coregistration point can be placed within an image area usable for
  coregistration due to the provided bad data mask.
* Fixed some wrong type hints.
* Added COREG_LOCAL.calculate_spatial_shifts() allowing to explicitly compute the shifts instead of implicitly
  running the getter properties. This improves API clarity and facilitates debugging.
* Added sphinx-autodoc-typehints to doc requirements.


1.0.2 (2020-10-12)
------------------

* Fixed linting.
* Fixed DeprecationWarning within CORE_LOCAL.view_CoRegPoints().
* Caught matplotlib warnings within tests.
* Added test/doc/lint/dev requirements as optional installation procedures to setup.py.


1.0.1 (2020-10-12)
------------------

* Excluded tests from being installed via 'pip install'. Set development status to 'stable'.
* Use SPDX license identifier and set all files to GLP3+ to be consistent with license headers in the source files.
* Improved installation instructions. Added conda-forge badge.


1.0.0 (2020-10-06)
------------------

* Revised COREG_LOCAL.view_CoRegPoints() and replaced basemap with cartopy.
* Revised environment_arosics.yml.
* Fixed issue #36.
* Closed issue #26.


0.9.26 (2020-10-02)
-------------------

* Fixed broken pip installation of basemap.


0.9.25 (2020-09-30)
-------------------

* Replaced requirement 'basemap' in setup.py and requirements.txt by ssh link to fix exception during 'pip install'.
* Updated the installation instructions as AROSICS is now on conda-forge.


0.9.24 (2020-09-28)
-------------------

* The 'pykrige', 'pyfftw' and 'basemap' requirements are no longer optional since they are easily installable from
  conda-forge now.
* Updated requirements and installation instructions.


0.9.23 (2020-09-25)
-------------------

* Moved all matplotlib imports from module level to function/class level to avoid static TLS ImportError.


0.9.22 (2020-10-02)
-------------------

* Moved all skimage imports from module level to function/class level to avoid static TLS ImportError.


0.9.21 (2020-10-15)
-------------------

* Replaced deprecated HTTP links.
* Fixed typo.
* arosics_ci.docker now inherits from ci_base_centos:0.1. Conda update now uses conda-forge channel.
* Don't inherit from gms_base.
* Re-added conda-forge::libgdal.
* Fixed syntax.
* Added pip to requirements.
* Updated CI setup files and .gitlab-ci.yml.
* Added some information about supported projections to the docs.


0.9.20 (2020-08-26)
-------------------

* AROSICS now uses pyproj>2.2 under the hood for projection transformations.
* Added minimal version of pyproj.


0.9.19 (2020-08-21)
-------------------

* Added tolerances to the window position validation to avoid float precision issues.
* Updated minimal version of geoarray.
* Fixed a bug which causes COREG.equalize_pixGrids() to run although the pixel grids of reference and target image
  are equal.
* Fixed ResourceWarning in COREG.show_matchWin() as well as in COREG.calculate_spatial_shifts().
* Fixed exception in COREG.view_CoRegPoints_folium() in case of a recent version of folium.


0.9.18 (2020-08-18)
-------------------

* Added geoarray update to CI config.
* Fixed DeprecationWarning coming from holoviews.


0.9.17 (2020-05-19)
-------------------

* Updated minimal version of py_tools_ds (fixes PyProj DeprecationWarning).


0.9.16 (2020-05-19)
-------------------

* Fixed create_github_release CI job.


0.9.15 (2020-04-09
-------------------

* Added create_release_from_gitlab_ci.sh and updated create_github_release CI job.


0.9.14 (2020-04-08)
-------------------

* Fixed create_github_release CI job.


0.9.13 (2020-04-08)
-------------------

* Fixed invalid yaml syntax.
* Added CI job 'create_github_release'.


0.9.12 (2020-04-08)
-------------------

* Revised .gitlab-ci.yml. Updated installation instructions
  (Python is now installed from conda-forge channel - fixes issue #35).
* Updated test_arosics_install CI job.
* Added funding information.


0.9.11 (2020-04-07)
-------------------

* Fixed typo.


0.9.10 (2020-04-07)
-------------------

* Added Zenodo badge and citation hint to README.rst.


0.9.9 (2020-04-07)
------------------

* Fixed line break.


0.9.8 (2020-04-07)
------------------

* Updated .zenodo.json.
* Added CITATION file.
* Updated copyright.
* Updated installation instructions and environment_arosics.yml.
* Added .zenodo.json file.
* Removed version pinnings from requirements_dev.txt.


0.9.7 (2020-04-06)
------------------

* Fix incompatibity with shapely 1.7.0
  (implies an update of the minimal version of py_tools_ds). Remove shapely version pinning.


0.9.6 (2020-02-11)
------------------

* Pinned shapely to versions older or equal than 1.6.4.


0.9.5 (2020-01-08)
------------------

* Updated minimal version of py_tools_ds.
* Updated conda environment file.


0.9.4 (2020-01-08)
------------------

* Removed pyresample dependency (not needed anymore).
* Fixed broken badge.
* Merge branch 'bugfix/adapt_to_geopandas_changes' into 'master'


0.9.3 (2019-11-27)
------------------

* Fixed issue #31 (ValueError: Unknown column geometry).
* Fixed issue #32 (NotImplementedError: fillna currently only supports filling with a scalar geometry).
* Added pandas to requirements.
* Changed badge target.
* Added downloads badge.


0.9.2 (2019-11-27)
------------------

* Removed deprecated PyPI upload code from .gitlab-ci.yml. Replaced relative links in README.rst by absolute ones.


0.9.1 (2019-07-26)
------------------

* Added title to README.rst. Try to disable title.
* Added pyresample to conda dependencies (might fix test_arosics_install). Replaced deprecated PyPI upload by twine.
* Changed description file in setup.cfg.
* Added missing cli_reference.rst content.
* Added missing cli_reference.rst.


0.9.0 (2019-11-27)
------------------

* Removed the deprecated README.md.
* Replaced HTML table by image.
* Added links and fixed typo.
* Revised about.rst, added Gitter badge.
* Revised README.rst.
* Resized images physically.
* Updated README.rst.
* Revised CONTRIBUTING.rst
* Improved code block style.
* Changed toc maxdepth.
* Added usage instructions.
* Updated api_cli_reference.rst and sub-sections.
* Updated usage.rst and sub-sections.
* Moved CLI reference to API reference subsection.
* Fix in installation.rst.
* Revised README.rst.
* Updated usage.rst.
* Updated installation.rst.
* Enabled TODOs to be rendered.
* Revised docstring style.
* Added caption.
* Added subsections to usage.rst.
* First empty version of usage.rst.
* Revised DESHIFTER.__doc__.
* Revised about.rst.
* Revised DESHIFTER.__doc__. Added sphinx type hints.
* Added about.rst. Updated index.rst and re-ordered HISTORY.rst.
* Test revision.
* Revised 'make docs' rule.
* Revised DESHIFTER.__doc__.
* Changed sphinx theme. Documentation now also includes __init__() methods.
* Increased sphinx documentation content width.


0.8.19 (2019-07-22)
-------------------

* Removed hardcoded test.
* Added license texts. Added funding note.


0.8.18 (2019-06-17)
-------------------

* Fixed issue #30 (Exception in case of non-quadratic pixels of the input images).


0.8.17 (2019-05-10)
-------------------

* Updated minimal version of geoarray.


0.8.16 (2019-02-27)
-------------------

* Updated minimal version of py_tools_ds (fixes issue #27).


0.8.15 (2019-02-19)
-------------------

* Fixed PyPi upload error (invalid value for classifiers within setup.py).).
* Updated minimal version of py_tools_ds.
* Added tests for ETRS/LAEA projection compatibility.
* Fixed some style issues.
* Added gitter badge. Added classifiers to setup.py.
* Added keywords.
* Code style improvements.


0.8.14 (2018-12-05)
-------------------

* Moved cmocean to conda requirements due to setup issue under Python 2.7.
* Removed '-y -q' from conda install commands contained in installation instructions in README files.
* Replaced 'importlib.util.find_spec' with 'pkgutil.find_loader' to ensure Python 2.7 compatibility.
* Updated minimal version of geoarray.
* Added Python 3.7 to classifiers in setup.py.


0.8.13 (2018-12-04)
-------------------

* Test fix.
* Fixed issue # 17 (Coregistration sometimes fails in case of floating point coordinates of the input images.)
* Fixed an issue causing SSIM computation to fail (due to float coordinates).


0.8.12 (2018-11-30)
-------------------

* Fixed issue #23 ('TypeError in case COREG or COREG_LOCAL is called with an in-memory reference or target image and
  path_out is set to 'auto'.').


0.8.11 (2018-11-28)
-------------------

* Fixed exception in case Tie_Point_Grid.to_PointShapefile() is called with skip_nodata=False.


0.8.9 (2018-11-27)
------------------

* Fixed figure of tie point grid broken due to matplotlib 3.0.0/basemap 1.2 incompatibility.


0.8.8 (2018-10-22)
------------------

* Fixed issue #21(pandas value error during dataframe merging).
* Fixed linting.
* Added folium and geojson to requirements. Fixed view_CoRegPoints_folium().
* CI setup now updates ci_env environment installed via docker_pyenvs instead of creating an independent environment.
* Fixed duplicate of pycodestyle in environment file.
* Fix.
* CI Python environment is not separate from the base env.
* Fixed mixed channels for gdal and libgdal causing libkea issue during CI.


0.8.7 (2018-08-10)
------------------

* Fix for incompatible version of pycodestyle during CI.
* Updated minimally required geoarray version.
* Added version.py.
* Bugfix.
* Implemented changes from the current branch of geoarray (feature/improve_metadata_handling).
* Updated docker runner build script.


0.8.6 (2018-07-20)
------------------

* Bugfix for issue #13 (ValueError related to pandas.merge).


0.8.5 (2018-04-25)
------------------

* Fixed documentation on output data format.
* Updated test_COREG.py.


0.8.4 (2018-03-08)
------------------

* Removed TestBandnames.
* Revised previous commit.


0.8.3 (2018-03-07)
------------------

* Fixed ValueError as reported in https://git.gfz-potsdam.de/EnMAP/sicor/issues/22.


0.8.2 (2018-01-23)
------------------

* Revised arosics_cli.py.
* Fixed issue #14.
* Added importlib (must be revised).


0.8.1 (2017-11-21)
------------------

* Added test for COREG_LOCAL.view_CoRegPoints_folium().


0.8.0 (2017-11-21)
------------------

* Added shift vector plots (COREG_LOCAL.view_CoRegPoints(shapes2plot='vectors') + tests.


0.7.0 (2017-11-20)
------------------

* Adapted docker installer to new external base image.
* Updated arosics_environment.yml.
* Updated docker installer workflow.
* Added environment_arosics.yml
* Updated minimal version py_tools_ds.
* Added geopandas to CI installer test.
* Updated minimum version of py_tools_ds in docker container setup.
* Added Test_Tie_Point_Grid.tearDown().
* Removed old functions for deshifting within COREG class:
* Moved several functions to py_tools_ds.
* Removed deprecated functions.
* Removed io and utilities modules.


0.6.8 (2017-11-16)
------------------

* Fixed Tie_Point_Grid.to_PointShapefile().
* Added tests for some functions within Tie_Point_Grid.
* Updated README files.
* Updated README files and installation.rst.
* Moved package geopandas to conda dependencies.


0.6.7 (2017-11-15)
------------------

* Fixed exceptions within Tie_Point_Grid.plot_shift_distribution(), calc_overall_mssim(), calc_rmse. Added test_tie_point_grid.py.


0.6.6 (2017-10-26)
------------------

* Updated minimal version of geoarray.
* Added requirements_pip.txt.


0.6.5 (2017-11-18)
------------------

* Bugfix for not checking validity of GeoArray_CoReg.footprint_poly.


0.6.4 (2017-10-12)
------------------

* Updated minimal versions of geoarray and py_tools_ds.


0.6.3 (2017-10-12)
------------------

* Excluded some funcs from coverage.
* test_arosics_install is now executed within latest Python.
* Updated docker setup.


0.6.2 (2017-10-11)
------------------

* Fixed pages.
* Updated .gitlab-ci.yml to make pages work again.


0.6.1 (2017-10-10)
------------------

* Simplified dependency checks.


0.6.0 (2017-10-10)
------------------

* Updated docker setup.
* Updated minimal versions of dependencies.
* Disabled coverage for deprecated funcs. Too small SCPS is now catched.
* Tie_Point_Grid.get_CoRegPoints_table(): local CS not rejectd anymore.
* Fixed test_shift_calculation_with_image_coords_only(). Fixed flake8 issues.
* SSIM now fails with a warning instead of raising an exception and forcing the whole coreg to fail.
* test_COREG.test_shift_calculation_with_image_coords_only: changed input gt.
* Revised COREG.show_matchWin().
* COREG.calculate_spatial_shifts(): removed deprecated function.
* Added test_shift_calculation_with_image_coords_only()


0.5.1 (2017-10-06)
------------------

* First attempt to implement autoclip to polygon to fix unequal matching window sizes in case of float coordinates.
* Updated test_COREG_LOCAL.
* Tie_Point_Grid: added type hints.
* DeShifter: cleaned up.
* Cleaned requirements.txt.


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
* Replaced all occurrences of "quality grid" with "tie point grid".


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
