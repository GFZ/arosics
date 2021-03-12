.. highlight:: shell

============
Contributing
============

Contributions are welcome, and they are greatly appreciated! Every
little bit helps, and credit will always be given.

You can contribute in many ways:

Types of Contributions
----------------------

Report Bugs
~~~~~~~~~~~

Report bugs at https://git.gfz-potsdam.de/danschef/arosics/issues

If you are reporting a bug, please include:

* Your operating system name and version.
* Any details about your local setup that might be helpful in troubleshooting.
* Detailed steps to reproduce the bug.

Fix Bugs
~~~~~~~~

Look through the GitHub issues for bugs. Anything tagged with "bug"
and "help wanted" is open to whoever wants to implement it.

Implement Features
~~~~~~~~~~~~~~~~~~

Look through the GitHub issues for features. Anything tagged with "enhancement"
and "help wanted" is open to whoever wants to implement it.

Write Documentation
~~~~~~~~~~~~~~~~~~~

arosics could always use more documentation, whether as part of the
official arosics docs, in docstrings, or even on the web in blog posts,
articles, and such.

Submit Feedback
~~~~~~~~~~~~~~~

The best way to send feedback is to file an issue at https://git.gfz-potsdam.de/danschef/arosics/issues.

If you are proposing a feature:

* Explain in detail how it would work.
* Keep the scope as narrow as possible, to make it easier to implement.
* Remember that this is a volunteer-driven project, and that contributions
  are welcome :)

You may also join our chat here: |Gitter|

.. |Gitter| image:: https://badges.gitter.im/Join%20Chat.svg
    :target: https://gitter.im/arosics/Lobby?utm_source=share-link&utm_medium=link&utm_campaign=share-link
    :alt: https://gitter.im/arosics/Lobby?utm_source=share-link&utm_medium=link&utm_campaign=share-link

Get Started!
------------

Ready to contribute? Here's how to set up `arosics` for local development.

#. Fork the `arosics` repo on GitHub.

#. Clone your fork locally:

   .. code-block:: bash

      $ git clone https://git.gfz-potsdam.de/danschef/arosics.git

#. Install your local copy into a virtualenv. Assuming you have virtualenvwrapper installed,
   this is how you set up your fork for local development:

   .. code-block:: bash

      $ mkvirtualenv arosics
      $ cd arosics/
      $ python setup.py develop

#. Create a branch for local development:

   .. code-block:: bash

      $ git checkout -b name-of-your-bugfix-or-feature

   Now you can make your changes locally.

#. When you're done making changes, check that your changes pass flake8 and the tests,
   including testing other Python versions with tox:

   .. code-block:: bash

      $ flake8 arosics tests
      $ python -m unittest discover
      $ tox

   To get flake8 and tox, just pip install them into your virtualenv.

#. Commit your changes and push your branch to GitHub:

   .. code-block:: bash

      $ git add .
      $ git commit -m "Your detailed description of your changes."
      $ git push origin name-of-your-bugfix-or-feature

#. Submit a pull request through the GitHub website.

Pull Request Guidelines
-----------------------

Before you submit a pull request, check that it meets these guidelines:

1. The pull request should include tests.
2. If the pull request adds functionality, the docs should be updated. Put
   your new functionality into a function with a docstring, and add the
   feature to the list in README.rst.
3. The pull request should work for Python 3.6+, and for PyPy. Check
   https://travis-ci.org/danschef/arosics/pull_requests
   and make sure that the tests pass for all supported Python versions.

Tips
----

To run a subset of tests:

.. code-block:: bash

    # e.g., to test if the COREG class can be properly initialized:
    $ python -m unittest tests.test_COREG.COREG_GLOBAL_init
