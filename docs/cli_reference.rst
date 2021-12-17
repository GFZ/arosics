Command line interface reference
********************************

arosics
-------

At the command line, arosics provides the **arosics** command:

.. argparse::
   :filename: ./../arosics/arosics_cli.py
   :func: get_arosics_argparser
   :prog: arosics


.. note::

  The verbose program mode gives some more output about the interim results,
  shows some figures and writes the used footprint and overlap polygons to disk.
  Maybe the figures must be manually closed in in order to continue the processing
  (depending on your Python configuration).
