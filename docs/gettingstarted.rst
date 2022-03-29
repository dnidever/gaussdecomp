***************
Getting Started
***************



Introduction
============

This is the Gaussian Decomposition software described in `Nidever et al. (2008) <https://ui.adsabs.harvard.edu/abs/2008ApJ...679..432N/abstract>`_ following the algorithm from `Haud (2000) <https://ui.adsabs.harvard.edu/abs/2000A%26A...364...83H>`_.  While it was designed to be used for HI spectra, it can be used for other types of data like Halpha spectra.

The software was originally run on the Leiden-Argentine-Bonn all-sky HI survey and there are hard-coded settings that need modification for certain datasets.  My plan is to modify the code at some point in the future to be more general and to allow these values to be more configurable.  Currently, I've been making dataset-specific copies of five of the programs and modifying the settings.

However, for most cubes the defaults should be okay.

Note that I wrote this software as a first-year graduate student back in 2005 and I haven't had much time to optimize it since then.  One of my goals is to speed it up.




Overview of Python code: :ref:`Python Overview`

Overview of IDL code: :ref:`IDL Overview`



.. _Python Overview:

Python Overview
===============

There are five main modules:

 - :mod:`~gaussdecomp.driver`:  Decomposes all of the spectra in a datacube.
 - :mod:`~gaussdecomp.fitter`:  Does the actual Gaussian Decomposition.
 - :mod:`~gaussdecomp.cube`:  Contains the :class:`~gaussdecomp.cube.Cube` class for a data cube.
 - :mod:`~gaussdecomp.spectrum`:  Contains the :class:`~gaussdecomp.spectrum.Spectrum` class for a single spectrum.
 - :mod:`~gaussdecomp.utils`:  Various utility functions.

Python code.

|gaussdecomp| Classes
---------------------

Every |gaussdecomp| model can do a few important things:

 - :meth:`~gaussdecomp.cube.Cube.train`: This trains the model on a training set.
 - :meth:`call <theborg.model.Model>`: Run the emulator for a set of data/labels, e.g., ``out = model(input)``.
 - :meth:`~gaussdecomp.spectrum.Spectrum.save`: Save the model to a file.
 - :meth:`~gaussdecomp.model.Model.load`: Load a model from a file.

   

.. _IDL Overview:

IDL Overview
============

The programs are:

 - gdriver.pro: Main driver program.  There are a few hard-coded defaults that you might want to change.
 - gloadspec.pro: Loads the cube.  You might want to make some modifications to how it loads the data.
 - gincrement.pro: Increments the current position.  Depending on the data the step-size will change.
 - parcheck.pro: Checks if Gaussians have "bad" parameters.  The thresholds for "bad" might need to be tweaked.
 - setlimits.pro: Set limits on all of the Gaussian parameters (height, velocity, width). 
 - hinoise.pro: The program that calculates the noise in each spectrum.  You should set which velocity range to use.



