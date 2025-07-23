sfit
====

Python extension for period finding using "least-squares periodogram"
with simultaneous linear external parameter decorrelation.  This is a
rewrite done for a student project based on a much older (roughly 2004
vintage) C program used for my thesis and some publications prior to
2014.

The computationally intensive part of the calculation (grid of
sin/cos) is done in parallel using pthreads.  The vectorized sincos
routines from the SLEEF library can optionally also be used to
further accelerate the calculation.  This is currently supported on
x86-64 machines but other architectures might be added if I have
access to suitable hardware and the vectorization is worthwhile.

Installation procedure should be standard setup.py as for most Python
modules.  Mac should be supported as well as Linux but I don't run a
Mac any more so it's not tested regularly.  To build with SLEEF this
needs to be visible to distutils when setup.py is executed; on modern
Debian or similar systems where there is a package available, this can
be installed using a command like "apt-get install libsleef-dev".

A simple demo script test.py is included for processing MEarth release
light curves.  Common mode and FWHM are decorrelated by default as
well as handling meridian flips and multiple light curves.  I usually
don't decorrelate FWHM, it's just there for demonstration purposes.
The script and plotting are very crude, apologies, this was the first
thing I wrote in Python.

