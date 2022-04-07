sfit
====

Python extension for period finding using "least-squares periodogram"
with simultaneous linear external parameter decorrelation.  This is a
rewrite done for a student project based on a much older (roughly 2004
vintage) C program used for my thesis and some publications prior to
2014.

The computationally intensive part of the calculation (grid of
sin/cos) is done in parallel using pthreads.

Installation procedure should be standard setup.py as for most Python
modules.  Mac should be supported as well as Linux but I don't run a
Mac any more so it's not tested regularly.

We use numpy's internal LAPACK and while this usually works I've had
trouble with some peculiar numpy installs built without LAPACK or on
systems such as Debian-based distributions where numpy was installed
from a package but LAPACK isn't available.  The latter problem might
be fixed by installing the appropriate LAPACK packages using a command
like "apt-get install liblapack-dev libatlas-base-dev".

It's not currently supported by the build system but it is possible to
use either my QR routines from "lib" also on this github account (by
removing -DUSE_LAPACK and altering the include and linker paths), or
another LAPACK if desired by playing with the build flags.

A simple demo script test.py is included for processing MEarth release
light curves.  Common mode and FWHM are decorrelated by default as
well as handling meridian flips and multiple light curves.  I usually
don't decorrelate FWHM, it's just there for demonstration purposes.
The script and plotting are very crude, apologies, this was the first
thing I wrote in Python.

