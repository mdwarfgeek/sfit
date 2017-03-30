#!/usr/bin/env python

import os
import distutils.core
import numpy.distutils.misc_util
import numpy.distutils.system_info

# Headers for numpy.
inc = numpy.distutils.misc_util.get_numpy_include_dirs()

# Use numpy's LAPACK setup.  This is not ideal, we really want single
# threaded given we're calling it out of a thread, whereas this seems
# to use the pthread versions when I've tested it.
opts = numpy.distutils.system_info.get_info("lapack_opt")

# Make sure keys we need to change exist.
for key in ["define_macros",
            "extra_compile_args",
            "extra_link_args",
            "include_dirs",
            "libraries",
            "sources"]:
  if key not in opts:
    opts[key] = []

# Add our setup to the LAPACK setup retrieved from numpy.
opts["define_macros"].extend([("USE_LAPACK", None)])
opts["extra_compile_args"].extend(["-pthread"])
opts["extra_compile_args"].extend(["-O3", "-ffast-math"])
opts["extra_link_args"].extend(["-pthread"])
opts["include_dirs"].extend(inc)
opts["libraries"].extend(["m"])
opts["sources"].extend(["sfit.c", "sysinfo.c", "wrap.c"])

mod = distutils.core.Extension("sfit", **opts)

distutils.core.setup(name="sfit",
                     version="0.01",
                     description="Least-squares periodogram",
                     ext_modules=[mod])

