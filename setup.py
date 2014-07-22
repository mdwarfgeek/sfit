#!/usr/bin/env python

import os
import distutils.core
import numpy.distutils.misc_util

inc = numpy.distutils.misc_util.get_numpy_include_dirs()

mod = distutils.core.Extension("sfit",
                               extra_compile_args=["-pthread", "-ffast-math"],
                               extra_link_args=["-pthread"],
                               include_dirs=inc,
                               library_dirs=["/usr/lib64/atlas"],
                               libraries=["lapack",
                                          "cblas", "f77blas", "atlas",
                                          "m"],
                               sources=["sfit.c", "sysinfo.c", "wrap.c"])

distutils.core.setup(name="sfit",
                     version="0.01",
                     description="Least-squares periodogram",
                     ext_modules=[mod])

