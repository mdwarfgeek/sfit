#!/usr/bin/env python

import os
import distutils.core
import numpy.distutils.misc_util

inc = numpy.distutils.misc_util.get_numpy_include_dirs()

# To use LAPACK, add -DWITH_LAPACK to compiler args and uncomment/modify
# libraries.

mod = distutils.core.Extension("sfit",
                               extra_compile_args=["-O3",
                                                   "-pthread",
                                                   "-ffast-math",
                                                   "-DWITH_LAPACK"],
                               extra_link_args=["-pthread"],
                               include_dirs=inc,
# With LAPACK
                               library_dirs=["/usr/lib64/atlas"],
                               libraries=["lapack",
                                          "cblas", "f77blas", "atlas",
                                          "m"],
## No LAPACK
#                               libraries=["m"],
                               sources=["linsolve.c", "sfit.c",
                                        "sysinfo.c", "wrap.c"])

distutils.core.setup(name="sfit",
                     version="0.01",
                     description="Least-squares periodogram",
                     ext_modules=[mod])

