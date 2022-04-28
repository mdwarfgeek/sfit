#!/usr/bin/env python

import os
import platform
import subprocess
import tempfile

import distutils.core
import distutils.ccompiler
import distutils.errors
import distutils.sysconfig

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
            "extra_objects",
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
opts["sources"].extend(["sfit.c", "sysinfo.c", "vsincos_gen.c", "wrap.c"])

# Feature test for SLEEF library on x86-64.
if platform.machine() == "x86_64":
  # Set up compiler.
  compiler = distutils.ccompiler.new_compiler()
  distutils.sysconfig.customize_compiler(compiler)

  for macro in opts["define_macros"]:
    compiler.define_macro(*macro)

  # Copy compiler args and libraries.
  base_compiler_args = list(opts["extra_compile_args"])

  libraries = list(opts["libraries"])
  libraries.append("sleef")

  need_sleef_lib = False

  features = [ "sse4", "avx" ]
  compargs = [ "-msse4.1", "-mavx" ]

  for feature, comparg in zip(features, compargs):
    # Copy compiler args and add the feature we're looking for.
    compiler_args = list(base_compiler_args)
    compiler_args.append(comparg)

    srcfile = "vsincos_sleef_"+feature+".c"
    tstprog = "test_vsincos_sleef_"+feature

    feat_objs = None

    try:
      # Compile the function.
      feat_objs = compiler.compile([srcfile], extra_preargs=compiler_args)

      # Compile and link test program.
      tst_objs = compiler.compile([tstprog+".c"], extra_preargs=compiler_args)

      objs = list(feat_objs)
      objs.extend(tst_objs)

      compiler.link_executable(objs, tstprog, libraries=["sleef", "m"])

    except distutils.errors.CompileError:
      print("No SLEEF " + feature)
    except distutils.errors.LinkError:
      print("No SLEEF " + feature)
    else:
      print("Found SLEEF " + feature)
      opts["define_macros"].extend([("HAVE_SLEEF_"+feature.upper(), None)])
      opts["extra_objects"].extend(feat_objs)
      need_sleef_lib = True

  if need_sleef_lib:
    opts["libraries"].extend(["sleef"])

mod = distutils.core.Extension("sfit", **opts)

distutils.core.setup(name="sfit",
                     version="0.01",
                     description="Least-squares periodogram",
                     ext_modules=[mod])

