#!/usr/bin/env python

import sys
import math
import numpy
import sfit
import matplotlib.pyplot as plt

argc = len(sys.argv)

if argc < 3:
  print "Usage:\t", sys.argv[0], "lcfile [...] pmin pmax"
  sys.exit(1)

(pmin, pmax) = map(float, sys.argv[argc-2:argc])
filelist = sys.argv[1:argc-2]

buf = [None]*len(filelist)

window=0

for ilc, lcfile in enumerate(filelist):
  # Read light curve
  lc = numpy.genfromtxt(lcfile,
                        dtype={"names": ("bjd", "mag", "e_mag", "texp",
                                       "dmag", "fwhm", "ellipt", "airmass",
                                         "xlc", "ylc", "angle",
                                         "skylev", "pkht",
                                         "s", "v", "r", "f",
                                         "cm", "corr_mag"),
                               "formats": ("f8", "f4", "f4", "f4",
                                           "f4", "f4", "f4", "f4",
                                           "f8", "f8", "f4",
                                           "f4", "f4",
                                           "i1", "i1", "i1", "i1",
                                           "f4", "f4") })

  # Make sure sorted on time.
  numpy.sort(lc, order="bjd")

  # 5-sigma clip.
  y = lc["mag"]
  ymed = numpy.median(y)
  ysig = 1.4826*numpy.median(numpy.absolute(y-ymed))

  lcclip = numpy.extract(numpy.absolute(y-ymed) < 5*ysig, lc)

  # After that, which segments are populated?
  segs = {s: -1 for s in lcclip["s"]}

  # Assign new indices.
  for i, s in enumerate(sorted(segs.keys())):
    segs[s] = i

  # Remap them to decide where we put the new DCs.
  idc = [segs[s] for s in lcclip["s"]]

  print "For", lcfile, "fitting", len(segs), "DC offsets"

  # Take off first timestamp.
  t = lcclip["bjd"] - lc[0]["bjd"]

  # Keep track of largest window for deciding frequency sampling.
  if t[-1] > window:
      window = t[-1]

  # Make matrix of external parameters (yuk).
  ep = numpy.empty([2, len(t)], numpy.double)
  ep[0] = lcclip["cm"]
  ep[1] = lcclip["fwhm"]

  # Add a tuple to the buffer.
  buf[ilc] = (t, lcclip["mag"], 1.0 / (lcclip["e_mag"]**2),
              ep,     # external parameters
              idc,    # DCs
              None)   # Sinusoids

# Sampling in frequency (with oversampling).
vsamp = 0.1 / window

# Period range.
pl = int(math.floor(1.0 / (vsamp*pmax)))
if pl < 1:
  pl = 1

ph = int(math.ceil(1.0 / (vsamp*pmin)))

nn = ph-pl+1

# Perform period search.
(pergrm, winfunc) = sfit.search(buf, pl, ph, vsamp)

# Best period.
p = numpy.argmin(pergrm)

# Parabolic interpolation for a better estimate.
if p > 0 and p < nn-1:
    aa = pergrm[p]
    bb = 0.5*(pergrm[p+1] - pergrm[p-1])
    cc = 0.5*(pergrm[p+1] + pergrm[p-1] - 2.0*aa)
    offset = -0.5*bb/cc
else:
    offset = 0.0

vbest = (pl+p+offset)*vsamp

print "Best period", 1.0/vbest, "days"

(chinull, bnull) = sfit.null(buf)
print "Null hypothesis", chinull, bnull

(chialt, balt) = sfit.single(buf, vbest)
print "Alternate hypothesis", chialt, balt

# Frequency grid for plot.
v = numpy.linspace(pl, ph, nn)
v *= vsamp

# "Amplitude spectrum" equivalent is square root of chi^2.
ampspec = numpy.sqrt(pergrm)

# Plots.
npanel = len(buf)+2

for ilc, lc in enumerate(buf):
    (t, y, wt, ep, idc, iamp) = lc

    b = balt[ilc]

    if ep is not None:
      nep = ep.shape[0]
    else:
      nep = 1

    if idc is not None:
      ndc = numpy.max(idc)+1
    else:
      ndc = 1

    if iamp is not None:
      namp = numpy.max(iamp)+1
    else:
      namp = 1

    # This is how the "b" array is packed.
    bdc = b[0:ndc]        # DCs
    bep = b[ndc:ndc+nep]  # external parameters
    bsc = b[ndc+nep:]     # sin, cos, sin, cos, ...

    plt.subplot(npanel, 1, ilc+1)
    plt.axis([0.0, 1.0, numpy.max(y), numpy.min(y)])

    phase = numpy.fmod(vbest*t, 1.0)
    blm = bdc[idc] + numpy.dot(bep, ep)

    plt.plot(phase, y-blm, ".")
    
    modx = numpy.linspace(0.0, 1.0, 1000)
    modp = 2*math.pi*modx

    mody = bsc[0] * numpy.sin(modp) + bsc[1] * numpy.cos(modp)

    print "Amplitude", math.sqrt(bsc[0]**2 + bsc[1]**2)

    plt.plot(modx, mody)

plt.subplot(npanel, 1, npanel-1)
plt.axis([numpy.min(v), numpy.max(v), numpy.max(ampspec), numpy.min(ampspec)])
plt.plot(v, ampspec)
plt.plot([vbest, vbest], plt.ylim(), linestyle='--')
plt.subplot(npanel, 1, npanel)
plt.axis([numpy.min(v), numpy.max(v), 0.0, 1.0])
plt.plot(v, winfunc)
plt.plot([vbest, vbest], plt.ylim(), linestyle='--')

plt.show()

sys.exit(0)