#!/usr/bin/env python

from __future__ import print_function

import sys
import math
import numpy
import sfit

import matplotlib.pyplot as plt
import matplotlib.gridspec

# Figure size based on pgplot.
figsize = (10.5, 7.8)  # inches

argc = len(sys.argv)

if argc < 3:
  print("Usage:\t", sys.argv[0], "lcfile [...] pmin pmax")
  sys.exit(1)

(pmin, pmax) = list(map(float, sys.argv[argc-2:argc]))
filelist = sys.argv[1:argc-2]

buf = [None]*len(filelist)

window=0

bjdbase = None

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
  lc = numpy.sort(lc, order="bjd")

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

  print("For", lcfile, "fitting", len(segs), "DC offsets")

  # Take off floor of first timestamp.
  if bjdbase is None:
    bjdbase = math.floor(lc[0]["bjd"])

  t = lcclip["bjd"] - bjdbase

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

# "Amplitude spectrum" equivalent is square root of chi^2.
ampspec = numpy.sqrt(pergrm)

# Best period.
p = numpy.argmin(ampspec)

# Parabolic interpolation for a better estimate.
if p > 0 and p < nn-1:
    aa = ampspec[p]
    bb = 0.5*(ampspec[p+1] - ampspec[p-1])
    cc = 0.5*(ampspec[p+1] + ampspec[p-1] - 2.0*aa)
    offset = -0.5*bb/cc
else:
    offset = 0.0

vbest = (pl+p+offset)*vsamp
pbest = 1.0 / vbest

print("Best period", pbest, "days")

(chinull, bnull, bcovnull) = sfit.null(buf)
#print "Null hypothesis", chinull, bnull, bcovnull

(chialt, balt, bcovalt) = sfit.single(buf, vbest)
#print "Alternate hypothesis", chialt, balt, bcovalt

# Frequency grid for plot.
v = numpy.linspace(pl, ph, nn)
v *= vsamp

# As period.
pspec = 1.0 / v

# Plots.
nbin = 100

for ilc, lc in enumerate(buf):
    (t, y, wt, ep, idc, iamp) = lc

    b = balt[ilc]
    bcov = bcovalt[ilc]

    if ep is not None:
      nep = ep.shape[0]
    else:
      nep = 0

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
    bcovsc = bcov[ndc+nep:,ndc+nep:]
    
    # "Corrected" y array
    blm = bdc[idc]
    if nep > 0:
      blm += numpy.dot(bep, ep)
    ycorr = y - blm

    # Uncertainty.
    yerr = numpy.sqrt(1.0 / wt)
    
    # T0
    phi = math.atan2(bsc[1], bsc[0]) / (2.0*math.pi)
    if phi < 0.0:
      phi += 1.0

    ioff = round(numpy.median(t) * vbest)  # move to middle of data set

    t0 = bjdbase+(ioff-phi)/vbest

    print("Dataset", ilc+1, "T0 =", t0)

    fig = plt.figure(figsize=figsize)
    
    gs = matplotlib.gridspec.GridSpec(2, 1)
    
    gsh = matplotlib.gridspec.GridSpecFromSubplotSpec(2, 1, hspace=0, subplot_spec=gs[0])
    
    axt = fig.add_subplot(gsh[0])

    axt.errorbar(t, ycorr, yerr, fmt="none", ecolor="#A0A0A0", capsize=2, alpha=0.5)

    bins = numpy.linspace(t[0], t[-1], nbin+1)

    xb = 0.5*(bins[0:nbin] + bins[1:nbin+1])
    
    ybn = numpy.histogram(t, bins=bins, weights=ycorr*wt)[0]
    ybd = numpy.histogram(t, bins=bins, weights=wt)[0]

    wb = ybd > 0
    
    axt.plot(xb[wb], ybn[wb]/ybd[wb], "o", color="black")
    
    modx = numpy.linspace(t[0], t[-1], 1000)
    modp = 2*math.pi*vbest*modx

    mody = bsc[0] * numpy.sin(modp) + bsc[1] * numpy.cos(modp)

    axt.plot(modx, mody, color="red")

    axt.set_xlim(math.floor(t[0]), math.ceil(t[-1]))
    axt.set_ylim(numpy.max(ycorr), numpy.min(ycorr))

    axt.xaxis.tick_top()
    axt.xaxis.set_label_position("top")
    
    axt.set_xlabel("BJD $-$ {0:.1f}".format(bjdbase))
    axt.set_ylabel("Delta mag")

    axp = fig.add_subplot(gsh[1], sharey=axt)

    phase = numpy.fmod(vbest*t, 1.0)

    axp.errorbar(phase, ycorr, yerr, fmt="none", ecolor="#A0A0A0", capsize=2, alpha=0.5)

    bins = numpy.linspace(0.0, 1.0, nbin+1)

    xb = 0.5*(bins[0:nbin] + bins[1:nbin+1])
    
    ybn = numpy.histogram(phase, bins=bins, weights=ycorr*wt)[0]
    ybd = numpy.histogram(phase, bins=bins, weights=wt)[0]

    wb = ybd > 0
    
    axp.plot(xb[wb], ybn[wb]/ybd[wb], "o", color="black")
    
    modx = numpy.linspace(0.0, 1.0, 1000)
    modp = 2*math.pi*modx

    mody = bsc[0] * numpy.sin(modp) + bsc[1] * numpy.cos(modp)

    # Evaluate model for this light curve only, to get error scaling.
    thischisq = sfit.single([lc], vbest)[0]
    thisndof = len(y) - len(b) - 1  # 1 extra for period
    if thisndof > 0:
      thiserrscl = numpy.sqrt(thischisq / thisndof)
    else:
      thiserrscl = 1.0
      
    # Error propagation from sin,cos terms to amplitude.
    aa = bsc[0]*bsc[0]
    bb = bsc[1]*bsc[1]
    ab = bsc[0]*bsc[1]

    ampsq = aa + bb
    vvv = bcovsc[0,0]*aa + bcovsc[1,1]*bb + (bcovsc[1,0] + bcovsc[0,1])*ab

    amp = numpy.sqrt(ampsq)
    e_amp = numpy.sqrt(vvv / ampsq) * thiserrscl
    
    print("Amplitude {0:.6f} +/- {1:.6f}".format(amp, e_amp))
    
    axp.plot(modx, mody, color="red")

    axp.set_xlim(0.0, 1.0)
    axp.set_ylim(numpy.max(ycorr), numpy.min(ycorr))
    
    axp.set_ylabel("Delta mag")
    axp.set_xlabel("Phase")

    gsl = matplotlib.gridspec.GridSpecFromSubplotSpec(2, 1, hspace=0, subplot_spec=gs[1])
    
    axp = fig.add_subplot(gsl[0])
    axp.plot(pspec, ampspec, color="black")
    axp.axvline(pbest, color="red", alpha=0.5)
    if pbest > 2:  # alias periods
      axp.axvline(1.0/(1.0-vbest), color="green", alpha=0.25)
      axp.axvline(1.0/(1.0+vbest), color="green", alpha=0.25)
    axp.set_xlim(numpy.min(pspec), numpy.max(pspec))
    axp.set_ylim(numpy.max(ampspec), numpy.min(ampspec))
    axp.semilogx()
    axp.get_xaxis().set_visible(False)
    
    axp.set_ylabel("sqrt($\chi^2$)")
    
    axw = fig.add_subplot(gsl[1], sharex=axp)
    axw.plot(pspec, winfunc, color="black")
    axw.axvline(pbest, color="red", alpha=0.5)
    if pbest > 2:  # alias periods
      axw.axvline(1.0/(1.0-vbest), color="green", alpha=0.25)
      axw.axvline(1.0/(1.0+vbest), color="green", alpha=0.25)
    axw.set_xlim(numpy.min(pspec), numpy.max(pspec))
    axw.set_ylim(0.0, 0.99)
    axw.get_xaxis().set_visible(True)
    
    axw.set_xlabel("Period (days), best P = {0:.3f} days".format(pbest))
    axw.set_ylabel("Window function")
    
    plt.tight_layout()
    plt.show()
