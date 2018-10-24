#include <Python.h>
#include <structmember.h>
#include <numpy/arrayobject.h>

#include <stdlib.h>
#include <signal.h>

#include "sfit.h"

struct lc_list {
  struct sfit_lc *lclist;
  int nlc;
  
  PyObject **tarr;
  PyObject **yarr;
  PyObject **wtarr;
  PyObject **eparr;
  PyObject **idcarr;
  PyObject **iamparr;
};

static void cleanup_lc (struct lc_list *l) {
  int ilc;

  if(l->lclist) {
    free((void *) l->lclist);
    l->lclist = NULL;
  }

  for(ilc = 0; ilc < l->nlc; ilc++) {
    if(l->tarr) {
      Py_XDECREF(l->tarr[ilc]);
    }
    if(l->yarr) {
      Py_XDECREF(l->yarr[ilc]);
    }
    if(l->wtarr) {
      Py_XDECREF(l->wtarr[ilc]);
    }
    if(l->eparr) {
      Py_XDECREF(l->eparr[ilc]);
    }
    if(l->idcarr) {
      Py_XDECREF(l->idcarr[ilc]);
    }
    if(l->iamparr) {
      Py_XDECREF(l->iamparr[ilc]);
    }
  }

  if(l->tarr) {
    free((void *) l->tarr);
    l->tarr = NULL;
  }
  if(l->yarr) {
    free((void *) l->yarr);
    l->yarr = NULL;
  }
  if(l->wtarr) {
    free((void *) l->wtarr);
    l->wtarr = NULL;
  }
  if(l->eparr) {
    free((void *) l->eparr);
    l->eparr = NULL;
  }
  if(l->idcarr) {
    free((void *) l->idcarr);
    l->idcarr = NULL;
  }
  if(l->iamparr) {
    free((void *) l->iamparr);
    l->iamparr = NULL;
  }
}

static int get_lc (struct lc_list *l,
                   PyObject *list) {
  int ilc, nlc;
  struct sfit_lc *lc;

  PyObject *lcargs;

  PyObject *targ, *yarg, *wtarg, *eparg, *idcarg, *iamparg;
  int maxdc, maxamp;

  int ndim, tdp;
  npy_intp *dim;

  int idp;

  /* Deal with light curves */
  nlc = PySequence_Length(list);

  l->lclist = (struct sfit_lc *) malloc(nlc * sizeof(struct sfit_lc));
  l->tarr = (PyObject **) malloc(nlc * sizeof(PyObject *));
  l->yarr = (PyObject **) malloc(nlc * sizeof(PyObject *));
  l->wtarr = (PyObject **) malloc(nlc * sizeof(PyObject *));
  l->eparr = (PyObject **) malloc(nlc * sizeof(PyObject *));
  l->idcarr = (PyObject **) malloc(nlc * sizeof(PyObject *));
  l->iamparr = (PyObject **) malloc(nlc * sizeof(PyObject *));
  if(!l->lclist ||
     !l->tarr || !l->yarr || !l->wtarr ||
     !l->eparr || !l->idcarr || !l->iamparr) {
    PyErr_SetString(PyExc_MemoryError,
                    "malloc");
    goto error;
  }

  l->nlc = nlc;

  for(ilc = 0; ilc < nlc; ilc++) {
    l->tarr[ilc] = NULL;
    l->yarr[ilc] = NULL;
    l->wtarr[ilc] = NULL;
    l->eparr[ilc] = NULL;
    l->idcarr[ilc] = NULL;
    l->iamparr[ilc] = NULL;
  }
  
  for(ilc = 0; ilc < nlc; ilc++) {
    lc = l->lclist+ilc;

    lcargs = PySequence_GetItem(list, ilc);
    if(!lcargs)
      goto error;

    /* Extract tuple */
    targ    = NULL;
    yarg    = NULL;
    wtarg   = NULL;
    eparg   = NULL;
    idcarg  = NULL;
    iamparg = NULL;

    PyArg_ParseTuple(lcargs, "OOO|OOO",
                     &targ, &yarg, &wtarg,
                     &eparg,
                     &idcarg,
                     &iamparg);
    
    /* Done with lcargs */
    Py_DECREF(lcargs);

#define GETNP(in, type, out)                                              \
    (out) = PyArray_FROM_OTF((in), (type), NPY_IN_ARRAY | NPY_FORCECAST); \
    if(!(out))                                                            \
      goto error;

    /* t, y, wt */
    GETNP(targ, NPY_DOUBLE, l->tarr[ilc])
    GETNP(yarg, NPY_DOUBLE, l->yarr[ilc])
    GETNP(wtarg, NPY_DOUBLE, l->wtarr[ilc])

    lc->ndp = PyArray_Size(l->tarr[ilc]);

    if(PyArray_Size(l->yarr[ilc]) != lc->ndp) {
      PyErr_SetString(PyExc_IndexError,
                      "array 'y' is not same length as 't'");
      goto error;
    }

    if(PyArray_Size(l->wtarr[ilc]) != lc->ndp) {
      PyErr_SetString(PyExc_IndexError,
                      "array 'wt' is not same length as 't'");
      goto error;
    }

    lc->t = PyArray_DATA(l->tarr[ilc]);
    lc->y = PyArray_DATA(l->yarr[ilc]);
    lc->wt = PyArray_DATA(l->wtarr[ilc]);

    if(eparg && eparg != Py_None) {
      /* External parameters */
      GETNP(eparg, NPY_DOUBLE, l->eparr[ilc]);

      ndim = PyArray_NDIM(l->eparr[ilc]);

      if(ndim > 0) {
        dim = PyArray_DIMS(l->eparr[ilc]);

        if(ndim > 1) {
          lc->nep = dim[0];
          tdp = dim[1];
        }
        else {
          lc->nep = 1;
          tdp = dim[0];
        }
        
        if(tdp != lc->ndp) {
          PyErr_SetString(PyExc_IndexError,
                          "array 'ep' is not same length as 't'");
          goto error;
        }
        
        lc->ep = PyArray_DATA(l->eparr[ilc]);
      }
      else {
        lc->ep = NULL;
        lc->nep = 0;
      }
    }
    else {
      lc->ep = NULL;
      lc->nep = 0;
    }

    if(idcarg && idcarg != Py_None) {
      /* DC */
      GETNP(idcarg, NPY_INT, l->idcarr[ilc]);

      if(PyArray_Size(l->idcarr[ilc]) != lc->ndp) {
        PyErr_SetString(PyExc_IndexError,
                        "array 'idc' is not same length as 't'");
        goto error;
      }

      lc->idc = PyArray_DATA(l->idcarr[ilc]);


      /* Figure out how many DCs */
      maxdc = 0;

      for(idp = 0; idp < lc->ndp; idp++)
        if(lc->idc[idp] > maxdc)
          maxdc = lc->idc[idp];

      lc->ndc = maxdc+1;
    }
    else {
      lc->idc = NULL;
      lc->ndc = 0;
    }

    if(iamparg && iamparg != Py_None) {
      /* Amplitudes */
      GETNP(iamparg, NPY_INT, l->iamparr[ilc]);

      if(PyArray_Size(l->iamparr[ilc]) != lc->ndp) {
        PyErr_SetString(PyExc_IndexError,
                        "array 'iamp' is not same length as 't'");
        goto error;
      }

      lc->iamp = PyArray_DATA(l->iamparr[ilc]);

      /* Figure out how many sin,cos */
      maxamp = 0;

      for(idp = 0; idp < lc->ndp; idp++)
        if(lc->iamp[idp] > maxamp)
          maxamp = lc->iamp[idp];

      lc->namp = maxamp+1;
    }
    else {
      lc->iamp = NULL;
      lc->namp = 1;
    }
  }

  return(0);

 error:
  cleanup_lc(l);

  return(1);
}

static PyObject *wrap_search (PyObject *self,
                              PyObject *args,
                              PyObject *keywds) {
  static char *kwlist[] = { "list", "pl", "ph", "vsamp", "nthr", NULL };
  PyObject *list = NULL;

  int pl, ph;
  double vsamp;
  int nthr = -1;

  struct lc_list l;

  npy_intp outdim[1];
  PyObject *chisq = NULL, *winfunc = NULL;

  struct sigaction sa, old_sa;
  int rv;

  /* Make sure lc_list is empty, all pointers null */
  memset(&l, 0, sizeof(l));

  /* Get arguments */
  if(!PyArg_ParseTupleAndKeywords(args, keywds, "Oiid|i", kwlist,
                                  &list, &pl, &ph, &vsamp, &nthr))
    goto error;

  /* Get LC arguments */
  if(get_lc(&l, list))
    goto error;

  /* Create output arrays */
  outdim[0] = ph-pl+1;

  chisq = PyArray_SimpleNew(1, outdim, NPY_DOUBLE);
  if(!chisq)
    goto error;

  winfunc = PyArray_SimpleNew(1, outdim, NPY_DOUBLE);
  if(!winfunc)
    goto error;

  /* Python catches SIGINT and throws an exception.  Here we restore the
     standard system handling while the search is running so it's possible
     to quit.  This was done to avoid having to dig into the nightmare
     of signal handling with pthreads. */
  memset(&sa, 0, sizeof(sa));
  sa.sa_handler = SIG_DFL;
  sigemptyset(&(sa.sa_mask));

  if(sigaction(SIGINT, &sa, &old_sa)) {
    PyErr_SetString(PyExc_RuntimeError,
                    "sigaction");
    goto error;
  }

  /* Compute */
  rv = sfit_search(l.lclist, l.nlc, pl, ph, vsamp, nthr,
                   (double *) PyArray_DATA(chisq),
                   (double *) PyArray_DATA(winfunc));

  /* Restore Python's signal handler */
  if(sigaction(SIGINT, &old_sa, NULL)) {
    PyErr_SetString(PyExc_RuntimeError,
                    "sigaction");
    goto error;
  }

  if(rv) {
    PyErr_SetString(PyExc_RuntimeError,
                    "sfit_search failed");
    goto error;
  }

  /* Done */
  cleanup_lc(&l);

  return(Py_BuildValue("NN",
                       PyArray_Return((PyArrayObject *) chisq),
                       PyArray_Return((PyArrayObject *) winfunc)));

 error:
  cleanup_lc(&l);

  PyArray_XDECREF_ERR((PyArrayObject *) chisq);
  PyArray_XDECREF_ERR((PyArrayObject *) winfunc);

  return(NULL);
}

static PyObject *wrap_null (PyObject *self,
                            PyObject *args,
                            PyObject *keywds) {
  static char *kwlist[] = { "list", NULL };
  PyObject *list = NULL;

  struct lc_list l;

  double *b = (double *) NULL;
  double *bcov = (double *) NULL;
  int bstride;

  double chisq;

  int rv;

  PyObject *bout = (PyObject *) NULL;
  PyObject *bcovout = (PyObject *) NULL;
  npy_intp outdim[3];

  int ilc, i, j, ncoeff;
  double *pin, *pout;

  /* Make sure lc_list is empty, all pointers null */
  memset(&l, 0, sizeof(l));

  /* Get arguments */
  if(!PyArg_ParseTupleAndKeywords(args, keywds, "O", kwlist,
                                  &list))
    goto error;

  /* Get LC arguments */
  if(get_lc(&l, list))
    goto error;

  /* Compute */
  rv = sfit_null(l.lclist, l.nlc,
                 &b, &bcov, &bstride,
                 &chisq);

  if(rv) {
    PyErr_SetString(PyExc_RuntimeError,
                    "sfit_null failed");
    goto error;
  }

  /* Create arrays */
  outdim[0] = l.nlc;
  outdim[1] = bstride;
  outdim[2] = bstride;

  bout = PyArray_SimpleNew(2, outdim, NPY_DOUBLE);
  if(!bout)
    goto error;

  bcovout = PyArray_SimpleNew(3, outdim, NPY_DOUBLE);
  if(!bcovout)
    goto error;

  memcpy(PyArray_DATA(bout), b,
         outdim[0]*outdim[1]*sizeof(double));

  pout = PyArray_DATA(bcovout);

  for(ilc = 0; ilc < l.nlc; ilc++) {
    ncoeff = l.lclist[ilc].ncoeff;

    pin = bcov + ilc*bstride*bstride;

    for(i = 0; i < ncoeff; i++) {
      for(j = 0; j < ncoeff; j++, pin++, pout++)
        *pout = *pin;

      for(; j < bstride; j++, pout++)
        *pout = 0;
    }

    for(; i < bstride; i++)
      for(j = 0; j < bstride; j++, pout++)
        *pout = 0;
  }

  /* Done */
  cleanup_lc(&l);

  free((void *) b);
  b = (double *) NULL;
  free((void *) bcov);
  bcov = (double *) NULL;

  return(Py_BuildValue("dNN",
                       chisq,
                       PyArray_Return((PyArrayObject *) bout),
                       PyArray_Return((PyArrayObject *) bcovout)));

 error:
  cleanup_lc(&l);

  if(b)
    free((void *) b);
  if(bcov)
    free((void *) bcov);

  PyArray_XDECREF_ERR((PyArrayObject *) bout);
  PyArray_XDECREF_ERR((PyArrayObject *) bcovout);

  return(NULL);
}

static PyObject *wrap_single (PyObject *self,
                              PyObject *args,
                              PyObject *keywds) {
  static char *kwlist[] = { "list", "v", NULL };
  PyObject *list = NULL;

  double v;

  struct lc_list l;

  double *b = (double *) NULL;
  double *bcov = (double *) NULL;
  int bstride;

  double chisq;

  int rv;

  PyObject *bout = (PyObject *) NULL;
  PyObject *bcovout = (PyObject *) NULL;
  npy_intp outdim[3];

  int ilc, i, j, ncoeff;
  double *pin, *pout;

  /* Make sure lc_list is empty, all pointers null */
  memset(&l, 0, sizeof(l));

  /* Get arguments */
  if(!PyArg_ParseTupleAndKeywords(args, keywds, "Od", kwlist,
                                  &list, &v))
    goto error;

  /* Get LC arguments */
  if(get_lc(&l, list))
    goto error;

  /* Compute */
  rv = sfit_single(l.lclist, l.nlc, v,
                   &b, &bcov, &bstride,
                   &chisq);

  if(rv) {
    PyErr_SetString(PyExc_RuntimeError,
                    "sfit_single failed");
    goto error;
  }

  /* Create array */
  outdim[0] = l.nlc;
  outdim[1] = bstride;
  outdim[2] = bstride;

  bout = PyArray_SimpleNew(2, outdim, NPY_DOUBLE);
  if(!bout)
    goto error;

  bcovout = PyArray_SimpleNew(3, outdim, NPY_DOUBLE);
  if(!bcovout)
    goto error;

  memcpy(PyArray_DATA(bout), b,
         outdim[0]*outdim[1]*sizeof(double));

  pout = PyArray_DATA(bcovout);

  for(ilc = 0; ilc < l.nlc; ilc++) {
    ncoeff = l.lclist[ilc].ncoeff + l.lclist[ilc].namp*2;

    pin = bcov + ilc*bstride*bstride;

    for(i = 0; i < ncoeff; i++) {
      for(j = 0; j < ncoeff; j++, pin++, pout++)
        *pout = *pin;

      for(; j < bstride; j++, pout++)
        *pout = 0;
    }

    for(; i < bstride; i++)
      for(j = 0; j < bstride; j++, pout++)
        *pout = 0;
  }
  
  /* Done */
  cleanup_lc(&l);

  free((void *) b);
  b = (double *) NULL;
  free((void *) bcov);
  bcov = (double *) NULL;

  return(Py_BuildValue("dNN",
                       chisq,
                       PyArray_Return((PyArrayObject *) bout),
                       PyArray_Return((PyArrayObject *) bcovout)));

 error:
  cleanup_lc(&l);

  if(b)
    free((void *) b);
  if(bcov)
    free((void *) bcov);

  PyArray_XDECREF_ERR((PyArrayObject *) bout);
  PyArray_XDECREF_ERR((PyArrayObject *) bcovout);

  return(NULL);
}

static PyMethodDef sfit_methods[] = {
  { "search", (PyCFunction) wrap_search, METH_VARARGS | METH_KEYWORDS,
    "(chisq, winfunc) = search(list, pl, ph, vsamp, nthr=-1)\n\n"
    "Main subroutine to compute the periodogram.\n\n"
    "Required arguments:\n"
    "list  -- Python list of light curves to fit.\n"
    "pl    -- Range of frequency sample index to compute.\n"
    "ph\n"
    "vsamp -- Sampling in frequency (1/time units).\n\n"
    "Optional arguments:\n"
    "nthr  -- Number of threads (-1 = number of CPUs).\n\n"
    "The list argument is a Python list of tuples containing the\n"
    "following arguments for each light curve:\n"
    " (t, y, wt, ep, idc, ndc, iamp, namp)\n"
    "t, y and wt are the time, data, and weight vectors.\n"
    "ep contains vectors of external parameters.\n"
    "idc and iamp are integer indices for fitting separate DC offsets\n"
    "or amplitudes.\n\n"
  },
  { "null", (PyCFunction) wrap_null, METH_VARARGS | METH_KEYWORDS,
    "(chisq, b, bcov) = null(list)\n\n"
    "Computes null hypothesis model.\n\n"
  },
  { "single", (PyCFunction) wrap_single, METH_VARARGS | METH_KEYWORDS,
    "(chisq, b, bcov) = single(list, v)\n\n"
    "Computes single fit at frequency v.\n\n"
  },
  { NULL, NULL, 0, NULL }
};

#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef sfit_mod = {
  PyModuleDef_HEAD_INIT,
  "sfit",       /* m_name */
  NULL,         /* m_doc */
  -1,           /* m_size */
  sfit_methods, /* m_methods */
  NULL,         /* m_reload */
  NULL,         /* m_traverse */
  NULL,         /* m_clear */
  NULL,         /* m_free */
};
#endif

#if PY_MAJOR_VERSION >= 3
PyMODINIT_FUNC PyInit_sfit (void) {
#else
PyMODINIT_FUNC initsfit (void) {
#endif
  PyObject *m;

  /* Init module */
#if PY_MAJOR_VERSION >= 3
  m = PyModule_Create(&sfit_mod);
#else
  m = Py_InitModule("sfit", sfit_methods);
#endif
  if(!m)
    goto error;

  /* Import numpy */
  import_array();

#if PY_MAJOR_VERSION >= 3
  return(m);
  
 error:
  return(NULL);
#else
 error:
  return;
#endif
}
