#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <assert.h>
#include <pthread.h>

#include "sfit.h"

#ifndef USE_LAPACK
#include <lfa.h>
#endif

#if defined(HAVE_SLEEF_SSE) || defined(HAVE_SLEEF_AVX)
#define HAVE_SLEEF
#include <cpuid.h>
#endif

/* Structure to pass between threads */

struct sfit_info {
  struct sfit_lc *lclist;
  int nlc;
  int pl;
  int ph;
  double vsamp;

  int ncoeffmax;
  int ncovmax;
  int nmeasmax;

  pthread_mutex_t mm;

  volatile int p;
  pthread_mutex_t mp;

  volatile double *chisqper;
  volatile double *winfunc;
  pthread_mutex_t mo;

  void (*vsincos_func) (double, double *, int, double *, double *);
};

/* Prototypes for external stuff */

int get_num_cpus (void);

#ifdef USE_LAPACK
void dgesvd_ (char *, char *,
	      int *, int *, double *, int *,
	      double *, double *, int *, double *, int *,
	      double *, int *, int *, unsigned int, unsigned int);
void dgelss_ (int *, int *, int *,
              double *, int *, double *, int *, double *, double *, int *,
              double *, int *, int *);

int svdsolcov (double *a, double *b, double *cov, int m, double scale);
#endif

#ifdef HAVE_SLEEF_SSE4
void vsincos_sleef_sse4 (double v, double *t, int ndp,
                         double *sinarg, double *cosarg);
#endif
#ifdef HAVE_SLEEF_AVX
void vsincos_sleef_avx (double v, double *t, int ndp,
                        double *sinarg, double *cosarg);
#endif
void vsincos_gen (double v, double *t, int ndp,
                  double *sinarg, double *cosarg);

/* Main calculation */

static inline int sfit_compute (struct sfit_info *inf,
                                double *sinarg,
                                double *cosarg,
                                double v,
                                int ncoeffmax,
                                int ncovmax,
                                int dosc,  /* 0 = null hyp., 1 = alt hyp. */
                                double *b_r,
                                double *bcov_r,
                                volatile double *chisq_r,
                                volatile double *winfunc_r) {
  int ilc, idp, iep, k, j, ncoeff;
  int pamp;

  struct sfit_lc *lc;
  double y, wt;

  double chisq, a0, a1, a2;
  double a[ncovmax], b[ncoeffmax], c[ncoeffmax];

#ifdef USE_LAPACK
  int lwork = 6*ncoeffmax;
  double s[ncoeffmax], dwork[lwork], rcond;
  int nrhs, rank, info;
#else
  double s[ncoeffmax], beta[ncoeffmax];
  int perm[ncoeffmax];
#endif

  double mod, err;

  int nerr;

  nerr = 0;

  chisq = 0;

  a0 = 0;
  a1 = 0;
  a2 = 0;

  for(ilc = 0; ilc < inf->nlc; ilc++) {
    lc = inf->lclist+ilc;

    /* Required by the optimized matrix calculation below. */
    assert(lc->off_dc == 0);

    if(dosc) {
      inf->vsincos_func(v, lc->t, lc->ndp, sinarg, cosarg);

      ncoeff = lc->ncoeff + 2*lc->namp;
    }
    else
      ncoeff = lc->ncoeff;

    /* Init matrices */
    for(k = 0; k < ncoeff; k++) {
      for(j = 0; j < ncoeff; j++)
        a[k*ncoeff+j] = 0;

      b[k] = 0;
    }

    /* Accumulate sums for fit */
    for(idp = 0; idp < lc->ndp; idp++) {
      y = lc->y[idp];
      wt = lc->wt[idp];

      /* External parameters */
      for(iep = 0; iep < lc->nep; iep++)
        c[lc->off_ep + iep] = lc->ep[iep*lc->ndp + idp];

      if(dosc) {
        /* sin, cos */
        if(lc->namp > 1)
          pamp = lc->off_amp + 2*lc->iamp[idp];
        else
          pamp = lc->off_amp;
        
        c[pamp]   = sinarg[idp];
        c[pamp+1] = cosarg[idp];
        
        a0 += wt;
        a1 += sinarg[idp] * wt;
        a2 += cosarg[idp] * wt;
      }

      /* There's only one DC for each light curve point, meaning the
         upper left square block of the matrix A is diagonal.  Here
         we take advantage of this by splitting the matrix into four:
         
         A = [ D   E^T ]
             [ E   R   ]

         where D is a diagonal matrix and we need to add 1 * wt to
         the diagonal for the appropriate DC.  The calculation of E
         can also be optimized, only one column corresponding to
         the DC for this point needs to be updated. */

      if(lc->ndc > 0) {
        /* Column we need to update */
        j = lc->off_dc + lc->idc[idp];

        /* D */
        a[j*ncoeff+j] += wt;

        /* and the corresponding element of b */
        b[j] += y * wt;

        /* E */
        for(k = lc->off_ep; k < ncoeff; k++)
          a[k*ncoeff+j] += c[k] * wt;
      }

      /* R */
      for(k = lc->off_ep; k < ncoeff; k++) {
        for(j = lc->off_ep; j <= k; j++)
          a[k*ncoeff+j] += c[k]*c[j] * wt;

        b[k] += y * c[k] * wt;
      }
    }

    /* Fill in duplicates */
    for(k = 0; k < ncoeff; k++)
      for(j = k+1; j < ncoeff; j++)
        a[k*ncoeff+j] = a[j*ncoeff+k];

    if(ncoeff > 0) {
#ifdef USE_LAPACK
      /* Do they want covariance? */
      if(bcov_r)
        info = svdsolcov(a, b,
                         &(bcov_r[ilc*ncoeffmax*ncoeffmax]),
                         ncoeff, -1.0);
      else {
        /* Solve */
        nrhs = 1;
        rcond = -1.0;
        
        dgelss_(&ncoeff, &ncoeff, &nrhs,
                a, &ncoeff, b, &ncoeff, s, &rcond, &rank,
                dwork, &lwork, &info);
      }
      
      /* XXX - error handling */
      
      if(info)
        fprintf(stderr, "dgelss: %d", info);
#else
      qr(a, s, beta, perm, ncoeff);
      qrsolve(a, s, beta, perm, b, ncoeff, -1.0);

      if(bcov_r)
        qrinvert(a, s, beta, perm,
                 &(bcov_r[ilc*ncoeffmax*ncoeffmax]),
                 ncoeff, -1.0);
#endif
    }

    /* Do they want b? */
    if(b_r) {
      for(k = 0; k < ncoeff; k++)
        b_r[ilc*ncoeffmax+k] = b[k];
    }

    /* Accumulate chi^2 */
    for(idp = 0; idp < lc->ndp; idp++) {
      if(lc->ndc > 0)
        mod = b[lc->off_dc + lc->idc[idp]];
      else
        mod = 0;

      if(lc->namp > 1)
        pamp = lc->off_amp + 2*lc->iamp[idp];
      else
        pamp = lc->off_amp;

      if(dosc)
        mod += b[pamp] * sinarg[idp] + b[pamp+1] * cosarg[idp];

      /* External parameters */
      for(iep = 0; iep < lc->nep; iep++)
        mod += b[lc->off_ep + iep] * lc->ep[iep*lc->ndp+idp];

      /* Error */
      err = lc->y[idp] - mod;

      /* Chi^2 */
      chisq += err*err * lc->wt[idp];
    }
  }

  if(chisq_r)
    *chisq_r = chisq;

  if(winfunc_r)
    *winfunc_r = sqrt(a1*a1 + a2*a2) / a0;

  return(nerr);
}

static int sfit_search_work_single (struct sfit_info *inf) {
  double *sinarg = (double *) NULL;
  double *cosarg = (double *) NULL;

  int p;
  int rv = 0;

  /* Allocate storage */
  sinarg = (double *) malloc(inf->nmeasmax * sizeof(double));
  cosarg = (double *) malloc(inf->nmeasmax * sizeof(double));
  if(!sinarg || !cosarg)
    return(1);
  
  for(p = inf->pl; p <= inf->ph; p++) {
    /* Calculation for this period */
    rv += sfit_compute(inf, sinarg, cosarg,
                       p * inf->vsamp,
                       inf->ncoeffmax,
                       inf->ncovmax,
                       1,     /* alt. hyp. */
                       NULL,  /* don't need b */
                       NULL,
                       inf->chisqper+(p-inf->pl),
                       inf->winfunc+(p-inf->pl));
  }

  /* Free storage */
  free((void *) sinarg);
  free((void *) cosarg);

  return(rv);
}

static void *sfit_search_work_thread (void *data) {
  struct sfit_info *inf = (struct sfit_info *) data;

  double *sinarg = (double *) NULL;
  double *cosarg = (double *) NULL;

  int p, havework;
  double chisq, winfunc;

  size_t rv = 0;

  /* Allocate thread-local storage */
  pthread_mutex_lock(&(inf->mm));

  sinarg = (double *) malloc(inf->nmeasmax * sizeof(double));
  cosarg = (double *) malloc(inf->nmeasmax * sizeof(double));

  pthread_mutex_unlock(&(inf->mm));

  if(!sinarg || !cosarg) {
    rv = 1;
    pthread_exit((void *) rv);
  }
  
  for(;;) {
    /* Try to get more work */
    pthread_mutex_lock(&(inf->mp));

    if(inf->p <= inf->ph) {
      /* Claim this work */
      p = inf->p;
      inf->p++;
      havework = 1;
    }
    else  /* done */
      havework = 0;

    pthread_mutex_unlock(&(inf->mp));

    if(!havework)
      break;

    /* Run main computation for this period */
    rv += sfit_compute(inf, sinarg, cosarg,
                       p * inf->vsamp,
                       inf->ncoeffmax,
                       inf->ncovmax,
                       1,     /* alt. hyp. */
                       NULL,  /* don't need b */
                       NULL,
                       &chisq, &winfunc);

    /* Store results */
    pthread_mutex_lock(&(inf->mo));

    inf->chisqper[p-inf->pl] = chisq;
    inf->winfunc[p-inf->pl] = winfunc;

    pthread_mutex_unlock(&(inf->mo));
  }

  /* Free storage */
  pthread_mutex_lock(&(inf->mm));

  free((void *) sinarg);
  free((void *) cosarg);

  pthread_mutex_unlock(&(inf->mm));

  return((void *) rv);
}

void sfit_info (struct sfit_info *inf,
                struct sfit_lc *lclist, int nlc) {
  int ilc, ncoeffsc, ncoeffmax, nmeasmax;
#ifdef HAVE_SLEEF
  unsigned int ia, ib, ic, id;
#endif

  /* Compute number of coefficients in each LC */
  ncoeffmax = 0;
  nmeasmax = 0;

  for(ilc = 0; ilc < nlc; ilc++) {
    /* Base number of coefficients (excluding sin,cos) */
    lclist[ilc].ncoeff = lclist[ilc].ndc + lclist[ilc].nep;

    /* Number including sin and cos */
    ncoeffsc = lclist[ilc].ncoeff + 2*lclist[ilc].namp;

    lclist[ilc].off_dc = 0;
    lclist[ilc].off_ep = lclist[ilc].ndc;
    lclist[ilc].off_amp = lclist[ilc].off_ep + lclist[ilc].nep;

    if(ncoeffsc > ncoeffmax)
      ncoeffmax = ncoeffsc;
    if(lclist[ilc].ndp > nmeasmax)
      nmeasmax = lclist[ilc].ndp;
  }

  inf->lclist = lclist;
  inf->nlc = nlc;

  inf->ncoeffmax = ncoeffmax;
  inf->ncovmax = ncoeffmax*ncoeffmax;
  inf->nmeasmax = nmeasmax;

  inf->vsincos_func = vsincos_gen;

#ifdef HAVE_SLEEF
  /* Runtime check for CPU features */
  __cpuid(1, ia, ib, ic, id);
#ifdef HAVE_SLEEF_SSE4
  if(ic & bit_SSE4_1) {
    inf->vsincos_func = vsincos_sleef_sse4;
  }
#endif
#ifdef HAVE_SLEEF_AVX
  if(ic & bit_AVX) {
    inf->vsincos_func = vsincos_sleef_avx;
  }
#endif
#endif
}

int sfit_search (struct sfit_lc *lclist, int nlc,
                 int pl, int ph, double vsamp,
                 int nthr,
                 double *chisqper, double *winfunc) {
  struct sfit_info inf;

  pthread_t *thrl = (pthread_t *) NULL;
  int *thrf = (int *) NULL;

  int ithr;

  size_t trv;
  int rv = 0;

  /* Fill in basic parts of structure (common to all routines) */
  sfit_info(&inf, lclist, nlc);

  /* How many threads? */
  if(nthr < 0)
    nthr = get_num_cpus();

  if(nthr < 0)  /* failed */
    nthr = 1;

  /* Fill in rest of structure */
  inf.pl = pl;
  inf.ph = ph;
  inf.vsamp = vsamp;

  inf.p = pl;

  inf.chisqper = chisqper;
  inf.winfunc = winfunc;

  /* Threaded? */
  if(nthr > 1) {
    if(pthread_mutex_init(&(inf.mm), (pthread_mutexattr_t *) NULL))
      return(1);

    if(pthread_mutex_init(&(inf.mp), (pthread_mutexattr_t *) NULL))
      return(1);

    if(pthread_mutex_init(&(inf.mo), (pthread_mutexattr_t *) NULL))
      return(1);

    thrl = (pthread_t *) malloc(nthr * sizeof(pthread_t));
    thrf = (int *) malloc(nthr * sizeof(int));
    if(!thrl || !thrf)
      return(1);

    for(ithr = 0; ithr < nthr; ithr++) {
      if(pthread_create(thrl+ithr, (pthread_attr_t *) NULL,
                        sfit_search_work_thread, (void *) &inf)) {
        /* error */
        rv++;
        thrf[ithr] = 0;
      }
      else
        thrf[ithr] = 1;
    }

    for(ithr = 0; ithr < nthr; ithr++) {
      if(thrf[ithr]) {
        if(pthread_join(thrl[ithr], (void **) &trv))
          /* error */
          rv++;
        else
          rv += trv;
      }
    }

    pthread_mutex_destroy(&(inf.mm));
    pthread_mutex_destroy(&(inf.mp));
    pthread_mutex_destroy(&(inf.mo));

    free((void *) thrl);
    free((void *) thrf);
  }
  else {
    /* Run here instead */
    rv = sfit_search_work_single(&inf);
  }

  return(rv);
}

int sfit_null (struct sfit_lc *lclist, int nlc,
               double **b, double **bcov, int *bstride,
               double *chisq) {
  struct sfit_info inf;
  int rv;

  /* Fill in basic parts of structure (common to all routines) */
  sfit_info(&inf, lclist, nlc);

  /* Allocate storage */
  *b = (double *) malloc(nlc * inf.ncoeffmax * sizeof(double));
  *bcov = (double *) malloc(nlc * inf.ncovmax * sizeof(double));
  if(!(*b) || !(*bcov))
    return(1);

  memset(*b, 0, nlc * inf.ncoeffmax * sizeof(double));
  memset(*bcov, 0, nlc * inf.ncovmax * sizeof(double));

  rv = sfit_compute(&inf, NULL, NULL,
                    0,
                    inf.ncoeffmax,
                    inf.ncovmax,
                    0,     /* null hyp. */
                    *b,
                    *bcov,
                    chisq,
                    NULL);

  *bstride = inf.ncoeffmax;

  return(rv);
}

int sfit_single (struct sfit_lc *lclist, int nlc,
                 double v,
                 double **b, double **bcov, int *bstride,
                 double *chisq) {
  struct sfit_info inf;

  double *sinarg = (double *) NULL;
  double *cosarg = (double *) NULL;

  int rv;

  /* Fill in basic parts of structure (common to all routines) */
  sfit_info(&inf, lclist, nlc);

  /* Allocate storage */
  sinarg = (double *) malloc(inf.nmeasmax * sizeof(double));
  cosarg = (double *) malloc(inf.nmeasmax * sizeof(double));
  *b = (double *) malloc(nlc * inf.ncoeffmax * sizeof(double));
  *bcov = (double *) malloc(nlc * inf.ncovmax * sizeof(double));
  if(!sinarg || !cosarg || !(*b) || !(*bcov))
    return(1);

  memset(*b, 0, nlc * inf.ncoeffmax * sizeof(double));
  memset(*bcov, 0, nlc * inf.ncovmax * sizeof(double));

  rv = sfit_compute(&inf, sinarg, cosarg,
                    v,
                    inf.ncoeffmax,
                    inf.ncovmax,
                    1,     /* alt. hyp. */
                    *b,
                    *bcov,
                    chisq,
                    NULL);

  *bstride = inf.ncoeffmax;

  /* Free internal storage */
  free((void *) sinarg);
  free((void *) cosarg);

  return(rv);
}

#ifdef USE_LAPACK

/* Computes minimum norm solution and covariance matrix for a linear
   least squares problem.  Borrowed from library. */

int svdsolcov (double *a, double *b, double *cov, int m, double scale) {
  int lwork = 6*m;
  double s[m], u[m*m], vt[m*m], work[lwork];
  int mthr, info;

  int i, j, k;
  double thresh, sum, sol;

  /* Perform decomposition */
  dgesvd_("A", "A",
	  &m, &m,
	  a, &m, s, u, &m, vt, &m,
	  work, &lwork, &info, 1, 1);

  if(scale < 0)
    /* Machine precision */
    scale = FLT_RADIX * DBL_EPSILON;

  /* Threshold for singular values */
  thresh = scale*s[0];
  if(thresh < DBL_MIN)  /* smallest positive number */
    thresh = DBL_MIN;

  /* Take reciprocal of elements of s above threshold */
  for(k = 0; k < m && s[k] > thresh; k++)
    s[k] = 1.0 / s[k];

  mthr = m;

  /* Compute inverse and solution vector */
  for(i = 0; i < m; i++) {
    sol = 0;

    for(j = 0; j < m; j++) {
      /* i'th row of v multiplied by j'th column of u transpose */
      sum = 0;
      for(k = 0; k < mthr; k++)
	sum += vt[i*m+k]*u[k*m+j] * s[k];

      cov[i*m+j] = sum;

      sol += b[j]*sum;
    }

    work[i] = sol;
  }    

  for(i = 0; i < m; i++)
    b[i] = work[i];

  return(info);
}

#endif
