#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <pthread.h>

#include "sfit.h"
#include "sincos.h"

#ifndef USE_LAPACK
#include <lfa.h>
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
  double phi, wt;

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

    if(dosc) {
      /* Compute sin, cos */
      for(idp = 0; idp < lc->ndp; idp++) {
        phi = fmod(v * lc->t[idp], 1.0);
        inline_bare_sincos(2*M_PI*phi, sinarg[idp], cosarg[idp]);
      }

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
      wt = lc->wt[idp];

      /* Blank vector */
      for(k = 0; k < ncoeff; k++)
        c[k] = 0;

      /* DC */
      if(lc->ndc > 0)
        c[lc->off_dc + lc->idc[idp]] = 1.0;

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

      /* Outer product, interesting parts only */
      for(k = 0; k < ncoeff; k++) {
        for(j = 0; j <= k; j++)
          a[k*ncoeff+j] += c[k]*c[j] * wt;

        b[k] += lc->y[idp] * c[k] * wt;
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
