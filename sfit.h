#ifndef SFIT_H
#define SFIT_H

struct sfit_lc {
  /* Data arrays */
  double *t;
  double *y;
  double *wt;
  int ndp;

  /* DC offsets */
  int *idc;
  int ndc;

  /* External parameters */
  double *ep;
  int nep;

  /* Amplitudes */
  int *iamp;
  int namp;

  /* For internal use */
  int ncoeff;
  int off_dc;
  int off_ep;
  int off_amp;
};

int sfit_search (struct sfit_lc *lclist, int nlc,
                 int pl, int ph, double vsamp,
                 int nthr,
                 double *chisqper, double *winfunc);

int sfit_null (struct sfit_lc *lclist, int nlc,
               double **b, int *bstride,
               double *chisq);

int sfit_single (struct sfit_lc *lclist, int nlc,
                 double v,
                 double **b, int *bstride,
                 double *chisq);

#endif  /* SFIT_H */
