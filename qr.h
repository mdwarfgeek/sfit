#ifndef QR_H
#define QR_H

/* -- qr.c: QR decomposition and solution of linear equations -- */

void qr (double *a, double *s, double *betaarr, int *perm, int n);
int qrsolve (double *a, double *s, double *betaarr, int *perm,
             double *b, int n, double rcond);
int qrinvert (double *a, double *s, double *betaarr, int *perm,
              double *ainv, int n, double rcond);

#endif  /* QR_H */
