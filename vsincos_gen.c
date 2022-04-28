#include <math.h>

#include "sincos.h"

void vsincos_gen (double v, double *t, int ndp,
                  double *sinarg, double *cosarg) {
  int idp;
  double phi;

  for(idp = 0; idp < ndp; idp++) {
    /* Range reduction */
    phi = v * t[idp];
    phi -= rint(phi);
    
    /* Compute sin, cos, preferably simultaneously.
       Results stored to arrays. */
    inline_bare_sincos(2*M_PI*phi, sinarg[idp], cosarg[idp]);
  }
}
