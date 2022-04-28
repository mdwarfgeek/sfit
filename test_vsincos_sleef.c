#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void FUNC (double v, double *t, int ndp,
           double *sinarg, double *cosarg);

int main (int argc, char *argv[]) {
  unsigned int ia, ib, ic, id;

  double t[4] = { 0, 0.25, 0.5, 0.75 };
  double sans[4] = { 0.0, 1.0,  0.0, -1.0 };
  double cans[4] = { 1.0, 0.0, -1.0, 0.0 };

  double s[4] = { -99.9, -99.9, -99.9, -99.9 };
  double c[4] = { -99.9, -99.9, -99.9, -99.9 };

  int i, rv;

  __cpuid(1, ia, ib, ic, id);
  if(ic & BIT) {
    printf("CPU supports " FEAT "\n");
  }
  else {
    printf("does not support " FEAT "\n");
    exit(1);
  }

  FUNC(1.0, t, 4, s, c);

  rv = 0;

  for(i = 0; i < 4; i++) {
    printf("sincos(2*PI*%.20le) = (%.20le, %.20le)\n", t[i], s[i], c[i]);

    if(fabs(s[i]-sans[i]) > 1.0e-15)
      rv = 1;
    if(fabs(c[i]-cans[i]) > 1.0e-15)
      rv = 1;
  }

  if(rv)
    printf("Answers are incorrect, failing test\n");
  else
    printf("Answers seem to be correct\n");

  return rv;
}
