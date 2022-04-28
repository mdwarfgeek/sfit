#include <math.h>
#include <emmintrin.h>

#include <sleef.h>

#undef FUNC
#define FUNC vsincos_sleef_sse4

#undef REGD
#define REGD __m128d

#undef REGD2
#define REGD2 Sleef___m128d_2

#undef VLEN
#define VLEN 2

#undef MM_PREFIX
#define MM_PREFIX(f) _mm_ ## f

#undef SLEEF_FUNC
#define SLEEF_FUNC(f, a) f ## d2_ ## a

#include "vsincos_sleef.c"
