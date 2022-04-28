#include <math.h>
#include <immintrin.h>

#include <sleef.h>

#undef FUNC
#define FUNC vsincos_sleef_avx

#undef REGD
#define REGD __m256d

#undef REGD2
#define REGD2 Sleef___m256d_2

#undef VLEN
#define VLEN 4

#undef MM_PREFIX
#define MM_PREFIX(f) _mm256_ ## f

#undef SLEEF_FUNC
#define SLEEF_FUNC(f, a) f ## d4_ ## a

#include "vsincos_sleef.c"
