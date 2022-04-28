#include <cpuid.h>

#undef FUNC
#define FUNC vsincos_sleef_sse4

#undef BIT
#define BIT bit_SSE4_1

#undef FEAT
#define FEAT "SSE"

#include "test_vsincos_sleef.c"
