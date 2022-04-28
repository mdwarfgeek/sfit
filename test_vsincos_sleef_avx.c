#include <cpuid.h>

#undef FUNC
#define FUNC vsincos_sleef_avx

#undef BIT
#define BIT bit_AVX

#undef FEAT
#define FEAT "AVX"

#include "test_vsincos_sleef.c"
