#ifndef SINCOS_H
#define SINCOS_H

#include <math.h>

#if defined(__GLIBC__)
#include <features.h>
#endif

#if defined(__GNUC__) && (defined(__i386) || defined(__amd64))

/* Inline x86 assembler routine, GNU C syntax */
#define inline_bare_sincos(a, s, c) \
  asm("fsincos" : "=t" (c), "=u" (s) : "0" (a))

#elif defined(__GLIBC__) && defined(_GNU_SOURCE)

/* Try to use glibc provided in math.h */
#define inline_bare_sincos(a, s, c) sincos((a), &(s), &(c))

#else

/* Generic implementation using standard math.h functions */
#define inline_bare_sincos(a, s, c) (s) = sin(a); (c) = cos(a)

#endif

#endif  /* SINCOS_H */
