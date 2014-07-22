#include <stdlib.h>
#include <stdint.h>

#ifdef _WIN32

#include <windows.h>

int get_num_cpus (void) {
  SYSTEM_INFO si;

  GetSystemInfo(&si);

  return(si.dwNumberOfProcessors);
}

#else  /* assume POSIX */

#include <sys/types.h>
#include <unistd.h>

#if defined(_SC_NPROCESSORS_ONLN)

int get_num_cpus (void) {
  return(sysconf(_SC_NPROCESSORS_ONLN));
}

#elif defined(_SC_NPROC_ONLN)

int get_num_cpus (void) {
  return(sysconf(_SC_NPROC_ONLN));
}

#elif defined(__FreeBSD__) || defined(__NetBSD__) || \
      defined(__OpenBSD__) || defined(__APPLE__)

#include <sys/param.h>
#include <sys/sysctl.h>

int get_num_cpus (void) {
  int mib[4] = {
    CTL_HW,
#if defined(HW_AVAILCPU)
    HW_AVAILCPU
#else
    HW_NCPU
#endif
  };
  int rv;
  size_t len = sizeof(rv);

  if(sysctl(mib, 2, &rv, &len, NULL, 0))
    return(-1);

  return(rv);
}

#else  /* don't know */

int get_num_cpus (void) {
  return(-1);
}

#endif

#endif
