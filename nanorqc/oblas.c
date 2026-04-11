#include "oblas.h"

// To check CPUID
#ifdef _MSC_VER
#include <intrin.h>
#else // gcc
#include <cpuid.h>
#endif

void check_cpuid(void) {
  cpu_oswaprow = oswaprow;
  cpu_oaxpy = oaxpy;
  cpu_oscal = oscal;
  cpu_oaxpy_b32 = oaxpy_b32;

#ifdef _MSC_VER
  int cpu_info[4], max_id;

  __cpuid(cpu_info, 0);
  max_id = cpu_info[0]; // EAX
  if (max_id >= 1) {
    __cpuid(cpu_info, 1); // CPUID 1
#ifdef OBLAS_SSE
    if (cpu_info[2] & (1 << 9)) { // 1<<9 in ECX = SSSE3
      cpu_oswaprow = oswaprow_sse;
      cpu_oaxpy = oaxpy_sse;
      cpu_oscal = oscal_sse;
      cpu_oaxpy_b32 = oaxpy_b32_sse;
    }
#endif
  }
  if (max_id >= 7) {
    __cpuidex(cpu_info, 7, 0); // CPUID 7 0
#ifdef OBLAS_AVX
    if (cpu_info[1] & (1 << 5)) { // 1<<5 in EBX = AVX2
      cpu_oswaprow = oswaprow_avx;
      cpu_oaxpy = oaxpy_avx;
      cpu_oscal = oscal_avx;
      cpu_oaxpy_b32 = oaxpy_b32_avx;
    }
#endif
  }

#else // gcc
  unsigned int max_id, eax, ebx, ecx, edx;

  __cpuid(0, eax, ebx, ecx, edx);
  max_id = eax; // EAX
  if (max_id >= 1) {
    __cpuid(1, eax, ebx, ecx, edx); // CPUID 1
#ifdef OBLAS_SSE
    if (ecx & (1 << 9)) { // 1<<9 in ECX = SSSE3
      cpu_oswaprow = oswaprow_sse;
      cpu_oaxpy = oaxpy_sse;
      cpu_oscal = oscal_sse;
      cpu_oaxpy_b32 = oaxpy_b32_sse;
    }
#endif
  }
  if (max_id >= 7) {
    __cpuid_count(7, 0, eax, ebx, ecx, edx); // CPUID 7 0
#ifdef OBLAS_AVX
    if (ebx & (1 << 5)) { // 1<<5 in EBX = AVX2
      cpu_oswaprow = oswaprow_avx;
      cpu_oaxpy = oaxpy_avx;
      cpu_oscal = oscal_avx;
      cpu_oaxpy_b32 = oaxpy_b32_avx;
    }
#endif
  }
#endif
}
