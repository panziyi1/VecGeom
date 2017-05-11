#ifndef VECCORE_SIMD_H
#define VECCORE_SIMD_H

#if defined(SSE) 
#include <x86intrin.h>
#else
#include <arm_neon.h>
#endif

#if defined(__MIC__)
#define VECCORE_SIMD_ALIGN 64
#elif defined(__AVX__)
#define VECCORE_SIMD_ALIGN 32
#else
#define VECCORE_SIMD_ALIGN 16
#endif

#endif
