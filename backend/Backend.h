/// \file vector/Backend.h
/// \author Johannes de Fine Licht (johannes.definelicht@cern.ch)

#ifndef VECGEOM_BACKEND_BACKEND_H_
#define VECGEOM_BACKEND_BACKEND_H_

#include "base/Global.h"

#ifdef VECGEOM_NVCC
#include "backend/cuda/Backend.h"
#elif defined(VECGEOM_VC)
#include "backend/vc/Backend.h"
#elif defined(VECGEOM_CILK)
#include "backend/cilk/Backend.h"
#elif defined(VECGEOM_MICVEC)
#include "backend/micvec/Backend.h"
#elif defined(VECGEOM_UMESIMD)
#include "backend/umesimd/Backend.h"
#else
#include "backend/scalar/Backend.h"
#endif

namespace vecgeom {

inline namespace VECGEOM_IMPL_NAMESPACE {

template <typename T>
VECGEOM_FORCE_INLINE
VECGEOM_CUDA_HEADER_BOTH
T NonZeroAbs(T const &x)
{
  return vecCore::math::Abs(x) + vecCore::NumericLimits<T>::Min();
}

template <typename T>
VECGEOM_FORCE_INLINE
VECGEOM_CUDA_HEADER_BOTH
T NonZero(T const &x)
{
  const auto tiny = vecCore::NumericLimits<T>::Min();
  return x + vecCore::Blend(x >= 0, tiny, -tiny);
}
}
}

#endif // VECGEOM_BACKEND_BACKEND_H_
