/// \file GenericKernels.h
/// \author Johannes de Fine Licht (johannes.definelicht@cern.ch)

#ifndef VECGEOM_VOLUMES_KERNEL_GENERICKERNELS_H_
#define VECGEOM_VOLUMES_KERNEL_GENERICKERNELS_H_

#include "base/Global.h"
#include "base/Transformation3D.h"
#include "base/Vector3D.h"

namespace VECGEOM_NAMESPACE {

template <class Backend>
struct GenericKernels {

  typedef typename Backend::precision_v Float_t;
  typedef typename Backend::int_v Int_t;
  typedef typename Backend::bool_v Bool_t;

}; // End struct GenericKernels

template<bool tolerant, typename T>
VECGEOM_CUDA_HEADER_BOTH
VECGEOM_INLINE
T MakePlusTolerant(T const &x)
{
  return (tolerant)? x+kSHalfTolerance*10 : x;
}

template<bool tolerant, typename T>
VECGEOM_CUDA_HEADER_BOTH
VECGEOM_INLINE
T MakeMinusTolerant(T const &x)
{
  return (tolerant)? x-kSHalfTolerance*10 : x;
}

template<bool tolerant, typename T>
VECGEOM_CUDA_HEADER_BOTH
VECGEOM_INLINE
T  MakePlusTolerantSquare(T const &x, T const& xsq)
{
  return (tolerant)? xsq+kSTolerance*10*x : xsq;
}

template<bool tolerant, typename T>
VECGEOM_CUDA_HEADER_BOTH
VECGEOM_INLINE
T MakeMinusTolerantSquare(T const &x, T const &xsq)
{
  return (tolerant)? xsq-kSTolerance*10*x : xsq;
}


//Replicated definitions with modified Tolerances appropriate for Sphere
template<bool tolerant, typename T>
VECGEOM_CUDA_HEADER_BOTH
VECGEOM_INLINE
T MakeSPlusTolerant(T const &x)
{
  return (tolerant)? x+kSHalfTolerance : x;
}

template<bool tolerant, typename T>
VECGEOM_CUDA_HEADER_BOTH
VECGEOM_INLINE
T MakeSMinusTolerant(T const &x)
{
  return (tolerant)? x-kSHalfTolerance : x;
}

template<bool tolerant, typename T>
VECGEOM_CUDA_HEADER_BOTH
VECGEOM_INLINE
T  MakeSPlusTolerantSquare(T const &x, T const& xsq)
{
  return (tolerant)? xsq+kSTolerance*x : xsq;
}

template<bool tolerant, typename T>
VECGEOM_CUDA_HEADER_BOTH
VECGEOM_INLINE
T MakeSMinusTolerantSquare(T const &x, T const &xsq)
{
  return (tolerant)? xsq-kSTolerance*x : xsq;
}

} // End global namespace

#endif // VECGEOM_VOLUMES_KERNEL_GENERICKERNELS_H_
