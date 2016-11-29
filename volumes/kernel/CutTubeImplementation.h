/*
 * CutTubeImplementation.h
 *
 *  Created on: 03.11.2016
 *      Author: mgheata
 */

#ifndef VECGEOM_VOLUMES_KERNEL_CUTTUBEIMPLEMENTATION_H_
#define VECGEOM_VOLUMES_KERNEL_CUTTUBEIMPLEMENTATION_H_

#include <cstdio>

#include "base/Vector3D.h"
#include "volumes/kernel/GenericKernels.h"
#include "volumes/CutTubeStruct.h"
#include "TubeImplementation.h"

namespace vecgeom {

VECGEOM_DEVICE_FORWARD_DECLARE(struct CutTubeImplementation;);
VECGEOM_DEVICE_DECLARE_CONV(struct, CutTubeImplementation);

inline namespace VECGEOM_IMPL_NAMESPACE {

class PlacedCutTube;
class UnplacedCutTube;

template <typename T>
struct CutTubeStruct;

struct CutTubeImplementation {

  using PlacedShape_t    = PlacedCutTube;
  using UnplacedStruct_t = CutTubeStruct<double>;
  using UnplacedVolume_t = UnplacedCutTube;

  VECGEOM_CUDA_HEADER_BOTH
  static void PrintType() {}

  template <typename Stream>
  static void PrintType(Stream &s, int transCodeT = translation::kGeneric, int rotCodeT = rotation::kGeneric)
  {
    s << "SpecializedCutTube<" << transCodeT << "," << rotCodeT << ">";
  }

  template <typename Stream>
  static void PrintImplementationType(Stream & /*s*/)
  {
  }

  template <typename Stream>
  static void PrintUnplacedType(Stream & /*s*/)
  {
  }

  template <typename Real_v, typename Bool_v>
  VECGEOM_FORCE_INLINE
  VECGEOM_CUDA_HEADER_BOTH
  static void Contains(UnplacedStruct_t const &unplaced, Vector3D<Real_v> const &point, Bool_v &inside);

  template <typename Real_v, typename Inside_t>
  VECGEOM_FORCE_INLINE
  VECGEOM_CUDA_HEADER_BOTH
  static void Inside(UnplacedStruct_t const &unplaced, Vector3D<Real_v> const &point, Inside_t &inside);

  template <typename Real_v>
  VECGEOM_FORCE_INLINE
  VECGEOM_CUDA_HEADER_BOTH
  static void DistanceToIn(UnplacedStruct_t const &unplaced, Vector3D<Real_v> const &point,
                           Vector3D<Real_v> const &direction, Real_v const &stepMax, Real_v &distance);

  template <typename Real_v>
  VECGEOM_FORCE_INLINE
  VECGEOM_CUDA_HEADER_BOTH
  static void DistanceToOut(UnplacedStruct_t const &unplaced, Vector3D<Real_v> const &point,
                            Vector3D<Real_v> const &direction, Real_v const &stepMax, Real_v &distance);

  template <typename Real_v>
  VECGEOM_FORCE_INLINE
  VECGEOM_CUDA_HEADER_BOTH
  static void SafetyToIn(UnplacedStruct_t const &unplaced, Vector3D<Real_v> const &point, Real_v &safety);

  template <typename Real_v>
  VECGEOM_FORCE_INLINE
  VECGEOM_CUDA_HEADER_BOTH
  static void SafetyToOut(UnplacedStruct_t const &unplaced, Vector3D<Real_v> const &point, Real_v &safety);

  template <typename Real_v, typename Bool_v>
  VECGEOM_FORCE_INLINE
  VECGEOM_CUDA_HEADER_BOTH
  static void NormalKernel(UnplacedStruct_t const &unplaced, Vector3D<Real_v> const &point, Vector3D<Real_v> &normal,
                           Bool_v &valid);
}; // End struct CutTubeImplementation

//********************************
//**** implementations start here
//********************************/

//______________________________________________________________________________
template <typename Real_v, typename Bool_v>
VECGEOM_FORCE_INLINE
VECGEOM_CUDA_HEADER_BOTH
void CutTubeImplementation::Contains(UnplacedStruct_t const &unplaced, Vector3D<Real_v> const &point, Bool_v &contains)
{
  contains                = Bool_v(false);
  Bool_v inside_cutplanes = Bool_v(false);
  // Check the cut planes first
  unplaced.GetCutPlanes().Contains<Real_v, Bool_v>(point, inside_cutplanes);
  if (vecCore::EarlyReturnAllowed()) {
    if (vecCore::MaskEmpty(inside_cutplanes)) return;
  }
  // Check the tube
  TubeImplementation<TubeTypes::UniversalTube>::Contains<Real_v>(unplaced.GetTubeStruct(), point, contains);
  contains &= inside_cutplanes;
}

//______________________________________________________________________________
template <typename Real_v, typename Inside_t>
VECGEOM_FORCE_INLINE
VECGEOM_CUDA_HEADER_BOTH
void CutTubeImplementation::Inside(UnplacedStruct_t const &unplaced, Vector3D<Real_v> const &point, Inside_t &inside)
{
  inside                    = Inside_t(EInside::kOutside);
  Inside_t inside_cutplanes = Inside_t(EInside::kOutside);
  // Check the cut planes first
  unplaced.GetCutPlanes().Inside<Real_v, Inside_t>(point, inside_cutplanes);
  if (vecCore::EarlyReturnAllowed()) {
    if (vecCore::MaskFull(inside_cutplanes == inside)) return;
  }
  // Check the tube
  TubeImplementation<TubeTypes::UniversalTube>::Inside<Real_v, Inside_t>(unplaced.GetTubeStruct(), point, inside);
  vecCore::MaskedAssign(inside, (inside_cutplanes == EInside::kOutside) ||
                                    (inside_cutplanes == EInside::kSurface && inside != EInside::kOutside),
                        inside_cutplanes);
}

//______________________________________________________________________________
template <typename Real_v>
VECGEOM_FORCE_INLINE
VECGEOM_CUDA_HEADER_BOTH
void CutTubeImplementation::DistanceToIn(UnplacedStruct_t const &unplaced, Vector3D<Real_v> const &point,
                                         Vector3D<Real_v> const &direction, Real_v const &stepMax, Real_v &distance)
{

#define USE_CONV_WRONG_SIDE 1
#define USE_CONV_FROM_BOUNDARY 1
  // Compute distance to cut planes
  using Bool_v = vecCore::Mask_v<Real_v>;
  distance     = InfinityLength<Real_v>();

#if USE_CONV_WRONG_SIDE == 1
  using Inside_v = vecCore::Index_v<Real_v>;
  // Check the cut planes first
  Inside_v inside_cutplanes;
  unplaced.GetCutPlanes().Inside<Real_v, Inside_v>(point, inside_cutplanes);
  Inside_v instart;
  // Check the tube
  TubeImplementation<TubeTypes::UniversalTube>::Inside<Real_v, Inside_v>(unplaced.GetTubeStruct(), point, instart);
  vecCore::MaskedAssign(instart, (inside_cutplanes == EInside::kOutside) ||
                                     (inside_cutplanes == EInside::kSurface && instart != EInside::kOutside),
                        inside_cutplanes);
  // Points already inside have to return negative distance
  vecCore::MaskedAssign(distance, instart == EInside::kInside, Real_v(-1.));
#endif

  Bool_v inside;
  Real_v dplanes;
  // Compute distance to cut planes
  unplaced.GetCutPlanes().DistanceToIn<Real_v>(point, direction, dplanes);
  // Mark tracks hitting the planes
  Bool_v hitplanes = dplanes < stepMax && dplanes > Real_v(-kTolerance);

#if USE_CONV_WRONG_SIDE == 1
  if (vecCore::EarlyReturnAllowed()) {
    if (vecCore::MaskFull((inside_cutplanes != EInside::kInside) && !hitplanes)) // No particles are hitting
      return;
  }
#endif

  // Propagate with dplanes only the particles that are hitting
  Vector3D<Real_v> propagated = point;
  vecCore::MaskedAssign(dplanes, !hitplanes, Real_v(0.));
  if (vecCore::EarlyReturnAllowed()) {
    if (!vecCore::MaskEmpty(hitplanes)) {
      propagated += dplanes * direction;
      // Hitting the planes does not guarantee that the hit point is between them
      unplaced.GetCutPlanes().Inside<Real_v, Inside_v>(propagated, inside_cutplanes);
    }
  } else {
    // In CUDA we have to propagate in all cases
    propagated += dplanes * direction;
    // Hitting the planes does not guarantee that the hit point is between them
    unplaced.GetCutPlanes().Inside<Real_v, Inside_v>(propagated, inside_cutplanes);
  }
  Bool_v done = inside_cutplanes == EInside::kOutside;
#if USE_CONV_WRONG_SIDE == 1
  done |= instart == EInside::kInside;
#endif
  // Tracks that cannot get between the cut planes cannot hit
  if (vecCore::EarlyReturnAllowed()) {
    if (vecCore::MaskFull(done)) // No particles are hitting
      return;
  }

  // All propagated points not yet marked as done should now lie between the cut planes
  // Check now which of the propagated points already entering the solid
  TubeImplementation<TubeTypes::UniversalTube>::Inside<Real_v>(unplaced.GetTubeStruct(), propagated, instart);
  inside = instart == EInside::kInside;
  vecCore::MaskedAssign(distance, inside && !done, dplanes);
  done |= inside;

  if (vecCore::EarlyReturnAllowed()) {
    if (vecCore::MaskFull(done)) { // Some particles hit top/bottom
      // The line below is needed for the convention
      vecCore::MaskedAssign(distance, vecCore::math::Abs(distance) < Real_v(kTolerance), Real_v(0.));
      return;
    }
  }

  // The limit distance for tube crossing cannot exceed the distance to
  // exiting the cut planes
  Real_v dexit = InfinityLength<Real_v>();
  unplaced.GetCutPlanes().DistanceToOut<Real_v>(propagated, direction, dexit);

  // Compute distance to tube
  Real_v dtube = InfinityLength<Real_v>();
  TubeImplementation<TubeTypes::UniversalTube>::DistanceToIn<Real_v>(unplaced.GetTubeStruct(), propagated, direction,
                                                                     stepMax, dtube);
  vecCore::MaskedAssign(dtube, dexit < dtube, InfinityLength<Real_v>());
  vecCore::MaskedAssign(distance, !done && (dtube + dplanes) < stepMax, dtube + dplanes);
// The line below is needed for the convention
#if USE_CONV_FROM_BOUNDARY == 1
//  vecCore::MaskedAssign(distance, vecCore::math::Abs(distance) < Real_v(kTolerance), Real_v(0.));
#endif
}

//______________________________________________________________________________
template <typename Real_v>
VECGEOM_FORCE_INLINE
VECGEOM_CUDA_HEADER_BOTH
void CutTubeImplementation::DistanceToOut(UnplacedStruct_t const &unplaced, Vector3D<Real_v> const &point,
                                          Vector3D<Real_v> const &direction, Real_v const &stepMax, Real_v &distance)
{
  // Compute distance to cut planes
  distance = InfinityLength<Real_v>();
  unplaced.GetCutPlanes().DistanceToOut<Real_v>(point, direction, distance);

  // Compute distance to tube
  Real_v dtube = InfinityLength<Real_v>();
  TubeImplementation<TubeTypes::UniversalTube>::DistanceToOut<Real_v>(unplaced.GetTubeStruct(), point, direction,
                                                                      stepMax, dtube);
  vecCore::MaskedAssign(distance, dtube < distance, dtube);
  // The line below is needed for the start on boundary convention
  //  vecCore::MaskedAssign(distance, vecCore::math::Abs(distance) < Real_v(kTolerance), Real_v(0.));
}

//______________________________________________________________________________
template <typename Real_v>
VECGEOM_FORCE_INLINE
VECGEOM_CUDA_HEADER_BOTH
void CutTubeImplementation::SafetyToIn(UnplacedStruct_t const &unplaced, Vector3D<Real_v> const &point, Real_v &safety)
{
  // Compute safety to cut planes. These will contain the sign, i.e. if on
  // the wrong side they will be negative
  unplaced.GetCutPlanes().SafetyToIn<Real_v>(point, safety);
  Real_v saftube;
  // Compute safety to tube
  TubeImplementation<TubeTypes::UniversalTube>::SafetyToIn<Real_v>(unplaced.GetTubeStruct(), point, saftube);
  // The safety is the maximum of the 2 values
  vecCore::MaskedAssign(safety, saftube > safety, saftube);
  // The line below is needed for the rounding convention
  //  vecCore::MaskedAssign(safety, vecCore::math::Abs(safety) < Real_v(kTolerance), Real_v(0.));
}

//______________________________________________________________________________
template <typename Real_v>
VECGEOM_CUDA_HEADER_BOTH
void CutTubeImplementation::SafetyToOut(UnplacedStruct_t const &unplaced, Vector3D<Real_v> const &point, Real_v &safety)
{
  // Compute safety to cut planes. These will contain the sign, i.e. if on
  // the wrong side they will be negative
  unplaced.GetCutPlanes().SafetyToOut<Real_v>(point, safety);
  Real_v saftube;
  // Compute safety to tube
  TubeImplementation<TubeTypes::UniversalTube>::SafetyToOut<Real_v>(unplaced.GetTubeStruct(), point, saftube);
  // The safety is the minimum of the 2 values
  vecCore::MaskedAssign(safety, saftube < safety, saftube);
  // The line below is needed for the rounding convention
  //  vecCore::MaskedAssign(safety, vecCore::math::Abs(safety) < Real_v(kTolerance), Real_v(0.));
}

//______________________________________________________________________________
template <typename Real_v, typename Bool_v>
VECGEOM_CUDA_HEADER_BOTH
void CutTubeImplementation::NormalKernel(UnplacedStruct_t const &unplaced, Vector3D<Real_v> const &point,
                                         Vector3D<Real_v> &normal, Bool_v &valid)
{
  // Compute safety to cut planes
  valid = Bool_v(true);
  Real_v safcut;
  unplaced.GetCutPlanes().SafetyToOut<Real_v>(point, safcut);
  // Compute safety to tube
  Real_v saftube;
  TubeImplementation<TubeTypes::UniversalTube>::SafetyToOut<Real_v>(unplaced.GetTubeStruct(), point, saftube);
  Bool_v istube = vecCore::math::Abs(saftube) < vecCore::math::Abs(safcut);
  // The statement below works only as long as Real_v is scalar
  if (istube) {
    TubeImplementation<TubeTypes::UniversalTube>::NormalKernel<Real_v, Bool_v>(unplaced.GetTubeStruct(), point, normal,
                                                                               valid);
    return;
  }
  // Select the correct cut plane
  if (point.z() < 0)
    normal = unplaced.fCutPlanes.GetNormal(0);
  else
    normal = unplaced.fCutPlanes.GetNormal(1);
}

//*****************************
//**** Implementations end here
//*****************************
}
} // End global namespace

#endif /* VECGEOM_VOLUMES_KERNEL_CUTTUBEIMPLEMENTATION_H_ */