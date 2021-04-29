/*
 * BooleanImplementation.h
 */

#ifndef BOOLEANINTERSECTIONIMPLEMENTATION_H_
#define BOOLEANINTERSECTIONIMPLEMENTATION_H_

#include "VecGeom/base/Global.h"
#include "VecGeom/base/Vector3D.h"
#include "VecGeom/volumes/BooleanStruct.h"
#include <VecCore/VecCore>

namespace vecgeom {

inline namespace VECGEOM_IMPL_NAMESPACE {

/**
 * partial template specialization for UNION implementation
 */
template <typename Dispatch>
struct BooleanImplementation<kIntersection, Dispatch> {
  using PlacedShape_t    = PlacedBooleanVolume<kIntersection>;
  using UnplacedVolume_t = UnplacedBooleanVolume<kIntersection>;
  using UnplacedStruct_t = BooleanStruct;

  VECCORE_ATT_HOST_DEVICE
  static void PrintType()
  { /* printf("SpecializedBooleanVolume<%i, %i, %i>", kIntersection, transCodeT, rotCodeT); */
  }

  template <typename Stream>
  static void PrintType(Stream &s)
  {
    // s << "SpecializedBooleanVolume<kIntersection"
    //  << "," << transCodeT << "," << rotCodeT << ">";
  }

  template <typename Stream>
  static void PrintType(Stream &st, int transCodeT = translation::kGeneric, int rotCodeT = rotation::kGeneric)
  {
    st << "SpecializedBooleanVolume<kIntersection" << transCodeT << "," << rotCodeT << ">";
  }

  template <typename Stream>
  static void PrintImplementationType(Stream &s)
  {
    // s << "BooleanImplementation<kIntersection"
    //  << "," << transCodeT << "," << rotCodeT << ">";
  }

  template <typename Stream>
  static void PrintUnplacedType(Stream &s)
  {
    // s << "UnplacedBooleanVolume";
  }

  template <typename Real_v, typename Bool_v>
  VECGEOM_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  static void Contains(BooleanStruct const &unplaced, Vector3D<Real_v> const &point, Bool_v &inside)
  {
    const auto insideA = Dispatch::Contains(unplaced.fLeftVolume, point);
    const auto insideB = Dispatch::Contains(unplaced.fRightVolume, point);
    inside             = insideA && insideB;
  }

  template <typename Real_v, typename Inside_t>
  VECGEOM_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  static void Inside(BooleanStruct const &unplaced, Vector3D<Real_v> const &point, Inside_t &inside)
  {
    // now use the Inside functionality of left and right components
    // algorithm taken from Geant4 implementation
    VPlacedVolume const *const fPtrSolidA = unplaced.fLeftVolume;
    VPlacedVolume const *const fPtrSolidB = unplaced.fRightVolume;

    const auto positionA = Dispatch::Inside(fPtrSolidA, point);

    if (positionA == EInside::kOutside) {
      inside = EInside::kOutside;
      return;
    }

    const auto positionB = Dispatch::Inside(fPtrSolidB, point);
    if (positionA == EInside::kInside && positionB == EInside::kInside) {
      inside = EInside::kInside;
      return;
    } else {
      if ((positionA == EInside::kInside && positionB == EInside::kSurface) ||
          (positionB == EInside::kInside && positionA == EInside::kSurface) ||
          (positionA == EInside::kSurface && positionB == EInside::kSurface)) {
        inside = EInside::kSurface;
        return;
      } else {
        inside = EInside::kOutside;
        return;
      }
    }
  }

  template <typename Real_v>
  VECGEOM_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  static void DistanceToIn(BooleanStruct const &unplaced, Vector3D<Real_v> const &point, Vector3D<Real_v> const &dir,
                           Real_v const &stepMax, Real_v &distance)
  {
    Vector3D<Real_v> hitpoint = point;

    auto inleft  = Dispatch::Contains(unplaced.fLeftVolume, hitpoint);
    auto inright = Dispatch::Contains(unplaced.fRightVolume, hitpoint);
    Real_v d1    = 0.;
    Real_v d2    = 0.;
    Real_v snext = 0.0;

    // just a pre-check before entering main algorithm
    if (inleft && inright) {
      d1 = Dispatch::PlacedDistanceToOut(unplaced.fLeftVolume, hitpoint, dir, stepMax);
      d2 = Dispatch::PlacedDistanceToOut(unplaced.fRightVolume, hitpoint, dir, stepMax);

      // if we are close to a boundary continue
      if (d1 < 2 * kTolerance) inleft = false;  // Backend::kFalse;
      if (d2 < 2 * kTolerance) inright = false; // Backend::kFalse;

      // otherwise exit
      if (inleft && inright) {
        // TODO: WE are inside both so should return a negative number
        distance = 0.0;
        return;
      }
    }

    // main loop
    while (1) {
      d1 = d2 = 0;
      if (!inleft) {
        d1 = Dispatch::DistanceToIn(unplaced.fLeftVolume, hitpoint, dir);
        d1 = Max(d1, kTolerance);
        if (d1 > 1E20) {
          distance = kInfLength;
          return;
        }
      }
      if (!inright) {
        d2 = Dispatch::DistanceToIn(unplaced.fRightVolume, hitpoint, dir);
        d2 = Max(d2, kTolerance);
        if (d2 > 1E20) {
          distance = kInfLength;
          return;
        }
      }

      if (d1 > d2) {
        // propagate to left shape
        snext += d1;
        inleft = true; // Backend::kTrue;
        hitpoint += d1 * dir;

        // check if propagated point is inside right shape
        // check is done with a little push
        inright = Dispatch::Contains(unplaced.fRightVolume, hitpoint + kTolerance * dir);
        if (inright) {
          distance = snext;
          return;
        }
        // here inleft=true, inright=false
      } else {
        // propagate to right shape
        snext += d2;
        inright = true; // Backend::kTrue;
        hitpoint += d2 * dir;

        // check if propagated point is inside left shape
        inleft = Dispatch::Contains(unplaced.fLeftVolume, hitpoint + kTolerance * dir);
        if (inleft) {
          distance = snext;
          return;
        }
      }
      // here inleft=false, inright=true
    } // end while loop
    distance = snext;
    return;
  }

  template <typename Real_v>
  VECGEOM_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  static void DistanceToOut(BooleanStruct const &unplaced, Vector3D<Real_v> const &point,
                            Vector3D<Real_v> const &direction, Real_v const &stepMax, Real_v &distance)
  {
    distance = Min(Dispatch::DistanceToOut(unplaced.fLeftVolume, point, direction),
                   Dispatch::DistanceToOut(unplaced.fRightVolume, point, direction));
  }

  template <typename Real_v>
  VECGEOM_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  static void SafetyToIn(BooleanStruct const &unplaced, Vector3D<Real_v> const &point, Real_v &safety)
  {
    // This is the Geant4 algorithm
    // TODO: ROOT seems to produce better safeties
    const auto insideA = Dispatch::Contains(unplaced.fLeftVolume, point);
    const auto insideB = Dispatch::Contains(unplaced.fRightVolume, point);

    if (!insideA && insideB) {
      safety = Dispatch::SafetyToIn(unplaced.fLeftVolume, point);
    } else {
      if (!insideB && insideA) {
        safety = Dispatch::SafetyToIn(unplaced.fRightVolume, point);
      } else {
        safety =
            Min(Dispatch::SafetyToIn(unplaced.fLeftVolume, point), Dispatch::SafetyToIn(unplaced.fRightVolume, point));
      }
    }
    return;
  }

  template <typename Real_v>
  VECGEOM_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  static void SafetyToOut(BooleanStruct const &unplaced, Vector3D<Real_v> const &point, Real_v &safety)
  {
    safety = Min(
        // TODO: could fail if left volume is placed shape
        Dispatch::SafetyToOut(unplaced.fLeftVolume, point),

        // TODO: consider introducing PlacedSafetyToOut function
        Dispatch::SafetyToOut(unplaced.fRightVolume, unplaced.fRightVolume->GetTransformation()->Transform(point)));
    vecCore::MaskedAssign(safety, safety < Real_v(0.), Real_v(0.));
  }

  template <typename Real_v, typename Bool_v>
  VECGEOM_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  static void NormalKernel(BooleanStruct const &unplaced, Vector3D<Real_v> const &point, Vector3D<Real_v> &normal,
                           Bool_v &valid)
  {
    Vector3D<Real_v> localNorm;
    Vector3D<Real_v> localPoint;
    valid = false; // Backend::kFalse;

    VPlacedVolume const *const fPtrSolidA = unplaced.fLeftVolume;
    VPlacedVolume const *const fPtrSolidB = unplaced.fRightVolume;
    Real_v safetyA, safetyB;

    if (Dispatch::Contains(fPtrSolidA, point)) {
      fPtrSolidA->GetTransformation()->Transform(point, localPoint);
      safetyA = Dispatch::SafetyToOut(fPtrSolidA, localPoint);
    } else {
      safetyA = Dispatch::SafetyToIn(fPtrSolidA, point);
    }

    if (Dispatch::Contains(fPtrSolidB, point)) {
      fPtrSolidB->GetTransformation()->Transform(point, localPoint);
      safetyB = Dispatch::SafetyToOut(fPtrSolidB, localPoint);
    } else {
      safetyB = Dispatch::SafetyToIn(fPtrSolidB, point);
    }
    const auto onA = safetyA < safetyB;
    if (vecCore::MaskFull(onA)) {
      fPtrSolidA->GetTransformation()->Transform(point, localPoint);
      valid = Dispatch::Normal(fPtrSolidA, localPoint, localNorm);
      fPtrSolidA->GetTransformation()->InverseTransformDirection(localNorm, normal);
      return;
    } else {
      //  if (vecCore::MaskEmpty(onA)) {  // to use real mask operation when supporting vectors
      fPtrSolidB->GetTransformation()->Transform(point, localPoint);
      valid = Dispatch::Normal(fPtrSolidB, localPoint, localNorm);
      fPtrSolidB->GetTransformation()->InverseTransformDirection(localNorm, normal);
      return;
    }
    // Some particles are on A, some on B. We never arrive here in the scalar case
    // If the interface to Normal will support the vector case, we have to write code here.
    return;
  }
}; // End struct BooleanImplementation

} // namespace VECGEOM_IMPL_NAMESPACE

} // namespace vecgeom

#endif /* BooleanImplementation_H_ */
