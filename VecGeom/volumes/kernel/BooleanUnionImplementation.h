/*
 * BooleanImplementation.h
 */

#ifndef BOOLEANUNIONIMPLEMENTATION_H_
#define BOOLEANUNIONIMPLEMENTATION_H_

#include "VecGeom/base/Global.h"
#include "VecGeom/base/Vector3D.h"
#include "VecGeom/volumes/BooleanStruct.h"

#ifdef VECGEOM_ROOT
#include "TGeoShape.h"
#endif

namespace vecgeom {

inline namespace VECGEOM_IMPL_NAMESPACE {

/**
 * partial template specialization for UNION implementation
 */
template <typename Dispatch>
struct BooleanImplementation<kUnion, Dispatch> {
  using PlacedShape_t    = PlacedBooleanVolume<kUnion>;
  using UnplacedVolume_t = UnplacedBooleanVolume<kUnion>;
  using UnplacedStruct_t = BooleanStruct;

  VECCORE_ATT_HOST_DEVICE
  static void PrintType()
  { /* printf("SpecializedBooleanVolume<%i, %i, %i>", kUnion, transCodeT, rotCodeT); */
  }

  template <typename Stream>
  static void PrintType(Stream &st, int transCodeT = translation::kGeneric, int rotCodeT = rotation::kGeneric)
  {
    st << "SpecializedBooleanVolume<kUnion" << transCodeT << "," << rotCodeT << ">";
  }

  template <typename Stream>
  static void PrintType(Stream &s)
  {
    //  s << "SpecializedBooleanVolume<kUnion"
    //    << "," << transCodeT << "," << rotCodeT << ">";
  }

  template <typename Stream>
  static void PrintImplementationType(Stream &s)
  {
    // s << "BooleanImplementation<kUnion"
    //   << "," << transCodeT << "," << rotCodeT << ">";
  }

  template <typename Stream>
  static void PrintUnplacedType(Stream &s)
  {
    s << "UnplacedBooleanVolume";
  }

  template <typename Real_v, typename Bool_v>
  VECGEOM_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  static void Contains(BooleanStruct const &unplaced, Vector3D<Real_v> const &point, Bool_v &inside)
  {
    inside = Dispatch::Contains(unplaced.fLeftVolume, point);
    if (vecCore::MaskFull(inside)) return;
    inside |= Dispatch::Contains(unplaced.fRightVolume, point);
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
    if (positionA == EInside::kInside) {
      inside = EInside::kInside;
      return;
    }

    const auto positionB = Dispatch::Inside(fPtrSolidB, point);
    if (positionB == EInside::kInside) {
      inside = EInside::kInside;
      return;
    }

    if ((positionA == EInside::kSurface) && (positionB == EInside::kSurface)) {
      Vector3D<Precision> normalA, normalB, localPoint, localNorm;
      fPtrSolidA->GetTransformation()->Transform(point, localPoint);
      Dispatch::Normal(fPtrSolidA, localPoint, localNorm);
      fPtrSolidA->GetTransformation()->InverseTransformDirection(localNorm, normalA);

      fPtrSolidB->GetTransformation()->Transform(point, localPoint);
      Dispatch::Normal(fPtrSolidB, localPoint, localNorm);
      fPtrSolidB->GetTransformation()->InverseTransformDirection(localNorm, normalB);

      if (normalA.Dot(normalB) < 0)
        inside = EInside::kInside; // touching solids -)(-
      else
        inside = EInside::kSurface; // overlapping solids =))
      return;
    } else {
      if ((positionB == EInside::kSurface) || (positionA == EInside::kSurface)) {
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
  static void DistanceToIn(BooleanStruct const &unplaced, Vector3D<Real_v> const &point,
                           Vector3D<Real_v> const &direction, Real_v const &stepMax, Real_v &distance)
  {
    const auto d1 = Dispatch::DistanceToIn(unplaced.fLeftVolume, point, direction, stepMax);
    const auto d2 = Dispatch::DistanceToIn(unplaced.fRightVolume, point, direction, stepMax);
    distance      = Min(d1, d2);
  }

  template <typename Real_v>
  VECGEOM_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  static void DistanceToOut(BooleanStruct const &unplaced, Vector3D<Real_v> const &point, Vector3D<Real_v> const &dir,
                            Real_v const &stepMax, Real_v &distance)
  {
    VPlacedVolume const *const ptrSolidA = unplaced.fLeftVolume;
    VPlacedVolume const *const ptrSolidB = unplaced.fRightVolume;

    Real_v dist = 0.;
    Real_v pushdist(1E-6);
    // size_t push          = 0;
    const auto positionA = Dispatch::Inside(ptrSolidA, point);
    Vector3D<Real_v> nextp(point);
    bool connectingstep(false);

    // reusable kernel as lambda
    auto kernel = [&](VPlacedVolume const *A, VPlacedVolume const *B) {
      do {
        connectingstep    = false;
        const auto disTmp = Dispatch::PlacedDistanceToOut(A, nextp, dir);
        dist += (disTmp >= 0. && disTmp < kInfLength) ? disTmp : 0;
        // give a push
        dist += pushdist;
        // push++;
        nextp = point + dist * dir;
        // B could be overlapping with A -- and/or connecting A to another part of A
        // if (B->Contains(nextp)) {
        if (Dispatch::Inside(B, nextp) != vecgeom::kOutside) {
          const auto disTmp = Dispatch::PlacedDistanceToOut(B, nextp, dir);
          dist += (disTmp >= 0. && disTmp < kInfLength) ? disTmp : 0;
          dist += pushdist;
          // push++;
          nextp          = point + dist * dir;
          connectingstep = true;
        }
      } while (connectingstep && (Dispatch::Inside(A, nextp) != kOutside));
    };

    if (positionA != kOutside) { // initially in A
      kernel(ptrSolidA, ptrSolidB);
    }
    // if( positionB != kOutside )
    else {
      kernel(ptrSolidB, ptrSolidA);
    }
    // At the end we need to subtract just one push distance, since intermediate distances
    // from pushed points are smaller than the real distance with the push value
    distance = dist - pushdist;
    if (distance < kTolerance && positionA == kOutside && Dispatch::Inside(ptrSolidB, point) == kOutside)
      distance = -kTolerance;
    return;
  }

  template <typename Real_v>
  VECGEOM_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  static void SafetyToIn(BooleanStruct const &unplaced, Vector3D<Real_v> const &point, Real_v &safety)
  {
    VPlacedVolume const *const fPtrSolidA = unplaced.fLeftVolume;
    VPlacedVolume const *const fPtrSolidB = unplaced.fRightVolume;
    const auto distA                      = Dispatch::SafetyToIn(fPtrSolidA, point);
    const auto distB                      = Dispatch::SafetyToIn(fPtrSolidB, point);
    safety                                = Min(distA, distB);
    // If safety is negative it should not be made 0 (convention)
    // vecCore::MaskedAssign(safety, safety < 0.0, 0.0);
  }

  template <typename Real_v>
  VECGEOM_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  static void SafetyToOut(BooleanStruct const &unplaced, Vector3D<Real_v> const &point, Real_v &safety)
  {

    safety                                = -kTolerance; // invalid side
    VPlacedVolume const *const fPtrSolidA = unplaced.fLeftVolume;
    VPlacedVolume const *const fPtrSolidB = unplaced.fRightVolume;

    const auto insideA = Dispatch::Inside(fPtrSolidA, point);
    const auto insideB = Dispatch::Inside(fPtrSolidB, point);

    // Is point already outside?
    if (insideA == kOutside && insideB == kOutside) return;

    if (insideA != kOutside && insideB != kOutside) /* in both */
    {
      safety = Max(Dispatch::SafetyToOut(fPtrSolidA, point),
                   Dispatch::SafetyToOut(fPtrSolidB, fPtrSolidB->GetTransformation()->Transform(point)));
    } else {
      if (insideA == kSurface || insideB == kSurface) return;
      /* only contained in B */
      if (insideA == kOutside) {
        safety = Dispatch::SafetyToOut(fPtrSolidB, fPtrSolidB->GetTransformation()->Transform(point));
      } else {
        safety = Dispatch::SafetyToOut(fPtrSolidA, point);
      }
    }
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

    // If point is inside A, then it must be on a surface of A (points on the
    // intersection between A and B cannot be on surface, or if they are they
    // are on a common surface and the normal can be computer for A or B)
    if (Dispatch::Contains(fPtrSolidA, point)) {
      fPtrSolidA->GetTransformation()->Transform(point, localPoint);
      valid = Dispatch::Normal(fPtrSolidA, localPoint, localNorm);
      fPtrSolidA->GetTransformation()->InverseTransformDirection(localNorm, normal);
      return;
    }
    // Same for points inside B
    if (Dispatch::Contains(fPtrSolidB, point)) {
      fPtrSolidB->GetTransformation()->Transform(point, localPoint);
      valid = Dispatch::Normal(fPtrSolidB, localPoint, localNorm);
      fPtrSolidB->GetTransformation()->InverseTransformDirection(localNorm, normal);
      return;
    }
    // Points outside both A and B can be on any surface. We use the safety.
    const auto safetyA = Dispatch::SafetyToIn(fPtrSolidA, point);
    const auto safetyB = Dispatch::SafetyToIn(fPtrSolidB, point);
    auto onA           = safetyA < safetyB;
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
