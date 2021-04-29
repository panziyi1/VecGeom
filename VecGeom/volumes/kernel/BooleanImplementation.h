/*
 * BooleanImplementation.h
 */

#ifndef BOOLEANIMPLEMENTATION_H_
#define BOOLEANIMPLEMENTATION_H_

#include "VecGeom/base/Vector3D.h"
#include "VecGeom/volumes/BooleanStruct.h"
#include <VecCore/VecCore>

namespace vecgeom {

VECGEOM_DEVICE_DECLARE_CONV_TEMPLATE_2v(struct, BooleanImplementation, BooleanOperation, Arg1, typename, Arg2);

/// Prototype for functions dispatching to matching PlacedVolume virtual interface
#define DEFINE_FUNC_DISPATCH_TO_PLACED_VOLUME(func)                                                             \
  template <typename... Args>                                                                                   \
  static VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE auto func(VPlacedVolume const *vol, const Args &... args) \
      ->decltype(vol->func(args...))                                                                            \
  {                                                                                                             \
    return vol->func(args...);                                                                                  \
  }

#define DEFINE_FUNC_DISPATCH_TO_PLACED_VOLUME_OUT(func)                                        \
  static VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE auto func(                               \
      VPlacedVolume const *vol, Vector3D<Precision> const &input, Vector3D<Precision> &output) \
      ->decltype(vol->func(input, output))                                                     \
  {                                                                                            \
    return vol->func(input, output);                                                           \
  }

struct DispatchVirtual {
  DEFINE_FUNC_DISPATCH_TO_PLACED_VOLUME(Contains);
  DEFINE_FUNC_DISPATCH_TO_PLACED_VOLUME(Inside);
  DEFINE_FUNC_DISPATCH_TO_PLACED_VOLUME(DistanceToIn);
  DEFINE_FUNC_DISPATCH_TO_PLACED_VOLUME(DistanceToOut);
  DEFINE_FUNC_DISPATCH_TO_PLACED_VOLUME(PlacedDistanceToOut);
  DEFINE_FUNC_DISPATCH_TO_PLACED_VOLUME(SafetyToIn);
  DEFINE_FUNC_DISPATCH_TO_PLACED_VOLUME(SafetyToOut);
  DEFINE_FUNC_DISPATCH_TO_PLACED_VOLUME_OUT(Normal);
};

inline namespace VECGEOM_IMPL_NAMESPACE {

template <BooleanOperation Op>
class PlacedBooleanVolume;
template <BooleanOperation Op>
class UnplacedBooleanVolume;

/// Boolean implementation, supporting user-defined dispatch to PlaceVolume components
template <BooleanOperation boolOp, typename Dispatch = DispatchVirtual>
struct BooleanImplementation {
  using PlacedShape_t    = PlacedBooleanVolume<boolOp>;
  using UnplacedVolume_t = UnplacedBooleanVolume<boolOp>;
  using UnplacedStruct_t = BooleanStruct;

  // empty since functionality will be implemented in
  // partially template specialized structs
};

/**
 * an ordinary (non-templated) implementation of a Boolean solid
 * using the virtual function interface of its constituents
 *
 * TEMPLATE SPECIALIZATION FOR SUBTRACTION
 */
template <typename Dispatch>
struct BooleanImplementation<kSubtraction, Dispatch> {
  using PlacedShape_t    = PlacedBooleanVolume<kSubtraction>;
  using UnplacedVolume_t = UnplacedBooleanVolume<kSubtraction>;
  using UnplacedStruct_t = BooleanStruct;

  VECCORE_ATT_HOST_DEVICE
  static void PrintType()
  { /*printf("SpecializedBooleanVolume<%i, %i, %i>", kSubtraction, transCodeT, rotCodeT);*/
  }

  template <typename Stream>
  static void PrintType(Stream &s)
  {
    /*s << "SpecializedBooleanVolume<kSubtraction"
      << "," << transCodeT << "," << rotCodeT << ">";
    */
  }

  template <typename Stream>
  static void PrintType(Stream &st, int transCodeT = translation::kGeneric, int rotCodeT = rotation::kGeneric)
  {
    st << "SpecializedBooleanVolume<kSubtraction" << transCodeT << "," << rotCodeT << ">";
  }

  template <typename Stream>
  static void PrintImplementationType(Stream &s)
  {
    /*
    s << "BooleanImplementation<kSubtraction"
      << "," << transCodeT << "," << rotCodeT << ">";
    */
  }

  template <typename Stream>
  static void PrintUnplacedType(Stream &s)
  {
    s << "BooleanStruct";
  }

  template <typename Real_v, typename Bool_v>
  VECGEOM_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  static void Contains(BooleanStruct const &unplaced, Vector3D<Real_v> const &point, Bool_v &inside)
  {
    Vector3D<Real_v> tmp;
    inside = Dispatch::Contains(unplaced.fLeftVolume, point);
    if (vecCore::MaskEmpty(inside)) return;

    auto rightInside = Dispatch::Contains(unplaced.fRightVolume, point);
    inside &= !rightInside;
  }

  template <typename Real_v, typename Inside_t>
  VECGEOM_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  static void Inside(BooleanStruct const &unplaced, Vector3D<Real_v> const &p, Inside_t &inside)
  {

    // now use the Inside functionality of left and right components
    // algorithm taken from Geant4 implementation
    VPlacedVolume const *const fPtrSolidA = unplaced.fLeftVolume;
    VPlacedVolume const *const fPtrSolidB = unplaced.fRightVolume;

    const auto positionA = Dispatch::Inside(fPtrSolidA, p);
    if (positionA == EInside::kOutside) {
      inside = EInside::kOutside;
      return;
    }

    const auto positionB = Dispatch::Inside(fPtrSolidB, p);

    if (positionA == EInside::kInside && positionB == EInside::kOutside) {
      inside = EInside::kInside;
      return;
    } else {
      if ((positionA == EInside::kInside && positionB == EInside::kSurface) ||
          (positionB == EInside::kOutside && positionA == EInside::kSurface)
          /*
           ||( positionA == EInside::kSurface && positionB == EInside::kSurface &&
             (   fPtrSolidA->Normal(p) -
               fPtrSolidB->Normal(p) ).mag2() >
             1000.0*G4GeometryTolerance::GetInstance()->GetRadialTolerance() ) )
          */) {
        inside = EInside::kSurface;
        return;
      } else {
        inside = EInside::kOutside;
        return;
      }
    }
    // going to be a bit more complicated due to Surface states
  }

  template <typename Real_v>
  VECGEOM_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  static void DistanceToIn(BooleanStruct const &unplaced, Vector3D<Real_v> const &p, Vector3D<Real_v> const &dir,
                           Real_v const &stepMax, Real_v &distance)
  {

    // TOBEDONE: ASK Andrei about the while loop
    // Compute distance from a given point outside to the shape.
    Real_v d1, d2, snxt = 0.;
    Vector3D<Real_v> hitpoint = p;
    // check if inside '-'
    auto insideRight = Dispatch::Contains(unplaced.fRightVolume, p);

    //  // epsilon is used to push across boundaries
    Precision epsil(1E-6);
    //
    //  // we should never subtract a volume such that B - A > 0
    //
    while (1) {
      if (insideRight) {
        //    // propagate to outside of '- / RightShape'
        d1 = Dispatch::PlacedDistanceToOut(unplaced.fRightVolume, hitpoint, dir, stepMax);
        snxt += (d1 >= 0. && d1 < kInfLength) ? (d1 + epsil) : 0.;
        hitpoint += (d1 >= 0. && d1 < kInfLength) ? (d1 + epsil) * dir : 0. * dir;

        // now master outside 'B'; check if inside 'A'
        if (Dispatch::Contains(unplaced.fLeftVolume, hitpoint)) {
          auto check = Dispatch::PlacedDistanceToOut(unplaced.fLeftVolume, hitpoint, dir);
          if (check > epsil) {
            distance = snxt;
            //	std::cerr << "hitting  " << distance << "\n";
            return;
          }
        }
      }

      // if outside of both we do a max operation
      // master outside '-' and outside '+' ;  find distances to both
      //        fLeftMat->MasterToLocal(&master[0], &local[0]);
      d2 = Dispatch::DistanceToIn(unplaced.fLeftVolume, hitpoint, dir, stepMax);
      d2 = Max(d2, Precision(0.));
      if (d2 == kInfLength) {
        distance = kInfLength;
        // std::cerr << "missing A " << d2 << "\n";
        return;
      }

      d1 = Dispatch::DistanceToIn(unplaced.fRightVolume, hitpoint, dir, stepMax);
      if (d2 < d1 - kTolerance) {
        snxt += d2 + epsil;
        // std::cerr << "returning  " << snxt << "\n";
        distance = snxt;
        return;
      }

      //        // propagate to '-'
      snxt += (d1 >= 0. && d1 < kInfLength) ? d1 + epsil : 0.;
      hitpoint += (d1 >= 0. && d1 < kInfLength) ? (d1 + epsil) * dir : epsil * dir;
      insideRight = true;
    } // end while
  }

  template <typename Real_v>
  VECGEOM_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  static void DistanceToOut(BooleanStruct const &unplaced, Vector3D<Real_v> const &point,
                            Vector3D<Real_v> const &direction, Real_v const &stepMax, Real_v &distance)
  {
    const auto distancel  = Dispatch::PlacedDistanceToOut(unplaced.fLeftVolume, point, direction, stepMax);
    const Real_v dinright = Dispatch::DistanceToIn(unplaced.fRightVolume, point, direction, stepMax);
    distance              = Min(distancel, dinright);
    return;
  }

  template <typename Real_v>
  VECGEOM_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  static void SafetyToIn(BooleanStruct const &unplaced, Vector3D<Real_v> const &point, Real_v &safety)
  {
    VPlacedVolume const *const fPtrSolidA = unplaced.fLeftVolume;
    VPlacedVolume const *const fPtrSolidB = unplaced.fRightVolume;

    // very approximate
    if (Dispatch::Contains(fPtrSolidA, point) && // case 1
        Dispatch::Contains(fPtrSolidB, point)) {
      safety = Dispatch::SafetyToOut(fPtrSolidB, fPtrSolidB->GetTransformation()->Transform(point));
    } else {
      // po
      safety = Dispatch::SafetyToIn(fPtrSolidA, point);
    }
  }

  template <typename Real_v>
  VECGEOM_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  static void SafetyToOut(BooleanStruct const &unplaced, Vector3D<Real_v> const &point, Real_v &safety)
  {
    const auto safetyleft  = Dispatch::SafetyToOut(unplaced.fLeftVolume, point);
    const auto safetyright = Dispatch::SafetyToIn(unplaced.fRightVolume, point);
    safety                 = Min(safetyleft, safetyright);
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

    // If point is inside B, then it must be on a surface of B
    if (Dispatch::Contains(fPtrSolidB, point)) {
      fPtrSolidB->GetTransformation()->Transform(point, localPoint);
      valid = Dispatch::Normal(fPtrSolidB, localPoint, localNorm);
      // The normal to the subtracted solid has to be inverted and transformed back
      localNorm *= -1.;
      fPtrSolidB->GetTransformation()->InverseTransformDirection(localNorm, normal);
      return;
    }

    // If point is outside A, then it must be on a surface of A
    if (Dispatch::Contains(fPtrSolidA, point)) {
      fPtrSolidA->GetTransformation()->Transform(point, localPoint);
      valid = Dispatch::Normal(fPtrSolidA, localPoint, localNorm);
      fPtrSolidA->GetTransformation()->InverseTransformDirection(localNorm, normal);
      return;
    }

    // Point is inside A and outside B, check safety
    fPtrSolidA->GetTransformation()->Transform(point, localPoint);
    Real_v safetyA = Dispatch::SafetyToOut(fPtrSolidA, localPoint);
    Real_v safetyB = Dispatch::SafetyToIn(fPtrSolidB, point);
    Bool_v onA     = safetyA < safetyB;
    if (vecCore::MaskFull(onA)) {
      valid = Dispatch::Normal(fPtrSolidA, localPoint, localNorm);
      fPtrSolidA->GetTransformation()->InverseTransformDirection(localNorm, normal);
      return;
    } else {
      //  if (vecCore::MaskEmpty(onA)) {  // to use real mask operation when supporting vectors
      fPtrSolidB->GetTransformation()->Transform(point, localPoint);
      valid = Dispatch::Normal(fPtrSolidB, localPoint, localNorm);
      // The normal to the subtracted solid has to be inverted and transformed back
      localNorm *= -1.;
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

// include stuff for boolean union
#include "BooleanUnionImplementation.h"

// include stuff for boolean intersection
#include "BooleanIntersectionImplementation.h"

#endif /* BOOLEANIMPLEMENTATION_H_ */
