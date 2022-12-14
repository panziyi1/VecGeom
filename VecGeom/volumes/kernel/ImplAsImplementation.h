#ifndef VECGEOM_VOLUMES_KERNEL_INDIRECTIMPLEMENTATION_H_
#define VECGEOM_VOLUMES_KERNEL_INDIRECTIMPLEMENTATION_H_

#include "VecGeom/base/Vector3D.h"
#include <VecCore/VecCore>
#include <VecGeom/volumes/PlacedImplAs.h>

namespace vecgeom {

VECGEOM_DEVICE_DECLARE_CONV_TEMPLATE_2t(struct, IndirectImplementation, typename, typename);

inline namespace VECGEOM_IMPL_NAMESPACE {

// structure that provides (indirect) implementation kernels
// for a certain unplaced volume which is implemented in terms of
// another one
template <typename UnplVol_t, typename DispatchingImplementation>
struct IndirectImplementation {

  using PlacedShape_t    = PlacedImplAs<UnplVol_t>; // typename DispatchingImplementation::PlacedShape_t;
  using UnplacedStruct_t = typename DispatchingImplementation::UnplacedStruct_t;
  using UnplacedVolume_t = UnplVol_t;

  VECCORE_ATT_HOST_DEVICE
  static void PrintType()
  {
    //  printf("SpecializedOrb<%i, %i>", transCodeT, rotCodeT);
  }

  template <typename Stream>
  static void PrintType(Stream &st, int transCodeT = translation::kGeneric, int rotCodeT = rotation::kGeneric)
  {
    DispatchingImplementation::PrintType(st, transCodeT, rotCodeT);
  }

  template <typename Stream>
  static void PrintImplementationType(Stream &st)
  {
    (void)st;
    // st << "IndirectImplementation<" << transCodeT << "," << rotCodeT << ">";
  }

  template <typename Stream>
  static void PrintUnplacedType(Stream &st)
  {
    (void)st;
  }

  template <typename Real_v, typename Bool_v>
  VECGEOM_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  static void Contains(UnplacedStruct_t const &s, Vector3D<Real_v> const &point, Bool_v &inside)
  {
    DispatchingImplementation::template Contains(s, point, inside);
  }

  template <typename Real_v, typename Inside_t>
  VECGEOM_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  static void Inside(UnplacedStruct_t const &s, Vector3D<Real_v> const &point, Inside_t &inside)
  {
    DispatchingImplementation::template Inside(s, point, inside);
  }

  template <typename Real_v>
  VECGEOM_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  static void DistanceToIn(UnplacedStruct_t const &s, Vector3D<Real_v> const &point, Vector3D<Real_v> const &direction,
                           Real_v const &stepMax, Real_v &distance)
  {
    DispatchingImplementation::template DistanceToIn(s, point, direction, stepMax, distance);
  }

  template <typename Real_v>
  VECGEOM_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  static void DistanceToOut(UnplacedStruct_t const &s, Vector3D<Real_v> const &point, Vector3D<Real_v> const &direction,
                            Real_v const &stepMax, Real_v &distance)
  {
    DispatchingImplementation::template DistanceToOut(s, point, direction, stepMax, distance);
  }

  template <typename Real_v>
  VECGEOM_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  static void SafetyToIn(UnplacedStruct_t const &s, Vector3D<Real_v> const &point, Real_v &safety)
  {
    DispatchingImplementation::template SafetyToIn(s, point, safety);
  }

  template <typename Real_v>
  VECGEOM_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  static void SafetyToOut(UnplacedStruct_t const &s, Vector3D<Real_v> const &point, Real_v &safety)
  {
    DispatchingImplementation::template SafetyToOut(s, point, safety);
  }

  template <typename Real_v>
  VECGEOM_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  static Vector3D<Real_v> NormalKernel(UnplacedStruct_t const &s, Vector3D<Real_v> const &point,
                                       typename vecCore::Mask_v<Real_v> &valid)
  {
    DispatchingImplementation::template NormalKernel(s, point, valid);
  }
};

} // end implementation namespace
} // end global vecgeom namespace

#endif // VECGEOM_VOLUMES_KERNEL_INDIRECTIMPLEMENTATION_H_
