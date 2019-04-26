/// \file SpecializedTorus.h

#ifndef VECGEOM_VOLUMES_SPECIALIZEDTORUS_H_
#define VECGEOM_VOLUMES_SPECIALIZEDTORUS_H_

#include "base/Global.h"

#include "volumes/kernel/TorusImplementation.h"
#include "volumes/PlacedTorus.h"
#include "volumes/SpecializedPlacedVolImplHelper.h"

namespace vecgeom {
inline namespace VECGEOM_IMPL_NAMESPACE {

// NOTE: we may want to specialize the torus like we do for the tube
// at the moment this is not done

template <TranslationCode transCodeT, RotationCode rotCodeT>
  //using SpecializedTorus = SIMDSpecializedVolImplHelper<TorusImplementation/*<transCodeT, rotCodeT>*/, transCodeT, rotCodeT>;
  using SpecializedTorus = LoopSpecializedVolImplHelper<TorusImplementation/*<transCodeT, rotCodeT>*/, transCodeT, rotCodeT>;

using SimpleTorus = SpecializedTorus<translation::kGeneric, rotation::kGeneric>;

// template <TranslationCode transCodeT, RotationCode rotCodeT>
// VECCORE_ATT_DEVICE
// VPlacedVolume *UnplacedTorus::Create(LogicalVolume const *const logical_volume,
//                                     Transformation3D const *const transformation,
// #ifdef VECCORE_CUDA
//                                     const int id,
// #endif
//                                     VPlacedVolume *const placement)
// {
//   (void)placement;
// #ifdef VECCORE_CUDA
//   return new SimpleTorus(logical_volume, transformation, id);
// #else
//   return new SimpleTorus(logical_volume, transformation);
// #endif
// }

} // namespace VECGEOM_IMPL_NAMESPACE
} // namespace vecgeom

#endif // VECGEOM_VOLUMES_SPECIALIZEDTORUS_H_
