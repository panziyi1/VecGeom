#ifndef VECGEOM_VOLUMES_SPECIALIZEDTUBE_H_
#define VECGEOM_VOLUMES_SPECIALIZEDTUBE_H_

#include "VecGeom/base/Global.h"

#include "VecGeom/volumes/kernel/TubeImplementation.h"
#include "VecGeom/volumes/PlacedTube.h"
#include "VecGeom/volumes/SpecializedPlacedVolImplHelper.h"

namespace vecgeom {
inline namespace VECGEOM_IMPL_NAMESPACE {

template <TranslationCode transCodeT, RotationCode rotCodeT, typename tubeTypeT>
using SpecializedTube = SpecializedVolImplHelper<TubeImplementation<tubeTypeT>, transCodeT, rotCodeT>;

using SimpleTube = SpecializedTube<translation::kGeneric, rotation::kGeneric, TubeTypes::UniversalTube>;

} // namespace VECGEOM_IMPL_NAMESPACE
} // namespace vecgeom

#endif // VECGEOM_VOLUMES_SPECIALIZEDTUBE_H_
