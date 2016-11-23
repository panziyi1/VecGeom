#ifndef VECGEOM_VOLUMES_SPECIALIZEDCUTTUBE_H_
#define VECGEOM_VOLUMES_SPECIALIZEDCUTTUBE_H_

#include "base/Global.h"

#include "volumes/kernel/CutTubeImplementation.h"
#include "volumes/PlacedCutTube.h"
#include "volumes/SpecializedPlacedVolImplHelper.h"
#include "volumes/UnplacedCutTube.h"

namespace vecgeom {
inline namespace VECGEOM_IMPL_NAMESPACE {

template <TranslationCode transCodeT, RotationCode rotCodeT>
using SpecializedCutTube = SIMDSpecializedVolImplHelper<CutTubeImplementation, transCodeT, rotCodeT>;

using SimpleCutTube = SpecializedCutTube<translation::kGeneric, rotation::kGeneric>;
}
} // End global namespace

#endif
