// This file is part of VecGeom and is distributed under the
// conditions in the file LICENSE.txt in the top directory.
// For the full list of authors see CONTRIBUTORS.txt and `git log`.

/// Declaration of the specialized Trd volume
/// @file volumes/SpecializedTrd.h
/// @author Georgios Bitzes

#ifndef VECGEOM_VOLUMES_SPECIALIZEDTRD_H_
#define VECGEOM_VOLUMES_SPECIALIZEDTRD_H_

#include "VecGeom/base/Global.h"

#include "VecGeom/volumes/kernel/TrdImplementation.h"
#include "VecGeom/volumes/PlacedTrd.h"
#include "VecGeom/volumes/SpecializedPlacedVolImplHelper.h"

namespace vecgeom {
inline namespace VECGEOM_IMPL_NAMESPACE {

template <TranslationCode transCodeT, RotationCode rotCodeT, typename TrdTypeT>
using SpecializedTrd = SpecializedVolImplHelper<TrdImplementation<TrdTypeT>, transCodeT, rotCodeT>;

using SimpleTrd = SpecializedTrd<translation::kGeneric, rotation::kGeneric, TrdTypes::UniversalTrd>;

} // namespace VECGEOM_IMPL_NAMESPACE
} // namespace vecgeom

#endif // VECGEOM_VOLUMES_SPECIALIZEDTRD_H_
