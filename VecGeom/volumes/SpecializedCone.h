/*
 * SpecializedCone.h
 *
 *  Created on: May 14, 2014
 *      Author: swenzel
 */

#ifndef VECGEOM_VOLUMES_SPECIALIZEDCONE_H_
#define VECGEOM_VOLUMES_SPECIALIZEDCONE_H_

#include "VecGeom/base/Global.h"
#include "VecGeom/volumes/kernel/ConeImplementation.h"
#include "VecGeom/volumes/PlacedCone.h"
#include "VecGeom/volumes/SpecializedPlacedVolImplHelper.h"

namespace vecgeom {
inline namespace VECGEOM_IMPL_NAMESPACE {

template <TranslationCode transCodeT, RotationCode rotCodeT, typename coneTypeT>
using SpecializedCone = SpecializedVolImplHelper<ConeImplementation<coneTypeT>, transCodeT, rotCodeT>;

using SimpleCone = SpecializedCone<translation::kGeneric, rotation::kGeneric, ConeTypes::UniversalCone>;

} // namespace VECGEOM_IMPL_NAMESPACE
} // namespace vecgeom

#endif // VECGEOM_VOLUMES_SPECIALIZEDCONE_H_
