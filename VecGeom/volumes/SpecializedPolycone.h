/*
 * SpecializedPolycone.h
 *
 *  Created on: Dec 8, 2014
 *      Author: swenzel
 */

#ifndef VECGEOM_VOLUMES_SPECIALIZEDPOLYCONE_H_
#define VECGEOM_VOLUMES_SPECIALIZEDPOLYCONE_H_

#include "VecGeom/base/Global.h"

#include "VecGeom/volumes/kernel/PolyconeImplementation.h"
#include "VecGeom/volumes/PlacedPolycone.h"
#include "VecGeom/volumes/SpecializedPlacedVolImplHelper.h"
#include "VecGeom/volumes/UnplacedPolycone.h"

#include <stdio.h>

namespace vecgeom {
inline namespace VECGEOM_IMPL_NAMESPACE {

template <TranslationCode transCodeT, RotationCode rotCodeT, typename polyconeTypeT>
using SpecializedPolycone = SpecializedVolImplHelper<PolyconeImplementation<polyconeTypeT>, transCodeT, rotCodeT>;

using SimplePolycone = SpecializedPolycone<translation::kGeneric, rotation::kGeneric, ConeTypes::UniversalCone>;

} // namespace VECGEOM_IMPL_NAMESPACE
} // namespace vecgeom

#endif /* VECGEOM_VOLUMES_SPECIALIZEDPOLYCONE_H_ */
