/// @file SpecializedHype.h
/// @author Raman Sehgal (raman.sehgal@cern.ch)

#ifndef VECGEOM_VOLUMES_SPECIALIZEDHYPE_H_
#define VECGEOM_VOLUMES_SPECIALIZEDHYPE_H_

#include "VecGeom/base/Global.h"
#include "VecGeom/volumes/kernel/HypeImplementation.h"
#include "VecGeom/volumes/PlacedHype.h"
#include "VecGeom/volumes/SpecializedPlacedVolImplHelper.h"

#include <stdio.h>

namespace vecgeom {
inline namespace VECGEOM_IMPL_NAMESPACE {

template <TranslationCode transCodeT, RotationCode rotCodeT, typename hypeTypeT>
using SpecializedHype = SpecializedVolImplHelper<HypeImplementation<hypeTypeT>, transCodeT, rotCodeT>;

using SimpleHype = SpecializedHype<translation::kGeneric, rotation::kGeneric, HypeTypes::UniversalHype>;

} // namespace VECGEOM_IMPL_NAMESPACE
} // namespace vecgeom

#endif // VECGEOM_VOLUMES_SPECIALIZEDHYPE_H_
