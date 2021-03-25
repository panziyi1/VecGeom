/// \file PlacedHype.cpp
/// \author Marilena Bandieramonte (marilena.bandieramonte@cern.ch)

#include "VecGeom/volumes/PlacedHype.h"
#include "VecGeom/volumes/SpecializedHype.h"
#include "VecGeom/volumes/Hype.h"
#include "VecGeom/base/Global.h"

//#if defined(VECGEOM_BENCHMARK) && defined(VECGEOM_ROOT)
#ifdef VECGEOM_ROOT
#include "TGeoHype.h"
#endif

#ifdef VECGEOM_GEANT4
#include "G4Hype.hh"
#endif

namespace vecgeom {
inline namespace VECGEOM_IMPL_NAMESPACE {

#ifndef VECCORE_CUDA

#ifdef VECGEOM_ROOT
TGeoShape const *PlacedHype::ConvertToRoot() const
{
  return GetUnplacedVolume()->ConvertToRoot(GetName());
}
#endif

#ifdef VECGEOM_GEANT4
G4VSolid const *PlacedHype::ConvertToGeant4() const
{
  return GetUnplacedVolume()->ConvertToGeant4(GetName());
}
#endif

#endif // VECGEOM_BENCHMARK

} // namespace VECGEOM_IMPL_NAMESPACE

#ifdef VECCORE_CUDA
VECGEOM_DEVICE_INST_PLACED_VOLUME_ALLSPEC_3(SpecializedHype, HypeTypes::UniversalHype)

#endif // VECCORE_CUDA
} // End namespace vecgeom
