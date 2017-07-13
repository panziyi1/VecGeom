/// \file PlacedTorus.cpp
#include "volumes/PlacedTorus.h"
#include "volumes/Torus.h"
#include "volumes/SpecializedTorus.h"
#include "base/Vector3D.h"

#ifdef VECGEOM_ROOT
#include "TGeoTorus.h"
#endif

#ifdef VECGEOM_GEANT4
#include "G4Torus.hh"
#endif

#ifdef VECGEOM_USOLIDS
#include "UTorus.hh"
#endif

namespace vecgeom {
inline namespace VECGEOM_IMPL_NAMESPACE {

#ifndef VECCORE_CUDA

VPlacedVolume const *PlacedTorus::ConvertToUnspecialized() const
{
  return new SimpleTorus(GetLabel().c_str(), GetLogicalVolume(), GetTransformation());
}

#ifdef VECGEOM_ROOT
TGeoShape const *PlacedTorus::ConvertToRoot() const
{
  const UnplacedTorus &ut = *(static_cast<UnplacedTorus const *>(GetUnplacedVolume()));
  return new TGeoTorus(GetLabel().c_str(), ut.rtor(), ut.rmin(), ut.rmax(), ut.sphi() * kRadToDeg,
                       ut.dphi() * kRadToDeg);
}
#endif

#ifdef VECGEOM_GEANT4
G4VSolid const *PlacedTorus::ConvertToGeant4() const
{
  const UnplacedTorus &ut = *(static_cast<UnplacedTorus const *>(GetUnplacedVolume()));
  return new G4Torus(GetLabel().c_str(), ut.rmin(), ut.rmax(), ut.rtor(), ut.sphi(), ut.dphi());
}
#endif

#endif // VECCORE_CUDA

} // End impl namespace

#ifdef VECCORE_CUDA

VECGEOM_DEVICE_INST_PLACED_VOLUME_ALLSPEC(SpecializedTorus)

#endif

} // namespace vecgeom
