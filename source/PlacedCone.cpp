/*
 * PlacedCone.cpp
 *
 *  Created on: Jun 13, 2014
 *      Author: swenzel
 */

#include "VecGeom/volumes/PlacedCone.h"
#include "VecGeom/volumes/Cone.h"

#if defined(VECGEOM_ROOT)
#include "TGeoCone.h"
#endif

#if defined(VECGEOM_GEANT4)
#include "G4Cons.hh"
#endif

namespace vecgeom {
inline namespace VECGEOM_IMPL_NAMESPACE {

VECCORE_ATT_HOST_DEVICE
void PlacedCone::PrintType() const
{
  printf("PlacedCone");
}

void PlacedCone::PrintType(std::ostream &s) const
{
  s << "PlacedCone";
}

#ifndef VECCORE_CUDA

#ifdef VECGEOM_ROOT
TGeoShape const *PlacedCone::ConvertToRoot() const
{
  return GetUnplacedVolume()->ConvertToRoot(GetName());
}
#endif

#ifdef VECGEOM_GEANT4
G4VSolid const *PlacedCone::ConvertToGeant4() const
{
  return GetUnplacedVolume()->ConvertToGeant4(GetName());
}
#endif

#endif

} // namespace VECGEOM_IMPL_NAMESPACE

#ifdef VECCORE_CUDA

VECGEOM_DEVICE_INST_PLACED_VOLUME(PlacedCone)

#endif // VECCORE_CUDA

} // End namespace vecgeom
