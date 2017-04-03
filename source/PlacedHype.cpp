/// \file PlacedHype.cpp
/// \author Marilena Bandieramonte (marilena.bandieramonte@cern.ch)

#include "volumes/PlacedHype.h"
#include "volumes/Hype.h"
#include "base/Global.h"

//#if defined(VECGEOM_BENCHMARK) && defined(VECGEOM_ROOT)
#ifdef VECGEOM_ROOT
#include "TGeoHype.h"
#endif

#ifdef VECGEOM_USOLIDS
#include "UBox.hh"
#endif

#ifdef VECGEOM_GEANT4
#include "G4Hype.hh"
#endif

namespace vecgeom {
inline namespace VECGEOM_IMPL_NAMESPACE {

#ifndef VECCORE_CUDA

VPlacedVolume const *PlacedHype::ConvertToUnspecialized() const
{
  std::cout << "Convert VEC*********\n";
  return new SimpleHype(GetLabel().c_str(), GetLogicalVolume(), GetTransformation());
}

#ifdef VECGEOM_ROOT
TGeoShape const *PlacedHype::ConvertToRoot() const
{
  std::cout << "Convert ROOT*********\n";
  return new TGeoHype(GetLabel().c_str(), GetRmin(), GetStIn() * kRadToDeg, GetRmax(), GetStOut() * kRadToDeg, GetDz());
}
#endif

#if defined(VECGEOM_USOLIDS) && !defined(VECGEOM_REPLACE_USOLIDS)
::VUSolid const *PlacedHype::ConvertToUSolids() const
{
  // assert(0 && "Hype unsupported for USolids.");
  // return NULL;
  // std::cerr << "**************************************************************\n";
  // std::cerr << "WARNING: Hyperboloid unsupported for USolids.; returning a box\n";
  // std::cerr << "**************************************************************\n";
  // return new UBox("",10,10,10);
  return NULL;
}
#endif

#ifdef VECGEOM_GEANT4
G4VSolid const *PlacedHype::ConvertToGeant4() const
{
  return new G4Hype(GetLabel().c_str(), GetRmin(), GetRmax(), GetStIn(), GetStOut(), GetDz());
}
#endif

#endif // VECGEOM_BENCHMARK

} // End global namespace

#ifdef VECCORE_CUDA

VECGEOM_DEVICE_INST_PLACED_VOLUME_ALLSPEC(SpecializedHype)

#endif // VECCORE_CUDA
} // End namespace vecgeom
