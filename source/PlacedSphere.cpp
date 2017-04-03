/// @file PlacedSphere.cpp
/// @author Raman Sehgal (raman.sehgal@cern.ch)

#include "volumes/PlacedSphere.h"
#include "volumes/Sphere.h"

#ifdef VECGEOM_USOLIDS
#include "USphere.hh"
#endif

#ifdef VECGEOM_ROOT
#include "TGeoSphere.h"
#endif

#ifdef VECGEOM_GEANT4
#include "G4Sphere.hh"
#endif

namespace vecgeom {
inline namespace VECGEOM_IMPL_NAMESPACE {

#ifndef VECCORE_CUDA

VPlacedVolume const *PlacedSphere::ConvertToUnspecialized() const
{
  return new SimpleSphere(GetLabel().c_str(), GetLogicalVolume(), GetTransformation());
}

#ifdef VECGEOM_ROOT
TGeoShape const *PlacedSphere::ConvertToRoot() const
{
  return new TGeoSphere(GetLabel().c_str(), GetInnerRadius(), GetOuterRadius(), GetStartThetaAngle() * kRadToDeg,
                        (GetStartThetaAngle() + GetDeltaThetaAngle()) * kRadToDeg, GetStartPhiAngle() * kRadToDeg,
                        (GetStartPhiAngle() + GetDeltaPhiAngle()) * kRadToDeg);
}
#endif

#if defined(VECGEOM_USOLIDS) && !defined(VECGEOM_REPLACE_USOLIDS)
::VUSolid const *PlacedSphere::ConvertToUSolids() const
{

  return new USphere(GetLabel().c_str(), GetInnerRadius(), GetOuterRadius(), GetStartPhiAngle(), GetDeltaPhiAngle(),
                     GetStartThetaAngle(), GetDeltaThetaAngle());
}
#endif

#ifdef VECGEOM_GEANT4
G4VSolid const *PlacedSphere::ConvertToGeant4() const
{
  return new G4Sphere(GetLabel().c_str(), GetInnerRadius(), GetOuterRadius(), GetStartPhiAngle(), GetDeltaPhiAngle(),
                      GetStartThetaAngle(), GetDeltaThetaAngle());
}
#endif

#endif // VECCORE_CUDA

} // End impl namespace

#ifdef VECCORE_CUDA

VECGEOM_DEVICE_INST_PLACED_VOLUME_ALLSPEC(SpecializedSphere)

#endif // VECCORE_CUDA

} // End global namespace
