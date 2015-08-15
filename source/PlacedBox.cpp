/// \file PlacedBox.cpp
/// \author Johannes de Fine Licht (johannes.definelicht@cern.ch)

#ifdef OFFLOAD_MODE
  #pragma offload_attribute(push,target(mic))
  #include <map>
#endif

#include "volumes/PlacedBox.h"

#include "base/AOS3D.h"
#include "base/SOA3D.h"
#include "volumes/SpecializedBox.h"
#ifdef VECGEOM_ROOT
#include "TGeoBBox.h"
#endif
#ifdef VECGEOM_USOLIDS
#include "UBox.hh"
#endif
#ifdef VECGEOM_GEANT4
#include "G4Box.hh"
#endif

#include <stdio.h>

namespace vecgeom {
inline namespace VECGEOM_IMPL_NAMESPACE {

VECGEOM_CUDA_HEADER_BOTH
void PlacedBox::PrintType() const {
  printf("PlacedBox");
}

#ifndef VECGEOM_NVCC

VPlacedVolume const* PlacedBox::ConvertToUnspecialized() const {
  return new SimpleBox(GetLabel().c_str(), logical_volume_, transformation_);
}

#ifdef VECGEOM_ROOT
TGeoShape const* PlacedBox::ConvertToRoot() const {
  return new TGeoBBox(GetLabel().c_str(), x(), y(), z());
}
#endif

#ifdef VECGEOM_USOLIDS
::VUSolid const* PlacedBox::ConvertToUSolids() const {
  return new UBox(GetLabel(), x(), y(), z());
}
#endif

#ifdef VECGEOM_GEANT4
G4VSolid const* PlacedBox::ConvertToGeant4() const {
  return new G4Box("", x(), y(), z());
}
#endif

#endif // VECGEOM_NVCC

#ifdef OFFLOAD_MODE

static
std::map<size_t, size_t> _vplaced_volumes;

size_t PlacedBox::CopyToXeonPhi() const {
  size_t addr = size_t(this);
  size_t ret;
  auto it = _vplaced_volumes.find(addr);
  if(it == _vplaced_volumes.end()) {
    size_t logical_volume = logical_volume_->CopyToXeonPhi();
    size_t transf = transformation_->CopyToXeonPhi();
#pragma offload target(mic) out(ret) in(addr,logical_volume,transf) nocopy(_vplaced_volumes)
{
    LogicalVolume *lv = (LogicalVolume*)logical_volume;
    VPlacedVolume *vpv = lv->Place((Transformation3D const *const)transf);
    _vplaced_volumes[addr] = size_t(vpv);
    ret = size_t(vpv);
}
    _vplaced_volumes[addr] = ret;
  }
  return _vplaced_volumes[addr];
}

#endif

} // End impl namespace

#ifdef VECGEOM_NVCC

VECGEOM_DEVICE_INST_PLACED_VOLUME_ALLSPEC( SpecializedBox )

#endif // VECGEOM_NVCC

} // End namespace vecgeom

#ifdef OFFLOAD_MODE
  #pragma offload_attribute(pop)
#endif
