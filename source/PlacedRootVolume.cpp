/// \file PlacedRootVolume.cpp
/// \author Johannes de Fine Licht (johannes.definelicht@cern.ch)

#include "volumes/PlacedRootVolume.h"

#include "base/AOS3D.h"
#include "base/SOA3D.h"

namespace vecgeom {

PlacedRootVolume::PlacedRootVolume(char const *const label,
                                   TGeoShape const *const rootShape,
                                   LogicalVolume const *const logicalVolume,
                                   Transformation3D const *const transformation)
    : VPlacedVolume(label, logicalVolume, transformation, NULL),
      fRootShape(rootShape) {}

PlacedRootVolume::PlacedRootVolume(TGeoShape const *const rootShape,
                                   LogicalVolume const *const logicalVolume,
                                   Transformation3D const *const transformation)
    : PlacedRootVolume(rootShape->GetName(), rootShape, logicalVolume,
                       transformation) {}

void PlacedRootVolume::PrintType() const {
  printf("PlacedRootVolume");
}

void PlacedRootVolume::Contains(SOA3D<Precision> const &points,
                                bool *const output) const {
  for (int i = 0; i < points.size(); ++i) output[i] = PlacedRootVolume::Contains(points[i]);
}

void PlacedRootVolume::Contains(AOS3D<Precision> const &points,
                                bool *const output) const {
  for (int i = 0; i < points.size(); ++i) output[i] = PlacedRootVolume::Contains(points[i]);
}

void PlacedRootVolume::Inside(SOA3D<Precision> const &points,
                              Inside_t *const output) const {
  for (int i = 0; i < points.size(); ++i) output[i] = PlacedRootVolume::Inside(points[i]);
}

void PlacedRootVolume::Inside(AOS3D<Precision> const &points,
                              Inside_t *const output) const {
  for (int i = 0; i < points.size(); ++i) output[i] = PlacedRootVolume::Inside(points[i]);
}

void PlacedRootVolume::DistanceToIn(SOA3D<Precision> const &position,
                                    SOA3D<Precision> const &direction,
                                    Precision const *const step_max,
                                    Precision *const output) const {
  for (int i = 0, iMax = position.size(); i < iMax; ++i) {
    output[i] = PlacedRootVolume::DistanceToIn(position[i], direction[i], step_max[i]);
  }
}

void PlacedRootVolume::DistanceToIn(AOS3D<Precision> const &position,
                                    AOS3D<Precision> const &direction,
                                    Precision const *const step_max,
                                    Precision *const output) const {
  for (int i = 0, iMax = position.size(); i < iMax; ++i) {
    output[i] = PlacedRootVolume::DistanceToIn(position[i], direction[i], step_max[i]);
  }
}

void PlacedRootVolume::DistanceToOut(SOA3D<Precision> const &position,
                                     SOA3D<Precision> const &direction,
                                     Precision const *const step_max,
                                     Precision *const output) const {
  for (int i = 0, iMax = position.size(); i < iMax; ++i) {
    output[i] = PlacedRootVolume::DistanceToOut(position[i], direction[i], step_max[i]);
  }
}

void PlacedRootVolume::DistanceToOut(AOS3D<Precision> const &position,
                                     AOS3D<Precision> const &direction,
                                     Precision const *const step_max,
                                     Precision *const output) const {
  for (int i = 0, iMax = position.size(); i < iMax; ++i) {
    output[i] = PlacedRootVolume::DistanceToOut(position[i], direction[i], step_max[i]);
  }
}

void PlacedRootVolume::SafetyToIn(SOA3D<Precision> const &position,
                                  Precision *const safeties) const  {
  for (int i = 0, iMax = position.size(); i < iMax; ++i) {
    safeties[i] = PlacedRootVolume::SafetyToIn(position[i]);
  }
}

void PlacedRootVolume::SafetyToIn(AOS3D<Precision> const &position,
                                  Precision *const safeties) const  {
  for (int i = 0, iMax = position.size(); i < iMax; ++i) {
    safeties[i] = PlacedRootVolume::SafetyToIn(position[i]);
  }
}

void PlacedRootVolume::SafetyToInMinimize(SOA3D<Precision> const &position,
                                          Precision *const safeties) const  {
  for (int i = 0, iMax = position.size(); i < iMax; ++i) {
    Precision result = PlacedRootVolume::SafetyToIn(position[i]);
    safeties[i] = (result < safeties[i]) ? result : safeties[i];
  }
}

void PlacedRootVolume::SafetyToOut(SOA3D<Precision> const &position,
                                   Precision *const safeties) const {
  for (int i = 0, iMax = position.size(); i < iMax; ++i) {
    safeties[i] = PlacedRootVolume::SafetyToOut(position[i]);
  }
}

void PlacedRootVolume::SafetyToOut(AOS3D<Precision> const &position,
                                   Precision *const safeties) const  {
  for (int i = 0, iMax = position.size(); i < iMax; ++i) {
    safeties[i] = PlacedRootVolume::SafetyToOut(position[i]);
  }
}

void PlacedRootVolume::SafetyToOutMinimize(SOA3D<Precision> const &position,
                                           Precision *const safeties) const {
  for (int i = 0, iMax = position.size(); i < iMax; ++i) {
    Precision result = PlacedRootVolume::SafetyToOut(position[i]);
    safeties[i] = (result < safeties[i]) ? result : safeties[i];
  }
}


void PlacedRootVolume::DistanceToInMinimize(SOA3D<Precision> const &positions,
                                            SOA3D<Precision> const &directions,
                                            int daughterindex,
                                            Precision *step,
                                            int * nextnodeids) const {
   for (int i = 0, iMax = positions.size(); i < iMax; ++i) {
     Precision result = PlacedRootVolume::DistanceToIn(positions[i], directions[i], step[i]);
     if(result < step[i]) {
           step[i]=result;
           nextnodeids[i]=daughterindex;
     }
   }
}

void PlacedRootVolume::DistanceToOut(SOA3D<Precision> const &positions,
                                     SOA3D<Precision> const &directions,
                                     Precision const * step_max,
                                     Precision * distance,
                                     int * nextnodeindex) const {
    for (int i = 0, iMax = positions.size(); i < iMax; ++i) {
        distance[i] = PlacedRootVolume::DistanceToOut(positions[i], directions[i], step_max[i]);
        nextnodeindex[i] = (distance[i] < step_max[i] ) ? -1 : -2;
      }
}


#ifdef VECGEOM_BENCHMARK
VPlacedVolume const* PlacedRootVolume::ConvertToUnspecialized() const {
  assert(0 && "Attempted to perform conversion on unsupported ROOT volume.");
  return NULL;
}
#ifdef VECGEOM_ROOT
TGeoShape const* PlacedRootVolume::ConvertToRoot() const {
  assert(0 && "Attempted to perform conversion on unsupported ROOT volume.");
  return NULL;
}
#endif
#ifdef VECGEOM_USOLIDS
::VUSolid const* PlacedRootVolume::ConvertToUSolids() const {
  assert(0 && "Attempted to perform conversion on unsupported ROOT volume.");
  return NULL;
}
#endif
#endif // VECGEOM_BENCHMARK

#ifdef VECGEOM_CUDA_INTERFACE
VPlacedVolume* PlacedRootVolume::CopyToGpu(
    LogicalVolume const *const logical_volume,
    Transformation3D const *const transformation,
    VPlacedVolume *const gpu_ptr) const {
  assert(0 && "Attempted to copy unsupported ROOT volume to GPU.");
}
VPlacedVolume* PlacedRootVolume::CopyToGpu(
    LogicalVolume const *const logical_volume,
    Transformation3D const *const transformation) const {
  assert(0 && "Attempted to copy unsupported ROOT volume to GPU.");
}
#endif

} // End namespace vecgeom
