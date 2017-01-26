/// \file PlacedPolyhedron.h
/// \author Johannes de Fine Licht (johannes.definelicht@cern.ch)

#ifndef VECGEOM_VOLUMES_PLACEDPOLYHEDRON_H_
#define VECGEOM_VOLUMES_PLACEDPOLYHEDRON_H_

#include "base/Global.h"

#include "volumes/PlacedVolume.h"
#include "volumes/kernel/PolyhedronImplementation.h"
#include "volumes/PlacedVolImplHelper.h"
#include "volumes/UnplacedPolyhedron.h"

class UPolyhedraHistorical;
struct PolyhedraSideRZ {   // Avoid clash in class name UPolyhedraSideRZ
  vecgeom::Precision r, z; // start of vector
};

namespace vecgeom {

VECGEOM_DEVICE_FORWARD_DECLARE(class PlacedPolyhedron;);
VECGEOM_DEVICE_DECLARE_CONV(class, PlacedPolyhedron);

inline namespace VECGEOM_IMPL_NAMESPACE {

class PlacedPolyhedron : public PlacedVolumeImplHelper<UnplacedPolyhedron, VPlacedVolume> {
  using Base = PlacedVolumeImplHelper<UnplacedPolyhedron, VPlacedVolume>;

public:
  typedef UnplacedPolyhedron UnplacedShape_t;

#ifndef VECGEOM_NVCC
  // constructor inheritance;
  using Base::Base;
  PlacedPolyhedron(char const *const label, LogicalVolume const *const logicalVolume,
                   Transformation3D const *const transformation, vecgeom::PlacedBox const *const boundingBox)
      : Base(label, logicalVolume, transformation, boundingBox)
  {
  }

  PlacedPolyhedron(LogicalVolume const *const logicalVolume, Transformation3D const *const transformation,
                   vecgeom::PlacedBox const *const boundingBox)
      : PlacedPolyhedron("", logicalVolume, transformation, boundingBox)
  {
  }
#else
  __device__ PlacedPolyhedron(LogicalVolume const *const logicalVolume, Transformation3D const *const transformation,
                              PlacedBox const *const boundingBox, const int id)
      : Base(logicalVolume, transformation, boundingBox, id)
  {
  }
#endif

  VECGEOM_CUDA_HEADER_BOTH
  virtual ~PlacedPolyhedron() {}

  // Accessors

  VECGEOM_CUDA_HEADER_BOTH
  UnplacedPolyhedron const *GetUnplacedVolume() const
  {
    return static_cast<UnplacedPolyhedron const *>(GetLogicalVolume()->GetUnplacedVolume());
  }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_FORCE_INLINE
  int GetSideCount() const { return GetUnplacedVolume()->GetSideCount(); }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_FORCE_INLINE
  int GetZSegmentCount() const { return GetUnplacedVolume()->GetZSegmentCount(); }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_FORCE_INLINE
  bool HasInnerRadii() const { return GetUnplacedVolume()->HasInnerRadii(); }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_FORCE_INLINE
  bool HasPhiCutout() const { return GetUnplacedVolume()->HasPhiCutout(); }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_FORCE_INLINE
  ZSegment const &GetZSegment(int index) const { return GetUnplacedVolume()->GetZSegment(index); }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_FORCE_INLINE
  Array<ZSegment> const &GetZSegments() const { return GetUnplacedVolume()->GetZSegments(); }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_FORCE_INLINE
  Array<Precision> const &GetZPlanes() const { return GetUnplacedVolume()->GetZPlanes(); }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_FORCE_INLINE
  Array<Precision> const &GetRMin() const { return GetUnplacedVolume()->GetRMin(); }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_FORCE_INLINE
  Array<Precision> const &GetRMax() const { return GetUnplacedVolume()->GetRMax(); }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_FORCE_INLINE
  Vector3D<Precision> GetPhiSection(int i) const { return GetUnplacedVolume()->GetPhiSection(i); }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_FORCE_INLINE
  SOA3D<Precision> const &GetPhiSections() const { return GetUnplacedVolume()->GetPhiSections(); }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_FORCE_INLINE
  Precision GetPhiStart() const { return GetUnplacedVolume()->GetPhiStart(); }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_FORCE_INLINE
  Precision GetPhiEnd() const { return GetUnplacedVolume()->GetPhiEnd(); }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_FORCE_INLINE
  Precision GetPhiDelta() const { return GetUnplacedVolume()->GetPhiDelta(); }

#ifndef VECGEOM_NVCC
  virtual std::ostream &StreamInfo(std::ostream &os) const { return GetUnplacedVolume()->StreamInfo(os); }

  virtual Precision Capacity() override { return GetUnplacedVolume()->Capacity(); }

  virtual Precision SurfaceArea() override { return GetUnplacedVolume()->SurfaceArea(); }

  void Extent(Vector3D<Precision> &aMin, Vector3D<Precision> &aMax) const override
  {
    GetUnplacedVolume()->Extent(aMin, aMax);
  }

  virtual Vector3D<Precision> GetPointOnSurface() const override { return GetUnplacedVolume()->GetPointOnSurface(); }

#if defined(VECGEOM_USOLIDS)
  std::string GetEntityType() const override { return GetUnplacedVolume()->GetEntityType(); }
#endif
#endif

  VECGEOM_CUDA_HEADER_BOTH
  int PhiSegmentIndex(Vector3D<Precision> const &point) const;

  bool IsOpen() const { return (GetUnplacedVolume()->GetPhiDelta() < kTwoPi); }
  bool IsGeneric() const { return false; }
  Precision GetStartPhi() const { return GetUnplacedVolume()->GetPhiStart(); }
  Precision GetEndPhi() const { return GetUnplacedVolume()->GetPhiEnd(); }
  int GetNumSide() const { return GetUnplacedVolume()->GetSideCount(); }
  int GetNumRZCorner() const { return GetZPlanes().size(); }

  // CUDA specific

  virtual int memory_size() const override { return sizeof(*this); }

  VECGEOM_CUDA_HEADER_BOTH
  virtual bool Normal(Vector3D<Precision> const &point, Vector3D<Precision> &normal) const override
  {
    return GetUnplacedVolume()->Normal(point, normal);
  }

// Comparison specific
#ifndef VECGEOM_NVCC
  virtual VPlacedVolume const *ConvertToUnspecialized() const override;
#ifdef VECGEOM_ROOT
  virtual TGeoShape const *ConvertToRoot() const override;
#endif
#if defined(VECGEOM_USOLIDS) && !defined(VECGEOM_REPLACE_USOLIDS)
  virtual ::VUSolid const *ConvertToUSolids() const override;
#endif
#ifdef VECGEOM_GEANT4
  virtual G4VSolid const *ConvertToGeant4() const override;
#endif
#endif // VECGEOM_NVCC
};

} // End inline namespace

} // End global namespace

#endif // VECGEOM_VOLUMES_PLACEDPOLYHEDRON_H_
