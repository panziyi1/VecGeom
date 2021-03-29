// This file is part of VecGeom and is distributed under the
// conditions in the file LICENSE.txt in the top directory.
// For the full list of authors see CONTRIBUTORS.txt and `git log`.

/// \brief Declaration of the unplaced Orb shape
/// \file volumes/UnplacedOrb.h
/// \author Raman Sehgal

#ifndef VECGEOM_VOLUMES_UNPLACEDORB_H_
#define VECGEOM_VOLUMES_UNPLACEDORB_H_

#include "VecGeom/base/Global.h"
#include "VecGeom/base/AlignedBase.h"
#include "VecGeom/base/Vector3D.h"
#include "VecGeom/volumes/UnplacedVolume.h"
#include "VecGeom/volumes/OrbStruct.h" // the pure Orb struct
#include "VecGeom/volumes/kernel/OrbImplementation.h"
#include "VecGeom/volumes/UnplacedVolumeImplHelper.h"
#ifdef VECGEOM_ROOT
class TGeoShape;
#endif
#ifdef VECGEOM_GEANT4
class G4VSolid;
#endif

namespace vecgeom {

VECGEOM_DEVICE_FORWARD_DECLARE(class UnplacedOrb;);
VECGEOM_DEVICE_DECLARE_CONV(class, UnplacedOrb);

inline namespace VECGEOM_IMPL_NAMESPACE {

class UnplacedOrb : public UnplacedVolumeImplHelper<OrbImplementation>, public AlignedBase {

private:
  OrbStruct<Precision> fOrb; ///< Structure holding the data for Orb

  Precision fCubicVolume, fSurfaceArea; ///< Variable to store the cached value of Volume and SurfaceArea

  Precision fEpsilon, fRTolerance; ///< Radial Tolerance

public:
  using Kernel = OrbImplementation;

  /// Default constructor for the unplaced orb.
  VECCORE_ATT_HOST_DEVICE
  UnplacedOrb();

  /// Constructor for the unplaced Orb.
  /** The constructor takes 1 parameter: radius of Orb
      @param r  Radius of Orb.
  */
  VECCORE_ATT_HOST_DEVICE
  UnplacedOrb(const Precision r);

  /// Setter for Radial tolerance
  VECCORE_ATT_HOST_DEVICE
  void SetRadialTolerance();

  /// Getter for Radial tolerance
  VECCORE_ATT_HOST_DEVICE
  Precision GetRadialTolerance() const { return fRTolerance; }

  /// Getter for the structure storing Orb data.
  VECCORE_ATT_HOST_DEVICE
  OrbStruct<Precision> const &GetStruct() const { return fOrb; }

  /// Getter for Radius
  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  Precision GetRadius() const { return fOrb.fR; }

  /// Setter for Radius
  VECCORE_ATT_HOST_DEVICE
  void SetRadius(Precision r);

  VECCORE_ATT_HOST_DEVICE
  void Extent(Vector3D<Precision> &aMin, Vector3D<Precision> &aMax) const override;

  Precision Capacity() const override { return fCubicVolume; }

  Precision SurfaceArea() const override { return fSurfaceArea; }

  virtual Vector3D<Precision> SamplePointOnSurface() const override;

  VECCORE_ATT_HOST_DEVICE
  virtual bool Normal(Vector3D<Precision> const &p, Vector3D<Precision> &normal) const override
  {
    bool valid;
    normal = OrbImplementation::NormalKernel(fOrb, p, valid);
    return valid;
  }

  /// Get the solid type as string.
  /** @return Name of the solid type.*/
  std::string GetEntityType() const;

  /// Get list of Orb parameters as an array.
  /** @param[in]  aNumber Not used.
      @param[out] aArray User array to be filled (rMin, stIn, rMax, stOut, dz)
  */
  VECCORE_ATT_HOST_DEVICE
  void GetParametersList(int aNumber, Precision *aArray) const;

  VECCORE_ATT_HOST_DEVICE
  UnplacedOrb *Clone() const;

  std::ostream &StreamInfo(std::ostream &os) const;

public:
  virtual int MemorySize() const final { return sizeof(*this); }

  VECCORE_ATT_HOST_DEVICE
  virtual void Print() const override;

  virtual void Print(std::ostream &os) const override;

#ifndef VECCORE_CUDA
  virtual SolidMesh *CreateMesh3D(Transformation3D const &trans, size_t nSegments) const override;
#endif

#ifdef VECGEOM_CUDA_INTERFACE
  virtual size_t DeviceSizeOf() const override { return DevicePtr<cuda::UnplacedOrb>::SizeOf(); }
  virtual DevicePtr<cuda::VUnplacedVolume> CopyToGpu() const override;
  virtual DevicePtr<cuda::VUnplacedVolume> CopyToGpu(DevicePtr<cuda::VUnplacedVolume> const gpu_ptr) const override;
#endif

  VECCORE_ATT_DEVICE
  virtual VPlacedVolume *PlaceVolume(LogicalVolume const *const volume, Transformation3D const *const transformation,
#ifdef VECCORE_CUDA
                                     const int id, const int copy_no, const int child_id,
#endif
                                     VPlacedVolume *const placement = NULL) const override;

/// Comparison specific conversion functions
#ifndef VECCORE_CUDA
#ifdef VECGEOM_ROOT
  TGeoShape const *ConvertToRoot(char const *label = "") const;
#endif

#ifdef VECGEOM_GEANT4
  G4VSolid const *ConvertToGeant4(char const *label = "") const;
#endif
#endif // VECCORE_CUDA
};
} // namespace VECGEOM_IMPL_NAMESPACE
} // namespace vecgeom

#endif
