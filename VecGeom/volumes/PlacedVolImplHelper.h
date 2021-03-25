#pragma once

#include "VecGeom/base/Global.h"
#include "VecGeom/volumes/PlacedVolume.h"

#include <algorithm>

#ifdef VECGEOM_DISTANCE_DEBUG
#include "VecGeom/volumes/utilities/ResultComparator.h"
#endif

namespace vecgeom {

VECGEOM_DEVICE_DECLARE_CONV_TEMPLATE_2v(struct, PlacedVolumeImplHelper, class, Arg1, class, Arg2);

inline namespace VECGEOM_IMPL_NAMESPACE {

// A helper template class to automatically implement (wire)
// interfaces from PlacedVolume using kernel functions and functionality
// from the unplaced shapes
template <class UnplacedShape_t, class BaseVol = VPlacedVolume>
struct PlacedVolumeImplHelper : public BaseVol {

  using Struct_t = typename UnplacedShape_t::UnplacedStruct_t;
  using Specialization = typename UnplacedShape_t::Implementation_t;

public:
  using BaseVol::BaseVol;

#ifndef VECCORE_CUDA
  PlacedVolumeImplHelper(VPlacedVolume const *other)
      : BaseVol(other->GetName(), other->GetLogicalVolume(), other->GetTransformation())
  {
  }

  PlacedVolumeImplHelper(char const *const label, LogicalVolume const *const logical_volume,
                                 Transformation3D const *const transformation)
      : BaseVol(label, logical_volume, transformation)
  {
  }

  PlacedVolumeImplHelper(char const *const label, LogicalVolume *const logical_volume,
                                 Transformation3D const *const transformation)
      : BaseVol(label, logical_volume, transformation)
  {
  }

  PlacedVolumeImplHelper(LogicalVolume const *const logical_volume,
                                 Transformation3D const *const transformation)
      : PlacedVolumeImplHelper("", logical_volume, transformation)
  {
  }

  // this constructor mimics the constructor from the Unplaced solid
  // it ensures that placed volumes can be constructed just like ordinary Geant4/ROOT solids
  template <typename... ArgTypes>
  PlacedVolumeImplHelper(char const *const label, ArgTypes... params)
      : PlacedVolumeImplHelper(label, new LogicalVolume(new UnplacedShape_t(params...)),
                                       &Transformation3D::kIdentity)
  {
  }
#else // Compiling for CUDA
  VECCORE_ATT_DEVICE PlacedVolumeImplHelper(LogicalVolume const *const logical_volume,
                                                    Transformation3D const *const transformation, const unsigned int id,
                                                    const int copy_no, const int child_id)
      : BaseVol(logical_volume, transformation, id, copy_no, child_id)
  {
  }
#endif

  using BaseVol::GetLogicalVolume;

  // destructor
  VECCORE_ATT_HOST_DEVICE
  virtual ~PlacedVolumeImplHelper() {}

  virtual int MemorySize() const override { return sizeof(*this); }

  VECCORE_ATT_HOST_DEVICE
  UnplacedShape_t const *GetUnplacedVolume() const
  {
    return static_cast<UnplacedShape_t const *>(GetLogicalVolume()->GetUnplacedVolume());
  }

  VECCORE_ATT_HOST_DEVICE
  Struct_t const *GetUnplacedStruct() const { return static_cast<Struct_t const *>(&GetUnplacedVolume()->GetStruct()); }

#if !defined(VECCORE_CUDA)
  virtual Precision Capacity() override { return const_cast<UnplacedShape_t *>(GetUnplacedVolume())->Capacity(); }

  virtual void Extent(Vector3D<Precision> &aMin, Vector3D<Precision> &aMax) const override
  {
    VPlacedVolume::Extent(aMin, aMax);
  }

  VECCORE_ATT_HOST_DEVICE
  virtual bool Normal(Vector3D<Precision> const &point, Vector3D<Precision> &normal) const override
  {
    return VPlacedVolume::Normal(point, normal);
  }
#endif

  virtual Precision SurfaceArea() const override { return GetUnplacedVolume()->SurfaceArea(); }

  VECCORE_ATT_HOST_DEVICE
  virtual bool UnplacedContains(Vector3D<Precision> const &point) const override
  {
    return GetUnplacedVolume()->UnplacedShape_t::Contains(point);
  }

  VECCORE_ATT_HOST_DEVICE
  virtual Precision DistanceToOut(Vector3D<Precision> const &point, Vector3D<Precision> const &direction,
                                  const Precision stepMax = kInfLength) const override
  {
    return GetUnplacedVolume()->UnplacedShape_t::DistanceToOut(point, direction, stepMax);
  }

  VECCORE_ATT_HOST_DEVICE
  virtual Precision SafetyToOut(Vector3D<Precision> const &point) const override
  {
    return GetUnplacedVolume()->UnplacedShape_t::SafetyToOut(point);
  }

  // Calls to the implementation.
  VECCORE_ATT_HOST_DEVICE
  virtual EnumInside Inside(Vector3D<Precision> const &point) const override
  {
    Inside_t output;
    Transformation3D const *tr = this->GetTransformation();
    Specialization::Inside(*this->GetUnplacedStruct(), tr->Transform<Precision>(point), output);
    return (EnumInside)output;
  }

  VECCORE_ATT_HOST_DEVICE
  virtual bool Contains(Vector3D<Precision> const &point) const override
  {
    bool output(false);
    Transformation3D const *tr = this->GetTransformation();
    Vector3D<Precision> lp     = tr->Transform(point);
    Specialization::Contains(*this->GetUnplacedStruct(), lp, output);
    return output;
  }

  VECCORE_ATT_HOST_DEVICE
  virtual bool Contains(Vector3D<Precision> const &point, Vector3D<Precision> &localPoint) const override
  {
    bool output(false);
    Vector3D<Precision> lp;
    Transformation3D const *tr = this->GetTransformation();
    localPoint                 = tr->Transform(point);
    Specialization::Contains(*this->GetUnplacedStruct(), localPoint, output);
#ifdef VECGEOM_DISTANCE_DEBUG
    DistanceComparator::CompareUnplacedContains(this, output, localPoint);
#endif
    return output;
  }

  VECCORE_ATT_HOST_DEVICE
  virtual Precision DistanceToIn(Vector3D<Precision> const &point, Vector3D<Precision> const &direction,
                                 const Precision stepMax = kInfLength) const override
  {
#ifndef VECCORE_CUDA
    assert(direction.IsNormalized() && " direction not normalized in call to DistanceToIn ");
#endif
    Precision output(kInfLength);
    Transformation3D const *tr = this->GetTransformation();
    Specialization::DistanceToIn(*this->GetUnplacedStruct(), tr->Transform(point),
                                 tr->TransformDirection(direction), stepMax, output);
#ifdef VECGEOM_DISTANCE_DEBUG
    DistanceComparator::CompareDistanceToIn(this, output, point, direction, stepMax);
#endif
    return output;
  }

  VECCORE_ATT_HOST_DEVICE
  virtual Precision PlacedDistanceToOut(Vector3D<Precision> const &point, Vector3D<Precision> const &direction,
                                        const Precision stepMax = kInfLength) const override
  {
#ifndef VECCORE_CUDA
    assert(direction.IsNormalized() && " direction not normalized in call to PlacedDistanceToOut ");
#endif
    Transformation3D const *tr = this->GetTransformation();
    Precision output(-1.);
    Specialization::template DistanceToOut(*this->GetUnplacedStruct(), tr->Transform(point),
                                           tr->TransformDirection(direction), stepMax, output);

#ifdef VECGEOM_DISTANCE_DEBUG
    DistanceComparator::CompareDistanceToOut(this, output, this->GetTransformation()->Transform(point),
                                             this->GetTransformation()->TransformDirection(direction), stepMax);
#endif
    return output;
  }

  VECCORE_ATT_HOST_DEVICE
  virtual Precision SafetyToIn(Vector3D<Precision> const &point) const override
  {
    Precision output(kInfLength);
    Transformation3D const *tr = this->GetTransformation();
    Specialization::SafetyToIn(*this->GetUnplacedStruct(), tr->Transform(point), output);
    return output;
  }

#ifdef VECGEOM_CUDA_INTERFACE
  using ThisClass_t = typename Specialization::PlacedShape_t;
  virtual size_t DeviceSizeOf() const override { return DevicePtr<CudaType_t<ThisClass_t>>::SizeOf(); }

  DevicePtr<cuda::VPlacedVolume> CopyToGpu(DevicePtr<cuda::LogicalVolume> const logical_volume,
                                           DevicePtr<cuda::Transformation3D> const transform,
                                           DevicePtr<cuda::VPlacedVolume> const in_gpu_ptr) const override
  {
    DevicePtr<CudaType_t<ThisClass_t>> gpu_ptr(in_gpu_ptr);
    gpu_ptr.Construct(logical_volume, transform, this->id(), this->GetCopyNo(), this->GetChildId());
    CudaAssertError();
    // Need to go via the void* because the regular c++ compilation
    // does not actually see the declaration for the cuda version
    // (and thus can not determine the inheritance).
    return DevicePtr<cuda::VPlacedVolume>((void *)gpu_ptr);
  }

  DevicePtr<cuda::VPlacedVolume> CopyToGpu(DevicePtr<cuda::LogicalVolume> const logical_volume,
                                           DevicePtr<cuda::Transformation3D> const transform) const override
  {
    DevicePtr<CudaType_t<ThisClass_t>> gpu_ptr;
    gpu_ptr.Allocate();
    return CopyToGpu(logical_volume, transform, DevicePtr<cuda::VPlacedVolume>((void *)gpu_ptr));
  }
#endif // VECGEOM_CUDA_INTERFACE


}; // End class PlacedVolumeImplHelper
} // namespace VECGEOM_IMPL_NAMESPACE
} // namespace vecgeom
