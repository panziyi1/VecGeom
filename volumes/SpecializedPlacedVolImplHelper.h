#pragma once

#include "base/Global.h"
#include "base/SOA3D.h"
#include "volumes/PlacedBox.h"

#include <algorithm>

#ifdef VECGEOM_DISTANCE_DEBUG
#include "volumes/utilities/ResultComparator.h"
#endif

namespace vecgeom {

// putting a forward declaration by hand
VECGEOM_DEVICE_DECLARE_CONV_TEMPLATE_1t_2v(class, CommonSpecializedVolImplHelper,typename,TranslationCode, translation::kGeneric, RotationCode, rotation::kGeneric)
VECGEOM_DEVICE_DECLARE_CONV_TEMPLATE_1t_2v(class, SIMDSpecializedVolImplHelper,class,TranslationCode, translation::kGeneric, RotationCode, rotation::kGeneric)

inline namespace VECGEOM_IMPL_NAMESPACE {

template <class Specialization, TranslationCode transC, RotationCode rotC>
class CommonSpecializedVolImplHelper : public Specialization::PlacedShape_t {

  using PlacedShape_t = typename Specialization::PlacedShape_t;
  using UnplacedVolume_t = typename Specialization::UnplacedVolume_t;

public:

#ifndef VECGEOM_NVCC
  CommonSpecializedVolImplHelper(char const *const label, LogicalVolume const *const logical_volume,
                                 Transformation3D const *const transformation,
                                 vecgeom::PlacedBox const *const boundingBox)
      : PlacedShape_t(label, logical_volume, transformation, boundingBox) {}

  CommonSpecializedVolImplHelper(char const *const label, LogicalVolume const *const logical_volume,
                                 Transformation3D const *const transformation)
      : CommonSpecializedVolImplHelper(label, logical_volume, transformation,
                                       details::UseIfSameType<PlacedShape_t, vecgeom::PlacedBox>::Get(this)) {}

  CommonSpecializedVolImplHelper(char const *const label, LogicalVolume *const logical_volume,
                                 Transformation3D const *const transformation,
                                 vecgeom::PlacedBox const *const boundingBox)
      : PlacedShape_t(label, logical_volume, transformation, boundingBox) {}

  CommonSpecializedVolImplHelper(char const *const label, LogicalVolume *const logical_volume,
                                 Transformation3D const *const transformation)
      : CommonSpecializedVolImplHelper(label, logical_volume, transformation,
                                       details::UseIfSameType<PlacedShape_t, vecgeom::PlacedBox>::Get(this)) {}

  CommonSpecializedVolImplHelper(LogicalVolume const *const logical_volume,
                                 Transformation3D const *const transformation,
                                 vecgeom::PlacedBox const *const boundingBox)
      : CommonSpecializedVolImplHelper("", logical_volume, transformation, boundingBox) {}

  CommonSpecializedVolImplHelper(LogicalVolume const *const logical_volume,
                                 Transformation3D const *const transformation)
      : CommonSpecializedVolImplHelper("", logical_volume, transformation) {}

  // this constructor mimics the constructor from the Unplaced solid
  // it ensures that placed volumes can be constructed just like ordinary Geant4/ROOT/USolids solids
  template <typename... ArgTypes>
  CommonSpecializedVolImplHelper(char const *const label, ArgTypes... params)
      : CommonSpecializedVolImplHelper(label, new LogicalVolume(new UnplacedVolume_t(params...)),
                                       &Transformation3D::kIdentity) {}

#else // Compiling for CUDA
  __device__ CommonSpecializedVolImplHelper(LogicalVolume const *const logical_volume,
                                            Transformation3D const *const transformation,
                                            PlacedBox const *const boundingBox, const int id)
      : PlacedShape_t(logical_volume, transformation, boundingBox, id) {}

  __device__ CommonSpecializedVolImplHelper(LogicalVolume const *const logical_volume,
                                            Transformation3D const *const transformation, const int id)
      : PlacedShape_t(logical_volume, transformation, details::UseIfSameType<PlacedShape_t, PlacedBox>::Get(this), id) {
  }
#endif
  using PlacedShape_t::SafetyToOut;
  using PlacedShape_t::DistanceToOut;
  using PlacedShape_t::UnplacedContains;
  using PlacedShape_t::Contains;
  using PlacedShape_t::SafetyToIn;
  using PlacedShape_t::DistanceToIn;
  using PlacedShape_t::Inside;
  using PlacedShape_t::PlacedShape_t;

  virtual int memory_size() const override { return sizeof(*this); }

  VECGEOM_CUDA_HEADER_BOTH
  virtual void PrintType() const override { Specialization::PrintType(); }

  virtual void PrintType(std::ostream &os) const override { Specialization::PrintType(os); }
  virtual void PrintImplementationType(std::ostream &os) const override { Specialization::PrintImplementationType(os); }
  virtual void PrintUnplacedType(std::ostream &os) const override { Specialization::PrintUnplacedType(os); }

  VECGEOM_CUDA_HEADER_BOTH
  virtual EnumInside Inside(Vector3D<Precision> const &point) const override {
    Inside_t output;
    Transformation3D const *tr = this->GetTransformation();
    Specialization::Inside(*this->GetUnplacedStruct(),
                          tr->Transform<transC, rotC, Precision>(point), output);
    return (EnumInside)output;
  }

  VECGEOM_CUDA_HEADER_BOTH
  virtual bool Contains(Vector3D<Precision> const &point) const override {
    bool output(false);
    Transformation3D const *tr = this->GetTransformation();
    Vector3D<Precision> lp = tr->Transform<transC, rotC, Precision>(point);
    Specialization::Contains(*this->GetUnplacedStruct(), lp, output);
    return output;
  }

  VECGEOM_CUDA_HEADER_BOTH
  virtual bool Contains(Vector3D<Precision> const &point,
                        Vector3D<Precision> &localPoint) const override {
    bool output(false);
    Vector3D<double> lp;
    Transformation3D const *tr = this->GetTransformation();
    localPoint = tr->Transform<transC, rotC,Precision>(point);
    Specialization::Contains(
      *this->GetUnplacedStruct(),
      localPoint,
      output
      );
#ifdef VECGEOM_DISTANCE_DEBUG
    DistanceComparator::CompareUnplacedContains( this, output, localPoint );
#endif
return output;
  }

  VECGEOM_CUDA_HEADER_BOTH
  virtual Precision DistanceToIn(Vector3D<Precision> const &point,
                                 Vector3D<Precision> const &direction,
                                 const Precision stepMax = kInfinity) const override {
#ifndef VECGEOM_NVCC
      assert( direction.IsNormalized() && " direction not normalized in call to DistanceToIn " );
#endif
      Precision output(kInfinity);
      Transformation3D const *tr = this->GetTransformation();
      Specialization::DistanceToIn(
              *this->GetUnplacedStruct(),
              tr->Transform<transC,rotC>(point),
      tr->TransformDirection<rotC>(direction),
      stepMax,
      output
      ); 
#ifdef VECGEOM_DISTANCE_DEBUG
    DistanceComparator::CompareDistanceToIn( this, output, point, direction, stepMax );
#endif
    return output;
  }

  VECGEOM_CUDA_HEADER_BOTH
  virtual Precision PlacedDistanceToOut(Vector3D<Precision> const &point,
                                        Vector3D<Precision> const &direction,
                                        const Precision stepMax = kInfinity) const override {
#ifndef VECGEOM_NVCC
      assert( direction.IsNormalized() && " direction not normalized in call to  PlacedDistanceToOut " );
#endif
      (void)point;
      (void)direction;
      (void)stepMax;
    // Transformation3D const *tr = this->GetTransformation();
     /*  Real_v output = Specialization<VecCore::VecGeomBackend::Scalar<double>>::DistanceToOut(
        *this->GetUnplacedStruct(),
        tr->Transform<transC, rotC>(point),
        tr->TransformDirection<rotC>(direction),
        stepMax);

  #ifdef VECGEOM_DISTANCE_DEBUG
      DistanceComparator::CompareDistanceToOut(
              this,
              output,
              this->GetTransformation()->Transform(point),
              this->GetTransformation()->TransformDirection(direction),
              stepMax );
  #endif
  return output; */
      return 0.;
    }

  VECGEOM_CUDA_HEADER_BOTH
  virtual Precision SafetyToIn(Vector3D<Precision> const &point) const override {
    Precision output(kInfinity);
    Transformation3D const *tr = this->GetTransformation();
    Specialization::SafetyToIn(
      *this->GetUnplacedStruct(),
      tr->Transform<transC,rotC>(point),
      output
    );
    return output;
  }



  virtual VECGEOM_BACKEND_PRECISION_TYPE
  DistanceToInVec(Vector3D<VECGEOM_BACKEND_PRECISION_TYPE> const &position,
               Vector3D<VECGEOM_BACKEND_PRECISION_TYPE> const &direction,
               const VECGEOM_BACKEND_PRECISION_TYPE step_max = kInfinity) const override {
    (void)position;
    (void)direction;
    (void)step_max;
    return -1.;
  }
  virtual VECGEOM_BACKEND_PRECISION_TYPE
  DistanceToOutVec(Vector3D<VECGEOM_BACKEND_PRECISION_TYPE> const &position,
                Vector3D<VECGEOM_BACKEND_PRECISION_TYPE> const &direction,
                VECGEOM_BACKEND_PRECISION_TYPE const step_max = kInfinity) const override {
    (void)position;
    (void)direction;
    (void)step_max;
    return -1.;
  }
  virtual VECGEOM_BACKEND_PRECISION_TYPE
  SafetyToInVec(Vector3D<VECGEOM_BACKEND_PRECISION_TYPE> const &position) const override {
    (void)position;
    return -1.;
  }
  virtual VECGEOM_BACKEND_PRECISION_TYPE
  SafetyToOutVec(Vector3D<VECGEOM_BACKEND_PRECISION_TYPE> const &position) const override {
    (void)position;
    return -1.;
  }

}; // End class CommonSpecializedVolImplHelper

// needs to be in the specializations
template <class Specialization, typename Real_v, int transC, int rotC>
static void ContainsLoopKernel(typename Specialization::UnplacedStruct_t const &shapestruct,
                               Transformation3D const &trans, const size_t offset, const size_t size,
                               SOA3D<Precision> const &points, bool *const output) {

  using Bool_v = typename vecCore::Mask_v<Real_v>;
  for (decltype(points.size()) i(offset); i < size; i += vecCore::VectorSize<Real_v>()) {
    Vector3D<Real_v> point(vecCore::FromPtr<Real_v>(points.x() + i), vecCore::FromPtr<Real_v>(points.y() + i),
                           vecCore::FromPtr<Real_v>(points.z() + i));
    Bool_v result(false);
    Specialization::template Contains<Real_v>(shapestruct, trans.Transform<transC, rotC>(point), result);
    // with bools we can't write back vectorized as masks have a different size than a pod bool
    for (unsigned j = 0; j < vecCore::VectorSize<Real_v>(); ++j) {
      output[j + i] = vecCore::MaskLaneAt(result, j);
    }
  }
}

template <class Specialization, typename Real_v, int transC, int rotC>
static void InsideLoopKernel(typename Specialization::UnplacedStruct_t const &shapestruct,
                               Transformation3D const &trans, const size_t offset, const size_t size,
                               SOA3D<Precision> const &points, Inside_t *const output) {
 using Index_t = vecCore::Index_v<Real_v>;
 for (decltype(points.size()) i(offset); i < size; i += vecCore::VectorSize<Real_v>()) {
    Vector3D<Real_v> point(vecCore::FromPtr<Real_v>(points.x() + i), vecCore::FromPtr<Real_v>(points.y() + i),
                           vecCore::FromPtr<Real_v>(points.z() + i));
    Index_t result;
    Specialization::template Inside<Real_v>(shapestruct, trans.Transform<transC, rotC>(point), result);
    // TODO: make a proper store here
    for(size_t j = 0; j < vecCore::VectorSize<Index_t>(); ++j)
      output[i+j] = vecCore::LaneAt<Index_t>(result, j);
  }
}

template <class Specialization, typename Real_v, int transC, int rotC>
static void SafetyToInLoopKernel(
    typename Specialization::UnplacedStruct_t const &shapestruct,
    Transformation3D const &trans, const size_t offset, const size_t size, SOA3D<Precision> const &points,
    double *const output) {

  for (decltype(points.size()) i(offset); i < size; i += vecCore::VectorSize<Real_v>()) {
    Vector3D<Real_v> point(vecCore::FromPtr<Real_v>(points.x() + i), vecCore::FromPtr<Real_v>(points.y() + i),
                           vecCore::FromPtr<Real_v>(points.z() + i));
    Real_v result(kInfinity);
    Specialization::template SafetyToIn<Real_v>(shapestruct, trans.Transform<transC, rotC>(point), result);
    vecCore::Store(result, output + i);
  }
}


 template <class Specialization, typename Real_v, int transC, int rotC>
 static void DistanceToInLoopKernel(typename Specialization::UnplacedStruct_t const &shapestruct,
                                    Transformation3D const &trans, const size_t offset, const size_t size,
                                    SOA3D<Precision> const &points, SOA3D<Precision> const &directions,
                                    Precision const *const stepMax, double *const output) {

   for (decltype(points.size()) i(offset); i < size; i += vecCore::VectorSize<Real_v>()) {
     Vector3D<Real_v> point(vecCore::FromPtr<Real_v>(points.x() + i), vecCore::FromPtr<Real_v>(points.y() + i),
                            vecCore::FromPtr<Real_v>(points.z() + i));
     Vector3D<Real_v> dir(vecCore::FromPtr<Real_v>(directions.x() + i), vecCore::FromPtr<Real_v>(directions.y() + i),
                          vecCore::FromPtr<Real_v>(directions.z() + i));
     Real_v step_max(vecCore::FromPtr<Real_v>(stepMax + i));
     Real_v result(kInfinity);
     Specialization::template DistanceToIn<Real_v>(shapestruct, trans.Transform<transC, rotC>(point),
                                           trans.TransformDirection<rotC>(dir), step_max, result);
     vecCore::Store(result, output + i);
   }
 }

  template <class Specialization, int transC, int rotC>
  class SIMDSpecializedVolImplHelper : public CommonSpecializedVolImplHelper<Specialization, transC, rotC> {
    using CommonHelper_t = CommonSpecializedVolImplHelper<Specialization, transC, rotC>;

  public:
    using CommonHelper_t::SafetyToOut;
    using CommonHelper_t::DistanceToOut;
    using CommonHelper_t::UnplacedContains;
    using CommonHelper_t::Contains;
    using CommonHelper_t::SafetyToIn;
    using CommonHelper_t::DistanceToIn;
    using CommonHelper_t::Inside;
    using CommonHelper_t::CommonHelper_t;

    VECGEOM_CUDA_HEADER_BOTH
    virtual ~SIMDSpecializedVolImplHelper() {}

    virtual void SafetyToIn(SOA3D<Precision> const &points, Precision *const output) const override {
      const auto kS = vecCore::VectorSize<VectorBackend::Real_v>();
      auto offset = points.size() - points.size() % kS;
      //   auto shape = ((UnplacedVolume_t *)this)->UnplacedVolume_t::GetUnplacedStruct();
      auto shape = this->GetUnplacedStruct();
      auto transf = this->GetTransformation();

      // vector loop treatment
      SafetyToInLoopKernel<Specialization, VectorBackend::Real_v, transC, rotC>(*shape, *transf, 0, offset, points,
                                                                                   output);
      // tail treatment
      SafetyToInLoopKernel<Specialization, ScalarBackend::Real_v, transC, rotC>(*shape, *transf, offset,
                                                                                   points.size(), points, output);
    }

    virtual void SafetyToInMinimize(SOA3D<Precision> const &/* points */, Precision *const /* safeties */) const override {
      // we do no longer need this (probably)
      // SafetyToInMinimizeTemplate(points, safeties);
    }

    virtual void DistanceToIn(SOA3D<Precision> const &points, SOA3D<Precision> const &directions,
                              Precision const *const stepMax, Precision *const output) const override {
      auto offset = points.size() - points.size() % vecCore::VectorSize<VectorBackend::Real_v>();
      auto shape = this->GetUnplacedStruct();
      auto transf = this->GetTransformation();
      // vector loop treatment
      DistanceToInLoopKernel<Specialization, VectorBackend::Real_v, transC, rotC>(*shape, *transf, 0, offset, points, directions,
                                                                    stepMax, output);
      // tail treatment
      DistanceToInLoopKernel<Specialization, ScalarBackend::Real_v, transC, rotC>(
          *shape, *transf, offset, points.size(), points, directions, stepMax, output);
    }

    virtual void DistanceToInMinimize(SOA3D<Precision> const &/* points */, SOA3D<Precision> const &/* directions */,
                                      int /* daughterindex */, Precision *const /* output */,
                                      int *const /* nextnodeids */) const override {
      // we do no longer need this (probably)
      // DistanceToInMinimizeTemplate(points, directions, daughterindex, output, nextnodeids);
    }

    virtual void Contains(SOA3D<Precision> const &points, bool *const output) const override {
      auto offset = points.size() - points.size() % vecCore::VectorSize<VectorBackend::Real_v>();
      auto shape = this->GetUnplacedStruct();
      auto transf = this->GetTransformation();
      // vector loop treatment
      ContainsLoopKernel<Specialization, VectorBackend::Real_v, transC, rotC>(*shape, *transf, 0, offset, points,
                                                                                   output);
      // tail treatment
      ContainsLoopKernel<Specialization, ScalarBackend::Real_v, transC, rotC>(*shape, *transf, offset, points.size(),
                                                                                 points, output);
    }

    virtual void Inside(SOA3D<Precision> const &points, Inside_t *const output) const override {
      // I would be in favor of getting rid of this interface (unless someone asks for it)
      // Inside is only provided for Geant4 which currently does not have a basket interface
      // InsideTemplate(points, output);
      auto offset = points.size() - points.size() % vecCore::VectorSize<VectorBackend::Real_v>();
      auto shape = this->GetUnplacedStruct();
      auto transf = this->GetTransformation();
      // vector loop treatment
      InsideLoopKernel<Specialization, VectorBackend::Real_v, transC, rotC>(*shape, *transf, 0, offset, points,
                                                                                       output);
      // tail treatment
      InsideLoopKernel<Specialization, ScalarBackend::Real_v, transC, rotC>(*shape, *transf, offset, points.size(),
                                                                                     points, output);
    }

#ifdef VECGEOM_CUDA_INTERFACE
    using ThisClass_t = SIMDSpecializedVolImplHelper<Specialization, transC, rotC>;

    virtual size_t DeviceSizeOf() const override { return DevicePtr<CudaType_t<ThisClass_t>>::SizeOf(); }

    DevicePtr<cuda::VPlacedVolume> CopyToGpu(DevicePtr<cuda::LogicalVolume> const logical_volume,
                                             DevicePtr<cuda::Transformation3D> const transform,
                                             DevicePtr<cuda::VPlacedVolume> const in_gpu_ptr) const override {
      DevicePtr<CudaType_t<ThisClass_t>> gpu_ptr(in_gpu_ptr);
      gpu_ptr.Construct(logical_volume, transform, DevicePtr<cuda::PlacedBox>(), this->id());
      CudaAssertError();
      // Need to go via the void* because the regular c++ compilation
      // does not actually see the declaration for the cuda version
      // (and thus can not determine the inheritance).
      return DevicePtr<cuda::VPlacedVolume>((void *)gpu_ptr);
    }

    DevicePtr<cuda::VPlacedVolume> CopyToGpu(DevicePtr<cuda::LogicalVolume> const logical_volume,
                                             DevicePtr<cuda::Transformation3D> const transform) const override {
    DevicePtr<CudaType_t<ThisClass_t>> gpu_ptr;
    gpu_ptr.Allocate();
    return CopyToGpu(logical_volume, transform, DevicePtr<cuda::VPlacedVolume>((void *)gpu_ptr));
  }
#endif // VECGEOM_CUDA_INTERFACE

  }; // end SIMD Helper

  template <class Specialization, int transC, int rotC>
  class LoopSpecializedVolImplHelper : public CommonSpecializedVolImplHelper<Specialization, transC, rotC> {
    using Helper_t = LoopSpecializedVolImplHelper<Specialization, transC, rotC>;

  public:
    virtual void SafetyToIn(SOA3D<Precision> const &points, Precision *const output) const override {
      (void)points;
      (void)output;
      // SafetyToInTemplate(points, output);
    }

    virtual void SafetyToInMinimize(SOA3D<Precision> const &points, Precision *const safeties) const override {
      (void)points;
      (void)safeties;
      // SafetyToInMinimizeTemplate(points, safeties);
    }

    virtual void Contains(SOA3D<Precision> const &points, bool *const output) const override {
      (void)points;
      (void)output;
      // ContainsTemplate(points, output);
    }

    virtual void Inside(SOA3D<Precision> const &points, Inside_t *const output) const override {
      (void)points;
      (void)output;
      // InsideTemplate(points, output);
    }

    virtual void DistanceToIn(SOA3D<Precision> const &points, SOA3D<Precision> const &directions,
                              Precision const *const stepMax, Precision *const output) const override {
      (void)points;
      (void)directions;
      (void)stepMax;
      (void)output;
      // DistanceToInTemplate(points, directions, stepMax, output);
    }

    virtual void DistanceToInMinimize(SOA3D<Precision> const &points, SOA3D<Precision> const &directions,
                                      int daughterindex, Precision *const output,
                                      int *const nextnodeids) const override {
      (void)points;
      (void)directions;
      (void)daughterindex;
      (void)output;
      (void)nextnodeids;
      // DistanceToInMinimizeTemplate(points, directions, daughterindex, output, nextnodeids);
    }
  };  // end Loop Helper
} } // End global namespace
