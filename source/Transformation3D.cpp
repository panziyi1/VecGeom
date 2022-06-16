/// \file Transformation3D.cpp
/// \author Johannes de Fine Licht (johannes.definelicht@cern.ch)
#include "VecGeom/base/Transformation3D.h"

#ifdef VECGEOM_CUDA_INTERFACE
#include "VecGeom/backend/cuda/Interface.h"
#endif
#include "VecGeom/base/SpecializedTransformation3D.h"

#ifdef VECGEOM_ROOT
#include "TGeoMatrix.h"
#endif

#include <sstream>
#include <stdio.h>

namespace vecgeom {
inline namespace VECGEOM_IMPL_NAMESPACE {

#ifdef VECGEOM_ROOT
// function to convert this transformation to a TGeo transformation
// mainly used for the benchmark comparisons with ROOT
TGeoMatrix *Transformation3D::ConvertToTGeoMatrix() const
{
  double rotd[9];
  if (fHasRotation) {
    for (auto i = 0; i < 9; ++i)
      rotd[i] = Rotation()[i];
  }

  if (fIdentity) {
    return new TGeoIdentity();
  }
  if (fHasTranslation && !fHasRotation) {
    return new TGeoTranslation(fTranslation[0], fTranslation[1], fTranslation[2]);
  }
  if (fHasRotation && !fHasTranslation) {
    TGeoRotation *tmp = new TGeoRotation();
    tmp->SetMatrix(rotd);
    return tmp;
  }
  if (fHasTranslation && fHasRotation) {
    TGeoRotation *tmp = new TGeoRotation();
    tmp->SetMatrix(rotd);
    return new TGeoCombiTrans(fTranslation[0], fTranslation[1], fTranslation[2], tmp);
  }
  return 0;
}
#endif

std::ostream &operator<<(std::ostream &os, Transformation3D const &transformation)
{
  os << "Transformation {" << transformation.Translation() << ", "
     << "(" << transformation.Rotation(0) << ", " << transformation.Rotation(1) << ", " << transformation.Rotation(2)
     << ", " << transformation.Rotation(3) << ", " << transformation.Rotation(4) << ", " << transformation.Rotation(5)
     << ", " << transformation.Rotation(6) << ", " << transformation.Rotation(7) << ", " << transformation.Rotation(8)
     << ")}"
     << "; identity(" << transformation.IsIdentity() << "); rotation(" << transformation.HasRotation() << ")";
  return os;
}

#ifdef VECGEOM_CUDA_INTERFACE

DevicePtr<cuda::Transformation3D> Transformation3D::CopyToGpu(DevicePtr<cuda::Transformation3D> const gpu_ptr) const
{

  gpu_ptr.Construct(fTranslation[0], fTranslation[1], fTranslation[2], fRotation[0], fRotation[1], fRotation[2],
                    fRotation[3], fRotation[4], fRotation[5], fRotation[6], fRotation[7], fRotation[8]);
  CudaAssertError();
  return gpu_ptr;
}

DevicePtr<cuda::Transformation3D> Transformation3D::CopyToGpu() const
{

  DevicePtr<cuda::Transformation3D> gpu_ptr;
  gpu_ptr.Allocate();
  return this->CopyToGpu(gpu_ptr);
}

/**
 * Copy a large number of transformation instances to the GPU.
 * \param trafos Host instances to copy.
 * \param gpu_ptrs Device pointers to indicate where the transformations should be placed.
 * The device memory must have been allocated before copying.
 */
void Transformation3D::CopyManyToGpu(const std::vector<Transformation3D const *>& trafos,
                                     const std::vector<DevicePtr<cuda::Transformation3D>>& gpu_ptrs)
{
  assert(trafos.size() == gpu_ptrs.size());

  // Memory for constructor data
  // Store it as
  // tx0, tx1, tx2, ...,
  // ty0, ty1, ty2, ...,
  // ...
  // rot0_0, rot0_1, rot0_2, ...,
  // ...
  std::vector<Precision> trafoData(12 * trafos.size());

  std::size_t trafoCounter = 0;
  for (Transformation3D const * trafo : trafos) {
    for (unsigned int i = 0; i < 3; ++i) trafoData[trafoCounter +  i    * trafos.size()] = trafo->Translation(i);
    for (unsigned int i = 0; i < 9; ++i) trafoData[trafoCounter + (i+3) * trafos.size()] = trafo->Rotation(i);
    ++trafoCounter;
  }

  ConstructManyOnGpu<cuda::Transformation3D>(trafos.size(), gpu_ptrs.data(),
      trafoData.data(),                     trafoData.data() +  1 * trafos.size(), trafoData.data() +  2 * trafos.size(), // translations
      trafoData.data() + 3 * trafos.size(), trafoData.data() +  4 * trafos.size(), trafoData.data() +  5 * trafos.size(), // rotations
      trafoData.data() + 6 * trafos.size(), trafoData.data() +  7 * trafos.size(), trafoData.data() +  8 * trafos.size(),
      trafoData.data() + 9 * trafos.size(), trafoData.data() + 10 * trafos.size(), trafoData.data() + 11 * trafos.size()
  );
}

#endif // VECGEOM_CUDA_INTERFACE

} // namespace VECGEOM_IMPL_NAMESPACE

#ifdef VECCORE_CUDA

namespace cxx {

template size_t DevicePtr<cuda::Transformation3D>::SizeOf();
template void DevicePtr<cuda::Transformation3D>::Construct(const Precision tx, const Precision ty, const Precision tz,
                                                           const Precision r0, const Precision r1, const Precision r2,
                                                           const Precision r3, const Precision r4, const Precision r5,
                                                           const Precision r6, const Precision r7,
                                                           const Precision r8) const;
template void ConstructManyOnGpu<Transformation3D>(std::size_t, DevicePtr<cuda::Transformation3D> const *,
                                                   Precision const * tx, Precision const * ty, Precision const * tz,
                                                   Precision const * r0, Precision const * r1, Precision const * r2,
                                                   Precision const * r3, Precision const * r4, Precision const * r5,
                                                   Precision const * r6, Precision const * r7, Precision const * r8);

} // namespace cxx

#endif // VECCORE_CUDA

} // namespace vecgeom
