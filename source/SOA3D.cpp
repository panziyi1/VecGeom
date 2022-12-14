/// \file SOA3D.cpp
/// \author Johannes de Fine Licht (johannes.definelicht@cern.ch)

#include "VecGeom/base/SOA3D.h"
#ifdef VECGEOM_ENABLE_CUDA
#include "VecGeom/backend/cuda/Interface.h"
#endif

namespace vecgeom {

#ifdef VECCORE_CUDA

namespace cxx {

template size_t DevicePtr<cuda::SOA3D<Precision>>::SizeOf();
template void DevicePtr<cuda::SOA3D<Precision>>::Construct(DevicePtr<Precision> x, DevicePtr<Precision> y,
                                                           DevicePtr<Precision> z, size_t size) const;

} // End cxx namespace

#endif

} // End namespace vecgeom
