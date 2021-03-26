/// \file CudaManager.cu
/// \author Johannes de Fine Licht (johannes.definelicht@cern.ch)

#include "VecGeom/management/CudaManager.h"

#include <stdio.h>

#include "VecGeom/backend/cuda/Backend.h"

namespace vecgeom {

// we put some global data into a separate namespace
// this is done since CUDA does not support static const members in class definitions
namespace globaldevicegeomdata {
static __device__ VPlacedVolume *gCompactPlacedVolBuffer = nullptr;

static __device__ NavIndex_t *gNavIndex = nullptr;

static __device__ int gMaxDepth = 0;

VECCORE_ATT_DEVICE
VPlacedVolume *&GetCompactPlacedVolBuffer()
{
  return gCompactPlacedVolBuffer;
}

VECCORE_ATT_DEVICE
NavIndex_t *&GetNavIndex()
{
  return gNavIndex;
}

VECCORE_ATT_DEVICE
int GetMaxDepth()
{
  return gMaxDepth;
}
}

inline namespace cuda {

__global__ void InitDeviceCompactPlacedVolBufferPtrCudaKernel(void *gpu_ptr)
{
  // gpu_ptr is some pointer on the device that was allocated by some other means
  globaldevicegeomdata::GetCompactPlacedVolBuffer() = (vecgeom::cuda::VPlacedVolume *)gpu_ptr;
}

void InitDeviceCompactPlacedVolBufferPtr(void *gpu_ptr)
{
  InitDeviceCompactPlacedVolBufferPtrCudaKernel<<<1, 1>>>(gpu_ptr);
}

__global__ void InitDeviceNavIndexPtrCudaKernel(void *gpu_ptr, int maxdepth)
{
  // gpu_ptr is some pointer on the device that was allocated by some other means
  globaldevicegeomdata::GetNavIndex() = (NavIndex_t *)gpu_ptr;
  globaldevicegeomdata::gMaxDepth = maxdepth;
}

void InitDeviceNavIndexPtr(void *gpu_ptr, int maxdepth)
{
  InitDeviceNavIndexPtrCudaKernel<<<1, 1>>>(gpu_ptr, maxdepth);
}

__global__ void CudaManagerPrintGeometryKernel(vecgeom::cuda::VPlacedVolume const *const world)
{
  printf("Geometry loaded on GPU:\n");
  world->PrintContent();
}

void CudaManagerPrintGeometry(vecgeom::cuda::VPlacedVolume const *const world)
{
  CudaManagerPrintGeometryKernel<<<1, 1>>>(world);
  cxx::CudaAssertError();
  cudaDeviceSynchronize();
}
}
} // End namespace vecgeom
