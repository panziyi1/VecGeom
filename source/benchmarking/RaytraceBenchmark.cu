/// \file Raytracer.cu
/// \author Guilherme Amadio. Rewritten to use navigation from common code by Andrei Gheata.

#include "VecGeom/benchmarking/Raytracer.h"

#include <VecGeom/base/Transformation3D.h>
#include <VecGeom/management/GeoManager.h>
#include <VecGeom/management/CudaManager.h>
#include <VecGeom/navigation/NavStateIndex.h>
#include <VecGeom/volumes/PlacedVolume.h>
#include <VecGeom/benchmarking/Raytracer.h>
#include <VecGeom/base/Stopwatch.h>

#include <cassert>
#include <cstdio>

using namespace vecgeom;

void check_cuda_err(cudaError_t result, char const *const func, const char *const file, int const line)
{
  if (result) {
    fprintf(stderr, "CUDA error = %s at %s:%d\n", cudaGetErrorString(result), file, line);
    cudaDeviceReset();
    exit(1);
  }
}

#define checkCudaErrors(val) check_cuda_err((val), #val, __FILE__, __LINE__)

__global__ void RenderKernel(RaytracerData_t rtdata, char *input_buffer, unsigned char *output_buffer)
{
  int px = threadIdx.x + blockIdx.x * blockDim.x;
  int py = threadIdx.y + blockIdx.y * blockDim.y;

  if ((px >= rtdata.fSize_px) || (py >= rtdata.fSize_py)) return;

  int ray_index = py * rtdata.fSize_px + px;

  Ray_t *ray          = (Ray_t *)(input_buffer + ray_index * sizeof(Ray_t));
  Color_t pixel_color = Raytracer::RaytraceOne(rtdata, *ray, px, py);

  int pixel_index = 4 * ray_index;

  output_buffer[pixel_index + 0] = pixel_color.fComp.red;
  output_buffer[pixel_index + 1] = pixel_color.fComp.green;
  output_buffer[pixel_index + 2] = pixel_color.fComp.blue;
  output_buffer[pixel_index + 3] = 255;
}

__global__ void RenderTile(RaytracerData_t rtdata, int offset_x, int offset_y, int tile_size_x, int tile_size_y,
                           unsigned char *tile)
{
  int local_px = threadIdx.x + blockIdx.x * blockDim.x;
  int local_py = threadIdx.y + blockIdx.y * blockDim.y;

  if (local_px >= tile_size_x || local_py >= tile_size_y) return;

  int pixel_index = 4 * (local_py * tile_size_x + local_px);

  int global_px = offset_x + local_px;
  int global_py = offset_y + local_py;

  Ray_t ray;
  Color_t pixel_color = Raytracer::RaytraceOne(rtdata, ray, global_px, global_py);

  tile[pixel_index + 0] = pixel_color.fComp.red;
  tile[pixel_index + 1] = pixel_color.fComp.green;
  tile[pixel_index + 2] = pixel_color.fComp.blue;
  tile[pixel_index + 3] = 255;
}

// subdivide image in 16 tiles and launch each tile on a separate CUDA stream
void RenderTiledImage(vecgeom::cuda::RaytracerData_t *rtdata, unsigned char *output_buffer)
{
  cudaStream_t streams[4];

  unsigned char *tile_host[16];
  unsigned char *tile_device[16];

  int tile_size_x = rtdata->fSize_px / 4 + 1;
  int tile_size_y = rtdata->fSize_py / 4 + 1;

  // subdivide each tile in 4x4 thread blocks
  dim3 blocks(tile_size_x / 4 + 1, tile_size_y / 4 + 1), threads(4, 4);

  for (int i = 0; i < 4; ++i)
    checkCudaErrors(cudaStreamCreate(&streams[i]));

  for (int i = 0; i < 16; ++i) {
    // allocate tile on host and device
    checkCudaErrors(cudaMalloc(&tile_device[i], 4 * tile_size_x * tile_size_y));
    // CUDA streams require host memory to be pinned
    checkCudaErrors(cudaMallocHost(&tile_host[i], 4 * tile_size_x * tile_size_y));
  }

  // wait for memory to reach GPU before launching kernels
  checkCudaErrors(cudaDeviceSynchronize());

  for (int ix = 0; ix < 4; ++ix) {
    for (int iy = 0; iy < 4; ++iy) {
      int idx      = 4 * ix + iy;
      int offset_x = ix * tile_size_x;
      int offset_y = iy * tile_size_y;

      RenderTile<<<blocks, threads, 0, streams[iy]>>>(*rtdata, offset_x, offset_y, tile_size_x, tile_size_y,
                                                      tile_device[idx]);
      // copy back rendered tile to system memory
      checkCudaErrors(cudaMemcpyAsync(tile_host[idx], tile_device[idx], (size_t)4 * tile_size_x * tile_size_y,
                                      cudaMemcpyDeviceToHost, streams[iy]));
      checkCudaErrors(cudaFree(tile_device[idx]));
    }
  }

  // ensure all tiles have been copied back
  checkCudaErrors(cudaDeviceSynchronize());
  checkCudaErrors(cudaGetLastError());

  for (int ix = 0; ix < 4; ++ix) {
    for (int iy = 0; iy < 4; ++iy) {
      int idx      = 4 * ix + iy;
      int offset_x = ix * tile_size_x;
      int offset_y = iy * tile_size_y;

      // copy each tile into the final destination
      for (int i = 0; i < tile_size_x; ++i) {
        for (int j = 0; j < tile_size_y; ++j) {
          int px          = offset_x + i;
          int py          = offset_y + j;
          int tile_index  = 4 * (j * tile_size_x + i);
          int pixel_index = 4 * (py * rtdata->fSize_px + px);

          if ((px >= rtdata->fSize_px) || (py >= rtdata->fSize_py)) continue;

          output_buffer[pixel_index + 0] = tile_host[idx][tile_index + 0];
          output_buffer[pixel_index + 1] = tile_host[idx][tile_index + 1];
          output_buffer[pixel_index + 2] = tile_host[idx][tile_index + 2];
          output_buffer[pixel_index + 3] = tile_host[idx][tile_index + 3];
        }
      }
      checkCudaErrors(cudaFreeHost(tile_host[idx]));
    }
  }
  checkCudaErrors(cudaGetLastError());
}

int RaytraceBenchmarkGPU(vecgeom::cuda::RaytracerData_t *rtdata, bool use_tiles)
{
  // Allocate ray data and output data on the device
  size_t statesize = NavStateIndex::SizeOfInstance(rtdata->fMaxDepth);
  size_t raysize   = Ray_t::SizeOfInstance();
  printf(" State size is %lu, ray size is %lu\n", statesize, raysize);

  printf("=== Allocating %.3f MB of ray data on the device\n", (float)rtdata->fNrays * raysize / 1048576);
  // char *input_buffer_gpu = nullptr;
  char *input_buffer = new char[rtdata->fNrays * raysize];
  checkCudaErrors(cudaMallocManaged((void **)&input_buffer, rtdata->fNrays * raysize));

  unsigned char *output_buffer = nullptr;
  checkCudaErrors(
      cudaMallocManaged((void **)&output_buffer, 4 * sizeof(unsigned char) * rtdata->fSize_px * rtdata->fSize_py));

  // Load and synchronize the geometry on the GPU
  auto &cudaManager = cxx::CudaManager::Instance();
  cudaManager.LoadGeometry((vecgeom::cxx::VPlacedVolume *)rtdata->fWorld);
  cudaManager.Synchronize();

  // CudaManager is altering the stack size... setting an appropriate value
  size_t def_stack_limit = 0, def_heap_limit = 0;
  cudaDeviceGetLimit(&def_stack_limit, cudaLimitStackSize);
  cudaDeviceGetLimit(&def_heap_limit, cudaLimitMallocHeapSize);
  std::cout << "=== cudaLimitStackSize = " << def_stack_limit << "  cudaLimitMallocHeapSize = " << def_heap_limit
            << std::endl;
  auto err = cudaDeviceSetLimit(cudaLimitStackSize, 8192);
  cudaDeviceGetLimit(&def_stack_limit, cudaLimitStackSize);
  std::cout << "=== CUDA thread stack size limit set now to: " << def_stack_limit << std::endl;

  auto gpu_world = cudaManager.world_gpu();
  assert(gpu_world && "GPU world volume is a null pointer");

  // Initialize the navigation state for the view point
  NavStateIndex vpstate;
  Raytracer::LocateGlobalPoint(rtdata->fWorld, rtdata->fStart, vpstate, true);
  rtdata->fVPstate = vpstate;
  rtdata->fWorld   = gpu_world;

  rtdata->Print();

  Stopwatch timer;
  timer.Start();

  if (use_tiles) {
    RenderTiledImage(rtdata, output_buffer);
  } else {
    // Construct rays in place
    for (int iray = 0; iray < rtdata->fNrays; ++iray)
      Ray_t::MakeInstanceAt(input_buffer + iray * raysize);

    dim3 blocks(rtdata->fSize_px / 8 + 1, rtdata->fSize_py / 8 + 1), threads(8, 8);
    RenderKernel<<<blocks, threads>>>(*rtdata, input_buffer, output_buffer);
  }

  checkCudaErrors(cudaGetLastError());
  checkCudaErrors(cudaDeviceSynchronize());

  auto time_gpu = timer.Stop();
  std::cout << "Time on GPU: " << time_gpu << "\n";

  write_ppm("output.ppm", output_buffer, rtdata->fSize_px, rtdata->fSize_py);

  checkCudaErrors(cudaFree(input_buffer));
  checkCudaErrors(cudaFree(output_buffer));
  return 0;
}
