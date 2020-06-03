/// \file Raytracer.cu
/// \author Andrei Gheata (andrei.gheata@cern.ch)

#include <VecGeom/base/Transformation3D.h>
#include <VecGeom/management/GeoManager.h>
#include <VecGeom/management/CudaManager.h>
#include <VecGeom/navigation/NavigationState.h>
#include <VecGeom/volumes/PlacedVolume.h>
#include <VecGeom/benchmarking/Raytracer.h>

#include <cassert>
#include <cstdio>

using namespace vecgeom;

void check_cuda(cudaError_t result, char const *const func, const char *const file, int const line)
{
  if (result) {
    fprintf(stderr, "CUDA error = %s at %s:%d\n", cudaGetErrorString(result), file, line);
    cudaDeviceReset();
    exit(1);
  }
}

#define checkCudaErrors(val) check_cuda((val), #val, __FILE__, __LINE__)

namespace visitorcuda {

class GlobalToLocalVisitor {
private:
  int fType = 0; ///< visit type: 0=no computation 1=master-to-local NavStatePath 2=master-to-local NavStateIndex
  Vector3D<Precision> fGlobalPoint; ///< point coordinates in the global frame
  Vector3D<Precision> fLocal;       ///< make sure we store the result somewhere

public:
  GlobalToLocalVisitor() {}

  VECCORE_ATT_HOST_DEVICE
  void SetType(int type) { fType = type; }

  VECCORE_ATT_HOST_DEVICE
  void SetGlobalPoint(Precision x, Precision y, Precision z) { fGlobalPoint.Set(x, y, z); }

  VECCORE_ATT_HOST_DEVICE
  void apply(NavStatePath *state, NavIndex_t nav_ind)
  {
    switch (fType) {
    case 0:
      break;
    case 1:
      fLocal = state->GlobalToLocal(fGlobalPoint);
      break;
    case 2:
      fLocal = NavStateIndex::GlobalToLocalImpl(nav_ind, fGlobalPoint);
      break;
    }
  }
};

/// Traverses the geometry tree keeping track of the state context (volume path or navigation state)
/// and applies the injected Visitor
VECCORE_ATT_HOST_DEVICE
template <typename Visitor>
void visitAllPlacedVolumesPassNavIndex(VPlacedVolume const *currentvolume, Visitor *visitor, NavStatePath *state,
                                       NavIndex_t nav_ind)
{
  if (currentvolume != NULL) {
    state->Push(currentvolume);
    visitor->apply(state, nav_ind);
    for (auto daughter : currentvolume->GetDaughters()) {
      auto nav_ind_d = NavStateIndex::PushImpl(nav_ind, daughter);
      visitAllPlacedVolumesPassNavIndex(daughter, visitor, state, nav_ind_d);
    }
    state->Pop();
  }
}

} // namespace visitorcuda

__global__
TestNavIndexGPUKernel(vecgeom::cuda::VPlacedVolume const* const gpu_world, vecgeom::cuda::NavStatePath * const state, int type, int npasses)
{
  using namespace visitorcuda;
  
  state->Clear();
  GlobalToLocalVisitor visitor;

  NavIndex_t nav_ind_top = 1; // The navigation index corresponding to the world
  visitor.SetType(type);

  for (auto i = 0; i < npasses; ++i)
    visitAllPlacedVolumesPassNavIndex(world, &visitor, state, nav_ind_top);
}



void TestNavIndexGPU(vecgeom::cuda::VPlacedVolume const* const world, int npasses)
{
  // Load and synchronize the geometry on the GPU
  int maxdepth = vecgeom::cxx::GeoManager::Instance().getMaxDepth();
  size_t statesize = NavigationState::SizeOfInstance(maxdepth);

  vecgeom::cxx::CudaManager::Instance().LoadGeometry((vecgeom::cxx::VPlacedVolume*) world);
  vecgeom::cxx::CudaManager::Instance().Synchronize();
  
  auto gpu_world = vecgeom::cxx::CudaManager::Instance().world_gpu();
  assert(gpu_world && "GPU world volume is a null pointer");

  char *input_buffer = nullptr;
  checkCudaErrors(cudaMallocManaged((void **)&input_buffer, statesize));
  auto state = NavStatePath::MakeInstanceAt(maxdepth, (void *)(input_buffer));


  Stopwatch timer;
  timer.Start();
  TestNavIndexGPU<<<1, 1>>>(gpu_world, state, 0, npasses);
  auto tbaseline = timer.Stop();

  timer.Start();
  TestNavIndexGPU<<<1, 1>>>(gpu_world, state, 1, npasses);
  auto tnavstate = timer.Stop();
  std::cout << "NavStatePath::GlobalToLocal took: " << tnavstate - tbaseline << " sec.\n";

  timer.Start();
  TestNavIndexGPU<<<1, 1>>>(gpu_world, state, 2, npasses);
  auto tnavindex = timer.Stop();
  std::cout << "NavStateIndex::GlobalToLocal took: " << tnavindex - tbaseline << " sec.\n";

  std::cout << "Speedup per GlobalToLocal averaged over all states: " << std::setprecision(3)
            << (tnavstate - tbaseline) / (tnavindex - tbaseline) << "\n";

  checkCudaErrors(cudaGetLastError());
  checkCudaErrors(cudaDeviceSynchronize());

  checkCudaErrors(cudaFree(input_buffer));
}
