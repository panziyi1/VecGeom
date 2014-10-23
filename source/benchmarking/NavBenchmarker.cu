
#include "base/Global.h"

#include "volumes/PlacedVolume.h"
#include "base/SOA3D.h"
#include "navigation/SimpleNavigator.h"
#include "navigation/NavigationState.h"
#include "management/CudaManager.h"
#include "backend/cuda/Backend.h"

#ifdef VECGEOM_CUDA
#include "base/Stopwatch.h"
#include "backend/cuda/Backend.h"
#include "management/CudaManager.h"
#endif

namespace vecgeom_cuda
{

__global__
void NavKernelTest(
    VPlacedVolume const *const volume,
    const SOA3D<Precision> positions,
    const SOA3D<Precision> directions,
    const int n,
    Precision *const steps) {

  unsigned tid = ThreadIndex(); 
NavigationState old(5), newstate(5);
SimpleNavigator nav;
nav.LocatePoint(volume,positions[tid],old,true);
double step;
nav.FindNextBoundaryAndStep(positions[tid],directions[tid], old, newstate, kInfinity, step);
steps[tid]=step;            
//distance[tid] = volume->SafetyToOut(positions[tid]);

}

} // end of namespace vecgeom_cuda

namespace vecgeom
{

  template <typename Type>
  Type* AllocateAligned(unsigned n) {
    return (Type*) _mm_malloc(n*sizeof(Type), kAlignmentBoundary);
  }

  template <typename Type>
  void FreeAligned(Type *const pmemory) {
    if (pmemory) _mm_free(pmemory);
  }

  void RunNavigationCuda( VPlacedVolume const* volume, unsigned npoints,
                          Precision *const posX, Precision *const posY, Precision *const posZ,
                          Precision *const dirX, Precision *const dirY, Precision *const dirZ,
                          Precision *const distances, Precision *const safeties) {

#ifdef VECGEOM_CUDA
   // transfer geometry to GPU
   typedef vecgeom_cuda::VPlacedVolume const* CudaVolume;
   typedef vecgeom_cuda::SOA3D<Precision> CudaSOA3D;

   CudaManager::Instance().LoadGeometry(volume);
   CudaManager::Instance().Synchronize();
   // std::list<CudaVolume> volumesGpu;
   // for (std::list<VolumePointers>::const_iterator v = fVolumes.begin();
   //      v != fVolumes.end(); ++v) {
   //   volumesGpu.push_back(
   //     reinterpret_cast<CudaVolume>(
   //       CudaManager::Instance().LookupPlaced(v->Specialized())
   //     )
   //   );
   // }

   // copy points to the GPU
   Precision *posXGpu = AllocateOnGpu<Precision>(sizeof(Precision)*npoints);
   Precision *posYGpu = AllocateOnGpu<Precision>(sizeof(Precision)*npoints);
   Precision *posZGpu = AllocateOnGpu<Precision>(sizeof(Precision)*npoints);
   CopyToGpu(posX, posXGpu, npoints*sizeof(Precision));
   CopyToGpu(posY, posYGpu, npoints*sizeof(Precision));
   CopyToGpu(posZ, posZGpu, npoints*sizeof(Precision));
   CudaSOA3D positionGpu = CudaSOA3D(posXGpu, posYGpu, posZGpu, npoints);

   // copy directions to the GPU
   Precision *dirXGpu = AllocateOnGpu<Precision>(sizeof(Precision)*npoints);
   Precision *dirYGpu = AllocateOnGpu<Precision>(sizeof(Precision)*npoints);
   Precision *dirZGpu = AllocateOnGpu<Precision>(sizeof(Precision)*npoints);
   CopyToGpu(dirX, dirXGpu, npoints*sizeof(Precision));
   CopyToGpu(dirY, dirYGpu, npoints*sizeof(Precision));
   CopyToGpu(dirZ, dirZGpu, npoints*sizeof(Precision));
   CudaSOA3D directionGpu = CudaSOA3D(dirXGpu, dirYGpu, dirZGpu, npoints);

   // allocate space for kernel output
   Precision *stepsGpu = AllocateOnGpu<Precision>(sizeof(Precision)*npoints);

   vecgeom_cuda::LaunchParameters launch =
     vecgeom_cuda::LaunchParameters(npoints);

   // launch kernel in GPU
   vecgeom_cuda::Stopwatch timer;
   timer.Start();
   // for (unsigned r = 0; r < fRepetitions; ++r) {
   vecgeom_cuda::NavKernelTest<<< launch.grid_size, launch.block_size>>>(
     *volume, positionGpu, directionGpu, npoints, stepsGpu );
   // }
   Precision elapsedCuda = timer.Stop();

   Precision *const stepsCpu = AllocateAligned<Precision>();
   Precision * steps = (Precision *) _mm_malloc(sizeof(Precision)*np,32);
   CopyFromGpu(stepsGpu, stepsCpu, npoints*sizeof(Precision));

   FreeFromGpu(stepsGpu);
   FreeFromGpu(posXGpu);
   FreeFromGpu(posYGpu);
   FreeFromGpu(posZGpu);
   FreeFromGpu(dirXGpu);
   FreeFromGpu(dirYGpu);
   FreeFromGpu(dirZGpu);

   // compare steps from navigator with the ones above
   std::cout<<"TODO: compare steps from navigator with the previous ones.\n";
#endif  // VECGEOM_CUDA
  }

}
