
#include "base/Global.h"

#include "volumes/PlacedVolume.h"
#include "base/SOA3D.h"
#include "navigation/SimpleNavigator.h"
#include "navigation/NavigationState.h"
#include "management/CudaManager.h"
#include "backend/cuda/Backend.h"


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

}
