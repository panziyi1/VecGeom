#include "base/Global.h"
#include "navigation/NavigationState.h"
#include "navigation/NavStatePool.h"
#include "management/CudaManager.h"
#include "navigation/SimpleNavigator.h"
#include "backend/cuda/Backend.h"

#include <stdio.h>

__global__ 
void ProcessNavStates( void* gpu_ptr /* a pointer to buffer of navigation states */, int depth, int n )
{
  const int i = vecgeom::cuda::ThreadIndex();
  if( i >= n ) return;

  if(i==0){
    printf("*** Size of NavigationState on the GPU: %d bytes\n", sizeof(vecgeom::cuda::NavigationState));
  }

  // get the navigationstate for this thread/lane
  using vecgeom::cuda::NavigationState;
  using vecgeom::cuda::NavStatePool;
  NavStatePool &states  = *reinterpret_cast<NavStatePool*>(gpu_ptr);
  NavigationState *state  = states[i];

  // vecgeom::cuda::NavigationState * state = reinterpret_cast<vecgeom::cuda::NavigationState*>( gpu_ptr +
  //       vecgeom::cuda::NavigationState::SizeOf(depth)*i );

  // actually do something to the states; here just popping off the top volume
  state->Pop();

   // state->SetPathPointer( oldpathpointer );
}

void LaunchNavigationKernel( void* gpu_ptr, int depth, int n )
{
  vecgeom::cuda::LaunchParameters launch =
      vecgeom::cuda::LaunchParameters(n);

  ProcessNavStates<<< launch.grid_size , launch.block_size >>>( gpu_ptr, depth, n );
}


