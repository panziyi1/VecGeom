#ifndef navgpu_h
#define navgpu_h

#include "base/Global.h"

namespace vecgeom {
  class VPlacedVolume;
 
 void RunNavigationCuda( VPlacedVolume const* volume, unsigned npoints,
                          Precision *const posX, Precision *const posY, Precision *const posZ,
                          Precision *const dirX, Precision *const dirY, Precision *const dirZ,
                         Precision const *const pSteps,      Precision *const steps );

}

#endif
