#include<iostream>
#include<dlfcn.h>
#include "volumes/PlacedVolume.h" 

using namespace vecgeom; 
int main()
{
    void *handle;
    handle = dlopen("libCMSGeom.so", RTLD_NOW);
    if (!handle)
    {
      std::cout << "The error is " << dlerror() << "\n";
    }

    // the Create Detector function type:
    typedef VPlacedVolume const * (*CreateFunc_t)();
  
    CreateFunc_t f = (CreateFunc_t) dlsym(handle,"_Z16generateDetectorv");

    // call the create function
    vecgeom::VPlacedVolume const * world = f();

    // check if world has your expected label
    std::cout << world->GetLabel() << "\n";

    // from here on the geometry is available via the world
    // GeoManager has to be closed ....

    return 0 ;
}
