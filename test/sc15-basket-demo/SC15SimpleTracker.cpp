/// author: Andrei Gheata (andrei.gheata@cern.ch) Oct 14, 2015
#include <sstream>
#include <string>
#include "volumes/Tube.h"
#include "volumes/LogicalVolume.h"
#include "management/GeoManager.h"
#include "SC15Navigators.h"
#include <iostream>
using namespace vecgeom;

#define SETINNERNAV(depth) \
  if(layer==depth) vol->SetUserExtensionPtr( (void*) InnerMostTubeNavigator<depth>::Instance() );

#define SETLAYERNAV(depth) \
  if(layer==depth) vol->SetUserExtensionPtr( (void*) LayerNavigator<depth>::Instance() );

#define SETWORLDNAV(depth) \
  if(layer==0) vol->SetUserExtensionPtr( (void*) WorldNavigator<true>::Instance() );

void AssignNavigatorToVolume( LogicalVolume *vol, int layer, int maxlayers ){
  assert( maxlayers <= 10 ); // we are only listing 10 template specializations up to depth 10 here
  if( layer == maxlayers - 1 ) // it is the last layer
  {
    SETINNERNAV(0);
    SETINNERNAV(1);
    SETINNERNAV(2);
    SETINNERNAV(3);
    SETINNERNAV(4);
    SETINNERNAV(5);
    SETINNERNAV(6);
    SETINNERNAV(7);
    SETINNERNAV(8);
    SETINNERNAV(9);
    SETINNERNAV(10);
  }
  if( layer < maxlayers - 1 ){
    SETLAYERNAV(0);
    SETLAYERNAV(1);
    SETLAYERNAV(2);
    SETLAYERNAV(3);
    SETLAYERNAV(4);
    SETLAYERNAV(5);
    SETLAYERNAV(6);
    SETLAYERNAV(7);
    SETLAYERNAV(8);
    SETLAYERNAV(9);
    SETLAYERNAV(10);
  }
}


VPlacedVolume *CreateSimpleTracker(int nlayers) {
  std::cout << "Creating SimpleTracker geometry having " << nlayers << " layers"
            << std::endl;
  // World size
  const double world_size = 50.;

  // Top volume
  UnplacedBox *uTop = new UnplacedBox(world_size, world_size, world_size);
  LogicalVolume *top = new LogicalVolume("world", uTop);
  top->SetUserExtensionPtr( (void *) WorldNavigator<true>::Instance() );
  
// Cylindrical layers
  double rmax;
  std::string layerBase = "layer_";
  double deltaR = world_size/nlayers;

  LogicalVolume *mother = top;
  for (auto layer=0; layer<nlayers; ++layer) {
    rmax = world_size - layer * deltaR;
    std::ostringstream layerName;
    layerName << layerBase << layer;
    UnplacedTube *uLayer = new UnplacedTube(0, rmax, world_size, 0, kTwoPi);
    LogicalVolume *layerVol = new LogicalVolume(layerName.str().c_str(), uLayer);
    
    AssignNavigatorToVolume(layerVol, layer, nlayers);

    // Place in mother
    mother->PlaceDaughter(layerName.str().c_str(), layerVol, &Transformation3D::kIdentity);
  }  

  VPlacedVolume *world = top->Place();
  GeoManager::Instance().SetWorld(world);
  GeoManager::Instance().CloseGeometry();
  return world;
}

int main(int argc, char* argv[]) {
  int nlayers = 10;
  if (argc > 1) {
    nlayers = std::atoi(argv[1]);
  }
  CreateSimpleTracker(nlayers);

  // loop over all logical volumes and print navigator
  std::vector<LogicalVolume *> lvols;
  GeoManager::Instance().GetAllLogicalVolumes(lvols);

  for( auto v : lvols ){
    auto nav = (VNavigator*) v->GetUserExtensionPtr();
    std::cerr << v->GetName() << " has navigator " << nav->GetName() << "\n";
  }

}
