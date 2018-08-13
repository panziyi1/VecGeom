#include "RootPersistencyTest.h"

using namespace std;
using namespace vecgeom;

int main(int argc, char *argv[])
{
  cout << "///// Running RootPersistencyTest /////" << endl << endl;

  double dx = 20.;
  double dy = 30.;
  double dz = 40.;

  UnplacedBox worldUnplaced = UnplacedBox(dx * 4, dy * 4, dz * 4);
  UnplacedBox boxUnplaced(dx, dy, dz);
  LogicalVolume world("w0rld", &boxUnplaced);
  LogicalVolume box("b0x", &boxUnplaced);
  Transformation3D placement = Transformation3D(5, 5, 5);
  world.PlaceDaughter(&box, &placement);

  VPlacedVolume *worldPlaced = world.Place();

  GeoManager::Instance().SetWorldAndClose(worldPlaced);

  GeoManager::Instance().GetWorld()->PrintContent();

  cout << "placed vol count: " << GeoManager::Instance().GetPlacedVolumesCount() << endl;
  cout << "registered vol count: " << GeoManager::Instance().GetRegisteredVolumesCount() << endl;

  cout << endl << "writing on vecgeom_export.root";

  RootGeoManager::Instance().Export("vecgeom_export.root");

  cout << endl << "---------------------" << endl;

  cout << "reading from vecgeom_export.root" << endl << endl;

  RootGeoManager::Instance().Import("vecgeom_export.root");

  GeoManager::Instance().GetWorld()->PrintContent();
  cout << endl << "placed vol count: " << GeoManager::Instance().GetPlacedVolumesCount() << endl;
  cout << "registered vol count: " << GeoManager::Instance().GetRegisteredVolumesCount() << endl << endl;

  cout << "printing all logical vol: " << endl;
  for (auto el : GeoManager::Instance().GetLogicalVolumesMap()) {
    cout << el.first << ") ";
    el.second->Print();
    cout << endl;
  }

  cout << endl << endl;
  return 0;
}
