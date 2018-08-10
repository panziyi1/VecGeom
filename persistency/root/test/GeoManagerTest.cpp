#include "GeoManagerTest.h"

using namespace std;
using namespace vecgeom;

void write()
{
  cout << "/////Running GeoManagerTest in write mode/////" << endl << endl;

  RootPersistencyProxy(); // calling the proxy

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

  cout << "writing on geom.root" << endl << endl;
  TFile fo("geom.root", "RECREATE");
  fo.WriteObject(&GeoManager::Instance(), "geom_saved");
  fo.Close();

  GeoManager::Instance().GetWorld()->Print();

  cout << endl << "---------------------" << endl;
  cout << "placed vol count: " << GeoManager::Instance().GetPlacedVolumesCount() << endl;
  cout << "registered vol count: " << GeoManager::Instance().GetRegisteredVolumesCount() << endl;

  cout << endl << endl;
}

void read()
{

  cout << endl << "/////Running GeoManagerTest in read mode/////" << endl << endl;

  cout << "reading from geom.root" << endl << endl;
  TFile fi("geom.root");

  GeoManager *geom;
  fi.GetObject("geom_saved", geom);

  cout << "comparing geom pointer with &GeoManager::Instance():  " << geom << " = " << &GeoManager::Instance() << endl
       << endl;

  GeoManager::Instance().GetWorld()->Print();

  cout << endl << "---------------------" << endl;
  cout << "placed vol count: " << GeoManager::Instance().GetPlacedVolumesCount() << endl;
  cout << "registered vol count: " << GeoManager::Instance().GetRegisteredVolumesCount() << endl << endl;

  cout << "printing all logical vol: " << endl;
  for (auto el : GeoManager::Instance().GetLogicalVolumesMap()) {
    cout << el.first << ") ";
    el.second->Print();
    cout << endl;
  }
  cout << endl << endl;
}
void usage()
{
  cout << "Usage:" << endl;
  cout << "\t./GeoManagerTest w && ./GeoManagerTest r" << endl << endl;
}

int main(int argc, char *argv[])
{
  if (argc != 2) {
    usage();
    return 0;
  }
  RootPersistencyProxy();
  switch (argv[1][0]) {
  case 'w':
    write();
    break;
  case 'r':
    read();
    break;
  default:
    usage();
  }
  return 0;
}
