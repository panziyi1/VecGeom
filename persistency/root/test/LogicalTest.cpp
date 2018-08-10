#include "RootPersistencyTest.h"

using namespace std;
using namespace vecgeom;

void logical_test()
{
  cout << "Running logical_test" << endl << endl;

  RootPersistencyProxy(); // calling the proxy

  UnplacedBox worldUnplaced = UnplacedBox(4, 8, 6);
  UnplacedBox boxUnplaced   = UnplacedBox(4, 8, 6);
  LogicalVolume world("world", &worldUnplaced);
  LogicalVolume box("box", &boxUnplaced);
  Transformation3D placement(0.1, 0, 0);

  char name[10];
  for (auto i = 0; i < 10; ++i) {
    sprintf(name, "box_%d", i);
    world.PlaceDaughter(name, &box, &placement);
  }

  VPlacedVolume *worldPlaced = world.Place();

  cout << "writing on logical.root" << endl << endl;

  TFile fo("logical.root", "RECREATE");
  fo.WriteObject(&world, "world_saved");
  fo.WriteObject(worldPlaced, "worldP_saved");
  fo.Close();

  cout << "reading from logical.root" << endl << endl;
  TFile fi("logical.root");

  LogicalVolume *rworld;
  VPlacedVolume *rworldPlaced;

  fi.GetObject("world_saved", rworld);
  fi.GetObject("worldP_saved", rworldPlaced);

  world.Print();
  cout << endl << "------------------------------" << endl;
  rworld->Print();
  cout << endl << endl;

  cout << "Number of daughters: " << rworld->GetDaughters().size() << ", expected: " << world.GetDaughters().size()
       << endl
       << endl;
  for (auto placed : rworld->GetDaughters())
    placed->Print();

  cout << endl << endl;

  worldPlaced->Print();
  cout << endl << "------------------------------" << endl;
  rworldPlaced->Print();
  cout << endl << endl;
}