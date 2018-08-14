#include "RootPersistencyComponentsTest.h"

using namespace std;
using namespace vecgeom;

bool logical_test()
{
  cout << "///// Running logical_test /////" << endl << endl;

  RootPersistencyProxy(); // calling the proxy

  UnplacedBox worldUnplaced = UnplacedBox(4, 8, 6);
  UnplacedBox boxUnplaced   = UnplacedBox(4, 18, 6);
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

  // testing
  bool all_test_ok = true;
  bool test_ok;

  test_ok = true;
  cout << "[1] comparing LogicalVolume *world\n\n"
       << ">> before\n";
  world.PrintContent();
  cout << "\n\n"
       << ">> after\n";
  rworld->PrintContent();
  cout << "\n" << endl;
  auto daughters  = world.GetDaughters();
  auto rdaughters = rworld->GetDaughters();
  if (rworld->GetDaughters().size() != world.GetDaughters().size()) {
    cout << "number of daughters don't match\n";
    test_ok = false;
  }
  for (auto vol = daughters.begin(), rvol = rdaughters.begin(), volEnd = daughters.end(), rvolEnd = rdaughters.end();
       vol != volEnd && rvol != rvolEnd; ++vol, ++rvol) {
    if ((*vol)->GetLabel() != (*rvol)->GetLabel()) {
      cout << "label doesn't match\n";
      test_ok = false;
    }
    if (!(*((*vol)->GetTransformation()) == *((*rvol)->GetTransformation()))) {
      cout << "transformation doesn't match\n";
      test_ok = false;
    }
    if (((UnplacedBox *)((*vol)->GetLogicalVolume()->GetUnplacedVolume()))->dimensions() !=
        ((UnplacedBox *)((*rvol)->GetLogicalVolume()->GetUnplacedVolume()))->dimensions()) {
      cout << "dimensions of unplaced doesn't match\n";
      test_ok = false;
    }
  }

  if (test_ok)
    cout << "test passed\n\n" << endl;
  else {
    cout << "! test not passed\n\n" << endl;
    all_test_ok = false;
  }

  cout << endl << endl;
  return all_test_ok;
}