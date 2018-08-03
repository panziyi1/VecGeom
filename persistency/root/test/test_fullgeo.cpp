#include "test_fullgeo.hh"

using namespace vecgeom;
using namespace std;

int main(int argc, char *argv[])
{
  proxy();
  if (argc > 1 && strcmp(argv[1], "w") == 0) { // if w -> write the file

    TFile fo("geo_full.root", "RECREATE");

    cout << "writing on geo_full.root:" << endl << endl;

    UnplacedBox worldUnplaced = UnplacedBox(4, 8, 6);
    UnplacedBox boxUnplaced   = UnplacedBox(4, 8, 6);

    LogicalVolume world("world", &worldUnplaced);
    LogicalVolume box("box", &boxUnplaced);
    // LogicalVolume box2("box2", &boxUnplaced);

    Transformation3D placement(0.1, 0, 0);
    world.PlaceDaughter("box", &box, &placement);
    // world.PlaceDaughter("box2", &box2, &placement);

    VPlacedVolume *worldPlaced = world.Place();

    world.Print();
    cout << endl << endl;
    worldPlaced->Print();

    fo.WriteObject(&world, "world_saved");
    fo.WriteObject(worldPlaced, "worldP_saved");

    // GeoManager::Instance().SetWorldAndClose(worldPlaced);

  } else {
    TFile fi("geo_full.root");

    LogicalVolume *world;
    VPlacedVolume *worldPlaced;

    cout << "reading from geo_full.root:" << endl << endl;

    fi.GetObject("world_saved", world);
    fi.GetObject("worldP_saved", worldPlaced);

    world->Print();
    worldPlaced->Print();
  }

  return 0;
}
