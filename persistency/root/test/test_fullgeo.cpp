#include "test_fullgeo.hh"

using namespace vecgeom;
using namespace std;

int main(int argc, char *argv[])
{
  if (argc > 1 && strcmp(argv[1], "w") == 0) { // if w -> write the file

    TFile fo("geo_full.root", "RECREATE");

    cout << "writing on geo_full.root:";

    UnplacedBox worldUnplaced = UnplacedBox(4, 8, 6);
    UnplacedBox boxUnplaced   = UnplacedBox(4, 8, 6);

    LogicalVolume world("world", &worldUnplaced);
    LogicalVolume box("box", &boxUnplaced);

    Transformation3D placement(0.1, 0, 0);
    world.PlaceDaughter("box", &box, &placement);

    VPlacedVolume *worldPlaced = world.Place();

    GeoManager::Instance().SetWorldAndClose(worldPlaced);

    world.Print();

    fo.WriteObject(&world, "world_saved");
  } else {
    TFile fi("geo_full.root");

    LogicalVolume *world;

    fi.GetObject("world_saved", world);

    cout << "\n\t";
    world->Print();
  }

  return 0;
}
