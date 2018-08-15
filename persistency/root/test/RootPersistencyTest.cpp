#include "RootPersistencyTest.h"

using namespace std;
using namespace vecgeom;

int main(int argc, char *argv[])
{
  cout << "\033[1;34m///// Running RootPersistencyTest /////\033[0m" << endl << endl;

  UnplacedBox worldUnplaced(80., 80., 80.);

  UnplacedBox boxUnplaced(40., 40., 40.);

  UnplacedParaboloid parUnplaced(6., 9., 10.);

  UnplacedParallelepiped palUnplaced(40., 40., 40., .2, .5, .1);

  UnplacedSphere sphUnplaced(.2, .5);

  auto trdUnplaced = GeoManager::MakeInstance<UnplacedTrd>(0.14999999999999999, 0.14999999999999999, 24.707000000000001,
                                                           24.707000000000001, 22.699999999999999);

  using vec3dd       = vecgeom::Vector3D<double>;
  vec3dd trapvert[8] = {vec3dd(-10.0, -20.0, -40.0), vec3dd(+10.0, -20.0, -40.0), vec3dd(-10.0, +20.0, -40.0),
                        vec3dd(+10.0, +20.0, -40.0), vec3dd(-30.0, -40.0, +40.0), vec3dd(+30.0, -40.0, +40.0),
                        vec3dd(-30.0, +40.0, +40.0), vec3dd(+30.0, +40.0, +40.0)};
  auto trapUnplaced  = GeoManager::MakeInstance<UnplacedTrapezoid>(trapvert);

  double verticesx1[8] = {-3, -3, 3, 3, -2, -2, 2, 2};
  double verticesy1[8] = {-3, 3, 3, -3, -2, 2, 2, -2};
  UnplacedGenTrap gentrpUnplaced(verticesx1, verticesy1, 5);

  auto hypeUnplaced = GeoManager::MakeInstance<UnplacedHype>(10, 20, kPi / 4, kPi / 3, 50);
  UnplacedOrb orbUnplaced(9.);

  double thb  = 3 * kPi / 4;
  double phib = kPi / 3;
  double tht  = kPi / 4;
  double phit = 2 * kPi / 3;
  Vector3D<double> nbottom(std::sin(thb) * std::cos(phib), std::sin(thb) * std::sin(phib), std::cos(thb));
  Vector3D<double> ntop(std::sin(tht) * std::cos(phit), std::sin(tht) * std::sin(phit), std::cos(tht));
  auto cuttUnplaced = GeoManager::MakeInstance<UnplacedCutTube>(0, 5., 10., 0., 2 * kPi, nbottom, ntop);

  LogicalVolume world("world", &boxUnplaced);
  LogicalVolume box("box", &boxUnplaced);
  LogicalVolume par("par", &parUnplaced);
  LogicalVolume pal("pal", &palUnplaced);
  LogicalVolume sph("sph", &sphUnplaced);
  LogicalVolume trd("trd", trdUnplaced);
  LogicalVolume trap("trap", trapUnplaced);
  LogicalVolume gentrp("gentrp", &gentrpUnplaced);
  LogicalVolume hype("hype", hypeUnplaced);
  LogicalVolume orb("orb", &orbUnplaced);
  LogicalVolume cutt("cutt", cuttUnplaced);

  Transformation3D placement = Transformation3D(5, 5, 5);
  Transformation3D placement2 = Transformation3D(0, 0, 0, 90, 0, 0);
  world.PlaceDaughter(&box, &placement);
  world.PlaceDaughter(&par, &placement);
  world.PlaceDaughter(&pal, &placement2);
  world.PlaceDaughter(&sph, &placement);
  world.PlaceDaughter(&trd, &placement);
  world.PlaceDaughter(&trap, &placement);
  world.PlaceDaughter(&gentrp, &placement);
  world.PlaceDaughter(&hype, &placement2);
  world.PlaceDaughter(&orb, &placement2);
  world.PlaceDaughter(&cutt, &placement2);

  VPlacedVolume *worldPlaced = world.Place();

  GeoManager::Instance().SetWorld(worldPlaced);

  GeoManager::Instance().GetWorld()->PrintContent();

  cout << endl << "placed vol count: " << GeoManager::Instance().GetPlacedVolumesCount() << endl;
  cout << "registered vol count: " << GeoManager::Instance().GetRegisteredVolumesCount() << endl;

  cout << endl << "\033[0;31mwriting on vecgeom_export.root\n" << endl;

  RootGeoManager::Instance().Export("vecgeom_export.root");

  cout << endl << "\033[1;30m---------------------\033[0;34m" << endl;

  cout << "reading from vecgeom_export.root\n" << endl;

  RootGeoManager::Instance().Import("vecgeom_export.root");
  cout << "\033[0m" << endl;
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
