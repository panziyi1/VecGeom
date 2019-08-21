#include "RootPersistencyTest.h"

using namespace std;
using namespace vecgeom;

void write()
{
  cout << "\033[1;34m///// Running RootPersistencyTest /////\033[0m" << endl << endl;

  UnplacedBox worldUnplaced(80., 80., 80.);

  UnplacedBox boxUnplaced(40., 40., 40.);

  UnplacedParaboloid parUnplaced(6., 9., 10.);

  UnplacedParallelepiped palUnplaced(40., 40., 40., .2, .5, .1);

  UnplacedSphere sphUnplaced(.2, .5);

  auto trdUnplaced = GeoManager::MakeInstance<UnplacedTrd>(0.14999999999999999, 0.14999999999999999, 24.707000000000001,
                                                           24.707000000000001, 22.699999999999999);

  auto hypeUnplaced = GeoManager::MakeInstance<UnplacedHype>(10, 20, kPi / 4, kPi / 3, 50);

  using vec3dd       = vecgeom::Vector3D<double>;
  vec3dd trapvert[8] = {vec3dd(-10.0, -20.0, -40.0), vec3dd(+10.0, -20.0, -40.0), vec3dd(-10.0, +20.0, -40.0),
                        vec3dd(+10.0, +20.0, -40.0), vec3dd(-30.0, -40.0, +40.0), vec3dd(+30.0, -40.0, +40.0),
                        vec3dd(-30.0, +40.0, +40.0), vec3dd(+30.0, +40.0, +40.0)};
  auto trapUnplaced  = GeoManager::MakeInstance<UnplacedTrapezoid>(trapvert);
  UnplacedOrb orbUnplaced(9.);

  double thb  = 3 * kPi / 4;
  double phib = kPi / 3;
  double tht  = kPi / 4;
  double phit = 2 * kPi / 3;
  Vector3D<double> nbottom(std::sin(thb) * std::cos(phib), std::sin(thb) * std::sin(phib), std::cos(thb));
  Vector3D<double> ntop(std::sin(tht) * std::cos(phit), std::sin(tht) * std::sin(phit), std::cos(tht));
  auto cuttUnplaced = GeoManager::MakeInstance<UnplacedCutTube>(0, 5., 10., 0., 2 * kPi, nbottom, ntop);

  LogicalVolume *lb     = new LogicalVolume("boxlv", new UnplacedBox(10., 10., 10.));
  UnplacedAssembly *ass = new UnplacedAssembly();
  LogicalVolume *assLv  = new LogicalVolume("assemblylv", ass);
  ass->AddVolume(lb->Place(new Transformation3D(-20., 0., 0.)));
  ass->AddVolume(lb->Place(new Transformation3D(20., 0., 0.)));
  VPlacedVolume *const assPlaced = assLv->Place();

  auto tubeUnplaced = GeoManager::MakeInstance<UnplacedTube>(45, 50, 50, 0, 2 * kPi);

  auto coneUnplaced = GeoManager::MakeInstance<UnplacedCone>(10, 20, 15, 25, 100, 0, 2. * M_PI);

  double L = 10.;
  LogicalVolume ltube("tube", new GenericUnplacedTube(0., 0.9 * L / 4., L, 0., vecgeom::kTwoPi));
  auto upperhole = ltube.Place(new Transformation3D(-L / 4, -L / 4, 0.));
  auto lowerhole = ltube.Place(new Transformation3D(L / 4, L / 4, 0.));
  UnplacedBooleanVolume<kUnion> unionUnplaced(kUnion, upperhole, lowerhole);
  UnplacedBooleanVolume<kIntersection> intersectionUnplaced(kIntersection, upperhole, lowerhole);
  UnplacedBooleanVolume<kSubtraction> subtractionUnplaced(kSubtraction, upperhole, lowerhole);

  double RMINVec0[2];
  RMINVec0[0] = 1;
  RMINVec0[1] = 1;
  double RMAXVec0[2];
  RMAXVec0[0] = 2;
  RMAXVec0[1] = 2;
  double Z_Values0[2];
  Z_Values0[0]       = -1;
  Z_Values0[1]       = 1;
  auto polyhUnplaced = GeoManager::MakeInstance<UnplacedPolyhedron>(0.0, kPi, 2, 2, Z_Values0, RMINVec0, RMAXVec0);
  
  double rmin[] = {0.1, 0., 0., 0.2};
  double rmax[] = {1., 2., 2., 1.5};
  double z[]    = {-1, -0.5, 0.5, 10};
  auto pconUnplaced = GeoManager::MakeInstance<UnplacedPolycone>(0, kTwoPi, 4, z, rmin, rmax);

  UnplacedScaledShape scaledUnplaced(tubeUnplaced, 0.5, 1.3, 1.);

  Vector3D<double> p0(0., 0., 2.), p1(0., 0., 0.), p2(2., 0., 0.), p3(0., 2., 0.);
  UnplacedTet tetUnplaced = UnplacedTet(p0, p1, p2, p3);
  UnplacedTorus2 torusUnplaced(1.2, 3.1, 5., 0, kTwoPi);
/*

    double verticesx1[8] = {-3, -3, 3, 3, -2, -2, 2, 2};
    double verticesy1[8] = {-3, 3, 3, -3, -2, 2, 2, -2};
    UnplacedGenTrap gentrpUnplaced(verticesx1, verticesy1, 5);

  

    UnplacedMultiUnion multiuUnplaced;
    double sized = 10. * std::pow(0.5 / 10, 1. / 3.);
    for (size_t i = 0; i < 10; ++i) {
      Vector3D<double> pos(RNG::Instance().uniform(-10., 10.), RNG::Instance().uniform(-10., 10.),
                           RNG::Instance().uniform(-10., 10.));
      double sizernd = RNG::Instance().uniform(0.8 * sized, 1.2 * sized);
      Transformation3D trans(pos.x(), pos.y(), pos.z(), RNG::Instance().uniform(-180, 180),
                             RNG::Instance().uniform(-180, 180), RNG::Instance().uniform(-180, 180));
      trans.SetProperties();
      UnplacedBox *box = new UnplacedBox(sizernd, sizernd, sizernd);
      multiuUnplaced.AddNode(box, trans);
    }
    multiuUnplaced.Close();



    double x[12], y[12];
    for (size_t i = 0; i < (size_t)12; ++i) {
      x[i] = 5 * std::sin(i * (2. * M_PI) / 12);
      y[i] = 5 * std::cos(i * (2. * M_PI) / 12);
    }
    UnplacedSExtruVolume sextruUnplaced(12, x, y, -5, 5);

    UnplacedTessellated tslUnplaced = UnplacedTessellated();
    TessellatedOrb(10., 10, tslUnplaced);
    tslUnplaced.Close();

  */
  LogicalVolume world("world", &boxUnplaced);
  LogicalVolume box("box", &boxUnplaced);
  LogicalVolume par("par", &parUnplaced);
  LogicalVolume pal("pal", &palUnplaced);
  LogicalVolume sph("sph", &sphUnplaced);
  LogicalVolume trd("trd", trdUnplaced);
  LogicalVolume hype("hype", hypeUnplaced);
  LogicalVolume trap("trap", trapUnplaced);
  LogicalVolume orb("orb", &orbUnplaced);
  LogicalVolume cutt("cutt", cuttUnplaced);
  LogicalVolume tube("tube", tubeUnplaced);
  LogicalVolume cone("cone", coneUnplaced);
  LogicalVolume lunion("union", &unionUnplaced);
  LogicalVolume lintersection("intersection", &intersectionUnplaced);
  LogicalVolume lsubtraction("subtraction", &subtractionUnplaced);
  LogicalVolume pcon("pcon", pconUnplaced);
  LogicalVolume polyh("polyh", polyhUnplaced);
  LogicalVolume scaled("scaled", &scaledUnplaced);
  LogicalVolume tet("tet", &tetUnplaced);
  LogicalVolume torus("torus", &torusUnplaced);
  /*
    LogicalVolume gentrp("gentrp", &gentrpUnplaced);
    LogicalVolume xtru("xtru", ExtrudedMultiLayer(false));
    LogicalVolume multiu("multiu", &multiuUnplaced);
    LogicalVolume sextru("sextru", &sextruUnplaced);
    LogicalVolume tsl("tsl", &tslUnplaced);
  */

  Transformation3D placement  = Transformation3D(5, 5, 5);
  Transformation3D placement2 = Transformation3D(0, 0, 0, 90, 0, 0);
  Transformation3D origin(0, 0, 0);

  world.PlaceDaughter(&box, &placement);
  world.PlaceDaughter(&par, &placement);
  world.PlaceDaughter(&pal, &placement2);
  world.PlaceDaughter(&sph, &placement);
  world.PlaceDaughter(&trd, &placement);
  world.PlaceDaughter(&hype, &placement2);
  world.PlaceDaughter(&trap, &placement);
  world.PlaceDaughter(&orb, &placement2);
  world.PlaceDaughter(&cutt, &placement2);
  world.PlaceDaughter(assPlaced);
  world.PlaceDaughter(&tube, &placement2);
  world.PlaceDaughter(&cone, &placement);
  world.PlaceDaughter(&lunion, new Transformation3D(-L / 2, -L / 2, 0.));
  world.PlaceDaughter(&lintersection, new Transformation3D(-L / 2, -L / 2, 0.));
  world.PlaceDaughter(&lsubtraction, new Transformation3D(-L / 2, -L / 2, 0.));
  world.PlaceDaughter(&pcon, &placement2);
  world.PlaceDaughter(&polyh, &placement2);
  world.PlaceDaughter(&scaled, &placement);
  world.PlaceDaughter(&tet, &placement);
  world.PlaceDaughter(&torus, &placement2);
  /*
    world.PlaceDaughter(&gentrp, &placement);
    world.PlaceDaughter(&xtru, &placement2);
    world.PlaceDaughter(&multiu, &placement);
    world.PlaceDaughter(&sextru, &placement);
    world.PlaceDaughter(&tsl, &origin);
  */
  VPlacedVolume *worldPlaced = world.Place();

  GeoManager::Instance().SetWorldAndClose(worldPlaced);

  /*
    for (auto &lvol : GeoManager::Instance().GetLogicalVolumesMap()) {
      if (lvol.second->GetDaughtersp()->size() < 4) {
        lvol.second->SetNavigator(NewSimpleNavigator<>::Instance());
      } else if (lvol.second->GetDaughtersp()->size() < 10) {
        lvol.second->SetNavigator(SimpleABBoxNavigator<>::Instance());
      } else {
        lvol.second->SetNavigator(HybridNavigator<>::Instance());
        HybridManager2::Instance().InitStructure((lvol.second));
      }

      if (lvol.second->ContainsAssembly()) {
        lvol.second->SetLevelLocator(SimpleAssemblyAwareABBoxLevelLocator::GetInstance());
      } else {
        lvol.second->SetLevelLocator(SimpleABBoxLevelLocator::GetInstance());
      }
    }
  */
  GeoManager::Instance().GetWorld()->PrintContent();

  cout << endl << "placed vol count: " << GeoManager::Instance().GetPlacedVolumesCount() << endl;
  cout << "registered vol count: " << GeoManager::Instance().GetRegisteredVolumesCount() << endl;

  cout << endl << "\033[0;31mwriting on vecgeom_export.root\n" << endl;

  // gDebug = 2;
  RootGeoManager::Instance().Export("vecgeom_export.root");

  cout << "\033[0m" << endl;
  GeoManager::Instance().Clear();
}

void read()
{
  cout << endl << "\033[0;31mreading from vecgeom_export.root\n" << endl;

  // gDebug = 2;
  RootGeoManager::Instance().Import("vecgeom_export.root");
  cout << "\033[0m" << endl;
  GeoManager::Instance().GetWorld()->PrintContent();
  cout << endl << "placed vol count: " << GeoManager::Instance().GetPlacedVolumesCount() << endl;
  cout << "registered vol count: " << GeoManager::Instance().GetRegisteredVolumesCount() << endl << endl;

  cout << "\033[0m" << endl;
  GeoManager::Instance().Clear();
}

void usage(std::string name)
{
  cerr << "Usage: " << name << " [mode]\n"
       << "Available modes:\n"
       << "\tw -> write geometry to root file\n"
       << "\tr -> read geometry from root file\n"
       << "\tno args -> executes write and read in sequence" << endl;
}

int main(int argc, char *argv[])
{
  if (argc == 1) {
    write();
    read();
    return 0;
  } else if (argc == 2) {
    switch (argv[1][0]) {
    case 'w':
      write();
      return 0;
    case 'r':
      read();
      return 0;
    default:
      usage(argv[0]);
    }

  } else {
    usage(argv[0]);
  }
  return 1;
}
