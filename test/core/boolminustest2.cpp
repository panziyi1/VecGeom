/*
 * boolminustest.cpp
 *
 *  Created on: Aug 13, 2014
 *      Author: swenzel
 */

#include "VecGeom/base/Transformation3D.h"
#include "VecGeom/base/Vector3D.h"
#include "VecGeom/volumes/Box.h"
#include "VecGeom/volumes/Tube.h"
#include "VecGeom/volumes/TBooleanMinusVolume.h"
#include "VecGeom/management/GeoManager.h"
#include "VecGeom/benchmarking/Benchmarker.h"

using namespace vecgeom;

typedef TUnplacedBooleanMinusVolume BoxMinusBox_t;

typedef TSpecializedBooleanMinusVolume<PlacedBox, PlacedBox, translation::kGeneric, rotation::kIdentity>
    SpecializedVol_t;

int main()
{
  UnplacedBox worldUnplaced(10., 10., 10.);
  LogicalVolume world("world", &worldUnplaced);

  // components for boolean solid
  UnplacedBox motherbox(5., 5., 5.);
  UnplacedBox subtractedbox(2., 2., 2);
  // translation for boolean solid right shape
  Transformation3D translation(-2.5, 0, 0);

  VPlacedVolume *worldPlaced = world.Place();
  GeoManager::Instance().SetWorld(worldPlaced);

  PlacedBox const *placedsubtractedbox = new PlacedBox(new LogicalVolume("", &subtractedbox), &translation);
  PlacedBox const *placedmotherbox = new PlacedBox(new LogicalVolume("", &motherbox), &Transformation3D::kIdentity);

  // now make the unplaced boolean solid
  BoxMinusBox_t booleansolid(placedmotherbox, placedsubtractedbox);
  LogicalVolume booleanlogical("booleanL", &booleansolid);
  // placement of boolean solid
  Transformation3D placement(5, 5, 5);

  // make a special solid by hand ( need to sort out factory later )
  SpecializedVol_t *sp = new SpecializedVol_t("booleanspecialized", &booleanlogical, &placement);
  // booleansolid.CreateSpecializedVolume()
  // world.PlaceDaughter("boolean", &booleanlogical, &placement);

  // add this boolean solid to the world
  world.PlaceDaughter(sp);

  Benchmarker tester(GeoManager::Instance().GetWorld());
  tester.SetVerbosity(3);
  tester.SetPoolMultiplier(1);
  tester.SetRepetitions(1024);
  tester.SetPointCount(1 << 10);
  tester.RunInsideBenchmark();
  tester.RunToOutBenchmark();

  return 0;
}
