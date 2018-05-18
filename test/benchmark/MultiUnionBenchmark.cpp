#ifndef VECGEOM_ENABLE_CUDA

#include "volumes/LogicalVolume.h"
#include "benchmarking/Benchmarker.h"
#include "management/GeoManager.h"
#include "ArgParser.h"
#include "base/Stopwatch.h"
#include <iostream>
#include "volumes/MultiUnion.h"
#include "base/RNG.h"
#include "base/Transformation3D.h"

#ifdef VECGEOM_ROOT
#include "utilities/Visualizer.h"
#include "TPolyMarker3D.h"
#include "TPolyLine3D.h"
#include "volumes/Box.h"
#endif

using namespace vecgeom;

#endif

#ifdef VECGEOM_ROOT
void AddComponentToVisualizer(MultiUnionStruct const &munion, size_t i, Visualizer &visualizer)
{
  VPlacedVolume *placedbox    = munion.fVolumes[i];
  Vector3D<double> dimensions = ((UnplacedBox *)placedbox->GetUnplacedVolume())->dimensions();
  TPolyLine3D pl(2);
  pl.SetLineColor(kBlue);
  Vector3D<double> vert[8];
  Vector3D<double> local;
  local[0] = -dimensions[0];
  local[1] = -dimensions[1];
  local[2] = -dimensions[2];
  placedbox->GetTransformation()->InverseTransform(local, vert[0]);
  local[0] = -dimensions[0];
  local[1] = dimensions[1];
  local[2] = -dimensions[2];
  placedbox->GetTransformation()->InverseTransform(local, vert[1]);
  local[0] = dimensions[0];
  local[1] = dimensions[1];
  local[2] = -dimensions[2];
  placedbox->GetTransformation()->InverseTransform(local, vert[2]);
  local[0] = dimensions[0];
  local[1] = -dimensions[1];
  local[2] = -dimensions[2];
  placedbox->GetTransformation()->InverseTransform(local, vert[3]);
  local[0] = -dimensions[0];
  local[1] = -dimensions[1];
  local[2] = dimensions[2];
  placedbox->GetTransformation()->InverseTransform(local, vert[4]);
  local[0] = -dimensions[0];
  local[1] = dimensions[1];
  local[2] = dimensions[2];
  placedbox->GetTransformation()->InverseTransform(local, vert[5]);
  local[0] = dimensions[0];
  local[1] = dimensions[1];
  local[2] = dimensions[2];
  placedbox->GetTransformation()->InverseTransform(local, vert[6]);
  local[0] = dimensions[0];
  local[1] = -dimensions[1];
  local[2] = dimensions[2];
  placedbox->GetTransformation()->InverseTransform(local, vert[7]);

  pl.SetPoint(0, vert[0].x(), vert[0].y(), vert[0].z());
  pl.SetPoint(1, vert[1].x(), vert[1].y(), vert[1].z());
  visualizer.AddLine(pl);
  pl.SetPoint(0, vert[1].x(), vert[1].y(), vert[1].z());
  pl.SetPoint(1, vert[2].x(), vert[2].y(), vert[2].z());
  visualizer.AddLine(pl);
  pl.SetPoint(0, vert[2].x(), vert[2].y(), vert[2].z());
  pl.SetPoint(1, vert[3].x(), vert[3].y(), vert[3].z());
  visualizer.AddLine(pl);
  pl.SetPoint(0, vert[3].x(), vert[3].y(), vert[3].z());
  pl.SetPoint(1, vert[0].x(), vert[0].y(), vert[0].z());
  visualizer.AddLine(pl);
  pl.SetPoint(0, vert[4].x(), vert[4].y(), vert[4].z());
  pl.SetPoint(1, vert[5].x(), vert[5].y(), vert[5].z());
  visualizer.AddLine(pl);
  pl.SetPoint(0, vert[5].x(), vert[5].y(), vert[5].z());
  pl.SetPoint(1, vert[6].x(), vert[6].y(), vert[6].z());
  visualizer.AddLine(pl);
  pl.SetPoint(0, vert[6].x(), vert[6].y(), vert[6].z());
  pl.SetPoint(1, vert[7].x(), vert[7].y(), vert[7].z());
  visualizer.AddLine(pl);
  pl.SetPoint(0, vert[7].x(), vert[7].y(), vert[7].z());
  pl.SetPoint(1, vert[4].x(), vert[4].y(), vert[4].z());
  visualizer.AddLine(pl);
  pl.SetPoint(0, vert[0].x(), vert[0].y(), vert[0].z());
  pl.SetPoint(1, vert[4].x(), vert[4].y(), vert[4].z());
  visualizer.AddLine(pl);
  pl.SetPoint(0, vert[1].x(), vert[1].y(), vert[1].z());
  pl.SetPoint(1, vert[5].x(), vert[5].y(), vert[5].z());
  visualizer.AddLine(pl);
  pl.SetPoint(0, vert[2].x(), vert[2].y(), vert[2].z());
  pl.SetPoint(1, vert[6].x(), vert[6].y(), vert[6].z());
  visualizer.AddLine(pl);
  pl.SetPoint(0, vert[3].x(), vert[3].y(), vert[3].z());
  pl.SetPoint(1, vert[7].x(), vert[7].y(), vert[7].z());
  visualizer.AddLine(pl);
}
#endif

int main(int argc, char *argv[])
{
#ifndef VECGEOM_ENABLE_CUDA
  OPTION_INT(npoints, 1000000);
  OPTION_INT(nrep, 4);
  OPTION_DOUBLE(nsolids, 100);
  constexpr double size = 10.;

  UnplacedBox worldUnplaced = UnplacedBox(2 * size, 2 * size, 2 * size);

  UnplacedMultiUnion multiunion = UnplacedMultiUnion();
  double sized                  = size * std::pow(0.5 / nsolids, 1. / 3.);
  for (size_t i = 0; i < nsolids; ++i) {
    Vector3D<double> pos(RNG::Instance().uniform(-size, size), RNG::Instance().uniform(-size, size),
                         RNG::Instance().uniform(-size, size));
    double sizernd = RNG::Instance().uniform(0.8 * sized, 1.2 * sized);
    Transformation3D *trans =
        new Transformation3D(pos.x(), pos.y(), pos.z(), RNG::Instance().uniform(-180, 180),
                             RNG::Instance().uniform(-180, 180), RNG::Instance().uniform(-180, 180));
    UnplacedBox *box = new UnplacedBox(sizernd, sizernd, sizernd);
    std::string name = "box";
    name += std::to_string(i);
    LogicalVolume *lbox = new LogicalVolume(name.c_str(), box);
    auto placed         = lbox->Place(trans);
    multiunion.AddNode(placed);
  }
  multiunion.Close();

  std::cout << "Benchmarking multi-union solid having " << nsolids << " random boxes\n";

#ifdef VECGEOM_ROOT
  Visualizer visualizer;
  // Visualize bounding box
  Vector3D<double> deltas = 0.5 * (multiunion.GetStruct().fMaxExtent - multiunion.GetStruct().fMinExtent);
  Vector3D<double> origin = 0.5 * (multiunion.GetStruct().fMaxExtent + multiunion.GetStruct().fMinExtent);
  SimpleBox box("bbox", deltas.x(), deltas.y(), deltas.z());
  visualizer.AddVolume(box, Transformation3D(origin.x(), origin.y(), origin.z()));
  for (size_t i = 0; i < multiunion.GetNumberOfSolids(); ++i) {
    AddComponentToVisualizer(multiunion.GetStruct(), i, visualizer);
  }
  TPolyMarker3D pm(npoints);
  pm.SetMarkerColor(kRed);
  pm.SetMarkerStyle(7);
  Vector3D<double> point;
  for (size_t i = 0; i < npoints; ++i) {
    point.Set(RNG::Instance().uniform(-2 * size, 2 * size), RNG::Instance().uniform(-2 * size, 2 * size),
              RNG::Instance().uniform(-2 * size, 2 * size));
    if (multiunion.Contains(point)) pm.SetNextPoint(point[0], point[1], point[2]);
  }

  visualizer.AddPoints(pm);
  visualizer.Show();
#endif

  return 0;

/*
LogicalVolume world("world", &worldUnplaced);
LogicalVolume lunion("multiunion", &multiunion);

Transformation3D placement(0, 0, 0);
const VPlacedVolume *placedMunion = world.PlaceDaughter("multi-union", &lunion, &placement);

VPlacedVolume *worldPlaced = world.Place();

GeoManager::Instance().SetWorldAndClose(worldPlaced);

Benchmarker tester(GeoManager::Instance().GetWorld());
tester.SetVerbosity(3);
tester.SetRepetitions(nrep);
tester.SetPointCount(npoints);
tester.SetToInBias(0.8);
tester.SetPoolMultiplier(1);
//  tester.RunToInBenchmark();
//  tester.RunToOutBenchmark();
return tester.RunBenchmark();
*/
#else
  return 0;
#endif
}