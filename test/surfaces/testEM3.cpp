#include <iostream>

#include "test/benchmark/ArgParser.h"
#include "VecGeom/volumes/LogicalVolume.h"
#include "VecGeom/volumes/Box.h"
#include "VecGeom/management/GeoManager.h"

#include <VecGeom/surfaces/Model.h>
#include <VecGeom/surfaces/BrepHelper.h>

using namespace vecgeom;
// Forwards
void CreateVecGeomWorld(int, int);

int main(int argc, char *argv[])
{
  using BrepHelper = vgbrep::BrepHelper<vecgeom::Precision>;

  OPTION_INT(layers, 50);
  OPTION_INT(absorbers, 2.);

  CreateVecGeomWorld(layers, absorbers);

  if (!BrepHelper::Instance().CreateLocalSurfaces())
    return 1;
  if (!BrepHelper::Instance().CreateGlobalSurfacesFlatTop())
    return 2;

  // Test clearing surface data
  BrepHelper::Instance().ClearData();
}

void CreateVecGeomWorld(int NbOfLayers, int NbOfAbsorbers)
{
  const double CalorSizeYZ       = 40;
  const double GapThickness      = 2.3;
  const double AbsorberThickness = 5.7;

  const double LayerThickness = GapThickness + AbsorberThickness;
  const double CalorThickness = NbOfLayers * LayerThickness;

  const double WorldSizeX  = 1.2 * CalorThickness;
  const double WorldSizeYZ = 1.2 * CalorSizeYZ;

  auto worldSolid = new vecgeom::UnplacedBox(0.5 * WorldSizeX, 0.5 * WorldSizeYZ, 0.5 * WorldSizeYZ);
  auto worldLogic = new vecgeom::LogicalVolume("World", worldSolid);
  vecgeom::VPlacedVolume *worldPlaced = worldLogic->Place();

  //
  // Calorimeter
  //
  auto calorSolid = new vecgeom::UnplacedBox(0.5 * CalorThickness, 0.5 * CalorSizeYZ, 0.5 * CalorSizeYZ);
  auto calorLogic = new vecgeom::LogicalVolume("Calorimeter", calorSolid);
  vecgeom::Transformation3D origin;
  worldLogic->PlaceDaughter(calorLogic, &origin);

  //
  // Layers
  //
  auto layerSolid = new vecgeom::UnplacedBox(0.5 * LayerThickness, 0.5 * CalorSizeYZ, 0.5 * CalorSizeYZ);

  //
  // Absorbers
  //
  auto gapSolid = new vecgeom::UnplacedBox(0.5 * GapThickness, 0.5 * CalorSizeYZ, 0.5 * CalorSizeYZ);
  auto gapLogic = new vecgeom::LogicalVolume("Gap", gapSolid);
  vecgeom::Transformation3D gapPlacement(-0.5 * LayerThickness + 0.5 * GapThickness, 0, 0);

  auto absorberSolid = new vecgeom::UnplacedBox(0.5 * AbsorberThickness, 0.5 * CalorSizeYZ, 0.5 * CalorSizeYZ);
  auto absorberLogic = new vecgeom::LogicalVolume("Absorber", absorberSolid);
  vecgeom::Transformation3D absorberPlacement(0.5 * LayerThickness - 0.5 * AbsorberThickness, 0, 0);

  // Create a new LogicalVolume per layer, we need unique IDs for scoring.
  double xCenter = -0.5 * CalorThickness + 0.5 * LayerThickness;
  for (int i = 0; i < NbOfLayers; i++) {
    auto layerLogic = new vecgeom::LogicalVolume("Layer", layerSolid);
    vecgeom::Transformation3D placement(xCenter, 0, 0);
    calorLogic->PlaceDaughter(layerLogic, &placement);

    layerLogic->PlaceDaughter(gapLogic, &gapPlacement);
    layerLogic->PlaceDaughter(absorberLogic, &absorberPlacement);

    xCenter += LayerThickness;
  }

  vecgeom::GeoManager::Instance().SetWorldAndClose(worldPlaced);
}
