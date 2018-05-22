/*
 * EmbreeManager.cpp
 *
 *  Created on: May 18, 2018
 *      Author: swenzel
 */

#include "management/EmbreeManager.h"
#include "volumes/LogicalVolume.h"
#include "volumes/PlacedVolume.h"
#include "management/GeoManager.h"
#include "management/ABBoxManager.h"
#include "volumes/PlacedVolume.h"
#include "base/Stopwatch.h"

#include <embree3/rtcore.h>

// TODO: study BACKFACE_CULLING OPTION

namespace vecgeom {
inline namespace VECGEOM_IMPL_NAMESPACE {

 double g_step;
 LogicalVolume const *g_lvol;
 VPlacedVolume const *g_pvol;
 VPlacedVolume const *g_lastexited;
 Vector3D<double> const* g_pos;
 Vector3D<double> const* g_dir;
 Vector3D<float> const* g_normals;
 int g_count;
 bool* g_geomIDs;
 EmbreeManager::BoxIdDistancePair_t* g_hitlist;

void EmbreeManager::InitStructure(LogicalVolume const *lvol)
{
  auto numregisteredlvols = GeoManager::Instance().GetRegisteredVolumesCount();
  if (fStructureHolder.size() != numregisteredlvols) {
    fStructureHolder.resize(numregisteredlvols, nullptr);
  }
  if (fStructureHolder[lvol->id()] != nullptr) {
    RemoveStructure(lvol);
  }
  BuildStructure(lvol);
}

void EmbreeManager::BuildStructure(LogicalVolume const *vol)
{
  // for a logical volume we are referring to the functions that builds everything giving just bounding
  // boxes
  int nDaughters{0};
  // get the boxes (and number of boxes), must be called before the BuildStructure
  // function call since otherwise nDaughters is not guaranteed to be initialized
  auto boxes                  = ABBoxManager::Instance().GetABBoxes(vol, nDaughters);
  auto structure              = BuildStructureFromBoundingBoxes(boxes, nDaughters);
  fStructureHolder[vol->id()] = structure;
  assert((int)vol->GetDaughters().size() == nDaughters);
}

EmbreeManager::EmbreeAccelerationStructure *EmbreeManager::BuildStructureFromBoundingBoxes(
    ABBoxManager::ABBoxContainer_t abboxes, size_t numberofdaughters) const
{
  Stopwatch timer;
  timer.Start();
  // init (device) + scene
  // make the scene
  auto device = rtcNewDevice("VecGeomDevice"); // --> could be a global device??
  auto scene = rtcNewScene(device);
  rtcSetSceneFlags(scene, RTC_SCENE_FLAG_CONTEXT_FILTER_FUNCTION);
  rtcSetSceneBuildQuality(scene, RTC_BUILD_QUALITY_HIGH);

  auto structure = new EmbreeManager::EmbreeAccelerationStructure();
  structure->fDevice = device;
  structure->fScene = scene;

  structure->fNormals = new Vector3D<float>[numberofdaughters * 12];
  structure->fNumberObjects = numberofdaughters;

  for (int d = 0; d < numberofdaughters; ++d) {
    auto lower = abboxes[2 * d];
    auto upper = abboxes[2 * d + 1];
    // add each box as different geometry to the scene (will be able to obtain the geometryID and hence the object id)
    AddBoxGeometryToScene(*structure, lower, upper);
  }

  // commit the scene
  rtcCommitScene(scene);

  auto elasped = timer.Stop();
  std::cout << "EMBREE SETUP TOOK " << elasped << "s \n";
  return structure;
}

void EmbreeManager::RemoveStructure(LogicalVolume const *lvol)
{
  // FIXME: take care of memory deletion within acceleration structure
  if (fStructureHolder[lvol->id()]) delete fStructureHolder[lvol->id()];
}

// it does not have to be a box? could by any polygonal hull
void EmbreeManager::AddArbitraryBBoxToScene(EmbreeAccelerationStructure& structure, VPlacedVolume const *pvol) const
{
  if (!pvol) {
    return;
  }
  const auto transf   = pvol->GetTransformation();
  const auto unplaced = pvol->GetLogicalVolume()->GetUnplacedVolume();

  // calc extend
  Vector3D<double> lower;
  Vector3D<double> upper;
  unplaced->Extent(lower, upper);

  AddBoxGeometryToScene(structure, lower, upper, *transf);
}


void EmbreeManager::AddBoxGeometryToScene(EmbreeAccelerationStructure& structure,
                                          Vector3D<Precision> const &lower_local,
                                          Vector3D<Precision> const &upper_local, Transformation3D const &transf) const
{
  auto embreeDevice = structure.fDevice;
  auto embreeScene = structure.fScene;

  /* create a triangulated cube with 12 triangles and 8 vertices */
  unsigned int meshID;
  RTCGeometry geom_0 = rtcNewGeometry(embreeDevice, RTC_GEOMETRY_TYPE_TRIANGLE);
  rtcSetGeometryBuildQuality(geom_0, RTC_BUILD_QUALITY_HIGH);
  rtcSetGeometryTimeStepCount(geom_0, 1);
  meshID = rtcAttachGeometry(embreeScene, geom_0);

  // to get vertices in the frame of the scene:
  const auto lower = transf.InverseTransform(lower_local);
  const auto upper = transf.InverseTransform(upper_local);

  struct Vertex {
    float x, y, z, r;
  };
  struct Triangle {
    int v0, v1, v2;
  };

  //
  /* set vertices and vertex colors */
  Vertex *vertices =
      (Vertex *)rtcSetNewGeometryBuffer(geom_0, RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT3, 4 * sizeof(float), 8);
  vertices[0].x = lower.x();
  vertices[0].y = lower.y();
  vertices[0].z = lower.z();
  vertices[1].x = lower.x();
  vertices[1].y = lower.y();
  vertices[1].z = upper.z();
  vertices[2].x = lower.x();
  vertices[2].y = upper.y();
  vertices[2].z = lower.z();
  vertices[3].x = lower.x();
  vertices[3].y = upper.y();
  vertices[3].z = upper.z();
  vertices[4].x = upper.x();
  vertices[4].y = lower.y();
  vertices[4].z = lower.z();
  vertices[5].x = upper.x();
  vertices[5].y = lower.y();
  vertices[5].z = upper.z();
  vertices[6].x = upper.x();
  vertices[6].y = upper.y();
  vertices[6].z = lower.z();
  vertices[7].x = upper.x();
  vertices[7].y = upper.y();
  vertices[7].z = upper.z();

  //  /* set triangles and face colors */
  int tri = 0;
  Triangle *triangles =
      (Triangle *)rtcSetNewGeometryBuffer(geom_0, RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT3, 3 * sizeof(int), 12);
  // build up triangles by using indices to vertices

  // left side
  triangles[tri].v0 = 0;
  triangles[tri].v1 = 2;
  triangles[tri].v2 = 1;
  tri++;
  triangles[tri].v0 = 1;
  triangles[tri].v1 = 2;
  triangles[tri].v2 = 3;
  tri++;
  // right side
  triangles[tri].v0 = 4;
  triangles[tri].v1 = 5;
  triangles[tri].v2 = 6;
  tri++;
  triangles[tri].v0 = 5;
  triangles[tri].v1 = 7;
  triangles[tri].v2 = 6;
  tri++;
  // bottom side
  triangles[tri].v0 = 0;
  triangles[tri].v1 = 1;
  triangles[tri].v2 = 4;
  tri++;
  triangles[tri].v0 = 1;
  triangles[tri].v1 = 5;
  triangles[tri].v2 = 4;
  tri++;
  // top side
  triangles[tri].v0 = 2;
  triangles[tri].v1 = 6;
  triangles[tri].v2 = 3;
  tri++;
  triangles[tri].v0 = 3;
  triangles[tri].v1 = 6;
  triangles[tri].v2 = 7;
  tri++;
  // front side
  triangles[tri].v0 = 0;
  triangles[tri].v1 = 4;
  triangles[tri].v2 = 2;
  tri++;
  triangles[tri].v0 = 2;
  triangles[tri].v1 = 4;
  triangles[tri].v2 = 6;
  tri++;

  // back side
  triangles[tri].v0 = 1;
  triangles[tri].v1 = 3;
  triangles[tri].v2 = 5;
  tri++;
  triangles[tri].v0 = 3;
  triangles[tri].v1 = 7;
  triangles[tri].v2 = 5;
  tri++;

  // calculate normals to detect on which side of triangle we are
  auto calcNormal = [triangles, vertices](int i) {
    const auto &vertex0 = vertices[triangles[i].v0];
    Vector3D<float> p0(vertex0.x, vertex0.y, vertex0.z);
    const auto &vertex1 = vertices[triangles[i].v1];
    Vector3D<float> p1(vertex1.x, vertex1.y, vertex1.z);
    const auto &vertex2 = vertices[triangles[i].v2];
    Vector3D<float> p2(vertex2.x, vertex2.y, vertex2.z);
    // the normal is (p1 - p0) x (p2 - p0);
    return (p1 - p0).Cross(p2 - p0).FixZeroes();
  };

  // show all normals
  for (int i = 0; i < 12; ++i) {
    // std::cerr << "normal " << i << " " << calcNormal(i) << "\n";
    structure.fNormals[meshID * 12 + i] = calcNormal(i);
  }

  rtcCommitGeometry(geom_0);
  rtcReleaseGeometry(geom_0);
}

//void EmbreeManager::AddBoxGeometryToScene(RTCDevice embreeDevice, RTCScene embreeScene, Vector3D<Precision> const &lower, Vector3D<Precision> const &upper) const
//{
//  // just convenience structs
//  struct Vertex {
//    float x, y, z, r;
//  };
//  struct Quad {
//    int v0, v1, v2, v3;
//  };
//
//  /* create a triangulated cube with 12 triangles and 8 vertices */
//  unsigned int meshID;
//  RTCGeometry geom_0 = rtcNewGeometry(embreeDevice, RTC_GEOMETRY_TYPE_QUAD);
//  rtcSetGeometryBuildQuality(geom_0, RTC_BUILD_QUALITY_HIGH);
//  rtcSetGeometryTimeStepCount(geom_0, 1);
//  meshID = rtcAttachGeometry(embreeScene, geom_0);
//
//  //
//  /* set vertices and vertex colors */
//  Vertex *vertices =
//      (Vertex *)rtcSetNewGeometryBuffer(geom_0, RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT3, 4 * sizeof(float), 8);
//  vertices[0].x = lower.x();
//  vertices[0].y = lower.y();
//  vertices[0].z = lower.z();
//  vertices[1].x = lower.x();
//  vertices[1].y = lower.y();
//  vertices[1].z = upper.z();
//  vertices[2].x = lower.x();
//  vertices[2].y = upper.y();
//  vertices[2].z = lower.z();
//  vertices[3].x = lower.x();
//  vertices[3].y = upper.y();
//  vertices[3].z = upper.z();
//  vertices[4].x = upper.x();
//  vertices[4].y = lower.y();
//  vertices[4].z = lower.z();
//  vertices[5].x = upper.x();
//  vertices[5].y = lower.y();
//  vertices[5].z = upper.z();
//  vertices[6].x = upper.x();
//  vertices[6].y = upper.y();
//  vertices[6].z = lower.z();
//  vertices[7].x = upper.x();
//  vertices[7].y = upper.y();
//  vertices[7].z = upper.z();
//
//  //  /* set triangles and face colors */
//  int tri = 0;
//  Quad *quads =
//      (Quad *)rtcSetNewGeometryBuffer(geom_0, RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT4, 4 * sizeof(int), 6);
//  // build up quads by using indices to vertices
//  // left side
//  quads[tri].v0 = 0;
//  quads[tri].v1 = 2;
//  quads[tri].v2 = 1;
//  quads[tri].v3 = 3;
//  tri++;
//  // right side
//  quads[tri].v0 = 4;
//  quads[tri].v1 = 5;
//  quads[tri].v2 = 6;
//  quads[tri].v3 = 7;
//  tri++;
//  // bottom side
//  quads[tri].v0 = 0;
//  quads[tri].v1 = 1;
//  quads[tri].v2 = 4;
//  quads[tri].v3 = 5;
//  tri++;
//  // top side
//  quads[tri].v0 = 2;
//  quads[tri].v1 = 6;
//  quads[tri].v2 = 3;
//  quads[tri].v3 = 7;
//  tri++;
//  // front side
//  quads[tri].v0 = 0;
//  quads[tri].v1 = 4;
//  quads[tri].v2 = 2;
//  quads[tri].v3 = 6;
//  tri++;
//
//  // back side
//  quads[tri].v0 = 1;
//  quads[tri].v1 = 3;
//  quads[tri].v2 = 5;
//  quads[tri].v3 = 7;
//
//  rtcCommitGeometry(geom_0);
//}



} // end inline namespace
} // end global namespace


