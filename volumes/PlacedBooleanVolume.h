#ifndef VECGEOM_VOLUMES_PLACEDTBOOLEAN_H_
#define VECGEOM_VOLUMES_PLACEDTBOOLEAN_H_

#include "base/Global.h"
#include "backend/Backend.h"
 
#include "volumes/PlacedVolume.h"
#include "volumes/UnplacedVolume.h"
#include "volumes/UnplacedBooleanVolume.h"
#ifdef VECGEOM_ROOT
#include "TGeoShape.h"
#include "TGeoVolume.h"
#include "TGeoCompositeShape.h"
#include "TGeoBoolNode.h"
#include "TGeoMatrix.h"
#include "TGeoManager.h"
#endif

#ifdef VECGEOM_USOLIDS
#include "UBox.hh"
#endif
#ifdef VECGEOM_GEANT4
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4ThreeVector.hh"
#endif


namespace VECGEOM_NAMESPACE {

class PlacedBooleanVolume : public VPlacedVolume {

    typedef UnplacedBooleanVolume UnplacedVol_t;

public:

#ifndef VECGEOM_NVCC
  PlacedBooleanVolume(char const *const label,
            LogicalVolume const *const logicalVolume,
            Transformation3D const *const transformation,
            PlacedBox const *const boundingBox)
      : VPlacedVolume(label, logicalVolume, transformation, boundingBox) {}

  PlacedBooleanVolume(LogicalVolume const *const logicalVolume,
            Transformation3D const *const transformation,
            PlacedBox const *const boundingBox)
      : PlacedBooleanVolume("", logicalVolume, transformation, boundingBox) {}
#else
  __device__
  PlacedBooleanVolume(LogicalVolume const *const logicalVolume,
            Transformation3D const *const transformation,
            PlacedBox const *const boundingBox,
            const int id)
      : VPlacedVolume(logicalVolume, transformation, boundingBox, id) {}
#endif

  VECGEOM_CUDA_HEADER_BOTH
  virtual ~PlacedBooleanVolume() {}


  VECGEOM_CUDA_HEADER_BOTH
  UnplacedVol_t const* GetUnplacedVolume() const {
    return static_cast<UnplacedVol_t const *>(
        logical_volume()->unplaced_volume());
  }


  VECGEOM_CUDA_HEADER_BOTH
  virtual void PrintType() const { };

  // CUDA specific
  virtual int memory_size() const { return sizeof(*this); }

#ifdef VECGEOM_CUDA_INTERFACE
  virtual VPlacedVolume* CopyToGpu(
      LogicalVolume const *const logical_volume,
      Transformation3D const *const transformation,
      VPlacedVolume *const gpu_ptr) const;
  virtual VPlacedVolume* CopyToGpu(
      LogicalVolume const *const logical_volume,
      Transformation3D const *const transformation) const;
#endif

  // Comparison specific

#ifndef VECGEOM_NVCC
  virtual VPlacedVolume const* ConvertToUnspecialized() const {
   return this;
  }
#ifdef VECGEOM_ROOT
 virtual TGeoShape const* ConvertToRoot() const {
      printf("Converting to ROOT\n");
      // what do we need?
      VPlacedVolume const * left = GetUnplacedVolume()->fLeftVolume;
      VPlacedVolume const * right = GetUnplacedVolume()->fRightVolume;
      Transformation3D const * leftm = left->transformation();
      Transformation3D const * rightm = right->transformation();

      TGeoShape *shape;
      if( GetUnplacedVolume()->GetOp() == kSubtraction ){
        TGeoSubtraction * node = new TGeoSubtraction(
              const_cast<TGeoShape*>(left->ConvertToRoot()),
              const_cast<TGeoShape*>(right->ConvertToRoot()),
              leftm->ConvertToTGeoMatrix(), rightm->ConvertToTGeoMatrix());
        shape = new TGeoCompositeShape("RootComposite",node);
      }
      if( GetUnplacedVolume()->GetOp() == kUnion ){
        TGeoUnion * node = new TGeoUnion(
              const_cast<TGeoShape*>(left->ConvertToRoot()),
              const_cast<TGeoShape*>(right->ConvertToRoot()),
              leftm->ConvertToTGeoMatrix(), rightm->ConvertToTGeoMatrix());
        shape = new TGeoCompositeShape("RootComposite",node);
      }
      if( GetUnplacedVolume()->GetOp() == kIntersection ){
        TGeoIntersection * node = new TGeoIntersection(
              const_cast<TGeoShape*>(left->ConvertToRoot()),
              const_cast<TGeoShape*>(right->ConvertToRoot()),
              leftm->ConvertToTGeoMatrix(), rightm->ConvertToTGeoMatrix());
        shape = new TGeoCompositeShape("RootComposite",node);
      }
      return shape;
  }
#endif
#ifdef VECGEOM_USOLIDS
  virtual ::VUSolid const* ConvertToUSolids() const {
      printf("Converting to USOLIDS\n");
      return new UBox("",10,10,10);
  }
#endif
#ifdef VECGEOM_GEANT4
  virtual G4VSolid const* ConvertToGeant4() const {
      printf("Converting to Geant4\n");
      VPlacedVolume const * left = GetUnplacedVolume()->fLeftVolume;
      VPlacedVolume const * right = GetUnplacedVolume()->fRightVolume;
      Transformation3D const * rightm = right->transformation();
      if( GetUnplacedVolume()->GetOp() == kSubtraction ){
        return new G4SubtractionSolid( GetLabel(),
              const_cast<G4VSolid*>(left->ConvertToGeant4()),
              const_cast<G4VSolid*>(right->ConvertToGeant4()),
              0, /* no rotation for the moment */
              G4ThreeVector(rightm->Translation(0),  rightm->Translation(1),  rightm->Translation(2))
            );
      }
      if( GetUnplacedVolume()->GetOp() == kUnion ){
             return new G4UnionSolid( GetLabel(),
                   const_cast<G4VSolid*>(left->ConvertToGeant4()),
                   const_cast<G4VSolid*>(right->ConvertToGeant4()),
                   0, /* no rotation for the moment */
                   G4ThreeVector(rightm->Translation(0),  rightm->Translation(1),  rightm->Translation(2))
                 );
           }
      if( GetUnplacedVolume()->GetOp() == kIntersection ){
                  return new G4IntersectionSolid( GetLabel(),
                        const_cast<G4VSolid*>(left->ConvertToGeant4()),
                        const_cast<G4VSolid*>(right->ConvertToGeant4()),
                        0, /* no rotation for the moment */
                        G4ThreeVector(rightm->Translation(0),  rightm->Translation(1),  rightm->Translation(2))
                      );
                }

  }
  #endif
#endif // VECGEOM_BENCHMARK

}; // end class declaration

} // End global namespace

#endif // VECGEOM_VOLUMES_PLACEDTBOOLEAN_H_
