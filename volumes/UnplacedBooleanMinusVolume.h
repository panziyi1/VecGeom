#ifndef UNPLACEDBOOLEANMINUSVOLUME_H_
#define UNPLACEDBOOLEANMINUSVOLUME_H_

#include "base/Global.h"
#include "base/AlignedBase.h"
#include "base/Vector3D.h"
#include "volumes/UnplacedVolume.h"
#include "volumes/PlacedVolume.h"

namespace VECGEOM_NAMESPACE {

/**
 * A class representing a simple UNPLACED substraction boolean volume A-B
 * It takes two template arguments:
 * 1.: the mother (or left) volume A in unplaced form
 * 2.: the subtraction (or right) volume B in placed form;
 * the placement is with respect to the left volume
 *
 *
 *
 * will be a boolean solid where two boxes are subtracted
 * and B is only translated (not rotated) with respect to A
 *
 */
class UnplacedBooleanMinusVolume : public VUnplacedVolume, public AlignedBase {

public:
    VPlacedVolume const* fLeftVolume;
    VPlacedVolume const* fRightVolume;

public:
  // need a constructor
  UnplacedBooleanMinusVolume( VPlacedVolume const* left,
                              VPlacedVolume const* right ) :
                                   fLeftVolume(left),
                                   fRightVolume(right) {}

  virtual int memory_size() const { return sizeof(*this); }

  #ifdef VECGEOM_CUDA_INTERFACE
  virtual VUnplacedVolume* CopyToGpu() const;
  virtual VUnplacedVolume* CopyToGpu(VUnplacedVolume *const gpu_ptr) const;
  #endif


  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  Precision Capacity() const {
    // TBDONE -- need some sampling
    return 0.;
  }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  Precision SurfaceArea() const {
    // TBDONE -- need some sampling
    return 0.;
  }


  VECGEOM_CUDA_HEADER_BOTH
  void Extent( Vector3D<Precision> &, Vector3D<Precision> &) const {
     // TBDONE
  };


  VECGEOM_CUDA_HEADER_BOTH
  Vector3D<Precision> GetPointOnSurface() const {
    // TBDONE
      return Vector3D<Precision>() ;
  }


  VECGEOM_CUDA_HEADER_BOTH
  virtual void Print() const {} ;

  virtual void Print(std::ostream &os) const {};

  template <TranslationCode transCodeT, RotationCode rotCodeT>
   VECGEOM_CUDA_HEADER_DEVICE
   static VPlacedVolume* Create(LogicalVolume const *const logical_volume,
                                Transformation3D const *const transformation,
 #ifdef VECGEOM_NVCC
                                const int id,
 #endif
                                VPlacedVolume *const placement = NULL);


 private:

   VECGEOM_CUDA_HEADER_DEVICE
   virtual VPlacedVolume* SpecializedVolume(
       LogicalVolume const *const volume,
       Transformation3D const *const transformation,
       const TranslationCode trans_code, const RotationCode rot_code,
 #ifdef VECGEOM_NVCC
       const int id,
 #endif
       VPlacedVolume *const placement = NULL) const;

}; // End class

} // End global namespace



#endif /* UNPLACEDBOOLEANMINUSVOLUME_H_ */
