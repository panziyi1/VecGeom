#ifndef VECGEOM_VOLUMES_BOOLEANSTRUCT_H_
#define VECGEOM_VOLUMES_BOOLEANSTRUCT_H_
#include "base/Global.h"

namespace vecgeom {

// Declare types shared by cxx and cuda.
namespace BooleanOp { const int kUnion = 0, kIntersection = 1, kSubtraction = 2; }

inline namespace VECGEOM_IMPL_NAMESPACE {

/**
 * A class representing a simple UNPLACED boolean volume A-B
 * It takes two template arguments:
 * 1.: the mother (or left) volume A in unplaced form
 * 2.: the (or right) volume B in placed form, acting on A with a boolean operation;
 * the placement is with respect to the left volume
 */

typedef int BooleanOperation;

struct BooleanStruct {
  VPlacedVolume const *fLeftVolume  = nullptr;
  VPlacedVolume const *fRightVolume = nullptr;
  BooleanOperation const fOp        = -1;
  mutable double fCapacity          = -1;
  mutable double fSurfaceArea       = -1;

  VECCORE_ATT_HOST_DEVICE
  BooleanStruct(BooleanOperation op, VPlacedVolume const *left, VPlacedVolume const *right)
      : fLeftVolume(left), fRightVolume(right), fOp(op)
  {
  }
  BooleanStruct(TRootIOCtor *) {}
}; // End struct

} // namespace VECGEOM_IMPL_NAMESPACE

} // namespace vecgeom

#endif /* VECGEOM_VOLUMES_BOOLEANSTRUCT_H_ */
