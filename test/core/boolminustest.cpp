/*
 * boolminustest.cpp
 *
 *  Created on: Aug 13, 2014
 *      Author: swenzel
 */

#include "VecGeom/volumes/TUnplacedBooleanMinusVolume.h"
#include "VecGeom/volumes/kernel/TBooleanMinusImplementation.h"
#include "VecGeom/volumes/Box.h"
#include "VecGeom/base/Transformation3D.h"
#include "VecGeom/base/Vector3D.h"
#include "VecGeom/volumes/PlacedVolume.h"
#include <iostream>

using namespace vecgeom;

// if we don't want to give all this information, we can also be very unprecise and construct something based on
// virtual functions
typedef TUnplacedBooleanMinusVolume<VPlacedVolume, VPlacedVolume> GenericSubtraction_t;

// let's try the vector interface
Real_v DistanceToOut(BoxMinusBox_t const &unplaced, Vector3D<Real_v> const &point, Vector3D<Real_v> const &dir)
{
  Real_v dist(0.);
  TBooleanMinusImplementation<translation::kIdentity, rotation::kIdentity>::DistanceToOut<BoxMinusBox_t, kVc>(
      unplaced, point, dir, kInfLength, dist);
  return dist;
}

double DistanceToOut(BoxMinusBox_t const &unplaced, Vector3D<Precision> const &point, Vector3D<Precision> const &dir)
{
  double dist(0.);
  TBooleanMinusImplementation<translation::kIdentity, rotation::kIdentity>::DistanceToOut<BoxMinusBox_t, kScalar>(
      unplaced, point, dir, kInfLength, dist);
  return dist;
}

double SafetyToOut(BoxMinusBox_t const &unplaced, Vector3D<Precision> const &point)
{
  double dist(0);
  TBooleanMinusImplementation<translation::kIdentity, rotation::kIdentity>::SafetyToOut<BoxMinusBox_t, kScalar>(
      unplaced, point, dist);
  return dist;
}
/*
double DistanceToOut2( GenericSubtraction_t const & unplaced, Vector3D<Precision> const & point, Vector3D<Precision>
const & dir)
{
    double dist(0);
    TBooleanMinusImplementation<translation::kIdentity,
rotation::kIdentity>::DistanceToOut<GenericSubtraction_t,kScalar>(
            unplaced, point, dir, kInfLength, dist);
    return dist;
}
*/

int main()
{
  UnplacedBox world(10., 10., 10.);
  UnplacedBox motherbox(5., 5., 5.);
  UnplacedBox subtractedbox(2., 2., 2);
  Transformation3D translation(-2.5, 0, 0);

  PlacedBox const *placedsubtractedbox = new PlacedBox(new LogicalVolume("", &subtractedbox), &translation);
  PlacedBox const *placedmotherbox = new PlacedBox(new LogicalVolume("", &motherbox), &Transformation3D::kIdentity);

  // now make the boolean solid
  BoxMinusBox_t complexsolid(placedmotherbox, placedsubtractedbox);

  // now calculate distance to out; here directly talking to the solid
  std::cerr << DistanceToOut(complexsolid, Vector3D<Precision>(0., 0., 0.), Vector3D<Precision>(0, 0, -1)) << "\n";
  std::cerr << DistanceToOut(complexsolid, Vector3D<Precision>(0., 0., 0.), Vector3D<Precision>(-1, 0, 0)) << "\n";
  std::cerr << SafetyToOut(complexsolid, Vector3D<Precision>(0., 0., 0.));
  std::cerr << DistanceToOut(complexsolid, Vector3D<Real_v>(0.), Vector3D<Real_v>(-1.));
}
