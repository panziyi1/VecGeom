/*
 * PolyconeHistorical.h
 *
 *  Created on: Apr 27, 2017
 *      Author: rsehgal
 */

#ifndef VOLUMES_POLYCONEHISTORICAL_H_
#define VOLUMES_POLYCONEHISTORICAL_H_

/*
#include "base/Global.h"
#include "base/AlignedBase.h"
#include "base/Vector3D.h"
#include "volumes/UnplacedVolume.h"
#include "volumes/UnplacedCone.h"
#include "base/Vector.h"
#include <vector>
#include "volumes/Wedge.h"
*/
namespace vecgeom {
inline namespace VECGEOM_IMPL_NAMESPACE {

class PolyconeHistorical {
public:
  double fHStart_angle   = 0.;
  double fHOpening_angle = 0.;
  int fHNum_z_planes     = 0;
  double *fHZ_values     = nullptr;
  double *fHRmin         = nullptr;
  double *fHRmax         = nullptr;

  VECCORE_ATT_HOST_DEVICE
  PolyconeHistorical() {}
  VECCORE_ATT_HOST_DEVICE
  PolyconeHistorical(int z_planes) : fHNum_z_planes(z_planes)
  {
    fHZ_values = new double[z_planes];
    fHRmin     = new double[z_planes];
    fHRmax     = new double[z_planes];
    for (int i = 0; i < z_planes; i++) {
      fHZ_values[i] = 0.0;
      fHRmin[i]     = 0.0;
      fHRmax[i]     = 0.0;
    }
  }
  VECCORE_ATT_HOST_DEVICE
  ~PolyconeHistorical()
  {
    delete[] fHZ_values;
    delete[] fHRmin;
    delete[] fHRmax;
  }
  VECCORE_ATT_HOST_DEVICE
  PolyconeHistorical(const PolyconeHistorical &source)
  {
    fHStart_angle   = source.fHStart_angle;
    fHOpening_angle = source.fHOpening_angle;
    fHNum_z_planes  = source.fHNum_z_planes;

    fHZ_values = new double[fHNum_z_planes];
    fHRmin     = new double[fHNum_z_planes];
    fHRmax     = new double[fHNum_z_planes];

    for (int i = 0; i < fHNum_z_planes; i++) {
      fHZ_values[i] = source.fHZ_values[i];
      fHRmin[i]     = source.fHRmin[i];
      fHRmax[i]     = source.fHRmax[i];
    }
  }
  VECCORE_ATT_HOST_DEVICE
  PolyconeHistorical &operator=(const PolyconeHistorical &right)
  {
    if (&right == this) return *this;

    fHStart_angle   = right.fHStart_angle;
    fHOpening_angle = right.fHOpening_angle;
    fHNum_z_planes  = right.fHNum_z_planes;

    delete[] fHZ_values;
    delete[] fHRmin;
    delete[] fHRmax;
    fHZ_values = new double[fHNum_z_planes];
    fHRmin     = new double[fHNum_z_planes];
    fHRmax     = new double[fHNum_z_planes];

    for (int i = 0; i < fHNum_z_planes; i++) {
      fHZ_values[i] = right.fHZ_values[i];
      fHRmin[i]     = right.fHRmin[i];
      fHRmax[i]     = right.fHRmax[i];
    }

    return *this;
  }
};
} // namespace VECGEOM_IMPL_NAMESPACE
} // namespace vecgeom

#endif /* VOLUMES_POLYCONEHISTORICAL_H_ */
