/// \file Raytracer.h
/// \author Andrei Gheata (andrei.gheata@cern.ch)

#ifndef VECGEOM_RAYTRACER_H_
#define VECGEOM_RAYTRACER_H_

#include <list>
#include <VecGeom/base/Global.h>
#include <VecGeom/base/Vector3D.h>
#include <VecGeom/benchmarking/Color.h>

#ifdef VECGEOM_ENABLE_CUDA
#include <VecGeom/backend/cuda/Interface.h>
#endif

namespace vecgeom {

enum ERTmodel { kRTxray = 0, kRTspecular, kRTtransparent, kRTdiffuse };
enum ERTView { kRTVparallel = 0, kRTVperspective };

VECGEOM_DEVICE_FORWARD_DECLARE(class VPlacedVolume;);

inline namespace VECGEOM_IMPL_NAMESPACE {

class NavigationState;
class VPlacedVolume;

struct Ray_t {
  using VPlacedVolumePtr_t = VPlacedVolume const *;

  Vector3D<double> fPos;
  Vector3D<double> fDir;
  NavigationState *fCrtState  = nullptr;
  NavigationState *fNextState = nullptr;
  VPlacedVolumePtr_t fVolume  = nullptr;
  int fNcrossed               = 0;
  Color_t fColor;
  bool fDone = false;
};

/// \brief Raytracing a logical volume content using a given model
///
/// In order to run a benchmark, a world volume must be provided to the
/// raytracer.
class Raytracer {

private:
  using VPlacedVolumePtr_t = VPlacedVolume const *;

  double fScale     = 0;                 ///< Scaling from pixels to world coordinates
  double fShininess = 1.;                ///< Shininess exponent in the specular model
  Vector3D<double> fSourceDir;           ///< Light source direction
  Vector3D<double> fStart;               ///< Eye position in perspectove mode
  Vector3D<double> fScreenPos;           ///< Screen position
  Vector3D<double> fDir;                 ///< Start direction of all rays in parallel view mode
  Vector3D<double> fUp;                  ///< Up vector in the shooting rectangle plane
  Vector3D<double> fRight;               ///< Right vector in the shooting rectangle plane
  Vector3D<double> fLeftC;               ///< left-down corner of the ray shooting rectangle
  int fVerbosity      = 0;               ///< Verbosity level
  int fNrays          = 0;               ///< Number of rays left to propagate
  int fSize_px        = 1024;            ///< Image pixel size in x
  int fSize_py        = 1024;            ///< Image pixel size in y
  int fVisDepth       = 1;               ///< Visible geometry depth
  Color_t fLightColor = 0xFFFFFFFF;      ///< Light color
  Color_t fObjColor   = 0x0000FFFF;      ///< Object color
  ERTmodel fModel     = kRTxray;         ///< Selected RT model
  ERTView fView       = kRTVperspective; ///< View type

  VPlacedVolumePtr_t fWorld = nullptr;    ///< World volume
  std::list<VPlacedVolumePtr_t> fVolumes; ///< List of physical volumes (to migrate to GPU)
  NavigationState *fVPstate = nullptr;    ///< Viewpoint state
  char *fNavStates          = nullptr;    ///< Actual location in memory of the navigation states
  Ray_t *fRays              = nullptr;    ///< Data block for rays

public:
  /// \brief dummy constructor, requiring to set the camera parameters and world volume after
  Raytracer() {}

  /// \param world Mother volume containing daughters that will be benchmarked.
  /// \param screen_position position of the screen in world coordinates, rays are starting from a plane passing through
  /// it, normal vector pointing to the origin of the world reference frame \param up_vector the projection of this
  /// vector on the camera plane determines the 'up' direction of the image \param img_size_px image size on X in pixels
  /// \param img_size_py image size on Y in pixels
  Raytracer(VPlacedVolumePtr_t world, Vector3D<double> const &screen_position, Vector3D<double> const &up_vector,
            int img_size_px, int img_size_py, ERTmodel model, ERTView view);

  ~Raytracer();

  /// \brief set image size in pixels
  void SetImageSize(int img_size_px, int img_size_py);

  /// \brief set raytracing model
  void SetRTModel(ERTmodel model) { fModel = model; }

  /// \brief set raytracing model
  void SetLightSourceDir(Vector3D<double> const &dir)
  {
    if (fSourceDir.Mag2() > 0) {
      fSourceDir = dir;
      fSourceDir.Normalize();
    }
  }

  /// \brief set light color
  void SetLightColor(unsigned int col) { fLightColor = col; }

  /// \brief set object color
  void SetObjColor(unsigned int col) { fObjColor = col; }

  /// \brief set visible depth
  void SetVisDepth(int depth) { fVisDepth = depth; }

  /// \brief Entry point to propagate all rays
  void PropagateRays();

  /// \brief Initialize a ray and propagate it inside the world volume if needed
  /// \param px x-axis pixel index
  /// \param py y-axis pixel index
  void StartRay(int iray);

  /// \brief Propagate a ray for a single step until next boundary
  /// \param px x-axis pixel index
  /// \param py y-axis pixel index
  /// \return 0 or 1, if ray has stopped or has to continue
  int PropagateOneStep(int iray);

  /// \brief Propagate all rays for a single step until reaching boundaries
  /// \return Number of surviving rays
  int PropagateAllOneStep();

  /// \return World whose daughters are benchmarked.
  VPlacedVolumePtr_t GetWorld() const { return fWorld; }

  /// \param World volume containing daughters to be benchmarked.
  void SetWorld(VPlacedVolumePtr_t world);

private:
  void ApplyRTmodel(Ray_t &ray, double step);
  void CreateNavigators();
  void GenerateVolumePointers(VPlacedVolumePtr_t vol);

#ifdef VECGEOM_ENABLE_CUDA
  /*
    void RunInsideCuda(double *posX, double *posY, double *posZ, bool *contains, Inside_t *inside);
    void RunToInCuda(double *posX, double *posY, double *posZ, double *dirX, double *dirY, double
    *dirZ, double *distances, double *safeties); void RunToOutCuda(double *posX, double *posY, double
    *posZ, double *dirX, double *dirY, double *dirZ, double *distances, double *safeties);
  */
  void GetVolumePointers(std::list<cxx::DevicePtr<cuda::VPlacedVolume>> &volumesGpu);
#endif
};

} // End namespace VECGEOM_IMPL_NAMESPACE
} // End namespace vecgeom

#endif // VECGEOM_RAYTRACER_H_