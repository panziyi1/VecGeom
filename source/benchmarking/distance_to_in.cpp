/**
 * @file distance_to_in.cpp
 * @author Johannes de Fine Licht (johannes.definelicht@cern.ch)
 */

#include "UBox.hh"
#include "TGeoBBox.h"
#include "base/iterator.h"
#include "base/soa3d.h"
#include "base/stopwatch.h"
#include "base/transformation_matrix.h"
#include "benchmarking/distance_to_in.h"
#include "volumes/logical_volume.h"
#include "volumes/placed_box.h"
#include "volumes/utilities/volume_utilities.h"
#include <random>

namespace VECGEOM_NAMESPACE {

DistanceToInBenchmarker::DistanceToInBenchmarker(VPlacedVolume const *const world)
    : Benchmark(world) {
	psteps_ = (Precision *) _mm_malloc( sizeof(Precision) * n_points_, kAlignmentBoundary );
	for(int i=0;i<n_points_;++i)
	{
		psteps_[i]=kInfinity;
	}
}

DistanceToInBenchmarker::~DistanceToInBenchmarker() {
  delete point_pool_;
  delete dir_pool_;
  _mm_free(psteps_);
}

void DistanceToInBenchmarker::set_pool_multiplier(const unsigned pool_multiplier) {
  if (pool_multiplier < 1) {
    std::cerr << "Pool multiplier must be an integral number >= 1.\n";
    return;
  }
  pool_multiplier_ = pool_multiplier;
}

void DistanceToInBenchmarker::GenerateVolumePointers(VPlacedVolume const *const vol) {

  volumes_.emplace(volumes_.end(), vol);

  for (Iterator<Daughter> i = vol->daughters().begin();
       i != vol->daughters().end(); ++i) {
    GenerateVolumePointers(*i);
  }

}

BenchmarkResult DistanceToInBenchmarker::GenerateBenchmarkResult(
    const Precision elapsed, const BenchmarkType type) const {
  const BenchmarkResult benchmark = {
    .elapsed = elapsed,
    .type = type,
    .repetitions = repetitions(),
    .volumes = static_cast<unsigned>(volumes_.size()),
    .points = n_points_,
    .bias = bias_
  };
  return benchmark;
}

void DistanceToInBenchmarker::PrepareBenchmark() {

  // Allocate memory
  if (point_pool_) delete point_pool_;
  if (dir_pool_) delete dir_pool_;
  point_pool_ = new SOA3D<Precision>(n_points_*pool_multiplier_);
  dir_pool_ = new SOA3D<Precision>(n_points_*pool_multiplier_);

  // Generate pointers to representations in each geometry
  volumes_.clear();
  GenerateVolumePointers(world_);

  // Generate points not contained in any daughters and set the fraction hitting
  // a daughter to the specified bias.
  volumeutilities::FillUncontainedPoints(*world_, *point_pool_);
  volumeutilities::FillBiasedDirections(*world_, *point_pool_, bias_, *dir_pool_);

  point_pool_->setfillsize(n_points_*pool_multiplier_);
  dir_pool_->setfillsize(n_points_*pool_multiplier_);
}

void DistanceToInBenchmarker::BenchmarkAll() {

  PrepareBenchmark();

  // Allocate output memory
  Precision *const distances_specialized   = AllocateDistance();
  Precision *const distances_specializedvec = AllocateDistance();
  Precision *const distances_unspecialized = AllocateDistance();
  Precision *const distances_usolids       = AllocateDistance();
  Precision *const distances_root          = AllocateDistance();

  // Run all four benchmarks
  results_.push_back(RunSpecialized(distances_specialized));
  results_.push_back(RunSpecializedVec(distances_specializedvec));
  results_.push_back(RunUnspecialized(distances_unspecialized));
  results_.push_back(RunUSolids(distances_usolids));
  results_.push_back(RunRoot(distances_root));

  // Compare results
  unsigned mismatches = 0;
  const Precision tolerance = 1e-12;
  for (unsigned i = 0; i < n_points_; ++i) {
    const bool root_mismatch =
        abs(distances_specializedvec[i] - distances_root[i]) > tolerance &&
        !(distances_specialized[i] == kInfinity &&
          distances_root[i] == 1e30);
    const bool usolids_mismatch =
        abs(distances_specializedvec[i] - distances_usolids[i]) > tolerance &&
        !(distances_specialized[i] == kInfinity &&
          distances_usolids[i] == UUtils::kInfinity);
    if (root_mismatch || usolids_mismatch) {
      if (verbose() > 1) {
        if (!mismatches) std::cout << "VecGeom / USolids / ROOT\n";
        std::cout << distances_specialized[i]  << " / "
                  << distances_usolids[i] << " / "
                  << distances_root[i]    << std::endl;
      }
      mismatches++;
    }
  }
  if (verbose()) {
    std::cout << mismatches << " / " << n_points_
              << " mismatches detected.\n";
  }

  // Clean up memory
  FreeDistance(distances_specialized);
  FreeDistance(distances_unspecialized);
  FreeDistance(distances_usolids);
  FreeDistance(distances_root);
}

void DistanceToInBenchmarker::BenchmarkSpecialized() {
  PrepareBenchmark();
  Precision *const distances = AllocateDistance();
  results_.push_back(RunSpecialized(distances));
  FreeDistance(distances);
}

void DistanceToInBenchmarker::BenchmarkUnspecialized() {
  PrepareBenchmark();
  Precision *const distances = AllocateDistance();
  results_.push_back(RunUnspecialized(distances));
  FreeDistance(distances);
}

void DistanceToInBenchmarker::BenchmarkUSolids() {
  PrepareBenchmark();
  Precision *const distances = AllocateDistance();
  results_.push_back(RunUSolids(distances));
  FreeDistance(distances);
}


void DistanceToInBenchmarker::BenchmarkRoot() {
  PrepareBenchmark();
  Precision *const distances = AllocateDistance();
  results_.push_back(RunRoot(distances));
  FreeDistance(distances);
}

BenchmarkResult DistanceToInBenchmarker::RunSpecialized(Precision *const distances) const {
  if (verbose()) std::cout << "Running specialized benchmark...";
  Stopwatch timer;
  timer.Start();
  for (unsigned r = 0; r < repetitions(); ++r) {
    const int index = (rand() % pool_multiplier_) * n_points_;
    for (std::vector<VolumePointers>::const_iterator d = volumes_.begin();
         d != volumes_.end(); ++d) {
      for (unsigned i = 0; i < n_points_; ++i) {
        const int p = index + i;
        distances[i] = d->specialized()->DistanceToIn(
          (*point_pool_)[p], (*dir_pool_)[p]
        );
      }
    }
  }
  const Precision elapsed = timer.Stop();
  if (verbose()) std::cout << " Finished in " << elapsed << "s.\n";
  return GenerateBenchmarkResult(elapsed, kSpecialized);
}

// QUESTION: HOW DO I GET A DIFFERENT VECTOR CONTAINER FROM THE POOL EACH TIME?
BenchmarkResult DistanceToInBenchmarker::RunSpecializedVec(Precision *const distances) const {
  if (verbose()) std::cout << "Running specialized benchmark with vector interface...";
  Stopwatch timer;
  timer.Start();
  for (unsigned r = 0; r < repetitions(); ++r) {
    const int index = (rand() % pool_multiplier_) * n_points_;
    for (std::vector<VolumePointers>::const_iterator d = volumes_.begin();
         d != volumes_.end(); ++d){
    		// call vector interface
    		d->specialized()->DistanceToIn(	(*point_pool_), (*dir_pool_), psteps_, distances );
    }
  }
  const Precision elapsed = timer.Stop();
  if (verbose()) std::cout << " Finished in " << elapsed << "s.\n";
  return GenerateBenchmarkResult(elapsed, kSpecialized);
}

BenchmarkResult DistanceToInBenchmarker::RunUnspecialized(Precision *const distances) const {
  if (verbose()) std::cout << "Running unspecialized benchmark...";
  Stopwatch timer;
  timer.Start();
  for (unsigned r = 0; r < repetitions(); ++r) {
    const int index = (rand() % pool_multiplier_) * n_points_;
    for (std::vector<VolumePointers>::const_iterator d = volumes_.begin();
         d != volumes_.end(); ++d) {
      for (unsigned i = 0; i < n_points_; ++i) {
        const int p = index + i;
        distances[i] = d->unspecialized()->DistanceToIn(
          (*point_pool_)[p], (*dir_pool_)[p]
        );
      }
    }
  }
  const Precision elapsed = timer.Stop();
  if (verbose()) std::cout << " Finished in " << elapsed << "s.\n";
  return GenerateBenchmarkResult(elapsed, kUnspecialized);
}

BenchmarkResult DistanceToInBenchmarker::RunUSolids(Precision *const distances) const {
  if (verbose()) std::cout << "Running USolids benchmark...";
  Stopwatch timer;
  timer.Start();
  for (unsigned r = 0; r < repetitions(); ++r) {
    const int index = (rand() % pool_multiplier_) * n_points_;
    for (std::vector<VolumePointers>::const_iterator v = volumes_.begin();
         v != volumes_.end(); ++v) {
      TransformationMatrix const *matrix = v->unspecialized()->matrix();
      for (unsigned i = 0; i < n_points_; ++i) {
        const int p = index + i;
        const Vector3D<Precision> point =
            matrix->Transform<1, 0>((*point_pool_)[p]);
        const Vector3D<Precision> dir =
            matrix->TransformRotation<0>((*dir_pool_)[p]);
        distances[i] = v->usolids()->DistanceToIn(
          UVector3(point[0], point[1], point[2]),
          UVector3(dir[0], dir[1], dir[2])
        );
      }
    }
  }
  const Precision elapsed = timer.Stop();
  if (verbose()) std::cout << " Finished in " << elapsed << "s.\n";
  return GenerateBenchmarkResult(elapsed, kUSolids);
}

BenchmarkResult DistanceToInBenchmarker::RunRoot(Precision *const distances) const {
  if (verbose()) std::cout << "Running ROOT benchmark...";
  Stopwatch timer;
  timer.Start();
  for (unsigned r = 0; r < repetitions(); ++r) {
    const int index = (rand() % pool_multiplier_) * n_points_;
    for (std::vector<VolumePointers>::const_iterator v = volumes_.begin();
         v != volumes_.end(); ++v) {
      TransformationMatrix const *matrix = v->unspecialized()->matrix();
      for (unsigned i = 0; i < n_points_; ++i) {
        const int p = index + i;
        Vector3D<Precision> point =
            matrix->Transform<1, 0>((*point_pool_)[p]);
        Vector3D<Precision> dir =
            matrix->TransformRotation<0>((*dir_pool_)[p]);
        distances[i] = v->root()->DistFromOutside(&point[0], &dir[0]);
      }
    }
  }
  const Precision elapsed = timer.Stop();
  if (verbose()) std::cout << " Finished in " << elapsed << "s.\n";
  return GenerateBenchmarkResult(elapsed, kRoot);
}

} // End global namespace
