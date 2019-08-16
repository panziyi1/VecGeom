#include "navigation/SimpleSafetyEstimator.h"

namespace vecgeom {
inline namespace VECGEOM_IMPL_NAMESPACE {

// SimpleSafetyEstimator *SimpleSafetyEstimator::fgInstance = nullptr;

#ifdef VECCORE_CUDA
VECCORE_ATT_DEVICE
SimpleSafetyEstimator *gSimpleSafetyEstimator = nullptr;

VECCORE_ATT_DEVICE
VSafetyEstimator *SimpleSafetyEstimator::Instance()
{
  if (gSimpleSafetyEstimator == nullptr) gSimpleSafetyEstimator = new SimpleSafetyEstimator();
  return gSimpleSafetyEstimator;
}
#endif

} // namespace VECGEOM_IMPL_NAMESPACE
} // namespace vecgeom
