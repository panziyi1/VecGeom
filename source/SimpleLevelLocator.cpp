#include "navigation/SimpleLevelLocator.h"

namespace vecgeom {
inline namespace VECGEOM_IMPL_NAMESPACE {

template <>
TSimpleLevelLocator<true> *TSimpleLevelLocator<true>::fgInstance = nullptr;

template <>
TSimpleLevelLocator<false> *TSimpleLevelLocator<false>::fgInstance = nullptr;

} // namespace VECGEOM_IMPL_NAMESPACE
} // namespace vecgeom