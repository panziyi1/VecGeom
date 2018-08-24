#include "navigation/SimpleABBoxLevelLocator.h"

namespace vecgeom {
inline namespace VECGEOM_IMPL_NAMESPACE {

template <>
TSimpleABBoxLevelLocator<true> *TSimpleABBoxLevelLocator<true>::fgInstance = nullptr;

template <>
TSimpleABBoxLevelLocator<false> *TSimpleABBoxLevelLocator<false>::fgInstance = nullptr;

} // namespace VECGEOM_IMPL_NAMESPACE
} // namespace vecgeom