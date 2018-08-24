#include "navigation/SimpleABBoxNavigator.h"

namespace vecgeom {
inline namespace VECGEOM_IMPL_NAMESPACE {

template <>
SimpleABBoxNavigator<true> *SimpleABBoxNavigator<true>::fgInstance = nullptr;

template <>
SimpleABBoxNavigator<false> *SimpleABBoxNavigator<false>::fgInstance = nullptr;

} // namespace VECGEOM_IMPL_NAMESPACE
} // namespace vecgeom