#include "navigation/NewSimpleNavigator.h"

namespace vecgeom {
inline namespace VECGEOM_IMPL_NAMESPACE {

template <>
NewSimpleNavigator<true> *NewSimpleNavigator<true>::fgInstance = nullptr;

template <>
NewSimpleNavigator<false> *NewSimpleNavigator<false>::fgInstance = nullptr;

#ifdef VECCORE_CUDA
VECCORE_ATT_DEVICE
VNavigator *gSimpleNavigator = nullptr;

template <>
VECCORE_ATT_DEVICE
VNavigator *NewSimpleNavigator<false>::Instance()
{
  if (gSimpleNavigator == nullptr) gSimpleNavigator = new NewSimpleNavigator();
  return gSimpleNavigator;
}
#endif

} // namespace VECGEOM_IMPL_NAMESPACE
} // namespace vecgeom