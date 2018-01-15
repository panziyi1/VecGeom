/// \file TessellatedHelpers.cpp
/// \author Mihaela Gheata (mihaela.gheata@cern.ch)

#include <ostream>
#include "volumes/TessellatedSection.h"

namespace vecgeom {
inline namespace VECGEOM_IMPL_NAMESPACE {

#include "volumes/TessellatedSection.h"

std::ostream &operator<<(std::ostream &os, vecgeom::TessellatedSection<double, true, false> const &ts)
{
#ifndef VECCORE_ENABLE_UMESIMD
  for (size_t i = 0; i < ts.GetNfacets(); ++i) {
    os << ts.GetFacet(i) << "\n";
  }
#endif
  return os;
}

} // End inline implementation namespace
} // End global namespace
