/// \file Iterator.h
/// \author Johannes de Fine Licht (johannes.definelicht@cern.ch)

#ifndef VECGEOM_BASE_ITERATOR_H_
#define VECGEOM_BASE_ITERATOR_H_

#include "base/Global.h"

#include <iterator>

namespace VECGEOM_NAMESPACE {

/// \brief Custom iterator class for use with container classes.
template <typename Type>
class Iterator : public std::iterator<std::forward_iterator_tag, Type> {

private:

  Type *fElement;

public:

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  Iterator(Type *const e) : fElement(e) {}

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  Iterator(Iterator const &other) : fElement(other.fElement) {}

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  Iterator& operator=(Iterator const &other) {
    fElement = other.fElement;
    return *this;
  }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  bool operator==(Iterator const &other) {
    return fElement == other.fElement;
  }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  bool operator!=(Iterator const &other) {
    return fElement != other.fElement;
  }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  Iterator& operator++() {
    this->fElement++;
    return *this;
  }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  Iterator& operator++(int) {
    Iterator temp(*this);
    ++(*this);
    return temp;
  }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  Type& operator*() {
    return *fElement;
  }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  Type* operator->() {
    return fElement;
  }

};

} // End global namespace

#endif // VECGEOM_BASE_ITERATOR_H_