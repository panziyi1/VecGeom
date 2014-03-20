/**
 * @file backend.cu
 * @author Johannes de Fine Licht (johannes.definelicht@cern.ch)
 */

#include <cassert>
#include <iostream>
 
#include "backend/cuda/interface.h"

namespace vecgeom {

cudaError_t CudaCheckError(const cudaError_t err) {
  if (err != cudaSuccess) {
    std::cout << "CUDA reported error with message: \""
              << cudaGetErrorString(err) << "\"\n";
  }
  return err;
}

cudaError_t CudaCheckError() {
  return CudaCheckError(cudaGetLastError());
}

void CudaAssertError(const cudaError_t err) {
  assert(CudaCheckError(err) == cudaSuccess);
}

void CudaAssertError() {
  CudaAssertError(cudaGetLastError());
}

} // End namespace vecgeom