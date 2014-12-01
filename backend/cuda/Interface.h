/// \file cuda/Backend.h
/// \author Johannes de Fine Licht (johannes.definelicht@cern.ch)

#ifndef VECGEOM_BACKEND_CUDA_INTERFACE_H_
#define VECGEOM_BACKEND_CUDA_INTERFACE_H_

#include "base/Global.h"

#ifdef VECGEOM_CUDA

#include "driver_types.h" // Required for cudaError_t type
#include "cuda_runtime.h"

#include <type_traits>

namespace vecgeom { 

#ifdef VECGEOM_NVCC

inline namespace cuda {

template <typename DataClass, typename... ArgsTypes>
__global__
void ConstructOnGpu(DataClass *const gpu_ptr, ArgsTypes... params) {
   new (gpu_ptr) DataClass(params...);
}

template <typename DataClass, typename... ArgsTypes>
void Generic_CopyToGpu(DataClass *const gpu_ptr, ArgsTypes... params)
{
   ConstructOnGpu<<<1, 1>>>(gpu_ptr, params...);
}

} // End cuda namespace

namespace cxx {

template <typename Type>
Type* AllocateOnDevice() {
  Type *ptr;
  vecgeom::cxx::CudaAssertError(vecgeom::cxx::CudaMalloc((void**)&ptr, sizeof(Type)));
  return ptr;
}

} // End cxx namespace

#else

namespace cuda {

template <typename Type> Type* AllocateOnDevice();
template <typename DataClass, typename... ArgsTypes>
void Generic_CopyToGpu(DataClass *const gpu_ptr, ArgsTypes... params);

} // End cuda namespace

#endif

namespace cxx {

template <typename Type> class SOA3D;
template <typename Type> class AOS3D;

cudaError_t CudaCheckError(const cudaError_t err);

cudaError_t CudaCheckError();

void CudaAssertError(const cudaError_t err);

void CudaAssertError();

cudaError_t CudaMalloc(void** ptr, unsigned size);

cudaError_t CudaCopyToDevice(void* tgt, void const* src, unsigned size);

cudaError_t CudaCopyFromDevice(void* tgt, void const* src, unsigned size);

cudaError_t CudaFree(void* ptr);

template <typename Type>
Type* AllocateOnGpu(const unsigned int size) {
  Type *ptr;
  CudaAssertError(CudaMalloc((void**)&ptr, size));
  return ptr;
}

template <typename Type>
Type* AllocateOnGpu() {
  return AllocateOnGpu<Type>(sizeof(Type));
}

template <typename Type>
void FreeFromGpu(Type *const ptr) {
  CudaAssertError(CudaFree(ptr));
}

template <typename Type>
void CopyToGpu(Type const *const src, Type *const tgt, const unsigned size) {
  CudaAssertError(
    CudaCopyToDevice(tgt, src, size)
  );
}

template <typename Type>
void CopyToGpu(Type const *const src, Type *const tgt) {
  CopyToGpu<Type>(src, tgt, sizeof(Type));
}

template <typename Type>
void CopyFromGpu(Type const *const src, Type *const tgt, const unsigned size) {
  CudaAssertError(
    CudaCopyFromDevice(tgt, src, size)
  );
}

class DevicePtrBase
{
   void *fPtr;

protected:
   DevicePtrBase(const DevicePtrBase&orig) : fPtr(orig.fPtr) {}
   DevicePtrBase &operator=(const DevicePtrBase&orig) { fPtr = orig.fPtr; return *this; }

   void MemcpyToDevice(const void* what, unsigned long nbytes)
   {
      CudaAssertError(cudaMemcpy(fPtr,what,nbytes,
                                 cudaMemcpyHostToDevice) );
   }

   void MemcpyToHostAsync(void* where, unsigned long nbytes, cudaStream_t stream)
   {
      CudaAssertError(cudaMemcpyAsync(where, fPtr, nbytes, cudaMemcpyDeviceToHost, stream));
   }
   void *GetPtr() const { return fPtr; }

   void Malloc(unsigned long size) {
      vecgeom::cxx::CudaAssertError(vecgeom::cxx::CudaMalloc((void**)&fPtr, size));
   }
   void Free() {
      vecgeom::cxx::CudaAssertError(vecgeom::cxx::CudaFree((void**)&fPtr));
   }
   void Increment(long add) {
      fPtr = (char*)fPtr + add;
   }
public:
   DevicePtrBase() : fPtr(0) {}
   explicit DevicePtrBase(void *input) : fPtr(input) {}
   ~DevicePtrBase() { /* does not own content per se */ }

};

template <typename Type>
class DevicePtr : private DevicePtrBase
{
public:
   DevicePtr() = default;
   DevicePtr(const DevicePtr&) = default;
   DevicePtr &operator=(const DevicePtr&orig) = default;

   // should be taking a DevicePtr<void*>
   explicit DevicePtr(void *input) : DevicePtrBase(input) {}

   // Need to go via the explicit route accepting all conversion
   // because the regular c++ compilation
   // does not actually see the declaration for the cuda version
   // (and thus can not determine the inheritance).
   template <typename inputType>
      explicit DevicePtr(DevicePtr<inputType> const &input) : DevicePtrBase((void*)input) {}

   // Disallow conversion from const to non-const.
   DevicePtr(DevicePtr<const Type> const &input, 
      typename std::enable_if<!std::is_const<Type>::value, Type>::type * = nullptr) = delete;

#ifdef VECGEOM_NVCC
   // Allows implicit conversion from DevicePtr<Derived> to DevicePtr<Base>
   template <typename inputType,
      typename std::enable_if<std::is_base_of<Type, inputType>::value>::type* = nullptr>
      DevicePtr(DevicePtr<inputType> const &input) : DevicePtrBase(input.GetPtr()) {}

   // Disallow conversion from const to non-const.
   template <typename inputType,
      typename std::enable_if<std::is_base_of<Type, inputType>::value>::type* = nullptr>
      DevicePtr(DevicePtr<const inputType> const &input) = delete;
#endif

#ifdef VECGEOM_NVCC
   template <typename DataClass, typename... ArgsTypes>
   void Construct(ArgsTypes... params) const
   {
      ConstructOnGpu<<<1, 1>>>(*(*this), params...);
   }

   static size_t SizeOf()
   {
      return sizeof<DataClass>;
   }

#else
   template <typename... ArgsTypes>
   void Construct(ArgsTypes... params) const;

   static size_t SizeOf();
#endif
};

template <typename Type>
class DevicePtr<const Type> : private DevicePtrBase
{
public:
   DevicePtr() = default;
   DevicePtr(const DevicePtr&) = default;
   DevicePtr &operator=(const DevicePtr&orig) = default;

   // should be taking a DevicePtr<void*>
   explicit DevicePtr(void *input) : DevicePtrBase(input) {}

   // Need to go via the explicit route accepting all conversion
   // because the regular c++ compilation
   // does not actually see the declaration for the cuda version
   // (and thus can not determine the inheritance).
   template <typename inputType>
      explicit DevicePtr(DevicePtr<inputType> const &input) : DevicePtrImpl<const Type>((void*)input) {}

   // Implicit conversion from non-const to const.
   DevicePtr(DevicePtr<typename std::remove_const<Type>::type > const &input) : DevicePtrImpl<const Type>((void*)input) {}
 
#ifdef VECGEOM_NVCC
   // Allows implicit conversion from DevicePtr<Derived> to DevicePtr<Base>
   template <typename inputType,
      typename std::enable_if<std::is_base_of<Type, inputType>::value>::type* = nullptr>
      DevicePtr(DevicePtr<inputType> const &input) : DevicePtrImpl<const Type>(input.GetPtr()) {}
#endif

#ifdef VECGEOM_NVCC
   template <typename DataClass, typename... ArgsTypes>
   void Construct(ArgsTypes... params) const
   {
      ConstructOnGpu<<<1, 1>>>(*(*this), params...);
   }

   static size_t SizeOf()
   {
      return sizeof<DataClass>;
   }

#else
   template <typename... ArgsTypes>
   void Construct(ArgsTypes... params) const;

   static size_t SizeOf();
#endif
};


} // End cxx namespace

} // End vecgeom namespace

#endif // VECGEOM_CUDA

#endif // VECGEOM_BACKEND_CUDA_INTERFACE_H_
