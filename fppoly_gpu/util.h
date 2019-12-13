#pragma once

#include <iostream>

#ifndef NO_CUDA

static void cuda_check_status(cudaError_t status) {
  if (status != cudaSuccess) {
    std::cerr << "error: CUDA API call : " << cudaGetErrorString(status)
              << std::endl;
    exit(1);
  }
}

static void cuda_check_last_kernel(const std::string &errstr) {
  auto status = cudaGetLastError();

  if (status != cudaSuccess) {
    std::cout << "error: CUDA kernel launch : " << errstr << " : "
              << cudaGetErrorString(status) << std::endl;
    exit(1);
  }
}

/// allocate space on GPU for n instances of type T
template <typename T> T *malloc_device(const size_t n) {
  void *p;
  auto status = cudaMalloc(&p, n * sizeof(T));
  cuda_check_status(status);

  return (T *)p;
}

/// copy n*T from host to device
template <typename T>
void copy_to_device(T *to, const T *from, const size_t n) {
  auto status = cudaMemcpy(to, from, n * sizeof(T), cudaMemcpyHostToDevice);
  cuda_check_status(status);
}

/// copy n*T from device to host
template <typename T> void copy_to_host(T *to, const T *from, const size_t n) {
  auto status = cudaMemcpy(to, from, n * sizeof(T), cudaMemcpyDeviceToHost);
  cuda_check_status(status);
}

#endif
