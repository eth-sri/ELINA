/*
 *
 *  This source file is part of ELINA (ETH LIbrary for Numerical Analysis).
 *  ELINA is Copyright © 2018 Department of Computer Science, ETH Zurich
 *  This software is distributed under GNU Lesser General Public License
 * Version 3.0. For more information, see the ELINA project website at:
 *  http://elina.ethz.ch
 *
 *  THE SOFTWARE IS PROVIDED "AS-IS" WITHOUT ANY WARRANTY OF ANY KIND, EITHER
 *  EXPRESS, IMPLIED OR STATUTORY, INCLUDING BUT NOT LIMITED TO ANY WARRANTY
 *  THAT THE SOFTWARE WILL CONFORM TO SPECIFICATIONS OR BE ERROR-FREE AND ANY
 *  IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE,
 *  TITLE, OR NON-INFRINGEMENT.  IN NO EVENT SHALL ETH ZURICH BE LIABLE FOR ANY
 *  DAMAGES, INCLUDING BUT NOT LIMITED TO DIRECT, INDIRECT,
 *  SPECIAL OR CONSEQUENTIAL DAMAGES, ARISING OUT OF, RESULTING FROM, OR IN
 *  ANY WAY CONNECTED WITH THIS SOFTWARE (WHETHER OR NOT BASED UPON WARRANTY,
 *  CONTRACT, TORT OR OTHERWISE).
 *
 *  @file cuda_util.h
 *  @author Christoph Müller
 *  @brief Provides convenient C++ wrappers around CUDA functions.
 */

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
