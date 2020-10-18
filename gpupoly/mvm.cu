/*
 *  GPUPoly library
 *  This source file is part of ELINA (ETH LIbrary for Numerical Analysis).
 *  ELINA is Copyright © 2020 Department of Computer Science, ETH Zurich
 *  This software is distributed under GNU Lesser General Public License Version 3.0.
 *  For more information, see the ELINA project website at:
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
 */


 /*! \file src/mvm.cu
	  \brief Implementation of Matrix-Vector multiplication
	  \author Fran&ccedil;ois Serre

	  Implementation of the Matrix-Vector multiplication members of the class Matrix, defined in src/matrix.h.
  */


#include <assert.h>
#include <algorithm>
#include <cmath>
#include <type_traits>
#include <thrust/count.h>
#include <thrust/functional.h>
#include "utils.h"
#include "intv.h"
#include "vector.h"
#include "matrix.h"
#include "cutlass/cutlass.h"
#include "cutlass/gemm/device/gemm.h"

/*
template <bool up, typename T> inline __device__ double& accessScalar(T& input);
template <> inline __device__ double& accessScalar<true>(double& input) { return input; }
template <> inline __device__ double& accessScalar<false>(double& input) { return input; }
template <> inline __device__ double& accessScalar<true>(Intv& input) { return input.high; }
template <> inline __device__ double& accessScalar<false>(Intv& input) { return input.low; }
template <bool up, typename T> inline __device__ double accessScalar(const T& input);
template <> inline __device__ double accessScalar<true>(const double& input) { return input; }
template <> inline __device__ double accessScalar<false>(const double& input) { return input; }
template <> inline __device__ double accessScalar<true>(const Intv& input) { return input.high; }
template <> inline __device__ double accessScalar<false>(const Intv& input) { return input.low; }*/

template <int lgBlockSize, bool up, typename TA, typename TB, typename TC, typename T>
__global__ void mvmKernel_dr(T* dest, const TA* A, const TB* b, const TC* c, const size_t n, const size_t N)
{
	__shared__ T red[1 << lgBlockSize];
	size_t i = blockIdx.y;
	T res = (threadIdx.x == 0 && c) ? Intv<T>::access_dr<up>(c[i]) : 0;
	for (size_t j = threadIdx.x; j < n; j += (1 << lgBlockSize))
		res = Intv<T>::fma_dr<up>(A[i * N + j], b[j], res);
#pragma unroll
	for (int j = 0; j < lgBlockSize; j++)
	{
		red[threadIdx.x] = res;
		__syncthreads();
		int k = threadIdx.x + (1 << j);
		if (k < (1 << lgBlockSize))
			res = Intv<T>::add_dr<up>(res, red[k]);
		__syncthreads();
	}
	if (threadIdx.x == 0)
	{
		dest[i] = res;
	}
}
template <typename T>
template<bool up, typename Tr>
void Matrix<T>::mvm_dr(Vector<T>& dest, const Vector<Tr>& rhs,std::shared_ptr<const Vector<T>> offset) const
{
	assert(rhs.size() == n_);
	assert(!offset || offset->size() == m_);
	assert(!offset || !offset->interval());
	dest.resize(m_, false);
	dim3 block(1024, 1, 1);
	dim3 grid(1, m_, 1);
		if (interval_)
		{
			if (rhs.interval())
				mvmKernel_dr<10, up, Intv<T>, Intv<Tr>, T,T> << <grid, block >> > (dest, *this, rhs, offset ? (const T*)*offset : nullptr, n_, pitch());
			else
				mvmKernel_dr<10, up,  Intv<T>, Tr, T,T> << <grid, block >> > (dest, *this, rhs, offset ? (const T*)*offset : nullptr, n_, pitch());
		}
		else
		{
			if (rhs.interval())
				mvmKernel_dr<10, up,  T, Intv<Tr>, T,T> << <grid, block >> > (dest, *this, rhs, offset ? (const T*)*offset : nullptr, n_, pitch());
			else
				mvmKernel_dr<10, up,  T, Tr, T,T> << <grid, block >> > (dest, *this, rhs, offset ? (const T*)*offset : nullptr, n_, pitch());
		}
	gpuChkKer();
}
template void Matrix<double>::mvm_dr<false>(Vector<double>& dest, const Vector<double>& rhs, std::shared_ptr<const Vector<double>> offset) const;
template void Matrix<double>::mvm_dr<true>(Vector<double>& dest, const Vector<double>& rhs, std::shared_ptr<const Vector<double>> offset) const;
template void Matrix<float>::mvm_dr<false>(Vector<float>& dest, const Vector<float>& rhs, std::shared_ptr<const Vector<float>> offset) const;
template void Matrix<float>::mvm_dr<true>(Vector<float>& dest, const Vector<float>& rhs, std::shared_ptr<const Vector<float>> offset) const;
template void Matrix<double>::mvm_dr<false>(Vector<double>& dest, const Vector<float>& rhs, std::shared_ptr<const Vector<double>> offset) const;
template void Matrix<double>::mvm_dr<true>(Vector<double>& dest, const Vector<float>& rhs, std::shared_ptr<const Vector<double>> offset) const;
template void Matrix<float>::mvm_dr<false>(Vector<float>& dest, const Vector<double>& rhs, std::shared_ptr<const Vector<float>> offset) const;
template void Matrix<float>::mvm_dr<true>(Vector<float>& dest, const Vector<double>& rhs, std::shared_ptr<const Vector<float>> offset) const;


template <int lgBlockSize, typename TA, typename TB, typename T>
__global__ void mvmKernel(Intv<T>* dest, const TA* A, const TB* b, const size_t n, const size_t N)
{
	extern __shared__ __align__(sizeof(Intv<T>)) unsigned char smem[];
	Intv<T>* red = reinterpret_cast<Intv<T>*>(smem);
	size_t i = blockIdx.y;
	Intv<T> res(0);
	for (size_t j = threadIdx.x; j < n; j += (1 << lgBlockSize))
		Intv<T>::fma(res, b[j],A[i * N + j], res);
#pragma unroll
	for (size_t j = 0; j < lgBlockSize; j++)
	{
		red[threadIdx.x] = res;
		__syncthreads();
		size_t k = threadIdx.x + (1 << j);
		if (k < (1 << lgBlockSize))
			res += red[k];
		__syncthreads();
	}
	if (threadIdx.x == 0)
		dest[i] = res;
}

template <typename T>
template <typename Tr>
void Matrix<T>::mvm(Vector<Tr>& dest, const Vector<Tr>& rhs) const
{
	const int lgBlockSize = 8;
	const int blockSize = 1<<lgBlockSize;
	const int sharedSize = blockSize * sizeof(Intv<T>);
	assert(rhs.size() == n_);

	dest.resize(m_, true);
	dim3 block(blockSize, 1, 1);
	dim3 grid(1, m_, 1);
	
		if (interval_)
		{
			if (rhs.interval())
				mvmKernel<lgBlockSize,Intv<T>, Intv<Tr>, Tr> << <grid, block, sharedSize >> > (dest, *this,rhs, n_, N_ / 2);
			else
				mvmKernel<lgBlockSize, Intv<T>,Tr, Tr> << <grid, block, sharedSize >> > (dest, *this, rhs,  n_, N_ / 2);
		}
		else
		{
			if (rhs.interval())
				mvmKernel<lgBlockSize,T, Intv<Tr>,Tr> << <grid, block, sharedSize >> > (dest, *this, rhs, n_, N_);
			else
				mvmKernel<lgBlockSize,T,Tr, Tr> << <grid, block, sharedSize >> > (dest, *this, rhs,  n_, N_);
		}

	gpuChkKer();
}


template void Matrix<double>::mvm(Vector<double>& dest, const Vector<double>& rhs) const;
template void Matrix<float>::mvm(Vector<float>& dest, const Vector<float>& rhs) const;
template void Matrix<double>::mvm(Vector<float>& dest, const Vector<float>& rhs) const;
template void Matrix<float>::mvm(Vector<double>& dest, const Vector<double>& rhs) const;