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

 /*! \file src/mmm.cu
	 \brief Implementation of Matrix-Matrix multiplication
	 \author Fran&ccedil;ois Serre

	 Implementation of the Matrix-Matrix multiplication member of the class Matrix, defined in src/matrix.h.
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



  /// Interval arithmetic implementation for the cutlass library
namespace cutlass {
	namespace arch {
		/*
		/// Specialisation of Fused Multiply-Add in the case of an interval-real multiplication
		template <typename LayoutA, typename LayoutB, typename LayoutC, typename T>
		struct Mma<gemm::GemmShape<1, 1, 1>, 1, Intv<T>, LayoutA, T, LayoutB, Intv<T>, LayoutC, OpMultiplyAdd> {
			CUTLASS_DEVICE
				void operator()(Array<Intv<T>, 1>& d, Array<Intv<T>, 1> const& a, Array<T, 1> const& b, Array<Intv<T>, 1> const& c) {
				Intv<T>::fma(d[0],a[0], b[0], c[0]);
			}
		};

		/// Specialisation of Fused Multiply-Add in the case of an real-interval multiplication
		template <typename LayoutA, typename LayoutB, typename LayoutC,typename T>
		struct Mma<gemm::GemmShape<1, 1, 1>, 1, T, LayoutA, Intv<T>, LayoutB, Intv<T>, LayoutC, OpMultiplyAdd> {
			CUTLASS_DEVICE
				void operator()(Array<Intv<T>, 1>& d, Array<T, 1> const& a, Array<Intv<T>, 1> const& b, Array<Intv<T>, 1> const& c) {
				Intv<T>::fma(d[0], a[0], b[0], c[0]);
			}
		};

		/// Specialisation of Fused Multiply-Add in the case of a real-real multiplication
		template <typename LayoutA, typename LayoutB, typename LayoutC,typename T>
		struct Mma<gemm::GemmShape<1, 1, 1>, 1, T, LayoutA, T, LayoutB, Intv<T>, LayoutC, OpMultiplyAdd> {
			CUTLASS_DEVICE
				void operator()(Array<Intv<T>, 1>& d, Array<T, 1> const& a, Array<T, 1> const& b, Array<Intv<T>, 1> const& c) {
				Intv<T>::fma(d[0], a[0], b[0], c[0]);
			}
		};
		*/
		template <typename LayoutA, typename LayoutB, typename LayoutC, typename Ta, typename Tb, typename Tc>
		struct Mma<gemm::GemmShape<1, 1, 1>, 1, Ta, LayoutA, Tb, LayoutB, Intv<Tc>, LayoutC, OpMultiplyAdd> {
			CUTLASS_DEVICE
				void operator()(Array<Intv<Tc>, 1>& d, Array<Ta, 1> const& a, Array<Tb, 1> const& b, Array<Intv<Tc>, 1> const& c) {
				Intv<Tc>::fma(d[0], a[0], b[0], c[0]);
			}
		};

	}
}


/*

__global__ void mmm_int_sim_old(Intv* dest, const Intv* lhs, const double* rhs, const size_t dest_N, const size_t lhs_N, const size_t rhs_N, const size_t m, const size_t n, const size_t K)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	int j = blockIdx.y * blockDim.y + threadIdx.y;
	if (i >= m || j >= n)
		return;
	Intv res;
	for (size_t k = 0; k < K; k++)
		Intv::fma(res, lhs[i * lhs_N + k], rhs[k * rhs_N + j], res);
	dest[i * dest_N + j] = res;
}

__global__ void mmm_sim_int(double* dest, const double* lhs, const double* rhs, const size_t dest_N, const size_t lhs_N, const size_t rhs_N, const size_t K)
{
	int i = blockIdx.x;
	int j = threadIdx.x;
	double low = 0;
	double high = 0;
	for (size_t k = 0; k < K; k++)
	{
		low = fmin(__fma_rd(lhs[i * lhs_N + k], rhs[k * rhs_N + 2 * j], low), __fma_rd(lhs[i * lhs_N + k], rhs[k * rhs_N + 2 * j + 1], low));
		high = fmax(__fma_ru(lhs[i * lhs_N + k], rhs[k * rhs_N + 2 * j], high), __fma_ru(lhs[i * lhs_N + k], rhs[k * rhs_N + 2 * j + 1], high));
	}
	dest[i * dest_N + 2 * j] = low;
	dest[i * dest_N + 2 * j + 1] = high;
}

Matrix::CublasHandle Matrix::cublasHandle;

__global__ void mmm_sim_sim(double* dest, const double* lhs, const double* rhs, const size_t dest_N, const size_t lhs_N, const size_t rhs_N, const size_t K)
{
	int i = blockIdx.x;
	int j = threadIdx.x;
	double low = 0;
	double high = 0;
	for (size_t k = 0; k < K; k++)
	{
		low = __fma_rd(lhs[i * lhs_N + k], rhs[k * rhs_N + j], low);
		high = __fma_ru(lhs[i * lhs_N + k], rhs[k * rhs_N + j], high);
	}
	dest[i * dest_N + 2 * j] = low;
	dest[i * dest_N + 2 * j + 1] = high;
}

__global__ void mmm_notsound(double* dest, const double* lhs, const double* rhs, const size_t dest_N, const size_t lhs_N, const size_t rhs_N, const size_t K)
{
	int i = blockIdx.x;
	int j = threadIdx.x;
	double res = 0;
	for (size_t k = 0; k < K; k++)
	{
		res += lhs[i * lhs_N + k] * rhs[k * rhs_N + j];
	}
	dest[i * dest_N + j] = res;
}
*/

template <typename TA, typename TB, typename TD>
void cutlassGemm(
	const int m, const int n, const int k,
	const TA* dataA, const size_t pitchA,
	const TB* dataB, const size_t pitchB,
	TD* dest, const size_t pitchDest, cublasHandle_t& handle
)
{
	cutlass::gemm::device::Gemm<
		TA, cutlass::layout::RowMajor, // lhs matrix
		TB, cutlass::layout::RowMajor,  // rhs matrix
		TD, cutlass::layout::RowMajor, TD, // result matrix
		cutlass::arch::OpClassSimt,
		cutlass::arch::Sm61,
		cutlass::gemm::GemmShape<64, 32, 8>, // Threadblock
		cutlass::gemm::GemmShape<8, 32, 8>, // WarpShape
		cutlass::gemm::GemmShape<1, 1, 1> // Instruction shape
	> gemm;

	gemm({ {m, n, k},  // Gemm Problem dimensions
		{ dataA, int(pitchA) },    // Tensor-ref for source matrix A
		{ dataB,  int(pitchB) },    // Tensor-ref for source matrix B
		{ NULL,  int(pitchDest) },    // Tensor-ref for source matrix C
		{ dest,  int(pitchDest) },   // Tensor-ref for destination matrix D (may be different memory than source C matrix)
		{ 1,0 } });
	gpuChkKer();
}
template <>
void cutlassGemm<double,double,double>(
	const int m, const int n, const int k,
	const double* dataA, const size_t pitchA,
	const double* dataB, const size_t pitchB,
	double* dest, const size_t pitchDest, cublasHandle_t& handle
)
{
	double alpha = 1;
	double beta = 0;
	//C=AB => Ct=BtAt
	cublasDgemm(handle,
		CUBLAS_OP_N, CUBLAS_OP_N,
		n, m, k,
		&alpha,
		dataB, pitchB,
		dataA, pitchA,
		&beta,
		dest, pitchDest
	);
}
template <>
void cutlassGemm<float, float, float>(
	const int m, const int n, const int k,
	const float* dataA, const size_t pitchA,
	const float* dataB, const size_t pitchB,
	float* dest, const size_t pitchDest, cublasHandle_t& handle
	)
{
	float alpha = 1;
	float beta = 0;
	//C=AB => Ct=BtAt
	cublasSgemm(handle,
		CUBLAS_OP_N, CUBLAS_OP_N,
		n, m, k,
		&alpha,
		dataB, pitchB,
		dataA, pitchA,
		&beta,
		dest, pitchDest
	);
}

template <typename T>
template <typename Tr>
void Matrix<T>::mmm(const Matrix<Tr>& rhs, Matrix<T>& dest, bool sound) const
{
	assert(rhs.m() == n());
	//std::cout << "MMM lhs:" << m() << "x" << n() << (interval() ? "int" : "") << " rhs:" << rhs.m() << "x" << rhs.n() << (rhs.interval() ? "int" : "") << std::endl;
	if (!sound && !interval() && !rhs.interval())
	{
		//std::cout << "oui!" << std::endl;
		dest.reshape(m(), rhs.n(), false);
		cutlassGemm<T, Tr, T>(m(), rhs.n(), n(), *this, pitch(), rhs, rhs.pitch(), dest, dest.pitch(),cublasHandle);
		/*T alpha = 1;
		T beta = 0;
		//C=AB => Ct=BtAt
		cublasDgemm(cublasHandle,
			CUBLAS_OP_N, CUBLAS_OP_N,
			rhs.n(), m(), n(),
			&alpha,
			rhs, rhs.pitch(),
			*this, pitch(),
			&beta,
			dest, dest.pitch()
		);*/
		return;
	}
	dest.reshape(m(), rhs.n(), true);
	if (interval())
	{
		assert(!rhs.interval());
		cutlassGemm<Intv<T>, Tr, Intv<T>>(m(), rhs.n(), n(), *this, pitch(), rhs, rhs.pitch(), dest, dest.pitch(), cublasHandle);
	}
	else
	{
		if (rhs.interval())
			cutlassGemm<T, Intv<Tr>, Intv<T>>(m(), rhs.n(), n(), *this, pitch(), rhs, rhs.pitch(), dest, dest.pitch(), cublasHandle);
		else
			cutlassGemm<T, Tr, Intv<T>>(m(), rhs.n(), n(), *this, pitch(), rhs, rhs.pitch(), dest, dest.pitch(), cublasHandle);
	}
}

template void Matrix<double>::mmm<double>(const Matrix<double>& rhs, Matrix<double>& dest, bool sound) const;
template void Matrix<float>::mmm<float>(const Matrix<float>& rhs, Matrix<float>& dest, bool sound) const;
template void Matrix<double>::mmm<float>(const Matrix<float>& rhs, Matrix<double>& dest, bool sound) const;
template void Matrix<float>::mmm<double>(const Matrix<double>& rhs, Matrix<float>& dest, bool sound) const;


/*template class Matrix<double>;
template class Matrix<float>;*/