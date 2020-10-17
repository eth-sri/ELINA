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


 /*!
   \file src/affineexpr.cu
   \brief AffineExpr implementation.
   \author Fran&ccedil;ois Serre

   Implementation of the methods of the class AffineExpr, that represents a set of affine expressions. Definition of this class can be found in src/layers/conv2d.h.
  */

#include "affineexpr.h"
template<typename T>
__global__ void makeIdMatrix(T* dest, size_t N, int outputSize, const int* annoyingList)
{
	int col = blockIdx.x * blockDim.x + threadIdx.x;
	int row = blockIdx.y;
	int output = annoyingList[row];
	if (col < outputSize)
		dest[row * N + col] = (output == col);
}

template <typename T>
std::shared_ptr<const Matrix<T>> AffineExpr<T>::getA() const
{
	if (!A)
	{
		const int blockSize = 256;
		auto res = std::make_shared<Matrix<T>>(m, n, false);
		dim3 block(blockSize, 1, 1);
		dim3 grid((n + blockSize - 1) / blockSize, m, 1);
		makeIdMatrix<T> << <grid, block >> > (*res, res->pitch(), n, rows);
		gpuErrchk(cudaPeekAtLastError());
		gpuErrchk(cudaDeviceSynchronize());
		return res;
	}
	else
		return A;
}
template <bool upper, typename TB,typename T>
__global__ void backPropConstantInit(T* dest, const size_t size, const int* annoyingList, const T* exprb, const TB* b)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx < size)
	{
		T lhs = Intv<T>::access_dr<upper>(b[annoyingList[idx]]);
		T rhs = exprb ? exprb[idx] : 0;
		T res = Intv<T>::template add_dr<upper>(lhs, rhs);
		dest[idx] = res;
	}
}

template <typename T>
template <typename Tr>
Vector<T> AffineExpr<T>::evaluate(const Vector<Tr>& bounds) const
{
	Vector<T> res(m, false);
	if (A)
		if (up)
			A->template mvm_dr<true>(res, bounds, b);
		else
			A->template mvm_dr<false>(res, bounds, b);
	else
		if (bounds.interval())
			if (up)
				backPropConstantInit<true, Intv<Tr>,T> << <(m + 255) / 256, 256 >> > (res, m, rows, b ? (const T*)*b : nullptr, bounds);
			else
				backPropConstantInit<false, Intv<Tr>,T> << <(m + 255) / 256, 256 >> > (res, m, rows, b ? (const T*)*b : nullptr, bounds);
		else
			if (up)
				backPropConstantInit<true, Tr,T> << <(m + 255) / 256, 256 >> > (res, m, rows, b ? (const T*)*b : nullptr, bounds);
			else
				backPropConstantInit<false, Tr,T> << <(m + 255) / 256, 256 >> > (res, m, rows, b ? (const T*)*b : nullptr, bounds);
	gpuChkKer();
	return res;
}

template Vector<double> AffineExpr<double>::evaluate(const Vector<double>& bounds) const;
template Vector<float> AffineExpr<float>::evaluate(const Vector<float>& bounds) const;
template Vector<double> AffineExpr<double>::evaluate(const Vector<float>& bounds) const;
template Vector<float> AffineExpr<float>::evaluate(const Vector<double>& bounds) const;

template <int lgBlockSize, bool up, typename T>
__global__ void putBack(Intv<T>* dest, const T* bound, const int* annoyingList, int size)
{
	int i = (blockIdx.x << lgBlockSize) + threadIdx.x;
	if (i < size)
	{
		int row = annoyingList?annoyingList[i]:i;
		Intv<T> cur = dest[row];
		T u = bound[i];
		if (up)
		{
			if (u < cur.high)
				cur.high = u;
		}
		else
		{
			if (u > cur.low)
				cur.low = u;
		}
		dest[row] = cur;
	}
}
template <typename T>
void AffineExpr<T>::evaluateAndUpdate(Vector<T>& dest, const Vector<T>& bounds) const
{
	assert(bounds.size() == n);
	auto res = evaluate(bounds);
	if(up)
	putBack<10, true,T> << <(m + 1023) / 1024, 1024 >> > (dest, res, rows, m);
	else
	putBack<10, false,T> << <(m + 1023) / 1024, 1024 >> > (dest, res, rows, m);
	gpuChkKer();
}

template <int lgBlockSize>
__global__ void AESelect(int* rows, const int* oldRows, const size_t n)
{
	size_t idx = (blockIdx.x << lgBlockSize) + threadIdx.x;
	if (idx < n)
	{
		rows[idx] = oldRows[rows[idx]];
	}
}

template <typename T>
void AffineExpr<T>::selectRows(size_t size, int* smallRows)
{
	A = A ? std::make_shared<const Matrix<T>>(A->template selectRows<T>(size, smallRows)) : nullptr;
	b = b ? std::make_shared<const Vector<T>>(b->select(size, smallRows)) : nullptr;
	dim3 block(1024, 1, 1);
	dim3 grid((size + 1023) / 1024, 1, 1);
	AESelect<10> << <grid, block >> > (smallRows, rows, size);
	cudaMemcpy(rows, smallRows, size * sizeof(int), cudaMemcpyDeviceToDevice);
	m = size;
}

template class AffineExpr<double>;
template class AffineExpr<float>;