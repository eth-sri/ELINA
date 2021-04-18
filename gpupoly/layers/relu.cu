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
   \file src/layers/relu.cu
   \brief Implementation of the methods of ReLU layer.
   \author Fran&ccedil;ois Serre

   Implementation of the methods of ReLU layer.
  */
#include "relu.h"
#include "../intv.h"
#include "../filters.h"
#include <type_traits>


template<> inline Vector<double>& ReLU::activCst<double>() { return activCstD; }
template<> inline Vector<double>& ReLU::activFac<double>() { return activFacD; }
template<> inline Vector<float>& ReLU::activCst<float>() { return activCstS; }
template<> inline Vector<float>& ReLU::activFac<float>() { return activFacS; }
template<> inline const Vector<double>& ReLU::activCst<double>()const { return activCstD; }
template<> inline const Vector<double>& ReLU::activFac<double>() const { return activFacD; }
template<> inline const Vector<float>& ReLU::activCst<float>() const { return activCstS; }
template<> inline const Vector<float>& ReLU::activFac<float>() const { return activFacS; }
void ReLU::eval(Vector<double>& dest, bool sound, bool precise) {eval<double>(dest, sound,precise);}
void ReLU::eval(Vector<float>& dest, bool sound, bool precise) {eval<float>(dest, sound,precise);}
void ReLU::backSubstitute(typename AffineExpr<double>::Queue& queue, const AffineExpr<double>& expr) const {backSubstitute<double>(queue, expr);}
void ReLU::backSubstitute(typename AffineExpr<float>::Queue& queue, const AffineExpr<float>& expr) const {backSubstitute<float>(queue, expr);}


template <typename T, bool useAreaHeuristic>
__global__ void evalReLU(Intv<T>* dest, Intv<T>* activCst, Intv<T>* activFac, const Intv<T>* inputs, const size_t size)
{
	size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx >= size)
		return;
	Intv<T> cur = inputs[idx];

	if (cur.high <= 0) // we're on the lower part of the ReLU
	{
		activCst[idx] = T(0);
		activFac[idx] = T(0);
		cur = T(0);
	}
	else if (cur.low >= 0) // we're on the upper part of the ReLU
	{
		activCst[idx] = T(0);
		activFac[idx] = T(1);
	}
	else // we're in between
	{
		T lambda = useAreaHeuristic ? (cur.high > -cur.low) : 0;
		activFac[idx].low = lambda;
		activCst[idx].low = 0;
		lambda = cur.high / (cur.high - cur.low);
		activFac[idx].high = lambda;
		activCst[idx].high = Intv<T>::template mul_dr<true>(-cur.low, lambda);
		cur.low = 0;
	}
	dest[idx] = cur;
}

template <typename T>
void ReLU::eval(Vector<T>& dest, bool sound, bool precise)
{
	assert(dest.size() == this->outputSize);
	assert(dest.interval());
	if (precise)
	{
		nn.reEvaluateLayer(parent, ContainsZero<T>(), true, sound);
		nn.reEvaluateLayer(parent, ContainsZero<T>(), false, sound);
	}
	if(useAreaHeuristic)
		evalReLU<T, true> << <(this->outputSize + 255) / 256, 256 >> > (dest, activCst<T>(), activFac<T>(), nn.getConcreteBounds<T>(parent), this->outputSize);
	else
		evalReLU<T, false> << <(this->outputSize + 255) / 256, 256 >> > (dest, activCst<T>(), activFac<T>(), nn.getConcreteBounds<T>(parent), this->outputSize);
}


template <int lgBlockSize, typename TA, typename Tdest, bool upper, typename T>
static __global__ void backSubstituteReLU(Tdest* destA, T* destb, const TA* exprA, const T* exprb, const Intv<T>* modelFac, const Intv<T>* modelCst, const size_t dest_N, const size_t expr_N, const size_t n)
{
	__shared__ T red[1 << lgBlockSize];
	size_t row = blockIdx.x;
	T res = (threadIdx.x == 0 && exprb) ? exprb[row] : T(0);
	for (size_t col = threadIdx.x; col < n; col += (1 << lgBlockSize))
	{
		T in1 = Intv<T>::access_dr<false>(exprA[row * expr_N + col]);
		T coef1 = (in1 > 0 == upper) ? modelFac[col].high : modelFac[col].low;
		T off1 = (in1 > 0 == upper) ? modelCst[col].high : modelCst[col].low;
		Intv<T>::access_dr<false>(destA[row * dest_N + col]) = Intv<T>::template mul_dr<false>(in1, coef1);
		T res1 = Intv<T>::template fma_dr<upper>(off1, in1, res);
		if (std::is_same<Intv<T>, TA>::value)
		{
			T in2 = Intv<T>::access_dr<true>(exprA[row * expr_N + col]);
			T coef2 = (in2 > 0 == upper) ? modelFac[col].high : modelFac[col].low;
			T off2 = (in2 > 0 == upper) ? modelCst[col].high : modelCst[col].low;
			Intv<T>::access_dr<true>(destA[row * dest_N + col]) = Intv<T>::template mul_dr<true>(in2, coef2);
			T res2 = Intv<T>::template fma_dr<upper>(off2, in2, res);
			res = upper ? Intv<T>::max(res1, res2) : Intv<T>::min(res1, res2);
		}
		else
		{
			if (std::is_same<Intv<T>, Tdest>::value)
				Intv<T>::access_dr<true>(destA[row * dest_N + col]) = Intv<T>::template mul_dr<true>(in1, coef1);
			res = res1;
		}
	}
#pragma unroll
	for (int j = 0; j < lgBlockSize; j++)
	{
		red[threadIdx.x] = res;
		__syncthreads();
		int k = threadIdx.x + (1 << j);
		if (k < (1 << lgBlockSize))
			res = Intv<T>::template add_dr<upper>(res, red[k]);
		__syncthreads();
	}
	if (threadIdx.x == 0)
		destb[row] = res;
}

template <int lgBlockSize, bool upper, typename T>
static __global__ void backSubstituteReLUInit(
	T* destA, const size_t N, const size_t n,
	const int* rows,
	const Intv<T>* modelFac)
{
	int row = blockIdx.y;
	int col = threadIdx.x + (blockIdx.x << lgBlockSize);
	if (col < n)
	{
		int realRow = rows[row];
		destA[row * N + col] = (realRow == col) ? Intv<T>::access_dr<upper>(modelFac[col]) : T(0);
	}
}

template <typename T>
void ReLU::backSubstitute(typename AffineExpr<T>::Queue& queue, const AffineExpr<T>& expr) const
{
	if (expr.A)
	{
		std::shared_ptr<Matrix<T>> A;
		auto b = std::make_shared<Vector<T>>(expr.m, false);
		dim3 block(1024, 1, 1);
		dim3 grid(expr.m, 1, 1);
		if (expr.sound)
		{
			A = std::make_shared<Matrix<T>>(expr.m, expr.n, true);
			if (expr.A->interval())
				if (expr.up)
					backSubstituteReLU<10, Intv<T>, Intv<T>, true,T> << <grid, block >> > (*A, *b, (const Intv<T>*)*expr.A, expr.b ? (const T*)*expr.b : nullptr, activFac<T>(), activCst<T>(), A->pitch(), expr.A->pitch(), expr.n);
				else
					backSubstituteReLU<10, Intv<T>, Intv<T>, false,T> << <grid, block >> > (*A, *b, (const Intv<T>*)*expr.A, expr.b ? (const T*)*expr.b : nullptr, activFac<T>(), activCst<T>(), A->pitch(), expr.A->pitch(), expr.n);
			else
				if (expr.up)
					backSubstituteReLU<10, T, Intv<T>, true,T> << <grid, block >> > (*A, *b, (const T*)*expr.A, expr.b ? (const T*)*expr.b : nullptr, activFac<T>(), activCst<T>(), A->pitch(), expr.A->pitch(), expr.n);
				else
					backSubstituteReLU<10, T, Intv<T>, false,T> << <grid, block >> > (*A, *b, (const T*)*expr.A, expr.b ? (const T*)*expr.b : nullptr, activFac<T>(), activCst<T>(), A->pitch(), expr.A->pitch(), expr.n);
		}
		else
		{
			assert(!expr.A->interval());
			A = std::make_shared<Matrix<T>>(expr.m, expr.n, false);
			if (expr.up)
				backSubstituteReLU<10, T, T, true,T> << <grid, block >> > (*A, *b, (const T*)*expr.A, expr.b ? (const T*)*expr.b : nullptr, activFac<T>(), activCst<T>(), A->pitch(), expr.A->pitch(), expr.n);
			else
				backSubstituteReLU<10, T, T, false,T> << <grid, block >> > (*A, *b, (const T*)*expr.A, expr.b ? (const T*)*expr.b : nullptr, activFac<T>(), activCst<T>(), A->pitch(), expr.A->pitch(), expr.n);

		}
		gpuChkKer();
		queue.emplace(expr.m, expr.n, parent, expr.up, expr.rows, A, b, expr.cs,expr.sound);
	}
	else
	{
		auto A = std::make_shared<Matrix<T>>(expr.m, expr.n, false);
		auto b = std::make_shared<Vector<T>>(expr.evaluate(activCst<T>()));
		dim3 block(1024, 1, 1);
		dim3 grid((expr.n + 1023) / 1024, expr.m, 1);
		if (expr.up)
			backSubstituteReLUInit<10, true,T> << <grid, block >> > (
				*A, A->pitch(), expr.n,
				expr.rows,
				activFac<T>());
		else
			backSubstituteReLUInit<10, false,T> << <grid, block >> > (
				*A, A->pitch(), expr.n,
				expr.rows,
				activFac<T>());
		queue.emplace(expr.m, expr.n, parent, expr.up, expr.rows, A, b,ConvShape::diagonal(expr.n),expr.sound);
	}
}

