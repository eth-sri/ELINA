/*
 *  GPUPoly library
 *  This source file is part of ELINA (ETH LIbrary for Numerical Analysis).
 *  ELINA is Copyright � 2020 Department of Computer Science, ETH Zurich
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
   \file src/layers/maxpool2d.cu
   \brief MaxPool2D implementation.
   \author Fran&ccedil;ois Serre

   Implementation of the methods of the class MaxPool2D. Definition of this class can be found in src/layers/maxpool2d.h.
  */

#include "maxpool2d.h"
#include "../filters.h"
#include <limits>
#include <cmath>

template<> inline Vector<double>& MaxPool2D::modelCst<double>() { return modelCstD; }
template<> inline Vector<double>& MaxPool2D::modelFac<double>() { return modelFacD; }
template<> inline Vector<float>& MaxPool2D::modelCst<float>() { return modelCstS; }
template<> inline Vector<float>& MaxPool2D::modelFac<float>() { return modelFacS; }
template<> inline const Vector<double>& MaxPool2D::modelCst<double>()const { return modelCstD; }
template<> inline const Vector<double>& MaxPool2D::modelFac<double>() const { return modelFacD; }
template<> inline const Vector<float>& MaxPool2D::modelCst<float>() const { return modelCstS; }
template<> inline const Vector<float>& MaxPool2D::modelFac<float>() const { return modelFacS; }
void MaxPool2D::eval(Vector<double>& dest, bool sound, bool precise) { eval<double>(dest, sound,precise); }
void MaxPool2D::eval(Vector<float>& dest, bool sound, bool precise) { eval<float>(dest, sound,precise); }
void MaxPool2D::backSubstitute(typename AffineExpr<double>::Queue& queue, const AffineExpr<double>& expr) const { backSubstitute<double>(queue, expr); }
void MaxPool2D::backSubstitute(typename AffineExpr<float>::Queue& queue, const AffineExpr<float>& expr) const { backSubstitute<float>(queue, expr); }


template <int lgBlockSize, typename T>
__global__ void evalMaxPool2D(
	Intv<T>* dest,
	Intv<T>* modelCst,
	Intv<T>* modelFac,
	size_t* modelNeuron,
	const Intv<T>* input,
	ConvShape cs)
{
	//blockDim.x: #output col
	//blockDim.y: #output row
	//blockDim.z: #output channel
	int tid = ((threadIdx.z * blockDim.y) + threadIdx.y) * blockDim.x + threadIdx.x;


	// First look into the current pool for the elements with the highest bounds

	int highestUpId; // Id of the elmenet with the highest upper bound of the pool
	T highestUp = -HUGE_VAL; // value of the upper bound of that element
	T highestUpDown; // value of the lower bound of that element
	
	T highestDown = -HUGE_VAL; // value of the highest lower bound in the pool
	for (int delta_row = threadIdx.y; delta_row < cs.kernel_size_rows; delta_row += blockDim.y)
	{
		int input_row = blockIdx.y * cs.stride_rows + delta_row - cs.padding_top;
		if (input_row >= 0 && input_row < cs.input_rows)
			for (int delta_col = threadIdx.x; delta_col < cs.kernel_size_cols; delta_col += blockDim.x)
			{
				int input_col = blockIdx.x * cs.stride_cols + delta_col - cs.padding_left;
				if (input_col >= 0 && input_col < cs.input_cols)
				{
					int curId = input_row * cs.input_cols + input_col;
					const Intv<T> inp = input[blockIdx.z * cs.input_rows * cs.input_cols + curId];
					highestDown = Intv<T>::max(highestDown, inp.low);
					if (highestUp < inp.high)
					{
						highestUp = inp.high;
						highestUpId = curId;
						highestUpDown = inp.low;
					}
				}
			}
	}

	// Reduce the result throughout the pool
	__shared__ T redHU[1 << lgBlockSize];
	__shared__ int redHUI[1 << lgBlockSize];
	__shared__ T redHUD[1 << lgBlockSize];
	__shared__ T redHD[1 << lgBlockSize];
	int max = blockDim.x * blockDim.y * blockDim.z;
#pragma unroll
	for (int j = 0; j < lgBlockSize; j++)
	{
		//red[tid] = res;
		redHU[tid] = highestUp;
		redHUI[tid] = highestUpId;
		redHUD[tid] = highestUpDown;
		redHD[tid] = highestDown;
		__syncthreads();
		int k = tid + (1 << j);
		if (k < max)
		{
			highestDown = Intv<T>::max(highestDown, redHD[k]);
			if (highestUp < redHU[k])
			{
				highestUp = redHU[k];
				highestUpId = redHUI[k];
				highestUpDown = redHUD[k];
			}
		}
		__syncthreads();
	}

	if (tid == 0)
	{
		redHU[0] = highestUp;
		redHUI[0] = highestUpId;
		redHUD[0] = highestUpDown;
		redHD[0] = highestDown;
	}
	__syncthreads();
	highestUp = redHU[0];
	highestUpId = redHUI[0];
	highestUpDown = redHUD[0];
	highestDown = redHD[0];

	// now all threads of the pool know the characteristics of the best neurons.

	// We now check if the neuron with the highest upper bound is uncontested, i.e. if there are no other neuron in the pool that has a higher upper bound than its lower bound. If it is contested, we maintain the best contestant upper bound.
	bool uncontested = true;
	T bestContestant = highestUpDown;
	for (int delta_row = threadIdx.y; delta_row < cs.kernel_size_rows; delta_row += blockDim.y)
	{
		int input_row = blockIdx.y * cs.stride_rows + delta_row - cs.padding_top;
		if (input_row >= 0 && input_row < cs.input_rows)
			for (int delta_col = threadIdx.x; delta_col < cs.kernel_size_cols; delta_col += blockDim.x)
			{
				int input_col = blockIdx.x * cs.stride_cols + delta_col - cs.padding_left;
				if (input_col >= 0 && input_col < cs.input_cols)
				{
					int curId = input_row * cs.input_cols + input_col;
					const Intv<T> inp = input[blockIdx.z * cs.input_rows * cs.input_cols + curId];
					if (curId != highestUpId && bestContestant < inp.high)
					{
						uncontested = false;
						bestContestant = inp.high;
					}
				}
			}
	}
	
	// reduction is easier in this case :)
	uncontested = __syncthreads_and(uncontested);

	if (uncontested)
	{
		// write the result
		if (tid == 0)
		{
			const auto destAdr = (blockIdx.z * gridDim.y + blockIdx.y) * gridDim.x + blockIdx.x;
			dest[destAdr] = Intv<T>(highestDown, highestUp);
			modelCst[destAdr] = T(0);
			modelFac[destAdr] = T(1);
			modelNeuron[destAdr] = highestUpId;
		}
	}
	else 
	{
		// we reduce the best contestant value
#pragma unroll
		for (int j = 0; j < lgBlockSize; j++)
		{
			redHUD[tid] = bestContestant;
			__syncthreads();
			int k = tid + (1 << j);
			if (k < max)
				bestContestant = Intv<T>::max(bestContestant, redHUD[k]);
			__syncthreads();
		}

		// now thread 0 knows the upper bound of the second highest value

		// write the result
		if (tid == 0)
		{
			const auto destAdr = (blockIdx.z * gridDim.y + blockIdx.y) * gridDim.x + blockIdx.x;
			dest[destAdr] = Intv<T>(highestDown, highestUp);
			T facUp = (highestUp - bestContestant) / (highestUp - highestUpDown);
			T cstUp = Intv<T>::template fma_dr<true>(facUp,-highestUpDown,bestContestant);
			bool facDown = (highestUp - highestDown) > (highestDown - highestUpDown);
			T cstDown = facDown ? 0 : highestDown;
			modelCst[destAdr] = Intv<T>(cstDown, cstUp);
			modelFac[destAdr] = Intv<T>(facDown, facUp);
			modelNeuron[destAdr] = highestUpId;
		}

	}
}

template<typename T>
void MaxPool2D::eval(Vector<T>& dest, bool sound, bool precise)
{
	AlwaysKeep<T> filter;
	if (precise)
	{
		nn.reEvaluateLayer<T>(parent, filter, false, sound);
		nn.reEvaluateLayer<T>(parent, filter, true, sound);
	}
	const Vector<T>& input = nn.getConcreteBounds<T>(parent);
	assert(input.interval());
	//assert(input.size() == inputSize);
	dest.resize(outputSize, true);
	dim3 block(std::min(16, cs.kernel_size_cols), std::min(16, cs.kernel_size_rows));
	dim3 grid(cs.output_cols, cs.output_rows, cs.filters);
	const int sharedSize = block.x * block.y * block.z * sizeof(Intv<T>);
	gpuChkKer();
		evalMaxPool2D<8, T> << <grid, block, sharedSize >> > (
			dest,
			modelCst<T>(),
			modelFac<T>(),
			modelNeuron,
			input,
			cs);
	gpuChkKer();
}



template <typename Texpr, typename Tdest, typename T>
__global__ void MaxPoolBackPropagate(
	const Texpr* expr, size_t expr_N,
	const Intv<T>* modelCst,
	const Intv<T>* modelFac,
	const size_t* modelNeuron,
	Tdest* dest, size_t dest_N, int dest_m, int dest_n,
	ConvShape cs,
	ConvShape prevShape,
	const int* rows,
	bool upper
)
{
	const int col = threadIdx.x + blockIdx.x * blockDim.x;
	const int row = threadIdx.y + blockIdx.y * blockDim.y;

	if (col >= cs.input_rows * cs.input_cols || row >= dest_m)
		return;

	const int in_ch = blockIdx.z;
	const int in_row = col / cs.input_cols;
	const int in_col = col % cs.input_cols;

	int out_row_min = max(0, (in_row + cs.padding_top - cs.kernel_size_rows + cs.stride_rows) / cs.stride_rows);
	int out_row_max = min(cs.output_rows, (in_row + cs.padding_top) / cs.stride_rows + 1);//TODO not sure if padd_top or padding_bottom
	int out_col_min = max(0, (in_col + cs.padding_left - cs.kernel_size_cols + cs.stride_cols) / cs.stride_cols);
	int out_col_max = min(cs.output_cols, (in_col + cs.padding_left) / cs.stride_cols + 1);//TODO not sure if padd_top or padding_bottom


	if (prevShape) {
		if (prevShape.isDiagonal())
		{
			const int realRow = rows[row];
			const int realOutputCol = realRow % cs.output_cols;
			const int realOutputRow =  (realRow / cs.output_cols) % cs.output_rows;
			out_row_min = max(out_row_min, realOutputRow);
			out_row_max = min(out_row_max, realOutputRow + 1);
			out_col_min = max(out_col_min, realOutputCol );
			out_col_max = min(out_col_max, realOutputCol + 1);
		}
		else
		{
			const int realRow = rows[row];
			const int realOutputCol =  realRow % prevShape.output_cols ;
			const int realOutputRow = (realRow / prevShape.output_cols) % prevShape.output_rows;
			out_row_min = max(out_row_min, realOutputRow * prevShape.stride_rows - prevShape.padding_top);
			out_row_max = min(out_row_max, realOutputRow * prevShape.stride_rows - prevShape.padding_top + prevShape.kernel_size_rows);//TODO not sure if padd_top or padding_bottom
			out_col_min = max(out_col_min, realOutputCol * prevShape.stride_cols - prevShape.padding_left);
			out_col_max = min(out_col_max, realOutputCol * prevShape.stride_cols - prevShape.padding_left + prevShape.kernel_size_cols);//TODO not sure if padd_top or padding_bottom
		}
	}

	Tdest res;
	for (int out_row = out_row_min;
		out_row < out_row_max;
		out_row++)
	{
		//const int delta_row = in_row + cs.padding_rows - cs.stride_rows * out_row;
		for (int out_col = out_col_min;
			out_col < out_col_max;
			out_col++)
		{
			//const int delta_col = in_col + cs.padding_cols - cs.stride_cols * out_col;
			//if (out_col >= 0 && out_col < cs.output_cols && out_row >= 0 && out_row < cs.output_rows)
			{
				size_t curId = (in_ch * cs.output_rows + out_row) * cs.output_cols + out_col;
				if (modelNeuron[curId] == col)
				{
					const Texpr in = expr[row * expr_N + curId];
					const Intv<T> coef =modelFac[curId];

					T in1 = Intv<T>::access_dr<false>(in);
					T coef1 = (in1 > 0 == upper) ? coef.high : coef.low;
					Intv<T>::template access_dr<false>(res) = Intv<T>::template fma_dr<false>(in1, coef1, Intv<T>::access_dr<false>(res));
					if (std::is_same<Intv<T>, Texpr>::value)
					{
						T in2 = Intv<T>::access_dr<true>(in);
						T coef2 = (in2 > 0 == upper) ? coef.high : coef.low;
						Intv<T>::template access_dr<true>(res) = Intv<T>::template fma_dr<true>(in2, coef2, Intv<T>::access_dr<true>(res));
					}
					else if (std::is_same<Intv<T>, Tdest>::value)
					{
						Intv<T>::template access_dr<true>(res) = Intv<T>::template fma_dr<true>(in1, coef1, Intv<T>::access_dr<true>(res));
					}
				}
			}
		}
	}


	dest[row * dest_N + in_ch * cs.input_rows * cs.input_cols + col] = res;
}
template <int lgBlockSize, typename TA, bool upper, typename T>
static __global__ void MaxPoolBackSubstituteCst(T* destb, const TA* exprA, const T* exprb, const Intv<T>* modelCst, const size_t expr_N, const size_t n)
{
	__shared__ T red[1 << lgBlockSize];
	size_t row = blockIdx.x;
	T res = (threadIdx.x == 0 && exprb) ? exprb[row] : 0;
	for (size_t col = threadIdx.x; col < n; col += (1 << lgBlockSize))
	{
		T in1 = Intv<T>::access_dr<false>(exprA[row * expr_N + col]);
		T off1 = (in1 > 0 == upper) ? modelCst[col].high : modelCst[col].low;
		T res1 = Intv<T>::template fma_dr<upper>(off1, in1, res);
		if (std::is_same<Intv<T>, TA>::value)
		{
			T in2 = Intv<T>::access_dr<true>(exprA[row * expr_N + col]);
			T off2 = (in2 > 0 == upper) ? modelCst[col].high : modelCst[col].low;
			T res2 = Intv<T>::template fma_dr<upper>(off2, in2, res);
			res = upper ? Intv<T>::max(res1, res2) : Intv<T>::min(res1, res2);
		}
		else
		{
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

template <typename T>
void MaxPool2D::backSubstitute(typename AffineExpr<T>::Queue& queue, const AffineExpr<T>& expr) const {
	assert(expr.A);
	//assert(expr.b);
	std::shared_ptr<Matrix<T>> A;
	ConvShape ncs = expr.cs * cs;
	dim3 block(16, 16, 1);
	dim3 grid((cs.inputSize() / cs.input_channels + block.x - 1) / block.x, (expr.m + block.y - 1) / block.y, cs.input_channels);
	if (expr.sound)
	{
		A = std::make_shared<Matrix<T>>(expr.m, cs.inputSize(), true);
		if (expr.A->interval())
		{
				MaxPoolBackPropagate<Intv<T>,Intv<T>, T> << <grid, block >> > (
					*expr.A, expr.A->pitch(),
					modelCst<T>(),
					modelFac<T>(),
					modelNeuron,
					*A, A->pitch(), A->m(), A->n(),
					cs,
					expr.cs,
					expr.rows,
					expr.up
					);
		}
		else
		{
				MaxPoolBackPropagate<T, Intv<T>,T> << <grid, block >> > (
					*expr.A, expr.A->pitch(),
					modelCst<T>(),
					modelFac<T>(),
					modelNeuron,
					*A, A->pitch(), A->m(), A->n(),
					cs,
					expr.cs,
					expr.rows,
					expr.up
					);
		}
	}
	else
	{
		assert(!expr.A->interval());
		A = std::make_shared<Matrix<T>>(expr.m, cs.inputSize(), false);
			MaxPoolBackPropagate<T,T, T> << <grid, block >> > (
				*expr.A, expr.A->pitch(),
				modelCst<T>(),
				modelFac<T>(),
				modelNeuron,
				*A, A->pitch(), A->m(), A->n(),
				cs,
				expr.cs,
				expr.rows,
				expr.up
				);
	}
	gpuChkKer();

	auto b = std::make_shared<Vector<T>>(expr.m, false);
	block=dim3(1024, 1, 1);
	grid=dim3(expr.m, 1, 1);
	
	if (expr.A->interval())
		if (expr.up)
			MaxPoolBackSubstituteCst<10, Intv<T>, true,T> << <grid, block >> > (*b, (const Intv<T>*)*expr.A, expr.b ? (const T*)*expr.b : nullptr, modelCst<T>(), expr.A->pitch(), expr.n);
		else
			MaxPoolBackSubstituteCst<10, Intv<T>, false,T> << <grid, block >> > (*b, (const Intv<T>*)*expr.A, expr.b ? (const T*)*expr.b : nullptr, modelCst<T>(), expr.A->pitch(), expr.n);
	else
		if (expr.up)
			MaxPoolBackSubstituteCst<10, T, true,T> << <grid, block >> > (*b, (const T*)*expr.A, expr.b ? (const T*)*expr.b : nullptr, modelCst<T>(), expr.A->pitch(), expr.n);
		else
			MaxPoolBackSubstituteCst<10, T, false,T> << <grid, block >> > (*b, (const T*)*expr.A, expr.b ? (const T*)*expr.b : nullptr, modelCst<T>(), expr.A->pitch(), expr.n);
	gpuChkKer();

	queue.emplace(expr.m, cs.inputSize(), parent, expr.up, expr.rows, A, b, ncs,expr.sound);
}

MaxPool2D::MaxPool2D(NeuralNetwork& nn, int pool_rows, int pool_cols, int input_rows, int input_cols, int input_channels, int stride_rows, int stride_cols, int padding_top, int padding_left, int padding_bottom, int padding_right, int parent) :
	//Layer(input_rows / pool_rows * input_cols / pool_cols * input_channels),
	NeuralNetwork::Layer(nn,((input_rows + padding_left + padding_bottom - pool_rows + stride_rows) / stride_rows)* ((input_cols + padding_left + padding_right - pool_cols + stride_cols) / stride_cols)* input_channels),
	parent(parent),
	cs(input_channels, pool_rows, pool_cols, input_rows, input_cols, input_channels, stride_rows, stride_cols, padding_top, padding_left, padding_bottom, padding_right),
	modelCstS(outputSize, true), modelCstD(outputSize, true),
	modelFacS(outputSize, true), modelFacD(outputSize, true)
{
	gpuErrchk(cudaMalloc((void**)&modelNeuron,outputSize * sizeof(size_t)));
	/*std::cout << "MaxPool2D:" << std::endl;
	cs.print();*/
}

MaxPool2D::~MaxPool2D() {
	cudaFree(modelNeuron);
}
