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
   \file src/layers/conv2d.cu
   \brief Conv2D implementation.
   \author Fran&ccedil;ois Serre

   Implementation of the methods of the class Conv2D. Definition of this class can be found in src/layers/conv2d.h.
  */

#include "conv2d.h"
#include <algorithm>


template<typename T> void Conv2D<T>::eval(Vector<double>& dest, bool sound, bool precise) { eval<double>(dest, sound); }
template<typename T> void Conv2D<T>::eval(Vector<float>& dest, bool sound, bool precise) { eval<float>(dest, sound); }



template<typename T>
Conv2D<T>::Conv2D(
	NeuralNetwork& nn,
	const bool channels_first,
	const int filters,
	const int kernel_size_rows, const int kernel_size_cols,
	const int input_rows, const int input_cols, const int input_channels,
	const int stride_rows, const int stride_cols,
	const int padding_rows, const int padding_cols,
	const Matrix<T>& A,
	int parent) :
	NeuralNetwork::Layer(nn, ((input_rows + 2 * padding_rows - kernel_size_rows + stride_rows) / stride_rows)* ((input_cols + 2 * padding_cols - kernel_size_cols + stride_cols) / stride_cols)* filters),
	cs(channels_first, filters, kernel_size_rows, kernel_size_cols, input_rows, input_cols, input_channels, stride_rows, stride_cols, padding_rows, padding_cols),
	inputSize(input_rows* input_cols* input_channels),
	parent(parent),
	conv(A), convt(A.transpose())
{
	/*std::cout << "Conv2D:" << std::endl;
	cs.print();*/
}




template <bool channels_first, typename Tin, typename Tout, typename Tconv, typename Tarith>
__global__ void backSubstituteReLU(
	const Tin* expr, size_t expr_N,
	const Tconv* conv, size_t convN,
	Tout* dest, size_t dest_N, int dest_m, int dest_n,
	ConvShape cs,
	ConvShape prevShape,
	ConvShape newShape,
	const int* rows
)
{
	const int in_ch = blockIdx.x;
	const int row = blockIdx.y;
	const int realRow = rows[row];
	//extern __shared__ double sm[];
	//Texpr* sm = reinterpret_cast<Texpr*>(smc);
	int in_row_min = 0;
	int in_row_max = cs.input_rows;
	int in_col_min = 0;
	int in_col_max = cs.input_cols;
	if (newShape)
	{
		const int realOutputCol = channels_first ? realRow % newShape.output_cols : (realRow / newShape.filters) % newShape.output_cols;
		const int realOutputRow = channels_first ? (realRow / newShape.output_cols) % newShape.output_rows : realRow / (newShape.output_cols * newShape.filters);
		in_row_min = max(0, realOutputRow * newShape.stride_rows - newShape.padding_rows);
		in_row_max = min(cs.input_rows, realOutputRow * newShape.stride_rows - newShape.padding_rows + newShape.kernel_size_rows);
		in_col_min = max(0, realOutputCol * newShape.stride_cols - newShape.padding_cols);
		in_col_max = min(cs.input_cols, realOutputCol * newShape.stride_cols - newShape.padding_cols + newShape.kernel_size_cols);
	}
	for (int in_row = threadIdx.y + in_row_min; in_row < in_row_max; in_row += blockDim.y)
		for (int in_col = threadIdx.x + in_col_min; in_col < in_col_max; in_col += blockDim.x)
		{
			int col = in_row * cs.input_cols + in_col;

			int out_row_min = max(0, (in_row + cs.padding_rows - cs.kernel_size_rows + cs.stride_rows) / cs.stride_rows);
			int out_row_max = min(cs.output_rows, (in_row + cs.padding_rows) / cs.stride_rows + 1);
			int out_col_min = max(0, (in_col + cs.padding_cols - cs.kernel_size_cols + cs.stride_cols) / cs.stride_cols);
			int out_col_max = min(cs.output_cols, (in_col + cs.padding_cols) / cs.stride_cols + 1);

			if (prevShape) {
				if (prevShape.isDiagonal())
				{
					const int realOutputCol = channels_first ? realRow % cs.output_cols : (realRow / cs.filters) % cs.output_cols;
					const int realOutputRow = channels_first ? (realRow / cs.output_cols) % cs.output_rows : realRow / (cs.output_cols * cs.filters);
					out_row_min = max(out_row_min, realOutputRow);
					out_row_max = min(out_row_max, realOutputRow + 1);
					out_col_min = max(out_col_min, realOutputCol);
					out_col_max = min(out_col_max, realOutputCol + 1);
				}
				else
				{
					const int realOutputCol = channels_first ? realRow % prevShape.output_cols : (realRow / prevShape.filters) % prevShape.output_cols;
					const int realOutputRow = channels_first ? (realRow / prevShape.output_cols) % prevShape.output_rows : realRow / (prevShape.output_cols * prevShape.filters);
					out_row_min = max(out_row_min, realOutputRow * prevShape.stride_rows - prevShape.padding_rows);
					out_row_max = min(out_row_max, realOutputRow * prevShape.stride_rows - prevShape.padding_rows + prevShape.kernel_size_rows);
					out_col_min = max(out_col_min, realOutputCol * prevShape.stride_cols - prevShape.padding_cols);
					out_col_max = min(out_col_max, realOutputCol * prevShape.stride_cols - prevShape.padding_cols + prevShape.kernel_size_cols);
				}
			}

			Tout res(0);
				for (int out_row = out_row_min; out_row < out_row_max; out_row++)
				{
					const int delta_row = in_row + cs.padding_rows - cs.stride_rows * out_row;
					for (int out_col = out_col_min; out_col < out_col_max; out_col++)
					{
						for (int out_f = 0; out_f < cs.filters; out_f++)
						{

						const int delta_col = in_col + cs.padding_cols - cs.stride_cols * out_col;

						const Tin a = expr[channels_first ?
							row * expr_N + (out_f * cs.output_rows + out_row) * cs.output_cols + out_col :
							row * expr_N + (out_row * cs.output_cols + out_col) * cs.filters + out_f
						];
						const Tconv b = conv[channels_first ?
							(in_ch * cs.filters + out_f) * convN + delta_row * cs.kernel_size_cols + delta_col :
							(delta_row * cs.kernel_size_cols + delta_col) * convN + in_ch * cs.filters + out_f
						];
						Intv<Tarith>::fma(res, a, b, res);
					}
				}
			}


			dest[channels_first ?
				row * dest_N + in_ch * cs.input_rows * cs.input_cols + col :
				row * dest_N + col * cs.input_channels + in_ch
			] = res;
		}
}

template <bool channels_first, typename T, typename Td>
__global__ void backSubstituteReLUInit(
	const T* conv, size_t convN,
	Td* dest, size_t dest_N, 
	ConvShape cs,
	const int* rows
)
{
	const unsigned int row = blockIdx.x;
	const unsigned int realRow = rows[row];
	const int realOutputFilter = channels_first ? realRow /( cs.output_cols *cs.output_rows): realRow % cs.filters;
	const int realOutputCol = channels_first ? realRow % cs.output_cols : (realRow / cs.filters) % cs.output_cols;
	const int realOutputRow = channels_first ? (realRow / cs.output_cols) % cs.output_rows : realRow / (cs.output_cols * cs.filters);

	for (int delta_row = threadIdx.z; delta_row < cs.kernel_size_rows; delta_row += blockDim.z)
	{
		int input_row = realOutputRow * cs.stride_rows + delta_row - cs.padding_rows;
		if (input_row >= 0 && input_row < cs.input_rows)
			for (int delta_col = threadIdx.y; delta_col < cs.kernel_size_cols; delta_col += blockDim.y)
			{
				int input_col = realOutputCol * cs.stride_cols + delta_col - cs.padding_cols;
				if (input_col >= 0 && input_col < cs.input_cols)
					for (int in_ch = threadIdx.x; in_ch < cs.input_channels; in_ch += blockDim.x)
					{
						const T a = conv[(delta_row * cs.kernel_size_cols + delta_col) * convN + in_ch * cs.filters + realOutputFilter];
						dest[channels_first ?
							row * dest_N + (in_ch * cs.input_rows + input_row)*cs.input_cols+input_col :
							row * dest_N + (input_row*cs.input_cols+input_col) * cs.input_channels + in_ch
						] = a;
					}
			}
	}
}

template <typename T, typename Te>
void convBackSubstitute(typename AffineExpr<Te>::Queue& queue, const AffineExpr<Te>& expr, const Matrix<T>& conv,const Matrix<T>& convt, const ConvShape& cs, const int parent)
{
	if (!expr.A)
	{
		bool intervalOut = expr.sound && (std::is_same<T, double>() && std::is_same<Te, float>());
		auto A = std::make_shared<Matrix<Te>>(expr.m, cs.inputSize(), intervalOut);
		A->zeroFill();
		dim3 block(std::min(16, cs.input_channels), std::min(4, cs.kernel_size_cols), std::min(4, cs.kernel_size_rows));
		dim3 grid(expr.m);
		if(intervalOut)
			if (cs.channels_first)
				backSubstituteReLUInit<true, T, Intv<Te>> << <grid, block >> > (
					conv, conv.pitch(),
					*A, A->pitch(),
					cs,
					expr.rows
					);
			else
				backSubstituteReLUInit<false, T, Intv<Te>> << <grid, block >> > (
					conv, conv.pitch(),
					*A, A->pitch(),
					cs,
					expr.rows
					);
		else
		if (cs.channels_first)
			backSubstituteReLUInit<true,T,Te> << <grid, block >> > (
				conv, conv.pitch(),
				*A, A->pitch(),
				cs,
				expr.rows
				);
		else
			backSubstituteReLUInit<false,T,Te> << <grid, block >> > (
				conv, conv.pitch(),
				*A, A->pitch(),
				cs,
				expr.rows
				);
		gpuChkKer();
		queue.emplace(expr.m, cs.inputSize(), parent, expr.up, expr.rows, A, expr.b, cs, expr.sound);
		return;//*/

	}
	ConvShape ncs = expr.cs * cs;

	auto A = std::make_shared<Matrix<Te>>(expr.m, cs.inputSize(), expr.sound);
	dim3 block(std::min(cs.input_cols, 16), std::min(cs.input_rows, 16), 1);
	if (ncs)
	{
		A->zeroFill();
		block.x = std::min(ncs.kernel_size_cols, 16);
		block.y = std::min(ncs.kernel_size_rows, 16);
	}
	dim3 grid(cs.input_channels, A->m(), 1);
	size_t sm = 0;
	if (expr.sound)
	{
		if (expr.A->interval())
		{
			if (cs.channels_first)
				backSubstituteReLU<true,Intv<Te>, Intv<Te>,T,Te> << <grid, block, sm >> > (
					*expr.A, expr.A->pitch(),
					convt, convt.pitch(),
					*A, A->pitch(), A->m(), A->n(),
					cs,
					expr.cs,
					ncs,
					expr.rows
					);
			else
				backSubstituteReLU<false, Intv<Te>, Intv<Te>, T, Te> << <grid, block, sm >> > (
					*expr.A, expr.A->pitch(),
					conv, conv.pitch(),
					*A, A->pitch(), A->m(), A->n(),
					cs,
					expr.cs,
					ncs,
					expr.rows
					);
		}
		else
		{
			if (cs.channels_first)
				backSubstituteReLU<true, Te, Intv<Te>, T, Te> << <grid, block, sm >> > (
					*expr.A, expr.A->pitch(),
					convt, convt.pitch(),
					*A, A->pitch(), A->m(), A->n(),
					cs,
					expr.cs,
					ncs,
					expr.rows
					);
			else
				backSubstituteReLU<false, Te, Intv<Te>, T, Te> << <grid, block, sm >> > (
					*expr.A, expr.A->pitch(),
					conv, conv.pitch(),
					*A, A->pitch(), A->m(), A->n(),
					cs,
					expr.cs,
					ncs,
					expr.rows
					);
		}
	}
	else
	{
		assert(!expr.A->interval());
		if (cs.channels_first)
			backSubstituteReLU<true, Te, Te, T, Te> << <grid, block, sm >> > (
				*expr.A, expr.A->pitch(),
				convt, convt.pitch(),
				*A, A->pitch(), A->m(), A->n(),
				cs,
				expr.cs,
				ncs,
				expr.rows
				);
		else
			backSubstituteReLU<false, Te, Te, T, Te> << <grid, block, sm >> > (
				*expr.A, expr.A->pitch(),
				conv, conv.pitch(),
				*A, A->pitch(), A->m(), A->n(),
				cs,
				expr.cs,
				ncs,
				expr.rows
				);
	}

	gpuChkKer();
	queue.emplace(expr.m, cs.inputSize(), parent, expr.up, expr.rows, A, expr.b, ncs, expr.sound);
}

template<typename T> void Conv2D<T>::backSubstitute(typename AffineExpr<double>::Queue& queue, const AffineExpr<double>& expr) const { convBackSubstitute(queue, expr, conv,convt,cs,parent); }
template<typename T> void Conv2D<T>::backSubstitute(typename AffineExpr<float>::Queue& queue, const AffineExpr<float>& expr) const { convBackSubstitute(queue, expr, conv, convt, cs,parent); }
#ifdef STRONG_FP_SOUNDNESS
template<> void Conv2D<double>::backSubstitute(typename AffineExpr<float>::Queue& queue, const AffineExpr<float>& expr) const { convBackSubstitute(queue, expr, *convf, *convtf, cs, parent); }
#endif

template <int lgBlockSize, bool channel_first, typename Tc, typename Td>
__global__ void convRun(
	Intv<Td>* dest,
	const Intv<Td>* input,
	const Tc* conv, size_t convN,
	ConvShape cs)
{
	extern __shared__ __align__(sizeof(Intv<Td>)) unsigned char smem[];
	Intv<Td>* red = reinterpret_cast<Intv<Td>*>(smem);

	//extern __shared__ Intv<T> red[];
	//blockDim.x: #output col
	//blockDim.y: #output row
	//blockDim.z: #output filter
	int tid = ((threadIdx.z * blockDim.y) + threadIdx.y) * blockDim.x + threadIdx.x;
	Intv<Td> res = 0;
	for (int delta_row = threadIdx.z; delta_row < cs.kernel_size_rows; delta_row += blockDim.z)
	{
		int input_row = blockIdx.y * cs.stride_rows + delta_row - cs.padding_rows;
		if (input_row >= 0 && input_row < cs.input_rows)
			for (int delta_col = threadIdx.y; delta_col < cs.kernel_size_cols; delta_col += blockDim.y)
			{
				int input_col = blockIdx.x * cs.stride_cols + delta_col - cs.padding_cols;
				if (input_col >= 0 && input_col < cs.input_cols)
					for (int ch = threadIdx.x; ch < cs.input_channels; ch += blockDim.x)
					{
						const Intv<Td> inp = input[channel_first ? ch * cs.input_rows * cs.input_cols + input_row * cs.input_cols + input_col : (input_row * cs.input_cols + input_col) * cs.input_channels + ch];
						const Tc a = conv[(delta_row * cs.kernel_size_cols + delta_col) * convN + ch * gridDim.z + blockIdx.z];
						Intv<Td>::fma(res, inp, a, res);
					}
			}
	}
	int max = blockDim.x * blockDim.y * blockDim.z;
#pragma unroll
	for (int j = 0; j < lgBlockSize; j++)
	{
		red[tid] = res;
		__syncthreads();
		int k = tid + (1 << j);
		if (k < max)
			res += red[k];
		__syncthreads();
	}
	if (tid == 0)
		dest[channel_first ? (blockIdx.z * gridDim.y + blockIdx.y) * gridDim.x + blockIdx.x : (blockIdx.y * gridDim.x + blockIdx.x) * gridDim.z + blockIdx.z] = res;
}

template <typename T>
template <typename Te>
void Conv2D<T>::eval(Vector<Te>& dest, bool sound)
{
	const Vector<Te>& input = nn.getConcreteBounds<Te>(parent);
	assert(input.interval());
	assert(input.size() == inputSize);
	dest.resize(outputSize, true);
	dim3 block(std::min(16, cs.input_channels), std::min(4, cs.kernel_size_cols), std::min(4, cs.kernel_size_rows));
	dim3 grid(cs.output_cols, cs.output_rows, cs.filters);
	const int sharedSize = block.x * block.y * block.z * sizeof(Intv<T>);
	gpuChkKer();
	if (cs.channels_first)
		convRun<8, true,T,Te> << <grid, block, sharedSize >> > (
			dest,
			input,
			conv, conv.pitch(),
			cs);
	else
		convRun<8, false,T,Te> << <grid, block, sharedSize >> > (
			dest,
			input,
			conv, conv.pitch(),
			cs);
	gpuChkKer();
}



#ifdef STRONG_FP_SOUNDNESS

template<> Conv2D<double>::Conv2D(
	NeuralNetwork& nn,
	const bool channels_first,
	const int filters,
	const int kernel_size_rows, const int kernel_size_cols,
	const int input_rows, const int input_cols, const int input_channels,
	const int stride_rows, const int stride_cols,
	const int padding_rows, const int padding_cols,
	const Matrix<double>& A,
	int parent) :
	NeuralNetwork::Layer(nn, ((input_rows + 2 * padding_rows - kernel_size_rows + stride_rows) / stride_rows)* ((input_cols + 2 * padding_cols - kernel_size_cols + stride_cols) / stride_cols)* filters),
	cs(channels_first, filters, kernel_size_rows, kernel_size_cols, input_rows, input_cols, input_channels, stride_rows, stride_cols, padding_rows, padding_cols),
	inputSize(input_rows* input_cols* input_channels),
	parent(parent),

	convf(std::make_shared<const Matrix<float>>(A, false)),
	convtf(std::make_shared<const Matrix<float>>(A.transpose(), false)),

	conv(A), convt(A.transpose())

{
	/*std::cout << "Conv2D:" << std::endl;
	cs.print();*/
}




#endif



template class Conv2D<double>;
template class Conv2D<float>;