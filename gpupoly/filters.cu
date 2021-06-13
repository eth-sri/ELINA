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
	\file src/filters.cu
	\brief Implementation of neuron filter functors.
	\author Fran&ccedil;ois Serre

	Implementation of the method used in neuron filter functors. Definition of these can be found in src/filters.cu.
   */

#include "filters.h"
#include <thrust/sequence.h>


template <int lgBlockSize, typename T>
__global__ void listCriticalNeuronsCZ(int* dest, int* nbCritical, const Intv<T>* v, const int* oldList, size_t m)
{
	__shared__ int red[1 << lgBlockSize];
	__shared__ int offset;
	size_t idx = (blockIdx.x << lgBlockSize) + threadIdx.x;
	bool res;
	if (idx < m)
	{
		Intv<T> cur = oldList ? v[oldList[idx]] : v[idx];
		//res = cur.high > 0 && cur.low < 0;
		res = cur.high > 0 && cur.low < 0;
	}
	else
		res = false;
	int sum = res;
#pragma unroll
	for (int j = 0; j < lgBlockSize; j++)
	{
		red[threadIdx.x] = sum;
		__syncthreads();
		int k = threadIdx.x - (1 << j);
		if (k >= 0)
			sum += red[k];
		__syncthreads();
	}
	if (threadIdx.x == (1 << lgBlockSize) - 1)
		offset = atomicAdd(nbCritical, sum) - 1;
	__syncthreads();
	if (res)
		dest[offset + sum] = idx;
}


template <typename T>
int ContainsZero<T>::listCriticalNeurons(int* dest, const Vector<T>& v, int* tmpInt, const int* oldList, const int oldNbCritical) const
{
	int an = 0;
	gpuErrchk(cudaMemcpy(tmpInt, &an, sizeof(int), cudaMemcpyHostToDevice)); //TODO: replace with cudamemset
	size_t m = oldList ? oldNbCritical : v.size();
	listCriticalNeuronsCZ<8, T> << < (m + 255) / 256, 256 >> > (dest, tmpInt, v, oldList, m);
	gpuChkKer();
	gpuErrchk(cudaMemcpy(&an, tmpInt, sizeof(int), cudaMemcpyDeviceToHost));
	//std::cout << an << "/" << m << std::endl;
	return an;
}

template <typename T>
int AlwaysKeep<T>::listCriticalNeurons(int* dest, const Vector<T>& v, int* tmpInt, const int* oldList, const int oldNbCritical) const
{
	thrust::sequence(thrust::device_pointer_cast<int>(dest), thrust::device_pointer_cast<int>(dest + v.size())); //TODO: seems that the old list is always this
	return v.size();
}

template class ContainsZero<double>;
template class AlwaysKeep<double>;
template class ContainsZero<float>;
template class AlwaysKeep<float>;