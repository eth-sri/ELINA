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

 /*! \file src/vector.cu
	 \brief Implementation of Vector.
	 \author Fran&ccedil;ois Serre
 */

#include "matrix.h"

template <typename T>
GPUMem<false> Vector<T>::gpuMem;

template <typename T>
void Vector<T>::resize(size_t n, bool interval)
{
	n_ = n;
	interval_ = interval;
	if (capacity < memSize())
	{
		if (data_)
			gpuMem.free(data_,capacity);
		capacity = memSize();
		data_=static_cast<T*>(gpuMem.alloc(capacity));
	}
}

template <typename T>
Vector<T>::Vector() :n_(0), interval_(false), capacity(0), data_(NULL) {}

template <typename T>
Vector<T>::Vector(const size_t n, const bool interval) : n_(n), interval_(interval), capacity(memSize())
{
	if (n > 0)
		data_ = static_cast<T*>(gpuMem.alloc(capacity));
	else
		data_ = NULL;
}

template <typename T>
Vector<T>::Vector(const std::vector<T>& v) : Vector(v.size(),false)//interval_(false), n_(v.size()), capacity(memSize())
{
	if(n_>0)
		gpuErrchk(cudaMemcpy(data_, v.data(), memSize(), cudaMemcpyHostToDevice));
}

template <typename T>
Vector<T>::Vector(size_t size, const T* v) : Vector(size, false)
{
	if(size>0)
	gpuErrchk(cudaMemcpy(data_, v, memSize(), cudaMemcpyHostToDevice));
}

template <typename T>
Vector<T>::Vector(const std::vector<Intv<T>>& v) : Vector(v.size(),true)// interval_(true), n_(v.size()), capacity(memSize())
{
	if (n_ > 0)
		gpuErrchk(cudaMemcpy(data_, v.data(), memSize(), cudaMemcpyHostToDevice));
}



template <typename T>
Vector<T>::Vector(Vector<T>&& other) :n_(other.n_), interval_(other.interval_), data_(other.data_), capacity(other.capacity)
{
	other.data_ = NULL;
	other.capacity = 0;
}



template <typename T>
Vector<T>& Vector<T>::operator= (const std::vector<T>& other)
{
	resize(other.size(), false);
	gpuErrchk(cudaMemcpy(data_, other.data(), memSize(), cudaMemcpyHostToDevice));
	return *this;
}

template <typename T>
Vector<T>& Vector<T>::operator= (const std::vector<Intv<T>>& other)
{
	resize(other.size(), true);
	gpuErrchk(cudaMemcpy(data_, other.data(), memSize(), cudaMemcpyHostToDevice));
	return *this;
}

template <typename T>
Vector<T>& Vector<T>::operator= (Vector<T>&& other)
{
	if (&other != this)
	{
		n_ = other.n_;
		interval_ = other.interval_;
		if (data_)
			gpuMem.free(data_,capacity);
		capacity = other.capacity;
		data_ = other.data_;
		other.data_ = NULL;
		other.capacity = 0;
	}
	return *this;
}

template <typename T>
Vector<T>::~Vector()
{
	if (data_)
		gpuMem.free(data_,capacity);
}

template <typename T>
struct isNegativeDouble
{
	__host__ __device__
		bool operator()(const T& x) const {
		return x < 0;
	}
};

template <typename T>
struct isNegativeIntv
{
	__host__ __device__
		bool operator()(const Intv<T>& x) const {
		return x.low < 0;
	}
};

template <typename T>
bool Vector<T>::isPositive() const {
	if (interval_)
	{
		isNegativeIntv<T> in;
		size_t res = thrust::count_if(beginIntv(), endIntv(), in);
		return res == 0;
	}
	else
	{
		isNegativeDouble<T> in;
		size_t res = thrust::count_if(begin(), end(), in);
		return res == 0;
	}

}

template <int lgBlockSize, typename Tin, typename Tout>
__global__ void VectorSelect(Tout* v, const int* annoyingList, const size_t n, const Tin* ov)
{
	size_t idx = (blockIdx.x << lgBlockSize) + threadIdx.x;
	if (idx < n)
	{
		size_t row = annoyingList? annoyingList[idx]:idx;
		v[idx] = ov[row];
	}
}

template <typename T>
Vector<T> Vector<T>::select(size_t size, const int* rows) const
{
	Vector<T> res(size, interval_);
	dim3 block(1024, 1, 1);
	dim3 grid((size + 1023) / 1024, 1, 1);
	if (interval_)
		VectorSelect<10, Intv<T>,Intv<T>> << <grid, block >> > (res, rows, size, *this);
	else
		VectorSelect<10, T, T> << <grid, block >> > (res, rows, size, *this);
	gpuChkKer();
	return res;
}

template <typename T>
Vector<T>::Vector(const Vector<T>& other) : Vector(other.size(), other.interval())//n_(other.n_), interval_(other.interval_), capacity(memSize())
{
	if (n_ > 0)
		gpuErrchk(cudaMemcpy(data_, other.data_, memSize(), cudaMemcpyDeviceToDevice));
}


template <typename T>
template <typename To>
Vector<T>::Vector(const Vector<To>& other) : Vector(other.size(), other.interval() || (std::is_same<T,float>() && std::is_same<To,double>()))//n_(other.n_), interval_(other.interval_), capacity(memSize())
{
	dim3 block(1024, 1, 1);
	dim3 grid((n_ + 1023) / 1024, 1, 1);
	if (other.interval())
		VectorSelect<10, Intv<To>, Intv<T>> << <grid, block >> > (*this, nullptr, n_, other);
	else
		if(interval_)
			VectorSelect<10, To, Intv<T>> << <grid, block >> > (*this, nullptr, n_, other);
		else
			VectorSelect<10, To, T> << <grid, block >> > (*this, nullptr, n_, other);
	gpuChkKer();
//	if (n_ > 0)
//		gpuErrchk(cudaMemcpy(data_, other.data_, memSize(), cudaMemcpyDeviceToDevice));
}

template Vector<double>::Vector(const Vector<float>&);
template Vector<float>::Vector(const Vector<double>&);



template <typename T>
Vector<T>& Vector<T>::operator= (const Vector<T>& other)
{
	if (this != &other)
	{
		resize(other.size(), other.interval());
		gpuErrchk(cudaMemcpy(data_, other.data_, memSize(), cudaMemcpyDeviceToDevice));
	}
	return *this;
}



template <typename T>
template<typename To>
Vector<T>& Vector<T>::operator= (const Vector<To>& other)
{
		resize(other.size(), other.interval() || (std::is_same<T, float>() && std::is_same<To, double>()));
		dim3 block(1024, 1, 1);
		dim3 grid((n_ + 1023) / 1024, 1, 1);
		if (other.interval())
			VectorSelect<10, Intv<To>, Intv<T>> << <grid, block >> > (*this, nullptr, n_, other);
		else
			if (interval_)
				VectorSelect<10, To, Intv<T>> << <grid, block >> > (*this, nullptr, n_, other);
			else
				VectorSelect<10, To, T> << <grid, block >> > (*this, nullptr, n_, other);
		gpuChkKer();
		//gpuErrchk(cudaMemcpy(data_, other.data_, memSize(), cudaMemcpyDeviceToDevice));
	return *this;
}

template Vector<double>& Vector<double>::operator= (const Vector<float>&);
template Vector<float>& Vector<float>::operator= (const Vector<double>&);



template class Vector<double>;
template class Vector<float>;