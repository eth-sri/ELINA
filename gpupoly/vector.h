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

/*! \file src/vector.h
	\brief Vectors stored in GPU memory containing real or interval values.
	\author Fran&ccedil;ois Serre
*/


#pragma once
#include "cuda_runtime.h"
#include "utils.h"
#include "intv.h"
#include "gpumem.h"
#include <thrust/device_ptr.h>
#include <thrust/copy.h>
#include <vector>
#include <iostream>
#include <limits>

/*! Vector of real or interval values stored in GPU memory
*/
template<typename T>
class Vector
{
	//friend class Matrix<T>;
public:
	static GPUMem<false> gpuMem;
private:
	/// True if elements are real intervals
	bool interval_;

	/// Size in number of elements
	size_t n_;

	/// Pointer to the data in GPU memory
	T* data_;

	/// Max. capacity of the underlying storage
	size_t capacity;

public:
	Vector();
	Vector(const size_t n, const bool interval = false);
	Vector(const std::vector<T>& v);
	Vector(size_t size, const T* v);
	Vector(const std::vector<Intv<T>>& v);
	Vector(const Vector<T>& other);
	template<typename To> 
	Vector(const Vector<To>& other);
	Vector(Vector<T>&& other);
	Vector& operator= (const std::vector<T>& other);
	Vector& operator= (const std::vector<Intv<T>>& other);
	Vector& operator= (const Vector<T>& other);
	template<typename To>
	Vector& operator= (const Vector<To>& other);
	Vector& operator= (Vector<T>&& other);
	operator Intv<T>* () {
		assert(interval_);
		return reinterpret_cast<Intv<T>*>(data_);
	}
	operator const Intv<T>* () const {
		assert(interval_);
		return reinterpret_cast<const Intv<T>*>(data_);
	}
	operator T* () {
		assert(!interval_);
		return data_;
	}
	operator const T* () const {
		assert(!interval_);
		return data_;
	}
	~Vector();
	inline size_t size() const
	{
		return n_;
	}
	void resize(size_t n, bool interval);
	inline thrust::device_ptr<T> begin() const
	{
		return thrust::device_pointer_cast<T>(data_);
	}
	inline thrust::device_ptr<T> end() const
	{
		return thrust::device_pointer_cast<T>(data_) + (interval_ ? 2 * n_ : n_);
	}
	inline thrust::device_ptr<const Intv<T>> beginIntv() const
	{
		return thrust::device_pointer_cast<const Intv<T>>(*this);
	}
	inline thrust::device_ptr<const Intv<T>> endIntv() const
	{
		return thrust::device_pointer_cast<const Intv<T>>(*this) + n_;
	}
	inline const T* data() const 
	{
		return data_;
	}
	/*inline T* data()
	{
		return data_;
	}*/
	inline bool interval() const
	{
		return interval_;
	}
	Vector<T> select(size_t size, const int* rows) const;
	bool isPositive() const;
	void print() const {
		if (interval_)
		{
			std::vector<Intv<T>> tmp(n_);
			cudaMemcpy(tmp.data(), data_, n_ * sizeof(Intv<T>), cudaMemcpyDeviceToHost);
			for (int i = 0; i < tmp.size(); i++)
			{
				std::cout << "[" << tmp[i].low << " " << tmp[i].high << "] ";
			}
			std::cout<<std::endl;
		}
		else
		{
			std::vector<T> tmp(n_);
			cudaMemcpy(tmp.data(), data_, n_ * sizeof(T), cudaMemcpyDeviceToHost);
			//double max = -std::numeric_limits<double>::infinity();
			//int argmax = -1;
			for (int i = 0; i < tmp.size(); i++)
			{
				std::cout << tmp[i] << " ";
				/*if (tmp[i].high > max)
				{
					max = tmp[i].high;
					argmax = i;
				}*/
			}
			std::cout << std::endl;
			//std::cout << std::endl << argmax << std::endl;
		}
	}
	void check() const {
		if (interval_)
		{
			std::vector<Intv<T>> tmp(n_);
			cudaMemcpy(tmp.data(), data_, n_ * sizeof(Intv<T>), cudaMemcpyDeviceToHost);
			for (auto i:tmp)
				if (i.high < i.low)
				{
					print();
					throw - 1;
				}
		}
	}
	template <typename Tr>
	static void add(Vector<T>& dest, const Vector<T>& lhs, const Vector<Tr>& rhs);
	template <bool up>
	static  std::shared_ptr<const Vector<T>> add_dr( std::shared_ptr<const Vector<T>> lhs, std::shared_ptr<const Vector<T>> rhs, std::shared_ptr<Vector<T>> dest=NULL);
	inline size_t memSize() const
	{
		return interval_ ? n_ * sizeof(Intv<T>) : n_ * sizeof(T);
	}
	
};

