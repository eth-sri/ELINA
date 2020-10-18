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


/*! \file src/matrix.h
	\brief Definition of the class Matrix
	\author Fran&ccedil;ois Serre
*/

#pragma once
#include <vector>
#include <iostream>
#include "cublas_v2.h"
#include "intv.h"
#include "vector.h"
#include "gpumem.h"

//! Class that represents a matrix stored in GPU memory containing either real or interval values. Methods are implemented in src/matrix.cu, src/mmm.cu and src/mvm.cu.
template<typename T>
class Matrix
{
	static class CublasHandle
	{
		cublasHandle_t* handle;
	public:
		CublasHandle() :handle(nullptr) {}
		~CublasHandle()
		{
			if (handle)
			{
				cublasDestroy(*handle);
				delete handle;
			}
		}
		operator cublasHandle_t& () {
			if (!handle)
			{
				handle = new cublasHandle_t;
				cublasCreate(handle);
			}
			return *handle;
		}
	} cublasHandle;
public:
	static GPUMem<false> gpuMem;
private:
	bool interval_;
	size_t m_;
	size_t n_;
	size_t N_; // gap between two rows in A (smallest multiple of alignment quantum bigger than n)
	T* data_;
	size_t capacity;
public:
	//! Size in bytes of the data (including alignment padding)
	inline size_t memSize() const
	{
		return m_ * N_ * sizeof(T);
	}
	//! Size in bytes of each row (including alignment padding)
	size_t pitchBytes() const
	{
		return N_ * sizeof(T);
	}
	//! Number of elements in each row, including padding
	size_t pitch() const
	{
		return interval_ ? N_ / 2 : N_;
	}
	//! Constructor
	/*!
	Creates a matrix of size 0x0
	*/
	Matrix();
	//! Constructor
	/*!
	Creates (and allocates) a new matrix
	\param m Number of rows
	\param n Number of columns
	\param interval If true, each element of the matrix is an interval (Intv), otherwise a scalar (T)
	*/
	Matrix(const size_t m, const size_t n, const bool interval = false);
	//! Constructor
	/*!
	Creates (and allocates) a new scalar matrix, and store the elements from the host memory.
	\param M Data of the matrix, stored as a std::vector of rows, each row being itself an std::vector.
	*/
	Matrix(const std::vector<std::vector<T>>& M);
	//! Constructor
	/*!
	Creates (and allocates) a new scalar matrix, and store the elements from the host memory.
	\param m Number of rows
	\param n Number of columns
	\param data (Host) Pointer to the data (in row major order, with no padding)
	*/
	Matrix(size_t m, size_t n, const T* data);
	//! Copy constructor
	Matrix(const Matrix<T>& other);
	//! Move constructor
	Matrix(Matrix<T>&& other);
	//! Assignment operator
	Matrix<T>& operator= (const Matrix<T>& other);
	//! Move assignment operator
	Matrix<T>& operator= (Matrix<T>&& other);
	//! Destructor
	~Matrix()
	{
		if (data_)
			gpuMem.free(data_,capacity);
	}
	//! Reshape a matrix
	/*!
	Change the size of a matrix
	\warning After this operation, the data contained in the matrix is no longer valid
	*/
	void reshape(const size_t m, const size_t n, const bool interval);
	//! Number of columns (without padding)
	inline size_t n() const
	{
		return n_;
	}
	//! Number of rows
	inline size_t m() const
	{
		return m_;
	}
	//! If true, the matrix contains interval values (Intv); otherwise scalar values (T)
	inline bool interval() const
	{
		return interval_;
	}
	/*inline size_t N() const
	{
		return N_;
	}*/
	//! Get a pointer to the data
	inline void* data() {
		return data_;
	}
	//! Get a const pointer to the data
	inline const void* data() const {
		return data_;
	}
	//! Get a pointer to the data
	operator Intv<T>* () {
		assert(interval_);
		return reinterpret_cast<Intv<T>*>(data_);
	}
	//! Get a const pointer to the data
	operator const Intv<T>* () const {
		assert(interval_);
		return reinterpret_cast<const Intv<T>*>(data_);
	}
	//! Get a pointer to the data
	operator T* () {
		assert(!interval_);
		return data_;
	}
	//! Get a const pointer to the data
	operator const T* () const {
		assert(!interval_);
		return data_;
	}
	//! Get the transposed matrix
	Matrix<T> transpose() const;

	void zeroFill()
	{
		cudaMemset2D(data(), pitchBytes(), 0, n() * sizeof(T) * (1 + interval()), m());
	}
	//! Select rows
	/*!
	Creates a new matrix that contains a selection of rows of the current one.
	\param size Number of selected rows
	\param rows An array containing the index of the row to select (pointer to a GPU array)
	*/
	template <typename Td>
	Matrix<Td> selectRows(size_t size, const int* rows) const;


	template <typename To>
	Matrix(const Matrix<To>& other, bool sound);
	//! Matrix vector multiplication with directed gathering
	/*!
	Computes a (sound) matrix-vector multiplication, and gather for each element an endpoint.
	\param up Whether the upper endpoint should be gathered (true), or the lower one (false)
	\param dest Destination Vector
	\param rhs Vector the current matrix has to be multiplied to
	\param offset An optional bias to add (must be a scalar Vector)
	*/
	template<bool up, typename Tr>
	void mvm_dr(Vector<T>& dest, const Vector<Tr>& rhs, std::shared_ptr<const Vector<T>> offset=nullptr) const;
	//! Matrix vector multiplication.
	/*!
	Computes a (sound) matrix-vector multiplication.
	\param dest Destination Vector
	\param rhs Vector the current matrix has to be multiplied to
	\param offset An optional bias to add (must be a scalar Vector)
	*/
	template<typename Tr>
	void mvm(Vector<Tr>& dest, const Vector<Tr>& rhs) const;
	//! Matrix-matrix multiplication.
	/*!
	Computes a (sound) matrix-matrix multiplication.
	\param rhs right-hand side Matrix
	\param dest Destination Matrix
	*/
	template <typename Tr>
	void mmm(const Matrix<Tr>& rhs, Matrix<T>& dest, bool sound=true) const;
//! Prints the content of the Matrix.
	void print() const {
		if (interval_)
		{
			std::vector<Intv<T>> tmp(m_ * n_);
			cudaMemcpy2D(
				tmp.data(), n_*sizeof(Intv<T>),
				data_, pitchBytes(),
				n_ * sizeof(Intv<T>),m_,
				cudaMemcpyDeviceToHost);
			for (int i = 0; i < m_; i++)
			{
				for (int j = 0; j < n_; j++)
					std::cout << "["<< tmp[i*n_+j].low << " " << tmp[i * n_ + j].high << "] ";
				std::cout << std::endl;
			}
			std::cout << std::endl;
		}
		else
		{
			std::vector<T> tmp(m_ * n_);
			cudaMemcpy2D(
				tmp.data(), n_ * sizeof(T),
				data_, pitchBytes(),
				n_ * sizeof(T), m_,
				cudaMemcpyDeviceToHost);
			for (int i = 0; i < m_; i++)
			{
				for (int j = 0; j < n_; j++)
					std::cout << tmp[i * n_ + j] << " ";
				std::cout << std::endl;
			}
			std::cout << std::endl;
		}
	}
	//! Matrix-matrix addition.
/*!
Computes the sum of two matrices
\param dest Destination Matrix
\param lhs left-hand side Matrix
\param rhs right-hand side Matrix
*/
	static void add(Matrix<T>& dest, const Matrix<T>& lhs, const Matrix<T>& rhs, bool sound=true);

};

