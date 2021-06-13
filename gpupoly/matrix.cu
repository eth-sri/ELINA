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

 /*! \file src/matrix.cu
	 \brief Implementation of Matrix
     \author Fran&ccedil;ois Serre

	 Implementation of non arithmetic members of the class Matrix, defiined in src/matrix.h.
 */

#include "matrix.h"

#define GLOBAL_ALIGNMENT_SIZE 128

size_t lowestMultiple(const size_t n, const size_t m)
{
	return (n / m + (n % m != 0)) * m;
}

template<typename T>
GPUMem<false> Matrix<T>::gpuMem;

template<typename T>
Matrix<T>::Matrix() :capacity(0), data_(NULL),m_(0),n_(0),N_(0){}

template<typename T>
Matrix<T>::Matrix(size_t m, size_t n, bool interval) :capacity(0),data_(NULL)
{
	reshape(m, n, interval);
}

template<typename T>
Matrix<T>::Matrix(const std::vector<std::vector<T>>& M) :capacity(0), data_(NULL)
{
	m_ = M.size();
	if (m_ == 0)
	{
		interval_ = false;
		n_ = 0;
		N_ = 0;
		data_ = NULL;
		capacity = 0;
	}
	else
	{
		reshape(m_, M.front().size(), false);
		T* dataHost;
		gpuErrchk(cudaMallocHost(&dataHost, m_ * n_ * sizeof(T)));
		for (size_t i = 0; i < m_; i++)
		{
			assert(M[i].size() == n_);
			for (size_t j = 0; j < n_; j++)
				dataHost[i * n_ + j] = M[i][j];
		}
		gpuErrchk(cudaMemcpy2D(data_, N_ * sizeof(T), dataHost, n_ * sizeof(T), n_ * sizeof(T), m_, cudaMemcpyHostToDevice));
		cudaFreeHost(dataHost);
	}
}

template<typename T>
Matrix<T>::Matrix(size_t m, size_t n, const T* data) :Matrix<T>(m, n, false)
{
	if (m * n > 0)
		gpuErrchk(cudaMemcpy2D(data_, pitchBytes(), data, n*sizeof(T), n*sizeof(T), m, cudaMemcpyHostToDevice));
}

template<typename T>
Matrix<T>::Matrix(const Matrix<T>& other):Matrix(other.m_,other.n_,other.interval_)
{
	if(m_ * n_>0)
		gpuErrchk(cudaMemcpy2D(data_, pitchBytes(), other.data_, other.pitchBytes(),(1+interval_)* n_ * sizeof(T), m_, cudaMemcpyDeviceToDevice));
}

template<typename T>
Matrix<T>::Matrix(Matrix<T>&& other)
{
	interval_ = other.interval_;
	n_ = other.n_;
	N_ = other.N_;
	m_ = other.m_;
	data_ = other.data_;
	capacity = other.capacity;
	other.capacity = 0;
	other.data_ = NULL;
}

template<typename T>
Matrix<T>& Matrix<T>::operator= (Matrix<T>&& other)
{
	if (&other != this)
	{
		n_ = other.n_;
		N_ = other.N_;
		m_ = other.m_;
		interval_ = other.interval_;
		if (data_)
			gpuMem.free(data_,capacity);
		data_ = other.data_;
		capacity = other.capacity;
		other.capacity = 0;
		other.data_ = NULL;
	}
	return *this;
}

template<typename T>
Matrix<T>& Matrix<T>::operator= (const Matrix<T>& other)
{
	if (&other != this)
	{
		reshape(other.m_, other.n_, other.interval_);
		if (m_ * n_ > 0)
		gpuErrchk(cudaMemcpy2D(data_, N_ * sizeof(T), other.data_, N_ * sizeof(T), (1+interval_)*n_ * sizeof(T), m_, cudaMemcpyDeviceToDevice));
	}
	return *this;
}

template<typename T>
void Matrix<T>::reshape(const size_t m, const size_t n, const bool interval)
{
	m_ = m;
	n_ = n;
	interval_ = interval;
	N_ = lowestMultiple(interval ? 2 * n : n, GLOBAL_ALIGNMENT_SIZE / sizeof(T));
	if (capacity < memSize())
	{
		if (data_)
			gpuMem.free(data_,capacity);
		capacity = memSize();
		data_=static_cast<T*>(gpuMem.alloc(capacity));
	}
}

template <int blockSize, typename TA, typename TB, typename TD,typename T>
__global__ void addMatrix(
	TD* dest, const size_t destN,
	const TA* lhs, const size_t lhsN,
	const TB* rhs, const size_t rhsN,
	int n)
{
	int i = blockIdx.y;
	int j = threadIdx.x + blockIdx.x * blockSize;
	if (j < n)
		Intv<T>::add(dest[i * destN + j], lhs[i * lhsN + j], rhs[i * rhsN + j]);
}

template<typename T>
void Matrix<T>::add(Matrix<T>& dest, const Matrix<T>& lhs, const Matrix<T>& rhs,bool sound)
{
	const int blockSize = 256;
	assert(lhs.n_ == rhs.n_);
	assert(lhs.m_ == rhs.m_);
	
	dim3 block(blockSize, 1, 1);
	dim3 grid((lhs.n_+blockSize-1)/blockSize, lhs.m_, 1);
	if (sound)
	{
		dest.reshape(lhs.m_, lhs.n_, true);
		if (lhs.interval_)
		{
			if (rhs.interval_)
				addMatrix<blockSize, Intv<T>, Intv<T>,Intv<T>,T> << <grid, block >> > (dest, dest.pitch(), lhs, lhs.pitch(), rhs, rhs.pitch(), lhs.n());
			else
				addMatrix<blockSize, Intv<T>, T, Intv<T>,T> << <grid, block >> > (dest, dest.pitch(), lhs, lhs.pitch(), rhs, rhs.pitch(), lhs.n());
		}
		else
		{
			if (rhs.interval_)
				addMatrix<blockSize, T, Intv<T>, Intv<T>, T> << <grid, block >> > (dest, dest.pitch(), lhs, lhs.pitch(), rhs, rhs.pitch(), lhs.n());
			else
				addMatrix<blockSize, T, T, Intv<T>, T> << <grid, block >> > (dest, dest.pitch(), lhs, lhs.pitch(), rhs, rhs.pitch(), lhs.n());
		}
	}
	else {
		assert(!lhs.interval());
		assert(!rhs.interval());
		dest.reshape(lhs.m_, lhs.n_, false);
#if false // the following lines use cublas for adding two scalar matrices... but it seems to be faster to use directly our code
		double alpha = 1;
		cublasDgeam(cublasHandle,
			CUBLAS_OP_N, CUBLAS_OP_N,
			lhs.n(), lhs.m(),
			&alpha,
			lhs, lhs.pitch(),
			&alpha,
			rhs, rhs.pitch(),
			dest, dest.pitch()
		);
#else
		addMatrix<blockSize, T, T, T,T> << <grid, block >> > (dest, dest.pitch(), lhs, lhs.pitch(), rhs, rhs.pitch(), lhs.n());
#endif
	}
}

template<typename T>
template<typename Tr>
void Vector<T>::add(Vector<T>& dest, const Vector<T>& lhs, const Vector<Tr>& rhs)
{
	const int blockSize = 256;
	assert(lhs.size() == rhs.size());
	dest.resize(lhs.n_, true);
	dim3 block(blockSize, 1, 1);
	dim3 grid((lhs.n_ + blockSize - 1) / blockSize,1, 1);

	if (lhs.interval())
	{
		if (rhs.interval())
			addMatrix<blockSize, Intv<T>, Intv<Tr>,Intv<T>,T> << <grid, block >> > (dest, dest.size(), lhs, lhs.size(), rhs, rhs.size(), lhs.size());
		else
			addMatrix<blockSize, Intv<T>, Tr, Intv<T>, T> << <grid, block >> > (dest, dest.size(), lhs, lhs.size(), rhs, rhs.size(), lhs.size());
	}
	else
	{
		if (rhs.interval())
			addMatrix<blockSize, T, Intv<Tr>,  Intv<T>, T> << <grid, block >> > (dest, dest.size(), lhs, lhs.size(), rhs, rhs.size(), lhs.size());
		else
			addMatrix<blockSize, T, Tr, Intv<T>, T> << <grid, block >> > (dest, dest.size(), lhs, lhs.size(), rhs, rhs.size(), lhs.size());
	}
}
template void Vector<double>::add(Vector<double>& dest, const Vector<double>& lhs, const Vector<double>& rhs);
template void Vector<double>::add(Vector<double>& dest, const Vector<double>& lhs, const Vector<float>& rhs);
template void Vector<float>::add(Vector<float>& dest, const Vector<float>& lhs, const Vector<double>& rhs);
template void Vector<float>::add(Vector<float>& dest, const Vector<float>& lhs, const Vector<float>& rhs);

template <int blockSize, bool up, typename T>
__global__ void addMatrix_dr(
	T* dest, const size_t destN,
	const T* lhs, const size_t lhsN,
	const T* rhs, const size_t rhsN,
	int n)
{
	int i = blockIdx.y;
	int j = threadIdx.x + blockIdx.x * blockSize;
	if (j < n)
		dest[i * destN + j] = Intv<T>::template add_dr<up>(lhs[i * lhsN + j], rhs[i * rhsN + j]);
}

template<typename T>
template <bool up>
std::shared_ptr<const Vector<T>> Vector<T>::add_dr(std::shared_ptr<const Vector<T>> lhs, std::shared_ptr<const Vector<T>> rhs, std::shared_ptr<Vector<T>> dest)
{
	const int blockSize = 256;
	assert(!lhs || !lhs->interval());
	assert(!rhs || !rhs->interval());
	if (!lhs)
		return rhs;
	if (!rhs)
		return lhs;
	assert(lhs->n_ == rhs->n_);
	if (dest)
		dest->resize(lhs->n_, false);
	else
		dest = std::make_shared<Vector<T>>(lhs->n_, false);
	dim3 block(blockSize, 1, 1);
	dim3 grid((lhs->n_ + blockSize - 1) / blockSize, 1, 1);
	addMatrix_dr<blockSize,up,T> << <grid, block >> > (*dest, dest->n_, *lhs, lhs->n_, *rhs, rhs->n_, lhs->n_);
	return dest;
}

template
std::shared_ptr<const Vector<double>> Vector<double>::add_dr<true>(std::shared_ptr<const Vector<double>> lhs, std::shared_ptr<const Vector<double>> rhs, std::shared_ptr<Vector<double>> dest);
template
std::shared_ptr<const Vector<double>> Vector<double>::add_dr<false>(std::shared_ptr<const Vector<double>> lhs, std::shared_ptr<const Vector<double>> rhs, std::shared_ptr<Vector<double>> dest);

template
std::shared_ptr<const Vector<float>> Vector<float>::add_dr<true>(std::shared_ptr<const Vector<float>> lhs, std::shared_ptr<const Vector<float>> rhs, std::shared_ptr<Vector<float>> dest);
template
std::shared_ptr<const Vector<float>> Vector<float>::add_dr<false>(std::shared_ptr<const Vector<float>> lhs, std::shared_ptr<const Vector<float>> rhs, std::shared_ptr<Vector<float>> dest);


template <int lgBlockSize,typename Tout, typename Tin>
__global__ void backSubstituteDense(Tout* A,const size_t pitchOut, const int* annoyingList, const size_t n,  const Tin* oA, const size_t pitchIn)
{
	size_t i = (blockIdx.y << lgBlockSize) + threadIdx.x;
	size_t row = annoyingList?annoyingList[blockIdx.x]: blockIdx.x;
	if (i < n)
		A[blockIdx.x * pitchOut + i] = oA[row * pitchIn + i];
}

template<typename T>
template<typename Td>
Matrix<Td> Matrix<T>::selectRows(size_t size, const int* rows, bool forceIntervalOut) const
{
	bool intervalOut = interval_ || forceIntervalOut;
	Matrix<Td> res(size, n_, intervalOut);
	dim3 block(1024, 1, 1);
	dim3 grid(size,(n_ + 1023) / 1024, 1);
	if(interval_)
		backSubstituteDense<10,Intv<Td>,Intv<T>> << <grid, block >> > (res,res.pitch(), rows, n_, *this,pitch() );
	else
		if(intervalOut)
			backSubstituteDense<10, Intv<Td>,T> << <grid, block >> > (res,res.pitch(), rows, n_,  *this,pitch());
		else
			backSubstituteDense<10, Td,T> << <grid, block >> > (res, res.pitch(), rows, n_,  *this,pitch());
	gpuChkKer();
	return res;
}

template Matrix<double> Matrix<double>::selectRows(size_t size, const int* rows, bool forceIntervalOut) const;
template Matrix<float> Matrix<double>::selectRows(size_t size, const int* rows, bool forceIntervalOut) const;
template Matrix<double> Matrix<float>::selectRows(size_t size, const int* rows, bool forceIntervalOut) const;
template Matrix<float> Matrix<float>::selectRows(size_t size, const int* rows, bool forceIntervalOut) const;


template<typename T>
template<typename To>
Matrix<T>::Matrix(const Matrix<To>& other, bool sound) :Matrix(other.m(), other.n(), other.interval() || (sound && std::is_same<T, float>() && std::is_same<To, double>()))
{
	dim3 block(1024, 1, 1);
	dim3 grid(m(),(n_ + 1023) / 1024, 1);
	gpuChkKer();
	if (other.interval())
		backSubstituteDense<10, Intv<T>, Intv<To>> << <grid, block >> > (*this,pitch(), nullptr, n_, other,other.pitch());
	else
		if (interval())
			backSubstituteDense<10, Intv<T>, To> << <grid, block >> > (*this, pitch(), nullptr, n_, other, other.pitch());
		else
			backSubstituteDense<10, T, To> << <grid, block >> > (*this, pitch(), nullptr, n_, other, other.pitch());
	gpuChkKer();
}
template Matrix<double>::Matrix(const Matrix<float>&,bool);
template Matrix<float>::Matrix(const Matrix<double>&,bool);
template Matrix<double>::Matrix(const Matrix<double>&, bool);
template Matrix<float>::Matrix(const Matrix<float>&, bool);


template<>
Matrix<double> Matrix<double>::transpose() const {
	assert(!interval());
	Matrix<double> res(n(), m(), false);
	double alpha = 1;
	double beta = 0;
	cublasDgeam(cublasHandle, CUBLAS_OP_T, CUBLAS_OP_T,
		m(), n(),
		&alpha,
		*this, pitch(),
		&beta,
		*this, pitch(),
		res, res.pitch()
	);
	return res;
}

template<>
Matrix<float> Matrix<float>::transpose() const {
	assert(!interval());
	Matrix<float> res(n(), m(), false);
	float alpha = 1;
	float beta = 0;
	cublasSgeam(cublasHandle, CUBLAS_OP_T, CUBLAS_OP_T,
		m(), n(),
		&alpha,
		*this, pitch(),
		&beta,
		*this, pitch(),
		res, res.pitch()
	);
	return res;
}

template <typename T> 
typename Matrix<T>::CublasHandle Matrix<T>::cublasHandle;

template class Matrix<double>;

template class Matrix<float>;
