/*
 * GPUPoly library
 * A deep neural network verifier library running on GPU
 * Copyright (C) 2020 Francois Serre
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA
 */

 /*! \file src/intv.h
	 \brief Interval representation and interval arithmetic on GPU.
 */
 /*!
   \file src/intv.h
   \brief Interval representation and basic arithmetic interval on GPU.
   \author Fran&ccedil;ois Serre

   Definition of the class Intv that represents an interval, and sound arithmetic operation on them, on GPU.
 */

#pragma once
#include <cuda.h>
#include "config.h"
#define CUDART_INF_F            __int_as_float(0x7f800000)
#define CUDART_INF              __longlong_as_double(0x7ff0000000000000ULL)
 //! Interval
 /*! Real interval represented by its endpoints, i.e. the set { x | low <= x <= high}. The arithmetic operations defined on them are sound, and exploit the directed rounding instructions on GPU.
 */
template<typename T>
struct Intv
{
	T low; //!< Lower endpoint

	T high; //!< Higher endpoint

	//! Constructor
	/*!
	Creates an interval only containing 0
	*/
	__device__ __host__ Intv() : low(0), high(0) {}

	//! Constructor
	/*!
	Creates an interval defined by its endpoints.
	\param low Lower endpoint
	\param high Higher endpoint
	*/
	__device__ __host__ Intv(const T& low, const T& high) :low(low), high(high) {}

	//! Constructor
	/*!
	Creates a singleton.
	\param i Element that the singleton contains.
	*/
	__device__ __host__ Intv(const T& i) : low(i), high(i) {}

	/*/// Converts an integer to a trivial interval... deleted because bugprone
	CUTLASS_HOST_DEVICE Intv(const int i) : low(i), high(i) {}*/

	/// Multiplication of two intervals
	__device__ Intv<T> operator*(const Intv<T>& rhs) const;

	/// Addition of two intervals
	__device__ Intv<T> operator+(const Intv<T>& rhs) const;

	/// Test if two intervals are different
	__device__ bool operator!=(const Intv<T>& rhs) const;

	/// Inplace addition of intervals
	__device__ Intv<T>& operator+=(const Intv<T>& rhs);

	/// Inplace multiplication of intervals
	__device__ Intv<T>& operator*=(const Intv<T>& rhs);

	Intv<T>& operator=(const double& src);
	Intv<T>& operator=(const float& src);
	Intv<T>& operator=(const Intv<double>& src);
	Intv<T>& operator=(const Intv<float>& src);
	

	/*! \name Fused multiply-add operations
	The following function computes d = a * b + c, with the terms a and b being either a scalar or an interval. These are floating point sound. If strong soundness is enabled, using their results as terms of a sum makes this sum sound even with an arbitrary order.
	*/
	template<typename Ta, typename Tb>
	static __device__ void fma(Intv<T>& d, const Ta& a, const Tb& b, const Intv<T>& c);
	/*! \name Fused multiply-add operations
	The following function computes d = a * b + c without floating point soundness.
	*/
	template<typename Tb>
	static __device__ void fma(T& d, const T& a, const Tb& b, const T& c);

	/*! \name Add operations
	The following functions compute d = a + b, with the terms a and b being either a scalar or an interval, and c a scalar.
	*/
	template<typename Ta, typename Tb>
	static __device__ void add(Intv<T>& d, const Ta& a, const Tb& b);
	template<typename Tb>
	static __device__ void add(T& d, const T& a, const Tb& b);

	/*! \name Fused Multiply Add operation with directed gathering
	The following function computes a * b + c, with the terms a and b being either a scalar or an interval, and c a scalar, and returns the highest (resp. lowest) bound of the result if up is set to true (resp. false).
	*/
	template <bool up, typename Ta, typename Tb>
	static __device__ T fma_dr(const Ta& ta, const Tb& tb, const T& tc);

	template <bool up>
	static __device__ T mul_dr(const T& ta, const T& tb);
	template <bool up>
	static __device__ T add_dr(const T& ta, const T& tb);
	static __device__ T max(const T& ta, const T& tb);
	static __device__ T min(const T& ta, const T& tb);


	
	template <bool up, typename Ts>
	static __device__ T access_dr(const Ts& input);
	template <bool up, typename Ts>
	static __device__ T access_dr(const Intv<Ts>& input);
	template <bool up>
	static __device__ T& access_dr(T& input);
	template <bool up>
	static __device__ T& access_dr(Intv<T>& input);

	/*template<typename Ts>
	static __device__ void convert(Intv<T>& dest, const Ts& source);
	template<typename Ts>
	static __device__ void convert(T& dest, const Ts& source);*/
};




template<> inline __device__ Intv<double> Intv<double>::operator*(const Intv<double>& rhs) const { return Intv(fmin(fmin(__dmul_rd(this->low, rhs.low), __dmul_rd(this->low, rhs.high)), fmin(__dmul_rd(this->high, rhs.low), __dmul_rd(this->high, rhs.high))), fmax(fmax(__dmul_ru(this->low, rhs.low), __dmul_ru(this->low, rhs.high)), fmax(__dmul_ru(this->high, rhs.low), __dmul_ru(this->high, rhs.high)))); }
template<> inline __device__ Intv<float> Intv<float>::operator*(const Intv<float>& rhs) const { return Intv(fminf(fminf(__fmul_rd(this->low, rhs.low), __fmul_rd(this->low, rhs.high)), fminf(__fmul_rd(this->high, rhs.low), __fmul_rd(this->high, rhs.high))), fmaxf(fmaxf(__fmul_ru(this->low, rhs.low), __fmul_ru(this->low, rhs.high)), fmaxf(__fmul_ru(this->high, rhs.low), __fmul_ru(this->high, rhs.high)))); }

template<typename T> inline __device__ Intv<T> Intv<T>::operator+(const Intv<T>& rhs) const { Intv<T> res; Intv<T>::add(res, *this, rhs);	return res; }

template<typename T> inline __device__ bool Intv<T>::operator!=(const Intv<T>& rhs) const { return (this->high != rhs.high) || (this->low != rhs.low); }

template<typename T> inline __device__ Intv<T>& Intv<T>::operator+=(const Intv<T>& rhs) { Intv<T>::add(*this, *this, rhs);	return *this; }

template<> inline __device__ Intv<double>& Intv<double>::operator*=(const Intv<double>& rhs) { double tmp = fmin(fmin(__dmul_rd(low, rhs.low), __dmul_rd(low, rhs.high)), fmin(__dmul_rd(high, rhs.low), __dmul_rd(high, rhs.high)));	high = fmax(fmax(__dmul_ru(low, rhs.low), __dmul_ru(low, rhs.high)), fmax(__dmul_ru(high, rhs.low), __dmul_ru(high, rhs.high)));	low = tmp;	return *this; }
template<> inline __device__ Intv<float>& Intv<float>::operator*=(const Intv<float>& rhs) { float tmp = fminf(fminf(__fmul_rd(low, rhs.low), __fmul_rd(low, rhs.high)), fminf(__fmul_rd(high, rhs.low), __fmul_rd(high, rhs.high))); high = fmaxf(fmaxf(__fmul_ru(low, rhs.low), __fmul_ru(low, rhs.high)), fmaxf(__fmul_ru(high, rhs.low), __fmul_ru(high, rhs.high)));	low = tmp;	return *this; }


template<> inline __device__ __host__ Intv<double>& Intv<double>::operator=(const Intv<double>& source) { low = source.low; high = source.high; return *this; }
template<> inline __device__ __host__ Intv<double>& Intv<double>::operator=(const double& source) { low = source; high = source; return *this; }
template<> inline __device__ __host__ Intv<double>& Intv<double>::operator=(const Intv<float>& source) { low = source.low; high = source.high; return *this; }
template<> inline __device__ __host__ Intv<double>& Intv<double>::operator=(const float& source) { low = source; high = source; return *this; }
template<> inline __device__ Intv<float>& Intv<float>::operator=(const Intv<double>& source) { low = __double2float_rd(source.low); high = __double2float_ru(source.high); return *this; }
template<> inline __device__ Intv<float>& Intv<float>::operator=(const double& source) { low = __double2float_rd(source); high = __double2float_ru(source); return *this; }
template<> inline __device__ __host__ Intv<float>& Intv<float>::operator=(const Intv<float>& source) { low = source.low; high = source.high; return *this; }
template<> inline __device__ __host__ Intv<float>& Intv<float>::operator=(const float& source) { low = source; high = source; return *this; }
/*template <> template<> inline __device__ void convert<Intv<double>>::operator()(Intv<double>& dest, const Intv<float>& source) { dest = Intv<double>(source.low, source.high); }
template <> template<> inline __device__ void convert<double>::operator()(double& dest, const float& source) { dest = source; }
template <> template<> inline __device__ void convert< Intv<float>>::operator()(Intv<float>& dest, const Intv<double>& source) { dest = Intv<float>(__double2float_rd(source.low), __double2float_ru(source.high)); }
template <> template<> inline __device__ void convert< Intv<float>>::operator()(Intv<float>& dest, const double& source) { dest = Intv<float>(__double2float_rd(source), __double2float_ru(source)); }
*/

#ifdef STRONG_FP_SOUNDNESS
// strong floating point soundness (i.e. user evaluates with precision T, in any order for additions)
template<> template<> inline __device__ void Intv<double>::fma(Intv<double>& d, const Intv<double>& a, const Intv<double>& b, const Intv<double>& c) { double prod1 = a.low * b.low;	double prod2 = a.low * b.high;	double prod3 = a.high * b.low;	double prod4 = a.high * b.high;	d.low = __dadd_rd(nextafter(fmin(fmin(prod1, prod2), fmin(prod3, prod4)), -CUDART_INF), c.low);	d.high = __dadd_ru(nextafter(fmax(fmax(prod1, prod2), fmax(prod3, prod4)), CUDART_INF), c.high); }
template<> template<> inline __device__ void Intv<double>::fma(Intv<double>& d, const Intv<double>& a, const double& b, const Intv<double>& c) { double prod1 = a.low * b;	double prod2 = a.high * b;	d.low = __dadd_rd(nextafter(fmin(prod1, prod2), -CUDART_INF), c.low);	d.high = __dadd_ru(nextafter(fmax(prod1, prod2), CUDART_INF), c.high); }
template<> template<> inline __device__ void Intv<double>::fma(Intv<double>& d, const double& a, const Intv<double>& b, const Intv<double>& c) { double prod1 = a * b.low;	double prod2 = a * b.high;	d.low = __dadd_rd(nextafter(fmin(prod1, prod2), -CUDART_INF), c.low);	d.high = __dadd_ru(nextafter(fmax(prod1, prod2), CUDART_INF), c.high); }
template<> template<> inline __device__ void Intv<double>::fma(Intv<double>& d, const double& a, const double& b, const Intv<double>& c) { double prod = a * b; d.low = __dadd_rd(nextafter(prod, -CUDART_INF), c.low);	d.high = __dadd_ru(nextafter(prod, CUDART_INF), c.high); }
template<> template<> inline __device__ void Intv<float>::fma(Intv<float>& d, const Intv<float>& a, const Intv<float>& b, const Intv<float>& c) { float prod1 = a.low * b.low;	float prod2 = a.low * b.high;	float prod3 = a.high * b.low;	float prod4 = a.high * b.high;	d.low = __fadd_rd(nextafterf(fminf(fminf(prod1, prod2), fminf(prod3, prod4)), -CUDART_INF_F), c.low);	d.high = __fadd_ru(nextafterf(fmaxf(fmaxf(prod1, prod2), fmaxf(prod3, prod4)), CUDART_INF_F), c.high); }
template<> template<> inline __device__ void Intv<float>::fma(Intv<float>& d, const Intv<float>& a, const float& b, const Intv<float>& c) { float prod1 = a.low * b;	float prod2 = a.high * b;	d.low = __fadd_rd(nextafterf(fminf(prod1, prod2), -CUDART_INF_F), c.low);	d.high = __fadd_ru(nextafterf(fmaxf(prod1, prod2), CUDART_INF_F), c.high); }
template<> template<> inline __device__ void Intv<float>::fma(Intv<float>& d, const float& a, const Intv<float>& b, const Intv<float>& c) { float prod1 = a * b.low;	float prod2 = a * b.high;	d.low = __fadd_rd(nextafterf(fminf(prod1, prod2), -CUDART_INF_F), c.low);	d.high = __fadd_ru(nextafterf(fmaxf(prod1, prod2), CUDART_INF_F), c.high); }
template<> template<> inline __device__ void Intv<float>::fma(Intv<float>& d, const float& a, const float& b, const Intv<float>& c) { float prod = a * b;	d.low = __fadd_rd(nextafterf(prod, -CUDART_INF_F), c.low); 	d.high = __fadd_ru(nextafterf(prod, CUDART_INF_F), c.high); }
#else
// weak floating point soundness (i.e. user evaluates with infinite precision)
template<> template<> inline __device__ void Intv<double>::fma(Intv<double>& d, const Intv<double>& a, const Intv<double>& b, const Intv<double>& c) { double tmp = fmin(fmin(__fma_rd(a.low, b.low, c.low), __fma_rd(a.low, b.high, c.low)), fmin(__fma_rd(a.high, b.low, c.low), __fma_rd(a.high, b.high, c.low)));	d.high = fmax(fmax(__fma_ru(a.low, b.low, c.high), __fma_ru(a.low, b.high, c.high)), fmax(__fma_ru(a.high, b.low, c.high), __fma_ru(a.high, b.high, c.high)));	d.low = tmp; }
template<> template<> inline __device__ void Intv<double>::fma(Intv<double>& d, const Intv<double>& a, const double& b, const Intv<double>& c) { double tmp = fmin(__fma_rd(a.low, b, c.low), __fma_rd(a.high, b, c.low));	d.high = fmax(__fma_ru(a.low, b, c.high), __fma_ru(a.high, b, c.high));	d.low = tmp; }
template<> template<> inline __device__ void Intv<double>::fma(Intv<double>& d, const double& a, const Intv<double>& b, const Intv<double>& c) { double tmp = fmin(__fma_rd(a, b.low, c.low), __fma_rd(a, b.high, c.low));	d.high = fmax(__fma_ru(a, b.low, c.high), __fma_ru(a, b.high, c.high));	d.low = tmp; }
template<> template<> inline __device__ void Intv<double>::fma(Intv<double>& d, const double& a, const double& b, const Intv<double>& c) { d.low = __fma_rd(a, b, c.low);	d.high = __fma_ru(a, b, c.high); }
template<> template<> inline __device__ void Intv<float>::fma(Intv<float>& d, const Intv<float>& a, const Intv<float>& b, const Intv<float>& c) { float tmp = fminf(fminf(__fmaf_rd(a.low, b.low, c.low), __fmaf_rd(a.low, b.high, c.low)), fminf(__fmaf_rd(a.high, b.low, c.low), __fmaf_rd(a.high, b.high, c.low)));	d.high = fmaxf(fmaxf(__fmaf_ru(a.low, b.low, c.high), __fmaf_ru(a.low, b.high, c.high)), fmaxf(__fmaf_ru(a.high, b.low, c.high), __fmaf_ru(a.high, b.high, c.high)));	d.low = tmp; }
template<> template<> inline __device__ void Intv<float>::fma(Intv<float>& d, const Intv<float>& a, const float& b, const Intv<float>& c) { float tmp = fminf(__fmaf_rd(a.low, b, c.low), __fmaf_rd(a.high, b, c.low)); d.high = fmaxf(__fmaf_ru(a.low, b, c.high), __fmaf_ru(a.high, b, c.high));	d.low = tmp; }
template<> template<> inline __device__ void Intv<float>::fma(Intv<float>& d, const float& a, const Intv<float>& b, const Intv<float>& c) { float tmp = fminf(__fmaf_rd(a, b.low, c.low), __fmaf_rd(a, b.high, c.low)); d.high = fmaxf(__fmaf_ru(a, b.low, c.high), __fmaf_ru(a, b.high, c.high));	d.low = tmp; }
template<> template<> inline __device__ void Intv<float>::fma(Intv<float>& d, const float& a, const float& b, const Intv<float>& c) { d.low = __fmaf_rd(a, b, c.low); d.high = __fmaf_ru(a, b, c.high); }
#endif

template<> template<> inline __device__ void Intv<double>::fma(Intv<double>& d, const Intv<double>& a, const Intv<float>& b, const Intv<double>& c) { const Intv<double> b_d(b.low, b.high); fma(d, a, b_d, c); }
template<> template<> inline __device__ void Intv<double>::fma(Intv<double>& d, const Intv<double>& a, const float& b, const Intv<double>& c) { const double b_d = b;	return fma(d, a, b_d, c); }
template<> template<> inline __device__ void Intv<double>::fma(Intv<double>& d, const double& a, const Intv<float>& b, const Intv<double>& c) { const Intv<double> b_d(b.low, b.high);	return fma(d, a, b_d, c); }
template<> template<> inline __device__ void Intv<double>::fma(Intv<double>& d, const double& a, const float& b, const Intv<double>& c) { const double b_d = b; return fma(d, a, b_d, c); }
template<> template<> inline __device__ void Intv<float>::fma(Intv<float>& d, const Intv<float>& a, const Intv<double>& b, const Intv<float>& c) { fma(d, a, Intv<float>(__double2float_rd(b.low), __double2float_ru(b.high)), c); }
template<> template<> inline __device__ void Intv<float>::fma(Intv<float>& d, const Intv<float>& a, const double& b, const Intv<float>& c) { fma(d, a, Intv<float>(__double2float_rd(b), __double2float_ru(b)), c); }
template<> template<> inline __device__ void Intv<float>::fma(Intv<float>& d, const float& a, const Intv<double>& b, const Intv<float>& c) { fma(d, a, Intv<float>(__double2float_rd(b.low), __double2float_ru(b.high)), c); }
template<> template<> inline __device__ void Intv<float>::fma(Intv<float>& d, const float& a, const double& b, const Intv<float>& c) { fma(d, a, Intv<float>(__double2float_rd(b), __double2float_ru(b)), c); }


template<typename T>
template<typename Tb>
inline __device__ void Intv<T>::fma(T& d, const T& a, const Tb& b, const T& c)
{
	d = a * b + c;
}

template<> template<> inline __device__ void Intv<double>::add(Intv<double>& d, const Intv<double>& a, const Intv<double>& b) { d.low = __dadd_rd(a.low, b.low);	d.high = __dadd_ru(a.high, b.high); }
template<> template<> inline __device__ void Intv<double>::add(Intv<double>& d, const Intv<double>& a, const Intv<float>& b) { Intv<double> b_d(b.low, b.high);	add(d, a, b_d); }
template<> template<> inline __device__ void Intv<double>::add(Intv<double>& d, const Intv<double>& a, const double& b) { d.low = __dadd_rd(a.low, b); d.high = __dadd_ru(a.high, b); }
template<> template<> inline __device__ void Intv<double>::add(Intv<double>& d, const Intv<double>& a, const float& b) { double b_d = b;	add(d, a, b_d); }
template<> template<> inline __device__ void Intv<double>::add(Intv<double>& d, const double& a, const Intv<double>& b) { d.low = __dadd_rd(a, b.low); d.high = __dadd_ru(a, b.high); }
template<> template<> inline __device__ void Intv<double>::add(Intv<double>& d, const double& a, const Intv<float>& b) { Intv<double> b_d(b.low, b.high);	add(d, a, b_d); }
template<> template<> inline __device__ void Intv<double>::add(Intv<double>& d, const double& a, const double& b) { d.low = __dadd_rd(a, b); d.high = __dadd_ru(a, b); }
template<> template<> inline __device__ void Intv<double>::add(Intv<double>& d, const double& a, const float& b) { double b_d = b; add(d, a, b_d); }
template<> template<> inline __device__ void Intv<float>::add(Intv<float>& d, const Intv<float>& a, const Intv<double>& b) { d.low = __fadd_rd(a.low, __double2float_rd(b.low));	d.high = __fadd_ru(a.high, __double2float_ru(b.high)); }
template<> template<> inline __device__ void Intv<float>::add(Intv<float>& d, const Intv<float>& a, const Intv<float>& b) { d.low = __fadd_rd(a.low, b.low); d.high = __fadd_ru(a.high, b.high); }
template<> template<> inline __device__ void Intv<float>::add(Intv<float>& d, const Intv<float>& a, const double& b) { d.low = __fadd_rd(a.low, __double2float_rd(b));	d.high = __fadd_ru(a.high, __double2float_ru(b)); }
template<> template<> inline __device__ void Intv<float>::add(Intv<float>& d, const Intv<float>& a, const float& b) { d.low = __fadd_rd(a.low, b); d.high = __fadd_ru(a.high, b); }
template<> template<> inline __device__ void Intv<float>::add(Intv<float>& d, const float& a, const Intv<double>& b) { d.low = __fadd_rd(a, __double2float_rd(b.low));	d.high = __fadd_ru(a, __double2float_ru(b.high)); }
template<> template<> inline __device__ void Intv<float>::add(Intv<float>& d, const float& a, const Intv<float>& b) { d.low = __fadd_rd(a, b.low); d.high = __fadd_ru(a, b.high); }
template<> template<> inline __device__ void Intv<float>::add(Intv<float>& d, const float& a, const double& b) { d.low = __fadd_rd(a, __double2float_rd(b));	d.high = __fadd_ru(a, __double2float_ru(b)); }
template<> template<> inline __device__ void Intv<float>::add(Intv<float>& d, const float& a, const float& b) { d.low = __fadd_rd(a, b);	d.high = __fadd_ru(a, b); }



template<typename T>
template<typename Tb>
inline __device__ void Intv<T>::add(T& d, const T& a, const Tb& b)
{
	d = a + b;
}

template <typename T>
template <bool up, typename Ta, typename Tb>
inline __device__ T Intv<T>::fma_dr(const Ta& a, const Tb& b, const T& c)
{
	Intv<T> res;
	fma(res, a, b, c);
	return up ? res.high : res.low;
}

template<> template <bool up> inline __device__ double Intv<double>::mul_dr(const double& ta, const double& tb) { return up ? __dmul_ru(ta, tb) : __dmul_rd(ta, tb); }
template<> template <bool up> inline __device__ float Intv<float>::mul_dr(const float& ta, const float& tb) { return up ? __fmul_ru(ta, tb) : __fmul_rd(ta, tb); }

template<> template <bool up>inline __device__ double Intv<double>::add_dr(const double& ta, const double& tb) { return up ? __dadd_ru(ta, tb) : __dadd_rd(ta, tb); }
template<> template <bool up> inline __device__ float Intv<float>::add_dr(const float& ta, const float& tb) { return up ? __fadd_ru(ta, tb) : __fadd_rd(ta, tb); }

template<> inline __device__ double Intv<double>::max(const double& ta, const double& tb) { return fmax(ta, tb); }
template<> inline __device__ float Intv<float>::max(const float& ta, const float& tb) { return fmaxf(ta, tb); }

template<> inline __device__ double Intv<double>::min(const double& ta, const double& tb) { return fmin(ta, tb); }
template<> inline __device__ float Intv<float>::min(const float& ta, const float& tb) { return fminf(ta, tb); }

template<> template <> inline __device__ double Intv<double>::access_dr<false,double>(const double& input) { return input; }
template<> template <> inline __device__ double Intv<double>::access_dr<true, double>(const double& input) { return input; }
template<> template <> inline __device__ double Intv<double>::access_dr<false, double>(const Intv<double>& input) { return input.low; }
template<> template <> inline __device__ double Intv<double>::access_dr<true, double>(const Intv<double>& input) { return input.high; }
template<> template <> inline __device__ float Intv<float>::access_dr<false, float>(const float& input) { return input; }
template<> template <> inline __device__ float Intv<float>::access_dr<true, float>(const float& input) { return input; }
template<> template <> inline __device__ float Intv<float>::access_dr<false, float>(const Intv<float>& input) { return input.low; }
template<> template <> inline __device__ float Intv<float>::access_dr<true, float>(const Intv<float>& input) { return input.high; }
template<> template <> inline __device__ double Intv<double>::access_dr<false, float>(const float& input) { return input; }
template<> template <> inline __device__ double Intv<double>::access_dr<true, float>(const float& input) { return input; }
template<> template <> inline __device__ double Intv<double>::access_dr<false, float>(const Intv<float>& input) { return input.low; }
template<> template <> inline __device__ double Intv<double>::access_dr<true, float>(const Intv<float>& input) { return input.high; }
template<> template <> inline __device__ float Intv<float>::access_dr<false, double>(const double& input) { return __double2float_rd(input); }
template<> template <> inline __device__ float Intv<float>::access_dr<true, double>(const double& input) { return __double2float_ru(input); }
template<> template <> inline __device__ float Intv<float>::access_dr<false, double>(const Intv<double>& input) { return __double2float_rd(input.low); }
template<> template <> inline __device__ float Intv<float>::access_dr<true, double>(const Intv<double>& input) { return __double2float_ru(input.high); }
template<typename T> template <bool up> inline __device__ T& Intv<T>::access_dr(T& input) { return input; }
template<typename T> template <bool up> inline __device__ T& Intv<T>::access_dr(Intv<T>& input) { return up ? input.high : input.low; }

/*template <typename T> template<typename Ts> inline __device__ void Intv<T>::convert(Intv<T>& dest, const Ts& source);
template <typename T> template<typename Ts> inline __device__ void Intv<T>::convert(T& dest, const Ts& source);*/

/*
namespace IntvOps
{
	template <typename Td>
	struct convert {
		template <typename Ts>
		static __device__ void operator()(Td& dest, const Ts& source);
	};



	template <> template<> inline __device__ void convert<Intv<double>>::operator()(Intv<double>& dest, const Intv<double>& source) { dest = source; }
	template <> template<> inline __device__ void convert<double>::operator()(double& dest, const double& source) { dest = source; }
	template <> template<> inline __device__ void convert<Intv<float>>::operator()(Intv<float>& dest, const Intv<float>& source) { dest = source; }
	template <> template<> inline __device__ void convert<float>::operator()(float& dest, const float& source) { dest = source; }
	template <> template<> inline __device__ void convert<Intv<double>>::operator()(Intv<double>& dest, const Intv<float>& source) { dest = Intv<double>(source.low, source.high); }
	template <> template<> inline __device__ void convert<double>::operator()(double& dest, const float& source) { dest = source; }
	template <> template<> inline __device__ void convert< Intv<float>>::operator()(Intv<float>& dest, const Intv<double>& source) { dest = Intv<float>(__double2float_rd(source.low), __double2float_ru(source.high)); }
	template <> template<> inline __device__ void convert< Intv<float>>::operator()(Intv<float>& dest, const double& source) { dest = Intv<float>(__double2float_rd(source), __double2float_ru(source)); }
}
*/