/*
 *
 *  This source file is part of ELINA (ETH LIbrary for Numerical Analysis).
 *  ELINA is Copyright © 2019 Department of Computer Science, ETH Zurich
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
 *
 *  @file rounding.h
 *  @author Christoph Müller
 *  @brief Provides a templated wrapper around the CUDA rounding intrinsics.
 */

#ifndef __ROUNDING_H_INCLUDED__
#define __ROUNDING_H_INCLUDED__


template <typename T>
__inline__ __device__
T add_ru(const T lhs, const T rhs)
{
    std::cerr << "This function should not be called for types other than float or double" << std::endl;

    exit(1);
};


template <>
__inline__ __device__
float add_ru<float>(const float lhs, const float rhs)
{
    return __fadd_ru(lhs, rhs);
}


template <>
__inline__ __device__
double add_ru<double>(const double lhs, const double rhs)
{
    return __dadd_ru(lhs, rhs);
}


template <typename T>
__inline__ __device__
T add_rd(const T lhs, const T rhs)
{
    std::cerr << "This function should not be called for types other than float or double" << std::endl;

    exit(1);
};


template <>
__inline__ __device__
float add_rd<float>(const float lhs, const float rhs)
{
    return __fadd_rd(lhs, rhs);
}


template <>
__inline__ __device__
double add_rd<double>(const double lhs, const double rhs)
{
    return __dadd_rd(lhs, rhs);
}


template <typename T>
__inline__ __device__
T mul_ru(const T lhs, const T rhs)
{
    std::cerr << "This function should not be called for types other than float or double" << std::endl;

    exit(1);
};


template <>
__inline__ __device__
float mul_ru<float>(const float lhs, const float rhs)
{
    return __fmul_ru(lhs, rhs);
}


template <>
__inline__ __device__
double mul_ru<double>(const double lhs, const double rhs)
{
    return __dmul_ru(lhs, rhs);
}


template <typename T>
__inline__ __device__
T mul_rd(const T lhs, const T rhs)
{
    std::cerr << "This function should not be called for types other than float or double" << std::endl;

    exit(1);
};


template <>
__inline__ __device__
float mul_rd<float>(const float lhs, const float rhs)
{
    return __fmul_rd(lhs, rhs);
}


template <>
__inline__ __device__
double mul_rd<double>(const double lhs, const double rhs)
{
    return __dmul_rd(lhs, rhs);
}


template <typename T>
__inline__ __device__
T div_ru(const T lhs, const T rhs)
{
    std::cerr << "This function should not be called for types other than float or double" << std::endl;

    exit(1);
};


template <>
__inline__ __device__
float div_ru<float>(const float lhs, const float rhs)
{
    return __fdiv_ru(lhs, rhs);
}


template <>
__inline__ __device__
double div_ru<double>(const double lhs, const double rhs)
{
    return __ddiv_ru(lhs, rhs);
}


template <typename T>
__inline__ __device__
T div_rd(const T lhs, const T rhs)
{
    std::cerr << "This function should not be called for types other than float or double" << std::endl;

    exit(1);
};


template <>
__inline__ __device__
float div_rd<float>(const float lhs, const float rhs)
{
    return __fdiv_rd(lhs, rhs);
}


template <>
__inline__ __device__
double div_rd<double>(const double lhs, const double rhs)
{
    return __ddiv_rd(lhs, rhs);
}

#endif //__ROUNDING_H_INCLUDED__
