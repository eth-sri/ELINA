#pragma once


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
