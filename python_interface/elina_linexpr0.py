from elina_linexpr0_h import *


# ====================================================================== #
# I. Memory management and printing
# ====================================================================== #

def elina_linexpr0_alloc(lin_discr, size):
    """
    Allocate a linear expressions with coefficients by default of type ElinaScalar and c_double.
    If sparse representation, corresponding new dimensions are initialized with ELINA_DIM_MAX.
    
    Parameters
    ----------
    lin_discr : c_uint
        Enum of type ElinaLinexprDiscr that defines the representation (sparse or dense).
    size : c_size_t
        Size of the internal array.

    Returns
    -------
    linexpr : ElinaLinexpr0Ptr
        Pointer to the newly allocated ElinaLinexpr0

    """

    linexpr = None
    try:
        elina_linexpr0_alloc_c = elina_auxiliary_api.elina_linexpr0_alloc
        elina_linexpr0_alloc_c.restype = ElinaLinexpr0Ptr
        elina_linexpr0_alloc_c.argtypes = [c_uint, c_size_t]
        linexpr = elina_linexpr0_alloc_c(lin_discr, size)
    except:
        print('Problem with loading/calling "elina_linexpr0_alloc" from "libelinaux.so"')
        print('Make sure you are passing c_uint, c_size_t to the function')

    return linexpr


def elina_linexpr0_realloc(linexpr, size):
    """
    Change the dimensions of the array in an ElinaLinexpr0.
    If new coefficients are added, their type is of type ElinaScalar and c_double.
    If sparse representation, corresponding new dimensions are initialized with ELINA_DIM_MAX.
    
    Parameters
    ----------
    linexpr : ElinaLinexpr0Ptr
        Pointer to the ElinaLinexpr0 that needs to be reallocated.
    size : c_size_t
        Size of the internal array.
        
    Returns
    -------
    None

    """

    try:
        elina_linexpr0_realloc_c = elina_auxiliary_api.elina_linexpr0_realloc
        elina_linexpr0_realloc_c.restype = None
        elina_linexpr0_realloc_c.argtypes = [ElinaLinexpr0Ptr, c_size_t]
        elina_linexpr0_realloc_c(linexpr, size)
    except:
        print('Problem with loading/calling "elina_linexpr0_realloc" from "libelinaux.so"')
        print('Make sure you are passing ElinaLinexpr0Ptr, c_size_t to the function')


def elina_linexpr0_minimize(linexpr):
    """
    Reduce the coefficients of an ElinaLinexpr0 (transform ElinaInterval-s into ElinaScalar-s when possible).
    In case of sparse representation, also remove zero coefficients.
    
    Parameters
    ----------
    linexpr : ElinaLinexpr0Ptr
        Pointer to the ElinaLinexpr0 that needs to be minimized.

    Returns
    -------
    None

    """

    try:
        elina_linexpr0_minimize_c = elina_auxiliary_api.elina_linexpr0_minimize
        elina_linexpr0_minimize_c.restype = None
        elina_linexpr0_minimize_c.argtypes = [ElinaLinexpr0Ptr]
        elina_linexpr0_minimize_c(linexpr)
    except:
        print('Problem with loading/calling "elina_linexpr0_minimize" from "libelinaux.so"')
        print('Make sure you are passing ElinaLinexpr0Ptr to the function')


def elina_linexpr0_free(linexpr):
    """
    Free an ElinaLinexpr0.
    
    Parameters
    ----------
    linexpr : ElinaLinexpr0Ptr
        Pointer to the ElinaLinexpr0 that needs to be freed.

    Returns
    -------
    None

    """

    try:
        elina_linexpr0_free_c = elina_auxiliary_api.elina_linexpr0_free
        elina_linexpr0_free_c.restype = None
        elina_linexpr0_free_c.argtypes = [ElinaLinexpr0Ptr]
        elina_linexpr0_free_c(linexpr)
    except:
        print('Problem with loading/calling "elina_linexpr0_free" from "libelinaux.so"')
        print('Make sure you are passing ElinaLinexpr0Ptr to the function')


def elina_linexpr0_copy(linexpr2):
    """
    Duplicate an ElinaLinexpr0.
    
    Parameters
    ----------
    linexpr2 : ElinaLinexpr0Ptr
        Pointer to the ElinaLinexpr0 that needs to be duplicated.

    Returns
    -------
    linexpr1 : ElinaLinexpr0Ptr
        Duplicate of an ElinaLinexpr0.

    """

    linexpr1 = None
    try:
        elina_linexpr0_copy_c = elina_auxiliary_api.elina_linexpr0_copy
        elina_linexpr0_copy_c.restype = ElinaLinexpr0Ptr
        elina_linexpr0_copy_c.argtypes = [ElinaLinexpr0Ptr]
        linexpr1 = elina_linexpr0_copy_c(linexpr2)
    except:
        print('Problem with loading/calling "elina_linexpr0_copy" from "libelinaux.so"')
        print('Make sure you are passing ElinaLinexpr0Ptr to the function')

    return linexpr1


def elina_linexpr0_fprint(stream, linexpr, name_of_dim):
    """
    Print an ElinaLinexpr to stream.
    
    Parameters
    ----------
    stream : c_void_p
        Stream onto which to print.
    linexpr : ElinaLinexpr0Ptr
        Pointer to the ElinaLinexpr0 that needs to be printed.
    name_of_dim : POINTER(c_char_p)
        Dimension names.

    Returns
    -------
    None

    """

    try:
        elina_linexpr0_fprint_c = elina_auxiliary_api.elina_linexpr0_fprint
        elina_linexpr0_fprint_c.restype = None
        elina_linexpr0_fprint_c.argtypes = [c_void_p, ElinaLinexpr0Ptr, POINTER(c_char_p)]
        elina_linexpr0_fprint_c(stream, linexpr, name_of_dim)
    except:
        print('Problem with loading/calling "elina_linexpr0_fprint" from "libelinaux.so"')
        print('Make sure you are passing c_void_p, ElinaLinexpr0Ptr, POINTER(c_char_p) to the function')


def elina_linexpr0_print(linexpr, name_of_dim):
    """
    Print an ElinaLinexpr0 to stdout.
    
    Parameters
    ----------
    linexpr : ElinaLinexpr0Ptr
        Pointer to the ElinaLinexpr0 that needs to be printed.
    name_of_dim : POINTER(c_char_p)
        Dimension names.

    Returns
    -------
    None

    """

    try:
        elina_linexpr0_print_c = elina_auxiliary_api.elina_linexpr0_print
        elina_linexpr0_print_c.restype = None
        elina_linexpr0_print_c.argtypes = [ElinaLinexpr0Ptr, POINTER(c_char_p)]
        elina_linexpr0_print_c(linexpr, name_of_dim)
    except:
        print('Problem with loading/calling "elina_linexpr0_print" from "libelinaux.so"')
        print('Make sure you are passing ElinaLinexpr0Ptr, POINTER(c_char_p) to the function')


# ====================================================================== #
# II. Tests
# ====================================================================== #

def elina_linexpr0_is_integer(linexpr, intdim):
    """
    Test if the expression depends only on integer variables, assuming that the first intdim dimensions are integer.
    
    Parameters
    ----------
    linexpr : ElinaLinexpr0Ptr
        Pointer to the ElinaLinexpr0 that needs to be tested.
    intdim : c_size_t
        Number of dimensions we assume to be integer.

    Returns
    -------
    result : c_bool
        Result of the test.

    """

    result = None
    try:
        elina_linexpr0_is_integer_c = elina_auxiliary_api.elina_linexpr0_is_integer
        elina_linexpr0_is_integer_c.restype = c_bool
        elina_linexpr0_is_integer_c.argtypes = [ElinaLinexpr0Ptr, c_size_t]
        result = elina_linexpr0_is_integer_c(linexpr, intdim)
    except:
        print('Problem with loading/calling "elina_linexpr0_is_integer" from "libelinaux.so"')
        print('Make sure you are passing ElinaLinexpr0Ptr, c_size_t to the function')

    return result


def elina_linexpr0_is_real(linexpr, intdim):
    """
    Test if the expression depends only on real variables, assuming that the first intdim dimensions are integer.

    Parameters
    ----------
    linexpr : ElinaLinexpr0Ptr
        Pointer to the ElinaLinexpr0 that needs to be tested.
    intdim : c_size_t
        Number of dimensions we assume to be integer.

    Returns
    -------
    result : c_bool
        Result of the test.

    """

    result = None
    try:
        elina_linexpr0_is_real_c = elina_auxiliary_api.elina_linexpr0_is_real
        elina_linexpr0_is_real_c.restype = c_bool
        elina_linexpr0_is_real_c.argtypes = [ElinaLinexpr0Ptr, c_size_t]
        result = elina_linexpr0_is_real_c(linexpr, intdim)
    except:
        print('Problem with loading/calling "elina_linexpr0_is_real" from "libelinaux.so"')
        print('Make sure you are passing ElinaLinexpr0Ptr, c_size_t to the function')

    return result


def elina_linexpr0_type(linexpr):
    """
    Check the type of an ElinaLinexpr0.
    
    Parameters
    ----------
    linexpr : ElinaLinexpr0Ptr
        Pointer to the ElinaLinexpr0 that needs to be checked for its type.

    Returns
    -------
    result : c_uint
        Type of the ElinaLinexpr0 according to ElinaLinexprType.

    """

    result = None
    try:
        elina_linexpr0_type_c = elina_auxiliary_api.elina_linexpr0_type
        elina_linexpr0_type_c.restype = c_uint
        elina_linexpr0_type_c.argtypes = [ElinaLinexpr0Ptr]
        result = elina_linexpr0_type_c(linexpr)
    except:
        print('Problem with loading/calling "elina_linexpr0_type" from "libelinaux.so"')
        print('Make sure you are passing ElinaLinexpr0Ptr to the function')

    return result


def elina_linexpr0_is_linear(linexpr):
    """
    Test if all involved coefficients are scalars.
    
    Parameters
    ----------
    linexpr : ElinaLinexpr0Ptr
        Pointer to the ElinaLinexpr0 that needs to be tested.

    Returns
    -------
    result : c_bool
        Result of the test.

    """

    result = None
    try:
        elina_linexpr0_is_linear_c = elina_auxiliary_api.elina_linexpr0_is_linear
        elina_linexpr0_is_linear_c.restype = c_bool
        elina_linexpr0_is_linear_c.argtypes = [ElinaLinexpr0Ptr]
        result = elina_linexpr0_is_linear_c(linexpr)
    except:
        print('Problem with loading/calling "elina_linexpr0_is_linear" from "libelinaux.so"')
        print('Make sure you are passing ElinaLinexpr0Ptr to the function')

    return result


def elina_linexpr0_is_quasilinear(linexpr):
    """
    Test if all involved coefficients apart from the constants are scalars.

    Parameters
    ----------
    linexpr : ElinaLinexpr0Ptr
        Pointer to the ElinaLinexpr0 that needs to be tested.

    Returns
    -------
    result : c_bool
        Result of the test.

    """

    result = None
    try:
        elina_linexpr0_is_quasilinear_c = elina_auxiliary_api.elina_linexpr0_is_quasilinear
        elina_linexpr0_is_quasilinear_c.restype = c_bool
        elina_linexpr0_is_quasilinear_c.argtypes = [ElinaLinexpr0Ptr]
        result = elina_linexpr0_is_quasilinear_c(linexpr)
    except:
        print('Problem with loading/calling "elina_linexpr0_is_quasilinear" from "libelinaux.so"')
        print('Make sure you are passing ElinaLinexpr0Ptr to the function')

    return result


def elina_linexpr0_array_type(texpr, size):
    """
    Check the type of an ElinaLinexpr0.

    Parameters
    ----------
    texpr : ElinaLinexpr0Array
        Array that needs to be checked for its type.
    size : c_size_t
        Size of the array.

    Returns
    -------
    result : c_uint
        Type of the ElinaLinexpr0Array according to ElinaLinexprType.

    """

    result = None
    try:
        elina_linexpr0_array_type_c = elina_auxiliary_api.elina_linexpr0_array_type
        elina_linexpr0_array_type_c.restype = c_uint
        elina_linexpr0_array_type_c.argtypes = [ElinaLinexpr0Array, c_size_t]
        result = elina_linexpr0_array_type_c(texpr, size)
    except:
        print('Problem with loading/calling "elina_linexpr0_array_type" from "libelinaux.so"')
        print('Make sure you are passing ElinaLinexpr0Array, c_size_t to the function')

    return result


def elina_linexpr0_array_is_linear(texpr, size):
    """
    Test if an ElinaLinexpr0Array is linear.

    Parameters
    ----------
    texpr : ElinaLinexpr0Array
        Array that needs to be tested.
    size : c_size_t
        Size of the array.

    Returns
    -------
    result : c_uint
        Result of the test.
        
    """

    result = None
    try:
        elina_linexpr0_array_is_linear_c = elina_auxiliary_api.elina_linexpr0_array_is_linear
        elina_linexpr0_array_is_linear_c.restype = c_bool
        elina_linexpr0_array_is_linear_c.argtypes = [ElinaLinexpr0Array, c_size_t]
        result = elina_linexpr0_array_is_linear_c(texpr, size)
    except:
        print('Problem with loading/calling "elina_linexpr0_array_is_linear" from "libelinaux.so"')
        print('Make sure you are passing ElinaLinexpr0Array, c_size_t to the function')

    return result


def elina_linexpr0_array_is_quasilinear(texpr, size):
    """
    Test if an ElinaLinexpr0Array is quasilinear.

    Parameters
    ----------
    texpr : ElinaLinexpr0Array
        Array that needs to be tested.
    size : c_size_t
        Size of the array.

    Returns
    -------
    result : c_uint
        Result of the test.

    """

    result = None
    try:
        elina_linexpr0_array_is_quasilinear_c = elina_auxiliary_api.elina_linexpr0_array_is_quasilinear
        elina_linexpr0_array_is_quasilinear_c.restype = c_bool
        elina_linexpr0_array_is_quasilinear_c.argtypes = [ElinaLinexpr0Array, c_size_t]
        result = elina_linexpr0_array_is_quasilinear_c(texpr, size)
    except:
        print('Problem with loading/calling "elina_linexpr0_array_is_quasilinear" from "libelinaux.so"')
        print('Make sure you are passing ElinaLinexpr0Array, c_size_t to the function')

    return result


# ====================================================================== #
# III. Access
# ====================================================================== #

def elina_linexpr0_size(linexpr):
    """
    Return the size of an ElinaLinexpr0.
    
    Parameters
    ----------
    linexpr : ElinaLinexpr0Ptr
        Pointer to the ElinaLinexpr0 that needs to be checked for its size.

    Returns
    -------
    size_linexpr = c_size_t
        Size of the ElinaLinexpr0.

    """

    size_linexpr = None
    try:
        elina_linexpr0_size_c = elina_auxiliary_api.elina_linexpr0_size
        elina_linexpr0_size_c.restype = c_size_t
        elina_linexpr0_size_c.argtypes = [ElinaLinexpr0Ptr]
        size_linexpr = elina_linexpr0_size_c(linexpr)
    except:
        print('Problem with loading/calling "elina_linexpr0_size" from "libelinaux.so"')
        print('Make sure you are passing ElinaLinexpr0Ptr to the function')

    return size_linexpr


def elina_linexpr0_cstref(linexpr):
    """
    Get the constant of an ElinaLinexpr0Ptr.
    
    Parameters
    ----------
    linexpr : ElinaLinexpr0Ptr
        Pointer to the ElinaLinexpr0 that we want to get the constant from.

    Returns
    -------
    coeff : ElinaCoeffPtr
        Pointer to the ElinaCoeff constant of an ElinaLinexpr0.

    """

    coeff = None
    try:
        elina_linexpr0_cstref_c = elina_auxiliary_api.elina_linexpr0_cstref
        elina_linexpr0_cstref_c.restype = ElinaCoeffPtr
        elina_linexpr0_cstref_c.argtypes = [ElinaLinexpr0Ptr]
        coeff = elina_linexpr0_cstref_c(linexpr)
    except:
        print('Problem with loading/calling "elina_linexpr0_cstref" from "libelinaux.so"')
        print('Make sure you are passing ElinaLinexpr0Ptr to the function')

    return coeff


def elina_linexpr0_coeffref(linexpr, dim):
    """
    Get the coefficient associated to a dimension of an ElinaLinexpr0Ptr.
    In case of sparse representation, possibly induce the addition of a new linear term.
    Return NULL if:
    - In case of dense representation, dim >= expr->size.
    - In case of sparse representation, dim == ELINA_DIM_MAX.
    
    Parameters
    ----------
    linexpr : ElinaLinexpr0Ptr
        Pointer to the ElinaLinexpr0 that we want to get the coefficient associated to a dimension from.
    dim : ElinaDim
        Dimension for which we want to get the coefficient.

    Returns
    -------
    coeff : ElinaCoeffPtr
        Pointer to the ElinaCoeff coefficient of an ElinaLinexpr0.


    """

    coeff = None
    try:
        elina_linexpr0_coeffref_c = elina_auxiliary_api.elina_linexpr0_coeffref
        elina_linexpr0_coeffref_c.restype = ElinaCoeffPtr
        elina_linexpr0_coeffref_c.argtypes = [ElinaLinexpr0Ptr, ElinaDim]
        coeff = elina_linexpr0_coeffref_c(linexpr, dim)
    except:
        print('Problem with loading/calling "elina_linexpr0_coeffref" from "libelinaux.so"')
        print('Make sure you are passing ElinaLinexpr0Ptr to the function')

    return coeff


def elina_linexpr0_get_cst(coeff, linexpr):
    """
    Get the constant of an ElinaLinexpr0 and assign it to an ElinaCoeff.
    
    Parameters
    ----------
    coeff : ElinaCoeffPtr
        Destination.
    linexpr : ElinaLinexpr0Ptr
        Source.

    Returns
    -------
    None

    """

    try:
        elina_linexpr0_get_cst_c = elina_auxiliary_api.elina_linexpr0_get_cst
        elina_linexpr0_get_cst_c.restype = None
        elina_linexpr0_get_cst_c.argtypes = [ElinaCoeffPtr, ElinaLinexpr0Ptr]
        elina_linexpr0_get_cst_c(coeff, linexpr)
    except:
        print('Problem with loading/calling "elina_linexpr0_get_cst" from "libelinaux.so"')
        print('Make sure you are passing ElinaCoeffPtr, ElinaLinexpr0Ptr to the function')


def elina_linexpr0_get_coeff(coeff, linexpr, dim):
    """
    Get the coefficient of the dim dimension of an ElinaLinexpr0 and assign it to an ElinaCoeff.

    Parameters
    ----------
    coeff : ElinaCoeffPtr
        Destination.
    linexpr : ElinaLinexpr0Ptr
        Source.
    dim : ElinaDim
        Dimension for which we want to get the coefficient.

    Returns
    -------
    result : c_bool
        Return true in case elina_linexpr0_coeffref returns NULL, false otherwise.

    """

    result = None
    try:
        elina_linexpr0_get_coeff_c = elina_auxiliary_api.elina_linexpr0_get_coeff
        elina_linexpr0_get_coeff_c.restype = c_bool
        elina_linexpr0_get_coeff_c.argtypes = [ElinaCoeffPtr, ElinaLinexpr0Ptr, ElinaDim]
        result = elina_linexpr0_get_coeff_c(coeff, linexpr, dim)
    except:
        print('Problem with loading/calling "elina_linexpr0_get_coeff" from "libelinaux.so"')
        print('Make sure you are passing ElinaCoeffPtr, ElinaLinexpr0Ptr, ElinaDim to the function')

    return result


def elina_linexpr0_set_cst(linexpr, cst):
    """
    Set the constant of an ElinaLinexpr0 by using an ElinaCoeff.
    
    Parameters
    ----------
    linexpr : ElinaLinexpr0Ptr
        Destination.
    cst : ElinaCoeffPtr
        Source.

    Returns
    -------
    None

    """

    try:
        elina_linexpr0_set_cst_c = elina_auxiliary_api.elina_linexpr0_set_cst
        elina_linexpr0_set_cst_c.restype = None
        elina_linexpr0_set_cst_c.argtypes = [ElinaLinexpr0Ptr, ElinaCoeffPtr]
        elina_linexpr0_set_cst_c(linexpr, cst)
    except:
        print('Problem with loading/calling "elina_linexpr0_set_cst" from "libelinaux.so"')
        print('Make sure you are passing ElinaLinexpr0Ptr, ElinaCoeffPtr to the function')


def elina_linexpr0_set_cst_scalar(linexpr, scalar):
    """
    Set the constant of an ElinaLinexpr0 by using an ElinaScalar.

    Parameters
    ----------
    linexpr : ElinaLinexpr0Ptr
        Destination.
    scalar : ElinaScalarPtr
        Source.

    Returns
    -------
    None

    """

    try:
        elina_linexpr0_set_cst_scalar_c = elina_auxiliary_api.elina_linexpr0_set_cst_scalar
        elina_linexpr0_set_cst_scalar_c.restype = None
        elina_linexpr0_set_cst_scalar_c.argtypes = [ElinaLinexpr0Ptr, ElinaScalarPtr]
        elina_linexpr0_set_cst_scalar_c(linexpr, scalar)
    except:
        print('Problem with loading/calling "elina_linexpr0_set_cst_scalar" from "libelinaux.so"')
        print('Make sure you are passing ElinaLinexpr0Ptr, ElinaScalarPtr to the function')


def elina_linexpr0_set_cst_scalar_int(linexpr, num):
    """
    Set the constant of an ElinaLinexpr0 by using a c_int.

    Parameters
    ----------
    linexpr : ElinaLinexpr0Ptr
        Destination.
    num : c_int
        Source.

    Returns
    -------
    None

    """

    try:
        elina_linexpr0_set_cst_scalar_int_c = elina_auxiliary_api.elina_linexpr0_set_cst_scalar_int
        elina_linexpr0_set_cst_scalar_int_c.restype = None
        elina_linexpr0_set_cst_scalar_int_c.argtypes = [ElinaLinexpr0Ptr, c_int]
        elina_linexpr0_set_cst_scalar_int_c(linexpr, num)
    except:
        print('Problem with loading/calling "elina_linexpr0_set_cst_scalar_int" from "libelinaux.so"')
        print('Make sure you are passing ElinaLinexpr0Ptr, c_int to the function')


def elina_linexpr0_set_cst_scalar_frac(linexpr, num, den):
    """
    Set the constant of an ElinaLinexpr0 by using an integer fraction.
    
    Parameters
    ----------
    linexpr : ElinaLinexpr0Ptr
        Destination.
    num : c_long
        Source.
    den : c_ulong
        Source.

    Returns
    -------
    None

    """

    try:
        elina_linexpr0_set_cst_scalar_frac_c = elina_auxiliary_api.elina_linexpr0_set_cst_scalar_frac
        elina_linexpr0_set_cst_scalar_frac_c.restype = None
        elina_linexpr0_set_cst_scalar_frac_c.argtypes = [ElinaLinexpr0Ptr, c_long, c_ulong]
        elina_linexpr0_set_cst_scalar_frac_c(linexpr, num, den)
    except:
        print('Problem with loading/calling "elina_linexpr0_set_cst_scalar_frac" from "libelinaux.so"')
        print('Make sure you are passing ElinaLinexpr0Ptr, c_long, c_ulong to the function')


def elina_linexpr0_set_cst_scalar_double(linexpr, num):
    """
    Set the constant of an ElinaLinexpr0 by using a c_double.

    Parameters
    ----------
    linexpr : ElinaLinexpr0Ptr
        Destination.
    num : c_double
        Source.

    Returns
    -------
    None

    """

    try:
        elina_linexpr0_set_cst_scalar_double_c = elina_auxiliary_api.elina_linexpr0_set_cst_scalar_double
        elina_linexpr0_set_cst_scalar_double_c.restype = None
        elina_linexpr0_set_cst_scalar_double_c.argtypes = [ElinaLinexpr0Ptr, c_double]
        elina_linexpr0_set_cst_scalar_double_c(linexpr, num)
    except:
        print('Problem with loading/calling "elina_linexpr0_set_cst_scalar_double" from "libelinaux.so"')
        print('Make sure you are passing ElinaLinexpr0Ptr, c_double to the function')


def elina_linexpr0_set_cst_interval(linexpr, interval):
    """
    Set the constant of an ElinaLinexpr0 by using an ElinaInterval.

    Parameters
    ----------
    linexpr : ElinaLinexpr0Ptr
        Destination.
    interval : ElinaIntervalPtr
        Source.

    Returns
    -------
    None

    """
    try:
        elina_linexpr0_set_cst_interval_c = elina_auxiliary_api.elina_linexpr0_set_cst_interval
        elina_linexpr0_set_cst_interval_c.restype = None
        elina_linexpr0_set_cst_interval_c.argtypes = [ElinaLinexpr0Ptr, ElinaIntervalPtr]
        elina_linexpr0_set_cst_interval_c(linexpr, interval)
    except:
        print('Problem with loading/calling "elina_linexpr0_set_cst_interval" from "libelinaux.so"')
        print('Make sure you are passing ElinaLinexpr0Ptr, ElinaIntervalPtr to the function')


def elina_linexpr0_set_cst_interval_scalar(linexpr, inf, sup):
    """
    Set the constant of an ElinaLinexpr0 by using an interval of two ElinaScalar-s.

    Parameters
    ----------
    linexpr : ElinaLinexpr0Ptr
        Destination.
    inf : ElinaScalarPtr
        Source.
    sup : ElinaScalarPtr
        Source.

    Returns
    -------
    None

    """

    try:
        elina_linexpr0_set_cst_interval_scalar_c = elina_auxiliary_api.elina_linexpr0_set_cst_interval_scalar
        elina_linexpr0_set_cst_interval_scalar_c.restype = None
        elina_linexpr0_set_cst_interval_scalar_c.argtypes = [ElinaLinexpr0Ptr, ElinaScalarPtr, ElinaScalarPtr]
        elina_linexpr0_set_cst_interval_scalar_c(linexpr, inf, sup)
    except:
        print('Problem with loading/calling "elina_linexpr0_set_cst_interval_scalar" from "libelinaux.so"')
        print('Make sure you are passing ElinaLinexpr0Ptr, ElinaScalarPtr, ElinaScalarPtr to the function')


def elina_linexpr0_set_cst_interval_int(linexpr, inf, sup):
    """
    Set the constant of an ElinaLinexpr0 by using an interval defined by two c_int-s.

    Parameters
    ----------
    linexpr : ElinaLinexpr0Ptr
        Destination.
    inf : c_int
        Source.
    sup : c_int
        Source.

    Returns
    -------
    None

    """

    try:
        elina_linexpr0_set_cst_interval_int_c = elina_auxiliary_api.elina_linexpr0_set_cst_interval_int
        elina_linexpr0_set_cst_interval_int_c.restype = None
        elina_linexpr0_set_cst_interval_int_c.argtypes = [ElinaLinexpr0Ptr, c_int, c_int]
        elina_linexpr0_set_cst_interval_int_c(linexpr, inf, sup)
    except:
        print('Problem with loading/calling "elina_linexpr0_set_cst_interval_int" from "libelinaux.so"')
        print('Make sure you are passing ElinaLinexpr0Ptr, c_int, c_int to the function')


def elina_linexpr0_set_cst_interval_frac(linexpr, numinf, deninf, numsup, densup):
    """
    Set the constant of an ElinaLinexpr0 by using an interval defined by two integer fractions.

    Parameters
    ----------
    linexpr : ElinaLinexpr0Ptr
        Destination.
    numinf : c_int
        Source.
    numsup : c_int
        Source.
    deninf : c_int
        Source.
    densup : c_int
        Source.

    Returns
    -------
    None

    """

    try:
        elina_linexpr0_set_cst_interval_frac_c = elina_auxiliary_api.elina_linexpr0_set_cst_interval_frac
        elina_linexpr0_set_cst_interval_frac_c.restype = None
        elina_linexpr0_set_cst_interval_frac_c.argtypes = [ElinaLinexpr0Ptr, c_long, c_ulong, c_long, c_ulong]
        elina_linexpr0_set_cst_interval_frac_c(linexpr, numinf, deninf, numsup, densup)
    except:
        print('Problem with loading/calling "elina_linexpr0_set_cst_interval_frac" from "libelinaux.so"')
        print('Make sure you are passing ElinaLinexpr0Ptr, c_long, c_ulong, c_long, c_ulong to the function')


def elina_linexpr0_set_cst_interval_double(linexpr, inf, sup):
    """
    Set the constant of an ElinaLinexpr0 by using an interval defined by two c_double-s.

    Parameters
    ----------
    linexpr : ElinaLinexpr0Ptr
        Destination.
    inf : c_double
        Source.
    sup : c_double
        Source.

    Returns
    -------
    None

    """

    try:
        elina_linexpr0_set_cst_interval_double_c = elina_auxiliary_api.elina_linexpr0_set_cst_interval_double
        elina_linexpr0_set_cst_interval_double_c.restype = None
        elina_linexpr0_set_cst_interval_double_c.argtypes = [ElinaLinexpr0Ptr, c_double, c_double]
        elina_linexpr0_set_cst_interval_double_c(linexpr, inf, sup)
    except:
        print('Problem with loading/calling "elina_linexpr0_set_cst_interval_double" from "libelinaux.so"')
        print('Make sure you are passing ElinaLinexpr0Ptr, c_double, c_double to the function')


def elina_linexpr0_set_coeff(linexpr, dim, coeff):
    """
    Set the dim dimension coefficient of an ElinaLinexpr0 by using an ElinaCoeff.

    Parameters
    ----------
    linexpr : ElinaLinexpr0Ptr
        Destination.
    dim : ElinaDim
        Destination.
    coeff : ElinaCoeffPtr
        Source.

    Returns
    -------
    None

    """

    try:
        elina_linexpr0_set_coeff_c = elina_auxiliary_api.elina_linexpr0_set_coeff
        elina_linexpr0_set_coeff_c.restype = None
        elina_linexpr0_set_coeff_c.argtypes = [ElinaLinexpr0Ptr, ElinaDim, ElinaCoeffPtr]
        elina_linexpr0_set_coeff_c(linexpr, dim, coeff)
    except:
        print('Problem with loading/calling "elina_linexpr0_set_coeff" from "libelinaux.so"')
        print('Make sure you are passing ElinaLinexpr0Ptr, ElinaDim, ElinaCoeffPtr to the function')


def elina_linexpr0_set_coeff_scalar(linexpr, dim, scalar):
    """
    Set the dim dimension coefficient of an ElinaLinexpr0 by using an ElinaScalar.

    Parameters
    ----------
    linexpr : ElinaLinexpr0Ptr
        Destination.
    dim : ElinaDim
        Destination.
    scalar : ElinaScalar
        Source.

    Returns
    -------
    None

    """

    try:
        elina_linexpr0_set_coeff_scalar_c = elina_auxiliary_api.elina_linexpr0_set_coeff_scalar
        elina_linexpr0_set_coeff_scalar_c.restype = None
        elina_linexpr0_set_coeff_scalar_c.argtypes = [ElinaLinexpr0Ptr, ElinaDim, ElinaScalarPtr]
        elina_linexpr0_set_coeff_scalar_c(linexpr, dim, scalar)
    except:
        print('Problem with loading/calling "elina_linexpr0_set_coeff_scalar" from "libelinaux.so"')
        print('Make sure you are passing ElinaLinexpr0Ptr, ElinaDim, ElinaScalarPtr to the function')


def elina_linexpr0_set_coeff_scalar_int(linexpr, dim, num):
    """
    Set the dim dimension coefficient of an ElinaLinexpr0 by using a c_int.

    Parameters
    ----------
    linexpr : ElinaLinexpr0Ptr
        Destination.
    dim : ElinaDim
        Destination.
    num : c_int
        Source.

    Returns
    -------
    None

    """

    try:
        elina_linexpr0_set_coeff_scalar_int_c = elina_auxiliary_api.elina_linexpr0_set_coeff_scalar_int
        elina_linexpr0_set_coeff_scalar_int_c.restype = None
        elina_linexpr0_set_coeff_scalar_int_c.argtypes = [ElinaLinexpr0Ptr, ElinaDim, c_int]
        elina_linexpr0_set_coeff_scalar_int_c(linexpr, dim, num)
    except:
        print('Problem with loading/calling "elina_linexpr0_set_coeff_scalar_int" from "libelinaux.so"')
        print('Make sure you are passing ElinaLinexpr0Ptr, ElinaDim, c_int to the function')


def elina_linexpr0_set_coeff_scalar_frac(linexpr, dim, num, den):
    """
    Set the dim dimension coefficient of an ElinaLinexpr0 by using an integer fraction.

    Parameters
    ----------
    linexpr : ElinaLinexpr0Ptr
        Destination.
    dim : ElinaDim
        Destination.
    num : c_long
        Source.
    den : c_long
        Source.

    Returns
    -------
    None

    """

    try:
        elina_linexpr0_set_coeff_scalar_frac_c = elina_auxiliary_api.elina_linexpr0_set_coeff_scalar_frac
        elina_linexpr0_set_coeff_scalar_frac_c.restype = None
        elina_linexpr0_set_coeff_scalar_frac_c.argtypes = [ElinaLinexpr0Ptr, ElinaDim, c_long, c_ulong]
        elina_linexpr0_set_coeff_scalar_frac_c(linexpr, dim, num, den)
    except:
        print('Problem with loading/calling "elina_linexpr0_set_coeff_scalar_frac" from "libelinaux.so"')
        print('Make sure you are passing ElinaLinexpr0Ptr, ElinaDim, c_long, c_ulong to the function')


def elina_linexpr0_set_coeff_scalar_double(linexpr, dim, num):
    """
    Set the dim dimension coefficient of an ElinaLinexpr0 by using a c_double.

    Parameters
    ----------
    linexpr : ElinaLinexpr0Ptr
        Destination.
    dim : ElinaDim
        Destination.
    num : c_double
        Source.

    Returns
    -------
    None

    """

    try:
        elina_linexpr0_set_coeff_scalar_double_c = elina_auxiliary_api.elina_linexpr0_set_coeff_scalar_double
        elina_linexpr0_set_coeff_scalar_double_c.restype = None
        elina_linexpr0_set_coeff_scalar_double_c.argtypes = [ElinaLinexpr0Ptr, ElinaDim, c_double]
        elina_linexpr0_set_coeff_scalar_double_c(linexpr, dim, num)
    except:
        print('Problem with loading/calling "elina_linexpr0_set_coeff_scalar_double" from "libelinaux.so"')
        print('Make sure you are passing ElinaLinexpr0Ptr, ElinaDim, c_double to the function')


def elina_linexpr0_set_coeff_interval(linexpr, dim, interval):
    """
    Set the dim dimension coefficient of an ElinaLinexpr0 by using an ElinaInterval.

    Parameters
    ----------
    linexpr : ElinaLinexpr0Ptr
        Destination.
    dim : ElinaDim
        Destination.
    interval : ElinaIntervalPtr
        Source.

    Returns
    -------
    None

    """

    try:
        elina_linexpr0_set_coeff_interval_c = elina_auxiliary_api.elina_linexpr0_set_coeff_interval
        elina_linexpr0_set_coeff_interval_c.restype = None
        elina_linexpr0_set_coeff_interval_c.argtypes = [ElinaLinexpr0Ptr, ElinaDim, ElinaIntervalPtr]
        elina_linexpr0_set_coeff_interval_c(linexpr, dim, interval)
    except:
        print('Problem with loading/calling "elina_linexpr0_set_coeff_interval" from "libelinaux.so"')
        print('Make sure you are passing ElinaLinexpr0Ptr, ElinaDim, ElinaIntervalPtr to the function')


def elina_linexpr0_set_coeff_interval_scalar(linexpr, dim, inf, sup):
    """
    Set the dim dimension coefficient of an ElinaLinexpr0 by using an interval defined by two ElinaScalar-s.

    Parameters
    ----------
    linexpr : ElinaLinexpr0Ptr
        Destination.
    dim : ElinaDim
        Destination.
    inf : ElinaScalar
        Source.
    sup : ElinaScalar
        Source.

    Returns
    -------
    None

    """

    try:
        elina_linexpr0_set_coeff_interval_scalar_c = elina_auxiliary_api.elina_linexpr0_set_coeff_interval_scalar
        elina_linexpr0_set_coeff_interval_scalar_c.restype = None
        elina_linexpr0_set_coeff_interval_scalar_c.argtypes = [ElinaLinexpr0Ptr, ElinaDim, ElinaScalarPtr, ElinaScalarPtr]
        elina_linexpr0_set_coeff_interval_scalar_c(linexpr, dim, inf, sup)
    except:
        print('Problem with loading/calling "elina_linexpr0_set_coeff_interval_scalar" from "libelinaux.so"')
        print('Make sure you are passing ElinaLinexpr0Ptr, ElinaDim, ElinaScalarPtr, ElinaScalarPtr to the function')


def elina_linexpr0_set_coeff_interval_int(linexpr, dim, inf, sup):
    """
    Set the dim dimension coefficient of an ElinaLinexpr0 by using an interval defined by two c_int-s.

    Parameters
    ----------
    linexpr : ElinaLinexpr0Ptr
        Destination.
    dim : ElinaDim
        Destination.
    inf : c_int
        Source.
    sup : c_int
        Source.

    Returns
    -------
    None

    """

    try:
        elina_linexpr0_set_coeff_interval_int_c = elina_auxiliary_api.elina_linexpr0_set_coeff_interval_int
        elina_linexpr0_set_coeff_interval_int_c.restype = None
        elina_linexpr0_set_coeff_interval_int_c.argtypes = [ElinaLinexpr0Ptr, ElinaDim, c_int, c_int]
        elina_linexpr0_set_coeff_interval_int_c(linexpr, dim, inf, sup)
    except:
        print('Problem with loading/calling "elina_linexpr0_set_coeff_interval_int" from "libelinaux.so"')
        print('Make sure you are passing ElinaLinexpr0Ptr, ElinaDim, c_int, c_int to the function')


def elina_linexpr0_set_coeff_interval_frac(linexpr, dim, numinf, numsup, deninf, densup):
    """
    Set the dim dimension coefficient of an ElinaLinexpr0 by using an interval defined by two integer fractions.

    Parameters
    ----------
    linexpr : ElinaLinexpr0Ptr
        Destination.
    dim : ElinaDim
        Destination.
    numinf : c_long
        Source.
    numsup : c_ulong
        Source.
    deninf : c_long
        Source.
    densup : c_ulong
        Source.

    Returns
    -------
    None

    """

    try:
        elina_linexpr0_set_coeff_interval_frac_c = elina_auxiliary_api.elina_linexpr0_set_coeff_interval_frac
        elina_linexpr0_set_coeff_interval_frac_c.restype = None
        elina_linexpr0_set_coeff_interval_frac_c.argtypes = [ElinaLinexpr0Ptr, ElinaDim, c_long, c_ulong, c_long, c_ulong]
        elina_linexpr0_set_coeff_interval_frac_c(linexpr, dim, numinf, numsup, deninf, densup)
    except:
        print('Problem with loading/calling "elina_linexpr0_set_coeff_interval_frac" from "libelinaux.so"')
        print('Make sure you are passing ElinaLinexpr0Ptr, ElinaDim, c_long, c_ulong, c_long, c_ulong to the function')


def elina_linexpr0_set_coeff_interval_double(linexpr, dim, inf, sup):
    """
    Set the dim dimension coefficient of an ElinaLinexpr0 by using an interval defined by two c_double-s.

    Parameters
    ----------
    linexpr : ElinaLinexpr0Ptr
        Destination.
    dim : ElinaDim
        Destination.
    inf : c_double
        Source.
    sup : c_double
        Source.

    Returns
    -------
    None

    """

    try:
        elina_linexpr0_set_coeff_interval_double_c = elina_auxiliary_api.elina_linexpr0_set_coeff_interval_double
        elina_linexpr0_set_coeff_interval_double_c.restype = None
        elina_linexpr0_set_coeff_interval_double_c.argtypes = [ElinaLinexpr0Ptr, ElinaDim, c_double, c_double]
        elina_linexpr0_set_coeff_interval_double_c(linexpr, dim, inf, sup)
    except:
        print('Problem with loading/calling "elina_linexpr0_set_coeff_interval_double" from "libelinaux.so"')
        print('Make sure you are passing ElinaLinexpr0Ptr, ElinaDim, c_double, c_double to the function')


def elina_linexpr0_set_list_generic():
    # TODO not implemented

    return -1


def elina_linexpr0_set_list():
    # TODO not implemented

    return -1


# ====================================================================== #
# IV. Change of dimensions and permutations
# ====================================================================== #

def elina_linexpr0_add_dimensions_with(linexpr, dimchange):
    """
    Add dimensions to an ElinaLinexpr0 by following the semantics of an ElinaDimchange.
    
    Parameters
    ----------
    linexpr : ElinaLinexpr0Ptr
        Pointer to the ElinaLinexpr0 to which we want to add dimensions.
    dimchange : ElinaDimchangePtr
        Pointer to the ElinaDimchange which semantics we need to follow.

    Returns
    -------
    None

    """

    try:
        elina_linexpr0_add_dimensions_with_c = elina_auxiliary_api.elina_linexpr0_add_dimensions_with
        elina_linexpr0_add_dimensions_with_c.restype = None
        elina_linexpr0_add_dimensions_with_c.argtypes = [ElinaLinexpr0Ptr, ElinaDimchangePtr]
        elina_linexpr0_add_dimensions_with_c(linexpr, dimchange)
    except:
        print('Problem with loading/calling "elina_linexpr0_add_dimensions_with" from "libelinaux.so"')
        print('Make sure you are passing ElinaLinexpr0Ptr, ElinaDimchangePtr to the function')


def elina_linexpr0_add_dimensions(linexpr2, dimchange):
    """
    Add dimensions to an ElinaLinexpr0 by following the semantics of an ElinaDimchange.

    Parameters
    ----------
    linexpr2 : ElinaLinexpr0Ptr
        Pointer to the ElinaLinexpr0 to which we want to add dimensions.
    dimchange : ElinaDimchangePtr
        Pointer to the ElinaDimchange which semantics we need to follow.

    Returns
    -------
    linexpr1 : ElinaLinexpr0Ptr
        Pointer to the newly created ElinaLinexpr0 with added dimensions.

    """

    linexpr1 = None
    try:
        elina_linexpr0_add_dimensions_c = elina_auxiliary_api.elina_linexpr0_add_dimensions
        elina_linexpr0_add_dimensions_c.restype = ElinaLinexpr0Ptr
        elina_linexpr0_add_dimensions_c.argtypes = [ElinaLinexpr0Ptr, ElinaDimchangePtr]
        linexpr1 = elina_linexpr0_add_dimensions_c(linexpr2, dimchange)
    except:
        print('Problem with loading/calling "elina_linexpr0_add_dimensions" from "libelinaux.so"')
        print('Make sure you are passing ElinaLinexpr0Ptr, ElinaDimchangePtr to the function')

    return linexpr1


def elina_linexpr0_permute_dimensions_with(linexpr, perm):
    """
    Apply given permutation to the dimensions of an ElinaLinexpr0.
    If dense, the size of the permutation should be the same as the size of the ElinaLinexpr0.
    If sparse, the dimensions present in the expression should just be less than the size of the permutation.
    
    Parameters
    ----------
    linexpr : ElinaLinexpr0Ptr
        Pointer to the ElinaLinexpr0 which dimensions we want to permute.
    perm : ElinaDimpermPtr
        Pointer to the ElinaDimpermPtr which describes the permutation.

    Returns
    -------
    None
    
    """

    try:
        elina_linexpr0_permute_dimensions_with_c = elina_auxiliary_api.elina_linexpr0_permute_dimensions_with
        elina_linexpr0_permute_dimensions_with_c.restype = None
        elina_linexpr0_permute_dimensions_with_c.argtypes = [ElinaLinexpr0Ptr, ElinaDimpermPtr]
        elina_linexpr0_permute_dimensions_with_c(linexpr, perm)
    except:
        print('Problem with loading/calling "elina_linexpr0_permute_dimensions_with" from "libelinaux.so"')
        print('Make sure you are passing ElinaLinexpr0Ptr, ElinaDimpermPtr to the function')


def elina_linexpr0_permute_dimensions(linexpr2, perm):
    """
    Apply given permutation to the dimensions of an ElinaLinexpr0.
    If dense, the size of the permutation should be the same as the size of the ElinaLinexpr0.
    If sparse, the dimensions present in the expression should just be less than the size of the permutation.

    Parameters
    ----------
    linexpr2 : ElinaLinexpr0Ptr
        Pointer to the ElinaLinexpr0 which dimensions we want to permute.
    perm : ElinaDimpermPtr
        Pointer to the ElinaDimpermPtr which describes the permutation.

    Returns
    -------
    linexpr1 : ElinaLinexpr0Ptr
        Pointer to the newly created ElinaLinexpr0 with permuted dimensions.

    """

    linexpr1 = None
    try:
        elina_linexpr0_permute_dimensions_c = elina_auxiliary_api.elina_linexpr0_permute_dimensions
        elina_linexpr0_permute_dimensions_c.restype = ElinaLinexpr0Ptr
        elina_linexpr0_permute_dimensions_c.argtypes = [ElinaLinexpr0Ptr, ElinaDimpermPtr]
        linexpr1 = elina_linexpr0_permute_dimensions_c(linexpr2, perm)
    except:
        print('Problem with loading/calling "elina_linexpr0_permute_dimensions" from "libelinaux.so"')
        print('Make sure you are passing ElinaLinexpr0Ptr, ElinaDimpermPtr to the function')

    return linexpr1


# ====================================================================== #
# V. Hashing, comparison
# ====================================================================== #

def elina_linexpr0_hash(linexpr):
    """
    Calculate the hash code of an ElinaLinexpr0.
    
    Parameters
    ----------
    linexpr : ElinaLinexpr0Ptr
        Pointer to the ElinaLinexpr0 that needs to be hashed.

    Returns
    -------
    result : c_long
        Resulting hash.

    """

    result = None
    try:
        elina_linexpr0_hash_c = elina_auxiliary_api.elina_linexpr0_hash
        elina_linexpr0_hash_c.restype = c_long
        elina_linexpr0_hash_c.argtypes = [ElinaLinexpr0Ptr]
        result = elina_linexpr0_hash_c(linexpr)
    except:
        print('Problem with loading/calling "elina_linexpr0_hash" from "libelinaux.so"')
        print('Make sure you are passing ElinaLinexpr0Ptr to the function')

    return result


def elina_linexpr0_equal(linexpr1, linexpr2):
    """
    Test if an ElinaCoeff is equal to an integer.

    Parameters
    ----------
    linexpr1 : ElinaLinexpr0Ptr
        Pointer to the ElinaLinexpr0 that needs to be tested for equality.
    linexpr2 : ElinaLinexpr0Ptr
        Pointer to the ElinaLinexpr0 that needs to be tested for equality.

    Returns
    -------
    result : c_bool
        Result of the equality test.

    """

    result = None
    try:
        elina_linexpr0_equal_c = elina_auxiliary_api.elina_linexpr0_equal
        elina_linexpr0_equal_c.restype = c_bool
        elina_linexpr0_equal_c.argtypes = [ElinaLinexpr0Ptr, ElinaLinexpr0Ptr]
        result = elina_linexpr0_equal_c(linexpr1, linexpr2)
    except:
        print('Problem with loading/calling "elina_linexpr0_equal" from "libelinaux.so"')
        print('Make sure you are passing ElinaLinexpr0Ptr, ElinaLinexpr0Ptr to the function')

    return result


def elina_linexpr0_compare(linexpr1, linexpr2):
    """
    Compare two ElinaLinexpr0 by lexicographic ordering, terminating by constant coefficients.
    
    Parameters
    ----------
    linexpr1 : ElinaLinexpr0Ptr
        Pointer to the ElinaLinexpr0 that needs to be compared.
    linexpr2 : ElinaLinexpr0Ptr
        Pointer to the ElinaLinexpr0 that needs to be compared.

    Returns
    -------
    result : c_int
        Result of the comparison.
    """

    result = None
    try:
        elina_linexpr0_compare_c = elina_auxiliary_api.elina_linexpr0_compare
        elina_linexpr0_compare_c.restype = c_int
        elina_linexpr0_compare_c.argtypes = [ElinaLinexpr0Ptr, ElinaLinexpr0Ptr]
        result = elina_linexpr0_compare_c(linexpr1, linexpr2)
    except:
        print('Problem with loading/calling "elina_linexpr0_compare" from "libelinaux.so"')
        print('Make sure you are passing ElinaLinexpr0Ptr, ElinaLinexpr0Ptr to the function')

    return result


# ====================================================================== #
# Vb. Array of expressions
# ====================================================================== #


def elina_linexpr0_array_free(texpr, size):
    """
    Free an ElinaLinexpr0Array.
    
    Parameters
    ----------
    texpr : ElinaLinexpr0Array
        ElinaLinexpr0Array that needs to be freed.
    size : c_size_t
        Size of the ElinaLinexpr0Array that needs to be freed.

    Returns
    -------
    None
    
    """

    try:
        elina_linexpr0_array_free_c = elina_auxiliary_api.elina_linexpr0_array_free
        elina_linexpr0_array_free_c.restype = None
        elina_linexpr0_array_free_c.argtypes = [ElinaLinexpr0Array, c_size_t]
        elina_linexpr0_array_free_c(texpr, size)
    except:
        print('Problem with loading/calling "elina_linexpr0_array_free" from "libelinaux.so"')
        print('Make sure you are passing ElinaLinexpr0Array, c_size_t to the function')
