from elina_interval_h import *

# ====================================================================== #
# Basics
# ====================================================================== #


def elina_interval_alloc():
    """
    Allocate a new ElinaInterval, using DOUBLE as default type for the scalars inside.
    
    Returns
    -------
    interval : ElinaIntervalPtr
        Pointer to the newly allocated ElinaInterval.
    
    """

    interval = None
    try:
        elina_interval_alloc_c = elina_auxiliary_api.elina_interval_alloc
        elina_interval_alloc_c.restype = ElinaIntervalPtr
        interval = elina_interval_alloc_c()
    except:
        print('Problem with loading/calling "elina_interval_alloc" from "libelinaux.so"')

    return interval


def elina_interval_reinit(interval, elina_scalar_discr):
    """
    Reinitialise a given ElinaInterval, according to the provided scalar type.
    
    Parameters
    ----------
    interval : ElinaIntervalPtr
        Pointer to the ElinaInterval that needs reinitialisation.
    elina_scalar_discr : c_uint
        Enum of type ElinaScalarDiscr that defines the core of the ElinaScalars (0 = double, 1 = mpq, 2 = mpfr) used in the interval.

    Returns
    -------
    None
    
    """

    try:
        elina_interval_reinit_c = elina_auxiliary_api.elina_interval_reinit
        elina_interval_reinit_c.restype = None
        elina_interval_reinit_c.argtypes = [ElinaIntervalPtr, c_uint]
        elina_interval_reinit_c(interval, elina_scalar_discr)
    except:
        print('Problem with loading/calling "elina_interval_reinit" from "libelinaux.so"')
        print('Make sure you are passing ElinaIntervalPtr and c_uint to the function')


def elina_interval_free(interval):
    """
    Free an ElinaInterval.
    
    Parameters
    ----------
    interval : ElinaIntervalPtr
        Pointer to the ElinaInterval that needs to be freed.

    Returns
    -------
    None
    
    """

    try:
        elina_interval_free_c = elina_auxiliary_api.elina_interval_free
        elina_interval_free_c.restype = None
        elina_interval_free_c.argtypes = [ElinaIntervalPtr]
        elina_interval_free_c(interval)
    except:
        print('Problem with loading/calling "elina_interval_free" from "libelinaux.so"')
        print('Make sure you are passing ElinaIntervalPtr to the function')


def elina_interval_fprint(stream, interval):
    """
    Print an ElinaInterval onto a given stream.
    
    Parameters
    ----------
    stream : c_void_p
        Stream on which to print.
    interval : ElinaIntervalPtr
        Pointer to the ElinaInterval that needs to be printed.

    Returns
    -------
    None

    """

    try:
        elina_interval_fprint_c = elina_auxiliary_api.elina_interval_fprint
        elina_interval_fprint_c.restype = None
        elina_interval_fprint_c.argtypes = [c_void_p, ElinaIntervalPtr]
        elina_interval_fprint_c(stream, interval)
    except:
        print('Problem with loading/calling "elina_interval_fprint" from "libelinaux.so"')
        print('Make sure you are passing c_void_p and ElinaIntervalPtr to the function')


# ====================================================================== #
# Assignments
# ====================================================================== #


def elina_interval_set(interval1, interval2):
    """
    Set the value of one ElinaInterval to the value of another ElinaInterval.
    
    Parameters
    ----------
    interval1 : ElinaIntervalPtr
        Destination.
    interval2 : ElinaIntervalPtr
        Source.

    Returns
    -------
    None

    """

    try:
        elina_interval_set_c = elina_auxiliary_api.elina_interval_set
        elina_interval_set_c.restype = None
        elina_interval_set_c.argtypes = [ElinaIntervalPtr, ElinaIntervalPtr]
        elina_interval_set_c(interval1, interval2)
    except:
        print('Problem with loading/calling "elina_interval_set" from "libelinaux.so"')
        print('Make sure you are passing ElinaIntervalPtr and ElinaIntervalPtr to the function')


def elina_interval_set_scalar(interval, inf, sup):
    """
    Set the value of an ElinaInterval by using two ElinaScalar-s.
    
    Parameters
    ----------
    interval : ElinaIntervalPtr
        Destination.
    inf : ElinaScalarPtr
        Lower bound.
    sup : ElinaScalarPtr
        Upper bound.

    Returns
    -------
    None
    
    """

    try:
        elina_interval_set_scalar_c = elina_auxiliary_api.elina_interval_set_scalar
        elina_interval_set_scalar_c.restype = None
        elina_interval_set_scalar_c.argtypes = [ElinaIntervalPtr, ElinaScalarPtr, ElinaScalarPtr]
        elina_interval_set_scalar_c(interval, inf, sup)
    except:
        print('Problem with loading/calling "elina_interval_set_scalar" from "libelinaux.so"')
        print('Make sure you are passing ElinaIntervalPtr, ElinaScalarPtr and ElinaScalarPtr to the function')


def elina_interval_set_mpq(interval, inf, sup):
    """
    Set the value of an ElinaInterval by using two Mpq-s.

    Parameters
    ----------
    interval : ElinaIntervalPtr
        Destination.
    inf : Mpq_t
        Lower bound.
    sup : Mpq_t
        Upper bound.

    Returns
    -------
    None

    """

    try:
        elina_interval_set_mpq_c = elina_auxiliary_api.elina_interval_set_mpq
        elina_interval_set_mpq_c.restype = None
        elina_interval_set_mpq_c.argtypes = [ElinaIntervalPtr, Mpq_t, Mpq_t]
        elina_interval_set_mpq_c(interval, inf, sup)
    except:
        print('Problem with loading/calling "elina_interval_set_mpq" from "libelinaux.so"')
        print('Make sure you are passing ElinaIntervalPtr, Mpq_t and Mpq_t to the function')


def elina_interval_set_int(interval, inf, sup):
    """
    Set the value of an ElinaInterval by using two long integers.
    
    Parameters
    ----------
    interval : ElinaIntervalPtr
        Destination.
    inf : c_long
        Lower bound.
    sup : c_long
        Upper bound.

    Returns
    -------
    None

    """

    try:
        elina_interval_set_int_c = elina_auxiliary_api.elina_interval_set_int
        elina_interval_set_int_c.restype = None
        elina_interval_set_int_c.argtypes = [ElinaIntervalPtr, c_long, c_long]
        elina_interval_set_int_c(interval, inf, sup)
    except:
        print('Problem with loading/calling "elina_interval_set_int" from "libelinaux.so"')
        print('Make sure you are passing ElinaIntervalPtr, c_long and c_long to the function')


def elina_interval_set_frac(interval, num_inf, den_inf, num_sup, den_sup):
    """
    Set the value of an ElinaInterval by using four long integers.
    
    Parameters
    ----------
    interval : ElinaIntervalPtr
        Destination.
    num_inf : c_long
        Numerator of lower bound.
    den_inf : c_long 
        Denominator of lower bound.
    num_sup : c_long
        Numerator of upper bound.
    den_sup : c_long
        Denominator of upper bound.
        
    Returns
    -------
    None
    
    """

    try:
        elina_interval_set_frac_c = elina_auxiliary_api.elina_interval_set_frac
        elina_interval_set_frac_c.restype = None
        elina_interval_set_frac_c.argtypes = [ElinaIntervalPtr, c_long, c_ulong, c_long, c_ulong]
        elina_interval_set_frac_c(interval, num_inf, den_inf, num_sup, den_sup)
    except:
        print('Problem with loading/calling "elina_interval_set_frac" from "libelinaux.so"')
        print('Make sure you are passing ElinaIntervalPtr, c_long, c_ulong, c_long and c_ulong to the function')


def elina_interval_set_double(interval, inf, sup):
    """
    Set the value of an ElinaInterval by using two doubles.
    
    Parameters
    ----------
    interval : ElinaIntervalPtr
        Destination.
    inf : c_double
        Lower bound.
    sup : c_double
        Upper bound.

    Returns
    -------
    None
    
    """

    try:
        elina_interval_set_double_c = elina_auxiliary_api.elina_interval_set_double
        elina_interval_set_double_c.restype = None
        elina_interval_set_double_c.argtypes = [ElinaIntervalPtr, c_double, c_double]
        elina_interval_set_double_c(interval, inf, sup)
    except:
        print('Problem with loading/calling "elina_interval_set_double" from "libelinaux.so"')
        print('Make sure you are passing ElinaIntervalPtr, c_double and c_double to the function')


def elina_interval_set_mpfr(interval, inf, sup):
    """
    Set the value of an ElinaInterval by using two Mpfr-s.
    
    Parameters
    ----------
    interval : ElinaIntervalPtr
        Destination.
    inf : Mpfr_t
        Lower bound.
    sup : Mpfr_t
        Upper bound.

    Returns
    -------
    None

    """

    try:
        elina_interval_set_mpfr_c = elina_auxiliary_api.elina_interval_set_mpfr
        elina_interval_set_mpfr_c.restype = None
        elina_interval_set_mpfr_c.argtypes = [ElinaIntervalPtr, Mpfr_t, Mpfr_t]
        elina_interval_set_mpfr_c(interval, inf, sup)
    except:
        print('Problem with loading/calling "elina_interval_set_mpfr" from "libelinaux.so"')
        print('Make sure you are passing ElinaIntervalPtr, Mpfr_t and Mpfr_t to the function')


def elina_interval_set_top(interval):
    """
    Set an ElinaInterval to the universe interval [-oo, +oo].
    
    Parameters
    ----------
    interval : ElinaIntervalPtr
        Pointer to the ElinaInterval that needs to be set.

    Returns
    -------
    None
    
    """

    try:
        elina_interval_set_top_c = elina_auxiliary_api.elina_interval_set_top
        elina_interval_set_top_c.restype = None
        elina_interval_set_top_c.argtypes = [ElinaIntervalPtr]
        elina_interval_set_top_c(interval)
    except:
        print('Problem with loading/calling "elina_interval_set_top" from "libelinaux.so"')
        print('Make sure you are passing ElinaIntervalPtr to the function')


def elina_interval_set_bottom(interval):
    """
    Set an ElinaInterval to the empty interval [1, -1].
    
    Parameters
    ----------
    interval : ElinaIntervalPtr
        Pointer to the ElinaInterval that needs to be set.

    Returns
    -------
    None

    """

    try:
        elina_interval_set_bottom_c = elina_auxiliary_api.elina_interval_set_bottom
        elina_interval_set_bottom_c.restype = None
        elina_interval_set_bottom_c.argtypes = [ElinaIntervalPtr]
        elina_interval_set_bottom_c(interval)
    except:
        print('Problem with loading/calling "elina_interval_set_bottom" from "libelinaux.so"')
        print('Make sure you are passing ElinaIntervalPtr to the function')


# ====================================================================== #
# Combined allocation and assignments
# ====================================================================== #

def elina_interval_alloc_set(interval2):
    """
    Allocate a new ElinaInterval and initialise it with another ElinaInterval.
    
    Parameters
    ----------
    interval2 : ElinaIntervalPtr
        Pointer to the ElinaInterval used for initialisation.

    Returns
    -------
    interval1 : ElinaIntervalPtr
        Pointer to the newly allocated and initialised ElinaInterval

    """

    interval1 = None
    try:
        elina_interval_alloc_set_c = elina_auxiliary_api.elina_interval_alloc_set
        elina_interval_alloc_set_c.restype = ElinaIntervalPtr
        elina_interval_alloc_set_c.argtypes = [ElinaIntervalPtr]
        interval1 = elina_interval_alloc_set_c(interval2)
    except:
        print('Problem with loading/calling "elina_interval_alloc_set" from "libelinaux.so"')
        print('Make sure you are passing ElinaIntervalPtr to the function')

    return interval1


# ======================================================================
# Tests
# ======================================================================


def elina_interval_is_top(interval):
    """
    Test if an ElinaInterval belongs to the universe interval.
    
    Parameters
    ----------
    interval : ElinaIntervalPtr.
        Pointer to the ElinaInterval that needs to be tested.

    Returns
    -------
    result : c_bool
        Result of the test.

    """

    result = None
    try:
        elina_interval_is_top_c = elina_auxiliary_api.elina_interval_is_top
        elina_interval_is_top_c.restype = c_bool
        elina_interval_is_top_c.argtypes = [ElinaIntervalPtr]
        result = elina_interval_is_top_c(interval)
    except:
        print('Problem with loading/calling "elina_interval_is_top" from "libelinaux.so"')
        print('Make sure you are passing ElinaIntervalPtr to the function')

    return result


def elina_interval_is_bottom(interval):
    """
    Test if an ElinaInterval belongs to the empty interval.

    Parameters
    ----------
    interval : ElinaIntervalPtr.
        Pointer to the ElinaInterval that needs to be tested.

    Returns
    -------
    result : c_bool
        Result of the test.

    """

    result = None
    try:
        elina_interval_is_bottom_c = elina_auxiliary_api.elina_interval_is_bottom
        elina_interval_is_bottom_c.restype = c_bool
        elina_interval_is_bottom_c.argtypes = [ElinaIntervalPtr]
        result = elina_interval_is_bottom_c(interval)
    except:
        print('Problem with loading/calling "elina_interval_is_top" from "libelinaux.so"')
        print('Make sure you are passing ElinaIntervalPtr to the function')

    return result


def elina_interval_is_leq(interval1, interval2):
    """
    Test if one ElinaInterval is included into another ElinaInterval.
    
    Parameters
    ----------
    interval1 : ElinaIntervalPtr
        Pointer to the ElinaInterval that needs to be tested for inclusion in interval2.
    interval2 : ElinaIntervalPtr
        Pointer to the ElinaInterval that needs to be tested if it includes interval1.

    Returns
    -------
    result : c_bool
        Result of the test.
        
    """

    result = None
    try:
        elina_interval_is_leq_c = elina_auxiliary_api.elina_interval_is_leq
        elina_interval_is_leq_c.restype = c_bool
        elina_interval_is_leq_c.argtypes = [ElinaIntervalPtr, ElinaIntervalPtr]
        result = elina_interval_is_leq_c(interval1, interval2)
    except:
        print('Problem with loading/calling "elina_interval_is_leq" from "libelinaux.so"')
        print('Make sure you are passing ElinaIntervalPtr and ElinaIntervalPtr to the function')

    return result


def elina_interval_cmp(interval1, interval2):
    """
    Compare an ElinaInterval with another ElinaInterval.
    
    Parameters
    ----------
    interval1 : ElinaIntervalPtr
        Pointer to the ElinaInterval that needs to be compared.
    interval2 : ElinaIntervalPtr
        Pointer to the ElinaInterval that needs to be compared.

    Returns
    -------
    result : c_int
        The result of the comparison.
        Return 
        -1 if interval1 is included in interval2
        +1 if interval2 is included in interval1
        -2 if interval1.inf < interval2.inf
        +2 if interval1.inf > interval2.inf

    """

    result = None
    try:
        elina_interval_cmp_c = elina_auxiliary_api.elina_interval_cmp
        elina_interval_cmp_c.restype = c_int
        elina_interval_cmp_c.argtypes = [ElinaIntervalPtr, ElinaIntervalPtr]
        result = elina_interval_cmp_c(interval1, interval2)
    except:
        print('Problem with loading/calling "elina_interval_cmp" from "libelinaux.so"')
        print('Make sure you are passing ElinaIntervalPtr and ElinaIntervalPtr to the function')

    return result


def elina_interval_equal(interval1, interval2):
    """
    Test if an ElinaInterval is equal to another ElinaInterval.
    
    Parameters
    ----------
    interval1 : ElinaIntervalPtr
        Pointer to the ElinaInterval to be tested for equality.
    interval2 : ElinaIntervalPtr
        Pointer to the ElinaInterval to be tested for equality.

    Returns
    -------
    result : c_bool
        Result of the equality test.

    """

    result = None
    try:
        elina_interval_equal_c = elina_auxiliary_api.elina_interval_equal
        elina_interval_equal_c.restype = c_bool
        elina_interval_equal_c.argypes = [ElinaIntervalPtr, ElinaIntervalPtr]
        result = elina_interval_equal_c(interval1, interval2)
    except:
        print('Problem with loading/calling "elina_interval_equal" from "libelinaux.so"')
        print('Make sure you are passing ElinaIntervalPtr and ElinaIntervalPtr to the function')

    return result


def elina_interval_equal_int(interval, b):
    """
    Test if an ElinaInterval is equal to an integer.
    
    Parameters
    ----------
    interval : ElinaIntervalPtr
        Pointer to the ElinaInterval to be tested for equality.
    b : c_int
        Integer to be tested for equality.

    Returns
    -------
    result : c_bool
        Result of the equality test.

    """

    result = None
    try:
        elina_interval_equal_int_c = elina_auxiliary_api.elina_interval_equal_int
        elina_interval_equal_int_c.restype = c_bool
        elina_interval_equal_int_c.argtypes = [ElinaIntervalPtr, c_int]
        result = elina_interval_equal_int_c(interval, b)
    except:
        print('Problem with loading/calling "elina_interval_equal_int" from "libelinaux.so"')
        print('Make sure you are passing ElinaIntervalPtr and c_int to the function')

    return result


# ====================================================================== #
# Other operations
# ====================================================================== #


def elina_interval_neg(interval1, interval2):
    """
    Set an ElinaInterval to the negative of another ElinaInterval.
    
    Parameters
    ----------
    interval1 : ElinaIntervalPtr
        Destination.
    interval2 : ElinaIntervalPtr
        Source.

    Returns
    -------
    None

    """

    try:
        elina_interval_neg_c = elina_auxiliary_api.elina_interval_neg
        elina_interval_neg_c.restype = None
        elina_interval_neg_c.argtypes = [ElinaIntervalPtr, ElinaIntervalPtr]
        elina_interval_neg_c(interval1, interval2)
    except:
        print('Problem with loading/calling "elina_interval_neg" from "libelinaux.so"')
        print('Make sure you are passing ElinaIntervalPtr and ElinaIntervalPtr to the function')


def elina_interval_hash(interval):
    """
    Calculate the hash code of an ElinaInterval.
    
    Parameters
    ----------
    interval : ElinaIntervalPtr
        Pointer to the ElinaInterval that needs to be hashed.

    Returns
    -------
    result : c_long
        The resulting hash.

    """

    result = None
    try:
        elina_interval_hash_c = elina_auxiliary_api.elina_interval_hash
        elina_interval_hash_c.restype = c_long
        elina_interval_hash_c.argtypes = [ElinaIntervalPtr]
        result = elina_interval_hash_c(interval)
    except:
        print('Problem with loading/calling "elina_interval_hash" from "libelinaux.so"')
        print('Make sure you are passing ElinaIntervalPtr to the function')

    return result


# ======================================================================
# Array of intervals */
# ======================================================================


def elina_interval_array_alloc(size):
    """
    Allocate a new ElinaIntervalArray.
    
    Parameters
    ----------
    size : c_size_t
        Size of the ElinaIntervalArray allocated.

    Returns
    -------
    interval_array : ElinaIntervalArray
        The newly allocated ElinaIntervalArray.

    """

    interval_array = None
    try:
        elina_interval_array_alloc_c = elina_auxiliary_api.elina_interval_array_alloc
        elina_interval_array_alloc_c.restype = ElinaIntervalArray
        elina_interval_array_alloc_c.argtypes = [c_size_t]
        interval_array = elina_interval_array_alloc_c(size)
    except:
        print('Problem with loading/calling "elina_interval_array_alloc" from "libelinaux.so"')
        print('Make sure you are passing c_size_t to the function')

    return interval_array


def elina_interval_array_free(interval_array, size):
    """
    Free an ElinaIntervalArray.
    
    Parameters
    ----------
    interval_array : ElinaIntervalArray
        ElinaIntervalArray to be freed.
    size : c_size_t
        Size of the ElinaIntervalArray.

    Returns
    -------
    None

    """

    try:
        elina_interval_array_free_c = elina_auxiliary_api.elina_interval_array_free
        elina_interval_array_free_c.restype = None
        elina_interval_array_free_c.argtypes = [ElinaIntervalArray, c_size_t]
        elina_interval_array_free_c(interval_array, size)
    except:
        print('Problem with loading/calling "elina_interval_array_free" from "libelinaux.so"')
        print('Make sure you are passing ElinaIntervalPtr and c_size_t to the function')
