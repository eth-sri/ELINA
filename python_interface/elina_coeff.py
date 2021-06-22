#
#
#  This source file is part of ELINA (ETH LIbrary for Numerical Analysis).
#  ELINA is Copyright © 2021 Department of Computer Science, ETH Zurich
#  This software is distributed under GNU Lesser General Public License Version 3.0.
#  For more information, see the ELINA project website at:
#  http://elina.ethz.ch
#
#  THE SOFTWARE IS PROVIDED "AS-IS" WITHOUT ANY WARRANTY OF ANY KIND, EITHER
#  EXPRESS, IMPLIED OR STATUTORY, INCLUDING BUT NOT LIMITED TO ANY WARRANTY
#  THAT THE SOFTWARE WILL CONFORM TO SPECIFICATIONS OR BE ERROR-FREE AND ANY
#  IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE,
#  TITLE, OR NON-INFRINGEMENT.  IN NO EVENT SHALL ETH ZURICH BE LIABLE FOR ANY     
#  DAMAGES, INCLUDING BUT NOT LIMITED TO DIRECT, INDIRECT,
#  SPECIAL OR CONSEQUENTIAL DAMAGES, ARISING OUT OF, RESULTING FROM, OR IN
#  ANY WAY CONNECTED WITH THIS SOFTWARE (WHETHER OR NOT BASED UPON WARRANTY,
#  CONTRACT, TORT OR OTHERWISE).
#
#


from elina_coeff_h import *

# ====================================================================== #
# Basics
# ====================================================================== #


def elina_coeff_alloc(discr):
    """
    Allocate a new ElinaCoeff, using a specific type as the core.
    
    Parameters
    ----------
    discr : c_int
        Discriminant specifying the type of the core of the ElinaCoeff.

    Returns
    -------
    coeff : ElinaCoeffPtr
        Pointer to the newly allocated ElinaCoeff.

    """

    coeff = None
    try:
        elina_coeff_alloc_c = elina_auxiliary_api.elina_coeff_alloc
        elina_coeff_alloc_c.restype = ElinaCoeffPtr
        elina_coeff_alloc_c.argtypes = [c_uint]
        coeff = elina_coeff_alloc_c(discr)
    except:
        print('Problem with loading/calling "elina_coeff_alloc" from "libelinaux.so"')
        print('Make sure you are passing c_uint to the function')

    return coeff


def elina_coeff_reinit(coeff, coeff_discr, scalar_discr):
    """
    Reinitialise a given ElinaCoeff, according to the provided types.
    
    Parameters
    ----------
    coeff : ElinaCoeffPtr
        Pointer to the ElinaCoeff that needs to be reinitiliased.
    coeff_discr : c_uint
        Enum of type ElinaCoeffDiscr that defines the core of the ElinaCoeff.
    scalar_discr : c_uint
        Enum of type ElinaScalarDiscr that defines the core of the ElinaScalar (0 = double, 1 = mpq, 2 = mpfr).

    Returns
    -------
    None

    """

    try:
        elina_coeff_reinit_c = elina_auxiliary_api.elina_coeff_reinit
        elina_coeff_reinit_c.restype = None
        elina_coeff_reinit_c.argtypes = [ElinaCoeffPtr, c_uint, c_uint]
        elina_coeff_reinit_c(coeff, coeff_discr, scalar_discr)
    except:
        print('Problem with loading/calling "elina_coeff_reinit" from "libelinaux.so"')
        print('Make sure you are passing ElinaCoeffPtr, c_uint, c_uint to the function')


def elina_coeff_free(coeff):
    """
    Free an ElinaCoeff.
    
    Parameters
    ----------
    coeff : ElinaCoeffPtr
        Pointer to the ElinaCoeff that needs to be freed.

    Returns
    -------

    """

    try:
        elina_coeff_free_c = elina_auxiliary_api.elina_coeff_free
        elina_coeff_free_c.restype = None
        elina_coeff_free_c.argtypes = [ElinaCoeffPtr]
        elina_coeff_free_c(coeff)
    except:
        print('Problem with loading/calling "elina_coeff_free" from "libelinaux.so"')
        print('Make sure you are passing ElinaCoeffPtr to the function')


def elina_coeff_fprint(stream, coeff):
    """
    Print an ElinaCoeff onto a given stream.
    
    Parameters
    ----------
    stream : c_void_p
        Stream on which to print.
    coeff : ElinaCoeffPtr
        Pointer to the ElinaCoeff that needs to be printed.

    Returns
    -------
    None

    """

    try:
        elina_coeff_fprint_c = elina_auxiliary_api.elina_coeff_fprint
        elina_coeff_fprint_c.restype = None
        elina_coeff_fprint_c.argtypes = [c_void_p, ElinaCoeffPtr]
        elina_coeff_fprint_c(stream, coeff)
    except:
        print('Problem with loading/calling "elina_coeff_fprint" from "libelinaux.so"')
        print('Make sure you are passing c_void_p, ElinaCoeffPtr to the function')


def elina_coeff_reduce(coeff):
    """
    Reduce an ElinaCoeff of core type ElinaInterval [a, a], to an ElinaScalar.
    
    Parameters
    ----------
    coeff : ElinaCoeffPtr
        Pointer to the ElinaCoeff that needs to be reduced.

    Returns
    -------
    None

    """

    try:
        elina_coeff_reduce_c = elina_auxiliary_api.elina_coeff_reduce
        elina_coeff_reduce_c.restype = None
        elina_coeff_reduce_c.argtypes = [ElinaCoeffPtr]
        elina_coeff_reduce_c(coeff)
    except:
        print('Problem with loading/calling "elina_coeff_reduce" from "libelinaux.so"')
        print('Make sure you are passing ElinaCoeffPtr to the function')


# ====================================================================== #
# Assignments
# ====================================================================== #


def elina_coeff_set(coeff1, coeff2):
    """
    Set the value of one ElinaCoeff to the value of another ElinaCoeff.
    
    Parameters
    ----------
    coeff1 : ElinaCoeffPtr
        Destination.
    coeff2 : ElinaCoeffPtr
        Source

    Returns
    -------
    None

    """

    try:
        elina_coeff_set_c = elina_auxiliary_api.elina_coeff_set
        elina_coeff_set_c.restype = None
        elina_coeff_set_c.argtypes = [ElinaCoeffPtr, ElinaCoeffPtr]
        elina_coeff_set_c(coeff1, coeff2)
    except:
        print('Problem with loading/calling "elina_coeff_set" from "libelinaux.so"')
        print('Make sure you are passing ElinaCoeffPtr, ElinaCoeffPtr to the function')


def elina_coeff_set_scalar(coeff, scalar):
    """
    Set the value of an ElinaCoeff with core ElinaScalar by using an ElinaScalar.
    
    Parameters
    ----------
    coeff : ElinaCoeffPtr
        Destination.
    scalar : ElinaScalarPtr
        Source.

    Returns
    -------

    """

    try:
        elina_coeff_set_scalar_c = elina_auxiliary_api.elina_coeff_set_scalar
        elina_coeff_set_scalar_c.restype = None
        elina_coeff_set_scalar_c.argtypes = [ElinaCoeffPtr, ElinaScalarPtr]
        elina_coeff_set_scalar_c(coeff, scalar)
    except:
        print('Problem with loading/calling "elina_coeff_set_scalar" from "libelinaux.so"')
        print('Make sure you are passing ElinaCoeffPtr, ElinaScalarPtr to the function')


def elina_coeff_set_scalar_mpq(coeff, mpq_t):
    """
    Set the value of an ElinaCoeff with core ElinaScalar by using a Mpq_t.
    
    Parameters
    ----------
    coeff : ElinaCoeffPtr
        Destination.
    mpq_t : Mpq_t
        Source.

    Returns
    -------
    None
    
    """

    try:
        elina_coeff_set_scalar_mpq_c = elina_auxiliary_api.elina_coeff_set_scalar_mpq
        elina_coeff_set_scalar_mpq_c.restype = None
        elina_coeff_set_scalar_mpq_c.argypes = [ElinaCoeffPtr, Mpq_t]
        elina_coeff_set_scalar_mpq_c(coeff, mpq_t)
    except:
        print('Problem with loading/calling "elina_coeff_set_scalar_mpq_c" from "libelinaux.so"')
        print('Make sure you are passing ElinaCoeffPtr, Mpq_t to the function')


def elina_coeff_set_scalar_int(coeff, num):
    """
    Set the value of an ElinaCoeff with core ElinaScalar by using a long integer.
    
    Parameters
    ----------
    coeff : ElinaCoefPtr
        Destination.
    num : c_long
        Source.

    Returns
    -------
    None

    """

    try:
        elina_coeff_set_scalar_int_c = elina_auxiliary_api.elina_coeff_set_scalar_int
        elina_coeff_set_scalar_int_c.restype = None
        elina_coeff_set_scalar_int_c.argypes = [ElinaCoeffPtr, c_long]
        elina_coeff_set_scalar_int_c(coeff, num)
    except:
        print('Problem with loading/calling "elina_coeff_set_scalar_int" from "libelinaux.so"')
        print('Make sure you are passing ElinaCoeffPtr, c_long to the function')


def elina_coeff_set_scalar_frac(coeff, num, den):
    """
    Set the value of an ElinaCoeff with core ElinaScalar by using fraction of two long integers.
    
    Parameters
    ----------
    coeff : ElinaCoeffPtr
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
        elina_coeff_set_scalar_frac_c = elina_auxiliary_api.elina_coeff_set_scalar_frac
        elina_coeff_set_scalar_frac_c.restype = None
        elina_coeff_set_scalar_frac_c.argypes = [ElinaCoeffPtr, c_long, c_ulong]
        elina_coeff_set_scalar_frac_c(coeff, num, den)
    except:
        print('Problem with loading/calling "elina_coeff_set_scalar_frac" from "libelinaux.so"')
        print('Make sure you are passing ElinaCoeffPtr, c_long, c_ulong to the function')


def elina_coeff_set_scalar_double(coeff, num):
    """
    Set the value of an ElinaCoeff with core ElinaScalar by using a double.
    
    Parameters
    ----------
    coeff : ElinaCoeffPtr
        Destination.
    num : c_double
        Source.

    Returns
    -------
    None
    
    """

    try:
        elina_coeff_set_scalar_double_c = elina_auxiliary_api.elina_coeff_set_scalar_double
        elina_coeff_set_scalar_double_c.restype = None
        elina_coeff_set_scalar_double_c.argypes = [ElinaCoeffPtr, c_double]
        elina_coeff_set_scalar_double_c(coeff, num)
    except:
        print('Problem with loading/calling "elina_coeff_set_scalar_double" from "libelinaux.so"')
        print('Make sure you are passing ElinaCoeffPtr, c_double to the function')


def elina_coeff_set_scalar_mpfr(coeff, mpfr_t):
    """
    Set the value of an ElinaCoeff with core ElinaScalar by using a Mpfr_t.
    
    Parameters
    ----------
    coeff : ElinaCoeffPtr
        Destination.
    mpfr_t : Mpfr_t
        Source.

    Returns
    -------
    None

    """

    try:
        elina_coeff_set_scalar_mpfr_c = elina_auxiliary_api.elina_coeff_set_scalar_mpfr
        elina_coeff_set_scalar_mpfr_c.restype = None
        elina_coeff_set_scalar_mpfr_c.argtypes = [ElinaCoeffPtr, Mpfr_t]
        elina_coeff_set_scalar_mpfr_c(coeff, mpfr_t)
    except:
        print('Problem with loading/calling "elina_coeff_set_scalar_mpfr" from "libelinaux.so"')
        print('Make sure you are passing ElinaCoeffPtr, Mpfr_t to the function')


def elina_coeff_set_interval(coeff, interval):
    """
    Set the value of an ElinaCoeff with core ElinaInterval by using an ElinaInterval
    
    Parameters
    ----------
    coeff : ElinaCoeffPtr
        Destination.
    interval : ElinaIntervalPtr
        Source.

    Returns
    -------
    None

    """

    try:
        elina_coeff_set_interval_c = elina_auxiliary_api.elina_coeff_set_interval
        elina_coeff_set_interval_c.restype = None
        elina_coeff_set_interval_c.argtypes = [ElinaCoeffPtr, ElinaIntervalPtr]
        elina_coeff_set_interval_c(coeff, interval)
    except:
        print('Problem with loading/calling "elina_coeff_set_interval" from "libelinaux.so"')
        print('Make sure you are passing ElinaCoeffPtr, ElinaIntervalPtr to the function')


def elina_coeff_set_interval_scalar(coeff, inf, sup):
    """
    Set the value of an ElinaCoeff with core ElinaInterval by using two ElinaScalar-s.
    
    Parameters
    ----------
    coeff : ElinaCoeffPtr
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
        elina_coeff_set_interval_scalar_c = elina_auxiliary_api.elina_coeff_set_interval_scalar
        elina_coeff_set_interval_scalar_c.restype = None
        elina_coeff_set_interval_scalar_c.argtypes = [ElinaCoeffPtr, ElinaScalarPtr, ElinaScalarPtr]
        elina_coeff_set_interval_scalar_c(coeff, inf, sup)
    except:
        print('Problem with loading/calling "elina_coeff_set_interval_scalar" from "libelinaux.so"')
        print('Make sure you are passing ElinaCoeffPtr, ElinaScalarPtr, ElinaScalarPtr to the function')


def elina_coeff_set_interval_mpq(coeff, inf, sup):
    """
    Set the value of an ElinaCoeff with core ElinaInterval by using two Mpq_t-s.
    
    Parameters
    ----------
    coeff : ElinaCoeffPtr
        Destination.
    inf : Mpq_t
        Source.
    sup : Mpq_t
        Source.

    Returns
    -------
    None

    """

    try:
        elina_coeff_set_interval_mpq_c = elina_auxiliary_api.elina_coeff_set_interval_mpq
        elina_coeff_set_interval_mpq_c.restype = None
        elina_coeff_set_interval_mpq_c.argypes = [ElinaCoeffPtr, Mpq_t, Mpq_t]
        elina_coeff_set_interval_mpq_c(coeff, inf, sup)
    except:
        print('Problem with loading/calling "elina_coeff_set_interval_mpq" from "libelinaux.so"')
        print('Make sure you are passing ElinaCoeffPtr, Mpq_t, Mpq_t to the function')


def elina_coeff_set_interval_int(coeff, inf, sup):
    """
    Set the value of an ElinaCoeff with core ElinaInterval by using two long integers.
    
    Parameters
    ----------
    coeff : ElinaCoeffPtr
        Destination.
    inf : c_long
        Source.
    sup : c_long
        Source.

    Returns
    -------
    None

    """

    try:
        elina_coeff_set_interval_int_c = elina_auxiliary_api.elina_coeff_set_interval_int
        elina_coeff_set_interval_int_c.restype = None
        elina_coeff_set_interval_int_c.argtypes = [ElinaCoeffPtr, c_long, c_long]
        elina_coeff_set_interval_int_c(coeff, inf, sup)
    except:
        print('Problem with loading/calling "elina_coeff_set_interval_int" from "libelinaux.so"')
        print('Make sure you are passing ElinaCoeffPtr, c_long, c_long to the function')


def elina_coeff_set_interval_frac(coeff, numinf, deninf, numsup, densup):
    """
    Set the value of an ElinaCoeff with core ElinaInterval by using two pairs of long integers as fractions.
    
    Parameters
    ----------
    coeff : ElinaCoeffPtr
        Destination.
    numinf : c_long
        Source.
    deninf : c_ulong
        Source.
    numsup : c_long
        Source.
    densup : c_ulong
        Source.

    Returns
    -------
    None

    """

    try:
        elina_coeff_set_interval_frac_c = elina_auxiliary_api.elina_coeff_set_interval_frac
        elina_coeff_set_interval_frac_c.restype = None
        elina_coeff_set_interval_frac_c.argtypes = [ElinaCoeffPtr, c_long, c_ulong, c_long, c_ulong]
        elina_coeff_set_interval_frac_c(coeff, numinf, deninf, numsup, densup)
    except:
        print('Problem with loading/calling "elina_coeff_set_interval_frac" from "libelinaux.so"')
        print('Make sure you are passing ElinaCoeffPtr, c_long, c_ulong, c_long, c_ulong to the function')


def elina_coeff_set_interval_double(coeff, inf, sup):
    """
    Set the value of an ElinaCoeff with core ElinaInterval by using two double-s.
    
    Parameters
    ----------
    coeff : ElinaCoeffPtr
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
        elina_coeff_set_interval_double_c = elina_auxiliary_api.elina_coeff_set_interval_double
        elina_coeff_set_interval_double_c.restype = None
        elina_coeff_set_interval_double_c.argtypes = [ElinaCoeffPtr, c_double, c_double]
        elina_coeff_set_interval_double_c(coeff, inf, sup)
    except:
        print('Problem with loading/calling "elina_coeff_set_interval_double" from "libelinaux.so"')
        print('Make sure you are passing ElinaCoeffPtr, c_double, c_double to the function')


def elina_coeff_set_interval_top(coeff):
    """
    Set the value of an ElinaCoeff with core ElinaInterval by using the universe interval [-oo, +oo].
    
    Parameters
    ----------
    coeff : ElinaCoeffPtr
        Destination.
        
    Returns
    -------
    None
    
    """

    try:
        elina_coeff_set_interval_top_c = elina_auxiliary_api.elina_coeff_set_interval_top
        elina_coeff_set_interval_top_c.restype = None
        elina_coeff_set_interval_top_c.argtypes = [ElinaCoeffPtr]
        elina_coeff_set_interval_top_c(coeff)
    except:
        print('Problem with loading/calling "elina_coeff_set_interval_top" from "libelinaux.so"')
        print('Make sure you are passing ElinaCoeffPtr to the function')


def elina_coeff_set_interval_mpfr(coeff, inf, sup):
    """
    Set the value of an ElinaCoeff with core ElinaInterval by using two Mpfr_t-s.
    
    Parameters
    ----------
    coeff : ElinaCoeffPtr
        Destination.
    inf : Mpfr_t
        Source.
    sup : Mpfr_t
        Source.

    Returns
    -------
    None

    """

    try:
        elina_coeff_set_interval_mpfr_c = elina_auxiliary_api.elina_coeff_set_interval_mpfr
        elina_coeff_set_interval_mpfr_c.restype = None
        elina_coeff_set_interval_mpfr_c.argtypes = [ElinaCoeffPtr, Mpfr_t, Mpfr_t]
        elina_coeff_set_interval_mpfr_c(coeff, inf, sup)
    except:
        print('Problem with loading/calling "elina_coeff_set_interval_mpfr" from "libelinaux.so"')
        print('Make sure you are passing ElinaCoeffPtr, Mpfr_t, Mpfr_t to the function')


# ====================================================================== #
# Combined allocation and assignment
# ====================================================================== #


def elina_coeff_alloc_set(coeff2):
    """
    Allocate a new ElinaCoeff and initialise it with another ElinaCoeff.
    
    Parameters
    ----------
    coeff2 : ElinaCoeffPtr
        Pointer to the ElinaCoeff used for initialisation.

    Returns
    -------
    coeff1: ElinaCoeffPtr
        Pointer to the newly allocated and initialised ElinaCoeff.

    """

    coeff1 = None
    try:
        elina_coeff_alloc_set_c = elina_auxiliary_api.elina_coeff_alloc_set
        elina_coeff_alloc_set_c.restype = ElinaCoeffPtr
        elina_coeff_alloc_set_c.argtypes = [ElinaCoeffPtr]
        coeff1 = elina_coeff_alloc_set_c(coeff2)
    except:
        print('Problem with loading/calling "elina_coeff_alloc_set" from "libelinaux.so"')
        print('Make sure you are passing ElinaCoeffPtr to the function')

    return coeff1


def elina_coeff_alloc_set_scalar(scalar):
    """
    Allocate a new ElinaCoeff and initialise it with an ElinaScalar.
    
    Parameters
    ----------
    scalar : ElinaScalarPtr
        Pointer to the ElinaScalar used for initialisation.

    Returns
    -------
    coeff : ElinaCoeffPtr
        Pointer to the newly allocated and initialised ElinaCoeff.

    """

    coeff = None
    try:
        elina_coeff_alloc_set_scalar_c = elina_auxiliary_api.elina_coeff_alloc_set_scalar
        elina_coeff_alloc_set_scalar_c.restype = None
        elina_coeff_alloc_set_scalar_c.argtypes = [ElinaScalarPtr]
        coeff = elina_coeff_alloc_set_scalar_c(scalar)
    except:
        print('Problem with loading/calling "elina_coeff_alloc_set_scalar" from "libelinaux.so"')
        print('Make sure you are passing ElinaScalarPtr to the function')

    return coeff


def elina_coeff_alloc_set_interval(interval):
    """
    Allocate a new ElinaCoeff and initialise it with an ElinaInterval.
    
    Parameters
    ----------
    interval : ElinaIntervalPtr
        Pointer to the ElinaInterval used for initialisation.

    Returns
    -------
    coeff : ElinaCoeffPtr
        Pointer to the newly allocated and initialised ElinaCoeff.

    """

    coeff = None
    try:
        elina_coeff_alloc_set_interval_c = elina_auxiliary_api.elina_coeff_alloc_set_interval
        elina_coeff_alloc_set_interval_c.restype = ElinaCoeffPtr
        elina_coeff_alloc_set_interval_c.argtypes = [ElinaIntervalPtr]
        elina_coeff_alloc_set_interval_c(interval)
    except:
        print('Problem with loading/calling "elina_coeff_alloc_set_interval" from "libelinaux.so"')
        print('Make sure you are passing ElinaIntervalPtr to the function')

    return coeff


# ====================================================================== #
# Tests
# ====================================================================== #


def elina_coeff_cmp(coeff1, coeff2):
    """
    Compare an ElinaCoeff with another ElinaCoeff.
    
    Parameters
    ----------
    coeff1 : ElinaCoeffPtr
        Pointer to the ElinaCoeff that needs to be compared.
    coeff2 : ElinaCoeffPtr
        Pointer to the ElinaCoeff that needs to be compared.

    Returns
    -------
    result : c_int
        The result of the comparison.
        Return 
        corresponding to elina_scalar_cmp if the two ElinaCoeff-s have ElinaScalar core
        corresponding to elina_interval_cmp if the two ElinaCoeff-s have ElinaInterval core
        -3 if the first ElinaCoeff has an ElinaScalar core
        +3 if the second ElinaCoeff has an ElinaScalar core

    """

    result = None
    try:
        elina_coeff_cmp_c = elina_auxiliary_api.elina_coeff_cmp
        elina_coeff_cmp_c.restype = c_int
        elina_coeff_cmp_c.argtypes = [ElinaCoeffPtr, ElinaCoeffPtr]
        result = elina_coeff_cmp_c(coeff1, coeff2)
    except:
        print('Problem with loading/calling "elina_coeff_cmp" from "libelinaux.so"')
        print('Make sure you are passing ElinaCoeffPtr, ElinaCoeffPtr to the function')

    return result


def elina_coeff_equal(coeff1, coeff2):
    """
    Test if an ElinaCoeff is equal to another ElinaCoeff.
    
    Parameters
    ----------
    coeff1 : ElinaCoeffPtr
        Pointer to the ElinaCoeff that needs to be tested for equality.
    coeff2 : ElinaCoefPtr
        Pointer to the ElinaCoeff that needs to be tested for equality.

    Returns
    -------
    result : c_bool
        Result of the equality test.

    """

    result = None
    try:
        elina_coeff_equal_c = elina_auxiliary_api.elina_coeff_equal
        elina_coeff_equal_c.restype = c_bool
        elina_coeff_equal_c.argtypes = [ElinaCoeffPtr, ElinaCoeffPtr]
        result = elina_coeff_equal_c(coeff1, coeff2)
    except:
        print('Problem with loading/calling "elina_coeff_equal" from "libelinaux.so"')
        print('Make sure you are passing ElinaCoeffPtr, ElinaCoeffPtr to the function')

    return result


def elina_coeff_zero(coeff):
    """
    Test if an ElinaCoeff is a zero ElinaScalar or an ElinaInterval with zero bounds.

    Parameters
    ----------
    coeff : ElinaCoefPtr
        Pointer to the ElinaCoeff that needs to be tested.

    Returns
    -------
    result : c_bool
        Result of the zero test.

    """

    result = None
    try:
        elina_coeff_zero_c = elina_auxiliary_api.elina_coeff_zero
        elina_coeff_zero_c.restype = c_int
        elina_coeff_zero_c.argtypes = [ElinaCoeffPtr]
        result = elina_coeff_zero_c(coeff)
    except:
        print('Problem with loading/calling "elina_coeff_zero" from "libelinaux.so"')
        print('Make sure you are passing ElinaCoeffPtr to the function')

    return result


def elina_coeff_equal_int(coeff, i):
    """
    Test if an ElinaCoeff is equal to an integer.

    Parameters
    ----------
    coeff : ElinaCoeffPtr
        Pointer to the ElinaCoeff that needs to be tested for equality.
    i : c_int
       Integer that needs to be tested for equality.

    Returns
    -------
    result : c_bool
        Result of the equality test.

    """

    result = None
    try:
        elina_coeff_equal_int_c = elina_auxiliary_api.elina_coeff_equal_int
        elina_coeff_equal_int_c.restype = c_bool
        elina_coeff_equal_int_c.argtypes = [ElinaCoeffPtr, c_int]
        result = elina_coeff_equal_int_c(coeff, i)
    except:
        print('Problem with loading/calling "elina_coeff_equal_int" from "libelinaux.so"')
        print('Make sure you are passing ElinaCoeffPtr, c_int to the function')

    return result


# ====================================================================== #
# Other operations
# ====================================================================== #


def elina_coeff_neg(coeff1, coeff2):
    """
    Set the value of an ElinaCoeff to the negative of another ElinaCoeff.
    
    Parameters
    ----------
    coeff1 : ElinaCoeffPtr
        Destination.
    coeff2 : ElinaCoeffPtr
        Source.
        
    Returns
    -------
    None
    
    """

    try:
        elina_coeff_neg_c = elina_auxiliary_api.elina_coeff_neg
        elina_coeff_neg_c.restype = None
        elina_coeff_neg_c.argtypes = [ElinaCoeffPtr, ElinaCoeffPtr]
        elina_coeff_neg_c(coeff1, coeff2)
    except:
        print('Problem with loading/calling "elina_coeff_neg" from "libelinaux.so"')
        print('Make sure you are passing ElinaCoeffPtr, ElinaCoeffPtr to the function')


def elina_coeff_hash(coeff):
    """
    Calculate the hash code of an ElinaCoeff.
    
    Parameters
    ----------
    coeff : ElinaCoeffPtr
        Pointer to the ElinaCoeff that needs to be hashed.

    Returns
    -------
    result : c_long
        The resulting hash.

    """

    result = None
    try:
        elina_coeff_hash_c = elina_auxiliary_api.elina_coeff_hash
        elina_coeff_hash_c.restype = c_long
        elina_coeff_hash_c.argtypes = [ElinaCoeffPtr]
        result = elina_coeff_hash_c(coeff)
    except:
        print('Problem with loading/calling "elina_coeff_hash" from "libelinaux.so"')
        print('Make sure you are passing c_long to the function')

    return result
