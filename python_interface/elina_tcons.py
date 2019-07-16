#
#
#  This source file is part of ELINA (ETH LIbrary for Numerical Analysis).
#  ELINA is Copyright © 2019 Department of Computer Science, ETH Zurich
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


from elina_tcons0_h import *
from elina_auxiliary_imports import *

# ********************************************************************** #
# I. elina_tcons0_t
# ********************************************************************** #


# ====================================================================== #
# I.1 Memory management and printing
# ====================================================================== #

def elina_tcons0_make(constyp, texpr, scalar):
    """
    Create an ElinaTcons0 of a type defined as ElinaConstyp given an ElinaTexpr0.
    The expression and the coefficient are not duplicated, just pointed to.
    
    Parameters
    ----------
    constyp : c_uint
       Enum compatible with ElinaConstyp defining the type of the constraint.
    texpr : ElinaTexpr0Ptr
       Pointer to the ElinaTexpr0 that is used for the constraint.
    scalar : ElinaScalarPtr
       Pointer to the ElinaScalar used as a modulo in case of constyp == EQMOD. Null otherwise.
    
    Returns
    -------
    tcons = ElinaTcons0
       Newly created ElinaTcons0.
    
    """

    tcons = None
    try:
        elina_tcons0_make_c = elina_auxiliary_api.elina_tcons0_make
        elina_tcons0_make_c.restype = ElinaTcons0
        elina_tcons0_make_c.argtypes = [c_uint, ElinaTexpr0Ptr, ElinaScalarPtr]
        tcons = elina_tcons0_make_c(constyp, texpr, scalar)
    except:
        print('Problem with loading/calling "elina_tcons0_make" from "libelinaux.so"')
        print('Make sure you are passing c_uint, ElinaTexpr0Ptr, ElinaScalarPtr to the function')

    return tcons


def elina_tcons0_make_unsat():
    """
    Create the constraint -1 >= 0.

    Returns
    -------
    tcons : ElinaTcons0
        Newly created ElinaTcons0.

    """

    tcons = None
    try:
        elina_tcons0_make_unsat_c = elina_auxiliary_api.elina_tcons0_make_unsat
        elina_tcons0_make_unsat_c.restype = ElinaTcons0
        tcons = elina_tcons0_make_unsat_c()
    except:
        print('Problem with loading/calling "elina_tcons0_make_unsat" from "libelinaux.so"')

    return tcons


def elina_tcons0_from_lincons0(lincons):
    """
    Create an ElinaTcons0 from an ElinaLincons0.
    
    Parameters
    ----------
    lincons : ElinaLincons0Ptr
        Pointer to the ElinaLincons0 that needs to be transformed to ElinaTcons0.

    Returns
    -------
    tcons : ElinaTcons0
        The newly created ElinaTcons0.

    """

    tcons = None
    try:
        elina_tcons0_from_lincons0_c = elina_auxiliary_api.elina_tcons0_from_lincons0
        elina_tcons0_from_lincons0_c.restype = ElinaTcons0
        elina_tcons0_from_lincons0_c.argtypes = [ElinaLincons0Ptr]
        tcons = elina_tcons0_from_lincons0_c(lincons)
    except:
        print('Problem with loading/calling "elina_tcons0_from_lincons0" from "libelinaux.so"')
        print('Make sure you are passing ElinaLincons0Ptr to the function')

    return tcons


def elina_tcons0_copy(tcons2):
    """
    Duplicate an ElinaTcons0.

    Parameters
    ----------
    tcons2 : ElinaTcons0Ptr
        Pointer to the ElinaTcons0 that needs to be duplicated.

    Returns
    -------
    tcons1 : ElinaTcons0
        The duplicate ElinaTcons0.

    """

    tcons1 = None
    try:
        elina_tcons0_copy_c = elina_auxiliary_api.elina_tcons0_copy
        elina_tcons0_copy_c.restype = ElinaTcons0
        elina_tcons0_copy_c.argtypes = [ElinaTcons0Ptr]
        tcons1 = elina_tcons0_copy_c(tcons2)
    except:
        print('Problem with loading/calling "elina_tcons0_copy" from "libelinaux.so"')
        print('Make sure you are passing ElinaTcons0Ptr to the function')

    return tcons1


def elina_tcons0_clear(tcons):
    """
    Free an ElinaTcons0 and set the pointer to NULL.

    Parameters
    ----------
    tcons : ElinaTcons0Ptr
        Pointer to the ElinaTcons0 that needs to be freed.

    Returns
    -------
    None

    """

    try:
        elina_tcons0_clear_c = elina_auxiliary_api.elina_tcons0_clear
        elina_tcons0_clear_c.restype = None
        elina_tcons0_clear_c.argtypes = [ElinaTcons0Ptr]
        elina_tcons0_clear_c(tcons)
    except:
        print('Problem with loading/calling "elina_tcons0_clear" from "libelinaux.so"')
        print('Make sure you are passing ElinaTcons0Ptr to the function')


def elina_tcons0_fprint(stream, tcons, name_of_dim):
    """
    Print an ElinaTcons0, having dimension names to stream.

    Parameters
    ----------
    stream : c_void_p
        Stream onto which to print.
    tcons : ElinaTcons0Ptr
        Pointer to the ElinaTcons0 that needs to be printed.
    name_of_dim : POINTER(c_char_p)
        Dimension names.

    Returns
    -------
    None

    """

    try:
        elina_tcons0_fprint_c = elina_auxiliary_api.elina_tcons0_fprint
        elina_tcons0_fprint_c.restype = None
        elina_tcons0_fprint_c.argtypes = [c_void_p, ElinaTcons0Ptr, POINTER(c_char_p)]
        elina_tcons0_fprint_c(stream, tcons, name_of_dim)
    except:
        print('Problem with loading/calling "elina_tcons0_fprint" from "libelinaux.so"')
        print('Make sure you are passing c_void_p, ElinaTcons0Ptr, POINTER(c_char_p) to the function')


# ====================================================================== #
# I.2 Tests
# ====================================================================== #

def elina_tcons0_is_interval_cst(tcons):
    """
    Test if an ElinaTcons0 contains only constant leaves, i.e. no variables.
    
    Parameters
    ----------
    tcons : ElinaTcons0Ptr
        Pointer to the ElinaTcons0 that needs to be tested.

    Returns
    -------
    result : c_bool
        Result of the test.

    """

    result = None
    try:
        elina_tcons0_is_interval_cst_c = elina_auxiliary_api.elina_tcons0_is_interval_cst
        elina_tcons0_is_interval_cst_c.restype = c_bool
        elina_tcons0_is_interval_cst_c.argtypes = [ElinaTcons0Ptr]
        result = elina_tcons0_is_interval_cst_c(tcons)
    except:
        print('Problem with loading/calling "elina_tcons0_is_interval_cst" from "libelinaux.so"')
        print('Make sure you are passing ElinaTcons0Ptr to the function')

    return result


def elina_tcons0_is_interval_linear(tcons):
    """
    Test if an ElinaTcons0 is linear with possibly interval coefficients.
    No rounding performed.

    Parameters
    ----------
    tcons : ElinaTcons0Ptr
        Pointer to the ElinaTcons0 that needs to be tested.

    Returns
    -------
    result : c_bool
        Result of the test.

    """

    result = None
    try:
        elina_tcons0_is_interval_linear_c = elina_auxiliary_api.elina_tcons0_is_interval_linear
        elina_tcons0_is_interval_linear_c.restype = c_bool
        elina_tcons0_is_interval_linear_c.argtypes = [ElinaTcons0Ptr]
        result = elina_tcons0_is_interval_linear_c(tcons)
    except:
        print('Problem with loading/calling "elina_tcons0_is_interval_linear" from "libelinaux.so"')
        print('Make sure you are passing ElinaTcons0Ptr to the function')

    return result


def elina_tcons0_is_interval_polynomial(tcons):
    """
    Test if an ElinaTcons0 is polynomial with possibly interval coefficients.
    No rounding performed.

    Parameters
    ----------
    tcons : ElinaTcons0Ptr
        Pointer to the ElinaTcons0 that needs to be tested.

    Returns
    -------
    result : c_bool
        Result of the test.

    """

    result = None
    try:
        elina_tcons0_is_interval_polynomial_c = elina_auxiliary_api.elina_tcons0_is_interval_polynomial
        elina_tcons0_is_interval_polynomial_c.restype = c_bool
        elina_tcons0_is_interval_polynomial_c.argtypes = [ElinaTcons0Ptr]
        result = elina_tcons0_is_interval_polynomial_c(tcons)
    except:
        print('Problem with loading/calling "elina_tcons0_is_interval_polynomial" from "libelinaux.so"')
        print('Make sure you are passing ElinaTcons0Ptr to the function')

    return result


def elina_tcons0_is_interval_polyfrac(tcons):
    """
    Test if an ElinaTcons0 is polynomial fraction with possibly interval coefficients.
    No rounding performed.

    Parameters
    ----------
    tcons : ElinaTcons0Ptr
        Pointer to the ElinaTcons0 that needs to be tested.

    Returns
    -------
    result : c_bool
        Result of the test.

    """

    result = None
    try:
        elina_tcons0_is_interval_polyfrac_c = elina_auxiliary_api.elina_tcons0_is_interval_polyfrac
        elina_tcons0_is_interval_polyfrac_c.restype = c_bool
        elina_tcons0_is_interval_polyfrac_c.argtypes = [ElinaTcons0Ptr]
        result = elina_tcons0_is_interval_polyfrac_c(tcons)
    except:
        print('Problem with loading/calling "elina_tcons0_is_interval_polyfrac" from "libelinaux.so"')
        print('Make sure you are passing ElinaTcons0Ptr to the function')

    return result


def elina_tcons0_is_scalar(tcons):
    """
    Test if an ElinaTcons0 contains only scalar coefficients.

    Parameters
    ----------
    tcons : ElinaTcons0Ptr
        Pointer to the ElinaTcons0 that needs to be tested.

    Returns
    -------
    result : c_bool
        Result of the test.

    """

    result = None
    try:
        elina_tcons0_is_scalar_c = elina_auxiliary_api.elina_tcons0_is_scalar
        elina_tcons0_is_scalar_c.restype = c_bool
        elina_tcons0_is_scalar_c.argtypes = [ElinaTcons0Ptr]
        result = elina_tcons0_is_scalar_c(tcons)
    except:
        print('Problem with loading/calling "elina_tcons0_is_scalar" from "libelinaux.so"')
        print('Make sure you are passing ElinaTcons0Ptr to the function')

    return result


# ====================================================================== #
# I.3 Change of dimensions and permutations
# ====================================================================== #

def elina_tcons0_add_dimensions_with(tcons, dimchange):
    """
    Add dimensions to an ElinaTcons0 by following the semantics of an ElinaDimchange.

    Parameters
    ----------
    tcons : ElinaTcons0Ptr
        Pointer to the ElinaTcons0 to which we want to add dimensions.
    dimchange : ElinaDimchangePtr
        Pointer to the ElinaDimchange which semantics we want to follow.

    Returns
    -------
    None

    """

    try:
        elina_tcons0_add_dimensions_with_c = elina_auxiliary_api.elina_tcons0_add_dimensions_with
        elina_tcons0_add_dimensions_with_c.restype = None
        elina_tcons0_add_dimensions_with_c.argtypes = [ElinaTcons0Ptr, ElinaDimchangePtr]
        elina_tcons0_add_dimensions_with_c(tcons, dimchange)
    except:
        print('Problem with loading/calling "elina_tcons0_add_dimensions_with" from "libelinaux.so"')
        print('Make sure you are passing ElinaTcons0Ptr, ElinaDimchangePtr to the function')


def elina_tcons0_add_dimensions(tcons2, dimchange):
    """
    Add dimensions to an ElinaTcons0 by following the semantics of an ElinaDimchange.

    Parameters
    ----------
    tcons2 : ElinaTconsPtr
        Pointer to the ElinaTcons0 to which we want to add dimensions.
    dimchange : ElinaDimchangePtr
        Pointer to the ElinaDimchange which semantics we want to follow.

    Returns
    -------
    tcons1 : ElinaTcons0Ptr
        Pointer to the newly created ElinaTcons0 with added dimensions.

    """

    tcons1 = None
    try:
        elina_tcons0_add_dimensions_c = elina_auxiliary_api.elina_tcons0_add_dimensions
        elina_tcons0_add_dimensions_c.restype = ElinaTcons0
        elina_tcons0_add_dimensions_c.argtypes = [ElinaTcons0Ptr, ElinaDimchangePtr]
        tcons1 = elina_tcons0_add_dimensions_c(tcons2, dimchange)
    except:
        print('Problem with loading/calling "elina_tcons0_add_dimensions" from "libelinaux.so"')
        print('Make sure you are passing ElinaTcons0Ptr, ElinaDimchangePtr to the function')

    return tcons1


def elina_tcons0_remove_dimensions_with(tcons, dimchange):
    """
    Remove dimensions from an ElinaTcons0 by following the semantics of an ElinaDimchange.

    Parameters
    ----------
    tcons : ElinaTcons0Ptr
        Pointer to the ElinaTcons0 from which we want to remove dimensions.
    dimchange : ElinaDimchangePtr
        Pointer to the ElinaDimchange which semantics we want to follow.

    Returns
    -------
    None

    """

    try:
        elina_tcons0_remove_dimensions_with_c = elina_auxiliary_api.elina_tcons0_remove_dimensions_with
        elina_tcons0_remove_dimensions_with_c.restype = None
        elina_tcons0_remove_dimensions_with_c.argtypes = [ElinaTcons0Ptr, ElinaDimchangePtr]
        elina_tcons0_remove_dimensions_with_c(tcons, dimchange)
    except:
        print('Problem with loading/calling "elina_tcons0_remove_dimensions_with" from "libelinaux.so"')
        print('Make sure you are passing ElinaTcons0Ptr, ElinaDimchangePtr to the function')


def elina_tcons0_remove_dimensions(tcons2, dimchange):
    """
    Remove dimensions to an ElinaTcons0 by following the semantics of an ElinaDimchange.

    Parameters
    ----------
    tcons2 : ElinaTconsPtr
        Pointer to the ElinaTcons0 from which we want to remove dimensions.
    dimchange : ElinaDimchangePtr
        Pointer to the ElinaDimchange which semantics we want to follow.

    Returns
    -------
    tcons1 : ElinaTcons0Ptr
        Pointer to the newly created ElinaTcons0 with removed dimensions.

    """

    tcons1 = None
    try:
        elina_tcons0_remove_dimensions_c = elina_auxiliary_api.elina_tcons0_remove_dimensions
        elina_tcons0_remove_dimensions_c.restype = ElinaTcons0
        elina_tcons0_remove_dimensions_c.argtypes = [ElinaTcons0Ptr, ElinaDimchangePtr]
        tcons1 = elina_tcons0_remove_dimensions_c(tcons2, dimchange)
    except:
        print('Problem with loading/calling "elina_tcons0_remove_dimensions" from "libelinaux.so"')
        print('Make sure you are passing ElinaTcons0Ptr, ElinaDimchangePtr to the function')

    return tcons1


def elina_tcons0_permute_dimensions_with(tcons, dimperm):
    """
    Apply given permutation to the dimensions of an ElinaLincons0.

    Parameters
    ----------
    tcons : ElinaTcons0Ptr
        Pointer to the ElinaTcons0 which dimensions we want to permute.
    dimperm : ElinaDimpermPtr
        Pointer to the ElinaDimperm which describes the permutation.

    Returns
    -------
    None

    """

    try:
        elina_tcons0_permute_dimensions_with_c = elina_auxiliary_api.elina_tcons0_permute_dimensions_with
        elina_tcons0_permute_dimensions_with_c.restype = None
        elina_tcons0_permute_dimensions_with_c.argtypes = [ElinaTcons0Ptr, ElinaDimpermPtr]
        elina_tcons0_permute_dimensions_with_c(tcons, dimperm)
    except:
        print('Problem with loading/calling "elina_tcons0_permute_dimensions_with" from "libelinaux.so"')
        print('Make sure you are passing ElinaTcons0Ptr, ElinaDimpermPtr to the function')


def elina_tcons0_permute_dimensions(tcons2, dimperm):
    """
    Apply given permutation to the dimensions of an ElinaLincons0.

    Parameters
    ----------
    tcons2 : ElinaTcons0Ptr
        Pointer to the ElinaTcons0 which dimensions we want to permute.
    dimperm : ElinaDimpermPtr
        Pointer to the ElinaDimperm which describes the permutation.

    Returns
    -------
    tcons1 : ElinaTcons0
        Pointer to the newly created ElinaLincons0 with permuted dimensions.

    """

    tcons1 = None
    try:
        elina_tcons0_permute_dimensions_c = elina_auxiliary_api.elina_tcons0_permute_dimensions
        elina_tcons0_permute_dimensions_c.restype = ElinaTcons0
        elina_tcons0_permute_dimensions_c.argtypes = [ElinaTcons0Ptr, ElinaDimpermPtr]
        tcons1 = elina_tcons0_permute_dimensions_c(tcons2, dimperm)
    except:
        print('Problem with loading/calling "elina_tcons0_permute_dimensions" from "libelinaux.so"')
        print('Make sure you are passing ElinaTcons0Ptr, ElinaDimpermPtr to the function')

    return tcons1


# ********************************************************************** #
# II. Array of linear constraints
# ********************************************************************** #

def elina_tcons0_array_make(size):
    """
    Allocate an ElinaLincons0Array of a given size.
    The constraints are initialized with NULL pointers.

    Parameters
    ----------
    size : c_size_t
        Size of the ElinaTcons0Array.

    Returns
    -------
    tcons_array : ElinaTcons0Array
        Newly allocated ElinaTcons0Array.

    """

    tcons_array = None
    try:
        elina_tcons0_array_make_c = elina_auxiliary_api.elina_tcons0_array_make
        elina_tcons0_array_make_c.restype = ElinaTcons0Array
        elina_tcons0_array_make_c.argtypes = [c_size_t]
        tcons_array = elina_tcons0_array_make_c(size)
    except:
        print('Problem with loading/calling "elina_tcons0_array_make" from "libelinaux.so"')
        print('Make sure you are passing c_size_t to the function')

    return tcons_array


def elina_tcons0_array_resize(tcons_array, size):
    """
    Resize an ElinaTcons0Array to a different size.
    New constraints are initialized with NULL pointers.
    Removed constraints with non-NULL pointers are freed.

    Parameters
    ----------
    tcons_array : ElinaTcons0ArrayPtr
        Pointer to the ElinaTcons0Array that needs to be resized.
    size : c_size_t
        New size of the ElinaTcons0Array.

    Returns
    -------
    None

    """

    try:
        elina_tcons0_array_resize_c = elina_auxiliary_api.elina_tcons0_array_resize
        elina_tcons0_array_resize_c.restype = None
        elina_tcons0_array_resize_c.argtypes = [ElinaTcons0ArrayPtr, c_size_t]
        elina_tcons0_array_resize_c(tcons_array, size)
    except:
        print('Problem with loading/calling "elina_tcons0_array_resize" from "libelinaux.so"')
        print('Make sure you are passing ElinaTcons0ArrayPtr, c_size_t to the function')


def elina_tcons0_array_clear(tcons_array):
    """
    Free an ElinaTcons0Array.

    Parameters
    ----------
    tcons_array : ElinaTcons0ArrayPtr
        Pointer the ElinaTcons0Array that needs to be freed.

    Returns
    -------
    None

    """

    try:
        elina_tcons0_array_clear_c = elina_auxiliary_api.elina_tcons0_array_clear
        elina_tcons0_array_clear_c.restype = None
        elina_tcons0_array_clear_c.argtypes = [ElinaTcons0ArrayPtr]
        elina_tcons0_array_clear_c(tcons_array)
    except:
        print('Problem with loading/calling "elina_tcons0_array_clear" from "libelinaux.so"')
        print('Make sure you are passing ElinaTcons0ArrayPtr, c_size_t to the function')


def elina_tcons0_array_fprint(stream, tcons_array, name_of_dim):
    """
    Print ElinaTcons0Array to stream, having dimension names.

    Parameters
    ----------
    stream : c_void_p
        Stream onto which to print.
    tcons_array : ElinaTcons0ArrayPtr
        Pointer to the ElinaTcons0Array that needs to be printed.
    name_of_dim : POINTER(c_char_p)
        Dimension names.

    Returns
    -------
    None

    """

    try:
        elina_tcons0_array_fprint_c = elina_auxiliary_api.elina_tcons0_array_fprint
        elina_tcons0_array_fprint_c.restype = None
        elina_tcons0_array_fprint_c.argtypes = [c_void_p, ElinaTcons0ArrayPtr, POINTER(c_char_p)]
        elina_tcons0_array_fprint_c(stream, tcons_array, name_of_dim)
    except:
        print('Problem with loading/calling "elina_tcons0_array_fprint" from "libelinaux.so"')
        print('Make sure you are passing c_void_p, ElinaTcons0ArrayPtr, POINTER(c_char_p) to the function')


def elina_tcons0_array_is_interval_linear(tcons_array):
    """
    Test if the expressions involved in an ElinaTcons0Array are interval linear.

    Parameters
    ----------
    tcons_array : ElinaTcons0ArrayPtr
        Pointer to the ElinaTcons0Array that needs to be tested.

    Returns
    -------
    result : c_bool
        Result of the test.

    """

    result = None
    try:
        elina_tcons0_array_is_interval_linear_c = elina_auxiliary_api.elina_tcons0_array_is_interval_linear
        elina_tcons0_array_is_interval_linear_c.restype = c_bool
        elina_tcons0_array_is_interval_linear_c.argtypes = [ElinaTcons0ArrayPtr]
        result = elina_tcons0_array_is_interval_linear_c(tcons_array)
    except:
        print('Problem with loading/calling "elina_tcons0_array_is_interval_linear" from "libelinaux.so"')
        print('Make sure you are passing ElinaTcons0ArrayPtr to the function')

    return result


# ====================================================================== #
# II.1 Change of dimensions and permutations
# ====================================================================== #

def elina_tcons0_array_add_dimensions_with(tcons_array, dimchange):
    """
    Add dimensions to an ElinaTcons0Array by following the semantics of an ElinaDimchange.

    Parameters
    ----------
    tcons_array : ElinaTcons0ArrayPtr
        Pointer to the ElinaTcons0Array to which we want to add dimensions.
    dimchange : ElinaDimchangePtr
        Pointer to the ElinaDimchange which semantics we want to follow.

    Returns
    -------
    None

    """

    try:
        elina_tcons0_array_add_dimensions_with_c = elina_auxiliary_api.elina_tcons0_array_add_dimensions_with
        elina_tcons0_array_add_dimensions_with_c.restype = None
        elina_tcons0_array_add_dimensions_with_c.argtypes = [ElinaTcons0ArrayPtr, ElinaDimchangePtr]
        elina_tcons0_array_add_dimensions_with_c(tcons_array, dimchange)
    except:
        print('Problem with loading/calling "elina_tcons0_array_add_dimensions_with" from "libelinaux.so"')
        print('Make sure you are passing ElinaTcons0ArrayPtr, ElinaDimchangePtr to the function')


def elina_tcons0_array_add_dimensions(tcons_array2, dimchange):
    """
    Add dimensions to an ElinaTcons0Array by following the semantics of an ElinaDimchange.

    Parameters
    ----------
    tcons_array2 : ElinaTcons0ArrayPtr
        Pointer to the ElinaTcons0Array to which we want to add dimensions.
    dimchange : ElinaDimchangePtr
        Pointer to the ElinaDimchange which semantics we want to follow.

    Returns
    -------
    tcons_array1 : ElinaTcons0Array
        Newly created ElinaTcons0Array with added dimensions.

    """

    tcons_array1 = None
    try:
        elina_tcons0_array_add_dimensions_c = elina_auxiliary_api.elina_tcons0_array_add_dimensions
        elina_tcons0_array_add_dimensions_c.restype = ElinaTcons0Array
        elina_tcons0_array_add_dimensions_c.argtypes = [ElinaTcons0ArrayPtr, ElinaDimchangePtr]
        tcons_array1 = elina_tcons0_array_add_dimensions_c(tcons_array2, dimchange)
    except:
        print('Problem with loading/calling "elina_tcons0_array_add_dimensions_c" from "libelinaux.so"')
        print('Make sure you are passing ElinaTcons0ArrayPtr, ElinaDimchangePtr to the function')

    return tcons_array1


def elina_tcons0_array_remove_dimensions_with(tcons_array, dimchange):
    """
    Remove dimensions from an ElinaTcons0Array by following the semantics of an ElinaDimchange.

    Parameters
    ----------
    tcons_array : ElinaTcons0ArrayPtr
        Pointer to the ElinaTcons0Array from which we want to remove dimensions.
    dimchange : ElinaDimchangePtr
        Pointer to the ElinaDimchange which semantics we want to follow.

    Returns
    -------
    None

    """

    try:
        elina_tcons0_array_remove_dimensions_with_c = elina_auxiliary_api.elina_tcons0_array_remove_dimensions_with
        elina_tcons0_array_remove_dimensions_with_c.restype = None
        elina_tcons0_array_remove_dimensions_with_c.argtypes = [ElinaTcons0ArrayPtr, ElinaDimchangePtr]
        elina_tcons0_array_remove_dimensions_with_c(tcons_array, dimchange)
    except:
        print('Problem with loading/calling "elina_tcons0_array_remove_dimensions_with" from "libelinaux.so"')
        print('Make sure you are passing ElinaTcons0ArrayPtr, ElinaDimchangePtr to the function')


def elina_tcons0_array_remove_dimensions(tcons_array2, dimchange):
    """
    Remove dimensions from an ElinaTcons0Array by following the semantics of an ElinaDimchange.

    Parameters
    ----------
    tcons_array2 : ElinaTcons0ArrayPtr
        Pointer to the ElinaTcons0Array from which we want to remove dimensions.
    dimchange : ElinaDimchangePtr
        Pointer to the ElinaDimchange which semantics we want to follow.

    Returns
    -------
    tcons_array1 : ElinaTcons0Array
        Newly created ElinaTcons0Array with removed dimensions.

    """

    tcons_array1 = None
    try:
        elina_tcons0_array_remove_dimensions_c = elina_auxiliary_api.elina_tcons0_array_remove_dimensions
        elina_tcons0_array_remove_dimensions_c.restype = ElinaTcons0Array
        elina_tcons0_array_remove_dimensions_c.argtypes = [ElinaTcons0ArrayPtr, ElinaDimchangePtr]
        tcons_array1 = elina_tcons0_array_remove_dimensions_c(tcons_array2, dimchange)
    except:
        print('Problem with loading/calling "elina_tcons0_array_remove_dimensions" from "libelinaux.so"')
        print('Make sure you are passing ElinaTcons0ArrayPtr, ElinaDimchangePtr to the function')

    return tcons_array1


def elina_tcons0_array_permute_dimensions_with(tcons_array, dimperm):
    """
    Apply given permutation to the dimensions of an ElinaLincons0Array.

    Parameters
    ----------
    tcons_array : ElinaTcons0ArrayPtr
        Pointer to the ElinaTcons0Array which dimensions we want to permute.
    dimperm : ElinaDimpermPtr
        Pointer to the ElinaDimpermPtr which describes the permutation.

    Returns
    -------
    None

    """

    try:
        elina_tcons0_array_permute_dimensions_with_c = elina_auxiliary_api.elina_tcons0_array_permute_dimensions_with
        elina_tcons0_array_permute_dimensions_with_c.restype = None
        elina_tcons0_array_permute_dimensions_with_c.argtypes = [ElinaTcons0ArrayPtr, ElinaDimpermPtr]
        elina_tcons0_array_permute_dimensions_with_c(tcons_array, dimperm)
    except:
        print('Problem with loading/calling "elina_tcons0_array_permute_dimensions_with" from "libelinaux.so"')
        print('Make sure you are passing ElinaTcons0ArrayPtr, ElinaDimpermPtr to the function')


def elina_tcons0_array_permute_dimensions(tcons_array2, dimperm):
    """
    Apply given permutation to the dimensions of an ElinaLincons0Array.

    Parameters
    ----------
    tcons_array2 : ElinaTcons0ArrayPtr
        Pointer to the ElinaLincons0Array which dimensions we want to permute.
    dimperm : ElinaDimpermPtr
        Pointer to the ElinaDimpermPtr which describes the permutation.

    Returns
    -------
    tcons_array1 : ElinaTcons0Array
        Newly created ElinaTcons0Array with permuted dimensions.

    """

    tcons_array1 = None
    try:
        elina_tcons0_array_permute_dimensions_c = elina_auxiliary_api.elina_tcons0_array_permute_dimensions
        elina_tcons0_array_permute_dimensions_c.restype = ElinaTcons0Array
        elina_tcons0_array_permute_dimensions_c.argtypes = [ElinaTcons0ArrayPtr, ElinaDimpermPtr]
        tcons_array1 = elina_tcons0_array_permute_dimensions_c(tcons_array2, dimperm)
    except:
        print('Problem with loading/calling "elina_tcons0_array_permute_dimensions" from "libelinaux.so"')
        print('Make sure you are passing ElinaTcons0ArrayPtr, ElinaDimpermPtr to the function')

    return tcons_array1
