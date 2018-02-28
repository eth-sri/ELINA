from elina_lincons0_h import *


# ********************************************************************** #
# I. elina_lincons0_t
# ********************************************************************** #

# ====================================================================== #
# I.1 Memory management and printing
# ====================================================================== #

def elina_lincons0_make(constyp, linexpr, scalar):
    """
    Create an ElinaLincons0 of a type defined as ElinaConstyp given an ElinaLinexpr0.
    The expression and the coefficient are not duplicated, just pointed to.
    
    Parameters
    ----------
    constyp : c_uint
        Enum compatible with ElinaConstyp defining the type of the constraint.
    linexpr : ElinaLinexpr0Ptr
        Pointer to the ElinaLinexpr0 that is used for the constraint.
    scalar : ElinaScalarPtr
        Pointer to the ElinaScalar used as a modulo in case of constyp == EQMOD. Null otherwise.

    Returns
    -------
    lincons = ElinaLincons0
        Newly created ElinaLincons0.

    """

    lincons = None
    try:
        elina_lincons0_make_c = elina_auxiliary_api.elina_lincons0_make
        elina_lincons0_make_c.restype = ElinaLincons0
        elina_lincons0_make_c.argtypes = [c_uint, ElinaLinexpr0Ptr, ElinaScalarPtr]
        lincons = elina_lincons0_make_c(constyp, linexpr, scalar)
    except:
        print('Problem with loading/calling "elina_lincons0_make" from "libelinaux.so"')
        print('Make sure you are passing c_uint, ElinaLinexpr0Ptr, ElinaScalarPtr to the function')

    return lincons


def elina_lincons0_make_unsat():
    """
    Create the constraint -1 >= 0.
    
    Returns
    -------
    lincons : ElinaLincons0
        Newly created ElinaLincons0.

    """

    lincons = None
    try:
        elina_lincons0_make_unsat_c = elina_auxiliary_api.elina_lincons0_make_unsat
        elina_lincons0_make_unsat_c.restype = ElinaLincons0
        lincons = elina_lincons0_make_unsat_c()
    except:
        print('Problem with loading/calling "elina_lincons0_make_unsat" from "libelinaux.so"')

    return lincons


def elina_lincons0_copy(lincons2):
    """
    Duplicate an ElinaLincons0.
    
    Parameters
    ----------
    lincons2 : ElinaLincons0Ptr
        Pointer to the ElinaLincons0 that needs to be duplicated.

    Returns
    -------
    lincons1 : ElinaLincons0
        The duplicate ElinaLincons0.

    """

    lincons1 = None
    try:
        elina_lincons0_copy_c = elina_auxiliary_api.elina_lincons0_copy
        elina_lincons0_copy_c.restype = ElinaLincons0
        elina_lincons0_copy_c.argtypes = [ElinaLincons0Ptr]
        lincons1 = elina_lincons0_copy_c(lincons2)
    except:
        print('Problem with loading/calling "elina_lincons0_copy" from "libelinaux.so"')
        print('Make sure you are passing ElinaLincons0Ptr to the function')

    return lincons1


def elina_lincons0_clear(lincons):
    """
    Free an ElinaLincons0 and set the pointer to NULL.
    
    Parameters
    ----------
    lincons : ElinaLincons0Ptr
        Pointer to the ElinaLincons0 that needs to be freed.

    Returns
    -------
    None

    """

    try:
        elina_lincons0_clear_c = elina_auxiliary_api.elina_lincons0_clear
        elina_lincons0_clear_c.restype = None
        elina_lincons0_clear_c.argtypes = [ElinaLincons0Ptr]
        elina_lincons0_clear_c(lincons)
    except:
        print('Problem with loading/calling "elina_lincons0_clear" from "libelinaux.so"')
        print('Make sure you are passing ElinaLincons0Ptr to the function')


def elina_lincons0_print(lincons, name_of_dim):
    """
    Print an ElinaLincons0, having dimension names to stdout.
    
    Parameters
    ----------
    lincons : ElinaLincons0Ptr
        Pointer to the ElinaLincons0 that needs to be printed.
    name_of_dim : POINTER(c_char_p)
        Dimension names.

    Returns
    -------
    None

    """

    try:
        elina_lincons0_print_c = elina_auxiliary_api.elina_lincons0_print
        elina_lincons0_print_c.restype = None
        elina_lincons0_print_c.argtypes = [ElinaLincons0Ptr, POINTER(c_char_p)]
        elina_lincons0_print_c(lincons, name_of_dim)
    except:
        print('Problem with loading/calling "elina_lincons0_print" from "libelinaux.so"')
        print('Make sure you are passing ElinaLincons0Ptr, POINTER(c_char_p) to the function')


def elina_lincons0_fprint(stream, lincons, name_of_dim):
    """
    Print an ElinaLincons0, having dimension names to stream.
    
    Parameters
    ----------
    stream : c_void_p
        Stream onto which to print.
    lincons : ElinaLincons0Ptr
        Pointer to the ElinaLincons0 that needs to be printed.
    name_of_dim : POINTER(c_char_p)
        Dimension names.
    
    Returns
    -------
    None

    """

    try:
        elina_lincons0_fprint_c = elina_auxiliary_api.elina_lincons0_fprint
        elina_lincons0_fprint_c.restype = None
        elina_lincons0_fprint_c.argtypes = [c_void_p, ElinaLincons0Ptr, POINTER(c_char_p)]
        elina_lincons0_fprint_c(stream, lincons, name_of_dim)
    except:
        print('Problem with loading/calling "elina_lincons0_fprint" from "libelinaux.so"')
        print('Make sure you are passing ElinaLincons0Ptr, POINTER(c_char_p) to the function')


# ====================================================================== #
# I.2 Tests
# ====================================================================== #

def elina_lincons0_is_unsat(lincons):
    """
    Test if the constraint is b >= 0 or [a, b] >= 0 with b < 0.
    
    Parameters
    ----------
    lincons : ElinaLincons0Ptr
        Pointer to the ElinaLincons0 that needs to be tested.

    Returns
    -------
    result : c_bool
        Result of the test.

    """

    result = None
    try:
        elina_lincons0_is_unsat_c = elina_auxiliary_api.elina_lincons0_is_unsat
        elina_lincons0_is_unsat_c.restype = c_bool
        elina_lincons0_is_unsat_c.argtypes = [ElinaLincons0Ptr]
        result = elina_lincons0_is_unsat_c(lincons)
    except:
        print('Problem with loading/calling "elina_lincons0_is_unsat" from "libelinaux.so"')
        print('Make sure you are passing ElinaLincons0Ptr to the function')

    return result


def elina_lincons0_is_sat(lincons):
    """
    Test if the constraint is trivially satisfiable, e.g. [a,b] >= 0 with a > 0.
    
    Parameters
    ----------
    lincons : ElinaLincons0Ptr
        Pointer to the ElinaLincons0 that needs to be tested.

    Returns
    -------
    result = c_bool
        Result of the test.

    """

    result = None
    try:
        elina_lincons0_is_sat_c = elina_auxiliary_api.elina_lincons0_is_sat
        elina_lincons0_is_sat_c.restype = c_bool
        elina_lincons0_is_sat_c.argtypes = [ElinaLincons0Ptr]
        result = elina_lincons0_is_sat_c(lincons)
    except:
        print('Problem with loading/calling "elina_lincons0_is_sat" from "libelinaux.so"')
        print('Make sure you are passing ElinaLincons0Ptr to the function')

    return result


# ====================================================================== #
# I.3 Change of dimensions and permutations
# ====================================================================== #

def elina_lincons0_add_dimensions_with(lincons, dimchange):
    """
    Add dimensions to an ElinaLincons0 by following the semantics of an ElinaDimchange.

    Parameters
    ----------
    lincons : ElinaLincons0Ptr
        Pointer to the ElinaLincons0 to which we want to add dimensions.
    dimchange : ElinaDimchangePtr
        Pointer to the ElinaDimchange which semantics we want to follow.

    Returns
    -------
    None

    """

    try:
        elina_lincons0_add_dimensions_with_c = elina_auxiliary_api.elina_lincons0_add_dimensions_with
        elina_lincons0_add_dimensions_with_c.restype = None
        elina_lincons0_add_dimensions_with_c.argtypes = [ElinaLincons0Ptr, ElinaDimchangePtr]
        elina_lincons0_add_dimensions_with_c(lincons, dimchange)
    except:
        print('Problem with loading/calling "elina_lincons0_add_dimensions_with" from "libelinaux.so"')
        print('Make sure you are passing ElinaLincons0Ptr, ElinaDimchangePtr to the function')


def elina_lincons0_add_dimensions(lincons2, dimchange):
    """
    Add dimensions to an ElinaLincons0 by following the semantics of an ElinaDimchange.

    Parameters
    ----------
    lincons2 : ElinaLincons0Ptr
        Pointer to the ElinaLincons0 to which we want to add dimensions.
    dimchange : ElinaDimchangePtr
        Pointer to the ElinaDimchange which semantics we want to follow.

    Returns
    -------
    lincons1 : ElinaLincons0Ptr
        Pointer to the newly created ElinaLincons0 with added dimensions.

    """

    lincons1 = None
    try:
        elina_lincons0_add_dimensions_c = elina_auxiliary_api.elina_lincons0_add_dimensions
        elina_lincons0_add_dimensions_c.restype = ElinaLincons0Ptr
        elina_lincons0_add_dimensions_c.argtypes = [ElinaLincons0Ptr, ElinaDimchangePtr]
        lincons1 = elina_lincons0_add_dimensions_c(lincons2, dimchange)
    except:
        print('Problem with loading/calling "elina_lincons0_add_dimensions" from "libelinaux.so"')
        print('Make sure you are passing ElinaLincons0Ptr, ElinaDimchangePtr to the function')

    return lincons1


def elina_lincons0_permute_dimensions_with(lincons, perm):
    """
    Apply given permutation to the dimensions of an ElinaLincons0.

    Parameters
    ----------
    lincons : ElinaLincons0Ptr
        Pointer to the ElinaLincons0 which dimensions we want to permute.
    perm : ElinaDimpermPtr
        Pointer to the ElinaDimperm which describes the permutation.

    Returns
    -------
    None

    """

    try:
        elina_lincons0_permute_dimensions_with_c = elina_auxiliary_api.elina_lincons0_permute_dimensions_with
        elina_lincons0_permute_dimensions_with_c.restype = None
        elina_lincons0_permute_dimensions_with_c.argtypes = [ElinaLincons0Ptr, ElinaDimpermPtr]
        elina_lincons0_permute_dimensions_with_c(lincons, perm)
    except:
        print('Problem with loading/calling "elina_lincons0_permute_dimensions_with" from "libelinaux.so"')
        print('Make sure you are passing ElinaLincons0Ptr, ElinaDimpermPtr to the function')


def elina_lincons0_permute_dimensions(lincons2, perm):
    """
    Apply given permutation to the dimensions of an ElinaLincons0.

    Parameters
    ----------
    lincons2 : ElinaLincons0Ptr
        Pointer to the ElinaLincons0 which dimensions we want to permute.
    perm : ElinaDimpermPtr
        Pointer to the ElinaDimperm which describes the permutation.

    Returns
    -------
    lincons1 : ElinaLincons0Ptr
        Pointer to the newly created ElinaLincons0 with permuted dimensions.

    """

    lincons1 = None
    try:
        elina_lincons0_permute_dimensions_c = elina_auxiliary_api.elina_lincons0_permute_dimensions
        elina_lincons0_permute_dimensions_c.restype = ElinaLincons0Ptr
        elina_lincons0_permute_dimensions_c.argtypes = [ElinaLincons0Ptr, ElinaDimpermPtr]
        lincons1 = elina_lincons0_permute_dimensions_c(lincons2, perm)
    except:
        print('Problem with loading/calling "elina_lincons0_permute_dimensions" from "libelinaux.so"')
        print('Make sure you are passing ElinaLincons0Ptr, ElinaDimpermPtr to the function')

    return lincons1


# ********************************************************************** #
# II. Array of linear constraints
# ********************************************************************** #

def elina_lincons0_array_make(size):
    """
    Allocate an ElinaLincons0Array of a given size.
    The constraints are initialized with NULL pointers.
    
    Parameters
    ----------
    size : c_size_t
        Size of the ElinaLincons0Array.

    Returns
    -------
    lincons_array : ElinaLincons0Array
        Newly allocated ElinaLincons0Array.

    """

    lincons_array = None
    try:
        elina_lincons0_array_make_c = elina_auxiliary_api.elina_lincons0_array_make
        elina_lincons0_array_make_c.restype = ElinaLincons0Array
        elina_lincons0_array_make_c.argtypes = [c_size_t]
        lincons_array = elina_lincons0_array_make_c(size)
    except:
        print('Problem with loading/calling "elina_lincons0_array_make" from "libelinaux.so"')
        print('Make sure you are passing c_size_t to the function')

    return lincons_array


def elina_lincons0_array_resize(lincons_array, size):
    """
    Resize an ElinaLincons0Array to a different size.
    New constraints are initialized with NULL pointers.
    Removed constraints with non-NULL pointers are freed.
    
    Parameters
    ----------
    lincons_array : ElinaLincons0ArrayPtr
        Pointer to the ElinaLincons0Array that needs to be resized.
    size : c_size_t
        New size of the ElinaLincons0Array.

    Returns
    -------
    None
    
    """

    try:
        elina_lincons0_array_resize_c = elina_auxiliary_api.elina_lincons0_array_resize
        elina_lincons0_array_resize_c.restype = None
        elina_lincons0_array_resize_c.argtypes = [ElinaLincons0ArrayPtr, c_size_t]
        elina_lincons0_array_resize_c(lincons_array, size)
    except:
        print('Problem with loading/calling "elina_lincons0_array_resize" from "libelinaux.so"')
        print('Make sure you are passing ElinaLincons0ArrayPtr, c_size_t to the function')


def elina_lincons0_array_clear(lincons_array):
    """
    Free an ElinaLincons0Array.
    
    Parameters
    ----------
    lincons_array : ElinaLincons0ArrayPtr
        Pointer the ElinaLincons0Array that needs to be freed.

    Returns
    -------
    None
    
    """

    try:
        elina_lincons0_array_clear_c = elina_auxiliary_api.elina_lincons0_array_clear
        elina_lincons0_array_clear_c.restype = None
        elina_lincons0_array_clear_c.argtypes = [ElinaLincons0ArrayPtr]
        elina_lincons0_array_clear_c(lincons_array)
    except:
        print('Problem with loading/calling "elina_lincons0_array_clear" from "libelinaux.so"')
        print('Make sure you are passing ElinaLincons0ArrayPtr to the function')


def elina_lincons0_array_print(lincons_array, name_of_dim):
    """
    Print an ElinaLincons0Array to stdout, having dimension names.
    
    Parameters
    ----------
    lincons_array : ElinaLincons0ArrayPtr
        Pointer to the ElinaLincons0Array that needs to be printed.
    name_of_dim : POINTER(c_char_p)
        Dimension names.
        
    Returns
    -------
    None

    """

    try:
        elina_lincons0_array_print_c = elina_auxiliary_api.elina_lincons0_array_print
        elina_lincons0_array_print_c.restype = None
        elina_lincons0_array_print_c.argtypes = [ElinaLincons0ArrayPtr, POINTER(c_char_p)]
        elina_lincons0_array_print_c(lincons_array, name_of_dim)
    except:
        print('Problem with loading/calling "elina_lincons0_array_print" from "libelinaux.so"')
        print('Make sure you are passing ElinaLincons0ArrayPtr, POINTER(c_char_p) to the function')


def elina_lincons0_array_fprint(stream, lincons_array, name_of_dim):
    """
    Print ElinaLincons0Array to stream, having dimension names.
    
    Parameters
    ----------
    stream : c_void_p
        Stream onto which to print.
    lincons_array : ElinaLincons0ArrayPtr
        Pointer to the ElinaLincons0Array that needs to be printed.
    name_of_dim : POINTER(c_char_p)
        Dimension names.

    Returns
    -------
    None

    """

    try:
        elina_lincons0_array_fprint_c = elina_auxiliary_api.elina_lincons0_array_fprint
        elina_lincons0_array_fprint_c.restype = None
        elina_lincons0_array_fprint_c.argtypes = [c_void_p, ElinaLincons0ArrayPtr, POINTER(c_char_p)]
        elina_lincons0_array_fprint_c(stream, lincons_array, name_of_dim)
    except:
        print('Problem with loading/calling "elina_lincons0_array_fprint" from "libelinaux.so"')
        print('Make sure you are passing c_void_p, ElinaLincons0ArrayPtr, POINTER(c_char_p) to the function')


def elina_lincons0_array_type(lincons_array):
    """
    Checks the type of an ElinaLincons0Array.
    
    Parameters
    ----------
    lincons_array : ElinaLincons0ArrayPtr
        Pointer to the ElinaLincons0Array that needs to be checked for its type.

    Returns
    -------
    array_type : c_uint
        Type of the ElinaLincons0Array according to ElinaLinconstyp.

    """

    array_type = None
    try:
        elina_lincons0_array_type_c = elina_auxiliary_api.elina_lincons0_array_type
        elina_lincons0_array_type_c.restype = c_uint
        elina_lincons0_array_type_c.argtypes = [ElinaLincons0ArrayPtr]
        array_type = elina_lincons0_array_type_c(lincons_array)
    except:
        print('Problem with loading/calling "elina_lincons0_array_type" from "libelinaux.so"')
        print('Make sure you are passing ElinaLincons0ArrayPtr to the function')

    return array_type


def elina_lincons0_array_is_linear(lincons_array):
    """
    Test if the expressions involved in an ElinaLincons0Array are linear.
    
    Parameters
    ----------
    lincons_array : ElinaLincons0ArrayPtr
        Pointer to the ElinaLincons0Array that needs to be tested.

    Returns
    -------
    result : c_bool
        Result of the test.
        
    """

    result = None
    try:
        elina_lincons0_array_is_linear_c = elina_auxiliary_api.elina_lincons0_array_is_linear
        elina_lincons0_array_is_linear_c.restype = c_bool
        elina_lincons0_array_is_linear_c.argtypes = [ElinaLincons0ArrayPtr]
        result = elina_lincons0_array_is_linear_c(lincons_array)
    except:
        print('Problem with loading/calling "elina_lincons0_array_is_linear" from "libelinaux.so"')
        print('Make sure you are passing ElinaLincons0ArrayPtr to the function')

    return result


def elina_lincons0_array_is_quasilinear(lincons_array):
    """
    Test if the expressions involved in an ElinaLincons0Array are quasilinear.
    
    Parameters
    ----------
    lincons_array : ElinaLincons0ArrayPtr
        Pointer to the ElinaLincons0Array that needs to be tested.

    Returns
    -------
    result : c_bool
        Result of the test.

    """

    result = None
    try:
        elina_lincons0_array_is_quasilinear_c = elina_auxiliary_api.elina_lincons0_array_is_quasilinear
        elina_lincons0_array_is_quasilinear_c.restype = c_bool
        elina_lincons0_array_is_quasilinear_c.argtypes = [ElinaLincons0ArrayPtr]
        result = elina_lincons0_array_is_quasilinear_c(lincons_array)
    except:
        print('Problem with loading/calling "elina_lincons0_array_is_quasilinear" from "libelinaux.so"')
        print('Make sure you are passing ElinaLincons0ArrayPtr to the function')

    return result


# ====================================================================== #
# II.1 Change of dimensions and permutations
# ====================================================================== #

def elina_lincons0_array_add_dimensions_with(lincons_array, dimchange):
    """
    Add dimensions to an ElinaLincons0Array by following the semantics of an ElinaDimchange.

    Parameters
    ----------
    lincons_array : ElinaLincons0ArrayPtr
        Pointer to the ElinaLincons0Array to which we want to add dimensions.
    dimchange : ElinaDimchangePtr
        Pointer to the ElinaDimchange which semantics we want to follow.

    Returns
    -------
    None

    """

    try:
        elina_lincons0_array_add_dimensions_with_c = elina_auxiliary_api.elina_lincons0_array_add_dimensions_with
        elina_lincons0_array_add_dimensions_with_c.restype = None
        elina_lincons0_array_add_dimensions_with_c.argtypes = [ElinaLincons0ArrayPtr, ElinaDimchangePtr]
        elina_lincons0_array_add_dimensions_with_c(lincons_array, dimchange)
    except:
        print('Problem with loading/calling "elina_lincons0_array_add_dimensions_with" from "libelinaux.so"')
        print('Make sure you are passing ElinaLincons0ArrayPtr, ElinaDimchangePtr to the function')


def elina_lincons0_array_add_dimensions(lincons_array2, dimchange):
    """
    Add dimensions to an ElinaLincons0Array by following the semantics of an ElinaDimchange.
    
    Parameters
    ----------
    lincons_array2 : ElinaLincons0ArrayPtr
        Pointer to the ElinaLincons0Array to which we want to add dimensions.
    dimchange : ElinaDimchangePtr
        Pointer to the ElinaDimchange which semantics we want to follow.

    Returns
    -------
    lincons_array1 : ElinaLincons0Array
        Newly created ElinaLincons0Array with added dimensions.

    """

    lincons_array1 = None
    try:
        elina_lincons0_array_add_dimensions_c = elina_auxiliary_api.elina_lincons0_array_add_dimensions
        elina_lincons0_array_add_dimensions_c.restype = ElinaLincons0Array
        elina_lincons0_array_add_dimensions_c.argtypes = [ElinaLincons0ArrayPtr, ElinaDimchangePtr]
        lincons_array1 = elina_lincons0_array_add_dimensions_c(lincons_array2, dimchange)
    except:
        print('Problem with loading/calling "elina_lincons0_array_add_dimensions" from "libelinaux.so"')
        print('Make sure you are passing ElinaLincons0ArrayPtr, ElinaDimchangePtr to the function')

    return lincons_array1


def elina_lincons0_array_permute_dimensions_with(lincons_array, perm):
    """
    Apply given permutation to the dimensions of an ElinaLincons0Array.

    Parameters
    ----------
    lincons_array : ElinaLincons0ArrayPtr
        Pointer to the ElinaLincons0Array which dimensions we want to permute.
    perm : ElinaDimpermPtr
        Pointer to the ElinaDimpermPtr which describes the permutation.

    Returns
    -------
    None
    
    """

    try:
        elina_lincons0_array_permute_dimensions_with_c = \
            elina_auxiliary_api.elina_lincons0_array_permute_dimensions_with
        elina_lincons0_array_permute_dimensions_with_c.restype = None
        elina_lincons0_array_permute_dimensions_with_c.argtypes = [ElinaLincons0ArrayPtr, ElinaDimpermPtr]
        elina_lincons0_array_permute_dimensions_with_c(lincons_array, perm)
    except:
        print('Problem with loading/calling "elina_lincons0_array_permute_dimensions_with" from "libelinaux.so"')
        print('Make sure you are passing ElinaLincons0ArrayPtr, ElinaDimpermPtr to the function')


def elina_lincons0_array_permute_dimensions(lincons_array2, perm):
    """
    Apply given permutation to the dimensions of an ElinaLincons0Array.

    Parameters
    ----------
    lincons_array2 : ElinaLincons0ArrayPtr
        Pointer to the ElinaLincons0Array which dimensions we want to permute.
    perm : ElinaDimpermPtr
        Pointer to the ElinaDimpermPtr which describes the permutation.

    Returns
    -------
    lincons_array1 : ElinaLincons0Array
        Newly created ElinaLincons0Array with permuted dimensions.

    """

    lincons_array1 = None
    try:
        elina_lincons0_array_permute_dimensions_c = elina_auxiliary_api.elina_lincons0_array_permute_dimensions
        elina_lincons0_array_permute_dimensions_c.restype = ElinaLincons0Array
        elina_lincons0_array_permute_dimensions_c.argtypes = [ElinaLincons0ArrayPtr, ElinaDimpermPtr]
        lincons_array1 = elina_lincons0_array_permute_dimensions_c(lincons_array2, perm)
    except:
        print('Problem with loading/calling "elina_lincons0_array_permute_dimensions" from "libelinaux.so"')
        print('Make sure you are passing ElinaLincons0ArrayPtr, ElinaDimpermPtr to the function')

    return lincons_array1
