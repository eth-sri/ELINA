from elina_abstract0_h import *


# ********************************************************************** #
# I. General management
# ********************************************************************** #

# ============================================================ #
# I.1 Memory
# ============================================================ #

def elina_abstract0_copy(man, a1):
    """
    Return a copy of an ElinaAbstract0.
    Destructive update does not affect the initial value.
    
    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    a1 : ElinaAbstract0Ptr
        Pointer to the ElinaAbstract0.
        
    Returns
    -------
    a0 : ElinaAbstract0Ptr
        Pointer to the newly created copy ElinaAbstract0.

    """

    a0 = None
    try:
        elina_abstract0_copy_c = elina_auxiliary_api.elina_abstract0_copy
        elina_abstract0_copy_c.restype = ElinaAbstract0Ptr
        elina_abstract0_copy_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr]
        a0 = elina_abstract0_copy_c(man, a1)
    except:
        print('Problem with loading/calling "elina_abstract0_copy" from "libelinaux.so"')
        print('Make sure you are passing ElinaManagerPtr, ElinaAbstract0Ptr to the function')

    return a0


def elina_abstract0_free(man, a):
    """
    Free all the memory used by the value of an ElinaAbstract0.
    
    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    a : ElinaAbstract0Ptr
        Pointer to the ElinaAbstract0 which memory needs to be freed.
        
    Returns
    -------
    None
    
    """

    try:
        elina_abstract0_free_c = elina_auxiliary_api.elina_abstract0_free
        elina_abstract0_free_c.restype = None
        elina_abstract0_free_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr]
        elina_abstract0_free_c(man, a)
    except:
        print('Problem with loading/calling "elina_abstract0_free" from "libelinaux.so"')
        print('Make sure you are passing ElinaManagerPtr, ElinaAbstract0Ptr to the function')


def elina_abstract0_size(man, a):
    """
    Return the abstract size of the value in an ElinaAbstract0.
    
    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    a : ElinaAbstract0Ptr
        Pointer to the ElinaAbstract0 for which we want to obtain the size of the value.

    Returns
    -------
    size : c_size_t
        The size of the value in an ElinaAbstract0.

    """

    size = None
    try:
        elina_abstract0_size_c = elina_auxiliary_api.elina_abstract0_size
        elina_abstract0_size_c.restype = c_size_t
        elina_abstract0_size_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr]
        size = elina_abstract0_size_c(man, a)
    except:
        print('Problem with loading/calling "elina_abstract0_size" from "libelinaux.so"')
        print('Make sure you are passing ElinaManagerPtr, ElinaAbstract0Ptr to the function')

    return size


# ============================================================ #
# I.2 Control of internal representation
# ============================================================ #

def elina_abstract0_minimize(man, a):
    """
    Minimize the size of an ElinaAbstract0.
    This may result in a later recomputation of internal information.
    
    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    a : ElinaAbstract0Ptr
        Pointer to the ElinaAbstract0 which size we want to minimize.

    Returns
    -------
    None

    """

    try:
        elina_abstract0_minimize_c = elina_auxiliary_api.elina_abstract0_minimize
        elina_abstract0_minimize_c.restype = None
        elina_abstract0_minimize_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr]
        elina_abstract0_minimize_c(man, a)
    except:
        print('Problem with loading/calling "elina_abstract0_minimize" from "libelinaux.so"')
        print('Make sure you are passing ElinaManagerPtr, ElinaAbstract0Ptr to the function')


def elina_abstract0_canonicalize(man, a):
    """
    Put an ElinaAbstract0 in canonical form. (not a clear definition yet)
    
    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    a : ElinaAbstract0Ptr
        Pointer to the ElinaAbstract0 we want to canonicalize.
    
    Returns
    -------
    None
    
    """

    try:
        elina_abstract0_canonicalize_c = elina_auxiliary_api.elina_abstract0_canonicalize
        elina_abstract0_canonicalize_c.restype = None
        elina_abstract0_canonicalize_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr]
        elina_abstract0_canonicalize_c(man, a)
    except:
        print('Problem with loading/calling "elina_abstract0_canonicalize" from "libelinaux.so"')
        print('Make sure you are passing ElinaManagerPtr, ElinaAbstract0Ptr to the function')


def elina_abstract0_hash(man, a):
    """
    Return a hash of an ELinaAbstract0.
    Two ElinaAbstract0 in canonical form, if considered as equal by elina_abstract0_is_eq, should receive the same hash.
    
    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    a : ElinaAbstract0Ptr
        Pointer to the ElinaAbstract0 we want to hash.

    Returns
    -------
    result : c_int
        The resulting hash.

    """

    result = None
    try:
        elina_abstract0_hash_c = elina_auxiliary_api.elina_abstract0_hash
        elina_abstract0_hash_c.restype = c_int
        elina_abstract0_hash_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr]
        result = elina_abstract0_hash_c(man, a)
    except:
        print('Problem with loading/calling "elina_abstract0_hash" from "libelinaux.so"')
        print('Make sure you are passing ElinaManagerPtr, ElinaAbstract0Ptr to the function')

    return result


def elina_abstract0_approximate(man, a, algorithm):
    """
    Perform some transformation on an ElinaAbstract0.
    The transformation may result in an information loss.
    
    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    a : ElinaAbstract0Ptr
        Pointer to the ElinaAbstract0 we want to transform.
    algorithm : c_int
        Integer specifying the transformation.

    Returns
    -------
    None
    
    """

    try:
        elina_abstract0_approximate_c = elina_auxiliary_api.elina_abstract0_approximate
        elina_abstract0_approximate_c.restype = None
        elina_abstract0_approximate_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr, c_int]
        elina_abstract0_approximate_c(man, a, algorithm)
    except:
        print('Problem with loading/calling "elina_abstract0_approximate" from "libelinaux.so"')
        print('Make sure you are passing ElinaManagerPtr, ElinaAbstract0Ptr, c_int to the function')


# ============================================================ #
# I.3 Printing
# ============================================================ #

def elina_abstract0_fprint(stream, man, a, name_of_dim):
    """
    Print an ElinaAbstract0 to stream having dimension names.
    
    Parameters
    ----------
    stream : c_void_p
        Stream onto which to print.
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    a : ElinaAbstract0Ptr
        Pointer to the ElinaAbstract0 we want to print.
    name_of_dim : POINTER(c_char_p)
        Dimension names.

    Returns
    -------
    None

    """

    try:
        elina_abstract0_fprint_c = elina_auxiliary_api.elina_abstract0_fprint
        elina_abstract0_fprint_c.restype = None
        elina_abstract0_fprint_c.argtypes = [c_void_p, ElinaManagerPtr, ElinaAbstract0Ptr, POINTER(c_char_p)]
        elina_abstract0_fprint_c(stream, man, a, name_of_dim)
    except:
        print('Problem with loading/calling "elina_abstract0_fprint" from "libelinaux.so"')
        print('Make sure you are passing c_void_p, ElinaManagerPtr, ElinaAbstract0Ptr, '
              'POINTER(c_char_p) to the function')


def elina_abstract0_fprintdiff(stream, man, a1, a2, name_of_dim):
    """
    Print the difference between two ElinaAbstract0-s to stream having dimension names.
    The meaning of difference is library dependent.

    Parameters
    ----------
    stream : c_void_p
        Stream onto which to print.
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    a1 : ElinaAbstract0Ptr
        Pointer to the first ElinaAbstract0.
    a2 : ElinaAbstract0Ptr
        Pointer to the second ElinaAbstract0.
    name_of_dim : POINTER(c_char_p)
        Dimension names.

    Returns
    -------
    None

    """

    try:
        elina_abstract0_fprintdiff_c = elina_auxiliary_api.elina_abstract0_fprintdiff
        elina_abstract0_fprintdiff_c.restype = None
        elina_abstract0_fprintdiff_c.argtypes = \
            [c_void_p, ElinaManagerPtr, ElinaAbstract0Ptr, ElinaAbstract0Ptr, POINTER(c_char_p)]
        elina_abstract0_fprintdiff_c(stream, man, a1, a2, name_of_dim)
    except:
        print('Problem with loading/calling "elina_abstract0_fprintdiff" from "libelinaux.so"')
        print('Make sure you are passing c_void_p, ElinaManagerPtr, ElinaAbstract0Ptr, ElinaAbstract0Ptr, '
              'POINTER(c_char_p) to the function')


def elina_abstract0_fdump(stream, man, a):
    """
    Print the internal representation of an ElinaAbstract0 to stream.
    Used for debugging purposes.

    Parameters
    ----------
    stream : c_void_p
        Stream onto which to print.
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    a : ElinaAbstract0Ptr
        Pointer to the ElinaAbstract0 which internal representation we want to print.

    Returns
    -------
    None

    """

    try:
        elina_abstract0_fdump_c = elina_auxiliary_api.elina_abstract0_fdump
        elina_abstract0_fdump_c.restype = None
        elina_abstract0_fdump_c.argtypes = [c_void_p, ElinaManagerPtr, ElinaAbstract0Ptr]
        elina_abstract0_fdump_c(stream, man, a)
    except:
        print('Problem with loading/calling "elina_abstract0_fdump" from "libelinaux.so"')
        print('Make sure you are passing c_void_p, ElinaManagerPtr, ElinaAbstract0Ptr to the function')


# ============================================================ #
# I.4 Serialization
# ============================================================ #

def elina_abstract0_serialize_raw(man, a):
    """
    Allocate an ElinaMembuf and output an ElinaAbstract0 in raw binary format to it.
    It is the user responsibility to free the memory afterwards (with free).
    
    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    a : ElinaAbstract0Ptr
        Pointer to the ElinaAbstract0 we want to output to an ElinaMembuf.

    Returns
    -------
    membuf : ElinaMembuf
        The resulting ElinaMembuf.

    """

    membuf = None
    try:
        elina_abstract0_serialize_raw_c = elina_auxiliary_api.elina_abstract0_serialize_raw
        elina_abstract0_serialize_raw_c.restype = ElinaMembuf
        elina_abstract0_serialize_raw_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr]
        membuf = elina_abstract0_serialize_raw_c(man, a)
    except:
        print('Problem with loading/calling "elina_abstract0_serialize_raw" from "libelinaux.so"')
        print('Make sure you are passing ElinaManagerPtr, ElinaAbstract0Ptr to the function')

    return membuf


def elina_abstract0_deserialize_raw(man, ptr, size):
    """
    Create an ElinaAbstract0 from raw binary format read from an input stream.
    
    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    ptr : c_void_p
        Stream from which to read.
    size : POINTER(c_size_t)
        Store the number of bytes read in size. 
    
    Returns
    -------
    a : ElinaAbstract0Ptr
        Pointer to the resulting ElinaAbstract0.
        
    """

    a = None
    try:
        elina_abstract0_deserialize_raw_c = elina_auxiliary_api.elina_abstract0_deserialize_raw
        elina_abstract0_deserialize_raw_c.restype = ElinaAbstract0Ptr
        elina_abstract0_deserialize_raw_c.argtypes = [ElinaManagerPtr, c_void_p, POINTER(c_size_t)]
        a = elina_abstract0_deserialize_raw_c(man, ptr, size)
    except:
        print('Problem with loading/calling "elina_abstract0_deserialize_raw" from "libelinaux.so"')
        print('Make sure you are passing ElinaManagerPtr, ElinaAbstract0Ptr to the function')

    return a


# ********************************************************************** #
# II. Constructor, accessors, tests and property extraction
# ********************************************************************** #

# ============================================================ #
# II.1 Basic constructors
# ============================================================ #

def elina_abstract0_bottom(man, intdim, realdim):
    """
    Create a bottom (empty) ElinaAbstract0.

    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    intdim : c_size_t
        Number of integer variables.
    realdim : c_size_t
        Number of real variables.
    
    Returns
    -------
    a : ElinaAbstract0Ptr
        Pointer to the resulting ElinaAbstract0.

    """

    a = None
    try:
        elina_abstract0_bottom_c = elina_auxiliary_api.elina_abstract0_bottom
        elina_abstract0_bottom_c.restype = ElinaAbstract0Ptr
        elina_abstract0_bottom_c.argtypes = [ElinaManagerPtr, c_size_t, c_size_t]
        a = elina_abstract0_bottom_c(man, intdim, realdim)
    except:
        print('Problem with loading/calling "elina_abstract0_bottom" from "libelinaux.so"')
        print('Make sure you are passing EElinaManagerPtr, c_size_t, c_size_t to the function')

    return a


def elina_abstract0_top(man, intdim, realdim):
    """
    Create a top (universe) ElinaAbstract0.

    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    intdim : c_size_t
        Number of integer variables.
    realdim : c_size_t
        Number of real variables.

    Returns
    -------
    a : ElinaAbstract0Ptr
        Pointer to the resulting ElinaAbstract0.

    """

    a = None
    try:
        elina_abstract0_top_c = elina_auxiliary_api.elina_abstract0_top
        elina_abstract0_top_c.restype = ElinaAbstract0Ptr
        elina_abstract0_top_c.argtypes = [ElinaManagerPtr, c_size_t, c_size_t]
        a = elina_abstract0_top_c(man, intdim, realdim)
    except:
        print('Problem with loading/calling "elina_abstract0_top" from "libelinaux.so"')
        print('Make sure you are passing ElinaManagerPtr, c_size_t, c_size_t to the function')

    return a


def elina_abstract0_of_box(man, intdim, realdim, tinterval):
    """
    Create a hypercube ElinaAbstract0, defined by an ElinaIntervalArray.

    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    intdim : c_size_t
        Number of integer variables.
    realdim : c_size_t
        Number of real variables.
    tinterval : ElinaIntervalArray
        ElinaIntervalArray defining the hypercube.

    Returns
    -------
    a : ElinaAbstract0Ptr
        Pointer to the resulting ElinaAbstract0.

    """

    a = None
    try:
        elina_abstract0_of_box_c = elina_auxiliary_api.elina_abstract0_of_box
        elina_abstract0_of_box_c.restype = ElinaAbstract0Ptr
        elina_abstract0_of_box_c.argtypes = [ElinaManagerPtr, c_size_t, c_size_t, ElinaIntervalArray]
        a = elina_abstract0_of_box_c(man, intdim, realdim, tinterval)
    except:
        print('Problem with loading/calling "elina_abstract0_of_box" from "libelinaux.so"')
        print('Make sure you are passing ElinaManagerPtr, c_size_t, c_size_t, ElinaIntervalArray to the function')

    return a


# ============================================================ #
# II.2 Accessors
# ============================================================ #

def elina_abstract0_dimension(man, a):
    """
    Return the dimensionality of an ElinaAbstract0.

    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    a : ElinaAbstract0Ptr
        Pointer to the ElinaAbstract0 which dimensionality we want to obtain.

    Returns
    -------
    dim : ElinaDim
        Dimensionality of the ElinaAbstract0.

    """

    dim = None
    try:
        elina_abstract0_dimension_c = elina_auxiliary_api.elina_abstract0_dimension
        elina_abstract0_dimension_c.restype = ElinaDim
        elina_abstract0_dimension_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr]
        dim = elina_abstract0_dimension_c(man, a)
    except:
        print('Problem with loading/calling "elina_abstract0_dimension" from "libelinaux.so"')
        print('Make sure you are passing ElinaManagerPtr, ElinaAbstract0Ptr to the function')

    return dim


# ============================================================ #
# II.3 Tests
# ============================================================ #

def elina_abstract0_is_bottom(man, a):
    """
    Test if an ElinaAbstract0 is bottom.
    
    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    a : ElinaAbstract0Ptr
        Pointer to the ElinaAbstract0 that needs to be tested.

    Returns
    -------
    result : c_bool
        Result of the test.
        True means that the predicate is certainly true.
        False by default means we don't know.
        However, if the flag exact in the manager is true, then false means that the predicate is false.
        
    """

    result = None
    try:
        elina_abstract0_is_bottom_c = elina_auxiliary_api.elina_abstract0_is_bottom
        elina_abstract0_is_bottom_c.restype = c_bool
        elina_abstract0_is_bottom_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr]
        result = elina_abstract0_is_bottom_c(man, a)
    except:
        print('Problem with loading/calling "elina_abstract0_is_bottom" from "libelinaux.so"')
        print('Make sure you are passing ElinaManagerPtr, ElinaAbstract0Ptr to the function')

    return result


def elina_abstract0_is_top(man, a):
    """
    Test if an ElinaAbstract0 is top.

    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    a : ElinaAbstract0Ptr
        Pointer to the ElinaAbstract0 that needs to be tested.

    Returns
    -------
    result : c_bool
        Result of the test.
        True means that the predicate is certainly true.
        False by default means we don't know.
        However, if the flag exact in the manager is true, then false means that the predicate is false.

    """

    result = None
    try:
        elina_abstract0_is_top_c = elina_auxiliary_api.elina_abstract0_is_top
        elina_abstract0_is_top_c.restype = c_bool
        elina_abstract0_is_top_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr]
        result = elina_abstract0_is_top_c(man, a)
    except:
        print('Problem with loading/calling "elina_abstract0_is_top" from "libelinaux.so"')
        print('Make sure you are passing ElinaManagerPtr, ElinaAbstract0Ptr to the function')

    return result


def elina_abstract0_is_leq(man, a1, a2):
    """
    Test if one ElinaAbstract0 is less or equal than another.

    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    a1 : ElinaAbstract0Ptr
        Pointer to the first ElinaAbstract0 that needs to be tested.
    a2 : ElinaAbstract0Ptr
        Pointer to the second ElinaAbstract0 that needs to be tested.

    Returns
    -------
    result : c_bool
        Result of the test.
        True means that the predicate is certainly true.
        False by default means we don't know.
        However, if the flag exact in the manager is true, then false means that the predicate is false.

    """

    result = None
    try:
        elina_abstract0_is_leq_c = elina_auxiliary_api.elina_abstract0_is_leq
        elina_abstract0_is_leq_c.restype = c_bool
        elina_abstract0_is_leq_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr, ElinaAbstract0Ptr]
        result = elina_abstract0_is_leq_c(man, a1, a2)
    except:
        print('Problem with loading/calling "elina_abstract0_is_leq" from "libelinaux.so"')
        print('Make sure you are passing ElinaManagerPtr, ElinaAbstract0Ptr, ElinaAbstract0Ptr to the function')

    return result


def elina_abstract0_is_eq(man, a1, a2):
    """
    Test if one ElinaAbstract0 is equal to another.

    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    a1 : ElinaAbstract0Ptr
        Pointer to the first ElinaAbstract0 that needs to be tested.
    a2 : ElinaAbstract0Ptr
        Pointer to the second ElinaAbstract0 that needs to be tested.

    Returns
    -------
    result : c_bool
        Result of the test.
        True means that the predicate is certainly true.
        False by default means we don't know.
        However, if the flag exact in the manager is true, then false means that the predicate is false.

    """

    result = None
    try:
        elina_abstract0_is_eq_c = elina_auxiliary_api.elina_abstract0_is_eq
        elina_abstract0_is_eq_c.restype = c_bool
        elina_abstract0_is_eq_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr, ElinaAbstract0Ptr]
        result = elina_abstract0_is_eq_c(man, a1, a2)
    except:
        print('Problem with loading/calling "elina_abstract0_is_eq" from "libelinaux.so"')
        print('Make sure you are passing ElinaManagerPtr, ElinaAbstract0Ptr, ElinaAbstract0Ptr to the function')

    return result


def elina_abstract0_sat_lincons(man, a, lincons):
    """
    Test if an ElinaAbstract0 satisfies an ElinaLincons0.

    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    a : ElinaAbstract0Ptr
        Pointer to the first ElinaAbstract0 that needs to be tested.
    lincons : ElinaLincons0
        Pointer to the ElinaLincons0P that needs to be used for the test.
        
    Returns
    -------
    result : c_bool
        Result of the test.
        True means that the predicate is certainly true.
        False by default means we don't know.
        However, if the flag exact in the manager is true, then false means that the predicate is false.

    """

    result = None
    try:
        elina_abstract0_sat_lincons_c = elina_auxiliary_api.elina_abstract0_sat_lincons
        elina_abstract0_sat_lincons_c.restype = c_bool
        elina_abstract0_sat_lincons_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr, ElinaLincons0Ptr]
        result = elina_abstract0_sat_lincons_c(man, a, lincons)
    except:
        print('Problem with loading/calling "elina_abstract0_sat_lincons" from "libelinaux.so"')
        print('Make sure you are passing ElinaManagerPtr, ElinaAbstract0Ptr, ElinaLincons0Ptr to the function')

    return result


def elina_abstract0_sat_tcons(man, a, tcons):
    """
    Test if an ElinaAbstract0 satisfies an ElinaTcons0.

    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    a : ElinaAbstract0Ptr
        Pointer to the first ElinaAbstract0 that needs to be tested.
    tcons : ElinaTcons0Ptr
        Pointer to the ElinaTcons0 that needs to be used for the test.

    Returns
    -------
    result : c_bool
        Result of the test.
        True means that the predicate is certainly true.
        False by default means we don't know.
        However, if the flag exact in the manager is true, then false means that the predicate is false.

    """

    result = None
    try:
        elina_abstract0_sat_tcons_c = elina_auxiliary_api.elina_abstract0_sat_tcons
        elina_abstract0_sat_tcons_c.restype = c_bool
        elina_abstract0_sat_tcons_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr, ElinaTcons0Ptr]
        result = elina_abstract0_sat_tcons_c(man, a, tcons)
    except:
        print('Problem with loading/calling "elina_abstract0_sat_tcons" from "libelinaux.so"')
        print('Make sure you are passing ElinaManagerPtr, ElinaAbstract0Ptr, ElinaTcons0Ptr to the function')

    return result


def elina_abstract0_sat_interval(man, a, dim, interval):
    """
    Test if an ElinaDim is included in the interval in an ElinaAbstract0.

    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    a : ElinaAbstract0Ptr
        Pointer to the first ElinaAbstract0 that needs to be tested.
    dim : ElinaDim
        ElinaDim that needs to be used for the test.
    interval : ElinaIntervalPtr
        Pointer to the interval that needs to be used for the test.
    
    Returns
    -------
    result : c_bool
        Result of the test.
        True means that the predicate is certainly true.
        False by default means we don't know.
        However, if the flag exact in the manager is true, then false means that the predicate is false.

    """

    result = None
    try:
        elina_abstract0_sat_interval_c = elina_auxiliary_api.elina_abstract0_sat_interval
        elina_abstract0_sat_interval_c.restype = c_bool
        elina_abstract0_sat_interval_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr, ElinaDim, ElinaIntervalPtr]
        result = elina_abstract0_sat_interval_c(man, a, dim, interval)
    except:
        print('Problem with loading/calling "elina_abstract0_sat_interval" from "libelinaux.so"')
        print('Make sure you are passing ElinaManagerPtr, ElinaAbstract0Ptr, ElinaDim, '
              'ElinaIntervalPtr to the function')

    return result


def elina_abstract0_is_dimension_unconstrained(man, a, dim):
    """
    Test if an ElinaDim is unconstrained in an ElinaAbstract0.
    
    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    a : ElinaAbstract0Ptr
        Pointer to the first ElinaAbstract0 that needs to be tested.
    dim : ElinaDim
        ElinaDim that needs to be used for the test.

    Returns
    -------
    result : c_bool
        Result of the test.
        True means that the predicate is certainly true.
        False by default means we don't know.
        However, if the flag exact in the manager is true, then false means that the predicate is false.

    """

    result = None
    try:
        elina_abstract0_is_dimension_unconstrained_c = elina_auxiliary_api.elina_abstract0_is_dimension_unconstrained
        elina_abstract0_is_dimension_unconstrained_c.restype = c_bool
        elina_abstract0_is_dimension_unconstrained_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr, ElinaDim]
        result = elina_abstract0_is_dimension_unconstrained_c(man, a, dim)
    except:
        print('Problem with loading/calling "elina_abstract0_is_dimension_unconstrained" from "libelinaux.so"')
        print('Make sure you are passing ElinaManagerPtr, ElinaAbstract0Ptr, ElinaDim to the function')

    return result


# ============================================================ #
# II.4 Extraction of properties
# ============================================================ #

def elina_abstract0_bound_linexpr(man, a, linexpr):
    """
    Returns the ElinaInterval taken by an ElinaLinexpr0 over an ElinaAbstract0.
    
    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    a : ElinaAbstract0Ptr
        Pointer to the ElinaAbstract0.
    linexpr : ElinaLinexpr0Ptr
        Pointer to the ElinaLinexpr0.

    Returns
    -------
    interval : ElinaIntervalPtr
        Pointer to the ElinaInterval.

    """

    interval = None
    try:
        elina_abstract0_bound_linexpr_c = elina_auxiliary_api.elina_abstract0_bound_linexpr
        elina_abstract0_bound_linexpr_c.restype = ElinaIntervalPtr
        elina_abstract0_bound_linexpr_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr, ElinaLinexpr0Ptr]
        interval = elina_abstract0_bound_linexpr_c(man, a, linexpr)
    except:
        print('Problem with loading/calling "elina_abstract0_bound_linexpr" from "libelinaux.so"')
        print('Make sure you are passing ElinaManagerPtr, ElinaAbstract0Ptr, ElinaLinexpr0Ptr to the function')

    return interval


def elina_abstract0_bound_texpr(man, a, texpr):
    """
    Returns the ElinaInterval taken by an ElinaTexpr0 over an ElinaAbstract0.

    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    a : ElinaAbstract0Ptr
        Pointer to the ElinaAbstract0.
    texpr : ElinaTexpr0Ptr
        Pointer to the ElinaTexpr0.

    Returns
    -------
    interval : ElinaIntervalPtr
        Pointer to the ElinaInterval.

    """

    interval = None
    try:
        elina_abstract0_bound_texpr_c = elina_auxiliary_api.elina_abstract0_bound_texpr
        elina_abstract0_bound_texpr_c.restype = ElinaIntervalPtr
        elina_abstract0_bound_texpr_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr, ElinaTexpr0Ptr]
        interval = elina_abstract0_bound_texpr_c(man, a, texpr)
    except:
        print('Problem with loading/calling "elina_abstract0_bound_texpr" from "libelinaux.so"')
        print('Make sure you are passing ElinaManagerPtr, ElinaAbstract0Ptr, ElinaTexpr0Ptr to the function')

    return interval


def elina_abstract0_bound_dimension(man, a, dim):
    """
    Returns the ElinaInterval taken by an ElinaDim over an ElinaAbstract0.

    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    a : ElinaAbstract0Ptr
        Pointer to the ElinaAbstract0.
    dim : ElinaDim
        Pointer to the ElinaDim.

    Returns
    -------
    interval : ElinaIntervalPtr
        Pointer to the ElinaInterval.

    """

    interval = None
    try:
        elina_abstract0_bound_dimension_c = elina_auxiliary_api.elina_abstract0_bound_dimension
        elina_abstract0_bound_dimension_c.restype = ElinaIntervalPtr
        elina_abstract0_bound_dimension_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr, ElinaDim]
        interval = elina_abstract0_bound_dimension_c(man, a, dim)
    except:
        print('Problem with loading/calling "elina_abstract0_bound_dimension" from "libelinaux.so"')
        print('Make sure you are passing ElinaManagerPtr, ElinaAbstract0Ptr, ElinaDim to the function')

    return interval


def elina_abstract0_to_lincons_array(man, a):
    """
    Converts an ElinaAbstract0 to a polyhedra represented as a conjunction of multiple ElinaLincons0.
    The constraints are normally guaranteed to be scalar (without intervals).

    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    a : ElinaAbstract0Ptr
        Pointer to the ElinaAbstract0.

    Returns
    -------
    lincons_array : ElinaLincons0Array
        ElinaLincons0Array representing the polyhedra.
    
    """

    lincons_array = None
    try:
        elina_abstract0_to_lincons_array_c = elina_auxiliary_api.elina_abstract0_to_lincons_array
        elina_abstract0_to_lincons_array_c.restype = ElinaLincons0Array
        elina_abstract0_to_lincons_array_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr]
        lincons_array = elina_abstract0_to_lincons_array_c(man, a)
    except:
        print('Problem with loading/calling "elina_abstract0_to_lincons_array" from "libelinaux.so"')
        print('Make sure you are passing ElinaManagerPtr, ElinaAbstract0Ptr to the function')

    return lincons_array


def elina_abstract0_to_tcons_array(man, a):
    """
    Converts an ElinaAbstract0 to a conjunction of multiple ElinaTcons0.
    The constraints are normally guaranteed to be scalar (without intervals).

    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    a : ElinaAbstract0Ptr
        Pointer to the ElinaAbstract0.

    Returns
    -------
    tcons_array : ElinaTcons0Array
        ElinaTcons0Array representing the tree expressions constraints.

    """

    tcons_array = None
    try:
        elina_abstract0_to_tcons_array_c = elina_auxiliary_api.elina_abstract0_to_tcons_array
        elina_abstract0_to_tcons_array_c.restype = ElinaTcons0Array
        elina_abstract0_to_tcons_array_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr]
        tcons_array = elina_abstract0_to_tcons_array_c(man, a)
    except:
        print('Problem with loading/calling "elina_abstract0_to_tcons_array" from "libelinaux.so"')
        print('Make sure you are passing ElinaManagerPtr, ElinaAbstract0Ptr to the function')

    return tcons_array


def elina_abstract0_to_box(man, a):
    """
    Converts an ElinaAbstract0 to an ElinaIntervalArray hypercube.
    The size of the resulting ElinaIntervalArray is elina_abstract0_dimension(man,a).
    This function can be reimplemented by using elina_abstract0_bound_linexpr.
    
    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    a : ElinaAbstract0Ptr
        Pointer to the ElinaAbstract0.

    Returns
    -------
    interval_array : ElinaIntervalArray
        ElinaIntervalArray representing the hypercube.

    """

    interval_array = None
    try:
        elina_abstract0_to_box_c = elina_auxiliary_api.elina_abstract0_to_box
        elina_abstract0_to_box_c.restype = ElinaIntervalArray
        elina_abstract0_to_box_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr]
        interval_array = elina_abstract0_to_box_c(man, a)
    except:
        print('Problem with loading/calling "elina_abstract0_to_box" from "libelinaux.so"')
        print('Make sure you are passing ElinaManagerPtr, ElinaAbstract0Ptr to the function')

    return interval_array


# ********************************************************************** #
# III. Operations
# ********************************************************************** #

# ============================================================ #
# III.1 Meet and Join
# ============================================================ #

def elina_abstract0_meet(man, destructive, a1, a2):
    """
    Meet of two ElinaAbstract0-s.

    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    destructive : c_bool
        Boolean flag.
    a1 : ElinaAbstract0Ptr
        Pointer to the first ElinaAbstract0.
    a2 : ElinaAbstract0Ptr
        Pointer to the second ElinaAbstract0.

    Returns
    -------
    a0 : ElinaAbstract0Ptr
        Pointer to the resulting ElinaAbstract0.

    """

    a0 = None
    try:
        elina_abstract0_meet_c = elina_auxiliary_api.elina_abstract0_meet
        elina_abstract0_meet_c.restype = ElinaAbstract0Ptr
        elina_abstract0_meet_c.argtypes = [ElinaManagerPtr, c_bool, ElinaAbstract0Ptr, ElinaAbstract0Ptr]
        a0 = elina_abstract0_meet_c(man, destructive, a1, a2)
    except:
        print('Problem with loading/calling "elina_abstract0_meet" from "libelinaux.so"')
        print('Make sure you are passing ElinaManagerPtr, c_bool, ElinaAbstract0Ptr, ElinaAbstract0Ptr to the function')

    return a0


def elina_abstract0_join(man, destructive, a1, a2):
    """
    Join of two ElinaAbstract0-s.

    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    destructive : c_bool
        Boolean flag.
    a1 : ElinaAbstract0Ptr
        Pointer to the first ElinaAbstract0.
    a2 : ElinaAbstract0Ptr
        Pointer to the second ElinaAbstract0.

    Returns
    -------
    a0 : ElinaAbstract0Ptr
        Pointer to the resulting ElinaAbstract0.

    """

    a0 = None
    try:
        elina_abstract0_join_c = elina_auxiliary_api.elina_abstract0_join
        elina_abstract0_join_c.restype = ElinaAbstract0Ptr
        elina_abstract0_join_c.argtypes = [ElinaManagerPtr, c_bool, ElinaAbstract0Ptr, ElinaAbstract0Ptr]
        a0 = elina_abstract0_join_c(man, destructive, a1, a2)
    except:
        print('Problem with loading/calling "elina_abstract0_join" from "libelinaux.so"')
        print('Make sure you are passing ElinaManagerPtr, c_bool, ElinaAbstract0Ptr, ElinaAbstract0Ptr to the function')

    return a0


def elina_abstract0_meet_array(man, tab, size):
    """
    Meet of an array of ElinaAbstract0-s.
    Raises an [[exc_invalid_argument]] exception if [[size==0]].

    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    tab : ElinaAbstract0Array
        Array of ElinaAbstract0Ptr-s for which we obtain the meet.
    size : c_size_t
        
    
    Returns
    -------
    a0 : ElinaAbstract0Ptr
        Pointer to the resulting ElinaAbstract0.

    """

    a = None
    try:
        elina_abstract0_meet_array_c = elina_auxiliary_api.elina_abstract0_meet_array
        elina_abstract0_meet_array_c.restype = ElinaAbstract0Ptr
        elina_abstract0_meet_array_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Array, c_size_t]
        a = elina_abstract0_meet_array_c(man, tab, size)
    except:
        print('Problem with loading/calling "elina_abstract0_meet_array" from "libelinaux.so"')
        print('Make sure you are passing ElinaManagerPtr, ElinaAbstract0Array, c_size_t to the function')

    return a


def elina_abstract0_join_array(man, tab, size):
    """
    Join of an array of ElinaAbstract0-s.
    Raises an [[exc_invalid_argument]] exception if [[size==0]].

    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    tab : ElinaAbstract0Array
        Array of ElinaAbstract0Ptr-s for which we obtain the meet.
    size : c_size_t
        Size of the ElinaAbstract0Array.

    Returns
    -------
    a0 : ElinaAbstract0Ptr
        Pointer to the resulting ElinaAbstract0.

    """

    a = None
    try:
        elina_abstract0_join_array_c = elina_auxiliary_api.elina_abstract0_join_array
        elina_abstract0_join_array_c.restype = ElinaAbstract0Ptr
        elina_abstract0_join_array_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Array, c_size_t]
        a = elina_abstract0_join_array_c(man, tab, size)
    except:
        print('Problem with loading/calling "elina_abstract0_join_array" from "libelinaux.so"')
        print('Make sure you are passing ElinaManagerPtr, ElinaAbstract0Array, c_size_t to the function')

    return a


def elina_abstract0_meet_lincons_array(man, destructive, a1, lincons_array):
    """
    Meet of an ElinaAbstract0 with an ElinaLincons0Array.

    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    destructive : c_bool
        Boolean flag.
    a1 : ElinaAbstract0
        Pointer to the ElinaAbstract0.
    lincons_array : ElinaLincons0ArrayPtr
        Pointer to the ElinaLincons0Array.
        
    Returns
    -------
    a0 : ElinaAbstract0Ptr
        Pointer to the resulting ElinaAbstract0.

    """

    a0 = None
    try:
        elina_abstract0_meet_lincons_array_c = elina_auxiliary_api.elina_abstract0_meet_lincons_array
        elina_abstract0_meet_lincons_array_c.restype = ElinaAbstract0Ptr
        elina_abstract0_meet_lincons_array_c.argtypes = [ElinaManagerPtr, c_bool,
                                                         ElinaAbstract0Ptr, ElinaLincons0ArrayPtr]
        a0 = elina_abstract0_meet_lincons_array_c(man, destructive, a1, lincons_array)
    except:
        print('Problem with loading/calling "elina_abstract0_meet_lincons_array" from "libelinaux.so"')
        print('Make sure you are passing ElinaManagerPtr, c_bool, ElinaAbstract0Ptr, '
              'ElinaLincons0ArrayPtr to the function')

    return a0


def elina_abstract0_meet_tcons_array(man, destructive, a1, tcons_array):
    """
    Meet of an ElinaAbstract0 with an ElinaTcons0Array.

    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    destructive : c_bool
        Boolean flag.
    a1 : ElinaAbstract0
        Pointer to the ElinaAbstract0.
    tcons_array : ElinaTcons0ArrayPtr
        Pointer to the ElinaTcons0Array.
        
    Returns
    -------
    a0 : ElinaAbstract0Ptr
        Pointer to the resulting ElinaAbstract0.

    """

    a0 = None
    try:
        elina_abstract0_meet_tcons_array_c = elina_auxiliary_api.elina_abstract0_meet_tcons_array
        elina_abstract0_meet_tcons_array_c.restype = ElinaAbstract0Ptr
        elina_abstract0_meet_tcons_array_c.argtypes = [ElinaManagerPtr, c_bool,
                                                       ElinaAbstract0Ptr, ElinaTcons0ArrayPtr]
        a0 = elina_abstract0_meet_tcons_array_c(man, destructive, a1, tcons_array)
    except:
        print('Problem with loading/calling "elina_abstract0_meet_tcons_array" from "libelinaux.so"')
        print('Make sure you are passing ElinaManagerPtr, c_bool, ElinaAbstract0Ptr, '
              'ElinaTcons0ArrayPtr to the function')

    return a0


# ============================================================ #
# III.2 Assignment and Substitutions
# ============================================================ #

def elina_abstract0_assign_linexpr_array(man, destructive, org, tdim, linexpr_array, size, dest):
    """
    Parallel assignment of several dimensions of an ElinaAbstract0 by using an ElinaLinexpr0Array.
    
    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    destructive : c_bool
        Boolean flag.
    org : ElinaAbstract0Ptr
        Pointer to the ElinaAbstract0 which dimensions need to be assigned.
    tdim : ElinaDimPtr
        Array of ElinaDim that need to be assigned.
    linexpr_array : ElinaLinexpr0Array
        Array of ElinaLinexpr0 used for assigning the dimensions.
    size : c_size_t
        Size of the ElinaLinexpr0Array.
    dest : ElinaAbstract0Ptr
        If not NULL the result of the transformation is intersected with this ElinaAbstract0. 
        This is useful for precise backward transformations in lattices like intervals or octagons.

    Returns
    -------
    a : ElinaAbstract0Ptr
        Pointer to the resulting ElinaAbstract0.

    """

    a = None
    try:
        elina_abstract0_assign_linexpr_array_c = elina_auxiliary_api.elina_abstract0_assign_linexpr_array
        elina_abstract0_assign_linexpr_array_c.restype = ElinaAbstract0Ptr
        elina_abstract0_assign_linexpr_array_c.argtypes = [ElinaManagerPtr, c_bool, ElinaAbstract0Ptr,
                                                           ElinaDimPtr, ElinaLinexpr0Array, c_size_t,
                                                           ElinaAbstract0Ptr]
        a = elina_abstract0_assign_linexpr_array_c(man, destructive, org, tdim, linexpr_array, size, dest)
    except:
        print('Problem with loading/calling "elina_abstract0_assign_linexpr_array" from "libelinaux.so"')
        print('Make sure you are passing ElinaManagerPtr, c_bool, ElinaAbstract0Ptr, ElinaDimPtr, ElinaLinexpr0Array, '
              'c_size_t, ElinaAbstract0Ptr to the function')

    return a


def elina_abstract0_assign_texpr_array(man, destructive, org, tdim, texpr_array, size, dest):
    """
    Parallel assignment of several dimensions of an ElinaAbstract0 by using an ElinaTexpr0Array.

    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    destructive : c_bool
        Boolean flag.
    org : ElinaAbstract0Ptr
        Pointer to the ElinaAbstract0 which dimensions need to be assigned.
    tdim : ElinaDimPtr
        Array of ElinaDim that need to be assigned.
    texpr_array : ElinaTexpr0Array
        Array of ElinaTexpr0 used for assigning the dimensions.
    size : c_size_t
        Size of the ElinaTexpr0Array.
    dest : ElinaAbstract0Ptr
        If not NULL the result of the transformation is intersected with this ElinaAbstract0. 
        This is useful for precise backward transformations in lattices like intervals or octagons.

    Returns
    -------
    a : ElinaAbstract0Ptr
        Pointer to the resulting ElinaAbstract0.

    """

    a = None
    try:
        elina_abstract0_assign_texpr_array_c = elina_auxiliary_api.elina_abstract0_assign_texpr_array
        elina_abstract0_assign_texpr_array_c.restype = ElinaAbstract0Ptr
        elina_abstract0_assign_texpr_array_c.argtypes = [ElinaManagerPtr, c_bool, ElinaAbstract0Ptr, ElinaDimPtr,
                                                         ElinaTexpr0Array, c_size_t, ElinaAbstract0Ptr]
        a = elina_abstract0_assign_texpr_array_c(man, destructive, org, tdim, texpr_array, size, dest)
    except:
        print('Problem with loading/calling "elina_abstract0_assign_texpr_array" from "libelinaux.so"')
        print('Make sure you are passing ElinaManagerPtr, c_bool, ElinaAbstract0Ptr, ElinaDimPtr, ElinaTexpr0Array, '
              'c_size_t, ElinaAbstract0Ptr to the function')

    return a


def elina_abstract0_substitute_linexpr_array(man, destructive, org, tdim, linexpr_array, size, dest):
    """
    Parallel substitution of several dimensions of an ElinaAbstract0 by using an ElinaLinexpr0Array.

    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    destructive : c_bool
        Boolean flag.
    org : ElinaAbstract0Ptr
        Pointer to the ElinaAbstract0 which dimensions need to be substituted.
    tdim : ElinaDimPtr
        Array of ElinaDim that need to be substituted.
    linexpr_array : ElinaLinexpr0Array
        Array of ElinaLinexpr0 used for substituting the dimensions.
    size : c_size_t
        Size of the ElinaLinexpr0Array.
    dest : ElinaAbstract0Ptr
        If not NULL the result of the transformation is intersected with this ElinaAbstract0. 
        This is useful for precise backward transformations in lattices like intervals or octagons.

    Returns
    -------
    a : ElinaAbstract0Ptr
        Pointer to the resulting ElinaAbstract0.
        
    """

    a = None
    try:
        elina_abstract0_substitute_linexpr_array_c = elina_auxiliary_api.elina_abstract0_substitute_linexpr_array
        elina_abstract0_substitute_linexpr_array_c.restype = ElinaAbstract0Ptr
        elina_abstract0_substitute_linexpr_array_c.argtypes = [ElinaManagerPtr, c_bool, ElinaAbstract0Ptr, ElinaDimPtr,
                                                               ElinaLinexpr0Array, c_size_t, ElinaAbstract0Ptr]
        a = elina_abstract0_substitute_linexpr_array_c(man, destructive, org, tdim, linexpr_array, size, dest)
    except:
        print('Problem with loading/calling "elina_abstract0_substitute_linexpr_array" from "libelinaux.so"')
        print('Make sure you are passing ElinaManagerPtr, c_bool, ElinaAbstract0Ptr, ElinaDimPtr, ElinaLinexpr0Array, '
              'c_size_t, ElinaAbstract0Ptr to the function')

    return a


def elina_abstract0_substitute_texpr_array(man, destructive, org, tdim, texpr_array, size, dest):
    """
    Parallel substitution of several dimensions of an ElinaAbstract0 by using an ElinaTexpr0Array.

    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    destructive : c_bool
        Boolean flag.
    org : ElinaAbstract0Ptr
        Pointer to the ElinaAbstract0 which dimensions need to be substituted.
    tdim : ElinaDimPtr
        Array of ElinaDim that need to be substituted.
    texpr_array : ElinaTexpr0Array
        Array of ElinaTexpr0A used for substituting the dimensions.
    size : c_size_t
        Size of the ElinaTexpr0Array.
    dest : ElinaAbstract0Ptr
        If not NULL the result of the transformation is intersected with this ElinaAbstract0. 
        This is useful for precise backward transformations in lattices like intervals or octagons.

    Returns
    -------
    a : ElinaAbstract0Ptr
        Pointer to the resulting ElinaAbstract0.

    """

    a = None
    try:
        elina_abstract0_substitute_texpr_array_c = elina_auxiliary_api.elina_abstract0_substitute_texpr_array
        elina_abstract0_substitute_texpr_array_c.restype = ElinaAbstract0Ptr
        elina_abstract0_substitute_texpr_array_c.argtypes = [ElinaManagerPtr, c_bool, ElinaAbstract0Ptr, ElinaDimPtr,
                                                             ElinaTexpr0Array, c_size_t, ElinaAbstract0Ptr]
        a = elina_abstract0_substitute_texpr_array_c(man, destructive, org, tdim, texpr_array, size, dest)
    except:
        print('Problem with loading/calling "elina_abstract0_substitute_texpr_array" from "libelinaux.so"')
        print('Make sure you are passing ElinaManagerPtr, c_bool, ElinaAbstract0Ptr, ElinaDimPtr, ElinaTexpr0Array, '
              'c_size_t, ElinaAbstract0Ptr to the function')

    return a


# ============================================================ #
# III.3 Projections
# ============================================================ #

def elina_abstract0_forget_array(man, destructive, a1, tdim, size):
    # TODO finish documentation for this function
    """
    !?What does this function do!?
    
    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    destructive : c_bool
        Boolean flag.
    a1 : ElinaAbstract0Ptr
    tdim : ElinaDimPtr
    size : c_size-t

    Returns
    -------
    a0 : ElinaAbstract0Ptr
        Pointer to the resulting ElinaAbstract0.

    """

    a0 = None
    try:
        elina_abstract0_forget_array_c = elina_auxiliary_api.elina_abstract0_forget_array
        elina_abstract0_forget_array_c.restype = ElinaAbstract0Ptr
        elina_abstract0_forget_array_c.argtypes = [ElinaManagerPtr, c_bool, ElinaAbstract0Ptr, ElinaDimPtr, c_size_t]
        a0 = elina_abstract0_forget_array_c(man, destructive, a1, tdim, size)
    except:
        print('Problem with loading/calling "elina_abstract0_forget_array" from "libelinaux.so"')
        print('Make sure you are passing ElinaManagerPtr, c_bool, ElinaAbstract0Ptr, ElinaDimPtr, c_size_t '
              'to the function')

    return a0


# ============================================================ #
# III.4 Change and permutation of dimensions
# ============================================================ #

def elina_abstract0_add_dimensions(man, destructive, a1, dimchange, project):
    """
    Add dimensions to an ElinaAbstract0 by following the semantics of an ElinaDimchange.

    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    destructive : c_bool
        Boolean flag.
    a1 : ElinaAbstract0Ptr
        Pointer to the ElinaAbstract0 to which we want to add dimensions.
    dimchange : ElinaDimchangePtr
        Pointer to the ElinaDimchange which semantics we want to follow.
    project : c_bool
        Projection boolean flag.

    Returns
    -------
    a0 : ElinaAbstract0Ptr
        Pointer to the resulting ElinaAbstract0.

    """

    a0 = None
    try:
        elina_abstract0_add_dimensions_c = elina_auxiliary_api.elina_abstract0_add_dimensions
        elina_abstract0_add_dimensions_c.restype = ElinaAbstract0Ptr
        elina_abstract0_add_dimensions_c.argtypes = [ElinaManagerPtr, c_bool, ElinaAbstract0Ptr,
                                                     ElinaDimchangePtr, c_bool]
        a0 = elina_abstract0_add_dimensions_c(man, destructive, a1, dimchange, project)
    except:
        print('Problem with loading/calling "elina_abstract0_add_dimensions" from "libelinaux.so"')
        print('Make sure you are passing ElinaManagerPtr, c_bool, ElinaAbstract0Ptr,ElinaManagerPtr, '
              'c_bool, ElinaAbstract0Ptr, ElinaDimchangePtr, c_bool to the function')

    return a0


def elina_abstract0_remove_dimensions(man, destructive, a1, dimchange):
    """
    Remove dimensions from an ElinaAbstract0 by following the semantics of an ElinaDimchange.

    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    destructive : c_bool
        Boolean flag.
    a1 : ElinaAbstract0Ptr
        Pointer to the ElinaAbstract0 from which we want to remove dimensions.
    dimchange : ElinaDimchangePtr
        Pointer to the ElinaDimchange which semantics we want to follow.

    Returns
    -------
    a0 : ElinaAbstract0Ptr
        Pointer to the resulting ElinaAbstract0.

    """

    a0 = None
    try:
        elina_abstract0_remove_dimensions_c = elina_auxiliary_api.elina_abstract0_remove_dimensions
        elina_abstract0_remove_dimensions_c.restype = ElinaAbstract0Ptr
        elina_abstract0_remove_dimensions_c.argtypes = [ElinaManagerPtr, c_bool, ElinaAbstract0Ptr, ElinaDimchangePtr]
        a0 = elina_abstract0_remove_dimensions_c(man, destructive, a1, dimchange)
    except:
        print('Problem with loading/calling "elina_abstract0_remove_dimensions" from "libelinaux.so"')
        print('Make sure you are passing ElinaManagerPtr, c_bool, ElinaAbstract0Ptr,ElinaManagerPtr, '
              'c_bool, ElinaAbstract0Ptr, ElinaDimchangePtr to the function')

    return a0


def elina_abstract0_permute_dimensions(man, destructive, a1, dimperm):
    """
    Permute dimensions of an ElinaAbstract0 by following the semantics of an ElinaDimperm.

    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    destructive : c_bool
        Boolean flag.
    a1 : ElinaAbstract0Ptr
        Pointer to the ElinaAbstract0 which dimensions we want to permute.
    dimperm : ElinaDimpermPtr
        Pointer to the ElinaDimperm which semantics we want to follow.

    Returns
    -------
    a0 : ElinaAbstract0Ptr
        Pointer to the resulting ElinaAbstract0.

    """

    a0 = None
    try:
        elina_abstract0_permute_dimensions_c = elina_auxiliary_api.elina_abstract0_permute_dimensions
        elina_abstract0_permute_dimensions_c.restype = ElinaAbstract0Ptr
        elina_abstract0_permute_dimensions_c.argtypes = [ElinaManagerPtr, c_bool, ElinaAbstract0Ptr, ElinaDimpermPtr]
        a0 = elina_abstract0_permute_dimensions_c(man, destructive, a1, dimperm)
    except:
        print('Problem with loading/calling "elina_abstract0_permute_dimensions" from "libelinaux.so"')
        print('Make sure you are passing ElinaManagerPtr, c_bool, ElinaAbstract0Ptr,ElinaManagerPtr, '
              'c_bool, ElinaAbstract0Ptr, ElinaDimpermPtr to the function')

    return a0


# ============================================================ #
# III.5 Expansion and folding of dimensions
# ============================================================ #

def elina_abstract0_expand(man, destructive, a1, dim, n):
    """
    Expand an ElinaDim of an ElinaAbstract0 into itself + n additional dimensions. 
    Resulting in (n+1) unrelated dimensions having same relations with other dimensions.
    The (n+1) dimensions are put as follows:
    - original dimension dim
    - if the dimension is integer, the n additional dimensions are put at the end of the integer dimensions; 
    - if it is real, at the end of the real dimensions.
    
    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    destructive : c_bool
        Boolean flag.
    a1 : ElinaAbstract0Ptr
        Pointer to the ElinaAbstract0 which dimension we want to expand.
    dim : ElinaDim
        ElinaDim that needs to be expanded.
    n : c_size_t
        Amount of additional dimensions.
    
    Returns
    -------
    a0 : ElinaAbstract0Ptr
        Pointer to the resulting ElinaAbstract0.

    """

    a0 = None
    try:
        elina_abstract0_expand_c = elina_auxiliary_api.elina_abstract0_expand
        elina_abstract0_expand_c.restype = ElinaAbstract0Ptr
        elina_abstract0_expand_c.argtypes = [ElinaManagerPtr, c_bool, ElinaAbstract0Ptr, ElinaDim, c_size_t]
        a0 = elina_abstract0_expand_c(man, destructive, a1, dim, n)
    except:
        print('Problem with loading/calling "elina_abstract0_expand" from "libelinaux.so"')
        print('Make sure you are passing ElinaManagerPtr, c_bool, ElinaAbstract0Ptr, ElinaDim, '
              'c_size_t to the function')

    return a0


def elina_abstract0_fold(man, destructive, a1, tdim, n):
    """
    Fold the dimensions of an ElinaAbstract0, specified by an array of ElinaDim-s of size >= 1.
    Put the result in the first dimension and remove other dimensions.

    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    destructive : c_bool
        Boolean flag.
    a1 : ElinaAbstract0Ptr
        Pointer to the ElinaAbstract0 which dimension we want to fold.
    tdim : ElinaDimPtr
        Array of ElinaDim-s that needs to be folded.
    n : c_size_t
        Size of the array of ElinaDim-s.

    Returns
    -------
    a0 : ElinaAbstract0Ptr
        Pointer to the resulting ElinaAbstract0.

    """

    a0 = None
    try:
        elina_abstract0_fold_c = elina_auxiliary_api.elina_abstract0_fold
        elina_abstract0_fold_c.restype = ElinaAbstract0Ptr
        elina_abstract0_fold_c.argtypes = [ElinaManagerPtr, c_bool, ElinaAbstract0Ptr, ElinaDimPtr, c_size_t]
        a0 = elina_abstract0_fold_c(man, destructive, a1, tdim, n)
    except:
        print('Problem with loading/calling "elina_abstract0_fold" from "libelinaux.so"')
        print('Make sure you are passing ElinaManagerPtr, c_bool, ElinaAbstract0Ptr, ElinaDimPtr, '
              'c_size_t to the function')

    return a0


# ============================================================ #
# III.6 Widening
# ============================================================ #

def elina_abstract0_widening(man, a1, a2):
    # TODO finish documentation for this function
    """
    !?What does this function do!?
    
    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    a1 : ElinaAbstract0Ptr
        Pointer to the first ElinaAbstract0.
    a2 : ElinaAbstract0Ptr
        Pointer to the second ElinaAbstract0.

    Returns
    -------
    a0 : ElinaAbstract0Ptr
        Pointer to the resulting ElinaAbstract0.

    """

    a0 = None
    try:
        elina_abstract0_widening_c = elina_auxiliary_api.elina_abstract0_widening
        elina_abstract0_widening_c.restype = ElinaAbstract0Ptr
        elina_abstract0_widening_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr, ElinaAbstract0Ptr]
        a0 = elina_abstract0_widening_c(man, a1, a2)
    except:
        print('Problem with loading/calling "elina_abstract0_widening" from "libelinaux.so"')
        print('Make sure you are passing ElinaManagerPtr, ElinaAbstract0Ptr, ElinaAbstract0Ptr to the function')

    return a0


# ============================================================ #
# III.7 Closure operation
# ============================================================ #

def elina_abstract0_closure(man, destructive, a1):
    """
    Return the topological closure of a possibly opened ElinaAbstract0.
    
    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    destructive : c_bool
        Boolean flag.
    a1 : ElinaAbstract0Ptr
        Pointer to the ElinaAbstract0.
        
    Returns
    -------
    a0 : ElinaAbstract0Ptr
        Pointer to the resulting ElinaAbstract0.

    """

    a0 = None
    try:
        elina_abstract0_closure_c = elina_auxiliary_api.elina_abstract0_closure
        elina_abstract0_closure_c.restype = ElinaAbstract0Ptr
        elina_abstract0_closure_c.argtypes = [ElinaManagerPtr, c_bool, ElinaAbstract0Ptr]
        a0 = elina_abstract0_closure_c(man, destructive, a1)
    except:
        print('Problem with loading/calling "elina_abstract0_closure" from "libelinaux.so"')
        print('Make sure you are passing ElinaManagerPtr, c_bool, ElinaAbstract0Ptr to the function')

    return a0


# ********************************************************************** #
# IV. Functions offered by the ELINA interface
# ********************************************************************** #

def elina_abstract0_manager(a):
    """
    Return a reference to the manager contained in an ElinaAbstract0.
    The reference should not be freed.
    
    Parameters
    ----------
    a : ElinaAbstract0Ptr
        Pointer to the ElinaAbstract0.

    Returns
    -------
    man : ElinaManagerPtr
        Pointer to the resulting ElinaManager.

    """

    man = None
    try:
        elina_abstract0_manager_c = elina_auxiliary_api.elina_abstract0_manager
        elina_abstract0_manager_c.restype = ElinaManagerPtr
        elina_abstract0_manager_c.argtypes = [ElinaAbstract0Ptr]
        man = elina_abstract0_manager_c(a)
    except:
        print('Problem with loading/calling "elina_abstract0_manager" from "libelinaux.so"')
        print('Make sure you are passing ElinaAbstract0Ptr to the function')

    return man


def elina_abstract0_of_lincons_array(man, intdim, realdim, lincons_array):
    """
    Transform an ElinaLincons0Array to an ElinaAbstract0.

    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    intdim : c_size_t
        Number of integer dimensions.
    realdim : c_size_t
        Number of real dimensions.
    lincons_array : ElinaLincons0ArrayPtr
        Pointer to the ElinaLincons0Array that needs to be transformed.

    Returns
    -------
    a : ElinaAbstract0Ptr
        Pointer to the resulting ElinaAbstract0.

    """

    a = None
    try:
        elina_abstract0_of_lincons_array_c = elina_auxiliary_api.elina_abstract0_of_lincons_array
        elina_abstract0_of_lincons_array_c.restype = ElinaAbstract0Ptr
        elina_abstract0_of_lincons_array_c.argtypes = [ElinaManagerPtr, c_size_t, c_size_t, ElinaLincons0ArrayPtr]
        a = elina_abstract0_of_lincons_array_c(man, intdim, realdim, lincons_array)
    except:
        print('Problem with loading/calling "elina_abstract0_of_lincons_array" from "libelinaux.so"')
        print('Make sure you are passing ElinaManagerPtr, c_size_t, c_size_t, ElinaLincons0ArrayPtr to the function')

    return a


def elina_abstract0_of_tcons_array(man, intdim, realdim, tcons_array):
    """
    Transform an ElinaTcons0Array to an ElinaAbstract0.

    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    intdim : c_size_t
        Number of integer dimensions.
    realdim : c_size_t
        Number of real dimensions.
    tcons_array : ElinaTcons0ArrayPtr
        Pointer to the ElinaTcons0Array that needs to be transformed.

    Returns
    -------
    a : ElinaAbstract0Ptr
        Pointer to the resulting ElinaAbstract0.

    """

    a = None
    try:
        elina_abstract0_of_lincons_array_c = elina_auxiliary_api.elina_abstract0_of_lincons_array
        elina_abstract0_of_lincons_array_c.restype = ElinaAbstract0Ptr
        elina_abstract0_of_lincons_array_c.argtypes = [ElinaManagerPtr, c_size_t, c_size_t, ElinaTcons0ArrayPtr]
        a = elina_abstract0_of_lincons_array_c(man, intdim, realdim, tcons_array)
    except:
        print('Problem with loading/calling "elina_abstract0_of_lincons_array" from "libelinaux.so"')
        print('Make sure you are passing ElinaManagerPtr, c_size_t, c_size_t, ElinaTcons0ArrayPtr to the function')

    return a


def elina_abstract0_assign_linexpr(man, destructive, org, dim, linexpr, dest):
    """
    Assign a single ElinaDim of an ElinaAbstract0 by using an ElinaLinexpr0.

    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    destructive : c_bool
        Boolean flag.
    org : ElinaAbstract0Ptr
        Pointer to the ElinaAbstract0 which dimension we need to assign.
    dim : ElinaDim
        ElinaDim we want to assign.
    linexpr : ElinaLinexpr0Ptr
        Pointer to the ElinaLinexpr0 tht needs to be used for assignment.
    dest : ElinaAbstract0Ptr
        Pointer to an ElinaAbstract0 used as an optional argument.
        If not NULL, the result of the transformation is intersected with dest.
        This is useful for precise backward transformations in lattices like intervals or octagons.
    
    Returns
    -------
    a : ElinaAbstract0Ptr
        Pointer to the resulting ElinaAbstract0.

    """

    a = None
    try:
        elina_abstract0_assign_linexpr_c = elina_auxiliary_api.elina_abstract0_assign_linexpr
        elina_abstract0_assign_linexpr_c.restype = ElinaAbstract0Ptr
        elina_abstract0_assign_linexpr_c.argtypes = [ElinaManagerPtr, c_bool, ElinaAbstract0Ptr,
                                                     ElinaDim, ElinaLinexpr0Ptr, ElinaAbstract0Ptr]
        a = elina_abstract0_assign_linexpr_c(man, destructive, org, dim, linexpr, dest)
    except:
        print('Problem with loading/calling "elina_abstract0_assign_linexpr" from "libelinaux.so"')
        print('Make sure you are passing ElinaManagerPtr, c_bool, ElinaAbstract0Ptr, ElinaDim, '
              'ElinaLinexpr0Ptr, ElinaAbstract0Ptr to the function')

    return a


def elina_abstract0_assign_texpr(man, destructive, org, dim, texpr, dest):
    """
    Assign a single ElinaDim of an ElinaAbstract0 by using an ElinaTexpr0P.

    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    destructive : c_bool
        Boolean flag.
    org : ElinaAbstract0Ptr
        Pointer to the ElinaAbstract0 which dimension we need to assign.
    dim : ElinaDim
        ElinaDim we want to assign.
    texpr : ElinaTexpr0Ptr
        Pointer to the ElinaTexpr0 tht needs to be used for assignment.
    dest : ElinaAbstract0Ptr
        Pointer to an ElinaAbstract0 used as an optional argument.
        If not NULL, the result of the transformation is intersected with dest.
        This is useful for precise backward transformations in lattices like intervals or octagons.

    Returns
    -------
    a : ElinaAbstract0Ptr
        Pointer to the resulting ElinaAbstract0.

    """

    a = None
    try:
        elina_abstract0_assign_texpr_c = elina_auxiliary_api.elina_abstract0_assign_texpr
        elina_abstract0_assign_texpr_c.restype = ElinaAbstract0Ptr
        elina_abstract0_assign_texpr_c.argtypes = [ElinaManagerPtr, c_bool, ElinaAbstract0Ptr,
                                                   ElinaDim, ElinaTexpr0Ptr, ElinaAbstract0Ptr]
        a = elina_abstract0_assign_texpr_c(man, destructive, org, dim, texpr, dest)
    except:
        print('Problem with loading/calling "elina_abstract0_assign_texpr" from "libelinaux.so"')
        print('Make sure you are passing ElinaManagerPtr, c_bool, ElinaAbstract0Ptr, ElinaDim, '
              'ElinaTexpr0Ptr, ElinaAbstract0Ptr to the function')

    return a


def elina_abstract0_substitute_linexpr(man, destructive, org, dim, linexpr, dest):
    """
    Substitute a single ElinaDim of an ElinaAbstract0 by using an ElinaLinexpr0.

    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    destructive : c_bool
        Boolean flag.
    org : ElinaAbstract0Ptr
        Pointer to the ElinaAbstract0 which dimension we need to substitute.
    dim : ElinaDim
        ElinaDim we need to substitute.
    linexpr : ElinaLinexpr0Ptr
        Pointer to the ElinaLinexpr0 tht needs to be used for substitution.
    dest : ElinaAbstract0Ptr
        Pointer to an ElinaAbstract0 used as an optional argument.
        If not NULL, the result of the transformation is intersected with dest.
        This is useful for precise backward transformations in lattices like intervals or octagons.

    Returns
    -------
    a : ElinaAbstract0Ptr
        Pointer to the resulting ElinaAbstract0.

    """

    a = None
    try:
        elina_abstract0_substitute_linexpr_c = elina_auxiliary_api.elina_abstract0_substitute_linexpr
        elina_abstract0_substitute_linexpr_c.restype = ElinaAbstract0Ptr
        elina_abstract0_substitute_linexpr_c.argtypes = [ElinaManagerPtr, c_bool, ElinaAbstract0Ptr,
                                                         ElinaDim, ElinaLinexpr0Ptr, ElinaAbstract0Ptr]
        a = elina_abstract0_substitute_linexpr_c(man, destructive, org, dim, linexpr, dest)
    except:
        print('Problem with loading/calling "elina_abstract0_substitute_linexpr" from "libelinaux.so"')
        print('Make sure you are passing ElinaManagerPtr, c_bool, ElinaAbstract0Ptr, ElinaDim, '
              'ElinaLinexpr0Ptr, ElinaAbstract0Ptr to the function')

    return a


def elina_abstract0_substitute_texpr(man, destructive, org, dim, texpr, dest):
    """
    Substitute a single ElinaDim of an ElinaAbstract0 by using an ElinaTexpr0.

    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    destructive : c_bool
        Boolean flag.
    org : ElinaAbstract0Ptr
        Pointer to the ElinaAbstract0 which dimension we need to substitute.
    dim : ElinaDim
        ElinaDim we need to substitute.
    texpr : ElinaTexpr0Ptr
        Pointer to the ElinaTexpr0 tht needs to be used for substitution.
    dest : ElinaAbstract0Ptr
        Pointer to an ElinaAbstract0 used as an optional argument.
        If not NULL, the result of the transformation is intersected with dest.
        This is useful for precise backward transformations in lattices like intervals or octagons.

    Returns
    -------
    a : ElinaAbstract0Ptr
        Pointer to the resulting ElinaAbstract0.

    """

    a = None
    try:
        elina_abstract0_substitute_texpr_c = elina_auxiliary_api.elina_abstract0_substitute_texpr
        elina_abstract0_substitute_texpr_c.restype = ElinaAbstract0Ptr
        elina_abstract0_substitute_texpr_c.argtypes = [ElinaManagerPtr, c_bool, ElinaAbstract0Ptr,
                                                       ElinaDim, ElinaTexpr0Ptr, ElinaAbstract0Ptr]
        a = elina_abstract0_substitute_texpr_c(man, destructive, org, dim, texpr, dest)
    except:
        print('Problem with loading/calling "elina_abstract0_substitute_texpr" from "libelinaux.so"')
        print('Make sure you are passing ElinaManagerPtr, c_bool, ElinaAbstract0Ptr, ElinaDim, '
              'ElinaTexpr0Ptr, ElinaAbstract0Ptr to the function')

    return a


def elina_abstract0_apply_dimchange2(man, destructive, a1, dimchange2, project):
    """
    Add then remove/project dimensions, from an ElinaAbstract0, by following the semantics of an ElinaDimchange2.
    If project is true, the newly added dimension are projected on their 0-hyperplane.

    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    destructive : c_bool
        Boolean flag.
    a1 : ElinaAbstract0Ptr
        Pointer to the ElinaAbstract0 which we need to transform.
    dimchange2 : ElinaDimchange2Ptr
        Pointer to the ElinaDimchange2 which semantics we need to follow.
    project : c_bool
        If true, the newly added dimension are projected on their 0-hyperplane.

    Returns
    -------
    a0 : ElinaAbstract0Ptr
        Pointer to the resulting ElinaAbstract0.

    """

    a0 = None
    try:
        elina_abstract0_apply_dimchange2_c = elina_auxiliary_api.elina_abstract0_apply_dimchange2
        elina_abstract0_apply_dimchange2_c.restype = ElinaAbstract0Ptr
        elina_abstract0_apply_dimchange2_c.argtypes = [ElinaManagerPtr, c_bool, ElinaAbstract0Ptr,
                                                       ElinaDimchange2Ptr, c_bool]
        a0 = elina_abstract0_apply_dimchange2_c(man, destructive, a1, dimchange2, project)
    except:
        print('Problem with loading/calling "elina_abstract0_apply_dimchange2" from "libelinaux.so"')
        print('Make sure you are passing ElinaManagerPtr, c_bool, ElinaAbstract0Ptr, ElinaDimchange2Ptr, '
              'c_bool to the function')

    return a0


def elina_abstract0_widening_threshold(man, a1, a2, lincons_array):
    # TODO finish documentation for this function
    """
    !? What does this function do !?
    
    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    a1 : ElinaAbstract0Ptr
        Pointer to the first ElinaAbstract0.
    a2 : ElinaAbstract0Ptr
        Pointer to the second ElinaAbstract0.
    lincons_array : ElinaLincons0ArrayPtr
        Pointer to the ElinaLincons0Array.

    Returns
    -------
    a0 : ElinaAbstract0Ptr
        Pointer to the resulting ElinaAbstract0.

    """

    a0 = None
    try:
        elina_abstract0_widening_threshold_c = elina_auxiliary_api.elina_abstract0_widening_threshold
        elina_abstract0_widening_threshold_c.restype = ElinaAbstract0Ptr
        elina_abstract0_widening_threshold_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr,
                                                         ElinaAbstract0Ptr, ElinaLincons0ArrayPtr]
        a0 = elina_abstract0_widening_threshold_c(man, a1, a2, lincons_array)
    except:
        print('Problem with loading/calling "elina_abstract0_widening_threshold" from "libelinaux.so"')
        print('Make sure you are passing ElinaManagerPtr, ElinaAbstract0Ptr, ElinaAbstract0Ptr to the function')

    return a0
