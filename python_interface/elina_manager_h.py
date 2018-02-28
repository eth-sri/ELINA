# ************************************************************************* #
# elina_manager.h: global manager passed to all functions
# ************************************************************************* #

from elina_coeff_h import *


# ********************************************************************** #
# I. Types */
# ********************************************************************** #


# ====================================================================== #
# I.O General usage
# ====================================================================== #

class ElinaMembuf(Structure):
    """
    ElinaMembuf ctype compatible with elina_membuf_t from elina_manager.h.

    Fields
    ------
    ptr : c_void_p
        Void pointer.
    size : c_size_t

    """

    _fields_ = [('ptr', c_void_p), ('size', c_size_t)]


# ====================================================================== #
# I.1 Identifying functions
# ====================================================================== #

class ElinaFunID(CtypesEnum):
    """
    Enum compatible with elina_funid_t from elina_manager.h.

    Options
    -------
    ELINA_FUNID_UNKNOWN
    ELINA_FUNID_COPY
    ELINA_FUNID_FREE
    ELINA_FUNID_ASIZE
    ELINA_FUNID_MINIMIZE
    ELINA_FUNID_CANONICALIZE
    ELINA_FUNID_HASH
    ELINA_FUNID_APPROXIMATE
    ELINA_FUNID_FPRINT
    ELINA_FUNID_FPRINTDIFF
    ELINA_FUNID_FDUMP
    ELINA_FUNID_SERIALIZE_RAW
    ELINA_FUNID_DESERIALIZE_RAW
    ELINA_FUNID_BOTTOM
    ELINA_FUNID_TOP
    ELINA_FUNID_OF_BOX
    ELINA_FUNID_DIMENSION
    ELINA_FUNID_IS_BOTTOM
    ELINA_FUNID_IS_TOP
    ELINA_FUNID_IS_LEQ
    ELINA_FUNID_IS_EQ
    ELINA_FUNID_IS_DIMENSION_UNCONSTRAINED
    ELINA_FUNID_SAT_INTERVAL
    ELINA_FUNID_SAT_LINCONS
    ELINA_FUNID_SAT_TCONS
    ELINA_FUNID_BOUND_DIMENSION
    ELINA_FUNID_BOUND_LINEXPR
    ELINA_FUNID_BOUND_TEXPR
    ELINA_FUNID_TO_BOX
    ELINA_FUNID_TO_LINCONS_ARRAY
    ELINA_FUNID_TO_TCONS_ARRAY
    ELINA_FUNID_TO_GENERATOR_ARRAY
    ELINA_FUNID_MEET
    ELINA_FUNID_MEET_ARRAY
    ELINA_FUNID_MEET_LINCONS_ARRAY
    ELINA_FUNID_MEET_TCONS_ARRAY
    ELINA_FUNID_JOIN
    ELINA_FUNID_JOIN_ARRAY
    ELINA_FUNID_ADD_RAY_ARRAY
    ELINA_FUNID_ASSIGN_LINEXPR_ARRAY
    ELINA_FUNID_SUBSTITUTE_LINEXPR_ARRAY
    ELINA_FUNID_ASSIGN_TEXPR_ARRAY
    ELINA_FUNID_SUBSTITUTE_TEXPR_ARRAY
    ELINA_FUNID_ADD_DIMENSIONS
    ELINA_FUNID_REMOVE_DIMENSIONS
    ELINA_FUNID_PERMUTE_DIMENSIONS
    ELINA_FUNID_FORGET_ARRAY
    ELINA_FUNID_EXPAND
    ELINA_FUNID_FOLD
    ELINA_FUNID_WIDENING
    ELINA_FUNID_CLOSURE
    ELINA_FUNID_SIZE
    ELINA_FUNID_CHANGE_ENVIRONMENT
    ELINA_FUNID_RENAME_ARRAY
    ELINA_FUNID_SIZE2

    """

    ELINA_FUNID_UNKNOWN = 0
    ELINA_FUNID_COPY = 1
    ELINA_FUNID_FREE = 2
    ELINA_FUNID_ASIZE = 3                   # For avoiding name conflict with ELINA_FUNID_SIZE
    ELINA_FUNID_MINIMIZE = 4
    ELINA_FUNID_CANONICALIZE = 5
    ELINA_FUNID_HASH = 6
    ELINA_FUNID_APPROXIMATE = 7
    ELINA_FUNID_FPRINT = 8
    ELINA_FUNID_FPRINTDIFF = 9
    ELINA_FUNID_FDUMP = 10
    ELINA_FUNID_SERIALIZE_RAW = 11
    ELINA_FUNID_DESERIALIZE_RAW = 12
    ELINA_FUNID_BOTTOM = 13
    ELINA_FUNID_TOP = 14
    ELINA_FUNID_OF_BOX = 15
    ELINA_FUNID_DIMENSION = 16
    ELINA_FUNID_IS_BOTTOM = 17
    ELINA_FUNID_IS_TOP = 18
    ELINA_FUNID_IS_LEQ = 19
    ELINA_FUNID_IS_EQ = 20
    ELINA_FUNID_IS_DIMENSION_UNCONSTRAINED = 21
    ELINA_FUNID_SAT_INTERVAL = 22
    ELINA_FUNID_SAT_LINCONS = 23
    ELINA_FUNID_SAT_TCONS = 24
    ELINA_FUNID_BOUND_DIMENSION = 25
    ELINA_FUNID_BOUND_LINEXPR = 26
    ELINA_FUNID_BOUND_TEXPR = 27
    ELINA_FUNID_TO_BOX = 28
    ELINA_FUNID_TO_LINCONS_ARRAY = 29
    ELINA_FUNID_TO_TCONS_ARRAY = 30
    ELINA_FUNID_TO_GENERATOR_ARRAY = 31
    ELINA_FUNID_MEET = 32
    ELINA_FUNID_MEET_ARRAY = 33
    ELINA_FUNID_MEET_LINCONS_ARRAY = 34
    ELINA_FUNID_MEET_TCONS_ARRAY = 35
    ELINA_FUNID_JOIN = 36
    ELINA_FUNID_JOIN_ARRAY = 37
    ELINA_FUNID_ADD_RAY_ARRAY = 38
    ELINA_FUNID_ASSIGN_LINEXPR_ARRAY = 39
    ELINA_FUNID_SUBSTITUTE_LINEXPR_ARRAY = 40
    ELINA_FUNID_ASSIGN_TEXPR_ARRAY = 41
    ELINA_FUNID_SUBSTITUTE_TEXPR_ARRAY = 42
    ELINA_FUNID_ADD_DIMENSIONS = 43
    ELINA_FUNID_REMOVE_DIMENSIONS = 44
    ELINA_FUNID_PERMUTE_DIMENSIONS = 45
    ELINA_FUNID_FORGET_ARRAY = 46
    ELINA_FUNID_EXPAND = 47
    ELINA_FUNID_FOLD = 48
    ELINA_FUNID_WIDENING = 49
    ELINA_FUNID_CLOSURE = 50
    ELINA_FUNID_SIZE = 51
    ELINA_FUNID_CHANGE_ENVIRONMENT = 52
    ELINA_FUNID_RENAME_ARRAY = 53
    ELINA_FUNID_SIZE2 = 54

elina_name_of_funid = (c_char_p * ElinaFunID.ELINA_FUNID_SIZE2)()


# ====================================================================== #
# I.2 Exceptions
# ====================================================================== #

class ElinaExc(CtypesEnum):
    """
    Enum compatible with elina_funid_t from elina_manager.h.

    Options
    -------
    ELINA_EXC_NONE :                    No exception detected
    ELINA_EXC_TIMEOUT :                 Timeout detected
    ELINA_EXC_OUT_OF_SPACE :            Out of space detected
    ELINA_EXC_OVERFLOW :                Magnitude overflow detected
    ELINA_EXC_INVALID_ARGUMENT :        Invalid arguments
    ELINA_EXC_NOT_IMPLEMENTED :         Not implemented
    ELINA_EXC_SIZE :                    Size    
    
    """

    ELINA_EXC_NONE = 0                  # no exception detected
    ELINA_EXC_TIMEOUT = 1               # timeout detected
    ELINA_EXC_OUT_OF_SPACE = 2          # out of space detected
    ELINA_EXC_OVERFLOW = 3              # magnitude overflow detected
    ELINA_EXC_INVALID_ARGUMENT = 4      # invalid arguments
    ELINA_EXC_NOT_IMPLEMENTED = 5       # not implemented
    ELINA_EXC_SIZE = 6

elina_name_of_exception = (c_char_p * ElinaExc.ELINA_EXC_SIZE)()


class ElinaExcLog(Structure):
    """
    ElinaExcLog ctype compatible with elina_exclog_t from elina_manager.h.
    Serves as an exception log.

    Fields
    ------
    exn : c_uint
    funid : c_uint
    msg : c_char_p
    tail : ElinaExcLogPtr
    
    """

    _fields_ = [('exn', c_uint), ('funid', c_uint), ('msg', c_char_p), ('tail', ElinaExcLogPtr)]

ElinaExcLogPtr = POINTER(ElinaExcLog)


class ElinaResult(Structure):
    """
    ElinaResult ctype compatible with elina_result_t from elina_manager.h.

    Fields
    ------
    exclog : ElinaExcLogPtr
        Pointer to the ElinaExcLog that represents the history of exceptions.
    exn : c_uint
        Enum specifying the exception for the last called function as defined in ElinaExc.
    flag_exact : c_bool
        Flag indicating whether the result is mathematically exact or undefined.
    flag_best : c_bool
        Flag indicating whether the result is best correct approximation or undefined.
    
    """

    _fields_ = [('exclog', ElinaExcLogPtr), ('exn', c_uint), ('flag_exact', c_bool), ('flag_best', c_bool)]


# ====================================================================== #
# I.2 Options
# ====================================================================== #

class ElinaFunOpt(Structure):
    """
    ElinaFunOpt ctype compatible with elina_funopt_t from elina_manager.h.

    Fields
    ------
    algorithm : c_int
        Algorithm selection:
        - 0 is default algorithm;
        - MAX_INT is most accurate available;
        - MIN_INT is most efficient available;
        - otherwise, no accuracy or speed meanin
    timeout : c_size_t
        Above the given computation time (unit !?), the function may abort with the exception flag flag_time_out on.    
    max_object_size : c_size_t
        Defined in abstract object size unit.
    flag_exact_wanted : c_bool
        Returns information whether exact result is wanted.
    flag_best_wanted : c_bool
        Returns information whether best correct approximation result is wanted.

    """

    _fields_ = [('algorithm', c_int), ('timeout', c_size_t), ('max_object_size', c_size_t),
                ('flag_exact_wanted', c_bool), ('flag_best_wanted', c_bool)]

ElinaFunOptPtr = POINTER(ElinaFunOpt)


class ElinaOption(Structure):
    """
    ElinaOption ctype compatible with elina_option_t from elina_manager.h
    
    Fields
    ------
    funopt : ElinaFunOpt * ElinaFunID.ELINA_FUNID_SIZE
        Array of ElinaFunOpt of size ElinaFunID.ELINA_FUNID_SIZE.
    abort_if_exception : c_bool * ElinaExc.ELINA_EXC_SIZE
        Array of bools of size ElinaExc.ELINA_EXC_SIZE determining for which types of exceptions we abort.
    scalar_discr : c_uint
        Enum specifying the preferred type for scalars as defined in ElinaScalarDiscr.
    
    """

    _fields_ = [('funopt', ElinaFunOpt * ElinaFunID.ELINA_FUNID_SIZE),
                ('abort_if_exception', c_bool * ElinaExc.ELINA_EXC_SIZE),
                ('scalar_discr', c_uint)]


# ====================================================================== #
# I.3 Manager
# ====================================================================== #

class ElinaManager(Structure):
    """
    ElinaManager ctype compatible with elina_manager_t from elina_manager.h
    
    Fields
    ------
    library : c_char_p
        Name of the effective library.
    version : c_char_p
        Version of the effective library.
    internal : c_void_p
        Library dependant, should be different for each thread (working space).
    funptr : c_void_p * ElinaFunID.ELINA_FUNID_SIZE
        Array of function pointers of size ElinaFunID.ELINA_FUNID_SIZE initialised by the effective library.
    option : ElinaOption
        Options (in).
    result : ElinaResult
        Exceptions and other indications (out).
    internal_free : CFUNCTYPE(None, c_void_p)
        Deallocation function for internal.
    count : c_size_t
        Reference counter.
    
    """

    _fields_ = [('library', c_char_p), ('version', c_char_p), ('internal', c_void_p),
                ('funptr', c_void_p * ElinaFunID.ELINA_FUNID_SIZE), ('option', ElinaOption),
                ('result', ElinaResult), ('internal_free', CFUNCTYPE(None, c_void_p)), ('count', c_size_t)]

ElinaManagerPtr = POINTER(ElinaManager)
