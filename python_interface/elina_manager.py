#
#
#  This source file is part of ELINA (ETH LIbrary for Numerical Analysis).
#  ELINA is Copyright © 2018 Department of Computer Science, ETH Zurich
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


from elina_manager_h import *
from elina_auxiliary_imports import *


# ********************************************************************** #
# II. User Functions
# ********************************************************************** #

def elina_manager_clear_exclog(man):
    """
    Erase the current exception log of an ElinaManager
    
    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ELinaManager which exception log needs to be erased.

    Returns
    -------
    None

    """

    try:
        elina_manager_clear_exclog_c = elina_auxiliary_api.elina_manager_clear_exclog
        elina_manager_clear_exclog_c.restype = None
        elina_manager_clear_exclog_c.argtypes = [ElinaManagerPtr]
        elina_manager_clear_exclog_c(man)
    except:
        print('Problem with loading/calling "elina_manager_clear_exclog" from "libelinaux.so"')
        print('Make sure you are passing ElinaManagerPtr to the function')


def elina_manager_free(man):
    """
    Dereference the counter of an ElinaManager and possibly free internal field if it is not yet put to NULL.

    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ELinaManager that needs to be freed.

    Returns
    -------
    None

    """

    try:
        elina_manager_free_c = elina_auxiliary_api.elina_manager_free
        elina_manager_free_c.restype = None
        elina_manager_free_c.argtypes = [ElinaManagerPtr]
        elina_manager_free_c(man)
    except:
        print('Problem with loading/calling "elina_manager_free" from "libelinaux.so"')
        print('Make sure you are passing ElinaManagerPtr to the function')


def elina_manager_get_library(man):
    """
    Get method for obtaining the library name out of an ElinaManager.

    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ELinaManager which library name we want to get.

    Returns
    -------
    library : c_char_p
        Returns the name of the library as string.

    """

    library = None
    try:
        elina_manager_get_library_c = elina_auxiliary_api.elina_manager_get_library
        elina_manager_get_library_c.restype = c_char_p
        elina_manager_get_library_c.argtypes = [ElinaManagerPtr]
        library = elina_manager_get_library_c(man)
    except:
        print('Problem with loading/calling "elina_manager_get_library" from "libelinaux.so"')
        print('Make sure you are passing ElinaManagerPtr to the function')

    return library


def elina_manager_get_version(man):
    """
    Get method for obtaining the library version out of an ElinaManager.
    
    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ELinaManager which library version we want to get.
    
    Returns
    -------
    library : c_char_p
        Returns the version of the library as string.
    
    """
    version = None
    try:
        elina_manager_get_version_c = elina_auxiliary_api.elina_manager_get_version
        elina_manager_get_version_c.restype = c_char_p
        elina_manager_get_version_c.argtypes = [ElinaManagerPtr]
        version = elina_manager_get_version_c(man)
    except:
        print('Problem with loading/calling "elina_manager_get_version" from "libelinaux.so"')
        print('Make sure you are passing ElinaManagerPtr to the function')

    return version


def elina_manager_get_funopt(man, funid):
    """
    Get method for obtaining the function options out of an ElinaManager given a function ID.

    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ELinaManager which function options we want to get.
    funid : c_uint
        ID of the function which options we want to get.

    Returns
    -------
    fun_opt : ElinaFunOpt
        Options of the function specified with the function id.

    """

    fun_opt = None
    try:
        elina_manager_get_funopt_c = elina_auxiliary_api.elina_manager_get_funopt
        elina_manager_get_funopt_c.restype = ElinaFunOpt
        elina_manager_get_funopt_c.argtypes = [ElinaManagerPtr, c_uint]
        fun_opt = elina_manager_get_funopt_c(man, funid)
    except:
        print('Problem with loading/calling "elina_manager_get_funopt" from "libelinaux.so"')
        print('Make sure you are passing ElinaManagerPtr, c_uint to the function')

    return fun_opt


def elina_manager_get_abort_if_exception(man, exn):
    """
    Get method for obtaining the abort_if_exception flag from an ElinaManager given an exception ID.

    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ELinaManager which abort_if_exception flag we want to get.
    exn : c_uint
        ID of the exception which options we want to get.

    Returns
    -------
    abort_if_exception : c_bool
        Value of the abort_if_exception flag for the given ElinaManager and exception ID.

    """

    abort_if_exception = None
    try:
        elina_manager_get_abort_if_exception_c = elina_auxiliary_api.elina_manager_get_abort_if_exception
        elina_manager_get_abort_if_exception_c.restype = c_bool
        elina_manager_get_abort_if_exception_c.argtypes = [ElinaManagerPtr, c_uint]
        abort_if_exception = elina_manager_get_abort_if_exception_c(man, exn)
    except:
        print('Problem with loading/calling "elina_manager_get_abort_if_exception" from "libelinaux.so"')
        print('Make sure you are passing ElinaManagerPtr, c_uint to the function')

    return abort_if_exception


def elina_manager_get_exception(man):
    """
    Get the last exception raised from an ElinaManager.

    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ELinaManager which last exception raised we want to get.

    Returns
    -------
    exn : c_uint
        ID of the last exception raised.

    """

    exn = None
    try:
        elina_manager_get_exception_c = elina_auxiliary_api.elina_manager_get_exception
        elina_manager_get_exception_c.restype = c_uint
        elina_manager_get_exception_c.argtypes = [ElinaManagerPtr]
        exn = elina_manager_get_exception_c(man)
    except:
        print('Problem with loading/calling "elina_manager_get_exception" from "libelinaux.so"')
        print('Make sure you are passing ElinaManagerPtr to the function')

    return exn


def elina_manager_get_exclog(man):
    """
    Get the full exception log of an ElinaManager.

    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ELinaManager which exception log we want to get.

    Returns
    -------
    exc_log : ElinaExcLogPtr
        Pointer to the ElinaExcLog.

    """

    exc_log = None
    try:
        elina_manager_get_exclog_c = elina_auxiliary_api.elina_manager_get_exclog
        elina_manager_get_exclog_c.restype = ElinaExcLogPtr
        elina_manager_get_exclog_c.argtypes = [ElinaManagerPtr]
        exc_log = elina_manager_get_exclog_c(man)
    except:
        print('Problem with loading/calling "elina_manager_get_exclog" from "libelinaux.so"')
        print('Make sure you are passing ElinaManagerPtr to the function')

    return exc_log


def elina_manager_get_flag_exact(man):
    """
    Get method for obtaining the flag_exact from an ElinaManager.

    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ELinaManager which flag_exact we want to get.

    Returns
    -------
    flag_exact : c_bool
        Value of the flag_exact for the given ElinaManager.

    """

    flag_exact = None
    try:
        elina_manager_get_flag_exact_c = elina_auxiliary_api.elina_manager_get_flag_exact
        elina_manager_get_flag_exact_c.restype = c_bool
        elina_manager_get_flag_exact_c.argtypes = [ElinaManagerPtr]
        flag_exact = elina_manager_get_flag_exact_c(man)
    except:
        print('Problem with loading/calling "elina_manager_get_flag_exact" from "libelinaux.so"')
        print('Make sure you are passing ElinaManagerPtr to the function')

    return flag_exact


def elina_manager_get_flag_best(man):
    """
    Get method for obtaining the flag_best from an ElinaManager.

    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ELinaManager which best flag_best we want to get.

    Returns
    -------
    flag_best : c_bool
        Value of the flag_best for the given ElinaManager.

    """

    flag_best = None
    try:
        elina_manager_get_flag_best_c = elina_auxiliary_api.elina_manager_get_flag_best
        elina_manager_get_flag_best_c.restype = c_bool
        elina_manager_get_flag_best_c.argtypes = [ElinaManagerPtr]
        flag_best = elina_manager_get_flag_best_c(man)
    except:
        print('Problem with loading/calling "elina_manager_get_flag_best" from "libelinaux.so"')
        print('Make sure you are passing ElinaManagerPtr to the function')

    return flag_best


def elina_funopt_init(fopt):
    """
    Initialise an ElinaFunOpt.

    Parameters
    ----------
    fopt : ElinaFunOptPtr
        Pointer to the ElinaFunOpt that needs to be initialised.

    Returns
    -------
    None
    
    """

    try:
        elina_manager_get_flag_best_c = elina_auxiliary_api.elina_manager_get_flag_best
        elina_manager_get_flag_best_c.restype = None
        elina_manager_get_flag_best_c.argtypes = [ElinaFunOptPtr]
        elina_manager_get_flag_best_c(fopt)
    except:
        print('Problem with loading/calling "elina_manager_get_flag_best" from "libelinaux.so"')
        print('Make sure you are passing ElinaFunOptPtr to the function')


def elina_manager_set_funopt(man, funid, fopt):
    """
    Set the ElinaFunOpt for a given function ID and an ElinaManager.

    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager which function options we want to set.
    funid : c_uint
        ID of the function which options we want to set.
    fopt : ElinaFunOptPtr
        Pointer to the ElinaFunOpt that needs to be used as new value.

    Returns
    -------
    None

    """

    try:
        elina_manager_set_funopt_c = elina_auxiliary_api.elina_manager_set_funopt
        elina_manager_set_funopt_c.restype = None
        elina_manager_set_funopt_c.argtypes = [ElinaManagerPtr, c_uint, ElinaFunOptPtr]
        elina_manager_set_funopt_c(man, funid, fopt)
    except:
        print('Problem with loading/calling "elina_manager_set_funopt" from "libelinaux.so"')
        print('Make sure you are passing ElinaManagerPtr, c_uint, ElinaFunOptPtr to the function')


def elina_manager_set_abort_if_exception(man, exn, flag):
    """
    Set the abort_if_exception flag for a given exception ID and an ElinaManager.

    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager which function options we want to set.
    exn : c_uint
        ID of the exception which abort_if_exception flag we want to set.
    flag : c_bool
        New value of the abort_if_exception flag.

    Returns
    -------
    None

    """

    try:
        elina_manager_set_abort_if_exception_c = elina_auxiliary_api.elina_manager_set_abort_if_exception
        elina_manager_set_abort_if_exception_c.restype = None
        elina_manager_set_abort_if_exception_c.argtypes = [ElinaManagerPtr, c_uint, c_bool]
        elina_manager_set_abort_if_exception_c(man, exn, flag)
    except:
        print('Problem with loading/calling "elina_manager_set_abort_if_exception" from "libelinaux.so"')
        print('Make sure you are passing ElinaManagerPtr, c_uint, c_bool to the function')


def elina_fpy_init():
    """
    Tries to set the FPU rounding-mode towards +oo, returns true if successful.
    
    Returns
    -------
    result : c_bool
        Returns if the operation was successful.

    """

    result = None
    try:
        elina_fpy_init_c = elina_auxiliary_api.elina_fpy_init
        elina_fpy_init_c.restype = c_bool
        result = elina_fpy_init_c()
    except:
        print('Problem with loading/calling "elina_fpy_init" from "libelinaux.so"')

    return result


# ********************************************************************** #
# III. Implementor Functions
# ********************************************************************** #

def elina_manager_alloc(library, version, internal, internal_free):
    """
    Allocates an ElinaManager.
    
    Parameters
    ----------
    library : c_char_p
        Name of the library.
    version : c_char_p
        Version of the library.
    internal : c_void_p
        Library dependant.
    internal_free : CFUNCTYPE(None, c_void_p)
        Deallocation function for internal.

    Returns
    -------
    man : ElinaManagerPtr
        Pointer to the newly allocated ElinaManager.

    """

    man = None
    try:
        elina_manager_alloc_c = elina_auxiliary_api.elina_manager_alloc
        elina_manager_alloc_c.restype = ElinaManagerPtr
        elina_manager_alloc_c.argtypes = [c_char_p, c_char_p, c_void_p, CFUNCTYPE(None, c_void_p)]
        man = elina_manager_alloc_c(library, version, internal, internal_free)
    except:
        print('Problem with loading/calling "elina_manager_alloc" from "libelinaux.so"')
        print('Make sure you are passing c_char_p, c_char_p, c_void_p, CFUNCTYPE(None, c_void_p) to the function')

    return man


def elina_manager_copy(man2):
    """
    Copy an ElinaManager.
    
    Parameters
    ----------
    man2 : ElinaManagerPtr
        Source.

    Returns
    -------
    man1 : ElinaManagerPtr
        Destination.
        
    """

    man1 = None
    try:
        elina_manager_copy_c = elina_auxiliary_api.elina_manager_copy
        elina_manager_copy_c.restype = ElinaManagerPtr
        elina_manager_copy_c.argtypes = [ElinaManagerPtr]
        man1 = elina_manager_copy_c(man2)
    except:
        print('Problem with loading/calling "elina_manager_copy" from "libelinaux.so"')
        print('Make sure you are passing ElinaManagerPtr to the function')

    return man1


def elina_manager_raise_exception(man, exn, funid, msg):
    """
    Raise an exception and put fields flag_exact and flag_best of the result field in the ElinaManager to false.

    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    exn : c_uint
        Exception ID.
    funid : c_uint
        Function ID.
    msg : c_char_p
        Message.

    Returns
    -------
    None
    
    """

    try:
        elina_manager_raise_exception_c = elina_auxiliary_api.elina_manager_raise_exception
        elina_manager_raise_exception_c.restype = None
        elina_manager_raise_exception_c.argtypes = [ElinaManagerPtr, c_uint, c_uint, c_char_p]
        elina_manager_raise_exception_c(man, exn, funid, msg)
    except:
        print('Problem with loading/calling "elina_manager_raise_exception" from "libelinaux.so"')
        print('Make sure you are passing ElinaManagerPtr, c_uint, c_uint, c_char_p to the function')


def elina_exc_cons(exn, funid, msg, tail):
    # TODO finish documentation for this function
    """
    !?What does this function do!?
    
    Parameters
    ----------
    exn : c_uint
        Exception ID.
    funid : c_uint
        Function ID.
    msg : c_char_p
        Message.
    tail : ElinaExcLogPtr
        Pointer to the ElinaExcLog.

    Returns
    -------
    None

    """

    try:
        elina_exc_cons_c = elina_auxiliary_api.elina_exc_cons
        elina_exc_cons_c.restype = None
        elina_exc_cons_c.argtypes = [c_uint, c_uint, c_char_p, ElinaExcLogPtr]
        elina_exc_cons_c(exn, funid, msg, tail)
    except:
        print('Problem with loading/calling "elina_exc_cons" from "libelinaux.so"')
        print('Make sure you are passing c_uint, c_uint, c_char_p, ElinaExcLogPtr to the function')


def elina_exclog_free(head):
    """
    Free an ElinaExcLog.
    
    Parameters
    ----------
    head : ElinaExcLogPtr
        Pointer to the ElinaExcLog that needs to be freed.

    Returns
    -------
    None

    """

    try:
        elina_exclog_free_c = elina_auxiliary_api.elina_exclog_free
        elina_exclog_free_c.restype = None
        elina_exclog_free_c.argtypes = [ElinaExcLogPtr]
        elina_exclog_free_c(head)
    except:
        print('Problem with loading/calling "elina_exclog_free" from "libelinaux.so"')
        print('Make sure you are passing ElinaExcLogPtr to the function')
