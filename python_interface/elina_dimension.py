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


from elina_dimension_h import *

# ====================================================================== #
# Functions
# ====================================================================== #


# ---------------------------------------------------------------------- #
# elina_dimchange_t
# ---------------------------------------------------------------------- #


def elina_dimchange_init(dimchange, intdim, realdim):
    """
    Initialise a given ElinaDimchange.
    
    Parameters
    ----------
    dimchange : ElinaDimchangePtr
        Pointer to the ElinaDimchange that needs to be initialised.
    intdim : c_size_t
        Number of integer dimensions.
    realdim : c_size_t
        Number of real dimensions.

    Returns
    -------
    None

    """

    try:
        elina_dimchange_init_c = elina_auxiliary_api.elina_dimchange_init
        elina_dimchange_init_c.restype = None
        elina_dimchange_init_c.argtypes = [ElinaDimchangePtr, c_size_t, c_size_t]
        elina_dimchange_init_c(dimchange, intdim, realdim)
    except:
        print('Problem with loading/calling "elina_dimchange_init" from "libelinaux.so"')
        print('Make sure you are passing ElinaDimchangePtr, c_size_t, c_size_t to the function')


def elina_dimchange_alloc(intdim, realdim):
    """
    Allocate a new ElinaDimchange.
    
    Parameters
    ----------
    intdim : c_size_t
        Number of integer dimensions.
    realdim : c_size_t
        Number of real dimensions.

    Returns
    -------
    dimchange : ElinaDimchangePtr
        Pointer to the newly allocated ElinaDimchange

    """

    dimchange = None
    try:
        elina_dimchange_alloc_c = elina_auxiliary_api.elina_dimchange_alloc
        elina_dimchange_alloc_c.restype = ElinaDimchangePtr
        elina_dimchange_alloc_c.argtypes = [c_size_t, c_size_t]
        dimchange = elina_dimchange_alloc_c(intdim, realdim)
    except:
        print('Problem with loading/calling "elina_dimchange_alloc" from "libelinaux.so"')
        print('Make sure you are passing c_size_t, c_size_t to the function')

    return dimchange


def elina_dimchange_fprint(stream, dimchange):
    """
    Print an ElinaDimchange onto a given stream.
    
    Parameters
    ----------
    stream : c_void_p
        Stream on which to print.
    dimchange : ElinaDimchangePtr
        Pointer to the ElinaDimchange that needs to be printed.

    Returns
    -------
    None

    """

    try:
        elina_dimchange_fprint_c = elina_auxiliary_api.elina_dimchange_fprint
        elina_dimchange_fprint_c.restype = None
        elina_dimchange_fprint_c.argtypes = [c_void_p, ElinaDimchangePtr]
        elina_dimchange_fprint_c(stream, dimchange)
    except:
        print('Problem with loading/calling "elina_dimchange_fprint" from "libelinaux.so"')
        print('Make sure you are passing c_void_p, ElinaDimchangePtr to the function')


def elina_dimchange_add_invert(dimchange):
    """
    Assuming that dimchange is a transformation for add_dimensions, invert it to obtain the inverse transformation using remove_dimensions.
    
    Parameters
    ----------
    dimchange : ElinaDimchangePtr
        Pointer to the ElinaDimchange that needs to be inverted.

    Returns
    -------
    None

    """

    try:
        elina_dimchange_add_invert_c = elina_auxiliary_api.elina_dimchange_add_invert
        elina_dimchange_add_invert_c.restype = None
        elina_dimchange_add_invert_c.argtypes = [ElinaDimchangePtr]
        elina_dimchange_add_invert_c(dimchange)
    except:
        print('Problem with loading/calling "elina_dimchange_add_invert" from "libelinaux.so"')
        print('Make sure you are passing ElinaDimchangePtr to the function')


# ---------------------------------------------------------------------- #
# elina_dimchange2_t
# ---------------------------------------------------------------------- #


def elina_dimchange2_clear(dimchange2):
    """
    Clear an ElinaDimchange2.
    
    Parameters
    ----------
    dimchange2 : ElinaDimchange2Ptr
        Pointer to the ElinaDimchange2 that needs to be cleared.

    Returns
    -------
    None

    """

    try:
        elina_dimchange2_clear_c = elina_auxiliary_api.elina_dimchange2_clear
        elina_dimchange2_clear_c.restype = None
        elina_dimchange2_clear_c.argtypes = [ElinaDimchange2Ptr]
        elina_dimchange2_clear_c(dimchange2)
    except:
        print('Problem with loading/calling "elina_dimchange2_clear" from "libelinaux.so"')
        print('Make sure you are passing ElinaDimchange2Ptr to the function')


def elina_dimchange_free(dimchange):
    """
    Free an ElinaDimchange.
    
    Parameters
    ----------
    dimchange : ElinaDimchangePtr
        Pointer to the ElinaDimchange that needs to be freed.

    Returns
    -------
    None

    """

    try:
        elina_dimchange_free_c = elina_auxiliary_api.elina_dimchange_free
        elina_dimchange_free_c.restype = None
        elina_dimchange_free_c.argtypes = [ElinaDimchangePtr]
        elina_dimchange_free_c(dimchange)
    except:
        print('Problem with loading/calling "elina_dimchange_free" from "libelinaux.so"')
        print('Make sure you are passing ElinaDimchangePtr to the function')

def elina_dimchange2_free(dimchange2):
    """
    Free an ElinaDimchange2.
    
    Parameters
    ----------
    dimchange2 : ElinaDimchange2Ptr
        Pointer to the ElinaDimchange2 that needs to be freed.

    Returns
    -------
    None

    """

    try:
        elina_dimchange2_free_c = elina_auxiliary_api.elina_dimchange2_free
        elina_dimchange2_free_c.restype = None
        elina_dimchange2_free_c.argtypes = [ElinaDimchange2Ptr]
        elina_dimchange2_free_c(dimchange2)
    except:
        print('Problem with loading/calling "elina_dimchange2_free" from "libelinaux.so"')
        print('Make sure you are passing ElinaDimchange2Ptr to the function')





def elina_dimchange2_fprint(stream, dimchange2):
    """
    Print an ElinaDimchange2 onto a given stream.
    
    Parameters
    ----------
    stream : c_void_p
        Stream onto which to print.
    dimchange2 : ElinaDimchange2Ptr
        Pointer to the ElinaDimchange2 that needs to be printed.

    Returns
    -------
    None

    """

    try:
        elina_dimchange2_fprint_c = elina_auxiliary_api.elina_dimchange2_fprint
        elina_dimchange2_fprint_c.restype = None
        elina_dimchange2_fprint_c.argtypes = [c_void_p, ElinaDimchange2Ptr]
        elina_dimchange2_fprint_c(stream, dimchange2)
    except:
        print('Problem with loading/calling "elina_dimchange2_fprint" from "libelinaux.so"')
        print('Make sure you are passing c_void_p, ElinaDimchange2Ptr to the function')


# ---------------------------------------------------------------------- #
# elina_dimperm_t
# ---------------------------------------------------------------------- #


def elina_dimperm_init(dimperm, size):
    """
    Initialise an ElinaDimperm of a given size.
    
    Parameters
    ----------
    dimperm : ElinaDimpermPtr
        Pointer to the ElinaDimperm that needs to be initialised.
    size : c_size_t
        Size.

    Returns
    -------
    None

    """

    try:
        elina_dimperm_init_c = elina_auxiliary_api.elina_dimperm_init
        elina_dimperm_init_c.restype = None
        elina_dimperm_init_c.argtypes = [ElinaDimpermPtr, c_size_t]
        elina_dimperm_init_c(dimperm, size)
    except:
        print('Problem with loading/calling "elina_dimperm_init" from "libelinaux.so"')
        print('Make sure you are passing ElinaDimpermPtr, c_size_t to the function')


def elina_dimperm_alloc(size):
    """
    Allocate an ElinaDimperm.
    
    Parameters
    ----------
    size : c_size_t
        Size.

    Returns
    -------
    dimperm : ElinaDimpermPtr
        Pointer to the newly allocated ElinaDimperm.

    """

    dimperm = None
    try:
        elina_dimperm_alloc_c = elina_auxiliary_api.elina_dimperm_alloc
        elina_dimperm_alloc_c.restype = ElinaDimpermPtr
        elina_dimperm_alloc_c.argtypes = [c_size_t]
        dimperm = elina_dimperm_alloc_c(size)
    except:
        print('Problem with loading/calling "elina_dimperm_alloc" from "libelinaux.so"')
        print('Make sure you are passing c_size_t to the function')

    return dimperm


def elina_dimperm_fprint(stream, dimperm):
    """
    Print an ElinaDimperm onto a given stream.
    
    Parameters
    ----------
    stream : c_void_p
        Stream onto which to print.
    dimperm : ElinaDimpermPtr
        Pointer to the ElinaDimperm that needs to be printed.

    Returns
    -------
    None

    """

    try:
        elina_dimperm_fprint_c = elina_auxiliary_api.elina_dimperm_fprint
        elina_dimperm_fprint_c.restype = None
        elina_dimperm_fprint_c.argtypes = [c_void_p, ElinaDimpermPtr]
        elina_dimperm_fprint_c(stream, dimperm)
    except:
        print('Problem with loading/calling "elina_dimperm_fprint" from "libelinaux.so"')
        print('Make sure you are passing c_void_p, ElinaDimpermPtr to the function')


def elina_dimperm_set_id(perm):
    """
    Set a given ElinaDimperm to the identity permutation.
    
    Parameters
    ----------
    perm : ElinaDimpermPtr
        Pointer to the ElinaDimperm that needs to be set to the identity permutation.

    Returns
    -------
    None

    """

    try:
        elina_dimperm_set_id_c = elina_auxiliary_api.elina_dimperm_set_id
        elina_dimperm_set_id_c.restype = None
        elina_dimperm_set_id_c.argtypes = [ElinaDimpermPtr]
        elina_dimperm_set_id_c(perm)
    except:
        print('Problem with loading/calling "elina_dimperm_set_id" from "libelinaux.so"')
        print('Make sure you are passing ElinaDimpermPtr to the function')


def elina_dimperm_compose(perm, perm1, perm2):
    """
    Compose 2 permutations and store the result into an already allocated permutation.
    The size of both permutations should be equal.
    At exit, we have perm.dim[i] = perm2.dim[perm1.dim[i]].
    
    Parameters
    ----------
    perm : ElinaDimpermPtr
        Pointer to the ElinaDimperm where we want to store the composed permutation
    perm1 : ElinaDimpermPtr
        Pointer to the first permutation that needs to be composed.
    perm2 : ElinaDimpermPtr
        Pointer to the second permutation that needs to be composed.

    Returns
    -------
    None

    """

    try:
        elina_dimperm_compose_c = elina_auxiliary_api.elina_dimperm_compose
        elina_dimperm_compose_c.restype = None
        elina_dimperm_compose_c.argtypes = [ElinaDimpermPtr, ElinaDimpermPtr, ElinaDimpermPtr]
        elina_dimperm_compose_c(perm, perm1, perm2)
    except:
        print('Problem with loading/calling "elina_dimperm_set_id" from "libelinaux.so"')
        print('Make sure you are passing ElinaDimpermPtr, ElinaDimpermPtr, ElinaDimpermPtr to the function')


def elina_dimperm_invert(nperm, perm):
    """
    Invert a permutation perm and store it in an already allocated permutation.
    The size of both permutations are supposed to be equal.
    
    Parameters
    ----------
    nperm : ElinaDimpermPtr
        Destination.
    perm : ElinaDimpermPtr
        Source.

    Returns
    -------
    None

    """

    try:
        elina_dimperm_invert_c = elina_auxiliary_api.elina_dimperm_invert_c
        elina_dimperm_invert_c.restype = None
        elina_dimperm_invert_c.argtypes = [ElinaDimpermPtr, ElinaDimpermPtr]
        elina_dimperm_invert_c(nperm, perm)
    except:
        print('Problem with loading/calling "elina_dimperm_invert" from "libelinaux.so"')
        print('Make sure you are passing ElinaDimpermPtr, ElinaDimpermPtr to the function')
