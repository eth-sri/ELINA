from elina_auxiliary_imports import *

# ====================================================================== #
# Datatypes
# ====================================================================== #


ElinaDim = c_uint
ElinaDimPtr = POINTER(ElinaDim)
ELINA_DIM_MAX = c_uint(-1)


class ElinaDimension(Structure):
    """
    ElinaDimension ctype compatible with elina_dimension_t from elina_dimension.h.
    Datatype for specifying the dimensionality of an abstract value.
    
    Fields:
    -------
    intdim : c_size_t
    realdim : c_size_t
    
    """

    _fields_ = [('intdim', c_size_t), ('realdim', c_size_t)]


class ElinaDimchange(Structure):
    """
    ElinaDimchange ctype compatible with elina_dimchange_t from elina_dimension.h.
    Datatype for specifying change of dimension (addition or removal).
    
    Fields:
    -------
    dim : ElinaDimPtr
        Assumed to be an array of size intdim+realdim.
    intdim : c_size_t
        Number of integer dimensions to add/remove.
    realdim : c_size_t
        Number of real dimensions to add/remove.

    """

    _fields_ = [('dim', ElinaDimPtr), ('intdim', c_size_t), ('realdim', c_size_t)]

ElinaDimchangePtr = POINTER(ElinaDimchange)


class ElinaDimchange2(Structure):
    """
    ElinaDimchange2 ctype compatible with elina_dimchange2_t from elina_dimension.h.
    Datatype for specifying double changes of dimensions (combination of addition and then removal).
    Used by level 1 function change_environment.
    
    Fields:
    -------
    add : ElinaDimchangePtr
        If not NULL, specifies the adding new dimensions.
    remove : ElinaDimchangePtr
        If not NULL, specifies the removal of dimensions.
        
    """

    _fields_ = [('add', ElinaDimchangePtr), ('remove', ElinaDimchangePtr)]

ElinaDimchange2Ptr = POINTER(ElinaDimchange2)


class ElinaDimperm(Structure):
    """
    ElinaDimperm ctype compatible with elina_dimperm_t from elina_dimension.h.
    Datatype for permutations, representint the permutation i -> dimperm.p[i] for 0<=i<dimperm.size.
    
    Fields:
    -------
    dim : ElinaDimPtr
        Array assumed to be of size size.
    size : c_size_t
    
    """

    _fields_ = [('dim', ElinaDimPtr), ('size', c_size_t)]

ElinaDimpermPtr = POINTER(ElinaDimperm)
