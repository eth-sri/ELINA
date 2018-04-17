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
