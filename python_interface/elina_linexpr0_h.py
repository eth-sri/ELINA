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
from elina_coeff_h import *


# ====================================================================== #
# Datatypes */
# ====================================================================== #


class ElinaLinexprDiscr(CtypesEnum):
    """
    Enum compatible with discr field in elina_linexpr0_t from elina_linexpr0.h.
    Discriminant for dense or sparse representation.
    
    """

    ELINA_LINEXPR_DENSE = 0
    ELINA_LINEXPR_SPARSE = 1


class ElinaLinterm(Structure):
    """
    ElinaLinterm ctype compatible with elina_linterm_t from elina_linexpr0.h.
    A term, for use in sparse representation, meant to be an abstract datatype.
    
    Fields
    ------
    dim : ElinaDim
    coeff : ElinaCoeff
    
    """

    _fields_ = [('dim', ElinaDim), ('coeff', ElinaCoeff)]

ElinaLintermPtr = POINTER(ElinaLinterm)


class ElinaLinexprUnion(Union):
    """
    Ctype representing the union field in elina_linexpr0_t from from elina_linexpr0.h.
    Important invariant if sparse representation:
    - Linear terms are sorted in increasing order with respect to their dimension.
    - ELINA_DIM_MAX dimensions are meaningless: they serve as free linterm when a new dimension is needed.
    
    Fields
    ------
    coeff : ElinaCoeffPtr
        Array of coefficients.
    linterm : ElinaLintermPtr
        Array of linear terms.

    """

    _fields_ = [('coeff', ElinaCoeffPtr), ('linterm', ElinaLintermPtr)]


class ElinaLinexpr0(Structure):
    """
    ElinaLinexpr0 ctype compatible with elina_linexpr0_t from elina_linexpr0.h.   
    A linear expression, meant to be an abstract datatype.
    
    Fields
    ------
    cst : ElinaCoeff
        Constant.
    discr : c_uint
        Discriminant for the array type in the Union, according to ElinaLinexprDiscr enum.
    size : c_size_t
        Size of the array 
    p : ElinaLinexprUnion
        Union field representing an array of coefficients or an array of linear terms.
    
    """

    _fields_ = [('cst', ElinaCoeff), ('discr', c_uint), ('size', c_size_t), ('p', ElinaLinexprUnion)]

ElinaLinexpr0Ptr = POINTER(ElinaLinexpr0)
ElinaLinexpr0Array = POINTER(ElinaLinexpr0Ptr)


class ElinaLinexprType(CtypesEnum):
    """
    Enum compatible with elina_linexpr_type_t from elina_linexpr0.h.
    
    - An interval linear expression is the more general form.
    - A quasilinear expression is such that the only non-scalar coefficient is the constant coefficient.
    - A linear expression contains no non-scalar coefficients.
    
    """

    ELINA_LINEXPR_INTLINEAR = 0
    ELINA_LINEXPR_QUASILINEAR = 1
    ELINA_LINEXPR_LINEAR = 2


class ElinaCoefftag(CtypesEnum):
    """
    Enum compatible with elina_coefftag_t from elina_linexpr0.h.
    
    Options
    -------
    ELINA_COEFF :               Waiting for an ElinaCoeffPtr and a dimension
    ELINA_COEFF_S :             Waiting for an ElinaScalarPtr and a dimension
    ELINA_COEFF_S_MPQ :         Waiting for a Mpq_t and a dimension
    ELINA_COEFF_S_MPFR :        Waiting for a Mpfr_t and a dimension
    ELINA_COEFF_S_INT :         Waiting for a c_int and a dimension
    ELINA_COEFF_S_FRAC :        Waiting for 2 c_int-s and a dimension
    ELINA_COEFF_S_DOUBLE :      Waiting for a c_double and a dimension
    ELINA_COEFF_I :             Waiting for an ElinaIntervalPtr and a dimension
    Elina_COEFF_SCALAR :        Waiting for 2 ElinaScalarPtr-s and a dimension
    ELINA_COEFF_I_MPQ :         Waiting for 2 Mpq_t-s and a dimension
    ELINA_COEFF_I_MPFR :        Waiting for 2 Mpfr_t-s and a dimension
    ELINA_COEFF_I_INT :         Waiting for 2 c_int-s and a dimension
    ELINA_COEFF_I_FRAC :        Waiting for 4 c_int-s and a dimension
    ELINA_COEFF_I_DOUBLE :      Waiting for 2 c_double-s and a dimension
    ELINA_CST :                 Waiting for an ElinaCoeffPtr
    ELINA_CST_S :               Waiting for an ElinaScalarPtr
    ELINA_CST_S_MPQ :           Waiting for a Mpq_t
    ELINA_CST_S_MPFR :          Waiting for a Mpfr_t
    ELINA_CST_S_INT :           Waiting for a c_int
    ELINA_CST_S_FRAC :          Waiting for 2 c_int-s
    ELINA_CST_S_DOUBLE :        Waiting for a c_double
    ELINA_CST_I :               Waiting for an ElinaIntervalPtr
    ELINA_CST_I_SCALAR :        Waiting for 2 ElinaScalarPtr-s
    ELINA_CST_I_MPQ :           Waiting for 2 Mpq_t-s
    ELINA_CST_I_MPFR :          Waiting for 2 Mpfr_t-s
    ELINA_CST_I_INT :           Waiting for 2 c_int-s
    ELINA_CST_I_FRAC :          Waiting for 4 c_int-s
    ELINA_CST_I_DOUBLE :        Waiting for 2 c_doub
    ELINA_END :                 END
    
    """

    ELINA_COEFF = 0
    ELINA_COEFF_S = 1
    ELINA_COEFF_S_MPQ = 2
    ELINA_COEFF_S_MPFR = 3
    ELINA_COEFF_S_INT = 4
    ELINA_COEFF_S_FRAC = 5
    ELINA_COEFF_S_DOUBLE = 6
    ELINA_COEFF_I = 7
    Elina_COEFF_SCALAR = 8
    ELINA_COEFF_I_MPQ = 9
    ELINA_COEFF_I_MPFR = 10
    ELINA_COEFF_I_INT = 11
    ELINA_COEFF_I_FRAC = 12
    ELINA_COEFF_I_DOUBLE = 13
    ELINA_CST = 14
    ELINA_CST_S = 15
    ELINA_CST_S_MPQ = 16
    ELINA_CST_S_MPFR = 17
    ELINA_CST_S_INT = 18
    ELINA_CST_S_FRAC = 19
    ELINA_CST_S_DOUBLE = 20
    ELINA_CST_I = 21
    ELINA_CST_I_SCALAR = 22
    ELINA_CST_I_MPQ = 23
    ELINA_CST_I_MPFR = 24
    ELINA_CST_I_INT = 25
    ELINA_CST_I_FRAC = 26
    ELINA_CST_I_DOUBLE = 27
    ELINA_END = 28
