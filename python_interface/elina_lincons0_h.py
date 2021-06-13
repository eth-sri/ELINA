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


from elina_linexpr0 import *


# ************************************************************************* #
# elina_lincons0.h: linear constraints and arrays
# ************************************************************************* #

# ====================================================================== #
# Datatypes */
# ====================================================================== #


class ElinaConstyp(CtypesEnum):
    """
    Enum compatible with elina_linconstyp_t from elina_lincons0.h.
    
    Options
    -------
    ELINA_CONS_EQ :         Equality constraint
    ELINA_CONS_SUPEQ :      >= constraint
    ELINA_CONS_SUP :        > constraint
    ELINA_CONS_EQMOD :      Congruence equality constraint
    ELINA_CONS_DISEQ :      Disequality constraint
    
    """

    ELINA_CONS_EQ = 0
    ELINA_CONS_SUPEQ = 1
    ELINA_CONS_SUP = 2
    ELINA_CONS_EQMOD = 3
    ELINA_CONS_DISEQ = 4


class ElinaLincons0(Structure):
    """
    ElinaLincons0 ctype compatible with elina_lincons0_t from elina_lincons0.h.
    
    Fields
    ------
    linexpr0 : ElinaLinexpr0Ptr
        Expression.
    constyp : c_uint
        Type of constraint.
    scalar : ElinaScalarPtr
        For EQMOD constraint indicates the modulo, otherwise NULL.

    """

    _fields_ = [('linexpr0', ElinaLinexpr0Ptr), ('constyp', c_uint), ('scalar', ElinaScalarPtr)]

ElinaLincons0Ptr = POINTER(ElinaLincons0)


class ElinaLincons0Array(Structure):
    """
    ElinaLinconsArray0 ctype compatible with elina_lincons0_array_t from elina_lincons0.h.

    Fields
    ------
    p : ElinaLincons0Ptr
        Expression.
    size : c_size_t
        Size of the array.
        
    """

    _fields_ = [('p', ElinaLincons0Ptr), ('size', c_size_t)]

ElinaLincons0ArrayPtr = POINTER(ElinaLincons0Array)
