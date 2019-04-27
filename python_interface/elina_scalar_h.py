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


from elina_auxiliary_imports import *


# ********************************************************************** #
# I. Datatypes
# ********************************************************************** #


class ElinaScalarDiscr(CtypesEnum):
    """ Enum compatible with discr field in elina_scalar_t from elina_scalar.h """

    ELINA_SCALAR_DOUBLE = 0
    ELINA_SCALAR_MPQ = 1
    ELINA_SCALAR_MPFR = 2


class ElinaScalarUnion(Union):
    """ Ctype representing the union field in elina_scalar_t from elina_scalar.h """

    _fields_ = [('dbl', c_double), ('mpq_ptr', POINTER(Mpq)), ('mpfr_ptr', POINTER(Mpfr))]


class ElinaScalar(Structure):
    """ ElinaScalar ctype compatible with elina_scalar_t from elina_scalar.h """

    _fields_ = [('discr', c_uint), ('val', ElinaScalarUnion)]


ElinaScalarPtr = POINTER(ElinaScalar)

elina_scalar_print_prec = c_int(20)
