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


from ctypes import *
from enum import IntEnum
import os
import ctypes

# elina_auxiliary_api = CDLL("libelinaux.so")
elina_auxiliary_api = ctypes.CDLL(os.path.join(os.path.dirname(os.path.abspath(__file__)),"../elina_auxiliary/libelinaux.so"))


class CtypesEnum(IntEnum):
    """ A ctypes compatible IntEnum superclass """

    @classmethod
    def from_param(cls, obj):
        return int(obj)


class Mpz(Structure):
    """ Mpz ctype compatible with mpz_t from gmp.h """

    _fields_ = [('_mp_alloc', c_int), ('_mp_size', c_int), ('_mp_d', POINTER(c_long))]

Mpz_t = Mpz * 1


class Mpq(Structure):
    """ Mpq ctype compatible with mpq_t from gmp.h """

    _fields_ = [('_mp_num', Mpz), ('_mp_den', Mpz)]

Mpq_t = Mpq * 1


class Mpfr(Structure):
    """ Mpfr ctype compatible with mpfr_t from mpfr.h """

    _fields_ = [('_mpfr_prec', c_long), ('_mpfr_sign', c_int), ('_mpfr_exp', c_longlong), ('_mpfr_d', POINTER(c_ulonglong))]

Mpfr_t = Mpfr * 1


class MpfrRnd(CtypesEnum):
    """ Rounding enums compatible with roundings from mpfr.h """

    MPFR_RNDN = 0  # round to nearest with ties to even
    MPFR_RNDZ = 1  # round toward zero
    MPFR_RNDU = 2  # round toward +Inf
    MPFR_RNDD = 3  # round toward -Inf
    MPFR_RNDA = 4  # round away from zero
    MPFR_RNDF = 5  # faithful rounding ( not implemented yet)
    MPFR_RNDNA = -1  # round to nearest, with ties away from zero (mpfr_round)
