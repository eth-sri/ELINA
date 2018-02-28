from ctypes import *
from enum import IntEnum

elina_auxiliary_api = CDLL("libelinaux.so")


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
