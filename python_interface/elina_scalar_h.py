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
