from elina_interval_h import *

# ************************************************************************* #
# elina_coeff.h: coefficients, that are either scalars or intervals
# ************************************************************************* #


class ElinaCoeffDiscr(CtypesEnum):
    """ Enum compatible with discr field in elina_coeff_t from elina_coeff.h """

    ELINA_COEFF_SCALAR = 0
    ELINA_COEFF_INTERVAL = 1


class ElinaCoeffUnion(Union):
    """ Ctype representing the union field in elina_coeff_t from elina_coeff.h """

    _fields_ = [('scalar', ElinaScalarPtr), ('interval', ElinaIntervalPtr)]


class ElinaCoeff(Structure):
    """ ElinaCoeff ctype compatible with elina_coeff_t from elina_coeff.h """

    _fields_ = [('discr', c_uint), ('val', ElinaCoeffUnion)]


ElinaCoeffPtr = POINTER(ElinaCoeff)
