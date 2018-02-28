from elina_scalar_h import *

# ************************************************************************* #
# elina_interval.h: intervals
# ************************************************************************* #


class ElinaInterval(Structure):
    """ ElinaInterval ctype compatible with elina_interval_t from elina_interval.h """

    _fields_ = [('inf', ElinaScalarPtr), ('sup', ElinaScalarPtr)]

ElinaIntervalPtr = POINTER(ElinaInterval)
ElinaIntervalArray = POINTER(ElinaIntervalPtr)
