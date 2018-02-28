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