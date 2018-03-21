# ************************************************************************* #
# elina_tcons0.c: tree expressions
# ************************************************************************* #

from elina_linexpr0_h import *

# ====================================================================== #
# Datatypes
# ====================================================================== #

# IMPORTANT NOTE
# --------------
# correct use of floating-point ELINA_RTYPE_xxx currently supposes that the
# FPU rounds towards +oo


class ElinaTexprOp(CtypesEnum):
    """
    Enum compatible with elina_texpr_op_t from elina_texpr0.h.
    Defines the available operators.
    
    Options
    -------
    ELINA_TEXPR_ADD :   Binary
    ELINA_TEXPR_SUB :   Binary
    ELINA_TEXPR_MUL :   Binary
    ELINA_TEXPR_DIV :   Binary
    ELINA_TEXPR_MOD :   Binary, Either integer or real, no rounding
    ELINA_TEXPR_POW :   Binary
    ELINA_TEXPR_NEG :   Unary, No rounding
    ELINA_TEXPR_CAST :  Unary
    ELINA_TEXPR_SQRT :  Unary
    
    """

    ELINA_TEXPR_ADD = 0
    ELINA_TEXPR_SUB = 1
    ELINA_TEXPR_MUL = 2
    ELINA_TEXPR_DIV = 3
    ELINA_TEXPR_MOD = 4
    ELINA_TEXPR_POW = 5

    ELINA_TEXPR_NEG = 6
    ELINA_TEXPR_CAST = 7
    ELINA_TEXPR_SQRT = 8


class ElinaTexprRtype(CtypesEnum):
    """
    Enum compatible with elina_texpr_rtype_t from elina_texpr0.h.
    Numerical type defining the destination of the rounding.

    Options
    -------
    ELINA_RTYPE_REAL :          Real, no rounding
    ELINA_RTYPE_INT :           Integer
    ELINA_RTYPE_SINGLE :        IEEE 754 32-bit single precision, e.g.: C's float
    ELINA_RTYPE_DOUBLE :        IEEE 754 64-bit double precision, e.g.: C's double
    ELINA_RTYPE_EXTENDED :      Non-standard 80-bit double extended, e.g.: Intel's long double
    ELINA_RTYPE_QUAD :          Non-standard 128-bit quadruple precision, e.g.: Motorola's long double
    ELINA_RTYPE_SIZE :          Not to be used!
    
    """

    ELINA_RTYPE_REAL = 0
    ELINA_RTYPE_INT = 1
    ELINA_RTYPE_SINGLE = 2
    ELINA_RTYPE_DOUBLE = 3
    ELINA_RTYPE_EXTENDED = 4
    ELINA_RTYPE_QUAD = 5
    ELINA_RTYPE_SIZE = 6


class ElinaTexprRdir(CtypesEnum):
    """
    Enum compatible with elina_texpr_rdir_t from elina_texpr0.h.
    Rounding direction.
    
    Options
    -------
    ELINA_RDIR_NEAREST :    Round to nearest with ties to even
    ELINA_RDIR_ZERO :       Round toward zero
    ELINA_RDIR_UP :         Round toward +inf
    ELINA_RDIR_DOWN :       Round toward -inf
    ELINA_RDIR_RND :        All possible modes, non deterministically
    ELINA_RDIR_SIZE :       Not to be used!
    
    """

    ELINA_RDIR_NEAREST = MpfrRnd.MPFR_RNDN.value
    ELINA_RDIR_ZERO = MpfrRnd.MPFR_RNDZ.value
    ELINA_RDIR_UP = MpfrRnd.MPFR_RNDU.value
    ELINA_RDIR_DOWN = MpfrRnd.MPFR_RNDD.value
    ELINA_RDIR_RND = 4
    ELINA_RDIR_SIZE = 5

class ElinaTexpr0Node(Structure):
    """
    ElinaTexpr0Node ctype compatible with elina_texpxr0_node_t from elina_texpr0.h.
    Internal (operator) node.
        
    Fields
    ------
    op : c_uint
        Enum that specifies the operation as defined in ElinaTexprOp.
    type : c_uint
        Enum that specifies the destination type of the rounding.
    dir : c_uint
        Enum that specifies the direction of the rounding.
    exprA : ElinaTexpr0Ptr
        Pointer to the first operand (expression) in the operation.
    exprB : ElinaTexpr0Ptr
        Pointer to the second operand (expression) in the operation.

    """
    pass

ElinaTexpr0NodePtr = POINTER(ElinaTexpr0Node)


class ElinaTexprDiscr(CtypesEnum):
    """
    Enum compatible with elina_texpr_discr_t from elina_texpr0.h.
    Discriminant for the union in ElinaTexpr0 (node types).

    Options
    -------
    ELINA_TEXPR_CST :   ElinaCoeff
    ELINA_TEXPR_DIM :   ElinaDim
    ELINA_TEXPR_NODE :  ElinaTexpr0NodePtr

    """
    ELINA_TEXPR_CST = 0
    ELINA_TEXPR_DIM = 1
    ELINA_TEXPR_NODE = 2


class ElinaTexpr0Union(Structure):
    """
    ElinaTexpr0Union ctype compatible with the union in elina_texpr0_t from elina_texpr0.h.

    Fields
    ------
    cst : ElinaCoeff
        Active in case of leaf node of type ElinaCoeff.
    dim : ElinaDim
        Active in case of leaf node of type ElinaDim.
    node : ElinaTexpr0NodePtr
        Active otherwise.
    
    """

    _fields_ = [('cst', ElinaCoeff), ('dim', ElinaDim), ('node', ElinaTexpr0NodePtr)]


class ElinaTexpr0(Structure):
    """
    ElinaTexpr0 ctype compatible with elina_texpr0_t from elina_texpr0.h.
    
    Fields
    ------
    discr : c_uint
        Discriminant for the union.
    val : ElinaTexpr0Union
        Union containing the core of the expression.
    
    """

    _fields_ = [('discr', c_uint), ('val', ElinaTexpr0Union)]

ElinaTexpr0Ptr = POINTER(ElinaTexpr0)
ElinaTexpr0Array = POINTER(ElinaTexpr0Ptr)

ElinaTexpr0Node._fields_ = [('op', c_uint), ('type', c_uint), ('dir', c_uint),
                ('exprA', ElinaTexpr0Ptr), ('exprB', ElinaTexpr0Ptr)]
