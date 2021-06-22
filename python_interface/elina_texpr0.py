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


from elina_auxiliary_imports import *
from elina_texpr0_h import *


# ======================================================================
# I. Constructors and Destructors
# ======================================================================

def elina_texpr0_cst(coeff):
    """
    Create a constant leaf ElinaTexpr0 by using an ElinaCoeff.
    
    Parameters
    ----------
    coeff : ElinaCoeffPtr
        Pointer to the ElinaCoeff used for creation.

    Returns
    -------
    texpr : ElinaTexpr0Ptr
        Pointer to the newly created ElinaTexpr0.

    """

    texpr = None
    try:
        elina_texpr0_cst_c = elina_auxiliary_api.elina_texpr0_cst
        elina_texpr0_cst_c.restype = ElinaTexpr0Ptr
        elina_texpr0_cst_c.argtypes = [ElinaCoeffPtr]
        texpr = elina_texpr0_cst_c(coeff)
    except:
        print('Problem with loading/calling "elina_texpr0_cst" from "libelinaux.so"')
        print('Make sure you are passing ElinaCoeffPtr to the function')

    return texpr


def elina_texpr0_cst_scalar(scalar):
    """
    Create a constant leaf ElinaTexpr0 by using an ElinaScalar.

    Parameters
    ----------
    scalar : ElinaScalarPtr
        Pointer to the ElinaScalar used for creation.

    Returns
    -------
    texpr : ElinaTexpr0Ptr
        Pointer to the newly created ElinaTexpr0.

    """

    texpr = None
    try:
        elina_texpr0_cst_scalar_c = elina_auxiliary_api.elina_texpr0_cst_scalar
        elina_texpr0_cst_scalar_c.restype = ElinaTexpr0Ptr
        elina_texpr0_cst_scalar_c.argtypes = [ElinaScalarPtr]
        texpr = elina_texpr0_cst_scalar_c(scalar)
    except:
        print('Problem with loading/calling "elina_texpr0_cst_scalar" from "libelinaux.so"')
        print('Make sure you are passing ElinaScalarPtr to the function')

    return texpr


def elina_texpr0_cst_scalar_mpq(mpq_t):
    """
    Create a constant leaf ElinaTexpr0 by using a scalar defined by a Mpq_t.

    Parameters
    ----------
    mpq_t : Mpq_t
        Mpq_t used for creation.

    Returns
    -------
    texpr : ElinaTexpr0Ptr
        Pointer to the newly created ElinaTexpr0.

    """

    texpr = None
    try:
        elina_texpr0_cst_scalar_mpq_c = elina_auxiliary_api.elina_texpr0_cst_scalar_mpq
        elina_texpr0_cst_scalar_mpq_c.restype = ElinaTexpr0Ptr
        elina_texpr0_cst_scalar_mpq_c.argtypes = [Mpq_t]
        texpr = elina_texpr0_cst_scalar_mpq_c(mpq_t)
    except:
        print('Problem with loading/calling "elina_texpr0_cst_scalar_mpq" from "libelinaux.so"')
        print('Make sure you are passing Mpq_t to the function')

    return texpr


def elina_texpr0_cst_scalar_mpfr(mpfr_t):
    """
    Create a constant leaf ElinaTexpr0 by using a scalar defined by a Mpfr_t.

    Parameters
    ----------
    mpfr_t : Mpfr_t
        Mpfr_t used for creation.

    Returns
    -------
    texpr : ElinaTexpr0Ptr
        Pointer to the newly created ElinaTexpr0.

    """

    texpr = None
    try:
        elina_texpr0_cst_scalar_mpfr_c = elina_auxiliary_api.elina_texpr0_cst_scalar_mpfr
        elina_texpr0_cst_scalar_mpfr_c.restype = ElinaTexpr0Ptr
        elina_texpr0_cst_scalar_mpfr_c.argtypes = [Mpfr_t]
        texpr = elina_texpr0_cst_scalar_mpfr_c(mpfr_t)
    except:
        print('Problem with loading/calling "elina_texpr0_cst_scalar_mpfr" from "libelinaux.so"')
        print('Make sure you are passing Mpfr_t to the function')

    return texpr


def elina_texpr0_cst_scalar_int(num):
    """
    Create a constant leaf ElinaTexpr0 by using a scalar defined by a c_long.

    Parameters
    ----------
    num : c_long
        c_long used for creation.

    Returns
    -------
    texpr : ElinaTexpr0Ptr
        Pointer to the newly created ElinaTexpr0.

    """

    texpr = None
    try:
        elina_texpr0_cst_scalar_int_c = elina_auxiliary_api.elina_texpr0_cst_scalar_int
        elina_texpr0_cst_scalar_int_c.restype = ElinaTexpr0Ptr
        elina_texpr0_cst_scalar_int_c.argtypes = [c_long]
        texpr = elina_texpr0_cst_scalar_int_c(num)
    except:
        print('Problem with loading/calling "elina_texpr0_cst_scalar_int" from "libelinaux.so"')
        print('Make sure you are passing c_long to the function')

    return texpr


def elina_texpr0_cst_scalar_frac(num, den):
    """
    Create a constant leaf ElinaTexpr0 by using a scalar defined by an integer fraction.

    Parameters
    ----------
    num : c_long
        c_long used as numerator for creation.
    den : c_ulong
        c_ulong used as denominator for creation.

    Returns
    -------
    texpr : ElinaTexpr0Ptr
        Pointer to the newly created ElinaTexpr0.

    """

    texpr = None
    try:
        elina_texpr0_cst_scalar_frac_c = elina_auxiliary_api.elina_texpr0_cst_scalar_frac
        elina_texpr0_cst_scalar_frac_c.restype = ElinaTexpr0Ptr
        elina_texpr0_cst_scalar_frac_c.argtypes = [c_long, c_ulong]
        texpr = elina_texpr0_cst_scalar_frac_c(num, den)
    except:
        print('Problem with loading/calling "elina_texpr0_cst_scalar_frac" from "libelinaux.so"')
        print('Make sure you are passing c_long, c_ulong to the function')

    return texpr


def elina_texpr0_cst_scalar_double(num):
    """
    Create a constant leaf ElinaTexpr0 by using a scalar defined by a c_double.

    Parameters
    ----------
    num : c_double
        c_double used for creation.

    Returns
    -------
    texpr : ElinaTexpr0Ptr
        Pointer to the newly created ElinaTexpr0.

    """

    texpr = None
    try:
        elina_texpr0_cst_scalar_double_c = elina_auxiliary_api.elina_texpr0_cst_scalar_double
        elina_texpr0_cst_scalar_double_c.restype = ElinaTexpr0Ptr
        elina_texpr0_cst_scalar_double_c.argtypes = [c_double]
        texpr = elina_texpr0_cst_scalar_double_c(num)
    except:
        print('Problem with loading/calling "elina_texpr0_cst_scalar_double" from "libelinaux.so"')
        print('Make sure you are passing c_double to the function')

    return texpr


def elina_texpr0_cst_interval(interval):
    """
    Create a constant leaf ElinaTexpr0 by using an ElinaInterval.

    Parameters
    ----------
    interval : ElinaIntervalPtr
        Pointer to the ELinaInterval used for creation.

    Returns
    -------
    texpr : ElinaTexpr0Ptr
        Pointer to the newly created ElinaTexpr0.

    """

    texpr = None
    try:
        elina_texpr0_cst_interval_c = elina_auxiliary_api.elina_texpr0_cst_interval
        elina_texpr0_cst_interval_c.restype = ElinaTexpr0Ptr
        elina_texpr0_cst_interval_c.argtypes = [ElinaIntervalPtr]
        texpr = elina_texpr0_cst_interval_c(interval)
    except:
        print('Problem with loading/calling "elina_texpr0_cst_interval" from "libelinaux.so"')
        print('Make sure you are passing ElinaIntervalPtr to the function')

    return texpr


def elina_texpr0_cst_interval_scalar(inf, sup):
    """
    Create a constant leaf ElinaTexpr0 by using an interval defined by two ElinaScalar-s.

    Parameters
    ----------
    inf : ElinaScalarPtr
        Pointer to the ElinaScalar used as a lower bound for creation.
    sup : ElinaScalarPtr
        Pointer to the ElinaScalar used as an upper bound for creation.


    Returns
    -------
    texpr : ElinaTexpr0Ptr
        Pointer to the newly created ElinaTexpr0.

    """

    texpr = None
    try:
        elina_texpr0_cst_interval_scalar_c = elina_auxiliary_api.elina_texpr0_cst_interval_scalar
        elina_texpr0_cst_interval_scalar_c.restype = ElinaTexpr0Ptr
        elina_texpr0_cst_interval_scalar_c.argtypes = [ElinaScalarPtr, ElinaScalarPtr]
        texpr = elina_texpr0_cst_interval_scalar_c(inf, sup)
    except:
        print('Problem with loading/calling "elina_texpr0_cst_interval_scalar" from "libelinaux.so"')
        print('Make sure you are passing ElinaScalarPtr, ElinaScalarPtr to the function')

    return texpr


def elina_texpr0_cst_interval_mpq(inf, sup):
    """
    Create a constant leaf ElinaTexpr0 by using an interval defined by two Mpq_t-s.

    Parameters
    ----------
    inf : Mpq_t
        Mpq_t used as a lower bound for creation.
    sup : Mpq_t
        Mpq_t used as an upper bound for creation.


    Returns
    -------
    texpr : ElinaTexpr0Ptr
        Pointer to the newly created ElinaTexpr0.

    """

    texpr = None
    try:
        elina_texpr0_cst_interval_mpq_c = elina_auxiliary_api.elina_texpr0_cst_interval_mpq
        elina_texpr0_cst_interval_mpq_c.restype = ElinaTexpr0Ptr
        elina_texpr0_cst_interval_mpq_c.argtypes = [Mpq_t, Mpq_t]
        texpr = elina_texpr0_cst_interval_mpq_c(inf, sup)
    except:
        print('Problem with loading/calling "elina_texpr0_cst_interval_mpq" from "libelinaux.so"')
        print('Make sure you are passing Mpq_t, Mpq_t to the function')

    return texpr


def elina_texpr0_cst_interval_mpfr(inf, sup):
    """
    Create a constant leaf ElinaTexpr0 by using an interval defined by two Mpfr_t-s.

    Parameters
    ----------
    inf : Mpfr_t
        Mpfr_t used as a lower bound for creation.
    sup : Mpfr_t
        Mpfr_t used as an upper bound for creation.


    Returns
    -------
    texpr : ElinaTexpr0Ptr
        Pointer to the newly created ElinaTexpr0.

    """

    texpr = None
    try:
        elina_texpr0_cst_interval_mpfr_c = elina_auxiliary_api.elina_texpr0_cst_interval_mpfr
        elina_texpr0_cst_interval_mpfr_c.restype = ElinaTexpr0Ptr
        elina_texpr0_cst_interval_mpfr_c.argtypes = [Mpfr_t, Mpfr_t]
        texpr = elina_texpr0_cst_interval_mpfr_c(inf, sup)
    except:
        print('Problem with loading/calling "elina_texpr0_cst_interval_mpfr" from "libelinaux.so"')
        print('Make sure you are passing Mpfr_t, Mpfr_t to the function')

    return texpr


def elina_texpr0_cst_interval_int(inf, sup):
    """
    Create a constant leaf ElinaTexpr0 by using an interval defined by two integers.

    Parameters
    ----------
    inf : c_long
        c_long used as a lower bound for creation.
    sup : c_long
        c_long used as an upper bound for creation.


    Returns
    -------
    texpr : ElinaTexpr0Ptr
        Pointer to the newly created ElinaTexpr0.

    """

    texpr = None
    try:
        elina_texpr0_cst_interval_int_c = elina_auxiliary_api.elina_texpr0_cst_interval_int
        elina_texpr0_cst_interval_int_c.restype = ElinaTexpr0Ptr
        elina_texpr0_cst_interval_int_c.argtypes = [c_long, c_long]
        texpr = elina_texpr0_cst_interval_int_c(inf, sup)
    except:
        print('Problem with loading/calling "elina_texpr0_cst_interval_int" from "libelinaux.so"')
        print('Make sure you are passing c_long, c_long to the function')

    return texpr


def elina_texpr0_cst_interval_frac(numinf, deninf, numsup, densup):
    """
    Create a constant leaf ElinaTexpr0 by using an interval defined by two integer fractions.

    Parameters
    ----------
    numinf : c_long
        c_long used as a numerator of lower bound for creation.
    deninf : c_ulong
        c_ulong used as a denominator of lower bound for creation.
    numsup : c_long
        c_long used as a numerator of upper bound for creation.
    densup : c_ulong
        c_ulong used as a denominator of upper bound for creation.


    Returns
    -------
    texpr : ElinaTexpr0Ptr
        Pointer to the newly created ElinaTexpr0.

    """

    texpr = None
    try:
        elina_texpr0_cst_interval_frac_c = elina_auxiliary_api.elina_texpr0_cst_interval_frac
        elina_texpr0_cst_interval_frac_c.restype = ElinaTexpr0Ptr
        elina_texpr0_cst_interval_frac_c.argtypes = [c_long, c_ulong, c_long, c_ulong]
        texpr = elina_texpr0_cst_interval_frac_c(numinf, deninf, numsup, densup)
    except:
        print('Problem with loading/calling "elina_texpr0_cst_interval_frac" from "libelinaux.so"')
        print('Make sure you are passing c_long, c_ulong, c_long, c_ulong to the function')

    return texpr


def elina_texpr0_cst_interval_double(inf, sup):
    """
    Create a constant leaf ElinaTexpr0 by using an interval defined by two doubles.

    Parameters
    ----------
    inf : c_double
        c_double used as a lower bound for creation.
    sup : c_double
        c_double used as an upper bound for creation.


    Returns
    -------
    texpr : ElinaTexpr0Ptr
        Pointer to the newly created ElinaTexpr0.

    """

    texpr = None
    try:
        elina_texpr0_cst_interval_double_c = elina_auxiliary_api.elina_texpr0_cst_interval_double
        elina_texpr0_cst_interval_double_c.restype = ElinaTexpr0Ptr
        elina_texpr0_cst_interval_double_c.argtypes = [c_double, c_double]
        texpr = elina_texpr0_cst_interval_double_c(inf, sup)
    except:
        print('Problem with loading/calling "elina_texpr0_cst_interval_double" from "libelinaux.so"')
        print('Make sure you are passing c_double, c_double to the function')

    return texpr


def elina_texpr0_cst_interval_top():
    """
    Create a constant leaf ElinaTexpr0 by using the universe interval [-oo, +oo].

    Returns
    -------
    texpr : ElinaTexpr0Ptr
        Pointer to the newly created ElinaTexpr0.

    """

    texpr = None
    try:
        elina_texpr0_cst_interval_top_c = elina_auxiliary_api.elina_texpr0_cst_interval_top
        elina_texpr0_cst_interval_top_c.restype = ElinaTexpr0Ptr
        texpr = elina_texpr0_cst_interval_top_c()
    except:
        print('Problem with loading/calling "elina_texpr0_cst_interval_top" from "libelinaux.so"')

    return texpr


def elina_texpr0_dim(dim):
    """
    Create a dimension leaf ElinaTexpr0 by using an ElinaDim.

    Parameters
    ----------
    dim : ElinaDim
        ElinaDim used for creation.

    Returns
    -------
    texpr : ElinaTexpr0Ptr
        Pointer to the newly created ElinaTexpr0.

    """

    texpr = None
    try:
        elina_texpr0_dim_c = elina_auxiliary_api.elina_texpr0_dim
        elina_texpr0_dim_c.restype = ElinaTexpr0Ptr
        elina_texpr0_dim_c.argtypes = [ElinaDim]
        texpr = elina_texpr0_dim_c(dim)
    except:
        print('Problem with loading/calling "elina_texpr0_dim" from "libelinaux.so"')
        print('Make sure you are passing ElinaDim to the function')

    return texpr


def elina_texpr0_unop(op, operand, type, dir):
    """
    Create an unary operator ElinaTexpr0.

    Parameters
    ----------
    op : c_uint
        Enum defining the operation as defined in ElinaTexprOp. 
    operand : ElinaTexpr0Ptr
        Operand on which to perform the operation.
    type : c_uint
        Enum defining the resulting type of the rounding as defined by ElinaTexprRtype.
    dir : c_uint
        Enum defining the direction of the rounding as defined by ElinaTexprRdir.
        
    Returns
    -------
    texpr : ElinaTexpr0Ptr
        Pointer to the newly created ElinaTexpr0.

    """

    texpr = None
    try:
        elina_texpr0_unop_c = elina_auxiliary_api.elina_texpr0_unop
        elina_texpr0_unop_c.restype = ElinaTexpr0Ptr
        elina_texpr0_unop_c.argtypes = [c_uint, ElinaTexpr0Ptr, c_uint, c_uint]
        texpr = elina_texpr0_unop_c(op, operand, type, dir)
    except:
        print('Problem with loading/calling "elina_texpr0_unop" from "libelinaux.so"')
        print('Make sure you are passing c_uint, ElinaTexpr0Ptr, c_uint, c_uint to the function')

    return texpr


def elina_texpr0_binop(op, operand_a, operand_b, type, dir):
    """
    Create a binary operator ElinaTexpr0.

    Parameters
    ----------
    op : c_uint
        Enum defining the operation as defined in ElinaTexprOp. 
    operand_a : ElinaTexpr0Ptr
        First operand on which to perform the operation.
    operand_b : ElinaTexpr0Ptr
        Second operand on which to perform the operation.
    type : c_uint
        Enum defining the resulting type of the rounding as defined by ElinaTexprRtype.
    dir : c_uint
        Enum defining the direction of the rounding as defined by ElinaTexprRdir.

    Returns
    -------
    texpr : ElinaTexpr0Ptr
        Pointer to the newly created ElinaTexpr0.

    """

    texpr = None
    try:
        elina_texpr0_binop_c = elina_auxiliary_api.elina_texpr0_binop
        elina_texpr0_binop_c.restype = ElinaTexpr0Ptr
        elina_texpr0_binop_c.argtypes = [c_uint, ElinaTexpr0Ptr, ElinaTexpr0Ptr, c_uint, c_uint]
        texpr = elina_texpr0_binop_c(op, operand_a, operand_b, type, dir)
    except:
        print('Problem with loading/calling "elina_texpr0_binop" from "libelinaux.so"')
        print('Make sure you are passing c_uint, ElinaTexpr0Ptr, ElinaTexpr0Ptr, c_uint, c_uint to the function')

    return texpr


def elina_texpr0_copy(texpr2):
    """
    Copy an ElinaTexpr0.

    Parameters
    ----------
    texpr2 : ElinaTexpr0Ptr
        Pointer to the ElinaTexpr0 that needs to be copied.
    
    Returns
    -------
    texpr1 : ElinaTexpr0Ptr
        Pointer to the newly created copy ElinaTexpr0.

    """


    texpr1 = None
    try:
        elina_texpr0_copy_c = elina_auxiliary_api.elina_texpr0_copy
        elina_texpr0_copy_c.restype = ElinaTexpr0Ptr
        elina_texpr0_copy_c.argtypes = [ElinaTexpr0Ptr]
        texpr1 = elina_texpr0_copy_c(texpr2)
    except:
        print('Problem with loading/calling "elina_texpr0_copy" from "libelinaux.so"')
        print('Make sure you are passing ElinaTexpr0Ptr to the function')

    return texpr1


def elina_texpr0_free(texpr):
    """
    Free an ElinaTexpr0.

    Parameters
    ----------
    texpr : ElinaTexpr0Ptr
        Pointer to the ElinaTexpr0 that needs to be freed.
    
    Returns
    -------
    None
    
    """

    try:
        elina_texpr0_free_c = elina_auxiliary_api.elina_texpr0_free
        elina_texpr0_free_c.restype = None
        elina_texpr0_free_c.argtypes = [ElinaTexpr0Ptr]
        elina_texpr0_free_c(texpr)
    except:
        print('Problem with loading/calling "elina_texpr0_free" from "libelinaux.so"')
        print('Make sure you are passing ElinaTexpr0Ptr to the function')


def elina_texpr0_from_linexpr0(linexpr):
    """
    Create a comb-like ElinaTexpr0 by using an ElinaLinexpr0.

    Parameters
    ----------
    linexpr : ELinaLinexpr0Ptr
        Pointer to the ELinaLinexpr0 that needs to be transformed to ElinaTexpr0.

    Returns
    -------
    texpr : ElinaTexpr0Ptr
        Pointer to the newly created ElinaTexpr0.

    """

    texpr = None
    try:
        elina_texpr0_from_linexpr0_c = elina_auxiliary_api.elina_texpr0_from_linexpr0
        elina_texpr0_from_linexpr0_c.restype = ElinaTexpr0Ptr
        elina_texpr0_from_linexpr0_c.argtypes = [ElinaLinexpr0Ptr]
        texpr = elina_texpr0_from_linexpr0_c(linexpr)
    except:
        print('Problem with loading/calling "elina_texpr0_from_linexpr0" from "libelinaux.so"')
        print('Make sure you are passing ElinaLinexpr0Ptr to the function')

    return texpr


# ====================================================================== #
# II. Printing
# ====================================================================== #

def elina_texpr0_fprint(stream, texpr, name_of_dim):
    """
    Print an ElinaTexpr0 to stream by using dimension names.

    Parameters
    ----------
    stream : c_void_p
        Stream onto which to print.
    texpr : ElinaTexpr0Ptr
        Pointer to the ElinaTexpr0 that needs to be printed.
    name_of_dim : POINTER(c_char_p)
        Dimension names.
        
    Returns
    -------
    None
    
    """

    try:
        elina_texpr0_fprint_c = elina_auxiliary_api.elina_texpr0_fprint
        elina_texpr0_fprint_c.restype = None
        elina_texpr0_fprint_c.argtypes = [c_void_p, ElinaTexpr0Ptr, POINTER(c_char_p)]
        elina_texpr0_fprint_c(stream, texpr, name_of_dim)
    except:
        print('Problem with loading/calling "elina_texpr0_fprint" from "libelinaux.so"')
        print('Make sure you are passing c_void_p, ElinaTexpr0Ptr, POINTER(c_char_p) to the function')


def elina_texpr0_print(texpr, name_of_dim):
    """
    Print an ElinaTexpr0 to stdout by using dimension names.

    Parameters
    ----------
    texpr : ElinaTexpr0Ptr
        Pointer to the ElinaTexpr0 that needs to be printed.
    name_of_dim : POINTER(c_char_p)
        Dimension names.

    Returns
    -------
    None
    
    """

    try:
        elina_texpr0_print_c = elina_auxiliary_api.elina_texpr0_print
        elina_texpr0_print_c.restype = None
        elina_texpr0_print_c.argtypes = [ElinaTexpr0Ptr, POINTER(c_char_p)]
        elina_texpr0_print_c(texpr, name_of_dim)
    except:
        print('Problem with loading/calling "elina_texpr0_print" from "libelinaux.so"')
        print('Make sure you are passing ElinaTexpr0Ptr, POINTER(c_char_p) to the function')


# ====================================================================== #
# III. Tests, size
# ====================================================================== #

def elina_texpr0_depth(texpr):
    """
    Return the depth of an ElinaTexpr0, in terms of operator nodes.

    Parameters
    ----------
    texpr : ElinaTexpr0Ptr
        Pointer to the ElinaTexpr0 for which we want to know the depth.

    Returns
    -------
    depth : c_size_t
        Depth of the ElinaTexpr0.

    """

    depth = None
    try:
        elina_texpr0_depth_c = elina_auxiliary_api.elina_texpr0_depth
        elina_texpr0_depth_c.restype = c_size_t
        elina_texpr0_depth_c.argtypes = [ElinaTexpr0Ptr]
        depth = elina_texpr0_depth_c(texpr)
    except:
        print('Problem with loading/calling "elina_texpr0_depth" from "libelinaux.so"')
        print('Make sure you are passing ElinaTexpr0Ptr to the function')

    return depth


def elina_texpr0_size(texpr):
    """
    Return the number of operator nodes in an ElinaTexpr0.

    Parameters
    ----------
    texpr : ElinaTexpr0Ptr
        Pointer to the ElinaTexpr0 for which we want to know the size.

    Returns
    -------
    size : c_size_t
        Size of the ElinaTexpr0.

    """

    size = None
    try:
        elina_texpr0_size_c = elina_auxiliary_api.elina_texpr0_size
        elina_texpr0_size_c.restype = c_size_t
        elina_texpr0_size_c.argtypes = [ElinaTexpr0Ptr]
        size = elina_texpr0_size_c(texpr)
    except:
        print('Problem with loading/calling "elina_texpr0_size" from "libelinaux.so"')
        print('Make sure you are passing ElinaTexpr0Ptr to the function')

    return size


def elina_texpr0_max_dim(texpr):
    """
    Return the maximum ElinaDim + 1 of all dimensions in the expression.
    Return 0 if there are no dimensions at all. 

    Parameters
    ----------
    texpr : ElinaTexpr0Ptr
        Pointer to the ElinaTexpr0 for which we want to know the maximum dimension.

    Returns
    -------
    max_dim : ElinaDim
        Maximum ElinaDim + 1 of the given ElinaTexpr0.

    """

    max_dim = None
    try:
        elina_texpr0_max_dim_c = elina_auxiliary_api.elina_texpr0_max_dim
        elina_texpr0_max_dim_c.restype = ElinaDim
        elina_texpr0_max_dim_c.argtypes = [ElinaTexpr0Ptr]
        max_dim = elina_texpr0_max_dim_c(texpr)
    except:
        print('Problem with loading/calling "elina_texpr0_max_dim" from "libelinaux.so"')
        print('Make sure you are passing ElinaTexpr0Ptr to the function')

    return max_dim


def elina_texpr0_has_dim(texpr, dim):
    """
    Check if an ElinaTexpr0 contains a certain ElinaDim. 

    Parameters
    ----------
    texpr : ElinaTexpr0Ptr
        Pointer to the ElinaTexpr0 which we want to check.
    dim : ElinaDim
        ElinaDim used for the check.

    Returns
    -------
    result : c_bool
        Result of the check.

    """

    result = None
    try:
        elina_texpr0_has_dim_c = elina_auxiliary_api.elina_texpr0_has_dim
        elina_texpr0_has_dim_c.restype = c_bool
        elina_texpr0_has_dim_c.argtypes = [ElinaTexpr0Ptr, ElinaDim]
        result = elina_texpr0_has_dim_c(texpr, dim)
    except:
        print('Problem with loading/calling "elina_texpr0_has_dim" from "libelinaux.so"')
        print('Make sure you are passing ElinaTexpr0Ptr to the function')

    return result


def elina_texpr0_dimlist(texpr):
    """
    Return an ordered, ELINA_DIM_MAX-terminated array of occurring dimensions in a given ElinaTexpr0.
    
    Parameters
    ----------
    texpr : ElinaTexpr0Ptr
        Pointer to the ElinaTexpr0.

    Returns
    -------
    dim_array : ElinaDimPtr
        Array of ordered dimensions for the ElinaTexpr0.

    """

    dim_array = None
    try:
        elina_texpr0_dimlist_c = elina_auxiliary_api.elina_texpr0_dimlist
        elina_texpr0_dimlist_c.restype = ElinaDimPtr
        elina_texpr0_dimlist_c.argtypes = [ElinaTexpr0Ptr]
        dim_array = elina_texpr0_dimlist_c(texpr)
    except:
        print('Problem with loading/calling "elina_texpr0_dimlist" from "libelinaux.so"')
        print('Make sure you are passing ElinaTexpr0Ptr to the function')

    return dim_array


def elina_texpr0_is_interval_cst(texpr):
    """
    Test if an ElinaTexpr0 is consisted of only constant leaves.
    
    Parameters
    ----------
    texpr : ElinaTexpr0Ptr
        Pointer to the ElinaTexpr0 that needs to be tested.

    Returns
    -------
    result : c_bool
        Result of the test.

    """

    result = None
    try:
        elina_texpr0_is_interval_cst_c = elina_auxiliary_api.elina_texpr0_is_interval_cst
        elina_texpr0_is_interval_cst_c.restype = c_bool
        elina_texpr0_is_interval_cst_c.argtypes = [ElinaTexpr0Ptr]
        result = elina_texpr0_is_interval_cst_c(texpr)
    except:
        print('Problem with loading/calling "elina_texpr0_is_interval_cst" from "libelinaux.so"')
        print('Make sure you are passing ElinaTexpr0Ptr to the function')

    return result


def elina_texpr0_is_interval_linear(texpr):
    """
    Test if an ElinaTexpr0 is linear with possibly interval coefficients, but no rounding.

    Parameters
    ----------
    texpr : ElinaTexpr0Ptr
        Pointer to the ElinaTexpr0 that needs to be tested.

    Returns
    -------
    result : c_bool
        Result of the test.

    """

    result = None
    try:
        elina_texpr0_is_interval_linear_c = elina_auxiliary_api.elina_texpr0_is_interval_linear
        elina_texpr0_is_interval_linear_c.restype = c_bool
        elina_texpr0_is_interval_linear_c.argtypes = [ElinaTexpr0Ptr]
        result = elina_texpr0_is_interval_linear_c(texpr)
    except:
        print('Problem with loading/calling "elina_texpr0_is_interval_linear" from "libelinaux.so"')
        print('Make sure you are passing ElinaTexpr0Ptr to the function')

    return result


def elina_texpr0_is_interval_polynomial(texpr):
    """
    Test if an ElinaTexpr0 is polynomial with possibly interval coefficients, but no rounding.

    Parameters
    ----------
    texpr : ElinaTexpr0Ptr
        Pointer to the ElinaTexpr0 that needs to be tested.

    Returns
    -------
    result : c_bool
        Result of the test.

    """

    result = None
    try:
        elina_texpr0_is_interval_polynomial_c = elina_auxiliary_api.elina_texpr0_is_interval_polynomial
        elina_texpr0_is_interval_polynomial_c.restype = c_bool
        elina_texpr0_is_interval_polynomial_c.argtypes = [ElinaTexpr0Ptr]
        result = elina_texpr0_is_interval_polynomial_c(texpr)
    except:
        print('Problem with loading/calling "elina_texpr0_is_interval_polynomial" from "libelinaux.so"')
        print('Make sure you are passing ElinaTexpr0Ptr to the function')

    return result


def elina_texpr0_is_interval_polyfrac(texpr):
    """
    Test if an ElinaTexpr0 is polynomial fraction with possibly interval coefficients, but no rounding.

    Parameters
    ----------
    texpr : ElinaTexpr0Ptr
        Pointer to the ElinaTexpr0 that needs to be tested.

    Returns
    -------
    result : c_bool
        Result of the test.

    """

    result = None
    try:
        elina_texpr0_is_interval_polyfrac_c = elina_auxiliary_api.elina_texpr0_is_interval_polyfrac
        elina_texpr0_is_interval_polyfrac_c.restype = c_bool
        elina_texpr0_is_interval_polyfrac_c.argtypes = [ElinaTexpr0Ptr]
        result = elina_texpr0_is_interval_polyfrac_c(texpr)
    except:
        print('Problem with loading/calling "elina_texpr0_is_interval_polyfrac" from "libelinaux.so"')
        print('Make sure you are passing ElinaTexpr0Ptr to the function')

    return result


def elina_texpr0_is_scalar(texpr):
    """
    Test if an ElinaTexpr0 has only scalar coefficients.

    Parameters
    ----------
    texpr : ElinaTexpr0Ptr
        Pointer to the ElinaTexpr0 that needs to be tested.

    Returns
    -------
    result : c_bool
        Result of the test.

    """

    result = None
    try:
        elina_texpr0_is_scalar_c = elina_auxiliary_api.elina_texpr0_is_scalar
        elina_texpr0_is_scalar_c.restype = c_bool
        elina_texpr0_is_scalar_c.argtypes = [ElinaTexpr0Ptr]
        result = elina_texpr0_is_scalar_c(texpr)
    except:
        print('Problem with loading/calling "elina_texpr0_is_scalar" from "libelinaux.so"')
        print('Make sure you are passing ElinaTexpr0Ptr to the function')

    return result


def elina_texpr0_array_is_interval_linear(texpr_array, size):
    """
    Test if an ElinaTexpr0Array is linear with possibly interval coefficients, but no rounding.

    Parameters
    ----------
    texpr_array : ElinaTexpr0Array
        ElinaTexpr0Array that needs to be tested.
    size : c_size_t
        Size of the ElinaTexpr0Array.

    Returns
    -------
    result : c_bool
        Result of the test.

    """


    result = None
    try:
        elina_texpr0_array_is_interval_linear_c = elina_auxiliary_api.elina_texpr0_array_is_interval_linear
        elina_texpr0_array_is_interval_linear_c.restype = c_bool
        elina_texpr0_array_is_interval_linear_c.argtypes = [ElinaTexpr0Array, c_size_t]
        result = elina_texpr0_array_is_interval_linear_c(texpr_array, size)
    except:
        print('Problem with loading/calling "elina_texpr0_array_is_interval_linear" from "libelinaux.so"')
        print('Make sure you are passing ElinaTexpr0Array, c_size_t to the function')

    return result


def elina_texpr0_array_is_interval_polynomial(texpr_array, size):
    """
    Test if an ElinaTexpr0Array is polynomial with possibly interval coefficients, but no rounding.

    Parameters
    ----------
    texpr_array : ElinaTexpr0Array
        ElinaTexpr0Array that needs to be tested.
    size : c_size_t
        Size of the ElinaTexpr0Array.

    Returns
    -------
    result : c_bool
        Result of the test.

    """

    result = None
    try:
        elina_texpr0_array_is_interval_polynomial_c = elina_auxiliary_api.elina_texpr0_array_is_interval_polynomial
        elina_texpr0_array_is_interval_polynomial_c.restype = c_bool
        elina_texpr0_array_is_interval_polynomial_c.argtypes = [ElinaTexpr0Array, c_size_t]
        result = elina_texpr0_array_is_interval_polynomial_c(texpr_array, size)
    except:
        print('Problem with loading/calling "elina_texpr0_array_is_interval_polynomial" from "libelinaux.so"')
        print('Make sure you are passing ElinaTexpr0Array, c_size_t to the function')

    return result


def elina_texpr0_array_is_interval_polyfrac(texpr_array, size):
    """
    Test if an ElinaTexpr0Array is polynomial fraction with possibly interval coefficients, but no rounding.

    Parameters
    ----------
    texpr_array : ElinaTexpr0Array
        ElinaTexpr0Array that needs to be tested.
    size : c_size_t
        Size of the ElinaTexpr0Array.

    Returns
    -------
    result : c_bool
        Result of the test.

    """

    result = None
    try:
        elina_texpr0_array_is_interval_polyfrac_c = elina_auxiliary_api.elina_texpr0_array_is_interval_polyfrac
        elina_texpr0_array_is_interval_polyfrac_c.restype = c_bool
        elina_texpr0_array_is_interval_polyfrac_c.argtypes = [ElinaTexpr0Array, c_size_t]
        result = elina_texpr0_array_is_interval_polyfrac_c(texpr_array, size)
    except:
        print('Problem with loading/calling "elina_texpr0_array_is_interval_polyfrac" from "libelinaux.so"')
        print('Make sure you are passing ElinaTexpr0Array, c_size_t to the function')

    return result


def elina_texpr0_array_is_scalar(texpr_array, size):
    """
    Test if an ElinaTexpr0Array has only scalar coefficients (non-interval).

    Parameters
    ----------
    texpr_array : ElinaTexpr0Array
        ElinaTexpr0Array that needs to be tested.
    size : c_size_t
        Size of the ElinaTexpr0Array.

    Returns
    -------
    result : c_bool
        Result of the test.

    """

    result = None
    try:
        elina_texpr0_array_is_interval_polyfrac_c = elina_auxiliary_api.elina_texpr0_array_is_interval_polyfrac
        elina_texpr0_array_is_interval_polyfrac_c.restype = c_bool
        elina_texpr0_array_is_interval_polyfrac_c.argtypes = [ElinaTexpr0Array, c_size_t]
        result = elina_texpr0_array_is_interval_polyfrac_c(texpr_array, size)
    except:
        print('Problem with loading/calling "elina_texpr0_array_is_interval_polyfrac" from "libelinaux.so"')
        print('Make sure you are passing ElinaTexpr0Array, c_size_t to the function')

    return result


# ====================================================================== #
# IV. Operations
# ====================================================================== #

def elina_texpr0_substitute(texpr2, dim, texpr3):
    """
    Substitute every occurrence of an ElinaDim in an ElinaTexpr0 with a copy of a different ElinaTexpr0.

    Parameters
    ----------
    texpr2 : ElinaTexpr0Ptr
        Pointer to the ElinaTexpr0 in which we look for occurences of dim.
    dim : ElinaDim
        ElinaDim for which we find occurrences.
    texpr3 : ElinaTexpr0Ptr
        Pointer to the ElinaTexpr0 which we copy to an occurring dimension dim in texpr2.

    Returns
    -------
    texpr1 : ElinaTexpr0Ptr
        Pointer to the newly created ElinaTexpr0 which is a modified texpr2 with substitions of texpr3.

    """

    texpr1 = None
    try:
        elina_texpr0_substitute_c = elina_auxiliary_api.elina_texpr0_substitute
        elina_texpr0_substitute_c.restype = ElinaTexpr0Ptr
        elina_texpr0_substitute_c.argtypes = [ElinaTexpr0Ptr, ElinaDim, ElinaTexpr0Ptr]
        texpr1 = elina_texpr0_substitute_c(texpr2, dim, texpr3)
    except:
        print('Problem with loading/calling "elina_texpr0_substitute" from "libelinaux.so"')
        print('Make sure you are passing ElinaTexpr0Ptr, ElinaDim, ElinaTexpr0Ptr to the function')

    return texpr1


def elina_texpr0_substitute_with(texpr1, dim, texpr2):
    """
    Substitute every occurrence of an ElinaDim in an ElinaTexpr0 with a copy of a different ElinaTexpr0.

    Parameters
    ----------
    texpr1 : ElinaTexpr0Ptr
        Pointer to the ElinaTexpr0 in which we look for occurences of dim.
    dim : ElinaDim
        ElinaDim for which we find occurrences.
    texpr2 : ElinaTexpr0Ptr
        Pointer to the ElinaTexpr0 which we copy to an occurring dimension dim in texpr1.

    Returns
    -------
    None
    
    """

    try:
        elina_texpr0_substitute_with_c = elina_auxiliary_api.elina_texpr0_substitute_with
        elina_texpr0_substitute_with_c.restype = None
        elina_texpr0_substitute_with_c.argtypes = [ElinaTexpr0Ptr, ElinaDim, ElinaTexpr0Ptr]
        elina_texpr0_substitute_with_c(texpr1, dim, texpr2)
    except:
        print('Problem with loading/calling "elina_texpr0_substitute_with" from "libelinaux.so"')
        print('Make sure you are passing ElinaTexpr0Ptr, ElinaDim, ElinaTexpr0Ptr to the function')


# ======================================================================
# V. Change of dimensions and permutations */
# ======================================================================

def elina_texpr0_add_dimensions(texpr2, dimchange):
    """
    Add dimensions to an ElinaTexpr0 following the semantics of an ElinaDimchange.
    
    Parameters
    ----------
    texpr2 : ElinaTexpr0Ptr
        Pointer to the ElinaTexpr0 to which we want to add dimensions.
    dimchange : ElinaDimchangePtr
        Pointer to the ElinaDimchange which semantics we want to follow.

    Returns
    -------
    texpr1 : ElinaTexpr0Ptr
        Pointer to the newly created ElinaTexpr0 with added dimensions.

    """

    texpr1 = None
    try:
        elina_texpr0_add_dimensions_c = elina_auxiliary_api.elina_texpr0_add_dimensions
        elina_texpr0_add_dimensions_c.restype = ElinaTexpr0Ptr
        elina_texpr0_add_dimensions_c.argtypes = [ElinaTexpr0Ptr, ElinaDimchangePtr]
        texpr1 = elina_texpr0_add_dimensions_c(texpr2, dimchange)
    except:
        print('Problem with loading/calling "elina_texpr0_add_dimensions" from "libelinaux.so"')
        print('Make sure you are passing ElinaTexpr0Ptr, ElinaDimchangePtr to the function')

    return texpr1


def elina_texpr0_remove_dimensions(texpr2, dimchange):
    """
    Remove dimensions of an ElinaTexpr0 following the semantics of an ElinaDimchange.

    Parameters
    ----------
    texpr2 : ElinaTexpr0Ptr
        Pointer to the ElinaTexpr0 from which we want to remove dimensions.
    dimchange : ElinaDimchangePtr
        Pointer to the ElinaDimchange which semantics we want to follow.

    Returns
    -------
    texpr1 : ElinaTexpr0Ptr
        Pointer to the newly created ElinaTexpr0 with removed dimensions.

    """

    texpr1 = None
    try:
        elina_texpr0_remove_dimensions_c = elina_auxiliary_api.elina_texpr0_remove_dimensions
        elina_texpr0_remove_dimensions_c.restype = ElinaTexpr0Ptr
        elina_texpr0_remove_dimensions_c.argtypes = [ElinaTexpr0Ptr, ElinaDimchangePtr]
        texpr1 = elina_texpr0_remove_dimensions_c(texpr2, dimchange)
    except:
        print('Problem with loading/calling "elina_texpr0_remove_dimensions" from "libelinaux.so"')
        print('Make sure you are passing ElinaTexpr0Ptr, ElinaDimchangePtr to the function')

    return texpr1


def elina_texpr0_permute_dimensions(texpr2, dimperm):
    """
    Permute dimensions of an ElinaTexpr0 following the semantics of an ElinaDimperm.

    Parameters
    ----------
    texpr2 : ElinaTexpr0Ptr
        Pointer to the ElinaTexpr0 which dimensions we want to permute.
    dimperm : ElinaDimpermPtr
        Pointer to the ElinaDimperm which semantics we want to follow.

    Returns
    -------
    texpr1 : ElinaTexpr0Ptr
        Pointer to the newly created ElinaTexpr0 with permuted dimensions.

    """

    texpr1 = None
    try:
        elina_texpr0_permute_dimensions_c = elina_auxiliary_api.elina_texpr0_permute_dimensions
        elina_texpr0_permute_dimensions_c.restype = ElinaTexpr0Ptr
        elina_texpr0_permute_dimensions_c.argtypes = [ElinaTexpr0Ptr, ElinaDimpermPtr]
        texpr1 = elina_texpr0_permute_dimensions_c(texpr2, dimperm)
    except:
        print('Problem with loading/calling "elina_texpr0_permute_dimensions" from "libelinaux.so"')
        print('Make sure you are passing ElinaTexpr0Ptr, ElinaDimpermPtr to the function')

    return texpr1


def elina_texpr0_add_dimensions_with(texpr, dimchange):
    """
    Add dimensions to an ElinaTexpr0 following the semantics of an ElinaDimchange.
    
    Parameters
    ----------
    texpr : ElinaTexpr0Ptr
        Pointer to the ElinaTexpr0 to which we want to add dimensions.
    dimchange : ElinaDimchangePtr
        Pointer to the ElinaDimchange which semantics we want to follow.
    
    Returns
    -------
    None
        
    """

    try:
        elina_texpr0_add_dimensions_with_c = elina_auxiliary_api.elina_texpr0_add_dimensions_with
        elina_texpr0_add_dimensions_with_c.restype = None
        elina_texpr0_add_dimensions_with_c.argtypes = [ElinaTexpr0Ptr, ElinaDimchangePtr]
        elina_texpr0_add_dimensions_with_c(texpr, dimchange)
    except:
        print('Problem with loading/calling "elina_texpr0_add_dimensions_with" from "libelinaux.so"')
        print('Make sure you are passing ElinaTexpr0Ptr, ElinaDimchangePtr to the function')


def elina_texpr0_remove_dimensions_with(texpr, dimchange):
    """
    Remove dimensions of an ElinaTexpr0 following the semantics of an ElinaDimchange.

    Parameters
    ----------
    texpr : ElinaTexpr0Ptr
        Pointer to the ElinaTexpr0 from which we want to remove dimensions.
    dimchange : ElinaDimchangePtr
        Pointer to the ElinaDimchange which semantics we want to follow.

    Returns
    -------
    None

    """

    try:
        elina_texpr0_remove_dimensions_with_c = elina_auxiliary_api.elina_texpr0_remove_dimensions_with
        elina_texpr0_remove_dimensions_with_c.restype = None
        elina_texpr0_remove_dimensions_with_c.argtypes = [ElinaTexpr0Ptr, ElinaDimchangePtr]
        elina_texpr0_remove_dimensions_with_c(texpr, dimchange)
    except:
        print('Problem with loading/calling "elina_texpr0_remove_dimensions_with" from "libelinaux.so"')
        print('Make sure you are passing ElinaTexpr0Ptr, ElinaDimchangePtr to the function')


def elina_texpr0_permute_dimensions_with(texpr, dimperm):
    """
    Permute dimensions of an ElinaTexpr0 following the semantics of an ElinaDimperm.

    Parameters
    ----------
    texpr : ElinaTexpr0Ptr
        Pointer to the ElinaTexpr0 which dimensions we want to permute.
    dimperm : ElinaDimpermPtr
        Pointer to the ElinaDimperm which semantics we want to follow.

    Returns
    -------
    None

    """

    try:
        elina_texpr0_permute_dimensions_with_c = elina_auxiliary_api.elina_texpr0_permute_dimensions_with
        elina_texpr0_permute_dimensions_with_c.restype = None
        elina_texpr0_permute_dimensions_with_c.argtypes = [ElinaTexpr0Ptr, ElinaDimpermPtr]
        elina_texpr0_permute_dimensions_with_c(texpr, dimperm)
    except:
        print('Problem with loading/calling "elina_texpr0_permute_dimensions_with" from "libelinaux.so"')
        print('Make sure you are passing ElinaTexpr0Ptr, ElinaDimpermPtr to the function')


# ====================================================================== #
# VI. Hashing, comparisons
# ====================================================================== #

def elina_texpr0_hash(texpr):
    """
    Hash an ElinaTexpr0.
    
    Parameters
    ----------
    texpr : ElinaTexpr0Ptr
        Pointer to the ElinaTexpr0 that needs to be hashed.

    Returns
    -------
    result : c_long
        The resulting hash.

    """

    result = None
    try:
        elina_texpr0_hash_c = elina_auxiliary_api.elina_texpr0_hash
        elina_texpr0_hash_c.restype = c_long
        elina_texpr0_hash_c.argtypes = [ElinaTexpr0Ptr, ElinaDimpermPtr]
        result = elina_texpr0_hash_c(texpr)
    except:
        print('Problem with loading/calling "elina_texpr0_hash" from "libelinaux.so"')
        print('Make sure you are passing ElinaTexpr0Ptr to the function')

    return result


def elina_texpr0_equal(texpr1, texpr2):
    """
    Test if two ElinaTexpr0 are equal.
    
    Parameters
    ----------
    texpr1 : ElinaTexpr0Ptr
        Pointer to the first ElinaTexpr0 that needs to be tested.
    texpr2: ElinaTexpr0Ptr
        Pointer to the second ElinaTexpr0 that needs to be tested.

    Returns
    -------
    result : c_bool
        Result of the test.
    

    """

    result = None
    try:
        elina_texpr0_equal_c = elina_auxiliary_api.elina_texpr0_equal
        elina_texpr0_equal_c.restype = c_bool
        elina_texpr0_equal_c.argtypes = [ElinaTexpr0Ptr, ElinaTexpr0Ptr]
        result = elina_texpr0_equal_c(texpr1, texpr2)
    except:
        print('Problem with loading/calling "elina_texpr0_equal" from "libelinaux.so"')
        print('Make sure you are passing ElinaTexpr0Ptr, ElinaTexpr0Ptr to the function')

    return result
