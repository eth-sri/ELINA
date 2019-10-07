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


from zonoml_imports import *
from elina_manager_h import *
from elina_abstract0_h import *
import numpy as np 
from numpy.ctypeslib import ndpointer
import ctypes


_doublepp = ndpointer(dtype=np.uintp, ndim=1, flags='C')

# ====================================================================== #
# Basics
# ====================================================================== #

def zonoml_manager_alloc():
    """
    Allocates an ElinaManager.

    Returns
    -------
    man : ElinaManagerPtr
        Pointer to the newly allocated ElinaManager.

    """

    man = None
    try:
        zonoml_manager_alloc_c = zonoml_api.zonoml_manager_alloc
        zonoml_manager_alloc_c.restype = ElinaManagerPtr
        zonoml_manager_alloc_c.argtypes = None
        man = zonoml_manager_alloc_c()
    except:
        print('Problem with loading/calling "zonoml_manager_alloc" from "libzonoml.so"')

    return man


def zonotope_from_network_input(man, intdim, realdim, inf_array, sup_array):
    """
    Create the perturbed zonotope from input
    
    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    intdim : c_size_t
        Number of integer variables.
    realdim: c_size_t
        Number of real variables
    inf_array: POINTER(double)
        lower bound array
    sup_array: POINTER(double)
        upper bound array
    
    Returns
    -------
    res: ElinaAbstract0Ptr
         Pointer to the new abstract object

    """

    try:
        zonotope_from_network_input_c = zonoml_api.zonotope_from_network_input
        zonotope_from_network_input_c.restype = ElinaAbstract0Ptr
        zonotope_from_network_input_c.argtypes = [ElinaManagerPtr, c_size_t, c_size_t, ndpointer(ctypes.c_double), ndpointer(ctypes.c_double)]
        res = zonotope_from_network_input_c(man,intdim, realdim, inf_array, sup_array)
    except Exception as inst:
        print('Problem with loading/calling "zonotope_from_network_input" from "libzonoml.so"')
        print(inst)
    return res


def elina_abstract0_from_zonotope(man, intdim, realdim, num_error_terms, zonotope):
    """
    Create the perturbed zonotope from input
    
    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    intdim : c_size_t
        Number of integer variables.
    realdim: c_size_t
        Number of real variables
    inf_array: POINTER(double)
        lower bound array
    sup_array: POINTER(double)
        upper bound array
    
    Returns
    -------
    res: ElinaAbstract0Ptr
         Pointer to the new abstract object

    """

    try:
        elina_abstract0_from_zonotope_c = zonoml_api.elina_abstract0_from_zonotope
        elina_abstract0_from_zonotope_c.restype = ElinaAbstract0Ptr
        elina_abstract0_from_zonotope_c.argtypes = [ElinaManagerPtr, c_size_t, c_size_t, c_size_t, _doublepp]
        res = elina_abstract0_from_zonotope_c(man,intdim, realdim, num_error_terms, zonotope)
    except Exception as inst:
        print('Problem with loading/calling "elina_abstract0_from_zonotope" from "libzonoml.so"')
        print(inst)
    return res



def ffn_matmult_zono(man, destructive, element, start_offset, weights, bias, num_var, expr_offset, expr_size):
    """
    FFN Matrix multiplication
    
    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    destructive: c_bool
        Boolean flag
    element : ElinaAbstract0Ptr
        Pointer to the ElinaAbstract0 which dimensions need to be assigned.
    start_offset: ElinaDim
        The start offset from which the dimensions should be assigned.
    weights: _doublepp
        weight matrix
    bias: POINTER(double)
        bias vector
    num_var: c_size_t
        number of neurons to be assigned
    expr_offset: c_size_t
        the offset of the first variable in the assignment expression
    expr_size: c_size_t  
        number of variables in an assignment expression
    Returns
    -------
    res: ElinaAbstract0Ptr
         Pointer to the new abstract object

    """

    try:
        ffn_matmult_zono_c = zonoml_api.ffn_matmult_zono
        ffn_matmult_zono_c.restype = ElinaAbstract0Ptr
        ffn_matmult_zono_c.argtypes = [ElinaManagerPtr, c_bool,  ElinaAbstract0Ptr, ElinaDim, _doublepp, ndpointer(ctypes.c_double), c_size_t, c_size_t, c_size_t]
        res = ffn_matmult_zono_c(man, destructive, element, start_offset, weights, bias, num_var, expr_offset, expr_size)
    except Exception as inst:
        print('Problem with loading/calling "ffn_matmult_zono" from "libzonoml.so"')
        print(inst)

    return res


def ffn_matmult_without_bias_zono(man, destructive, element, start_offset, weights, num_var, expr_offset, expr_size):
    """
    FFN Matrix multiplication without bias
    
    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    destructive: c_bool
        Boolean flag
    element : ElinaAbstract0Ptr
        Pointer to the ElinaAbstract0 which dimensions need to be assigned.
    start_offset: ElinaDim
        The start offset from which the dimensions should be assigned.
    weights: _doublepp
        weight matrix
    num_var: c_size_t
        number of neurons to be assigned
    expr_offset: c_size_t
        the offset of the first variable in the assignment expression
    expr_size: c_size_t  
        number of variables in an assignment expression
    Returns
    -------
    res: ElinaAbstract0Ptr
         Pointer to the new abstract object

    """

    try:
        ffn_matmult_without_bias_zono_c = zonoml_api.ffn_matmult_without_bias_zono
        ffn_matmult_without_bias_zono_c.restype = ElinaAbstract0Ptr
        ffn_matmult_without_bias_zono_c.argtypes = [ElinaManagerPtr, c_bool,  ElinaAbstract0Ptr, ElinaDim, _doublepp, c_size_t, c_size_t, c_size_t]
        res = ffn_matmult_without_bias_zono_c(man, destructive, element, start_offset, weights,  num_var, expr_offset, expr_size)
    except Exception as inst:
        print('Problem with loading/calling "ffn_matmult_without_bias_zono" from "libzonoml.so"')
        print(inst)

    return res


def ffn_add_bias_zono(man, destructive, element, start_offset,  bias, num_var):
    """
    FFN bias add
    
    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    destructive: c_bool
        Boolean flag
    element : ElinaAbstract0Ptr
        Pointer to the ElinaAbstract0 which dimensions need to be assigned.
    start_offset: ElinaDim
        The start offset from which the dimensions should be assigned.
    bias: POINTER(double)
        bias vector
    num_var: c_size_t
        number of neurons to be assigned
    Returns
    -------
    res: ElinaAbstract0Ptr
         Pointer to the new abstract object

    """

    try:
        ffn_add_bias_zono_c = zonoml_api.ffn_add_bias_zono
        ffn_add_bias_zono_c.restype = ElinaAbstract0Ptr
        ffn_add_bias_zono_c.argtypes = [ElinaManagerPtr, c_bool,  ElinaAbstract0Ptr, ElinaDim, ndpointer(ctypes.c_double), c_size_t]
        res = ffn_add_bias_zono_c(man, destructive, element, start_offset, bias, num_var)
    except Exception as inst:
        print('Problem with loading/calling "ffn_add_bias_zono" from "libzonoml.so"')
        print(inst)

    return res   

def ffn_sub_bias_zono(man, destructive, element, start_offset,  bias, is_minuend, num_var):
    """
    FFN bias add
    
    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    destructive: c_bool
        Boolean flag
    element : ElinaAbstract0Ptr
        Pointer to the ElinaAbstract0 which dimensions need to be assigned.
    start_offset: ElinaDim
        The start offset from which the dimensions should be assigned.
    bias: POINTER(double)
        bias vector
    num_var: c_size_t
        number of neurons to be assigned
    Returns
    -------
    res: ElinaAbstract0Ptr
         Pointer to the new abstract object

    """

    try:
        ffn_sub_bias_zono_c = zonoml_api.ffn_sub_bias_zono
        ffn_sub_bias_zono_c.restype = ElinaAbstract0Ptr
        ffn_sub_bias_zono_c.argtypes = [ElinaManagerPtr, c_bool,  ElinaAbstract0Ptr, ElinaDim, ndpointer(ctypes.c_double), c_bool, c_size_t]
        res = ffn_sub_bias_zono_c(man, destructive, element, start_offset, bias, is_minuend, num_var)
    except Exception as inst:
        print('Problem with loading/calling "ffn_sub_bias_zono" from "libzonoml.so"')
        print(inst)

    return res 


def ffn_mul_bias_zono(man, destructive, element, start_offset,  bias, num_var):
    """
    FFN bias add
    
    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    destructive: c_bool
        Boolean flag
    element : ElinaAbstract0Ptr
        Pointer to the ElinaAbstract0 which dimensions need to be assigned.
    start_offset: ElinaDim
        The start offset from which the dimensions should be assigned.
    bias: POINTER(double)
        bias vector
    num_var: c_size_t
        number of neurons to be assigned
    Returns
    -------
    res: ElinaAbstract0Ptr
         Pointer to the new abstract object

    """

    try:
        ffn_mul_bias_zono_c = zonoml_api.ffn_mul_bias_zono
        ffn_mul_bias_zono_c.restype = ElinaAbstract0Ptr
        ffn_mul_bias_zono_c.argtypes = [ElinaManagerPtr, c_bool,  ElinaAbstract0Ptr, ElinaDim, ndpointer(ctypes.c_double), c_size_t]
        res = ffn_mul_bias_zono_c(man, destructive, element, start_offset, bias, num_var)
    except Exception as inst:
        print('Problem with loading/calling "ffn_mul_bias_zono" from "libzonoml.so"')
        print(inst)

    return res 


def conv_matmult_zono(man, destructive, element, start_offset, filter_weights, filter_bias, input_size, expr_offset, filter_size, num_filters, strides, output_size, pad_top, pad_left, has_bias):
    """
    Convolutional Matrix multiplication
    
    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    destructive: c_bool
        Boolean flag
    element : ElinaAbstract0Ptr
        Pointer to the ElinaAbstract0 which dimensions need to be assigned.
    start_offset: ElinaDim
        The start offset from which the dimensions should be assigned.
    filter_weights: POINTER(double)
        filter weights
    filter_bias: POINTER(double)
        filter biases
    input_size: POINTER(c_size_t)
        size of the input
    expr_offset: c_size_t
        the offset of the first variable in the assignment expression
    filter_size: POINTER(c_size_t)  
        size of the filters
    num_filters: c_size_t
        number of filters
    strides: POINTER(c_size_t)
       size of the strides
    is_valid_padding: c_bool
       if the padding is valid
    has_bias: c_bool
       if the filter has bias
    Returns
    -------
    res: ElinaAbstract0Ptr
         Pointer to the new abstract object

    """
    try:
        conv_matmult_zono_c = zonoml_api.conv_matmult_zono
        conv_matmult_zono_c.restype = ElinaAbstract0Ptr
        conv_matmult_zono_c.argtypes = [ElinaManagerPtr, c_bool, ElinaAbstract0Ptr, ElinaDim, ndpointer(ctypes.c_double), ndpointer(ctypes.c_double), POINTER(c_size_t), c_size_t, POINTER(c_size_t), c_size_t, POINTER(c_size_t), POINTER(c_size_t), c_size_t, c_size_t, c_bool]
        res = conv_matmult_zono_c(man, destructive, element, start_offset, filter_weights, filter_bias, input_size, expr_offset, filter_size, num_filters, strides, output_size, pad_top, pad_left, has_bias)
    except Exception as inst:
        print('Problem with loading/calling "conv_matmult_zono" from "libzonoml.so"')
        print(inst)
    return res


def relu_zono(man,destructive,elem,x):
    """
    Performs the ReLU operation
    
    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    destructive : c_bool
        Boolean flag.
    elem : ElinaAbstract0Ptr
        Pointer to the ElinaAbstract0 which dimensions need to be assigned.
    tdim : ElinaDim
        The dimension to be assigned.
    Returns
    -------
    res : ElinaAbstract0Ptr
        Pointer to the new abstract object.

    """

    res = None
    try:
        relu_zono_c = zonoml_api.relu_zono
        relu_zono_c.restype = ElinaAbstract0Ptr
        relu_zono_c.argtypes = [ElinaManagerPtr, c_bool, ElinaAbstract0Ptr, ElinaDim]
        res = relu_zono_c(man,destructive,elem,x)
    except:
        print('Problem with loading/calling "relu_zono" from "libzonoml.so"')

    return res

def relu_zono_refined(man,destructive,elem,x, new_inf, new_sup):
    """
    Performs the ReLU operation refined
    
    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    destructive : c_bool
        Boolean flag.
    elem : ElinaAbstract0Ptr
        Pointer to the ElinaAbstract0 which dimensions need to be assigned.
    tdim : ElinaDim
        The dimension to be assigned.
    new_inf: c_double
        The modified lower bound
    new_sup: c_double
        The modified upper bound
    Returns
    -------
    res : ElinaAbstract0Ptr
        Pointer to the new abstract object.

    """

    res = None
    try:
        relu_zono_refined_c = zonoml_api.relu_zono_refined
        relu_zono_refined_c.restype = ElinaAbstract0Ptr
        relu_zono_refined_c.argtypes = [ElinaManagerPtr, c_bool, ElinaAbstract0Ptr, ElinaDim, c_double, c_double]
        res = relu_zono_refined_c(man,destructive,elem,x, new_inf, new_sup)
    except:
        print('Problem with loading/calling "relu_zono_refined" from "libzonoml.so"')

    return res


def maxpool_zono_refined(man,destructive,elem,x, new_inf, new_sup):
    """
    Performs the Maxpool operation refined
    
    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    destructive : c_bool
        Boolean flag.
    elem : ElinaAbstract0Ptr
        Pointer to the ElinaAbstract0 which dimensions need to be assigned.
    tdim : ElinaDim
        The dimension to be assigned.
    new_inf: c_double
        The modified lower bound
    new_sup: c_double
        The modified upper bound
    Returns
    -------
    res : ElinaAbstract0Ptr
        Pointer to the new abstract object.

    """

    res = None
    try:
        maxpool_zono_refined_c = zonoml_api.maxpool_zono_refined
        maxpool_zono_refined_c.restype = ElinaAbstract0Ptr
        maxpool_zono_refined_c.argtypes = [ElinaManagerPtr, c_bool, ElinaAbstract0Ptr, ElinaDim, c_double, c_double]
        res = maxpool_zono_refined_c(man,destructive,elem,x, new_inf, new_sup)
    except:
        print('Problem with loading/calling "maxpool_zono_refined" from "libzonoml.so"')

    return res


def relu_zono_layerwise(man,destructive,elem,start_offset, num_dim):
    """
    Performs the ReLU operation
    
    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    destructive : c_bool
        Boolean flag.
    elem : ElinaAbstract0Ptr
        Pointer to the ElinaAbstract0 which dimensions need to be assigned.
    start_offset : ElinaDim
        The starting dimension.
    num_dim : ElinaDim
        The number of variables on which relu should be applied

    Returns
    -------
    res : ElinaAbstract0Ptr
        Pointer to the new abstract object.

    """

    res = None
    try:
        relu_zono_layerwise_c = zonoml_api.relu_zono_layerwise
        relu_zono_layerwise_c.restype = ElinaAbstract0Ptr
        relu_zono_layerwise_c.argtypes = [ElinaManagerPtr, c_bool, ElinaAbstract0Ptr, ElinaDim, ElinaDim]
        res = relu_zono_layerwise_c(man,destructive,elem,start_offset, num_dim)
    except:
        print('Problem with loading/calling "relu_zono_layerwise" from "libzonoml.so"')

    return res

def sigmoid_zono(man,destructive,elem,x):
    """
    Performs the Sigmoid operation
    
    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    destructive : c_bool
        Boolean flag.
    elem : ElinaAbstract0Ptr
        Pointer to the ElinaAbstract0 which dimensions need to be assigned.
    tdim : ElinaDim
        The dimension to be assigned.
    Returns
    -------
    res : ElinaAbstract0Ptr
        Pointer to the new abstract object.

    """

    res = None
    try:
        sigmoid_zono_c = zonoml_api.sigmoid_zono
        sigmoid_zono_c.restype = ElinaAbstract0Ptr
        sigmoid_zono_c.argtypes = [ElinaManagerPtr, c_bool, ElinaAbstract0Ptr, ElinaDim]
        res = sigmoid_zono_c(man,destructive,elem,x)
    except:
        print('Problem with loading/calling "sigmoid_zono" from "libzonoml.so"')

    return res

def sigmoid_zono_layerwise(man,destructive,elem,start_offset, num_dim):
    """
    Performs the Sigmoid operation layerwise
    
    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    destructive : c_bool
        Boolean flag.
    elem : ElinaAbstract0Ptr
        Pointer to the ElinaAbstract0 which dimensions need to be assigned.
    start_offset : ElinaDim
        The starting dimension.
    num_dim : ElinaDim
        The number of variables on which sigmoid should be applied

    Returns
    -------
    res : ElinaAbstract0Ptr
        Pointer to the new abstract object.

    """

    res = None
    try:
        sigmoid_zono_layerwise_c = zonoml_api.sigmoid_zono_layerwise
        sigmoid_zono_layerwise_c.restype = ElinaAbstract0Ptr
        sigmoid_zono_layerwise_c.argtypes = [ElinaManagerPtr, c_bool, ElinaAbstract0Ptr, ElinaDim, ElinaDim]
        res = sigmoid_zono_layerwise_c(man,destructive,elem,start_offset, num_dim)
    except:
        print('Problem with loading/calling "sigmoid_zono_layerwise" from "libzonoml.so"')

    return res

def tanh_zono(man,destructive,elem,x):
    """
    Performs the Tanh operation
    
    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    destructive : c_bool
        Boolean flag.
    elem : ElinaAbstract0Ptr
        Pointer to the ElinaAbstract0 which dimensions need to be assigned.
    tdim : ElinaDim
        The dimension to be assigned.
    Returns
    -------
    res : ElinaAbstract0Ptr
        Pointer to the new abstract object.

    """

    res = None
    try:
        tanh_zono_c = zonoml_api.tanh_zono
        tanh_zono_c.restype = ElinaAbstract0Ptr
        tanh_zono_c.argtypes = [ElinaManagerPtr, c_bool, ElinaAbstract0Ptr, ElinaDim]
        res = tanh_zono_c(man,destructive,elem,x)
    except:
        print('Problem with loading/calling "tanh_zono" from "libzonoml.so"')

    return res


def tanh_zono_layerwise(man,destructive,elem,start_offset, num_dim):
    """
    Performs the Tanh operation layerwise
    
    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    destructive : c_bool
        Boolean flag.
    elem : ElinaAbstract0Ptr
        Pointer to the ElinaAbstract0 which dimensions need to be assigned.
    start_offset : ElinaDim
        The starting dimension.
    num_dim : ElinaDim
        The number of variables on which tanh should be applied

    Returns
    -------
    res : ElinaAbstract0Ptr
        Pointer to the new abstract object.

    """

    res = None
    try:
        tanh_zono_layerwise_c = zonoml_api.tanh_zono_layerwise
        tanh_zono_layerwise_c.restype = ElinaAbstract0Ptr
        tanh_zono_layerwise_c.argtypes = [ElinaManagerPtr, c_bool, ElinaAbstract0Ptr, ElinaDim, ElinaDim]
        res = tanh_zono_layerwise_c(man,destructive,elem,start_offset, num_dim)
    except:
        print('Problem with loading/calling "tanh_zono_layerwise" from "libzonoml.so"')

    return res

def maxpool_zono(man, destructive, elem, pool_size, input_size, src_offset, strides, dimensionality, dst_offset, is_valid_padding):
    """
    Performs the Maxpool operation
    
    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    destructive : c_bool
        Boolean flag.
    elem : ElinaAbstract0Ptr
        Pointer to the ElinaAbstract0 which dimensions need to be assigned.
    pool_size: POINTER(c_size_t)
        The size of the Maxpool filter 
    input_size : POINTER(c_size_t)
        The number of variables on which Maxpool will be applied.
    src_offset: c_size_t
        The source offset in the abstract element for Maxpool
    strides: POINTER(c_size_t)
	The size of the sliding window
    dimensionality: c_size_t
        The number of the dimensions in the input and the Maxpool filter
    dst_offset: c_size_t
        The destination offset in the abstract element for Maxpool
    is_valid_padding: c_bool
        whether the padding is valid or same
    Returns
    -------
    res : ElinaAbstract0Ptr
        Pointer to the new abstract object.

    """

    res = None
    try:
        maxpool_zono_c = zonoml_api.maxpool_zono
        maxpool_zono_c.restype = ElinaAbstract0Ptr
        maxpool_zono_c.argtypes = [ElinaManagerPtr, c_bool, ElinaAbstract0Ptr, POINTER(c_size_t), POINTER(c_size_t), c_size_t, POINTER(c_size_t), c_size_t, c_size_t,c_bool]
        res = maxpool_zono_c(man,destructive,elem,pool_size,input_size, src_offset, strides, dimensionality, dst_offset,is_valid_padding)
    except Exception as inst:
        print('Problem with loading/calling "maxpool_zono" from "libzonoml.so"')
        print(inst)

    return res


def is_greater_zono(man, element, y, x):
    """
    Check if y is strictly greater than x in the abstract element 
    
    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    destructive : c_bool
        Boolean flag.
    y : ElinaDim
        The dimension y in the constraint y-x>0.
    x: ElinaDim
	The dimension x in the constraint y-x>0.
    
    Returns
    -------
    res = boolean

    """
    res= None
    try:
        is_greater_c = zonoml_api.is_greater
        is_greater_c.restype = c_bool
        is_greater_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr, ElinaDim, ElinaDim]
        res = is_greater_c(man,element,y, x)
    except Exception as inst:
        print('Problem with loading/calling "is_greater" from "libzonoml.so"')
        print(inst)
    return res


def affine_form_is_box(man, element, x):
    """
    Check if the affine form for x in the abstract element is a box 
    
    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    destructive : c_bool
        Boolean flag.
    x: ElinaDim
	The dimension x in the constraint y-x>0.
    
    Returns
    -------
    res = boolean

    """
    res= None
    try:
        affine_form_is_box_c = zonoml_api.affine_form_is_box
        affine_form_is_box_c.restype = c_bool
        affine_form_is_box_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr, ElinaDim]
        res = affine_form_is_box_c(man,element, x)
    except Exception as inst:
        print('Problem with loading/calling "affine_form_is_box" from "libzonoml.so"')
        print(inst)
    return res

def zono_add(man, element, dst_offset, src_offset, num_var):
    """
    Add the affine forms (y:=y+x) in different sections of the abstract element
    
    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    element : ElinaAbstract0Ptr
        Abstract element.
    dst_offset : c_size_t
        The destination offset (y)
    src_offset: c_size_t
	The source offset (x)
    num_var: c_size_t
        number of variables
    Returns
    -------
    None

    """
    
    try:
        zono_add_c = zonoml_api.zono_add
        zono_add_c.restype = None
        zono_add_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr, c_size_t, c_size_t, c_size_t]
        zono_add_c(man,element,dst_offset, src_offset, num_var)
    except Exception as inst:
        print('Problem with loading/calling "zono_add" from "libzonoml.so"')
        print(inst)


def zono_copy_section(man, element, dst_offset, src_offset, num_var):
    """
    copy affine forms from one section of the abstract element to another 
    
    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    element : ElinaAbstract0Ptr
        Abstract element.
    dst_offset : c_size_t
        The destination offset
    src_offset: c_size_t
	The source offset
    num_var: c_size_t
        number of variables
    Returns
    -------
    None

    """
    
    try:
        zono_copy_section_c = zonoml_api.zono_copy_section
        zono_copy_section_c.restype = None
        zono_copy_section_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr, c_size_t, c_size_t, c_size_t]
        zono_copy_section_c(man,element,dst_offset, src_offset, num_var)
    except Exception as inst:
        print('Problem with loading/calling "zono_copy_section" from "libzonoml.so"')
        print(inst)


def get_interval_width_var_zono(man, element, i):
    """
    get interval width for an affine form
    
    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    element : ElinaAbstract0Ptr
        Abstract element.
    i : c_size_t
        The index for the affine form
    Returns
    -------
    width = c_double

    """
    
    try:
        get_interval_width_var_zono_c = zonoml_api.get_interval_width_var_zono
        get_interval_width_var_zono_c.restype = c_double
        get_interval_width_var_zono_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr, c_size_t]
        width = get_interval_width_var_zono_c(man,element,i)
    except Exception as inst:
        print('Problem with loading/calling "get_interval_width_var_zono" from "libzonoml.so"')
        print(inst)

    return width


def handle_gather_layer(man, destructive, element, indexes):
    """
    handle the gather operation
    
    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    destructive: c_bool
        whether the returned element is the same as input
    element : ElinaAbstract0Ptr
        Abstract element.
    indexes : POINTER(c_size_t)
        The indexes for the gather
    Returns
    -------
    width = c_double

    """
    
    try:
        handle_gather_layer_c = zonoml_api.handle_gather_layer
        handle_gather_layer_c.restype = ElinaAbstract0Ptr
        handle_gather_layer_c.argtypes = [ElinaManagerPtr, c_bool, ElinaAbstract0Ptr, ndpointer(c_size_t)]
        res = handle_gather_layer_c(man,destructive, element, indexes)
    except Exception as inst:
        print('Problem with loading/calling "handle_gather_layer" from "libzonoml.so"')
        print(inst)

    return res


def get_affine_form_for_dim(man,element,dim):
    """
    get the list of coefficients for the given affine form
    
    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    element : ElinaAbstract0Ptr
        Abstract element.
    dim : c_size_t
        The dimension of the affine form
    Returns
    -------
    array = POINTER(c_double)

    """
    
    try:
        get_affine_form_for_dim_c = zonoml_api.get_affine_form_for_dim
        get_affine_form_for_dim_c.restype = POINTER(c_double)
        get_affine_form_for_dim_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr, c_size_t]
        array = get_affine_form_for_dim_c(man, element, dim)
    except Exception as inst:
        print('Problem with loading/calling "get_affine_form_for_dim" from "libzonoml.so"')
        print(inst)

    return array   
