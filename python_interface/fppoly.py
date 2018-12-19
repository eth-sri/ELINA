#
#
#  This source file is part of ELINA (ETH LIbrary for Numerical Analysis).
#  ELINA is Copyright  2018 Department of Computer Science, ETH Zurich
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


from fppoly_imports import *
from elina_manager_h import *
from elina_abstract0_h import *
from elina_interval_h import *
import numpy as np 
from numpy.ctypeslib import ndpointer
import ctypes

_doublepp = ndpointer(dtype=np.uintp, ndim=1, flags='C')
# ====================================================================== #
# Basics
# ====================================================================== #

def fppoly_manager_alloc():
    """
    Allocates an ElinaManager.

    Returns
    -------
    man : ElinaManagerPtr
        Pointer to the newly allocated ElinaManager.

    """

    man = None
    try:
        fppoly_manager_alloc_c = fppoly_api.fppoly_manager_alloc
        fppoly_manager_alloc_c.restype = ElinaManagerPtr
        fppoly_manager_alloc_c.argtypes = None
        man = fppoly_manager_alloc_c()
    except:
        print('Problem with loading/calling "fppoly_manager_alloc" from "libfppoly.so"')

    return man

def fppoly_from_network_input(man, intdim, realdim, inf_array, sup_array):
    """
    Create an abstract element from perturbed input
    
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

    res = None
    try:
        fppoly_from_network_input_c = fppoly_api.fppoly_from_network_input
        fppoly_from_network_input_c.restype = ElinaAbstract0Ptr
        fppoly_from_network_input_c.argtypes = [ElinaManagerPtr, c_size_t, c_size_t,ndpointer(ctypes.c_double),ndpointer(ctypes.c_double)]
        res = fppoly_from_network_input_c(man,intdim, realdim, inf_array,sup_array)
    except Exception as inst:
        print('Problem with loading/calling "fppoly_from_network_input" from "libfppoly.so"')
        print(inst)	

    return res

def fppoly_from_network_input_poly(man, intdim, realdim, inf_array, sup_array, lexpr_weights, lexpr_cst, lexpr_dim,  uexpr_weights, uexpr_cst, uexpr_dim, expr_size):
    """
    Create an abstract element from perturbed input
    
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
    lexpr_weights: POINTER(double)
        coefficients of the lower polyhedra constraints
    lexpr_cst: POINTER(double)
        constants of the lower polyhedra constraints
    lexpr_dim: POINTER(c_size_t)
        the indexes of the variables in the lower polyhedra constraints
    uexpr_weights: POINTER(double)
        coefficients of the upper polyhedra constraints
    uexpr_cst: POINTER(double)
        constants of the upper polyhedra constraints
    uexpr_dim: POINTER(c_size_t)
        the indexes of the variables in the upper polyhedra constraints
    expr_size: c_size_t
        size of the polyhedra constraints
    Returns
    -------
    res: ElinaAbstract0Ptr
         Pointer to the new abstract object

    """

    res = None
    try:
        fppoly_from_network_input_poly_c = fppoly_api.fppoly_from_network_input_poly
        fppoly_from_network_input_poly_c.restype = ElinaAbstract0Ptr
        fppoly_from_network_input_poly_c.argtypes = [ElinaManagerPtr, c_size_t, c_size_t,ndpointer(ctypes.c_double),ndpointer(ctypes.c_double),ndpointer(ctypes.c_double),ndpointer(ctypes.c_double),ndpointer(ctypes.c_size_t),ndpointer(ctypes.c_double),ndpointer(ctypes.c_double),ndpointer(ctypes.c_size_t), c_size_t]
        res = fppoly_from_network_input_poly_c(man,intdim, realdim, inf_array,sup_array, lexpr_weights, lexpr_cst, lexpr_dim, uexpr_weights, uexpr_cst, uexpr_dim ,expr_size)
    except Exception as inst:
        print('Problem with loading/calling "fppoly_from_network_input_poly" from "libfppoly.so"')
        print(inst)	

    return res





def ffn_handle_first_relu_layer(man, element,weights, bias,  size, num_pixels):
    """
    handle the first FFN ReLU layer
    
    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    weights : POINTER(POINTER(c_double))
        The weight matrix.
    bias : POINTER(c_double)
        The bias vector
    size: c_size_t
	Number of neurons in the first layer
    num_pixels:
        Number of pixels in the input
    Returns
    -------
    res : ElinaAbstract0Ptr
        Pointer to the new abstract object.

    """

    try:
        ffn_handle_first_relu_layer_c = fppoly_api.ffn_handle_first_relu_layer
        ffn_handle_first_relu_layer_c.restype = None
        ffn_handle_first_relu_layer_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr, _doublepp, ndpointer(ctypes.c_double),  c_size_t, c_size_t]
        ffn_handle_first_relu_layer_c(man,element,weights, bias,  size, num_pixels)
    except Exception as inst:
        print('Problem with loading/calling "ffn_handle_first_relu_layer" from "libfppoly.so"')
        print(inst)	

    return

def ffn_handle_first_sigmoid_layer(man, element,weights, bias,  size, num_pixels):
    """
    handle the FFN first Sigmoid layer
    
    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    weights : POINTER(POINTER(c_double))
        The weight matrix.
    bias : POINTER(c_double)
        The bias vector
    size: c_size_t
	Number of neurons in the first layer
    num_pixels:
        Number of pixels in the input
    Returns
    -------
    res : ElinaAbstract0Ptr
        Pointer to the new abstract object.

    """
    
    try:
        ffn_handle_first_sigmoid_layer_c = fppoly_api.ffn_handle_first_sigmoid_layer
        ffn_handle_first_sigmoid_layer_c.restype = None
        ffn_handle_first_sigmoid_layer_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr, _doublepp, ndpointer(ctypes.c_double),  c_size_t, c_size_t]
        ffn_handle_first_sigmoid_layer_c(man,element,weights, bias,  size, num_pixels)
    except Exception as inst:
        print('Problem with loading/calling "ffn_handle_first_sigmoid_layer" from "libfppoly.so"')
        print(inst)	
    
    return

def ffn_handle_first_tanh_layer(man, element,weights, bias,  size, num_pixels):
    """
    handle the first FFN Tanh layer
    
    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    weights : POINTER(POINTER(c_double))
        The weight matrix
    bias : POINTER(c_double)
        The bias vector
    size: c_size_t
	Number of neurons in the first layer
    num_pixels:
        Number of pixels in the input
    Returns
    -------
    res : ElinaAbstract0Ptr
        Pointer to the new abstract object.

    """
    
    try:
        ffn_handle_first_tanh_layer_c = fppoly_api.ffn_handle_first_tanh_layer
        ffn_handle_first_tanh_layer_c.restype = None
        ffn_handle_first_tanh_layer_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr, _doublepp, ndpointer(ctypes.c_double),  c_size_t, c_size_t]
        ffn_handle_first_tanh_layer_c(man,element,weights, bias,  size, num_pixels)
    except Exception as inst:
        print('Problem with loading/calling "ffn_handle_first_tanh_layer" from "libfppoly.so"')
        print(inst)	
    
    return




def ffn_handle_intermediate_relu_layer(man, element, weights, bias, num_out_neurons, num_in_neurons):
    """
    handle the intermediate FFN ReLU layer
    
    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    element : ElinaAbstract0Ptr
        Pointer to the abstract element
    weights: POINTER(POINTER(c_double))
        The weight matrix.
    bias: POINTER(c_size_t)
        The bias vector
    num_out_neurons: c_size_t
        number of output neurons
    num_in_neurons: c_size_t
	number of input neurons
    
    Returns
    -------
    None

    """

    try:
        ffn_handle_intermediate_relu_layer_c = fppoly_api.ffn_handle_intermediate_relu_layer
        ffn_handle_intermediate_relu_layer_c.restype = None
        ffn_handle_intermediate_relu_layer_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr, _doublepp, ndpointer(ctypes.c_double), c_size_t, c_size_t]
        ffn_handle_intermediate_relu_layer_c(man,element,weights,bias, num_out_neurons, num_in_neurons)
    except Exception as inst:
        print('Problem with loading/calling "affine_transform_handle_relu_layer" from "libfppoly.so"')
        print(inst)




def ffn_handle_intermediate_sigmoid_layer(man, element, weights, bias, num_out_neurons, num_in_neurons):
    """
    handle the intermediate FFN Sigmoid layer
    
    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    element : ElinaAbstract0Ptr
        Pointer to the abstract element
    weights: POINTER(POINTER(c_double))
        The weight matrix.
    bias: POINTER(c_size_t)
        The bias vector
    num_out_neurons: c_size_t
        number of output neurons
    num_in_neurons: c_size_t
	number of input neurons
    
    Returns
    -------
    None

    """
    
    try:
        ffn_handle_intermediate_sigmoid_layer_c = fppoly_api.ffn_handle_intermediate_sigmoid_layer
        ffn_handle_intermediate_sigmoid_layer_c.restype = None
        ffn_handle_intermediate_sigmoid_layer_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr, _doublepp, ndpointer(ctypes.c_double), c_size_t, c_size_t]
        ffn_handle_intermediate_sigmoid_layer_c(man,element,weights,bias, num_out_neurons, num_in_neurons)
    except Exception as inst:
        print('Problem with loading/calling "ffn_handle_sigmoid_layer" from "libfppoly.so"')
        print(inst)


def ffn_handle_intermediate_tanh_layer(man, element, weights, bias, num_out_neurons, num_in_neurons):
    """
    handle the intermediate FFN Tanh layer
    
    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    element : ElinaAbstract0Ptr
        Pointer to the abstract element
    weights: POINTER(POINTER(c_double))
        The weight matrix.
    bias: POINTER(c_size_t)
        The bias vector
    num_out_neurons: c_size_t
        number of output neurons
    num_in_neurons: c_size_t
	number of input neurons
    
    Returns
    -------
    None

    """
    
    try:
        ffn_handle_intermediate_tanh_layer_c = fppoly_api.ffn_handle_intermediate_tanh_layer
        ffn_handle_intermediate_tanh_layer_c.restype = None
        ffn_handle_intermediate_tanh_layer_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr, _doublepp, ndpointer(ctypes.c_double), c_size_t, c_size_t]
        ffn_handle_intermediate_tanh_layer_c(man,element,weights,bias, num_out_neurons, num_in_neurons)
    except Exception as inst:
        print('Problem with loading/calling "ffn_handle_tanh_layer" from "libfppoly.so"')
        print(inst)



def ffn_handle_last_relu_layer(man, element, weights, bias, num_out_neurons, num_in_neurons, has_relu):
    """
    handle the last FFN ReLU layer
    
    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    element : ElinaAbstract0Ptr
        Pointer to the ElinaAbstract0 abstract element
    weights: POINTER(POINTER(c_double))
        The weight matrix 
    bias : POINTER(c_size_t)
        The bias vector
    num_out_neurons: c_size_t
        The number of output neurons
    num_in_neurons: c_size_t
	The number of input_neurons
    has_relu: c_bool
        if the last layer has a ReLU activation
    Returns
    -------
    None

    """

    try:
        ffn_handle_last_relu_layer_c = fppoly_api.ffn_handle_last_relu_layer
        ffn_handle_last_relu_layer_c.restype = None
        ffn_handle_last_relu_layer_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr, _doublepp, ndpointer(ctypes.c_double), c_size_t, c_size_t, c_bool]
        ffn_handle_last_relu_layer_c(man,element,weights,bias, num_out_neurons, num_in_neurons, has_relu)
    except Exception as inst:
        print('Problem with loading/calling "ffn_handle_last_relu_layer" from "libfppoly.so"')
        print(inst)

def ffn_handle_last_sigmoid_layer(man, element, weights, bias, num_out_neurons, num_in_neurons, has_sigmoid):
    """
    handle the last FFN ReLU layer
    
    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    element : ElinaAbstract0Ptr
        Pointer to the ElinaAbstract0 abstract element
    weights: POINTER(POINTER(c_double))
        The weight matrix 
    bias : POINTER(c_size_t)
        The bias vector
    num_out_neurons: c_size_t
        The number of output neurons
    num_in_neurons: c_size_t
	The number of input_neurons
    has_sigmoid: c_bool
        if the last layer has a Sigmoid activation
    Returns
    -------
    None

    """
    
    try:
        ffn_handle_last_sigmoid_layer_c = fppoly_api.ffn_handle_last_sigmoid_layer
        ffn_handle_last_sigmoid_layer_c.restype = None
        ffn_handle_last_sigmoid_layer_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr, _doublepp, ndpointer(ctypes.c_double), c_size_t, c_size_t, c_bool]
        ffn_handle_last_sigmoid_layer_c(man,element,weights,bias, num_out_neurons, num_in_neurons, has_sigmoid)
    except Exception as inst:
        print('Problem with loading/calling "ffn_handle_last_sigmoid_layer" from "libfppoly.so"')
        print(inst)


def ffn_handle_last_tanh_layer(man, element, weights, bias, num_out_neurons, num_in_neurons, has_tanh):
    """
    handle the last FFN ReLU layer
    
    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    element : ElinaAbstract0Ptr
        Pointer to the ElinaAbstract0 abstract element
    weights: POINTER(POINTER(c_double))
        The weight matrix 
    bias : POINTER(c_size_t)
        The bias vector
    num_out_neurons: c_size_t
        The number of output neurons
    num_in_neurons: c_size_t
	The number of input_neurons
    has_tanh: c_bool
        if the last layer has a Tanh activation
    Returns
    -------
    None

    """
    
    try:
        ffn_handle_last_tanh_layer_c = fppoly_api.ffn_handle_last_tanh_layer
        ffn_handle_last_tanh_layer_c.restype = None
        ffn_handle_last_tanh_layer_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr, _doublepp, ndpointer(ctypes.c_double), c_size_t, c_size_t, c_bool]
        ffn_handle_last_tanh_layer_c(man,element,weights,bias, num_out_neurons, num_in_neurons, has_tanh)
    except Exception as inst:
        print('Problem with loading/calling "ffn_handle_last_tanh_layer" from "libfppoly.so"')
        print(inst)


def is_greater(man, element, y, x):
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
        is_greater_c = fppoly_api.is_greater
        is_greater_c.restype = c_bool
        is_greater_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr, ElinaDim, ElinaDim]
        res = is_greater_c(man,element,y, x)
    except Exception as inst:
        print('Problem with loading/calling "is_greater" from "libfppoly.so"')
        print(inst)
    return res

def conv_handle_first_layer(man, element, filter_weights, filter_bias,  input_size, filter_size, num_filters, strides, is_valid_padding, has_bias):
    """
    Convolutional Matrix multiplication in the first layer
    
    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    element : ElinaAbstract0Ptr
        Pointer to the ElinaAbstract0 abstract element.
    filter_weights: POINTER(double)
        filter weights
    filter_bias: POINTER(double)
        filter biases
    input_size: POINTER(c_size_t)
        size of the input
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
    None

    """
    try:
        conv_handle_first_layer_c = fppoly_api.conv_handle_first_layer
        conv_handle_first_layer_c.restype = None
        conv_handle_first_layer_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr, ndpointer(ctypes.c_double), ndpointer(ctypes.c_double), ndpointer(ctypes.c_size_t), POINTER(c_size_t), c_size_t, POINTER(c_size_t), c_bool, c_bool]
        conv_handle_first_layer_c(man,element, filter_weights, filter_bias, input_size, filter_size, num_filters, strides, is_valid_padding, has_bias)
    except Exception as inst:
        print('Problem with loading/calling "conv_handle_first_layer" from "libfppoly.so"')
        print(inst)
    return

def conv_handle_intermediate_relu_layer(man, element, filter_weights, filter_bias, input_size, filter_size, num_filters, strides, is_valid_padding, has_bias):
    """
    Convolutional Matrix multiplication in an Intermediate layer
    
    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    element : ElinaAbstract0Ptr
        Pointer to the ElinaAbstract0 abstract element.
    filter_weights: POINTER(double)
        filter weights
    filter_bias: POINTER(double)
        filter biases
    input_size: POINTER(c_size_t)
        size of the input
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
    None

    """
    try:
        conv_handle_intermediate_relu_layer_c = fppoly_api.conv_handle_intermediate_relu_layer
        conv_handle_intermediate_relu_layer_c.restype = None
        conv_handle_intermediate_relu_layer_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr, ndpointer(ctypes.c_double), ndpointer(ctypes.c_double), ndpointer(ctypes.c_size_t), POINTER(c_size_t), c_size_t, POINTER(c_size_t), c_bool, c_bool]
        conv_handle_intermediate_relu_layer_c(man, element, filter_weights, filter_bias, input_size, filter_size, num_filters, strides, is_valid_padding, has_bias)
    except Exception as inst:
        print('Problem with loading/calling "conv_handle_intermediate_relu_layer" from "libfppoly.so"')
        print(inst)


def handle_maxpool_layer(man, element, pool_size, input_size):
    """
    handle the Maxpool layer
    
    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    element : ElinaAbstract0Ptr
        Pointer to the ElinaAbstract0 abstract element.
    pool_size: POINTER(c_size_t)
        The size of the Maxpool filter 
    input_size : POINTER(c_size_t)
        The number of variables on which Maxpool will be applied.
    
    Returns
    -------
    res : c_size_t
        Number of neurons in the last layer

    """
    res=None
    try:
        handle_maxpool_layer_c = fppoly_api.handle_maxpool_layer
        handle_maxpool_layer_c.restype = c_size_t
        handle_maxpool_layer_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr, ndpointer(ctypes.c_size_t), ndpointer(ctypes.c_size_t)]
        res = handle_maxpool_layer_c(man, element, pool_size, input_size)
    except Exception as inst:
        print('Problem with loading/calling "handle_maxpool_layer" from "libfppoly.so"')
        print(inst)
    return res


def box_for_neuron(element,layerno, neuron_no):
    """
    returns bounds for a neuron in a layer
    
    Parameters
    ----------
    element : ElinaAbstract0Ptr
        Pointer to the ElinaAbstract0.
    layerno: c_size_t
        the layer number
    neuron_no: c_size_t
        the neuron number in the layer
    Returns
    -------
    interval_array : ElinaIntervalPtr
        ElinaIntervalArray representing the hypercube.

    """

    interval = None
    try:
        box_for_neuron_c = fppoly_api.box_for_neuron
        box_for_neuron_c.restype = ElinaIntervalPtr
        box_for_neuron_c.argtypes = [ElinaAbstract0Ptr, c_size_t, c_size_t]
        interval = box_for_neuron_c(element,layerno, neuron_no)
    except:
        print('Problem with loading/calling "box_for_neuron" from "fppoly.so"')
        print('Make sure you are passing ElinaAbstract0Ptr, c_size_t, c_size_t to the function')

    return interval

def box_for_layer(element,layerno):
    """
    returns bounds for all neurons in a layer
    
    Parameters
    ----------
    element : ElinaAbstract0Ptr
        Pointer to the ElinaAbstract0.
    layerno: c_size_t
        the layer number
    Returns
    -------
    interval_array : ElinaIntervalArray
        ElinaIntervalArray representing the hypercube.

    """

    interval_array = None
    try:
        box_for_layer_c = fppoly_api.box_for_layer
        box_for_layer_c.restype = ElinaIntervalArray
        box_for_layer_c.argtypes = [ElinaAbstract0Ptr, c_size_t]
        interval_array = box_for_layer_c(element,layerno)
    except:
        print('Problem with loading/calling "box_for_layer" from "fppoly.so"')
        print('Make sure you are passing ElinaAbstract0Ptr, c_size_t to the function')

    return interval_array
