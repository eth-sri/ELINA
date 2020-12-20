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


from fppoly_imports import *
from elina_manager_h import *
from elina_abstract0_h import *
from elina_interval_h import *
from elina_linexpr0_h import *
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

def fppoly_set_network_input_box(man, element, intdim, realdim, inf_array, sup_array):
    """
        Create an abstract element from perturbed input
        
        Parameters
        ----------
        man : ElinaManagerPtr
        Pointer to the ElinaManager.
        element: ElinaAbstract0Ptr
        Pointer to the abstract object
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
        fppoly_set_network_input_box_c = fppoly_api.fppoly_set_network_input_box
        fppoly_set_network_input_box_c.restype = None
        fppoly_set_network_input_box_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr, c_size_t, c_size_t,ndpointer(ctypes.c_double),ndpointer(ctypes.c_double)]
        res = fppoly_set_network_input_box_c(man,element, intdim, realdim, inf_array,sup_array)
    except Exception as inst:
        print('Problem with loading/calling "fppoly_set_network_input_box" from "libfppoly.so"')
        print(inst)

    return res

def fppoly_from_network_input_poly(man, intdim, realdim, inf_array, sup_array,
        lexpr_weights, lexpr_cst, lexpr_dim, uexpr_weights, uexpr_cst,
        uexpr_dim, expr_size, spatial_indices, spatial_neighbors, spatial_size,
        spatial_gamma):
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
    spatial_indices: POINTER(c_size_t)
        vector field indices
    spatial_neighbors: POINTER(c_size_t)
        neighboring vector field indices
    spatial_size: c_size_t
        number of spatial constraints
    spatial_gamma: double
        flow constraint parameter
    Returns
    -------
    res: ElinaAbstract0Ptr
         Pointer to the new abstract object

    """

    res = None
    try:
        fppoly_from_network_input_poly_c = fppoly_api.fppoly_from_network_input_poly
        fppoly_from_network_input_poly_c.restype = ElinaAbstract0Ptr
        fppoly_from_network_input_poly_c.argtypes = [
            ElinaManagerPtr, c_size_t, c_size_t, ndpointer(ctypes.c_double),
            ndpointer(ctypes.c_double), ndpointer(ctypes.c_double),
            ndpointer(ctypes.c_double), ndpointer(ctypes.c_size_t),
            ndpointer(ctypes.c_double), ndpointer(ctypes.c_double),
            ndpointer(ctypes.c_size_t), c_size_t, ndpointer(ctypes.c_size_t),
            ndpointer(ctypes.c_size_t), c_size_t, c_double
        ]
        res = fppoly_from_network_input_poly_c(
            man, intdim, realdim, inf_array, sup_array, lexpr_weights,
            lexpr_cst, lexpr_dim, uexpr_weights, uexpr_cst, uexpr_dim,
            expr_size, spatial_indices, spatial_neighbors, spatial_size,
            spatial_gamma
        )

    except Exception as inst:
        print('Problem with loading/calling "fppoly_from_network_input_poly" from "libfppoly.so"')
        print(inst)	

    return res





def handle_fully_connected_layer(man, element,weights, bias,  size, num_pixels, predecessors, num_predecessors):
    """
    handle the first FFN ReLU layer
    
    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    element : ElinaAbstract0Ptr
        Pointer to the ElinaAbstract0.
    weights : POINTER(POINTER(c_double))
        The weight matrix.
    bias : POINTER(c_double)
        The bias vector
    size: c_size_t
	Number of neurons in the first layer
    num_pixels: c_size_t
        Number of pixels in the input
    predecessors: POINTER(c_size_t)
        the layers before the current layer
    num_predecessors: c_size_t
        the number of predecessors of the current layer
    Returns
    -------
    res : ElinaAbstract0Ptr
        Pointer to the new abstract object.

    """

    try:
        handle_fully_connected_layer_c = fppoly_api.handle_fully_connected_layer
        handle_fully_connected_layer_c.restype = None
        handle_fully_connected_layer_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr, _doublepp, ndpointer(ctypes.c_double),  c_size_t, c_size_t, POINTER(c_size_t), c_size_t]
        handle_fully_connected_layer_c(man,element,weights, bias,  size, num_pixels, predecessors, num_predecessors)
    except Exception as inst:
        print('Problem with loading/calling "handle_fully_connected_layer" from "libfppoly.so"')
        print(inst)	

    return


def handle_fully_connected_layer_no_alloc(man, element,weights, bias,  size, num_pixels, predecessors, num_predecessors):
    """
        handle the first FFN ReLU layer
        
        Parameters
        ----------
        man : ElinaManagerPtr
        Pointer to the ElinaManager.
        element : ElinaAbstract0Ptr
        Pointer to the ElinaAbstract0.
        weights : POINTER(POINTER(c_double))
        The weight matrix.
        bias : POINTER(c_double)
        The bias vector
        size: c_size_t
        Number of neurons in the first layer
        num_pixels: c_size_t
        Number of pixels in the input
        predecessors: POINTER(c_size_t)
        the layers before the current layer
        num_predecessors: c_size_t
        the number of predecessors of the current layer
        Returns
        -------
        res : ElinaAbstract0Ptr
        Pointer to the new abstract object.
        
        """
    
    try:
        handle_fully_connected_layer_no_alloc_c = fppoly_api.handle_fully_connected_layer_no_alloc
        handle_fully_connected_layer_no_alloc_c.restype = None
        handle_fully_connected_layer_no_alloc_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr, _doublepp, ndpointer(ctypes.c_double),  c_size_t, c_size_t, POINTER(c_size_t), c_size_t]
        handle_fully_connected_layer_no_alloc_c(man,element,weights, bias,  size, num_pixels, predecessors, num_predecessors)
    except Exception as inst:
        print('Problem with loading/calling "handle_fully_connected_layer_no_alloc" from "libfppoly.so"')
        print(inst)

    return

def handle_concatenation_layer(man, element, predecessors, num_predecessors, C):
    """
        handle the first FFN ReLU layer

        Parameters
        ----------
        man : ElinaManagerPtr
        Pointer to the ElinaManager.
        element : ElinaAbstract0Ptr
        Pointer to the ElinaAbstract0.
    predecessors: POINTER(c_size_t)
        the layers before the current layer
    num_predecessors: c_size_t
        the number of predecessors of the current layer
    C: POINTER(c_size_t)
        the number of channels in each predecessor
    Returns
    -------
    res : None

    """

    try:
        handle_concatenation_layer_c = fppoly_api.handle_concatenation_layer
        handle_concatenation_layer_c.restype = None
        handle_concatenation_layer_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr, POINTER(c_size_t), c_size_t, POINTER(c_size_t)]
        handle_concatenation_layer_c(man, element, predecessors, num_predecessors, C)
    except Exception as inst:
        print('Problem with loading/calling "handle_concatenation_layer" from "libfppoly.so"')
        print(inst)

    return

def handle_tiling_layer(man, element, predecessors, num_predecessors, repeat):
    """
        handle the first FFN ReLU layer

        Parameters
        ----------
        man : ElinaManagerPtr
        Pointer to the ElinaManager.
        element : ElinaAbstract0Ptr
        Pointer to the ElinaAbstract0.
    predecessors: POINTER(c_size_t)
        the layers before the current layer
    num_predecessors: c_size_t
        the number of predecessors of the current layer
    repeat: c_size_t
        the number of tiles
    Returns
    -------
    res : None

    """

    try:
        handle_tiling_layer_c = fppoly_api.handle_tiling_layer
        handle_tiling_layer_c.restype = None
        handle_tiling_layer_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr, POINTER(c_size_t), c_size_t, c_size_t]
        handle_tiling_layer_c(man, element, predecessors, num_predecessors, repeat)
    except Exception as inst:
        print('Problem with loading/calling "handle_tiling_layer" from "libfppoly.so"')
        print(inst)

    return


def handle_sub_layer(man, element, cst, is_minuend, size, predecessors, num_predecessors):
    """
    handle the first FFN log layer
    
    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    cst : POINTER(c_double)
        The cst vector
    is_minuend: c_bool
        whether the assignment is y=x-cst or y=cst-x
    size: c_size_t
	Number of neurons in the first layer
    predecessors: POINTER(c_size_t)
        the layers before the current layer
    num_predecessors:  c_size_t
        the number of predecessors of the current layer
    Returns
    -------
    res : None

    """
    
    try:
        handle_sub_layer_c = fppoly_api.handle_sub_layer
        handle_sub_layer_c.restype = None
        handle_sub_layer_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr, ndpointer(ctypes.c_double),  c_bool, c_size_t, POINTER(c_size_t), c_size_t]
        handle_sub_layer_c(man, element, cst, is_minuend, size, predecessors, num_predecessors)
    except Exception as inst:
        print('Problem with loading/calling "handle_sub_layer" from "libfppoly.so"')
        print(inst)	
    
    return

def handle_mul_layer(man, element, cst,  size, predecessors, num_predecessors):
    """
    handle the first FFN log layer
    
    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    cst : POINTER(c_double)
        The cst vector
    size: c_size_t
	Number of neurons in the first layer
    predecessors: POINTER(c_size_t)
        the layers before the current layer
    num_predecessors: c_size_t
        the number of predecessors of the current layer
    Returns
    -------
    res : None

    """
    
    try:
        handle_mul_layer_c = fppoly_api.handle_mul_layer
        handle_mul_layer_c.restype = None
        handle_mul_layer_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr, ndpointer(ctypes.c_double),  c_size_t, POINTER(c_size_t), c_size_t]
        handle_mul_layer_c(man, element, cst, size, predecessors, num_predecessors)
    except Exception as inst:
        print('Problem with loading/calling "handle_mul_layer" from "libfppoly.so"')
        print(inst)	
    
    return


def handle_relu_layer(man, element, num_neurons, predecessors, num_predecessors, use_default_heuristics):
    """
    handle ReLU layer
    
    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    element : ElinaAbstract0Ptr
        Pointer to the abstract element
    num_neurons: c_size_t
        number of neurons
    predecessors: POINTER(c_size_t)
        the layers before the current layer
    num_predecessors: c_size_t
        the number of predecessors of the current layer
    use_default_heuristics: c_bool
        whether to use area heuristic
    Returns
    -------
    None

    """

    try:
        handle_relu_layer_c = fppoly_api.handle_relu_layer
        handle_relu_layer_c.restype = None
        handle_relu_layer_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr, c_size_t, POINTER(c_size_t), c_size_t, c_bool]
        handle_relu_layer_c(man, element, num_neurons, predecessors, num_predecessors, use_default_heuristics)
    except Exception as inst:
        print('Problem with loading/calling "handle_relu_layer" from "libfppoly.so"')
        print(inst)
        

def handle_sign_layer(man, element, num_neurons, predecessors, num_predecessors):
    """
    handle Sign layer
    
    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    element : ElinaAbstract0Ptr
        Pointer to the abstract element
    num_neurons: c_size_t
        number of neurons
    predecessors: POINTER(c_size_t)
        the layers before the current layer
    num_predecessors: c_size_t
        the number of predecessors of the current layer
    Returns
    -------
    None

    """

    try:
        handle_sign_layer_c = fppoly_api.handle_sign_layer
        handle_sign_layer_c.restype = None
        handle_sign_layer_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr, c_size_t, POINTER(c_size_t), c_size_t]
        handle_sign_layer_c(man, element, num_neurons, predecessors, num_predecessors)
    except Exception as inst:
        print('Problem with loading/calling "handle_sign_layer" from "libfppoly.so"')
        print(inst)


def handle_sigmoid_layer(man, element, num_neurons, predecessors, num_predecessors):
    """
    handle Sigmoid layer
    
    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    element : ElinaAbstract0Ptr
        Pointer to the abstract element
    num_neurons: c_size_t
        number of neurons
    predecessors: POINTER(c_size_t)
        the layers before the current layer
    num_predecessors: c_size_t
        the number of predecessors of the current layer
    Returns
    -------
    None

    """

    try:
        handle_sigmoid_layer_c = fppoly_api.handle_sigmoid_layer
        handle_sigmoid_layer_c.restype = None
        handle_sigmoid_layer_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr, c_size_t, POINTER(c_size_t), c_size_t]
        handle_sigmoid_layer_c(man, element, num_neurons, predecessors, num_predecessors)
    except Exception as inst:
        print('Problem with loading/calling "handle_sigmoid_layer" from "libfppoly.so"')
        print(inst)
        
def handle_tanh_layer(man, element, num_neurons, predecessors, num_predecessors):
    """
    handle Tanh layer
    
    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    element : ElinaAbstract0Ptr
        Pointer to the abstract element
    num_neurons: c_size_t
        number of neurons
    predecessors: POINTER(c_size_t)
        the layers before the current layer
    num_predecessors: c_size_t
        the number of predecessors of the current layer
    Returns
    -------
    None

    """

    try:
        handle_tanh_layer_c = fppoly_api.handle_tanh_layer
        handle_tanh_layer_c.restype = None
        handle_tanh_layer_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr, c_size_t, POINTER(c_size_t), c_size_t]
        handle_tanh_layer_c(man, element, num_neurons, predecessors, num_predecessors)
    except Exception as inst:
        print('Problem with loading/calling "handle_tanh_layer" from "libfppoly.so"')
        print(inst)
        
        
def handle_parabola_layer(man, element, num_neurons, predecessors, num_predecessors):
    """
    handle Parabola layer
    
    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    element : ElinaAbstract0Ptr
        Pointer to the abstract element
    num_neurons: c_size_t
        number of neurons
    predecessors: POINTER(c_size_t)
        the layers before the current layer
    num_predecessors: c_size_t
        the number of predecessors of the current layer
    Returns
    -------
    None

    """

    try:
        handle_parabola_layer_c = fppoly_api.handle_parabola_layer
        handle_parabola_layer_c.restype = None
        handle_parabola_layer_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr, c_size_t, POINTER(c_size_t), c_size_t]
        handle_parabola_layer_c(man, element, num_neurons, predecessors, num_predecessors)
    except Exception as inst:
        print('Problem with loading/calling "handle_parabola_layer" from "libfppoly.so"')
        print(inst)
        
        
def handle_log_layer(man, element, num_neurons, predecessors, num_predecessors):
    """
    handle Log layer
    
    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    element : ElinaAbstract0Ptr
        Pointer to the abstract element
    num_neurons: c_size_t
        number of neurons
    predecessors: POINTER(c_size_t)
        the layers before the current layer
    num_predecessors: c_size_t
        the number of predecessors of the current layer
    Returns
    -------
    None

    """

    try:
        handle_log_layer_c = fppoly_api.handle_log_layer
        handle_log_layer_c.restype = None
        handle_log_layer_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr, c_size_t, POINTER(c_size_t), c_size_t]
        handle_log_layer_c(man, element, num_neurons, predecessors, num_predecessors)
    except Exception as inst:
        print('Problem with loading/calling "handle_log_layer" from "libfppoly.so"')
        print(inst)



def is_greater(man, element, y, x, use_area_heuristic):
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
    use_area_heuristic: c_bool
        whether to use area heuristic
    Returns
    -------
    res = boolean

    """
    res= None
    try:
        is_greater_c = fppoly_api.is_greater
        is_greater_c.restype = c_bool
        is_greater_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr, ElinaDim, ElinaDim, c_bool]
        res = is_greater_c(man,element,y, x, use_area_heuristic)
    except Exception as inst:
        print('Problem with loading/calling "is_greater" from "libfppoly.so"')
        print(inst)
    return res

def handle_convolutional_layer(man, element, filter_weights, filter_bias,  input_size, filter_size, num_filters, strides, output_size, pad_top, pad_left, has_bias, predecessors, num_predecessors):
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
    predecessors: POINTER(c_size_t)
        the layers before the current layer
    num_predecessors: c_size_t
        the number of predecessors of the current layer
    Returns
    -------
    None

    """
    try:
        handle_convolutional_layer_c = fppoly_api.handle_convolutional_layer
        handle_convolutional_layer_c.restype = None
        handle_convolutional_layer_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr, ndpointer(ctypes.c_double), ndpointer(ctypes.c_double), ndpointer(ctypes.c_size_t), POINTER(c_size_t), c_size_t, POINTER(c_size_t), POINTER(c_size_t),c_size_t, c_size_t, c_bool, POINTER(c_size_t), c_size_t]
        handle_convolutional_layer_c(man,element, filter_weights, filter_bias, input_size, filter_size, num_filters, strides, output_size, pad_top, pad_left, has_bias, predecessors, num_predecessors)
    except Exception as inst:
        print('Problem with loading/calling "handle_convolutional_layer" from "libfppoly.so"')
        print(inst)
    return


def handle_pool_layer(man, element, pool_size, input_size, strides, pad_top, pad_left, output_size, predecessors, num_predecessors, is_maxpool):
    """
    handle the pooling layer
    
    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    element : ElinaAbstract0Ptr
        Pointer to the ElinaAbstract0 abstract element.
    pool_size: POINTER(c_size_t)
        The size of the pooling filter
    input_size: POINTER(c_size_t)
        the size of the input
    strides: POINTER(c_size_t)
        Strides for the pooling layer
    pad_top: c_size_t
        padding at top
    pad_left: c_size_t
        padding at left
    output_size : POINTER(c_size_t)
        The size of the output
    predecessors: POINTER(c_size_t)
        the layers before the current layer
    num_predecessors: c_size_t
        the number of predecessors of the current layer
    is_maxpool : c_bool
        whether it is maxpool or averagepool    
    Returns
    -------
    res : c_size_t
        Number of neurons in the last layer

    """
    res=None
    try:
        handle_pool_layer_c = fppoly_api.handle_pool_layer
        handle_pool_layer_c.restype = c_size_t
        handle_pool_layer_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr, POINTER(c_size_t), POINTER(c_size_t), POINTER(c_size_t), c_size_t, c_size_t, POINTER(c_size_t), POINTER(c_size_t), c_size_t , c_bool]
        res = handle_pool_layer_c(man, element, pool_size, input_size, strides, pad_left, pad_top, output_size, predecessors, num_predecessors, is_maxpool)
    except Exception as inst:
        print('Problem with loading/calling "handle_pool_layer" from "libfppoly.so"')
        print(inst)
    return res

def handle_residual_layer(man, element, num_neurons, predecessors, num_predecessors):
    """
    handle the Residual layer
    
    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    element : ElinaAbstract0Ptr
        Pointer to the ElinaAbstract0 abstract element.
    num_neurons: c_size_t
        The number of neurons in the residual layer 
    predecessors: POINTER(c_size_t)
        the layers before the current layer
    num_predecessors: c_size_t
        the number of predecessors of the current layer
    Returns
    -------
    None

    """

    try:
        handle_residual_layer_c = fppoly_api.handle_residual_layer
        handle_residual_layer_c.restype = None
        handle_residual_layer_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr, c_size_t, POINTER(c_size_t), c_size_t]
        handle_residual_layer_c(man, element, num_neurons, predecessors, num_predecessors)
    except Exception as inst:
        print('Problem with loading/calling "handle_residual_layer" from "libfppoly.so"')
        print(inst)



def box_for_neuron(man, element,layerno, neuron_no):
    """
    returns bounds for a neuron in a layer
    
    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
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
        box_for_neuron_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr, c_size_t, c_size_t]
        interval = box_for_neuron_c(man, element,layerno, neuron_no)
    except:
        print('Problem with loading/calling "box_for_neuron" from "fppoly.so"')
        print('Make sure you are passing ElinaManagerPtr, ElinaAbstract0Ptr, c_size_t, c_size_t to the function')

    return interval

def box_for_layer(man, element,layerno):
    """
    returns bounds for all neurons in a layer
    
    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
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
        box_for_layer_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr, c_size_t]
        interval_array = box_for_layer_c(man, element,layerno)
    except:
        print('Problem with loading/calling "box_for_layer" from "fppoly.so"')
        print('Make sure you are passing ElinaManagerPtr, ElinaAbstract0Ptr, c_size_t to the function')

    return interval_array

def get_num_neurons_in_layer(man, element,layerno):
    """
    returns the number of neurons in a layer
    
    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    element : ElinaAbstract0Ptr
        Pointer to the ElinaAbstract0.
    layerno: c_size_t
        the layer number
    Returns
    -------
    interval_array : ElinaIntervalArray
        ElinaIntervalArray representing the hypercube.

    """

    res = 0
    try:
        get_num_neurons_in_layer_c = fppoly_api.get_num_neurons_in_layer
        get_num_neurons_in_layer_c.restype = c_size_t
        get_num_neurons_in_layer_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr, c_size_t]
        res = get_num_neurons_in_layer_c(man, element,layerno)
    except:
        print('Problem with loading/calling "get_num_neurons_in_layer" from "fppoly.so"')
        print('Make sure you are passing ElinaManagerPtr, ElinaAbstract0Ptr, c_size_t to the function')

    return res


def update_bounds_for_neuron(man, element,layerno, neuron_no, lb, ub):
    """
    returns bounds for a neuron in a layer
    
    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    element : ElinaAbstract0Ptr
        Pointer to the ElinaAbstract0.
    layerno: c_size_t
        the layer number
    neuron_no: c_size_t
        the neuron number in the layer
    lb: c_double
        the updated lower bound
    ub: c_double
        the updated upper bound
    Returns
    -------
    None

    """

    
    try:
        update_bounds_for_neuron_c = fppoly_api.update_bounds_for_neuron
        update_bounds_for_neuron_c.restype = None
        update_bounds_for_neuron_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr, c_size_t, c_size_t, c_double, c_double]
        update_bounds_for_neuron_c(man, element,layerno, neuron_no, lb, ub)
    except:
        print('Problem with loading/calling "update_bounds_for_neuron" from "fppoly.so"')
        print('Make sure you are passing ElinaManagerPtr, ElinaAbstract0Ptr, c_size_t, c_size_t, c_double, c_double to the function')



def get_upper_bound_for_linexpr0(man,element,linexpr0, size, layerno):
    """
    returns bounds for a linexpr0 over neurons in "layerno"
    
    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    element : ElinaAbstract0Ptr
        Pointer to the ElinaAbstract0.
    linexpr0 : POINTER(ElinaLinexpr0Ptr)
        Pointer to the Elinalinexpr0
    size: c_size_t
        Size of the linexpr0 array
    layerno: c_size_t
        the layer number
    Returns
    -------
    ub : POINTER(c_double)
        array of upper bounds

    """

    res = None
    try:
        get_upper_bound_for_linexpr0_c = fppoly_api.get_upper_bound_for_linexpr0
        get_upper_bound_for_linexpr0_c.restype = POINTER(c_double)
        get_upper_bound_for_linexpr0_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr, ElinaLinexpr0Array, c_size_t, c_size_t]
        res = get_upper_bound_for_linexpr0_c(man, element, linexpr0, size, layerno)
    except:
        print('Problem with loading/calling "get_upper_bound_for_linexpr0" from "fppoly.so"')
        print('Make sure you are passing ElinaManagerPtr, ElinaAbstract0Ptr, ElinaLinexpr0Ptr, c_size_t to the function')

    return res
    

def get_lexpr_for_output_neuron(man,element,i):
    """
    returns lower polyhedra constraint for the i-th output neuron in terms of the input neurons
    
    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    element : ElinaAbstract0Ptr
        Pointer to the ElinaAbstract0.
    i: c_size_t
        output neuron number

    Returns
    -------
    expr :     ElinaLinexpr0Ptr
        The lower polyhedra expression for the output neuron in terms of input parameters and pixels

    """

    linexpr0 = None
    try:
        get_lexpr_for_output_neuron_c = fppoly_api.get_lexpr_for_output_neuron
        get_lexpr_for_output_neuron_c.restype = ElinaLinexpr0Ptr
        get_lexpr_for_output_neuron_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr, c_size_t]
        linexpr0 = get_lexpr_for_output_neuron_c(man,element,i)
    except:
        print('Problem with loading/calling "get_lexpr_for_output_neuron" from "fppoly.so"')
        print('Make sure you are passing ElinaManagerPtr, ElinaAbstract0Ptr, c_size_t to the function')

    return linexpr0


def get_uexpr_for_output_neuron(man,element,i):
    """
    returns lower polyhedra constraint for the i-th output neuron in terms of the input neurons
    
    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    element : ElinaAbstract0Ptr
        Pointer to the ElinaAbstract0.
    i: c_size_t
        output neuron number

    Returns
    -------
    expr :     ElinaLinexpr0Ptr
        The upper polyhedra expression for the output neuron in terms of input parameters and pixels

    """

    linexpr0 = None
    try:
        get_uexpr_for_output_neuron_c = fppoly_api.get_uexpr_for_output_neuron
        get_uexpr_for_output_neuron_c.restype = ElinaLinexpr0Ptr
        get_uexpr_for_output_neuron_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr, c_size_t]
        linexpr0 = get_uexpr_for_output_neuron_c(man,element,i)
    except:
        print('Problem with loading/calling "get_uexpr_for_output_neuron" from "fppoly.so"')
        print('Make sure you are passing ElinaManagerPtr, ElinaAbstract0Ptr, c_size_t to the function')

    return linexpr0


def create_lstm_layer(man, element,h, predecessors, num_predecessors):
    """
    creates an lstm layer for the neural network, this should be called only once per each lstm layer

    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    element : ElinaAbstract0Ptr
        Pointer to the ElinaAbstract0.
    h: c_size_t
        size of h_t
    predecessors: POINTER(c_size_t)
        the layers before the current layer
    num_predecessors: c_size_t
        the number of predecessors of the current layer
    Returns
    --------
    None
    """
    try:
        create_lstm_layer_c = fppoly_api.create_lstm_layer
        create_lstm_layer_c.restype = None
        create_lstm_layer_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr, c_size_t, POINTER(c_size_t), c_size_t]
        create_lstm_layer_c(man,element,h, predecessors, num_predecessors)
    except:
        print('Problem with loading/calling "create_lstm_layer" from "fppoly.so"')
        print('Make sure you are passing ElinaManagerPtr, ElinaAbstract0Ptr, c_size_t to the function')

    return
   

def handle_lstm_layer(man, element, weights, bias, d, h, predecessors, num_predecessors):
    """
    computes the hidden states and output vectors of the lstm unit, to be called at each time step after creating an LSTM unit

    Parameters
    -----------    
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    element : ElinaAbstract0Ptr
        Pointer to the ElinaAbstract0.
    weights : POINTER(POINTER(c_double))
        The weight matrix of size 4*h \times d+h, with h rows each for f_t, i_t, o_t, and c_t in order, 
        columnwise the first d entries correspond to x_t and the remaining correspond to h_t
    bias : POINTER(c_double)
        The bias vector of size 4*h, in the same format as weights
    d: c_size_t
       size of x_t
    h: c_size_t
       size of h_t
    predecessors: POINTER(c_size_t)
        the layers before the current layer
    num_predecessors: c_size_t
        the number of predecessors of the current layer
    Returns
    --------
    None
    """
    try:
        handle_lstm_layer_c = fppoly_api.handle_lstm_layer
        handle_lstm_layer_c.restype = None
        handle_lstm_layer_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr, _doublepp, ndpointer(ctypes.c_double), c_size_t, c_size_t, POINTER(c_size_t), c_size_t]
        handle_lstm_layer_c(man,element,weights,bias,d,h, predecessors, num_predecessors)
    except:
        print('Problem with loading/calling "handle_lstm_layer" from "fppoly.so"')
        print('Make sure you are passing ElinaManagerPtr, ElinaAbstract0Ptr, _doublepp, ndpointer(ctypes.c_double), c_size_t, c_size_t, c_bool to the function')

    return


def free_non_lstm_layer_expr(man,element,layerno):
    """
        returns bounds for a linexpr0 over neurons in "layerno"
        
        Parameters
        ----------
        man : ElinaManagerPtr
        Pointer to the ElinaManager.
        element : ElinaAbstract0Ptr
        Pointer to the ElinaAbstract0.
        layerno: c_size_t
        the layer number
        Returns
        -------
        None
        
        """
    
    try:
        free_non_lstm_layer_expr_c = fppoly_api.free_non_lstm_layer_expr
        free_non_lstm_layer_expr_c.restype = None
        free_non_lstm_layer_expr_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr, c_size_t]
        free_non_lstm_layer_expr_c(man, element, layerno)
    except:
        print('Problem with loading/calling "free_non_lstm_layer_expr" from "fppoly.so"')
        print('Make sure you are passing ElinaManagerPtr, ElinaAbstract0Ptr, c_size_t to the function')



def update_activation_upper_bound_for_neuron(man, element,layerno, neuron_no, coeff, dim, size):
    """
    returns bounds for a neuron in a layer
    
    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    element : ElinaAbstract0Ptr
        Pointer to the ElinaAbstract0.
    layerno: c_size_t
        the layer number
    neuron_no: c_size_t
        the neuron number in the layer
    lb: c_double
        the updated lower bound
    ub: c_double
        the updated upper bound
    Returns
    -------
    None

    """

    
    try:
        update_activation_upper_bound_for_neuron_c = fppoly_api.update_activation_upper_bound_for_neuron
        update_activation_upper_bound_for_neuron_c.restype = None
        update_activation_upper_bound_for_neuron_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr, c_size_t, c_size_t, ndpointer(ctypes.c_double), ndpointer(ctypes.c_size_t), c_size_t]
        update_activation_upper_bound_for_neuron_c(man, element,layerno, neuron_no, coeff, dim, size)
    except Exception as inst:
        print(inst)
        print('Problem with loading/calling "update_activation_upper_bound_for_neuron" from "fppoly.so"')
        print('Make sure you are passing ElinaManagerPtr, ElinaAbstract0Ptr, c_size_t, c_size_t, POINTER(c_double), POINTER(c_size_t), c_size_t to the function')
        
        
def update_activation_lower_bound_for_neuron(man, element,layerno, neuron_no, coeff, dim, size):
    """
    update bounds for a neuron in a layer
    
    Parameters
    ----------
    man : ElinaManagerPtr
        Pointer to the ElinaManager.
    element : ElinaAbstract0Ptr
        Pointer to the ElinaAbstract0.
    layerno: c_size_t
        the layer number
    neuron_no: c_size_t
        the neuron number in the layer
    lb: c_double
        the updated lower bound
    ub: c_double
        the updated upper bound
    Returns
    -------
    None

    """

    
    try:
        update_activation_lower_bound_for_neuron_c = fppoly_api.update_activation_lower_bound_for_neuron
        update_activation_lower_bound_for_neuron_c.restype = None
        update_activation_lower_bound_for_neuron_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr, c_size_t, c_size_t, ndpointer(ctypes.c_double), ndpointer(ctypes.c_size_t), c_size_t]
        update_activation_lower_bound_for_neuron_c(man, element,layerno, neuron_no, coeff, dim, size)
    except Exception as inst:
        print(inst)
        print('Problem with loading/calling "update_activation_lower_bound_for_neuron" from "fppoly.so"')
        print('Make sure you are passing ElinaManagerPtr, ElinaAbstract0Ptr, c_size_t, c_size_t, POINTER(c_double), POINTER(c_size_t), c_size_t to the function')
