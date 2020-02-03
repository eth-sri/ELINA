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


import numpy as np
from numpy.ctypeslib import ndpointer
import ctypes

from fppoly_gpu_imports import *
from elina_manager_h import *
from elina_abstract0_h import *
from elina_interval_h import *
from elina_linexpr0_h import *

c_float_type = ctypes.c_double
float_type = np.double

_doublepp = ndpointer(dtype=np.uintp, ndim=1, flags='C')


def fppoly_manager_alloc(gpu_number, backstep_depth):

    """
    Allocates an ElinaManager.

    Returns
    -------
        man : ElinaManagerPtr
            Pointer to the newly allocated ElinaManager

    """

    man = None

    try:
        fppoly_manager_alloc_c = fppoly_gpu_api.fppoly_manager_alloc
        fppoly_manager_alloc_c.restype = ElinaManagerPtr
        fppoly_manager_alloc_c.argtypes = [c_size_t, c_size_t]
        man = fppoly_manager_alloc_c(gpu_number, backstep_depth)
    except:
        print('Problem with loading/calling "fppoly_manager_alloc" from "libfppoly_gpu.so".')

    return man


def fppoly_from_network_input(man, intdim, realdim, inf_array, sup_array):

    """
    Create an abstract element from perturbed input.

    Parameters
    ----------
        man : ElinaManagerPtr
            Pointer to the ElinaManager
        intdim : c_size_t
            Number of integer variables
        realdim: c_size_t
            Number of real variables
        inf_array: POINTER(c_float_type)
            lower bound array
        sup_array: POINTER(c_float_type)
            upper bound array

    Returns
    -------
        res: ElinaAbstract0Ptr
             Pointer to the new abstract object

    """

    res = None

    try:
        fppoly_from_network_input_c = fppoly_gpu_api.fppoly_from_network_input
        fppoly_from_network_input_c.restype = ElinaAbstract0Ptr
        fppoly_from_network_input_c.argtypes = [ElinaManagerPtr, c_size_t, c_size_t, ndpointer(c_float_type), ndpointer(c_float_type)]
        res = fppoly_from_network_input_c(man, intdim, realdim, inf_array, sup_array)
    except Exception as inst:
        print('Problem with loading/calling "fppoly_from_network_input" from "libfppoly_gpu.so".')
        print(inst)

    return res


def fppoly_set_network_input_box(man, element, intdim, realdim, inf_array, sup_array):

    """
    Create an abstract element from perturbed input.

    Parameters
    ----------
        man : ElinaManagerPtr
            Pointer to the ElinaManager
        element: ElinaAbstract0Ptr
            Pointer to the abstract object
        intdim : c_size_t
            Number of integer variables
        realdim: c_size_t
            Number of real variables
        inf_array: POINTER(c_float_type)
            lower bound array
        sup_array: POINTER(c_float_type)
            upper bound array

    Returns
    -------
        res: ElinaAbstract0Ptr
            Pointer to the new abstract object

    """

    res = None

    try:
        fppoly_set_network_input_box_c = fppoly_gpu_api.fppoly_set_network_input_box
        fppoly_set_network_input_box_c.restype = None
        fppoly_set_network_input_box_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr, c_size_t, c_size_t, ndpointer(c_float_type), ndpointer(c_float_type)]
        res = fppoly_set_network_input_box_c(man, element, intdim, realdim, inf_array, sup_array)
    except Exception as inst:
        print('Problem with loading/calling "fppoly_set_network_input_box" from "libfppoly_gpu.so".')
        print(inst)

    return res


def fppoly_from_network_input_poly(man, intdim, realdim, inf_array, sup_array, lexpr_weights, lexpr_cst, lexpr_dim, uexpr_weights, uexpr_cst, uexpr_dim, expr_size):

    """
    Create an abstract element from perturbed input.

    Parameters
    ----------
        man : ElinaManagerPtr
            Pointer to the ElinaManager
        intdim : c_size_t
            Number of integer variables
        realdim: c_size_t
            Number of real variables
        inf_array: POINTER(c_float_type)
            lower bound array
        sup_array: POINTER(c_float_type)
            upper bound array
        lexpr_weights: POINTER(c_float_type)
            coefficients of the lower polyhedra constraints
        lexpr_cst: POINTER(c_float_type)
            constants of the lower polyhedra constraints
        lexpr_dim: POINTER(c_size_t)
            the indexes of the variables in the lower polyhedra constraints
        uexpr_weights: POINTER(c_float_type)
            coefficients of the upper polyhedra constraints
        uexpr_cst: POINTER(c_float_type)
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
        fppoly_from_network_input_poly_c = fppoly_gpu_api.fppoly_from_network_input_poly
        fppoly_from_network_input_poly_c.restype = ElinaAbstract0Ptr
        fppoly_from_network_input_poly_c.argtypes = [ElinaManagerPtr, c_size_t, c_size_t, ndpointer(c_float_type), ndpointer(c_float_type), ndpointer(c_float_type), ndpointer(c_float_type), ndpointer(c_size_t), ndpointer(c_float_type), ndpointer(c_float_type), ndpointer(c_size_t), c_size_t]
        res = fppoly_from_network_input_poly_c(man, intdim, realdim, inf_array, sup_array, lexpr_weights, lexpr_cst, lexpr_dim, uexpr_weights, uexpr_cst, uexpr_dim, expr_size)
    except Exception as inst:
        print('Problem with loading/calling "fppoly_from_network_input_poly" from "libfppoly_gpu.so".')
        print(inst)

    return res


def ffn_handle_first_relu_layer(man, element, weights, bias, size, num_pixels, predecessors):

    """
    Handle the first FFN ReLU layer.

    Parameters
    ----------
        man : ElinaManagerPtr
            Pointer to the ElinaManager
        element : ElinaAbstract0Ptr
            Pointer to the ElinaAbstract0
        weights : POINTER(POINTER(c_float_type))
            The weight matrix.
        bias : POINTER(c_float_type)
            The bias vector
        size: c_size_t
            Number of neurons in the first layer
        num_pixels:
            Number of pixels in the input
        predecessors:
            the layers before the current layer

    Returns
    -------
        res : ElinaAbstract0Ptr
            Pointer to the new abstract object

    """

    try:
        ffn_handle_first_relu_layer_c = fppoly_gpu_api.ffn_handle_first_relu_layer
        ffn_handle_first_relu_layer_c.restype = None
        ffn_handle_first_relu_layer_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr, _doublepp, ndpointer(c_float_type), c_size_t, c_size_t, POINTER(c_size_t)]
        ffn_handle_first_relu_layer_c(man, element, weights, bias, size, num_pixels, predecessors)
    except Exception as inst:
        print('Problem with loading/calling "ffn_handle_first_relu_layer" from "libfppoly_gpu.so".')
        print(inst)

    return


def ffn_handle_first_relu_layer_no_alloc(man, element, weights, bias, size, num_pixels, predecessors):

    """
    Handle the first FFN ReLU layer.

    Parameters
    ----------
        man : ElinaManagerPtr
            Pointer to the ElinaManager
        element : ElinaAbstract0Ptr
            Pointer to the ElinaAbstract0
        weights : POINTER(POINTER(c_float_type))
            The weight matrix.
        bias : POINTER(c_float_type)
            The bias vector
        size: c_size_t
            Number of neurons in the first layer
        num_pixels:
            Number of pixels in the input
        predecessors:
            the layers before the current layer

    Returns
    -------
        res : ElinaAbstract0Ptr
            Pointer to the new abstract object

    """

    print('"ffn_handle_first_relu_layer_no_alloc" is not implemented in "libfppoly_gpu.so". It only exists in the CPU version of fppoly.')

    return


def ffn_handle_first_sigmoid_layer(man, element, weights, bias, size, num_pixels, predecessors):

    """
    Handle the FFN first Sigmoid layer

    Parameters
    ----------
        man : ElinaManagerPtr
            Pointer to the ElinaManager.
        weights : POINTER(POINTER(c_float_type))
            The weight matrix.
        bias : POINTER(c_float_type)
            The bias vector
        size: c_size_t
            Number of neurons in the first layer
        num_pixels:
            Number of pixels in the input
        predecessors:
            the layers before the current layer

    Returns
    -------
        res : ElinaAbstract0Ptr
            Pointer to the new abstract object.

    """

    print('"ffn_handle_first_sigmoid_layer" is not implemented in "libfppoly_gpu.so". It only exists in the CPU version of fppoly.')

    return


def ffn_handle_first_sigmoid_layer_no_alloc(man, element, weights, bias, size, num_pixels, predecessors):

    """
    Handle the FFN first Sigmoid layer

    Parameters
    ----------
        man : ElinaManagerPtr
            Pointer to the ElinaManager.
        weights : POINTER(POINTER(c_float_type))
            The weight matrix.
        bias : POINTER(c_float_type)
            The bias vector
        size: c_size_t
            Number of neurons in the first layer
        num_pixels:
            Number of pixels in the input
        predecessors:
            the layers before the current layer

    Returns
    -------
        res : ElinaAbstract0Ptr
            Pointer to the new abstract object.

    """

    print('"ffn_handle_first_sigmoid_layer_no_alloc" is not implemented in "libfppoly_gpu.so". It only exists in the CPU version of fppoly.')

    return


def ffn_handle_first_tanh_layer(man, element, weights, bias, size, num_pixels, predecessors):

    """
    Handle the first FFN Tanh layer

    Parameters
    ----------
        man : ElinaManagerPtr
            Pointer to the ElinaManager.
        weights : POINTER(POINTER(c_float_type))
            The weight matrix
        bias : POINTER(c_float_type)
            The bias vector
        size: c_size_t
            Number of neurons in the first layer
        num_pixels:
            Number of pixels in the input
        predecessors:
            the layers before the current layer

    Returns
    -------
        res : ElinaAbstract0Ptr
            Pointer to the new abstract object.

    """

    print('"ffn_handle_first_tanh_layer" is not implemented in "libfppoly_gpu.so". It only exists in the CPU version of fppoly.')

    return


def ffn_handle_first_tanh_layer_no_alloc(man, element, weights, bias, size, num_pixels, predecessors):

    """
    Handle the first FFN Tanh layer

    Parameters
    ----------
        man : ElinaManagerPtr
            Pointer to the ElinaManager.
        weights : POINTER(POINTER(c_float_type))
            The weight matrix
        bias : POINTER(c_float_type)
            The bias vector
        size: c_size_t
            Number of neurons in the first layer
        num_pixels:
            Number of pixels in the input
        predecessors:
            the layers before the current layer

    Returns
    -------
        res : ElinaAbstract0Ptr
            Pointer to the new abstract object.

    """

    print('"ffn_handle_first_tanh_layer_no_alloc" is not implemented in "libfppoly_gpu.so". It only exists in the CPU version of fppoly.')

    return


def ffn_handle_first_parabola_layer(man, element, weights, bias, size, num_pixels, predecessors):

    """
    Handle the first FFN Parabolic layer

    Parameters
    ----------
        man : ElinaManagerPtr
            Pointer to the ElinaManager.
        weights : POINTER(POINTER(c_float_type))
            The weight matrix
        bias : POINTER(c_float_type)
            The bias vector
        size: c_size_t
            Number of neurons in the first layer
        num_pixels:
            Number of pixels in the input
        predecessors:
            the layers before the current layer

    Returns
    -------
        res : ElinaAbstract0Ptr
            Pointer to the new abstract object.

    """

    print('"ffn_handle_first_parabola_layer" is not implemented in "libfppoly_gpu.so". It only exists in the CPU version of fppoly.')

    return


def ffn_handle_first_parabola_layer_no_alloc(man, element, weights, bias, size, num_pixels, predecessors):

    """
    Handle the first FFN Parabolic layer

    Parameters
    ----------
        man : ElinaManagerPtr
            Pointer to the ElinaManager.
        weights : POINTER(POINTER(c_float_type))
            The weight matrix
        bias : POINTER(c_float_type)
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

    print('"ffn_handle_first_parabola_layer_no_alloc" is not implemented in "libfppoly_gpu.so". It only exists in the CPU version of fppoly.')

    return


def ffn_handle_first_log_layer(man, element, weights, bias, size, num_pixels, predecessors):

    """
    Handle the first FFN log layer

    Parameters
    ----------
        man : ElinaManagerPtr
            Pointer to the ElinaManager.
        weights : POINTER(POINTER(c_float_type))
            The weight matrix
        bias : POINTER(c_float_type)
            The bias vector
        size: c_size_t
            Number of neurons in the first layer
        num_pixels:
            Number of pixels in the input
        predecessors:
            the layers before the current layer

    Returns
    -------
        res : ElinaAbstract0Ptr
            Pointer to the new abstract object.

    """

    print('"ffn_handle_first_log_layer" is not implemented in "libfppoly_gpu.so". It only exists in the CPU version of fppoly.')

    return


def ffn_handle_first_log_layer_no_alloc(man, element, weights, bias, size, num_pixels, predecessors):

    """
    Handle the first FFN log layer

    Parameters
    ----------
        man : ElinaManagerPtr
            Pointer to the ElinaManager.
        weights : POINTER(POINTER(c_float_type))
            The weight matrix
        bias : POINTER(c_float_type)
            The bias vector
        size: c_size_t
            Number of neurons in the first layer
        num_pixels:
            Number of pixels in the input
        predecessors:
            the layers before the current layer

    Returns
    -------
        res : ElinaAbstract0Ptr
            Pointer to the new abstract object.

    """

    print('"ffn_handle_first_log_layer_no_alloc" is not implemented in "libfppoly_gpu.so". It only exists in the CPU version of fppoly.')

    return


def ffn_handle_first_sub_layer(man, element, cst, is_minuend, size, predecessors):

    """
    Handle the first FFN sub layer

    Parameters
    ----------
        man : ElinaManagerPtr
            Pointer to the ElinaManager
        cst : POINTER(c_float_type)
            The cst vector
        is_minuend: c_bool
            Whether the assignment is y=x-cst or y=cst-x
        size: c_size_t
            Number of neurons in the first layer
        predecessors:
            The layers before the current layer

    Returns
    -------
        res : None

    """

    try:
        ffn_handle_first_sub_layer_c = fppoly_gpu_api.ffn_handle_first_sub_layer
        ffn_handle_first_sub_layer_c.restype = None
        ffn_handle_first_sub_layer_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr, ndpointer(c_float_type),  c_bool, c_size_t, POINTER(c_size_t)]
        ffn_handle_first_sub_layer_c(man, element, cst, is_minuend, size, predecessors)
    except Exception as inst:
        print('Problem with loading/calling "ffn_handle_first_sub_layer" from "libfppoly_gpu.so"')
        print(inst)

    return


def ffn_handle_first_mul_layer(man, element, cst, size, predecessors):

    """
    Handle the first FFN mul layer

    Parameters
    ----------
        man : ElinaManagerPtr
            Pointer to the ElinaManager
        cst : POINTER(c_float_type)
            The cst vector
        size: c_size_t
            Number of neurons in the first layer
        predecessors:
            The layers before the current layer

    Returns
    -------
        res : None

    """

    try:
        ffn_handle_first_mul_layer_c = fppoly_gpu_api.ffn_handle_first_mul_layer
        ffn_handle_first_mul_layer_c.restype = None
        ffn_handle_first_mul_layer_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr, ndpointer(c_float_type),  c_size_t, POINTER(c_size_t)]
        ffn_handle_first_mul_layer_c(man, element, cst, size, predecessors)
    except Exception as inst:
        print('Problem with loading/calling "ffn_handle_first_mul_layer" from "libfppoly_gpu.so"')
        print(inst)

    return


def ffn_handle_intermediate_affine_layer(man, element, weights, bias, num_out_neurons, num_in_neurons, predecessors, use_area_heuristic):

    """
    Handle the intermediate FFN ReLU layer.

    Parameters
    ----------
        man : ElinaManagerPtr
            Pointer to the ElinaManager
        element : ElinaAbstract0Ptr
            Pointer to the abstract element
        weights: POINTER(POINTER(c_float_type))
            The weight matrix
        bias: POINTER(c_size_t)
            The bias vector
        num_out_neurons: c_size_t
            number of output neurons
        num_in_neurons: c_size_t
            number of input neurons
        predecessors:
            the layers before the current layer
        use_area_heuristic: c_bool
            whether to use area heuristic

    Returns
    -------
        None

    """

    try:
        ffn_handle_intermediate_affine_layer_c = fppoly_gpu_api.ffn_handle_intermediate_affine_layer
        ffn_handle_intermediate_affine_layer_c.restype = None
        ffn_handle_intermediate_affine_layer_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr, _doublepp, ndpointer(c_float_type), c_size_t, c_size_t, POINTER(c_size_t), c_bool]
        ffn_handle_intermediate_affine_layer_c(man, element, weights, bias, num_out_neurons, num_in_neurons, predecessors, use_area_heuristic)
    except Exception as inst:
        print('Problem with loading/calling "ffn_handle_intermediate_affine_layer" from "libfppoly_gpu.so".')
        print(inst)

    return


def ffn_handle_intermediate_affine_layer_no_alloc(man, element, weights, bias, num_out_neurons, num_in_neurons, predecessors, use_area_heuristic):

    """
    Handle the intermediate FFN ReLU layer.

    Parameters
    ----------
        man : ElinaManagerPtr
            Pointer to the ElinaManager
        element : ElinaAbstract0Ptr
            Pointer to the abstract element
        weights: POINTER(POINTER(c_float_type))
            The weight matrix
        bias: POINTER(c_size_t)
            The bias vector
        num_out_neurons: c_size_t
            number of output neurons
        num_in_neurons: c_size_t
            number of input neurons
        predecessors:
            the layers before the current layer
        use_area_heuristic: c_bool
            whether to use area heuristic

    Returns
    -------
        None

    """

    print('"ffn_handle_intermediate_affine_layer_no_alloc" is not implemented in "libfppoly_gpu.so". It only exists in the CPU version of fppoly.')

    return


def ffn_handle_intermediate_relu_layer(man, element, weights, bias, num_out_neurons, num_in_neurons, predecessors, use_area_heuristic):

    """
    Handle the intermediate FFN ReLU layer.

    Parameters
    ----------
        man: ElinaManagerPtr
            Pointer to the ElinaManager
        element: ElinaAbstract0Ptr
            Pointer to the abstract element
        weights: POINTER(POINTER(c_float_type))
            The weight matrix
        bias: POINTER(c_size_t)
            The bias vector
        num_out_neurons: c_size_t
            number of output neurons
        num_in_neurons: c_size_t
            number of input neurons
        predecessors:
            the layers before the current layer
        use_area_heuristic: c_bool
            whether to use area heuristic

    Returns
    -------
        None

    """

    try:
        ffn_handle_intermediate_relu_layer_c = fppoly_gpu_api.ffn_handle_intermediate_relu_layer
        ffn_handle_intermediate_relu_layer_c.restype = None
        ffn_handle_intermediate_relu_layer_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr, _doublepp, ndpointer(c_float_type), c_size_t, c_size_t, POINTER(c_size_t), c_bool]
        ffn_handle_intermediate_relu_layer_c(man, element, weights, bias, num_out_neurons, num_in_neurons, predecessors, use_area_heuristic)
    except Exception as inst:
        print('Problem with loading/calling "ffn_handle_intermediate_relu_layer" from "libfppoly_gpu.so".')
        print(inst)

    return


def ffn_handle_intermediate_relu_layer_no_alloc(man, element, weights, bias, num_out_neurons, num_in_neurons, predecessors, use_area_heuristic):

    """
    Handle the intermediate FFN ReLU layer.

    Parameters
    ----------
        man: ElinaManagerPtr
            Pointer to the ElinaManager
        element: ElinaAbstract0Ptr
            Pointer to the abstract element
        weights: POINTER(POINTER(c_float_type))
            The weight matrix
        bias: POINTER(c_size_t)
            The bias vector
        num_out_neurons: c_size_t
            number of output neurons
        num_in_neurons: c_size_t
            number of input neurons
        predecessors:
            the layers before the current layer
        use_area_heuristic: c_bool
            whether to use area heuristic

    Returns
    -------
        None

    """

    print('"ffn_handle_intermediate_relu_layer_no_alloc" is not implemented in "libfppoly_gpu.so". It only exists in the CPU version of fppoly.')

    return


def ffn_handle_intermediate_sigmoid_layer(man, element, weights, bias, num_out_neurons, num_in_neurons, predecessors, use_area_heuristic):

    """
    Handle the intermediate FFN Sigmoid layer

    Parameters
    ----------
        man : ElinaManagerPtr
            Pointer to the ElinaManager.
        element : ElinaAbstract0Ptr
            Pointer to the abstract element
        weights: POINTER(POINTER(c_float_type))
            The weight matrix.
        bias: POINTER(c_size_t)
            The bias vector
        num_out_neurons: c_size_t
            number of output neurons
        num_in_neurons: c_size_t
            number of input neurons
        predecessors:
            the layers before the current layer
        use_area_heuristic: c_bool
            whether to use area heuristic

    Returns
    -------
        None

    """

    print('"ffn_handle_intermediate_sigmoid_layer" is not implemented in "libfppoly_gpu.so". It only exists in the CPU version of fppoly.')

    return


def ffn_handle_intermediate_sigmoid_layer_no_alloc(man, element, weights, bias, num_out_neurons, num_in_neurons, predecessors, use_area_heuristic):

    """
    Handle the intermediate FFN Sigmoid layer

    Parameters
    ----------
        man : ElinaManagerPtr
            Pointer to the ElinaManager.
        element : ElinaAbstract0Ptr
            Pointer to the abstract element
        weights: POINTER(POINTER(c_float_type))
            The weight matrix.
        bias: POINTER(c_size_t)
            The bias vector
        num_out_neurons: c_size_t
            number of output neurons
        num_in_neurons: c_size_t
            number of input neurons
        predecessors:
            the layers before the current layer
        use_area_heuristic: c_bool
            whether to use area heuristic

        Returns
        -------
            None

    """

    print('"ffn_handle_intermediate_sigmoid_layer_no_alloc" is not implemented in "libfppoly_gpu.so". It only exists in the CPU version of fppoly.')

    return


def ffn_handle_intermediate_tanh_layer(man, element, weights, bias, num_out_neurons, num_in_neurons, predecessors, use_area_heuristic):

    """
    Handle the intermediate FFN Tanh layer

    Parameters
    ----------
        man : ElinaManagerPtr
            Pointer to the ElinaManager.
        element : ElinaAbstract0Ptr
            Pointer to the abstract element
        weights: POINTER(POINTER(c_float_type))
            The weight matrix.
        bias: POINTER(c_size_t)
            The bias vector
        num_out_neurons: c_size_t
            number of output neurons
        num_in_neurons: c_size_t
            number of input neurons
        predecessors:
            the layers before the current layer
        use_area_heuristic: c_bool
            whether to use area heuristic
    Returns
    -------
        None

    """

    print('"ffn_handle_intermediate_tanh_layer" is not implemented in "libfppoly_gpu.so". It only exists in the CPU version of fppoly.')

    return


def ffn_handle_intermediate_tanh_layer_no_alloc(man, element, weights, bias, num_out_neurons, num_in_neurons, predecessors, use_area_heuristic):

    """
    Handle the intermediate FFN Tanh layer

    Parameters
    ----------
        man : ElinaManagerPtr
            Pointer to the ElinaManager.
        element : ElinaAbstract0Ptr
            Pointer to the abstract element
        weights: POINTER(POINTER(c_float_type))
            The weight matrix.
        bias: POINTER(c_size_t)
            The bias vector
        num_out_neurons: c_size_t
            number of output neurons
        num_in_neurons: c_size_t
            number of input neurons
        predecessors:
            the layers before the current layer
        use_area_heuristic: c_bool
             whether to use area heuristic

    Returns
    -------
        None

    """

    print('"ffn_handle_intermediate_tanh_layer_no_alloc" is not implemented in "libfppoly_gpu.so". It only exists in the CPU version of fppoly.')

    return


def ffn_handle_intermediate_parabola_layer(man, element, weights, bias, num_out_neurons, num_in_neurons, predecessors, use_area_heuristic):

    """
    Handle the intermediate FFN Parabolic layer

    Parameters
    ----------
        man : ElinaManagerPtr
            Pointer to the ElinaManager.
        element : ElinaAbstract0Ptr
            Pointer to the abstract element
        weights: POINTER(POINTER(c_float_type))
            The weight matrix.
        bias: POINTER(c_size_t)
            The bias vector
        num_out_neurons: c_size_t
            number of output neurons
        num_in_neurons: c_size_t
            number of input neurons
        predecessors:
            the layers before the current layer
        use_area_heuristic: c_bool
            whether to use area heuristic

    Returns
    -------
        None

    """

    print('"ffn_handle_intermediate_parabola_layer" is not implemented in "libfppoly_gpu.so". It only exists in the CPU version of fppoly.')

    return


def ffn_handle_intermediate_parabola_layer_no_alloc(man, element, weights, bias, num_out_neurons, num_in_neurons, predecessors, use_area_heuristic):
    """
    Handle the intermediate FFN Parabolic layer

    Parameters
    ----------
        man : ElinaManagerPtr
            Pointer to the ElinaManager.
        element : ElinaAbstract0Ptr
            Pointer to the abstract element
        weights: POINTER(POINTER(c_float_type))
            The weight matrix.
        bias: POINTER(c_size_t)
            The bias vector
        num_out_neurons: c_size_t
            number of output neurons
        num_in_neurons: c_size_t
            number of input neurons
        predecessors:
            the layers before the current layer
        use_area_heuristic: c_bool
            whether to use area heuristic

    Returns
    -------
        None

    """

    print('"ffn_handle_intermediate_parabola_layer_no_alloc" is not implemented in "libfppoly_gpu.so". It only exists in the CPU version of fppoly.')

    return


def ffn_handle_intermediate_log_layer(man, element, weights, bias, num_out_neurons, num_in_neurons, predecessors, use_area_heuristic):

    """
    Handle the intermediate FFN Log layer

    Parameters
    ----------
        man : ElinaManagerPtr
            Pointer to the ElinaManager.
        element : ElinaAbstract0Ptr
            Pointer to the abstract element
        weights: POINTER(POINTER(c_float_type))
            The weight matrix.
        bias: POINTER(c_size_t)
            The bias vector
        num_out_neurons: c_size_t
            number of output neurons
        num_in_neurons: c_size_t
            number of input neurons
        predecessors:
            the layers before the current layer
        use_area_heuristic: c_bool
            whether to use area heuristic

    Returns
    -------
    None

    """

    print('"ffn_handle_intermediate_log_layer" is not implemented in "libfppoly_gpu.so". It only exists in the CPU version of fppoly.')

    return


def ffn_handle_intermediate_log_layer_no_alloc(man, element, weights, bias, num_out_neurons, num_in_neurons, predecessors, use_area_heuristic):

    """
    Handle the intermediate FFN Log layer

    Parameters
    ----------
        man : ElinaManagerPtr
            Pointer to the ElinaManager.
        element : ElinaAbstract0Ptr
            Pointer to the abstract element
        weights: POINTER(POINTER(c_float_type))
            The weight matrix.
        bias: POINTER(c_size_t)
            The bias vector
        num_out_neurons: c_size_t
            number of output neurons
        num_in_neurons: c_size_t
            number of input neurons
        predecessors:
            the layers before the current layer
        use_area_heuristic: c_bool
            whether to use area heuristic

    Returns
    -------
        None

    """

    print('"ffn_handle_intermediate_log_layer_no_alloc" is not implemented in "libfppoly_gpu.so". It only exists in the CPU version of fppoly.')

    return


def ffn_handle_intermediate_sub_layer(man, element, bias, is_minuend, num_in_neurons, predecessors, use_area_heuristic):

    """
    Handle the intermediate FFN Sub layer

    Parameters
    ----------
        man : ElinaManagerPtr
            Pointer to the ElinaManager
        element : ElinaAbstract0Ptr
            Pointer to the abstract element
        bias: POINTER(c_size_t)
            The bias vector
        is_minuend: c_bool
            Whether it is y = x-b pr b-x
        num_in_neurons: c_size_t
            Number of input neurons
        predecessors:
            The layers before the current layer
        use_area_heuristic: c_bool
            Whether to use area heuristic

    Returns
    -------
        None

    """

    try:
        ffn_handle_intermediate_sub_layer_c = fppoly_gpu_api.ffn_handle_intermediate_sub_layer
        ffn_handle_intermediate_sub_layer_c.restype = None
        ffn_handle_intermediate_sub_layer_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr, ndpointer(c_float_type), c_bool, c_size_t, POINTER(c_size_t),c_bool]
        ffn_handle_intermediate_sub_layer_c(man,element, bias, is_minuend, num_in_neurons, predecessors, use_area_heuristic)
    except Exception as inst:
        print('Problem with loading/calling "ffn_handle_intermediate_sub_layer" from "libfppoly_gpu.so"')
        print(inst)

    return


def ffn_handle_intermediate_mul_layer(man, element, bias, num_in_neurons, predecessors, use_area_heuristic):

    """
    Handle the intermediate FFN Mul layer

    Parameters
    ----------
        man : ElinaManagerPtr
            Pointer to the ElinaManager
        element : ElinaAbstract0Ptr
            Pointer to the abstract element
        bias: POINTER(c_size_t)
            The bias vector
        num_in_neurons: c_size_t
            Number of input neurons
        predecessors:
            The layers before the current layer
        use_area_heuristic: c_bool
            Whether to use area heuristic

    Returns
    -------
        None

    """

    try:
        ffn_handle_intermediate_mul_layer_c = fppoly_gpu_api.ffn_handle_intermediate_mul_layer
        ffn_handle_intermediate_mul_layer_c.restype = None
        ffn_handle_intermediate_mul_layer_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr, ndpointer(c_float_type), c_size_t, POINTER(c_size_t),c_bool]
        ffn_handle_intermediate_mul_layer_c(man,element, bias, num_in_neurons, predecessors, use_area_heuristic)
    except Exception as inst:
        print('Problem with loading/calling "ffn_handle_intermediate_mul_layer" from "libfppoly_gpu.so"')
        print(inst)

    return


def ffn_handle_last_relu_layer(man, element, weights, bias, num_out_neurons, num_in_neurons, predecessors, has_relu, use_area_heuristic):

    """
    Handle the last FFN ReLU layer.

    Parameters
    ----------
        man: ElinaManagerPtr
            Pointer to the ElinaManager
        element: ElinaAbstract0Ptr
            Pointer to the ElinaAbstract0 abstract element
        weights: POINTER(POINTER(c_float_type))
            The weight matrix
        bias: POINTER(c_size_t)
            The bias vector
        num_out_neurons: c_size_t
            The number of output neurons
        num_in_neurons: c_size_t
            The number of input_neurons
        predecessors:
            the layers before the current layer
        has_relu: c_bool
            if the last layer has a ReLU activation
        use_area_heuristic: c_bool
            whether to use area heuristic

    Returns
    -------
        None

    """

    try:
        ffn_handle_last_relu_layer_c = fppoly_gpu_api.ffn_handle_last_relu_layer
        ffn_handle_last_relu_layer_c.restype = None
        ffn_handle_last_relu_layer_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr, _doublepp, ndpointer(c_float_type), c_size_t, c_size_t, POINTER(c_size_t), c_bool, c_bool]
        ffn_handle_last_relu_layer_c(man, element, weights, bias, num_out_neurons, num_in_neurons, predecessors, has_relu, use_area_heuristic)
    except Exception as inst:
        print('Problem with loading/calling "ffn_handle_last_relu_layer" from "libfppoly_gpu.so".')
        print(inst)

    return


def ffn_handle_last_relu_layer_no_alloc(man, element, weights, bias, num_out_neurons, num_in_neurons, predecessors, has_relu, use_area_heuristic):

    """
    Handle the last FFN ReLU layer.

    Parameters
    ----------
        man: ElinaManagerPtr
            Pointer to the ElinaManager
        element: ElinaAbstract0Ptr
            Pointer to the ElinaAbstract0 abstract element
        weights: POINTER(POINTER(c_float_type))
            The weight matrix
        bias: POINTER(c_size_t)
            The bias vector
        num_out_neurons: c_size_t
            The number of output neurons
        num_in_neurons: c_size_t
            The number of input_neurons
        has_relu: c_bool
            if the last layer has a ReLU activation
        predecessors:
            the layers before the current layer
        use_area_heuristic: c_bool
            whether to use area heuristic

    Returns
    -------
        None

    """

    print('"ffn_handle_last_relu_layer_no_alloc" is not implemented in "libfppoly_gpu.so". It only exists in the CPU version of fppoly.')

    return


def ffn_handle_last_sigmoid_layer(man, element, weights, bias, num_out_neurons, num_in_neurons, predecessors, has_sigmoid, use_area_heuristic):

    """
    Handle the last FFN Sigmoid layer

    Parameters
    ----------
        man : ElinaManagerPtr
            Pointer to the ElinaManager.
        element : ElinaAbstract0Ptr
            Pointer to the ElinaAbstract0 abstract element
        weights: POINTER(POINTER(c_float_type))
            The weight matrix
        bias : POINTER(c_size_t)
            The bias vector
        num_out_neurons: c_size_t
            The number of output neurons
        num_in_neurons: c_size_t
            The number of input_neurons
        predecessors:
            the layers before the current layer
        has_sigmoid: c_bool
            if the last layer has a Sigmoid activation
        use_area_heuristic: c_bool
            whether to use area heuristic

    Returns
    -------
    None

    """

    print('"ffn_handle_last_sigmoid_layer" is not implemented in "libfppoly_gpu.so". It only exists in the CPU version of fppoly.')

    return


def ffn_handle_last_sigmoid_layer_no_alloc(man, element, weights, bias, num_out_neurons, num_in_neurons, predecessors, has_sigmoid, use_area_heuristic):

    """
    Handle the last FFN Sigmoid layer

    Parameters
    ----------
        man : ElinaManagerPtr
            Pointer to the ElinaManager.
        element : ElinaAbstract0Ptr
            Pointer to the ElinaAbstract0 abstract element
        weights: POINTER(POINTER(c_float_type))
            The weight matrix
        bias : POINTER(c_size_t)
            The bias vector
        num_out_neurons: c_size_t
            The number of output neurons
        num_in_neurons: c_size_t
            The number of input_neurons
        predecessors:
            the layers before the current layer
        has_sigmoid: c_bool
            if the last layer has a Sigmoid activation
        use_area_heuristic: c_bool
            whether to use area heuristic

    Returns
    -------
        None

    """

    print('"ffn_handle_last_sigmoid_layer_no_alloc" is not implemented in "libfppoly_gpu.so". It only exists in the CPU version of fppoly.')

    return


def ffn_handle_last_tanh_layer(man, element, weights, bias, num_out_neurons, num_in_neurons, predecessors, has_tanh, use_area_heuristic):

    """
    Handle the last FFN Tanh layer

    Parameters
    ----------
        man : ElinaManagerPtr
            Pointer to the ElinaManager.
        element : ElinaAbstract0Ptr
            Pointer to the ElinaAbstract0 abstract element
        weights: POINTER(POINTER(c_float_type))
            The weight matrix
        bias : POINTER(c_size_t)
            The bias vector
        num_out_neurons: c_size_t
            The number of output neurons
        num_in_neurons: c_size_t
            The number of input_neurons
        predecessors:
            the layers before the current layer
        has_tanh: c_bool
            if the last layer has a Tanh activation
        use_area_heuristic: c_bool
            whether to use area heuristic

    Returns
    -------
        None

    """

    print('"ffn_handle_last_tanh_layer" is not implemented in "libfppoly_gpu.so". It only exists in the CPU version of fppoly.')

    return


def ffn_handle_last_tanh_layer_no_alloc(man, element, weights, bias, num_out_neurons, num_in_neurons, predecessors, has_tanh, use_area_heuristic):

    """
    Handle the last FFN Tanh layer

    Parameters
    ----------
        man : ElinaManagerPtr
            Pointer to the ElinaManager.
        element : ElinaAbstract0Ptr
            Pointer to the ElinaAbstract0 abstract element
        weights: POINTER(POINTER(c_float_type))
            The weight matrix
        bias : POINTER(c_size_t)
            The bias vector
        num_out_neurons: c_size_t
            The number of output neurons
        num_in_neurons: c_size_t
            The number of input_neurons
        predecessors:
            the layers before the current layer
        has_tanh: c_bool
            if the last layer has a Tanh activation
        use_area_heuristic: c_bool
            whether to use area heuristic

    Returns
    -------
        None

    """

    print('"ffn_handle_last_tanh_layer_no_alloc" is not implemented in "libfppoly_gpu.so". It only exists in the CPU version of fppoly.')

    return


def ffn_handle_last_parabola_layer(man, element, weights, bias, num_out_neurons, num_in_neurons, predecessors, has_parabola, use_area_heuristic):

    """
    Handle the last FFN Parabolic layer

    Parameters
    ----------
        man : ElinaManagerPtr
            Pointer to the ElinaManager.
        element : ElinaAbstract0Ptr
            Pointer to the ElinaAbstract0 abstract element
        weights: POINTER(POINTER(c_float_type))
            The weight matrix
        bias : POINTER(c_size_t)
            The bias vector
        num_out_neurons: c_size_t
            The number of output neurons
        num_in_neurons: c_size_t
            The number of input_neurons
        predecessors:
            the layers before the current layer
        has_parabola: c_bool
            if the last layer has a Parabola activation
        use_area_heuristic: c_bool
            whether to use area heuristic
    Returns
    -------
        None

    """

    print('"ffn_handle_last_parabola_layer" is not implemented in "libfppoly_gpu.so". It only exists in the CPU version of fppoly.')

    return


def ffn_handle_last_parabola_layer_no_alloc(man, element, weights, bias, num_out_neurons, num_in_neurons, predecessors, has_parabola, use_area_heuristic):

    """
    Handle the last FFN Parabolic layer

    Parameters
    ----------
        man : ElinaManagerPtr
            Pointer to the ElinaManager.
        element : ElinaAbstract0Ptr
            Pointer to the ElinaAbstract0 abstract element
        weights: POINTER(POINTER(c_float_type))
            The weight matrix
        bias : POINTER(c_size_t)
            The bias vector
        num_out_neurons: c_size_t
            The number of output neurons
        num_in_neurons: c_size_t
            The number of input_neurons
        predecessors:
            the layers before the current layer
        has_parabola: c_bool
            if the last layer has a Parabola activation
        use_area_heuristic: c_bool
            whether to use area heuristic

    Returns
    -------
        None

    """

    print('"ffn_handle_last_parabola_layer_no_alloc" is not implemented in "libfppoly_gpu.so". It only exists in the CPU version of fppoly.')

    return


def ffn_handle_last_log_layer(man, element, weights, bias, num_out_neurons, num_in_neurons, predecessors, has_log, use_area_heuristic):

    """
    Handle the last FFN Log layer

    Parameters
    ----------
        man : ElinaManagerPtr
            Pointer to the ElinaManager.
        element : ElinaAbstract0Ptr
            Pointer to the ElinaAbstract0 abstract element
        weights: POINTER(POINTER(c_float_type))
            The weight matrix
        bias : POINTER(c_size_t)
            The bias vector
        num_out_neurons: c_size_t
            The number of output neurons
        num_in_neurons: c_size_t
            The number of input_neurons
        predecessors:
            the layers before the current layer
        has_log: c_bool
            if the last layer has a Log activation
        use_area_heuristic: c_bool
            whether to use area heuristic

    Returns
    -------
        None

    """

    print('"ffn_handle_last_log_layer" is not implemented in "libfppoly_gpu.so". It only exists in the CPU version of fppoly.')

    return


def ffn_handle_last_log_layer_no_alloc(man, element, weights, bias, num_out_neurons, num_in_neurons, predecessors, has_log, use_area_heuristic):

    """
    Handle the last FFN Log layer

    Parameters
    ----------
        man : ElinaManagerPtr
            Pointer to the ElinaManager.
        element : ElinaAbstract0Ptr
            Pointer to the ElinaAbstract0 abstract element
        weights: POINTER(POINTER(c_float_type))
            The weight matrix
        bias : POINTER(c_size_t)
            The bias vector
        num_out_neurons: c_size_t
            The number of output neurons
        num_in_neurons: c_size_t
            The number of input_neurons
        has_log: c_bool
            if the last layer has a Log activation
        use_area_heuristic: c_bool
            whether to use area heuristic

    Returns
    -------
        None

    """

    print('"ffn_handle_last_log_layer_no_alloc" is not implemented in "libfppoly_gpu.so". It only exists in the CPU version of fppoly.')

    return


def is_greater(man, element, y, x, use_area_heuristic):

    """
    Check if y is strictly greater than x in the abstract element.

    Parameters
    ----------
        man: ElinaManagerPtr
            Pointer to the ElinaManager.
        destructive : c_bool
            Boolean flag.
        y: ElinaDim
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
        is_greater_c = fppoly_gpu_api.is_greater
        is_greater_c.restype = c_bool
        is_greater_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr, ElinaDim, ElinaDim, c_bool]
        res = is_greater_c(man, element, y, x, use_area_heuristic)
    except Exception as inst:
        print('Problem with loading/calling "is_greater" from "libfppoly_gpu.so".')
        print(inst)

    return res


def conv_handle_first_layer(man, element, filter_weights, filter_bias, input_size, filter_size, num_filters, strides, output_size, pad_top, pad_left, has_bias, predecessors):

    """
    Convolutional Matrix multiplication in the first layer.

    Parameters
    ----------
        man : ElinaManagerPtr
            Pointer to the ElinaManager
        element : ElinaAbstract0Ptr
            Pointer to the ElinaAbstract0 abstract element
        filter_weights: POINTER(c_float_type)
            filter weights
        filter_bias: POINTER(c_float_type)
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
        predecessors:
            the layers before the current layer

    Returns
    -------
        None

    """

    try:
        conv_handle_first_layer_c = fppoly_gpu_api.conv_handle_first_layer
        conv_handle_first_layer_c.restype = None
        conv_handle_first_layer_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr, ndpointer(c_float_type), ndpointer(c_float_type), ndpointer(c_size_t), POINTER(c_size_t), c_size_t, POINTER(c_size_t), POINTER(c_size_t), c_size_t, c_size_t, c_bool, POINTER(c_size_t)]
        conv_handle_first_layer_c(man, element, filter_weights, filter_bias, input_size, filter_size, num_filters, strides, output_size, pad_top, pad_left, has_bias, predecessors)
    except Exception as inst:
        print('Problem with loading/calling "conv_handle_first_layer" from "libfppoly_gpu.so".')
        print(inst)

    return


def conv_handle_intermediate_relu_layer(man, element, filter_weights, filter_bias, input_size, filter_size, num_filters, strides, output_size, pad_top, pad_left, has_bias, predecessors, use_area_heuristic, retain_training_data = None, gradient = np.empty(0, dtype = float_type)):

    """
    Convolutional Matrix multiplication in an Intermediate layer.

    Parameters
    ----------
        man : ElinaManagerPtr
            Pointer to the ElinaManager
        element : ElinaAbstract0Ptr
            Pointer to the ElinaAbstract0 abstract element
        filter_weights: POINTER(c_float_type)
            filter weights
        filter_bias: POINTER(c_float_type)
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
        predecessors:
            the layers before the current layer
        use_area_heuristic: c_bool
            whether to use area heuristic
        retain_training_data: c_bool
            if true, method will retain training data
        gradient: POINTER(c_float_type)
            gradient for DeepPoly training

    Returns
    -------
        None

    """

    try:
        conv_handle_intermediate_relu_layer_c = fppoly_gpu_api.conv_handle_intermediate_relu_layer
        conv_handle_intermediate_relu_layer_c.restype = None
        conv_handle_intermediate_relu_layer_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr, ndpointer(c_float_type), ndpointer(c_float_type), ndpointer(c_size_t), POINTER(c_size_t), c_size_t, POINTER(c_size_t), POINTER(c_size_t), c_size_t, c_size_t, c_bool, POINTER(c_size_t), c_bool, c_bool, ndpointer(c_float_type)]
        conv_handle_intermediate_relu_layer_c(man, element, filter_weights, filter_bias, input_size, filter_size, num_filters, strides, output_size, pad_top, pad_left, has_bias, predecessors, use_area_heuristic, retain_training_data, gradient)
    except Exception as inst:
        print('Problem with loading/calling "conv_handle_intermediate_relu_layer" from "libfppoly_gpu.so".')
        print(inst)

    return


def conv_handle_intermediate_affine_layer(man, element, filter_weights, filter_bias, input_size, filter_size, num_filters, strides, output_size, pad_top, pad_left, has_bias, predecessors, use_area_heuristic, retain_training_data = None, gradient = np.empty(0, dtype = float_type)):

    """
    Convolutional Matrix multiplication in an Intermediate layer.

    Parameters
    ----------
        man: ElinaManagerPtr
            Pointer to the ElinaManager
        element: ElinaAbstract0Ptr
            Pointer to the ElinaAbstract0 abstract element
        filter_weights: POINTER(c_float_type)
            filter weights
        filter_bias: POINTER(c_float_type)
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
        predecessors:
            the layers before the current layer
        use_area_heuristic: c_bool
            whether to use area heuristic
        retain_training_data: c_bool
            if true, method will retain training data
        gradient: POINTER(c_float_type)
            gradient for DeepPoly training

    Returns
    -------
        None

    """

    try:
        conv_handle_intermediate_affine_layer_c = fppoly_gpu_api.conv_handle_intermediate_affine_layer
        conv_handle_intermediate_affine_layer_c.restype = None
        conv_handle_intermediate_affine_layer_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr, ndpointer(c_float_type), ndpointer(c_float_type), ndpointer(c_size_t), POINTER(c_size_t), c_size_t, POINTER(c_size_t), POINTER(c_size_t), c_size_t, c_size_t, c_bool, POINTER(c_size_t), c_bool, c_bool, ndpointer(c_float_type)]
        conv_handle_intermediate_affine_layer_c(man, element, filter_weights, filter_bias, input_size, filter_size, num_filters, strides, output_size, pad_top, pad_left, has_bias, predecessors, use_area_heuristic, retain_training_data, gradient)
    except Exception as inst:
        print('Problem with loading/calling "conv_handle_intermediate_affine_layer" from "libfppoly_gpu.so".')
        print(inst)

    return


def handle_maxpool_layer(man, element, pool_size, input_size, predecessors):

    """
    Handle the Maxpool layer.

    Parameters
    ----------
        man: ElinaManagerPtr
            Pointer to the ElinaManager.
        element: ElinaAbstract0Ptr
            Pointer to the ElinaAbstract0 abstract element.
        pool_size: POINTER(c_size_t)
            The size of the Maxpool filter
        input_size : POINTER(c_size_t)
            The number of variables on which Maxpool will be applied.
        predecessors:
            the layers before the current layer

    Returns
    -------
        res : c_size_t
            Number of neurons in the last layer

    """

    print('"handle_maxpool_layer" is not implemented in "libfppoly_gpu.so". It only exists in the CPU version of fppoly.')

    return 0


def handle_residual_relu_layer(man, element, num_neurons, predecessors, use_area_heuristic):

    """
    Handle the Residual layer.

    Parameters
    ----------
        man: ElinaManagerPtr
            Pointer to the ElinaManager.
        element: ElinaAbstract0Ptr
            Pointer to the ElinaAbstract0 abstract element.
        num_neurons: c_size_t
            The number of neurons in the residual layer
        predecessors:
            the layers before the current layer

    Returns
    -------
        None

    """

    try:
        handle_residual_relu_layer_c = fppoly_gpu_api.handle_residual_relu_layer
        handle_residual_relu_layer_c.restype = None
        handle_residual_relu_layer_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr, c_size_t, POINTER(c_size_t), c_bool]
        handle_residual_relu_layer_c(man, element, num_neurons, predecessors, use_area_heuristic)
    except Exception as inst:
        print('Problem with loading/calling "handle_residual_relu_layer" from "libfppoly_gpu.so".')
        print(inst)

    return


def handle_residual_affine_layer(man, element, num_neurons, predecessors, use_area_heuristic):

    """
    Handle the Residual layer.

    Parameters
    ----------
        man: ElinaManagerPtr
            Pointer to the ElinaManager.
        element: ElinaAbstract0Ptr
            Pointer to the ElinaAbstract0 abstract element.
        num_neurons: c_size_t
            The number of neurons in the residual layer
        predecessors:
            the layers before the current layer

    Returns
    -------
        None

    """

    try:
        handle_residual_affine_layer_c = fppoly_gpu_api.handle_residual_affine_layer
        handle_residual_affine_layer_c.restype = None
        handle_residual_affine_layer_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr, c_size_t, POINTER(c_size_t), c_bool]
        handle_residual_affine_layer_c(man, element, num_neurons, predecessors, use_area_heuristic)
    except Exception as inst:
        print('Problem with loading/calling "handle_residual_affine_layer" from "libfppoly_gpu.so".')
        print(inst)

    return


def box_for_neuron(man, element, layerno, neuron_no):

    """
    Returns bounds for a neuron in a layer.

    Parameters
    ----------
        man: ElinaManagerPtr
            Pointer to the ElinaManager
        element: ElinaAbstract0Ptr
            Pointer to the ElinaAbstract0
        layerno: c_size_t
            the layer number
        neuron_no: c_size_t
            the neuron number in the layer

    Returns
    -------
        interval_array: ElinaIntervalPtr
            ElinaIntervalArray representing the hypercube.

    """

    interval = None

    try:
        box_for_neuron_c = fppoly_gpu_api.box_for_neuron
        box_for_neuron_c.restype = ElinaIntervalPtr
        box_for_neuron_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr, c_size_t, c_size_t]
        interval = box_for_neuron_c(man, element, layerno, neuron_no)
    except:
        print('Problem with loading/calling "box_for_neuron" from "fppoly_gpu.so".')
        print('Make sure you are passing ElinaManagerPtr, ElinaAbstract0Ptr, c_size_t, c_size_t to the function.')

    return interval


def box_for_layer(man, element, layerno):

    """
    Returns bounds for all neurons in a layer.

    Parameters
    ----------
        man: ElinaManagerPtr
            Pointer to the ElinaManager
        element: ElinaAbstract0Ptr
            Pointer to the ElinaAbstract0
        layerno: c_size_t
            the layer number

    Returns
    -------
        interval_array: ElinaIntervalArray
            ElinaIntervalArray representing the hypercube.

    """

    interval_array = None

    try:
        box_for_layer_c = fppoly_gpu_api.box_for_layer
        box_for_layer_c.restype = ElinaIntervalArray
        box_for_layer_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr, c_size_t]
        interval_array = box_for_layer_c(man, element, layerno)
    except:
        print('Problem with loading/calling "box_for_layer" from "fppoly_gpu.so".')
        print('Make sure you are passing ElinaManagerPtr, ElinaAbstract0Ptr, c_size_t to the function.')

    return interval_array


def get_num_neurons_in_layer(man, element, layerno):

    """
    Returns the number of neurons in a layer.

    Parameters
    ----------
        man: ElinaManagerPtr
            Pointer to the ElinaManager
        element: ElinaAbstract0Ptr
            Pointer to the ElinaAbstract0
        layerno: c_size_t
            the layer number

    Returns
    -------
        interval_array: ElinaIntervalArray
            ElinaIntervalArray representing the hypercube

    """

    res = 0

    try:
        get_num_neurons_in_layer_c = fppoly_gpu_api.get_num_neurons_in_layer
        get_num_neurons_in_layer_c.restype = c_size_t
        get_num_neurons_in_layer_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr, c_size_t]
        res = get_num_neurons_in_layer_c(man, element, layerno)
    except:
        print('Problem with loading/calling "get_num_neurons_in_layer" from "fppoly_gpu.so".')
        print('Make sure you are passing ElinaManagerPtr, ElinaAbstract0Ptr, c_size_t to the function.')

    return res


def update_bounds_for_neuron(man, element, layerno, neuron_no, lb, ub):

    """
    Returns bounds for a neuron in a layer.

    Parameters
    ----------
        man: ElinaManagerPtr
            Pointer to the ElinaManager
        element: ElinaAbstract0Ptr
            Pointer to the ElinaAbstract0
        layerno: c_size_t
            the layer number
        neuron_no: c_size_t
            the neuron number in the layer
        lb: c_float_type
            the updated lower bound
        ub: c_float_type
            the updated upper bound

    Returns
    -------
        None

    """

    try:
        update_bounds_for_neuron_c = fppoly_gpu_api.update_bounds_for_neuron
        update_bounds_for_neuron_c.restype = None
        update_bounds_for_neuron_c.argtypes = [ElinaManagerPtr, ElinaAbstract0Ptr, c_size_t, c_size_t, c_float_type, c_float_type]
        update_bounds_for_neuron_c(man, element, layerno, neuron_no, lb, ub)
    except:
        print('Problem with loading/calling "update_bounds_for_neuron" from "fppoly_gpu.so".')
        print('Make sure you are passing ElinaManagerPtr, ElinaAbstract0Ptr, c_size_t, c_size_t, c_float_type, c_float_type to the function.')

    return


def get_bounds_for_linexpr0(man, element, linexpr0, layerno):

    """
    Returns bounds for a linexpr0 over neurons in "layerno".

    Parameters
    ----------
        man: ElinaManagerPtr
            Pointer to the ElinaManager
        element: ElinaAbstract0Ptr
            Pointer to the ElinaAbstract0
        linexpr0: ElinaLinexpr0Ptr
            Pointer to the Elinalinexpr0
        layerno: c_size_t
            the layer number

    Returns
    -------
        interval: ElinaIntervalPtr
            Poiner to the Elinainterval

    """

    interval = None

    print('"get_bounds_for_linexpr0" is not implemented in "libfppoly_gpu.so". It only exists in the CPU version of fppoly.')

    return interval


def get_lexpr_for_output_neuron(man, element, i):

    """
    Returns lower polyhedra constraint for the i-th output neuron in terms of the input neurons.

    Parameters
    ----------
        man: ElinaManagerPtr
            Pointer to the ElinaManager
        element: ElinaAbstract0Ptr
            Pointer to the ElinaAbstract0
        i: c_size_t
            Output neuron number

    Returns
    -------
        expr: ElinaLinexpr0Ptr
            The lower polyhedra expression for the output neuron in terms of input parameters and pixels

    """

    linexpr0 = None

    print('"get_lexpr_for_output_neuron" is not implemented in "libfppoly_gpu.so". It only exists in the CPU version of fppoly.')

    return linexpr0


def get_uexpr_for_output_neuron(man, element, i):

    """
    Returns lower polyhedra constraint for the i-th output neuron in terms of the input neurons.

    Parameters
    ----------
        man: ElinaManagerPtr
            Pointer to the ElinaManager.
        element: ElinaAbstract0Ptr
            Pointer to the ElinaAbstract0.
        i: c_size_t
            Output neuron number

    Returns
    -------
        expr: ElinaLinexpr0Ptr
            The upper polyhedra expression for the output neuron in terms of input parameters and pixels

    """

    linexpr0 = None

    print('"get_uexpr_for_output_neuron" is not implemented in "libfppoly_gpu.so". It only exists in the CPU version of fppoly.')

    return linexpr0


def create_lstm_layer(man, element, h, predecessors):

    """
    Creates an lstm layer for the neural network, this should be called only once per each lstm layer

    Parameters
    ----------
        man : ElinaManagerPtr
            Pointer to the ElinaManager.
        element : ElinaAbstract0Ptr
            Pointer to the ElinaAbstract0.
        h: c_size_t
            size of h_t
        predecessors:
            the layers before the current layer

    Returns
    --------
        None

    """

    print('"create_lstm_layer" is not implemented in "libfppoly_gpu.so". It only exists in the CPU version of fppoly.')

    return


def handle_lstm_layer(man, element, weights, bias, d, h, predecessors, use_area_heuristic):

    """
    Computes the hidden states and output vectors of the lstm unit, to be called at each time step after creating an LSTM unit

    Parameters
    -----------
        man : ElinaManagerPtr
            Pointer to the ElinaManager.
        element : ElinaAbstract0Ptr
            Pointer to the ElinaAbstract0.
        weights : POINTER(POINTER(c_float_type))
            The weight matrix of size 4*h \times d+h, with h rows each for f_t, i_t, o_t, and c_t in order,
            columnwise the first d entries correspond to x_t and the remaining correspond to h_t
        bias : POINTER(c_float_type)
            The bias vector of size 4*h, in the same format as weights
        d: c_size_t
           size of x_t
        h: c_size_t
           size of h_t
        predecessors:
            the layers before the current layer
        use_area_heuristic: c_bool
            whether to use area heuristic

    Returns
    --------
        None

    """

    print('"handle_lstm_layer" is not implemented in "libfppoly_gpu.so". It only exists in the CPU version of fppoly.')

    return


def free_non_lstm_layer_expr(man, element, layerno):

    """
    Returns bounds for a linexpr0 over neurons in "layerno"

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

    print('"free_non_lstm_layer_expr" is not implemented in "libfppoly_gpu.so". It only exists in the CPU version of fppoly.')

    return


def clean_training_data():

    """
    Frees all training data arrays.

    Parameters
    ----------
        None

    Returns
    -------
        num_neurons: c_int
            Number of neurons
    """

    res = None

    try:
        clean_training_data_c = fppoly_gpu_api.clean_training_data
        clean_training_data_c.restype = None
        clean_training_data_c.argtypes = None
        res = clean_training_data_c()
    except Exception as inst:
        print('Problem with loading/calling "clean_training_data" from "libfppoly_gpu.so".')
        print(inst)

    return res


def get_num_neurons_training_layer():

    """
    Receive the number of neurons of the layer currently trained.

    Parameters
    ----------
        None

    Returns
    -------
        num_neurons: c_int
            Number of neurons
    """

    res = 0

    try:
        get_num_neurons_training_layer_c = fppoly_gpu_api.get_num_neurons_training_layer
        get_num_neurons_training_layer_c.restype = c_int
        get_num_neurons_training_layer_c.argtypes = None
        res = get_num_neurons_training_layer_c()
    except Exception as inst:
        print('Problem with loading/calling "get_num_neurons_training_layer" from "libfppoly_gpu.so".')
        print(inst)

    return res


def get_adv():

    """
    Receive pointer to adv vector

    Parameters
    ----------
        None

    Returns
    -------
        lcst: POINTER(c_float_type)
            Pointer to the array
    """

    res = None

    try:
        get_adv_c = fppoly_gpu_api.get_adv
        get_adv_c.restype = POINTER(c_float_type)
        get_adv_c.argtypes = None
        res = get_adv_c()
    except Exception as inst:
        print('Problem with loading/calling "get_adv" from "libfppoly_gpu.so"')
        print(inst)

    return res
