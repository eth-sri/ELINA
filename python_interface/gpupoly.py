#  GPUPoly library
#  This source file is part of ELINA (ETH LIbrary for Numerical Analysis).
#  ELINA is Copyright � 2020 Department of Computer Science, ETH Zurich
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

## @file python/gpupoly.py
## @author Fran&ccedil;ois Serre
## @brief Python bindings for the GPUPoly library.
#
#  Defines the class Network, an interface to use the GPUPoly GPU library from Python.
#

import ctypes.util
import os
import numpy as np

## Python friendly interface to the GPUPoly GPU library.
class Network:
    if os.name == 'nt':
        os.add_dll_directory("${CUDAToolkit_BIN_DIR}")
        os.add_dll_directory("${GPUPoly_BINARY_DIR}")
        _lib = ctypes.cdll.LoadLibrary(ctypes.util.find_library('gpupoly'))
    else:
        #_lib=ctypes.cdll.LoadLibrary('${GPUPoly_BINARY_DIR}/dpGPUlib.so')
        _lib=ctypes.cdll.LoadLibrary('/local/home/francois/vs/SALSa/out/libgpupoly.so')

    _lib.create.argtypes = [ctypes.c_int]
    _lib.create.restype = ctypes.c_void_p
    _lib.test_d.argtypes = [
        ctypes.c_void_p,
        np.ctypeslib.ndpointer(dtype=np.float64, ndim=1),
        np.ctypeslib.ndpointer(dtype=np.float64, ndim=1),
        ctypes.c_int,
        ctypes.c_bool
    ]
    _lib.test_d.restype = ctypes.c_int
    _lib.test_s.argtypes = [
        ctypes.c_void_p,
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1),
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1),
        ctypes.c_int,
        ctypes.c_bool
    ]
    _lib.test_s.restype = ctypes.c_int
    _lib.clean.argtypes=[ctypes.c_void_p]
    _lib.clean.restype = None
    _lib.addLinear_d.argtypes = [
        ctypes.c_void_p,
        ctypes.c_int,
        ctypes.c_int,
        np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='C_CONTIGUOUS')
    ]
    _lib.addLinear_s.argtypes = [
        ctypes.c_void_p,
        ctypes.c_int,
        ctypes.c_int,
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS')
    ]
    _lib.addConv2D_d.argtypes = [
        ctypes.c_void_p,
        ctypes.c_int,
        ctypes.c_bool,
        ctypes.c_int,
        ctypes.c_int * 2,
        ctypes.c_int * 4,
        ctypes.c_int * 2,
        ctypes.c_int * 2,
        np.ctypeslib.ndpointer(dtype=np.float64, ndim=4, flags='C_CONTIGUOUS')
    ]
    _lib.addConv2D_s.argtypes = [
        ctypes.c_void_p,
        ctypes.c_int,
        ctypes.c_bool,
        ctypes.c_int,
        ctypes.c_int * 2,
        ctypes.c_int * 4,
        ctypes.c_int * 2,
        ctypes.c_int * 2,
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=4, flags='C_CONTIGUOUS')
    ]
    _lib.addBias_d.argtypes = [
        ctypes.c_void_p,
        ctypes.c_int,
        np.ctypeslib.ndpointer(dtype=np.float64, ndim=1)
    ]
    _lib.addBias_s.argtypes = [
        ctypes.c_void_p,
        ctypes.c_int,
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1)
    ]
    _lib.addReLU.argtypes = [
        ctypes.c_void_p,
        ctypes.c_int
    ]
    _lib.addMaxPool2D.argtypes = [
        ctypes.c_void_p,
        ctypes.c_int,
        ctypes.c_bool,
        ctypes.c_int * 2,
        ctypes.c_int * 4,
        ctypes.c_int * 2,
        ctypes.c_int * 2
    ]
    _lib.addParSum.argtypes = [
        ctypes.c_void_p,
        ctypes.c_int,
        ctypes.c_int
    ]
    _lib.addConcat.argtypes = [
        ctypes.c_void_p,
        ctypes.c_int,
        ctypes.c_int
    ]

    ## Constructor
    #
    #  Construct a new neural network to verify.
    #  @param inputSize Number of elements of the input layer.
    def __init__(self, input_size):
        self._nn = self._lib.create(input_size)
        self.input_size = input_size
        self._last_layer_id = 0

    def clean(self):
        if self._nn is not None:
            self._lib.clean(self._nn)
            self._nn=None
    ## Destructor
    def __del__(self):
        self.clean()

    ## Verifies an image
    #
    #  Checks if an image (or more precisely, an input box) entirely classifies as a given label
    #  Once this function has been called with a network, no additional layer should be added (i.e. it is necessary to clean and recreate a network)
    #  It is of course possible to call this function multiple times on the same network.
    #
    #  \param down A numpy array of inputSize doubles that represent the lower bound of the box
    #  \param up A numpy array of inputSize doubles that represent the upper bound of the box
    #  \param label Label in which the image is supposed to classify
    #  \param soundness Whether to use sound arithmetic.
    #  \returns whether the image could be certified.
    def test(self, down, up, label, soundness=True):
        down=down.flatten()
        up=up.flatten()
        assert down.shape == (self.input_size,)
        assert up.shape == (self.input_size,)
        assert up.dtype == down.dtype
        if up.dtype==np.float64:
            return self._lib.test_d(self._nn, down, up, label, soundness)
        return self._lib.test_s(self._nn, down, up, label, soundness)

    ## Fully connected linear layer
    #
    #  Adds a dense linear layer (without activation nor bias). If x is a column vector representing the input of the layer, the output is A*x.
    #
    #  \param A Numpy matrix that represents A. A has outputSize rows, and its number of columns equals the parent outputSize.
    #  \param parent Index of the parent layer (or 0 for the input layer). It can be None, in which case the parent is the last added layer.
    #  \returns the index of the newly created layer.
    def add_linear(self, a, parent=None):
        if parent is None:
            parent = self._last_layer_id
        assert a.ndim == 2
        output_size = a.shape[0]
        if a.dtype==np.float64:
            self._last_layer_id = self._lib.addLinear_d(self._nn, parent, output_size, np.ascontiguousarray(a))
        else:
            self._last_layer_id = self._lib.addLinear_s(self._nn, parent, output_size, np.ascontiguousarray(a))
        return self._last_layer_id

    ## ReLU layer
    #
    #  Adds a ReLU layer to the network.
    #
    #  \param parent Index of the parent layer (or 0 for the input layer). It can be None, in which case the parent is the last added layer.
    #  \returns the index of the newly created layer.
    def add_relu(self, parent=None):
        if parent is None:
            parent = self._last_layer_id
        self._last_layer_id = self._lib.addReLU(self._nn, parent)
        return self._last_layer_id

    ## Bias layer
    #
    #  Adds a Bias layer to the network, i.e. a layer that adds a constant vector to its input.
    #
    #  \param parent Index of the parent layer (or 0 for the input layer). It can be None, in which case the parent is the last added layer.
    #  \returns the index of the newly created layer.
    def add_bias(self, b, parent=None):
        if parent is None:
            parent = self._last_layer_id
        assert b.ndim==1
        if b.dtype==np.float64:
            self._last_layer_id = self._lib.addBias_d(self._nn, parent, b)
        else:
            self._last_layer_id = self._lib.addBias_s(self._nn, parent, b)
        return self._last_layer_id


    ## Convolution layer
    #
    #  Adds a convolution layer, without activation nor bias.
    #
    #  \param input_rows Dimension of the input.
    #  \param input_cols Dimension of the input.
    #  \param conv Convolution coefficients, given as a 4 dimensions numpy array  (in order [row, column, channel, filter].
    #  \param channels_first If true, the layer expects input with the shape [batch, channel, row, col], and its output has the shape [batch, filter, row, col]. If false, the layer expects input with the shape [batch, row, col, channel], and its output has the shape [batch, row, col, filter].
    #  \param batches Number of batches.
    #  \param strides An integer or a list of two integers, indicating the stride shape (respectively the number of rows and columns if a list).
    #  \param padding An integer or a list of two integers, indicating the padding (respectively the number of pixels to add at the top and bottom, and the number of pixels to add on the left and right).
    #  \param parent Index of the parent layer (or 0 for the input layer). It can be None, in which case the parent is the last added layer.
    #  \returns the index of the newly created layer.
    def add_conv_2d(self, input_rows, input_cols, conv, channel_first=True, batches=1, strides=1, padding=0,
                    parent=None):
        if parent is None:
            parent = self._last_layer_id
        assert conv.ndim == 4
        kernel = conv.shape[0:2]
        input_shape = [batches, input_rows, input_cols, conv.shape[2]]
        filters = conv.shape[3]
        if not isinstance(strides, list):
            strides = [strides, strides]
        if not isinstance(padding, list):
            padding = [padding, padding]
        if conv.dtype==np.float64:
            self._last_layer_id = self._lib.addConv2D_d(
                self._nn,
                parent,
                channel_first,
                filters,
                (ctypes.c_int * 2)(*kernel),
                (ctypes.c_int * 4)(*input_shape),
                (ctypes.c_int * 2)(*strides),
                (ctypes.c_int * 2)(*padding),
                np.ascontiguousarray(conv))
        else:
            self._last_layer_id = self._lib.addConv2D_s(
                self._nn,
                parent,
                channel_first,
                filters,
                (ctypes.c_int * 2)(*kernel),
                (ctypes.c_int * 4)(*input_shape),
                (ctypes.c_int * 2)(*strides),
                (ctypes.c_int * 2)(*padding),
                np.ascontiguousarray(conv))
        return self._last_layer_id

    ## MaxPool2D layer.
    #
    #  Adds a max pooling layer.
    #
    #  \param pool Pool shape (Python list of 2 integers, respectively the number of rows and columns).
    #  \param input_rows Dimension of the input.
    #  \param input_cols Dimension of the input.
    #  \param channels Number of channels.
    #  \param channels_first If true, the layer expects input with the shape [batch, channel, row, col], and its output has the shape [batch, filter, row, col]. If false, the layer expects input with the shape [batch, row, col, channel], and its output has the shape [batch, row, col, filter].
    #  \param batches Number of batches.
    #  \param strides An integer or a list of two integers, indicating the stride shape (respectively the number of rows and columns if a list).
    #  \param padding An integer or a list of two integers, indicating the padding (respectively the number of pixels to add at the top and bottom, and the number of pixels to add on the left and right).
    #  \param parent Index of the parent layer (or 0 for the input layer). It can be None, in which case the parent is the last added layer.
    #  \returns the index of the newly created layer.
    def add_maxpool_2d(self, pool, input_rows, input_cols, channels, channel_first=True, batches=1, strides=None,
                       padding=0,
                       parent=None):
        if parent is None:
            parent = self._last_layer_id
        if not isinstance(pool, list):
            pool = [pool, pool]
        if strides is None:
            strides=pool
        elif not isinstance(strides, list):
            strides = [strides, strides]
        if not isinstance(padding, list):
            padding = [padding, padding]
        input_shape = [batches, input_rows, input_cols, channels]
        self._last_layer_id = self._lib.addMaxPool2D(
            self._nn,
            parent,
            channel_first,
            (ctypes.c_int * 2)(*pool),
            (ctypes.c_int * 4)(*input_shape),
            (ctypes.c_int * 2)(*strides),
            (ctypes.c_int * 2)(*padding))
        return self._last_layer_id

    ## ParSum layer.
    #
    #  Adds a "ParSum" layer, i.e. a layer that sums up the result of two previous layers.
    #  \param parent1 Index of the first parent layer.
    #  \param parent2 Index of the second parent layer.
    #  \returns the index of the newly created layer.
    def add_parsum(self, parent1, parent2):
        self._last_layer_id = self._lib.addParSum(self._nn, parent1, parent2)
        return self._last_layer_id

    ## Concatenation layer.
    #
    #  Adds a concatenation layer (like the one we can find in skipnets). It concatenates the result of two previous layers.
    #  \param parent1 Index of the first parent layer.
    #  \param parent2 Index of the second parent layer.
    #  \returns the index of the newly created layer.
    def add_concat(self, parent1, parent2):
        self._last_layer_id = self._lib.addConcat(self._nn, parent1, parent2)
        return self._last_layer_id

