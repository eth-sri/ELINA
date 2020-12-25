#  GPUPoly library
#  This source file is part of ELINA (ETH LIbrary for Numerical Analysis).
#  ELINA is Copyright © 2020 Department of Computer Science, ETH Zurich
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

## @file python_interface/gpupoly.py
## @author Fran&ccedil;ois Serre
## @brief Python bindings for the GPUPoly library.
#
#  Defines the class Network, an interface to use the GPUPoly library from Python.
#

import ctypes.util
import os
import numpy as np


## Python friendly interface to the GPUPoly library.
class Network:
    if os.name == 'nt':
        #os.add_dll_directory("${CUDAToolkit_BIN_DIR}")
        #os.add_dll_directory("${GPUPoly_BINARY_DIR}")
        os.add_dll_directory("${CUDAToolkit_BIN_DIR}")
        os.add_dll_directory("${GPUPoly_BINARY_DIR}")
        _lib = ctypes.cdll.LoadLibrary(ctypes.util.find_library('gpupoly'))
    else:
        # _lib=ctypes.cdll.LoadLibrary('${GPUPoly_BINARY_DIR}/dpGPUlib.so.0.10')
        _lib = ctypes.cdll.LoadLibrary(os.path.join(os.path.dirname(os.path.abspath(__file__)),"../gpupoly/libgpupoly.so"))

    def _nullable_ndptr(*args, **kwargs):
        base = np.ctypeslib.ndpointer(*args, **kwargs)
        def from_param(cls, obj):
            if obj is None:
                return obj
            return base.from_param(obj)
        return type(base.__name__, (base,), {'from_param': classmethod(from_param)})

    _lib.create.argtypes = [ctypes.c_int]
    _lib.create.restype = ctypes.c_void_p

    _lib.test_d.argtypes = [
        ctypes.c_void_p,
        np.ctypeslib.ndpointer(dtype=np.float64, ndim=1),
        np.ctypeslib.ndpointer(dtype=np.float64, ndim=1),
        ctypes.c_int,
        ctypes.c_bool
    ]
    _lib.test_d.restype = ctypes.c_bool
    _lib.test_s.argtypes = [
        ctypes.c_void_p,
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1),
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1),
        ctypes.c_int,
        ctypes.c_bool
    ]
    _lib.test_s.restype = ctypes.c_bool

    _lib.setLayerBox_d.argtypes = [
        ctypes.c_void_p,
        np.ctypeslib.ndpointer(dtype=np.float64, ndim=1),
        np.ctypeslib.ndpointer(dtype=np.float64, ndim=1),
        ctypes.c_int
    ]
    _lib.setLayerBox_d.restype = None
    _lib.setLayerBox_s.argtypes = [
        ctypes.c_void_p,
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1),
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1),
        ctypes.c_int
    ]
    _lib.setLayerBox_s.restype = None

    _lib.relax_s.argtypes = [
        ctypes.c_void_p,
        ctypes.c_int,
        ctypes.c_bool,
        ctypes.c_bool
    ]
    _lib.relax_s.restype = None
    _lib.relax_d.argtypes = [
        ctypes.c_void_p,
        ctypes.c_int,
        ctypes.c_bool,
        ctypes.c_bool
    ]
    _lib.relax_d.restype = None

    _lib.evalAffineExpr_d.argtypes = [
        ctypes.c_void_p,
        np.ctypeslib.ndpointer(dtype=np.float64, ndim=2),
        ctypes.c_int,
        ctypes.c_int,
        _nullable_ndptr(dtype=np.float64, ndim=2, flags='C_CONTIGUOUS'),
        _nullable_ndptr(dtype=np.float64, ndim=1),
        ctypes.c_int,
        ctypes.c_bool
    ]
    _lib.evalAffineExpr_d.restype = None
    _lib.evalAffineExpr_s.argtypes = [
        ctypes.c_void_p,
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=2),
        ctypes.c_int,
        ctypes.c_int,
        _nullable_ndptr(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
        _nullable_ndptr(dtype=np.float32, ndim=1),
        ctypes.c_int,
        ctypes.c_bool
    ]
    _lib.evalAffineExpr_s.restype = None

    _lib.getOutputSize.argtypes = [ctypes.c_void_p, ctypes.c_int]
    _lib.getOutputSize.restype = ctypes.c_int

    _lib.clean.argtypes = [ctypes.c_void_p]
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
        ctypes.c_int,
        ctypes.c_int * 2,
        ctypes.c_int * 3,
        ctypes.c_int * 2,
        ctypes.c_int * 2,
        np.ctypeslib.ndpointer(dtype=np.float64, ndim=4, flags='C_CONTIGUOUS')
    ]
    _lib.addConv2D_s.argtypes = [
        ctypes.c_void_p,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.c_int * 2,
        ctypes.c_int * 3,
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
        ctypes.c_int * 2,
        ctypes.c_int * 3,
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

    NO_BACKSUBSTITUTION = 0
    FULL_BACKSUBSTITUTION = 1
    BACKSUBSTITUTION_WHILE_CONTAINS_ZERO = 2

    ## Constructor
    #
    #  Construct a new neural network to verify.
    #  @param inputSize Number of elements of the input layer.
    def __init__(self, input_size):
        self._nn = self._lib.create(input_size)
        self.input_size = input_size
        self._last_layer_id = 0

    ## Removes all layers of the network.
    def clean(self):
        if self._nn is not None:
            self._lib.clean(self._nn)
            self._nn = None

    ## Destructor. Note that python may delay the call to this method, a manual call to "clean" is therefore required in case another network would need to be loaded.
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
    #  \returns whether the difference between the output logits wrt L_oo ball.
    def test(self, down, up, label, soundness=True):
        down = down.flatten()
        up = up.flatten()
        self._last_dtype=up.dtype
        assert down.shape == (self.input_size,)
        assert up.shape == (self.input_size,)
        assert up.dtype == down.dtype
        if up.dtype == np.float64:
            return self._lib.test_d(self._nn, down, up, label, soundness)
        return self._lib.test_s(self._nn, down, up, label, soundness)


    def test_old(self, down, up, label, soundness=True):
        self.setLayerBox(down, up)  # Set input layer concrete bounds

        # creates a matrix that computes the difference with the expected layer.
        diffMatrix = np.delete(-np.eye(self._lib.getOutputSize(self._nn, self._last_layer_id)), label, 0)
        diffMatrix[:, label] = 1
        diffMatrix = diffMatrix.astype(self._last_dtype)

        # relax all layers, using simple interval analysis first.
        for i in range(self._last_layer_id):
            self.relax(i + 1, soundness=soundness, refineActivationsInput=False)

        # Evaluates an expression that computes the difference between expected label and each element of the output layer
        res = self.evalAffineExpr(diffMatrix, back_substitute=self.BACKSUBSTITUTION_WHILE_CONTAINS_ZERO, sound=soundness)
        #print("res1 ", res)
        if (res > 0).all(): # Expected layer is higher than all others
            return True

        # We failed to verify, so we redo the analysis with backsubstitution before activation layer.
        for i in range(self._last_layer_id):
            self.relax(i + 1, soundness=soundness)
        res = self.evalAffineExpr(diffMatrix, back_substitute=self.BACKSUBSTITUTION_WHILE_CONTAINS_ZERO, sound=soundness)
        #print("res2 ", res)
        return (res>0).all()


    def eval(self, x):
        self.setLayerBox(x, x)  # Set input layer concrete bounds
        # relax all layers, using simple interval analysis first.
        for i in range(self._last_layer_id):
            self.relax(i + 1, soundness=False, refineActivationsInput=False)
        return self.evalAffineExpr()


    ## Sets the concrete bounds of a layer.
    #
    #  \param down A numpy array represents the lower bound of the box. Must have the size of the layer.
    #  \param up A numpy array represents the upper bound of the box. Must have the size of the layer.
    #  \param layer Index of the layer (by default, sets the input box).
    def setLayerBox(self, down, up, layer=0):
        down = down.flatten()
        up = up.flatten()
        assert up.dtype == down.dtype
        self._last_dtype = up.dtype
        if up.dtype == np.float64:
            return self._lib.setLayerBox_d(self._nn, down, up, layer)
        return self._lib.setLayerBox_s(self._nn, down, up, layer)

    ## Propagates forward the concrete bounds by interval analysis. Activation layers have their approximation models computed.
    #
    # \param layer Index of the layer
    # \param refineActivationsInput If true, and layer is an activation layer, then the input is first refined first via the appropriate back-substitution.
    # \param soundness Whether to use sound (but slower) arithmetic.
    # \param dtype Datatype of the concrete bounds to be used. Can be np.float32, np.float64, or None, in which case it uses the same type as the last setLayerBox.
    def relax(self, layer, refineActivationsInput=True, soundness=True, dtype=None):
        if dtype is None:
            dtype = self._last_dtype
        if dtype == np.float64:
            self._lib.relax_d(self._nn, layer, refineActivationsInput, soundness)
        else:
            self._lib.relax_s(self._nn, layer, refineActivationsInput, soundness)

    ## Evaluate the concrete bounds of a list of affine expressions.
    #
    # Evaluate the concrete bounds of a list of m affine expressions of the neurons of a given layer via back-substitution.
    # The affine expressions have the form Ax+b, where A is a m*n matrix, b a vector of size m, and x represents the n neurons of the layer layerId.
    # \param a A numpy array of dimension [m,n]. If None, equivalent to the identity.
    # \param b A numpy array of dimension [m]. If None, equivalent to a zero vector.
    # \param layer Index of the layer
    # \param back_substitute If set to FULL_BACKSUBSITUTION, back-substitution is always performed back to the inputs. If set to BACKSUBSTITUTION_WHILE_CONTAINS_ZERO, the backsubstitution is stopped as soon as 0 is not within the concrete bounds anymore. If set to NO_BACKSUBSTITUTION, only use the current concrete bounds of the layer.
    # \param sound If True, use floating-point sound arithmetic (slower).
    # \param dtype Datatype of the concrete bounds to be used. Can be np.float32, np.float64, or None, in which case it uses the same type as the last setLayerBox.
    # \returns A numpy array of size [m,2] containing the concrete bounds.
    def evalAffineExpr(self, a=None, b=None, layer=None, back_substitute=NO_BACKSUBSTITUTION, sound=True, dtype=None):
        if layer is None:
            layer = self._last_layer_id
        if dtype is None:
            dtype = self._last_dtype
        n = self._lib.getOutputSize(self._nn, layer)
        if a is None:
            m = n
        else:
            assert a.ndim == 2
            assert a.dtype == dtype
            m = a.shape[0]
            a = np.ascontiguousarray(a)
        if b is not None:
            assert b.shape == (m,)
            assert b.dtype == dtype
        res = np.ascontiguousarray(np.ndarray((m, 2), dtype=dtype))
        if dtype == np.float64:
            self._lib.evalAffineExpr_d(self._nn, res, layer, m, a, b, back_substitute, sound)
        else:
            self._lib.evalAffineExpr_s(self._nn, res, layer, m, a, b, back_substitute, sound)
        return np.reshape(res, (m, 2))


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
        if a.dtype == np.float64:
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
        assert b.ndim == 1
        if b.dtype == np.float64:
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
    def add_conv_2d(self, input_rows, input_cols, conv, strides=1, padding=0,parent=None):
        if parent is None:
            parent = self._last_layer_id
        assert conv.ndim == 4
        kernel = conv.shape[0:2]
        input_shape = [input_rows, input_cols, conv.shape[2]]
        filters = conv.shape[3]
        if not isinstance(strides, list):
            strides = [strides, strides]
        if not isinstance(padding, list):
            padding = [padding, padding]
        if conv.dtype == np.float64:
            self._last_layer_id = self._lib.addConv2D_d(
                self._nn,
                parent,
                filters,
                (ctypes.c_int * 2)(*kernel),
                (ctypes.c_int * 3)(*input_shape),
                (ctypes.c_int * 2)(*strides),
                (ctypes.c_int * 2)(*padding),
                np.ascontiguousarray(conv))
        else:
            self._last_layer_id = self._lib.addConv2D_s(
                self._nn,
                parent,
                filters,
                (ctypes.c_int * 2)(*kernel),
                (ctypes.c_int * 3)(*input_shape),
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
    def add_maxpool_2d(self, pool, input_rows, input_cols, channels, strides=None,
                       padding=0,
                       parent=None):
        if parent is None:
            parent = self._last_layer_id
        if not isinstance(pool, list):
            pool = [pool, pool]
        if strides is None:
            strides = pool
        elif not isinstance(strides, list):
            strides = [strides, strides]
        if not isinstance(padding, list):
            padding = [padding, padding]
        input_shape = [input_rows, input_cols, channels]
        self._last_layer_id = self._lib.addMaxPool2D(
            self._nn,
            parent,
            (ctypes.c_int * 2)(*pool),
            (ctypes.c_int * 3)(*input_shape),
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
