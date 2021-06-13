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


## @file python_interface/onnx2gpupoly.py
## @author Fran&ccedil;ois Serre
## @brief Onnx parser for the GPUPoly library.
#
#  Defines a function that parses an ONNX network, and returns a GPUPoly network handle that can be used to certify images.
#


import onnx
import onnx.utils
import onnx.numpy_helper
import numpy as np

from gpupoly import Network

## Load a pre-trained Onnx network into the GPUPoly library.
def onnx2gpupoly(graph, useAreaHeuristic=True):
    # dictionary that will containing the layer outputs/values that have already been processed.
    # Each value is a tuple containing an integer as a first element.
    # If the integer is -1, the value is a constant, and this constant is given as a numpy ndarray as a second element.
    # otherwise, the first value is the index of the corresponding layer in GPUPoly, and the second element is a list containing its shape.
    implemented = {}

    # Process the initializers of the graph
    for constant in graph.initializer:
        implemented[constant.name] = [-1, onnx.numpy_helper.to_array(constant)]

    # Process the input of the graph, and creates a new GPUPoly network with the corresponding size
    inputName = None
    for i in graph.input:
        if not i.name in implemented:
            assert inputName is None, "Graph must contain only one input"
            inputName = i.name
            inputShape = []
            inputType = i.type.tensor_type.elem_type
            for j in i.type.tensor_type.shape.dim:
                inputShape.append(j.dim_value)
            implemented[inputName] = [0, inputShape] # Input layer has the index 0 in GPUPoly.

    def product(list):
        res=1
        for i in list:
            res=res*i
        return res
    nn=Network(product(inputShape))

    # Stores the layer of the graph in a dictionary indexed by their output name.
    layers = {}
    for layer in graph.node:
        layers[layer.output[0]] = layer

    # Check if the name of a value has already been processed; otherwise process the corresponding operator.
    def getLayer(name):
        if name in implemented:
            return implemented[name]
        res = implement(layers[name])
        implemented[name] = res
        return res

    # Actual processing of Onnx operators
    def implement(layer):
        if (layer.op_type == "Add"):
            lhsIndex, lhs = getLayer(layer.input[0])
            rhsIndex, rhs = getLayer(layer.input[1])
            # Both inputs are constants
            if lhsIndex == -1 and rhsIndex == -1:
                return -1, lhs + rhs
            # Both inputs are actual layers
            if lhsIndex > -1 and rhsIndex > -1:
                assert lhs == rhs, "Incompatible dimensions (multidirectional broadcasting is not yet fully supported)"
                return nn.add_parsum(lhsIndex,rhsIndex),lhs
            # From now, one of the input is a layer, the other a constant.
            if rhsIndex == -1:
                assert product(lhs) == rhs.size, "Incompatible dimensions (multidirectional broadcasting is not yet fully supported)"
                return nn.add_bias(rhs,lhsIndex), lhs
            assert product(rhs) == lhs.size, "Incompatible dimensions (multidirectional broadcasting is not yet fully supported)"
            return nn.add_bias(lhs, rhsIndex), rhs
        if (layer.op_type=="BatchNormalization"):
            return getLayer(layer.input[0])
        if (layer.op_type=="Cast"):
            lhsIndex, lhs = getLayer(layer.input[0])
            assert lhsIndex==-1,"Cast only supports constant input"
            supportedTypes={0:"float_",1:"single",2:"uint8",3:"int8",4:"uint16",5:"int16",6:"int32",7:"int64",9:"bool_",10:"float16",11:"double",12:"uint32",13:"uint64"}
            assert layer.attribute[0].name == "to", "Wrong attribute for cast operator"
            assert layer.attribute[0].i in supportedTypes, "Unknown type number"
            return -1,lhs.astype(supportedTypes[layer.attribute[0].i])
        if (layer.op_type == "Concat"):
            inputs = []
            for inp in layer.input:
                inputIndex, inputShape = getLayer(inp)
                assert inputIndex == -1,"Only concatenation of constants is currently supported."
                inputs.append(inputShape)
            axis = layer.attribute[0].i
            return -1, np.concatenate(inputs, axis)
        if (layer.op_type == "Constant"):
            return -1, onnx.numpy_helper.to_array(layer.attribute[0].t)
        if (layer.op_type == "Conv"):
            inputIndex, inputShape = getLayer(layer.input[0])
            dataIndex, data = getLayer(layer.input[1])
            padding_top = 0
            padding_left = 0
            padding_bottom = 0
            padding_right = 0
            stride_rows = 1
            stride_cols = 1
            group=1
            for attribute in layer.attribute:
                if attribute.name == "pads":
                    padding_top = attribute.ints[0]
                    padding_left = attribute.ints[1]
                    padding_bottom = attribute.ints[2]
                    padding_right = attribute.ints[3]
                if attribute.name == "strides":
                    stride_rows = attribute.ints[0]
                    stride_cols = attribute.ints[1]
                if attribute.name == "group":
                    group = attribute.i
            assert dataIndex == -1
            assert data.ndim == 4
            print(inputShape)
            print(data.shape)
            if group==2:
                ndata=np.zeros((data.shape[0],2*data.shape[1],data.shape[2],data.shape[3]),data.dtype)
                ndata[:(data.shape[0]//2), :data.shape[1], :, :] = data[:(data.shape[0]//2), :, :, :]
                ndata[(data.shape[0] // 2):, data.shape[1]:, :, :] = data[(data.shape[0] // 2):, :, :, :]
                data=ndata
            assert inputShape[0] == 1
            assert data.shape[1] == inputShape[-3]
            outputShape = inputShape[:-3]
            outputShape.append(data.shape[0])
            outputShape.append((inputShape[-2] + padding_top + padding_bottom - data.shape[-2] + stride_rows) // stride_rows)
            outputShape.append((inputShape[-1] + padding_right + padding_left - data.shape[-1] + stride_cols) // stride_cols)
            outputIndex = nn.add_conv_2d(inputShape[-2],inputShape[-1],np.transpose(data,[2,3,1,0]),[stride_rows,stride_cols],[padding_top,padding_left,padding_bottom,padding_right],inputIndex)
            if len(layer.input)==3:
                dataIndex, data = getLayer(layer.input[2])
                assert dataIndex==-1
                assert data.shape==(outputShape[-3],)
                data=data.repeat(outputShape[-2]*outputShape[-1])
                assert(data.size==product(outputShape))
                outputIndex=nn.add_bias(data,outputIndex)
            return outputIndex, outputShape
        if (layer.op_type=="Dropout"):
            return getLayer(layer.input[0])
        if (layer.op_type=="Flatten"):
            inputIndex, inputShape = getLayer(layer.input[0])
            assert inputIndex>=0
            outputShape=[inputShape[0],product(inputShape[1:])]
            return inputIndex,outputShape
        if (layer.op_type == "Gather"):
            
            lhsIndex, lhs = getLayer(layer.input[0])
            rhsIndex, rhs = getLayer(layer.input[1])
            assert lhsIndex == -1 and rhsIndex == -1
            axis = layer.attribute[0].i
            return -1, np.take(lhs, rhs, axis)

        if (layer.op_type == "Gemm"):
            inputIndex, inputShape = getLayer(layer.input[0])
            dataIndex, data = getLayer(layer.input[1])
            assert dataIndex == -1
            assert data.ndim == 2
            assert inputShape[-1] == data.shape[1]
            outputShape = inputShape[:-1]
            outputShape.append(data.shape[0])
            outputIndex = nn.add_linear(data,inputIndex)
            #print("Gemm Matrix", data)
            if len(layer.input)==3:
                dataIndex, data = getLayer(layer.input[2])
                assert dataIndex==-1
                assert data.shape==(outputShape[-1],)
                #print("Gemm bias ", data)
                outputIndex=nn.add_bias(data,outputIndex)
            return outputIndex, outputShape
        if (layer.op_type=="Identity"):
            return getLayer(layer.input[0])
        if (layer.op_type == "MatMul"):
            inputIndex, inputShape = getLayer(layer.input[0])
            dataIndex, data = getLayer(layer.input[1])
            assert dataIndex == -1
            assert data.ndim == 2
            assert inputShape[-1] == data.shape[0]
            outputShape = inputShape[:-1]
            outputShape.append(data.shape[1])
            outputIndex = len(implemented)
            print(f"({outputIndex}) <- MatMul({inputIndex}, data {data.shape}) : {inputShape} -> {outputShape}")
            return outputIndex, outputShape
        if (layer.op_type == "MaxPool"):
            inputIndex, inputShape = getLayer(layer.input[0])
            padding_rows = 0
            padding_cols = 0
            stride_rows = 1
            stride_cols = 1
            kernel_rows = 1
            kernel_cols = 1

            for attribute in layer.attribute:
                if attribute.name == "pads":
                    padding_rows = attribute.ints[0]
                    assert padding_rows == attribute.ints[1]
                    padding_cols = attribute.ints[2]
                    assert padding_cols == attribute.ints[3]
                if attribute.name == "strides":
                    stride_rows = attribute.ints[0]
                    stride_cols = attribute.ints[1]
                if attribute.name == "kernel_shape":
                    kernel_rows = attribute.ints[0]
                    kernel_cols = attribute.ints[1]
            print(inputShape)
            assert inputShape[0] == 1
            outputShape = inputShape[:-2]
            outputShape.append((inputShape[-2] + 2 * padding_rows - kernel_rows + stride_rows) // stride_rows)
            outputShape.append((inputShape[-1] + 2 * padding_cols - kernel_cols + stride_cols) // stride_cols)
            outputIndex = nn.add_maxpool_2d([kernel_rows,kernel_cols],inputShape[-2],inputShape[-1],inputShape[-3],[stride_rows,stride_cols],[padding_rows,padding_cols],inputIndex)
            return outputIndex, outputShape
        if (layer.op_type == "Relu"):
            inputIndex, inputShape = getLayer(layer.input[0])
            return nn.add_relu(inputIndex, useAreaHeuristic),inputShape

        if (layer.op_type == "Reshape"):
            inputIndex, inputShape = getLayer(layer.input[0])
            dataIndex, data = getLayer(layer.input[1])
            assert dataIndex == -1
            outputShape = data.tolist()
            if inputIndex == -1:
                return -1, np.reshape(inputShape, outputShape)
            inProd = 1
            for i in inputShape:
                inProd = inProd * i
            outProd = 1
            index = None
            for i in range(len(outputShape)):
                if outputShape[i] == -1:
                    assert index is None
                    index = i
                else:
                    outProd = outProd * outputShape[i]
            if not index is None:
                data[index] = inProd // outProd
            return inputIndex, data.tolist()
        if (layer.op_type == "Shape"):
            inputIndex, inputShape = getLayer(layer.input[0])
            assert inputIndex > 0
            return -1, np.array(inputShape)
        if (layer.op_type == "Unsqueeze"):
            inputIndex, input = getLayer(layer.input[0])
            assert inputIndex == -1
            axes = list(layer.attribute[0].ints)
            return -1, np.expand_dims(input, axes)



        if (layer.op_type == "Transpose"):
            #inputIndex, inputShape=getLayer(layer.input[0])
            #print(f"({inputIndex}) <- Transpose({inputIndex}) : {inputShape}")
            #return inputIndex, inputShape
            raise NotImplementedError("Transpose is not supported.")
        raise NotImplementedError(f"Unknown layer: {layer.op_type}")

        # raise


    #assert (len(graph.output) == 1)

    getLayer(graph.output[0].name)
    return nn
