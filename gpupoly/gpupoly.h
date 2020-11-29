/*
 *  GPUPoly library
 *  This source file is part of ELINA (ETH LIbrary for Numerical Analysis).
 *  ELINA is Copyright © 2020 Department of Computer Science, ETH Zurich
 *  This software is distributed under GNU Lesser General Public License Version 3.0.
 *  For more information, see the ELINA project website at:
 *  http://elina.ethz.ch
 *
 *  THE SOFTWARE IS PROVIDED "AS-IS" WITHOUT ANY WARRANTY OF ANY KIND, EITHER
 *  EXPRESS, IMPLIED OR STATUTORY, INCLUDING BUT NOT LIMITED TO ANY WARRANTY
 *  THAT THE SOFTWARE WILL CONFORM TO SPECIFICATIONS OR BE ERROR-FREE AND ANY
 *  IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE,
 *  TITLE, OR NON-INFRINGEMENT.  IN NO EVENT SHALL ETH ZURICH BE LIABLE FOR ANY
 *  DAMAGES, INCLUDING BUT NOT LIMITED TO DIRECT, INDIRECT,
 *  SPECIAL OR CONSEQUENTIAL DAMAGES, ARISING OUT OF, RESULTING FROM, OR IN
 *  ANY WAY CONNECTED WITH THIS SOFTWARE (WHETHER OR NOT BASED UPON WARRANTY,
 *  CONTRACT, TORT OR OTHERWISE).
 */

 /*!
   \file include/gpupoly.h
   \brief API for the GPUPoly library.
   \author Fran&ccedil;ois Serre

   This file contains definitions of the functions to be used to interact with the library. Their implementation can be found in src/bindings.cu.
   An example program that uses this library is provided in test/example.cpp.
  */
#pragma once

#include "gpupoly_export.h"


  // Forward declaration of the class NeuralNetwork, wich represents a pre-trained neural network. In this API, a pointer to it is used as a handle.
#ifdef __cplusplus
	class NeuralNetwork; 
#else
typedef struct NeuralNetwork NeuralNetwork; // pure C version.
#endif

#ifdef __cplusplus
extern "C" {
#endif
	//! Creates a new network.
	/*!
	  \param inputSize Number of elements of the inputs (e.g. 784 for MNIST dataset)
	  \returns A handle to be used by all other API calls to refer to this network.
	*/
	GPUPOLY_EXPORT NeuralNetwork* create(
		int inputSize
	);
	/*! \name Verifies an image
	  Checks if an image (or more precisely, an input box) entirely classifies as a given label
	  Once this function has been called with a network, no additionnal layer should be added (i.e. it is necessary to clean and recreate a network)
	  It is of course possible to call this function multiple times on the same network.

	  Depending on GPUPoly configuration, this function can deliver different certification levels, listed here in increasing strength order:
		0. No certification: GPUPoly was not able to prove that the whole input box classifies as the given label.
		1. Certification without floating-point soundness. This yields no formal guarantee, but is in practice likely to mean 2.
		2. Certified for an infinitly precise evaluation. If a box is certified with this level, any element of the box will classify to the given label, as long as the inference is computed using infinitly precise arithmetic.
		3. Certified for an evaluation that uses double precision IEEE754 arithmetic, with any summation order and any rounding mode.
		4. Certified for an evaluation that uses single precision IEEE754 arithmetic, with any summation order and any rounding mode.

	  By default, GPUPoly tries to deliver level 2 certifications, and will therefore output a result that will be either 0 or 2. The parameter soundness allows to disable floating-point soundness to speed up the computation (which then becomes roughly twice or three times faster). This mode doesn't offer any formal guarantee anymore, but yields in practice the same results. Without floating point soundness, GPUPoly outputs certification levels that can be 0 or 1.
	  Conversely, it is possible to make GPUPoly search for stronger certification levels using the `STRONG_FP_SOUNDNESS` compile option. If GPUPoly is compiled with this option, the letter `S` will appear after its version number (instead of `W`), and it will output certification levels of 0, 3 or 4. If the floating-point soundness parameter is disabled, GPUPoly has the same behaviour with or without the `STRONG_FP_SOUNDNESS` option.
	  */

      //!@{

	//! Tries to certify a box.
	/*!
	  \param nn Handle to the network
	  \param dataDown An array of inputSize floats that represent the lower bound of the box
	  \param dataUp An array of inputSize floats that represent the upper bound of the box
	  \param expectedLabel Label in which the image is supposed to classify
	  \param soundness Whether to use sound arithmetic.
	  \returns if certified.
	 */
	GPUPOLY_EXPORT bool test_s(
		NeuralNetwork* nn,
		const float* dataDown,
		const float* dataUp,
		int expectedLabel,
		bool soundness = true
	);
	//! Creates a new network.
	/*!
	  \param nn Handle to the network
	  \param dataDown An array of inputSize doubles that represent the lower bound of the box
	  \param dataUp An array of inputSize doubles that represent the upper bound of the box
	  \param expectedLabel Label in which the image is supposed to classify
	  \param soundness Whether to use sound arithmetic.
	  \returns if certified.
	 */
	GPUPOLY_EXPORT bool test_d(
		NeuralNetwork* nn,
		const double* dataDown,
		const double* dataUp,
		int expectedLabel,
		bool soundness = true
	);
	//!\}


	//! Sets the concrete bounds of a layer.
	//!@{
	//! Sets the single precision concrete bounds of a layer.
	/*!
	  \param nn Handle to the network
	  \param dataDown An array of inputSize floats that represent the lower bound of the box
	  \param dataUp An array of inputSize floats that represent the upper bound of the box
	  \param layer Index of the layer
	 */
	GPUPOLY_EXPORT void setLayerBox_s(
		NeuralNetwork* nn,
		const float* dataDown,
		const float* dataUp,
		int layer
	);
	//! Sets the double precision concrete bounds of a layer.
	/*!
	  \param nn Handle to the network
	  \param dataDown An array of inputSize doubles that represent the lower bound of the box
	  \param dataUp An array of inputSize doubles that represent the upper bound of the box
	  \param layer Index of the layer
	 */
	GPUPOLY_EXPORT void setLayerBox_d(
		NeuralNetwork* nn,
		const double* dataDown,
		const double* dataUp,
		int layer
	);
	//!\}


	//! Propagates forward the concrete bounds by interval analysis. Activation layers have their approximation models computed.
	//!@{
	//! Propagates forward single precision bounds.
	/*!
	  \param nn Handle to the network
	  \param layer Index of the layer
	  \param refineActivationsInput If true, and layer is an activation layer, then the input is first refined first via the appropriate back-substitution.
	  \param soundness Whether to use sound (but slower) arithmetic.
	 */
	GPUPOLY_EXPORT void relax_s(
		NeuralNetwork* nn, 
		int layer, 
		bool refineActivationsInput,
		bool soundness
	);
	//! Propagates forward double precision bounds.
	/*!
	  \param nn Handle to the network
	  \param layer Index of the layer
	  \param refineActivationsInput If true, and layer is an activation layer, then the input is first refined first via the appropriate back-substitution.
	  \param soundness Whether to use sound (but slower) arithmetic.
	*/
	GPUPOLY_EXPORT void relax_d(
		NeuralNetwork* nn,
		int layer,
		bool refineActivationsInput,
		bool soundness
	);
	//!\}

	//! Evaluate the concrete bounds of a list of affine expressions.
	/*!
		Evaluate the concrete bounds of a list of m affine expressions of the neurons of a given layer via back-substitution, using the concrete bounds of the last box that was tested.
		The	affine expressions have the form Ax+b, where A is a m*n matrix, b a vector of size m, and x represents the n neurons of the layer layerId.
	*/
	//!@{
	//! Evaluate single precision expressions.
	/*!
	  \param nn Handle to the network
	  \param dest An array of 2*m floats where the (interleaved) result will be written.
	  \param layer Index of the layer
	  \param m Number of expressions
	  \param A Content of the matrix A in row major order. A has m rows, and its number of columns equals the outputSize of the layer.
	  \param b Content of the vector b. b has m elements.
	  \param backsubstitute If 0, only evaluate with the concrete bounds of layer. If 1, always backsubstitute back to the inputs. If 2, backsubstitute until 0 is not strictly included within the bounds.
	  \param soundness Whether to use sound (but slower) arithmetic.
	 */
	GPUPOLY_EXPORT void evalAffineExpr_s(
		NeuralNetwork* nn,
		float* dest,
		int layer,
		int m,
		const float* A,
		const float* b,
		int backsubstitute,
		bool soundness
	);
	//! Evaluate double precision expressions.
	/*!
	  \param nn Handle to the network
	  \param dest An array of 2*m doubles where the (interleaved) result will be written.
	  \param layer Index of the layer
	  \param m Number of expressions
	  \param A Content of the matrix A in row major order. A has m rows, and its number of columns equals the outputSize of the layer.
	  \param b Content of the vector b. b has m elements.
	  \param backsubstitute If 0, only evaluate with the concrete bounds of layer. If 1, always backsubstitute back to the inputs. If 2, backsubstitute until 0 is not strictly included within the bounds.
	  \param soundness Whether to use sound (but slower) arithmetic.
	 */
	GPUPOLY_EXPORT void evalAffineExpr_d(
		NeuralNetwork* nn,
		double* dest,
		int layer,
		int m,
		const double* A,
		const double* b,
		int backsubstitute,
		bool soundness
	);
	//!\}

	//! Gets the dimention of a given layer.
	GPUPOLY_EXPORT int getOutputSize(NeuralNetwork* nn, int layer);

	//! Cleans the network
	/*!
	  Cleans the network from memory.
	  This function has to be called at the end of the program, or before specifying a new network (otherwise create will return false).
	*/
	GPUPOLY_EXPORT void clean(
		NeuralNetwork* nn
	);




	/*! \name Adding layers to a network
	  The following functions add an additional layer to the network
	  They all return an internal index of the layer.
	  They all have one (or two) "parent" parameter, where the index of a previous layer should be specified (or 0 for the input layer).
	  Note that layer added successively might not have consecutive indexes.
	*/
	//!@{

	// Affine layers:


	//! Linear layer
	/*!
	  Adds a dense linear layer (without activation nor bias). If x is a column vector representing the input of the layer, the output is A*x.

	  \param parent index of the parent layer (or 0 for the input layer)
	  \param outputSize Size of the output (a.k.a. the number of neurons of the layer)
	  \param data Content of the matrix A in row major order. A has outputSize rows, and its number of columns equals the parent outputSize.
	  \returns the index of the newly created layer.
	 */
	GPUPOLY_EXPORT int addLinear_d(
		NeuralNetwork* nn,
		int parent,
		int outputSize,
		const double* data
	);
	GPUPOLY_EXPORT int addLinear_s(
		NeuralNetwork* nn,
		int parent,
		int outputSize,
		const float* data
	);
	
	//! Bias layer
	/*!
	  Adds a bias layer, i.e. a layer that adds a constant vector to its input.

	  \param parent index of the parent layer (or 0 for the input layer)
	  \param data Constant vector to be added to the input. Contains as many elements as the input.
	  \returns the index of the newly created layer.
	 */
	GPUPOLY_EXPORT int addBias_s(
		NeuralNetwork* nn,
		int parent,
		const float* data
	);
	GPUPOLY_EXPORT int addBias_d(
		NeuralNetwork* nn,
		int parent,
		const double* data
	);

	//! Convolution layer
	/*!
	  Adds a convolution layer, without activation.

	  \param parent index of the parent layer (or 0 for the input layer).
	  \param channels_first if true, the layer expects input with the shape [batch, channel, row, col], and its output has the shape [batch, filter, row, col]. If false, the layer expects input with the shape [batch, row, col, channel], and its output has the shape [batch, row, col, filter].
	  \param filters number of filters.
	  \param kernel_shape dimentions of the kernel (expects an array of 2 ints, respectively the number of rows and columns.
	  \param input_shape dimentions of the input (expects an array of 4 ints, respectively the number of batches, rows, columns and channels).
	  \param stride_shape stride shape (expects an array of 2 ints, respectively the number of rows and columns.
	  \param padding padding (expects an array of 2 ints, respectively the number of pixels to add at the top and bottom, and the number of pixels to add on the left and right).
	  \param data convolution coefficients (given in row major, then column, then channel and filter minor order). Contains filters*kernel_size_rows*kernel_size_cols*input_shape_channels elements.
	  \returns the index of the newly created layer.
	 */
	GPUPOLY_EXPORT int addConv2D_d(
		NeuralNetwork* nn,
		int parent,
		bool channels_first,
		int filters,
		int* kernel_shape,
		int* input_shape,
		int* stride_shape,
		int* padding,
		const double* data
	);

	GPUPOLY_EXPORT int addConv2D_s(
		NeuralNetwork* nn,
		int parent,
		bool channels_first,
		int filters,
		int* kernel_shape,
		int* input_shape,
		int* stride_shape,
		int* padding,
		const float* data
	);

	// Activation layers:

	//! ReLU layer
	/*!
	  Adds a ReLU layer to the network.

	  \param parent index of the parent layer (or 0 for the input layer)
	  \returns the index of the newly created layer.
	 */
	GPUPOLY_EXPORT int addReLU(
		NeuralNetwork* nn,
		int parent
	);

	//! MaxPool2D layer.
	/*!
	  Adds a max pooling layer.

	 \param parent index of the parent layer (or 0 for the input layer)
	 \param channels_first if true, the layer expects input with the shape [batch, channel, row, col], and its output has the shape [batch, channel, row, col]. If false, the layer expects input with the shape [batch, row, col, channel], and its output has the shape [batch, row, col, channel].
	 \param pool_shape pool shape (expects an array of 2 ints, respectively the number of rows and columns
	 \param input_shape dimentions of the input (expects an array of 4 ints, respectively the number of batches, rows, columns and channels).
	 \param stride_shape stride shape (expects an array of 2 ints, respectively the number of rows and columns
	 \param padding padding (expects an array of 2 ints, respectively the number of pixels to add at the top and bottom, and the number of pixels to add on the left and right)
	 \returns the index of the newly created layer.
	 */
	GPUPOLY_EXPORT int addMaxPool2D(
		NeuralNetwork* nn,
		int parent,
		bool channels_first,
		int* pool_shape,
		int* input_shape,
		int* stride_shape,
		int* padding
	);

	// Residual layers:

	//! ParSum layer.
	/*!
	  Adds a "ParSum" layer, i.e. a layer that sums up the result of two previous layers.
	  \param parent1 Index of the first parent layer.
	  \param parent2 Index of the second parent layer.
	  \returns the index of the newly created layer.
	 */
	GPUPOLY_EXPORT int addParSum(
		NeuralNetwork* nn,
		int parent1,
		int parent2
	);

	// 
	//! Concatenation layer.
	/*!
	  Adds a concatenation layer (like the one we can find in skipnets). It concatenates the result of two previous layers.
	  \param parent1 Index of the first parent layer.
	  \param parent2 Index of the second parent layer.
	  \returns the index of the newly created layer.
	 */
	GPUPOLY_EXPORT int addConcat(
		NeuralNetwork* nn,
		int parent1,
		int parent2
	);

	//!\}

#ifdef __cplusplus
}
#endif