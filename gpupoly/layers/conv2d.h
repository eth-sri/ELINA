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
   \file src/layers/conv2d.h
   \brief Convolution layer.
   \author Fran&ccedil;ois Serre

   Neural network layer that performs a 2D convolution of its input.
  */



#pragma once
#include "../network.h"

//! 2D - Convolution layer. 
/*!
   Neural network layer that performs a 2D convolution of its input.
*/
template<typename T>
class Conv2D :public NeuralNetwork::Layer
{
	const int parent; //!< Index of parent layer.
	const Matrix<T> conv; //!< Convolution data (coefficients).
	const Matrix<T> convt;
#ifdef STRONG_FP_SOUNDNESS
	std::shared_ptr<const Matrix<float>> convf;
	std::shared_ptr<const Matrix<float>> convtf;
#endif

	const ConvShape cs; //!< Shape of the convolution

	const size_t inputSize; //!< input size (product of input dimentions)
	

public:
	//! Constructor
/*!
  Creates a new convolution layer, without activation nor bias.

  \param filters number of filters
  \param kernel_size_rows dimention of the kernel (rows)
  \param kernel_size_cols dimention of the kernel (cols)
  \param input_rows dimention of the input (number of rows)
  \param input_cols dimention of the input (number of cols)
  \param input_channels dimention of the input (number of channels)
  \param stride_rows stride shape (number of rows)
  \param stride_cols stride shape (number of cols)
  \param padding_rows padding (number of pixels to add at the top and bottom)
  \param padding_cols padding (number of pixels to add on the left and right)
  \param A Convolution data (Vector of coefficients given in row major, then column, then channel and filter minor order). Contains filters*kernel_size_rows*kernel_size_cols*input_shape_channels elements.)
  \param parent index of the parent layer (or 0 for the input layer)
 */
	Conv2D(
		NeuralNetwork& nn,
		const int filters,
		const int kernel_size_rows, const int kernel_size_cols,
		const int input_rows, const int input_cols, const int input_channels,
		const int stride_rows, const int stride_cols,
		const int padding_rows, const int padding_cols,
		const Matrix<T>& A,
		int parent);
	template <typename Te> void eval(Vector<Te>& dest, bool sound);
	//template <typename Te> void backSubstitute(typename AffineExpr<Te>::Queue& queue, const AffineExpr<Te>& expr) const;
	virtual void eval(Vector<double>& dest, bool sound, bool precise);
	virtual void eval(Vector<float>& dest, bool sound, bool precise);
	virtual void backSubstitute(typename AffineExpr<double>::Queue& queue, const AffineExpr<double>& expr) const;
	virtual void backSubstitute(typename AffineExpr<float>::Queue& queue, const AffineExpr<float>& expr) const;
};
