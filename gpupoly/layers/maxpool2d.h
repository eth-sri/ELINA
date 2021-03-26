/*
 *  GPUPoly library
 *  This source file is part of ELINA (ETH LIbrary for Numerical Analysis).
 *  ELINA is Copyright � 2020 Department of Computer Science, ETH Zurich
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
	\file src/layers/maxpool2d.h
	\brief MaxPool2D layer.
	\author Fran&ccedil;ois Serre

	Neural network layer that performs a 2D maxpool of its input.
   */
#pragma once
#include "../network.h"

class MaxPool2D :public NeuralNetwork::Layer {
	const int parent; //!< Index of parent layer.
	const ConvShape cs; //!< Shape of the transform
	Vector<float> modelCstS; //!< Model constants
	Vector<double> modelCstD; //!< Model constants
	Vector<float> modelFacS; //!< Model factors
	Vector<double> modelFacD; //!< Model factors
	template<typename T> Vector<T>& modelCst();
	template<typename T> Vector<T>& modelFac();
	template<typename T> const Vector<T>& modelCst() const;
	template<typename T> const Vector<T>& modelFac() const;
	size_t* modelNeuron; //!< Which input neuron is activated
public:
	//! Constructor
/*!
  Creates a new maxpool layer.
  \param pool_rows dimention of the kernel (rows)
  \param pool_cols dimention of the kernel (cols)
  \param input_rows dimention of the input (number of rows)
  \param input_cols dimention of the input (number of cols)
  \param input_channels dimention of the input (number of channels)
  \param stride_rows stride shape (number of rows)
  \param stride_cols stride shape (number of cols)
  \param padding_rows padding (number of pixels to add at the top and bottom)
  \param padding_cols padding (number of pixels to add on the left and right)
  \param parent index of the parent layer (or 0 for the input layer)
 */
	MaxPool2D(
		NeuralNetwork& nn,
		int pool_rows, int pool_cols,
		int input_rows, int input_cols, int input_channels,
		int stride_rows, int stride_cols,
		int padding_top, int padding_left,
        int padding_bottom, int padding_right,
        int parent);
	virtual ~MaxPool2D();
	template <typename T>
	void eval(Vector<T>& dest, bool sound, bool precise);
	template <typename T>
	void backSubstitute(typename AffineExpr<T>::Queue& queue, const AffineExpr<T>& expr) const;
	virtual void eval(Vector<double>& dest, bool sound, bool precise);
	virtual void eval(Vector<float>& dest, bool sound, bool precise);
	virtual void backSubstitute(typename AffineExpr<double>::Queue& queue, const AffineExpr<double>& expr) const;
	virtual void backSubstitute(typename AffineExpr<float>::Queue& queue, const AffineExpr<float>& expr) const;
	virtual bool hasInternalModel() const { return true; }
};