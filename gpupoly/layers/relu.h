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
   \file src/layers/relu.h
   \brief ReLU layer.
   \author Fran&ccedil;ois Serre

   Definition of the class ReLU.
  */


#pragma once
#include "../network.h"

//! ReLU layer
/*!
ReLU layer: perform max(., 0) elementwise on the result.
Back-substition model consists of one upper elementwise affine expression (activFac[i].high * x[i] + activCst[i].high), and one lower (activFac[i].low * x[i] + activCst[i].low).
*/
class ReLU : public NeuralNetwork::Layer {
	const int parent; //!< Index of the parent layer
	Vector<double> activCstD; //!< Model constants
	Vector<float> activCstS; //!< Model constants
	Vector<double> activFacD; //!< Model factors
	Vector<float> activFacS; //!< Model factors
	template<typename T> Vector<T>& activCst();
	template<typename T> Vector<T>& activFac();
	template<typename T> const Vector<T>& activCst() const;
	template<typename T> const Vector<T>& activFac() const;

public:
	//! Constructor
/*!
  Constructs a new bias layer.

  \param size Number of elements in the input (and therefore output).
  \param parent Index of the first parent layer (or 0 for the input layer).
*/
	ReLU(NeuralNetwork& nn, int size, const int parent) :NeuralNetwork::Layer(nn, size), activCstD(size, true), activCstS(size,true), activFacS(size, true), activFacD(size, true), parent(parent) {}
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

