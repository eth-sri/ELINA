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
   \file src/layers/bias.h
   \brief Bias layer.
   \author Fran&ccedil;ois Serre

   Neural network layer that adds a vector a constants (bias) to its inputs.
  */

#pragma once

#include "../network.h"

//! Bias layer.
/*!
 Layer that adds a constant vector to its inputs.
  */
template<typename T>
class Bias :public NeuralNetwork::Layer
{
	const int parent; //!< index of the parent layer
	std::shared_ptr<const Vector<T>> b; //!< Constant Vector to be added
public:
	//! Constructor
	/*!
	  Constructs a new bias layer.

	  \param b Pointer to a Vector containing the bias coefficients. This Vector may contain intervals.
	  \param parent index of the parent layer (or 0 for the input layer).
	*/
	Bias(NeuralNetwork& nn,std::shared_ptr<const Vector<T>> b, const int parent) :NeuralNetwork::Layer(nn, b->size()), b(b), parent(parent) {}

	virtual void eval(Vector<double>& dest, bool sound, bool precise)
	{
		Vector<double>::add(dest, nn.template getConcreteBounds<double>(parent), *b);
	}
	virtual void eval(Vector<float>& dest, bool sound, bool precise)
	{
		Vector<float>::add(dest, nn.template getConcreteBounds<float>(parent), *b);
	}


	virtual void backSubstitute(typename AffineExpr<double>::Queue& queue, const AffineExpr<double>& expr) const
	{
		auto nb = std::make_shared<Vector<double>>(expr.evaluate(*b));
		queue.emplace(expr.m, expr.n, parent, expr.up, expr.rows, expr.A, nb, expr.cs,expr.sound);
	}
	virtual void backSubstitute(typename AffineExpr<float>::Queue& queue, const AffineExpr<float>& expr) const
	{
		auto nb = std::make_shared<Vector<float>>(expr.evaluate(*b));
		queue.emplace(expr.m, expr.n, parent, expr.up, expr.rows, expr.A, nb, expr.cs, expr.sound);
	}

};