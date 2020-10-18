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
   \file src/layers/parsum.h
   \brief Concatenation layer.
   \author Fran&ccedil;ois Serre

   Neural network layer that sums its two inputs.
  */

#pragma once
#include "../network.h"

// Layer that sums the result of two previous layers.
//! Partial sum layer.
/*!
  Neural network layer that sums its two inputs.
*/
class ParSum :public NeuralNetwork::Layer
{
	const int parentL; // The two 
	const int parentR; // parent layers

public:
	//! Constructor
/*!
  Constructs a new bias layer.

  \param size Size of the first input.
  \param parentUp Index of the first parent layer (or 0 for the input layer).
  \param sizeDown Size of the second input.
  \param parentDown Index of the first parent layer (or 0 for the input layer).
*/
	ParSum(NeuralNetwork& nn, int size, int parentL, int parentR) :NeuralNetwork::Layer(nn, size), parentL(parentL), parentR(parentR) {}
	virtual void eval(Vector<double>& dest, bool sound, bool precise)
	{
		assert(nn.getConcreteBounds<double>(parentL).size() == outputSize);
		assert(nn.getConcreteBounds<double>(parentR).size() == outputSize);
		Vector<double>::add(dest, nn.getConcreteBounds<double>(parentL), nn.getConcreteBounds<double>(parentR));
	}
	virtual void eval(Vector<float>& dest, bool sound, bool precise)
	{
		assert(nn.getConcreteBounds<float>(parentL).size() == outputSize);
		assert(nn.getConcreteBounds<float>(parentR).size() == outputSize);
		Vector<float>::add(dest, nn.getConcreteBounds<float>(parentL), nn.getConcreteBounds<float>(parentR));
	}
	virtual void backSubstitute(typename AffineExpr<double>::Queue& queue, const AffineExpr<double>& expr) const
	{
		queue.emplace(expr.m, expr.n, parentL, expr.up, expr.rows, expr.A, expr.b, expr.cs, expr.sound);
		queue.emplace(expr.m, expr.n, parentR, expr.up, expr.rows, expr.A, nullptr, expr.cs, expr.sound);
	}
	virtual void backSubstitute(typename AffineExpr<float>::Queue& queue, const AffineExpr<float>& expr) const
	{
		queue.emplace(expr.m, expr.n, parentL, expr.up, expr.rows, expr.A, expr.b, expr.cs, expr.sound);
		queue.emplace(expr.m, expr.n, parentR, expr.up, expr.rows, expr.A, nullptr, expr.cs, expr.sound);
	}
};