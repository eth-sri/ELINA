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
  \file src/layers/input.h
  \brief Conv2D implementation.
  \author Fran&ccedil;ois Serre

  Implementation of the methods of the class Conv2D. Definition of this class can be found in src/layers/conv2d.h.
*/
#pragma once
#include "../network.h"

//! Input layer.
/*!
  Input layer... Does nothing, but is convenient for the main algorithm.
*/
class Input :public NeuralNetwork::Layer
{
public:
	//! Constructor
/*!
  Creates a new input layer.

  \param size Number of input element of the network.
 */
	Input(NeuralNetwork& nn, int size) :NeuralNetwork::Layer(nn,size) {}
	virtual void eval(Vector<double>& dest, bool sound, bool precise)
	{
		return;
	}
	virtual void eval(Vector<float>& dest, bool sound, bool precise)
	{
		return;
	}

	virtual void backSubstitute(typename AffineExpr<double>::Queue& queue, const AffineExpr<double>& expr) const
	{
		return;
	}
	virtual void backSubstitute(typename AffineExpr<float>::Queue& queue, const AffineExpr<float>& expr) const
	{
		return;
	}

};