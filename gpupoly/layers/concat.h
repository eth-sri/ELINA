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
   \file src/layers/concat.h
   \brief Concatenation layer.
   \author Fran&ccedil;ois Serre

   Neural network layer that concatenate its two inputs.
  */


#pragma once
#include "../network.h"

//! Concatenation layer.
/*!
  Neural network layer that concatenate its two inputs.
*/
class Concat :public NeuralNetwork::Layer
{
	int parentUp; //!< First parent
	int sizeUp; //!< Size of the first parent
	int parentDown; //!< Second parent
	int sizeDown; //!< Size of the second parent
public:
	//! Constructor
	/*!
	  Constructs a new bias layer.

	  \param sizeUp Size of the first input.
	  \param parentUp Index of the first parent layer (or 0 for the input layer).
	  \param sizeDown Size of the second input.
	  \param parentDown Index of the first parent layer (or 0 for the input layer).
	*/
	Concat(NeuralNetwork& nn, int sizeUp, int parentUp, int sizeDown, int parentDown) :NeuralNetwork::Layer(nn, sizeUp+sizeDown),sizeUp(sizeUp),parentUp(parentUp),sizeDown(sizeDown),parentDown(parentDown) {}

	template<typename T>
	void eval(Vector<T>& dest)
	{
		const Vector<T>& up = nn.getConcreteBounds<T>(parentUp);
		const Vector<T>& down = nn.getConcreteBounds<T>(parentDown);

		gpuErrchk(cudaMemcpy(
			(Intv<T>*)dest,
			(const Intv<T>*)up,
			up.memSize(),
			cudaMemcpyDeviceToDevice));
		gpuErrchk(cudaMemcpy(
			(Intv<T>*)dest + sizeUp,
			(const Intv<T>*)down,
			down.memSize(),
			cudaMemcpyDeviceToDevice));
	}

	virtual void eval(Vector<float>& dest, bool sound, bool precise)
	{
		eval(dest);
	}
	virtual void eval(Vector<double>& dest, bool sound, bool precise)
	{
		eval(dest);
	}
		
	template<typename T>
	void backSubstitute(typename AffineExpr<T>::Queue& queue, const AffineExpr<T>& expr) const
	{
		auto AUp = std::make_shared<Matrix<T>>(expr.A->m(), sizeUp, expr.A->interval());
		cudaMemcpy2D(
			AUp->data(), AUp->pitchBytes(),
			expr.A->data(), expr.A->pitchBytes(),
			sizeUp * sizeof(T) * (1 + expr.A->interval()), expr.A->m(),
			cudaMemcpyDeviceToDevice
		);
		queue.emplace(expr.m, sizeUp, parentUp, expr.up, expr.rows, AUp, expr.b,ConvShape(),expr.sound);
		auto ADown = std::make_shared<Matrix<T>>(expr.A->m(), sizeDown, expr.A->interval());
		cudaMemcpy2D(
			ADown->data(), ADown->pitchBytes(),
			(T*)expr.A->data()+sizeUp* (1 + expr.A->interval()), expr.A->pitchBytes(),
			sizeDown * sizeof(T) * (1 + expr.A->interval()), expr.A->m(),
			cudaMemcpyDeviceToDevice
		);
		queue.emplace(expr.m, sizeDown, parentDown, expr.up, expr.rows, ADown, nullptr, ConvShape(), expr.sound);
	}
	virtual void backSubstitute(typename AffineExpr<double>::Queue& queue, const AffineExpr<double>& expr) const
	{
		backSubstitute(queue, expr);
	}
	virtual void backSubstitute(typename AffineExpr<float>::Queue& queue, const AffineExpr<float>& expr) const
	{
		backSubstitute(queue, expr);
	}

};