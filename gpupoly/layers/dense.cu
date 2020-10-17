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
   \file src/layers/dense.cu
   \brief Dense linear layer implementation.
   \author Fran&ccedil;ois Serre

   Implementation of class Dense, a neural network layer that performs a linear transform (a.k.a. a matrix-vector multiplication) defined in src/dense.h.
*/


#include "dense.h"


template <typename T>
Dense<T>::Dense(NeuralNetwork& nn, const Matrix<T>& A, const int parent) :NeuralNetwork::Layer(nn, A.m()), parent(parent), A(A) {}

template <typename T>
void Dense<T>::eval(Vector<double>& dest, bool sound, bool precise)
{
	A.mvm(dest, nn.template getConcreteBounds<double>(parent)); // a simple matrix-vector multiplication
}
template <typename T>
void Dense<T>::eval(Vector<float>& dest, bool sound, bool precise)
{
	A.mvm(dest, nn.template getConcreteBounds<float>(parent)); // a simple matrix-vector multiplication
}

template <typename T>
void Dense<T>::backSubstitute(typename AffineExpr<double>::Queue& queue, const AffineExpr<double>& expr) const
{
	std::shared_ptr<Matrix<double>> nA;
	if (expr.A)
	{
		nA = std::make_shared<Matrix<double>>();
		expr.A->mmm(A, *nA, expr.sound);
	}
	else
		nA = std::make_shared<Matrix<double>>(A.template selectRows<double>(expr.m, expr.rows));
	queue.emplace(expr.m, A.n(), parent, expr.up, expr.rows, nA, expr.b, ConvShape(), expr.sound);
}
template <typename T>
void Dense<T>::backSubstitute(typename AffineExpr<float>::Queue& queue, const AffineExpr<float>& expr) const
{
	std::shared_ptr<Matrix<float>> nA;
	if (expr.A)
	{
		nA = std::make_shared<Matrix<float>>();
		expr.A->mmm(A, *nA, expr.sound);
	}
	else
		nA = std::make_shared<Matrix<float>>(A.template selectRows<float>(expr.m, expr.rows));
	queue.emplace(expr.m, A.n(), parent, expr.up, expr.rows, nA, expr.b, ConvShape(), expr.sound);
}

#ifdef STRONG_FP_SOUNDNESS
template <> Dense<double>::Dense(NeuralNetwork& nn, const Matrix<double>& A, const int parent) :NeuralNetwork::Layer(nn, A.m()), parent(parent), A(A),Af(std::make_shared<const Matrix<float>>(A,false)) {}
template <> void Dense<double>::eval(Vector<float>& dest, bool sound, bool precise)
{
	Af->mvm(dest, nn.template getConcreteBounds<float>(parent)); // a simple matrix-vector multiplication
}
template <>
void Dense<double>::backSubstitute(typename AffineExpr<float>::Queue& queue, const AffineExpr<float>& expr) const
{
	std::shared_ptr<Matrix<float>> nA;
	if (expr.A)
	{
		nA = std::make_shared<Matrix<float>>();
		expr.A->mmm(*Af, *nA, expr.sound);
	}
	else
		nA = std::make_shared<Matrix<float>>(Af->template selectRows<float>(expr.m, expr.rows));
	queue.emplace(expr.m, Af->n(), parent, expr.up, expr.rows, nA, expr.b, ConvShape(), expr.sound);
}
#endif


template class Dense<double>;
template class Dense<float>;
