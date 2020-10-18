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
  \file src/filters.h
  \brief Neuron filter functors.
  \author Fran&ccedil;ois Serre

  Definition of the class Filter that represents a neuron filter functor, and definition of some of them. Implementation of their method can be found in src/filters.cu.
 */

#pragma once
#include "vector.h"



 //! Functor to filter out neurons during backsubstitution.
 /*!
 Functor to filter out neurons during backsubstitution.
 */
template<typename T>
class NeuronFilter
{
public:
	//! Dtor
	virtual ~NeuronFilter() {}

	//! Provide a list of neurons of the previous layer which precision will impact the quality of the model, and that should therefore have their interval evaluated thoroughly via backsubstitution.
	virtual int listCriticalNeurons(int* dest, const Vector<T>& v, int* tmpInt, const int* oldList = nullptr, const int oldNbCritical = 0) const = 0;
};


//! Filter that selects neurons that contain 0 in their interval.
/*!
  Filter that selects the neurons that contain 0 within their concrete bounds.
*/
template<typename T>
class ContainsZero :public NeuronFilter<T> {
public:
	// Dtor
	virtual ~ContainsZero() {}

	// Provide a list of neurons of the previous layer which precision will impact the quality of the model, and that should therefore have their interval evaluated thoroughly via backsubstitution.
	virtual int listCriticalNeurons(int* dest, const Vector<T>& v, int* tmpInt, const int* oldList = nullptr, const int oldNbCritical = 0) const;
};

//! Filter that selects all neurons.
/*!
  Filter that selects all neurons.
*/
template<typename T>
class AlwaysKeep :public NeuronFilter<T> {
public:
	// Dtor
	virtual ~AlwaysKeep() {}

	// Provide a list of neurons of the previous layer which precision will impact the quality of the model, and that should therefore have their interval evaluated thoroughly via backsubstitution.
	virtual int listCriticalNeurons(int* dest, const Vector<T>& v, int* tmpInt, const int* oldList = nullptr, const int oldNbCritical = 0) const;
};