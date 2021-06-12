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
   \file src/network.h
   \brief NeuralNetwork definition.
   \author Fran&ccedil;ois Serre

   Definition of the class NeuralNetwork that represents a neural network. It contains the main verification algorithm.
  */

#pragma once
#include "affineexpr.h"
#include "matrix.h"
#include "vector.h"
#include "filters.h"
#include <vector>
#include <memory>

class NeuralNetwork
{
public:


	//! Neural network layer
/*!
Interface for the layers of a neural network
*/
	class Layer
	{
	protected:
		NeuralNetwork& nn;
	public:
		const size_t outputSize; 	//!< number of neurons, i.e. number of elements of the output vector

		//! Constructor
		/*!
		  Constructs a new bias layer.
			\param nn Reference to the network.
		  \param outputSize Number of neurons (values) that this layer outputs.
		*/
		Layer(NeuralNetwork& nn, const size_t outputSize) : nn(nn), outputSize(outputSize) {}

		//! Computes the output value using the results of parent layers.
		/*!
		  Computes the output interval of a layer via forward evaluation.
		  Activation layers may order a backsubstitution of their input first to increase accuracy, and have to compute a model (stored internally) for future back-substitutions.

		  \param dest Vector where to put the result.
		*/
		virtual void eval(Vector<float>& dest, bool sound, bool precise) = 0;
		virtual void eval(Vector<double>& dest, bool sound, bool precise) = 0;

		//! Performs a backsubstitution 
		/*!
		  Given an affine expression (expr) of the neurons of this layer, express it as a sum of affine expressions of parent (previous) layers, and store these expressions in queue.

		  \param queue Reference to a queue where the backsubstituted expression(s) are to be pushed.
		  \param expr An affine expression of the neurons of this layer.
		*/
		virtual void backSubstitute(typename AffineExpr<double>::Queue& queue, const AffineExpr<double>& expr) const = 0;
		virtual void backSubstitute(typename AffineExpr<float>::Queue& queue, const AffineExpr<float>& expr) const = 0;

		//! Destructor
		virtual ~Layer() {}

		virtual bool hasInternalModel() const { return false; }
	};


	//! Constructor
	/*!
	  Constructs a new (pre-trained) network to be checked

	  \param inputSize Number of elements of the input layer
	*/
	NeuralNetwork(const size_t inputSize);

	NeuralNetwork(const NeuralNetwork&) = delete;

	~NeuralNetwork();

	
	//! Adds a new layer in the network
	/*!
	  Constructs a new neural network to be checked

	  \param layer A layer, as prototyped earlier in this file.
	  \returns The layer id that was assigned to this layer.
	*/
	int addLayer(Layer* layer);

	Layer* operator[](const size_t i) const {
		return layers[i];
	}

	//! Tries to verify a box.
	/*!
	  \param input A vector containing the interval for each input element
	  \param label Label in which the image is supposed to classify
	  \param sound Whether to use sound arithmetic.
	  \returns An integer indicating the certification level (see gpupoly.h for details).
	 */
	template <typename T>
	bool run(const Vector<T>& input, const int label, bool sound);


	//! Get a reference to the current concrete bounds of a layer
	template <typename T>
	Vector<T>& getConcreteBounds(int layer);


	//! Evaluates an affine expression
	/*!
	  Evaluates affine combinations of neurons at a layer of a neural network by backsubstitution to the inputs.
	  These expressions have the form Ax+b, where x is the value of the neurons at layer fromLayer
	*/
	template<typename T>
	void evaluateAffine(Vector<T>& dest, const NeuronFilter<T>& al, int layer, bool up, bool sound, const std::shared_ptr<const Matrix<T>>& A = nullptr, const std::shared_ptr<const Vector<T>>& b = nullptr);

	template<typename T>
	AffineExpr<T> getSensitivityExpr(T* const destA, T* const destb, int layer, bool up, bool sound, const std::shared_ptr<const Matrix<T>>& A = nullptr, const std::shared_ptr<const Vector<T>>& b = nullptr);
	
	//! Orders a reevaluation (via back-substitution) of the concrete bounds of a layer
	/*!
	  \param layer Id of the layer to be reevaluated
	  \param al A NeuronFilter indicating a stopping criteria for the back-substitution
	  \param up If true, reevaluates the upper bound, otherwise the lower bound.
	  \param sound Whether to use sound arithmetic.
	 */
	template <typename T>
	void reEvaluateLayer(int layer, const NeuronFilter<T>& al, bool up, bool sound)
	{
		evaluateAffine(getConcreteBounds<T>(layer), al, layer, up, sound);
	}
private:



	std::vector<Layer*> layers;
	size_t maxLayerSize;

	std::vector<std::shared_ptr<Vector<double>>> concreteBoundsD;
	std::vector<std::shared_ptr<Vector<float>>> concreteBoundsS;



	int* annoyingNeurons;
	int* annoyingNeuronList;
	int* annoyingNeuronList2;




	
};

