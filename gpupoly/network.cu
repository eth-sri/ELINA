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


 /*! \file src/network.cu
	  \brief Implementation of NeuralNetwork
	  \author Fran&ccedil;ois Serre

	  Implementation of the members of the class NeuralNetwork, defined in src/network.h.
  */

#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <vector>
#include <assert.h>
#include <cmath>
#include <queue>
#include <thrust/fill.h>
#include "layers/input.h"
#include "filters.h"
#include "network.h"
#include "cuda.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
using namespace std;


template <>
Vector<float>& NeuralNetwork::getConcreteBounds(int layer)
{
	return *concreteBoundsS[layer];
}

template <>
Vector<double>& NeuralNetwork::getConcreteBounds(int layer)
{
	return *concreteBoundsD[layer];
}

template <typename T>
__global__ void setFinal(T* A, int* rows, int label, size_t outputSize, size_t N)
{
	size_t i = threadIdx.x;
	rows[i] = i;
	for (int j = 0; j < outputSize; j++)
		if (j == i + (i >= label))
			A[i * N + j] = -1;
		else if (j == label)
			A[i * N + j] = 1;
		else
			A[i * N + j] = 0;
}


template <typename T>
void NeuralNetwork::evaluateAffine(Vector<T>& dest, const NeuronFilter<T>& al, int layer, bool up, bool sound, const std::shared_ptr<const Matrix<T>>& A, const std::shared_ptr<const Vector<T>>& b)
{
	if (!annoyingNeuronList)
		gpuErrchk(cudaMalloc((void**)&annoyingNeuronList, maxLayerSize * sizeof(int)));
	if (!annoyingNeuronList2)
		gpuErrchk(cudaMalloc((void**)&annoyingNeuronList2, maxLayerSize * sizeof(int)));

	int an = al.listCriticalNeurons(annoyingNeuronList, dest, annoyingNeurons);
	if (an)
	{
		int maxNeurBP = (1 << 30) / (maxLayerSize * sizeof(Intv<T>));
		for (int start = 0; start < an; start += maxNeurBP)
		{
			int length = std::min(maxNeurBP, an - start);
			auto partialA = A ? std::make_shared<const Matrix<T>>(A->template selectRows<T> (length, annoyingNeuronList + start)) : nullptr;
			auto partialb = b ? std::make_shared<const Vector<T>>(b->select(length, annoyingNeuronList+start)) : nullptr;
			auto inExpr = AffineExpr<T>(length, layers[layer]->outputSize, layer, up, annoyingNeuronList + start,partialA,partialb,ConvShape(),sound);
			typename AffineExpr<T>::Queue exprs;
			exprs.push(inExpr);
			int nbEval = 0;
			while (!exprs.empty())
			{
				AffineExpr<T> tmp = exprs.top();
				assert(tmp.sound == sound);
				exprs.pop();
				//assert(tmp.m == size);
				assert(tmp.n == layers[tmp.layer]->outputSize);
				if (exprs.empty())
				{
					//concreteBounds[tmp.layer]->check();
					tmp.evaluateAndUpdate(dest, getConcreteBounds<T>(tmp.layer));
					//dest.check();
					nbEval++;

					if (nbEval > 1)
					{
						int an = al.listCriticalNeurons(annoyingNeuronList2, dest, annoyingNeurons, tmp.rows, tmp.m);
						if (an < tmp.m)
						{
							if (an == 0)
								return;
							tmp.selectRows(an, annoyingNeuronList2);
						}
					}

				}
				if (!exprs.empty() && exprs.top().layer == tmp.layer)
				{
					AffineExpr<T> tmp2 = exprs.top();
					exprs.pop();
					assert(tmp.sound == tmp2.sound);
					auto A = std::make_shared<Matrix<T>>();
					Matrix<T>::add(*A, *tmp.getA(), *tmp2.getA(),tmp.sound);
					std::shared_ptr<const Vector<T>> b;
					if (tmp.up)
						b = Vector<T> ::template add_dr<true> (tmp.b, tmp2.b);
					else
						b = Vector<T> ::template add_dr<false> (tmp.b, tmp2.b);
					ConvShape cs;
					if (tmp.cs && tmp2.cs)
					{
						if (
							tmp.cs.filters == tmp2.cs.filters &&
							tmp.cs.output_rows == tmp2.cs.output_rows &&
							tmp.cs.output_cols == tmp2.cs.output_cols &&
							tmp.cs.input_rows == tmp2.cs.input_rows &&
							tmp.cs.input_cols == tmp2.cs.input_cols &&
							tmp.cs.input_channels == tmp2.cs.input_channels
							)
						{
							if (tmp.cs.kernel_size_cols > tmp2.cs.kernel_size_cols)
								cs = tmp.cs;
							else
								cs = tmp2.cs;
						}
						else
						{
							std::cout << "Error merging:" << std::endl;
							tmp.cs.print();
							tmp2.cs.print();
							std::cout << std::endl;
						}

					}
					exprs.emplace(tmp.m, tmp.n, tmp.layer, tmp.up, tmp.rows, A, b, cs,tmp.sound);
				}
				else
					layers[tmp.layer]->backSubstitute(exprs, tmp);
			}
		}
	}
}

NeuralNetwork::NeuralNetwork(const size_t inputSize) :
	layers(), maxLayerSize(0),
	concreteBoundsS(), concreteBoundsD(),
	annoyingNeuronList(nullptr),
	annoyingNeuronList2(nullptr)
{
	gpuErrchk(cudaMalloc((void**)&annoyingNeurons, sizeof(int)));
	addLayer(new Input(*this,inputSize));
}

NeuralNetwork::~NeuralNetwork() {
	if (annoyingNeurons) // we were not moved
	{
		cudaFree(annoyingNeurons);
		for (auto l : layers)
			delete l;
	}
	if (annoyingNeuronList)
		cudaFree(annoyingNeuronList);
	if (annoyingNeuronList2)
		cudaFree(annoyingNeuronList2);
}

int NeuralNetwork::addLayer(Layer* layer)
{
	assert(!annoyingNeuronList);
	size_t m = layer->outputSize;
	layers.push_back(layer);
	if (maxLayerSize < m)
		maxLayerSize = m;
	concreteBoundsS.push_back(std::make_shared<Vector<float>>(m, true));
	concreteBoundsD.push_back(std::make_shared<Vector<double>>(m, true));
	return layers.size() - 1;
}

template <typename T>
bool NeuralNetwork::run(const Vector<T>& input, const int label, bool sound)
{
	if (!annoyingNeuronList)
		gpuErrchk(cudaMalloc((void**)&annoyingNeuronList, maxLayerSize * sizeof(int)));
	if (!annoyingNeuronList2)
		gpuErrchk(cudaMalloc((void**)&annoyingNeuronList2, maxLayerSize * sizeof(int)));

	size_t outputSize = layers.back()->outputSize; // size of the output of the neural network
	getConcreteBounds<T>(0) = input;


	// creates an additional "layer" to check if the neuron "label" is the bigger one
	auto finalA = std::make_shared<Matrix<T>>(outputSize - 1, outputSize, false);
	setFinal<T> << <1, outputSize - 1 >> > (*finalA, annoyingNeuronList, label, outputSize, finalA->pitch());

	Vector<T> res;

	// computes the coefficients for each layer
	for (int p = 1; p < layers.size(); p++)
		layers[p]->eval(getConcreteBounds<T>(p), sound, false);
	finalA->mvm(res, getConcreteBounds<T>(layers.size() - 1));
	evaluateAffine<T>(res, ContainsZero<T>(), layers.size() - 1, false, sound, finalA);
	if (res.isPositive())
		return true;


	// computes the coefficients for each layer
	for (int p = 1; p < layers.size(); p++)
		layers[p]->eval(getConcreteBounds<T>(p), sound, true);
	
	finalA->mvm(res, getConcreteBounds<T>(layers.size() - 1));
	evaluateAffine<T>(res, ContainsZero<T>(), layers.size() - 1, false,sound, finalA);
	return res.isPositive();
}



template bool NeuralNetwork::run(const Vector<float>& input, const int label, bool sound);
template bool NeuralNetwork::run(const Vector<double>& input, const int label, bool sound);
