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
   \file benchmark.cu
   \brief Benchmark for mlsys for non-ONNX networks. 
   \author Fran&ccedil;ois Serre

   This program tests the non-ONNX networks submitted in mlsys, avoiding the performance penalty of ERAN.
  */

#include "gpupoly.h"
#include <iostream>
#include <vector>
#include <stack>
#include <fstream>
#include <algorithm>
#include <tuple>
#include <chrono>
#include <assert.h>
using namespace std;


// Parse a vector formated as "[a, b, c]" from an input stream
template <typename T>
vector <T> readVector(istream& is)
{
	vector<T> l;
	char c;
	is >> c;
	assert(c == '[');
	while ((c = is.peek()) != ']')
	{
		if (c == ' ' || c == ',' || c == '\n' || c == '\r')
		{
			is.get();
			continue;
		}
		T tmp;
		is >> tmp;
		l.push_back(tmp);
	}
	is.get();
	return l;
}
vector <double> readVector(istream& is, size_t n)
{
	vector<double> l(n);
	for (size_t i = 0; i < n; i++)
	{
		char c = is.peek();
		while (c == ' ' || c == ',' || c == '\n' || c == '\r' || c == '[' || c == ']')
		{
			is.get();
			c = is.peek();
		}
		is >> l[i];
	}
	char c = is.peek();
	while (c == ' ' || c == ',' || c == '\n' || c == '\r' || c == ']')
	{
		is.get();
		c = is.peek();
	}
	return l;
}
// Parse a matrix formated as "[[a, b], [c,d]]" from an input stream
static vector<vector<double>> readMatrix(istream& is)
{
	vector<vector<double>> l;
	char c;
	is >> c;
	assert(c == '[');
	while ((c = is.peek()) != ']')
	{
		if (c == ' ' || c == ',' || c == '\n' || c == '\r')
		{
			is.get();
			continue;
		}
		l.push_back(readVector<double>(is));

	}
	is.get();
	return l;
}

vector<tuple<vector<double>, vector<double>, int>> readCifar10(const double epsilon, const vector<double>& mean, const vector<double>& std, const bool shuffle) {

#ifdef _MSC_VER
	ifstream f("C:\\Users\\Francois\\source\\cifar10_test.csv");
#else
	ifstream f("/local/home/gagandeep/eran/data/cifar10_test_full.csv");
#endif
	vector<tuple<vector<double>, vector<double>, int>> res;
	char c;
	while ((c = f.peek()) != EOF)
	{
		if (c == ' ' || c == '\n' || c == '\r')
		{
			f.get();
			continue;
		}
		int label;
		f >> label;
		f.get();
		vector<double> l(3072);
		vector<double> u(3072);
		for (int i = 0; i < 3072; i++)
		{
			int e;
			f >> e;
			f.get();
			double down = std::min(std::max(e / 255.0 - epsilon, 0.), 1.);
			double up = std::min(std::max(e / 255.0 + epsilon, 0.), 1.);
			int adr = shuffle ? (i % 3) * 1024 + i / 3 : i;
			l[adr] = (down - mean[i % 3]) / std[i % 3];
			u[adr] = (up - mean[i % 3]) / std[i % 3];
		}
		//res.push_back(pair<int, vector<Intv>>(label, tmp));
		res.emplace_back(l, u, label);
	}

	f.close();
	//cout << "Read " << res.size() << " inputs." << endl;
	return res;
}

vector<tuple<vector<double>, vector<double>, int>> readMNIST(const double epsilon, const vector<double>& mean, const vector<double>& std) {
#ifdef _MSC_VER
	ifstream f("C:\\Users\\Francois\\source\\mnist_test_full.csv");
#else
	ifstream f("/local/home/gagandeep/eran/data/mnist_test_full.csv");
#endif
	vector<tuple<vector<double>, vector<double>, int>> res;
	char c;
	while ((c = f.peek()) != EOF)
	{
		if (c == ' ' || c == '\n' || c == '\r')
		{
			f.get();
			continue;
		}
		int label;
		f >> label;
		f.get();
		vector<double> l(784);
		vector<double> u(784);
		for (int i = 0; i < 784; i++)
		{
			int e;
			f >> e;
			f.get();
			double down = std::min(std::max(e / 255.0 - epsilon, 0.), 1.);
			double up = std::min(std::max(e / 255.0 + epsilon, 0.), 1.);
			l[i] = (down - mean[0]) / std[0];
			u[i] = (up - mean[0]) / std[0];
		}
		//res.push_back(pair<int, vector<Intv>>(label, tmp));
		res.emplace_back(l, u, label);
	}

	f.close();
	//cout << "Read " << res.size() << " inputs." << endl;
	return res;
}



NeuralNetwork* readNetwork(const string filename, const size_t inputSize, const bool channels_first, std::vector<double>& mean, std::vector<double>& std)
{
	ifstream f(filename);
	char c;
	//std::vector<Layer*> layers;
	auto nn = create(inputSize);
	
	int branches;
	//std::list<int> prev = { -1 };

	//layers.emplace_back(new Input(inputSize));
	int prev = 0;
	bool meanSet = false;
	while ((c = f.peek()) != EOF)
	{
		if (c == ' ' || c == '\n' || c == '\r')
		{
			f.get();
			continue;
		}
		string tmp;
		f >> tmp;
		//		Layer::ActivationType type;
		if (tmp == "Normalize")
		{
			f.ignore(6);
			mean = readVector<double>(f);
			f.ignore(5);
			std = readVector<double>(f);
			meanSet = true;
			continue;
		}
		else if (tmp == "Conv2D")
		{
			f >> tmp;

			f.ignore(9); // " filters="
			int filters;
			f >> filters;
			f.ignore(14); // ", kernel_size="
			auto kernel_size = readVector<int>(f);
			f.ignore(14); // ", input_shape="
			auto input_shape = readVector<int>(f);
			std::vector<int> stride = { 1,1 };
			int padding = 0;
			if (f.peek() == ',')
			{
				f.ignore(9); // ", stride="
				stride = readVector<int>(f);
				f.ignore(10); // ", padding="
				f >> padding;
			}
			auto A = readVector(f, filters * kernel_size[0] * kernel_size[1] * input_shape[2]);
			auto b = readVector<double>(f);
			//std::cout << "Conv2D " << tmp << endl;
			//std::cout << "filters:" << filters << " kernel:[" << kernel_size[0] << "," << kernel_size[1] << "] input:[" << input_shape[0] << "," << input_shape[1] << "," << input_shape[2] << "] stride:[" << stride[0] << "," << stride[1] << "] padding:" << padding << endl;
			//std::cout << b.size() << endl;

			/*Layer* conv = new Conv2D<true>(filters, kernel_size[0], kernel_size[1], input_shape[0], input_shape[1], input_shape[2], stride[0], stride[1], padding, padding, Vector(A), prev);
			prev = layers.addLayer(conv);*/
			int pad[2] = { padding,padding };
			prev = addConv2D_d(nn, prev, filters, kernel_size.data(), input_shape.data(), stride.data(), pad, A.data());

			size_t m = ((input_shape[0] + 2 * pad[0] - kernel_size[0] + stride[0]) / stride[0]) * ((input_shape[1] + 2 * pad[1] - kernel_size[1] + stride[1]) / stride[1]) * filters;
			std::vector<double> bias(m);
			size_t k = 0;
			if (channels_first)
				for (int i = 0; i < filters; i++)
					for (int j = 0; j < m / filters; j++)
						bias[k++] = b[i];
			else
				for (int j = 0; j < m / filters; j++)
					for (int i = 0; i < filters; i++)
						bias[k++] = b[i];

			prev = addBias_d(nn, prev, bias.data());

			if (tmp == "ReLU,")
				//prev = layers.addLayer(new ReLU(b2.size(), nullptr, prev));
				prev = addReLU(nn, prev);
			else if (tmp == "Affine,")
				;
			else
				throw - 1;
		}
		else if (tmp == "ReLU")
		{
			auto A = readMatrix(f);
			int m = A.size();
			int n = A.front().size();
			double* dataA = new double[m * n];
			size_t k = 0;
			for (auto row : A)
				for (auto e : row)
					dataA[k++] = e;
			auto b = readVector<double>(f);
			prev = addLinear_d(nn, prev, m, dataA);
			delete[] dataA;
			prev = addBias_d(nn, prev, b.data());
			prev = addReLU(nn, prev);
			//std::cout << "ReLu" << endl;
		}
		/*else if (tmp == "Tanh")
		{
			type = Layer::Tanh;
			auto A = readMatrix(f);
			auto b = readVector<double>(f);
			int id = layers.size();
			layers.emplace_back(new OldDense(type, Matrix(A), std::make_shared<const Vector>(b),id, prev));
			//std::cout << "Tanh" << endl;
			//std::cout << "In:" << A[0].size() << " Out:" << A.size() << endl;
			prev = id;
		}*/
		else if (tmp == "Affine")
		{
			auto A = readMatrix(f);
			int m = A.size();
			int n = A.front().size();
			double* dataA = new double[m * n];
			size_t k = 0;
			for (auto row : A)
				for (auto e : row)
					dataA[k++] = e;
			auto b = readVector<double>(f);
			prev = addLinear_d(nn, prev, m, dataA);
			delete[] dataA;
			prev = addBias_d(nn, prev, b.data());
			//std::cout << "Affine" << endl;
		}
		else if (tmp == "MaxPooling2D")
		{
			f.ignore(11);// "pool_size="
			auto pool_size = readVector<int>(f);
			f.ignore(14);// ", input_shape="
			auto input_shape = readVector<int>(f);
			int padding[2] = { 0,0 };
			input_shape.insert(input_shape.begin(), 1);
			prev = addMaxPool2D(nn, prev,pool_size.data(),input_shape.data(),pool_size.data(),padding);
		}
		else if (tmp == "SkipNet1")
		{
			branches = prev;
			//cout << "SN1" << endl;
		}
		else if (tmp == "SkipNet2")
		{
			std::swap(prev, branches);
			//cout << "SN2" << endl;
		}
		else if (tmp == "SkipCat")
		{
			//cout << "SC" << endl;
			string t2;
			f >> t2;
			//cout << t2 << endl;
			//prev = layers.addLayer(new Concat(layers[branches]->outputSize, branches, layers[prev]->outputSize, prev));
			prev = addConcat(nn, branches, prev);
		}
		else if (tmp == "ParSum1")
		{
			branches = prev;
			//cout << "PS1" << endl;
		}
		else if (tmp == "ParSum2")
		{
			std::swap(prev, branches);
			//cout << "PS2" << endl;
		}
		else if (tmp == "ParSumReLU")
		{
			//cout << "PSRELU" << endl;
			//prev = layers.addLayer(new ParSum(layers[branches]->outputSize, branches, prev));
			//prev = layers.addLayer(new ReLU(layers[branches]->outputSize, nullptr, prev));
			prev = addParSum(nn, branches, prev);
			prev = addReLU(nn, prev);
		}
		else
		{
			cout << "Unknown layer type: " << tmp << endl;
			throw 1;
		}
	}
	f.close();
	if (!meanSet)
	{
		if (inputSize == 784)
		{
			mean.push_back(0);
			std.push_back(1);
		}
		else
		{
			mean = std::vector<double>({ 0.5,0.5,0.5 });
			std = std::vector<double>({ 1,1,1 });
		}
	}
	cout << "Read " << prev << " layers." << endl;
	//return layers;
	return nn;
}

void testNetwork(const string& network, const size_t inputSize, const bool channels_first, const double epsilon, const int nbTests=250,bool verbose=false,string log="")
{
	std::vector<double> mean;
	std::vector<double> std;

	cout << network << " e=" << epsilon << endl;
	auto nn=readNetwork(network, inputSize,channels_first, mean, std);
	//NeuralNetwork test(nn);
	auto inputs = inputSize == 784 ? readMNIST(0, mean, std) : readCifar10(0, mean, std, true); // loads the inputs


	vector<int> candidates;
	for (int i = 0; i < inputs.size() && i < nbTests; i++) // ... we count the number of inputs for which the network passes the certification
		if (test_d(nn, get<0>(inputs[i]).data(), get<1>(inputs[i]).data(), get<2>(inputs[i]),true))
			candidates.push_back(i);
	cout << "Candidates:" << candidates.size() << endl;

	vector<int> results(candidates.size());
	vector<double> times(candidates.size());


	inputs = inputSize == 784 ? readMNIST(epsilon, mean, std) : readCifar10(epsilon, mean, std, true); // loads the inputs
	
	
	//warm up
	test_d(nn, get<0>(inputs[0]).data(), get<1>(inputs[0]).data(), get<2>(inputs[0]), true);
	
	auto startTimer = std::chrono::high_resolution_clock::now();
	auto prev = startTimer;
	int good = 0;
	int total = 0;
	for (int j = 0; j < candidates.size() && j < nbTests; j++) // ... we count the number of inputs for which the network passes the certification
	{
		int i = candidates[j];
		bool res = test_d(nn, get<0>(inputs[i]).data(), get<1>(inputs[i]).data(), get<2>(inputs[i]));
		//res = test(get<0>(inputs[i]).data(), get<1>(inputs[i]).data(), get<2>(inputs[i]),true);
		good += res;
		total++;
		auto stopTimer = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> diff = stopTimer - prev;
		times[j] = diff.count();
		results[j] = res;
		prev = stopTimer;

		if (verbose)
		{

			std::chrono::duration<double> diff = stopTimer - startTimer;
			cout << i << ":	" << res << "	so far " << good * 100.f / (total) << "% in "<<(diff.count()/total) <<"s per image. \r" << endl << flush;
		}
	}
	auto stopTimer = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> diff = stopTimer - startTimer;
	cout << "Verified:" << good << endl;
	cout << (diff.count()*1000 / total) << "ms per image" << endl << endl;
	clean(nn);

	if (log.length())
	{
		ofstream f(log);
		for (int j = 0; j < candidates.size() && j < nbTests; j++)
			f << candidates[j] << "," << results[j] << "," << (int)(times[j] * 1000000) << endl;
	}
}

int main()
{
	cudaSetDevice(1);
#ifdef _MSC_VER
	testNetwork("C:\\Users\\Francois\\source\\FFNN_cifar10_new.pyt", 3072, true, 1.0 / 500, 500, false,"testtime.txt");
	//testNetwork("C:\\Users\\Francois\\source\\ffnnRELU__Point_6_500.pyt", 3072, true, 1.0 / 500, 500, false);
	
#else
		testNetwork("/local/home/christoph/ERAN/nets/mnist/FFNN_mnist_new.pyt", 784, true,8.0 /255,10000, false, "FFNN_mnist_new.txt");
		testNetwork("/local/home/christoph/ERAN/nets/mnist/ConvBig_mnist_new.pyt", 784, true,3.0/10,  10000, false,"ConvBig_mnist_new.txt");
		testNetwork("/local/home/christoph/ERAN/nets/mnist/convSuper_mnist_new.pyt", 784, true, 8.0/255,  10000, false,"convSuper_mnist_new.txt");

		testNetwork("/local/home/christoph/ERAN/nets/cifar10/FFNN_cifar10_new.pyt", 3072, true, 1.0 / 500, 10000, false, "FFNN_cifar10_new.txt");
		testNetwork("/local/home/christoph/ERAN/nets/cifar10/ConvBig_cifar10_new.pyt", 3072, true, 8.0 / 255, 10000, false, "ConvBig_cifar10_new.txt");
		testNetwork("/local/home/christoph/ERAN/nets/cifar10/ConvLargeIBP.pyt", 3072, true, 8.0 / 255, 10000, false, "ConvLargeIBP.txt");
		
		testNetwork("/local/home/christoph/ERAN/nets/cifar10/ResNetTiny_PGD.pyt", 3072, true, 1.0 / 500, 500, false, "ResNetTiny_PGD.txt");
		testNetwork("/local/home/christoph/ERAN/nets/cifar10/ResNet18_PGD_new.pyt", 3072, true, 1.0 / 500, 500, false, "ResNet18_PGD_new.txt");
				
		testNetwork("/local/home/christoph/ERAN/nets/cifar10/ResNetTiny_DiffAI.pyt", 3072, true, 8.0 / 255, 500, false, "ResNetTiny_DiffAI.txt");
		testNetwork("/local/home/christoph/ERAN/nets/cifar10/SkipNet18_DiffAI.pyt", 3072, true, 8.0 / 255, 500, false, "SkipNet18_DiffAI.txt");
		testNetwork("/local/home/christoph/ERAN/nets/cifar10/ResNet18_DiffAI.pyt", 3072, true, 8.0 / 255, 500, false, "ResNet18_DiffAI.txt");
		testNetwork("/local/home/christoph/ERAN/nets/cifar10/ResNet34_DiffAI.pyt", 3072, true, 8.0 / 255, 500, false, "ResNet34_DiffAI.txt");

		testNetwork("/local/home/christoph/ERAN/nets/cifar10/ResNet18_PGD.pyt", 3072, true, 1.0 / 500, 500, false, "ResNet18_PGD.txt");
#endif
		return 0;
}