#include <iostream>

#include "gpupoly.h"

int main()
{
	NeuralNetwork* nn = create(2); // creates a network that has 2 input elements
	if (!nn) 
	{
		std::cout << "Error, it seems that the library is already used." << std::endl;
		return -1;
	}

	const double dataA[] = {
		1,1,
		1,-1
	}; // coefficients of a 2*2 matrix

	const double datab[] = { 0., 0.5 }; // coefficient of a bias vector

	// adds a linear layer
	int prev = addLinear_d(
		nn,
		0, // input of this layer is the input of the network, i.e. layer index 0
		2, // output of this layer has size 2
		dataA
	);

	// adds a bias
	prev = addBias_d(
		nn,
		prev, // input of this layer is the output of the linear layer
		datab
	);

	// adds a ReLU layer
	prev = addReLU(nn, prev, true);

	const double lowerBound1[] = { 1,2 };
	const double upperBound1[] = { 1,2 }; // First image: two pixels, no epsilon

	// Test it
	std::cout << "The first box ";
	if (test_d(nn,lowerBound1, upperBound1, 0))
		std::cout << "evaluates to 0!" << std::endl;
	else
		std::cout << "did not satisfy the property." << std::endl;

	const double lowerBound2[] = { -0.5,1.5 };
	const double upperBound2[] = { 0.5,2.5 }; // Second image: two pixels that both have a radius of epsilon=0.5

	// Test it
	std::cout << "The second box ";
	if (test_d(nn,lowerBound2, upperBound2, 0))
		std::cout << "evaluates to 0!" << std::endl;
	else
		std::cout << "did not satisfy the property." << std::endl;

	// Done with th network => clean it
	clean(nn);


	return 0;
}