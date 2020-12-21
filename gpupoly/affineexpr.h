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
   \file src/affineexpr.h
   \brief AffineExpr definition.
   \author Fran&ccedil;ois Serre

   Definition of the class AffineExpr, that represents a set of affine expressions. Implementation of the methods of this class can be found in src/layers/conv2d.cu.
  */

#pragma once
#include "matrix.h"
#include "vector.h"
#include <queue>


  //! Structure containing the parameters of a convolution.
struct ConvShape {
	int filters; //!< number of filters

	int kernel_size_rows; //!< kernel 
	int kernel_size_cols; //!< dimentions

	int input_rows; //!< input dimention (number of rows)
	int input_cols; //!< input dimention ((number of cols)
	int input_channels; //!< input dimention (number of channels)

	int stride_rows; //!< stride in rows
	int stride_cols; //!< stride in columns
	int padding_rows; //!< padding added on top and bottom (number of pixels)
	int padding_cols; //!< padding added on the left and right (number of pixels)

	int output_rows; //!< output dimention (number of rows)
	int output_cols; //!< output dimention (number of cols)

//! Constructor
/*!
  Creates a new convolution parameter structure.

  \param filters number of filters
  \param kernel_size_rows dimention of the kernel (rows)
  \param kernel_size_cols dimention of the kernel (cols)
  \param input_rows dimention of the input (number of rows)
  \param input_cols dimention of the input (number of cols)
  \param input_channels dimention of the input (number of channels)
  \param stride_rows stride shape (number of rows)
  \param stride_cols stride shape (number of cols)
  \param padding_rows padding (number of pixels to add at the top and bottom)
  \param padding_cols padding (number of pixels to add on the left and right)
 */
	ConvShape(
		const int filters,
		const int kernel_size_rows, const int kernel_size_cols,
		const int input_rows, const int input_cols, const int input_channels,
		const int stride_rows, const int stride_cols,
		const int padding_rows, const int padding_cols) :
		filters(filters),
		kernel_size_rows(kernel_size_rows), kernel_size_cols(kernel_size_cols),
		input_rows(input_rows), input_cols(input_cols), input_channels(input_channels),
		stride_rows(stride_rows), stride_cols(stride_cols),
		padding_rows(padding_rows), padding_cols(padding_cols),
		output_rows((input_rows + 2 * padding_rows - kernel_size_rows + stride_rows) / stride_rows), output_cols((input_cols + 2 * padding_cols - kernel_size_cols + stride_cols) / stride_cols) {}
	ConvShape() : kernel_size_rows(0) {}

	//! Indicates whether this is a valid convolution shape.
	__device__ __host__ operator bool() const {
		return kernel_size_rows > 0;
	}
	//! Size of the inputs
	__device__ __host__ int inputSize() const
	{
		return input_rows * input_cols * input_channels;
	}
	//! Size of the output
	__device__ __host__ int outputSize() const
	{
		return output_rows * output_cols * filters;
	}
	//! Indicates whether the matrix of a convolution with these parameters would be diagonal
	__device__ __host__ bool isDiagonal() const
	{
		return
			input_channels == filters &&
			kernel_size_cols == 1 &&
			kernel_size_rows == 1 &&
			stride_cols == 1 &&
			stride_rows == 1 &&
			padding_rows == 0 &&
			padding_cols == 0;
	}
	//! Creates the parameters of a diagonal matrix
	//! \param size Size of the diagonal matrix
	static ConvShape diagonal(int size)
	{
		return ConvShape(
			size,
			1, 1,
			1, 1, size,
			1, 1, 0, 0);
	}
	//! Computes the parameters of a product of convolution
	/*!
	  Computes the parameters that the composition of a convolution parametered by the current structure with a convolution parametered by rhs would have.
	  \param rhs Parameters of the right-hand side convolution
	*/
	ConvShape operator*(const ConvShape& rhs) const
	{
		if (!*this || !rhs)
			return ConvShape();
		assert(rhs.outputSize() == inputSize());
		if (rhs.isDiagonal())
			return *this;
		if (isDiagonal())
			return rhs;
		if (kernel_size_rows > 0 &&
			input_channels == rhs.filters &&
			input_cols == rhs.output_cols &&
			input_rows == rhs.output_rows)
		{
			const int nsr = (kernel_size_rows - 1) * rhs.stride_rows + rhs.kernel_size_rows;
			const int nsc = (kernel_size_cols - 1) * rhs.stride_cols + rhs.kernel_size_cols;
			if (nsr >= rhs.input_rows || nsc >= rhs.input_cols)
				return ConvShape();
			return ConvShape(
				filters,
				nsr,
				nsc,
				rhs.input_rows,
				rhs.input_cols,
				rhs.input_channels,
				stride_rows * rhs.stride_rows,
				stride_cols * rhs.stride_cols,
				padding_rows * rhs.stride_rows + rhs.padding_rows,
				padding_cols * rhs.stride_cols + rhs.padding_cols
			);
		}
		return ConvShape();
	}
	//! Print the current parameters.
	void print() const {
		std::cout << "Inputs: [" << input_rows << "," << input_cols << "," << input_channels << "] kernel: [" << kernel_size_rows << "," << kernel_size_cols << "] stride: [" << stride_rows << "," << stride_cols << "] filters: " << filters << " padding: [" << padding_rows << "," << padding_cols << "] Outputs: [" << output_rows << "," << output_cols << "," << filters << "]" << std::endl;
	}
};
//! Structure representing a set of affine expressions.
/*!

*/
template <typename T>
struct AffineExpr
{
	size_t m; //!< number of expressions (i.e. number of rows of the matrix)
	size_t n; //!< number of variables in each expression (i.e. number of columns in the matrix)
	int layer; //!< index of the layer the affine expressions refers to
	bool up; //!< true if these expressions are upper bounds. 
	int* rows; //!< integer array in global memory containing the original position of each row of the matrix: the ith element indicates at what position the ith row was originally
	std::shared_ptr<const Matrix<T>> A; //!< pointer to a Matrix containing the coefficients of the expressions (the ith row at the jth column is the jth coefficient of the ith affine expression). If null, this means that the matrix is the identity.
	std::shared_ptr<const Vector<T>> b; //!< pointer to a Vector containing the constants of each expression (the ith one is the constant of the ith expression). If null, it means that the coefficiesnts are 0.
	ConvShape cs; //!< If evaluates to false, the expression is dense. Otherwise, its shape corresponds to a convolution having these parameters.
	bool sound; //!< If floating point sound artithmetic should be used even if slower.
	//! Constructor
/*!
  Creates a new convolution parameter structure.

  \param m number of expressions (i.e. number of rows of the matrix)
  \param n number of variables in each expression (i.e. number of columns in the matrix)
  \param layer index of the layer the affine expressions refers to.
  \param up true if these expressions are upper bounds.
  \param rows integer array in global memory containing the original position of each row of the matrix: the ith element indicates at what position the ith row was originally
  \param A pointer to a Matrix containing the coefficients of the expressions (the ith row at the jth column is the jth coefficient of the ith affine expression). If null, this means that the matrix is the identity.
  \param b pointer to a Vector containing the constants of each expression (the ith one is the constant of the ith expression). If null, it means that the coefficiesnts are 0.
  \param cs If evaluates to false, the expression is dense. Otherwise, its shape corresponds to a convolution having these parameters.
 */
	AffineExpr(
		const size_t m, const size_t n,
		const int layer, bool up, int* rows,
		const std::shared_ptr<const Matrix<T>>& A/* = nullptr*/,
		const std::shared_ptr<const Vector<T>>& b/* = nullptr*/,
		const ConvShape& cs /*= ConvShape()*/,
		bool sound
		) :
		m(m), n(n), layer(layer), up(up), rows(rows),
		A(A),
		b(b), cs(cs),sound(sound)
	{
		assert(!A || (A->m() == m && A->n() == n));
		assert(!b || (b->size() == m && !b->interval()));
	}
	//! Returns a pointer to the matrix. If A is null, the matrix is created.
	std::shared_ptr<const Matrix<T>> getA() const;

	//! Returns Evaluates the expressions given the value of the neurons of layer.
	template <typename Tr>
	Vector<T> evaluate(const Vector<Tr>& bounds) const;

	//! Returns Evaluates the expressions given the value of the neurons of layer, and updates dest if the new bounds are better.
	void evaluateAndUpdate(Vector<T>& dest, const Vector<T>& bounds) const;

	//! Removes the row of the expression that do not appear in rows.
	void selectRows(size_t size, int* rows);

	//! Functor that orders affine expressions by the index of the layer they refer to.
	struct SmallerLayer {
		bool operator() (const AffineExpr<T>& lhs, const AffineExpr<T>& rhs) {
			return lhs.layer < rhs.layer;
		}
	};

	//! Type of a priority queue of AffineExpr that classifies them according to their layer.
	typedef std::priority_queue < AffineExpr<T>, std::vector<AffineExpr<T>>, AffineExpr<T>::SmallerLayer > Queue;
};