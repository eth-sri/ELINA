/*
 * GPUPoly library
 * A deep neural network verifier library running on GPU
 * Copyright (C) 2020 Francois Serre
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA
 */

/*!
  \file src/gpumem.h
  \brief Memory management on GPU
  \author Fran&ccedil;ois Serre
  Definition and implementation of a memory allocator on GPU
*/

#pragma once
#include "cuda_runtime.h"
#include "utils.h"
#include <stack>
#include<utility>
#include <iostream>
#include <assert.h>
template <bool verbose=false>

//! Memory allocator
/*!
  Memory allocator that tries to reuse the objects that are allocated on GPU to avoid the cost of actual allocations.
*/
class GPUMem {
	std::stack<std::pair<void*, size_t>> allocated;
	int inWild = 0;
	
public:
	void print()
	{
		std::cout << "wild:" << inWild << " allocated:" << allocated.size();
	}
	GPUMem() {

	}
	~GPUMem() {
		assert(inWild == 0);
		assert(allocated.empty());
	}
	//! Allocate a new memory space
	/*!
	  Allocate a new memory space, with a capacity equal or higher than size.
	  \param size Requested size (in bytes). After the allocation, size contains the actual allocated size.
	  \return A pointer to the allocated area.
	*/
	void* alloc(size_t& size)
	{
		inWild++;
		if (!allocated.empty())
		{
			auto e = allocated.top();
			allocated.pop();
			if (e.second >= size)
			{
				if (verbose)
					std::cout << "reuse alloc:" << (size / 1024. / 1024.) << "MB		Granted " <<( e.second/1024./1024.)<< "MB. In wild:"<<inWild << std::endl;
					
				size = e.second;

				return e.first;
			}
			cudaFree(e.first);
		}
		void* res;
		if (verbose)
			std::cout << "new alloc:" << (size / 1024. / 1024.) << "MB		Free pool contains " << allocated.size() << " elements."<<std::endl;
		gpuErrchk(cudaMalloc(&res,size));
		return res;
	}
	//! Free memory space
	/*!
	  After use, this function must be called to allow reuse or release of allocated memory.
	  \param adr Address of the allocated area to free.
	  \param size Size that was provided by alloc.
	*/
	void free(void* adr, const size_t size)
	{
		inWild--;
		if (inWild > allocated.size())
			allocated.push(std::make_pair(adr, size));
		else
		{
			cudaFree(adr);
			if (inWild < allocated.size())
			{
				auto e = allocated.top();
				allocated.pop();
				cudaFree(e.first);
			}
		}
	}
};