#  GPUPoly library
#  This source file is part of ELINA (ETH LIbrary for Numerical Analysis).
#  ELINA is Copyright � 2020 Department of Computer Science, ETH Zurich
#  This software is distributed under GNU Lesser General Public License Version 3.0.
#  For more information, see the ELINA project website at:
#  http://elina.ethz.ch
#
#  THE SOFTWARE IS PROVIDED "AS-IS" WITHOUT ANY WARRANTY OF ANY KIND, EITHER
#  EXPRESS, IMPLIED OR STATUTORY, INCLUDING BUT NOT LIMITED TO ANY WARRANTY
#  THAT THE SOFTWARE WILL CONFORM TO SPECIFICATIONS OR BE ERROR-FREE AND ANY
#  IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE,
#  TITLE, OR NON-INFRINGEMENT.  IN NO EVENT SHALL ETH ZURICH BE LIABLE FOR ANY
#  DAMAGES, INCLUDING BUT NOT LIMITED TO DIRECT, INDIRECT,
#  SPECIAL OR CONSEQUENTIAL DAMAGES, ARISING OUT OF, RESULTING FROM, OR IN
#  ANY WAY CONNECTED WITH THIS SOFTWARE (WHETHER OR NOT BASED UPON WARRANTY,
#  CONTRACT, TORT OR OTHERWISE).

## @file CMakeLists.txt
## @author Fran&ccedil;ois Serre
## @brief Project definition for CMake
#
#  CMake project to be used to compile GPUPoly on platforms not handled by ELINA.
#

cmake_minimum_required (VERSION 3.17.1)
project ("GPUPoly" VERSION 0.10.0 LANGUAGES CXX CUDA)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CUDA_STANDARD 14)
set(CMAKE_CUDA_STANDARD_REQUIRED True)

find_package(CUDAToolkit)
include(GenerateExportHeader)

set(GPUPOLYCore "layers/maxpool2d.cu" "affineexpr.cu" "layers/conv2d.cu" "layers/dense.cu" "mmm.cu" "mvm.cu" "network.cu" "layers/relu.cu" "matrix.cu" "vector.cu" "filters.cu")

option(STRONG_FP_SOUNDNESS "If off, an element in a certified box is guaranteed to classify correctly if infinite precision arithmetic is used during the evaluation (weak certification). If on, an input in a certified box is guaranteed to classify correctly as long as IEEE754 arithmetic is used for evaluation, irrespective of the ordering of sums (strong certification)." ON)
set(CMAKE_CXX_VISIBILITY_PRESET hidden)
set(CMAKE_VISIBILITY_INLINES_HIDDEN yes)
configure_file("config.h.in" "config.h")
add_library (gpupoly SHARED ${GPUPOLYCore} "bindings.cu")
set_target_properties(gpupoly PROPERTIES VERSION 0.10.0 SOVERSION 0.10 CUDA_ARCHITECTURES 61)
generate_export_header(gpupoly)
target_include_directories (gpupoly PRIVATE .)
target_include_directories (gpupoly PRIVATE ${PROJECT_BINARY_DIR})
target_link_libraries(gpupoly CUDA::cublas) 
install(TARGETS gpupoly)


configure_file("../python_interface/gpupoly.py" "gpupoly.py")
configure_file("../python_interface/onnx2gpupoly.py" "onnx2gpupoly.py" COPYONLY)
configure_file("../python_interface/gpupoly_test.py" "gpupoly_test.py" COPYONLY)

add_executable (benchmark "benchmark.cu")
set_target_properties(benchmark PROPERTIES CUDA_ARCHITECTURES 61)
target_include_directories (benchmark PRIVATE ${PROJECT_BINARY_DIR})
target_include_directories(benchmark PRIVATE .)
target_link_libraries(benchmark gpupoly)