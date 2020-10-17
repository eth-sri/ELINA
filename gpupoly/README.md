# GPUPoly
GPUPoly is a deep neural network verifier library running on GPU. Given the structure and the weigths of a trained neural network, it checks whether an input box (given as a set of intervals) is guaranteed to fully classify as a given label.

## Strength of certification and floating-point soundness
Depending on its configuration, given a network and an input, SALSa can deliver different certification levels, listed here in increasing strength order:

0. No certification: GPUPoly was not able to prove that the whole input box classifies as the given label.
1. Certification without floating-point soundness. This yields no formal guarantee, but is in practice likely to mean 2.
2. Certified for an infinitly precise evaluation. If a box is certified with this level, any element of the box will classify to the given label, as long as the inference is computed using infinitly precise arithmetic.
3. Certified for an evaluation that uses double precision IEEE754 arithmetic, with any summation order and any rounding mode.
4. Certified for an evaluation that uses single precision IEEE754 arithmetic, with any summation order and any rounding mode.

By default, GPUPoly tries to deliver level 2 certifications, and will therefore output a result that will be either 0 or 2. A parameter of GPUPoly allows to disable floating-point soundness to speed up the computation (which then becomes roughly twice or three times faster). This mode doesn't offer any formal guarantee anymore, but yields in practice the same results. Without floating point soundness, GPUPoly outputs certification levels that can be 0 or 1.

Conversely, it is possible to make GPUPoly search for stronger certification levels using the `STRONG_FP_SOUNDNESS` compile option. If GPUPoly is compiled with this option, the letter `S` will appear after its version number (instead of `W`), and it will output certification levels of 0, 3 or 4. If the floating-point soundness parameter is disabled, SALSa has the same behaviour with or without the `STRONG_FP_SOUNDNESS` option. 

## Building the library
The library can be built using CMake on both Windows and linux architectures. In both cases, the project `gpupolyExample` allows to check that the library is built and works correctly.

### Windows
The library can be built using Visual Studio 2019 and the plugins `GitHub.VisualStudio` and `Visual Studio Tools for CMake`. With this configuration, it is possible to clone the repository from the start window of the IDE, and to compile and run it directly.

For older versions of Visual Studio, the CMake GUI can be used to create a Visual Studio project that can be used to compile and use the library.

The compilation creates a library `gpupoly.dll` and an executable `gpupolyExample.exe` that can be run to test the library.

### Linux
The library can be built from the command line:

```
git clone https://github.com/fserre/gpupoly.git
cd gpupoly
cmake .
make
make install
./gpupolyExample
```

## Using the library
### C and C++
The header file `include/gpupoly.h` contains the binding functions for the library. An example of use is provided in `test/example.cpp`.

### Python
After a succesful build, the build directory should contain the file `gpupoly.py` that contains python bindings. The file `python/example_tf.py` contains an example that trains a network using TensorFlow 2, and that uses SALSa to test its robustness.

### ONNX
The (python) script onnx2gpupoly.py allows to load an onnx network for use within GPUPoly.

# Changelog
## [0.9.0] - 2020.10.16
Renamed to GPUPoly for integration within ELINA.

## [0.10.0] - 2020.10.06
API is changed to support the use of different floating point precisions.

## [0.9.0] - 2020.09.13
Project is renamed SALSa.

## 2020.06.18
Initial commit.