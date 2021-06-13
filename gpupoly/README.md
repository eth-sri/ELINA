# GPUPoly
GPUPoly is a deep neural network verifier library running on GPU. Given the structure and the weigths of a trained neural network, it checks whether an input box (given as a set of intervals) is guaranteed to fully classify as a given label.

A short introduction to GPUPoly can be seen in the following video:

[![Introduction video](https://img.youtube.com/vi/rM95uwgFFVw/hqdefault.jpg)](https://www.youtube.com/watch?v=rM95uwgFFVw)

## Strength of certification and floating-point soundness
Depending on its configuration, given a network and an input, GPUPoly can deliver different certification levels, listed here in increasing strength order:

0. No certification: GPUPoly was not able to prove that the whole input box classifies as the given label.
1. Certification without floating-point soundness. This yields no formal guarantee, but is in practice likely to mean 2.
2. Certified for an infinitly precise evaluation. If a box is certified with this level, any element of the box will classify to the given label, as long as the inference is computed using infinitly precise arithmetic.
3. Certified for an evaluation that uses double precision IEEE754 arithmetic, with any summation order and any rounding mode.
4. Certified for an evaluation that uses single precision IEEE754 arithmetic, with any summation order and any rounding mode.

By default, GPUPoly tries to deliver level 3 or 4 certifications. A parameter of GPUPoly allows to disable floating-point soundness to speed up the computation (which then becomes roughly twice or three times faster). This mode doesn't offer any formal guarantee anymore, but yields in practice the same results. Without floating point soundness, GPUPoly outputs certification of 1.

Conversely, it is possible to make GPUPoly search for weaker certification levels by switching off the `STRONG_FP_SOUNDNESS` compile option. If GPUPoly is compiled with this option, the letter `S` will appear after its version number (instead of `W`), and it will output certification levels of 0, 3 or 4. If the floating-point soundness parameter is disabled, GPUPoly has the same behaviour with or without the `STRONG_FP_SOUNDNESS` option. 

## Building the library
The library is built with ELINA, but can also be built using CMake on both Windows and linux architectures. 

### Windows
The library can be built using Visual Studio 2021 and the plugins `GitHub.VisualStudio` and `Visual Studio Tools for CMake`. With this configuration, it is possible to clone the repository from the start window of the IDE, and to compile and run it directly.

For older versions of Visual Studio, the CMake GUI can be used to create a Visual Studio project that can be used to compile and use the library.

The compilation creates a library `gpupoly.dll`

### Linux
The library can be built from the command line:

```
git clone https://github.com/eth-sri/ELINA
cd ELINA/gpupoly
cmake .
make
make install
```

## Using the library
### C and C++
The header file `gpupoly.h` contains the binding functions for the library. An example of use is provided in `example.cpp`.

### Python
After a succesful build, the build directory should contain the file `gpupoly.py` that contains python bindings.

### ONNX
The (python) script onnx2gpupoly.py allows to load an onnx network for use within GPUPoly.

# Citing GPUPoly
You can cite GPUPoly using the following bibtex:
```
@InProceedings{gpupoly,
  author    = {Fran{\c c}ois Serre and Christoph M{\"u}ller and Gagandeep Singh and Markus P{\"u}schel and Martin Vechev},
  title     = {Scaling Polyhedral Neural Network Verification on {GPU}s},
  booktitle = {Proc. Machine Learning and Systems (MLSys'21)},
  year      = {2021}
}
```

# Changelog
## [0.12.0] - 2021.04.15
Add the option to disable area heuristic in the ReLU, add support for asymetric padding in Conv2d. Revert to CMake.

## [0.11.0] - 2020.10.16
Renamed to GPUPoly for integration within ELINA.

## [0.10.0] - 2020.10.06
API is changed to support the use of different floating point precisions.

## [0.9.0] - 2020.09.13
Project is renamed SALSa.

## 2020.06.18
Initial commit.