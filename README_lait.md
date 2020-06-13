# Lait

Lait is a learned abstract transformer for speeding up numerical program analysis. It selectively drops potentially redundant constraints from the outputs of the join transformer. The decision on which constraints to dropped is made by a graph neural netwrok. For more details on how Lait works, please read [Lait PLDI'20 paper](https://files.sri.inf.ethz.ch/website/papers/pldi20-lait.pdf). Lait now supports the Polyhedra and the Octagon domains.


## Setup
Apart from the dependencies required by ELINA, Lait requires `python3` (we use python3.6 but one can modify the version) and python libraries:
```
$ sudo apt install python3.6 python3.6-dev
$ curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py
$ python3.6 get-pip.py
$ pip install torch torchvision scikit-learn
```

After installing the python dependencies, one can compile and install Lait by:
```
$ cd elina_poly
$ make LAIT=1 && sudo make install
$ cd ../elina_oct
$ make LAIT=1 && sudo make install
```
With the flag `LAIT=1`, the Lait transformer functions `opt_pk_lait` and `opt_oct_lait` are added in the ELINA library.

## Example Usage
The examples usages of Lait can be found in `elina_poly/elina_test_poly_lait.c` and `elina_oct/elina_test_oct_lait.c`.