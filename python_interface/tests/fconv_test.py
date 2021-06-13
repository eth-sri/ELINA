#
#
#  This source file is part of ELINA (ETH LIbrary for Numerical Analysis).
#  ELINA is Copyright 2021 Department of Computer Science, ETH Zurich
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
#
#

import sys
sys.path.insert(0, '../')
import cdd
import scipy.spatial
import numpy as np
from fconv import *
from time import time


activation2fast = {
    "relu": fkrelu,
    "pool": fkpool,
    "tanh": fktanh,
    "sigm": fksigm
}

activation2cdd = {
    "relu": krelu_with_cdd,
    "pool": kpool_with_cdd,
    "tanh": ktanh_with_cdd,
    "sigm": ksigm_with_cdd
}


def read_matrix(filename) -> np.ndarray:
    f = open("../../fconv/octahedron_hrep/" + filename, "r")
    arr = [list(map(float, l.strip().split())) for l in f.readlines()]
    return np.array(arr)


def hrep2vrep(hrep: np.ndarray, number_type: str, check_bounded=True) -> np.ndarray:
    """Note that first element of the vertex will be 1."""
    hrep = cdd.Matrix(hrep, number_type=number_type)
    hrep.rep_type = cdd.RepType.INEQUALITY
    vrep = cdd.Polyhedron(hrep).get_generators()
    vrep = np.array(vrep)
    if check_bounded and vrep.size > 0:
        # Making sure that the resulting polyhedron is bounded
        assert (vrep[:, 0] == 1).all()
    if vrep.size == 0:
        return None
    return vrep


def hrep2volume(hrep: np.ndarray, number_type: str) -> float:
    vrep = hrep2vrep(hrep, number_type)
    assert (vrep[:, 0] == 1).all()
    vertices = vrep[:, 1:]
    return scipy.spatial.ConvexHull(vertices).volume


def compute_k1_relaxation(hrep, activation):
    # This function is not applicable to Pool relaxation
    assert activation in ["relu", "tanh", "sigm"]
    K = hrep.shape[1] - 1
    final = []
    for xi in range(K):
        inp_1d = []
        for h in hrep:
            skip = False
            for j in range(K):
                if j == xi:
                    continue
                if h[j + 1] != 0:
                    skip = True
            if not skip:
                inp_1d.append([h[0], h[xi + 1]])
        inp_1d = np.array(inp_1d)
        assert inp_1d.shape == (2, 2)
        out_1d = activation2cdd[activation](inp_1d)
        assert out_1d.shape[1] == 3
        out = []
        for h_1d in out_1d:
            h = [0] * (2 * K + 1)
            h[0] = h_1d[0]
            h[xi + 1] = h_1d[1]
            h[xi + 1 + K] = h_1d[2]
            out.append(h)
        out = np.array(out)
        final.append(out)
    final = np.vstack(final)
    return final


def test_volume_k1(filename, activation, precision="float"):
    assert activation in ["relu", "tanh", "sigm"]
    print("k1 volume test", activation, filename)
    hrep = read_matrix(filename)

    res_fast = activation2cdd[activation](hrep)
    res_1d = compute_k1_relaxation(hrep, activation)

    volume_fast = hrep2volume(res_fast, precision)
    volume_1d = hrep2volume(res_1d, precision)

    over_approximation = round(volume_1d / volume_fast, 4)

    print("\tComparison with k1 volume (higher is better)", over_approximation)


def test_volume_optimal(filename, activation, precision="float"):
    assert activation in ["relu", "pool", "tanh", "sigm"]
    print("optimal volume test", activation, filename)
    hrep = read_matrix(filename)

    time_fast = time()
    res_fast = activation2fast[activation](hrep)
    time_fast = time() - time_fast

    time_cdd = time()
    res_cdd = activation2cdd[activation](hrep)
    time_cdd = time() - time_cdd

    res_concat = np.concatenate([res_cdd, res_fast], axis=0)

    volume_cdd = hrep2volume(res_cdd, precision)
    volume_fast = hrep2volume(res_fast, precision)
    volume_concat = hrep2volume(res_concat, precision)

    over_approximation = round(volume_fast / volume_cdd, 4)
    soundness = round(volume_concat / volume_cdd, 4)
    speedup = round(time_cdd / time_fast)
    # The sound over-approximation is >= 1 and correctness check should = 1.
    # However, due to numerical errors in computing volume / converting to V
    # slight deviations are possible - it's okay.
    print("\tOver approximation", over_approximation)
    print("\tSoundness check (should be close to 1)", soundness)
    print("\tSpeedup", speedup)


# Unfortunately library for computing volume often fails due to numerical errors.
# Thus it is possible to perform this test only for some of the inputs.
test_volume_optimal("k1/1.txt", "relu")
test_volume_optimal("k1/2.txt", "relu")
test_volume_optimal("k3/1.txt", "relu")
test_volume_optimal("k3/5.txt", "relu")
test_volume_optimal("k3/1.txt", "pool")
test_volume_optimal("k3/2.txt", "pool")

for activation in ["tanh", "sigm"]:
    for filename in ["k1/1.txt", "k1/2.txt", "k2/1.txt", "k2/2.txt"]:
        test_volume_optimal(filename, activation, "fraction")

for activation in ["relu", "tanh", "sigm"]:
    for filename in ["k2/1.txt", "k2/2.txt"]:
        test_volume_k1(filename, activation, "fraction")
