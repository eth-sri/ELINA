#
#
#  This source file is part of ELINA (ETH LIbrary for Numerical Analysis).
#  ELINA is Copyright 2019 Department of Computer Science, ETH Zurich
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


def test_volume_difference(filename, activation):
    assert activation in ["relu", "pool", "tanh", "sigm"]
    print("volume test for:", activation, filename)
    hrep = read_matrix(filename)

    time_fast = time()
    if activation == "relu":
        res_fast = fkrelu(hrep)
    elif activation == "pool":
        res_fast = fkpool(hrep)
    elif activation == "tanh":
        res_fast = fktanh(hrep)
    elif activation == "sigm":
        res_fast = fksigm(hrep)
    time_fast = time() - time_fast

    time_cdd = time()
    if activation == "relu":
        res_cdd = krelu_with_cdd(hrep)
    elif activation == "pool":
        res_cdd = kpool_with_cdd(hrep)
    elif activation == "tanh":
        res_cdd = ktanh_with_cdd(hrep)
    elif activation == "sigm":
        res_cdd = ksigm_with_cdd(hrep)
    time_cdd = time() - time_cdd

    res_concat = np.concatenate([res_cdd, res_fast], axis=0)

    volume_cdd = hrep2volume(res_cdd, "float")
    volume_fast = hrep2volume(res_fast, "float")
    volume_concat = hrep2volume(res_concat, "float")

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
test_volume_difference("k1/1.txt", "relu")
test_volume_difference("k1/2.txt", "relu")
test_volume_difference("k3/1.txt", "relu")
test_volume_difference("k3/5.txt", "relu")
test_volume_difference("k3/1.txt", "pool")
test_volume_difference("k3/2.txt", "pool")
for activation in ["tanh", "sigm"]:
    test_volume_difference("k1/1.txt", activation)
    test_volume_difference("k1/2.txt", activation)
    test_volume_difference("k2/1.txt", activation)
    test_volume_difference("k2/2.txt", activation)
