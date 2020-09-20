#
#
#  This source file is part of ELINA (ETH LIbrary for Numerical Analysis).
#  ELINA is Copyright © 2019 Department of Computer Science, ETH Zurich
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


from fppoly_imports import *

import numpy as np
from ctypes import *


class MatDouble_c(Structure):
    _fields_ = [
        ('rows', c_int),
        ('cols', c_int),
        ('data', POINTER(c_double))
    ]


new_MatDouble_c = fconv_api.new_MatDouble
new_MatDouble_c.argtype = [c_int, c_int, POINTER(c_double)]
new_MatDouble_c.restype = MatDouble_c

free_MatDouble_c = fconv_api.free_MatDouble
free_MatDouble_c.argtype = MatDouble_c
free_MatDouble_c.restype = None

fkrelu_c = fconv_api.fkrelu
fkrelu_c.argtype = [MatDouble_c]
fkrelu_c.restype = MatDouble_c


def fkrelu(inp_mat: np.ndarray) -> np.ndarray:
    """
    Input in format b + Ax >= 0. The input has to be octahedron in a certain format.
    An example of possible inp is:
        [0.4, 1, 0],
        [0.5, -1, 0],
        [0.25, 0, 1],
        [0.75, 0, -1]
    Which describes a system:
        x1  <= 0.4
        -x1 <= 0.5
        x2  <= 0.25
        -x2 <= 0.75
    """
    rows, cols = inp_mat.shape
    k = cols - 1
    assert k >= 1
    inp_mat = inp_mat.flatten().tolist()
    data_c = (c_double * (rows * cols))(*inp_mat)

    inp_hrep = new_MatDouble_c(rows, cols, data_c)
    out_hrep = fkrelu_c(inp_hrep)
    assert out_hrep.cols == 2 * k + 1

    out = [0] * (out_hrep.rows * out_hrep.cols)
    for i in range(out_hrep.rows * out_hrep.cols):
        out[i] = out_hrep.data[i]
    out = np.array(out)
    out = out.reshape(out_hrep.rows, out_hrep.cols)

    free_MatDouble_c(inp_hrep)
    free_MatDouble_c(out_hrep)

    return out
