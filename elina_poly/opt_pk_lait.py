#
#
#  This source file is part of ELINA (ETH LIbrary for Numerical Analysis).
#  ELINA is Copyright © 2021 Department of Computer Science, ETH Zurich
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


import torch
from opt_pk_lait_model import PolicyGCN

model = None


def init(model_path):
    model_path = model_path.decode('ascii')
    global model
    model = PolicyGCN.load(model_path)


def lait(block_lens, features, edges):
    results = [0] * sum(block_lens)

    feature_blocks, index_blocks, reverse_index_blocks, edge_blocks = [], [], dict(), []
    i = 0
    for block_len in block_lens:
        if block_len < 20:
            i += block_len
            continue

        for _ in range(block_len):
            if features[i][4] > 2:
                feature_blocks.append(features[i])
                index_blocks.append(i)
                reverse_index_blocks[i] = len(index_blocks) - 1
            i += 1
    for edge in edges:
        if edge[0] in reverse_index_blocks and edge[1] in reverse_index_blocks:
            edge_blocks.append((reverse_index_blocks[edge[0]], reverse_index_blocks[edge[1]], edge[2]))

    if len(feature_blocks) == 0:
        return results

    Y = model.predict(feature_blocks, edge_blocks)
    for i, y in enumerate(Y):
        results[index_blocks[i]] = int(y)

    return results
