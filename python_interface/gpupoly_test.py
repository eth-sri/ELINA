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


## @file gpupoly_test.py
## @author Fran&ccedil;ois Serre
## @brief Example of use of the ONNX frontend for GPUPoly.
#
#  This example loads an ONNX network, and uses the GPUPoly library to check its robustness.
#


import tensorflow as tf
import onnxruntime as rt
import numpy as np
import onnx
import time
import argparse
from onnx2gpupoly import onnx2gpupoly


parser=argparse.ArgumentParser("GPUPoly ONNX frontend")
parser.add_argument("network")
parser.add_argument("eps", type=float)
parser.add_argument("-limit", type=int, default=0)
args=parser.parse_args()


print("Reading "+args.network+"...")
nn=onnx2gpupoly(onnx.load(args.network).graph)

if(nn.input_size==784):
    print("Using MNIST dataset.")
    _, (x_test, y_test) = tf.keras.datasets.mnist.load_data()  # download dataset
    x_test=x_test/255.0
    std=1
    mean=0
elif (nn.input_size==3072):
    print("Using CIFAR10 dataset.")
    _, (x_test, y_test) = tf.keras.datasets.cifar10.load_data()  # download dataset
    x_test = np.transpose(x_test / 255.0, [0, 3, 1, 2])
    std = [0.2023, 0.1994, 0.2010]
    mean = [0.4914, 0.4822, 0.4465]
    std=np.array(std).reshape([1,3,1,1])
    mean = np.array(mean).reshape([1, 3, 1, 1])
else:
    print("Unknown dataset size:"+str(nn.input_size))
    exit(1)


print("Network and dataset loaded; running inference...")
candidates=[]
i=0
if(args.limit==0):
    nbImages=len(x_test)
else:
    nbImages = args.limit

x = ((np.clip(x_test, 0, 1)-mean)/std).astype("float32")
while i<nbImages:
   if (nn.test(x[i], x[i], y_test[i], False)>0):
       candidates.append(i)
   i+=1
print(f"{len(candidates)} candidates among {nbImages} images. ({len(candidates)/nbImages*100}%)")


epsilon=args.eps
x_down = ((np.clip(x_test - epsilon, 0, 1)-mean)/std).astype("float32")
x_up = ((np.clip(x_test + epsilon, 0, 1)-mean)/std).astype("float32")
success4 = success3 = success2 = success1 = 0
start_time = time.perf_counter()
for i in candidates:
    res = nn.test(x_down[i], x_up[i], y_test[i], True)
    if res >= 1:
        success1 += 1
    if res >= 2:
        success2 += 1
    if res >= 3:
        success3 += 1
    if res == 4:
        success4 += 1
elapsed = time.perf_counter() - start_time
print(f"For epsilon = {epsilon}, {success1} image(s) certified over {len(candidates)} candidates ({success1 / len(candidates) * 100}%) in {elapsed / len(candidates) * 1000}ms per image.")
