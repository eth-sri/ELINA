import numpy as np
from fppoly import *
from elina_abstract0 import *
from elina_manager import *

from ctypes import *
from ctypes.util import *

libc = CDLL(find_library('c'))
cstdout = c_void_p.in_dll(libc, 'stdout')

#encode
# 2p1 + p2 + phi <= p1' <= p3+p4+2phi
# p1 + 2phi <= p2' <= p3 + p4 + phi
# p3 + 2p4 + phi <= p3' <= 4p4 
# p2 + p4 <= p4' <= p4 + 2phi

lexpr_weights = np.array([[2.0,1.0,1.0],[1.0,0.0,2.0],[1.0,2.0,1.0],[1.0,1.0,0.0]])
np.ascontiguousarray(lexpr_weights, dtype=np.double)
print("shape ",lexpr_weights.shape[0])
print("before ", lexpr_weights)
lpp = (lexpr_weights.__array_interface__['data'][0]+ np.arange(lexpr_weights.shape[0])*lexpr_weights.strides[0]).astype(np.uintp) 

print("weights ", lexpr_weights[0])

uexpr_weights = np.array([[1.0,1.0,2.0],[1.0,2.0,1.0],[0.0,0.0,4.0],[0.0,1.0,2.0]])
upp = (uexpr_weights.__array_interface__['data'][0]+ np.arange(uexpr_weights.shape[0])*uexpr_weights.strides[0]).astype(np.uintp)

lexpr_cst = np.array([0.0,0.0,0.0,0.0])
np.ascontiguousarray(lexpr_cst, dtype=np.double)

uexpr_cst = np.array([0.0,0.0,0.0,0.0])
np.ascontiguousarray(uexpr_cst, dtype=np.double)

lexpr_dim = np.array([[0,1,4],[0,1,4],[2,3,4],[1,3,4]])
lexpr_dim = (lexpr_dim.__array_interface__['data'][0]+ np.arange(lexpr_dim.shape[0])*lexpr_dim.strides[0]).astype(np.uintp)

uexpr_dim = np.array([[2,3,4],[2,3,4],[0,1,3],[1,3,4]]) 
uexpr_dim = (uexpr_dim.__array_interface__['data'][0]+ np.arange(uexpr_dim.shape[0])*uexpr_dim.strides[0]).astype(np.uintp)

lexpr_size = np.array([3,3,3,3])
lexpr_size = lexpr_size.astype(np.uintp)
np.ascontiguousarray(lexpr_size, dtype=np.uintp)

uexpr_size = np.array([3,3,3,3])
uexpr_size = uexpr_size.astype(np.uintp)
np.ascontiguousarray(uexpr_size, dtype=np.uintp)


#print("type ", type(lexpr_size))

#encode
# 0.1 <= p1 <= 0.2
# 0.2 <= p2 <= 0.3
# 0.3 <= p3 <= 0.4
# 0.4 <= p4 <= 0.5
# 0.001 <= phi <= 0.002

inf = np.array([0.1,0.2,0.3,0.4,0.001])
sup = np.array([0.2,0.3,0.4,0.5,0.002])


man = fppoly_manager_alloc() 
element = fppoly_from_network_input_poly(man,0,4,inf,sup, lpp, lexpr_cst, lexpr_dim, lexpr_size, upp, uexpr_cst, uexpr_dim ,uexpr_size) 


#handle first layer
weights = np.array([[1.0,-1.0,1.0,-1.0],[-1.0,1.0,-1.0,1.0],[1.0,1.0,-1.0,-1.0],[-1.0,-1.0,1.0,1.0]])
np.ascontiguousarray(weights, dtype=np.double)
weights = (weights.__array_interface__['data'][0]+ np.arange(weights.shape[0])*weights.strides[0]).astype(np.uintp) 


biases = np.array([1.0,-1.0,-1.0,1.0])
np.ascontiguousarray(biases, dtype=np.double)


ffn_handle_first_relu_layer(man,element,weights,biases,4, 4)
elina_abstract0_fprint(cstdout,man,element,None)
