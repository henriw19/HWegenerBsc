import numpy as np
from hwbsc.stabilizers_utils import *

linear = np.array([[1,1,0],[0,1,1]])

hypgraph = hypergraph_prod(linear)

print(hypgraph[0:6,13:26])