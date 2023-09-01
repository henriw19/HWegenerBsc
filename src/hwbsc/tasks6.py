import numpy as np
from hwbsc.stab_codes import *
from hwbsc.codes import *
from hwbsc.stabilizers_utils import *

H_perfect = np.array([ 'XZZXI', 
                       'IXZZX', 
                       'XIXZZ', 
                       'ZXIXZ' ])

H_perfect_bin = np.zeros((4,10),dtype=int)
for i, stabilizers in enumerate(H_perfect):
    H_perfect_bin[i] = pauli2binary(stabilizers)


# Exercise 6.01
#print(CssCode(repetition_code(3),repetition_code(3)))

# Exercise 6.02

#print("H_x * H_z.T = 2* kron(H,H.T) = 0")

# Exercise 6.03

hamming_hypergraph = hypergraph_prod(hamming_code(3)).H
#print("The code parameters of the hypergraph product of the hamming code are:\n",StabCode(hamming_hypergraph))

# Exercise 6.04
for i in range(3,9,2):
    code = repetition_code(3)
    hypergraph = hypergraph_prod(code)
    #print(code.shape[0]**2+ code.shape[1]**2, StabCode(hypergraph).n)

#print("The relation is N = n^2 + m^2")

# Exercise 6.05

for i in range(3,9,2):
    code = hamming_code(3)
    k = code.shape[1] - rank(code)
    n = code.shape[1]
    hypergraph = hypergraph_prod(code).H
    N = hypergraph.shape[1]//2

    print(N - (2*code.shape[0]*n - 2 * (code.shape[0] - rank(code))), StabCode(hypergraph).k)

print("K = N - 2*m*n - 2k")


