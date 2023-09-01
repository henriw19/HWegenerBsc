from hwbsc.pauli import pauli2binary
from hwbsc.stabilizers_utils import check_H_commutativity, compute_distance, find_log_op, log_op_basis, check_css, stabmat2css
from hwbsc.stab_codes import StabCode, SteaneCode, CssCode
from ldpc.mod2 import rank, nullspace, row_basis, row_span, row_echelon
import numpy as np
from hwbsc.codes import steane_code, hamming_code
from itertools import product


# Exercise 5.03
H_perfect = np.array([ 'XZZXI', 
                       'IXZZX', 
                       'XIXZZ', 
                       'ZXIXZ' ])

H_perfect_bin = np.zeros((4,10),dtype=int)
for i, stabilizers in enumerate(H_perfect):
    H_perfect_bin[i] = pauli2binary(stabilizers)
print("The five-qubit code in GF2:\n",H_perfect_bin)

print("Check wether the stabilizers in H commute:\n",check_H_commutativity(H_perfect_bin))


# Exercise 5.06

print("Stabilizer commutativtiy check:\n",check_H_commutativity(steane_code()))

print("The number of logical qubits encoded by the Steane code:\n",SteaneCode().k)


# Exercise 5.07

print("The Steane code has a distance of:\n",SteaneCode().d)

fttbin = np.zeros((2,8),dtype=int)
for i, stabilizers in enumerate(['XXXX',
                                 'ZZZZ']):
    fttbin[i] = pauli2binary(stabilizers)
print("The code parameters of the [[4,2,2]] code are indeed:\n", StabCode(fttbin ,distance = True).code_parameters())


# Exercise 5.08

H = np.array([[0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
       [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
       [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1, 0],
       [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1],
       [0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0],
       [0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1],
       [1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
       [0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
       [0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
       [0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
       [0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
       [0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]).astype(int)

print("The CSS representation of H is given by:\n",stabmat2css(H))
print("The code parameters are:\n",StabCode(stabmat2css(H),name = 'H',distance = True))


# Exercise 5.09

print("The simplification is that the problem of decoding is seperated into 2 classical decoding problems")

