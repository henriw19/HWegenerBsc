from hwbsc.stabilizers_utils import *
from hwbsc.stab_codes import *
from hwbsc.codes import *
from ldpc.mod2 import rank, nullspace, row_basis, row_span

H_perfect = np.array([ 'XZZXI', 
                       'IXZZX', 
                       'XIXZZ', 
                       'ZXIXZ' ])

H_perfect_bin = np.zeros((4,10),dtype=int)
for i, stabilizers in enumerate(H_perfect):
    H_perfect_bin[i] = pauli2binary(stabilizers)



logicals = find_log_op(steane_code())

print(logicals)


