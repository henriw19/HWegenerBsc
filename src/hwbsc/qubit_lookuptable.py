import numpy as np
import itertools

from hwbsc.code_utils import compute_code_parameters
from hwbsc.binary import decimal2binary
from hwbsc.codes import hamming_code
from hwbsc.pauli import paulis_commute

# test PCM
H_0 = np.array([ ['X','X','X','X'], ['Z','Z','Z','Z']])
H_perfect = np.array([ ['X','Z','Z','X','I'], ['I','X','Z','Z','X'] , ['X','I','X','Z','Z'], ['Z','X','I','X','Z']])

test_error = np.array(['X','I','I','I'])


def compute_syndrome(H:np.ndarray, error:np.ndarray):
    s = np.zeros(H.shape[0],dtype=int)
    for i in range(H.shape[0]):
        
        if paulis_commute(H[i],error) == False:
            s[i] = 1
    return s

class LookupTableDecoder():

    def __init__(self, H:np.ndarray):
        """
        ...

        """

        #Create the lookup table in the constructor. This will prevent the lookup table from
        # being re-created every time the decode method is called.
        self.H = H
        self.n = H.shape[1]
        self.k = H.shape[1] - H.shape[0]
        self.d = NotImplemented
        combinations = list(itertools.product([0, 1], repeat=8))
        paulis = ['I','X','Y','Z']
        errors = list(itertools.product(paulis, repeat=self.n))
        self.lookup_table = {}
        for i in range(4**self.n):

            syndrome = tuple(compute_syndrome(H, errors[i]))
            if syndrome not in self.lookup_table:    
                self.lookup_table[syndrome] = [errors[i]]
            else:
                self.lookup_table[syndrome].append(errors[i])

   

    def decode(self, syndrome:np.ndarray)->np.ndarray:
        """
        ...
        """
        decoding = self.lookup_table[tuple(syndrome)]
        return np.array(decoding)
    
def weight(p:np.ndarray):
    #calculate the weight of a pauli operator
    weight = 0
    for i in range(len(p)):
        if p[i] != 'I':
            weight += 1
    return weight
        


def find_log_op(H:np.ndarray):
    # condition 1, map to zero-syndrome   
    log_op = Lookuptable.lookup_table[tuple([0,0])]

    # condition 2, not stabilizers
    for i in range(len(log_op)-2):
        if (np.array_equal(log_op[i],np.array(['X','X','X','X'])) or np.array_equal(log_op[i],np.array(['Y','Y','Y','Y'])) or np.array_equal(log_op[i],np.array(['Z','Z','Z','Z']))):
            del log_op[i]

    # condition 3, anti-commutation
    filtered_log_op = [x for x in log_op if any(y for y in log_op if paulis_commute(x,y) == False)]

    return filtered_log_op

Lookuptable = LookupTableDecoder(H_0)
print(Lookuptable.lookup_table)




