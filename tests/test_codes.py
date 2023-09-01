import numpy as np
from hwbsc.codes import repetition_code, hamming_code

def test_parity_check():
    assert np.array_equal(repetition_code(3), np.array([[1,1,0],[0,1,1],[1,0,1]], dtype=int))
    assert np.array_equal(repetition_code(5), np.array([[1,1,0,0,0],[0,1,1,0,0],[0,0,1,1,0],[0,0,0,1,1],[1,0,0,0,1]], dtype=int))

print(hamming_code(4))