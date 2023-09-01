import numpy as np
from hwbsc.code_utils import lookup_table, find_codewords, compute_distance, codeword_basis, compute_code_parameters, construct_generator_matrix,encode_permutations
from hwbsc.codes import hamming_code, repetition_code
from ldpc.mod2 import row_basis

def test_lookup_table():
    assert np.array_equal(lookup_table(hamming_code())[0,0,0][0], np.zeros(7))
    assert np.array_equal(lookup_table(hamming_code())[0,0,0][-1], np.ones(7))
  

def test_find_codewords():
    assert np.array_equal(find_codewords(hamming_code()), np.array([(0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 1, 1, 1), (0, 0, 1, 1, 0, 1, 0), (0, 0, 1, 1, 1, 0, 1), (0, 1, 0, 1, 0, 1, 1), (0, 1, 0, 1, 1, 0, 0), (0, 1, 1, 0, 0, 0, 1), (0, 1, 1, 0, 1, 1, 0), (1, 0, 0, 1, 0, 0, 1), (1, 0, 0, 1, 1, 1, 0), (1, 0, 1, 0, 0, 1, 1), (1, 0, 1, 0, 1, 0, 0), (1, 1, 0, 0, 0, 1, 0), (1, 1, 0, 0, 1, 0, 1), (1, 1, 1, 1, 0, 0, 0), (1, 1, 1, 1, 1, 1, 1)]))
    assert np.array_equal(find_codewords(repetition_code(3)), np.array([[0,0,0],[1,1,1]]))
    
def test_compute_distance():
    assert compute_distance(hamming_code()) == 3
    assert compute_distance(repetition_code(3)) == 3
    assert compute_distance(repetition_code(5)) == 5

def test_codeword_basis():
    assert np.array_equal(codeword_basis(find_codewords(repetition_code(3)))[0], [1, 1, 1])
    assert np.array_equal(codeword_basis(find_codewords(hamming_code())), [[1,0,0,1,0,0,1],[0,1,0,1,0,1,1],[0,0,1,1,0,1,0],[0,0,0,0,1,1,1]])

    assert np.array_equal(np.sort(codeword_basis(find_codewords(hamming_code())), axis=0),np.sort(row_basis(find_codewords(hamming_code())),axis=0))

def test_compute_code_parameters():
    assert np.array_equal(compute_code_parameters(repetition_code(3)), [3,1,3])
    assert np.array_equal(compute_code_parameters(hamming_code()), [7,4,3])

def test_construct_generator_matrix():
    assert np.array_equal(construct_generator_matrix(repetition_code(3)), [[1],[1],[1]])
    assert np.array_equal(np.mod(np.dot(repetition_code(3), construct_generator_matrix(repetition_code(3))),2), [[0],[0],[0]])                      
    assert np.array_equal(np.mod(np.dot(hamming_code(), construct_generator_matrix(hamming_code())),2), [[0,0,0,0],[0,0,0,0],[0,0,0,0]])                      

def test_encode_permutations():
    assert np.array_equal(find_codewords(hamming_code()),encode_permutations(construct_generator_matrix(hamming_code()),4))