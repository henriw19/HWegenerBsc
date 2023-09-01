import numpy as np
from itertools import product
from typing import List
from ldpc.mod2 import row_echelon, rank

def lookup_table(H:np.ndarray) -> dict:
    """
    Generates a lookup table connecting the syndromes to the corresponding errors for a given parity check matrix.

    Parameters
    ----------
    parity_check_matrix : ndarray of shape (m, n)
        The parity check matrix of the code, where `n` is the total number of bits and `m` is the number of
        parity checks.

    Returns
    -------
    lookup table : dict
        A dictionary with keys corresponding to the possible syndromes and values corresponding to the
        corresponding error vectors.

    """


    lookup_table = {}
    for error in product([0,1],repeat = len(H[0])):
        syndrome = tuple(H @ error %2)
        if syndrome not in lookup_table:    
            lookup_table[syndrome] = [tuple(error)]
        else:
            lookup_table[syndrome].append(tuple(error))

    return lookup_table

def find_codewords(H:np.ndarray) -> np.ndarray:
    """
    Returns the codewords of a given parity check matrix.
    
    Parameters
    ----------
    parity_check_matrix : ndarray of shape (m, n)
        The parity check matrix of the code, where `n` is the total number of bits and `m` is the number of
        parity checks.

    Returns
    -------
    codewords : ndarray of shape (c, n)
        The codewords of the parity check matrix, where `n` is the total number of bits and `c` is the number of
        codewords.

    """

    codewords = np.array(lookup_table(H)[tuple(np.zeros(len(H)))])
    return codewords

def codeword_basis(codewords: np.ndarray) -> np.ndarray:
    """
    Returns the basis of a given set of codewords
    
    Parameters
    ----------
    codewords : ndarray of shape (c, n)
        The codewords of a parity check matrix, where `n` is the total number of bits and `c` is the number of
        codewords.

    Returns
    -------
    codeword_basis : ndarray of shape (k, n)
        The basis of the given codewords, where `n` is the total number of bits and `k` is the length of the basis

    """

    basis = list()
    i = 0
    while np.sum(row_echelon(codewords)[0][i]) > 0:
        basis.append(row_echelon(codewords)[0][i])
        i += 1
    return basis
    

def compute_dimension(H:np.ndarray)->int:
    """
    Returns the dimension of a given parity check matrix
    
    Parameters
    ----------
    parity_check_matrix : ndarray of shape (m, n)
        The parity check matrix of the code, where `n` is the total number of bits and `m` is the number of
        parity checks.

    Returns
    -------
    dimension : int
        The dimension of the parity check matrix.   

    """

    k = len(H[0]) - rank(H) 

    return k

def compute_distance(H:np.ndarray)->int:
    """
    Returns the distance of a given parity check matrix
    
    Parameters
    ----------
    parity_check_matrix : ndarray of shape (m, n)
        The parity check matrix of the code, where `n` is the total number of bits and `m` is the number of
        parity checks.

    Returns
    -------
    distance : int
        The distance of the parity check matrix.   

    """

    smallest = len(H[0])
    zero_syndrome = tuple(np.zeros(len(H)))
    for i in range(len(lookup_table(H)[zero_syndrome])):
        weight = np.sum(lookup_table(H)[zero_syndrome][i])
        if ((weight < smallest ) &(weight != 0)):
            smallest = weight
    d = smallest
    return d

def compute_code_parameters(H:np.ndarray)->List[int]:
    """
    Returns the (n, k, d) code parameters of a given parity check matrix.
    
    Parameters
    ----------
    parity_check_matrix : ndarray of shape (m, n)
        The parity check matrix of the code, where `n` is the total number of bits and `m` is the number of
        parity checks.

    Returns
    -------
    code parameters : ndarray of shape (1, 3)
        The code parameters of the parity check matrix.
    
    """

    n = H.shape[1]

    k = compute_dimension(H)
    d = compute_distance(H)

    return [n,k, d]

def construct_generator_matrix(H:np.ndarray)->np.ndarray:
    """
    Returns the generator matrix for a given parity check matrix.

    Parameters
    ----------
    parity_check_matrix : ndarray of shape (m, n)
        The parity check matrix of the code, where `n` is the total number of bits and `m` is the number of
        parity checks.

    Returns
    -------
    generator_matrix : ndarray of shape (m, n)
        The generator matrix of the code.
    
    """

    G = np.array(codeword_basis(find_codewords(H))).T

    return G

def encode_permutations(G:np.ndarray, n:int) -> np.ndarray:
    codewords = []
    for bin in product([0,1],repeat = n):
        codewords.append(G @ bin %2)

    return codewords