import numpy as np
from hwbsc.binary import decimal2binary

from ldpc.mod2 import  row_basis, row_span


def repetition_code(length: int)->np.ndarray:
    """
    Returns the parity check matrix for a repetition code of length n.

    Parameters
    ----------
    n : int
        The length of the repetition code.

    Returns
    -------
    parity_check_matrix : ndarray of shape (n, n)
        The parity check matrix of the repetition code.

    """

    parity_check_matrix = np.zeros((length,length),dtype=int)

    ##write function here
    for i in range(length):
        parity_check_matrix[i][i] = 1
        parity_check_matrix[i][(i+1)%length] = 1

    return parity_check_matrix


def hamming_code(k):
    """

    Returns the parity check matrix for the (k, N) Hamming Code

    Parameters
    ----------
    k : int
        The number of parity bits

    Returns
    -------
    parity_check_matrix : ndarray of shape (k, N)
        The parity check matrix of the Hamming Code.

    """

    n = 2**k - k - 1  
    N = n+k 
    H = np.zeros((k, N),dtype=int)

    for i in range(1,k+1):
        for j in range(2**(i-1),N+1):
            if decimal2binary(j)[-i] == 1:
                H[i-1,j-1] = 1
        
    return H

def steane_code():
    hamming = np.array([[1,1,1,1,0,0,0],
                       [0,1,1,0,1,1,0],
                       [0,0,1,1,1,0,1]])
    steane_pcm = np.block([[ np.zeros(hamming.shape,dtype=int), hamming],
                         [hamming, np.zeros(hamming.shape,dtype=int) ]])
    
    return steane_pcm

def hamming_code15():
    """
    Returns the parity check matrix for the (7, 4) Hamming Code

    Returns
    -------
    parity_check_matrix : ndarray of shape (3, 7)
        The parity check matrix of the Hamming Code.

    """

    hamming_code_pcm = np.array([[1,0,0,0,1,1,1,0,0,0,0,1,1,1,1],[0,1,0,0,1,0,0,1,1,0,1,0,1,1,1],[0,0,1,0,0,1,0,1,0,1,1,1,0,1,1],[0,0,0,1,0,0,1,0,1,1,1,1,1,0,1]])
    return hamming_code_pcm


