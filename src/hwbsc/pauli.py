import numpy as np

def paulis_commute(p1: np.ndarray,p2: np.ndarray) -> bool:
    """
    Check if two Pauli operators commute.

    Parameters
    ----------
    n : int
        The length of the repetition code.

    Returns
    -------
    parity_check_matrix : ndarray of shape (n, n)
        The parity check matrix of the repetition code.

    """
    
    counter = 0
    for i in range(len(p1)):       
            if (p1[i] != p2[i]) and p1[i] != 'I' and p2[i] != 'I':
                counter = counter +1 
          
    if counter%2 == 0:
        return True
    else:
        return False
    

# def check_commutivity(H:np.ndarray):
#     """
#     check commutivtiy relations for H
#     ...
#     """
#     for i in range(H.shape[0]):
#         for j in range(H.shape[0]):
#             if paulis_commute(H[i],H[j]) == False:
#                 return False
#     return True

def pauli2binary(pauli:np.ndarray):

    paulibin = np.zeros((len(pauli)*2), dtype=int)
    for i in range(len(pauli)):
        if pauli[i] == 'X':
            paulibin[i] = 1
        if pauli[i] == 'Y':
            paulibin[i]= 1
            paulibin[i+len(pauli)] = 1
        if pauli[i] == 'Z':
            paulibin[i+len(pauli)] = 1
    return paulibin


        
   
    
    

