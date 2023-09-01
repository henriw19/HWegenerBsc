import numpy as np
from hwbsc.pauli import pauli2binary
from ldpc.mod2 import rank, nullspace, row_basis, row_span, row_echelon, row_basis
from itertools import product
import itertools


def check_commutativity(p1:np.ndarray,p2:np.ndarray):
    m = int(len(p1))

    n = int(m/2)
    unit = np.eye(n,dtype=int)
    I_Q = np.zeros((m,m),dtype=int)
    
    I_Q[:n, n:] = unit
    I_Q[n:, :n] = unit
    return int((p1 @ I_Q )@ p2.T %2)


def compute_syndrome(H:np.ndarray, errors:np.ndarray):
    """
    Compute the syndrome for a list of errors and a given check matrix.

    Parameters:
    -----------
    H : np.ndarray
        check matrix
    
    errors : np.ndarray
        errors
    
    Returns:
    --------
    syndrome : np.ndarray
        syndrome of length number of stabilizers

    """
    syndrome = np.zeros(H.shape[0],dtype=int)
    for error in errors:
        for i,stabilizer in enumerate(H):
            syndrome[i] = stabilizer @ I_Q(H.shape[1]) @ error.T %2

    return syndrome


def check_H_commutativity(H:np.ndarray):
    for stabilizer1 in H:
        for stabilizer2 in H:
            if check_commutativity(stabilizer1,stabilizer2) == 1:
                return False
    return True

def I_Q(m:int):
    """
    Create I_Q matrix with size mxm

    Parameters:
    -----------
    m : int
        size of the matrix
    
    Returns:
    --------
    I_Q : np.ndarray
        I_Q matrix
    
    """
    n = int(m/2)
    unit = np.eye(n,dtype=int)
    I_Q = np.zeros((m,m),dtype=int)
    I_Q[:n, n:] = unit
    I_Q[n:, :n] = unit
    return I_Q

def find_logs2(H):
    iterations = np.array(list(itertools.product([0,1], repeat = H.shape[1])))
    print(len(iterations))
    rowspan = row_span(H)
    print(rowspan)
    allcomb = []
    for bin_vector in iterations:
        allcomb.append(bin_vector @ rowspan) % 2
    x1 = np.hstack((np.kron(nullspace(H),allcomb)),np.zeros((H.shape[0],H.shape[0])))
    return x1



def find_log_op(H:np.ndarray):
    """
    Find the logical operators for a given check matrix

    Parameters:
    -----------
    H : np.ndarray
        check matrix
    
    Returns:
    --------
    log_op : np.ndarray
        Logical operators for the given check matrix
    
    """

    # condition 1, map to zero-syndrome 
    hx = H[H.shape[0]//2:,:H.shape[1]//2]
    hz = H[:H.shape[0]//2,H.shape[1]//2:]
    nspace = np.block([[np.zeros(nullspace(hz).shape),nullspace(hz)],[nullspace(hx),np.zeros(nullspace(hx).shape)]])

    I_Q1 = I_Q(nspace.shape[1])
    nspaceswapped = (I_Q1 @ nspace.T).T %2
    log_op = set()

    # condition 2, not a stabilizer
    logicals = []
    lamda = np.copy(H)
    ranklamda = rank(lamda)
    for potential_log in nspaceswapped:
        temp = np.vstack((lamda,potential_log))
        if rank(temp) > rank(lamda):
            lamda = temp
            logicals.append(potential_log)

    # convert set back to array
    log_op = np.array(list(logicals),dtype=int)
 
    return log_op

def find_log_op_color(H:np.ndarray):
    """
    Find the logical operators for a given check matrix

    Parameters:
    -----------
    H : np.ndarray
        check matrix
    
    Returns:
    --------
    log_op : np.ndarray
        Logical operators for the given check matrix
    
    """

    # condition 1, map to zero-syndrome 
    nspace = nullspace(H)

    log_op = set()

    # condition 2, not a stabilizer
    logicals = []
    lamda = np.copy(H)
    ranklamda = rank(lamda)
    for potential_log in nspace:
        temp = np.vstack((lamda,potential_log))
        if rank(temp) > rank(lamda):
            lamda = temp
            logicals.append(potential_log)

    # convert set back to array
    log_op = np.array(list(logicals),dtype=int)
 
    return log_op

def log_op_basis(H:np.ndarray):
    """
    Find the logical operator basis for a given check matrix

    Parameters:
    -----------
    H : np.ndarray
        check matrix
    
    Returns:
    --------
    basis : np.ndarray
        basis for the logical operators for the given check matrix
    
    """

    log_op = find_log_op(H)
    basis = row_basis(log_op)
    return basis

def stabmat2css(H:np.ndarray):
    """
    Find the logical operator basis for a given check matrix, returns error if the matrix cannot be transformed

    Parameters:
    -----------
    H : np.ndarray
        check matrix
    
    Returns:
    --------
    basis : np.ndarray
        basis for the logical operators for the given check matrix
    
    """
    #swap upper half of the matrix with the lower half

    row_echelon_H = row_echelon(H)[0]
    rows, cols = row_echelon_H.shape
    row_echelon_H[:rows//2], row_echelon_H[rows//2:] = row_echelon_H[rows//2:], row_echelon_H[:rows//2].copy()
    row_echelon_H[:, :cols//2], row_echelon_H[:, cols//2:] = row_echelon_H[:, cols//2:], row_echelon_H[:, :cols//2].copy()
    double_echelon = row_echelon(row_echelon_H)[0]
    swapped = np.zeros(double_echelon.shape,dtype=int)
    swapped[:rows//2,cols//2:cols] = double_echelon[:rows//2,:cols//2]
    swapped[rows//2:,:cols//2] = double_echelon[rows//2:,cols//2:]

    if check_css(swapped) == False:
        raise AssertionError("The matrix cannot be transformed to a CSS matrix.")

    return swapped


def check_css(H:np.ndarray):
    """
    Find the logical operator basis for a given check matrix, returns error if the matrix cannot be transformed

    Parameters:
    -----------
    H : np.ndarray
        check matrix
    
    Returns:
    --------
    basis : np.ndarray
        basis for the logical operators for the given check matrix
    
    """

    n = H.shape[1]
    for i in range(H.shape[0]):
        for j in range(int(n/2)):
            if sum(H[i,:int(n/2)]) != 0 and sum(H[i,int(n/2):]) != 0:
                return False
    return True


def compute_distance(H:np.ndarray):
    """
    Compute the distance for a given check matrix

    Parameters:
    -----------
    H : np.ndarray
        check matrix
    
    Returns:
    --------
    min_weight : np.ndarray
        weight of minimum weight stabilizer or logical operator which corresponds to the distance of H
    
    """
    span = row_span(np.concatenate((H,log_op_basis(H)),axis=0))
    
    min_weight = len(H)
    for i in range(span.shape[0]):
        if np.sum(span[i]) < min_weight and np.sum(span[i]) != 0:
            min_weight = np.sum(span[i])
        
    return min_weight


def hypergraph_product(classical_code:np.ndarray):
    """
    Compute the hypergraph product for a classical parity check matrix

    Parameters:
    -----------
    H : np.ndarray
        parity check matrix
    
    Returns:
    --------
    hypergraph_product : np.ndarray
        the hypergraph product of the classical parity check matrix
    
    """
    m = classical_code.shape[0]
    n = classical_code.shape[1]

    # hx = np.hstack((np.kron(classical_code,np.eye(n)), np.kron(np.eye(m), classical_code.T) ))
    # hz = np.hstack((np.kron(np.eye(n),classical_code),np.kron(classical_code.T,np.eye(m))) )
    hx = np.hstack((np.kron(np.eye(n), classical_code), np.kron(classical_code.T, np.eye(m)) ))
    hz = np.hstack((np.kron(classical_code, np.eye(n),),np.kron(np.eye(m), classical_code.T)) )
    
    
    return np.vstack((np.hstack((np.zeros(hz.shape),hz)),np.hstack((hx,np.zeros(hx.shape))))).astype(int)


def hexa_color_pcm_generator(m,n):
    """
    n is number of horizontal faces
    m is number of vertical faces
    m must be a int divisible by 6
    n must be smaller than m and even
    """
    if m%6 != 0:
        raise ValueError("m must be a int divisible by 6")
    
    if n%2 != 0 or n>= m:
        raise ValueError("n must be smaller than m and even")
    
    num_faces = n*m
    num_nodes = 2* num_faces
    facelist = []
    pcm = []
    for horface in range(n):
        offset = 0
        offset2 = 0
        if horface%2 != 0:
            offset = n
        
        if horface%(n-1) == 0 and horface != 0:
            offset2 = -n

        for verface in range(m):
            offset3 = 0
            if verface%(m-1) == 0 and verface != 0:
                offset3 = -2*(m * n)
            
            stab = np.zeros(6)
            stab[0] = 2*n*verface + horface + offset        
            stab[1] = stab[0] + 1 + offset2  
            if horface%2 != 0:
                stab[2] = stab[0] + n + offset3                    
                stab[3] = stab[0] + n + 1 + offset2 + offset3
            else:
                stab[2] = stab[0] + n                    
                stab[3] = stab[0] + n + 1 + offset2 
            stab[4] = stab[0] + 2*n + offset3
            stab[5] = stab[0] + 2*n + 1 + offset2 + offset3

            binstab = np.zeros(num_nodes).astype(int)
            binstab[stab.astype(int)] = 1
            
            facelist.append(stab[0].astype(int))
            pcm.append(binstab)
            
    return np.array(pcm), facelist