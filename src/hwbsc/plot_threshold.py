import numpy as np
import random
import matplotlib.pyplot as plt
from hwbsc.decoders import LookupTableDecoder
from hwbsc.codes import hamming_code, repetition_code
from hwbsc.code_utils import compute_dimension

def sample_error(H:np.ndarray, p: float) -> np.ndarray:
    """
    Generate an array of zeros of the length of a given code with some elements
    randomly flipped to ones, representing a random error.

    Parameters
    ----------
    H : numpy.ndarray 
        The parity check matrix of the code
    p : float
        The probability of flipping a zero to a one. Must be a value between 0 and 1.

    Returns
    -------
    error : numpy.ndarray
        An array of zeros, with some elements randomly flipped to ones.
    
    """
 
    error = np.zeros(H.shape[1])
    for i in range(len(H[0])):
        if random.random() < p:
            error[i] = 1
    return error

def plot(code, x:int, n:int):
    """
    Plot the evaluation of different decoders to determine the pseudo-threshold

    Parameters
    ----------
    code : func 
        The decoding to be evaluated
    x : int
        The number of steps between 0 and 1 for the physical error rate
    n : int
        The number of evaluations for each physical error rate
    
    Returns
    -------
    0

    """

    physical_error = np.linspace(0,1,x)
    log_err_rate = np.zeros(x)
    Lookuptable = LookupTableDecoder(code)
    k= compute_dimension(code)

    for i in range(x):
        logical_error_counter = 0
        for j in range(n):
            error = sample_error(code,physical_error[i])
            syndrome = code @ error %2
            recovery = Lookuptable.decode(syndrome)
            if not np.array_equal(np.array(recovery), np.array(error)):
                logical_error_counter += 1
        log_err_rate[i] = logical_error_counter/n
    
    log_err_rate = 1-(1-log_err_rate)**(1/k)

    plt.semilogy(physical_error,log_err_rate, label = 'decoding plot')
    plt.grid(True)
    plt.plot(physical_error,physical_error, label = 'physical err. rate = log. err. rate')
    plt.xlabel('physical error rate')
    plt.ylabel('word error rate')
    plt.savefig("threshold")
    return 0



            


