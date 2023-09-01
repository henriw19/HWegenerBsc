import numpy as np

def decimal2binary(decimal_number: int)->np.ndarray:
    """
    Converts a decimal number to a binary number.
    
    Parameters
    ----------
    decimal_number : int or float
        The decimal number to be converted to binary.
    
    Returns
    -------
    binary_number : numpy.ndarray
        A 1D NumPy array representing the binary number. The elements
        of the array are either 0 or 1
    """

    # check that the input is an integer
    if(isinstance(decimal_number, int) == False):
        raise TypeError("Input must be an integer.")

    # check that the input is non-negative
    if(decimal_number < 0):
        raise ValueError("Input must be a non-negative integer.")

    if(decimal_number == 0):
        # strictly speaking, we need zero bits to represent zero. However, by
        # convention, we use a single length-0 bit string to represent 0.
        length_of_binary_number = 1
    else:
        # for all other non-zero numbers, we us np.ceil(np.log2(decimal_number+1)) to
        # compute the length of the bit string.
        length_of_binary_number = np.ceil(np.log2(decimal_number+1)).astype(int)
    
    # create a numpy array of zeros of the appropriate length.
    binary_number = np.zeros(length_of_binary_number, dtype = int)

    ### write the rest of the function here.
    
    decimal_number_mod = decimal_number
    for i in range(0,length_of_binary_number):
        if int(decimal_number_mod) - 2**(length_of_binary_number-i-1) < 0:
            binary_number[i] = 0
            
        else:
            binary_number[i] = 1
            decimal_number_mod -= 2**(length_of_binary_number-i-1)
            


    return binary_number

def binary2decimal(binary_number: np.ndarray) -> int:
    """
    Converts a binary number to a decimal number.
    
    Parameters
    ----------
    binary_number : array_like
        A 1D array-like object containing the binary number to be converted.
        The elements of the array must be either 0 or 1.
    
    Returns
    -------
    decimal_number : int
        The decimal representation of the binary number.
    """
    # check that the input is a numpy array
    if not isinstance(binary_number, np.ndarray):
        raise TypeError("Input must be a numpy array.")

    # check that the input contains only 0s and 1s
    for bit in binary_number:
        if(bit != 0 and bit != 1):
            raise ValueError("Input must contain only 0s and 1s.")

    ### write the rest of the function here.
    
    l = len(binary_number)
    
    binary_number_mod = np.empty(l)
    for i in range(0,l):
        binary_number_mod[i] = binary_number[i]*2**(l-i-1)

    decimal_number = int(np.sum(binary_number_mod))

    return decimal_number

def errorweight(binary_error: np.ndarray) -> int:
    return int(np.sum(binary_error))