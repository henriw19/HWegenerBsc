import numpy as np

def repetition_code_encoder(repetition_code_length:int, binary_string: np.ndarray) -> np.ndarray:
    """
    Encodes a binary string into an n-bit repetition code

    Parameters
    ----------
    repetition_code_length : int
        The number of times each bit should be repeated in the encoding.
    
    binary_string : array_like
        A 1D array-like object containing the binary array to be encoded.
        The elements of the array must be either 0 or 1.
    
    Returns
    -------
    encoding_output : numpy.ndarray
        A 1D NumPy array representing the encoded binary array. The elements
        of the array are either 0 or 1, with each original bit repeated n times.
    
    """

    ### write code to check the input types are valid (see previous skeleton code for examples).
    ## Also write tests to check whether the input types are correctly validated.

    # check that the input is a numpy array
    if not isinstance(binary_string, np.ndarray):
        raise TypeError("Input must be a numpy array.")

    # check that the input contains only 0s and 1s
    for bit in binary_string:
        if(bit != 0 and bit != 1):
            raise ValueError("Input must contain only 0s and 1s.")
        
    #check that repetition_code_length is an integer larger than 1
    if repetition_code_length < 2:
        raise TypeError("Input must be an integer larger than 1")
    
    if not isinstance(repetition_code_length, int):
        raise TypeError("Input must be an integer")

    binary_string_length = len(binary_string)
    encoding_output = np.zeros(shape=(binary_string_length, repetition_code_length), dtype=int)

    #write the rest of the function here.
    for i in range(0,len(binary_string)):
        encoding_output[i,:] = binary_string[i] 

    return encoding_output

def repetition_code_majority_vote_decoder(binary_encoding: np.ndarray) -> np.ndarray:
    """
    Decodes a repetition code of length n using the majority vote decoding algorithm.

    Parameters
    ----------
    binary_encoding : ndarray
        A numpy array of length n containing the encoded binary values of the repetition code.

    Returns
    -------
    decoded_message : ndarray
        A numpy array of length n containing the decoded binary value of the repetition code.

    """

    ### write code to check the input types are valid (see previous skeleton code for examples).
    ## Also write tests to check whether the input types are correctly validated.

    # check that the input is a numpy array
    if not isinstance(binary_encoding, np.ndarray):
        raise TypeError("Input must be a numpy array.")

    # check that the input contains only 0s and 1s
    for m in range(len(binary_encoding)):
        for n in range(len(binary_encoding[m])):
            if(binary_encoding[m,n] != 0 and binary_encoding[m,n] != 1):
                raise ValueError("Input must contain only 0s and 1s.")

    decoded_output = binary_encoding
    #write the rest of the function here.
    
    l = len(binary_encoding)
                
    for i in range(l):
        
        if  np.sum(binary_encoding[i,:]) > len(binary_encoding[0])/2:
            decoded_output[i,:] = 1
        else:
            decoded_output[i,:] = 0    
    return decoded_output