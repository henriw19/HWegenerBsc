import numpy as np
from ldpc.mod2 import rank
from hwbsc.binary import decimal2binary


class LookupTableDecoder():

    def __init__(self, H:np.ndarray):
        """
        Generate a lookup table for decoding a binary code with a given parity check matrix.

        Parameters
        ----------
        H : numpy.ndarray
            A 2D numpy array representing the parity check matrix of the code. Each row of the matrix
            represents a parity check equation in binary form.

        Returns
        -------
        dict
            A dictionary representing the lookup table for decoding the code. The keys to the dictionary 
            are the syndromes and the values are the corresponding errors that are most likely to have
            occurred.

        """

        #Create the lookup table in the constructor. This will prevent the lookup table from
        # being re-created every time the decode method is called.
        self.H = H
        self.n = H.shape[1]
        self.k = H.shape[1] - rank(H)
        self.d = NotImplemented

        self.lookup_table = {}
        for i in range(2**self.n):
            # Create all possible binary errors with length n using the decimal2binary function
            error = np.pad(decimal2binary(i), (self.n-len(decimal2binary(i)),0))

            # Calculate corresponding syndromes and fill the lookup table with the errors of least weight
            syndrome = tuple(H @ error %2)
            if syndrome not in self.lookup_table:    
                self.lookup_table[syndrome] = tuple(error)
            else:
                if np.sum(error) < np.sum(self.lookup_table[syndrome]):
                    self.lookup_table[syndrome] = tuple(error)


    def decode(self, syndrome:np.ndarray)->np.ndarray:
        """
        Decodes a given syndrome using the already generated lookup table

        Parameters
        ----------
        syndrome : np.ndarray
            A 1D numpy array representing the syndrome to be decoded
        
        Returns
        -------
        decoding : np.ndarray
            A 1D numpy array representing the proposed recovery operation for the error associated
            with the given syndrome.
        
        """
        decoding = self.lookup_table[tuple(syndrome)]
        return np.array(decoding[0])