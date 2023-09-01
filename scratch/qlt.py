import numpy as np

from hwbsc.qubit_lookuptable import LookupTableDecoder

H = np.array([ np.array(['X','X','X','X']), np.array(['Z','Z','Z','Z'])])

decoder = LookupTableDecoder(H)

syndrome = np.array([1,1])

decoding = decoder.decode(syndrome)

print(decoding)