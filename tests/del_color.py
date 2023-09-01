
import sympy as sp
import numpy as np
from ldpc.mod2 import row_echelon, rank, nullspace
import scipy
from hwbsc.stab_codes import *
from hwbsc.surface_code import *
from hwbsc.color_code import hexa_color_pcm_generator, color_graph
from hwbsc.stabilizers_utils import compute_distance,check_H_commutativity, check_commutativity, find_log_op, check_css, hexa_color_pcm_generator,find_log_op_color

import galois

mat =  [[1, 1, 0, 0, 1, 1, 0, 0, 1, 1 ,0 ,0 ,0 ,0 ,0 ,0 ,0,0, 0, 0 ,0 ,0 ,0 ,0],
        [0 ,1 ,1 ,0 ,0 ,1 ,1 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0, 0, 0 ,0 ,0,1 ,1, 0],
        [0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 1, 0, 0, 0 ,1 ,0 ,0 ,0 ,0 ,0, 0, 0, 0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0],
        [0, 0 ,0, 0 ,1 ,1, 0 ,0 ,1 ,1, 0, 0, 0 ,0 ,0, 0, 1, 1 ,0 ,0, 1 ,1 ,0, 0],
        [0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0 ,0 ,0 ,0 ,0, 0, 0 ,0 ,0 ,0 ,0 ,0],
        [0, 0, 0, 0, 0, 0, 1, 0, 0 ,0 ,0 ,0 ,0, 0 ,0 ,0 ,0 ,0 ,1 ,0 ,0, 0, 1, 0,],
        [0, 0, 0, 0, 0, 0, 0 ,1 ,0 ,0 ,0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1],
        [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0 ,0 ,0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0 ,1 ,1 ,0, 0, 0 ,0 ,0 ,0],
        [0, 0 ,0 ,0, 0, 0, 0, 0, 0, 0 ,0 ,1 ,0, 1, 1 ,0 ,0 ,0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0 ,0 ,0 ,0 ,0 ,0 ,1 ,0, 0 ,1 ,1, 0, 0, 1, 1, 0, 0, 1],
        [0, 0, 0, 0, 0, 0 ,0, 0, 0, 0, 0, 0 ,0 ,1, 1 ,0 ,0 ,1, 1, 0 ,0, 1, 1, 0],
        [0 ,0, 0, 0, 0, 0 ,0, 0, 0, 0, 0, 0 ,0 ,0, 1 ,1, 1, 1, 0, 0 ,0, 0 ,0, 0],
        [0 ,0, 0 ,0 ,0, 0 ,0, 0, 0, 0, 0, 0 ,0 ,0, 0, 1 ,0, 0, 0, 1 ,0, 0 ,0 ,0],
        [0 ,0, 0 ,0 ,0, 0 ,0, 0, 0, 0, 0, 0 ,0 ,0, 0, 0 ,1, 1, 1, 1 ,1, 1 ,1 ,1]]
        # [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0]]
subpcm =   [[1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0 ,0 ,0 ,0 ,0 ,0 ,1 ,1 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0, 0],
            [0 ,0 ,0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0 ,0 ,0 ,0 ,0 ,0 ,0 ,1 ,1 ,1 ,1 ,0 ,0 ,1 ,1 ,0 ,0 ,0 ,0 ,0 ,0],
            [1 ,1 ,0 ,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0],
            [0, 0, 0, 0, 0 ,1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,1 ,0 ,1 ,1 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0],
            [0 ,1 ,1 ,0 ,0 ,1 ,1 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,1 ,1 ,0],
            [0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,1 ,0 ,0 ,0 ,1 ,0 ,0 ,0 ,0 ,0],
            [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0],
            [0 ,0 ,0 ,1 ,0 ,0 ,0 ,1 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0],
            [0 ,0 ,0 ,1 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,1 ,0 ,0 ,0 ,1],
            [0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,1 ,0 ,1 ,0 ,0 ,1 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ,0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1],
            [1 ,0 ,0 ,1 ,1 ,0 ,0 ,1 ,0 ,0 ,0, 0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,1 ,0 ,0 ,1]]

# GF = galois.GF(2)

# H = np.array([[1,1,0,0,0,1,0],
# [0,0,1,1,1,0,1],[0,0,0,0,1,1,0],[0,0,0,0,1,1,1],[0,0,0,0,0,0,1],[1,0,1,1,1,0,1]])
# H = GF(H)

# syn = [1,1,1,0,1,0]

# P, L, U = H.plu_decompose()
# print(np.array(P)@np.array(syn)%2)


# # print(SurfCode(hyprep).logicals)
array = [0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0
 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
 0 0 0 0 0 0 1 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0
 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
 0 1]