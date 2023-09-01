import pytest
import numpy as np
from hwbsc.repetition_code import repetition_code_encoder
from hwbsc.repetition_code import repetition_code_majority_vote_decoder

def test_repetition_code_encoder_input_type():
    ###input type check types tests here.

    NotImplemented

def test_repetition_code_encoder():
    ###write tests here for the function.
    assert np.array_equal(repetition_code_encoder(3,np.array([1,0,1,0,1,0])), np.array([[1,1,1],[0,0,0],[1,1,1],[0,0,0],[1,1,1],[0,0,0]]))
    assert np.array_equal(repetition_code_encoder(4,np.array([1,0,1,0])), np.array([[1,1,1,1],[0,0,0,0],[1,1,1,1],[0,0,0,0]]))

    NotImplemented

def test_repetition_code_majority_vote_decoder():
    assert np.array_equal(repetition_code_majority_vote_decoder(np.array([[1,0,1],[0,0,0],[1,1,1],[0,1,0],[1,1,1],[0,1,0]])),np.array([[1,1,1],[0,0,0],[1,1,1],[0,0,0],[1,1,1],[0,0,0]]))
    assert np.array_equal(repetition_code_majority_vote_decoder(np.array([[1,1,1,0,0],[0,1,0,0,0],[1,1,1,0,0],[0,1,0,1,1],[1,1,1,0,1],[0,1,1,1,0]])),np.array([[1,1,1,1,1],[0,0,0,0,0],[1,1,1,1,1],[1,1,1,1,1],[1,1,1,1,1],[1,1,1,1,1]]))