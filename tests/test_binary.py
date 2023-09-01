import numpy as np
import pytest
from hwbsc.binary import decimal2binary
from hwbsc.binary import binary2decimal

def test_decimal2binary_input_type():
    with pytest.raises(TypeError):
        decimal2binary("string")  # string input should raise TypeError
    with pytest.raises(TypeError):
        decimal2binary(1.23)  # float input should raise TypeError
    with pytest.raises(TypeError):
        decimal2binary(np.array([1]))  # numpy array input should raise TypeError
    with pytest.raises(TypeError):
        decimal2binary(None)  # None input should raise TypeError

def test_decimal2binary_first_five_numbers():
    assert np.array_equal(decimal2binary(0), np.array([0], dtype=int))
    assert np.array_equal(decimal2binary(1), np.array([1], dtype=int))
    assert np.array_equal(decimal2binary(2), np.array([1, 0], dtype=int))
    assert np.array_equal(decimal2binary(3), np.array([1, 1], dtype=int))
    assert np.array_equal(decimal2binary(4), np.array([1, 0, 0], dtype=int))

def test_decimal2binary_joschka():
    # Test the first 10000 binary numbers by comparing against the numpy implementation
    for i in range(10000):
        expected = np.binary_repr(i)
        actual = decimal2binary(i)
        assert np.array_equal(actual, np.array([int(d) for d in expected]))


def test_binary2decimal_input_type_checks():
    with pytest.raises(TypeError):
        binary2decimal("string")  # string input should raise TypeError
    with pytest.raises(TypeError):
        binary2decimal(1.23)  # float input should raise TypeError
    with pytest.raises(TypeError):
        binary2decimal(None)  # None input should raise TypeError
    with pytest.raises(ValueError):
        binary2decimal(np.array([1, 2, 0]))  # input with values other than 0 and 1 should raise ValueError

def test_binary2decimal():
    
    assert binary2decimal(np.array([0]).astype(int)) == 0
    assert binary2decimal(np.array([1,0]).astype(int)) == 2
    assert binary2decimal(np.array([1,1,1]).astype(int)) == 7
    assert binary2decimal(np.array([1,1,0,0,0,0,1,1,0,1,0]).astype(int)) == 1562

def test_binary2decimal_joschka():

    # Test the first 10000 binary numbers
    for i in range(10000):
        binary = np.binary_repr(i)
        expected = int(binary, 2)
        actual = binary2decimal(np.array([int(d) for d in binary]))
        assert actual == expected


def test_decimal2binary():
    
  assert np.array_equal(decimal2binary(2), np.array([1,0]))
  assert np.array_equal(decimal2binary(7), np.array([1,1,1]))
  assert np.array_equal(decimal2binary(27), np.array([1,1,0,1,1]))
  assert np.array_equal(decimal2binary(1562), np.array([1,1,0,0,0,0,1,1,0,1,0]))
    
   