import pytest
from hwbsc.demo import square
from hwbsc.demo import cube

def test_square_positive():
    assert square(3) == 9
    assert square(5) == 25
    assert square(10) == 100

def test_square_negative():
    assert square(-3) == 9
    assert square(-5) == 25
    assert square(-10) == 100

def test_square_zero():
    assert square(0) == 0

def test_cube_positive():
    assert cube(3) == 27
    assert cube(5) == 125
    assert cube(10) == 1000

def test_cube_negative():
    assert cube(-3) == -27
    assert cube(-5) == -125
    assert cube(-10) == -1000

def test_cube_zero():
    assert cube(0) == 0