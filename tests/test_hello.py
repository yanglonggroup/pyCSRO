import pytest
from pkgname import hello

def test_add():
    x=5
    y=6
    assert 11 == hello.add(x,y)