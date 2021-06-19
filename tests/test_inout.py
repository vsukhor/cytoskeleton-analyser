import numpy as np
import pytest
from pathlib import Path

from cytoskeleton_analyser.inout import flatten_dict
from cytoskeleton_analyser.inout import read_to_dtype
from cytoskeleton_analyser.inout import read_to_ndarray
from cytoskeleton_analyser.inout import Reader


class TestReader:

    @pytest.fixture
    def text_data_path(self):

        return Path(__file__).parent/'data'/'reader_data.txt'

    @pytest.fixture
    def data(self, text_data_path):

        with open(text_data_path, 'r') as f:
            yield f

    def test_skip_lines(self, text_data_path):

        with Reader(text_data_path) as r:
            r.skip_lines(2)
            assert r.f.readline() == "third line\n"

    def test_value(self, text_data_path):

        with Reader(text_data_path) as r:
            r.skip_lines(3)
            assert r.value(2) == "data3"

    def test_array(self, text_data_path):

        ee = [6.e3, 7.e3, 8.e3]
        with Reader(text_data_path) as r:
            r.skip_lines(5)
            aa = r.array(2)
        for e, a in zip(ee, aa):
            assert e == a

    def test_array_float(self, text_data_path):

        ee = [6.e3, 7.e3, 8.e3]
        with Reader(text_data_path) as r:
            r.skip_lines(5)
            aa = r.array_float(2)
        for e, a in zip(ee, aa):
            assert e == a

    def test_array_int(self, text_data_path):

        ee = [6000, 7000, 8000]
        with Reader(text_data_path) as r:
            r.skip_lines(5)
            aa = r.array_float(2)
        for e, a in zip(ee, aa):
            assert e == a


binary_data = [3.54, 95.3, 74.12, 9, 8, 47]
binary_data_path = Path(__file__).parent/'data'/'read_to_data'


def test_read_to_dtype():

    with open(binary_data_path, 'rb') as f:
        for i in range(3):
            a = read_to_dtype(np.float32(), f, 1)
            assert a == pytest.approx(binary_data[i])
        for i in range(3):
            b = read_to_dtype(np.int32(), f, 1)
            assert b == binary_data[i+3]


def test_read_to_ndarray():

    with open(binary_data_path, 'rb') as f:
        a = np.empty(3, dtype=np.float32)
        for i in range(3):
            read_to_ndarray(a, f, i)
            assert a[i] == pytest.approx(binary_data[i])
        b = np.empty(3, dtype=np.int32())
        for i in range(3):
            read_to_ndarray(b, f, i)
            assert b[i] == pytest.approx(binary_data[i+3])


def test_flatten_dict():

    u = {'a': 1, 'b': {'c': 2, 'd': 3}, 'e': 4}
    v = flatten_dict(u)
    assert len(v) == 4
    assert v['a'] == 1
    assert v['b_c'] == 2
    assert v['b_d'] == 3
    assert v['e'] == 4
