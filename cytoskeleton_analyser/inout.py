# Copyright (c) 2021 Valerii Sukhorukov. All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# 1. Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the distribution.
# 3. Neither the name of the copyright holder nor the
#    names of its contributors may be used to endorse or promote products
#    derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDER ''AS IS'' AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
# USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# -----------------------------------------------------------------------------

""" Input-output functions and classes.
"""

import csv
import gzip
import logging
import pathlib as pl
import struct
import sys
import zipfile
from typing import BinaryIO, Final, Optional, Sequence, TextIO, Union

import numpy as np

from .cells import CellType


class Paths:
    """A set of project-specific paths for in/out data storage.
    """

    def __init__(
            self,
            data_in: pl.PurePath,
            data_out: pl.PurePath,
            cell: CellType,
            rind: Optional[int] = None,
    ):
        """
        :param data_in: Full path to input data folder (common for all
            cell types).
        :param data_out: Cell type- and run-specific path to output data
            folder.
        :param cell: Instance of CellType.
        :param rind: Simulation ran index.
        """

        self.inner = pl.PurePath(cell.typename, f'plm{cell.plmind}')
        if rind is not None:
            self.inner /= f"csk{rind}"
        self.data_in = pl.Path(data_in)
        self.celltype_in = self.data_in / f"{cell.typename}"
        self.plasma_in = self.celltype_in / f"plm{cell.plmind}"
        self.run = self.plasma_in / f"csk{rind}" \
            if rind is not None else None
        self.data_out = pl.Path(data_out / self.inner)

    @staticmethod
    def ensure(p: pl.Path) -> None:
        """Ensure that path ``p`` exists.

        :param p: Path name.
        """

        pl.Path(p).mkdir(parents=True, exist_ok=True)


class Reader:
    """Read in numerical configuration data from a text file.
    """

    def __init__(self,
                 fname: pl.Path):

        self.f: Final[TextIO] = open(fname, 'r')

    def __enter__(self):

        return self

    def __exit__(self, exc_type, exc_value, traceback):

        if self.f is not None:
            self.f.close()

    def skip_lines(self, n: int):
        """ Skip 'n' lines.
        """
        last_line = ''
        for _ in range(n):
            last_line = self.f.readline()

        return last_line

    def value(self, i: int):

        u = self.f.readline()
        return u.split()[i]

    def array(self, i: int):

        u = self.f.readline()
        return [float(x) for x in u.split()[i:]]

    def array_float(self, i: int) -> list[float]:

        return [float(x) for x in self.array(i)]

    def array_int(self, i: int) -> list[int]:

        return [int(x) for x in self.array(i)]


def read_to_dtype(
        v: Union[np.number, np.ndarray],
        f: BinaryIO,
        n: int = 1,
) -> Union[np.number, np.ndarray]:
    """Read in (out of place) one or many numbers.

    Data are read from the descriptor 'f' of a binary file.
    :param v: Data item of the desired type.
    :param f: Binary file.
    :param n: How many numbers to rread in.
    :return: Data extracted.
    """

    nda = isinstance(v, np.ndarray)
    tp = v.dtype if nda else np.dtype(v)
    u = struct.unpack(str(n) + tp.char, f.read(n * tp.itemsize))

    return u[:n] if nda else u[0]


def read_to_ndarray(
        v: np.ndarray,
        f :BinaryIO,
        i: int,
) -> None:
    """Read in (in place) an numpy.ndarray of numbers.

    Data are read from the descriptor 'f' of a binary file.
    :param v: Data placeholder.
    :param f: Binary file.
    :param i: Array index.
    """
    n = 1 if v.ndim == 1 else v.shape[1]
    v[i] = read_to_dtype(v[i], f, n)


def set_logger(
        name: str,
        file: pl.Path
) -> logging.Logger:
    """Convenience routine creating loggers used throughout the project.

    :param name: Logger name.
    :param: file: Full name (including path) for the log file.
    """

    import datetime

    logger = logging.getLogger(name)

    while len(logger.handlers):
        logger.removeHandler(logger.handlers[-1])

    logger.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(name)s - %(message)s')

    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setLevel(logging.DEBUG)
    stream_handler.setFormatter(formatter)

    file_handler = logging.FileHandler(file)
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(formatter)

    logger.addHandler(stream_handler)
    logger.addHandler(file_handler)

    logger.info('Starting time: ' +
                datetime.datetime.now().strftime("%d %b %Y %H:%M:%S"))
    logger.info('')

    return logger


def save_as_csv(
        fname: pl.Path,
        caption: Sequence[str],
        data: Sequence[np.ndarray],
) -> None:
    """Save a sequence of ndarrays as .csv data table.

    :param fname: Full path/file to use.
    :param caption: csv table caption.
    :param data: Data to be saved.
    """

    assert len(caption) == len(data)
    assert len(data) == 1 or \
           all([len(data[0]) == len(d) for d in data[1:]])

    opn = gzip.open if fname.suffix == 'gz' else open
    with opn(fname, 'wt', newline='') as f:
        w = csv.writer(f, delimiter=',')
        w.writerow(caption)
        for i in range(len(data[0])):
            w.writerow([d[i] for d in data])


def rmdir(d: pl.Path) -> None:
    """Remove directtory 'd'.
    """

    import shutil

    try:
        shutil.rmtree(d)
        print ('Deleted: '+ str(d))
    except OSError as e:
        print ("Error: %s - %s." % (e.filename, e.strerror))


def underscored(s: str) -> str:
    """Turn spaces in the string ``str`` into underscores.

    Used primarely for filename formatting.
    """

    return '_'.join(s.split(' '))


def unzip(path: pl.Path) -> bool:
    """Uncompress folder containing simulation raw data.
    """

    if not path.is_dir():
        zf = path.with_suffix('.zip')
        if zf.is_file():
            with zipfile.ZipFile(zf, 'r') as zip:
                zip.extractall(zf.parent)
            return True

    return False


def flatten_dict(
        d: dict,
        s: str = '',
        exc: Optional[list] = None
) -> dict:
    """Reformat a multi-layer dictionary into a flat one.

    :param d: Input dict.
    :param s: Prefix to be added to keys.
    :param exc: List of keys to exclude from the resulting dict.
    :return: Reformatted dictionary.
    """

    if exc is None:
        exc = []
    u = {}
    for k, v in d.items():
        if k in exc:
            u[k] = v
        else:
            if len(s):
                k = s + '_' + k
            if isinstance(v, dict):
                u |= flatten_dict(v, k)
            elif isinstance(v, list):
                u |= {k+'_'+str(i): vv for i, vv in enumerate(v)}
            else:
                u[k] = v
                print('')

    return u
