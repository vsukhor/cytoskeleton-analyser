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

"""Exposes main entry point for processing the configuration settings.
"""

from pathlib import Path

from . import Config
from ..cells import CellType
from ..inout import Paths
from ..inout import Reader
from ..inout import set_logger


def process(
        new: bool,
        paths: Paths,
        rind: int,
        cell: CellType,
) -> Config:
    """Main entry point for processing the configuration settings.

    :param new: Switch between:
        (True) generating a new analysis, and
        (False) importing results of an earlier analysis from a file.
    :param paths: Paths for input/output operations.
    :param rind: Simulation run index.
    :param cell: CellType object.
    :return: Instance of configuration structure
    """

    paths.data_out.mkdir(parents=True, exist_ok=True)

    logger = set_logger(
        f"config {cell.typename} Plm_{cell.plmind} run {rind}",
        paths.data_out / 'config.log')
    Config.logger = logger

    if new:
        filename = paths.celltype_in / \
                   f"mainLog_{rind}_CskGen_{rind}_Plm_{cell.plmind}.txt"
        _preprocess(filename)
        with Reader(filename) as r:
            cf = Config().read(r)

        assert cf.run['index'] == rind

        cf.to_json(paths.data_out)

    else:
        cf = Config().from_json(paths.data_out)

    return cf


def _preprocess(
        filename: Path
) -> None:
    """ Auxiliary function for preprocessing the configuration file.
    :param filename: File name.
    """

    with open(filename, 'r') as f:
        txt = ''
        for r in f.readlines():
            rr = r.split(' ')
            if rr[0] == 'origin[]:':
                if len(rr) < 10:
                    f.close()
                    return
                txt += ' '.join(rr[:6]) + '\n'
                txt += ' '.join(rr[6:12]) + '\n'
                txt += ' '.join(rr[12:18]) + '\n'
                txt += ' '.join(rr[18:])
            elif rr[0] == 'origin_Centrosome[]:' or \
                    rr[0] == 'origin_Golgi[]:':
                txt += ' '.join(rr[:6]) + '\n'
                txt += ' '.join(rr[6:12]) + '\n'
                txt += ' '.join(rr[12:])
            else:
                txt += r
    with open(filename, 'w') as out:
        out.write(txt)
