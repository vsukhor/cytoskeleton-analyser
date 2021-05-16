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

"""Class for project-wide structure of data summary.
"""

from __future__ import annotations
import json
from logging import Logger
from pathlib import Path
from typing import Optional, Sequence

from .inout import underscored


class Summary:
    """Data summary used in reports of extracted cytoskeleton features.

    The summary is may hold several related features in a ``data`` list
    container. Its elements are feature-specific dictionaries structured
    as follows:
    'name': feature name
    'model': model structure desctibing the best fit available or None
    'avg': (optional) data average
    'std': (optional) data standdard deviation
    """

    logger: Logger

    def __init__(
            self,
            logger: Optional[Logger] = None
    ):

        self.data = None
        if logger is not None:
            self.logger = logger

    def init(
            self,
            names: list[str],
            ibest: list[int],
            fits: list[Sequence],
            avgs: Optional[list[float]] = None,
            stds: Optional[list[float]] = None,
    ) -> Summary:

        bms = [fits[i][ibest[i]].summary() if len(fits[i]) else None
               for i in range(len(names))]

        self.data = [{'name': n, 'model': m} for n, m in zip(names, bms)]
        if avgs is not None:
            self.data = [{**d, 'avg': a} for a, d in zip(avgs, self.data)]
        if stds is not None:
            self.data = [{**d, 'std': s} for s, d in zip(stds, self.data)]

        return self

    def add_items(
            self,
            d: Sequence[dict],
    ) -> None:
        """Include additional elements to the current structure.

        :param d: Data to be included.
        """

        k = len(self.data)
        inds = range(k)
        d *= k

        for i in inds:
            for k, v in d[i].items():
                if k not in self.data[i]:
                    self.data[i] = self.data[i] | {k: v}

    def save(
            self,
            path: Path,
            name: str
    ) -> None:
        """Save the summary structure to a .json data file.

        :param path: File position on the data storage device.
        :param name: File name without suffix.
        """

        fname = (path / underscored(name)).with_suffix('.json')
        self.logger.info('Saving summary: ' + str(fname) + '\n')

        with open(fname, 'w') as f:
            json.dump(self.data, f)

    def read(
            self,
            path: Path,
            name: str
    ) -> Optional[Summary]:
        """Read the summary structure from a .json data file.

        :param path: File position on the data storage device.
        :param name: File name without suffix.
        :return: Initialized nstance of of this class.
        """

        fname = (path / underscored(name)).with_suffix('.json')
        if not Path(fname).is_file():
            self.logger.info('Skipping Fit Summary: file ' +
                             str(fname) + ' is not available.\n')
            return None
        self.logger.info('Reading summary: ' + str(fname) + '\n')
        with open(fname, 'r') as f:
            self.data = json.load(f)

        return self
