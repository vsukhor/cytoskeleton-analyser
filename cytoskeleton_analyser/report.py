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

"""Base class for reports of microtubule characteristics.
"""

from abc import ABC
from collections import namedtuple
from logging import Logger
from pathlib import Path
from typing import Optional

import cytoskeleton_analyser.fitting as fit
from .histograms import Histogram
from .summary import Summary

Stats = namedtuple('Stats', 'avg std units')


class Report(ABC):

    """Abstract base class for creating attribute-specific reports.
    """

    name: str
    logger: Logger
    signature: str
    path_out: Path
    units: str
    is_polar: bool = False
    is_halfpolar: bool = False
    fits_sim: list[tuple] = []
    fits_exp: list[tuple] = []
    show: bool

    @classmethod
    def __create(
            cls,
            logger: Logger,
            path_out: Path,
            name: str,
            show: bool = True,
    ):

        cls.logger = logger
        fit.base.set_logger(logger)
        Histogram.logger = logger
        Summary.logger = logger

        cls.path_out = path_out
        cls.name = name
        cls.show = show

    @classmethod
    def summarize(
            cls,
            has: tuple[Histogram,
                       Optional[list[float]],
                       Optional[list[float]]],
            ibest: list[int]
    ) -> None:
        """Create and saave a report summary.

        :param has: Input data: a an instance of Histogram class and
            optionally data average and its standard deviation.
        :param ibest: Index of the best fitting model out of the set
            available in the Histogram class.
        """
        h, avgs, stds = has
        s = Summary()
        if h.experimental is None:
            s.init(['simulated'], ibest, [h.simulated.fits], avgs, stds)
        else:
            s.init(['simulated', 'experimental'], ibest,
                   [h.simulated.fits, h.experimental.fits], avgs, stds)
        s.save(cls.path_out, cls.name)

        cls.logger.info('')

    @classmethod
    def restore(
            cls,
            path: Path
    ) -> dict:
        """Restore full set of characteristics for the reported feature.

        Upon readig from the file provided by ```path``, instance
        of the Histogram class containing the data and eventual models
        is instantiated.

        :param path: Full (dir + specific name) name of the input file.
        :return: Summary of the data restored.
        """

        cls.logger.info('Restoring ' + cls.name + ':')

        summary = cls.log_stats(path)
        h = Histogram(cls.name).from_csv(path)
        h.simulated.polar = cls.is_polar
        h.simulated.halfpolar = cls.is_halfpolar
        cls._restore_fits(path, h)
        cls.plot(h)

        cls.logger.info('')

        return {'summary': summary,
                'h_sim_avg': h.simulated.avg}

    @classmethod
    def _restore_fits(
            cls,
            path: Path,
            hh: Histogram
    ) -> None:
        """Restore fitted models for the reported feature from a file.

         :param path: Full (dir + specific name) name of the input file.
         :param hh: Histogram object to store the models.
        """

        smr = Summary().read(path, cls.name)
        if smr is None:
            return

        smr = smr.data
        fc = [fit.class_from_classname(fit, s['model']['name'])
              if s['model'] is not None and s['model']['name'] is not None
              else None for s in smr]

        ee = [e.create() if hasattr(e, 'create') else e for e in fc]
        cc = [fit.subtype_from_classname(e, s['model']['name'])
              if e is not None else None for e, s in zip(ee, smr)]

        for c, d, h in zip(cc, smr, hh.elements()):
            if c is not None:
                cls.logger.info(f"{d['name']}: ")
                h.fits = [
                    fit.restore(c, d['model']['p'], cls.units, h.bc, h.h,
                                fit.method_from_classname(c,
                                                          d['model']['name']))
                ]

    @classmethod
    def log_stats(
            cls,
            path: Path
    ) -> Optional[Summary]:
        """Add a log record describing the restored attribute.

        :param path: Path to the input file of the input file.
        """

        smr = Summary().read(path, cls.name)
        if smr is None:
            return
        smr = smr.data
        for s in smr:
            if 'avg' in s:
                val = f"{s['avg']} ± {s['std']}" if 'std' in s else \
                      f"{s['avg']}"
                cls.logger.info(f"{cls.name} {s['name']}: {val} {cls.units}")

        return smr

    @classmethod
    def plot(
            cls,
            h: Histogram
    ):
        """Method prototype for creating subclass-specific plots.
        """

        pass
