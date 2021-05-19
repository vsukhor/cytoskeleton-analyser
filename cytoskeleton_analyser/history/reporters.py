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

"""Classes for analysis and visualization of microtubule dynamic.

Extracted parameters of dynamic instability are based on data from
'state_sequences' module.
"""

from __future__ import annotations
from typing import Any

import numpy as np

import cytoskeleton_analyser.fitting as fit
from ..histograms import Histogram
from ..histograms import Simulated
from ..history.state_sequences import StateSequence
from ..report import Report
from ..summary import Summary


class _StateSequenceReport(Report):
    """Report extended to make StateSequence attributs more convenient.

    Includes plotting and saving data histogram and optional fits.
    """

    #: State sequence to report about.
    ph: StateSequence

    #: Human-readable designation of the reported feature.
    feature: str

    #: Attribute.
    attr: Any

    @classmethod
    def create(
            cls,
            ph: StateSequence,
            region: str,
    ) -> type[_StateSequenceReport]:
        """Initialise class attributes from record sequences.

        :param ph: Input record sequences.
        :param region: Cell Region.
        :return: This class.
        """

        cls.ph = ph
        cls.logger = ph.logger
        fit.base.set_logger(ph.logger)
        Histogram.logger = ph.logger
        Summary.logger = ph.logger
        cls.path_out = ph.path_out
        cls.attr = getattr(cls.ph, cls.feature)
        cls.name = region + ' ' + ph.name() + ' ' + \
            cls.feature + f' at end {ph.END}'

        return cls

    @classmethod
    def is_valid(cls) -> bool:
        """Attribute valid iff it is number other than zero.
        """

        return \
            not np.isnan(cls.attr.avg) and \
            not np.isclose(cls.attr.avg, 0.)

    @classmethod
    def plot_sequence_property(
            cls,
            h: Histogram,
            xlim: list,
            yscale: str = 'linear',
            show: bool = True,
    ) -> None:
        """State sequence-specific plots of data histograms.

        :param h: Histogram object to polot.
        :param xlim: Limits on x axis.
        :param yscale: Scaling: 'linear' or 'log'.
        :param show: If True, display the plot.
        """

        h.plot(
            cls.name,
            xlabel=f'{cls.feature} ({cls.units})',
            xlim=xlim,
            yscale=yscale,
            save_path=cls.path_out,
            show=show,
        )

    @classmethod
    def pipeline(
            cls,
            show: bool = True
    ) -> dict:
        """Prototype for microtubule attribute-specific data pipelines.
        """

        pass


class Elongation(_StateSequenceReport):
    """Class reporting microtubule elongation during the state
    sequence lifetime.
    """

    @classmethod
    def create(
            cls,
            ph: StateSequence,
            region: str,
    ) -> type[_StateSequenceReport]:
        """Initialise class attributes from record sequences.

        :param ph: Input record sequences.
        :param region: Cell Region.
        :return: This class.
        """

        cls.feature = __class__.__name__.lower()
        return super().create(ph, region)

    @classmethod
    def report(
            cls,
            show: bool = True
    ):

        res = cls.ph.report(cls.feature)
        if cls.is_valid():
            cls.units = res[cls.feature]['units']
            data = [np.abs(cls.ph.dlen_um)]

            h = Histogram(
                cls.name,
                Simulated()
                .initialise(data, cls.fits_sim, dx=0.12, density=True)
            )
            h.to_csv(cls.path_out)
            cls.plot_sequence_property(h, xlim=[0., 40.], show=show)
            cls.logger.info('')
            return res, h

        return res, None

    @classmethod
    def pipeline(
            cls,
            show: bool = True
    ) -> dict:

        e = fit.Exponential.create()
        p = fit.Exponential.Pars
        cls.fits_sim = [
            (e.d_h, p(a=np.max(np.abs(cls.ph.dlen_um)),
                      tau1=np.abs(cls.ph.elongation.avg)),
             cls.ph.elongation.units),
            # (fit.Gamma.loc0, [3., 0.2, 0.]),
        ]

        res, h = cls.report(show)
        if cls.is_valid():
            cls.summarize(
                (h, [res[cls.feature]['avg']], [res[cls.feature]['std']]),
                [0]
            )

        return res


class Velocity(_StateSequenceReport):
    """Velocity report class.
    """

    @classmethod
    def create(
            cls,
            ph: StateSequence,
            region: str,
    ):

        cls.feature = __class__.__name__.lower()
        return super().create(ph, region)

    @classmethod
    def report(
            cls,
            show: bool = True
    ):

        res = cls.ph.report(cls.feature)
        cls.units = res[cls.feature]['units']
        if cls.is_valid():
            data = [np.abs(cls.ph.vel[np.abs(cls.ph.vel) < 200.])]

            h = Histogram(
                cls.name,
                Simulated()
                .initialise(data, cls.fits_sim, dx=0.5, density=True)
            )
            h.to_csv(cls.path_out)
            cls.plot_sequence_property(h, xlim=[0., 50.], show=show)
            cls.logger.info('')
            return res, h
        return res, None

    @classmethod
    def pipeline(
            cls,
            show: bool = True
    ) -> dict:

        cls.fits_sim = [
        ]

        res, h = cls.report(show)
        if cls.is_valid():
            cls.summarize(
                (h, [res[cls.feature]['avg']], [res[cls.feature]['std']]),
                [0]
            )

        return res
