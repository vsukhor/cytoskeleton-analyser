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

"""Holds abstract base class for sequences of filament dynamic states.
"""

from __future__ import annotations
import logging
from abc import ABC
from pathlib import Path
from typing import Final

import numpy as np

from ...histograms import Histogram
from ...histograms import Simulated
from ...history.filament_history import FilamentHistory
from ...report import Stats


class _StateSequence(ABC):
    """Abstract base class for a family of classes holding sequences of
    filament end dynamic states.

    Provides storage of state records and their basic analysis.
    More extended analysis could come from 'reporters' module.
    Should be only used for length-specific subclassing.
    The subclasses represent sequences various lengths and types.
    """

    PARAMS: dict            #: General-purpose parameters.
    logger: logging.Logger  #: Current logger.
    END: int                #: Index of microtubule end: 0 or 1.
    path_out: Path          #: File system path for data output.
    WIDTH: int              #: Sequence width in the event records.

    #: Slice of vars(self).keys() for automatic update.
    DATA_INDS: slice

    def __init__(
            self,
    ):

        self.fi = np.empty(0, dtype=int)
        self.inds = np.empty(0, dtype=int)
        if self.PARAMS['history_version'] == 2:
            self.ngrw = np.empty(0, dtype=np.uint64)
            self.nshr = np.empty(0, dtype=np.uint64)
            self.dngrw = np.empty(0, dtype=np.int32)
            self.dnshr = np.empty(0, dtype=np.int32)
            self.dornt = np.empty(0, dtype=np.float32)
        self.dlen = np.empty(0, dtype=np.int32)
        self.dlen_um = np.empty(0, dtype=np.float32)
        self.dtime = np.empty(0, dtype=np.float64)
        self.vel = np.empty(0, dtype=np.float32)
        self.length = np.empty(0, dtype=np.uint32)
        self.time = np.empty(0, dtype=np.float64)
        self.pos_fr = np.empty([0, 3], dtype=np.float32)
        self.pos_to = np.empty([0, 3], dtype=np.float32)
        self.dist0_fr = np.empty(0, dtype=np.float32)
        self.dist0_to = np.empty(0, dtype=np.float32)
        self.dist_plm = np.empty(0, dtype=np.float32)
        self.dist_nuc = np.empty(0, dtype=np.float32)
        self.cas = np.empty([0, 6], dtype=np.float32)

        self._duration = Stats(np.nan, np.nan, '')
        self._elongation = Stats(np.nan, np.nan, '')
        self._reorientation = Stats(np.nan, np.nan, '')
        self._velocity = Stats(np.nan, np.nan, '')

    def append(
            self,
            r: FilamentHistory,
            ii: np.ndarray,
            j: int,
    ) -> None:
        """Append new events to the current data set.

        :param r: A track of filament-specific history.
        :param ii: Indexes in 'r' appropriate for the current
            state sequence.
        :param j: Index of the filament-specific history record the new
            events are taken from.
        """

        s: Final = self.WIDTH
        self.inds = np.hstack((self.inds, ii))

        self.fi = np.hstack((self.fi, np.array([j] * len(ii), dtype=int)))

        if self.PARAMS['history_version'] == 2:
            #: Number of growth events.
            self.ngrw = np.hstack((self.ngrw, r.ngrw[ii]))
            #: Number of shrink events.
            self.nshr = np.hstack((self.nshr, r.nshr[ii]))

            dngrw = np.int32(r.ngrw[ii + s] - r.ngrw[ii])
            self.dngrw = np.hstack((self.dngrw, dngrw))
            dnshr = np.int32(r.nshr[ii + s] - r.nshr[ii])
            self.dnshr = np.hstack((self.dnshr, dnshr))
            dlen = dngrw - dnshr

            # 3d angle between two unit vectors:
            dornt = np.arccos(
                np.clip([np.dot(r.ornt[i + s, :], r.ornt[i, :]) for i in ii],
                        -1.0, 1.0)
            ) / np.pi * 180.
            #: 3d angle between two unit vectors:
            self.dornt = np.hstack((self.dornt, dornt))
        else:
            dlen = np.int32(r.length[ii + s]) - np.int32(r.length[ii])
        #
        self.dlen = np.hstack((self.dlen, dlen))

        dlen_um = dlen * self.PARAMS['edge_len']
        self.dlen_um = np.hstack((self.dlen_um, dlen_um))

        dtime = r.time[ii + s] - r.time[ii]
        self.dtime = np.hstack((self.dtime, dtime))

        vel = np.divide(dlen_um, dtime,
                        out=np.zeros_like(dlen_um),
                        where=dtime != 0.) * 60.
        self.vel = np.hstack((self.vel, vel))

        self.length = np.hstack((self.length, r.length[ii]))
        self.time = np.hstack((self.time, r.time[ii]))
        self.pos_fr = np.vstack((self.pos_fr, r.pos[ii, :]))
        self.pos_to = np.vstack((self.pos_to, r.pos[ii + s, :]))
        self.dist0_fr = np.hstack((self.dist0_fr, r.dist0[ii]))
        self.dist0_to = np.hstack((self.dist0_to, r.dist0[ii + s]))
        self.dist_plm = np.hstack((self.dist_plm, r.dist_plm[ii]))
        self.dist_nuc = np.hstack((self.dist_nuc, r.dist_nuc[ii]))
        self.cas = np.vstack((self.cas, r.cas[ii, :]))

    def remove_instantaneous(self) -> None:
        """Remove states (if present) with duration comparable to
        the machine precision to reduce noise.
        """

        inds_not_null = np.where(
            ~np.isclose(self.dtime, np.zeros_like(self.dtime), atol=0.)
        )[0]
        for k in list(vars(self).keys())[self.DATA_INDS]:
            setattr(self, k, getattr(self, k)[inds_not_null])

    def stack_attributes(
            self,
            other,
            new,
    ):
        """Concatenate elements of attribute data sets.

        This instance is stacked with ``other`` into a ``new`` one.
        """

        for k in list(vars(self).keys())[self.DATA_INDS]:
            this = getattr(self, k)
            that = getattr(other, k)
            if this.ndim == 1:
                setattr(new, k, np.hstack((this, that)))
            else:
                setattr(new, k, np.vstack((this, that)))

        return new

    def num(self) -> int:
        """Number of simulation events.

        The events are from the state sequence trail.
        """

        return len(self.time)

    def __or__(
            self,
            other: _StateSequence,
    ):
        """Abstract 'or' operator.

        To be overwritten in subclasses.
        """
        pass

    def name(self) -> str:
        """Abstract method for the sequence name.

        To be overwritten in subclasses.
        """

        return ''

    def xy_distribution(
            self,
            edges: np.ndarray,
    ) -> np.ndarray:
        """2d histogram of event positions projected onto xy plane.

        :param edges: Nx2 array of N bin edges for space discretization
            in x- and y-directions respectively.
        :return: two-dim histogram of event positions in xy plane.
        """

        return np.histogram2d(
            self.pos_fr[:, 0],
            self.pos_fr[:, 1],
            bins=[edges[0, :], edges[1, :]]
        )[0]

    def report(
            self,
            attr: str,
    ) -> dict:
        """A very simple report.

        More extended versions are given by separate classes derived
        from 'reporters._StateSequenceReport'.

        :param attr: Class property to report about.
        :return: Average, standard deviation and measurement units.
        """

        assert self.logger is not None

        v = getattr(self, attr)
        self.logger.info(self.name() + f" {attr}: "
                                       f"{v.avg} ± "
                                       f"{v.std} {v.units}")
        return {attr: {'avg': v.avg,
                       'std': v.std,
                       'units': v.units}}

    # Major kinetic characteristics of the state sequence
    # over its lifetime.

    # Reorientation of the the microtubule end.

    @property
    def reorientation(self) -> Stats:

        if np.isnan(self._reorientation.avg):
            setattr(self, 'reorientation', None)
        return self._reorientation

    @reorientation.setter
    def reorientation(self, _) -> None:

        if self.dornt.size:
            self._reorientation = Stats(
                np.float64(self.dornt.mean()),
                np.float64(self.dornt.std()),
                'grad',
            )

    def report_reorientation(
            self,
            region: str,
            show: bool = True
    ) -> dict:
        """Extended version of reorientation report.

        Includes plotting histogram with eventual fits.

        :param region: Cell region.
        :param show: If True, display the histogram plot.
        :return: Summary of the reorientation stats.
        """

        res = self.report('reorientation')
        if self.dlen_um.shape[0]:
            name = region + ' ' + self.name() + \
                   f" reorientation at end {self.END}"
            fits = [
                # (fit.Gamma.loc0, [3., 0.1, 0.]),
            ]
            h = Histogram(
                name,
                Simulated().initialise([np.abs(self.dornt)],
                                       fits, dx=1., density=True)
            )
            h.to_csv(self.path_out)

            h.plot(
                name,
                xlabel='angle (grad)',
                xlim=[0., 180.],
                save_path=self.path_out,
                show=show,
            )

        return res

    # State sequence duration.

    @property
    def duration(self) -> Stats:

        if np.isnan(self._duration.avg):
            setattr(self, 'duration', None)
        return self._duration

    @duration.setter
    def duration(self, _) -> None:

        if self.dtime.size:
            self._duration = Stats(self.dtime.mean(),
                                   self.dtime.std(),
                                   'sec')

    # Microtubule length change during the state sequence duration.

    @property
    def elongation(self) -> Stats:

        if np.isnan(self._elongation.avg):
            setattr(self, 'elongation', None)
        return self._elongation

    @elongation.setter
    def elongation(self, _) -> None:

        if self.dlen_um.size:
            self._elongation = Stats(self.dlen_um.mean(),
                                     self.dlen_um.std(),
                                     'μm')

    # Average velocity of the mictrotubule end (de)polymerisation.

    @property
    def velocity(self) -> Stats:

        if np.isnan(self._velocity.avg):
            setattr(self, 'velocity', None)
        return self._velocity

    @velocity.setter
    def velocity(self, _) -> None:

        if self.vel.size:
            self._velocity = Stats(self.vel.mean(),
                                   self.vel.std(),
                                   'μm/min')
