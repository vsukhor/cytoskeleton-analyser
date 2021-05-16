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

"""Project-specific data histograms for visualisation and storage.
"""

from __future__ import annotations
import csv
import itertools
import math
import gzip
from logging import Logger
from pathlib import Path
from typing import Optional, Sequence

import matplotlib.pyplot as plt
import numpy as np

import cytoskeleton_analyser.visualization as vis
import cytoskeleton_analyser.fitting as fts
from .inout import underscored


class Experimental:
    """Histogram derived from experimental data.

    Optioanlly include fits to the frequency array.
    """

    def __init__(self):

        self.bc = np.empty(0, dtype=float)  #: Bin centers.
        self.h = np.empty(0, dtype=float)   #: Frequency values.
        self.fits = []                      #: Fits to 'h' over 'bc'.
        self.extent = []                    #: Plot x axis extent.
        self.avg = None                     #: Average of 'h' over 'bc'.

    def initialise(
            self,
            raw: tuple,
            fits: list
    ):
        """Create the histogram along with the accompanying data fields.
        """

        self.bc = np.array(raw[0])  # bin centers
        self.h = np.array(raw[1])
        self.extent = extent(self.bc)
        self.fits = fts.fit(fits, self.bc, self.h)
        self.avg = avg(self.bc, self.h)

        return self

    @classmethod
    def restore(
            cls,
            bc: Sequence,
            h: Sequence
    ) -> Experimental:
        """Restore the class from bin centers and frequency arrays.

        :param bc: Bin centers.
        :param h: Frequencies.
        :return: This class.
        """

        c = cls()
        c.bc = np.array(bc)
        c.h = np.array(h)
        c.extent = extent(bc)
        c.avg = avg(c.bc, c.h)

        return c


class Simulated:
    """Histogram derived from simulated data.

    Optioanlly include fits to the frequency array.
    """

    def __init__(self):

        self.bc: np.ndarray = np.empty(0, dtype=float) #: Bin centers.
        self.h: np.ndarray = np.empty(0, dtype=float)  #: Frequency values.
        self.fits: list = []           #: Fits to 'h' over 'bc'.
        self.extent: list = []         #: x axis extent for plotting.
        self.avg = None                #: Average of 'h' over 'bc'.
        self.polar: bool = False       #: Plot as a polar graph if True.
        self.halfpolar: bool = False   #: Use 0:pi domain of polar plots

    def initialise(
            self,
            data: Sequence,
            fits: list,
            dx: float,
            density: bool = False,
            exper_bc: Optional[np.ndarray] = None,
            polar: bool = False,
            halfpolar: bool = False,
    ):
        """Create the histogram along with the accompanying data fields.
        """

        data = np.concatenate(data)
        self.extent, dx = self._extent(data, dx, exper_bc, polar)
        bins = round((self.extent[1] - self.extent[0]) / dx)
        self.h, bin_edges = np.histogram(data, bins=bins,
                                         range=self.extent, density=density)
        self.bc = bin_edges[:-1] + dx / 2 \
            if exper_bc is None \
            else exper_bc
        self.fits = fts.fit(fits, self.bc, self.h)
        self.avg = avg(self.bc, self.h)
        self.polar = polar
        self.halfpolar = halfpolar

        return self

    @classmethod
    def restore(
            cls,
            bc: Sequence,
            h: Sequence,
    ) -> Simulated:
        """Restore the class from bin centers and frequency arrays.

        :param bc: Bin centers.
        :param h: Frequencies.
        :return: This class.
        """

        c = cls()
        c.bc = np.array(bc)
        c.h = np.array(h)
        c.avg = avg(c.bc, c.h)
        c.extent = extent(c.bc)

        return c

    @staticmethod
    def _extent(
            data: Sequence,
            dx: float,
            exper_bc: Optional[np.ndarray] = None,
            polar: bool = False
    ):
        """Auxiliary method.
        """

        if exper_bc is None:
            # Set the bins to correspond to the simulated histogram.
            limits = [math.floor(min(data) / dx) * dx,
                      math.ceil(max(data) / dx) * dx]
            if not polar:
                limits[0] -= dx / 2.
                limits[1] += dx / 2.
            return limits, dx

        else:
            dx = exper_bc[1] - exper_bc[0]  # discard dx given in args
            return extent(exper_bc), dx


class Histogram:
    """Container for histograms.

     The histograms may be derived from simulated and/or experimental
     data.
    """

    logger: Logger

    def __init__(
            self,
            name: str,
            simulated: Optional[Simulated] = None,
            experimental: Optional[Experimental] = None,
    ):
        """
        :param name: Name of the cytoskeleton attribute.
        :param simulated: Histogram derived from simulated data.
        :param experimental: Histogram derived from experimental data.
        """

        #: Name of the cytoskeleton attribute.
        self.name = name

        #: Histogram object derived from simulated data.
        self.simulated = simulated

        #: Histogram object derived from experimental data.
        self.experimental = experimental

    def plot(
            self,
            title: str,
            xlabel: str,
            xlim: Optional[list] = None,
            yscale: str = 'linear',
            save_path: Optional[Path] = None,
            gzipped: bool = False,
            show: bool = True,
    ) -> None:
        """Plot the child histograms.

        Either single or mirrored double plot are produced
        depending on whether the experimental data is available.
        Optionally, save the figure in .svg file format.

        :param title: Plot title.
        :param xlabel: X-axis label.
        :param xlim: X-axis limits.
        :param yscale: Scaling: 'linear' or 'log'.
        :param save_path: If not None, the plot will be saved here
            in .svg format.
        :param gzipped: If True, apply .gz compression on top of
            the .svg file.
        :param show: If True, display the histogram plot.
        """

        if self.experimental is not None:
            fig = plot_double(
                x=self.simulated.bc,
                h=[self.experimental.h, self.simulated.h],
                fits=[self.experimental.fits, self.simulated.fits],
                title=title,
                yscale=yscale)
        else:
            fig = plot_single(
                x=self.simulated.bc,
                h=self.simulated.h,
                fits=self.simulated.fits,
                title=title,
                xlabel=xlabel,
                xlim=xlim,
                yscale=yscale,
                polar=self.simulated.polar,
                halfpolar=self.simulated.halfpolar)

        if save_path is not None:
            vis.save_plt(save_path / underscored(self.name), gzipped)
        if show:
            plt.show()
        else:
            plt.close(fig)

    def elements(self):
        """Return whatever is available.
        """

        return self.simulated, self.experimental \
            if self.experimental is not None else \
            self.simulated,

    def to_csv(
            self,
            path: Path,
            gzipped: bool = False,
    ) -> None:
        """Save the histogram to a .csv data table file.

        :param path: Directory to place the resulting file.
        :param gzipped: If True, apply .gz compression.
        """

        def row(j):
            r = [float(self.simulated.bc[j]),
                 float(self.simulated.h[j])]
            return r if self.experimental is None else \
                r + [float(self.experimental.h[j])]

        fname = (path / underscored(self.name)).with_suffix('.csv')
        if gzipped:
            fname = fname.with_suffix('.gz')
            opn = gzip.open
        else:
            opn = open
        with opn(fname, 'wt', newline='') as f:
            self.logger.info("exporting to: " + str(fname))
            w = csv.writer(f, delimiter=',')
            caption = ['x', 'sim']
            if self.experimental is not None:
                caption += ['emp']
            w.writerow(caption)
            n = len(self.simulated.h)
            for i in range(n):
                w.writerow(row(i))
        self.logger.info(f"exported: {n} records\n")

    def from_csv(
            self,
            path: Path,
            gzipped: bool = False
    ) -> Histogram:
        """Import the histogram from a .csv data table file.

        :param path: Directory read the resulting file from.
        :param gzipped: If True, the file is .gz-compressed.
        """

        fname = (path / underscored(self.name)).with_suffix('.csv')
        if gzipped:
            fname = fname.with_suffix('.gz')
            opn = gzip.open
        else:
            opn = open
        self.logger.info("importing from: " + str(fname))
        with opn(fname, 'rt', newline='') as f:
            r = csv.reader(f, delimiter=',')
            cap = r.__next__()  # caption
            with_emp = len(cap) > 2

            bc, hs, he = [], [], []
            for a in r:
                bc.append(float(a[0]))
                hs.append(float(a[1]))
                if with_emp:
                    he.append(float(a[2]))

            self.simulated = Simulated.restore(bc, hs)
            if with_emp:
                self.experimental = Experimental.restore(bc, he)

        self.logger.info(f"imported: {len(self.simulated.bc)} records\n")

        return self


def plot_double(
        x: Sequence,
        h: list[np.ndarray],
        fits: Sequence,
        title: str,
        yscale: str,
) -> plt.Figure:
    """A plot with dual mirrored y-axes and common vertical x-axis.

    Is convenient for visual comparison of two histograms,
    e.g. that of empirical and simulated data.

    :param x: Bin centers.
    :param h: Frequencies.
    :param fits: Fits to the frequencies.
    :param title: Figure title.
    :param yscale: Scaling: 'linear' or 'log'.
    :return: The figure created here.
    """

    hmax = max(h[0].max(), h[1].max())
    lmax = max(max([[f.prediction.max() for f in ff] if len(ff) else [0.]
                    for ff in fits]))
    maxval = max(hmax, lmax)
    clrs = [vis.Colors('Set1').set_all(len(ll)) for ll in fits]
    fig, ax = plt.subplots(1, 2)
    for a, hh in zip(ax, h):
        a.barh(y=x, width=hh, height=x[1] - x[0])
    for a, c, ff in zip(ax, clrs, fits):
        for i, f in enumerate(ff):
            a.plot(f.prediction, x, c=c.all[i], lw=1., label=f.name)
    # Add some text for labels, title and custom x-axis tick labels, etc.
    #  ax[0].set_xlabel('frequency')
    ax[0].set_xlim([maxval, 0])
    ax[1].set_xlim([0, maxval])
    ax[0].set_yscale(yscale)
    ax[1].set_yscale(yscale)
    ax[0].spines['left'].set_visible(False)
    ax[0].spines['top'].set_visible(False)
    ax[0].spines['bottom'].set_visible(False)
    ax[1].spines['right'].set_visible(False)
    ax[1].spines['top'].set_visible(False)
    ax[1].spines['bottom'].set_visible(False)
    ax[0].yaxis.set_ticks_position('right')
    [a.set_xticks([]) for a in ax]
    ax[1].set_yticklabels([])
    ax[0].set_title('experimental')
    ax[1].set_title('simulated')
    #    ax.legend()
    fig.suptitle(title)
    fig.text(
        x=0.5,
        y=0.85,
        s='(1/μm)',
        ha='center',
        va='bottom',
    )
    for a, f in zip(ax, fits):
        if len(f):
            a.legend()
    #    fig.legend()
    fig.tight_layout()

    return fig


def plot_single(
        x: Sequence,
        h: Sequence,
        fits: list,
        title: str,
        xlabel: str,
        xlim: Optional[list] = None,
        yscale: str = 'linear',
        polar: bool = False,
        halfpolar: bool = False,
) -> plt.Figure:
    """A project-wide standard for plotting histograms.

    Optionally include a set of fitted functions.

    :param x: Bin centers.
    :param h: Frequencies.
    :param fits: Fits to the frequencies.
    :param title: Figure title.
    :param xlabel: X label.
    :param xlim: Limits on x axis.
    :param yscale: Scaling: 'linear' or 'log'.
    :param polar: Plot as a polar graph if True.
    :param halfpolar: Use only the first half of the polar domain.
    :return: The figure created here.
    """

    subplot_kw = {'projection': 'polar'} if polar else None
    fig, ax = plt.subplots(subplot_kw=subplot_kw)
    ax.bar(x=x, height=h, width=abs(x[1] - x[0]))
    if len(fits):
        c = vis.Colors('Set1').set_all(len(fits))
        ls = itertools.cycle(['--', '-.', ':'])
        [ax.plot(x, f.prediction, c=c.all[i], lw=1., label=f.name) and
         [ax.plot(x, cc, c=c.all[i], ls=next(ls), lw=1., label=f"  comp {j}")
          for j, cc in enumerate(f.cc)]
         for i, f in enumerate(fits)]
    # Add some text for labels, title and custom x-axis tick labels, etc.
    if polar:
        if halfpolar:
            ax.set_thetamin(0)
            ax.set_thetamax(180)
    else:
        ax.set_ylabel('frequency')
        ax.set_xlabel(xlabel)
        if xlim is not None:
            ax.set_xlim(left=xlim[0], right=xlim[1])
        ax.set_yscale(yscale)
    ax.set_title(title)
    if len(fits):
        ax.legend()
    #    fig.title(title)
    plt.tight_layout()

    return fig


def avg(
        bc: np.ndarray,
        h: np.ndarray,
) -> float:
    """Average value of a frequency sequence data.
    """

    return np.mean(h * bc) / (np.sum(h) * (bc[1] - bc[0]))


def extent(
        bc: Sequence[np.number]
) -> list[float]:
    """Axes extent form a data array.
    """

    dx = bc[1] - bc[0]
    if dx > 0:
        return [bc[0] - dx / 2,
                bc[-1] + dx / 2]
    else:
        return [bc[-1] + dx / 2,
                bc[0] - dx / 2]
