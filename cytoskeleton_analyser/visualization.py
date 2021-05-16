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

"""Customized plotting of project-specific data.
"""

from __future__ import annotations
import gzip
import itertools
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Sequence

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np

import cytoskeleton_analyser.fitting as fitting


# Turn interactive plotting off.
plt.ioff()  # to enable plotting to disk only


class Colors:
    """Customize colors from a matplotlib.pyplot colormap.
    """

    def __init__(
            self,
            cm: str,
    ):
        """:param cm: Name of matplotlib.pyplot colormap.
        """

        self.cmap = plt.get_cmap(cm)
        self.n = None
        self.strong: Optional[list] = None
        self.weak: Optional[list] = None
        self.all: Optional[list] = None

    def set_strong_weak(
            self,
            n: int,
    ) -> Colors:
        """Sets 'strong' and 'weak' class attributes from full color set.

        Sould be only used when the applied colormap discriminates
        colors into strong/weak pairs.
        """

        inds = np.arange(0, len(self.cmap.colors), 2)
        self.n = n
        _strong = self.cmap(inds)
        _weak = self.cmap(inds + 1)
        self.strong = [Colors.cycler(_strong, i) for i in range(self.n)]
        self.weak = [Colors.cycler(_weak, i) for i in range(self.n)]

        return self

    def set_all(
            self,
            n: int,
    ) -> Colors:
        """Sets 'n' and colors from the applied colorset.

        If desired number of colors exceeds that provided by the colorset,
        the missing colors are retrieved by cycling operation.

        :param n: Number of colors to be arranged.
        """

        inds = np.arange(0, len(self.cmap.colors))
        self.n = n
        _all = self.cmap(inds)
        self.all = [Colors.cycler(_all, i) for i in range(self.n)]

        return self

    @staticmethod
    def cycler(
            a: Sequence,
            i: int,
    ):
        """A minimalistic cycler over finite-length sequencies.
           :param a: Sequence of colors values in the colorset.
           :param i: color index.
        """

        return a[i % len(a)]


class Multiplot:
    """A multi-axis figure showing related data.

    Each axis displays a single data array along with an optional set of
    numerical fits as well as the expected saturation level
    if appropriate. All axes are equally scaled and ranged.
    """

    @dataclass
    class Data:
        figtitle: str
        axtitles: Optional[str]
        x: list[np.ndarray]
        x_label: str
        y: list[np.ndarray]
        y_label: str
        fit: list
        colormap: str

    def __init__(
            self,
            data: Data,
            exportfile: Optional[Path] = None,
            gzipped: bool = False,
            show: bool = True,
    ):
        """Produce the multiplot.

        No other method of this class needs to be called.
        """

        n = len(data.y)
        self.fig, ax = plt.subplots(
            ncols=n,
            sharey=True,
            figsize=(n*5., 5.),   # default: (6.4, 4.8)
        )
        if n == 1:
            ax = [ax]
        ax[0].set_ylabel(data.y_label)

        colors = Colors(data.colormap).set_strong_weak(n)

        for i in range(n):
            d = self._extract_axis_data(data, i)
            self._set_axis(ax[i], d, {'weak': colors.weak[i],
                                      'strong': colors.strong[i]})
        self.fig.suptitle(data.figtitle)
        #    lines_labels = [ax.get_legend_handles_labels() for ax in fig.axes]
        #    handles, labels = [sum(lol, []) for lol in zip(*lines_labels)]
        #    fig.legend(handles, labels, loc='upper center')
        self.fig.tight_layout()

        if exportfile is not None:
            save_plt(exportfile, gzipped)

        if show:
            plt.show()
        else:
            plt.close(self.fig)

    @staticmethod
    def _set_axis(
            ax: plt.Axes,
            d: Data,
            color: dict,
    ) -> None:
        """Set an instance of figure axes.

        :param ax: Axes to set.
        :param d: Data container.
        :param color: dict containing weak and strong color tints.
        """

        ls = itertools.cycle(['--', '-.', ':'])
        ax.set_xlabel(d.x_label)
        if d.axtitles is not None:
            ax.set_title(d.axtitles)

        cr = color['weak'] if d.fit is not None and len(d.fit) else \
            color['strong']
        ax.scatter(d.x, d.y, color=cr, s=0.2, marker='.', label='_nolegend_')

        if d.fit is not None:
            for f in d.fit:
                if f is None or fitting.is_mock(f):
                    continue
                if f.isbest:
                    ax.plot(d.x, f.prediction, c=color['strong'],
                            lw=1., label='_nolegend_')
                    ax.plot([d.x[0], d.x[-1]], [f.saturation, f.saturation],
                            c=[0.3]*3, lw=0.25, label='_nolegend_')
                else:
                    ax.plot(d.x, f.prediction, c='k', lw=0.5,
                            ls=next(ls), label='_nolegend_')
        # ax.grid(True)

    @staticmethod
    def _extract_axis_data(
            d: Data,
            i: int,
    ) -> Data:
        """Extract a Data object corresponding to index 'i'.

        :param d: Data object to extract from.
        :param i: Index of the data array.
        :return: Data object containing single data for plotting.
        """

        return Multiplot.Data(
            figtitle=d.figtitle,
            axtitles=d.axtitles[i]
            if d.axtitles is not None and len(d.axtitles) else None,
            x=d.x,
            x_label=d.x_label,
            y=d.y[i],
            y_label=d.y_label,
            fit=d.fit[i] if d.fit is not None and len(d.fit) else None,
            colormap=d.colormap,
        )


def pie(
        sizes: list,
        labels: list[str],
        title: Optional[str] = None,
        explode: Optional[Sequence] = None,
        exportfile: Optional[Path] = None,
        gzipped: bool = False,
        show: bool = True,
) -> None:
    """Produce, show and optionally save a pie plot.

    :param sizes: Data to plot.
    :param labels: Section labels.
    :param title: Figure title.
    :param explode: Sections to highlite.
    :param exportfile: If not None, the plot will be saved here
       in .svg format.
    :param gzipped: If True, apply .gz compression on top of
        the .svg file.
    :param show: If True, display the plots.
    """

    fig, ax = plt.subplots()
    ax.pie(sizes,
           explode=explode,
           labels=labels,
           autopct='%1.1f%%',
           startangle=90)
    ax.axis('equal')  # ensures that pie is drawn as a circle.
    if title is not None:
        ax.set_title(title)

    if exportfile is not None:
        save_plt(exportfile, gzipped)

    if show:
        plt.show()
    else:
        plt.close(fig)


def colormesh_plot(
        h: np.ndarray,
        title: str,
        cmap: str,
        norm: mcolors.Normalize,
        axis_off: bool = True,
        exportfile: Optional[Path] = None,
        show: bool = True,
) -> None:
    """Produce, show and optionally save a colormesh plot.

    :param h: Data to plot.
    :param title: Figure title.
    :param cmap: Colormap.
    :param norm: Normalization.
    :param axis_off: If true, hide axis.
    :param exportfile: If not None, the plot will be saved here
        in .tiff format.
    :param show: If True, display the plots.
    """

    fig, (ax) = plt.subplots(1, 1)
    im = ax.pcolormesh(h, cmap=cmap, norm=norm)
    fig.colorbar(im, ax=ax)
    ax.set_aspect('equal')
    ax.set_title(title)
    if axis_off:
        ax.set_axis_off()
    if exportfile is not None:
        save_plt(exportfile, frmt='tiff', gzipped=False)
    if show:
        plt.show()
    else:
        plt.close(fig)


def save_plt(
        exportfile: Path,
        gzipped: bool = True,
        frmt: str = 'svg',
) -> None:
    """Save a matplotlib.pyplot plot as an image file.

    :param exportfile: Full path including the file
        but without extension.
    :param gzipped: If True, apply .gz compression on top of the
        image format extension.
    :param frmt: Image file format: default: '.svg'.
    """

    assert len(str(exportfile).split()) == 1

#    exportfile = underscored(exportfile)
    kw = {"compression": "tiff_lzw"} if frmt == 'tiff' else {}
    if gzipped:
        with gzip.open(exportfile.with_suffix(f'.{frmt}.gz'), 'wb') as f:
            plt.savefig(f, format=frmt)
    else:
        with open(exportfile.with_suffix(f'.{frmt}'), 'wb') as f:
            plt.savefig(f, format=frmt, pil_kwargs=kw)
