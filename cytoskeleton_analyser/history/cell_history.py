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

"""
History of simulated microtubule state changes over the whole cell.
"""

import logging
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from .. import fitting as fit
from .. import visualization as vis
from ..inout import Paths
from ..inout import read_to_dtype
from ..inout import save_as_csv
from ..histograms import plot_single
from ..position.spatial_systems import FullDepth
from .event_collections import EventCollections
from .filament_history import FilamentHistory

# Turn interactive plotting off.
plt.ioff()


class CellHistory:
    """Encapsulates simulation-wide history data over all microtubules.

    ``CellHistory`` class is a container and analyser of dynamic
    instability-related characteristics of reconstructed microtubules.
    Data are recorded over an extended interval of simulation time for
    all microtubules present in the system, discriminated for
    microtubule plus- or minus-end. It holds an arrray of filament-
    specific histories (encapsilated by 'FilamentHistory' class)
    used for in-depth state type-specific analysis by a 'states.States'
    object.
    """

    def __init__(
            self,
            params: dict,
            rind: int,
            end: int,
            paths: Paths,
            logger: logging.Logger,
    ):
        """:param params: Main parameters.
           :param rind: Index of the simulation rum.
           :param end: Microtubule end: 0 or 1.
           :param paths: directory names for in/out data storage.
           :param logger: Logger
        """

        self.logger = logger
        fit.set_logger(self.logger)
        self.logger.info('')

        #: Common parameters.
        self.params = params
        if 'history_version' not in self.params:
            self.params['history_version'] = 2

        #: Format version of the recording file.
        self.version = self.params['history_version']

        #: Filament end index: 0 or 1.
        self.end = end

        #: Simulation run index.
        self.rind = rind

        #: Paths for data input/output.
        self.paths = paths

        #: Number of filament-specific history sequencies.
        self.fnum = np.uint64()

        #: Simulation system time.
        self.time = np.float64()
        if self.version == 1:
            #: Iteration index.
            self.iteration = np.uint64()
            #: Cell radius.
            self.cellrad = np.float32()
            #: Total mass of the microtubules (in polymerization units).
            self.mtmass = np.uint64()

        #: Container for the recorded data.
        self.records = self.read(paths.run /
                                 f'history_e{self.end}_{self.rind}.dat')
        # Time characteristics.
        self.recording_time = self.set_recording_time()
        #: Filament lengths.
        self.lengths = np.empty(len(self.records))

        #: 3d cell representation.
        self.cell = self.spatial_cell()

        #: Events and microtubule states categorised by type.
        self.collections = {}
        s = EventCollections.\
            initialise(self.params, self.end, self.paths.data_out, self.logger)
        for n in params['cell'].regions._fields:
            self.collections[n] = s(self.records, n)

    def read(
            self,
            fname: Path,
    ) -> list[FilamentHistory]:
        """Read the raw data into a list of (one per filament) records.

        Each one is a time-ordered sequence of kinetic characterisitics
        encapsulated by HistoryPerFilament class.

        :param fname: File containig raw data.
        :return: List of HistoryPerFilament instances.
        """

        self.logger.info('Importing history from ' + str(fname))

        with open(fname, 'rb') as f:
            if self.version == 1:
                self.iteration = read_to_dtype(self.iteration, f)
            self.time = read_to_dtype(self.time, f)
            if self.version == 1:
                self.cellrad = read_to_dtype(self.cellrad, f)
                self.mtmass = read_to_dtype(self.mtmass, f)
            self.fnum = read_to_dtype(self.fnum, f)

            records = []
            for _ in range(self.fnum):
                nrc = read_to_dtype(np.int64(), f)
                records.append(FilamentHistory(nrc, self.version).read(f))

        return records

    def report_microtubule_lengths(
            self
    ) -> dict:
        """Report average length of microtubules.

        Lengths are averaged over the history monitor period.
        :return: Dictoinary item 'MT_lengths'.
        """

        self.lengths = np.array([f.length[0] for f in self.records])
        a = self.params['edge_len'] * np.average(self.lengths)
        self.logger.info(f"Average MT length: "
                         f"{a} {self.cell.len_units}")

        return {'MT_lengths': a}

    def set_recording_time(
            self
    ) -> dict:
        """Set time boundaries and duration.

        Set time boundaries of the history monitor period and its
        duration as a dict.

        :return: dict items 'initial', 'final', 'duration'
        """

        res = {
            'initial': np.amin(np.array([np.amin(r.time)
                                         for r in self.records])),
            'final': np.amax(np.array([np.amax(r.time)
                                       for r in self.records])),
        }
        res['duration'] = res['final'] - res['initial']

        return res

    def report_duration(
            self,
    ) -> dict:
        """Make items of 'recording_time' dict available outside.
        """

        self.logger.info(f"Records are "
                         f"over {self.recording_time['duration']} sec "
                         f"(from {self.recording_time['initial']} "
                         f"to {self.recording_time['final']} sec)")

        return {'recording_time': self.recording_time}

    def report_state_frequencies(
            self,
            show: bool = True,
    ) -> dict:
        """Analysis of state types.

        Calculate and return absolute and relative frequencies of major
        state types as a dict. These includes microtubule end-specific
        frequency of catastrophes (shrinks), recoveries(growths) and
        pauses.

        :param show: If True, display the pie plot.
       """

        res = {}
        for kk, sts in self.collections.items():

            # Transition frequencies to each of the state types:
            dur = self.recording_time['duration']
            frequency = {
                'catastrophes': sts.shrinks.num() / dur,
                'recoveries':   sts.grows.num() / dur,
                'pauses':       sts.pauses.num() / dur,
            }
            for k, v in frequency.items():
                self.logger.info(f"Frequency of {kk} {k}: {v} 1/sec")

            f = list(frequency.values())
            if sum(f) > 0.:
                exportfile = self.paths.data_out / \
                             f'{kk}_frequencies_e{self.end}'
                vis.pie([a / sum(f) for a in f],
                        list(frequency.keys()),
                        title=f'{kk} frequencies for end {self.end}',
                        exportfile=exportfile,
                        show=show)

            # Fraction of all catastrophes (transitions to shrinking
            # state) that are not interrupted by pause:
            fraction_spontaneous_catastrophes = \
                sts.pure_s.gs.num() / \
                sts.shrinks.num() \
                if sts.shrinks.num() else np.nan
            self.logger.info(f"Fraction of spontaneous {kk} catastrophes: "
                             f"{fraction_spontaneous_catastrophes}")

            # Ratio of recovery to catastrophe events:
            ratio_rec2cat = \
                frequency['recoveries'] / \
                frequency['catastrophes'] \
                if frequency['catastrophes'] else np.nan
            self.logger.info(f"Ratio of {kk} recoveries to catastrophes: "
                             f"{ratio_rec2cat}")

            res['state_'+kk] = {'frequencies': frequency,
                                'fraction_spontaneous_catastrophes':
                                    fraction_spontaneous_catastrophes,
                                'ratio_rec2cat': ratio_rec2cat,
                                }
        return res

    def create_spatial_maps(
            self,
            show: bool = True,
    ) -> None:
        """Spatial maps of state intensities in cell xy projection.

        :param show: If True, display the generated maps.
        """

        self.collections['cytosol'].xy_shrink_growth_ratio(
            self.cell.discretize_xy(nbins=200)[0],
            show
        )

    def report_comet_radial_distribution(
            self,
            show: bool = True,
    ) -> None:
        """Distribution of growing microtubule ends.

        Positions of growing ends is discretised as a function of
        distance to cell center in xy plane.

        :param show: If True, display the histogrm plot.
        """

        fname = self.paths.data_out /\
                f'comet_radial_distribution_e{self.end}'

        edges, bincenters = self.cell.discretize_radius(nbins=300)
        h = self.collections['cytosol'] \
                .radial_distribution_of_growing_ends(edges)
        save_as_csv(fname.with_suffix('.csv'),
                    ['bincenters', 'freq'],
                    [bincenters, h])

        fig = plot_single(
            bincenters, h,
            fits=[],
            title=f'Growing ends for end {self.end}',
            xlabel=f"distance to center ({self.cell.len_units})",
            xlim=None)
        vis.save_plt(fname, gzipped=False)

        if show:
            plt.show()
        else:
            plt.close(fig)

    def spatial_cell(self) -> type[FullDepth]:
        """Create spatial representation of the cell.

        Only mesh representing plasma membrane is of interest for the
        history analysis here.
        """

        return FullDepth\
            .create(self.paths, self.params, self.rind, self.logger,
                    init=False)
