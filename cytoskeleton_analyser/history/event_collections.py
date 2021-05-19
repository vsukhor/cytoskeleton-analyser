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

"""Collections of kinetic states the microtubule end can adopt.

   Formalized as EventCollections class.
   Microtubule ends are characterized by a kinetic state which persists
   until an instantenious transition to a new state occurs. Hence, the
   transition event is characterized by a pair of state types before and
   after the transition, all together forming a discrete Markov chain.
"""

from __future__ import annotations
import logging
from collections import namedtuple
from pathlib import Path
from typing import Final, Sequence

import meshio
import numpy as np

import cytoskeleton_analyser.visualization as vis
from ..histograms import Histogram
from ..history.filament_history import FilamentHistory
from ..history.state_sequences import Single
from ..history.state_sequences import Double
from ..history.state_sequences import Triple
from ..report import Stats
from .reporters import Elongation
from .reporters import Velocity


class EventCollections:
    """Encapsulates kinetic states a microtubule end can adopt.

    Extracts these from transition records stored in HistoryPerFilament
    andclassifies according to the state that the transitions generate.
    Includes methods for reporting and saving state characteristics.
    """

    PARAMS: dict            #: Dictionary of history-wide parameters.
    END: int                #: Filament end: 0 or 1.
    path_out: Path          #: Path to data output.
    logger: logging.Logger

    #: Collection defining single-state sequence is a pair of
    #: consecutive states. The state notation is a low-case version of
    #: state_sequences.
    PureS = namedtuple('PureS', 'sg pg gs ps gp sp cs cg gd sd pd')

    #: Collection defining double-state sequence is a triples of
    #: consecutive states. The state notation is a low-case version of
    #: state_sequences.
    PureD = namedtuple('PureD', 'gsg psg csg gps cgs cgp sgs pgs sgp')

    #: Collection defining triple-state sequences is a
    #: quadruples of consecutive states. The state notation is a
    #: low-case version of state_sequences.
    PureT = namedtuple('PureT', 'sgps cgps')

    def __init__(
            self,
            records: Sequence[FilamentHistory],
            cell_region: str,
    ):
        """
        :param records: List of HistoryPerFilament records.
        :param cell_region: Only records from this cell region are in.
        """

        self.include_ornt: bool = self.PARAMS['history_version'] > 1

        self.region: str = cell_region
        cr = getattr(self.PARAMS['cell'].regions, self.region)
        lim = [cr.rad_min, cr.rad_max]
        self.pure_s = self.PureS(*[Single(f).set_from_records(records, lim)
                                   for f in self.PureS._fields])
        self.pure_d = self.PureD(*[Double(f).set_from_records(records, lim)
                                   for f in self.PureD._fields])
        self.pure_t = self.PureT(*[Triple(f).set_from_records(records, lim)
                                   for f in self.PureT._fields])

        # Microtubule ends adopt one of three states: a shrinking,
        # growing or paused. Each state can arise from multitple
        # transition types.

        # Tracks of single states:
        t1 = self.pure_s
        self.shrinks = t1.gs | t1.ps | t1.cs  #: Shrink state.
        self.grows = t1.sg | t1.pg | t1.cg    #: Growth state.
        self.pauses = t1.gp | t1.sp           #: Pause state.

        # Combinations of two consecutive states:
        t2 = self.pure_d
        #: Shrink followed by growth
        self.sg = t2.gsg | t2.psg | t2.csg
        #: Growth followed by shrink
        self.gs = t2.sgs | t2.pgs | t2.cgs

        # Combinations of three consecutive states:
        t3 = self.pure_t
        #: Growth followed by pause then followed by shrink.
        self.gps = t3.sgps | t3.cgps

    @classmethod
    def initialise(
            cls,
            params: dict,
            end: int,
            path_out: Path,
            logger: logging.Logger,
    ) -> type[EventCollections]:
        """Initialise class attributes for this and dependent classes.

        Should be called before the instance creation.

        :param params: Input parameters.
        :param end: Microtubule end: 0 or 1.
        :param path_out: Data output path.
        :param logger: Logger
        :return: This class with class attributes initialized.
        """

        cls.logger = logger
        cls.END = end
        cls.path_out = path_out
        cls.PARAMS = params

        Single.logger = cls.logger
        Single.END = cls.END
        Single.path_out = path_out
        Single.PARAMS = params

        Double.logger = cls.logger
        Double.END = cls.END
        Double.path_out = path_out
        Double.PARAMS = params

        Triple.logger = cls.logger
        Triple.END = cls.END
        Triple.path_out = path_out
        Triple.PARAMS = params

        Histogram.logger = cls.logger

        return cls

    def report_cycle_uncorrelated(
            self,
            show: bool,
    ) -> dict:
        """Prepare and report the microtubule lifecycle data.

        The lifecycle is defined as a sum of shrink growth and pause
        assumed to be independent phases, i.e without having specific
        order. (For the ordered case see 'report_cycle_correlated'.)
        The relative phase contributions are visualised and exported as
        a pie plot.

        :param show: If True, display the plots, otherwise keep hidden.
        :return: dict
            'duration' - overall duration,
            'time_fraction' - relative contributions of growth, pause
            and shrink phases
        """

        cycle = np.array([self.shrinks.duration.avg,
                          self.grows.duration.avg,
                          self.pauses.duration.avg])
        cycle_duration = Stats(cycle.sum(),
                               cycle.std(), 'sec')
        self.logger.info(f"Uncorrelated {self.region} cycle_duration: "
                         f"{cycle_duration.avg} ± "
                         f"{cycle_duration.std} sec")

        growth_dur_frac = self.grows.duration.avg / cycle_duration.avg \
            if not np.isclose(cycle_duration.avg, 0.) \
            else np.nan
        self.logger.info(f"Uncorrelated {self.region} "
                         f"growth time fraction: {growth_dur_frac}")

        shrink_dur_frac = self.shrinks.duration.avg / cycle_duration.avg \
            if not np.isclose(cycle_duration.avg, 0.) \
            else np.nan

        self.logger.info(f"Uncorrelated {self.region} "
                         f"shrink time fraction: {shrink_dur_frac}")

        pause_dur_frac = self.pauses.duration.avg / cycle_duration.avg \
            if not np.isclose(cycle_duration.avg, 0.) \
            else np.nan
        self.logger.info(f"Uncorrelated {self.region} "
                         f"pause time fraction: {pause_dur_frac}")

        fractions = [growth_dur_frac, shrink_dur_frac, pause_dur_frac]
        if not np.isnan(fractions).all():
            exportfile = self.path_out / \
                         f'uncorrelated_{self.region}_time_fracts_e{self.END}'
            vis.pie(fractions,
                    ['growth', 'shrink', 'pause'],
                    title=f"Uncorrelated {self.region} time fractions",
                    exportfile=exportfile,
                    show=show)
        else:
            self.logger.info(f'No growth/shrink cycle is detected '
                             f'for end {self.END}: '
                             f'omitting Time Fractions plot.')

        return {
            'duration': {'avg': cycle_duration.avg,
                         'std': cycle_duration.std},
            'time_fraction': {'growth': growth_dur_frac,
                              'shrink': shrink_dur_frac,
                              'pause': pause_dur_frac}
        }

    def report_cycle_correlated(
            self,
            show: bool,
    ) -> dict:
        """Prepare and report the microtubule lifecycle data.

        The lifecycle is defined as a sequence of consecutive shrink,
        growth and, when present, pause phases following one another as
        triples of history tracks. The relative phase contributions are
        visualised and exported as a pie plot.

        :param show: If True, display the plots, otherwise keep hidden.
        :return: dict
            'duration' - overall duration,
            'time_fraction' - relative contributions of growth, pause
            and shrink phases.
        """

        # shrink -> growth:
        sg_dur = self.sg.report('duration')
        sgda = sg_dur['duration']['avg']

        sg_g_dur_frac = self.grows.duration.avg / sgda \
            if not np.isclose(sgda, 0.) \
            else np.nan
        self.logger.info(f"Correlated {self.region} "
                         f"growth time fraction sg cycle: {sg_g_dur_frac}")

        sg_s_dur_frac = self.shrinks.duration.avg / sgda \
            if not np.isclose(sgda, 0.) \
            else np.nan
        self.logger.info(f"Correlated {self.region} "
                         f"shrink time fraction sg cycle: {sg_s_dur_frac}")

        # growth -> shrink:
        gs: Final = self.gs
        gs_dur = gs.report('duration')
        gsda = gs_dur['duration']['avg']

        gs_g_dur_frac = self.grows.duration.avg / gsda \
            if not np.isclose(gsda, 0.) \
            else np.nan
        self.logger.info(f"Correlated {self.region} "
                         f"growth time fraction gs cycle: {gs_g_dur_frac}")

        gs_s_dur_frac = self.shrinks.duration.avg / gsda \
            if not np.isclose(gsda, 0.) \
            else np.nan
        self.logger.info(f"Correlated {self.region} "
                         f"shrink time fraction gs cycle: {gs_s_dur_frac}")

        # growth -> shrink and growth -> pause -> shrink:
        gps: Final = self.gps
        gps_dur = gps.report('duration')
        gsda = gs_dur['duration']['avg']
        gpsda = gps_dur['duration']['avg']
        g0s_num = gs.num() + gps.num()

        if g0s_num and \
           not np.isnan(gsda) and \
           not np.isnan(gpsda):
            g0s_dur = Stats(
                (gs.num() * gsda +
                 gps.num() * gpsda) / g0s_num,
                np.sqrt((gs.num() * gs_dur['duration']['std'])**2 +
                        (gps.num() * gps_dur['duration']['std'])**2) / g0s_num,
                'sec'
            )
        else:
            g0s_dur = Stats(np.nan, np.nan, '')
        self.logger.info(f"Correlated {self.region} cycle g0s duration: "
                         f"{g0s_dur.avg} ± "
                         f"{g0s_dur.std} "
                         f"{g0s_dur.units}")

        g0s_g_dur_frac = self.grows.duration.avg / g0s_dur.avg \
            if not np.isclose(g0s_dur.avg, 0.) \
            and not np.isnan(g0s_dur.avg) \
            else np.nan
        self.logger.info(f"Correlated {self.region} "
                         f"growth time fraction g0s cycle: {g0s_g_dur_frac}")

        g0s_s_dur_frac = self.shrinks.duration.avg / g0s_dur.avg \
            if not np.isclose(g0s_dur.avg, 0.) \
            and not np.isnan(g0s_dur.avg) \
            else np.nan
        self.logger.info(f"Correlated {self.region} "
                         f"shrink time fraction g0s cycle: {g0s_s_dur_frac}")

        g0s_p_dur_frac = self.pauses.duration.avg / g0s_dur.avg \
            if not np.isclose(g0s_dur.avg, 0.) \
            and not np.isnan(g0s_dur.avg) \
            else np.nan
        self.logger.info(f"Correlated {self.region} "
                         f"pause time fraction g0s cycle: {g0s_p_dur_frac}")

        fractions = [g0s_g_dur_frac, g0s_s_dur_frac, g0s_p_dur_frac]
        if not np.isnan(fractions).all():
            exportfile = self.path_out / f'Correlated_{self.region}_' \
                                         f'time_fractions_g0s_e{self.END}'
            vis.pie(fractions,
                    ['growth', 'shrink', 'pause'],
                    title=f"Correlated {self.region} time fractions ",
                    exportfile=exportfile,
                    show=show)
        else:
            self.logger.info(f'No growth/shrink cycle is detected '
                             f'for {self.region} end {self.END}: '
                             f'omitting Time Fractions plot.')

        return {
            'sg': {
                'duration': sg_dur['duration'],
                'time_fraction': {'growth': sg_g_dur_frac,
                                  'shrink': sg_s_dur_frac}
            },
            'gs': {
                'duration': gs_dur['duration'],
                'time_fraction': {'growth': gs_g_dur_frac,
                                  'shrink': gs_s_dur_frac}
            },
            'g0s': {
                'duration': {'avg': g0s_dur.avg,
                             'std': g0s_dur.std},
                'time_fraction': {'growth': g0s_g_dur_frac,
                                  'shrink': g0s_s_dur_frac,
                                  'pause': g0s_p_dur_frac}
            },
        }

    def report_by_type(
            self,
            show: bool,
    ) -> dict:
        """Characteristics of the microtubule end states.

        Analyses:
        - duration: time until the transition to a new state
        - elongation (μm, may be positive or negative): microtubule
        length state incurred at this end.
        - reorientation: change in the end direction over the state
        lifetime cycle characteristics.

        :param show: If True, display the plots, otherwise keep hidden.
        :return: Summary of the state cheracteristics.
        """

        self.logger.info(f'End {self.END} {self.region} states by type: ')

        res = dict()

        # Pure event types.
        # These are single events and sequences of two or tree ones:
        if sum(t.num() for t in self.pure_s) == 0:
            self.logger.warning(f'WARNING: No event of any type detected in '
                                f'cell {self.region}: are region borders '
                                f'meaningful?.')

        for tt in (self.pure_s, self.pure_d, self.pure_t):
            for t in tt:
                name = t.name()
                res[name] = {}
                res[name] |= {'num': t.num()}
                if t.num():
                    res[name] |= t.report('duration')
                    res[name] |= t.report('elongation')
                    res[name] |= t.report('velocity')
                    if self.include_ornt:
                        res[name] |= t.report('reorientation')
                else:
                    self.logger.info(f'No events detected for type {name} '
                                     f'in {self.region}: skipping report.')
        # Compositions of events:
        # Growth:
        res['growth'] = {}
        res['growth'] |= {'num': self.grows.num()}
        if self.grows.num():
            res['growth'] |= self.grows.report('duration')
            res['growth'] |= Elongation.create(self.grows, self.region) \
                                       .pipeline(show)
            res['growth'] |= Velocity.create(self.grows, self.region)\
                                     .pipeline(show)
            if self.include_ornt:
                res['growth'] |= self.grows.\
                    report_reorientation(self.region, show)

        # Shrink:
        res['shrink'] = {}
        res['shrink'] |= {'num': self.shrinks.num()}
        if self.shrinks.num():
            res['shrink'] |= self.shrinks.report('duration')
            res['shrink'] |= Elongation.create(self.shrinks, self.region)\
                                       .pipeline(show)
            res['shrink'] |= Velocity.create(self.shrinks, self.region)\
                                     .pipeline(show)
            if self.include_ornt:
                res['shrink'] |= self.shrinks.\
                    report_reorientation(self.region, show)

        # Pause:
        res['pause'] = {}
        res['pause'] |= {'num': self.pauses.num()}
        if self.pauses.num():
            res['pause'] |= self.pauses.report('duration')
            res['pause'] |= self.pauses.report('elongation')
            res['pause'] |= self.pauses.report('velocity')
            if self.include_ornt:
                res['pause'] |= self.pauses.report('reorientation')

        # Compositions of event sequences:
        res['cycle_uncorrelated'] = \
            self.report_cycle_uncorrelated(show) if self.END else None
        res['cycle_correlated'] = \
            self.report_cycle_correlated(show) if self.END else None

        return {self.region: res}

    def plot_spatial(
            self,
            mesh: meshio.Mesh,
            start: float,
            duration: float,
            show: bool,
    ) -> None:
        """Plot cell map showing positions of pause and rescue events.

        :param mesh: Plasma membrane mesh.
        :param start: Initial time after the start of history recording
            (sec) for the events to show.
        :param duration: Duration of the interval over which the events
            are plotted.
        :param show: If True, display the plot, otherwise keep hidden.
        """

        from matplotlib import pyplot as plt
        from matplotlib.patches import Polygon
        from matplotlib.collections import PatchCollection

        patches = [Polygon(mesh.points[c, :2])
                   for c in mesh.cells_dict['triangle']]
        p = PatchCollection(patches,
                            edgecolor=None,
                            facecolor=(0.9, 0.9, 0.9, 0.2))

        fig, ax = plt.subplots()
        ax.add_collection(p)

        pauses_inds = np.where(
            (self.pauses.time >= start) &
            (self.pauses.time < start + duration)
        )[0]
        ax.scatter(self.pauses.pos_fr[pauses_inds, 0],
                   self.pauses.pos_fr[pauses_inds, 1],
                   s=1, marker='.', c='y')
        recoveries_inds = np.where(
            (self.grows.time >= start) &
            (self.grows.time < start + duration)
        )[0]
        ax.scatter(self.grows.pos_fr[recoveries_inds, 0],
                   self.grows.pos_fr[recoveries_inds, 1],
                   s=1, marker='.', c='g')
        plt.axis('square')
        plt.axis('off')
        vis.save_plt(self.path_out / f'transitions_e{self.END}',
                     gzipped=True)
        if show:
            plt.show()
        else:
            plt.close(fig)

    def radial_distribution_of_growing_ends(
            self,
            edges: np.ndarray
    ) -> np.ndarray:
        """Density distribution of growing microtubule ends.

        Density of growing microtubule ends as a function of distance to
        cell center in xy plane.

        :param edges: Distance histogram bins.
        :return: Frequency histogram of microtubule ends.
        """

        h = np.zeros_like(edges[:-1], dtype=float)
        for f, t in zip(self.grows.dist0_fr,
                        self.grows.dist0_to):
            h[np.where((f < edges) & (edges <= t))[0]] += 1

        return h

    def xy_shrink_growth_ratio(
            self,
            edges: np.ndarray,
            show: bool,
    ) -> None:
        """Generate spatial maps of event positions and their ratios.

        :param edges: x,y histogram bins.
        :param show: If True, display the maps.
        """

        import matplotlib.colors as colors

        def norm(v):
            return v / max(abs(v.min().min()), abs(v.max().max()))

        def message(name, reason):
            self.logger.info(f'Plotting {name} was skipped for '
                             f'end {self.END} because of {reason} values.')

        c = self.shrinks.xy_distribution(edges).T
        r = self.grows.xy_distribution(edges).T

        norm11 = colors.Normalize(vmin=-1., vmax=1.)
        norm01 = colors.Normalize(vmin=0., vmax=1.)

        # Map catastrophes:
        if ~np.isclose(c, 0.).all():
            vis.colormesh_plot(
                norm(c),
                title=f'catastrophes at end {self.END}',
                cmap='gray',
                norm=norm01,
                exportfile=self.path_out / f'catastrophes_e{self.END}',
                show=show,
            )
        else:
            message('catastrophes', 'all zero')

        # Map recoveries:
        if ~np.isclose(r, 0.).all():
            vis.colormesh_plot(
                norm(r),
                title=f'recoveries at end {self.END}',
                cmap='gray',
                norm=norm01,
                exportfile=self.path_out / f'recoveries_e{self.END}',
                show=show,
            )
        else:
            message('recoveries', 'all zero')

        # Map relations between catastrophe and recovery frequensies:
        cc = c * abs(self.shrinks.velocity.avg)
        rr = r * abs(self.grows.velocity.avg)
        h = np.divide(rr - cc, rr + cc, out=np.zeros_like(rr),
                      where=(~np.isnan(cc) * ~np.isnan(rr)*((rr + cc) != 0.)))
        if ~np.isclose(h, 0.).all():
            exportfile = self.path_out / \
                         f'relative_shrinks_vs_growth_e{self.END}'
            vis.colormesh_plot(
                norm(h),
                title=f'relative shrinks vs growth at end {self.END}',
                cmap='RdBu',
                norm=norm11,
                exportfile=exportfile,
                show=show,
            )
        else:
            message('shrinks relative to growth', 'all zero')
        if (~np.isnan(cc) * ~np.isnan(rr)).all():
            vis.colormesh_plot(
                norm(rr - cc),
                title=f'diff shrinks and growth at end {self.END}',
                cmap='RdBu',
                norm=norm11,
                exportfile=self.path_out/f'diff_shrinks_growth_e{self.END}',
                show=show,
            )
            vis.colormesh_plot(
                norm(rr + cc),
                title=f'sum shrinks and growth at end {self.END}',
                cmap='binary',
                norm=colors.LogNorm(),
                exportfile=self.path_out/f'sum_shrinks_growth_e{self.END}',
                show=show,
            )
        else:
            message('diff and sum of shrinks and growth', 'nan')
