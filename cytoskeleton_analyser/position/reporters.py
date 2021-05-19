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

"""Reporter classes.

Classes implementing presentation, storage and retrieval of major geometric
characteristics of the spatial systems. Provides a unified set of
class methods imtended to be used in a similar way.
Two operational modes are possible:
- initial analysis, visualization and export of the results:
R.create(tp, s).pipeline(sp)
- import and visualization of stored analysis results:
R.create(tp, s).restore(p)
where:
'R': the report class name,
'tp': type of the spatial system,
's': flag of interactive visualization,
'sp': list of spatial system instances to be analysed,
'p': path to the stored analysis results.
"""

from __future__ import annotations
import json
from pathlib import Path
from typing import Final, Optional

import numpy as np

import cytoskeleton_analyser.fitting as fit
from ..histograms import Experimental
from ..histograms import Histogram
from ..histograms import Simulated
from ..report import Report
from .spatial_systems import FullDepth
from .spatial_systems import ListOfSpatialSystems


class Features:
    """Classification of reported features
    """

    #: Features applicable to both full and sliced cell representations.
    common: Final[list[str]] = [
        'Lengths3d',
        'Lengths2d',
        'Curvature3d',
        'RadialMass',
        'RadialEnds',
        'AnglesToRad',
        'SegmentNumbers',
    ]

    #: Features applicable to full cell representation alone.
    only_full: Final[list[str]] = [
        'AgesByNode',
        'AgesByFilament',
    ]

    #: Features applicable to sliced cell representation alone.
    only_slice: Final[list[str]] = [
        'Curvature2dConv',
        'Curvature2dMboc17',
    ]

    #: Complete set of implemented features.
    all: Final[list[str]] = common + only_full + only_slice

    @staticmethod
    def is_common(f: str) -> bool:
        """True if ``f`` belongs to set common to both representations.

        :param f: Feature class name.
        """

        return any(f == g for g in Features.common)

    @staticmethod
    def is_full(f: str) -> bool:
        """True if feature ``f`` is applicable to full representation.

        :param f: Feature class name.
        """

        return any(f == g for g in Features.only_full)

    @staticmethod
    def is_slice(f: str) -> bool:
        """True if feature ``f`` is applicable to slice representation.

        :param f: Feature class name.
        """

        return any(f == g for g in Features.only_slice)

    @staticmethod
    def is_any(f: str) -> bool:
        """True if feature ``f`` is applicable to any representation.

        :param f: Feature class name.
        """

        return any(f == g for g in Features.all)

    @staticmethod
    def is_applicable(
            f: str,
            tp: type[FullDepth],
    ) -> bool:
        """True if feature ``f`` is applicable to representation ``tp``.

        :param f: Feature class name.
        :param tp: Feature class name.
        """

        return Features.is_common(f) or \
               tp.type == 'slice' and Features.is_slice(f) or \
               tp.type == 'full' and Features.is_full(f)

    @staticmethod
    def reporter(f: str) :
        """Convert feature name ``f`` to corresponding reporter type.

        :param f: Feature name.
        """

        if not Features.is_any(f):
            raise ValueError(f"{f} is not a valid position_feature.")

        return globals()[f]


class _Report(Report):
    """Adaptation of report.Report class for spatial systems.

    For subclassing specific to reported cytoskeleton attributes.
    """

    tp: type[FullDepth]  #: Type of the spatial system.

    @classmethod
    def _create(
            cls,
            tp: type[FullDepth],
            name: str,
            show: bool = True,
    ) -> None:

        cls.tp = tp
        cls.__create(
            tp.logger,
            tp.paths.data_out,
            name + '_' + tp.type,
            show
        )


class Lengths3d(_Report):
    """Reports of  microtubule lengths in .
    """

    #: Part of figure and report titles.
    LABEL: Final = '3d filament length in'

    @classmethod
    def create(
            cls,
            tp: type[FullDepth],
            show: bool = True,
    ) -> type[Lengths3d]:

        super()._create(tp, __class__.__name__, show)
        cls.units = tp.len_units
        return cls

    @classmethod
    def report(
            cls,
            sp: ListOfSpatialSystems,
    ) -> tuple[Histogram, list, list]:

        data = [s.len_total3d for s in sp]
        avg, std = cls.tp.print_avgstd(cls.LABEL, data, cls.units)

        h = Histogram(
            cls.name,
            Simulated().initialise(data, cls.fits_sim, dx=0.4, density=True)
        )
        h.to_csv(cls.path_out)
        cls.plot(h)
        cls.logger.info('')
        return h, [avg], [std]

    @classmethod
    def plot(
            cls,
            h: Histogram
    ) -> None:

        h.plot(
            cls.LABEL + ' ' + cls.tp.type,
            xlabel=f'length ({cls.units})',
            xlim=[0., 60.],
            save_path=cls.path_out,
            show=cls.show
        )

    @classmethod
    def pipeline(
            cls,
            sp: ListOfSpatialSystems,
    ) -> None:

        cls.fits_sim = [
            # (fit.Gamma.loc0, [2., 0.1, 0.]),
            # (fit.Weibull.full, [2., 1.]),
            # (fit.Rayleigh.f, [1.]),
        ]

        rep = cls.report(sp)
        cls.summarize(rep, [0])


class Lengths2d(_Report):
    """Reports lengths of xy projections of the microtubules.

    Examines apparent lengths of simulated microtubules and compares
    them with experimental data obtained using optical microscopy.
    For this purpose, implements experimental sets from microtubule
    length measurements with superresolution methods
    by Zhang et al. (MBoC 2017)
    """

    #: Part of figure and report titles.
    LABEL: Final = 'length of filament 2d projections in'

    @classmethod
    def create(
            cls,
            tp: type[FullDepth],
            show: bool = True,
    ) -> type[Lengths2d]:

        super()._create(tp, __class__.__name__, show)
        cls.units = tp.len_units
        return cls

    @classmethod
    def _experimental(
            cls,
            cell_type: str,
    ):

        import cytoskeleton_analyser.position.empirical_data.mboc17 as mboc17

        bc, (contr, ca_ras) = mboc17.length(density=True)
        if cell_type == 'RW_Protr':
            h = ca_ras
        elif cell_type == 'SpreRou':
            h = contr
        avg = mboc17.avg(bc, h)
        cls.logger.info('\nEmpirical length of filament 2d projections in ' +
                        f'{cls.tp.type}: {avg} {cls.units}')

        return bc, h, avg

    @classmethod
    def report(
            cls,
            sp: ListOfSpatialSystems,
    ) -> tuple[Histogram, list, list]:

        data = [s.len_total2d for s in sp]
        avg, std = cls.tp.print_avgstd(cls.LABEL, data, cls.units)

        ct = cls.tp.params['cell'].typename
        if cls.tp.type == 'slice' and \
                (ct == 'RW_Protr' or
                 ct == 'SpreRou'):
            bc, l, l_avg = cls._experimental(ct)
            e = Experimental().initialise((bc, l), cls.fits_exp)
            h = Histogram(
                cls.name,
                Simulated().initialise(
                    data, cls.fits_sim, dx=0.4, exper_bc=e.bc, density=True),
                experimental=e
            )
            avg, std = [avg, l_avg], [std, np.nan]
        else:
            h = Histogram(
                cls.name,
                Simulated().initialise(
                    data, cls.fits_sim, dx=0.4, density=True),
            )
            avg, std = [avg], [std]
        h.to_csv(cls.path_out)
        cls.plot(h)
        cls.logger.info('')
        return h, avg, std

    @classmethod
    def plot(
            cls,
            h: Histogram
    ) -> None:

        h.plot(
            cls.LABEL + ' ' + cls.tp.type,
            xlabel=f'length ({cls.units})',
            xlim=[0., 60.],
            save_path=cls.path_out,
            show=cls.show
        )

    @classmethod
    def pipeline(
            cls,
            sp: ListOfSpatialSystems,
    ) -> None:

        best = []
        if cls.tp.type == 'full':
            cls.fits_sim = [
#                (fit.Gamma.loc0, [1., 1, 0.]),
#                (fit.Weibull.full, [2., 3.]),
            ]
            best = [0]

        if cls.tp.type == 'slice':
            e = fit.Exponential.create()
            p = fit.Exponential.Pars
            tu = cls.units

            cls.fits_exp = [
#                (e.d_h, p(a=1., tau1=2.), tu),
#                (fit.Gamma.loc0, [1., 1, 0.]),
#                (fit.Weibull.full, [2., 3.]),
            ]
            cls.fits_sim = [
#                (e.d_h, p(a=1., tau1=2.), tu),
#                (fit.Gamma.loc0, [1., 1, 0.]),
#                (fit.Weibull.full, [2., 3.]),
            ]
            best = [1, 1]

        rep = cls.report(sp)
        cls.summarize(rep, best)


class RadialMass(_Report):
    """Reports distribution of microtubule masss.

    Microtubule mass is analysed as a function of distance to cell
    center in xy plane.
    """

    #: Part of figure and report titles.
    LABEL: Final = 'mass vs distance to center '

    @classmethod
    def create(
            cls,
            tp: type[FullDepth],
            show: bool = True,
    ) -> type[RadialMass]:

        super()._create(tp, __class__.__name__, show)
        cls.units = tp.len_units
        return cls

    @classmethod
    def report(
            cls,
            sp: ListOfSpatialSystems
    ) -> tuple[Histogram, list, list]:

        data = [np.concatenate(s.center_dist_2d) for s in sp]
        avg, std = cls.tp.print_avgstd(cls.LABEL, data, cls.units)

        h = Histogram(
            cls.name,
            Simulated().initialise(
                data, fits=cls.fits_sim, dx=0.25, density=True
            ),
        )
        h.to_csv(cls.path_out)
        cls.plot(h)
        cls.logger.info('')
        return h, [avg], [std]

    @classmethod
    def plot(
            cls,
            h: Histogram
    ):

        h.plot(
            cls.LABEL + cls.tp.type,
            xlabel=f'length ({cls.units})',
            xlim=[0., 30.],
            save_path=cls.path_out,
            show=cls.show,
        )

    @classmethod
    def pipeline(
            cls,
            sp: ListOfSpatialSystems,
    ) -> None:

        cls.report(sp)


class RadialEnds(_Report):
    """Reports positions of of microtubule plus ends.

    Analyse the distribution of microtubule plus ends as a function
    of distance to cell center in xy plane.
    """

    #: Part of figure and report titles.
    LABEL: Final = 'plus-ends vs distance to center '

    @classmethod
    def create(
            cls,
            tp: type[FullDepth],
            show: bool = True,
    ) -> type[RadialEnds]:

        super()._create(tp, __class__.__name__, show)
        cls.units = tp.len_units
        return cls

    @classmethod
    def report(
            cls,
            sp: ListOfSpatialSystems,
    ) -> tuple[Histogram, list, list]:

        data = [np.concatenate(s.center_dist_2d_ends) for s in sp]
        avg, std = cls.tp.print_avgstd(cls.LABEL, data, cls.units)

        h = Histogram(
            cls.name,
            Simulated().initialise(
                data, fits=cls.fits_sim, dx=0.25, density=True)
        )
        h.to_csv(cls.path_out)
        cls.plot(h)
        cls.logger.info('')
        return h, [avg], [std]

    @classmethod
    def plot(
            cls,
            h: Histogram
    ):

        h.plot(
            cls.LABEL + cls.tp.type,
            xlabel=f'length ({cls.units})',
            xlim=[0., 30.],
            save_path=cls.path_out,
            show=cls.show,
            )

    @classmethod
    def pipeline(
            cls,
            sp: ListOfSpatialSystems,
    ) -> None:

        cls.report(sp)


class SegmentNumbers(_Report):
    """Reports apparent number of microtubules.

    Measure statistics on apparent number of microtubules in full system
    and as visible in TIRF microscopy observations.
    """

    @classmethod
    def create(
            cls,
            tp: type[FullDepth],
            _ = None,
    ) -> type[SegmentNumbers]:

        super()._create(tp, __class__.__name__)
        cls.units = 'segments'
        return cls

    @classmethod
    def report(
            cls,
            sp: ListOfSpatialSystems,
    ) -> tuple[Optional[Histogram], list, list]:

        cls.logger.info(f"Number of filaments in {cls.tp.type}:")
        data = np.array([len(s.len_total3d) for s in sp])
        [cls.logger.info(f'\t {n}') for n in data]
        avg = np.mean(data)
        std = np.std(data)
        cls.logger.info(f"overall: {avg} ± {std} slices\n")

        fname = cls.path_out / f"{cls.name}.json"
        with open(fname, 'w') as f:
            json.dump({'num': {'avg': avg, 'std': std}}, f)

        cls.logger.info('')
        return None, [avg], [std]

    @classmethod
    def restore(
            cls,
            path: Path
    ) -> dict:

        cls.logger.info('Restoring ' + cls.name + ':')
        summary = cls.log_stats(path)

#        fname = f"{cls.path_out}{cls.name}.json"
#        with open(fname, 'w') as f:
#            json.dump({'num': {'avg': avg, 'std': std}}, f)
        cls.logger.info('')

        return {'summary': summary,
                'h_sim_avg': None,
                }

    @classmethod
    def pipeline(
            cls,
            sp: ListOfSpatialSystems,
    ) -> None:

        cls.report(sp)


class Curvature3d(_Report):
    """Reports 3d curvature of microtubule fibers.
    """

    #: Part of figure and report titles.
    LABEL: Final = '3d curvature of filaments in'

    @classmethod
    def create(
            cls,
            tp: type[FullDepth],
            show: bool = True,
    ) -> type[Curvature3d]:

        super()._create(tp, __class__.__name__, show)
        cls.units = '1/' + tp.len_units
        return cls

    @classmethod
    def report(
            cls,
            sp: ListOfSpatialSystems,
    ) -> tuple[Histogram, list, list]:

        data = [s.curv3d for s in sp]
        avg, std = cls.tp.print_avgstd(cls.LABEL, data, cls.units)

        h = Histogram(
            cls.name,
            Simulated().initialise(data, cls.fits_sim, dx=0.02, density=True)
        )
        h.to_csv(cls.path_out)
        cls.plot(h)
        return h, [avg], [std]

    @classmethod
    def plot(
            cls,
            h: Histogram
    ) -> None:

        h.plot(
            cls.LABEL + ' ' + cls.tp.type,
            xlabel=f'curvature ({cls.units})',
            xlim=[0., 1.5],
            save_path=cls.path_out,
            show=cls.show
        )

    @classmethod
    def pipeline(
            cls,
            sp: ListOfSpatialSystems,
    ) -> None:

        cls.fits_sim = [
            (fit.Rayleigh.f, [1.]),
        ]

        rep = cls.report(sp)
        cls.summarize(rep, [0])


class Curvature2dConv(_Report):
    """Reports of apparent curvature of microtubule projections
    to xy plane.
    """

    #: Part of figure and report titles.
    LABEL: Final = '2d curvature of projected filaments in'

    @classmethod
    def create(
            cls,
            tp: type[FullDepth],
            show: bool = True,
    ) -> type[Curvature2dConv]:

        super()._create(tp, __class__.__name__, show)
        cls.units = '1/' + tp.len_units
        return cls

    @classmethod
    def report(
            cls,
            sp: ListOfSpatialSystems,
    ) -> tuple[Histogram, list, list]:

        data = [s.curv2d for s in sp]
        avg, std = cls.tp.print_avgstd(cls.LABEL, data, cls.units)

        h = Histogram(
            cls.name,
            Simulated().initialise(data, cls.fits_sim, dx=0.02, density=True)
        )
        h.to_csv(cls.path_out)
        cls.plot(h)
        return h, [avg], [std]

    @classmethod
    def plot(
            cls,
            h: Histogram
    ) -> None:

        h.plot(
            cls.LABEL + ' ' + cls.tp.type,
            xlabel=f'curvature ({cls.units})',
            xlim=[0., 1.5],
            save_path=cls.path_out,
            show=cls.show
        )

    @classmethod
    def pipeline(
            cls,
            sp: ListOfSpatialSystems,
    ) -> None:

        cls.fits_sim = [
#            (fit.Gamma.loc0, [2., 0.1, 0.]),
#            (fit.Weibull.full, [2., 1.]),
#            (fit.Rayleigh.f, [.5]),
        ]

        rep = cls.report(sp)
        cls.summarize(rep, [0])


# ======================================================================================================================
class Curvature2dMboc17(_Report):
    """Reports microtubule curvatures using formulas applied in the
       processing of superresolution images by Zhang et al. MBoC 2017
    """

    #: Part of figure and report titles.
    LABEL: Final = 'empirical curvature of filament 2d projections in'

    @classmethod
    def create(
            cls,
            tp: type[FullDepth],
            show: bool = True,
    ) -> type[Curvature2dMboc17]:

        name = __class__.__name__
        if tp.type != 'slice':
            tp.logger.warning(f'WARNING: Analysis of {name} makes only '
                               f'sence for Tirf slices.')

        ct = tp.params['cell'].typename
        if ct != 'RW_Protr' and ct != 'SpreRou':
            tp.logger.warning(f'WARNING: Analysis of {name} makes only '
                               f'sence for specific cell types.')

        super()._create(tp, name, show)
        cls.units = '1/' + tp.len_units

        return cls

    @classmethod
    def experimental(
            cls,
            cell_type: str
    ) -> tuple[list, Histogram, float]:

        import cytoskeleton_analyser.position.empirical_data.mboc17 as zh

        bc, (contr, ca_ras) = zh.curvature(density=True)
        if cell_type == 'RW_Protr':
            h = ca_ras
        elif cell_type == 'SpreRou':
            h = contr
        else:
            assert False, 'Wrong Cell Type'
        avg = zh.avg(bc, h)
        cls.logger.info('\nEmpirical curvature of filament 2d projections in '
                        + f'{cls.tp.type}: {avg} {cls.units}')

        return bc, h, avg

    @classmethod
    def report(
            cls,
            sp: ListOfSpatialSystems,
        ) -> tuple[Histogram, list, list]:

        data = [s.curv2d_mboc17 for s in sp]
        avg, std = cls.tp.print_avgstd(cls.LABEL, data, cls.units)

        ct = cls.tp.params['cell'].typename
        if cls.tp.type == 'slice' and \
                (ct == 'RW_Protr' or
                 ct == 'SpreRou'):
            bc, c, c_avg = cls.experimental(ct)
            e = Experimental().initialise((bc, c), cls.fits_exp)
            h = Histogram(
                cls.name,
                Simulated().initialise(
                    data, cls.fits_sim, dx=0.02, exper_bc=e.bc, density=True),
                e,
            )
            avg, std = [avg, c_avg], [std, np.nan]
        else:
            h = Histogram(
                cls.name,
                Simulated().initialise(
                    data, cls.fits_sim, dx=0.02, density=True),
            )
            avg, std = [avg], [std]
        h.to_csv(cls.path_out)
        cls.plot(h)

        return h, avg, std

    @classmethod
    def plot(
            cls,
            h: Histogram
    ) -> None:

        h.plot(
            cls.LABEL + ' ' + cls.tp.type,
            xlabel=f'curvature ({cls.units})',
            xlim=[0., 1.5],
            save_path=cls.path_out,
            show=cls.show
        )

    @classmethod
    def pipeline(
            cls,
            sp: ListOfSpatialSystems,
    ) -> None:

        cls.fits_sim = [
#            (fit.Gamma.loc0, [2., 0.1, 0.]),
#            (fit.Weibull.full, [2., 1.]),
#            (fit.Rayleigh.f, [1.]),
        ]
        cls.fits_exp = [
#            (fit.Gamma.loc0, [2., 0.1, 0.]),
#            (fit.Weibull.full, [2., 1.]),
#            (fit.Rayleigh.f, [1.]),
        ]

        rep = cls.report(sp)
        cls.summarize(rep, [0, 0])


class AnglesToRad(_Report):
    """Reports distribution of angles between points on the microtubule
    and local radial direction.
    """

    #: Part of figure and report titles.
    LABEL: Final = 'angle to radial direction in'

    @classmethod
    def create(
            cls,
            tp: type[FullDepth],
            show: bool = True,
    ) -> type[AnglesToRad]:

        super()._create(tp, __class__.__name__, show)
        cls.units = 'grad'
        cls.is_polar = True
        cls.is_halfpolar = True
        return cls

    @classmethod
    def report(
            cls,
            sp: ListOfSpatialSystems,
    ) -> tuple[Histogram, list, list]:

        data = [s.threshold_radial_dev(cls.tp.params['cell'].
                                       regions.lamella.is_inside)
                for s in sp]
        data = [np.array([2.*np.pi - d if d >= np.pi else d for d in dd])
                for dd in data]
        avg, std = cls.tp.print_avgstd(cls.LABEL,
                                       [d/np.pi*180. for d in data], cls.units)

        h = Histogram(
            cls.name,
            Simulated()
                .initialise(data, cls.fits_sim, dx=2.*np.pi/180., density=True,
                            polar=cls.is_polar, halfpolar=cls.is_halfpolar),
        )
        h.to_csv(cls.path_out)
        cls.plot(h)
        return h, [avg], [std]

    @classmethod
    def plot(
            cls,
            h: Histogram
    ):

        h.plot(
            cls.LABEL + ' ' + cls.tp.type,
            xlabel=f'angle ({cls.units})',
            xlim=[0., 2. * np.pi],
            save_path=cls.path_out,
            show=cls.show
        )

    @classmethod
    def pipeline(
            cls,
            sp: ListOfSpatialSystems,
    ):

        cls.fits_sim = [
            (fit.VonMisesDouble.full, [0.6, 0.3, 32., 5., np.pi / 3.])
        ]

        rep = cls.report(sp)
        cls.summarize(rep, [0])


class AgesByNode(_Report):
    """Report ages (time after the polymerization event) of microtubules
    by filament nodes: individual node ages.
    """

    #: Part of figure and report titles.
    LABEL: Final = 'node ages in'

    @classmethod
    def create(
            cls,
            tp: type[FullDepth],
            show: bool = True,
    ) -> type[AgesByNode]:

        super()._create(tp, __class__.__name__, show)
        cls.units = 'sec'
        return cls

    @classmethod
    def report(
            cls,
            sp: ListOfSpatialSystems,
    ) -> tuple[Histogram, list, list]:

        data = [s.ages_cumulative for s in sp]
        avg, std = cls.tp.print_avgstd(cls.LABEL, data, cls.units)

        h = Histogram(
            cls.name,
            Simulated().initialise(data, cls.fits_sim, dx=10., density=True),
        )
        h.to_csv(cls.path_out)
        cls.plot(h)
        return h, [avg], [std]

    @classmethod
    def plot(
            cls,
            h: Histogram
    ) -> None:

        h.plot(
            cls.LABEL + ' ' + cls.tp.type,
            xlabel=f'age ({cls.units})',
            xlim=[0., 10000],
            yscale='log',
            save_path=cls.path_out,
            show=cls.show
        )

    @classmethod
    def pipeline(
            cls,
            sp: ListOfSpatialSystems,
    ) -> None:

        e = fit.Exponential.create()
        p = fit.Exponential.Pars
        tu = cls.units
        a = 1.
        tau1 = 10
        tau2 = 1000
        cls.fits_sim = [
            (e.d_h, p(a=a, tau1=tau1), tu),
            (e.d_d_h, p(a=a, tau1=tau1, b=0.9, tau2=tau2), tu),
        ]

        rep = cls.report(sp)
        cls.summarize(rep, [1])


class AgesByFilament(_Report):
    """Report ages of microtubules.

    Age is defined as time passed after the polymerization event.
    Here node age averages over all filament nodes are considered.
    """

    #: Part of figure and report titles.
    LABEL: Final = 'filament ages in'

    @classmethod
    def create(
            cls,
            tp: type[FullDepth],
            show: bool = True,
    ) -> type[AgesByFilament]:

        super()._create(tp, __class__.__name__, show)
        cls.units = 'sec'

        return cls

    @classmethod
    def report(
            cls,
            sp: ListOfSpatialSystems,
    ) -> tuple[Histogram, list, list]:

        data = [s.ages_by_filament for s in sp]
        avg, std = cls.tp.print_avgstd(cls.LABEL, data, cls.units)

        h = Histogram(
            cls.name,
            Simulated().initialise(data, cls.fits_sim, dx=10., density=True),
        )
        h.to_csv(cls.path_out)
        cls.plot(h)
        return h, [avg], [std]

    @classmethod
    def plot(
            cls,
            h: Histogram
    ) -> None:

        h.plot(
            cls.LABEL + ' ' + cls.tp.type,
            xlabel=f'age ({cls.units})',
            xlim=[0., 10000],
            yscale='log',
            save_path=cls.path_out,
            show=cls.show,
        )

    @classmethod
    def pipeline(
            cls,
            sp: ListOfSpatialSystems,
    ) -> None:

        e = fit.Exponential.create()
        p = fit.Exponential.Pars
        tu = cls.units
        a = 1.
        tau1 = 10
        tau2 = 1000
        cls.fits_sim = [
            (e.d_h, p(a=a, tau1=tau1), tu),
            (e.d_d_h, p(a=a, tau1=tau1, b=0.9, tau2=tau2), tu),
        ]

        rep = cls.report(sp)
        cls.summarize(rep, [1])
