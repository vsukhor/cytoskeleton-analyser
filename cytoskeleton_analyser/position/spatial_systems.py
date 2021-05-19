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

"""Spatially resolved 3d representations of specific cell systems.
"""

from __future__ import annotations
import logging
import struct

from pathlib import Path
from typing import Callable, Final, Optional, Union

import numpy as np
from scipy.spatial.distance import cdist

from ..inout import Paths
from ..inout import read_to_dtype
from ..plasma_membrane import PlasmaMembrane


class FullDepth:
    """Instantaneous snapshots of complete cell microtubule subsystem
       along with the bounding plasma membrane.

       The snapshots are uncorrelated data sets recording steady states
       run of the Monte Carlo system evolving in the coarse of a
       simulation at intervals ensuring absence of correlations between
       consecutive recordings.
       The representation is assembled from raw simulation data files.
       It uses files named 'positions_*' and 'csk_ages_*'.
    """

    #: Parameters of general applicability.
    params: dict

    #: Short desription.
    type: str

    #: Indexes of the time snapshots available for the current run.
    snaps: list[int]

    #: Collection of file system paths relevant for data input-output.
    paths: Paths

    #: Plasma membrane delimiting cell volume.
    plasma_membrane: PlasmaMembrane

    #: Holds instances of the system representtions at different times.
    all_items: ListOfSpatialSystems

    #: Python standard logger.
    logger: logging.Logger

    #: Units of length.
    len_units: Final[str] = 'μm'

    def __init__(
            self,
            snap_ind: Optional[int],
            rind: int
    ):
        """:param snap_ind: Index of the current snapshot.
           :param rind: Index of the simulation rum.
        """

        #: Snap index.
        self.snap_ind: int = snap_ind

        #: Simulation run index.
        self.rind: int = rind

        #: Monte Carlo iteration.
        self.iteration: int = np.uint64()

        #: Time of the snapshot after simulation start.
        self.time: np.float64 = np.float64()

        #: Origin position of cell coordinate system.
        self.origin: np.ndarray = np.zeros(3)     #: Position of cell center.

        #: Total number of polymerised tubulin units.
        self.mtmass: np.uint64 = np.uint64()

        #: Number of filaments.
        self.nfilaments: np.uint64 = np.uint64()

        #: Nimber of tubulin nodes per filament.
        self.nnodes = None

        cof = 'coarse_' if self.params['iscoarse'] else 'fine_'
        si = str(self.snap_ind) + '_' if self.snap_ind is not None else ''
        #: Cell and run human readable.
        self.signature: str = \
            self.params['cell'].typename + \
            f" plm {self.params['cell'].plmind}" \
            f" csk {self.rind}\n"

        #: File name of the file containing node positions.
        self.fname_pos: Path = \
            self.paths.run / f"positions_{cof}{si}{self.rind}"
        self.fname_ages: Path = \
            self.paths.run / f"csk_ages_{cof}{si}{self.rind}"
        self.figtitle3d: str = 'MT ' + self.type + ': \n' + self.signature

        # Positional coordinates and measures:
        #: 3d positions of filament nodes.
        self.pos = []
        #: Lengths in 3d.
        self.len_total3d = np.empty(0)
        #: Lengths of xy projections.
        self.len_total2d = np.empty(0)

        # Filament curvatures:
        self.curv3d = None
        self.curv2d = None
        self.curv2d_mboc17 = None

        # Distances to cell center.
        self.center_dist_2d = None
        self.center_dist_2d_ends = None
        self.angles_radius = None
        self.ages = []
        self.ages_cumulative = np.empty(0, dtype=float)
        self.ages_by_filament = np.empty(0, dtype=float)

    def read_positions(self) -> None:
        """Read in the 3d spatial coordinates for microtubule nodes.

        Data a read from a raw binary file produced in the simulation
        (standard naming starting with ``positions_``).
        This populates 'self.pos' list, storing the node coordinates
        per filament. Elements of this list are np.ndarrays containing
        xyz coordinates of the nodes forming the microtubule.
        The imported simulation snapshot contains also accompanying data
        (iteration index, time, total microtubule mass, etc.)
        """

        self.logger.info('Importing positions from ' +
                         str(self.fname_pos) + ' ...')
        with open(self.fname_pos, 'rb') as f:
            self.iteration = read_to_dtype(self.iteration, f)
            self.time = read_to_dtype(self.time, f)
            read_to_dtype(np.float32(), f)  # cell_rad
            self.mtmass = read_to_dtype(self.mtmass, f)
            self.nfilaments = read_to_dtype(self.nfilaments, f)

            self.nnodes = np.zeros(self.nfilaments, dtype=np.uint32)
            for i in range(self.nfilaments):
                self.nnodes[i] = read_to_dtype(self.nnodes[i], f)
                p = np.array(struct.unpack('fff' * self.nnodes[i],
                                           f.read(3 * 4 * self.nnodes[i])))
                self.pos.append(p.reshape((self.nnodes[i], 3)))
        self.logger.info('Positions import finished.')

        self.logger.info(f'Snapshot at time: {self.time} sec:')
        self.logger.info(f'\titeration {self.iteration}')
        self.logger.info(f'\ttotal mass: {self.mtmass} nodes')
        self.logger.info(f'\tnumber of filaments: {self.nfilaments}')
        self.logger.info('')

    def read_ages(self) -> None:
        """Read in ages of microtubule nodes at a current time snapshot.

        Data are read from the raw binary file produced in the
        simulation (standard naming starting with ``csk_ages_``) and
        store the in 'self.ages' list on a per-microtubule basis.
        """

        self.logger.info(f'Importing ages from {self.fname_ages} ...')
        with open(self.fname_ages, 'rb') as f:
            time = np.float64
            time = read_to_dtype(time, f)
            assert time == self.time
            nfilaments = np.uint32
            nfilaments = read_to_dtype(nfilaments, f)
            assert nfilaments == self.nfilaments
            nnodes = np.zeros(nfilaments, dtype=np.uint32)
            self.ages = [np.empty(nn, dtype=np.float32) for nn in self.nnodes]
            for i in range(self.nfilaments):
                nnodes[i] = read_to_dtype(nnodes[i], f)
                assert nnodes[i] == self.nnodes[i]
                self.ages[i] = np.array(read_to_dtype(self.ages[i],
                                                      f, nnodes[i]))

        self.logger.info('Ages import finished.')
        self.logger.info('')

        self.ages_cumulative = np.concatenate(self.ages)
        self.ages_by_filament = np.array([f.mean() for f in self.ages])

    def set_lengths(self) -> None:
        """Initialise class attributes related to microtubule dimensions.
        """

        def norm_per_element(pp, d) -> list:
            """Euclidean norm of dimension d."""
            return [np.linalg.norm(np.diff(p[:, :d], axis=0), axis=1)
                    for p in pp]

        def len_total(nrm):
            return np.array([np.sum(n) for n in nrm])

        self.len_total3d = len_total(norm_per_element(self.pos, 3))
        self.len_total2d = len_total(norm_per_element(self.pos, 2))

    def set_curvatures(self) -> None:
        """Initialise filament curvature-related attributes.

        Specify curvature over filament internal nodes:
        curv3d - curvature in 3d space
        curv2d - apparent curvature of filament projections onto xy plane
        curv2d_mboc17 - curvature according specific algorithm adopted
        for image proceccing by Zhang et al. 2017
        """

        from cytoskeleton_analyser.position.curvature import Curvature, D2, D3

        self.curv3d, _ = Curvature.ordinary(self.pos, D3)
        self.curv2d, _ = Curvature.ordinary(self.pos, D2)
        self.curv2d_mboc17 = Curvature.mboc17(self.pos)

    def initialise(self) -> None:
        """Initialization of derivative descriptors of the microtubule system.
        """

        self.read_positions()
        self.set_lengths()
        self.set_center_distance_2d()
        self.set_curvatures()
        self.set_radial_dev()
        self.read_ages()

    def set_center_distance_2d(self) -> None:
        """Initialise an array of edge distances to cell center in xy plane.
        """

        self.center_dist_2d = \
            [cdist(f[:, :2], np.array([self.origin[:2]])).T[0]
             for f in self.pos if f.shape[0]]
        end_length = 10
        self.center_dist_2d_ends = [f[-min(f.shape[0], end_length):]
                                    for f in self.center_dist_2d]

    def set_radial_dev(self) -> None:
        """Initialise angles between filaments edges and
           radial direction in xy plane.
        """

        assert self.center_dist_2d is not None

        nf = len(self.pos)
        self.angles_radius = [np.empty(0)] * nf
        for p, cd, k in zip(self.pos,
                            self.center_dist_2d,
                            list(range(nf))):
            if p.shape[0]:
                raddir = np.subtract(p[:, :2], self.origin[:2]).T / cd
                segm = np.array([np.diff(p[:, 0]),
                                 np.diff(p[:, 1])])
                ornt = segm / np.linalg.norm(segm, axis=0)
                self.angles_radius[k] = \
                    np.arctan2(ornt[1, :], ornt[0, :]) - \
                    np.arctan2(raddir[1, :-1], raddir[0, :-1])
                self.angles_radius[k] = \
                    np.where(self.angles_radius[k] < 0.,
                             self.angles_radius[k] + 2. * np.pi,
                             self.angles_radius[k])

    def threshold_radial_dev(
            self,
            is_within_range: Callable[[list[float]], list[bool]],
    ) -> np.ndarray:
        """Filter the vector of xy radial deviations.

        The filter includes only filament edges from within the
        specified interval of distances to cell center.

        :param is_within_range: Selection function returning
            a boolean map to the input list.
        :return: Filtered radial deviations.
        """

        res = np.empty(0)
        for c, a in zip(self.center_dist_2d, self.angles_radius):
            ii = is_within_range(c)
            if ii.size > 1:
                res = np.hstack((res, a[ii[:-1]]))

        return res

    def threshold_tangential_dev(
            self,
            is_within_range: Callable[[list[float]], list[bool]],
    ) -> np.ndarray:
        """Filtered the vector of xy tangential deviations.

        The result includes only the filament edges from within
        the specified interval of distances to cell center.

        :param is_within_range: Selection function returning a boolean
            map to the input list.
        :return: filtered tangential deviations.
        """

        res = np.empty(0)
        for c, a in zip(self.center_dist_2d, self.angles_radius):
            ii = is_within_range(c)
            ii[-1] = False
            aa = a[ii[:-1]]
            over = aa > 90
            aa[over] -= 90
            aa[~over] = 90 - aa[~over]
            res = np.hstack((res, aa))

        return res

    @classmethod
    def create(
            cls,
            paths: Paths,
            params: dict,
            rind: int,
            logger: Optional[logging.Logger] = None,
            snap_inds: Optional[list[int]] = None,
            init: bool = True,
    ) -> type[FullDepth]:
        """Main method to create the class and all its instances.

        :param paths: Collection of file system paths to retrieve and
            store the data.
        :param params: General-purpose parameters.
        :param rind: Index of the simulation run to process.
        :param logger: logging.Logger object
        :param snap_inds: Indexes of data snapshots avalilable in
            simulation output.
        :param init: Initialise the attribute arrays (in some use cases
            this is not necessary).
        :return: This class along with all its instances.
        """

        cls.params = params
        cls.type = 'full'
        cls.paths = paths

        cls.logger = logger

        cls.plasma_membrane = PlasmaMembrane(paths.plasma_in, params['cell'])

        if init:
            if snap_inds is None:
                snap_inds = cls.list_snaps(cls.params, paths.run)
            if len(snap_inds) == 0:
                cls.params['use_final'] = True

            if cls.params['use_final']:
                cls.all_items = [cls(None, rind)]
            else:
                cls.all_items = [cls(i, rind) for i in snap_inds]

            [s.initialise() for s in cls.all_items]

        return cls

    @classmethod
    def discretize_radius(
            cls,
            nbins: int
    ) -> tuple[np.ndarray, np.ndarray]:
        """Discrretize xy positions in radial direction.

        :param nbins: Desired number of bins.
        :return: Arrays containing positions of bin edges and centers.
        """

        maxdist = np.ceil(cls.plasma_membrane.radial_extent())
        edges = np.linspace(0, maxdist, nbins)

        bincenters = edges[:-1] + (edges[1] - edges[0]) / 2.

        return edges, bincenters

    @classmethod
    def discretize_xy(
            cls,
            nbins: int
    ) -> tuple[np.ndarray, np.ndarray]:
        """Discrretize positions in xy plane as a square matrix.

        :param nbins: Desired number of bins.
        :return: Arrays containing positions of bin edges and centers.
        """

        minxy = cls.plasma_membrane.min_[:2]
        maxxy = cls.plasma_membrane.max_[:2]
        size = np.ceil(max(maxxy - minxy))
        l0 = 1.1 * min(np.floor(minxy))
        l1 = 1.1 * max(np.floor(minxy) + size)

        edges = np.array([np.linspace(l0, l1, nbins)]*2)
        bincenters = np.array([e[:-1] + (e[1] - e[0]) / 2. for e in edges])

        return edges, bincenters

    @staticmethod
    def list_snaps(
            params: dict,
            path: str
    ) -> list[int]:
        """Return snapshot indexes avalilable in simulation output.
        """

        from os import listdir
        from os.path import isfile, join
        import re

        cof = 'coarse_' if params['iscoarse'] else 'fine_'
        files = [f for f in listdir(path)
                 if isfile(join(path, f)) and 'positions_' + cof in f]
        rr = [re.search(r"_[0-9]+_", f) for f in files]

        return sorted([int(r.group()[1:-1]) for r in rr if r is not None])

    @classmethod
    def print_avgstd(
            cls,
            name: str,
            v: list,
            units: str = ''
    ) -> tuple[float, float]:
        """Average and standard deviation.

        Calculate, print and return average and standard deviation
        of the data given in list of ndarrays 'v'. Whenever 'v'
        consists of multiple arrays, average and standard deviation
        values are reported also for the 'v' elements independently.

        :param name: Name of the attribute processed.
        :param v: List of data arrays.
        :param units: Data units.
        :return: Average and standard deviation over the data arrays.
        """

        assert len(v) == 1 or \
               [isinstance(a, type(v[0])) for a in v[1:]]

        cls.logger.info(name + ' ' + cls.type + ': ')

        m = np.array([np.mean(vv) for vv in v])
        s = np.array([np.std(vv) for vv in v])
        [cls.logger.info(f"    snap {sn if sn is not None else 'final'}: "
                         f"distr. mean {mm}  distr. std {ss} {units}")
         for sn, mm, ss, in zip(cls.snaps, m, s)]

        mm = np.mean(m)
        sm = np.std(m)
        cls.logger.info(f"avg. distr. mean: {mm} ± {sm} {units}")

        ms = np.mean(s)
        ss = np.std(s)
        cls.logger.info(f"avg. distr. std:  {ms} ± {ss} {units}")

        return mm, ms

    def plot3d(
            self,
            color_mode: Optional[str] = None,
            with_mesh: bool = False,
            mesh_flattened: bool = False,
            axes_visible: bool = True,
            export: bool = False,
            show: bool = True,
    ) -> None:
        """A matplotlib-based 3d visualization.

        :param color_mode: Colormapping of specific attributes:
            'by_height' - color filament edges accorging to z-position
            'by_age' - color filament edges according to the node age.
            other - red.
        :param with_mesh: include plasma membrane mesh to denote
            internal volume of the cell.
        :param mesh_flattened: if True, show only xy projection of
            the mesh as a background at z = 0.
        :param axes_visible: Show or hide the figure axes and title.
        :param export: If True, export the figure in svg format.
        :param show: If True, display the figure.
        """

        import matplotlib.pyplot as plt
        import mpl_toolkits.mplot3d.art3d as art3d

        # Turn interactive plotting off.
        plt.ioff()

        fig = plt.figure(figsize=(10, 10))
        fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
        ax = fig.gca(projection='3d', proj_type='ortho')

        if axes_visible:
            fig.suptitle(self.figtitle3d)
            labels = self.len_units
            ax.set_xlabel('x (' + labels + ')')
            ax.set_ylabel('y (' + labels + ')')
            ax.set_zlabel('z (' + labels + ')')
        else:
            ax.set_axis_off()
        ax.xaxis.pane.fill = False
        ax.yaxis.pane.fill = False
        ax.zaxis.pane.fill = False

        ax.grid(False)

        axlim = [1.1 * self.plasma_membrane.min_.min(),
                 1.1 * self.plasma_membrane.max_.max()]
        ax.set_xlim3d(axlim[0], axlim[1])
        ax.set_ylim3d(axlim[0], axlim[1])
        ax.set_zlim3d(axlim[0], axlim[1])  # 0., 2. * axlim)

        ax.view_init(azim=0, elev=90)

        if with_mesh:
            if mesh_flattened:
                mvs = np.copy(self.plasma_membrane.mesh.points)
                mvs[:, 2] = 0.
            else:
                mvs = self.plasma_membrane.mesh.points
            mvs = [mvs[c, :]
                   for c in self.plasma_membrane.mesh.cells_dict['triangle']]
            p = art3d.Poly3DCollection(mvs, zsort='min', edgecolor=None,
                                       facecolor=(0.9, 0.9, 0.9, 0.2))
            ax.add_collection3d(p)

        if color_mode == 'by_height':
            pp = [m.reshape(-1, 1, 3) for m in self.pos]
            segs = np.concatenate([np.concatenate([m[:-1], m[1:]], axis=1)
                                   for m in pp])
            cc = np.concatenate([m[:-1, 2] for m in self.pos])
            norm = plt.Normalize(self.plasma_membrane.min_[2],
                                 self.plasma_membrane.max_[2])
            c = [(n, 0., 1. - n) for n in norm(cc)]
            coll = art3d.Line3DCollection(segs, colors=c, lw=0.3)
            ax.add_collection(coll)
            # fig.colorbar(coll, ax=ax)

        elif color_mode == 'by_age':
            pp = [m.reshape(-1, 1, 3) for m in self.pos]
            segs = np.concatenate([np.concatenate([m[:-1], m[1:]], axis=1)
                                   for m in pp])
            cc = np.concatenate([a[:-1] for a in self.ages])
            c = plt.cm.jet(cc / cc.max())
            coll = art3d.Line3DCollection(segs, colors=c, lw=0.3)
            ax.add_collection(coll)

        else:
            for m in self.pos:
                ax.plot(m[:, 0], m[:, 1], m[:, 2], c='r', lw=0.3)

        if export:
            self.export_to_svg(color_mode)

        if show:
            plt.show()
        else:
            plt.close(fig)

    def export_to_svg(
            self,
            color_mode: str
    ) -> None:
        """ Save the system as a figure in svg format.

        Limitation: ignores age- and z position-specific coloring.

        :param color_mode: Colormap indication.
        """

        import gzip
        import copy

        magn = 30.
        axlim = max([1.1 * self.plasma_membrane.min_.min(),
                     1.1 * self.plasma_membrane.max_.max()])
        w = magn * axlim
        h = magn * axlim
        fname = self.paths.data_out / f"cell_{self.type}_{color_mode}.svg.gz"

        with gzip.GzipFile(fname, 'w') as o:
            o.write(
                f"<svg width='{w}' height='{h}' preserveAspectRatio='meet | "
                f"slice' xmlns='http://www.w3.org/2000/svg'>"
                .encode())

            if self.plasma_membrane is not None:
                for f in self.plasma_membrane.mesh.cells_dict['triangle']:
                    p = self.plasma_membrane.mesh.points[f, :2] * magn/2 + w/2
                    o.write("<polygon points='".encode())
                    o.write(f"{p[0, 0]} {p[0, 1]}, "
                            f"{p[1, 0]} {p[1, 1]}, "
                            f"{p[2, 0]} {p[2, 1]}".encode())
                    o.write("' stroke='none' "
                            "fill='rgba(0.9, 0.9, 0.9, 0.07)'/>\n".encode())

            for f in copy.deepcopy(self.pos):
                f *= magn / 2
                f += w / 2
                o.write("<polyline points='".encode())
                o.write(f"{f[0, 0]} {f[0, 1]} ".encode())
                for i in range(1, f.shape[0]):
                    o.write(f"{f[i, 0]} {f[i, 1]} ".encode())
                o.write("' stroke='blue' fill='none' stroke-width='1'/>\n"
                        .encode())
            o.write("</svg>".encode())


class Slice(FullDepth):
    """Encapsulates a system resulting from optical volume sectioning.

       Optical extraction of a cell subvolume between xy-parallel planes.
       This roughly emulates the subsectioning produced by confpcal or
       Total Internal Reflection Fluorescence (TIRF) Microscope:
       a cell content inside a volume slice proximal to basal cell
       surface up to a predefined thickness.
    """

    def __init__(
            self,
            full: FullDepth,
    ):
        """:param full: The original system.
        """

        super().__init__(full.snap_ind, full.rind)

        #: (μm) z-positions of lower and upper limiting planes.
        self.zlimits = self.params['slice_limits']

        self.fi = []     #: filament indexes
        self.ni = []     #: node indexes

        self._extract_from_full(full.pos)
        self.iteration = full.iteration
        self.time = full.time
        self.figtitle3d = f"MT {self.type}: {self.zlimits['bottom']} " \
                          f"to {self.zlimits['top']}" + \
                          self.len_units + ': \n' + self.signature

    def _extract_from_full(
            self,
            pos: list[np.ndarray]
    ) -> None:
        """Extract a subvolume from a full system representation.

        Because the process may involve splitting of the original
        filaments into apparently independent segments or exclusion
        of original filaments, 'self.fi' and 'self.ni' map filament
        and node indexes respectively to the original ones.
        Positions of filaments belonging to the new system are stored
        in 'self.pos'.

        :param pos: list of positions of original (complete) filaments
            in 3d space.
        """

        assert len(pos)
        self.pos = []
        self.fi = []
        self.ni = []
        sfi = -1
        sni = int
        for p in pos:
            # Assume that previous node is not valid.
            # Hence, start new filament slice.
            pnv = False
            for j in range(p.shape[0]):
                if self.zlimits['bottom'] < p[j, 2] < self.zlimits['top']:
                    if not pnv:  # start new filament slice
                        sfi += 1
                        sni = 1
                        self.pos.append([])
                        self.fi.append([])
                        self.ni.append([])
                        pnv = True
                    self.pos[sfi].append(p[j, :])  # node positions
                    self.fi[sfi].append(sfi)  # filament slice indexes
                    self.ni[sfi].append(sni)  # node index
                    sni += 1
                else:
                    if pnv:
                        pnv = False

        self.pos = [np.array(p) for p in self.pos]
        self.nfilaments = len(self.pos)
        self.nnodes = np.array([p.shape[0] for p in self.pos])
        self.mtmass = np.sum(self.nnodes)

    def initialise(self) -> None:
        """Populate the attribute arrays.
        """

        self.set_lengths()
        self.set_center_distance_2d()
        self.set_curvatures()
        self.set_radial_dev()

    @classmethod
    def derive_from(
            cls,
            full: type[FullDepth],
            init: bool = True,
    ) -> type[Slice]:
        """Main method to create the class and all its instances.

        :param full: The original system.
        :param init: Initialise the attribute arrays (in some use cases
            this is not necessary).
        :return: This class along with all its instances.
        """

        cls.params = full.params
        cls.type = 'slice'

        if init:
            cls.all_items = [cls(f) for f in full.all_items]
            [s.initialise() for s in cls.all_items]

        return cls


ListOfSpatialSystems = Union[
    list[FullDepth],
    list[Slice],
]
