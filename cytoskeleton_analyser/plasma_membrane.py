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

""" Plasma membrane class. Encapsulates a 3d mesh of cell boundary.

It defines volume available for the microtubules.
"""

from pathlib import Path

import meshio
import numpy as np
from scipy.spatial.distance import cdist

from .cells import CellType


class PlasmaMembrane:
    """A minimalistic cell plasma membrane.

    I is used mostly as a 3d cell boundary to limit the volume
    available for microtubules.
    """

    def __init__(
            self,
            path: Path,
            cell: CellType,
            origin: np.ndarray = np.zeros(3)
    ):
        """
        :param path: Path to the mesh file in .stl format.
        :param cell: Cell type.
        :param origin: Cell geometric origin point.
        """

        #: meshio.Mesh object containing mesh representing the membrane.
        self.mesh: meshio.Mesh = self.load(path, cell)

        #: Minimal position of mesh points.
        self.min_ = self.mesh.points.min(0)

        #: Maxiimal position of mesh points.
        self.max_ = self.mesh.points.max(0)

        #: Point of ell geometric origin.
        self.origin: np.ndarray = origin

    @staticmethod
    def load(
            path: Path,
            cell: CellType,
    ) -> meshio.Mesh:
        """Read in cell membrane from file into a meshio.Mesh object.

        :param path: Path to the mesh file in .stl format.
        :param cell: Cell type.
        :return: Initialized mesh.
        """

        fname = path / f"plasmaMesh_{cell.plmind}.stl"
        return meshio.read(fname)

    def radial_extent(self) -> float:
        """Max extents of the membrane mesh in xy plane.

        :return: Distance to the furthest mesh node in xy plane.
        """

        return max(cdist(self.mesh.points[:, :2],
                         np.array([self.origin[:2]])).T[0])
