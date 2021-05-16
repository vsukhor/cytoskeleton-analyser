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

""" Cell and microtubule organizing center configurations.

For each cell type class, the following is defined:
'typename' and 'plmind' - colloquial name and index of specific
membrane configuration used in the simulation respectively.
These form also a part of file path to access and store the cell
the cell data.

"""

from typing import Final, NamedTuple

import numpy as np

#: Configuration names of the microtubule organizing celters.
organizing_centers: Final[list[str]] = [
    'total',
    'anchored at centrosome',
    'anchored at nucleus',
    'anchored at golgi',
    'MTOC-free'
]


class Region:
    """Specific cell subcompartment.
    """
    def __init__(
            self,
            rad_min: float,
            rad_max: float,
    ):
        """
        :param rad_min: Min radial distance from cell center in xy plane.
        :param rad_max: Max radial distance from cell center in xy plane.
        """

        self.rad_min = rad_min
        self.rad_max = rad_max

    def is_inside(
            self,
            c: np.ndarray,
    ) -> np.ndarray:
        """
        :param c: Radial distance from cell center in xy plane.
        :return: True iff point at 'c' is within the region borders.
        """

        return np.where((self.rad_min <= c) & (c < self.rad_max))[0]


class Regions(NamedTuple):
    """Definition of cell regions.
    """

    cytosol: Region
    soma: Region
    lamella: Region
    edge: Region


class CellType(NamedTuple):
    """Prototype class for specifying cells of specific types.

    'typename' and 'plmind' attributed form also a part of file path to
    access and store the cell the cell data.
    """

    typename: str        #: Type name.
    plmind: int          #: Index of the plasma membrane mesh.
    regions: Regions     #: Subcellular compartments.
