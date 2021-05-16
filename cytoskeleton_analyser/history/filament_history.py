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

"""Filament-specific histories of microtubule state changes.

These are characteristics of reconstructed microtubules recorded for a
particular microtubule end. Filament-specific histories (encapsilated by
'HistoryPerFilament' class) used for state type-specific analysis by a
'states.States' object.
"""

from __future__ import annotations
from typing import BinaryIO

import numpy as np

from ..inout import read_to_ndarray


class FilamentHistory:
    """Encapsulates a time-ordered sequence of events.

    The events must have occured at specific end of the same filament.
    Type of an event is formed by a pair (state_from, state_to) of
    filament end states before and after the event respectively.
    Besides, the recording contains other event characteristics, such as
    time and the filament end spatial position.
    """

    def __init__(
            self,
            numrec: int,
            version: int
    ):
        """
        :param numrec: Nnumber of event records.
        :param version: Version of the record format.
        """

        #: Number of event records forming the sequence.
        self.numrec = numrec

        #: Time of the event.
        self.time = np.empty(numrec, dtype=np.float64)

        #: State of the filament end before the event.
        self.state_fr = np.empty(numrec, dtype=np.uint32)

        #: State of the filament end after the event.
        self.state_to = np.empty(numrec, dtype=np.uint32)

        #: Spatial position of the end node.
        self.pos = np.empty([numrec, 3], dtype=np.float32)

        if version == 2:

            #: End node direction (unit vector).
            self.ornt = np.empty([numrec, 3], dtype=np.float32)

        #: Filament length.
        self.length = np.empty(numrec, dtype=np.uint32)

        #: Filament age.
        self.age = np.empty(numrec, dtype=np.float64)

        if version == 2:

            #: Number of growth events.
            self.ngrw = np.empty(numrec, dtype=np.uint64)

            #: Number of shrink events.
            self.nshr = np.empty(numrec, dtype=np.uint64)

        #: cas field intensities.
        self.cas = np.empty([numrec, 6], dtype=np.float32)

        #: Distance to the plasma membrane.
        self.dist_plm = np.empty(numrec, dtype=np.float32)

        #: Distance to cell nuclear membrane.
        self.dist_nuc = np.empty(numrec, dtype=np.float32)

        #: Distance to cell center.
        self.dist0 = np.empty(numrec, dtype=np.float32)

    def read(
            self,
            f: BinaryIO,
    ) -> FilamentHistory:
        """Set the class fields by reading them from the file ``f``.
        """

        for i in range(self.numrec):
            # By taking [1:-1], 'numrec' and 'dist0' are excluded.
            for k in list(vars(self).keys())[1:-1]:
                read_to_ndarray(getattr(self, k), f, i)
        self.dist0 = np.linalg.norm(self.pos[:, :2], axis=1)

        return self
