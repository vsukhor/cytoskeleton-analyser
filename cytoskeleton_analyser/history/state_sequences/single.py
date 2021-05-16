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

"""State a microtubule end adopts in the progression of dynamic
instability.
"""

from __future__ import annotations
from typing import Final, Optional, Sequence

import numpy as np

from ...history.filament_history import FilamentHistory
from ..state_types import StateTypes
from ._state_sequence import _StateSequence


class Single(_StateSequence):
    """A specific state of the microtubule end.

    Class encapsulates a state the microtubule end adopts in the
    progression of its lifetime. Forms part of analysis of the
    microtubule dynamic instability.
    """

    #: Sequence width in the event records.
    WIDTH: Final = 1

    #: Slice of vars(self).keys() for automatic update.
    #: Exclude 'fr', 'to', 'params', 'duration', 'elongation',
    #: 'reorientation', 'velocity'.
    DATA_INDS: Final = slice(2, -4, None)

    def __init__(
            self,
            s: Optional[str] = None,
    ):

        #: 'fr' and 'to' are states of the filament end dynamics before
        #: the sequence start and during lifetime of its only state
        #: respectively.
        self.fr, self.to = StateTypes.from_str(s) \
            if s is not None else (None,) * (self.WIDTH + 1)

        super().__init__()

    def set_from_records(
            self,
            records: Sequence[FilamentHistory],
            dist_lim: Optional[Sequence[float]] = None,
    ) -> Single:
        """Populate the class datasets from a history-wide records.

        Select records specific for this states of this class
        and append them to own data holders for further analysis.

        :param records: The input records to select from.
        :param dist_lim: Limiting values of event distance from cell
            center defining cell regions.
        :return: Instance of this class initialised.
        """

        if dist_lim is None:
            dist_lim = [0., np.Inf]

        for j, f in enumerate(records):
            isg = np.where(
                (f.state_fr[:-self.WIDTH] == self.fr) &
                (f.state_to[:-self.WIDTH] == self.to) &
                (f.dist0[:-self.WIDTH] >= dist_lim[0]) &
                (f.dist0[:-self.WIDTH] < dist_lim[1])
            )[0]
            if isg.size:
                self.append(f, isg, j)
        self.remove_instantaneous()

        return self

    def __or__(
            self,
            other: Single,
    ) -> Single:

        return self.stack_attributes(
            other,
            Single('o' +
                   StateTypes.name_short(self.to)
                   )
        )

    def name(self) -> str:

        return StateTypes.name_short(self.fr) + \
               StateTypes.name_short(self.to) \
               if self.fr >= 0 else \
               StateTypes.name(self.to)
