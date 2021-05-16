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

"""Types of dynamic states the microtubule ends adopt.

Enumeration of state types of a microtubule in the coarse of
its simulated lifetime.
"""

from __future__ import annotations
from typing import Final


class StateTypes:
    """Enumeration of kinetic states a microtubule end can adopt.
    """

    O: Final[int] = -1  # undefined
    G: Final[int] = 0   # growing
    S: Final[int] = 1   # shrinking
    P: Final[int] = 2   # paused
    C: Final[int] = 3   # connected
    D: Final[int] = 4   # depolymerized

    @classmethod
    def from_str(
            cls,
            ss: str
    ) -> (int, int):

        return (cls.__dict__[s.upper()] for s in ss)

    @classmethod
    def name(
            cls,
            k: int
    ):

        return \
            'growth' if k == cls.G else \
            'shrink' if k == cls.S else \
            'pause' if k == cls.P else \
            'connect' if k == cls.C else \
            'depol' if k == cls.D else \
            'undef'

    @classmethod
    def name_short(
            cls,
            k: int
    ) -> str:

        return \
            'g' if k == cls.G else \
            's' if k == cls.S else \
            'p' if k == cls.P else \
            'c' if k == cls.C else \
            'd' if k == cls.D else \
            'o'
