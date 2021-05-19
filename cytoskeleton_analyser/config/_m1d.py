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

"""Specifics of m1d microtubue model encapsulated in M1d class.

This class is part of simulation configuration parameters.
"""
from __future__ import annotations

from ..inout import Reader


class M1d:
    """Parameters of m1d microtubule model.

    The constitutive part of simulation configuration parameters.
    """

    def __init__(self):

        # Names of attributes and dictionary keys correspond to names of
        # simulation configuraion parameters.
        self.velocity = {}      # sh0, sh1
        self.state_change = {}  # gs, sg, ps
        self.release = 0.
        self.cut = 0.
        self.cap = {}         # on, off
        self.pl_membr = {}    # orient_strength, interact_thresh
        self.nuc_membr = {}   # orient_strength, interact_thresh
        # thickness, orientation, growth_slowdown,
        # catastr_inhibition, catastr_activation,
        # rescue_inhibition, rescue_activation
        self.cas = {}

    def read(
            self,
            r: Reader
    ):
        """ Initialise this instance using the configuration Reader.
        """

        self.release = float(r.value(2))
        self.cut = float(r.value(2))
        self.cap['on'] = float(r.value(2))
        self.cap['off'] = float(r.value(2))

        self.state_change['gs'] = float(r.value(2))
        self.state_change['sg'] = float(r.value(2))
        self.state_change['ps'] = float(r.value(2))

        self.cas['thickness'] = float(r.value(2))
        self.cas['orientation'] = {}
        self.cas['orientation']['extent'] = float(r.value(2))
        self.cas['orientation']['intensity'] = float(r.value(2))
        self.cas['catastr_activation'] = {}
        self.cas['catastr_activation']['extent'] = float(r.value(2))
        self.cas['catastr_activation']['intensity'] = float(r.value(2))
        self.cas['catastr_inhibition'] = {}
        self.cas['catastr_inhibition']['extent'] = float(r.value(2))
        self.cas['catastr_inhibition']['intensity'] = float(r.value(2))
        self.cas['rescue_activation'] = {}
        self.cas['rescue_activation']['extent'] = float(r.value(2))
        self.cas['rescue_activation']['intensity'] = float(r.value(2))
        self.cas['rescue_inhibition'] = {}
        self.cas['rescue_inhibition']['extent'] = float(r.value(2))
        self.cas['rescue_inhibition']['intensity'] = float(r.value(2))
        self.cas['growth_slowdown'] = {}
        self.cas['growth_slowdown']['extent'] = float(r.value(2))
        self.cas['growth_slowdown']['intensity'] = float(r.value(2))

        self.pl_membr['orient_strength'] = float(r.value(2))
        self.pl_membr['interact_thresh'] = float(r.value(2))

        self.nuc_membr['orient_strength'] = float(r.value(2))
        self.nuc_membr['interact_thresh'] = float(r.value(2))

        self.velocity['sh0'] = float(r.value(2))
        self.velocity['sh1'] = float(r.value(2))

    def __eq__(self, other):
        """Equality operator.
        """

        if other.__class__ is not self.__class__:
            return NotImplemented

        return all([getattr(self, k) == getattr(other, k)
                    for k in list(vars(self).keys())])

    def as_dict(self) -> dict:
        """Translate this instance into a python dict.
        """

        return {k: getattr(self, k) for k in list(vars(self).keys())}

    def from_dict(
            self,
            d: dict
    ) -> M1d:
        """Initialise this instance from a python dict.
        """

        for k, v in d.items():
            setattr(self, k, v)

        return self
