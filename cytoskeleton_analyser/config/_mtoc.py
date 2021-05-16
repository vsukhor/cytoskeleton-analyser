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

"""Configuration of icrotubule organizing centers.
"""
from __future__ import annotations

from ..inout import Reader


class Mtoc:
    """
    Configuration parameters of microtubule organizing centers (MTOCs).

    Creates and stores sets of configuration parameters
    for describtion of four major MTOC types
    implemented in raw microtubule reconstructions:
    - inspace: microtubule minus-ends are not attached to any organizing center
    - golgi: microtubules minus ends are attached inside cell Golgi apparatus
    - centrosome: microtubule minus ends are seeded in cell centrosome
    - nucleus: microtubules are attached on the surface of nuclear membrane
    """
    def __init__(self):

        self.inspace = {}
        self.golgi = {}
        self.centrosome = {}
        self.nucleus = {}

    def read(
            self,
            r: Reader):
        """ Initialise this instance using configuration Reader.
        """

        # Parameter names used in the simulation output file are given
        # in comments.
        self.centrosome['total_frac'] = float(r.value(2))  # totalFract_Centrosome
        self.centrosome['minus_end_cargo_free'] = float(r.value(2))  # minusEnd_cargoFree_Centrosome
        self.centrosome['origin'] = r.array_float(1)  # origin_Centrosome
        self.centrosome['size'] = r.array_float(1)  # size_Centrosome
        self.nucleus['total_frac'] = float(r.value(2))  # totalFract_Nucleus
        self.nucleus['minus_end_cargo_free'] = float(r.value(2))  # minusEnd_cargoFree_Nucleus
        self.nucleus['polar_bias'] = int(r.value(2))  # polarBias_Nucleus
        self.nucleus['angular_spread'] = float(r.value(2))  # angularSpread_Nucleus
        self.golgi['total_frac'] = float(r.value(2))  # totalFract_Golgi
        self.golgi['minus_end_cargo_free'] = float(r.value(2))  # minusEnd_cargoFree_Golgi
        self.golgi['origin'] = r.array_float(1)  # origin_Golgi
        self.golgi['size'] = r.array_float(1)  # size_Golgi
        self.golgi['polar_bias'] = int(r.value(2))  # polarBias_Golgi
        self.golgi['angular_spread'] = float(r.value(2))  # angularSpread_Golgi
        self.inspace['total_frac'] = float(r.value(2))  # totalFract_InSpace
        self.inspace['minus_end_cargo_free'] = float(r.value(2))  # minusEnd_cargoFree_InSpace
        self.inspace['use_pbs'] = int(r.value(2))  # usePBC_InSpace
        self.inspace['par_int'] = int(r.value(2))  # parInt_InSpace
        self.inspace['par_real'] = float(r.value(2))  # parReal_InSpace

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
    ) -> Mtoc:
        """Initialise this instance from a python dict.
        """

        for k, v in d.items():
            setattr(self, k, v)

        return self
