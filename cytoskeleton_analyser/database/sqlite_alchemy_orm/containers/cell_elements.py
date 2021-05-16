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

"""Subclasses of container base class for tailored to cell components.
"""

from collections import namedtuple

from . import container
from .. import models


class MembraneNucleus(container.Optional):
    """Container class for configuration of cell nuclear membrane.
    """

    model = models.ConfigNucleus


class MembranePlasma(container.Optional):
    """Container class for configuration of cell plasma membrane.
    """

    model = models.ConfigPlasma


class InSpace(container.Base):
    """Container class for configuration of unanchored microtubule MTOC.
    """

    model = models.ConfigMtocInSpace


class Golgi(container.Base):
    """Container class for configuration of Golgi-type MTOC.
    """

    model = models.ConfigMtocGolgi


class Centrosome(container.Base):
    """Container class for configuration of centrosome-type MTOC.
    """
    model = models.ConfigMtocCentrosome


class Nucleus(container.Base):
    """Container class for configuration of Nucleus-type MTOC.
    """

    model = models.ConfigMtocNucleus

#: Types of Microtubule Organizing Centers (MTOCs).
Mtoc = namedtuple('Mtoc', 'InSpace Golgi Centrosome Nucleus')

# Microtubule Organizing Centers (MTOCs).
mtoc = Mtoc(InSpace, Golgi, Centrosome, Nucleus)
