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

"""
Full set of configuration parameters used in microtubule simulations.
"""
from __future__ import annotations
import json
from logging import Logger
from pathlib import Path
from typing import Final

from ..inout import Reader
from ._m1d import M1d
from ._mtoc import Mtoc


class Config:
    """Container for full set of configuration parameters.

    Includes all parameter sets sufficient to define a simulation run.
    """

    logger: Logger

    #: slice to omit self.mtoc and self.m1d, treated separately.
    DATA_INDS: Final = slice(0, -2, None)

    def __init__(self):

        #: Simulation run settings:
        #: 'index', 'seed'.
        self.run: dict = {}

        #: Nicleus settings:
        #: 'use', 'origin', 'size', 'orientation'.
        self.nucleus: dict = {}

        #: Plasma settings:
        #: 'use', 'type', 'run' index.
        self.plasma: dict = {}

        #: Cytoskeleton settings:
        #: 'number', 'run', 'generateNew'
        self.cytoskeleton: dict = {}

        #: Cytosol settings:
        #: 'volume', 'conc_coef', 'tubulin_total'
        self.cytosol: dict = {}

        #: Number of filaments seeded at simulation start.
        self.num_filaments: int = 0

        #: Total simulation time.
        self.time_total: float = 0.

        #: Parameters governing simulation data binary output:
        #: granularity: 'fine'  'coarse'
        #: 'record_multiple_frames'
        #: time between consecutive recordings: 'interval'
        self.save: dict = {}

        #: Parameters for output in pdb format:
        #: 'fine', 'coarse', 'record_multiple_frames',
        # 'interval', 'scaling'
        self.export_pdb: dict = {}

        #: record, size
        self.history: dict = {}

        #: threshold, omit_own_nodes, step
        self.proximity: dict = {}

        #: Internal parameters of the microtubule fiber:
        #: 'persist_length', 'step'
        self.filament: dict = {}

        #: Microtubule organizing center:
        self.mtoc: Mtoc = Mtoc()

        #: Parameters of dynamic instability:
        self.m1d: M1d = M1d()

    def read(
            self,
            r: Reader
    ) -> Config:
        """Read the project-specific configuration data.

        :param r: Instance of Reader to perform the task.
        """

        r.skip_lines(2)            # time, ''

        self.plasma['type'] = r.value(2).split(sep='/')[-2]  # workingDirIn

        r.skip_lines(3)            # workingDirOut, runIni, runEnd

        self.export_pdb['scaling'] = float(r.value(2))  # pdbScaling
        self.cytoskeleton['number'] = int(r.value(2))   # numCytoskeletons
        self.cytoskeleton['run'] = int(r.value(2))      # cytoskeletonRun
        self.cytoskeleton['generateNew'] = int(r.value(2))  # cskGenerateNew
        self.plasma['use'] = int(r.value(2))            # usePlasma
        self.plasma['run'] = int(r.value(2))            # plasmaMembraneRun
        self.nucleus['use'] = int(r.value(2))           # useNucleus

        r.skip_lines(1)                                 # Reading 'seeds'

        self.run['index'] = int(r.value(2))             # RUN
        self.run['seed'] = int(r.value(2))              # SEED

        r.skip_lines(4)    # '', 'Creating nucleus:', '', 'Reading configNuc'

        self.nucleus['origin'] = r.array_float(1)       # nucleus.origin
        self.nucleus['size'] = r.array_float(1)         # nucleus.size
        self.nucleus['orientation'] = r.array_float(1)  # nucleus.orientation

        r.skip_lines(4)

        # self.plasma['mesh_file'] = '/'.join(r.value(5).split(sep='/')[-4:])
        r.value(5)

        last_line = r.skip_lines(9)

        if last_line.split(' ')[0] == 'Unable':
            r.skip_lines(5)

        r.skip_lines(2)

        self.cytosol['volume'] = float(r.value(3))  # cytosol_volume, μm^3
        self.cytosol['conc_coef'] = float(r.value(3))  # cytosol_conc_coef

        r.skip_lines(2)        # '', Reading configCsk

        self.num_filaments = int(r.value(2))             # numFilaments
        # cTubulinTotal
        self.cytosol['tubulin_total'] = float(r.value(2))
        self.time_total = float(r.value(2))              # timeTotal

        r.skip_lines(1)                                  # logFrequency

        self.save['fine'] = int(r.value(2))              # saveFine
        self.save['coarse'] = int(r.value(2))            # saveCoarse
        # recordMultipleSaveFrames
        self.save['recordMultipleSaveFrames'] = int(r.value(2))
        self.save['interval'] = float(r.value(2))        # saveInterval

        self.export_pdb['fine'] = int(r.value(2))        # exportFine
        self.export_pdb['coarse'] = int(r.value(2))      # exportCoarse
        # recordMultiplePdbFrames
        self.export_pdb['record_multiple_frames'] = int(r.value(2))
        self.export_pdb['interval'] = float(r.value(2))  # pdbInterval

        self.history['record'] = int(r.value(2))     # recordHistory
        self.history['size'] = float(r.value(2))     # historySize

        r.skip_lines(1)                              # numProxThreads
        self.proximity['threshold'] = float(r.value(2))   # proximityThreshold
        self.proximity['omit_own_nodes'] = int(r.value(2))  # omitProximOwnNode
        self.proximity['step'] = int(r.value(2))     # proximityStep

        r.skip_lines(4)     # '', Reading configFTp, '', Reading configCskM1d
        self.filament['step'] = float(r.value(2))            # step
        self.filament['persist_length'] = float(r.value(2))  # persistLength

        r.skip_lines(5)
        self.mtoc.read(r)

        r.skip_lines(9)
        self.m1d.read(r)

        return self

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

        s = {k: getattr(self, k)
             for k in list(vars(self).keys())[self.DATA_INDS]}
        s['mtoc'] = self.mtoc.as_dict()
        s['m1d'] = self.m1d.as_dict()

        return s

    def from_dict(
            self,
            d: dict
    ) -> None:
        """Initialise this instance from a python dict.
        """

        for h in list(vars(self).keys())[self.DATA_INDS]:
            setattr(self, h, d[h])
        self.mtoc = Mtoc().from_dict(d['mtoc'])
        self.m1d = M1d().from_dict(d['m1d'])

    def to_json(
            self,
            path: Path
    ) -> None:
        """Serialise content of this instance to 'config.json' file.
        """

        s = self.as_dict()

        fname = path / 'config.json'
        self.logger.info("Dumping to " + str(fname))
        with open(fname, 'w') as f:
            json.dump(s, f)

    def from_json(
            self,
            path: Path
    ) -> Config:
        """Initialise this instance from 'config.json'.
        """

        fname = path / 'config.json'
        self.logger.info("Importing from " + str(fname))
        with open(fname, 'r') as f:
            d = json.load(f)
            self.from_dict(d)

        return self
