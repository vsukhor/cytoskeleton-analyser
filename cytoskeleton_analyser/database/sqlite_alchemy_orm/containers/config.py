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

"""Exposes Config class: database interactor for the configuration.
"""

from __future__ import annotations
from typing import Optional

from sqlalchemy.orm.session import Session

from . import container
from .. import models
from .cell_elements import MembranePlasma
from .cell_elements import MembraneNucleus
from .cell_elements import mtoc
from cytoskeleton_analyser.config.config import Config as ConfigData
from cytoskeleton_analyser.inout import flatten_dict


class Config(container.Base):
    """Database interactor for all simulation configuration parameters.
    """

    model = models.Config

    @staticmethod
    def set_db_fields(
            on: ConfigData,
    ) -> dict:
        """Set database fields to ConfigData ``òn``.
        """

        # These attributes are treated separately.
        omit = ['mtoc', 'm1d', 'plasma', 'nucleus']
        res = flatten_dict(on.as_dict(), exc=omit)
        core = [
            # Cytosol:
            'cytosol_volume',
            'cytosol_conc_coef',
            'cytosol_tubulin_total',
            # Proximity:
            'proximity_threshold',
            'proximity_omit_own_nodes',
            'proximity_step',
            # Filament
            'filament_persist_length',
            'filament_step',
            # Misc
            'num_filaments',
            'time_total'
        ]
        res['core'] = {k: v for k, v in res.items() if k in core}
        res['mtoc'] = {k: flatten_dict(v) for k, v in res['mtoc'].items()}
        res['m1d'] = flatten_dict(res['m1d'])

        res['plasma'] = flatten_dict(res['plasma'])
        MembranePlasma.is_used = res['plasma']['use'] > 0
        del res['plasma']['use']

        res['nucleus'] = flatten_dict(res['nucleus'])
        MembraneNucleus.is_used = res['nucleus']['use'] > 0
        del res['nucleus']['use']

        return res

    @classmethod
    def existing_db_rows(
            cls,
            dbf: dict,
            session: Session,
    ) -> list[models.Config]:
        """Query the database to get rows corresponding to dict ``dbf``.
        """

        return session \
            .query(cls.model) \
            .join(models.ConfigNucleus) \
            .join(models.ConfigPlasma) \
            .join(models.ConfigMtocInSpace) \
            .join(models.ConfigMtocGolgi) \
            .join(models.ConfigMtocCentrosome) \
            .join(models.ConfigMtocNucleus) \
            .filter(Config.columns(dbf['core'])) \
            .filter(Config.columns(dbf['m1d'])) \
            .filter(MembranePlasma.columns(dbf['plasma'])) \
            .filter(MembraneNucleus.columns(dbf['nucleus'])) \
            .filter(mtoc.InSpace.columns(dbf['mtoc']['inspace'])) \
            .filter(mtoc.Golgi.columns(dbf['mtoc']['golgi'])) \
            .filter(mtoc.Centrosome.columns(dbf['mtoc']['centrosome'])) \
            .filter(mtoc.Nucleus.columns(dbf['mtoc']['nucleus'])) \
            .all()

    @classmethod
    def is_in_db(
            cls,
            cfs: ConfigData,
            session: Session,
    ) -> tuple[bool, Optional[models.Config]]:
        """Check if raw corresponding to ``cfs`` is in the database.

        :param cfs: Data structure to query about.
        :param session: Database session.
        :return: tuple: If found: True and the raw found, otherwise
            False and None.
        """

        dbf = cls.set_db_fields(cfs)
        row = cls.existing_db_rows(dbf, session)
        assert len(row) < 2
        if len(row):
            return True, row[0]
        return False, None

    @classmethod
    def add_to_db(
            cls,
            on: ConfigData,
            session: Session,
    ) -> models.Config:
        """Append data from ConfigData instance ``òn`` to the database.
        """

        dbf = cls.set_db_fields(on)
        plasma = \
            MembranePlasma.add_to_db(dbf['plasma'], session)
        nucleus = \
            MembraneNucleus.add_to_db(dbf['nucleus'], session)
        mtoc_inspace = \
            mtoc.InSpace.add_to_db(dbf['mtoc']['inspace'], session)
        mtoc_golgi = \
            mtoc.Golgi.add_to_db(dbf['mtoc']['golgi'], session)
        mtoc_centrosome = \
            mtoc.Centrosome.add_to_db(dbf['mtoc']['centrosome'], session)
        mtoc_nucleus = \
            mtoc.Nucleus.add_to_db(dbf['mtoc']['nucleus'], session)

        row = cls.existing_db_rows(dbf, session)
        assert len(row) < 2
        if len(row):
            row = row[0]
        else:
            row = cls.model(
                mtoc_inspace_id=mtoc_inspace.id,
                mtoc_inspace=mtoc_inspace,
                mtoc_golgi_id=mtoc_golgi.id,
                mtoc_golgi=mtoc_golgi,
                mtoc_centrosome_id=mtoc_centrosome.id,
                mtoc_centrosome=mtoc_centrosome,
                mtoc_nucleus_id=mtoc_nucleus.id,
                mtoc_nucleus=mtoc_nucleus,
                run_inds='',    # will be set by append_run_ind later on
            )
            cls.assign(dbf['core'], row)
            cls.assign(dbf['m1d'], row)
            if plasma is not None:
                row.plasma_id = plasma.id
                row.plasma = plasma
            if nucleus is not None:
                row.nucleus_id = nucleus.id
                row.nucleus = nucleus

            session.add(row)
            session.commit()
        return row

    @staticmethod
    def append_run(
            rind: int,
            row: Config,
            session: Session,
    ) -> None:
        """Appends settings specific simulation run index to the db.

        :param rind: SImulation run index.
        :param row: Formatted database raw.
        :param session: Database session.
        """

        u = sorted([int(x) for x in row.run_inds.split()] + [rind])
        row.run_inds = ' '.join(str(x) for x in u)
        session.commit()

    @classmethod
    def update(
            cls,
            cf: ConfigData,
            session: Session,
    ) -> None:
        """Updates the database with the data from ConfigData object.
        """

        is_included, row = cls.is_in_db(cf, session)
        if not is_included:
            row = cls.add_to_db(cf, session)

        if cf.run['index'] not in [int(x) for x in row.run_inds.split()]:
            cls.append_run(cf.run['index'], row, session)
