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

"""Base classes for setting containers.

The containers are responsible for interaction with the database models.
"""
from __future__ import annotations
from typing import Any

from sqlalchemy import and_
from sqlalchemy.orm.session import Session


class Base:
    """Base class for implementing configuration setting containers.

    It should be subclassed to implement cell-component specific
    interactors with the database tables.
    """

    model: Any

    @classmethod
    def columns(
            cls,
            d: dict,
    ):
        """Database columns corresponding to keys of dictionary ``d``.
        """

        return and_(*(getattr(cls.model, k) == v for k, v in d.items()))

    @classmethod
    def assign(
            cls,
            d: dict,
            row: Base,
    ):

        for k, v in d.items():
            setattr(row, k, v)

    @classmethod
    def existing_db_rows(
            cls,
            d,
            session: Session,
    ) -> list:
        """SQL query to retrieve rows corresponding to dictionary ``d``.
        """

        return session \
            .query(cls.model) \
            .filter(cls.columns(d)) \
            .all()

    @classmethod
    def add_to_db(
            cls,
            d: dict,
            session: Session,
    ) -> Base:
        """Append dict items of dictionary ``d`` to the database.
        """

        row = cls.existing_db_rows(d, session)
        assert len(row) < 2
        if len(row):
            row = row[0]
        else:
            row = cls.model()
            for k, v in d.items():
                setattr(row, k, v)
            session.add(row)
            session.commit()

        return row


class Optional(Base):

    """Intermediate implementation for setting containers.

    Tailored to optional settings.
    """
    is_used: bool

    @classmethod
    def columns(
            cls,
            d: dict,
    ):
        """Database columns corresponding to keys of dictionary ``d``.
        """

        return \
            super().columns(d) \
            if cls.is_used \
            else False

    @classmethod
    def add_to_db(
            cls,
            d: dict,
            session: Session,
    ) -> Base:

        """Append dict items of dictionary ``d`` to the database.
        """

        return \
            super().add_to_db(d, session) \
            if cls.is_used \
            else None
