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

import logging
from pathlib import Path
from typing import Sequence

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from .containers.config import Config as ConfigContainer
from .models import Base
from cytoskeleton_analyser.cells import CellType
from cytoskeleton_analyser.inout import Paths
from cytoskeleton_analyser.config import Config as ConfigRecord


class Db:
    """Wrapper class for a database engine and session.

    Sets a session with newly created or existing database.
    If the database with the specified ``name`` does not exist in
    the directory ``path``, a new database is created.
    """

    #: Database-specific logger.
    logger: logging.Logger

    def __init__(
            self,
            path: Path,
            name: str,
    ):
        """
        :param path: Directory for the database file.
        :param name: Name of the database file
        """

        # Valid SQLite URL forms are:
        #   sqlite:///:memory: (or, sqlite://),
        #   sqlite:///relative/path/to/file.db,
        #   sqlite:////absolute/path/to/file.db
        df_fname = 'sqlite:///' + str(path) + '/' + name + '.sqlite3'
        engine = create_engine(df_fname, echo=True)

        Base.metadata.create_all(engine)

        session = sessionmaker(bind=engine)

        #: Database session.
        self.session = session()

        self.set_logger(path)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.session.close()

    @staticmethod
    def set_logger(
            path: Path,
    ) -> None:
        """Sets logger to redirect database log to a file.

        :param path: Path to the database-specific log file.
        """

        file_h = logging.FileHandler(path / 'db.log')
        file_h.setLevel(logging.INFO)

        lgr = logging.getLogger('sqlalchemy')
        lgr.addHandler(file_h)
        lgr.setLevel(logging.INFO)


def update(
        cfs: Sequence[ConfigRecord],
        cell: CellType,
        path: Path,
) -> None:
    """Update/upgrade the databbase with new data records.

    If the database does not exist, it is created in the directory
    ``path``.

    :param cfs: Configuration data sets (one per simulation run) to
        insert into the database.
    :param cell: Cell data.
    :param path: File system path to the database.
    """

    Paths.ensure(path)

    with Db(path, cell.typename) as db:
        for cf in cfs:
            ConfigContainer.update(cf, db.session)
