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

"""SqlAlchemy database models of MTOC configuration parameters.

Configuration parameters for Microtubule Organizing Centers (MTOCs) are
arranged in type-specific database tables.
"""

from sqlalchemy import Column
from sqlalchemy import Float
from sqlalchemy import Integer
from sqlalchemy.orm import relationship

from . import Base


class ConfigMtocInSpace(Base):
    """Database model of microtubules having free minus-end.
    """

    __tablename__ = 'config_mtoc_inspace'

    id = Column(Integer, primary_key=True)

    total_frac = Column(Float)
    minus_end_cargo_free = Column(Float)
    use_pbs = Column(Integer)
    par_int = Column(Integer)
    par_real = Column(Float)

    main_configs = relationship('Config',
                                backref='config_mtoc_inspace')


class ConfigMtocGolgi(Base):
    """Database model of MTOC settings for Golgi-anchored microtubules.
    """

    __tablename__ = 'config_mtoc_golgi'

    id = Column(Integer, primary_key=True)

    total_frac = Column(Float)
    minus_end_cargo_free = Column(Float)
    origin_0 = Column(Float)
    origin_1 = Column(Float)
    origin_2 = Column(Float)
    size_0 = Column(Float)
    size_1 = Column(Float)
    size_2 = Column(Float)
    polar_bias = Column(Integer)
    angular_spread = Column(Float)

    main_configs = relationship('Config',
                                backref='config_mtoc_golgi')


class ConfigMtocCentrosome(Base):
    """Database model of MTOC settings for centrosomal microtubules.
    """

    __tablename__ = 'config_mtoc_centrosome'

    id = Column(Integer, primary_key=True)

    total_frac = Column(Float)
    minus_end_cargo_free = Column(Float)
    origin_0 = Column(Float)
    origin_1 = Column(Float)
    origin_2 = Column(Float)
    size_0 = Column(Float)
    size_1 = Column(Float)
    size_2 = Column(Float)

    main_configs = relationship('Config',
                                backref='config_mtoc_centrosome')


class ConfigMtocNucleus(Base):
    """Database model of MTOC settings for Nucleus-anchored microtubules.
    """

    __tablename__ = 'config_mtoc_nucleus'

    id = Column(Integer, primary_key=True)

    total_frac = Column(Float)
    minus_end_cargo_free = Column(Float)
    polar_bias = Column(Integer)
    angular_spread = Column(Float)

    main_configs = relationship('Config',
                                backref='config_mtoc_nucleus')
