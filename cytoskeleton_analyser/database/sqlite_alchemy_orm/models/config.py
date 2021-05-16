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

"""SqlAlchemy database model of full configuration parameter set.
"""

from sqlalchemy import Column
from sqlalchemy import Float
from sqlalchemy import ForeignKey
from sqlalchemy import Integer
from sqlalchemy import String
from sqlalchemy.orm import relationship

from . import Base


class Config(Base):

    """SqlAlchemy database model of full configuration parameter set.
    """

    __tablename__ = 'config'

    id = Column(Integer, primary_key=True)

    # Plasma:
    plasma_id = Column(Integer,
                       ForeignKey('config_plasma.id'))
    plasma = relationship('ConfigPlasma',
                          backref='config',
                          cascade='all, delete')

    # Nucleus:
    nucleus_id = Column(Integer,
                        ForeignKey('config_nucleus.id'))
    nucleus = relationship('ConfigNucleus',
                           backref='config',
                           cascade='all, delete')

    # Mtoc:
    mtoc_inspace_id = Column(Integer,
                             ForeignKey('config_mtoc_inspace.id'))
    mtoc_inspace = relationship('ConfigMtocInSpace',
                                backref='config')

    mtoc_golgi_id = Column(Integer,
                           ForeignKey('config_mtoc_golgi.id'))
    mtoc_golgi = relationship('ConfigMtocGolgi',
                              backref='config')

    mtoc_centrosome_id = Column(Integer,
                                ForeignKey('config_mtoc_centrosome.id'))
    mtoc_centrosome = relationship('ConfigMtocCentrosome',
                                   backref='config')

    mtoc_nucleus_id = Column(Integer,
                             ForeignKey('config_mtoc_nucleus.id'))
    mtoc_nucleus = relationship('ConfigMtocNucleus',
                                backref='config')

    # M1d:
    release = Column(Float)
    cut = Column(Float)
    cap_on = Column(Float)
    cap_off = Column(Float)
    state_change_gs = Column(Float)
    state_change_sg = Column(Float)
    state_change_ps = Column(Float)
    cas_thickness = Column(Float)
    cas_orientation_extent = Column(Float)
    cas_orientation_intensity = Column(Float)
    cas_growth_slowdown_extent = Column(Float)
    cas_growth_slowdown_intensity = Column(Float)
    cas_catastr_inhibition_extent = Column(Float)
    cas_catastr_inhibition_intensity = Column(Float)
    cas_catastr_activation_extent = Column(Float)
    cas_catastr_activation_intensity = Column(Float)
    cas_rescue_inhibition_extent = Column(Float)
    cas_rescue_inhibition_intensity = Column(Float)
    cas_rescue_activation_extent = Column(Float)
    cas_rescue_activation_intensity = Column(Float)
    pl_membr_orient_strength = Column(Float)
    pl_membr_interact_thresh = Column(Float)
    nuc_membr_orient_strength = Column(Float)
    nuc_membr_interact_thresh = Column(Float)
    velocity_sh0 = Column(Float)
    velocity_sh1 = Column(Float)

    # Cytosol:
    cytosol_volume = Column(Float)
    cytosol_conc_coef = Column(Float)
    cytosol_tubulin_total = Column(Float)

    # Proximity:
    proximity_threshold = Column(Float)
    proximity_omit_own_nodes = Column(Integer)
    proximity_step = Column(Integer)

    # Filament:
    filament_persist_length = Column(Float)
    filament_step = Column(Float)

    # Misc:
    num_filaments = Column(Integer)
    time_total = Column(Float)

    # Runs:
#    runs_id = Column(Integer, ForeignKey('runs.id'))
#    runs = relationship('Runs', backref='config')
#    history = relationship('History', backref='config')
    run_inds = Column(String(1024))
