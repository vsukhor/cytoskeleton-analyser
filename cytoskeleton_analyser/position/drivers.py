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

from logging import Logger
from typing import Optional

from ..inout import Paths
from ..inout import set_logger
from .reporters import Features
from .spatial_systems import ListOfSpatialSystems
from .spatial_systems import FullDepth
from .spatial_systems import Slice


def process(
        new: bool,
        paths: Paths,
        rind: int,
        pars: dict,
        show: bool,
        features: Optional[list[str]] = None,
) -> Optional[dict]:
    """Main entry point of the analysis pipeline.

    Can be used for both de novo analysis and for retrieval of
    prepared reports.

    :param new: Boolean switch choosing new analysis vs the retrieval.
    :param paths: Paths for input/output operations.
    :param rind: Index of the simulation run to be processed.
    :param pars: General-purpose parameters.
    :param show: If True, display the generated plots interactively.
    :param features: Features to report. If is None, all are reported.
    :return: If used for retrieval of prepared reports, summary of
        major results, othervise None.
    """

    paths.data_out = paths.data_out / 'position'
    Paths.ensure(paths.data_out)

    logger = set_logger(f'position ind {rind}',
                        paths.data_out / 'processing.log')

    if new:
        return produce(paths, pars, rind, logger, show, features)
    else:
        return restore(paths, pars, rind, logger, show, features)


def produce(
        paths: Paths,
        pars: dict,
        rind: int,
        logger: Logger,
        show: bool,
        features: Optional[list[str]] = None,
) -> None:
    """Produce ab initio analysis of major geometric charactersistics.

    Is applicable to all types of reconstructed spatial systems.

    :param paths: Collection of file system paths for data in-out.
    :param pars: General-purpose parameters.
    :param rind: Index of the simulation run to be processed.
    :param logger: logging.Logger object
    :param show: If True, display the generated plots interactively.
    :param features: Features to report. If is None, all are reported.
    """

    # Create full-depth reconstruction.
    FullDepth.create(paths, pars, rind, logger)
    produce_for_type(FullDepth.all_items, show, features)

    if 'slice_limits' in pars:
        # Extract slices.
        Slice.derive_from(FullDepth)
        produce_for_type(Slice.all_items, show, features)


def produce_for_type(
        sp: ListOfSpatialSystems,
        show: bool,
        features: Optional[list[str]] = None,
) -> None:
    """ Produce ab initio analysis of major geometric charactersistics.

    Is applicable to a specific class of reconstructed spatial system.

    :param sp: List of instances of the reconstructed spatial systems.
    :param show: If True, display the generated plots interactively.
    :param features: Features to report. If is None, all are reported.
    """

    assert len(sp) == 1 or [isinstance(a, type(sp[0])) for a in sp[1:]]
    tp = type(sp[0])

    tp.logger.info('-- Processing ' + tp.type + ' --')
    tp.logger.info('')

    tp.snaps = [s.snap_ind for s in sp]

    if tp.type == 'full':

        sp[0].plot3d('', with_mesh=True, mesh_flattened=True,
                     axes_visible=False, export=True, show=show)
        sp[0].plot3d('by_age', with_mesh=True, mesh_flattened=True,
                     axes_visible=False, export=True, show=show)
        sp[0].plot3d('by_height', with_mesh=True, mesh_flattened=True,
                     axes_visible=False, export=True, show=show)

    if tp.type == 'slice':

        sp[0].plot3d('', with_mesh=True, mesh_flattened=True,
                     axes_visible=False, export=True, show=show)

    if features is None:
        features = Features.all

    for f in features:
        if Features.is_applicable(f, tp):
            Features.reporter(f)\
                    .create(tp, show)\
                    .pipeline(sp)


def restore(
        paths: Paths,
        pars: dict,
        rind: int,
        logger: Logger,
        show: bool,
        features: Optional[list[str]] = None,
) -> dict:
    """Import and display pre-assembled reports from a storage location.

    :param paths: Collection of file system paths for data in-out.
    :param pars: General-purpose parameters.
    :param rind: Index of the simulation run to be processed.
    :param logger: logging.Logger object
    :param show: If True, display the generated plots interactively.
    :param features: Features to report. If is None, all are reported.
    :return: Summary of major results for both the full and slice
        reconstructions.
    """

    # Full-depth reconstruction.
    FullDepth.create(paths, pars, rind, logger, init=False)
    res_fd = import_reports(FullDepth, show, features)

    if 'slice_limits' in pars:
        # Extract TIRF-like slices.
        Slice.derive_from(FullDepth, init=False)
        res_ts = import_reports(Slice, show, features)
    else:
        res_ts = None

#    emulate_microcsope_image(Slice.all_items[1].pos)

    return {'full': res_fd,
            'slice': res_ts}


def import_reports(
        tp: type[FullDepth],
        show: bool,
        features: Optional[list[str]] = None,
) -> dict:
    """System type-specific restoration of preassembled reports.

    :param tp: Type of the reconstructed cell system.
    :param show: If True, display the generated plots interactively.
    :param show: If True, display the generated plots interactively.
    :param features: Features to report. If is None, all are reported.
    :return: Summary of major results relevant for the system type.
    """

    tp.logger.info('-- Restoring ' + tp.type + ' --')
    tp.logger.info('')

    if features is None:
        features = Features.all

    return {
        f: Features.reporter(f)
                   .create(tp, show)
                   .restore(tp.paths.data_out)
        for f in features if Features.is_applicable(f, tp)
    }
