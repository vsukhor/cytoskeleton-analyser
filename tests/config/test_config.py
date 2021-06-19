import numpy as np
import pytest

from pathlib import Path

from cytoskeleton_analyser.config.config import Config
from cytoskeleton_analyser.inout import Reader


class TestConfig:

    data_path = Path(__file__).parent.parent/'data'/'config_data.txt'

    @pytest.fixture()
    def init(self):
        return Config().read(Reader(self.data_path))

    def test_read(self, init):
        cf = init
        assert cf.plasma['type'] == 'Protr_test'
        assert cf.export_pdb['scaling'] == pytest.approx(0.1)
        assert cf.cytoskeleton['number'] == 1
        assert cf.cytoskeleton['run'] == -1
        assert cf.cytoskeleton['generateNew'] == 1
        assert cf.plasma['use'] == 1
        assert cf.plasma['run'] == 3021
        assert cf.nucleus['use'] == 1
        assert cf.run['index'] == 428
        assert cf.run['seed'] == 797746458
        assert all([a == b for a, b in zip(cf.nucleus['origin'],
                                           [0., 0., 6.])])
        assert all([a == b for a, b in zip(cf.nucleus['size'],
                                           [6., 6., 5.])])
        assert all([a == 0. for a in cf.nucleus['orientation']])
        assert cf.cytosol['volume'] == pytest.approx(3267.27)
        assert cf.cytosol['conc_coef'] == pytest.approx(1.9676e+06)
        assert cf.num_filaments == 500
        assert cf.cytosol['tubulin_total'] == pytest.approx(18.)
        assert cf.time_total == pytest.approx(1000000.)
        assert cf.save['fine'] == 1
        assert cf.save['coarse'] == 1
        assert cf.save['recordMultipleSaveFrames'] == 0
        assert cf.save['interval'] == pytest.approx(1000.)
        assert cf.export_pdb['fine'] == 1
        assert cf.export_pdb['coarse'] == 1
        assert cf.export_pdb['record_multiple_frames'] == 0
        assert cf.export_pdb['interval'] == pytest.approx(1000.)
        assert cf.history['record'] == 0
        assert cf.history['size'] == pytest.approx(100.)
        assert cf.proximity['threshold'] == pytest.approx(0.3)
        assert cf.proximity['omit_own_nodes'] == 1
        assert cf.proximity['step'] == 1
        assert cf.filament['step'] == pytest.approx(0.2)
        assert cf.filament['persist_length'] == pytest.approx(20.)

