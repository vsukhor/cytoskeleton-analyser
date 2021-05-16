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

"""A set of data models based on combinations of exponentials.
"""

from __future__ import annotations
from collections import namedtuple
from typing import Final, NamedTuple, Optional, Sequence

import numpy as np

from .base import Base


class _ExpBase:
    """Base class for exponential models.
    """

    outer = None

    @classmethod
    def create(cls, outer) -> type(_ExpBase):
        cls.outer = outer
        return cls


class _SaturatedGrowthHomogenous(_ExpBase):
    """Homogenous saturated growth.
    """

    sl: Final = slice(0, 2)

    @classmethod
    def f(cls, x, a, tau1):
        return a * (1. - np.exp(-x / tau1))

    @classmethod
    def saturation(cls, p):
        return p[0]


class _SaturatedGrowthInhomogenous(_ExpBase):
    """Inhomogenous saturated growth.
    """
    sl: Final = np.array([0, 1, 6])

    @classmethod
    def f(cls, x, a, tau1, d):
        return a * (1. - np.exp(-x / tau1)) + d

    @classmethod
    def saturation(cls, p):
        return p[0] + p[6]


class _DecayHomogenous(_ExpBase):
    """Inhomogenous decay.
    """

    sl: Final = slice(0, 2)

    @classmethod
    def f(cls, x, a, tau1):
        return a * np.exp(-x / tau1)

    @classmethod
    def saturation(cls, _):
        return 0.


class _DecayInhomogenous(_ExpBase):
    """Homogenous decay.
    """

    sl: Final = np.array([0, 1, 6])

    @classmethod
    def f(cls, x, a, tau1, d):
        return a * np.exp(-x / tau1) + d

    @classmethod
    def saturation(cls, p):
        return p[6]


class _DecayDecayHomogenous(_ExpBase):
    """Homogenous dual decay.
    """

    sl: Final = slice(0, 4)

    @classmethod
    def f(cls, x, a, tau1, b, tau2):
        return a * (b * np.exp(-x / tau1) +
                    (1. - b) * np.exp(-x / tau2))

    @classmethod
    def saturation(cls, _):
        return 0.


class _DecayDecayInhomogenous(_ExpBase):
    """Inhomogenous dual decay.
    """

    sl: Final = np.array([0, 1, 2, 3, 6])

    @classmethod
    def f(cls, x, a, tau1, b, tau2, d):
        return a * (b * np.exp(-x / tau1) +
                    (1. - b) * np.exp(-x / tau2)) + d

    @classmethod
    def saturation(cls, p):
        return p[6]


class _DecaySaturatedGrowthHomogenous(_ExpBase):
    """Homogenous sum of decay and saturated growth.
    """

    sl: Final = np.array([0, 1, 4, 5])

    @classmethod
    def f(cls, x, a, tau1, c, tau3):
        return a * np.exp(-x / tau1) + \
               c * (1. - np.exp(-x / tau3))

    @classmethod
    def saturation(cls, p):
        return p[4]


class _SaturatedGrowthSaturatedGrowthHomogenous(_ExpBase):
    """Homogenous dual saturated growth.
    """

    sl: Final = slice(0, 4)

    @classmethod
    def f(cls, x, a, tau1, b, tau2):
        return a * (1. - b * np.exp(-x / tau1) -
                    (1. - b) * np.exp(-x / tau2))

    @classmethod
    def saturation(cls, p):
        return p[0]


class _SaturatedGrowthSaturatedGrowthInhomogenous(_ExpBase):
    """Inhomogenous dual saturated growth.
    """

    sl: Final = np.array([0, 1, 2, 3, 6])

    @classmethod
    def f(cls, x, a, tau1, b, tau2, d):
        return a * (1. - b * np.exp(-x / tau1) -
                    (1. - b) * np.exp(-x / tau2)) + d

    @classmethod
    def saturation(cls, p):
        return p[0] + p[6]


class _SaturatedGrowthSaturatedGrowthDecayHomogenous(_ExpBase):
    """Homogenous sum of decay and dual saturated growth.
    """

    sl: Final = slice(0, 6)

    @classmethod
    def f(cls, x, a, tau1, b, tau2, c, tau3):
        return a * (1. - b * np.exp(-x / tau1) -
                         (1. - b) * np.exp(-x / tau2)) + \
               c * np.exp(-x / tau3)

    @classmethod
    def saturation(cls, p):
        return p[0]


class _DecayDecaySaturatedGrowthHomogenous(_ExpBase):
    """Homogenous sum of dual decay and saturated growth.
    """

    sl: Final = slice(0, 6)

    @classmethod
    def f(cls, x, a, tau1, b, tau2, c, tau3):
        return a * (b * np.exp(-x / tau1) +
                    (1. - b) * np.exp(-x / tau2)) + \
               c * (1. - np.exp(-x / tau3))

    @classmethod
    def saturation(cls, p):
        return p[4]


class _DecayDecaySaturatedGrowthInhomogenous(_ExpBase):
    """Inhomogenous sum of dual decay and saturated growth.
    """

    sl: Final = slice(0, 7)

    @classmethod
    def f(cls, x, a, tau1, b, tau2, c, tau3, d):
        return a * (b * np.exp(-x / tau1) +
                    (1. - b) * np.exp(-x / tau2)) + \
               c * (1. - np.exp(-x / tau3)) + d

    @classmethod
    def saturation(cls, p):
        return p[4] + p[6]


class _SaturatedGrowthLine(_ExpBase):
    """Sum of saturated growth and a line pi/4.
    """

    sl: Final = np.array([0, 1, 4, 6])

    @classmethod
    def f(cls, x, a, tau1, c, d):
        return a * (1. - np.exp(-x / tau1)) + c * (x - d)

    @classmethod
    def saturation(cls, _):
        return None


# ======================================================================
class Exponential(Base):
    """Data models using combinations of exponential functions.
    """

    class Pars(NamedTuple):
        a: float = np.nan
        tau1: float = np.nan
        b: float = np.nan
        tau2: float = np.nan
        c: float = np.nan
        tau3: float = np.nan
        d: float = np.nan

    #: Max number of fitted parameters.
    numpars = len(Pars._fields)

    sg_h: type(_SaturatedGrowthHomogenous)
    sg_i: type(_SaturatedGrowthInhomogenous)
    d_h: type(_DecayHomogenous)
    d_i: type(_DecayInhomogenous)
    d_d_h: type(_DecayDecayHomogenous)
    d_d_i: type(_DecayDecayInhomogenous)
    d_sg_h: type(_DecaySaturatedGrowthHomogenous)
    sg_sg_h: type(_SaturatedGrowthSaturatedGrowthHomogenous)
    sg_sg_i: type(_SaturatedGrowthSaturatedGrowthInhomogenous)
    sg_sg_d_h: type(_SaturatedGrowthSaturatedGrowthDecayHomogenous)
    d_d_sg_h: type(_DecayDecaySaturatedGrowthHomogenous)
    d_d_sg_i: type(_DecayDecaySaturatedGrowthInhomogenous)
    sg_line: type(_SaturatedGrowthLine)

    def __init__(
            self,
            x: Sequence,
            p0: Sequence,
            x_units: Optional[str] = None,
    ):
        """
        :param x: Data domain.
        :param p0: Initial parameter values.
        :param x_units: Units of x values.
        """

        super().__init__(x, list(p0), x_units)
        maxtau = 10. * np.max(x)
        self.bounds = (np.array([0., 0., 0., 0., 0., 0., 0.]),
                       np.array([+np.inf, maxtau, 1., maxtau,
                                 +np.inf, maxtau, +np.inf]))
        assert all(len(b) == self.numpars for b in self.bounds)
        self.saturation_ind: Optional[int] = None
        self.name: str = ''

    def fit(
            self,
            f: type(_ExpBase),
            y: np.ndarray,
    ) -> Exponential:
        """Fitting is done here.

        :param f: Function to fit.
        :param y: Data y values.
        :return: Instance of this class with the fit results set.
        """

        self.p = np.full_like(self.p0, np.nan)
        self.name = self.__class__.__name__ + '.' + f.__name__
        super()._fit(f.f, y, f.sl)
        self.par = {k: v for k, v in zip(self.Pars._fields, self.p)}
        self.aik = self._aik(y)
        self.saturation = f.saturation(self.p)
        self.set_equilibration(y)
        self.report()

        return self

    def restore(
            self,
            f: type(_ExpBase),
            y: np.ndarray,
    ) -> Exponential:
        """Restore the modelled function and fit quality indicators.

        Applies optimal parameters determined previously to restore
        fitted function and the fit attributes.

        :param f: Function to appy.
        :param y: Fitted data.
        :return: Instance of this class with the attributes
            reinitialised.
        """

        self.p = np.array(self.p0)
        self.name = self.__class__.__name__ + '.' + f.__name__
        self.isbest = True
        self.prediction = self.predict(f.f, f.sl)
        self.set_quality_indicators(y)         # residnorm, chi2 , fano
        self.par = {k: v for k, v in zip(self.Pars._fields, self.p)}
        self.aik = self._aik(y)
        self.saturation = f.saturation(self.p)
        self.set_equilibration(y)
        self.report()

        return self

    @classmethod
    def create(cls):
        """Create all models.
        """

        cls.sg_h = _SaturatedGrowthHomogenous.create(cls)
        cls.sg_i = _SaturatedGrowthInhomogenous.create(cls)
        cls.d_h = _DecayHomogenous.create(cls)
        cls.d_i = _DecayInhomogenous.create(cls)
        cls.d_d_h = _DecayDecayHomogenous.create(cls)
        cls.d_d_i = _DecayDecayInhomogenous.create(cls)
        cls.d_sg_h = _DecaySaturatedGrowthHomogenous.create(cls)
        cls.sg_sg_h = _SaturatedGrowthSaturatedGrowthHomogenous.create(cls)
        cls.sg_sg_i = _SaturatedGrowthSaturatedGrowthInhomogenous.create(cls)
        cls.sg_sg_d_h = _SaturatedGrowthSaturatedGrowthDecayHomogenous.\
            create(cls)
        cls.d_d_sg_h = _DecayDecaySaturatedGrowthHomogenous.create(cls)
        cls.d_d_sg_i = _DecayDecaySaturatedGrowthInhomogenous.create(cls)
        cls.sg_line = _SaturatedGrowthLine.create(cls)

        return cls

    @classmethod
    def subtype(cls, name: str):
        """Subclass object from ``name``.
        """

        for _, v in vars(cls)['__annotations__'].items():
            u = v[5:-1]
            if u == name:
                return globals()[u]
