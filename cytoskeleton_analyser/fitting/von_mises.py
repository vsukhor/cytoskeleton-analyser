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

"""Class modelling data using VonMises distribution.
"""

from __future__ import annotations
from typing import Callable, Final, Sequence

import numpy as np

from .base import Base


class VonMisesTriple(Base):
    """Data model using three-component von Mises distribution.
    """


    nump: Final = 6              #: Number of fitted parameters.
    sl: Final = slice(0, nump)   #: Slice selecting parameters to fit.

    def __init__(
            self,
            x: Sequence,
            p0: Sequence,
    ):
        """
        :param x: Data domain.
        :param p0: Initial parameter values.
        """

        assert len(p0) == self.nump
        super().__init__(x, p0)
        self._class = self.__class__
        self.bounds = (np.array([0., 0., 0., 0., 0., 0.]),
                       np.array([+np.inf, +np.inf, 1., 1., +np.inf, 2*np.pi]))

    def _par(self) -> dict:
        """Parameters as a dictionary {'name': value}.
        """

        return {'k1': self.p[0],
                'k2': self.p[1],
                'a': self.p[2],
                'b': self.p[3],
                'k3': self.p[4],
                'mi3': self.p[5]}

    def fit(
            self,
            f: Callable,
            y: np.ndarray,
    ) -> VonMisesTriple:
        """Fitting is done here.

        :param f: Function to fit.
        :param y: Data y values.
        :return: Instance of this class with the fit results set.
        """

        self.p = np.full_like(self.p0, np.nan)
        self.name = self.__class__.__name__ + '.' + f.__name__
        super()._fit(f, y, self.sl)
        self.par = self._par()
        self.set_cc(self.par)
        self.aik = self._aik(y)
        self.report()

        return self

    def restore(
            self,
            f: Callable,
            y: np.ndarray,
    ) -> VonMisesTriple:
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
        self.prediction = self.predict(f, self.sl)
        self.cc = []
        self.set_quality_indicators(y)         # residnorm, chi2 , fano
        self.par = self._par()
        self.set_cc(self.par)
        self.aik = self._aik(y)
        self.report()

        return self

    def set_cc(
            self,
            par: dict
    ) -> None:

        self.cc = [
            self.c1(self.x, par['k1'], par['a'], par['b']),
            self.c2(self.x, par['k2'], par['a'], par['b']),
            self.c3(self.x, par['a'], par['k3'], par['mi3']),
        ]


    @classmethod
    def c1(cls, th, k1, a, b):
        """First component of the fitted function.
        """

        return a * b * np.exp(k1 * np.cos(th)) / (2.*np.pi * np.i0(k1))

    @classmethod
    def c2(cls, th, k2, a, b):
        """Second component of the fitted function.
        """

        return a * (1. - b) * (np.exp(k2 * np.cos(th - np.pi/2.)) +
                               np.exp(k2 * np.cos(th - 3.*np.pi/2.))) \
               / (2.*np.pi * np.i0(k2))

    @classmethod
    def c3(cls, th, a, k3, mi3):
        """Third component of the fitted function.
        """

        return (1. - a) * (np.exp(k3 * np.cos(th - mi3)) +
                           np.exp(k3 * np.cos(th - (2.*np.pi - mi3)))) \
               / (2.*np.pi * np.i0(k3))

    @classmethod
    def full(cls, th, k1, k2, a, b, k3, mi3):
        """Sum of component function.
        """

        return cls.c1(th, k1, a, b) + \
               cls.c2(th, k2, a, b) + \
               cls.c3(th, a, k3, mi3)


class VonMisesDouble(Base):
    """Implements data model using two-component von Mises distribution.
    """

    nump: Final = 5              #: Number of fitted parameters.
    sl: Final = slice(0, nump)   #: Slice selecting parameters to fit.

    def __init__(
            self,
            x: Sequence,
            p0: Sequence,
    ):
        """
        :param x: Data domain.
        :param p0: Initial parameter values.
        """

        assert len(p0) == self.nump
        super().__init__(x, p0)
        self._class = self.__class__
        self.bounds = ([0., 0., 0., 0., 0.],
                       [+np.inf, 1., +np.inf, +np.inf, 2*np.pi])

    def _par(self) -> dict:
        """Parameters as a dictionary {'name': value}.
        """

        return {'a': self.p[0],
                'b': self.p[1],
                'k1': self.p[2],
                'k2': self.p[3],
                'mi2': self.p[4]}

    def set_cc(
            self,
            par: dict
    ) -> None:

        self.cc = [
            self.c1(self.x, par['a'], par['b'], par['k1']),
            self.c2(self.x, par['a'], par['b'], par['k2'], par['mi2']),
        ]

    def fit(
            self,
            f: Callable,
            y: np.ndarray,
    ) -> VonMisesDouble:
        """Fitting is done here.

        :param f: Function to fit.
        :param y: Data y values.
        :return: Instance of this class with the fit results set.
        """

        self.p = np.full_like(self.p0, np.nan)
        self.name = self.__class__.__name__ + '.' + f.__name__
        super()._fit(f, y, self.sl)
        self.par = self._par()
        self.set_cc(self.par)
        self.aik = self._aik(y)
        self.report()

        return self

    def restore(
            self,
            f: Callable,
            y: np.ndarray,
    ) -> VonMisesDouble:
        """Restore modelled function and fit quality indicators.

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
        self.prediction = self.predict(f, self.sl)
        self.set_quality_indicators(y)         # residnorm, chi2 , fano
        self.par = self._par()
        self.set_cc(self.par)
        self.aik = self._aik(y)
        self.report()

        return self

    @classmethod
    def c1(cls, th, a, b, k1):
        """First component of the fitted function.
        """

        return a * b * (np.exp(k1 * np.cos(th - np.pi/2.)) +
                        np.exp(k1 * np.cos(th - 3.*np.pi/2.))) / \
               (2.*np.pi * np.i0(k1))

    @classmethod
    def c2(cls, th, a, b, k2, mi2):
        """Second component of the fitted function.
        """

        return a * (1. - b) * (np.exp(k2 * np.cos(th - mi2)) +
                               np.exp(k2 * np.cos(th - (2.*np.pi - mi2)))) / \
               (2.*np.pi * np.i0(k2))

    @classmethod
    def full(cls, th, a, b, k1, k2, mi2):
        """Sum of component function.
        """

        return cls.c1(th, a, b, k1) + \
               cls.c2(th, a, b, k2, mi2)
