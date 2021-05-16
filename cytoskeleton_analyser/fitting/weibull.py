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

"""Class modelling data using Weibull distribution.
"""

from __future__ import annotations
import math
import sys
from typing import Callable, Final, Sequence

import numpy as np
from scipy.special import gamma as gammafun

from .base import Base


class Weibull(Base):
    """Implements data model using Weibull distribution.
    """

    #: Number of fitted parameters.
    nump: Final = 2    # k, lambda

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
        self.bounds = (np.array([0., 0.]),
                       np.array([+np.inf, +np.inf]))

    def sl(
            self,
            f: Callable,
    ) -> slice:
        """Select parameters depending on fitted function.

        :param f: Fitted function.
        :return: Slice selecting parameters to fit.
        """

        if str(f) == str(self.full):
            return slice(0, self.nump)
        else:  # f == FitWeibull.lambda1:
            self.p[-1] = 1.
            return slice(0, self.nump - 1)

    def _par(self) -> dict:
        """Parameters as a dictionary {'name': value}.
        """

        return {'k': self.p[0],
                'lambda': self.p[1]}

    def _mean(self) -> float:
        """Distribution mean.
        """

        return self.p[1] * gammafun(1. + 1. / self.p[0])

    def _var(self) -> float:
        """Disrtribution variance.
        """

        return self.p[1] * self.p[1] * (
                gammafun(1. + 2. / self.p[0]) -
                gammafun(1. + 1. / self.p[0])**2
        )

    def fit(
            self,
            f: Callable,
            y: np.ndarray,
    ) -> Weibull:
        """Fitting is done here.

        :param f: Function to fit.
        :param y: Data y values.
        :return: Instance of this class with the fit results set.
        """

        self.p = np.full_like(self.p0, np.nan)
        self.name = self.__class__.__name__ + '.' + f.__name__
        sl = self.sl(f)
        super()._fit(f, y, sl)
        self.par = self._par()
        self.mean = self._mean()
        self.var = self._var()
        self.aik = self._aik(y)
        self.report()
        return self

    def restore(
            self,
            f: Callable,
            y: np.ndarray,
    ) -> Weibull:
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
        self.prediction = self.predict(f, self.sl(f))
        self.set_quality_indicators(y)     # residnorm, chi2 , fano
        self.par = self._par()
        self.mean = self._mean()
        self.var = self._var()
        self.aik = self._aik(y)
        self.report()

        return self

    @classmethod
    def full(
            cls,
            x: np.ndarray,
            k: float,
            la: float,
    ) -> np.ndarray:
        """Fitted function: pdf of Weibull distribution.

        :param x: Domain.
        :param k: Fitted parameter.
        :param la: Fitted parameter.
        """

        if math.isclose(la, 0.):
            return sys.float_info.max

        return (x / la)**(k - 1.) * np.exp(-(x / la)**k) * k / la

    @classmethod
    def lambda1(
            cls,
            x: np.ndarray,
            k: float,
    ) -> np.ndarray:
        """Fitted function: pdf of Weibull distribution with lambda = 1.

        :param x: Domain.
        :param k: Fitted parameter.
        """

        return x**(k - 1.) * np.exp(-x**k) * k
