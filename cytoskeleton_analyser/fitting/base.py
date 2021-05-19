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

"""Exposes base class for subclassing in specific data models.
"""

import logging
from collections import namedtuple
from typing import Callable, Optional, Sequence

import numpy as np
from scipy.optimize import curve_fit


class Base:
    """Base class for subclassing in specific data models.
    """

    logger: logging.Logger = None
    Summary = namedtuple('Summary', 'name p chi2, mean saturation')

    def __init__(
            self,
            x: Sequence,
            p0: Sequence[float],
            x_units: Optional[str] = None
    ):
        """
        :param x: Data domain.
        :param p0: Initial parameter values.
        :param x_units: Units of x values.
        """

        #: Model name.
        self.name: Optional[str]
        #: Data x-values.
        self.x: np.ndarray = np.array(x)
        #: Units for data x-values.
        self.x_units: Optional[str] = x_units

        #: Initial values of model parameters.
        self.p0: np.ndarray = np.array(p0)

        #: Optilized values of model parameters.
        self.p: np.ndarray = np.full_like(self.p0, np.nan)

        #: Model parameters as a dictionary {'name': value}
        self.par = {}

        #: Bounds on the model parameters.
        self.bounds = (np.full_like(p0, -np.inf),
                       np.full_like(p0, np.inf))

        self.prediction = np.empty(len(x))  #: Model prediction.
        self.cc = []            #: Components of the model prediction.
        self.mean = np.nan      #: Model mean.
        self.var = np.nan       #: Model variance.
        self.saturation = None          #: Predicted saturation value.
        self.equilibrium_ind = None     #: Data index at equilibration.

        # Indicators of model quality.
        self.cov = None   #: Covariamce matrix of the fitted parameters.
        self.residnorm = np.nan         #: Norm of residuals.
        self.chi2 = np.nan              #: Chi square.
        self.fano = np.nan              #: Fano factor
        self.aik = None                 #: AIK indicator.
        self.isbest = False

    def _fit(
            self,
            f: Callable,
            y: np.ndarray,
            sl: slice
    ):
        """Fitting procedure common to all models.

        :param f: Function to fit.
        :param y: Data y values.
        :param sl: Slice selecting model-specific parameters.
        """

        self.p[sl], self.cov = \
            curve_fit(
                f=f,
                xdata=self.x,
                ydata=y,
                p0=self.p0[sl],
                maxfev=1000000,
                bounds=(self.bounds[0][sl], self.bounds[1][sl]),
            )
        self.prediction = self.predict(f, sl)
        self.set_quality_indicators(y)

    def predict(
            self,
            f: Callable,
            sl: slice,
    ) -> np.ndarray:
        """Model prediction.

        :param f: Model function.
        :param sl: Slice selecting model-specific parameters.
        """

        return f(self.x, *self.p[sl])

    def set_quality_indicators(
            self,
            y: np.ndarray,
    ) -> None:
        """Calculate quality indicators of the model after minimization.

        :param y: Fitted data array.
        """

        residuals = y - self.prediction
        inds_nonzero = \
            np.where(np.logical_not(np.isclose(self.prediction, 0.)))[0]
        self.residnorm = np.linalg.norm(residuals[inds_nonzero])
        if inds_nonzero.size:
            self.chi2 = np.sum(residuals[inds_nonzero]**2 /
                               self.prediction[inds_nonzero])
            self.fano = np.mean(residuals[inds_nonzero].var() /
                                self.prediction[inds_nonzero])

    def report(self):
        """Dump a report summarizing the model to the logger.
        """

        assert self.logger is not None

        def is_reportable(a):
            return a is not None and not np.isnan(a)

        Base.logger.info(self.name + ":")
        for key, val in self.par.items():
            if is_reportable(val):
                Base.logger.info(f"    {key} = {val}")

        Base.logger.info(f"        norm of residuals = {self.residnorm:.3f}")
        if is_reportable(self.aik):
            Base.logger.info(f"    aik: {self.aik}")
        if is_reportable(self.chi2):
            Base.logger.info(f"    chi2: {self.chi2}")
        if is_reportable(self.fano):
            Base.logger.info(f"    fano: {self.fano}")

        if is_reportable(self.mean):
            var = f", var = {self.var:.3f}" if is_reportable(self.var) else ''
            Base.logger.info(f"    mean = {self.mean:.3f}{var}")

        if is_reportable(self.saturation):
            Base.logger.info(f"    saturation at: {self.saturation}")

        if is_reportable(self.equilibrium_ind):
            Base.logger.info(f"    equilibration time to 1 std. dev.: "
                             f"{self.x[self.equilibrium_ind]:.4f} " +
                             self.x_units)
        elif is_reportable(self.saturation):
            Base.logger.info("    no equilibration is achieved "
                             "before the simulation end")

    def set_equilibration(
            self,
            y: np.ndarray,
    ) -> None:
        """Determines predicted time to equilibration if such exists.

        Is only applicable if the model achieves saturation.
        Assumes that equilibration point is reached whenever
        monotonously decreasing difference between saturation value and
        the model prediction becomes less than standard deviation.
        Sets ``equilibrium_ind``: data array index at which
        equilibration is assumed to be achieved.

        :param y: Fitted data array.
        """

        if self.saturation is None:
            return

        sd = (y - self.prediction).std()
        eqi = len(self.prediction[np.abs(self.saturation -
                                         self.prediction) > sd]) + 1
        if eqi < y.shape[0]:
            self.equilibrium_ind = eqi

    def _aik(self,
             y: np.ndarray) -> float:
        """Calculate AIK indicator of the fited model.

        :param y: Fitted data array.
        """

        return \
            2 * len(self.par) + \
            y.size * np.log(np.sum((y - self.prediction)**2))

    def summary(self) -> dict:
        """Create a dictionary sumarizing the model.
        """

        return {
            'name': self.name,
            'p': self.p.tolist(),
            'chi2': self.chi2,
            'mean': self.mean,
            'satur': self.saturation,
            'eqind': self.equilibrium_ind,
        }


def set_logger(logger):

    Base.logger = logger
