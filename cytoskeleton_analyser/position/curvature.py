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

"""Calculation of microtubule curvatures.
"""

from typing import Final, Union

import numpy as np


class D2:
    """Approximation of 2d curvature.

    For a sequence of points on a 2d plane.
    """
    d: Final = int(2)   #: Dimensionality.

    @staticmethod
    def c(
            dp: np.ndarray,
            ddp: np.ndarray,
    ) -> np.ndarray:
        """Calculate 2d curvature for a sequence of points.

        :param dp: 2d gradient.
        :param ddp: 2d laplacian.
        """

        return \
            np.abs(ddp[0, :]*dp[1, :] - dp[0, :]*ddp[1, :]) / \
            (dp[0, :]**2 + dp[1, :]**2)**1.5


class D3:
    """Approximation of 3d curvature.

    For a sequence of points in 3d space.
    """

    d: Final = int(3)   #: Dimensionality.

    @staticmethod
    def c(
            dp: np.ndarray,
            ddp: np.ndarray,
    ) -> np.ndarray:
        """Calculate 3d curvature for a sequence of points.

        :param dp: 2d gradient.
        :param ddp: 2d laplacian.
        """

        return \
            np.sqrt((ddp[2, :] * dp[1, :] - dp[2, :] * ddp[1, :])**2 +
                    (ddp[0, :] * dp[2, :] - dp[0, :] * ddp[2, :])**2 +
                    (ddp[1, :] * dp[0, :] - dp[1, :] * ddp[0, :])**2) / \
            (dp[0, :]**2 + dp[1, :]**2 + dp[2, :]**2)**1.5


class Curvature:
    """Generic class to encapsulate various curvature types.
    """

    def __init__(
            self,
            p: list[np.ndarray],
    ):
        """
        :param p: List of independent segmented lines. Each line is a
            sequence of points in 3d space.
        """

        self.dim2 = Curvature.ordinary(p, D2)
        self.dim3 = Curvature.ordinary(p, D3)
        self.mboc17 = Curvature.mboc17(p)

    @staticmethod
    def ordinary(
            p3d: list[np.ndarray],
            crv: Union[type[D2], type[D3]],
    ) -> tuple[np.ndarray, list[np.ndarray]]:
        """Calculate conventional curvature values at nodes
        of segmented lines.

        :param p3d: List of independent segmented lines. Each line is a
            sequence of points in 3d space.
        :param crv: Curvature class to apply: D2 or D3.
        :return: tuple of concatenated and line-by-line curvature values
            at node points.
        """

        c = [np.empty(0, dtype=float) for _ in range(len(p3d))]
        for ip, pp in enumerate(p3d):
            nnodes = pp.shape[0]
            if nnodes < 2:
                continue
            dp = np.empty([crv.d, nnodes])
            ddp = np.empty([crv.d, nnodes])
            for kk in range(crv.d):
                dp[kk, :] = np.gradient(pp[:, kk])
                ddp[kk, :] = np.gradient(dp[kk, :])
            if dp.shape[1]:
                c[ip] = np.array(crv.c(dp, ddp))

        return np.concatenate(c), c

    @staticmethod
    def mboc17(
            p3d: list[np.ndarray],
    ) -> np.ndarray:
        """Calculate 2d curvatures of segmented lines.

        In contrast to ``Curvature.ordinary``, this method implements
        algorithm adapted to filaments extracted from microscopy images
        by:
        Zhen Zhang Yukako Nishimura, and Pakorn Kanchanawonga
        ``Extracting microtubule networks from superresolution single-
        molecule localization microscopy data' MBoC, 2017.``

        :param p3d:list of segmented lines. Each line is represented by
            a sequence of points in 3d space.
        """

        c1d = []

        for pp in p3d:
            p = pp[::2, :]       # to get a 400 nm step as in the paper

            x1 = p[0::3, 0]
            x0 = p[1::3, 0]
            x2 = p[2::3, 0]
            y1 = p[0::3, 1]
            y0 = p[1::3, 1]
            y2 = p[2::3, 1]
            s = x2.size
            x10 = x1[:s] - x0[:s]
            x20 = x2[:s] - x0[:s]
            x12 = x1[:s] - x2[:s]
            y10 = y1[:s] - y0[:s]
            y20 = y2[:s] - y0[:s]
            y12 = y1[:s] - y2[:s]

            c1d.extend(2*np.abs(x10 * y20 - x20 * y10) /
                       np.sqrt((x10**2 + y10**2) *
                               (x20**2 + y20**2) *
                               (x12**2 + y12**2)))

        return np.array(c1d)
