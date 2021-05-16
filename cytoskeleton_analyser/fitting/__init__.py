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

""" Model classes and auxiliary functions for data fitting.
"""

import sys
from typing import Callable, Optional, Sequence, Union
from types import ModuleType

import numpy as np

from .base import set_logger
from .rayleigh import Rayleigh
from .weibull import Weibull
from .exponential import Exponential
from .von_mises import VonMisesDouble, VonMisesTriple
from .gamma import Gamma
from .normal import Normal
from .lognorm import Lognorm
from .const import MockConst


def is_const(model) -> bool:
    """Check if ``model`` is a mock wrapping a constant.

    :return: True if the ``model`` is a mock model representing a
        constant number.
    """

    return isinstance(model, MockConst)


def is_mock(model) -> bool:
    """Check if ``model`` is a constant number.

    :return: True if the ``model`` is a mock model representing a
        constant value.
    """

    return \
        isinstance(model, np.float64) or \
        isinstance(model, np.float32)


def fit(
        ft: Sequence,
        x: np.ndarray,
        y: np.ndarray,
) -> list:
    """Fit alternative models to the data.

    For each element of alternative fitters ``ft`` create instance of
    fitting the class and fit specific function to the data.

    :param ft: Collection of alternative fitters to apply to the data.
    :param x: Data x values.
    :param y: Data y values.
    :return: Parametrised models resulting from the fitting.
    """

    res = []
    for f in ft:
        if hasattr(f[0], '__self__'):
            res.append(f[0].__self__(x, f[1]).fit(f[0], y))
        else:
            res.append(f[0].outer(x, f[1], f[2]).fit(f[0], y))

    return res


def restore(
        f: Union[Callable, type],
        p: np.array,
        tu: str,
        x: np.ndarray,
        y: np.ndarray,
        fun: Optional[Callable] = None,
):
    """Generic function for restoring fitted models.
    """

    if hasattr(f, '__self__'):
        return f.__self__(x, p).restore(f, p, x, y)
    elif hasattr(f, 'outer'):
        return f.outer(x, p, tu).restore(f, y)
    else:
        return f(x, p).restore(fun, y)


def class_from_classname(
        module: ModuleType,
        classname: str,
) -> type:
    """Class object by name and module.
    """

    return getattr(sys.modules[module.__name__],
                   classname.split('.')[0])


def subtype_from_classname(
        c: type,
        classname: str,
) -> type:
    """Child class (if exists) od a class.

    Uses compound name ``classname`` containng child class name to
    return child class object of a given class.
    """

    if hasattr(c, 'subtype'):
        return c.subtype(classname.split('.')[1])
    else:
        return c


def method_from_classname(
        c: type,
        classname: str,
) -> Optional[Callable]:
    """Method od a class (if exists)..

    Uses compound name ``classname`` containng method name to return
    method of a given class.
    """

    a = classname.split('.')[1]
    return getattr(c, a) if hasattr(c, a) else None
