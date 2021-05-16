"""State sequences the microtubule end adopts over its lifetime.

Part of the analysis of microtubule dynamic instability.
"""

from typing import Final, Union

from .single import Single
from .double import Double
from .triple import Triple


#: State sequences available.
StateSequence: Final = Union[
    type[Single],
    type[Double],
    type[Triple],
]
