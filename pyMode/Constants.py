"""Constant definitions to be used by other code.

Alec Hammond. 2019-02-16.
"""

from collections import namedtuple
import pyMode as pm


FreqRange = namedtuple('FreqRange', ['min', 'max'])

AIR = pm.Medium(index=1)
