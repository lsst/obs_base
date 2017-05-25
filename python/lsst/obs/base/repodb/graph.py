from __future__ import print_function, division, absolute_import

from . import base


class RepoGraph:

    def __init__(self, units=None):
        if units is None:
            units = {}
        self.units = units
