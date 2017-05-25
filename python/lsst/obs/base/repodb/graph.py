from __future__ import print_function, division, absolute_import

__all__ = ("RepoGraph",)


class RepoGraph:

    def __init__(self, units=None, datasets=None):
        if units is None:
            units = {}
        self.units = units
        if datasets is None:
            datasets = {}
        self.datasets = datasets
