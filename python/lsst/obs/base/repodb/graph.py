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

    def addDataset(self, DatasetClass, **units):
        storage = {}
        for k, p in DatasetClass.properties.items():
            try:
                v = units.pop(k)
            except KeyError:
                raise ValueError(
                    "No value provided for {}.{}".format(
                        type(self).__name__, k
                    )
                )
            if not isinstance(v, p.UnitClass):
                raise TypeError(
                    "Invalid type (expected {}, got {}) for {}.{}".format(
                        p.UnitClass.__name__,
                        type(v),
                        type(self).__name__,
                        k
                    )
                )
            storage[k] = v
        if units:
            raise ValueError(
                "Unused values when constructing {}: {}".format(
                    type(self).__name__,
                    units.keys()
                )
            )
        # TODO
        raise NotImplementedError()