from __future__ import print_function, division, absolute_import

__all__ = ("RepoGraph",)


class RepoGraph:
    """A view of the content in a data repository as a graph of `Unit` and
    `Dataset` instances.

    The connections between `Unit`s (via `ForeignKey` descriptors) and
    `Datasets` (via `UnitProperty` descriptors) naturally forms a graph-like
    data structure.  A `RepoGraph` is the conceptual owner of a collection of
    related `Unit` and `Dataset` instances, representing a (possibly
    restricted) view of the content in a data repository.

    `RepoGraphs` should generally only be constructed by the `makeGraph`
    method of `RepoDatabase`, which ensures all `Unit`s and `Dataset`s have
    the full set of links that relate them.  After construction, new
    `Datasets` instances may be added to the graph using `addDataset`
    (typically this is done to represent a `Dataset` a `SuperTask` is expected
    to produce), but the set of `Unit` instances (as well as the types of
    `Unit`s and `Dataset`s) is considered fixed.

    A `RepoGraph` may be restricted relative to its parent `RepoDatabase` in
    three ways:

     - The set of `Unit` types may be a subset of the `Unit` types in the
       `RepoDatabase`.

     - The set of `Dataset` types may be a subset of the `Dataset` types in the
       `RepoDatabase`.

     - The set of `Unit` instances may be a subset of the `Unit` instances in
       the `RepoDatabaes`, which naturally restricts the `Dataset` instances
       in the graph to only those with `Unit`s in the restricted set.

    All of these restrictions must be consistent: all `Unit` types (instances)
    used by any `Dataset` type (instance) in the graph must be present in the
    graph, as well as any `Unit` type (instance) referred to by a `ForeignKey`
    field by a `Unit` type (instance) otherwise required.
    """

    def __init__(self, units=None, datasets=None):
        if units is None:
            units = {}
        self.units = units
        if datasets is None:
            datasets = {}
        self.datasets = datasets

    def addDataset(self, DatasetClass, **units):
        """Create a new `Dataset` instance and add it to the graph.

        If a `Dataset` instance with the given name and units already exists
        in the graph, an equivalent `Dataset` instance will be returned that
        and the graph will not be modified (the returned instance may not be
        the same in-memory object that is present in the graph, but the two
        will compare as equal and should otherwise be indistinguishable).

        Parameters
        ----------
        DatasetClass : `type` inherited from `Dataset`.
            A type object for the instance about to be constructed.

        Additional keyword parameters contain the names (keys) and `Unit`
        instances (values) that label the `Dataset`.

        Returns
        -------
        A new instance of `DatasetClass`.
        """
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
        dataset = DatasetClass(**storage)
        for k, p in DatasetClass.properties.items():
            unit = p.__get__(dataset)
            unit.datasets.setdefault(DatasetClass, set()).add(dataset)
        self.datasets[DatasetClass].add(dataset)
        return dataset
