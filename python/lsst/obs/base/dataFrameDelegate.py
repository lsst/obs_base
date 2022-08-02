# This file is part of obs_base.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (http://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""Support for reading DataFrames."""

import logging
import pandas
from typing import Any, Mapping, Optional
from lsst.daf.butler import StorageClassDelegate

log = logging.getLogger(__name__)


class DataFrameDelegate(StorageClassDelegate):

    def getComponent(self, composite: pandas.DataFrame, componentName: str) -> Any:
        """Get a component from a DataFrame.

        Parameters
        ----------
        composite : `pandas.DataFrame`
            `DataFrame` to access component.
        componentName : `str`
            Name of component to retrieve.

        Returns
        -------
        component : `object`
            The component.

        Raises
        ------
        AttributeError
            The component can not be found.
        """
        if componentName == 'columns':
            return composite.columns
        else:
            raise AttributeError(
                "Do not now how to retrieve component {} from {}".format(componentName, type(composite))
            )

    def handleParameters(self, inMemoryDataset: Any, parameters: Optional[Mapping[str, Any]] = None) -> any:
        """Modify the in-memory dataset using the supplied parameters,
        returning a possibly new object.

        Parameters
        ----------
        inMemoryDataset : `object`
            Object to modify based on the parameters.
        parameters : `dict`, optional
            Parameters to apply. Values are specific to the parameter.
            Supported parameters are defined in the associated
            `StorageClass`.  If no relevant parameters are specified the
            inMemoryDataset will be return unchanged.

        Returns
        -------
        inMemoryDataset : `object`
            Updated form of supplied in-memory dataset, after parameters
            have been used.
        """
        if 'columns' in parameters:
            allColumns = list(inMemoryDataset.columns)
            allColumns.extend(list(inMemoryDataset.index.names))

            for column in parameters['columns']:
                if not isinstance(column, str):
                    raise NotImplementedError(
                        "InMemoryDataset of a DataFrame only supports string column names."
                    )
                if column not in allColumns:
                    raise ValueError(f"Unrecognized column name {column}.")

            # Exclude index columns from the subset
            readColumns = [
                name for name in parameters['columns']
                if name not in inMemoryDataset.index.names
            ]

            return inMemoryDataset[readColumns]
        else:
            return inMemoryDataset
